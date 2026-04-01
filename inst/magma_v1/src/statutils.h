/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef STATUTILS_H
#define STATUTILS_H

#include "_global.h"
#include "utils.h"
#include "buffer.h"
#include "mathutils.h"
#include "exceptions.h"
#include "output.h"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <Eigen/Eigen>

#include <limits>

#define CHISQ_ASYM_FUDGE 1.025

namespace StatType {///weak ordering: ill-defined / p-values in front, 'real' statistics from LogPvalue / Z onwards
  enum Value {None, Unknown, Pvalue, PermPvalue, LogPvalue, Z, ChiSquare, ChisqSumDf1, ChisqSumDf2, T, F, Beta, Gamma, Pareto, PowerSum};
}


class StatConverter;
class StatUtils {
public:
  static StatConverter* get_converter(StatType::Value from = StatType::None, StatType::Value to = StatType::None, bool do_throw=true);

  template<typename T>
  static void add_noise(T* data, long len, double sd, bool normalize=false, bool compute_scale=false) {
    double scale = 1;
    if (normalize && compute_scale) {
      scale = MathUtils::get_sd(data, len);
      if (scale <= 0) throw StatException("normalization", "the variance for the input variable is negative or 0");
      sd *= scale;
    }

    generate_normal(data, len, 0, sd, true);
    if (normalize) MathUtils::normalize(data, len, scale);
  }
  
  template<typename T>
  static void generate_normal(T* data, long len, bool add=false) {generate_normal(data, len, 0, 1, add);}

  template<typename T>
  static void generate_normal(T* data, long len, double mean, double sd, bool add=false) {
    boost::normal_distribution<T> distn(mean,sd);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<T> > norm_gen(_RNG, distn);

    if (add) for (T* end = data+len; data < end; data++) *data += norm_gen();  
    else for (T* end = data+len; data < end; data++) *data = norm_gen();  
    _RNG = norm_gen.engine();    
  }
  
  template<typename T>
  static void generate_chisq(T* data, long len, double df, bool add=false) {
    boost::random::chi_squared_distribution<T> distcs(df);
    boost::variate_generator<boost::mt19937, boost::random::chi_squared_distribution<T> > chisq_gen(_RNG, distcs);

    if (add) for (T* end = data+len; data < end; data++) *data += chisq_gen();  
    else for (T* end = data+len; data < end; data++) *data = chisq_gen();  
    _RNG = chisq_gen.engine();    
  }

};

class StatConverter {
  StatType::Value from_type;
  StatType::Value to_type;

public:
  StatConverter(StatType::Value from = StatType::None, StatType::Value to = StatType::None) : from_type(from), to_type(to) {}
  virtual ~StatConverter(){}

  virtual void set_param(const string& name, double value) {throw StatException("stat conversion", "trying to set unknown parameter for StatConverter"); UNUSED(name); UNUSED(value);}

  virtual double convert(double value) = 0;
  virtual void convert(double* buff, long len) = 0;
  virtual void convert_obs(double* buff, long len, double& obs) {
    obs = convert(obs);
    convert(buff, len);
  }    
 
  virtual void convert_obs_buff(Buffer<double>& obs_buff, Buffer<double>& perm_buff, long no_obs, long no_perm) {
    for (int i = 0; i < no_obs; i++) convert_obs(perm_buff[i], no_perm, obs_buff(i));
  }
  
  virtual StatConverter* clone() = 0;
};

class FromP : public StatConverter {
protected:
  double min_value;
  double max_value;

public:
  FromP(StatType::Value to) : StatConverter(StatType::Pvalue, to) {
    min_value = numeric_limits<double>::min();
    max_value = 1-numeric_limits<double>::epsilon();
  }
  virtual ~FromP() {}
};

class NullConverter : public StatConverter {
public:
  virtual ~NullConverter() {}
  virtual void convert(double* buff, long len) {UNUSED(buff); UNUSED(len);}
  virtual double convert(double value) {return value;}
  
  virtual NullConverter* clone() {return new NullConverter();}
};

class FudgePval : public FromP {
  double fudge;
  double fudge_min;
  StatConverter* internal;
  bool owner;  

public:
  FudgePval(double f=1, StatConverter* pval=0, bool owner=false) : FromP(StatType::Pvalue), fudge(f), fudge_min(0.5), internal(pval), owner(owner) {}
  virtual ~FudgePval() {if (owner) delete internal;}
  
  void set_internal(StatConverter* pval=0, bool owned=false) {
    if (internal && owner) delete internal; 
    internal = pval; owner = owned;
  }
  virtual void set_param(const string& name, double value) {
    if (name == "fudge") {
      fudge = value > 0 ? value : 1;
    } else if (name == "fudge_min") {
      fudge_min = value > 0 ? value : 0;
    } else {
      if (internal) internal->set_param(name, value); 
      else StatConverter::set_param(name, value);
    }
  }
  virtual void convert(double* buff, long len);
  virtual double convert(double value);

  virtual FudgePval* clone() {
    FudgePval* out = new FudgePval(fudge, internal->clone(), true);
    out->set_param("fudge_min", fudge_min);
    return out;
  }
};

class CompoundConverter : public StatConverter {
  StatConverter* conv_in;
  StatConverter* conv_out;  
  bool owner;
  
public:
  CompoundConverter(StatConverter* in, StatConverter* out, bool owner) : conv_in(in), conv_out(out), owner(owner) {}
  CompoundConverter(StatType::Value in, StatType::Value out) : owner(true) {
    conv_in = StatUtils::get_converter(in, StatType::Pvalue);     
    conv_out = StatUtils::get_converter(StatType::Pvalue, out);         
  }
  virtual ~CompoundConverter() {
    if (owner) {delete conv_in; delete conv_out;}
  }
  void set_param_in(const string& name, double value) {conv_in->set_param(name, value);}
  void set_param_out(const string& name, double value) {conv_out->set_param(name, value);}  

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual CompoundConverter* clone() {return new CompoundConverter(conv_in->clone(), conv_out->clone(), true);}  
};

class ConvertPvalToLogPval : public FromP {
public:
  ConvertPvalToLogPval() : FromP(StatType::LogPvalue) {}
  virtual ~ConvertPvalToLogPval() {}
  
  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertPvalToLogPval* clone() {return new ConvertPvalToLogPval();}
};

class ConvertLogPvalToPval : public StatConverter {
public:
  ConvertLogPvalToPval() {}
  virtual ~ConvertLogPvalToPval(){}
  
  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertLogPvalToPval* clone() {return new ConvertLogPvalToPval();}
};


class ConvertNormToPval : public StatConverter {
  double sd;
public:
  ConvertNormToPval(double sd=1) : StatConverter(StatType::Z, StatType::Pvalue), sd(sd) {}
  virtual ~ConvertNormToPval(){}
 
  double get_sd() {return sd;}
 
  virtual void set_param(const string& name, double value) {
    if (name == "sd") sd = value > 0 ? value : 1;
    else if (name == "var") sd = value > 0 ? sqrt(value) : 1;
    else StatConverter::set_param(name, value);
  } 
  
  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertNormToPval* clone() {return new ConvertNormToPval(sd);}
};

class ConvertPvalToNorm : public FromP {
public:
  ConvertPvalToNorm() : FromP(StatType::Z) {}
  virtual ~ConvertPvalToNorm(){}
  
  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertPvalToNorm* clone() {return new ConvertPvalToNorm();}
};


class ConvertChisqToPval : public StatConverter {
  double min_value;

public:
  ConvertChisqToPval() : StatConverter(StatType::ChiSquare, StatType::Pvalue) {min_value = numeric_limits<double>::min();}
  virtual ~ConvertChisqToPval(){}

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertChisqToPval* clone() {return new ConvertChisqToPval();}
};

class ConvertPvalToChisq : public FromP {
public:
  ConvertPvalToChisq() : FromP(StatType::ChiSquare) {}
  virtual ~ConvertPvalToChisq(){}
  
  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertPvalToChisq* clone() {return new ConvertPvalToChisq();}
};


class ConvertChisqSumDf1ToPval : public StatConverter { 
  ConvertChisqToPval uni_converter; double min_value; 

protected:
  BaseBuffer<double> lambda; int ncomp;
  
  virtual double convert_internal(double value) = 0;
  virtual void process_lambda() {}
  
public:
  ConvertChisqSumDf1ToPval() : StatConverter(StatType::ChisqSumDf1, StatType::Pvalue), min_value(numeric_limits<double>::min()), ncomp(0) {}
  virtual ~ConvertChisqSumDf1ToPval(){}

  double get_ncomp() {return ncomp;}
  BaseBuffer<double>& get_lambda() {return lambda;}

  template<typename T>
  int set_corrs(Buffer<T>& cor, double cutoff=1e-4) {
    int nvar = cor.ncol(); lambda.resize(nvar);
    SelfAdjointEigenSolver<Matrix<T,Dynamic,Dynamic> > eig(Map<Matrix<T,Dynamic,Dynamic> >(cor[0], nvar, nvar));
    const T* eval = eig.eigenvalues().data(); double total = 0; ncomp = 0;
    for (const T* curr = eval; curr < eval + nvar; curr++) {if (*curr > 0) total += *curr;}
    cutoff *= total / nvar;        
    
    for (const T* curr = eval; curr < eval + nvar; curr++) {if (*curr > cutoff) lambda[ncomp++] = *curr;}
    if (ncomp == 0) throw StatException("pval computation", "could not extract principal components");       
    process_lambda();

    return ncomp;
  }
  
  void set_lambda(BaseBuffer<double>& lam, int use=0);

  virtual void convert(double* buff, long len) {for (long i = 0; i < len; i++) buff[i] = convert(buff[i]);}
  virtual double convert(double value);
};

class ConvertChisqSumDf1ToPvalBrown : public ConvertChisqSumDf1ToPval {
  FudgePval* fudge;
  double df, inv_scale;

  double convert_internal(double value); 
  void process_lambda();
  
public:
  ConvertChisqSumDf1ToPvalBrown(double fudge_factor=CHISQ_ASYM_FUDGE) : df(1), inv_scale(1) {fudge = new FudgePval(fudge_factor);}  
  virtual ~ConvertChisqSumDf1ToPvalBrown() {delete fudge;}

  double get_df() {return df;}
  double get_scale() {return 1/inv_scale;}

  virtual ConvertChisqSumDf1ToPvalBrown* clone() {
    ConvertChisqSumDf1ToPvalBrown* out = new ConvertChisqSumDf1ToPvalBrown();
    out->fudge = fudge->clone();
    out->set_lambda(lambda, ncomp);
    return out;
  }
};

class ConvertChisqSumDf1ToPvalEmpirical : public ConvertChisqSumDf1ToPval {
protected:
  int min_iter, max_iter, max_comp, adap_thresh;  
  int block_iter, max_blocks;
  Persistent<BaseBuffer<double> >* perms;

  BaseBuffer<double> perm_stats;
  Buffer<double> sub_stats;
  int avail_iter, avail_sub;

  void init_perms();
  void init_subblocks();
  int primary_stats(int total);
  void secondary_stats(double* write, int offset);
  
  double convert_internal(double value); 
  void process_lambda() {
    if (ncomp > max_comp) throw StatException("pval computation", "number of principal components exceeds maximum");       
    avail_iter = 0; avail_sub = 0;
  }
  
public:
  ConvertChisqSumDf1ToPvalEmpirical(int min_iter, int max_iter, int max_comp, int adap_thresh, Persistent<BaseBuffer<double> >* input = 0) 
    : min_iter(min_iter), max_iter(max_iter), max_comp(max_comp), adap_thresh(adap_thresh), block_iter(1e6), max_blocks(10), perms(0), avail_iter(0), avail_sub(0) {
    if (input) perms = input->copy(); 
    else init_perms();
  }
  virtual ~ConvertChisqSumDf1ToPvalEmpirical() {delete perms;}

  double minimum_pval() {return 0.5/max_iter;}

  virtual ConvertChisqSumDf1ToPvalEmpirical* clone() {
    ConvertChisqSumDf1ToPvalEmpirical* out = new ConvertChisqSumDf1ToPvalEmpirical(min_iter, max_iter, max_comp, adap_thresh, perms);
    out->set_lambda(lambda, ncomp);
    out->perm_stats.assign(perm_stats); out->sub_stats.assign(sub_stats);
    out->avail_iter = avail_iter; out->avail_sub = avail_sub;
    return out;
  }
};  


class ConvertChisqSumDf1ToPvalImhof : public ConvertChisqSumDf1ToPval {
  double epsilon, pval_cutoff; pair<double,double> bound_range;
  ConvertChisqSumDf1ToPvalBrown bounds;
  ConvertChisqSumDf1ToPvalEmpirical* backstop;  
  BaseBuffer<double> lambda_sq; 
  
  class ImhofFunc {
    double obs, *lambda, *lambda_sq; int ncomp;
  public:
    ImhofFunc(double obs, double* lambda, double* lambda_sq, int ncomp) : obs(obs), lambda(lambda), lambda_sq(lambda_sq), ncomp(ncomp) {}

    double operator() (double u) const;  
  };

  double convert_internal(double value); 
  void process_lambda();
  
public:
  ConvertChisqSumDf1ToPvalImhof(double epsilon_factor=1e5, double cutoff=1e-50) : pval_cutoff(cutoff), bound_range(-5,2.5), backstop(0) {
    epsilon = numeric_limits<double>::epsilon() * epsilon_factor;
    if (epsilon < numeric_limits<double>::epsilon()) epsilon = numeric_limits<double>::epsilon();
    if (epsilon > sqrt(numeric_limits<double>::epsilon())) epsilon = sqrt(numeric_limits<double>::epsilon()); 
  }
  virtual ~ConvertChisqSumDf1ToPvalImhof() {delete backstop;}

  void set_backstop(int min_iter, int max_iter, int max_comp, int adap_thresh) {
    delete backstop;
    backstop = new ConvertChisqSumDf1ToPvalEmpirical(min_iter, max_iter, max_comp, adap_thresh);
  }

  virtual ConvertChisqSumDf1ToPvalImhof* clone() {
    ConvertChisqSumDf1ToPvalImhof* out = new ConvertChisqSumDf1ToPvalImhof();
    out->set_lambda(lambda, ncomp);
    out->epsilon = epsilon; 
    if (backstop) out->backstop = backstop->clone();
    return out;
  }
};




/// Brown approach
class ConvertChisqSumDf2ToPval : public StatConverter {
  double df;
  double scale;
  double min_value;
  
public:
  ConvertChisqSumDf2ToPval(double df=1, double scale=1) : StatConverter(StatType::ChisqSumDf2, StatType::Pvalue), df(df), scale(scale) {min_value = numeric_limits<double>::min();}
  virtual ~ConvertChisqSumDf2ToPval(){}

  double get_df() {return df;}
  double get_scale() {return 1/scale;}

  template<typename T>
  static void set_brown(ConvertChisqSumDf2ToPval& converter, Buffer<T>& cor) {
    int dim = cor.ncol(); T exp = 2*dim, var = 0;
    for (int i = 0; i < dim; i++) {
      for (int j = i+1 ; j < dim; j++) {
        T& rho = cor(j,i);
        var += rho >= 0 ? rho * (3.25 + 0.75*rho) : rho * (3.27 + 0.71*rho);
      }
    }
    converter.set_moments(exp, 4*dim + 2*var);
  }
  
  template<typename T>
  static void set_perm(ConvertChisqSumDf2ToPval& converter, T* perms, int nperm) {
    double mean, var; MathUtils::get_stats(perms, nperm, mean, var, false);
    converter.set_moments(mean, var);        
  }
  
  void set_moments(double mean, double var) {
    if (var <= 0) throw StatException("pval computation", "the variance for the chi-square statistic is negative or 0");
    double df = 2*mean*mean/var, scale = var / (2*mean);
    set_param("df", df); set_param("scale", 1/scale);
  }

  virtual void set_param(const string& name, double value) {
    if (name == "df") df = value > 0 ? value : 1;
    else if (name == "scale") scale = value > 0 ? value : 1;
    else StatConverter::set_param(name, value);
  }

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertChisqSumDf2ToPval* clone() {return new ConvertChisqSumDf2ToPval(df, scale);}
};


class ConvertTToPval : public StatConverter {
  double df;

public:
  ConvertTToPval() : StatConverter(StatType::T, StatType::Pvalue), df(1) {}
  ConvertTToPval(double degree) : StatConverter(StatType::T, StatType::Pvalue), df(degree) {}
  virtual ~ConvertTToPval(){}

  virtual void set_param(const string& name, double value) {
    if (name == "df") {
      df = value > 0 ? value : 1;
    } else StatConverter::set_param(name, value);
  }

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertTToPval* clone() {return new ConvertTToPval(df);}  
};

class ConvertPvalToT : public FromP {
  double df;

public:
  ConvertPvalToT() : FromP(StatType::T), df(1) {}
  ConvertPvalToT(double degree) : FromP(StatType::T), df(degree) {}
  virtual ~ConvertPvalToT(){}

  virtual void set_param(const string& name, double value) {
    if (name == "df") {
      df = value > 0 ? value : 1;
    } else StatConverter::set_param(name, value);
  }

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertPvalToT* clone() {return new ConvertPvalToT(df);}
};


class ConvertFToPval : public StatConverter {
  double df1;
  double df2;
  double min_value;

public:
  ConvertFToPval() : StatConverter(StatType::F, StatType::Pvalue), df1(1), df2(10000) {min_value = numeric_limits<double>::min();}
  ConvertFToPval(double df1, double df2) : StatConverter(StatType::F, StatType::Pvalue), df1(df1), df2(df2) {min_value = numeric_limits<double>::min();}
  virtual ~ConvertFToPval(){}

  virtual void set_param(const string& name, double value) {
    if (name == "df1") {
      df1 = value > 0 ? value : 1;
    } else if (name == "df2") {
      df2 = value > 0 ? value : 1;
    } else StatConverter::set_param(name, value);
  }

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertFToPval* clone() {return new ConvertFToPval(df1, df2);}
};

class ConvertBetaToPval : public StatConverter {
protected:
  double shape1;
  double shape2;
  double min_value;
  double max_value;

public:
  ConvertBetaToPval(double shape1=1, double shape2=1) : StatConverter(StatType::Beta, StatType::Pvalue), shape1(shape1), shape2(shape2) {
    min_value = numeric_limits<double>::min();
    max_value = 1-numeric_limits<double>::epsilon();    
  }
  virtual ~ConvertBetaToPval(){}

  virtual void set_param(const string& name, double value) {
    if (name == "shape1") {
      shape1 = value > 0 ? value : 1;
    } else if (name == "shape2") {
      shape2 = value > 0 ? value : 1;
    } else StatConverter::set_param(name, value);
  }
  virtual void set_minp(int nvar) {shape1 = 1; shape2 = nvar;}

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertBetaToPval* clone() {return new ConvertBetaToPval(shape1, shape2);}
};

class ConvertFittedBetaToPval : public ConvertBetaToPval {
  ConvertPvalToLogPval to_logp;
  Buffer<double> lin_buffer;

  bool use_linear;
  double intercept;
  double slope;

public:
  ConvertFittedBetaToPval() : ConvertBetaToPval(), use_linear(false), intercept(0), slope(0) {}
  virtual ~ConvertFittedBetaToPval(){}

  virtual void set_param(const string& name, double value) {ConvertBetaToPval::set_param(name, value); use_linear = false;}
  virtual void set_minp(int nvar) {ConvertBetaToPval::set_minp(nvar); use_linear = false;}
  virtual double fit_param(double* pval, int len);

  virtual void convert(double* buff, long len);
  virtual double convert(double value);

  virtual ConvertFittedBetaToPval* clone() {
    ConvertFittedBetaToPval* fbeta = new ConvertFittedBetaToPval();
    fbeta->shape1 = shape1; fbeta->shape2 = shape2;    
    fbeta->intercept = intercept; fbeta->slope = slope;
    fbeta->use_linear = use_linear;
    return fbeta;    
  }  
};

class ConvertGammaToPval : public StatConverter {
  double shape;
  double scale;
  double min_value;

public:
  ConvertGammaToPval(double shape=1, double scale=1) : StatConverter(StatType::Gamma, StatType::Pvalue), shape(shape), scale(scale) {min_value = numeric_limits<double>::min();}
  virtual ~ConvertGammaToPval(){}

  virtual void set_param(const string& name, double value) {
    if (name == "shape") {
      shape = value > 0 ? value : 1;
    } else if (name == "scale") {
      scale = value > 0 ? value : 1;
    } else StatConverter::set_param(name, value);
  }
  virtual void fit_param(double* logp, int len);

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertGammaToPval* clone() {return new ConvertGammaToPval(shape, scale);}
};

class ConvertParetoToPval : public StatConverter {
  double shape;
  double scale;
  double cutoff;
  double proportion;
  
  double (ConvertParetoToPval::*conv)(double&);
    
  double max_value;
  double max_default;
  
  void set_pfunc() {
    max_value = 0.999*((shape < 0) ? cutoff - scale / shape : max_default);
    conv = (shape != 0) ? &ConvertParetoToPval::conv_gpd : &ConvertParetoToPval::conv_exp;
  }
  
  double conv_exp(double& value) {return value > cutoff ? exp(-(value-cutoff)/scale) : 1;}
  double conv_gpd(double& value) {return value > cutoff ? pow(1+(value-cutoff)*shape/scale, -1/shape) : 1;}

public:
  ConvertParetoToPval(double shape=0, double scale=1, double cutoff=0, double prop=1) : StatConverter(StatType::Pareto, StatType::Pvalue), shape(shape), scale(scale), cutoff(cutoff) {
    proportion = MathUtils::clamp01(prop);
    max_default = numeric_limits<double>::max();
    set_pfunc();
  }
  virtual ~ConvertParetoToPval(){}

  virtual void set_param(const string& name, double value) {
    if (name == "prop") proportion = MathUtils::clamp01(value);
    else if (name == "cutoff") cutoff = value;
    else if (name == "shape") shape = value;
    else if (name == "scale") {
      scale = value > 0 ? value : 1;
    } else StatConverter::set_param(name, value);
  }

  virtual double fit_param(double* perm, int len, double tail=1);

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertParetoToPval* clone() {
    ConvertParetoToPval* out = new ConvertParetoToPval(shape, scale, cutoff, proportion);
    out->max_value = max_value;
    return out;  
  }
};

class ConvertPermutationToPval : public StatConverter {
  ConvertGammaToPval* fgamma;
  ConvertParetoToPval* fgpd;
  ConvertLogPvalToPval noop;
  Buffer<double> fit_buffer;

  double tail_prop;
  int tail_count;  
  double threshold;
  int min_perm;
  
  double thresh_max;
  bool use_noop;
    
public:
  ConvertPermutationToPval(double prop=0.05, int min_count=1000, int min_perm=0) : StatConverter(StatType::Unknown, StatType::Pvalue), tail_count(min_count), min_perm(min_perm), thresh_max(numeric_limits<double>::max()), use_noop(false) {
    fgamma = new ConvertGammaToPval(); fgpd = new ConvertParetoToPval();    
    tail_prop = MathUtils::clamp01(prop);
    threshold = thresh_max;
  }
  virtual ~ConvertPermutationToPval() {delete fgamma; delete fgpd;}

  void toggle_noop(bool to=true) {use_noop = to;}
  virtual void set_param(const string& name, double value) {
    if (name == "prop") tail_prop = MathUtils::clamp01(value);
    else if (name == "min_count") tail_count = value;
    else if (name == "min_perm") min_perm = value;
    else if (name == "noop") toggle_noop(value);
    else StatConverter::set_param(name, value);
  }

  virtual double fit_param(double* perm, int len, bool keep_ordered=false);

  virtual void convert(double* buff, long len);
  virtual double convert(double value);
  
  virtual ConvertPermutationToPval* clone() {
    ConvertPermutationToPval* out = new ConvertPermutationToPval(tail_prop, tail_count, min_perm);
    delete out->fgamma; out->fgamma = fgamma->clone(); 
    delete out->fgpd; out->fgpd = fgpd->clone(); 
    out->threshold = threshold; out->use_noop = use_noop;
    return out;
  }
}; 

class BoundedMinP : public StatConverter {
  Buffer<double> minp_buffer;
  int no_tests;
  int bounded_max;
  bool log_mode;
  
  ConvertPermutationToPval* internal;
  bool owner;  

  double constrain(double pval, double minp);

public:
  BoundedMinP(ConvertPermutationToPval* to_pval=0, bool log_mode=false, bool owner=false) : no_tests(1), bounded_max(100), log_mode(log_mode), internal(to_pval), owner(owner) {}
  virtual ~BoundedMinP() {if (owner) delete internal;}

  ConvertPermutationToPval* get_internal() {return internal;}
  void set_internal() {set_internal(0,false,false);}
  void set_internal(ConvertPermutationToPval* pval, bool log, bool owned=false) {
    if (internal && owner) delete internal; 
    internal = pval; log_mode = log; owner = owned;
  }
  virtual void set_param(const string& name, double value) {
    if (name == "ntest") no_tests = value > 1 ? value : 1;
    else if (name == "max") bounded_max = value > 1 ? value : 1;
    else {
      if (internal) internal->set_param(name, value); 
      else StatConverter::set_param(name, value);
    }
  }

  virtual void convert(double* buff, long len);
  virtual double convert(double value);

  virtual BoundedMinP* clone() {
    BoundedMinP* out = new BoundedMinP(internal->clone(), log_mode, true);
    out->set_param("ntest", no_tests); out->set_param("max", bounded_max);    
    return out;
  }  
};

class ConvertFromUnknown : public StatConverter {
  StatConverter* internal;
  bool owner;
  BuffSorter<double> sorter;
  
  double rank(double* buff, long len, double obs=0, bool comp_obs=false);
  
public:
  ConvertFromUnknown(StatConverter* to, bool owner) : internal(to), owner(owner) {}
  ConvertFromUnknown(StatType::Value to) : StatConverter(StatType::Unknown, to), internal(0), owner(true) {
    internal = StatUtils::get_converter(StatType::Pvalue, to);  
  }
  virtual ~ConvertFromUnknown() {
    if (owner) delete internal;
  }
  virtual void convert(double* buff, long len);
  virtual double convert(double value) {throw StatException("stat conversion", "unable to convert statistic without permutations"); UNUSED(value);}
  virtual void convert_obs(double* buff, long len, double& obs);
  
  virtual ConvertFromUnknown* clone() {return new ConvertFromUnknown(internal->clone(), true);}
};

class ConvertFromPermPval : public StatConverter {
  StatConverter* internal;
  bool owner;
  double min_value;
  double max_value;  
  
  void set_min(double min) {min_value = (min <= 0 || min >= 1) ? 0.5/1000 : min; max_value = 1 - min_value;}
public:
  ConvertFromPermPval(StatConverter* to, bool owner) : internal(to), owner(owner) {}
  ConvertFromPermPval(StatType::Value to) : StatConverter(StatType::Unknown, to), internal(0) {
    set_min(0.5/1000);
    internal = StatUtils::get_converter(StatType::Pvalue, to);  
  }
  virtual ~ConvertFromPermPval() {delete internal;}
  
  virtual void set_param(const string& name, double value) {
    if (name == "N") set_min(0.5/value);
    else  if (name == "min") set_min(value);
    else StatConverter::set_param(name, value);
  }
  virtual void convert(double* buff, long len) {
    for (long i = 0; i < len; i++) {buff[i] = buff[i] > min_value ? (buff[i] < max_value ? buff[i] : max_value) : min_value;}
    internal->convert(buff, len);
  }
  
  virtual double convert(double value) {
    return internal->convert(value > min_value ? (value < max_value ? value : max_value) : min_value);
  }
  
  virtual ConvertFromPermPval* clone() {
    ConvertFromPermPval* out = new ConvertFromPermPval(internal->clone(), true); 
    out->set_param("min", min_value);
    return out;
  }
};

#endif /**STATUTILS_H*/

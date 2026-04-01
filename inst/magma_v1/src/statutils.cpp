/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "statutils.h"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>


#define FITTEDBETA_THRESH 1.8
#define PARETO_MAX_TAIL 0.1

StatConverter* StatUtils::get_converter(StatType::Value from, StatType::Value to, bool do_throw) {
  if (from == to || to == StatType::Unknown) return new NullConverter();
  if (from == StatType::Unknown) return new ConvertFromUnknown(to);
  if (from == StatType::PermPvalue) return new ConvertFromPermPval(to);  
  if (from == StatType::Z && to == StatType::Pvalue) return new ConvertNormToPval();
  if (from == StatType::T && to == StatType::Pvalue) return new ConvertTToPval();
  if (from == StatType::F && to == StatType::Pvalue) return new ConvertFToPval();
  if (from == StatType::ChiSquare && to == StatType::Pvalue) return new ConvertChisqToPval();  
  if (from == StatType::Beta && to == StatType::Pvalue) return new ConvertBetaToPval();    
  if (from == StatType::Gamma && to == StatType::Pvalue) return new ConvertGammaToPval();      
  if (from == StatType::LogPvalue && to == StatType::Pvalue) return new ConvertLogPvalToPval();  
  if (from == StatType::Pvalue) {
    if (to == StatType::Z) return new ConvertPvalToNorm();
    if (to == StatType::T) return new ConvertPvalToT();
    if (to == StatType::ChiSquare) return new ConvertPvalToChisq();
    if (to == StatType::LogPvalue) return new ConvertPvalToLogPval();
  }
  
  if (from > StatType::PermPvalue && from > StatType::PermPvalue) {
    StatConverter* in = get_converter(from, StatType::Pvalue);
    StatConverter* out = get_converter(StatType::Pvalue, to);
    if (in && out) return new CompoundConverter(in, out, true);
    else {delete in; delete out;}
  }
  
  if (do_throw) throw StatException("stat computation", "requested conversion between statistics is not defined");
  return 0;
}


using boost::math::normal;
void ConvertNormToPval::convert(double* buff, long len) { try {
  static normal distn(0,1);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the Z-statistic is not finite");
    buff[i] = cdf(distn, -buff[i]/sd);
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from normal", e.what());}}

using boost::math::normal;
double ConvertNormToPval::convert(double value) { try {
  static normal distn(0,1);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the Z-statistic is not finite");
  return cdf(distn, -value/sd);
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from normal", e.what());}}




using boost::math::normal;
void ConvertPvalToNorm::convert(double* buff, long len) { try {
  static normal distn(0,1);
  for (long i = 0; i < len; i++) buff[i] = quantile(complement(distn, buff[i] >= max_value ? max_value : (buff[i] <= min_value ? min_value : buff[i])));
} catch (const exception& e) {throw StatException("stat computation", "an error occurred when trying to convert p-value to normal", e.what());}}

using boost::math::normal;
double ConvertPvalToNorm::convert(double value) { try {
  static normal distn(0,1);
  return quantile(complement(distn, value >= max_value ? max_value : (value <= min_value ? min_value : value)));
} catch (const exception& e) {throw StatException("stat computation", "an error occurred when trying to convert p-value to normal", e.what());}}




using boost::math::chi_squared;
void ConvertChisqToPval::convert(double* buff, long len) { try {
  static chi_squared distcs(1);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the chi-squared statistic is not finite");
    buff[i] = cdf(complement(distcs, (buff[i] <= min_value ? min_value : buff[i])));
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from chi-squared", e.what());}}

using boost::math::chi_squared;
double ConvertChisqToPval::convert(double value) { try {
  static chi_squared distcs(1);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the chi-squared statistic is not finite");
  return cdf(complement(distcs, (value <= min_value ? min_value : value)));
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from chi-squared", e.what());}}



using boost::math::chi_squared;
void ConvertPvalToChisq::convert(double* buff, long len) { try {
  static chi_squared distcs(1);
  for (long i = 0; i < len; i++) buff[i] = quantile(complement(distcs, buff[i] >= max_value ? max_value : (buff[i] <= min_value ? min_value : buff[i])));
} catch (const exception& e) {throw StatException("stat computation", "an error occurred when trying to convert p-value to chi-square", e.what());}}


using boost::math::chi_squared;
double ConvertPvalToChisq::convert(double value) { try {
  static chi_squared distcs(1);
  return quantile(complement(distcs, value >= max_value ? max_value : (value <= min_value ? min_value : value)));
} catch (const exception& e) {throw StatException("stat computation", "an error occurred when trying to convert p-value to chi-square", e.what());}}



void ConvertChisqSumDf1ToPval::set_lambda(BaseBuffer<double>& lam, int use) {
  lambda.assign(lam); ncomp = use > 0 ? use : lam.size();
  process_lambda();  
}

double ConvertChisqSumDf1ToPval::convert(double value) { try {
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the chi-square sum statistic is not finite");
  if (value < min_value) value = min_value;
  if (ncomp > 1) return convert_internal(value);
  else if (ncomp == 1) return uni_converter.convert(value/lambda[0]);
  else throw StatException("pval computation", "could not extract principal components");
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from chi-square sum", e.what());}}


void ConvertChisqSumDf1ToPvalBrown::process_lambda() {
  if (ncomp > 1) {
    double mean = 0, var = 0;
    for (int i = 0; i < ncomp; i++) {mean += lambda[i]; var += lambda[i]*lambda[i];}
    var *= 2;

    df = 2*mean*mean / var;
    inv_scale = 2*mean / var;
  } else {df = 1; inv_scale = 1;}
}

using boost::math::chi_squared;
double ConvertChisqSumDf1ToPvalBrown::convert_internal(double value) {return fudge->convert(cdf(complement(chi_squared(df), inv_scale*value)));}


void ConvertChisqSumDf1ToPvalEmpirical::init_perms() {
  if (perms) delete perms;
  int required = block_iter+2*max_comp;
  perms = new Persistent<BaseBuffer<double> >(new BaseBuffer<double>(required));
  StatUtils::generate_chisq(perms->get_content()->data(), required, 1);
  perm_stats.resize(2*block_iter); 
}

void ConvertChisqSumDf1ToPvalEmpirical::init_subblocks() {
  int margin = max_blocks * int(ceil(double(max_iter)/block_iter));
  if (margin > block_iter) throw StatException("pval computation", "maximum number of permutations exceeds subblock margin");       
  if (sub_stats.empty()) sub_stats.set_size(block_iter+margin, max_blocks);  

  avail_sub = min(max_blocks, ncomp);
  double* read = perms->get_content()->data();
  for (int i = 0; i < avail_sub; i++) MathUtils::vec_scale(read+2*i, sub_stats[i], lambda[i], block_iter);
  for (int i = avail_sub; i < ncomp; i++) MathUtils::vec_sum(read+2*i, sub_stats[i%avail_sub], lambda[i], block_iter);
  for (int i = 0; i < avail_sub; i++) memcpy(sub_stats[i]+block_iter, sub_stats[i], sizeof(double)*margin);
}

int ConvertChisqSumDf1ToPvalEmpirical::primary_stats(int total) {
  if (total > block_iter) total = block_iter;
  if (avail_iter < total) {
    double *read = perms->get_content()->data()+avail_iter, *write = perm_stats.data()+avail_iter;
    int length = total - avail_iter; 

    MathUtils::vec_scale(read, write, lambda[0], length);
    if (ncomp > 1) {for (int i = 1; i < ncomp; i++) MathUtils::vec_sum(read+i, write, lambda[i], length);}
    avail_iter = total;
  }
  return total;
}

void ConvertChisqSumDf1ToPvalEmpirical::secondary_stats(double* write, int offset) {
  if (avail_sub > 1) {
    MathUtils::vec_sum(sub_stats[0], sub_stats[1]+offset, write, block_iter);
    for (int i = 2; i < avail_sub; i++) MathUtils::vec_sum(sub_stats[i]+i*offset, write, block_iter);
  } else if (avail_sub == 1) memcpy(write, sub_stats[0], sizeof(double)*block_iter);
}


double ConvertChisqSumDf1ToPvalEmpirical::convert_internal(double value) {
  int curr_iter = 0, curr_count = 0; double* stats = perm_stats.data();
  for (int target = min_iter; curr_count < adap_thresh && curr_iter < block_iter; target *= 10) {
    if (avail_iter < target) target = primary_stats(target);
    for (; curr_iter < target; curr_iter++) {if (stats[curr_iter] >= value) curr_count++;}
  }
  
  if (curr_iter == block_iter && curr_count < adap_thresh) {
    if (!avail_sub) init_subblocks();
    double* buffer = perm_stats.data() + block_iter; int offset = 0;
    for (int target = curr_iter*10; curr_count < adap_thresh && curr_iter < max_iter; target *= 10) {    
      secondary_stats(buffer, offset++); curr_iter += block_iter;             
      for (int i = 0; i < block_iter; i++) {if (buffer[i] >= value) curr_count++;}    
    }
  }

  if (curr_count > 0) return curr_count / double(curr_iter);
  else return 0.5 / max_iter;
}

void ConvertChisqSumDf1ToPvalImhof::process_lambda() {
  lambda_sq.resize(ncomp);
  for (int i = 0; i < ncomp; i++) lambda_sq[i] = lambda[i]*lambda[i];
  
  bounds.set_lambda(lambda, ncomp);
  if (backstop) backstop->set_lambda(lambda, ncomp);
}


double ConvertChisqSumDf1ToPvalImhof::ImhofFunc::operator() (double u) const {
  double theta = 0, rho = 1, usq = u*u;
  for (int i = 0; i < ncomp; i++) {
    theta += atan(u*lambda[i]);
    rho *= 1 + lambda_sq[i]*usq;
  }
  theta = (theta - obs*u)/2;
  rho = pow(rho, 0.25);
  return sin(theta) / (u*rho);
}

using boost::math::quadrature::gauss_kronrod;
double ConvertChisqSumDf1ToPvalImhof::convert_internal(double value) { 
  double ref_pval = bounds.convert(value), pi = boost::math::constants::pi<double>();
  if (ref_pval <= 0) ref_pval = numeric_limits<double>::min();
  double prob = -1, log_diff = 0, integral, stored_prob = 2;
  ImhofFunc imhof(value, lambda.data(), lambda_sq.data(), ncomp);  

  for (int mode = 0; prob < 0 && mode < 3; mode++) {
    switch (mode) {
      case 0: integral = gauss_kronrod<double, 61>::integrate(imhof, 0, std::numeric_limits<double>::infinity(), 15, epsilon); break;
      case 1: integral = gauss_kronrod<double, 61>::integrate(imhof, 0, std::numeric_limits<double>::infinity(), 25, epsilon); break;        
      case 2: integral = gauss_kronrod<double, 201>::integrate(imhof, 0, std::numeric_limits<double>::infinity(), 25, epsilon); break;        
    };    

    prob = 0.5 + 1/pi * integral;   
    if (prob > 1) prob = -1;
    if (prob == 0) prob = numeric_limits<double>::min();

    if (prob > 0) {
      log_diff = log10(prob) - log10(ref_pval);
      if (log_diff < bound_range.first) prob = -1;
      if (log_diff > bound_range.second) {stored_prob = min(stored_prob, prob); prob = -1;}
    }

    if (ref_pval < pval_cutoff) break; /// restrict to only first gauss_kronrod call if p-value very low, to avoid extreme runtimes
  }
  if (prob < 0 && stored_prob <= 1) prob = stored_prob;

  if (backstop) {
    log_diff = log10(prob) - log10(ref_pval);
    if (log_diff > bound_range.second && prob > backstop->minimum_pval()) prob = -1;
    if (prob < 0) prob = backstop->convert(value);
  }
  
  if (prob < 0) throw StatException("pval computation", "Imhof algorithm failed to converge");

  return prob;
} 



using boost::math::chi_squared;
void ConvertChisqSumDf2ToPval::convert(double* buff, long len) { try {
  chi_squared distcs(df);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the chi-square statistic is negative or not finite");
    buff[i] = cdf(complement(distcs, scale*(buff[i] <= min_value ? min_value : buff[i])));
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from chi-square sum", e.what());}}

using boost::math::chi_squared;
double ConvertChisqSumDf2ToPval::convert(double value) { try {
  chi_squared distcs(df);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the chi-square statistic is negative or not finite");
  return cdf(complement(distcs, scale*(value <= min_value ? min_value : value)));
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from chi-square sum", e.what());}}


using boost::math::students_t;
void ConvertTToPval::convert(double* buff, long len) { try {
  students_t distt(df);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the t-statistic is not finite");
    buff[i] = cdf(distt, -buff[i]);
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from t-statistic", e.what());}}

using boost::math::students_t;
double ConvertTToPval::convert(double value) { try {
  students_t distt(df);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the t-statistic is not finite");
  return cdf(distt, -value);
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from t-statistic", e.what());}}



using boost::math::students_t;
void ConvertPvalToT::convert(double* buff, long len) { try {
  students_t distt(df);
  for (long i = 0; i < len; i++) buff[i] = quantile(complement(distt, buff[i] >= max_value ? max_value : (buff[i] <= min_value ? min_value : buff[i])));
} catch (const exception& e) {throw StatException("stat computation", "an error occurred when trying to convert p-value to t-statistic", e.what());}}

using boost::math::students_t;
double ConvertPvalToT::convert(double value) { try {
  students_t distt(df);
  return quantile(complement(distt, value >= max_value ? max_value : (value <= min_value ? min_value : value)));
} catch (const exception& e) {throw StatException("stat computation", "an error occurred when trying to convert p-value to t-statistic", e.what());}}


using boost::math::fisher_f;
void ConvertFToPval::convert(double* buff, long len) { try {
  fisher_f distf(df1, df2);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the F-statistic is not finite");
    buff[i] = cdf(complement(distf, (buff[i] <= min_value ? min_value : buff[i])));
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from F-statistic", e.what());}}

using boost::math::fisher_f;
double ConvertFToPval::convert(double value) { try {
  fisher_f distf(df1, df2);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the F-statistic is not finite");
  return cdf(complement(distf, (value <= min_value ? min_value : value)));
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from F-statistic", e.what());}}


void ConvertBetaToPval::convert(double* buff, long len) { try {
  boost::math::beta_distribution<> distb(shape1, shape2);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the beta statistic is not finite");
    if (buff[i] < min_value) buff[i] = min_value;
    else if (buff[i] > max_value) buff[i] = max_value;
    buff[i] = cdf(distb, buff[i]);
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from beta statistic", e.what());}}

double ConvertBetaToPval::convert(double value) { try {
  boost::math::beta_distribution<> distb(shape1, shape2);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the beta statistic is not finite");
  if (value < min_value) value = min_value;
  else if (value > max_value) value = max_value;
  return cdf(distb, value);
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from beta statistic", e.what());}}

double ConvertFittedBetaToPval::fit_param(double* pval, int len) {double mean, var;
  MathUtils::get_stats(pval, len, mean, var, false);

  if (var > 0) { 
    float sub = mean * (1-mean) / var - 1;
    shape1 = mean * sub; shape2 = (1-mean) * sub;
  } else shape1 = shape2 = 0;
  if (shape1 <= 0 || shape2 <= 0) throw StatException("pval computation", "invalid parameter values for fitted beta distribution");

  double dim = shape2 / shape1;
  
  use_linear = (dim < FITTEDBETA_THRESH);
  if (use_linear) {
    lin_buffer.set_size(len,2);
    memcpy(lin_buffer[0], pval, len*sizeof(double));
    to_logp.convert(lin_buffer[0], len);
    sort(lin_buffer[0], lin_buffer[1]);
    
    double offset = log(len+1);
    for (int i = 0; i < len; i++) lin_buffer(i,1) = log(len-i) - offset;
    double* coeff = MathUtils::simple_regression(lin_buffer[0], lin_buffer[1], len);
    intercept = coeff[0]; slope = coeff[1]; delete[] coeff;
  }
  
  return dim;
}

void ConvertFittedBetaToPval::convert(double* buff, long len) { 
  if (use_linear) {
    to_logp.convert(buff, len);
    for (double* end = buff + len; buff < end; buff++) *buff = exp(intercept + (*buff)*slope);
  } else ConvertBetaToPval::convert(buff, len);
}

double ConvertFittedBetaToPval::convert(double value) {
  if (use_linear) return exp(intercept + to_logp.convert(value)*slope);
  else return ConvertBetaToPval::convert(value);
}

void ConvertGammaToPval::convert(double* buff, long len) { try {
  boost::math::gamma_distribution<> distg(shape, scale);
  for (long i = 0; i < len; i++) {
    if (!(boost::math::isfinite)(buff[i])) throw StatException("pval computation", "the gamma statistic is not finite");
    if (buff[i] < min_value) buff[i] = min_value;
    buff[i] = cdf(complement(distg, buff[i]));
  }
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from gamma statistic", e.what());}}

double ConvertGammaToPval::convert(double value) { try {
  boost::math::gamma_distribution<> distg(shape, scale);
  if (!(boost::math::isfinite)(value)) throw StatException("pval computation", "the gamma statistic is not finite");
  if (value < min_value) value = min_value;
  return cdf(complement(distg, value));
} catch (const exception& e) {throw StatException("pval computation", "an error occurred when trying to compute p-value from gamma statistic", e.what());}}

void ConvertGammaToPval::fit_param(double* logp, int len) {double mean, var;
  MathUtils::get_stats(logp, len, mean, var, false);
  if (mean > 0 && var > 0) { 
    shape = mean*mean / var;
    scale = var / mean;
  } else shape = scale = 0;
  if (shape <= 0 || scale <= 0) throw StatException("pval computation", "invalid parameter values for fitted gamma distribution");
}

double ConvertParetoToPval::convert(double value) {
  if (value > max_value) value = max_value;
  return proportion*(this->*conv)(value);
}
void ConvertParetoToPval::convert(double* buff, long len) {
  for (long i = 0; i < len; i++) { 
    if (buff[i] > max_value) buff[i] = max_value;
    buff[i] = proportion*(this->*conv)(buff[i]);
  }
}

double ConvertParetoToPval::fit_param(double* perm, int len, double tail) {
  int use = min((tail > 1) ? tail : len*tail, len*PARETO_MAX_TAIL);
  proportion = use / double(len);
  
  bool failed = false;
  if (use > 100) {
    sort(perm, perm+len);
    if (use == len) cutoff = perm[0];
    else {
      perm += len-use-1;
      cutoff = (perm[0] + perm[1])/2; perm++;
    }
    double mean, var; 
    MathUtils::get_stats(perm, use, mean, var, false); mean -= cutoff;
    if (mean > 0 && var > 0) {
      shape = 0.5*(1-mean*mean/var);
      if (shape < 0.5) {
        if (abs(shape) < 1e-5) shape = 0;
        scale = mean*(1-shape);
      } else failed = true;    
    } else failed = true;
  } else failed = true;

  if (failed) throw StatException("pval computation", "could not fit generalized pareto distribution to permutation data");

  set_pfunc();
  return cutoff; 
}
          
double ConvertPermutationToPval::convert(double value) {return !use_noop ? ((value > threshold) ? fgpd->convert(value) : fgamma->convert(value)) : noop.convert(value);}
void ConvertPermutationToPval::convert(double* buff, long len) {
  if (!use_noop) {
    for (long i = 0; i < len; i++) buff[i] = (buff[i] > threshold) ? fgpd->convert(buff[i]) : fgamma->convert(buff[i]);
  } else noop.convert(buff, len);
}

double ConvertPermutationToPval::fit_param(double* perm, int len, bool keep_ordered) {
  use_noop = false;
  fgamma->fit_param(perm, len);
  if (len >= min_perm) {
    int use = max(int(tail_prop*len), tail_count);
    if (keep_ordered) {fit_buffer.assign(perm, len); perm = fit_buffer.begin();}    
    try {threshold = fgpd->fit_param(perm, len, use);}
    catch (const StatException& se) {threshold = thresh_max;}
  } else threshold = thresh_max;
  return threshold;
}


void ConvertPvalToLogPval::convert(double* buff, long len) {
  for (long i = 0; i < len; i++) {
    buff[i] = -2*log(buff[i] >= max_value ? max_value : (buff[i] <= min_value ? min_value : buff[i]));
  }
}
double ConvertPvalToLogPval::convert(double value) {
  return -2*log(value >= max_value ? max_value : (value <= min_value ? min_value : value));
}

void ConvertLogPvalToPval::convert(double* buff, long len) {
  for (long i = 0; i < len; i++) buff[i] = buff[i] > 0 ? exp(-0.5*buff[i]) : 1;
}
double ConvertLogPvalToPval::convert(double value) {return value > 0 ? exp(-0.5*value) : 1;}


void FudgePval::convert(double* buff, long len) {
  if (internal) internal->convert(buff, len);
  double corr;
  for (long i = 0; i < len; i++) {
    buff[i] = buff[i] >= max_value ? max_value : (buff[i] <= min_value ? min_value : buff[i]);
    corr = pow(fudge, log10(buff[i]));
    buff[i] = pow(buff[i], max(fudge_min,corr));
  }
}

double FudgePval::convert(double value) {
  if (internal) value = internal->convert(value);
  value = value >= max_value ? max_value : (value <= min_value ? min_value : value);
  double corr = pow(fudge, log10(value));
  return pow(value, max(fudge_min,corr));
}

double BoundedMinP::constrain(double pval, double minp) {
  if (no_tests == 1) return pval;
  if (log_mode) minp = exp(-0.5*minp);
  
  if (pval < minp) return minp;
  pval = min(pval, minp*no_tests); ///Bonferroni upper bound
  if (minp > numeric_limits<double>::epsilon()) { /// use more exact bound if possible
    double upper_exact = 1 - pow(1 - max(minp, numeric_limits<double>::epsilon()), no_tests);
    if (upper_exact > minp) pval = min(pval, upper_exact);
  }

  return pval;
} 

void BoundedMinP::convert(double* buff, long len) {if (!internal) return;
  if (len <= bounded_max && no_tests > 1) {
    minp_buffer.assign(buff, len);
    internal->convert(buff, len);
    for (int i = 0; i < len; i++) buff[i] = constrain(buff[i], minp_buffer(i));
  } else internal->convert(buff, len);
}

double BoundedMinP::convert(double value) {
  if (internal) value = constrain(internal->convert(value), value);
  return value;
}

void CompoundConverter::convert(double* buff, long len) {
  conv_in->convert(buff, len);
  conv_out->convert(buff, len);  
}

double CompoundConverter::convert(double value) {
  return conv_out->convert(conv_in->convert(value));
}



double ConvertFromUnknown::rank(double* buff, long len, double obs, bool comp_obs) {
  vector<long>& index = sorter.run(buff, len);

  if (comp_obs) {
    if (obs >= buff[index[0]]) obs = 0.5/len;
    else if (buff[index[len-1]] >= obs) obs = 1 - 0.5/len;
    else {
      long from = 0, to = len, curr;
      while (from < to-1) {
        curr = (from+to)/2;
        if (obs > buff[index[curr]]) to = curr;
        else from = curr; 
      }
      obs = (from+1) / double(len);
    }
  }

  double div = len+1;  
  for (long i = 0; i < len; i++) buff[index[i]] = (i+1) / div; 

  return obs;
}

void ConvertFromUnknown::convert(double* buff, long len) {
  rank(buff, len);
  internal->convert(buff, len);  
}

void ConvertFromUnknown::convert_obs(double* buff, long len, double& obs) {
  obs = rank(buff, len, obs, true);
  internal->convert_obs(buff, len, obs);
}

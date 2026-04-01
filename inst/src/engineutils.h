/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef ENGINEUTILS_H
#define ENGINEUTILS_H

#include "permutation.h"
#include "statutils.h"

#define EU_NULL_PVAL -1

class EngineParam {
public:  
  bool auto_clear;

protected:
  double default_value; 
  double dummy_value;
  int no_elem;

public:
  EngineParam() :  auto_clear(false), default_value(0) {}
  virtual ~EngineParam() {}

  virtual void clear() = 0;
  virtual bool in_range(int& index) {return (index >= 0 && index < no_elem);}

  virtual double& element(int index=0) = 0;
  double get(int index=0) {return in_range(index) ? element(index) : default_value;}
  virtual Buffer<double>& get_all(Buffer<double>& buff) {buff.set_size(no_elem,1); for (int i = 0; i < no_elem; i++) buff(i) = element(i); return buff;}
  void set(double val, int index=0) {element(index) = val;}  
  virtual void set(EngineParam* input) {if (!input) return;
    int count = min(no_elem, input->size());
    for (int i = 0; i < count; i++) element(i) = input->element(i);
  }

  int size() {return no_elem;}  
  void set_default(double val, bool reset=true) {default_value = val; if (reset) clear();}
  
  virtual EngineParam* clone() = 0;

}; 

inline ostream& operator<<(ostream& os, EngineParam& par) {
  os << par.get(0);
  for (int i = 1; i < par.size(); i++) os << "\t" << par.get(i);
  return os;
}

class MonoParam : public EngineParam {
  double value;

public:
  MonoParam() {value = default_value; no_elem = 1;}
  virtual ~MonoParam(){}
  virtual void clear() {value = default_value;} 

  virtual bool in_range(int& index) {return index == 0;}
  virtual double& element(int index=0) {return index==0 ? value : dummy_value;}
  
  MonoParam* clone() {return new MonoParam(*this);}
};

class BiParam : public EngineParam {
  double values[2];

public:
  BiParam() {values[0] = values[1] = default_value; no_elem = 2;}
  virtual ~BiParam(){}
  virtual void clear() {values[0] = values[1] = default_value;}  

  virtual double& element(int index=0) {return in_range(index) ? values[index] : dummy_value;}

  BiParam* clone() {return new BiParam(*this);}
};

class MultiParam : public EngineParam {
  vector<double> values;

public:
  MultiParam(int len=3) {resize(len);}
  virtual ~MultiParam(){} 
  virtual void clear() {values.assign(no_elem, default_value);}  

  virtual double& element(int index=0) {return in_range(index) ? values[index] : dummy_value;}
  virtual void set(EngineParam* input) {if (!input) return;
    resize(input->size());
    EngineParam::set(input);
  }
  
  void resize(int len) {values.resize(len, default_value); no_elem = len;}
  MultiParam* clone() {return new MultiParam(*this);}
};


class AggregateParam {
public: enum AggregateType{Sum, Mean, Min, Max, Product};
  double default_value;

protected:
  AggregateType type;
  int tot_elem;
  int tot_values;

  Buffer<double> values;
  vector<int> partitions;

  int curr_elem;
  bool computed;

  double compute_base_core(double* data, int len);
  double compute_base(int index, double scale);

public:
  AggregateParam(AggregateType type, int tot_elem=0) : default_value(0), type(type), tot_elem(tot_elem), tot_values(1) {
    if (type == Product) default_value = 1;
    clear();
  }

  void add(double value) {if (curr_elem < tot_elem) values(curr_elem++) = value;}
  void add(EngineParam* par) {if (!par) return;
    if (curr_elem >= tot_elem) return;
    if (curr_elem == 0) {tot_values = par->size(); clear();}  
    
    int count = min(tot_values, par->size()); 
    for (int i = 0; i < count; i++) values(curr_elem,i) = par->get(i);;
    curr_elem++;
  }
  AggregateParam& operator<< (double value) {add(value); return *this;}
  AggregateParam& operator<< (EngineParam* par) {add(par); return *this;}
  
  double compute(double scale=1, bool recompute=false);
  void extract(EngineParam* target, double scale=1, bool recompute=false);
  void set_partition(vector<int>& sizes) {partitions = sizes;}

  void set_size(int size) {tot_elem = size;}
  void rewind() {curr_elem = 0; computed = false;}
  void clear() {
    rewind();
    values.set_size(tot_elem+1, tot_values, true);
    if (default_value) values.assign_value(default_value);
    partitions.clear();
  }
};

class AggregateParamBlock {
  vector<AggregateParam*> params;

public:
  AggregateParamBlock() {}
  ~AggregateParamBlock() {for (int i = 0; i < params.size(); i++) delete params[i];}

  AggregateParam* add_param(AggregateParam::AggregateType type=AggregateParam::Sum) {
    AggregateParam* out = new AggregateParam(type);
    params.push_back(out);
    return out;
  }
  
  void set_size(int size) {for (int i = 0; i < params.size(); i++) params[i]->set_size(size);}        
  void clear_param() {for (int i = 0; i < params.size(); i++) params[i]->clear();}
};


class Aggregator {
public: enum AggregateType{Mean, Top};                                 
protected:
  enum MultiState {ms_Invalid, ms_SepCorr, ms_Count};
  enum Processed {pr_None, pr_Checked, pr_Computed};

  ConvertPvalToLogPval logp;
  StatConverter* to_pval;
  vector<StatConverter*> multi_pval;

  int no_pval;
  int no_tests;
  
  int pval_added;  
  int expected;
  Processed processed;

  double scale[2]; 
  Buffer<double> stats;
  Buffer<double> pval;
   
  Buffer<double> base_correlations;  
  vector<Buffer<double>* > multi_correlations;
  Buffer<double> corr_buffer;
  Buffer<double>& correlations(int index=0, bool create=false);
  Buffer<double>& sub_correlations(int index);
  

  Buffer<short> multi_status;
  vector<int> subset;
  
  void process();
  virtual void process_core() = 0;

  bool has_subset() {return !subset.empty() && (subset.size() < no_pval);}
  int no_subset() {return has_subset() ? subset.size() : no_pval;}
  virtual void check_size(Buffer<double>& corr, bool fail=true);

public:
  Aggregator() : to_pval(0) {multi_correlations.push_back(&base_correlations); rewind();}
  virtual ~Aggregator() {
    for (int i = 1; i < multi_correlations.size(); i++) delete multi_correlations[i];
    for (int i = 1; i < multi_pval.size(); i++) delete multi_pval[i];    
    multi_correlations.clear(); multi_pval.clear();
  }

  virtual void rewind(int size_override=0, int ntest=1);
  void clear_corrs();
  void add_pval(double p) {if (pval_added < expected) pval(pval_added) = p; pval_added++;}

  Buffer<double>& get_corrs(int index=0) {return correlations(index);}
  void set_corrs(Buffer<double>& cov, bool is_corr=false, int index=0);
  int set_perm_corrs(Buffer<float>& perm, int use=0, int index=0);
  void set_subset(vector<int>* sub=0);
    
  int get_size() {return no_subset();}
  double get_stat(int index=0) {if (processed < pr_Computed) process(); return stats(index);}
  double get_pval(int index=0);
  Buffer<double>& get_all_pval() {return pval;}
  double get_scale(bool root=false);
  StatConverter* get_converter() {return to_pval;}
};   

class MeanAggregator : public Aggregator {
  ConvertChisqSumDf2ToPval brown; 
  virtual void process_core();
public:
  MeanAggregator() {to_pval = &brown;}
  virtual ~MeanAggregator(){}
};   

class FittedAggregator : public Aggregator {
protected:
  CompoundConverter z_to_logp;  
  ConvertPermutationToPval gpd;

  SimulationMVN* fit_perm;
  Buffer<double> perm_buffer;
  Buffer<float> stat_buffer;

  int min_perm;
  int max_perm;
  float tail_thresh;
  float tail_margin_scale;
  int tail_min_perm;   
  
  int adap_count; 
  int block_size;
  
  int curr_perm; 
  vector<double> zstat;

  virtual void process_core();
  virtual void compute_stats() = 0;

  virtual ConvertPermutationToPval* gpd_converter(int index);
  void prep_converter(int index);
  virtual void process_perms(Buffer<float>& perm, int total) = 0;
  void advance_perm(int curr_sign);
  
public:
  FittedAggregator(Settings& settings) : z_to_logp(StatType::Z, StatType::LogPvalue), block_size(10000) {
    min_perm = settings.geti("gene_gpd_min_sim"); 
    max_perm = settings.geti("gene_gpd_nsim"); 
    adap_count = settings.geti("gene_gpd_adap_count");
    tail_thresh = settings.getn("gene_gpd_tail");
    tail_margin_scale = settings.getn("gene_gpd_tail_margin") + 1;    
    tail_min_perm = min(int(settings.geti("gene_gpd_tail_min")), int(0.75*max_perm)); 

    fit_perm = new SimulationMVN(max_perm, 1000, 0, settings.getn("gene_gpd_prune"), false);
    
    gpd.set_param("prop", settings.getn("gene_gpd_tail_fit"));
    gpd.set_param("min_count", 1000);    
    gpd.set_param("min_perm", tail_min_perm); 

    to_pval = &gpd;  
  }
  virtual ~FittedAggregator(){delete fit_perm;}

  virtual void rewind(int size_override=0, int ntest=1) {Aggregator::rewind(size_override, ntest); perm_buffer.set_size(0,0); curr_perm = 0;}
};

class TopAggregator : public FittedAggregator {
  ConvertPvalToNorm to_z;
  BoundedMinP bounded_gpd;

  virtual void compute_stats();
  virtual ConvertPermutationToPval* gpd_converter(int index);
  virtual void process_perms(Buffer<float>& perm, int total);
public:
  TopAggregator(Settings& settings) : FittedAggregator(settings) {
    bounded_gpd.set_internal(&gpd, true);
    to_pval = &bounded_gpd;
  }
  virtual ~TopAggregator() {}
};    


#endif /** ENGINEUTILS_H*/
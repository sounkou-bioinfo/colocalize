/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENECORRELATIONS_H
#define GENECORRELATIONS_H

#include "baseinput.h"
#include "geneutils.h"
#include "engineutils.h"
#include "genestats.h"
#include "genemodelengine_base.h"
#include "buffer.h"

using namespace std;
class CorrelationData;

class GeneCorrelations {
  Buffer<float> mat_buffer;

  GeneStats* stats;
  deque<CorrelationData*> genes;  

  bool inactive;
  long long max_range;

  void discard(CorrelationData* data);

public:
  GeneCorrelations() : stats(0), inactive(true) {}
  ~GeneCorrelations() {clear();}

  void clear();
  void set_range(long long range) {inactive = false; max_range = range;}
  void set_data(GeneStats* gs) {stats = gs;}

  void compute(BaseBuffer<double>& corrs, CorrelationData* data);  
};

class CorrelationUtils {
  static int top_power;
  static double top_corr_weight;
  static double top_corr_pow[3];

  static const double powersum_coefficients[3][5];   
  static double powersum_sum(double* rsq, int pow);

public:
  static void set_top_param(Settings& settings);  

  template<typename T>
  static double variance_chisq(Buffer<T>& xtx, bool pre_squared) {
    if (pre_squared) {
      double out = 0; int size = xtx.nrow();
      for (int i = 0; i < size; i++) {
        out += xtx(i,i);
        for (int j = i+1; j < size; j++) out += 2*xtx(j,i);
      }    
      return out;
    } else return MathUtils::corr_trace(xtx);
  }  
  
  template<typename T>
  static double covariance_chisq(Buffer<T>& corrs, bool pre_squared, double scale_to_corr=1) {
    if (pre_squared) return scale_to_corr*MathUtils::sum(corrs.begin(), corrs.length());
    else return scale_to_corr*scale_to_corr*MathUtils::prod_trace(corrs);
  }

  template<typename T>
  static pair<double,double> variance_powersum(Buffer<T>& xtx, bool pre_squared, int power) {double sum[] = {0,0,0,0}; int dim = xtx.nrow();
    for (int i = 0; i < dim; i++) {
      T &var = xtx(i,i), *xcol = xtx[i];
      if (pre_squared) {
        for (int j = i+1; j < dim; j++) {double curr = 1;
          double rsq = xcol[j] / sqrt(var / xtx(j,j));
          for (int j = 0; j < power; j++) {curr *= rsq; sum[j] += curr;}      
        }
      } else {
        for (int j = i+1; j < dim; j++) {double curr = 1;
          double rsq = xcol[j] * xcol[j] / var / xtx(j,j);
          for (int j = 0; j < power; j++) {curr *= rsq; sum[j] += curr;}      
        }
      }
    }
  
    pair<double,double> out;
    out.first = 2*sum[0] + dim;
    out.second = 2*powersum_sum(sum, power) + dim;    
    return out;
  }

  template<typename T>  
  static pair<double,double> covariance_powersum(Buffer<T>& corrs, bool pre_squared, int power, double scale_to_corr=1) {double sum[] = {0,0,0,0}; 
    if (!pre_squared) scale_to_corr *= scale_to_corr;
    for (int i = 0; i < corrs.length(); i++) {double curr = 1;
      double rsq = pre_squared ? corrs(i)*scale_to_corr : corrs(i)*corrs(i)*scale_to_corr;
      for (int j = 0; j < power; j++) {curr *= rsq; sum[j] += curr;}
    }    
  
    double cov = powersum_sum(sum, power);
    return pair<double,double>(sum[0], cov);
  }

  template<typename T>
  static bool variance_top_invert(Buffer<T>& xtx, double* target, bool pre_squared, double scale_full_block=1) {
    pair<double,double> var_parts = variance_powersum(xtx, pre_squared, top_power);
    double scale = 1/sqrt(scale_full_block);
    target[0] = (var_parts.first > 0) ? scale/sqrt(var_parts.first) : 0;
    target[1] = (var_parts.second > 0) ? scale/sqrt(var_parts.second) : 0;  
    return target[0] > 0 && target[1] > 0;
  }  
 
  template<typename T>
  static double covariance_top(Buffer<T>& corrs, double* inv_sd_i, double* inv_sd_j, bool pre_squared, double scale_to_corr=1) {
    pair<double,double> covar_parts = CorrelationUtils::covariance_powersum(corrs, pre_squared, top_power, scale_to_corr); 
    double plain_corr = pow(covar_parts.first * inv_sd_i[0] * inv_sd_j[0], top_corr_pow[0]);
    double power_corr = pow(covar_parts.second * inv_sd_i[1] * inv_sd_j[1], top_corr_pow[1]);
    return pow((1-top_corr_weight)*plain_corr + top_corr_weight*power_corr, top_corr_pow[2]);
  }
};                                



class CorrelationData {
friend class BlockCorrelationData;
protected:  
  double sd; 
  int id;
  
  template<typename T>
  T* convert(CorrelationData* other) {return (typeid(*this) == typeid(*other)) ? static_cast<T*>(other) : 0;}
  virtual double get_covariance_core(CorrelationData* other) = 0;
  virtual bool compatible_settings(CorrelationData* other) {return other != 0;}
public:    
  CorrelationData() : sd(-1), id(-1) {}
  virtual ~CorrelationData() {}

  virtual GeneLocation get_location() = 0; 

  bool check_stdev() {return sd > 0;}
  double get_stdev() {return sd;}
  double get_variance() {return sd >= 0 ? sd*sd : -1;}  
  double get_correlation(CorrelationData* other); 
  double get_covariance(CorrelationData* other);

  void set_stdev(double var) {sd = (var > 0) ? sqrt(var) : -1;}

  void check_readiness(CorrelationData* other);

  int& gene_id() {return id;}
  virtual bool is_compatible(CorrelationData* other) {return typeid(*this) == typeid(*other) && compatible_settings(other);}
  virtual bool has_data() = 0; 
  virtual void load_data() = 0;
  virtual GeneStats::DataBlock* eject_data() = 0;
}; 


class CoreCorrelationData : public CorrelationData {
friend class BlockCorrelationData;
protected:
  GeneStats::DataBlock* data;
  static Buffer<float> corr_buffer;
  
  int no_blocks;
  int curr_block;
  int curr_offset;
  int curr_size;
  Buffer<double> block_covar;
  Buffer<int> block_index;
  
  bool set_block(int id);  
  bool check_index(int i) {return i >= 0 && i < no_blocks;}
  bool check_index(int& i, int& j) {
    if (j < i) swap(i,j);
    return (check_index(i) && check_index(j) && i != j);
  }
  
  virtual Buffer<int>& get_partition() = 0;
  void set_partition();
  
public:
  CoreCorrelationData(GeneStats* stats) : data(0), curr_block(0) {no_blocks = stats->nblock(); block_covar.set_square(no_blocks);}
  virtual ~CoreCorrelationData() {delete data;}

  virtual void add_variance(double var, int block=0) {
    if (no_blocks == 1) {set_stdev(var); block_covar(0) = var;}
    else if (check_index(block)) block_covar(block,block) = var;    
  }
  virtual void compute_variance(Buffer<float>& xtx, int block=0, double scale_full_block=1);
  virtual void compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j) = 0;

  virtual GeneLocation get_location() {return data->location;}
  virtual GeneStats::DataBlock* eject_data() {
    GeneStats::DataBlock* out = data; data = 0;
    return out;    
  }  
};


class RawCorrelationData : public CoreCorrelationData {
protected:
  DefaultStats* def_stats;
  
  int min_step_dim;
  int max_step_size;
 
  int get_step_size(int dim);
  virtual Buffer<int>& get_partition() {return def_stats->get_partition(true);}

  double prep_covariance(CorrelationData* other_base, Buffer<float>& target); 
  virtual double get_covariance_core(CorrelationData* other);
public:
  RawCorrelationData(DefaultStats* stats, int min_step_dim, int max_step_size) : CoreCorrelationData(stats), def_stats(stats), min_step_dim(min_step_dim), max_step_size(max_step_size) {}
  virtual ~RawCorrelationData(){}

  virtual void compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j);

  virtual bool has_data() {return data != 0;}
  virtual void load_data() {data = def_stats->eject_gene(); set_block(-1);}
};


class PowerRawCorrelationData : public RawCorrelationData {
protected:
  Buffer<double> inverse_sds;

  virtual double get_covariance_core(CorrelationData* other);  
public:
  PowerRawCorrelationData(DefaultStats* stats, int min_step_dim, int max_step_size) : RawCorrelationData(stats, min_step_dim, max_step_size) {
    inverse_sds.set_size(2,no_blocks,true); 
  }
  virtual ~PowerRawCorrelationData(){}
  
  virtual void compute_variance(Buffer<float>& xtx, int block=0, double scale_full_block=1);
  virtual void compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j);
};      



class PCCorrelationData : public CoreCorrelationData {
protected:
  DefaultStats* def_stats;
  int pcid;

  virtual double covariance_function(Buffer<float>& corrs, double scale_to_corr=1);
  virtual double get_covariance_core(CorrelationData* other);
  virtual Buffer<int>& get_partition() {return def_stats->get_partition(false, pcid);}
  BufferWindow<float> get_buffer() {return data->get_buffer()->window(curr_offset, curr_size);}

public:
  PCCorrelationData(DefaultStats* stats, int pcid) : CoreCorrelationData(stats), def_stats(stats), pcid(pcid) {}
  virtual ~PCCorrelationData(){}

  virtual void compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j);

  virtual bool has_data() {return data != 0;}  
  virtual void load_data() {data = def_stats->eject_pcs(pcid); set_block(-1);}
};


class PowerPCCorrelationData : public PCCorrelationData {
protected:
  int power;
  double corr_weight;

  virtual double covariance_function(Buffer<float>& corrs, double scale_to_corr=1);

  virtual bool compatible_settings(CorrelationData* other);
public:
  PowerPCCorrelationData(DefaultStats* stats, int pcid, int power, double corr_weight=1) : PCCorrelationData(stats, pcid), power(power), corr_weight(corr_weight) {}
  virtual ~PowerPCCorrelationData(){}
};

class CompoundCorrelationData : public CorrelationData {
protected:
  Aggregator::AggregateType aggregate_type;
  double inverse_sds[2]; /// 1/sd for the two sub-components for Top; or general 1/sd for Mean (second value not used);

  virtual bool prep_variance() = 0;
  virtual void compute_variance(Buffer<double>& corrs) {compute_variance(corrs, prep_variance());}
  virtual void compute_variance(Buffer<double>& corrs, bool process_blocks);  
  double compute_covariance(Buffer<double>& corrs, double* inv_sd_i, double* inv_sd_j, bool to_corr=false) {return compute_covariance(corrs, inv_sd_i, inv_sd_j, aggregate_type, to_corr);}  
  double compute_covariance(Buffer<double>& corrs, double* inv_sd_i, double* inv_sd_j, Aggregator::AggregateType type, bool to_corr);    

  virtual bool compatible_settings(CorrelationData* other);  
public:
  CompoundCorrelationData(Aggregator::AggregateType type) : aggregate_type(type) {}
  virtual ~CompoundCorrelationData() {}
  
  Aggregator::AggregateType type() {return aggregate_type;}
};


class BlockCorrelationData : public CompoundCorrelationData {
protected:
  static Buffer<double> corr_buffer;
  CoreCorrelationData* data;

  int no_blocks;
  Buffer<double>* block_corr;
  Buffer<double> block_sds;  
  Buffer<int>* block_index;

  void block_init();
  virtual bool prep_variance();
  virtual double prep_covariance(Buffer<double>& buffer, BlockCorrelationData* other); 
  virtual double get_covariance_core(CorrelationData* other);
  
public:
  BlockCorrelationData(CorrelationData* cd, Aggregator::AggregateType type, bool do_init=true) : CompoundCorrelationData(type) {
    data = dynamic_cast<CoreCorrelationData*>(cd);
    if (!data) {
      delete cd;
      throw DataException("initializing block correlation data", "could not cast to CoreCorrelationData");
    }
    block_init();
    if (do_init) compute_variance(*block_corr);  
  }
  virtual ~BlockCorrelationData() {delete data;}

  virtual GeneLocation get_location() {return data->get_location();}
  CoreCorrelationData* get_internal() {return data;}

  Buffer<double>& get_corrs() {return *block_corr;}

  virtual bool is_compatible(CorrelationData* other_base) {
    BlockCorrelationData* other_block = convert<BlockCorrelationData>(other_base);
    if (other_block && aggregate_type == other_block->type()) return data->is_compatible(other_block->data);
    else return false;
  }
  virtual bool has_data() {return data->has_data();}  
  virtual void load_data() {data->load_data();}  
  virtual GeneStats::DataBlock* eject_data() {return data->eject_data();}
};

class MultiCorrelationData : public CompoundCorrelationData {
  vector<CorrelationData*> data;
  static Buffer<double> corr_buffer;
  Buffer<double> engine_corrs;
  Buffer<float> perm_buffer;
  int no_perm;
  bool fixed_mode;

  virtual bool prep_variance();
  virtual double get_covariance_core(CorrelationData* other);  

public:  
  MultiCorrelationData(Aggregator::AggregateType type, bool fixed_mode=true) : CompoundCorrelationData(type), no_perm(0), fixed_mode(fixed_mode) {}
  virtual ~MultiCorrelationData() {for (int i = 0; i < data.size(); i++) delete data[i];}

  GeneLocation get_location() {return (!data.empty()) ? data[0]->get_location() : GeneLocation();}

  void set_data(CorrelationData* input) {data.push_back(input);}
  void set_corrs(Buffer<float>& perms, int use=0) {set_corrs(engine_corrs, perms, use);}
  void set_corrs(Buffer<double>& corrs, Buffer<float>& perms, int use=0);  

  bool check_stdev() {
    if (sd <= 0) return false;
    for (int i = 0; i < data.size(); i++) {if (!data[i]->check_stdev()) return false;}
    return true;
  }
  
  virtual bool is_compatible(CorrelationData* other_base);
  virtual bool has_data();  
  virtual void load_data() {for (int i = 0; i < data.size(); i++) data[i]->load_data();}
  virtual GeneStats::DataBlock* eject_data();
}; 

#endif /** GENECORRELATIONS_H */


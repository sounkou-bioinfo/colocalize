/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "genestats.h"
#include "statutils.h"

using namespace std;
class PermutationEngine;
class PermutationBlock;

class Permutation {
protected: 
  PermutationGenerator* generator;
  DefaultStats* data;
  GeneStats::GSID data_id;

  double obs_stat;
  short direction;
  int sign_count[2];
  
  int no_failed;
  void check_failure(int total);  
  
public: 
  Permutation() : generator(0), data(0), direction(1) {}
  virtual ~Permutation() {delete generator;}

  virtual bool set_data(GeneStats* gs);
  virtual void rewind() {sign_count[0] = 0; sign_count[1] = 0; no_failed = 0;}
  virtual void set_obs(double stat, StatConverter* conv=0, short dir=1, bool rew=true) {
    if (rew) rewind(); else sign_count[0] = sign_count[1] = 0;
    obs_stat = stat; direction = dir; UNUSED(conv);
  }

  virtual bool next() = 0;
  virtual bool is_empty() = 0;
  virtual void failed() = 0;
  virtual void process(float value) {sign_count[1]++; if (direction*value >= direction*obs_stat) sign_count[0]++;}
  template<typename T>
  void process_bulk(T* values, int len) {for (T* end = values + len; values < end; values++) process(*values);}    

  virtual Buffer<float>& get_product(Buffer<float>& target, Buffer<float>& data) = 0;
  virtual BufferWindow<float> get_covar() = 0;

  virtual double pval() {return sign_count[1] > 0 ? sign_count[0] / double(sign_count[1]) : 1;}
  virtual int nsign() {return sign_count[0];}
  virtual int nperm() {return sign_count[1];}
  virtual int next_nperm() = 0;
  virtual int max_next_nperm() = 0;
  virtual int max_nperm() = 0;

  virtual PermutationBlock* get_block(bool from_start=false) {UNUSED(from_start); return 0;}
  virtual PermutationEngine* get_engine() {return 0;}
};


class PermutationBlock : public Permutation {
  BufferWindow<float> permutations;
  BufferWindow<float> covar;
  Buffer<float> perm_stats;
  float* curr_perm;
  Buffer<float>* attached;

  int no_perm;
  bool empty;    

  ConvertPvalToLogPval logp;
  StatConverter* to_pval;
  double convert(double stat);

  void compute_sign() {
    if (sign_count[0] >= 0) return; sign_count[0] = 0;
    for (int i = 0; i < sign_count[1]; i++) {if (perm_stats(i) >= obs_stat) sign_count[0]++;}
  }

  void set_converter(StatConverter* in) {
    delete to_pval;
    to_pval = in ? in->clone() : 0;
  }

  struct Process {
    virtual void obs(double value, double& stat) = 0;
    virtual void perm(float value, float*& buff) = 0;
  };
  
  struct ProcessMean : public Process {
    void obs(double value, double& stat) {stat += value;}  
    void perm(float value, float*& buff) {*(buff++) += value;}      
  } proc_mean;
  
  struct ProcessTop : public Process {
    void obs(double value, double& stat) {if (value > stat) stat = value;}  
    void perm(float value, float*& buff) {if (value > *buff) *buff = value; buff++;}
  } proc_top;

  Process* proc;

public: 
  PermutationBlock() : no_perm(0), empty(false), to_pval(0) {proc = &proc_mean;}
  virtual ~PermutationBlock() {delete to_pval;}

  void set_perm(BufferWindow<float> perm) {
    permutations = perm;
    no_perm = permutations.view().ncol();    
    perm_stats.set_size(no_perm, 1, true); ///for 'Top' processing, assumes input values are non-negative
    rewind(); obs_stat = 0;
  }
  void set_perm(BufferWindow<float> perm, BufferWindow<float> cov) {set_perm(perm); covar = cov;}
  void set_perm(Buffer<float>& perm, int offset=0, int total=0) {set_perm(perm.window(offset, total));}    
  void set_perm(Buffer<float>& perm, Buffer<float>& cov, int offset=0, int total=0) {set_perm(perm.window(offset, total), cov.window(offset, total));}  

  virtual void rewind() {Permutation::rewind(); sign_count[0] = -1; curr_perm = perm_stats.begin(); attached = 0; empty = false; set_converter(0);}      
  virtual void set_obs(double stat, StatConverter* conv=0, short dir=1, bool rew=true) {
    if (rew) rewind(); else sign_count[0] = -1; sign_count[1] = 0; 
    set_converter(conv); UNUSED(dir); proc->obs(convert(stat), obs_stat);
  }
  
  virtual bool next() {no_failed = 0; if (!empty) {empty = true; return true;} else return false;}
  virtual bool is_empty() {return empty;}
  virtual int next_nperm() {return !empty ? no_perm : 0;}
  virtual int max_next_nperm() {return next_nperm();} 
  virtual int max_nperm() {return no_perm;}
  virtual void failed() {no_failed++; check_failure(no_perm);}
  
  void proc_mode(bool use_top) {
    if (use_top) proc = &proc_top;
    else proc = &proc_mean;
  }
    
  virtual void process(float value) {proc->perm(convert(value), curr_perm); sign_count[1]++;}
  virtual Buffer<float>& get_product(Buffer<float>& target, Buffer<float>& data);
  virtual BufferWindow<float> get_covar();

  virtual double pval() {compute_sign(); return Permutation::pval();}
  virtual int nsign() {compute_sign(); return sign_count[0];}

  virtual PermutationBlock* get_block(bool from_start=false) {if (from_start) rewind(); return this;}

  void attach_perm(Buffer<float>* buffer) {attached = buffer;}
  Buffer<float>* get_attached() {return attached;}
  virtual void transform(StatConverter* conv);  
  virtual void extract(float* buffer, int len) {memcpy(buffer, perm_stats.begin(), len*sizeof(float));}
  virtual void merge(PermutationBlock* input);
  PermutationBlock* copy_block();
    
  int size(bool total) {return max((total ? no_perm : sign_count[1]), int(0));}
  
  void rescale();
};


class PermutationEngine : public Permutation {
  static const int core_max = 10000;
  static const int transient_max = 10000;

  deque<int> partitioning;
 
  Buffer<float> core_perm;
  Buffer<float> transient_perm;
  Buffer<float> core_covar;  
  Buffer<float> transient_covar;    

  PermutationBlock* main_block;
  PermutationBlock* add_block;
  int block_curr;  

  bool adaptive;
  bool empty;
  int min_perm;
  int max_perm;
  int adap_count;

  int curr_index;
  int curr_nperm;

  int get_needed(int curr);

public:
  PermutationEngine(int min) : main_block(0), add_block(0), block_curr(0), adaptive(false), empty(false), min_perm(min), max_perm(min), adap_count(0) {partitioning.push_back(min_perm);}
  PermutationEngine(int min, int max, int count) : main_block(0), add_block(0), block_curr(0), empty(false), min_perm(min), max_perm(max), adap_count(count) {
    if (max_perm < min_perm) max_perm = min_perm;
    adaptive = (max_perm > min_perm) && (adap_count > 0);
    partitioning.push_back(min_perm); 
    
    if (adaptive) {
      if (min_perm < 5000 && max_perm > 5000) partitioning.push_back(5000);
      if (min_perm < 10000 && max_perm > 10000) partitioning.push_back(10000);    
      partitioning.push_back(max_perm);
    }
  }  
  virtual ~PermutationEngine() {delete main_block; delete add_block;}

  bool set_data(GeneStats* gs);
  bool add_partition(int offset);
  void rewind() {Permutation::rewind(); curr_index = -1; curr_nperm = 0; block_curr = 0; empty = false;}  

  bool next();
  bool is_empty() {return empty || next_nperm() == 0;}
  void set_empty() {empty = true;}
  int next_nperm() {
    int index = curr_index >= 0 ? curr_index+curr_nperm : 0;
    if ((index >= max_perm) || (adaptive && sign_count[0] >= adap_count)) return 0;
    if (index < core_max) return min(get_needed(index), core_max) - index;
    else return min(get_needed(index)-index, transient_max);
  }
  int max_next_nperm() {
    if (adaptive) return min(max_perm, max(core_max, transient_max));
    else return next_nperm();
  }
  int max_nperm() {return max_perm;}
  
  void failed() {no_failed++; check_failure(curr_nperm);}
  Buffer<float>& get_product(Buffer<float>& target, Buffer<float>& data);
  BufferWindow<float> get_covar();

  PermutationBlock* get_block(bool from_start=false);
  PermutationEngine* get_engine() {return this;}

//  void insert(PermutationBlock* input) {if (input != main_block) main_block->merge(input);}
  void insert(PermutationBlock* input);
  void process_block(PermutationBlock* add=0);
};



class SimulationMVN : public Permutation {
  Buffer<float> sim_buffer;
  float prune_perc;

  int tot_perm;
  int max_dim;
  int avail_dim;
  int chunk_size;

  Buffer<float> projection;
  Buffer<float> means;
  Buffer<float> buffer;

  int curr_dim;
  int curr_nvar;
  int curr_index;
  int curr_nperm;
  
  int samp_size;
  int no_covar;
  bool load_covar;
  double resid_scale;
  
  Buffer<float> covariates;
  Buffer<float> covar_sims;
  Buffer<float> pheno_fitted;

  void simulate(int needed);
  int invert_covariance(Buffer<float>& covar);
  void shift(float* data, float mean, int step, int tot_step) {for (float* end = data+step*tot_step; data < end; data += step) *data += mean;}
public:
  SimulationMVN(int tot_perm, int max_dim, int init_dim=-1, float prune=0.95, bool load_covar=true) : prune_perc(prune), tot_perm(tot_perm), max_dim(max_dim), avail_dim(0), no_covar(0), load_covar(load_covar) {
    if (init_dim < 0) init_dim = max_dim/4.0;
    set_block(tot_perm/4);
    simulate(init_dim);
  }
  virtual ~SimulationMVN(){}

  bool set_data(GeneStats* gs);
  void rewind() {Permutation::rewind(); curr_dim = 0; curr_nvar = 0; curr_index = -1; curr_nperm = 0;}  

  int set_input_covariance(Buffer<double>& covar); 
  int set_input_covariance(Buffer<float>& covar); 
  int set_input_covariance(Buffer<float>& input, Buffer<float>& covar); 
  void set_block(int size, int dim=0) {chunk_size = max(size, 1); if (dim > 0) max_dim = dim;}
  void set_empty() {curr_index = tot_perm;}
  
  bool next();
  bool is_empty() {return curr_index >= tot_perm || chunk_size <= 0;} 
  void failed() {no_failed++; check_failure(tot_perm);}
  
  Buffer<float>& get_sims(Buffer<float>& target);
  Buffer<float>& get_product(Buffer<float>& target, Buffer<float>& data) {UNUSED(&data); return get_sims(target);}
  BufferWindow<float> get_covar();

  int next_nperm() {return !is_empty() ? min(chunk_size, tot_perm - curr_index) : 0;}
  int max_next_nperm() {return next_nperm();}
  int max_nperm() {return tot_perm;}
};

class MultiPhenoFWER {
  DefaultStats& data;
  
  Buffer<float> perm_buffer;
  Buffer<float> min_p;

  void set_values(int no_values);

public:
  MultiPhenoFWER(DefaultStats& data) : data(data) {}
  
  void generate(string outfile, int no_perm);

};


#endif /** PERMUTATION_H */


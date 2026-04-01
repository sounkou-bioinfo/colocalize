/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEMODELENGINE_BASE_H
#define GENEMODELENGINE_BASE_H

#include "utils.h"
#include "mathutils.h"
#include "statutils.h"
#include "genestats.h"
#include "permutation.h"
#include "engineutils.h"

#define EIGEN_CUTOFF 0.0001

class CorrelationData; class CoreCorrelationData; class PCCorrelationData;

class GeneModelEngine {
protected:
  Settings& settings;
  int gene_id;

  bool has_perm;
  Permutation* permutation;

  bool data_changed;
  bool gene_changed;
  bool needs_var;

  int sample_size;
  int no_snps;
  int no_param;
  int no_covar;
  bool has_covar;
  bool has_interactor;
  
  int no_tests;
  Buffer<double> pval;
  map<string,EngineParam*> parameters;

  template<typename T>
  T* register_param(const string& name, T* par, bool clear=true) {par->auto_clear = clear;
    if (has_param(name)) delete parameters[name];
    parameters[name] = par; 
    return par;
  }
  void delete_params() {
    for (map<string,EngineParam*>::iterator it = parameters.begin(); it != parameters.end(); ++it) delete it->second;
  }

  virtual GeneStats* get_data() = 0;
  virtual void set_data(GeneStats* input) = 0;
  virtual void init_gene() {
    no_tests = has_interactor ? 3 : 1;     
    pval.set_size(1,no_tests);      
    for (map<string,EngineParam*>::iterator it = parameters.begin(); it != parameters.end(); ++it) {if (it->second->auto_clear) it->second->clear();}
  }

  virtual double run_core() = 0;
public:
  GeneModelEngine(Settings& settings) : settings(settings), gene_id(-1), no_tests(1) {}
  virtual ~GeneModelEngine() {delete_params();}    

  virtual void clear() {}
  virtual double run(GeneStats* input, bool var, Permutation* perm = 0); 
  virtual void set(string setting, string value="") {UNUSED(setting); UNUSED(value);}

  virtual CorrelationData* corr_data() = 0;
  virtual ModelInfoBlock model_info(GeneStats* gs) = 0;
    
  virtual int ntest() {return no_tests;}
  virtual int get_nparam() {return no_param;}
  virtual double get_pval(int index=0) {return index < no_tests ? pval(index) : -1;}
  virtual void load_pval(Buffer<double>& target) {target.assign(pval);}

  bool has_param(const string& name) {return get_param(name) != 0;}
  EngineParam* get_param(const string& name, bool copy=false) {
    map<string,EngineParam*>::iterator found = parameters.find(name);
    if (found == parameters.end()) return 0;
    else return (copy) ? found->second->clone() : found->second;
  }
  double get_param_value(const string& name, int index=0, double def_value=0) {
    EngineParam* par = get_param(name);
    return par ? par->get(index) : def_value;
  }
 
  virtual bool perm_only() {return false;}
  virtual bool equals(GeneModelEngine* other) {return this == other;}  
};


class CoreEngine : public GeneModelEngine {
protected:
  int block_id; 
  int block_count;

  Buffer<float> data_buffer;
  Buffer<float> mat_buffer;
  Buffer<float> stat_buffer;
  
  virtual double prep_block_corr(Buffer<float>& buffer, int block_i, int block_j) = 0;
  void set_variance(CoreCorrelationData* cd, double var);
  void compute_variance(CoreCorrelationData* cd, Buffer<float>& xtx, double scale=1);
  void compute_covariance(CoreCorrelationData* cd);  

  int eigen_count(Buffer<float>& corr, double prune_cutoff=0.99);    

  virtual void init_gene();

  ModelInfoBlock make_info(GeneStats* gs, ModelState::Model model);
public:
  CoreEngine(Settings& settings) : GeneModelEngine(settings), block_id(0), block_count(1) {}
  virtual ~CoreEngine() {}
  virtual void clear() {GeneModelEngine::clear(); data_buffer.reset(); mat_buffer.reset(); stat_buffer.reset();}
};

class PrinCompEngine : public CoreEngine {
protected:
  DefaultStats* data;
  PCCorrelationData* covar;

  float prune_var_prop;
  float prune_min_val;
  int prune_drop;
  int pcid;
  
  virtual void prenormalize(Buffer<float>& mat, Buffer<float>& scale);
  virtual int prune_invert(Buffer<float>& mat, Buffer<float>& buffer, bool root);

  virtual GeneStats* get_data() {return data;}
  virtual void set_data(GeneStats* input) {data = input->convert<DefaultStats>();}
  virtual void init_gene();
 
  virtual double prep_block_corr(Buffer<float>& buffer, int block_i, int block_j);
 
public:  
  PrinCompEngine(Settings& settings) : CoreEngine(settings), covar(0), prune_var_prop(0.9999), prune_min_val(0), prune_drop(0), pcid(0) {} 
  virtual ~PrinCompEngine();

  virtual CorrelationData* corr_data();
  void set_pcid(int id) {pcid = id;}
  virtual void clear();
};

#include "genecorrelations.h"

#endif /** GENEMODELENGINE_BASE_H */

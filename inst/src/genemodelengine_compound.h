/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEMODELENGINE_COMPOUND_H
#define GENEMODELENGINE_COMPOUND_H

#include <map>
#include "genemodelengine.h"

class CompoundEngine : public GeneModelEngine {
protected:
  int no_sub_tests;
  int dimension;

  GeneStats* data;
  Aggregator* aggregator;
  Aggregator::AggregateType aggregate_type;
  AggregateParamBlock aggr_param;

  AggregateParam* no_param_aggr;
  AggregateParam* register_param_aggr(const string& name, AggregateParam::AggregateType type=AggregateParam::Sum);
  void update_dimension(int dim);

  virtual GeneStats* get_data() {return data;}
  virtual void set_data(GeneStats* input) {data = input;}
  void set_aggregator(Aggregator::AggregateType type) {set_aggregator(type, aggregator);}
  void set_aggregator(Aggregator::AggregateType type, Aggregator*& target);
  virtual void init_gene();

public:
  CompoundEngine(Settings& settings, Aggregator::AggregateType aggregate_type=Aggregator::Mean) : GeneModelEngine(settings), dimension(0), aggregator(0), aggregate_type(aggregate_type) {}
  virtual ~CompoundEngine() {delete aggregator;}
};

class PartitionedEngine : public CompoundEngine {
protected:
  int block_max;
  float block_prop;

  int no_blocks;
  int block_size;
  bool block_mode;

  AggregateParam* eff_size_aggr;
  AggregateParam* eff_size_adjusted_aggr;

  virtual void init_gene();
  void load_param(const string& name, GeneModelEngine* en);  
  int compute_blocks(vector<int>& blocks, int total);
  vector<int> compute_blocks(int total);
  virtual void set_blocks();

public:
  PartitionedEngine(Settings& settings, Aggregator::AggregateType aggregate_type=Aggregator::Mean) : CompoundEngine(settings, aggregate_type) {
    block_max = settings.geti("gene_chunked_size", 0);
    block_prop = settings.getn("gene_chunked_fraction", 0);

    no_param_aggr = aggr_param.add_param();    
    eff_size_aggr = register_param_aggr("eff_size"); 
    eff_size_adjusted_aggr = register_param_aggr("eff_size_adjusted"); 
  }
  virtual ~PartitionedEngine() {}
  virtual void clear() {GeneModelEngine::clear();}
};

class BlockEngine : public PartitionedEngine {
protected:
  Buffer<float> mat_buffer;

  CoreEngine* engine;
  bool owner; 

  PermutationBlock* curr_perm;
  BlockCorrelationData* covar;

  virtual void run_single();
  virtual void run_blocks();

  virtual double run_core();

  virtual void set_covar();
  virtual void set_perm(Permutation* source, int block_id, bool rewind);
  virtual PermutationBlock* process_perm(StatConverter* convert = 0);
  virtual void process_pval();
  
  virtual void init_gene() {
    PartitionedEngine::init_gene();
    delete covar; covar = 0;
    curr_perm = 0;
  }
                                  
public:
  BlockEngine(Settings& settings, CoreEngine* en, Aggregator::AggregateType aggregate_type=Aggregator::Mean) : PartitionedEngine(settings, aggregate_type), engine(0), owner(false), covar(0) {
    if (en) set_engine(en);
    set_aggregator(aggregate_type);
  }
  virtual ~BlockEngine() {if (engine && owner) delete engine; mat_buffer.reset(); delete covar;}

  virtual void set_engine(CoreEngine* en, bool transfer=true) {
    if (engine && owner) delete engine;
    engine = en; owner = transfer;
  }

  virtual CorrelationData* corr_data() {
    CorrelationData* out = covar; covar = 0;
    return out;
  }

  virtual bool equals(GeneModelEngine* other) {return GeneModelEngine::equals(other) && engine->equals(static_cast<BlockEngine*>(other)->engine);}
  virtual void clear() {PartitionedEngine::clear(); mat_buffer.reset(); delete covar; covar = 0;}  
  virtual ModelInfoBlock model_info(GeneStats* gs) {return engine->model_info(gs);}    
};

class MultiEngine : public CompoundEngine {
  Buffer<float> mat_buffer;

  int no_engine;
  vector<GeneModelEngine*> engine_bay;
  vector<bool> engine_owner;

  AggregateParam* snp_counts_aggr;

  MultiCorrelationData* covar;

  PermutationEngine* corr_perm;
  Buffer<float> perm_buffer;
  vector<Buffer<float>*> block_perm;
  
  int no_cperm;
  int used_perm;

  virtual void init_gene();    
  void set_covar();  
  void set_nsub(int num) {
    no_sub_tests = num; no_tests = no_sub_tests * (no_engine+1);
    pval.set_size(no_sub_tests, no_engine+1); 
    pval.assign_value(EU_NULL_PVAL); 
    aggregator->rewind(dimension, no_sub_tests);    
  }
  void store_block_perm(Buffer<float>& perms, int index);
  
  double run_core();    
public:
  MultiEngine(Settings& settings, Aggregator::AggregateType aggregate_type=Aggregator::Mean) : CompoundEngine(settings, aggregate_type), no_engine(0), covar(0), corr_perm(0) {
    set_aggregator(aggregate_type);
    no_param_aggr = aggr_param.add_param(AggregateParam::Mean); 
    snp_counts_aggr = register_param_aggr("snp_counts", AggregateParam::Max); 

    no_cperm = settings.geti("gene_multi_cperm", 250);    
  }
  virtual ~MultiEngine() {
    for (int i = 0; i < engine_bay.size(); i++) {if (engine_owner[i]) delete engine_bay[i];}
    for (int i = 0; i < block_perm.size(); i++) delete block_perm[i];
    delete corr_perm; delete covar;
  }
  virtual void clear() {
    GeneModelEngine::clear(); mat_buffer.reset(); 
    delete corr_perm; corr_perm = 0;
    delete covar; covar = 0;  
  }

  virtual void add_engine(GeneModelEngine* engine, bool owner=false);
  virtual CorrelationData* corr_data() {
    CorrelationData* out = covar; covar = 0;
    return out;
  }
  ModelInfoBlock model_info(GeneStats* gs);
  
  virtual bool equals(GeneModelEngine* other);
};


#endif /** GENEMODELENGINE_COMPOUND_H */

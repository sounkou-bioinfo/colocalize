/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genemodelengine_compound.h"

#define MULTI_MIN_PERM 50

AggregateParam* CompoundEngine::register_param_aggr(const string& name, AggregateParam::AggregateType type) {
  register_param(name, new MultiParam(1), true);
  return aggr_param.add_param(type);
}
                                                          
void CompoundEngine::update_dimension(int dim) {
  if (dimension != dim) {dimension = dim;
    aggregator->rewind(dimension);
    aggr_param.set_size(dimension);
  } 
}

void CompoundEngine::set_aggregator(Aggregator::AggregateType type, Aggregator*& target) {
  delete target;
  switch (type) {
    case Aggregator::Mean: target = new MeanAggregator(); break;
    case Aggregator::Top: target = new TopAggregator(settings); break;    
    default: target = 0;
  }
}

void CompoundEngine::init_gene() {
  GeneModelEngine::init_gene();
  no_sub_tests = no_tests;
  data_changed = data->set_current(gene_id);
}

void PartitionedEngine::init_gene() {
  CompoundEngine::init_gene();
  block_size = block_max ? block_max : 100000000;
  if (block_prop) block_size = max(min(int(block_prop*data->get_sampsize()), block_max), 1);
  set_blocks();
}

void PartitionedEngine::load_param(const string& name, GeneModelEngine* en) {
  EngineParam* par = get_param(name);
  if (par) par->set(en->get_param(name));
}  

int PartitionedEngine::compute_blocks(vector<int>& blocks, int total) {
  int nblock = total / block_size + (total % block_size > 0);
  int base_size = total / nblock, add = total % nblock;
  
  for (int i = 0; i < nblock; i++) {
    int size = (add-- > 0) ? base_size + 1 : base_size;
    blocks.push_back(size);
  }
  return nblock;
}

vector<int> PartitionedEngine::compute_blocks(int total) {
  vector<int> blocks; compute_blocks(blocks, total);
  return blocks;
}

void PartitionedEngine::set_blocks() {
  block_mode = no_snps > block_size;
  if (!block_mode) {no_blocks = 1; return;}
  
  vector<int> blocks = compute_blocks(no_snps); 
  no_blocks = blocks.size();
  data->partition(blocks);  
}


void BlockEngine::set_covar() {
  delete covar; 
  try {covar = new BlockCorrelationData(engine->corr_data(), aggregate_type);} 
  catch (const DataException& de) {covar = 0;}
  if (!covar) throw StatException("preparing correlations", "block correlation data object failed to initialize");
  
  if (data_changed && block_mode) aggregator->set_corrs(covar->get_corrs(), true);      
}

void BlockEngine::set_perm(Permutation* source, int block_id, bool rewind) {
  if (has_perm && source) {
    if (block_id == 0 || !curr_perm) curr_perm = source->get_block(rewind);
    if (curr_perm) curr_perm->proc_mode(aggregate_type == Aggregator::Top);
  } else curr_perm = 0;
}

PermutationBlock* BlockEngine::process_perm(StatConverter* convert) {
  PermutationBlock* out = curr_perm; curr_perm = 0;
  if (convert) out->transform(convert);
  return out;
}

void BlockEngine::process_pval() {for (int j = 0; j < no_sub_tests; j++) pval(j) = aggregator->get_pval(j);}

double BlockEngine::run_core() {
  no_tests = no_sub_tests = engine->model_info(data).size();
  pval.set_size(1, no_tests);
  if (block_mode) run_blocks();
  else run_single();
  return pval(0);
}

void BlockEngine::run_single() {
  engine->run(data, needs_var, permutation);  
  no_param = engine->get_nparam();
  engine->load_pval(pval);
  
  load_param("eff_size", engine); load_param("eff_size_adjusted", engine);
  if (needs_var) set_covar();
}

void BlockEngine::run_blocks() {
  if (data_changed) update_dimension(no_blocks);
  aggregator->rewind(dimension, no_sub_tests);
  aggr_param.clear_param();

  int max_par = 0; 
  bool load_effsize = engine->has_param("eff_size");
  bool load_var = needs_var || data_changed;

  for (int i = 0; i < no_blocks; i++) { 
    data->set_block(i); set_perm(permutation, i, true);
    engine->run(data, load_var, curr_perm);

    if (engine->ntest() != no_sub_tests) throw GeneException("analyzing gene", "unexpected number of p-values in blocked gene model"); 
    for (int j = 0; j < no_sub_tests; j++) aggregator->add_pval(engine->get_pval(j));
    *no_param_aggr << engine->get_nparam(); 
    max_par = max(max_par, engine->get_nparam());

    if (load_effsize) {
      *eff_size_aggr << engine->get_param("eff_size");
      *eff_size_adjusted_aggr << engine->get_param("eff_size_adjusted");    
    }
  } 

  if (load_var) set_covar();
  process_pval();
  no_param = max( max_par, int(ceil(no_param_aggr->compute() * aggregator->get_scale())) ); 

  if (load_effsize) {                            
    double scale = aggregator->get_scale(true);
    eff_size_aggr->extract(get_param("eff_size"), scale);
    eff_size_adjusted_aggr->extract(get_param("eff_size_adjusted"), scale);    
  }

  if (has_perm) {
    PermutationEngine* add_perm = permutation->get_engine();
    if (add_perm) {
      add_perm->process_block(process_perm());
      while (add_perm->next()) {
        for (int i = 0; i < no_blocks; i++) {
          data->set_block(i); set_perm(add_perm, i, false);
          engine->run(data, false, curr_perm);
        }
        add_perm->process_block(process_perm());  
      } 
    } else process_perm(aggregator->get_converter());
  }
  data->set_block();
}

ModelInfoBlock MultiEngine::model_info(GeneStats* gs) {
  ModelInfoBlock out = (engine_bay[0]->model_info(gs)).set_model(ModelState::MultiModel);
  for (int i = 0; i < no_engine; i++) out.append((engine_bay[i]->model_info(gs)).shift_index(out.size()));
  out.detect_state(); return out;
}    

void MultiEngine::init_gene() {
  CompoundEngine::init_gene();
  no_tests *= no_engine + 1;
  aggr_param.clear_param();
  if (!has_perm) {
    if (!corr_perm) corr_perm = new PermutationEngine(no_cperm);
    corr_perm->set_data(data);
    permutation = corr_perm; 
  }    
  for (int i = 0; i< block_perm.size(); i++) block_perm[i]->set_empty();
  delete covar; covar = 0;    
}
 
void MultiEngine::set_covar() {
  delete covar; 
  try {covar = new MultiCorrelationData(aggregate_type);} 
  catch (const DataException& de) {covar = 0;}
  if (!covar) throw StatException("preparing correlations", "multi-model correlation data object failed to initialize");
  for (int i = 0; i < no_engine; i++) covar->set_data(engine_bay[i]->corr_data());
}

double MultiEngine::run_core() {
  if (no_engine == 0) {_LOG.error("initializing multi-model") << "no gene analysis models have been set" << endl; die();}
  if (!permutation || !permutation->get_engine()) {_LOG.error("initializing multi-model") << "multi-model requires PermutationEngine as input" << endl; die();}

  bool use_top = aggregate_type == Aggregator::Top;
  int max_par = 0; used_perm = 0; no_param = 0; 
  for (int i = 0; i < no_engine; i++) {
    PermutationBlock* perm_block = permutation->get_block(i==0); 
    if ( i < (no_engine - 1) || (has_perm && !permutation->is_empty()) ) data->store_stats(true); else data->store_stats(false);
    engine_bay[i]->set("output_block_perm", data_changed ? "true" : "false");
    engine_bay[i]->run(data, needs_var, perm_block); 

    if (i == 0) set_nsub(engine_bay[i]->ntest()); 
    else if (engine_bay[i]->ntest() != no_sub_tests) throw GeneException("analyzing gene", "unexpected number of p-values in gene multi-model"); 
    for (int j = 0; j < no_sub_tests; j++) {
      pval(j,i+1) = engine_bay[i]->get_pval(j);
      aggregator->add_pval(pval(j,i+1));
    }

    *no_param_aggr << engine_bay[i]->get_nparam(); 
    max_par = max(max_par, engine_bay[i]->get_nparam());  
    if (engine_bay[i]->has_param("snp_counts")) *snp_counts_aggr << engine_bay[i]->get_param("snp_counts");

    if (i == 0) {
      perm_buffer.set_size(perm_block->size(true), no_engine);
      used_perm = perm_block->size(false);        
    } else used_perm = min(used_perm, perm_block->size(false));
    perm_block->extract(perm_buffer[i], used_perm);

    if (data_changed && perm_block->get_attached()) store_block_perm(*perm_block->get_attached(), i);
    if (permutation->get_engine()) {perm_block->proc_mode(use_top); permutation->get_engine()->insert(perm_block);}
  }

  if (used_perm < MULTI_MIN_PERM) throw GeneException("analyzing gene", "number of successful permutations was too low to compute between-model correlations");
  if (data_changed) {
    snp_counts_aggr->extract(get_param("snp_counts"));
    aggregator->set_perm_corrs(perm_buffer, used_perm);

    for (int i = 0; i < (no_sub_tests - 1); i++) {
      if (i >= block_perm.size()) break;

      bool has_valid = false;
      for (int j = 0; j < no_engine; j++) {if (pval(i+1,j+1) >= 0 && pval(i+1,j+1) <= 1) has_valid = true; break;}
      if (!has_valid || !block_perm[i] || block_perm[i]->empty() || block_perm[i]->nrow() < MULTI_MIN_PERM) continue;

      aggregator->set_perm_corrs(*block_perm[i], -1, i+1);
    }
  }

  if (needs_var) {
    set_covar();
    covar->set_corrs(aggregator->get_corrs(), perm_buffer, used_perm);
  }

  for (int j = 0; j < no_sub_tests; j++) pval(j) = aggregator->get_pval(j);
  no_param = min(int(round(no_param_aggr->compute())), max_par); 

  if (has_perm) {   
    PermutationEngine* add_perm = permutation->get_engine();
    if (add_perm) {
      add_perm->process_block();
      while (add_perm->next()) {
        for (int i = 0; i < no_engine; i++) {
          PermutationBlock* perm_block = add_perm->get_block(); 
          if ( (i < (no_engine - 1) || !add_perm->is_empty()) ) data->store_stats(true); else data->store_stats(false);        
          engine_bay[i]->set("output_block_perm", "false");
          engine_bay[i]->run(data, false, perm_block);
          perm_block->proc_mode(use_top); add_perm->insert(perm_block);
        }
        add_perm->process_block();  
      } 
    }
  }
  return pval(0);
}

void MultiEngine::store_block_perm(Buffer<float>& perms, int index) {int npval = no_sub_tests - 1;
  if (perms.ncol() != npval) throw GeneException("analyzing gene", "number of p-values and number of permutation sets do not match");
  while (block_perm.size() < npval) block_perm.push_back(new Buffer<float>());

  for (int i = 0; i < npval; i++) {Buffer<float>& buffer = *block_perm[i];
    if (pval(i+1, index+1) < 0 || pval(i+1, index+1) > 1) buffer.set_empty();
    else if (index == 0) buffer.set_size(perms.nrow(), no_engine, true);
    if (buffer.empty()) continue;
    if (buffer.nrow() > perms.nrow()) buffer.shrink_rows(perms.nrow());
    buffer.assign_col(perms[i], index);
  }
}

void MultiEngine::add_engine(GeneModelEngine* engine, bool owner) {
  if (engine->perm_only()) {_LOG.error("initializing multi-model") << "cannot use permutation-only gene analysis model in multi-model" << endl; die();}

  engine_bay.push_back(engine);
  engine_owner.push_back(owner);    
  no_engine = engine_bay.size();
  update_dimension(no_engine);
}

bool MultiEngine::equals(GeneModelEngine* other) {
  if (!GeneModelEngine::equals(other)) return false;
  vector<GeneModelEngine*>& other_bay = static_cast<MultiEngine*>(other)->engine_bay;
  for (int i = 0; i < no_engine; i++) {if (!engine_bay[i]->equals(other_bay[i])) return false;}
  return true;    
}

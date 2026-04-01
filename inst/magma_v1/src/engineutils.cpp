/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "engineutils.h"

double AggregateParam::compute_base_core(double* data, int len) {
  switch (type) {
    case Sum:     return MathUtils::sum(data, len); 
    case Mean:    return MathUtils::get_mean(data, len); 
    case Min:     return MathUtils::get_min(data, len); 
    case Max:     return MathUtils::get_max(data, len); 
    case Product: return MathUtils::product(data, len); 
    default:      return 0;
  } 
}

double AggregateParam::compute_base(int index, double scale) {double& out = values(tot_elem,index); 
  if (!partitions.empty()) {
    int npart = partitions.size(), offset = 0;
    Buffer<double> buff(npart,1); double* data = values[index];
    for (int i = 0; i < npart; i++) {
      if (offset >= curr_elem) {npart = i; break;}
      int size = min(partitions[i], curr_elem-offset);
      if (size > 1) {
        buff(i) = compute_base_core(data+offset, size);       
      } else buff(i) = data[offset];
      offset += partitions[i];
    }
    out = compute_base_core(buff.begin(), npart);    
  } else out = compute_base_core(values[index], curr_elem);  
  out *= scale; return out;
}

double AggregateParam::compute(double scale, bool recompute) {if (recompute) computed = false;
  if (!computed && curr_elem > 0) {
    for (int i = 0; i < tot_values; i++) compute_base(i, scale);
    computed = true;
  }
  return computed ? values(tot_elem) : default_value;
}

void AggregateParam::extract(EngineParam* target, double scale, bool recompute) {compute(scale, recompute);
  MultiParam* mtarget = dynamic_cast<MultiParam*>(target);
  if (mtarget) mtarget->resize(tot_values);

  if (computed) {
    int count = min(tot_values, target->size()); 
    for (int i = 0; i < count; i++) target->set(values(tot_elem,i),i); 
  } else {
    for (int i = 0; i < target->size(); i++) target->set(default_value,i); 
  }
}

void Aggregator::check_size(Buffer<double>& corr, bool fail) {
  if (corr.ncol() != no_pval) {
    if (fail) throw GeneException("aggregating p-values", "block correlation matrix does not match number of gene blocks");
    else corr.set_empty();
  }
}

void Aggregator::rewind(int size_override, int ntest) {
  no_pval = size_override; no_tests = ntest; expected = no_pval * no_tests;
  scale[0] = scale[1] = 0; stats.set_size(no_tests,1,true);
  pval.set_size(no_tests,no_pval); pval_added = 0; 
  multi_status.set_size(no_tests, ms_Count);
  for (int i = 1; i < multi_pval.size(); i++) {delete multi_pval[i]; multi_pval[i] = 0;}
  if (multi_pval.size() < no_tests) multi_pval.resize(no_tests, 0);
  subset.clear();
  processed = pr_None;
}

void Aggregator::clear_corrs() {
  base_correlations.set_empty();
  for (int i = 1; i < multi_correlations.size(); i++) {if (multi_correlations[i]) multi_correlations[i]->set_empty();}
}

Buffer<double>& Aggregator::correlations(int index, bool create) {
  if (index > 0 && index < no_tests) {
    if (create) {
      if (index >= multi_correlations.size()) multi_correlations.resize(index+1,0);
      if (!multi_correlations[index]) multi_correlations[index] = new Buffer<double>();
      return *multi_correlations[index];
    } else if (index < multi_correlations.size() && multi_correlations[index] && !multi_correlations[index]->empty()) return *(multi_correlations[index]); 
  }
  return base_correlations;
}

Buffer<double>& Aggregator::sub_correlations(int index) {Buffer<double>& corrs = correlations(index);
  if (has_subset()) {int nsub = no_subset();
    corr_buffer.set_square(nsub);
    if (nsub > 1) {
      for (int i = 0; i < nsub; i++) {
        corr_buffer(i,i) = 1;
        for (int j = 0; j < i; j++) {
          corr_buffer(i,j) = corrs(subset[i],subset[j]);
          corr_buffer(j,i) = corr_buffer(i,j);
        }      
      }    
    } else corr_buffer(0) = 1;
    return corr_buffer;
  } else return corrs;
}

void Aggregator::set_corrs(Buffer<double>& cov, bool is_corr, int index) {
  if (index == 0) clear_corrs();
  Buffer<double>& target = correlations(index, true); target.assign(cov); 
  if (!is_corr) MathUtils::cov_to_cor(target);
  check_size(target, index==0); processed = pr_None;
}

int Aggregator::set_perm_corrs(Buffer<float>& perm, int use, int index) {
  if (index == 0) clear_corrs();
  Buffer<double>& target = correlations(index, true);
  if (use <= 0 || use > perm.nrow()) use = perm.nrow();

  try {MathUtils::cov_matrix(perm, target, true, use);}
  catch (const MathException& me) {
    if (index == 0) throw GeneException("analyzing gene", "unable to aggregate p-values", me.what());
    else target.set_empty();  
  } 

  check_size(target, index==0); processed = pr_None;
  return use;
}

void Aggregator::set_subset(vector<int>* sub) {
  if (sub) subset = *sub; else subset.clear();
  if (processed > pr_Checked) processed = pr_Checked;
}

double Aggregator::get_pval(int index) {
  if (index >= no_tests) return -1;
  if (processed != pr_Computed) process();
  if (multi_status(index, ms_Invalid)) return -1;

  StatConverter* conv = multi_status(index, ms_SepCorr) ? multi_pval[index] : to_pval;
  try {return conv->convert(get_stat(index));} 
  catch (const StatException& se) {throw GeneException("analyzing gene", "unable to convert test statistic to p-value", se.what());} 
}
  
double Aggregator::get_scale(bool root) {int type = root;
  if (scale[type] == 0) {
    if (root) {for (int i = 0; i < base_correlations.length(); i++) scale[type] += sqrt(base_correlations(i));}
    else scale[type] = MathUtils::sum(base_correlations.begin(), base_correlations.length());
    scale[type] = no_pval/scale[type];
  }
  return scale[type];  
}  

void Aggregator::process() {
  if (processed == pr_Computed) return;
  if (processed < pr_Checked) {
    if (no_pval == 0 || pval_added != expected) throw GeneException("aggregating p-values", "number of p-values does not match number of gene blocks");

    multi_status.set_size(no_tests, ms_Count, true); 
    for (int i_test = 0; i_test < no_tests; i_test++) {for (int i_pval = 0; i_pval < no_pval; i_pval++) {double& p = pval(i_test,i_pval);
      if ((p < 0 || p > 1)) {
        if (p != EU_NULL_PVAL) throw GeneException("aggregating p-values", "not all gene blocks have valid p-values");
        else multi_status(i_test, ms_Invalid) = true; 
      }
      if (&(correlations(i_test)) != &base_correlations && !multi_status(i_test, ms_Invalid)) multi_status(i_test, ms_SepCorr) = true;
    }}      
  }

  process_core(); processed = pr_Computed;
}

void MeanAggregator::process_core() {
  ConvertChisqSumDf2ToPval::set_brown(brown, sub_correlations(0));
  for (int i_test = 0; i_test < no_tests; i_test++) {
    if (multi_status(i_test, ms_Invalid)) continue;
    if (multi_status(i_test, ms_SepCorr)) {
      if (!multi_pval[i_test]) multi_pval[i_test] = new ConvertChisqSumDf2ToPval();
      ConvertChisqSumDf2ToPval::set_brown(*(static_cast<ConvertChisqSumDf2ToPval*>(multi_pval[i_test])), sub_correlations(i_test));      
    }

    double& curr = stats(i_test); curr = 0;
    if (has_subset()) {for (int i_sub = 0; i_sub < subset.size(); i_sub++) curr += logp.convert(pval(i_test,subset[i_sub]));}
    else {for (int i_pval = 0; i_pval < no_pval; i_pval++) curr += logp.convert(pval(i_test,i_pval));}
  }
}
  
void FittedAggregator::process_core() {
  compute_stats();
  
  prep_converter(0);
  for (int i_test = 1; i_test < no_tests; i_test++) {
    if (multi_status(i_test, ms_SepCorr)) prep_converter(i_test);
  }
}

ConvertPermutationToPval* FittedAggregator::gpd_converter(int index) {
  if (index > 0 && index < no_tests) {
    if (!multi_pval[index]) multi_pval[index] = gpd.clone();
    return static_cast<ConvertPermutationToPval*>(multi_pval[index]);
  }
  return &gpd;
}

void FittedAggregator::prep_converter(int index) {
  if (no_pval > 1) {
    fit_perm->set_obs(zstat[index]);
    fit_perm->set_block(1000, no_subset());
    fit_perm->set_input_covariance(sub_correlations(index));

    perm_buffer.set_size(max_perm,1); curr_perm = 0;
    while (fit_perm->next()) {
      fit_perm->get_sims(stat_buffer);
      int total = min(max_perm - curr_perm, (int) stat_buffer.ncol());      
      process_perms(stat_buffer, total);
      advance_perm(total);
    }    
  
    z_to_logp.convert(perm_buffer.begin(), curr_perm);
    gpd_converter(index)->fit_param(perm_buffer.begin(), curr_perm);
  } else gpd_converter(index)->toggle_noop();
}


void FittedAggregator::advance_perm(int added) {
  int offset = max(curr_perm - added, 0); added = curr_perm - offset;
  if (added > 0) fit_perm->process_bulk(perm_buffer.begin() + offset, added);
  if (fit_perm->is_empty()) return;

  int curr_sign = fit_perm->nsign();
  bool adap_complete = curr_sign >= adap_count && curr_perm >= min_perm;
  bool fit_complete = (curr_sign >= tail_thresh * curr_perm * tail_margin_scale) || curr_perm >= tail_min_perm;
  
  if (curr_perm < max_perm && !(adap_complete && fit_complete)) {
    int next_perm = block_size - curr_perm % block_size;
    if (curr_perm < block_size / 10) next_perm = block_size / 10 - curr_perm;
    else if (curr_perm < block_size / 2) next_perm = block_size / 2 - curr_perm;        
    fit_perm->set_block(min(next_perm, max_perm-curr_perm));
  } else fit_perm->set_empty();
}

ConvertPermutationToPval* TopAggregator::gpd_converter(int index) {BoundedMinP* minp;
  if (index > 0 && index < no_tests) {
    if (!multi_pval[index]) multi_pval[index] = bounded_gpd.clone();
    minp = static_cast<BoundedMinP*>(multi_pval[index]);
  } else minp = &bounded_gpd;
  minp->set_param("ntest", no_subset()); 
  return minp->get_internal();
}

void TopAggregator::compute_stats() {
  zstat.assign(no_tests, -999); zstat[0] = 0;
  for (int i_test = 0; i_test < no_tests; i_test++) {
    if (multi_status(i_test, ms_Invalid)) {stats(i_test) = 0; continue;}
    double minp = 1; 
    if (has_subset()) {for (int i_sub = 0; i_sub < subset.size(); i_sub++) {double &p = pval(i_test,subset[i_sub]); if (p < minp) minp = p;}}
    else {for (int i_pval = 0; i_pval < no_pval; i_pval++) {double &p = pval(i_test,i_pval); if (p < minp) minp = p;}}
    stats(i_test) = logp.convert(minp);

    double &global_z = zstat[multi_status(i_test, ms_SepCorr) ? i_test : 0], curr_z = to_z.convert(minp);
    if (curr_z > global_z) global_z = curr_z;    
  }
}

void TopAggregator::process_perms(Buffer<float>& perm, int total) {
  int nvar = perm.nrow(); float max;
  for (float *curr = perm.begin(), *end = curr+total*nvar; curr < end; curr += nvar) {max = *curr;
    for (int j = 1; j < nvar; j++) {if (curr[j] > max) max = curr[j];}    
    perm_buffer(curr_perm++) = max;
  }
}


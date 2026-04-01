/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genemodelengine_snpwise.h"

double SNPwise::run_core() {
  if (has_covar) {
    if (has_interactor) run_interaction();
    else run_covar();
  } else run_plain();

  for (int i = 0; i < pval.length(); i++) {
    try {pval(i) = (to_gene) ? (no_snps == 1 ? uni_pval(stat(i)) : to_gene->convert(stat(i))) : -1;}
    catch (const exception& e) {
      if (i == 0) throw;
      else pval(i) = -1;
    } 
  }

  store_converter();
  return pval(0);
}

void SNPwise::run_plain() {
  BufferWindow<float> wrap_var = data->get_var(); Buffer<float>& var = wrap_var.view();
  float df_e = sample_size - 2; from_snp_p.set_param("df2", df_e);  
  if (data->has_pval()) {
    snp_pval.assign(data->get_pval());
    stat(0) = stat_from_pval(snp_pval.begin());     
  } else { 
    snp_pval.set_size(no_snps,1);
    BufferWindow<float> wrap_xty = data->get_xty(); Buffer<float>& xty = wrap_xty.view();
    for (int i = 0; i < no_snps; i++) snp_pval(i) = xty(i)*xty(i) / var(i) / (sample_size-1);
    stat(0) = stat_from_ssm(snp_pval.begin(), 0, df_e);
  }

  if (has_perm) {                                                                   
    data->load_data(data_buffer, true);
    init_perm(data_buffer);
    while (permutation->next()) { 
      permutation->get_product(stat_buffer, data_buffer);
      for (int i_perm = 0; i_perm < stat_buffer.ncol(); i_perm++) { try {
        float* curr = stat_buffer[i_perm];
        for (int i_snp = 0; i_snp < no_snps; i_snp++) snp_pval(i_snp) = curr[i_snp]*curr[i_snp] / (sample_size-1);
        process_perm(stat_from_ssm(snp_pval[0], 0, df_e));        
      } catch (const exception& e) {permutation->failed();} }
      flush_perm();
    }
  }
}

void SNPwise::run_covar() {
  BufferWindow<float> wrap_var = data->get_var(); Buffer<float>& var = wrap_var.view();
  BufferWindow<float> wrap_xty = data->get_xty(); Buffer<float>& xty = wrap_xty.view();
  BufferWindow<float> wrap_xtz = data->get_xtz(); Buffer<float>& xtz = wrap_xtz.view();  

  float div = 1/float(sample_size - 1), df_e = sample_size - no_covar - 2; from_snp_p.set_param("df2", df_e);  
  float ssm_base = data->ytz_prod*div; 
  mat_buffer.set_size(no_snps,1); snp_pval.set_size(no_snps,1); Buffer<float>& ytz = data->ytz; int failed = 0;
  for (int i = 0; i < no_snps; i++) {
    mat_buffer(i) = (sample_size-1)*var(i) - MathUtils::in_product(xtz[i], no_covar)*div;  
    if (mat_buffer(i) >= EIGEN_CUTOFF) { 
      mat_buffer(i) = 1/mat_buffer(i);
      float diff = MathUtils::in_product(ytz[0], xtz[i], no_covar)*div - xty(i);
      snp_pval(i) = ssm_base + diff*diff * mat_buffer(i);
    } else {failed++; snp_pval(i) = null_ssm(0.45, ssm_base, df_e);}      
  }
  if (failed && failed/float(no_snps) > 0.25) throw GeneException("running snpwise analysis", "Too many SNPs in gene are fully collinear with covariates");   

  stat(0) = stat_from_ssm(snp_pval.begin(), ssm_base, df_e);  

  if (has_perm) {
    BufferWindow<float> covar_perm;    
    data->load_data(data_buffer);
    init_perm(data_buffer);
    while (permutation->next()) {
      permutation->get_product(stat_buffer, data_buffer);
      covar_perm = permutation->get_covar(); Buffer<float>& covar = covar_perm.view();
      for (int i_perm = 0; i_perm < stat_buffer.ncol(); i_perm++) { try {
        float* curr = stat_buffer[i_perm];
        ssm_base = MathUtils::in_product(covar[i_perm], no_covar)*div;
        for (int i_snp = 0; i_snp < no_snps; i_snp++) {
          float diff = MathUtils::in_product(covar[i_perm], xtz[i_snp], no_covar)*div - curr[i_snp];
          snp_pval(i_snp) = ssm_base + diff*diff * mat_buffer(i_snp);
        }
        process_perm(stat_from_ssm(snp_pval[0], ssm_base, df_e));              
      } catch (const exception& e) {permutation->failed();} }
      flush_perm();
    }
  }
}

/// assumes no_covar = 1, only the interactor
void SNPwise::run_interaction() {
  Buffer<float> var; var.assign(data->get_var());
  MathUtils::matrix_scale(var, sample_size - 1);
  Buffer<float>& pheno = data->get_pheno();  
  Buffer<float>& covar = data->get_covar();

  data->load_interaction(data_buffer); 
  data->load_pairwise_product(stat_buffer, data_buffer);

  mat_buffer.set_size(no_snps, 6);
  snp_pval.set_size(no_snps, 3);

  int failed = 0;
  float div = 1/float(sample_size - 1), df_e[2]; df_e[0] = sample_size - 3; df_e[1] = sample_size - 4;
  BufferWindow<float> wrap_xty = data->get_xty(); Buffer<float>& xty = wrap_xty.view();  
  BufferWindow<float> wrap_xtz = data->get_xtz(); Buffer<float>& xtz = wrap_xtz.view();    
  float &ytz = data->ytz(0), &ytz_sq = data->ytz_prod, ssm_base = ytz_sq*div; 
  for (int i = 0; i < no_snps; i++) {
    float det = var(i)*(sample_size-1) - xtz(i)*xtz(i);
    if (det < EIGEN_CUTOFF) {
      snp_pval(i,0) = null_ssm(0.45, ssm_base, df_e[0]); 
      snp_pval(i,1) = null_ssm(0.45, snp_pval(i,0), df_e[1]); 
      snp_pval(i,2) = null_ssm(0.69, ssm_base, df_e[1]/2); 
      failed++; continue;
    }

    mat_buffer(i,0) = var(i)/det; mat_buffer(i,1) = -xtz(i)/det; mat_buffer(i,2) = (sample_size-1)/det;    
    snp_pval(i,0) = ytz_sq*mat_buffer(i,0) + 2*mat_buffer(i,1)*ytz*xty(i) + mat_buffer(i,2)*xty(i)*xty(i);
    float itz = MathUtils::in_product(covar[0], data_buffer[i], sample_size);
    mat_buffer(i,3) = mat_buffer(i,0)*itz + mat_buffer(i,1)*stat_buffer(i);
    mat_buffer(i,4) = mat_buffer(i,1)*itz + mat_buffer(i,2)*stat_buffer(i);
    det = MathUtils::in_product(data_buffer[i], sample_size) - itz*mat_buffer(i,3) - stat_buffer(i)*mat_buffer(i,4);
   
    if (det < EIGEN_CUTOFF) {
      snp_pval(i,1) = null_ssm(0.45, snp_pval(i,0), df_e[1]); 
      snp_pval(i,2) = snp_pval(i,1);
      failed++; continue;
    }      
   
    mat_buffer(i,5) = 1/det; 
    double ssm_full = ytz*mat_buffer(i,3) + xty(i)*mat_buffer(i,4) - MathUtils::in_product(pheno[0], data_buffer[i], sample_size);
    snp_pval(i,1) = snp_pval(i,2) = ssm_full*ssm_full * mat_buffer(i,5) + snp_pval(i,0);
  }
  if (failed && failed/float(no_snps) > 0.25) throw GeneException("running snpwise analysis", "Too many SNPs/interaction terms in gene are fully collinear with covariates");     

  ///interaction
  from_snp_p.set_param("df1", 1); from_snp_p.set_param("df2", df_e[1]);   
  stat(1) = stat_from_ssm_diff(snp_pval[1], snp_pval[0], df_e[1]);  ///!invalidates snp_pval[1]

  ///main effect
  from_snp_p.set_param("df1", 1); from_snp_p.set_param("df2", df_e[0]);   
  stat(2) = stat_from_ssm(snp_pval[0], ssm_base, df_e[0]);  ///!invalidates snp_pval[0]
  
  ///full model 
  from_snp_p.set_param("df1", 2); from_snp_p.set_param("df2", df_e[1]);     
  stat(0) = stat_from_ssm(snp_pval[2], ssm_base, df_e[1]/2);  ///!invalidates snp_pval[2]

  if (has_perm) {
    BufferWindow<float> covar_perm;
    data->load_data(stat_buffer);
    data_buffer.append(stat_buffer);

    init_perm(data_buffer);
    while (permutation->next()) {
      permutation->get_product(stat_buffer, data_buffer);
      covar_perm = permutation->get_covar(); Buffer<float>& covar = covar_perm.view();
      for (int i_perm = 0; i_perm < stat_buffer.ncol(); i_perm++) { try {
        float *xty_curr = stat_buffer[i_perm]+no_snps, *ity_curr = stat_buffer[i_perm], &ytz_curr = covar(i_perm), ytz_sq_curr = ytz_curr*ytz_curr;
        ssm_base = MathUtils::in_product(covar[i_perm], no_covar)*div;
        for (int i_snp = 0; i_snp < no_snps; i_snp++) {
          double ssm = ytz_sq_curr*mat_buffer(i_snp,0) + 2*mat_buffer(i_snp,1)*ytz_curr*xty_curr[i_snp] + mat_buffer(i_snp,2)*xty_curr[i_snp]*xty_curr[i_snp];
          double ssm_full = ytz_curr*mat_buffer(i_snp,3) + xty_curr[i_snp]*mat_buffer(i_snp,4) - ity_curr[i_snp];
          snp_pval(i_snp) = ssm_full*ssm_full * mat_buffer(i_snp,5) + ssm;
        }
        process_perm(stat_from_ssm(snp_pval.begin(), ssm_base, df_e[1]/2));
      } catch (const exception& e) {permutation->failed();} }
      flush_perm();
    }
  }
}

double SNPwise::null_ssm(double F, double ssm0, float df_ratio) {
  return ((sample_size-1)*F/df_ratio + ssm0) / (F/df_ratio+1);
}

double SNPwise::stat_from_ssm_diff(double* ssm1, double* ssm0, float df_ratio) {double sst = sample_size - 1;
  for (int i = 0; i < no_snps; i++) ssm1[i] = (ssm1[i] - ssm0[i]) / (sst - ssm1[i]) * df_ratio;
  from_snp_p.convert(ssm1, no_snps);
  return stat_from_pval(ssm1);  
}

double SNPwise::stat_from_ssm(double* ssm1, double ssm0, float df_ratio) {double sst = sample_size - 1;
  if (ssm0 == 0) {for (int i = 0; i < no_snps; i++) ssm1[i] = ssm1[i] / (sst - ssm1[i]) * df_ratio;}
  else {for (int i = 0; i < no_snps; i++) ssm1[i] = (ssm1[i] - ssm0) / (sst - ssm1[i]) * df_ratio;}
  from_snp_p.convert(ssm1, no_snps);
  return stat_from_pval(ssm1);  
}

double SNPwise::stat_from_stats(double* stats) {
  from_snp_p.convert(stats, no_snps);
  return stat_from_pval(stats);
}

void SNPwise::init_gene() {
  CoreEngine::init_gene();
  stat.set_size(1, no_tests); 

  if (gene_changed) clear_converters();
  converter_ready = restore_converter();
  if (!converter_ready) to_gene = get_converter();

  if (needs_var && block_id == 0) {
    delete covar;
    covar = get_covar();
  }
} 

void SNPwise::store_converter() {
  if (!permutation || block_count <= 1 || !to_gene) return;
  if (stored_converters.size() > block_count) clear_converters();
  if (stored_converters.size() < block_count) stored_converters.resize(block_count, 0);

  if (stored_converters[block_id] != to_gene) {
    delete stored_converters[block_id];
    stored_converters[block_id] = to_gene->clone();
  }
}

bool SNPwise::restore_converter() {
  if (stored_converters.empty()) return false;
  if (stored_converters.size() != block_count) {clear_converters(); return false;} 
  if (stored_converters[block_id] == 0) return false;
  to_gene = stored_converters[block_id]; return true;
}

void SNPwise::clear_converters() {
  for (int i = 0; i < stored_converters.size(); i++) delete stored_converters[i];
  stored_converters.clear();  
}


double SNPwise::prep_block_corr(Buffer<float>& target, int block_i, int block_j) {
  if (block_i > block_j) swap(block_i, block_j);
  return data->load_block_cor(target, block_i, block_j, true);
}   

int SNPwise::step_size(int no_snps) {
  if (no_snps < 2*min_snps) return 1;
  return min(no_snps/min_snps, max_step);
}

double SNPwiseMeanBrown::run_core() {
  if (covar || !converter_ready) {
    BufferWindow<float> wrap_xtx = data->get_xtx(true); Buffer<float>& xtx = wrap_xtx.view();
    if (covar) compute_variance(covar, xtx);
    
    if (!converter_ready) {
      MathUtils::cov_to_cor(xtx); MathUtils::matrix_square(xtx);
      ConvertChisqSumDf2ToPval::set_brown(brown, xtx);  
    
      double fudge_level = 1;                                             
      if (no_snps > 1) { 
        double df = brown.get_df(), scale = brown.get_scale();
        double sdi_curr = 1 / sqrt(2 * df * scale * scale);
        
        double sdi_mono = sqrt(0.25) / no_snps, sdi_indep = sqrt(0.25/no_snps);
        double prop = MathUtils::clamp01((sdi_curr - sdi_mono) / sdi_indep);
        fudge_level = CHISQ_ASYM_FUDGE + pow(1-prop, pcorr_power)*pcorr_base;
        
        no_param = round(1 + prop * (no_snps-1)); 
      } else no_param = 1;

      fudge.set_param("fudge", fudge_level);
      if (has_interactor) no_param *= 2;
    }
  }
  
  return SNPwise::run_core();
}

double SNPwiseMeanBrown::stat_from_pval(double* pval) {
  to_snp_stat.convert(pval, no_snps);
  return MathUtils::sum(pval, no_snps);
}

double SNPwiseMeanImhof::run_core() {
  if (covar || !converter_ready) {
    BufferWindow<float> wrap_xtx = data->get_xtx(true); Buffer<float>& xtx = wrap_xtx.view();
    if (covar) compute_variance(covar, xtx);
    
    if (!converter_ready) {
      MathUtils::cov_to_cor(xtx); 
      int ncomp = imhof->set_corrs(xtx, eigen_thresh); 
      BaseBuffer<double>& eval = imhof->get_lambda();

      double sdi_curr = 0;
      for (int i = 0; i < ncomp; i++) sdi_curr += eval[i] * eval[i];
      sdi_curr = 1 / sqrt(2*sdi_curr);
      double sdi_mono = sqrt(0.5) / no_snps, sdi_indep = sqrt(0.5/no_snps);

      double prop = MathUtils::clamp01((sdi_curr - sdi_mono) / sdi_indep);
      no_param = round(1 + prop * (no_snps-1)); 
      if (has_interactor) no_param *= 2;
    }
  }
  
  return SNPwise::run_core();
}

double SNPwiseMeanImhof::stat_from_pval(double* pval) {
  to_snp_stat.convert(pval, no_snps);
  return MathUtils::sum(pval, no_snps);
}

double SNPwiseTop::run_core() {
  top_used = top_fraction > 0 ? max(1.0, double(round(no_snps*top_fraction))) : top_count;
  use_all = top_count >= no_snps;        
  log_stat = (top_used != 1);

  double rescale = xtx_max < no_snps ? double(no_snps) / xtx_max : 1;
  BufferWindow<float> wrap_xtx = data->get_capped_xtx(xtx_max, true); Buffer<float>& xtx = wrap_xtx.view();     

  if (covar) compute_variance(covar, xtx, rescale); 
  MathUtils::matrix_square(MathUtils::cov_to_cor(xtx));  
  no_param = round(eigen_count(xtx)*rescale*(1+has_interactor));  

  return SNPwise::run_core();
}

double SNPwiseTop::stat_from_pval(double* pval) {
  if (!log_stat) { ///-> top_used == 1
    double out = pval[0];
    for (double *end = pval+no_snps; pval < end; pval++) if (*pval < out) out = *pval;
    return out;    
  } else to_snp_stat.convert(pval, no_snps);

  if (use_all) return MathUtils::sum(pval, no_snps);
  
  for (int i = top_used; i < no_snps; i++) {
    double curr = pval[i]; int swap = i;
    for (int j = 0; j < top_used; j++) {
      if (curr > pval[j]) {curr = pval[j]; swap = j;}
    }
    if (swap != i) {pval[swap] = pval[i]; pval[i] = curr;}
  }
  
  double out = 0;
  for (int i = 0; i < top_used; i++) out += pval[i];
  return out;
}

double SNPwiseTopOne::run_core() {

  if (covar) compute_variance(covar, data->get_xtx().view());  
  double out = SNPwise::run_core();
  has_perm = external_perm;

  if (no_snps > 1) {  
    double pcorr = bounded_gpd.convert(to_snp_stat.convert(0.05));
    double est_param = log(1-pcorr) / log(0.95);
    no_param = MathUtils::clamp(round(est_param), 1, no_snps)*(1+has_interactor);
  } else no_param = 1+has_interactor;
  return out;
}

double SNPwiseTopOne::stat_from_pval(double* pval) {
  double out = pval[0];
  for (double *end = pval+no_snps; pval < end; pval++) if (*pval < out) out = *pval;
  return to_snp_stat.convert(out);  
}

void SNPwiseTopOne::init_gene() {
  SNPwise::init_gene();
  max_dim = min(sample_size, 1000*(1+has_interactor) + no_covar);

  external_perm = permutation;
  fit_distribution = !converter_ready && data_changed && no_snps > 1;
  proc_perm = external_perm && !fit_distribution;

  if (external_perm) external_perm->rewind();
  if (!fit_perm) {
    fit_perm = new SimulationMVN(max_perm, max_dim, 100, prune);
    perm_buffer.set_size(1,max_perm);
  } 

  if (fit_distribution) {
    fit_perm->set_data(data);
    if (!external_perm) permutation = fit_perm;
    has_perm = true;      
  } else if (no_snps == 1) {
    gpd.toggle_noop(true);
    bounded_gpd.set_param("ntest", 1);
  }
} 

void SNPwiseTopOne::init_perm(Buffer<float>& data_block, int index) {
  if (fit_distribution) {
    curr_perm = 0; curr_stat = index;
    fit_perm->set_obs(stat(index));
    fit_perm->set_block(1000, max_dim);

    if (data_block.ncol() > no_snps) {///SNP data assumed to be at end of buffer
      int no_var = data_block.ncol(), no_other = no_var - no_snps;
      BufferWindow<float> other_wrap = data_block.window(0,no_other); Buffer<float>& other = other_wrap.view();

      process_buffer.set_square(no_var);
      process_buffer.matrix_sub(0, no_other).noalias() = other.matrix().transpose() * other.matrix();
      process_buffer.matrix_sub(no_other, no_snps) = data->get_xtx().view().matrix();

      data->load_product(cross_buffer, other, ProductType::Raw, false);
      process_buffer.matrix_sub(no_other, no_snps, 0, no_other) = cross_buffer.matrix();

      fit_perm->set_input_covariance(data_block, process_buffer);      
    } else fit_perm->set_input_covariance(data_block, data->get_xtx().view());
  } else if (proc_perm) permutation->set_obs(stat(index), to_gene);     
}

void SNPwiseTopOne::process_perm(double perm_stat) {
  if (fit_distribution) {
    perm_buffer(curr_perm++) = perm_stat;
    fit_perm->process(perm_stat);
  } else if (proc_perm) permutation->process(perm_stat);
}  

void SNPwiseTopOne::flush_perm() {
  if (!fit_distribution) return;
  if (permutation == external_perm) external_cutoff = curr_perm;

  bool adap_complete = fit_perm->nsign() >= adap_count && curr_perm >= min_perm;
  bool fit_complete = (fit_perm->nsign() >= tail_thresh * curr_perm * tail_margin_scale) || curr_perm >= tail_min_perm;
  if (curr_perm < max_perm && !(adap_complete && fit_complete)) {
    if (permutation == external_perm && permutation->is_empty() && fit_perm) permutation = fit_perm;
    if (!permutation->is_empty()) {
      if (permutation == fit_perm) {
        int next_perm = block_size - curr_perm % block_size;
        if (curr_perm < block_size / 10) next_perm = block_size / 10 - curr_perm;
        else if (curr_perm < block_size / 2) next_perm = block_size / 2 - curr_perm;        
        fit_perm->set_block(min(next_perm, max_perm-curr_perm));
      }
      perm_buffer.total_cols(curr_perm + permutation->next_nperm());
      return;
    }                     
  }

  if (fit_perm) fit_perm->set_empty();
  fit_distribution = false;
  gpd.fit_param(perm_buffer.begin(), curr_perm, external_perm);
  bounded_gpd.set_param("ntest", no_snps);

  if (external_perm) {
    external_perm->set_obs(stat(curr_stat), to_gene, 1, false);      
    for (int i = 0; i < external_cutoff; i++) external_perm->process(perm_buffer(i));
    proc_perm = true;                        
  }
}





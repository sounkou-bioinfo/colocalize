/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genemodelengine_pcreg.h"
                                                 
#include <Eigen/Eigen>


double LinearRegression::run_core() {
  if (no_pheno > 1) {
    if (has_perm) throw GeneException("running regression", "invalid state (has_perm) for analysis with multiple phenotypes");
    if (has_covar) throw GeneException("running regression", "invalid state (has_covar) for analysis with multiple phenotypes");  
  }

  if (has_covar) {
    if (has_interactor) run_interaction();
    else run_covar();
  } else run_plain();

  if (covar) {
    SubsetStats* subset = dynamic_cast<SubsetStats*>(data);
    if (subset) compute_variance(covar, subset->get_wtw().view());    
    else set_variance(covar, data->npcs());
  }
  return pval(0);
}

void LinearRegression::run_plain() {
  BufferWindow<float> wrap_xtx = data->get_xtx(true); Buffer<float>& xtx = wrap_xtx.view();
  no_param = prune_invert(xtx, stat_buffer, true); data->prep_pcs(xtx);

  float df_m = no_param, df_e = sample_size - no_param - 1, scale = df_e / df_m, ssm, sse;
  ssm = MathUtils::in_product(data->get_wty().view().begin(), no_param) / (sample_size - 1);
  sse = sample_size - ssm - 1;

  double stat = (ssm / sse) * scale; 
  converter.set_param("df1", df_m); converter.set_param("df2", df_e);
  pval(0) = converter.convert(stat);

  eff_size->element() = MathUtils::clamp01(ssm / (sample_size - 1));
  eff_size_adjusted->element() = MathUtils::clamp01(1 - sse / df_e);

  if (no_pheno > 1) {
    BufferWindow<float> wrap_wty = data->get_wty_multi(); Buffer<float>& wty = wrap_wty.view();    
 
    for (int p = 0; p < (no_pheno-1); p++) {
      ssm = 0;
      for (int i = 0; i < no_param; i++) {float& v = wty(p,i); ssm += v*v;}
      ssm /= (sample_size - 1); sse = sample_size - ssm - 1;
      pval(p+1) = converter.convert((ssm / sse) * scale);      
    }
  } 
  
  if (has_perm) {
    permutation->set_obs(stat, &converter);
    BufferWindow<float> geno = data->get_pcs();
    while (permutation->next()) {
      permutation->get_product(stat_buffer, geno.view());
      for (int i = 0; i < stat_buffer.ncol(); i++) { try {
        ssm = MathUtils::in_product(stat_buffer[i], no_param) / (sample_size - 1);
        stat = ssm / (sample_size - ssm - 1) * scale;        
        permutation->process(stat);      
      } catch (const exception& e) {permutation->failed();} }
    }
  } 
}

void LinearRegression::run_covar() {      
  BufferWindow<float> wrap_xtx = data->get_xtx(true); Buffer<float>& xtx = wrap_xtx.view();
  int eff_snps = prune_invert(xtx, stat_buffer, true); data->prep_pcs(xtx);
  no_param = covar_invert(mat_buffer, true); ///decomposes XtX, not covariance  
  
  data_buffer.set_size(eff_snps+no_covar, 1);
  data_buffer.window(0,eff_snps,true).view().assign(data->get_wty());
  data_buffer.window(eff_snps,no_covar,true).view().assign(data->ytz);
  MathUtils::matrix_prod(data_buffer, mat_buffer, stat_buffer, MathUtils::TransposeFirst);

  float scale = 1/float(sample_size - 1), df_m, df_e;
  float ssm_base = data->ytz_prod*scale;  
  float ssm = MathUtils::in_product(stat_buffer.begin(), stat_buffer.length());

  pval(0) = compute_pval(ssm, ssm_base, no_param, df_m, df_e);
  covar_rsq(ssm, ssm_base, df_e, sample_size - no_covar - 1, 0);

  if (has_perm) {
    float scale[2]; scale[0] = 1/float(sample_size-1); scale[1] = df_e / df_m;  
    int no_var = mat_buffer.ncol();

    BufferWindow<float> geno = data->get_pcs(); 
    Buffer<float>& covar = data->get_covar();
  
    data_buffer.set_size(sample_size, no_var);
    data_buffer.matrix().noalias() = geno.view().matrix() * mat_buffer.matrix_rowskip(0,eff_snps) + covar.matrix() * mat_buffer.matrix_rowskip(eff_snps,no_covar); 
    MathUtils::matrix_prod(data_buffer, data->get_pheno(), stat_buffer, MathUtils::TransposeFirst);

    ssm = MathUtils::in_product(stat_buffer.begin(), no_var);
    double stat = (ssm - ssm_base) / (sample_size - ssm - 1) * scale[1]; 
    permutation->set_obs(stat, &converter);
  
    BufferWindow<float> covar_perm;
    while (permutation->next()) { 
      permutation->get_product(stat_buffer, data_buffer);
      covar_perm = permutation->get_covar();
      for (int i = 0; i < stat_buffer.ncol(); i++) { try {
        ssm_base = MathUtils::in_product(covar_perm.view()[i], no_covar) * scale[0];
        ssm = MathUtils::in_product(stat_buffer[i], no_var);
        stat = (ssm - ssm_base) / (sample_size - ssm - 1) * scale[1]; 
        permutation->process(stat); 
      } catch (const exception& e) {permutation->failed();} }
    }
  }
}

void LinearRegression::run_interaction() {
  BufferWindow<float> wrap_xtx = data->get_xtx(true); Buffer<float>& xtx = wrap_xtx.view(); Buffer<float> effects;
  int eff_snps = prune_invert(xtx, stat_buffer, true); data->prep_pcs(xtx);
  int main_param = covar_invert(mat_buffer, true); ///decomposes XtX, not covariance  

  data_buffer.set_size(eff_snps+no_covar, 1);
  data_buffer.window(0,eff_snps,true).view().assign(data->get_wty());
  data_buffer.window(eff_snps,no_covar,true).view().assign(data->ytz);
  MathUtils::matrix_prod(data_buffer, mat_buffer, effects, MathUtils::TransposeFirst);

  float scale = 1/float(sample_size - 1);
  float ssm_base = data->ytz_prod*scale;  
  float ssm = MathUtils::in_product(effects.begin(), effects.length());

  data->load_interaction(data_buffer, ProcessType::Centered);
  MathUtils::matrix_prod(data_buffer, mat_buffer, MathUtils::TransposeFirst);
  int int_snps = prune_invert(mat_buffer, stat_buffer, true);  
  MathUtils::matrix_prod(data_buffer, mat_buffer, stat_buffer); ///stat_buffer => project interaction component
  MathUtils::normalize(stat_buffer);
  no_param = interaction_invert(mat_buffer, stat_buffer, true);

  data_buffer.set_size(eff_snps+int_snps+no_covar, 1); 
  data_buffer.window(0,eff_snps,true).view().assign(data->get_wty());
  data_buffer.matrix_rowskip(eff_snps, int_snps).noalias() = stat_buffer.matrix().transpose() * data->get_pheno().matrix();
  data_buffer.window(eff_snps+int_snps,no_covar,true).view().assign(data->ytz);
  MathUtils::matrix_prod(data_buffer, mat_buffer, effects, MathUtils::TransposeFirst);

  float ssm_full = MathUtils::in_product(effects.begin(), effects.length());

  float df0 = no_param - main_param, df_m, df_e;
  try {pval(2) = compute_pval(ssm, ssm_base, main_param, df_m, df_e);} catch (const exception& e) {pval(2) = -1;} ///main effect
  try {pval(1) = compute_pval(ssm_full, ssm, df0, df_m, df_e);} catch (const exception& e) {pval(1) = -1;} ///interaction
  pval(0) = compute_pval(ssm_full, ssm_base, no_param, df_m, df_e); ///full model
  covar_rsq(ssm_full, ssm_base, df_e, df_e + no_param, 0);

  if (has_perm) {
    float scale[2]; scale[0] = 1/float(sample_size-1); scale[1] = df_e / df_m;  
    int no_var = mat_buffer.ncol();

    BufferWindow<float> geno = data->get_pcs(); 
    Buffer<float>& covar = data->get_covar();
      
    data_buffer.set_size(sample_size, no_var);
    data_buffer.matrix().noalias() = geno.view().matrix() * mat_buffer.matrix_rowskip(0,eff_snps)
                                      + stat_buffer.matrix() * mat_buffer.matrix_rowskip(eff_snps,int_snps)
                                      + covar.matrix() * mat_buffer.matrix_rowskip(eff_snps+int_snps,no_covar); 
    MathUtils::matrix_prod(data_buffer, data->get_pheno(), effects, MathUtils::TransposeFirst);
    ssm_full = MathUtils::in_product(effects.begin(), no_var);

    double stat = (ssm_full - ssm_base) / (sample_size - ssm_full - 1) * scale[1]; 
    permutation->set_obs(stat, &converter);  
    BufferWindow<float> covar_perm;                     
    while (permutation->next()) { 
      permutation->get_product(stat_buffer, data_buffer);
      covar_perm = permutation->get_covar(); Buffer<float>& covar = covar_perm.view();
      for (int i = 0; i < stat_buffer.ncol(); i++) { try {
        ssm_full = MathUtils::in_product(stat_buffer[i], no_var);
        ssm_base = MathUtils::in_product(covar[i], no_covar) * scale[0];
        stat = (ssm_full - ssm_base) / (sample_size - ssm_full - 1) * scale[1]; 
        permutation->process(stat); 
      } catch (const exception& e) {permutation->failed();} }
    } 
  }     
}

void LinearRegression::covar_rsq(float ssm1, float ssm0, float df_e1, float df_e0, int index) {
  float rsq0 = ssm0 / df_e0, rsq1 = ssm1 / df_e1, par1 = df_e0 - df_e1; 
  float rsq1_adj = rsq1 - (1-rsq1) * par1 / df_e1;

  eff_size->element(index) = MathUtils::clamp01((rsq1 - rsq0) / (1 - rsq0));
  eff_size_adjusted->element(index)= MathUtils::clamp01((rsq1_adj - rsq0) / (1 - rsq0));
}

double LinearRegression::compute_pval(float ssm_model, float ssm_null, int npar, float& df_m, float& df_e, double& stat) {
  df_m = npar; df_e = sample_size - npar - no_covar - 1;
  stat = (ssm_model - ssm_null) / float(sample_size - ssm_model - 1) * df_e / df_m; 

  converter.set_param("df1", df_m); converter.set_param("df2", df_e);
  return converter.convert(stat);
}

int LinearRegression::covar_invert(Buffer<float>& target, bool root) {
  BufferWindow<float> wrap_wtz = data->get_wtz(); Buffer<float>& wtz = wrap_wtz.view(); 
  int eff_snps = wtz.ncol(), nvar = eff_snps + no_covar; target.set_square(nvar, true);
 
  for (int i = 0; i < nvar; i++) target(i,i) = sample_size - 1;
  for (int i = 0; i < eff_snps; i++) memcpy(target[i]+eff_snps, wtz[i], sizeof(float)*no_covar);

  int dim = MathUtils::invert_matrix(target, EIGEN_CUTOFF, root);
  if (dim <= no_covar) throw GeneException("running regression", "SNPs in gene are fully collinear with covariates");

  return dim - no_covar;
}

int LinearRegression::interaction_invert(Buffer<float>& target, Buffer<float>& int_data, bool root) {
  BufferWindow<float> wrap_wtz = data->get_wtz(); Buffer<float>& wtz = wrap_wtz.view(); 
  int eff_snps = wtz.ncol(), int_snps = int_data.ncol(), nvar = eff_snps + int_snps + no_covar; target.set_square(nvar, true);

  target.matrix_sub(eff_snps, int_snps, 0, eff_snps).noalias() = int_data.matrix().transpose() * data->get_pcs().view().matrix();
  target.matrix_sub(eff_snps+int_snps, no_covar, eff_snps, int_snps).noalias() = data->get_covar().matrix().transpose() * int_data.matrix();
  for (int i = 0; i < eff_snps; i++) memcpy(target[i]+eff_snps+int_snps, wtz[i], sizeof(float)*no_covar);
  for (int i = 0; i < nvar; i++) target(i,i) = sample_size - 1;

  int dim = MathUtils::invert_matrix(target, EIGEN_CUTOFF, root);
  if (dim <= no_covar) throw GeneException("running regression", "SNPs and interaction terms in gene are fully collinear with covariates");
  return dim - no_covar;
}


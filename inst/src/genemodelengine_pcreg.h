/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEMODELENGINE_PCREG_H
#define GENEMODELENGINE_PCREG_H

#include "genemodelengine_base.h"

class LinearRegression : public PrinCompEngine {
  ConvertFToPval converter;

  MultiParam* eff_size;
  MultiParam* eff_size_adjusted;  

  int prune_count;
  float prune_perc;

  int no_pheno;

  int covar_invert(Buffer<float>& target, bool root=false);
  int interaction_invert(Buffer<float>& target, Buffer<float>& int_data, bool root=false);

  double run_core();
  void run_plain();
  void run_covar();
  void run_interaction();
  
  double compute_pval(float ssm_model, float ssm_null, int npar) {float df_m, df_e; double stat; return compute_pval(ssm_model, ssm_null, npar, df_m, df_e, stat);}
  double compute_pval(float ssm_model, float ssm_null, int npar, float& df_m, float& df_e) {double stat; return compute_pval(ssm_model, ssm_null, npar, df_m, df_e, stat);}  
  double compute_pval(float ssm_model, float ssm_null, int npar, double& stat) {float df_m, df_e; return compute_pval(ssm_model, ssm_null, npar, df_m, df_e, stat);}
  double compute_pval(float ssm_model, float ssm_null, int npar, float& df_m, float& df_e, double& stat);  
  void covar_rsq(float ssm1, float ssm0, float df_e1, float df_e0, int index);  

  virtual void init_gene() {
    PrinCompEngine::init_gene();
    if (no_pheno != data->npheno()) throw GeneException("running regression", "unexpected number of phenotypes");
    if (no_pheno > 1) {
      if (no_tests > 1) throw GeneException("running regression", "invalid state (no_tests) for analysis with multiple phenotypes");
      no_tests = no_pheno;
      pval.set_size(1,no_tests); 
    }
    eff_size->resize(no_tests); eff_size_adjusted->resize(no_tests);    
    int prune_keep = prune_perc*no_snps; if (prune_count > 0) prune_keep = min(prune_keep, prune_count);
    prune_drop = no_snps - prune_keep;    
    if (needs_var && block_id == 0) covar = new PCCorrelationData(data, pcid);
  }

public:
  LinearRegression(Settings& settings, int no_pheno=1) : PrinCompEngine(settings), no_pheno(no_pheno) {
    prune_var_prop = settings.getn("prune_cutoff", 0.99);
    prune_min_val = settings.getn("prune_value", 0);
    prune_count = settings.geti("prune_count", 0);
    prune_perc = settings.getn("prune_prop", 1);

    eff_size = register_param("eff_size", new MultiParam(1));    
    eff_size_adjusted = register_param("eff_size_adjusted", new MultiParam(1));        
  }
  virtual ~LinearRegression(){}
  
  void clear() {PrinCompEngine::clear();}
  ModelInfoBlock model_info(GeneStats* gs) {
    ModelInfoBlock out = make_info(gs, ModelState::PrinCompReg);
    out.no_repeat = no_pheno;
    return out;
  }
};


#endif /** GENEMODELENGINE_PCREG_H */

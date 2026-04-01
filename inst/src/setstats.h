/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SETSTATS_H
#define SETSTATS_H

#include "data_genedata.h"
#include "setstats_utils.h"

#define DETERMINANT_THRESH 0.025
#define PC_PRUNE_VARS 1e-5

class SetVariables {
  friend class GeneVariableMap;
  friend class SelfContainedComputation;
  friend class RegressionComputation;
  
  double gc_correction;
  double pc_prune;

  GeneData& gene_info;
  GeneCovarData<double> outcome_info;
  vector<GeneSetData*> set_info;
  vector<GeneCovarData<float>*> covar_info;
 
  SetStatsUtils::GeneMatrix* inv_corrs;

  int mean_N;

  void init_internal();
  template<typename VAR> void set_internal(VAR& source, GeneCovarData<float>& var_info, SetStatsUtils::InternalVariable::Type var_id);
  
  void register_block(VariableIndex* index, int block, SetStatsUtils::VariableType type);
  void register_interactions(VariableIndex* index, int block, SetStatsUtils::InteractionType type);
  
  void prep_sets(GeneSetData& set_block);
  template<typename T> void prep_covar(GeneCovarData<T>& covar_block, double trunc_low, double trunc_high, bool zstat=false);
 
public:
  GeneVariableMap variables;
  int internal_block;
 
  SetVariables(GeneData& info, double gc_correction, double pc_prune) : gc_correction(gc_correction), pc_prune(pc_prune), gene_info(info), outcome_info("OUTCOME_ZSTAT"), inv_corrs(0), mean_N(0), variables(*this), internal_block(-1) {DataFrame::set_owner(&gene_info, this); init_internal();}
  ~SetVariables() {
    DataFrame::delete_owned(&gene_info, this);
    for (int i = 0; i < set_info.size(); i++) delete set_info[i];
    for (int i = 0; i < covar_info.size(); i++) delete covar_info[i];    
    delete inv_corrs;
  }

  int add_sets(GeneSetData* sets) {set_info.push_back(sets); return set_info.size() - 1;}
  int add_covar(GeneCovarData<float>* covar) {covar_info.push_back(covar); return covar_info.size() - 1;}

  int add_interactions(GeneSetData* ss);
  int add_interactions(GeneCovarData<float>* sc);

  bool lock_variables();  
  pair<int,int> process_missing(bool drop_genes, const string& impute_type, float impute_value=0);

  void prep_outcome(double trunc_low, double trunc_high);
  void prep_sets() {for (int s = 0; s < set_info.size(); s++) prep_sets(*set_info[s]);}
  void prep_covar(double trunc_sd) {for (int c = 0; c < covar_info.size(); c++) prep_covar(*covar_info[c], trunc_sd, trunc_sd);}

  void init_corrs();
     
  int no_set_blocks() {return set_info.size();}
  int no_covar_blocks() {return covar_info.size();}

  GeneData& get_genes() {return gene_info;}
  GeneSetData& get_sets(int block) {return *set_info[block];}
  GeneCovarData<float>& get_covar(int block) {return *covar_info[block];}

  GeneCovarData<double>& get_outcome() {return outcome_info;}  
  double* get_outcome_data(int index=0); 

  SetStatsUtils::VariableInfo& get_variable(int var_id, bool check_data=false); 

  int no_genes() {return gene_info.data_used(false);}
  int mean_sampsize();
  
  SetStatsUtils::GeneMatrix& get_corrs(bool verbose=false); 
};


class SelfContainedComputation {
  SetVariables& variable_data;
  ConvertNormToPval converter;
  BaseBuffer<double> buffer;

  void error(const string& msg) {error_msg.clear(); throw SetException("using stat computation", msg);}
  void failure(const string& msg) {error_msg = msg; throw SetException("using stat computation", msg);}

public:
  string error_msg;  
  SetStatsUtils::VariableInfo* variable; 
  
  double mu;
  double pval;
  int no_genes;
  
  SelfContainedComputation(SetVariables& vars) : variable_data(vars), variable(&variable_data.get_variable(-1)) {}
  
  void set_variable(int var_id);
  void run();
};


class RegressionComputation {
public:
  enum ParameterType {pt_Any, pt_Core, pt_PruneCore, pt_Intercept, pt_Covariate, pt_InternalCovariate, pt_HiddenCovariate, pt_UsedCovariate, pt_Variable, pt_MainVariable, pt_COUNT};
  enum ModelType {mt_Intercept, mt_PruneModel, mt_BaseModel, mt_FullModel, mt_COUNT};
private:
  SetVariables& variable_data;
  bool initialised, resid_core, allow_collinear;
    
  int no_param;
  int no_rows;
  
  int eff_size;
  int eff_core;
  
  SetStatsUtils::CollinearType collinear_core;
  SetStatsUtils::CollinearType collinear_mode;
 
  ConvertTToPval converter;
  SetStatsUtils::BlockInverter inverter;

  enum BufferType {bt_XtX, bt_XtXinvStore, bt_XtXinvUse, bt_XtZ, bt_Coefficients, bt_StdErrors, bt_Temporary, bt_COUNT};
  enum ComputedType {ct_Inverse, ct_Coefficients, ct_Variance, ct_StdErrors, ct_Fitted, ct_COUNT};

  vector<short> computed;
  
  vector<Buffer<double>*> buffers;  
  vector<BaseBuffer<double>*> residuals;

  vector<double> ssr; 
  vector<int> df;
  double rsq_multiplier;
  
  Buffer<double> predictors; ///projected
  BaseBuffer<double> outcome; ///projected
  BaseBuffer<double> zstat; ///source
  
  vector<SetStatsUtils::VariableInfo*> variables;
  vector<ParameterType> variable_types;
  vector<pair<int,int> > type_range; 
  
  bool check_init(bool has_init);
  int convert_index(int& param_index, ParameterType type);
  pair<int,int> range_sum(pair<int,int> A, pair<int,int> B);

  void init_model(int no_vars, bool prune_only);
  void init_params(int no_vars);
  void init_buffers();
  void set_degrees();
  
  void error(const string& msg, const string& sub_msg="");
  void failure(const string& msg, const string& sub_msg="");
  
  double product(double* var, int len); 
  double product(float* var_A, double* var_B, int len); 
  double product(double* var_A, double* var_B, int len);   
  double product(set<int>& var_A, double* var_B); 
  void variable_sum(BaseBuffer<double>& target, SetStatsUtils::VariableInfo& var, double scale);

  void invert_conditional(bool prune_only);
  void residualize();
  void compute_core(int offset);

  Buffer<double>& fit_coefficients(ModelType type) {return fit_coefficients(type, *buffers[type == mt_FullModel ? bt_Coefficients : bt_Temporary]);}
  Buffer<double>& fit_coefficients(ModelType type, Buffer<double>& target);
  void fit_model(ModelType type);
  void fit_stderrors();  
  
  Buffer<double>& xtx() {return *buffers[bt_XtX];}
  Buffer<double>& xtx_inv();
  Buffer<double>& xtz() {return *buffers[bt_XtZ];}

  double scaling(int index);

  Buffer<double>& coefficients();
  Buffer<double>& std_errors();
  double residual_variance(ModelType type=mt_FullModel, bool scale=true);
  BaseBuffer<double>& resid_values(ModelType type=mt_FullModel);  
  
public:
  string error_msg;
  
  RegressionComputation(SetVariables& vars, bool residualize, bool allow_collinear) : variable_data(vars), initialised(false), resid_core(residualize), allow_collinear(allow_collinear), no_rows(0), inverter(DETERMINANT_THRESH, PC_PRUNE_VARS, allow_collinear), rsq_multiplier(1) {variables.push_back(&variable_data.get_variable(-1)); variable_types.push_back(pt_Intercept);}
  ~RegressionComputation() {Utils::delete_vector(buffers); Utils::delete_vector(residuals);}

  void init(int nvar=1);
  
  void set_conditional(set<int>& ids, ParameterType type);
  void set_variable(int var_id);
  void set_variables(set<int>* var_ids);
  int update_nvar(int no_vars, bool is_init=false);
  
  void run();  

  int get_size(bool total=true) {return total ? no_rows : eff_size;}
  int get_df() {return df[mt_FullModel];};

  SetStatsUtils::VariableInfo& get_variable(int index, ParameterType type=pt_Variable) {return *variables[convert_index(index, type)];}
  vector<SetStatsUtils::VariableInfo*> get_variables(ParameterType type=pt_Variable);

  int get_ngenes(int index, ParameterType type=pt_Variable); 
  double get_coeff(int index, ParameterType type=pt_Variable);
  double get_coeff_std(int index, ParameterType type=pt_Variable);  
  double get_se(int index, ParameterType type=pt_Variable);
  double get_pval(int index, short direction, ParameterType type=pt_Variable);
	bool get_collinear(int index, ParameterType type=pt_Variable);
	bool has_type(int index, ParameterType type, ParameterType select_type=pt_Variable);

  double get_rsq(bool base_rsq=false, bool rescale=true);
  void load_fitted(Buffer<double>& target, ModelType type);
  BaseBuffer<double>& get_outcome() {return zstat;}

  int param_offset(ParameterType type) {return type_range[type].first;}  
  int param_count(ParameterType type) {return type_range[type].second;}
};

#include "setstats.tpp"


#endif /** SETSTATS_H */

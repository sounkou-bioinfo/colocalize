/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEMODELENGINE_SNPWISE_H                                
#define GENEMODELENGINE_SNPWISE_H

#include "genemodelengine_base.h"

class SNPwise : public CoreEngine {
protected:
  DefaultStats* data;
  RawCorrelationData* covar;

  ConvertFToPval from_snp_p;
  StatConverter& to_snp_stat;
  StatConverter* to_gene;
  vector<StatConverter*> stored_converters;
  bool converter_ready;
  
  Buffer<double> snp_pval;
  Buffer<double> stat;
  
  int min_snps;
  int max_step;
  double pcorr_base;
  double pcorr_power;

  virtual double run_core();
  void run_plain();
  void run_covar();
  void run_interaction();  

  int step_size(int no_snps);
  
  double null_ssm(double F, double ssm0, float df_ratio);
  
  double stat_from_ssm_diff(double* ssm1, double* ssm0, float df_ratio);
  double stat_from_ssm(double* ssm1, double ssm0, float df_ratio);
  double stat_from_stats(double* stats);
  virtual double stat_from_pval(double* pval) = 0;
  virtual double uni_pval(double curr_stat) {return exp(-0.5*curr_stat);}

  virtual void init_perm(Buffer<float>& data_block, int index=0) {permutation->set_obs(stat(index), to_gene); UNUSED(data_block);}
  virtual void process_perm(double perm_stat) {permutation->process(perm_stat);}
  virtual void flush_perm() {}

  virtual double prep_block_corr(Buffer<float>& buffer, int block_i, int block_j);

  virtual GeneStats* get_data() {return data;}
  virtual void set_data(GeneStats* input) {data = input->convert<DefaultStats>();}
  virtual RawCorrelationData* get_covar() {return new RawCorrelationData(data, min_snps, max_step);}

  virtual void init_gene();
  
  virtual StatConverter* get_converter() = 0;
  void store_converter();
  bool restore_converter();
  void clear_converters();

public:
  SNPwise(Settings& settings, StatConverter& to_stat) : CoreEngine(settings), to_snp_stat(to_stat), covar(0), to_gene(0) {
    from_snp_p.set_param("df1", 1);
    min_snps = settings.geti("corr_min_snps");
    max_step = settings.getn("corr_max_step"); 
    pcorr_base = settings.getn("gene_snpwise_pcorr", 0) / 100;     
    pcorr_power = settings.getn("gene_snpwise_pcorr_power", 4);
  }
  virtual ~SNPwise() {clear_converters(); delete covar;}

  virtual CorrelationData* corr_data() {
    CoreCorrelationData* out = covar; covar = 0;
    return out;    
  }
  void clear() {
    CoreEngine::clear();
    delete covar; covar = 0;
    clear_converters();
  }  
};


class SNPwiseMeanBrown : public SNPwise {
  ConvertChisqSumDf2ToPval brown;
  FudgePval fudge; 
  ConvertPvalToLogPval to_logp;

  double run_core();
  double stat_from_pval(double* pval);
    
  virtual StatConverter* get_converter() {return &fudge;}
public:
  SNPwiseMeanBrown(Settings& settings) : SNPwise(settings, to_logp), fudge(CHISQ_ASYM_FUDGE) {fudge.set_internal(&brown);}
  virtual ~SNPwiseMeanBrown(){}

  void clear() {SNPwise::clear();}
  ModelInfoBlock model_info(GeneStats* gs) {return make_info(gs, ModelState::SNPwiseModelMeanBrown);}
};

class SNPwiseMeanImhof : public SNPwise {
  ConvertChisqSumDf1ToPval* imhof; 
  ConvertPvalToChisq to_chisq;
  ConvertChisqToPval from_chisq;

  double eigen_thresh;

  double run_core();
  double stat_from_pval(double* pval);
  double uni_pval(double curr_stat) {return from_chisq.convert(curr_stat);}  

  virtual void init_perm(Buffer<float>& data_block, int index=0) {permutation->set_obs(stat(index)); UNUSED(data_block);}      
  virtual StatConverter* get_converter() {return imhof;}
public:
  SNPwiseMeanImhof(Settings& settings) : SNPwise(settings, to_chisq), eigen_thresh(settings.getn("gene_imhof_eigen_thresh")) {
    string subtype = settings.gets("gene_snpwise_subtype", "");
    int max_comp = settings.geti("gene_chunked_size");
    
    if (subtype == "") {
      ConvertChisqSumDf1ToPvalImhof* converter = new ConvertChisqSumDf1ToPvalImhof(settings.getn("gene_imhof_epsilon_factor"), settings.getn("gene_imhof_pval_cutoff"));
      converter->set_backstop(1e3, (int) settings.getn("gene_imhof_perm_max"), max_comp, settings.geti("gene_imhof_adap_thresh"));
      imhof = converter;      
    } else if (subtype == "brown") imhof = new ConvertChisqSumDf1ToPvalBrown();
    else if (subtype == "empirical") imhof = new ConvertChisqSumDf1ToPvalEmpirical(1e3, (int) settings.getn("gene_imhof_perm_max"), max_comp, settings.geti("gene_imhof_adap_thresh"));    
    else throw EngineException("initializing Imhof model", "Invalid subtype specified");   
  }
  virtual ~SNPwiseMeanImhof() {delete imhof;}

  void clear() {SNPwise::clear();}
  ModelInfoBlock model_info(GeneStats* gs) {return make_info(gs, ModelState::SNPwiseModelMeanImhof);}
};

class SNPwiseTop : public SNPwise {
  ConvertPvalToLogPval to_logp;

  int top_count;
  float top_fraction;
  int xtx_max;

  int top_used;  
  bool use_all;
  bool log_stat;
  
  double run_core();
  double stat_from_pval(double* pval);
  double uni_pval(double curr_stat) {return log_stat ? exp(-0.5*curr_stat) : curr_stat;}  

  virtual void init_perm(Buffer<float>& data_block, int index=0) {UNUSED(data_block);
    if (log_stat) permutation->set_obs(stat(index));
    else permutation->set_obs(stat(index), &to_snp_stat, -1);
  }

  virtual StatConverter* get_converter() {return 0;}
public:
  SNPwiseTop(Settings& settings) : SNPwise(settings, to_logp) {
    top_count = settings.geti("gene_snpwise_top_count", 1);
    top_fraction = settings.getn("gene_snpwise_top_fraction", 0);
    xtx_max = settings.geti("gene_snpwise_top_xtx_max", 1000);
  }
  virtual ~SNPwiseTop(){}

  void clear() {SNPwise::clear();}
  ModelInfoBlock model_info(GeneStats* gs) {return make_info(gs, ModelState::SNPwiseModelTop);}  
  bool perm_only() {return true;}
};

class SNPwiseTopOne : public SNPwise {
  ConvertPermutationToPval gpd;
  BoundedMinP bounded_gpd;
  ConvertPvalToLogPval to_logp;  
  
  Buffer<float> process_buffer;
  Buffer<float> cross_buffer;

  int min_perm;
  int max_perm;
  int max_dim;
  float prune;
  float tail_thresh;
  float tail_margin_scale;
  int tail_min_perm; 
  
  int block_size;
  int adap_count;

  bool fit_distribution;
  bool proc_perm;

  SimulationMVN* fit_perm;
  Permutation* external_perm;  
  Buffer<double> perm_buffer;

  int curr_stat;
  int curr_perm;
  int external_cutoff;

  double run_core();
  double stat_from_pval(double* pval);

  void compute_var();

  virtual void init_perm(Buffer<float>& data_block, int index=0);
  virtual void process_perm(double perm_stat);
  virtual void flush_perm();

  virtual RawCorrelationData* get_covar() {return new PowerRawCorrelationData(data, min_snps, max_step);}
  virtual void init_gene();
  
  virtual StatConverter* get_converter() {return &bounded_gpd;}
public:
  SNPwiseTopOne(Settings& settings) : SNPwise(settings, to_logp), max_dim(1000), block_size(10000), fit_perm(0) {
    min_perm = settings.geti("gene_gpd_min_sim"); 
    max_perm = settings.geti("gene_gpd_nsim"); 
    adap_count = settings.geti("gene_gpd_adap_count");
    prune = settings.getn("gene_gpd_prune");
    tail_thresh = settings.getn("gene_gpd_tail");
    tail_margin_scale = settings.getn("gene_gpd_tail_margin") + 1;
    tail_min_perm = min(int(settings.geti("gene_gpd_tail_min")), int(0.75*max_perm)); 

    gpd.set_param("prop", settings.getn("gene_gpd_tail_fit"));
    gpd.set_param("min_count", 1000);    
    gpd.set_param("min_perm", tail_min_perm);  
    bounded_gpd.set_internal(&gpd, true);
  }
  virtual ~SNPwiseTopOne() {delete fit_perm;}

  void clear() {SNPwise::clear(); process_buffer.reset(); cross_buffer.reset(); perm_buffer.reset(); delete fit_perm; fit_perm = 0;}
  ModelInfoBlock model_info(GeneStats* gs) {return make_info(gs, ModelState::SNPwiseModelTopOne);}    
};


#endif /** GENEMODELENGINE_SNPWISE_H */

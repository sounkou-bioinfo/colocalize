/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENE_H
#define GENE_H

#include "parse.h"
#include "output.h"
#include "buffer.h"
#include "utils.h"
#include "mathutils.h"
#include "statutils.h"


using namespace std;
class GeneStats;

class Gene {
  friend class LinkedGene;
protected:
  static Buffer<float> temp_buffer;     
  enum BaseStatType {bs_Missing, bs_Freq, bs_MAF, bs_MAC, bs_Mean, bs_Variance, bs_Rescale, bs_Weights, bs_Count};
  enum DropType {dt_Retained=0, dt_Missing, dt_Mono, dt_LowFreq, dt_HighFreq, dt_LowVar, dt_DiffMiss, dt_InfoMiss, dt_NonRare};
  int bs_Extra;
  
  const Settings& settings;

  Buffer<float> data;
  BitBuffer missing;

  Buffer<float> base_stats;
  Buffer<double> snp_pval;
  Buffer<int> snp_N; 
  Buffer<int> snp_type; 
  Buffer<long> snp_id;
  Buffer<float> snp_weight;  

  int samp_size;
  int no_snps;
  int no_snps_orig;
  int no_rare;
  float miss_code;

  ProcessType::Type state;
  bool normalize_geno; 

  BufferWindow<float> get_datablock(int& offset, int& total);
  void check_block(int& offset, int& total);
  BufferWindow<float> get_basestats(BaseStatType type, int offset, int total) {check_block(offset, total); return BufferWindow<float>(base_stats[type]+offset, total);}

  void snp_center();
  void snp_normalize();
  void clear_missing();
  void prep_product();                                                  

  void set_rescale();
  virtual Buffer<float>& get_rescale_block(int& offset, int& total, bool norm=false);
  virtual float get_rescale_pair(Buffer<float>& bs_I, BitBuffer& miss_I, int i_snp, Buffer<float>& bs_J, BitBuffer& miss_J, int j_snp, bool samp_scale);

  virtual float sum_rescale(int snp, bool norm=false);
  void do_rescale_pair(float& target, int i_snp, int j_snp, ProductType::Type mode, Gene* other=0);

  void rescale_product(Buffer<float>& target, int offset_i, int total_i, int offset_j, int total_j, ProductType::Type mode, Gene* other=0);
  void rescale_product_skip(Buffer<float>& target, int offset_i, int step_i, int total_i, int offset_j, int step_j, int total_j, ProductType::Type mode, Gene* other=0);
  
  virtual float product_var_core(float* var, int snp, ProcessType::Type mode);
  
  virtual void init() = 0;
public:
  Gene(const Settings& s, int N=0, float miss=-999) : bs_Extra(0), settings(s), samp_size(N), miss_code(miss) {normalize_geno = settings["geno_normalize"];}
  virtual ~Gene() {}

  void set_N(int N) {samp_size = N;}
  bool is_raw() {return state == ProcessType::Raw;}
  
  int nsnp(bool orig=false) {return orig ? no_snps_orig : no_snps;}
  int nrare() {return no_rare;}  
  int nmax() {return data.maxcol();}
  
  float get_variance (int snp) {return base_stats(snp, bs_Variance);}
  BufferWindow<float> get_variance_block (int offset=0, int total=0) {return get_basestats(bs_Variance, offset, total);}
  BufferWindow<float> get_mac (int offset=0, int total=0) {return get_basestats(bs_MAC, offset, total);} 
  BufferWindow<double> get_pval(int offset=0, int total=0) {check_block(offset, total); return BufferWindow<double>(snp_pval.begin()+offset, total);}
  BufferWindow<int> get_snpN (int offset=0, int total=0) {check_block(offset, total); return BufferWindow<int>(snp_N.begin()+offset, total);}
  BufferWindow<int> get_snptype (int offset=0, int total=0) {check_block(offset, total); return BufferWindow<int>(snp_type.begin()+offset, total);}  
  BufferWindow<float> get_snpweight (int offset=0, int total=0) {check_block(offset, total); return BufferWindow<float>(snp_weight.begin()+offset, total);}    

  float product_var(float* var, int snp, ProcessType::Type mode=ProcessType::Centered); /// SNP x var (single)
  void product_var(Buffer<float>& target, Buffer<float>& var, int offset=0, int total=0, bool transpose=false, ProcessType::Type mode=ProcessType::Centered); /// SNP x var (block)

  virtual void product(Buffer<float>& target, int offset=0, int total=0, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (internal)
  void product(Buffer<float>& target, int offset_i, int total_i, int offset_j, int total_j, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (internal, off-diagonal)
  void product(Buffer<float>& target, Gene* other, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (external)
  void product(Buffer<float>& target, Gene* other, int offset_i, int total_i, int offset_j, int total_j, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (external)  
 
  void product_skip(Buffer<float>& target, int offset_i, int step_i, int total_i, int offset_j, int step_j, int total_j, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (internal, off-diagonal)  
  void product_skip(Buffer<float>& target, Gene* other, int step_i, int step_j, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (external)
  void product_skip(Buffer<float>& target, Gene* other, int offset_i, int step_i, int total_i, int offset_j, int step_j, int total_j, ProductType::Type mode=ProductType::Centered); /// SNP x SNP (external)

  virtual Buffer<float>& data_raw() {return data;}
  virtual void data_processed(Buffer<float>& target, int offset=0, int total=0, bool norm=false);
  virtual void data_product(Buffer<float>& target, Buffer<float>& var, int offset=0, int total=0);
  virtual void data_projection(Buffer<float>& target, Buffer<float>& rotate, int offset=0, int total=0);

  template<typename T>
  T* convert() {
    if (typeid(*this) != typeid(T)) {_LOG.error("processing gene") << "gene objects are incompatible" << endl; die();} 
    return static_cast<T*>(this);
  }  
};

class LoadedGene : public Gene {
protected:
  int gene_id;
  string gene_name;

  float max_miss;
  float min_maf;
  float max_maf;
  int min_mac;  
  int max_mac;
  double max_diff_corr;
    
  bool has_rare;
  bool rare_only;
  bool rare_normalize;
  float rare_maf_thresh;
  int rare_mac_thresh;  
  int rare_multi_max;
    
  bool full_qc;
  bool dummy_pheno;
  bool is_genotype;
  
  bool multi_partition;
  bool snp_report;
  static SharedPtr<DelimitedOutput> report_ptr; 

  
  vector<int> dropped_snps;
  vector<int> rare_snps; 
  vector<int> snp_status;
  
  virtual void prep_gene();
  virtual void init();

  void process_data(float* pheno);  
  void filter_snps();

  void process_rare() {process_rare(rare_snps, !snp_type.empty(), false); rare_snps.clear();} 
  void process_rare(vector<int>& index, bool by_type, bool no_max);
  void process_rare_core(vector<int>& index);  
  void process_rare_single(int index);

  int count_unique(vector<int>& index, bool flip=false);
public:
  LoadedGene(const Settings& s, int N, float miss) : Gene(s, N, miss), snp_report(false) {init();}
  virtual ~LoadedGene() {}
  
  void set_id(int id, string& name) {gene_id = id; gene_name = name;}
  Buffer<float>& set_buffer(int nsnp); 
  double* get_buffer_pval() {snp_pval.set_size(no_snps,1); return snp_pval.begin();}
  int* get_buffer_N() {snp_N.set_size(no_snps,1); return snp_N.begin();}  
  int* get_buffer_type() {snp_type.set_size(no_snps,1,true); return snp_type.begin();}    
  long* get_buffer_id(bool force) {
    if (force || multi_partition || snp_report) {snp_id.set_size(no_snps,1,true); return snp_id.begin();}
    else return 0;
  }      
  float* get_buffer_weight() {snp_weight.set_size(no_snps,1,true); return snp_weight.begin();}    

  virtual void process(float* pheno) = 0; 
};

class PlainBinaryGene : public LoadedGene {
public:
  PlainBinaryGene(const Settings& s, int N, float miss) : LoadedGene(s, N, miss) {}
  virtual ~PlainBinaryGene() {}

  virtual void process(float* pheno);

};

class ProcessedGene : public LoadedGene {
public:
  ProcessedGene(const Settings& s, int N, float miss) : LoadedGene(s, N, miss) {}
  virtual ~ProcessedGene() {}

  virtual void process(float* pheno);
};

class LinkedGene : public Gene {
  int outlier_thresh;

  Gene* parent;
  Buffer<float> parent_base_stats;
  ProcessType::Type parent_state;
  int parent_samp_size;
  double samp_size_rescale;
  
  int subset_step;
  int subset_offset;
  bool allow_prod;

  virtual Buffer<float>& get_rescale_block(int& offset, int& total, bool norm=false);
  virtual float get_rescale_pair(Buffer<float>& bs_I, BitBuffer& miss_I, int i_snp, Buffer<float>& bs_J, BitBuffer& miss_J, int j_snp, bool samp_scale);  
  virtual float sum_rescale(int snp, bool norm=false);  

  virtual float product_var_core(float* var, int snp, ProcessType::Type mode);

  virtual void init();
  
  void block_prod() {_LOG.error("processing gene") << "cannot compute product from doubly-linked LinkedGene object" << endl; die();} 
  void block_data() {_LOG.error("processing gene") << "data cannot be accessed through LinkedGene object" << endl; die();} 
public:
  LinkedGene(const Settings& s, int step=2, int offset=0) : Gene(s), parent(0), subset_step(step), subset_offset(offset) {
    outlier_thresh = settings.geti("gene_linked_outlier_thresh", 1);
  }
  virtual ~LinkedGene() {}

  void set_subset(int step, int offset) {subset_step = step; subset_offset = offset;}
  void set_parent(LoadedGene* par) {
    samp_size_rescale = 1; 
    allow_prod = true; parent = par; init();
  }
  void set_parent(LinkedGene* par) {
    samp_size_rescale = par->samp_size_rescale;
    allow_prod = false; parent = par; init();
  }
  
  using Gene::product;
  virtual void product(Buffer<float>& target, int offset=0, int total=0, ProductType::Type mode=ProductType::Centered);
  
  virtual Buffer<float>& data_raw() {block_data(); return data;}
  virtual void data_processed(Buffer<float>& target, int offset=0, int total=0, bool norm=false) {block_data(); UNUSED(target); UNUSED(offset); UNUSED(total); UNUSED(norm);}
  virtual void data_product(Buffer<float>& target, Buffer<float>& var, int offset=0, int total=0) {block_data(); UNUSED(target); UNUSED(var); UNUSED(offset); UNUSED(total);}
  virtual void data_projection(Buffer<float>& target, Buffer<float>& rotate, int offset=0, int total=0) {block_data(); UNUSED(target); UNUSED(rotate); UNUSED(offset); UNUSED(total);}

};


#endif /** GENE_H */
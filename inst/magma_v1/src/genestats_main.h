/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENESTATS_MAIN_H
#define GENESTATS_MAIN_H

#include <deque>

#include "genestats_base.h"
#include "parse.h"
#include "baseinput.h"
#include "utils.h"
#include "geneutils.h"

using namespace std;

class DefaultStats : public GeneStats {
protected:
  class PrinComp;
  
  bool sex_covar;
  bool chrX_only_sex_covar;
  bool sex_covar_all;
  bool sex_covar_chrX;
  bool use_residuals;
  bool use_interactor;
  bool pval_based;
  bool read_N;

  int snp_covar_index;  
  int sex_covar_index;

  Buffer<float> DUMMY_BUFFER;

  Buffer<float> phenotype;
  Buffer<float> phenotype_multi;
  Buffer<short> gender;  
  Buffer<float> covar;
  Buffer<float> interactor;

  Buffer<float> xty; /// #SNP x 1
  Buffer<float> xty_multi; /// #pheno x SNP
  Buffer<float> xtz; /// #COV x #SNP
  Buffer<float> xtx;
  Buffer<float> xtx_buffer;  

  bool has_xtx;
  bool has_xty;
  bool has_xtz;
  bool buffer_xtx;

  void compute_xty();
  void compute_xtz();
  virtual float compute_xtv(float* var, int snp, ProcessType::Type mode);
  virtual void compute_xtv(Buffer<float>& buffer, Buffer<float>& var, bool transpose, ProcessType::Type mode);
  
  void compute_xtx();
  virtual void compute_xtx(Buffer<float>& target, int offset, int total);

  virtual void compute_xtx_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode);
  virtual void compute_xtx_block_skip(Buffer<float>& buffer, int block_i, int step_i, int block_j, int step_j, ProductType::Type mode);
  virtual void compute_wtw_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode);

  void compute_wty(Buffer<float>& rotate);    
  void compute_wtz(Buffer<float>& rotate);    
  void compute_pcs(Buffer<float>& rotate);


  PrinComp* princomp;
  vector<PrinComp*> pc_data;
  int pc_id;

  void set_pc_buff(int size) {set_buffer(princomp->data, samp_size, size);}

  virtual void rewind();

  void process_covar();
  virtual void process_snps();
  void load_snp_info();  

public:
  DefaultStats(Settings& s, BaseInput& bi) : GeneStats(s, bi), sex_covar(false) {
    use_interactor = baseinput.interact;
    use_residuals = settings["use_residuals"];
    pval_based = settings["has_pval"];
    read_N = pval_based && settings["pval_col_N"];
    
    sex_covar_all = settings["covar_use_sex"];
    sex_covar_chrX = settings["chrX_sex_covar"]; 
  }
  virtual ~DefaultStats() {clear();}  
  virtual void clear();
  
  virtual void init();
  void set_chromosome(int chr);  

  virtual BufferDataBlock* eject_pcs(int pcid=-1); 
  
  void set_pcs(int index=0, bool init=false);
  virtual void partition(vector<int>& blocks);
  virtual void set_block(int index=-1);  
  virtual void store_stats(bool set) {buffer_xtx = set;}

  PermutationGenerator* perm_generator();
  
  Buffer<float> ytz;
  float ytz_prod;
 
  BufferWindow<float> get_xty();
  BufferWindow<float> get_xty_multi();  
  BufferWindow<float> get_xtz();
  BufferWindow<float> get_wty(int pcid=-1);
  BufferWindow<float> get_wty_multi(int pcid=-1);  
  BufferWindow<float> get_wty(Buffer<float>& rotate);
  BufferWindow<float> get_wty_multi(Buffer<float>& rotate);  
  BufferWindow<float> get_wtz(int pcid=-1);  
  BufferWindow<float> get_wtz(Buffer<float>& rotate);
  BufferWindow<float> get_xtx(bool invalidate=false);
  BufferWindow<float> get_capped_xtx(int max_size, bool invalidate=false);  
  
  
   
  BufferWindow<float> get_pcs(int pcid=-1); 
  BufferWindow<float> get_pcs(Buffer<float>& rotate); 
  void prep_pcs(Buffer<float>& rotate); 

  float load_block_cor(Buffer<float>& buffer, int block_i, int block_j, bool raw) {return load_block_xtx(buffer, block_i, block_j, raw, true, ProductType::Correlation);}
  float load_block_xtx(Buffer<float>& buffer, int block_i, int block_j, bool raw, ProductType::Type mode) {return load_block_xtx(buffer, block_i, block_j, raw, true, mode);}
  float load_block_xtx(Buffer<float>& buffer, int block_i, int block_j, bool raw, bool skip, ProductType::Type mode);

  void load_product(Buffer<float>& buffer, Buffer<float>& var, ProductType::Type mode=ProductType::Centered, bool col_snp=true); /// IF col_snp=true: #VAR x #SNP
  void load_pairwise_product(Buffer<float>& buffer, Buffer<float>& var, ProductType::Type mode=ProductType::Centered); 

   
  Buffer<float>& get_pheno() {return phenotype;}
  Buffer<float>& get_pheno_multi() {return phenotype_multi;}
  Buffer<float>& get_covar() {return covar;}
  Buffer<float>& get_interactor() {return interactor;}  
  
  void load_pheno(Buffer<float>& buffer, bool blur=false);
  void load_fitted(Buffer<float>& buffer, bool blur=false);  
  void load_data(Buffer<float>& buffer, bool norm=false);  
  void load_interaction(Buffer<float>& buffer, ProcessType::Type process=ProcessType::Raw);    

  BufferWindow<double> get_pval();
  void load_pval_corr(Buffer<float>& buffer, ProductType::Type mode=ProductType::Correlation, bool block=true);  


  int npcs(bool total=false) {return (total || !has_partition) ? princomp->npcs : princomp->index(curr_block,1);}      
  Buffer<int>& get_partition(bool raw, int index=0) {
    if (raw) return block_index;
    else {
      set_pcs(index);
      return princomp->index;   
    }
  }
    
  bool has_pval() {return pval_based;}
  bool has_interactor() {return use_interactor;}
  
  void write_vars(string prefix);  
  void write_gene(string prefix);

protected:
 struct PrinComp {                      
    Buffer<float> data;
    Buffer<float> wty; /// #PC x 1
    Buffer<float> wty_multi; /// #pheno x PC
    Buffer<float> wtz; /// #COV x #PC
    
    Buffer<int> index;
    bool partitioned;
   
    int npcs;   
    Buffer<bool> computed; ///By row, 0: PCs, 1: WtY, 2: WtZ
  
    PrinComp() {rewind();}
    void clear() {data.reset(); wty.reset(); wtz.reset(); index.reset(); computed.reset(); wty_multi.reset();}
    void rewind() {computed.set_size(1,3,true); index.set_size(0,2); npcs = 0; partitioned = false;}
   
    void set_partition(int nblocks) {computed.set_size(nblocks,3,true); partitioned = true;}
   
    BufferWindow<float> pc_block(Buffer<float>& source, int block, bool is_vector=false) {return block >= 0 ? source.window(index(block,0), index(block,1), is_vector) : source.window();}
    BufferWindow<float> get_pcs(int block=-1) {return pc_block(data, block);}
    BufferWindow<float> get_wty(int block=-1) {return pc_block(wty, block, true);}
    BufferWindow<float> get_wty_multi(int block=-1) {return pc_block(wty_multi, block);}    
    BufferWindow<float> get_wtz(int block=-1) {return pc_block(wtz, block);}    
   
    Buffer<float>& get_data() { 
      if (npcs <= 0) throw GeneException("processing gene data", "principal components for gene are not available");
      data.shrink(npcs); 
      return data;
    }
  };
}; 

class DefaultPermGen : public PermutationGenerator {
  Buffer<float> phenotype;
  Buffer<float> covariates;
  Buffer<float> perm_buffer;  
  
  int no_covar;

  void generate_plain(float* target);
  void generate_covar(float* target);

public:
  DefaultPermGen(DefaultStats* data, int N) : PermutationGenerator(N) {
    no_covar = data->ncovar();
    if (no_covar > 0) {
      data->load_fitted(phenotype, true);
      covariates.assign(data->get_covar());
    } else data->load_pheno(phenotype, true);
  }
  virtual ~DefaultPermGen(){}
  
  void run(BufferMode mode, int nperm);
  
  bool has_covar() {return no_covar > 0;}
};


#endif /** GENESTATS_MAIN_H */
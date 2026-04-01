/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENESTATS_SUBSET_H
#define GENESTATS_SUBSET_H

#include <deque>

#include "genestats_main.h"

using namespace std;
                                
class SubsetStats : public DefaultStats {
  class LinkedSubset;
  class Subset {
    friend class SubsetStats::LinkedSubset;
  protected:
    Settings& settings;
    LinkedGene* data; Gene* parent; bool loaded;
    void load_data();
  public:
    int step; int offset; int size;  
    int min_size; int max_step;    
    
    Subset(Settings& settings) : settings(settings), data(0), loaded(false), min_size(0), max_step(1) {}
    virtual ~Subset() {clear();}
    
    void clear() {delete data; data = 0;}
    virtual bool is_empty() {return data == 0;}
    
    void set(int samp_size, int step_init=0);
    virtual void set_gene(LinkedGene* d) {clear(); data = d;}    
    virtual void set_data(Gene* snp_data);
    
    virtual Gene* get_data();
    virtual Gene* eject_data();    
  };
  
  class LinkedSubset : public Subset {
    enum Mode {Indep, Cloned, Linked} mode;
    Subset* proxy; Subset* link; 
  public:
    LinkedSubset(Settings& settings) : Subset(settings), mode(Indep), proxy(0) {}  
    virtual ~LinkedSubset() {if (mode == Linked) delete proxy;}  

    virtual bool is_empty() {return (mode == Indep || !proxy) ? data == 0 : proxy->is_empty();}

    using Subset::set;
    void set(Subset& source, int samp_size);    
    virtual void set_gene(LinkedGene* d); 
    virtual void set_data(Gene* snp_data);    
    
    virtual Gene* get_data();
    virtual Gene* eject_data();    
  };
  
  struct SubsetPC {
    int processed;
    Buffer<float>* data;
   
    SubsetPC() : processed(0), data(0) {}
  };

  Subset sub_xtx;
  LinkedSubset sub_corrs;  
  
  bool auto_subset;
  int auto_size;
  
  vector<SubsetPC> pc_sub;
  Buffer<float> wtw;
  
  virtual LoadedGene* set_gene(int size);

  using DefaultStats::compute_xtv;  
  virtual void compute_xtv(Buffer<float>& buffer, Buffer<float>& var, bool transpose, ProcessType::Type mode);
    
  using DefaultStats::compute_xtx;
  virtual void compute_xtx(Buffer<float>& target, int offset, int total);
  virtual void compute_xtx_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode);
  virtual void compute_xtx_block_skip(Buffer<float>& buffer, int block_i, int step_i, int block_j, int step_j, ProductType::Type mode);
  virtual void compute_wtw(Buffer<float>& buffer, int block);
  virtual void compute_wtw_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode);

  void subset_pcs();

  virtual void clear();
  virtual void rewind() {DefaultStats::rewind(); for (int i = 0; i < pc_sub.size(); i++) pc_sub[i].processed = 0;}

  bool check_size(int size);
  virtual void process_snps();
public:
  SubsetStats(Settings& s, BaseInput& bi) : DefaultStats(s, bi), sub_xtx(s), sub_corrs(s), pc_sub(0) {
    int min_size = settings.geti("gene_model_bigdata_minsize", 1000);
    int max_step = max((int) settings.geti("gene_model_bigdata_maxstep", 100), 2); 
        
    sub_xtx.min_size = settings.geti("gene_model_bigdata_minsize_xtx", min_size);       
    sub_xtx.max_step = settings.geti("gene_model_bigdata_maxstep_xtx", max_step);        
    sub_corrs.min_size = min(min_size, sub_xtx.min_size);
    sub_corrs.max_step = max(max_step, sub_xtx.max_step);
    
    auto_subset = settings.gets("gene_model_bigdata") == "auto";
    auto_size = settings.geti("gene_model_bigdata_minsize_auto");       
  }
  virtual ~SubsetStats() {clear();}
  virtual void init();
  
  BufferWindow<float> get_wtw();
  
  virtual GeneDataBlock* eject_gene();
  virtual BufferDataBlock* eject_pcs(int pcid=-1);   
}; 


#endif /** GENESTATS_SUBSET_H */
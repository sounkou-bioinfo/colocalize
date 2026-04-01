/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENESTATS_BASE_H
#define GENESTATS_BASE_H

#include <deque>

#include "parse.h"
#include "baseinput.h"
#include "utils.h"
#include "geneutils.h"

using namespace std;
class PermutationGenerator;


class GeneStats {
public: 
  struct DataBlock {
    GeneLocation location;    

    DataBlock(GeneLocation loc) : location(loc) {}      
    virtual ~DataBlock() {}

    virtual string type() = 0;  
    virtual int size() = 0;
    
    virtual Gene* get_gene() {return 0;}    
    virtual Buffer<float>* get_buffer() {return 0;}
    virtual Gene* eject_gene() {return 0;}
    virtual Buffer<float>* eject_buffer() {return 0;}
  };
  
  struct GeneDataBlock : public DataBlock {
    Gene* data;
    bool owns_data;
    
    GeneDataBlock(GeneLocation loc, Gene* gene, bool owner=true) : DataBlock(loc), data(gene), owns_data(owner) {}      
    virtual ~GeneDataBlock() {if (owns_data) delete data;}

    virtual string type() {return "gene";}
    virtual int size() {return data != 0 ? data->nsnp() : 0;}
    virtual Gene* get_gene() {return data;}    
    virtual Gene* eject_gene() {
      Gene* out = data; data = 0;
      return owns_data ? out : 0;
    }
  };

  struct BufferDataBlock : public DataBlock {
    Buffer<float>* data;

    BufferDataBlock(GeneLocation loc, Buffer<float>* buff) : DataBlock(loc), data(buff) {} 
    BufferDataBlock(GeneLocation loc, Buffer<float>& buff) : DataBlock(loc) {
      data = new Buffer<float>;
      data->swap(buff);
    }      
    virtual ~BufferDataBlock() {delete data;}

    virtual string type() {return "buffer";}
    virtual int size() {return data != 0 ? data->ncol() : 0;}    
    virtual Buffer<float>* get_buffer() {return data;}
    virtual Buffer<float>* eject_buffer() {
      Buffer<float>* out = data; data = 0;
      return out;
    }
  };

  class GSID {
    void* signature;
    int instance;
    
  public:
    GSID(GeneStats* gs=0, int id=-1) : signature(gs), instance(id) {}
    bool operator==(const GSID& other) {return (signature == other.signature) && (instance == other.instance);}
    bool operator!=(const GSID& other) {return !(*this == other);}
  };

  enum Notice {ModelChange, DroppedCovar, DroppedCovarChrX, NoSubset, AutoSubset};
  bool has_message(Notice type) {return messages.find(type) != messages.end();}
  string get_message(Notice type) {return has_message(type) ? messages[type] : "";}    

protected:
  static int dump_count;
  static int instance_count;  
  map<Notice,string> messages;  

  int gs_id;  
  Settings& settings;
  BaseInput& baseinput;

  Buffer<float> dat_buffer;
  deque<Gene*> gene_storage;
  deque<Buffer<float>*> buffer_storage;  

  LoadedGene* snp_data;
  short* filter;

  int gene_id;
  int chromosome;
  
  int data_size;
  int samp_size;
  int orig_size;
  int no_pheno;
  int no_covar;
  int no_var;
  int tot_snps;
  int tot_rare;  
  float gene_mac;

  bool has_partition;
  Buffer<int> block_index;
  int block_offset;
  int block_size;
  int curr_block;
  int block_min_snps;
  int block_max_step;
  
  Gene* ejected_gene;
  BaseInput::SubsetMode selection;
  bool initialized;

  virtual void rewind() {block_index.set_size(0,2); has_partition = false; set_block(); ejected_gene = 0;}
  
  void insert_genes(DataBlock* block);
  void insert_buffers(DataBlock* block);
  void insert_buffer(Buffer<float>* curr_buff);

  Gene* fetch_gene(int size);
  Buffer<float>* fetch_buffer(int nrow, int ncol, bool zero=false);
  virtual LoadedGene* set_gene(int size);
  void set_buffer(Buffer<float>& target, int nrow, int ncol, bool zero=false);
  virtual LoadedGene* make_gene() {return new ProcessedGene(settings, samp_size, baseinput.misscode());}

  virtual void set_filter() {filter = baseinput.get_filter(selection);}

  virtual void process_snps() = 0;
  virtual void load_snp_info();

  int block_stepsize(int no_snps);

public:
  GeneStats(Settings& s, BaseInput& bi) : gs_id(-1), settings(s), baseinput(bi), snp_data(0), gene_id(-1), orig_size(-1), no_pheno(1), no_covar(0), ejected_gene(0), selection(BaseInput::sm_AllKnown), initialized(false) {
    data_size = baseinput.indiv_total;
       
    block_min_snps = settings.geti("block_corr_min_snps");
    block_max_step = settings.geti("block_corr_max_step");  
  }
  virtual ~GeneStats() {clear();}
  virtual void clear();

  virtual void init() = 0;
  virtual void set_chromosome(int chr);
  virtual void load_snps(int geneid);

  void set_message(Notice type, string msg="") {messages[type] = msg;}

  void insert_data(DataBlock* block);  
  virtual GeneDataBlock* eject_gene();

  virtual void partition(vector<int>& blocks);
  virtual int get_block() {return (curr_block > 0 && curr_block < block_index.nrow()) ? curr_block : 0;}
  virtual void set_block(int index=-1);
  virtual void store_stats(bool set) {UNUSED(set);}

  virtual PermutationGenerator* perm_generator() = 0;

  BufferWindow<int> get_snptype() {return snp_data->get_snptype(block_offset, block_size);}
  BufferWindow<float> get_var();
  virtual int get_sampsize(bool data_size=true) {return (!data_size && orig_size >= 0) ? orig_size : samp_size;}
  float get_mac() {return gene_mac;}
  int npheno() {return no_pheno;}
  int ncovar() {return no_covar;}
  int nvar(bool total=false) {return total ? no_var : block_size;}
  int nsnps() {return tot_snps;}
  int nrare() {return tot_rare;}
  int nblock() {return max(int(block_index.nrow()), 1);}
  
  virtual bool has_pval() {return false;}
  virtual bool has_interactor() {return false;}

  virtual bool equals(GSID& other) {return other == get_id();}
  GSID get_id() {return GSID(this, gs_id);}
  int get_geneid() {return gene_id;}  
  bool set_current(int& id) {
    bool changed = id != gene_id; id = gene_id; 
    return changed;
  }
  int set_current(int& id, int& block, int& block_tot) {int changed = 0;
    if (id != gene_id) changed = 2;
    else if ( (block != curr_block) || (block_tot != block_index.nrow()) ) changed = 1;
    id = gene_id; block = get_block(); block_tot = nblock();
    return changed;
  }

  template<typename T>
  T* convert() {
    T* out = dynamic_cast<T*>(this);
    if (out == 0) {_LOG.error("processing gene data") << "gene data objects are incompatible" << endl; die();}
    return out;
  }

  virtual void write_vars(string prefix) = 0;  
  virtual void write_gene(string prefix) = 0;
};

class PermutationGenerator {
public:
  enum BufferMode {Core, Transient};
protected:
  enum BufferType {Main, Covar};

  Buffer<float>* buffers[2][2];
  int samp_size;  
  
  bool prep_buffer(BufferMode mode, int nperm);
  
public:
  PermutationGenerator(int N) : samp_size(N) {
    buffers[0][0] = buffers[0][1] = buffers[1][0] = buffers[1][1] = 0;
  }
  virtual ~PermutationGenerator() {}
  
  virtual void run(BufferMode mode, int nperm) = 0;
  void set_buffer(Buffer<float>& main, Buffer<float>& covar, BufferMode mode) {buffers[mode][Main] = &main; buffers[mode][Covar] = &covar;}
  
  virtual bool has_covar() {return false;}
};

#endif /** GENESTATS_BASE_H */

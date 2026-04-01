/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SETSTATS_UTILS_H
#define SETSTATS_UTILS_H

#include <vector>

#include "data_genedata.h"

class GeneSetData;
template<typename T> class GeneCovarData;

namespace SetStatsUtils {

enum VariableType {vt_GeneSet, vt_GeneCovar}; 
enum InteractionType {it_InteractionSS, it_InteractionSC}; 
enum VariableClass {vc_AnyExternal, vc_GeneSet, vc_GeneCovar, vc_GeneCovarExternal, vc_Internal, vc_Interaction, vc_InteractionSS, vc_InteractionSC};
enum CollinearType {ct_None, ct_Variables, ct_Covariates, ct_All};

class InternalVariable {
public:
  enum Type {iv_GeneSize, iv_GeneDensity, iv_SampleSize, iv_InverseMac, iv_COUNT, iv_UNKNOWN}; ///aligns with name vector in names()  

  static vector<string>& names(bool printable=false);
  static string name(int ivar, bool printable=false) {return (ivar >= 0 && ivar < iv_COUNT) ? names(printable)[ivar] : "UNKNOWN";}  
  static Type type(const string& name, bool strip);
  static string format(const string& name);
};

class GeneMatrix {
protected:
  int dimension; 
  int rank;
  bool inverted;  
  
public:
  GeneMatrix() : inverted(false) {}
  virtual ~GeneMatrix() {}
  virtual void store_matrix(Buffer<double>& target) = 0;
  virtual int invert_matrix(bool verbose) = 0; 

  virtual void product(double* target) = 0;
  virtual void product(float* var, double* target) = 0;  
  virtual void product(double* var, double* target) = 0;  
  virtual void product(set<int>* var, double* target) = 0;  

  virtual double square(float* var) = 0;  
  virtual double square(double* var) = 0;  
  virtual double square(set<int>* var) = 0;  
  
  bool is_inverted() {return inverted;}
  int size(bool total=true) {return total ? dimension : rank;}
};

class CompoundGeneMatrix : public GeneMatrix {
  struct Block {
    int from; int size; bool inverted; 
    Block(int from, int size) : from(from), size(size), inverted(false) {}  
    virtual ~Block() {}
    virtual int invert() = 0;      
    virtual void store(Buffer<double>& target) = 0;    
    
    virtual void product(double*& target) = 0;
    virtual void product(float*& src, double*& target) = 0;  
    virtual void product(double*& src, double*& target) = 0;      
    virtual void product(set<int>::iterator& src, set<int>::iterator end, double*& target) = 0;  
    
    virtual double square(float*& var) = 0;  
    virtual double square(double*& var) = 0;  
    virtual double square(set<int>::iterator& var, set<int>::iterator end) = 0;  
  };
  
  struct CorrBlock : public Block {
    Buffer<double> corrs; double correction; double prune; bool mirrored;
    CorrBlock(int from, int size, double correction, double prune) : Block(from, size), correction(correction), prune(prune), mirrored(false) {init();}

    void init();  
    void add_value(int row, int col, double value);
    void add_row(int row, BaseBuffer<double>& data);
    
    void mirror();
    int invert();
    void store(Buffer<double>& target) {target.matrix_sub(from, size) = corrs.matrix();}    

    void product(double*& target);
    void product(float*& src, double*& target);  
    void product(double*& src, double*& target);      
    void product(set<int>::iterator& src, set<int>::iterator end, double*& target);  

    double square(float*& var);  
    double square(double*& var);  
    double square(set<int>::iterator& var, set<int>::iterator end); 
  };
  
  struct IndepBlock : public Block {
    IndepBlock(int from, int size) : Block(from, size) {}
    int invert() {inverted = true; return size;}
    void store(Buffer<double>& target) {for (int i = from; i < from+size; i++) target(i,i) = 1;}

    void product(double*& target);
    void product(float*& src, double*& target);  
    void product(double*& src, double*& target);      
    void product(set<int>::iterator& src, set<int>::iterator end, double*& target);  
    
    double square(float*& var);  
    double square(double*& var);  
    double square(set<int>::iterator& var, set<int>::iterator end);  
  };

  double correction;
  double prune;
  vector<Block*> blocks;

  void init(GeneData& data);
  void init(GeneData& data, set<int>& filter);  
public:
  CompoundGeneMatrix(GeneData& data, double corr_factor, double prune) : correction(corr_factor), prune(prune) {init(data);}
  CompoundGeneMatrix(GeneData& data, set<int>& filter, double corr_factor, double prune) : correction(corr_factor), prune(prune) {init(data, filter);}
  ~CompoundGeneMatrix() {Utils::delete_vector(blocks);}

  int invert_matrix(bool verbose);
  void store_matrix(Buffer<double>& target);

  void product(double* target) {for (int i = 0; i < blocks.size(); i++) blocks[i]->product(target);}
  void product(float* var, double* target) {for (int i = 0; i < blocks.size(); i++) blocks[i]->product(var, target);}
  void product(double* var, double* target) {for (int i = 0; i < blocks.size(); i++) blocks[i]->product(var, target);}  
  void product(set<int>* var, double* target);  
  
  double square(float* var);  
  double square(double* var);  
  double square(set<int>* var);  
};

class SetInfo; class CovarInfo; class InternalCovarInfo; class InteractionPair; 

class VariableInfo {
protected:
  int var_id; int iid; 
public:
  VariableInfo(int var_id, int iid) : var_id(var_id), iid(iid) {}
  virtual ~VariableInfo() {}

  virtual bool exists() {return true;}
  virtual bool is_valid() {return get_status() == VariableIndex::Valid;}

  virtual bool is_set() {return false;}
  virtual bool is_covar() {return false;}
  virtual bool is_interaction() {return false;}
  virtual bool has_type(VariableClass type) = 0;

  int get_id() {return var_id;}  
  virtual char* get_name() = 0;
  virtual string get_label() {return "";}
  virtual string get_name_label() {return string(get_name()) + " " + get_label();}
  virtual short get_status() = 0;
  virtual int get_column() = 0;
  
  virtual SetInfo* as_set() {return 0;}
  virtual CovarInfo* as_covar() {return 0;}
  virtual InteractionPair* get_interaction() {return 0;}
};

class DummyInfo : public VariableInfo {
  static char dummy;
public:
  DummyInfo() : VariableInfo(-1,-1) {}

  bool exists() {return false;}
  bool has_type(VariableClass type) {UNUSED(type); return false;}

  char* get_name() {return &dummy;}
  string get_name_label() {return string(&dummy);}
  short get_status() {return VariableIndex::NotAvailable;}
  int get_column() {return -1;}   
};

struct InteractionPair {
  virtual InteractionType type() = 0;
  virtual VariableInfo* var(int i) = 0;
  set<int> ids() {set<int> out; out.insert(var(0)->get_id()); out.insert(var(1)->get_id()); return out;}
  string get_name() {return string("INTERACT::") + var(0)->get_name() + "::" + var(1)->get_name();}
};

struct SetSetPair : public InteractionPair {
  SetInfo *set1, *set2;
  SetSetPair(SetInfo* set1, SetInfo* set2) : set1(set1), set2(set2) {}
  InteractionType type() {return it_InteractionSS;}
  VariableInfo* var(int i); 
};

struct SetCovPair : public InteractionPair {
  SetInfo* set; CovarInfo* cov;
  SetCovPair(SetInfo* set, CovarInfo* cov) : set(set), cov(cov) {}
  InteractionType type() {return it_InteractionSC;}
  VariableInfo* var(int i); 
};


class SetInfo : public VariableInfo {
  GeneSetData* index;
public:
  SetInfo(GeneSetData* index, int var_id, int iid) : VariableInfo(var_id, iid), index(index) {}

  bool is_set() {return true;}
  virtual bool has_type(VariableClass type) {return type == vc_GeneSet || type == vc_AnyExternal;}
  
  char* get_name();
  virtual string get_label() {return "(set)";}  
  short get_status();
  int get_column();
  int get_ngenes();
  set<int>* get_data();
  double compute_sd(int tot_genes);
  GeneSetData* get_storage() {return index;}
  
  SetInfo* as_set() {return this;}
  virtual SetSetPair* get_interaction() {return 0;}
};

class CovarInfo : public VariableInfo {
protected:
  GeneCovarData<float>* index;
public:
  CovarInfo(GeneCovarData<float>* index, int var_id, int iid) : VariableInfo(var_id, iid), index(index) {}
  virtual ~CovarInfo() {}

  bool is_covar() {return true;}
  virtual bool has_type(VariableClass type) = 0; 
  
  char* get_name();
  virtual string get_label() {return "(covar)";}
  short get_status();
  int get_column();
  int get_obs_genes();
  double get_orig_mean();
  double get_orig_sd();
  float* get_data();
  GeneCovarData<float>* get_storage() {return index;}
  
  CovarInfo* as_covar() {return this;}
  virtual InternalCovarInfo* as_internal() {return 0;}
  SetCovPair* get_interaction() {return 0;}
};

class ExternalCovarInfo : public CovarInfo {
public:
  ExternalCovarInfo(GeneCovarData<float>* index, int var_id, int iid) : CovarInfo(index, var_id, iid) {}  

  bool has_type(VariableClass type) {return (type == vc_GeneCovar || type == vc_GeneCovarExternal || type == vc_AnyExternal);}
};

class InternalCovarInfo : public CovarInfo {
  InternalVariable::Type subtype;
public:
  InternalCovarInfo(GeneCovarData<float>* index, int var_id, int iid, short type) : CovarInfo(index, var_id, iid) {subtype = (type >= 0 && type < InternalVariable::iv_COUNT) ? (InternalVariable::Type) type : InternalVariable::iv_UNKNOWN;}  

  bool has_type(VariableClass type) {return (type == vc_GeneCovar || type == vc_Internal);}
  InternalVariable::Type get_subtype() {return subtype;}
  string format_name();
  
  InternalCovarInfo* as_internal() {return this;}
};


class SetSetInteractionInfo : public SetInfo {
  SetSetPair info;
public:
  SetSetInteractionInfo(GeneSetData* index, int var_id, int iid, SetInfo* main1, SetInfo* main2) : SetInfo(index, var_id, iid), info(main1, main2) {}  

  string get_label() {return "(set x set)";}  
  bool is_interaction() {return true;}
  bool has_type(VariableClass type) {return SetInfo::has_type(type) || (type == vc_Interaction) || (type == vc_InteractionSS);}
  SetSetPair* get_interaction() {return &info;}
};


class SetCovInteractionInfo : public ExternalCovarInfo {
  SetCovPair info;
public:
  SetCovInteractionInfo(GeneCovarData<float>* index, int var_id, int iid, SetInfo* main1, CovarInfo* main2) : ExternalCovarInfo(index, var_id, iid), info(main1, main2) {}  

  string get_label() {return "(set x cov)";}  
  bool is_interaction() {return true;}
  bool has_type(VariableClass type) {return ExternalCovarInfo::has_type(type) || (type == vc_Interaction) || (type == vc_InteractionSC);}
  SetCovPair* get_interaction() {return &info;}
};


template<typename T>
class EigenInverter {
  Buffer<T>& data;
  SelfAdjointEigenSolver<Matrix<T,Dynamic,Dynamic> > eigen;
  
  BaseBuffer<T> weights;
  bool consecutive;

  double cutoff;

  int size;
  int offset;
  int rank;

  void init(int use) {
    size = (use > 0) ? use : data.nrow(); rank = size; 
    if (data.nrow() != data.ncol()) throw MathException("matrix inversion", "trying to invert a non-square matrix");
    if (data.nrow() == 0) throw MathException("matrix inversion", "trying to invert empty matrix");
    if (data.nrow() < offset+size || data.ncol() < offset+size) throw MathException("matrix inversion", "input matrix for inversion is smaller than specified dimension");

    try { 
      eigen = SelfAdjointEigenSolver<Matrix<T,Dynamic,Dynamic> >(data.matrix_sub(offset, size));
    } catch (const exception& e) {check_mem_error(e, "inverting matrix"); throw MathException("matrix inversion", "an error occurred when trying to invert matrix", e.what());}      
  }
  
  bool has_weights() {
    if (weights.size() != size) throw MathException("matrix inversion", "incorrect number of weights");
    return true;
  }
  
  void compute_weights(bool root, bool invert) {
    const T* values = eigen.eigenvalues().data(); /// should be in ascending order
    weights.resize(size); int dropped = 0; consecutive = true;
    for (int i = 0; i < size; i++) {
      if (values[i] <= cutoff) {
        if (dropped < i) consecutive = false;        
        weights[i] = 0; dropped++;
      } else {
        if (invert) weights[i] = root ? 1/sqrt(values[i]) : 1/values[i]; 
        else weights[i] = root ? sqrt(values[i]) : values[i]; 
      }        
    }    
    rank = size - dropped;
  }

  void compute_root(Buffer<T>& target, bool clear_weights=true) { try { 
    if (has_weights()) {
      target.set_size(size, rank);
       
      if (!consecutive) {
        for (int r = 0, w = 0; r < size; r++) {
          if (weights[r] > 0) {
            Map<Matrix<T,Dynamic,1> >(target[w++], size, 1).noalias() = eigen.eigenvectors().col(r) * weights[r];
          }
        }
      } else {
        int offset = size - rank;
        Map<Matrix<T,Dynamic,Dynamic> > out_map(target[0], size, rank);              
        out_map.noalias() = eigen.eigenvectors().block(0, offset, size, rank) * Map<Matrix<T,Dynamic,1> >(weights.data() + offset, rank, 1).asDiagonal();      
      }
    } 
    
    if (clear_weights) {weights.clear(); consecutive = true;}
  } catch (const exception& e) {check_mem_error(e, "inverting matrix"); throw MathException("matrix inversion", "an error occurred when trying to invert matrix", e.what());} }     
  
  void compute_inverse(Buffer<T>& target, bool emplace=false, int offset=0, bool clear_weights=true) { try { 
    if (has_weights()) {
      if (emplace) {
        if (target.ncol() < size+offset || target.nrow() < size+offset) throw MathException("matrix inversion", "output buffer for inversion is too small");
      } else target.set_size(size, size); 
      Map<Matrix<T,Dynamic,Dynamic>, 0, OuterStride<> > out_map(target[offset]+offset, size, size, target.nrow());
      if (consecutive) {
        int offset = size - rank;
        out_map.noalias() = eigen.eigenvectors().block(0, offset, size, rank) * Map<Matrix<T,Dynamic,1> >(weights.data() + offset, rank, 1).asDiagonal() * eigen.eigenvectors().block(0, offset, size, rank).transpose();
      } else {
        out_map.noalias() = eigen.eigenvectors() * Map<Matrix<T,Dynamic,1> >(weights.data(), size, 1).asDiagonal() * eigen.eigenvectors().transpose();
      }  
    } 

    if (clear_weights) {weights.clear(); consecutive = true;}  
  } catch (const MathException& me) {throw me;}
    catch (const exception& e) {check_mem_error(e, "inverting matrix"); throw MathException("matrix inversion", "an error occurred when trying to invert matrix", e.what());} 
  }   
  
  double compute_logdeterminant() {
    double log_det = 0; const T* values = eigen.eigenvalues().data();  
    for (int i = 0; i < size; i++) log_det += (values[i] > 0) ? log(values[i]) : log(0);
    return log_det;    
  }  
  
public:
  EigenInverter(Buffer<T>& input, int use=0, int offset=0) : data(input), cutoff(0.1), offset(offset) {init(use);}
  
  double set_cutoff(double value, bool estimate=false) {
    if (estimate) {
      const T* eval = eigen.eigenvalues().data(); 
      if (estimate) {double sum = 0;
        for (int i = 0; i < size; i++) sum += eval[i];
        cutoff = value * sum / size;
      }
    } else cutoff = value;
    return cutoff;
  }

  void root(Buffer<T>& target) {compute_weights(true, false); compute_root(target);}
  void root_inverse(Buffer<T>& target) {compute_weights(true, true); compute_root(target);}  
  void inverse(Buffer<T>& target, bool emplace=false, int offset=0) {compute_weights(false, true); compute_inverse(target, emplace, offset);}

  double determinant(bool log_scale) {
    double log_det = compute_logdeterminant();
    return !log_scale ? exp(log_det) : log_det;
  }
  
  pair<double,double> determinant_ratio() {
    double log_det = compute_logdeterminant(), log_trace = 0;
    for (int i = offset, end = offset+size; i < end; i++) {
      T& curr = data(i,i); log_trace += (curr >= 0) ? log(curr) : log(0);
    }
    return pair<double,double>(log_det, log_ratio(log_det, log_trace));
  }

  int get_rank() {return rank;}
  
  static double log_ratio(double log_a, double log_b) {return (boost::math::isfinite)(log_a) && (boost::math::isfinite)(log_b) ? exp(log_a - log_b) : 0;}
  
  static int invert(Buffer<T>& input, double cutoff, bool estimate=false) {
    EigenInverter eigen(input); eigen.set_cutoff(cutoff, estimate);    
    eigen.inverse(input);
    return eigen.get_rank();
  }
};

/// For inversion of block-matrix of below structure, with A previously inverted
/// matrix X:   A    B
///             B^T  D
class BlockInverter {
  typedef Map<Matrix<double,Dynamic,Dynamic>, 0, OuterStride<> > SubMap;

	double ratio_thresh, eigen_prune;
	bool allow_collinear;

	Buffer<double> M;
	Buffer<double> W;
	Buffer<double> WM;

public:
	BlockInverter(double thresh, double prune, bool collinear=false) : ratio_thresh(thresh), eigen_prune(prune), allow_collinear(collinear) {}

  CollinearType invert(Buffer<double>* X, Buffer<double>* A_inv, Buffer<double>* target);
};

}; /// End of namespace


class SetVariables;

class GeneSetData : public VariableIndex {
public:
  PositiveNumber<int> mapped_genes; 
  NonNegativeNumber<int> used_genes;   
 
  SparseBooleanVariable data; 

  GeneSetData(const string& id="", const string& file="") : VariableIndex(id, file) {
    register_variable(mapped_genes, "MAPPED_GENES", true);  
    register_variable(used_genes, "USED_GENES", true); 
  }
  ~GeneSetData() {}
  
  void link_data(DataFrame& df) {df.link_variable(data, data_id);}
  void unlink_data(bool do_delete);
};

template<typename T>
class GeneCovarData : public CovariateIndex<NumericVariable<T> > {
public:
  static const short ProcessBlock;
  static const short StatsBlock;

  NumericVariable<double> orig_mean;
  NumericVariable<double> orig_sd;
  NonNegativeNumber<int> truncation_count;

  GeneCovarData(string id="", string file="") : CovariateIndex<NumericVariable<T> >(id, file) {
    this->register_variable(orig_mean, "ORIG_MEAN", true, StatsBlock);     
    this->register_variable(orig_sd, "ORIG_SD", true, StatsBlock);     
    this->register_variable(truncation_count, "TRUNC_COUNT", true, ProcessBlock);     
  }
  ~GeneCovarData() {}
  
  void block_init(short block) {if (block > StatsBlock) block_init(StatsBlock); CovariateIndex<NumericVariable<T> >::block_init(block);}
}; 

template<typename T> const short GeneCovarData<T>::StatsBlock = 10;
template<typename T> const short GeneCovarData<T>::ProcessBlock = 11;

class GeneVariableMap : public IndexedDataFrame {
  SetVariables& variable_data;
  SetStatsUtils::DummyInfo dummy;
  
public:
  PointerVariable<SetStatsUtils::VariableInfo> variable;
  
  GeneVariableMap(SetVariables& vars) : IndexedDataFrame(false), variable_data(vars), variable(true) {
    id_index->set_param("drop_duplicates", true);
    register_variable(variable, "VAR_INFO", true);              
  }

  int add_variable(const string& prefix, const string& name, int block_index, int offset, SetStatsUtils::VariableType type);
  int set_interaction(int iid, int main1_id, int main2_id, int block_index, int offset, SetStatsUtils::InteractionType type);
  
  set<int> find_var(const string& name, SetStatsUtils::VariableClass type=SetStatsUtils::vc_AnyExternal);
  bool var_exists(const string& name, SetStatsUtils::VariableClass type=SetStatsUtils::vc_AnyExternal) {return !find_var(name, type).empty();}

  SetStatsUtils::VariableInfo& var_info(int var_id) {return (var_id >= 0 && active_id(var_id) && variable.get(var_id)) ? *variable.get(var_id) : dummy;} 

  set<int> var_list(SetStatsUtils::VariableClass type) {set<int> dummy; return var_list(type, dummy);}
  set<int> var_list(SetStatsUtils::VariableClass type, set<int>& filter);  
  
  pair<int,int> no_sets();
  pair<int,int> no_covar(bool skip_internal);  
  pair<int,int> no_vars(bool skip_internal);
};


#endif /** SETSTATS_UTILS_H */


/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef UTILS_DATAFRAME_H
#define UTILS_DATAFRAME_H

#include <algorithm>
#include <map>

#include "buffer.h"
#include "mathutils.h"
#include "utils.h"
#include "utils_dataframe_ids.h"
#include "utils_dataframe_variables.h"


class DataFrame {
  static map<DataFrame*,void*> data_owner;
  static bool is_owner(DataFrame* df, void* obj) {return (data_owner.find(df) != data_owner.end() && data_owner[df] == obj);}
    
public:
  static void set_owner(DataFrame* df, void* obj) {data_owner[df] = obj;}
  static void unset_owner(DataFrame* df, void* obj) {if (is_owner(df, obj)) data_owner.erase(df);} 
  static void delete_owned(DataFrame* df, void* obj) {if (is_owner(df, obj)) {delete df; data_owner.erase(df);}} 

protected: 
  static const short GlobalBlock;
  
  void* owner;
  BaseBuffer<int> local_map;  /// internal ID to storage position
  vector<Variable*> variables;

  vector<bool> variable_init;
  vector<short> variable_block;
  map<string,int> variable_index;
  vector<Variable*> attached_variables;

  int index_max;   /// internal ID range (fixed)
  int index_used;  /// valid entries
  int storage_req;

  BaseBuffer<int> output_mask;
  
  virtual bool valid_id(int iid) {return (iid >= 0 && iid < index_max && local_map[iid]);}
  virtual int validate_id(int iid) {return (iid >= 0 && iid < index_max && local_map[iid]) ? iid : -1;}

  void clear();
  void init_storage(int size);
  void init_variables(int size, short block=GlobalBlock);  
  
  void count_used();
  int max_storage();

  void reorder(const vector<int>& order); ///yields consecutive and iid-aligned storage, invalidates old iids
  virtual void remap_external(vector<int>& remap) {UNUSED(remap);}

  bool make_consecutive(); ///yields consecutive and iid-aligned storage, invalidates old iids
  template<typename T> void data_sort(UniVariable<T>& var); ///yields consecutive and iid-aligned storage, invalidates old iids
  
  void register_variable(Variable& var, const string& name="", bool auto_init=false, short init_block=GlobalBlock);
  int registered(Variable& var);    

  void drop_id(int iid) {if (valid_id(iid)) {local_map[iid] = 0; index_used--;}}
  void filter_ids(vector<short>& filter, bool include);

  DataFrame() : owner(0), index_max(0), index_used(0), storage_req(0) {} 
public:
  virtual ~DataFrame() {for (int i = 0; i < attached_variables.size(); i++) delete attached_variables[i];}
 
  bool active_id(int iid) {return valid_id(iid);}

  Variable* get_variable(const string& name) {return (variable_index.find(name) != variable_index.end()) ? variables[variable_index[name]] : 0;}
  template<typename VAR> VAR* typed_variable(const string& name) {return dynamic_cast<VAR*>(get_variable(name));}
  bool using_variable(const string& name) {Variable* var = get_variable(name); return var ? var->initialised() : false;}

  virtual void block_init(short block) {init_variables(0, block);}

  long data_used(bool update) {if (update) count_used(); return index_used;} 
  long data_total() {return index_max;}  

  void shrink_storage(float threshold=1); ///shrinks storage size, does not invalidate iids
  bool is_consecutive(); ///checks if storage is consecutive and iids map directly onto storage
  
  void load_map(BaseBuffer<int>& target) {target.assign(local_map);} 
  template<typename T> void store_rows(NumericID<T>& target);

  ///only for writing output; mask will become invalid if underlying storage is changed
  ///input specified in storage ids
  int set_mask(set<int>& sids); 
  bool unset_mask(bool force=false);

  void push_message(Variable::MessageType type, double value=0);

  void align_variable(Variable& var, bool do_init=false); ///align storage with DataFrame
  void link_variable(Variable& var) {link_variable(var, "", true);} ///align() and register with DataFrame
  void link_variable(Variable& var, const string& name, bool do_init=true);
  void attach_variable(Variable* var, const string& name, bool do_init=true); ///link() and transfer ownership to DataFrame
  void unlink_variable(Variable& var);
};

class ExpandingDataFrame : public DataFrame {
protected:
  vector<int> id_index;
  
  bool index_built;

  void remap_external(vector<int>& remap) {for (int i = 0; i < index_max; i++) id_index[i] = remap[id_index[i]];}

public:
  NumericID<int> input_id;

  ExpandingDataFrame() : index_built(false) {
    register_variable(input_id, "INPUT_ID");
  } 
  virtual ~ExpandingDataFrame() {}

  bool init(int exp_amount);
  void clear_data() {clear();}
  void lock_index(bool store_eids=false);    

  virtual int add_row();
  virtual int add_block(int size);

  int operator[](int eid) {return index_built && eid < id_index.size() ? validate_id(id_index[eid]) : -1;}
};

class IndexedDataFrame : public DataFrame {
protected:
  IDIndex* id_index;  /// external ID to internal ID
  vector<char*> duplicates;
  int added; 
  bool name_var;
  
  bool index_built;
  bool index_available;
  
  bool ensure_index();
  void remap_external(vector<int>& remap) {id_index->update_order(remap);}

  virtual void add_iid(const string& eid) {if (!index_built) id_index->add_id(eid);}
  virtual int get_iid(const string& eid) {if (!index_built) {id_index->add_id(eid); return added++;} else return added;} 

  void remap_id(char* eid, int iid) {if (index_available) id_index->update_id(eid, iid);}

public:
  CStringVariable name;
  void store_idnames();

  IndexedDataFrame(bool keep_names=false) : added(0), name_var(keep_names), index_built(false), index_available(false) {
    id_index = new IDIndexCore();
    register_variable(name, "NAME");
  } 
  virtual ~IndexedDataFrame() {delete id_index;}

  virtual bool init(int exp_amount, int exp_length=0);

  virtual void clear_data(bool clear_ids=true);
  virtual void clear_index(bool keep_names=false);  
  virtual void build_index();
  virtual int extend_index(int amount);    

  void align_map(IndexedDataFrame& master, int* target, int size, bool keep_mid=false); ///yields consecutive and iid-aligned storage, invalidates old iids
  
  int id_size() {return id_index->id_size();} 
  vector<char*>& get_duplicates() {return duplicates;}

  int operator[](char* eid) {return index_available ? validate_id(id_index->get_index(eid)) : -1;}
  int operator[](const string& eid) {return index_available ? validate_id(id_index->get_index(const_cast<char*>(eid.c_str()))) : -1;}
};

class VariableIndex : public IndexedDataFrame {
public:
  enum VarStatus {Valid=0, Invalid, HighMissing, LowVariance, Monotonic, NotAvailable};

  string data_id;
  string data_file;

  NonNegativeNumber<short> status;
  NumericID<int> storage_column; ///column ID in data storage variable
  NumericID<int> variable_id; ///optional ID from external reference object 

  
  VariableIndex(string id="", string file="") : IndexedDataFrame(true), data_id(id), data_file(file) {
    id_index->set_param("drop_duplicates", true);
    register_variable(status, "STATUS", true);          
    register_variable(storage_column, "COL_ID", true);          
    register_variable(variable_id, "VAR_ID", false);          
  }
  virtual ~VariableIndex() {}
  using IndexedDataFrame::get_iid;
  
  virtual void link_data(DataFrame& df) = 0;
  virtual void unlink_data(bool do_delete) = 0;
  
  bool is_valid(int iid) {return status.get(iid) == Valid;}
  int count_valid();
};

template<typename VAR>
class CovariateIndex : public VariableIndex {
public:
  /// per-covariate 
  NumericID<int> file_column; ///column ID from input file
  BooleanVariable complete_column;
  NonNegativeNumber<int> values_observed;  

  /// per-gene
  BlockMultiVariable<VAR> data;  
  BooleanVariable in_source;
  NonNegativeNumber<int> miss_count;
 
  CovariateIndex(string id="", string file="") : VariableIndex(id, file) {
    register_variable(file_column, "FILE_COL", file != "");          
    register_variable(complete_column, "FILE_OBS", file != "");              
    register_variable(values_observed, "VALUE_OBS", true);
  }
  virtual ~CovariateIndex() {}
  
  void build_index();
  int extend_index(int amount);   
  void link_data(DataFrame& df) {link_data(df, false);}
  void link_data(DataFrame& df, bool init_zero);
  void unlink_data(bool do_delete);
};

#include "utils_dataframe.tpp"

#endif /** UTILS_DATAFRAME_H */


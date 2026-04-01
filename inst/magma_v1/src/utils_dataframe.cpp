/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "_global.h"
#include "utils_dataframe.h"
                                                                                
#define SHRINK_THRESHOLD 0.90

map<DataFrame*,void*> DataFrame::data_owner;
const short DataFrame::GlobalBlock = 1;

void DataFrame::push_message(Variable::MessageType type, double value) {
  if (type == Variable::mt_Filtered && value >= 0) count_used();
}

int DataFrame::set_mask(set<int>& sids) {
  if (local_map.is_empty()) {_LOG.error("processing DataFrame object") << "trying to mask unmapped DataFrame object" << endl; die();}
  if (!is_consecutive()) {_LOG.error("processing DataFrame object") << "trying to mask non-consecutive DataFrame object" << endl; die();}

  output_mask.resize(sids.size());
  if (!sids.empty()) {int curr = 0;
    for (set<int>::iterator it = sids.begin(); it != sids.end(); ++it) {
      if (*it >= index_used) {output_mask.resize(curr); break;}
      output_mask[curr++] = *it + 1;
    }
  }
  
  if (!output_mask.is_empty()) {
    for (int i = 0; i < variables.size(); i++) variables[i]->set_mask(&output_mask);  
  } else unset_mask(true);
  
  return output_mask.size();
}

bool DataFrame::unset_mask(bool force) {
  if (!output_mask.is_empty() || force) {
    for (int i = 0; i < variables.size(); i++) variables[i]->set_mask(0);
    output_mask.clear(); return true;
  } else return false;
}

void DataFrame::align_variable(Variable& var, bool do_init) {
  if (!local_map.is_empty()) var.set_map(storage_req, &local_map);  
  else var.set_size(storage_req);
  if (!output_mask.is_empty()) var.set_mask(&output_mask);
  if (do_init) var.init();    
}

void DataFrame::link_variable(Variable& var, const string& name, bool do_init) {
  register_variable(var, name);
  align_variable(var, do_init);
}

void DataFrame::attach_variable(Variable* var, const string& name, bool do_init) {
  link_variable(*var, name, do_init);
  attached_variables.push_back(var);
}

void DataFrame::unlink_variable(Variable& var) {
  int var_id = registered(var);
  if (var_id >= 0) {
    variables.erase(variables.begin() + var_id);      
    variable_init.erase(variable_init.begin() + var_id);      
    variable_block.erase(variable_block.begin() + var_id);      
    unlink_variable(var); 
  }
  var.set_mask(0);
}
 
void DataFrame::register_variable(Variable& var, const string& name, bool auto_init, short init_block) {
  var.set_owner(this); variables.push_back(&var); 
  variable_init.push_back(auto_init); variable_block.push_back(init_block);
  if (name != "") {
    if (variable_index.find(name) != variable_index.end()) {_LOG.error("processing DataFrame object") << "variable ID '" << name << "' is already in use" << endl; die();}
    variable_index[name] = variables.size()-1;
  }
}

int DataFrame::registered(Variable& var) {
  for (int i = 0; i < variables.size(); i++) {if (variables[i] == &var) return i;}
  return -1;
}

void DataFrame::clear() {
  unset_mask();
  index_max = index_used = storage_req = 0;
  for (int i = 0; i < variables.size(); i++) variables[i]->clear_data();
  local_map.clear();
}

void DataFrame::init_storage(int size) {
  if (!local_map.is_empty()) return;

  storage_req = size; 
  local_map.resize(index_max); 
  for (int i = 0; i < index_max; i++) local_map[i] = i+1;
  for (int i = 0; i < variables.size(); i++) variables[i]->set_map(storage_req, &local_map);
}

void DataFrame::init_variables(int size, short block) {
  for (int i = 0; i < variables.size(); i++) {
    if (size) variables[i]->set_size(size);
    if (variable_block[i] == block && variable_init[i]) variables[i]->init();
  }
}

///counts number of iids still in use (index_used) and the minimum storage size needed
void DataFrame::count_used() {
  int *map = local_map.data(), size = local_map.size(); index_used = 0; storage_req = 0;
  for (int i = 0; i < size; i++) {
    if (map[i] > 0) {
      index_used++;
      if (map[i] > storage_req) storage_req = map[i];
    }
  }
}

int DataFrame::max_storage() {int max = 0;
  for (int i = 0; i < variables.size(); i++) {
    int curr = (variables[i]->initialised()) ? variables[i]->data_size() : 0;
    if (curr > max) max = curr;
  }
  return max;
}

///rearranges the data according to the order input vector, defining an ordering based on stored values
///that is, vector order is of data stored in variables, not iid (order[new pos.] = curr. pos.; skips if curr. pos. entry is no longer in use)
///results in consecutive storage in the given order, and consecutive and ordered local_map (up to index_used, 0s thereafter)
///therefore, iid and storage id are identical, and storage buffers can be accessed directly with iid (note: old iids are invalidated)
void DataFrame::reorder(const vector<int>& order) {
  if (index_max <= 0) return;
  if (local_map.is_empty()) {_LOG.error("processing DataFrame object") << "unable to sort index" << endl; die();}

  vector<int> inverse_map(index_max, -1); ///map from storage position to local_map entry; not in use = -1
  for (int i = 0; i < index_max; i++) {if (local_map[i]) inverse_map[local_map[i]-1] = i;}
  vector<int> filtered_order; ///order vector with unused entries skipped
  for (int i = 0; i < order.size(); i++) {if (inverse_map[order[i]] >= 0) filtered_order.push_back(order[i]);}  

  BaseBuffer<int> new_map(index_max, true); vector<int> remap(index_max, -1); 
  for (int i = 0; i < filtered_order.size(); i++) {
    remap[inverse_map[filtered_order[i]]] = i;
    new_map[i] = i+1;
  }
  local_map.swap(new_map);

  ///ensure that external IDs (if any) are remapped to correct iid; unused entries are all set to -1
  if (filtered_order.size() < index_max) {for (int i = 0; i < index_max; i++) {if (remap[i] < 0) remap[i] = -1;}}
  remap_external(remap);

  count_used();
  for (int i = 0; i < variables.size(); i++) { ///put the stored data in the correct order, filter out unused
    if (variables[i]->initialised()) variables[i]->reorder(filtered_order);  
    else variables[i]->set_size(storage_req);  
  }
}

///remove unused entries from storage (unless reduction in needed storage is smaller than 1 - threshold)
///storage is now consecutive, iids remain valid but do not necessarily map directly onto storage
void DataFrame::shrink_storage(float threshold) {count_used();
  if (local_map.is_empty() || index_used >= index_max || (threshold < 1 && storage_req > threshold*max_storage())) return;

  bool do_remap = (index_used != storage_req); vector<int> remap; 
  if (do_remap) {remap.reserve(index_used); 
    for (int iid = 0, sid = 1; iid < index_max; iid++) {
      if (local_map[iid]) {
        remap.push_back(local_map[iid]-1); 
        local_map[iid] = sid++;
      } 
    }
    storage_req = remap.size();
  }
  
  for (int i = 0; i < variables.size(); i++) {
    if (do_remap && variables[i]->initialised()) variables[i]->shrink(remap);  
    else variables[i]->set_size(storage_req);  
  }    
}

///checks if storage is consecutive and iids map directly onto storage
bool DataFrame::is_consecutive() {
  if (!local_map.is_empty()) {
    for (int i = 0; i < index_max; i++) {
      if (local_map[i] > 0 && (local_map[i] != i+1 || (i > 0 && local_map[i-1] == 0))) return false;
    }
  }
  return true;
}          

///ensures that storage is consecutive and iids map directly onto storage (invalidates iids)
bool DataFrame::make_consecutive() {bool do_remap = false; 
  if (!local_map.is_empty()) {
    vector<short> index(index_max, false);  
    for (int i = 0; i < index_max; i++) {
      if (local_map[i] > 0) {
        if (local_map[i] != i+1 || (i > 0 && local_map[i-1] == 0)) do_remap = true;
        index[local_map[i]-1] = true;
      }  
    }
    
    if (do_remap) {
      vector<int> remap; remap.reserve(index_max);
      for (int i = 0; i < index_max; i++) {if (index[i]) remap.push_back(i);}
      reorder(remap);
    }
  }
  return do_remap;
}

void DataFrame::filter_ids(vector<short>& filter, bool include) {
  if (filter.size() < index_max) filter.resize(index_max, 0);
  
  if (include) for (int i = 0; i < index_max; i++) {if (!filter[i]) local_map[i] = 0;}
  else for (int i = 0; i < index_max; i++) {if (filter[i]) local_map[i] = 0;}  
  count_used();
}



int ExpandingDataFrame::add_row() {
  if (!index_built) {
    id_index.push_back(id_index.size());
    return id_index.size() - 1;

  } else return -1;
}

int ExpandingDataFrame::add_block(int size) {
  if (!index_built) {
    int offset = id_index.size();
    for (int i = 0; i < size; i++) id_index.push_back(offset + i); 
    return offset;
  } else return -1;
}

bool ExpandingDataFrame::init(int exp_amount) {
  if (!index_built) {
    id_index.reserve(exp_amount);
    init_variables(exp_amount);
    return true;
  }
  else return false;
}

void ExpandingDataFrame::lock_index(bool store_eids) {
  if (!index_built) {
    index_used = index_max = id_index.size();
    init_storage(index_max);
    index_built = true;

    if (store_eids) {
      input_id.init();
      for (int i = 0; i < id_index.size(); i++) input_id.set(id_index[i], i);
    }
  }
}



bool IndexedDataFrame::init(int exp_amount, int exp_length) {
  if (!index_built && !index_available) {
    id_index->init(exp_amount, exp_length); 
    init_variables(exp_amount);
    return true;
  }
  else return false;
}

void IndexedDataFrame::clear_data(bool clear_ids) {
  clear();
  id_index->clear(clear_ids);
  index_available = false;
}

void IndexedDataFrame::clear_index(bool keep_names) {
  if (index_built) {  
    id_index->clear_tree();
    shrink_storage(SHRINK_THRESHOLD);
    if (keep_names) store_idnames();
    id_index->clear(!name.initialised());    
  }
  index_available = false;
}

void IndexedDataFrame::build_index() {
  if (!index_built) {duplicates.clear();
    index_used = id_index->build_index(&duplicates);
    id_index->build_map();
    index_max = id_index->size(true);
    init_storage(index_max);
    index_built = index_available = true;
    if (name_var) store_idnames();
  }
}

int IndexedDataFrame::extend_index(int amount) {
  if (!index_built) {_LOG.error("processing IndexedDataFrame object") << "cannot extend DataFrame until index has been built" << endl; die();}
  if (amount > 0) {
    int from = index_max, to = index_max + amount;
    count_used(); local_map.resize(to); 
    
    for (int i = from; i < to; i++) local_map[i] = storage_req++ + 1;
    for (int i = 0; i < variables.size(); i++) variables[i]->set_map(storage_req, &local_map); 
    index_max += amount;
    return from;
  } else return 0;
}

bool IndexedDataFrame::ensure_index() {
  if (!index_built) build_index();
  return index_available;  
}

void IndexedDataFrame::store_idnames() {
  if (name.initialised()) return;

  if (!index_available) {_LOG.error("processing IndexedDataFrame object") << "ID name variable is no longer available" << endl; die();}
  BaseBuffer<char*> id_names(index_max, true); id_index->load_ids(id_names);
  name.init(); name.get_data().assign_data(id_names);
  name.filter();
}


///aligns the data with the master DataFrame, such that storage is consecutive and in the same order as in master; entries not in master are discarded
///a map of master iids directly onto storage is created in target
///invalidates iids; if retained, master_iid variable may be invalidated if it is reordered
void IndexedDataFrame::align_map(IndexedDataFrame& master, int* target, int size, bool keep_mid) {
  if (!index_available || !master.index_available || local_map.is_empty() || !name.initialised()) {_LOG.error("processing IndexedDataFrame object") << "unable to align with external map" << endl; die();} 

  NumericID<int>* master_iid = new NumericID<int>(); ///iid in master corresponding to entries in this data frame
  if (keep_mid) attach_variable(master_iid, "MASTER_IID"); else link_variable(*master_iid);
  for (int iid = 0; iid < index_max; iid++) {
    if (valid_id(iid)) {
      int m_iid = master[name.get(iid)];
      if (m_iid >= 0) master_iid->set(iid, m_iid);        
      else drop_id(iid); 
    }
  }
  data_sort(*master_iid); ///reorder data to coincide with master, drop_id calls ensure removal entries not in master

  Utils::fill_value(target, size, -1);
  BaseBuffer<int>& ids = master_iid->get_buffer();  
  for (int iid = 0; iid < index_max; iid++) {
    int sid = local_map[iid] - 1;
    if (valid_id(iid)) target[ids[sid]] = sid;    
  } 

  if (!keep_mid) {unlink_variable(*master_iid); delete master_iid;}
}

int VariableIndex::count_valid() {int count = 0;
  for (int iid = 0; iid < index_max; iid++) {
    if (valid_id(iid) && status.get(iid) == Valid) count++;
  }
  return count;
}


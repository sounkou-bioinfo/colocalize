/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

template<typename T>
void DataFrame::data_sort(UniVariable<T>& var) {
  if (registered(var) < 0) {_LOG.error("processing DataFrame object") << "cannot sort on unregistered variable" << endl; die();} 
  if (var.initialised()) {
    BaseBuffer<T>& data = var.get_buffer();
    reorder(GenericSorter<T>().run(data.data(), data.size()));
  }
}

template<typename T> 
void DataFrame::store_rows(NumericID<T>& target) {
  if (local_map.is_empty()) {_LOG.error("processing DataFrame object") << "unable to store DataFrame row indices without valid local_map" << endl; die();}
  if (registered(target) < 0) link_variable(target);
  for (int i = 0; i < index_max; i++) {if (local_map[i]) target.set(i,i);}
}

template<typename VAR>
void CovariateIndex<VAR>::build_index() {
  if (!index_built) {
    VariableIndex::build_index();
    data.set_width(index_max); 
  }
}

template<typename VAR>
int CovariateIndex<VAR>::extend_index(int amount) {
  int offset = VariableIndex::extend_index(amount);
  if (offset > 0) data.set_width(index_max);    
  return offset;
}    


template<typename VAR>
void CovariateIndex<VAR>::link_data(DataFrame& df, bool init_zero) {
  df.link_variable(data, data_id); 
  df.link_variable(in_source); 
  df.link_variable(miss_count);
  if (init_zero) data.get_buffer().assign_zero();
}

template<typename VAR>
void CovariateIndex<VAR>::unlink_data(bool do_delete) {
  storage_column.assign_missing();
  if (data.get_owner()) {
    DataFrame& df = *data.get_owner();
    df.unlink_variable(data);     
    df.unlink_variable(in_source);
    df.unlink_variable(miss_count);
  }
  if (do_delete) {
    data.clear_data(); 
    in_source.clear_data();
    miss_count.clear_data();
  }
}





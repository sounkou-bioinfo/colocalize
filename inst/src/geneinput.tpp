/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

template<typename VAR>
VAR* GeneInputCombine::get_variable(const string& name, int index) {
  return input_data[index]->using_variable(name) ? input_data[index]->typed_variable<VAR>(name) : 0;
}  

template<typename T>
bool GeneInputCombine::check_variable(UniVariable<T>&, const string& name) {
  for (int i = 0; i < input_data.size(); i++) {if (get_variable<UniVariable<T> >(name, i) == 0) return false;}
  return true;
}

template<typename VAR>
bool GeneInputCombine::check_variable(BlockMultiVariable<VAR>&, const string& name, bool check_width) {int width = 0;
  for (int i = 0; i < input_data.size(); i++) {
    BlockMultiVariable<VAR>* curr = get_variable<BlockMultiVariable<VAR> >(name, i);
    if (!curr) return false;
    if (check_width) {
      if (width == 0) width = curr->get_width();
      else if (curr->get_width() != width) return false;
    }
  }  
  return true;
}

template<typename VAR>
bool GeneInputCombine::check_variable(JaggedMultiVariable<VAR>&, const string& name) {
  for (int i = 0; i < input_data.size(); i++) {if (get_variable<JaggedMultiVariable<VAR> >(name, i) == 0) return false;}
  return true;
}

template<typename T>
vector<BaseBuffer<T>*> GeneInputCombine::get_data(UniVariable<T>&, const string& name) {vector<BaseBuffer<T>*> data;
  for (int i = 0; i < input_data.size(); i++) {
    UniVariable<T>* curr = get_variable<UniVariable<T> >(name, i);
    if (curr) data.push_back(&curr->get_buffer());
  }
  if (data.size() != input_data.size()) data.clear();
  return data;
}

template<typename VAR>
vector<Buffer<typename VAR::VTYPE>*> GeneInputCombine::get_data(BlockMultiVariable<VAR>&, const string& name, bool check_width) {vector<Buffer<typename VAR::VTYPE>*> data;
  for (int i = 0; i < input_data.size(); i++) {
    BlockMultiVariable<VAR>* curr = get_variable<BlockMultiVariable<VAR> >(name, i);
    if (curr && (!check_width || data.empty() || data.back()->ncol() == curr->get_buffer().ncol())) data.push_back(&curr->get_buffer());
  }
  if (data.size() != input_data.size()) data.clear();
  return data;
}

template<typename VAR>
vector<BaseBuffer<BaseBuffer<typename VAR::VTYPE>*>*> GeneInputCombine::get_data(JaggedMultiVariable<VAR>&, const string& name) {vector<BaseBuffer<BaseBuffer<typename VAR::VTYPE>*>*> data;
  for (int i = 0; i < input_data.size(); i++) {
    JaggedMultiVariable<VAR>* curr = get_variable<JaggedMultiVariable<VAR> >(name, i);
    if (curr) data.push_back(&curr->get_buffer());
  }
  if (data.size() != input_data.size()) data.clear();
  return data;
}

template<typename D>
bool GeneInput::Aggregator<D>::set_data(vector<D*>& input) {data.clear(); 
  if (input.size() == size) {
    for (int i = 0; i < size; i++) {if (input[i]) data.push_back(input[i]);}
    eff_size = data.size(); return eff_size > 0;
  } else return false;
}

template<typename D>
bool GeneInput::Aggregator<D>::set_data(vector<BaseBuffer<D>*>& input) {vector<D*> buffers;
  for (int i = 0; i < input.size(); i++) {buffers.push_back((input[i] && !input[i]->is_empty()) ? input[i]->data() : 0);}
  return set_data(buffers);
}

template<typename D>
bool GeneInput::WeightedAggregator<D>::set_data(vector<D*>& input) {this->data.clear(); weights.clear(); 
  if (input.size() == this->size && !weight_base.empty()) {
    for (int i = 0; i < this->size; i++) {
      if (input[i]) {this->data.push_back(input[i]); weights.push_back(weight_base[i]);}
    }
    this->eff_size = this->data.size(); 
    if (this->eff_size > 0) {
      if (corr_base) {
        corrs.set_size(this->eff_size, this->eff_size, true);
        vector<int> index;
        for (int i = 0; i < this->size; i++) {if (input[i]) index.push_back(i);}
        for (int i = 0; i < this->eff_size; i++) {
          corrs(i,i) = (*corr_base)(index[i],index[i]);
          for (int j = i+1; j < this->eff_size; j++) corrs(i,j) = corrs(j,i) = (*corr_base)(index[i],index[j]);
        }
      }
      return true;
    } else return false;
  } else return false;
}

template<typename D>
bool GeneInput::WeightedAggregator<D>::set_weights(vector<double*>& input) {weight_base.clear(); 
  if (input.size() == this->size) {weight_base = input; return true;}
  else return false;
}

template<typename D>
bool GeneInput::WeightedAggregator<D>::set_corrs(Buffer<double>* corrs) {corr_base = 0;
  if (corrs && !corrs->empty()) {
    if (corrs->ncol() != this->size || corrs->nrow() != this->size) return false;
    corr_base = corrs;
  }
  return true;
}

template<typename D>
bool GeneInput::MergeAggregator<D>::compute(int row) {
  for (int col = 0; col < this->eff_size; col++) {
    if (this->access_map(row,col) >= 0) {this->value = this->data[col][this->access_map(row,col)]; return true;}
  }
  return false;    
}

template<typename D>
bool GeneInput::SumAggregator<D>::compute(int row) {this->value = 0; int used = 0;
  for (int col = 0; col < this->eff_size; col++) {
    if (this->access_map(row,col) >= 0) {this->value += this->data[col][this->access_map(row,col)]; used++;}
  }
  return used > 0;    
}

template<typename D>
bool GeneInput::MeanAggregator<D>::compute(int row) {this->value = 0; int used = 0;
  for (int col = 0; col < this->eff_size; col++) {
    if (this->access_map(row,col) >= 0) {this->value += this->data[col][this->access_map(row,col)]; used++;}
  }
  if (used > 1) {
    if (numeric_limits<D>::is_integer) this->value = integer_round(this->value/double(used));
    else this->value /= used;
  }
  return used > 0;    
}

template<typename D>
double GeneInput::MeanAggregator<D>::integer_round(double value) {
  if (round_type == 1) return ceil(value);
  else if (round_type == -1) return floor(value);
  else return round(value);
}

template<typename D>
bool GeneInput::WeightedZstatAggregator<D>::compute(int row) {this->value = 0; double weight_sum = 0; 
  Buffer<double> wts(this->eff_size, 1, true);
  for (int col = 0; col < this->eff_size; col++) {
    if (this->access_map(row,col) >= 0) {
      int& index = this->access_map(row,col); double& wt = this->weights[col][index];
      if (wt > 0) {
        this->value += this->truncate(this->data[col][index]) * wt; weight_sum += wt*wt;
        if (this->corr_base) {
          wts(col) = wt;
          for (int col2 = 0; col2 < col; col2++) weight_sum += 2*wt*wts(col2)*this->corrs(col, col2);
        }
      }
    }
  }
  if (weight_sum > 0) {this->value /= sqrt(weight_sum); return true;}
  else return false;
}

template<typename D>
bool GeneInput::WeightedPvalAggregator<D>::compute(int row) {try {this->value = 0; double weight_sum = 0; 
  Buffer<double> wts(this->eff_size, 1, true);
  for (int col = 0; col < this->eff_size; col++) {
    if (this->access_map(row,col) >= 0) {
      int& index = this->access_map(row,col); double& wt = this->weights[col][index]; D& p = this->data[col][index];
      if (wt > 0 && p >= 0 && p <= 1) {
        this->value += this->truncate(from_pval.convert(p)) * wt; weight_sum += wt*wt;
        if (this->corr_base) {
          wts(col) = wt;
          for (int col2 = 0; col2 < col; col2++) weight_sum += 2*wt*wts(col2)*this->corrs(col, col2);
        }
      }
    }
  }
  if (weight_sum > 0) {this->value = to_pval.convert(this->value / sqrt(weight_sum)); return true;}
  else return false;
} catch (const StatException& se) {return false;}}

template<typename D>
GeneInput::Aggregator<D>* GeneInputCombine::make_aggregator(Buffer<int>& access_map, GeneInput::Mode mode) {using namespace GeneInput;
  if (GeneInput::is_weighted(mode)) return make_aggregator<D>(access_map, get_weights(), get_data_corrs(), mode);
  switch (mode) {
    case Merge: return new MergeAggregator<D>(access_map);
    case Sum: return new SumAggregator<D>(access_map);
    case Mean: return new MeanAggregator<D>(access_map);    
    case CeilingMean: return new MeanAggregator<D>(access_map, 1);        
    case FloorMean: return new MeanAggregator<D>(access_map, -1);            
    default: return 0;  
  }
}

template<typename D>
GeneInput::WeightedAggregator<D>* GeneInputCombine::make_aggregator(Buffer<int>& access_map, vector<double*> weights, Buffer<double>* corrs, GeneInput::Mode mode) {using namespace GeneInput;
  WeightedAggregator<D>* out = 0;
  switch (mode) {
    case WeightedZstat: out = new WeightedZstatAggregator<D>(access_map, get_truncate()); break;
    case WeightedPval: out = new WeightedPvalAggregator<D>(access_map, get_truncate()); break;    
  }
  if (out && out->set_weights(weights) && out->set_corrs(corrs)) return out;
  else {delete out; return 0;}
}

template<typename D, typename DATA>
GeneInput::Aggregator<D>* GeneInputCombine::make_aggregator(Buffer<int>& access_map, DATA input, GeneInput::Mode mode) {
  GeneInput::Aggregator<D>* aggregator = make_aggregator<D>(access_map, mode);
  if (aggregator && aggregator->set_data(input)) return aggregator;
  else {delete aggregator; return 0;}
}

template<typename T> 
bool GeneInputCombine::process_variable(UniVariable<T>& target, const string& name, GeneInput::Mode mode) {
  GeneInput::Aggregator<T>* aggregator = make_aggregator<T>(access_map, get_data(target, name), mode);
  if (aggregator) {
    target.init();
    for (int i = 0; i < gene_info.data_used(false); i++) {if (aggregator->compute(i)) target.set(i, aggregator->get_value(), true);}
    delete aggregator; return true;  
  } else return false; 
}

template<typename VAR> 
bool GeneInputCombine::process_variable(BlockMultiVariable<VAR>& target, const string& name, GeneInput::Mode mode) {typedef typename VAR::VTYPE CTYPE;
  vector<Buffer<CTYPE>*> buffers = get_data(target, name);
  if (!buffers.empty()) {
    int size = buffers[0]->ncol(); if (size == 0) return false;
    vector<vector<CTYPE*> > data(size);
    for (int i = 0; i < buffers.size(); i++) {Buffer<CTYPE>& curr = *buffers[i];  
      if (curr.ncol() != size) return false;
      for (int j = 0; j < curr.ncol(); j++) data[j].push_back(curr[j]);
    }
    return process_variable(target, data, mode);
  } else return false; 
}

template<typename VAR> 
bool GeneInputCombine::process_aligned_variable(BlockMultiVariable<VAR>& target, const string& name, GeneInput::Mode mode, bool align_main, bool skip_root) {typedef typename VAR::VTYPE CTYPE;
  vector<Buffer<CTYPE>*> buffers = get_data(target, name, false);
  if (!buffers.empty()) {
    ModelInfoBlock& info = gene_info.get_model_info(); 
    int size = align_main ? info.main_index().size() : info.partition_index(skip_root).size(); if (size == 0) return false; 
    vector<vector<CTYPE*> > data(size); 
    
    for (int i = 0; i < buffers.size(); i++) {vector<int> alignment;  
      if (align_main) alignment = input_data[i]->get_model_info().copy().align(info).main_index();
      else alignment = input_data[i]->get_model_info().copy().align(info).partition_index(skip_root);
      if (alignment.size() != size) return false;

      Buffer<CTYPE>& curr = *buffers[i];      
      for (int j = 0; j < size; j++) data[j].push_back(alignment[j] >= 0 ? curr[alignment[j]] : 0);
    }
    return process_variable(target, data, mode);
  } else return false; 
}


template<typename VAR> 
bool GeneInputCombine::process_variable(BlockMultiVariable<VAR>& target, vector<vector<typename VAR::VTYPE*> >& data, GeneInput::Mode mode) {using namespace GeneInput; typedef typename VAR::VTYPE CTYPE;
  Aggregator<CTYPE>* aggregator = make_aggregator<CTYPE>(access_map, mode);  
  if (aggregator) {
    int no_var = data.size(); target.init(no_var);
    for (int var = 0; var < no_var; var++) {
      if (aggregator->set_data(data[var])) {
        for (int i = 0; i < gene_info.data_used(false); i++) {
          if (aggregator->compute(i)) target.set(i, var, aggregator->get_value(), true);
        }
      }    
    }
    delete aggregator; return true;  
  } else return false;
}


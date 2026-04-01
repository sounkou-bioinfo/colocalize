/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "geneutils.h"
#include <map>

string ModelState::model_name(Model model) {
  switch (model) {
    case MultiModel: return "MULTI";
    case SNPwiseModelMeanBrown: return "SNPWISE_MEAN_BROWN";
    case SNPwiseModelMeanImhof: return "SNPWISE_MEAN";    
    case SNPwiseModelTop: return "SNPWISE_TOP";
    case SNPwiseModelTopOne: return "SNPWISE_TOP1";
    case PrinCompReg: return "PCREG";
    default: return "#UNKNOWN#";
  }
}

string ModelState::analysis_name(Analysis analysis) {
  switch (analysis) {
    case FullAnalysis: return "FULL";
    case InteractionAnalysis: return "INTERACT";
    case MainAnalysis: return "MAIN";
    default: return "#UNKNOWN#";
  }
}

ModelState::Model ModelState::model_type(const string& name) {string NAME = Utils::uppercase(name);
  if (NAME == "MULTI" || NAME == "JOINT") return MultiModel;
  if (NAME == "SNPWISE_MEAN") return SNPwiseModelMeanImhof;   
  if (NAME == "SNPWISE_TOP") return SNPwiseModelTop;
  if (NAME == "SNPWISE_TOP1") return SNPwiseModelTopOne;  
  if (NAME == "PCREG") return PrinCompReg;
  return UnknownModel;
}

ModelState::Analysis ModelState::analysis_type(const string& name) {string NAME = Utils::uppercase(name);
  if (NAME == "FULL") return FullAnalysis; 
  if (NAME == "INTERACT") return InteractionAnalysis;
  if (NAME == "MAIN") return MainAnalysis;   
  return UnknownAnalysis;
}

ModelInfoBlock& ModelInfoBlock::operator=(const ModelInfoBlock& other) {
  if (&other == this) return *this;
  clear(); state = other.state; partition_names = other.partition_names; no_repeat = other.no_repeat;
  for (int i = 0; i < other.info.size(); i++) info.push_back(other.info[i]->copy());
  return *this;
}

bool ModelInfoBlock::operator==(const ModelInfoBlock& other) {
  if (!Utils::vec_equals(state, other.state) || !Utils::vec_equals(partition_names, other.partition_names) || info.size() != other.info.size() || no_repeat != other.no_repeat) return false;
  for (int i = 0; i < info.size(); i++) {if (*info[i] != *other.info[i]) return false;}
  return true;
}

set<ModelState::Model> ModelInfoBlock::get_models() {set<ModelState::Model> out;
  for (int i = 0; i < info.size(); i++) out.insert(info[i]->model);
  return out;
}

set<ModelState::Analysis> ModelInfoBlock::get_analyses() {set<ModelState::Analysis> out;
  for (int i = 0; i < info.size(); i++) out.insert(info[i]->analysis);
  return out;
}

ModelInfoBlock& ModelInfoBlock::align(const ModelInfoBlock& other) {
  if (!Utils::vec_equals(partition_names, other.partition_names)) update_partitions(other.partition_names);
  map<int,int> index; vector<ModelOutputInfo*> new_info; 
  for (int i = 0; i < info.size(); i++) index[info[i]->hash()] = i;
  for (int i = 0; i < other.info.size(); i++) {
    int hash = other.info[i]->hash(); map<int,int>::iterator found = index.find(hash);
    if (found != index.end() && info[found->second] != 0) {new_info.push_back(info[found->second]); info[found->second] = 0;}
    else new_info.push_back(other.info[i]->copy(-1));
  }
  clear(); info.swap(new_info); detect_state(); return *this;  
  return *this;
}

vector<int> ModelInfoBlock::main_index() {
  vector<int> out; vector<pair<int,string> > index = main_index("");
  for (int i = 0; i < index.size(); i++) out.push_back(index[i].first);

  return out;  
}

vector<pair<int,string> > ModelInfoBlock::main_index(const string& prefix) {detect_state(); 
  vector<pair<int,string> > out;
  for (int i = 0; i < info.size(); i++) out.push_back(pair<int,string>(info[i]->index, make_name(prefix, i)));
  if (no_repeat > 1) {
    vector<pair<int,string> > tpl; out.swap(tpl);
    for (int i = 0; i < no_repeat; i++) {
      string num = Utils::num_string(i+1);
      for (int j = 0; j < tpl.size(); j++) {
        pair<int,string> curr = tpl[j];
        curr.first += i*tpl.size();
        curr.second += num; 
        out.push_back(curr);        
      }
    }
  }
  
  return out;  
}

vector<int> ModelInfoBlock::partition_index(bool skip_root) {detect_state();
  vector<int> rank(partition_names.size(), -2), rerank(state[MaxIndex]+1, -1); set<int> index_used; int index_count = 0;
  for (int i = 0; i < info.size(); i++) {if (rank[info[i]->partition] < -1 || rank[info[i]->partition] > info[i]->index) rank[info[i]->partition] = info[i]->index;}
  for (int i = skip_root ? 1 : 0; i < rank.size(); i++) {if (rank[i] >= 0) index_used.insert(rank[i]);}
  for (set<int>::iterator it = index_used.begin(); it != index_used.end(); ++it) rerank[*it] = index_count++;

  set<int> used; vector<int> out; if (skip_root) used.insert(0);
  for (int i = 0; i < info.size(); i++) {int part = info[i]->partition;
    if (used.find(part) != used.end()) continue; used.insert(part); 
    if (rank[part] > -2) out.push_back(rank[part] >= 0 ? rerank[rank[part]] : -1);
  }    

  return out;
}

void ModelInfoBlock::detect_state() {
  state.assign(_STATECOUNT_, 0); int& max_index = state[MaxIndex]; max_index = -1;
  set<ModelState::Model> models, printable_models; 
  set<ModelState::Analysis> analyses; set<int> partitions, hash;
  for (int i = 0; i < info.size(); i++) {
    models.insert(info[i]->model); analyses.insert(info[i]->analysis); 
    partitions.insert(info[i]->partition); hash.insert(info[i]->hash());
    if (ModelState::is_printable(info[i]->model)) printable_models.insert(info[i]->model);        
    if (info[i]->analysis == ModelState::InteractionAnalysis) state[HasInteraction] = true;
    if (info[i]->partition > 0) state[HasPartitions] = true;
    if (info[i]->model == ModelState::MultiModel) state[HasMultiModel] = true;
    if (info[i]->index > max_index) max_index = info[i]->index;
  }  
  state[ModelCount] = models.size(); state[PrintableModelCount] = printable_models.size(); 
  state[AnalysisCount] = analyses.size(); state[PartitionCount] = partitions.size();
  state[IsUnbalanced] = (state[ModelCount] * state[AnalysisCount] * state[PartitionCount] != hash.size()); 
}

bool ModelInfoBlock::is_consistent() {
  if (state[ModelInfoBlock::IsUnbalanced]) return false;
  if (state[ModelInfoBlock::ModelCount] > 1 && !state[ModelInfoBlock::HasMultiModel]) return false;
  if (state[ModelInfoBlock::AnalysisCount] != 1 && state[ModelInfoBlock::AnalysisCount] != 3) return false;
  if (state[ModelInfoBlock::AnalysisCount] == 1 && info[0]->analysis != ModelState::MainAnalysis) return false;
  return true;
}

ModelInfoBlock& ModelInfoBlock::filter(ModelState::Analysis analysis, bool keep) {vector<ModelOutputInfo*> new_info; 
  for (int i = 0; i < info.size(); i++) {if ((info[i]->analysis == analysis) == keep) {new_info.push_back(info[i]); info[i] = 0;}}
  clear(); info.swap(new_info); detect_state(); return *this;
}

ModelInfoBlock& ModelInfoBlock::filter(ModelState::Model model, bool keep) {vector<ModelOutputInfo*> new_info; 
  for (int i = 0; i < info.size(); i++) {if ((info[i]->model == model) == keep) {new_info.push_back(info[i]); info[i] = 0;}}
  clear(); info.swap(new_info); detect_state(); return *this;
}

ModelInfoBlock& ModelInfoBlock::truncate(int max_index) {
  if (info.size() > max_index) {
    for (int i = max_index; i < info.size(); i++) delete info[i];
    info.resize(max_index);
  }
  detect_state(); return *this;
}

ModelInfoBlock& ModelInfoBlock::update_partitions(const vector<string>& partitions) {
  map<string,int> index; vector<ModelOutputInfo*> new_info; 
  for (int i = 1; i < partitions.size(); i++) index[partitions[i]] = i;
  for (int i = 0; i < info.size(); i++) {ModelOutputInfo*& curr = info[i];
    if (curr->partition > 0) {string& name = partition_names[curr->partition];
      if (index.find(name) != index.end()) curr->partition = index[name];     
      else continue;
    }
    new_info.push_back(curr); curr = 0;
  }  
  partition_names = partitions;
  clear(); info.swap(new_info); detect_state(); return *this;
}

ModelInfoBlock& ModelInfoBlock::strip_partitions() {vector<ModelOutputInfo*> new_info; 
  for (int i = 0; i < info.size(); i++) {if (info[i]->partition == 0) {new_info.push_back(info[i]); info[i] = 0;}}
  partition_names.clear(); partition_names.push_back("#UNDEFINED#");
  clear(); info.swap(new_info); detect_state(); return *this;
}


ModelInfoBlock& ModelInfoBlock::set_model(ModelState::Model model) {
  for (int i = 0; i < info.size(); i++) info[i]->model = model;
  detect_state(); return *this;
}

ModelInfoBlock& ModelInfoBlock::set_partition(int partition) {
  for (int i = 0; i < info.size(); i++) info[i]->partition = partition;
  detect_state(); return *this;
}

ModelInfoBlock& ModelInfoBlock::prune_models() {
  ModelState::Model use = state[HasMultiModel] ? ModelState::MultiModel : info[0]->model;
  return filter(use, true).set_model(ModelState::UnknownModel);
}

ModelInfoBlock& ModelInfoBlock::shift_index(int add) {
  for (int i = 0; i < info.size(); i++) info[i]->index += add;
  return *this;
}

ModelInfoBlock& ModelInfoBlock::append(ModelInfoBlock& other) {
  for (int i = 0; i < other.info.size(); i++) info.push_back(other.info[i]->copy());
  detect_state(); return *this;   
}

ModelInfoBlock& ModelInfoBlock::sort(bool reindex) {int max_partition = 0; vector<int> rank;
  for (int i = 0; i < info.size(); i++) {if (info[i]->partition > max_partition) max_partition = info[i]->partition;}
  for (int i = 0; i < max_partition+1; i++) rank.push_back(i+1);
  return sort(rank, reindex);
}                          

ModelInfoBlock& ModelInfoBlock::sort(vector<int> partition_rank, bool reindex) {vector<ModelOutputInfo*> new_info; 
  for (int i = 0; i < info.size(); i++) {if (partition_rank[info[i]->partition] > 0) {new_info.push_back(info[i]); info[i] = 0;}}
  clear(); info.swap(new_info);
  ::sort(info.begin(), info.end(), InfoSort(partition_rank));
  if (reindex) {for (int i = 0; i < info.size(); i++) info[i]->index = i;}
  return *this;
}

string ModelInfoBlock::make_name(const string& prefix, int index, bool print_full) {
  string name = prefix; int added = 0; ModelOutputInfo* curr = info[index];
  if (print_full || state[ModelInfoBlock::HasPartitions]) {added++;
    if (curr->partition > 0) name += "_PART-" + get_partition(curr->partition);
    else name += "_ALL-PART";
  }
  if (print_full || state[ModelInfoBlock::HasInteraction]) {
    if (added++) name += "_"; name += "_";
    name += ModelState::analysis_name(curr->analysis);
  }
  if (print_full || state[ModelInfoBlock::PrintableModelCount] > 1 || state[ModelInfoBlock::HasMultiModel]) {
    if (added++) name += "_"; name += "_";
    name += ModelState::model_name(curr->model);
  }
  return name;
}

bool ModelInfoBlock::InfoSort::operator () (ModelOutputInfo* i, ModelOutputInfo* j) {
  if (i->model != j->model) return i->model < j->model;
  if (i->analysis != j->analysis) return i->analysis < j->analysis;
  if (i->partition != j->partition) return partition_rank[i->partition] < partition_rank[j->partition];
  return i->index <= j->index;
}


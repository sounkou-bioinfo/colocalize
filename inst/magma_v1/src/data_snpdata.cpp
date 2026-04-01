/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "data_snpdata.h"

string SNPData::make_dummy(char* str, char* str2) {string out;
  if (strcmp(str, str2) > 0) str = str2;
  do out.push_back(127);
  while (strcmp(str, out.c_str()) >= 0);
  return out;
}

void SNPData::load_synonyms(IDIndex* synonyms) {
  if (ensure_index()) {
    id_index = new IDIndexSequence(id_index, synonyms); 
    has_synonyms = true;
  }
}

map<int,set<char*> > SNPData::get_synonyms(set<string>& ids) {
  if (ensure_index()) {
    vector<char> index(index_max, 0);    
    for (set<string>::iterator it = ids.begin(); it != ids.end(); ++it) {
      int iid = (*this)[*it];
      if (iid >= 0) index[iid] = 1;
    }
    return id_index->reverse_map(index);  
  } else return map<int,set<char*> >();
}

bool SNPSynonyms::init(int exp_amount, int exp_length) {
  if (IndexedDataFrame::init(exp_amount, exp_length)) {
    if (exp_amount < 1000) exp_amount = 1000;
    id_index->init(exp_amount, exp_length);  
    snpid.set_size(exp_amount); synid.set_size(exp_amount);
    synonyms.reserve(exp_amount/2);
    return true;
  } else return false;
}

void SNPSynonyms::clear_data(bool clear_ids) {
  IndexedDataFrame::clear_data(clear_ids);
  for (int i = 0; i < synonyms.size(); i++) delete synonyms[i];
  synonyms.clear();  
}

int SNPSynonyms::get_count(const string& type) {
  if (type == "aliases") return counts[0];
  else if (type == "targets") return counts[1];  
  else if (type == "synonyms") return counts[2];
  else return -1;
}

void SNPSynonyms::set_mode(const string& mode) {
  if (mode == "drop") duplicate_mode = DupDropAll; 
  else if (mode == "drop-dup") duplicate_mode = DupDropDup; 
  else if (mode == "skip") duplicate_mode = DupSkipLine; 
  else if (mode == "skip-dup") duplicate_mode = DupSkipDup; 
  else duplicate_mode = DupError;  
}

int SNPSynonyms::resolve_id(vector<int>& index, int id) {
  while (index[id] != id) id = index[id];
  return id;
}

void SNPSynonyms::build_index() {
  if (!index_built) {
    index_max = index_used = storage_req = id_index->build_index();
    for (int i = 0; i < variables.size(); i++) variables[i]->set_map(storage_req, 0);
    index_built = index_available = true;
  }
}

vector<vector<char*> > SNPSynonyms::process(SNPData& snp_info) {vector<vector<char*> > duplicates;
  if (offset < added) process_synonym();

  IDIndexCore* index_core = dynamic_cast<IDIndexCore*>(id_index);
  BaseBuffer<char*> snp_names; snp_names.assign(index_core->get_names());
  build_index();

  if (index_max == 0) return duplicates;
  BaseBuffer<int>& map_to_iid = index_core->get_iid_map();
  BaseBuffer<char*>& map_to_names = index_core->get_name_map();  

  VariableStorageAccess::UniAccessor<int>* syn_ids = synid.get_data().get();  
  vector<int> syn_index; syn_index.reserve(synonyms.size());
  for (int i = 0; i < synonyms.size(); i++) syn_index[i] = i;
  for (int i = 0; i < index_max; i++) {
    if (i && IDSorter::idcmp(map_to_names[i], map_to_names[i-1])) {
      int iid1 = map_to_iid[i-1], iid2 = map_to_iid[i], sid1 = resolve_id(syn_index, syn_ids->elem(iid1)), sid2 = resolve_id(syn_index, syn_ids->elem(iid2));
      if (iid2 < iid1) {swap(iid1, iid2); swap(sid1, sid2);}
      synonyms[sid2]->delete_elem(iid2);
      if (sid1 != sid2) {
        synonyms[sid1] = synonyms[sid1]->merge(synonyms[sid2]);
        synonyms[sid2] = 0; syn_index[syn_ids->elem(iid2)] = sid1;      
      } 
    }
  }
  synid.clear_data(); syn_index.clear();

  int alias_count = 0, target_count = 0;
  VariableStorageAccess::UniAccessor<int>* snp_ids = snpid.get_data().get();  
  BaseBuffer<int> syn_target(index_max); Utils::fill_value(syn_target.data(), index_max, -1);
  for (int i = 0; i < synonyms.size(); i++) {SynonymRangeIterator syn(synonyms[i]); 
    if (!syn.empty()) {vector<int> mapped, unmapped;
      while(syn.next()) {
        if (snp_ids->elem(syn.value) >= 0) mapped.push_back(syn.value);
        else unmapped.push_back(syn.value);        
      }

      if (!mapped.empty()) {
        alias_count += unmapped.size(); target_count++;
        int target = snp_ids->elem(mapped[0]);    
        if (mapped.size() > 1) {duplicates.push_back(vector<char*>());
          for (int i = 0; i < mapped.size(); i++) duplicates.back().push_back(snp_names[mapped[i]]);
          if (duplicate_mode == DupError || duplicate_mode == DupSkipLine) continue;
          else if (duplicate_mode == DupDropAll) {for (int i = 0; i < mapped.size(); i++) snp_info.drop_id(snp_ids->elem(mapped[i])); continue;}
          else if (duplicate_mode == DupDropDup) for (int i = 1; i < mapped.size(); i++) snp_info.remap_id(snp_names[mapped[i]], target);
        } 
        for (int i = 0; i < unmapped.size(); i++) syn_target[unmapped[i]] = target;
      } 
    }
  } 
  for (int i = 0; i < index_max; i++) map_to_iid[i] = syn_target[map_to_iid[i]];

  counts[0] = alias_count; counts[1] = target_count; counts[2] = synonyms.size();
  snpid.clear_data(); snp_names.clear(); syn_target.clear(); 

  IDIndexExternal* syn_remap = new IDIndexExternal(map_to_names.data(), map_to_iid.data(), index_max);
  syn_remap->init(alias_count, id_index->id_size());
  syn_remap->process(true, false);
  clear_data(false);

  if (!(duplicate_mode == DupError && !duplicates.empty())) snp_info.load_synonyms(syn_remap); 
  return duplicates;
}

void SNPSynonyms::process_synonym() {
  if (added - offset > 1) {
    int current_syn = synonyms.size();  
    synonyms.push_back(new SynonymRange(offset, added));
    for (; offset < added; offset++) synid.set(offset, current_syn);
  } else offset = added;
}

SNPSynonyms::SynonymRange* SNPSynonyms::SynonymRange::merge(SynonymRange* other) {if (!other) return this;
  if (from <= other->from) {
    if (other->from < to) other->from = to + 1;
    if (other->to <= other->from) {
      SynonymRange* other_child = other->child; other->child = 0; delete other; 
      return (other_child) ? merge(other_child) : this;
    }
    if (child) child = child->merge(other);
    else child = other;
    return this;
  } else return other->merge(this);
}

void SNPSynonyms::SynonymRange::delete_elem(int pos) {
  if (pos < from) return;
  else if (pos == from) from++;
  else if (pos == to-1) to--; 
  else if (pos >= to) {if (child) child->delete_elem(pos);}
  else {child = new SynonymRange(pos+1, to, child); to = pos;} 
}

bool SNPSynonyms::SynonymRangeIterator::empty() {
  if (!source) return false;
  int range = source->to - source->from;
  if (range > 1) return false;
  
  for (SynonymRange* curr = source->child; range < 2 && curr; curr = curr->child) range += curr->to - curr->from;
  return range < 2;
}

bool SNPSynonyms::SynonymRangeIterator::next() {
  if (source) {
    if (++value < source->to) return true;
    while ( (source = source->child) ) {
      if (source->from < source->to) {value = source->from; return true;}
    }
  }
  return false;
}



/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "_global.h"
#include "utils_dataframe_ids.h"

#define NODE_SPLIT 10
#define LINEAR_CUTOFF 10
#define HASH_MAX 10

bool IDSorter::idcmp(char* s1, char* s2) {char l1 = *(s1-1), l2 = *(s2-1);
  if (l1 != l2) return false;
  if (l1 == 127) return strcmp(s1, s2) == 0;
  for (char *r1 = s1+l1-1, *r2 = s2+l2-1; r1 >= s1; r1--, r2--) {if (*r1 != *r2) return false;}
  return true;
}

int IDSorter::run(BaseBuffer<char*>& id_buff, BaseBuffer<int>& index_buff, int size, bool drop_duplicates, vector<char*>* store_duplicates) {
  vector<pair<char*,int> > pairs; pairs.reserve(size);
  for (int i = 0; i < size; i++) pairs.push_back(pair<char*,int>(id_buff[i], index_buff[i]));
  sort(pairs.begin(), pairs.end(), SortObj());

  if (drop_duplicates) {int write = 0;
    char** id_array = id_buff.data(); int* index_array = index_buff.data();
    for (int i = 0; i < size; i++) {
      if (i && idcmp(pairs[i].first, pairs[i-1].first)) {
        if (store_duplicates) {if (store_duplicates->empty() || !idcmp(store_duplicates->back(), pairs[i].first)) store_duplicates->push_back(pairs[i].first);}
        continue;
      }
      id_array[write] = pairs[i].first;
      index_array[write++] = pairs[i].second;      
    }
    size = write;
  } else for (int i = 0; i < size; i++) {id_buff[i] = pairs[i].first; index_buff[i] = pairs[i].second;}
  id_buff.resize(size); index_buff.resize(size);

  return size;  
}

void IDTree::clear() {
  delete root; root = 0;
  leaf_index.clear(); leaf_depth.clear();
  id_count = 0; tree_status = false; ids = 0;
}

void IDTree::build(BaseBuffer<char*>& id_buff) {
  ids = id_buff.data(); id_count = id_buff.size(); tree_status = true;

  if (id_count > 0) {
    char *first = ids[0], *last = ids[id_count-1]; int offset = 0;
    while(first[offset] && first[offset] == last[offset]) offset++; 
    
    leaf_insert(0); 
    root = partition(0, id_count, offset);
    if (offset > 0) root = new RootNode(string(first, offset), root); 
    leaf_insert(id_count);
  } else root = new DeadNode();
}

int IDTree::operator[](char* id) {int index = root->run(id);
  if (index == 0) return -1;
  int offset = leaf_index[index], length = leaf_index[index+1] - offset;
  if (length == 0) return -1;
  int depth = leaf_depth[index]; id += depth; int comp = strcmp(id, ids[offset]+depth);
  if (comp == 0) return offset; if (comp < 0 || length == 1) return -1; 
  int last = offset + length - 1;
  comp = strcmp(id, ids[last]+depth);
  if (comp == 0) return last; if (comp > 0 || length == 2) return -1;   

  offset++; length -= 2;
  while (length > LINEAR_CUTOFF) {
    int target = offset + length/2;
    comp = strcmp(id, ids[target]+depth);
    if (comp == 0) return target;
    if (comp > 0) {offset = target + 1; length -= length/2 + 1;}
    else length /= 2;
  }

  last = offset + length;
  for (int i = offset; i < last; i++) {
    comp = strcmp(id, ids[i]+depth);
    if (comp < 0) return -1;
    if (comp == 0) return i;
  }

  return -1;
}

IDNode* IDTree::partition(int from, int to, unsigned char level) {
  if (from < to) {
    vector<int> location; vector<char> value; char curr_char = ids[from][level]; 
    for (char **current = ids + from, **previous = current, **end = ids + to; current <= end; current++) {
      if (current == end || (*current)[level] != curr_char) {
        value.push_back(curr_char); location.push_back(previous - ids); 
        if (current < end) {previous = current; curr_char = (*current)[level];}
      }
    }    
    location.push_back(to);
  
    return make_node(location, value, level+1);
  } else return new DeadNode();
}

IDNode* IDTree::make_node(vector<int>& location, vector<char>& value, unsigned char level) {char first = value.front(), last = value.back();
  if (value.size() > 1 && !(first < 0 && last < 0) && !(first > 0 && last > 0))  {
    IDNode *pos_node = 0, *neg_node = 0; int null_index = 0, prev_index = 0, prev_sign = (first >= 0 ? (first > 0 ? 1 : 0) : -1), curr_sign = prev_sign;
    for (int i = 1; i <= value.size(); i++) {
      if (i < value.size()) curr_sign = (value[i] >= 0 ? (value[i] > 0 ? 1 : 0) : -1);      
      if (i == value.size() || curr_sign != prev_sign) {
        if (prev_sign != 0) {
          IDNode* node = make_node(location, value, level, prev_index, i - prev_index);
          if (prev_sign == -1) neg_node = node; else pos_node = node;
        } else null_index = leaf_insert(location[prev_index], '\0', level);
        prev_sign = curr_sign; prev_index = i;     
      }
    }
    if (neg_node && pos_node) return new PosNegNode(neg_node, pos_node, null_index);
    else if (neg_node) return new NegNode(neg_node, null_index);
    else return new PosNode(pos_node, null_index);
  } else if (value.size() == 1 && first == '\0') return new StubLeafNode('\0', leaf_insert(location[0], '\0', level));
  else return make_node(location, value, level, 0, value.size());
}

IDNode* IDTree::make_node(vector<int>& location, vector<char>& value, unsigned char level, int offset, int total) {
  int max_count = 0, last = offset + total - 1;
  for (int i = offset; i <= last; i++) {if (location[i+1] - location[i] > leaf_max) max_count++;} 
  char lower = value[offset], upper = value[last];

  if (total > 1) {
    if (upper - lower > total + NODE_SPLIT - 2) {int split = 0;
      for (int i = offset+1; i <= last; i++) {if (value[i] - value[i-1] >= NODE_SPLIT) {split = i; break;}}
      if (split) {
        IDNode* node1 = make_node(location, value, level, offset, split - offset); 
        IDNode* node2 = make_node(location, value, level, split, last - split + 1);
        return new SplitNode(value[split], node1, node2);          
      }
    }
      
    if (max_count > 0) {
      BranchNode* branch = new BranchNode(lower, upper);
      for (int i = offset; i <= last; i++) tree_status = branch->add_child(value[i], partition(location[i], location[i+1], level));
      return branch;
    } else return new LeafNode(lower, upper, leaf_insert(location, value, level, offset, total));          
  } else {
    if (max_count > 0) return new TrunkNode(value[offset], partition(location[offset], location[offset+1], level));
    else return new StubLeafNode(value[offset], leaf_insert(location[offset], value[offset], level));  
  }
}

int IDTree::leaf_insert(int location, char value, unsigned char depth) {
  leaf_index.push_back(location); leaf_depth.push_back(value == '\0' ? depth-1 : depth);
  return leaf_index.size() - 1;  
}

int IDTree::leaf_insert(vector<int>& location, vector<char>& value, unsigned char depth, int offset, int total) {
  int index = leaf_insert(location[offset], value[offset], depth);  
  for (int i = offset+1, end = offset+total; i < end; i++) {
    unsigned char d = (value[i] == '\0') ? depth-1 : depth;
    for (char c = value[i-1]; c < value[i]; c++) {
      leaf_index.push_back(location[i]); 
      leaf_depth.push_back(d);
    }
  }
  return index;
}

int BranchNode::run(char* id) {
  if (*id < lower || *id > upper || !(children[*id - lower])) return 0;
  return children[*id - lower]->run(id+1);
}

bool BranchNode::add_child(char value, IDNode* child) {
  if (value < lower || value > upper) return false;
  if (children[value-lower]) delete children[value-lower];
  children[value-lower] = child;
  return true;
}


void IDIndex::init(long amount, int exp_length) {
  expected_elem = amount;
  if (exp_length > 0) default_id_size = exp_length;
}
void IDIndex::clear(bool clear_ids) {UNUSED(clear_ids); has_index = false; has_map = false;}
void IDIndex::clear_tree() {has_map = false;}

void IDIndex::set_param(const string& name, long value) {
  if (name == "drop_duplicates") drop_duplicates = value;
  else if (name == "id_size") default_id_size = value;
  else if (name == "tree_depth") tree_depth = value;
  else {_LOG.error("processing IDIndex object") << "trying to set unknown parameter '" << name << "'" << endl; die();}
}

void IDIndex::merge_maps(map<int,set<char*> >& base, map<int,set<char*> >& add) {map<int,set<char*> >::iterator found;
  for (map<int,set<char*> >::iterator it = add.begin(); it != add.end(); ++it) {set<char*>& ids = it->second;
    if ((found = base.find(it->first)) != base.end()) found->second.insert(ids.begin(), ids.end());
    else base[it->first] = ids;
  }
}

void IDIndexSequence::init(long amount, int exp_length) {primary->init(amount, exp_length); if (secondary) secondary->init(amount, exp_length);}
void IDIndexSequence::clear(bool clear_ids) {IDIndex::clear(clear_ids); primary->clear(clear_ids); if (secondary) secondary->clear(clear_ids);}
void IDIndexSequence::clear_tree() {IDIndex::clear_tree(); primary->clear_tree(); if (secondary) secondary->clear_tree();}
void IDIndexSequence::set_param(const string& name, long value) {primary->set_param(name, value); if (secondary) secondary->set_param(name, value);}

int IDIndexSequence::build_index(vector<char*>* store_duplicates) {has_index = true;
  if (secondary) secondary->build_index();
  return primary->build_index(store_duplicates);
}

void IDIndexSequence::build_map() {has_map = true;
  primary->build_map();
  if (secondary) secondary->build_map();
}

map<int,set<char*> > IDIndexSequence::reverse_map(vector<char>& used) {
  map<int,set<char*> > out = primary->reverse_map(used);
  if (secondary) {
    map<int,set<char*> > add = secondary->reverse_map(used);
    merge_maps(out, add);  
  }
  return out;
}

void IDIndexSequence::load_ids(BaseBuffer<char*>& target) {
  primary->load_ids(target);
  if (secondary) secondary->load_ids(target);
}

bool IDIndexSequence::update_id(char* id, int target) {
  bool result = primary->update_id(id, target);
  return (!result && secondary) ? secondary->update_id(id, target) : result;
}

void IDIndexSequence::update_order(vector<int>& iid_remap) {
  primary->update_order(iid_remap);
  if (secondary) secondary->update_order(iid_remap);  
}

int IDIndexSequence::get_index(char* id) {
  int index = primary->get_index(id);
  return (secondary && index < 0) ? secondary->get_index(id) : index; 
}

void IDIndexSequence::set_secondary(IDIndex* sec) {
  if (secondary) delete secondary; secondary = sec;  
  if (primary->has_index) build_index();
  if (primary->has_map) build_map();
}

void IDIndexCore::init(long amount, int exp_length) {
  IDIndex::init(amount, exp_length);
  input_ids.init(amount, exp_length);  
  id_map.expansion_reserve(amount);
}

void IDIndexCore::clear(bool clear_ids) {
  IDIndex::clear(clear_ids);
  if (clear_ids) input_ids.clear(); 
  id_index.clear(); id_map.clear(); id_tree.clear();
  index_max = 0;
}
void IDIndexCore::clear_tree() {IDIndex::clear_tree(); id_tree.clear();}

int IDIndexCore::build_index(vector<char*>* store_duplicates) {if (has_index) return id_index.size();
  input_ids.unload_index(id_index); 
  int used = IDSorter().run(id_index, id_map, index_max, drop_duplicates, store_duplicates);

  has_index = true;
  return used;
}

void IDIndexCore::build_map() {if (has_map) return;
  if (!has_index) build_index();
  if (tree_depth > 0) id_tree.set_max(tree_depth);
  id_tree.build(id_index);

  if (!id_tree.get_status()) {_LOG.error("initialising DataFrame object") << "construction of ID tree failed" << endl; die();}
  has_map = true;  
}

bool IDIndexCore::update_id(char* id, int target) {
  if (!has_map) build_map();
  int index = id_tree[id];
  if (index >= 0) {id_map[index] = target; return true;}
  else return false;
}

void IDIndexCore::update_order(vector<int>& iid_remap) {
  for (int i = 0; i < index_max; i++) id_map[i] = iid_remap[id_map[i]];
}

int IDIndexCore::get_index(char* id) {
  int index = has_map ? id_tree[id] : -1;
  if (index >= 0) index = id_map[index];
  return index;
}

map<int,set<char*> > IDIndexCore::reverse_map(vector<char>& used) {
  map<int,set<char*> > out; map<int,set<char*> >::iterator found; 
  int idx_total = id_index.size(), max = used.size();
  for (int i = 0; i < idx_total; i++) {int iid = id_map[i];
    if (iid >= 0 && iid < max && used[iid]) {
      if ( (found = out.find(iid)) == out.end() ) {
        set<char*> add; add.insert(id_index[i]);
        out[iid] = add;      
      } else found->second.insert(id_index[i]);
    }
  }
  return out;
}

void IDIndexCore::load_ids(BaseBuffer<char*>& target) {int idx_total = id_index.size();
  for (int i = 0; i < idx_total; i++) {
    if (id_map[i] >= 0) {
      char*& curr_id = target[id_map[i]];
      if (!curr_id) curr_id = id_index[i];
    }
  }
}

int IDIndexExternal::process(bool is_sorted, bool is_filtered) {
  if (is_filtered) {
    id_map.assign(external_map, external_elem);
    for (int i = 0; i < external_elem; i++) input_ids.add(external_index[i], true);
  } else {
    id_map.resize(expected_elem); int write = 0, used_elem = expected_elem, incr_size = max((external_elem - used_elem) / 10, 100);
    for (int i = 0; i < external_elem; i++) {
      if (external_map[i] >= 0) {
        input_ids.add(external_index[i], true);        
        id_map[write++] = external_map[i];
        if (write >= used_elem) {used_elem += incr_size; id_map.resize(used_elem);}
      }
    }
    id_map.resize(write);
  }
 
  input_ids.unload_index(id_index); index_max = id_index.size();
  if (!is_sorted) IDSorter().run(id_index, id_map, index_max, drop_duplicates);
  has_index = true;
  return id_index.size();
}

void IDIndexPartition::init(long amount, int exp_length) {
  IDIndex::init(amount, exp_length);
  no_blocks = amount >= 5*block_size ? amount / block_size : 1;
  block_size = no_blocks > 1 ? 1.05*block_size : amount;

  for (int i = 0; i < no_blocks; i++) {
    blocks.push_back(new IDIndexCore());
    for (map<string,long>::iterator it = stored_params.begin(); it != stored_params.end(); ++it) blocks.back()->set_param(it->first, it->second);
    blocks.back()->init(block_size, exp_length);  
  }
}

void IDIndexPartition::clear(bool clear_ids) {
  IDIndex::clear(clear_ids);
  stored_params.clear();
  for (int i = 0; i < no_blocks; i++) blocks[i]->clear(clear_ids);
}
void IDIndexPartition::clear_tree() {
  IDIndex::clear_tree();
  for (int i = 0; i < no_blocks; i++) blocks[i]->clear_tree();
}

void IDIndexPartition::set_param(const string& name, long value) {
  IDIndex::set_param(name, value);
  for (int i = 0; i < no_blocks; i++) blocks[i]->set_param(name, value);
  stored_params[name] = value;
}

int IDIndexPartition::build_index(vector<char*>* store_duplicates) {
  if (!has_index) {
    index_used = 0; index_ids = 0; has_index = true; 
    for (int i = 0; i < no_blocks; i++) {IDIndexCore* curr = blocks[i];
      index_used += curr->build_index(store_duplicates);
      index_ids += curr->id_size();
    }    
    index_ids = round(index_ids / float(no_blocks));
  }
  return index_used;
}

void IDIndexPartition::build_map() {
  if (!has_map) {
    has_map = true;
    for (int i = 0; i < no_blocks; i++) blocks[i]->build_map();
  }
}

map<int,set<char*> > IDIndexPartition::reverse_map(vector<char>& used) {
  map<int,set<char*> > out = blocks[0]->reverse_map(used);
  for (int i = 1; i < no_blocks; i++) {
    map<int,set<char*> > add = blocks[i]->reverse_map(used);
    merge_maps(out, add);  
  }
  return out;
}

void IDIndexPartition::load_ids(BaseBuffer<char*>& target) {for (int i = 0; i < no_blocks; i++) blocks[i]->load_ids(target);}
void IDIndexPartition::update_order(vector<int>& iid_remap) {for (int i = 0; i < no_blocks; i++) blocks[i]->update_order(iid_remap);}

int IDIndexPartition::hash(const string& id) {
  if (no_blocks > 1) {
    unsigned long long hash = 5381; int len = min(int(id.length()), HASH_MAX);
    for (int i = 0; i < len; i++) hash = ((hash << 5) + hash) + (unsigned char) id[i];  
    return hash % no_blocks;  
  } else return 0;
}

int IDIndexPartition::hash(char* id) {
  if (no_blocks > 1) {
    unsigned long long hash = 5381; 
    for (char *curr = id, *end = id + HASH_MAX; *curr && curr < end; curr++) hash = ((hash << 5) + hash) + (unsigned char) *curr;  
    return hash % no_blocks;  
  } else return 0;
}


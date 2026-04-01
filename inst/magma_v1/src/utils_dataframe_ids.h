/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef UTILS_DATAFRAME_IDS_H
#define UTILS_DATAFRAME_IDS_H

#include <map>
#include "buffer.h"
#include "utils.h"


class IDSorter {
  struct SortObj {
    bool operator() (const pair<char*,int> &i, const pair<char*,int> &j) {return strcmp(i.first, j.first) < 0;}
  };

public:
  static bool idcmp(char* s1, char* s2);  

  int run(BaseBuffer<char*>& id_buff, BaseBuffer<int>& index_buff, int size, bool drop_duplicates, vector<char*>* store_duplicates=0);
};

class IDNode; 
class IDTree {
  vector<int> leaf_index;
  vector<unsigned char> leaf_depth;
  char** ids;

  IDNode* root;
  long id_count;
  bool tree_status;
  
  int leaf_max;

  void leaf_insert(int location) {leaf_index.push_back(location); leaf_depth.push_back(0);}
  int leaf_insert(int location, char value, unsigned char depth);
  int leaf_insert(vector<int>& location, vector<char>& value, unsigned char depth, int offset, int total);
  
  IDNode* make_node(vector<int>& location, vector<char>& value, unsigned char level);
  IDNode* make_node(vector<int>& location, vector<char>& value, unsigned char level, int offset, int total);   
  IDNode* partition(int from, int to, unsigned char level);
  
public:
  IDTree() : root(0), tree_status(false), leaf_max(100) {}
  ~IDTree() {clear();}
  
  void set_max(int value) {leaf_max = value;}  
  void build(BaseBuffer<char*>& id_buff);
  void clear();

  bool get_status() {return tree_status;}
  int operator[](char* id);
};
 
class IDIndex {
protected:
  int expected_elem;
  int default_id_size;
  
  int tree_depth;
  bool drop_duplicates;
  
  void merge_maps(map<int,set<char*> >& base, map<int,set<char*> >& add);

public:
  IDIndex() : expected_elem(1000), default_id_size(10), tree_depth(10), drop_duplicates(false), has_index(false), has_map(false) {}  
  virtual ~IDIndex() {}

  bool has_index;
  bool has_map;

  virtual void init(long amount, int exp_length=0);
  virtual void clear(bool clear_ids=true);
  virtual void clear_tree();  
  virtual void set_param(const string& name, long value);  

  virtual int build_index(vector<char*>* store_duplicates=0) = 0; 
  virtual void build_map() = 0;

  virtual map<int,set<char*> > reverse_map(vector<char>& used) = 0;
  virtual void load_ids(BaseBuffer<char*>& target) = 0;

  virtual void add_id(const string& id) = 0;
  virtual void add_id(const string& id, int iid) = 0;  
  virtual bool update_id(char* id, int target) = 0;  
  virtual void update_order(vector<int>& iid_remap) = 0;
  virtual int get_index(char* id) = 0;
  virtual int size(bool max) = 0;
  virtual int id_size() = 0;
}; 

class IDIndexSequence : public IDIndex {
  IDIndex* primary;
  IDIndex* secondary;
  
public:
  IDIndexSequence(IDIndex* pri, IDIndex* sec=0) : primary(pri), secondary(0) {set_secondary(sec);}
  ~IDIndexSequence() {delete primary; delete secondary;}

  void init(long amount, int exp_length=0);
  void clear(bool clear_ids=true);
  void clear_tree();  
  void set_param(const string& name, long value);  

  int build_index(vector<char*>* store_duplicates=0);
  void build_map();
  
  map<int,set<char*> > reverse_map(vector<char>& used); 
  void load_ids(BaseBuffer<char*>& target); 
  
  IDIndex* get_primary() {return primary;}
  IDIndex* get_secondary() {return secondary;}  
  void set_secondary(IDIndex* sec);

  void add_id(const string& id) {primary->add_id(id);}
  void add_id(const string& id, int iid) {primary->add_id(id, iid);}  
  bool update_id(char* id, int target);
  void update_order(vector<int>& iid_remap);
  int get_index(char* id);
  int size(bool max) {return primary->size(max);}
  int id_size() {return primary->id_size();}
};

class IDIndexCore : public IDIndex {
protected:
  IndexedArchiveBuffer input_ids;
  BaseBuffer<char*> id_index;   ///mapping ID to ID
  ExpandingBuffer<int> id_map;  ///mapping ID to internal ID
  IDTree id_tree;               ///ID to mapping ID

  int index_max;

public:
  IDIndexCore() : index_max(0) {}
  virtual ~IDIndexCore() {}

  virtual void init(long amount, int exp_length=0);
  void clear(bool clear_ids=true);
  void clear_tree();  

  int build_index(vector<char*>* store_duplicates=0); 
  void build_map();
  
  map<int,set<char*> > reverse_map(vector<char>& used); 
  void load_ids(BaseBuffer<char*>& target);    

  void add_id(const string& id) {input_ids.add(id); id_map[index_max] = index_max; index_max++;}
  void add_id(const string& id, int iid) {input_ids.add(id); id_map[index_max++] = iid;}  
  bool update_id(char* id, int target);
  void update_order(vector<int>& iid_remap);
  int get_index(char* id);
  int size(bool max) {return max ? index_max : id_index.size();}
  int id_size() {return input_ids.size(false);}
  
  BaseBuffer<int>& get_iid_map() {return id_map;} 
  BaseBuffer<char*>& get_name_map() {return id_index;}
  BaseBuffer<char*>& get_names() {return input_ids.get_index();} 
}; 

class IDIndexExternal : public IDIndexCore {
  char** external_index;
  int* external_map;
  int external_elem;
  
public:
  IDIndexExternal(char** ext_index, int* ext_map, int ext_elem) : external_index(ext_index), external_map(ext_map), external_elem(ext_elem) {}
  
  int process(bool is_sorted, bool is_filtered);
};
 
class IDIndexPartition : public IDIndex {
  map<string,long> stored_params;

  vector<IDIndexCore*> blocks;
  int block_size;
  int no_blocks;
  
  int index_max;
  int index_used;
  int index_ids;

  int hash(const string& id);
  int hash(char* id); 
  
public:
  IDIndexPartition(int block_size) : block_size(block_size), no_blocks(0), index_max(0), index_used(0), index_ids(0) {}
  ~IDIndexPartition() {for (int i = 0; i < no_blocks; i++) delete blocks[i];}

  void init(long amount, int exp_length=0);
  void clear(bool clear_ids=true);
  void clear_tree();  
  void set_param(const string& name, long value);  

  int build_index(vector<char*>* store_duplicates=0); 
  void build_map();

  map<int,set<char*> > reverse_map(vector<char>& used);
  void load_ids(BaseBuffer<char*>& target); 

  void add_id(const string& id) {blocks[hash(id)]->add_id(id, index_max++);}
  void add_id(const string& id, int iid) {blocks[hash(id)]->add_id(id, iid); index_max++;}
  bool update_id(char* id, int target) {return blocks[hash(id)]->update_id(id, target);}
  void update_order(vector<int>& iid_remap);  
  int get_index(char* id) {return blocks[hash(id)]->get_index(id);}
  int size(bool max) {return max ? index_max : index_used;}
  int id_size() {return index_ids;}
};
 
 
class IDNode {
public:
  virtual ~IDNode() {}
  virtual int run(char* id) = 0;
};

class RootNode : public IDNode {
  IDNode* child;
  char* prefix;
  int length;
  
public:
  RootNode(string prefix_str, IDNode* child) : child(child) {
    length = prefix_str.size(); prefix = new char[length];
    for (int i = 0; i < length; i++) prefix[i] = prefix_str[i]; 
  }
  ~RootNode() {delete[] prefix; delete child;}

  int run(char* id) {return strncmp(id, prefix, length) == 0 ? child->run(id+length) : 0;}      
};
 
class TrunkNode : public IDNode {
  IDNode* child;
  char value;

public: 
  TrunkNode(char value, IDNode* child) : child(child), value(value) {}
  ~TrunkNode() {delete child;}
  
  int run(char* id) {return *id == value ? child->run(id+1) : 0;}
};
  
class BranchNode : public IDNode {
  char lower;
  char upper;
  IDNode** children;

public:
  BranchNode(char lower, char upper) : lower(lower), upper(upper) {int range = int(upper)-lower+1; children = new IDNode*[range]; Utils::set_zero(children, range);} 
  ~BranchNode() {
    if (children) {
      for (int i = 0; i <= int(upper)-lower; i++) delete children[i];
      delete[] children;
    } 
  }
  
  int run(char* id);
  bool add_child(char value, IDNode* child);    
};
   
class LeafNode : public IDNode {
protected:
  char lower;
  char upper;
  int base_index;
  
public:
  LeafNode(char lower, char upper, int offset) : lower(lower), upper(upper), base_index(offset - lower) {}
  ~LeafNode() {}
  
  int run(char* id) {return (*id < lower || *id > upper) ? 0 : base_index + *id;}
};

class StubLeafNode : public IDNode {
  char value;
  int index;
  
public:
  StubLeafNode(char value, int offset) : value(value), index(offset) {}
  
  int run(char* id) {return *id == value ? index : 0;}
};

class DeadNode : public IDNode {
  public:
  int run(char* id) {UNUSED(id); return 0;}
};
  
class SplitNode : public IDNode {
  char pivot;
  IDNode* lower;
  IDNode* upper;

public:
  SplitNode(char value, IDNode* lower, IDNode* upper) : pivot(value), lower(lower), upper(upper) {}
  ~SplitNode() {delete lower; delete upper;}
  
  int run(char* id) {return *id < pivot ? lower->run(id) : upper->run(id);}
};

class ValenceNode : public IDNode {
protected:
  IDNode* child;
  int null_index;

public:
  ValenceNode(IDNode* child, int null_index) : child(child), null_index(null_index) {}
  virtual ~ValenceNode() {delete child;}
  
  virtual int run(char* id) = 0;
};

class PosNode : public ValenceNode {
public:
  PosNode(IDNode* child, int null_index) : ValenceNode(child, null_index) {}

  int run(char* id) {return *id > 0 ? child->run(id) : (*id == 0 ? null_index : 0);}  
};

class NegNode : public ValenceNode {
public:
  NegNode(IDNode* child, int null_index) : ValenceNode(child, null_index) {}

  int run(char* id) {return *id < 0 ? child->run(id) : (*id == 0 ? null_index : 0);}  
};

class PosNegNode : public IDNode {
  IDNode* negative;
  IDNode* positive;
  int null_index;  
  
public:
  PosNegNode(IDNode* neg, IDNode* pos, int null_index=0) : negative(neg), positive(pos), null_index(null_index) {}
  ~PosNegNode() {delete negative; delete positive;}
  
  int run(char* id) {return *id < 0 ? negative->run(id) : (*id > 0 ? positive->run(id) : null_index);}  
};


#endif /** UTILS_DATAFRAME_IDS_H */

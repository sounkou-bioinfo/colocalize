/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef DATA_SNPDATA_H
#define DATA_SNPDATA_H

#include "utils_dataframe.h"

class SNPData : public IndexedDataFrame {
  bool has_synonyms;
  
  string make_dummy(char* str, char* str2);
  
public:
  PvalueVariable<double> pval;
  PositiveNumber<int> nsamp;
  PositiveWeightVariable<float> weights;    

  SNPData() : has_synonyms(false) {
    delete id_index; id_index = new IDIndexPartition(100000);
    id_index->set_param("drop_duplicates", true);
    id_index->set_param("tree_depth", 25);
    register_variable(pval, "PVAL"); register_variable(nsamp, "NSAMP"); register_variable(weights, "WEIGHT");
  }
 
  void load_synonyms(IDIndex* synonyms);
  map<int,set<char*> > get_synonyms(set<string>& ids);
  
  using IndexedDataFrame::add_iid;
  using IndexedDataFrame::drop_id;
  using IndexedDataFrame::remap_id;
  using IndexedDataFrame::filter_ids;
};

class SNPSynonyms : public IndexedDataFrame {
private:
  struct SynonymRange {
    int from;
    int to;
    SynonymRange* child;

    SynonymRange(int from, int to) : from(from), to(to), child(0) {}    
    SynonymRange(int from, int to, SynonymRange* child) : from(from), to(to), child(child) {}
    ~SynonymRange() {delete child;}
    
    SynonymRange* merge(SynonymRange* other);
    void delete_elem(int pos);
  };
  
  struct SynonymRangeIterator {
    SynonymRange* source;
    int value;
    
    SynonymRangeIterator(SynonymRange* source) : source(source) {value = source ? source->from-1 : 0;}
    
    bool empty();
    bool next();
  };

  vector<SynonymRange*> synonyms;
  
  int offset;
  int counts[3];

  void set_mode(const string& mode);
  int resolve_id(vector<int>& index, int id);

  bool valid_id(int iid) {UNUSED(iid); return false;}
  int validate_id(int iid) {UNUSED(iid); return -1;}

public:
  NumericID<int> snpid;
  NumericID<int> synid;

  enum DupMode {DupError, DupDropAll, DupDropDup, DupSkipLine, DupSkipDup} duplicate_mode;
  SNPSynonyms(const string& mode) : offset(0) {
    set_mode(mode); counts[0] = counts[1] = counts[2] = 0;
    id_index->set_param("tree_depth", 25);
    register_variable(snpid, "SNPID", true); register_variable(synid, "SYNID", true);
  }  

  bool init(int exp_amount, int exp_length=0);
  void clear_data(bool clear_ids=true);  
  void build_index();  

  vector<vector<char*> > process(SNPData& snp_info);

  using IndexedDataFrame::get_iid;
  void process_synonym();

  int get_count(const string& type);
};


#endif /** DATA_SNPDATA_H */


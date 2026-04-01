/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEUTILS_H
#define GENEUTILS_H

#include <algorithm>
#include <vector>
#include <ostream>
#include <vector>

#include "utils.h"

using namespace std;

struct GeneLocation {
  int chromosome;
  unsigned long begin;
  unsigned long end;

  GeneLocation() : chromosome(-1), begin(0), end(0) {}
  GeneLocation(int) : chromosome(-1), begin(0), end(0) {} /// conversion constructor
  GeneLocation(int chr, unsigned long from, unsigned long to) : chromosome(chr), begin(from), end(to) {}

  GeneLocation operator+(const GeneLocation& other) const {
    return GeneLocation(chromosome == other.chromosome ? chromosome : -1, min(begin, other.begin), max(end, other.end));
  }

  bool operator==(const GeneLocation& other) const {
    return chromosome == other.chromosome && begin == other.begin && end == other.end;
  }
  
  bool operator!=(const GeneLocation& other) const {return !(*this == other);}


  ///NOTE: < and > operators are NOT symmetrical (and thus, same for <= and >=)
  ///If gene A is fully contained in gene B, both B < A and B > A return true
  bool operator<(const GeneLocation& other) const {
    return (chromosome < other.chromosome) ||
              (chromosome == other.chromosome && begin < other.begin) ||
                 (chromosome == other.chromosome && begin == other.begin && end < other.end);
  }
  bool operator<=(const GeneLocation& other) const {return *this == other || *this < other;}
  
  bool operator>(const GeneLocation& other) const {
    return (chromosome > other.chromosome) ||
              (chromosome == other.chromosome && end > other.end) ||
                 (chromosome == other.chromosome && end == other.end && begin > other.begin);
  }
  bool operator>=(const GeneLocation& other) const {return *this == other || *this > other;}

  ///Assumes other gene is on same chromosome
  bool before(const GeneLocation& other) const {return end < other.begin;}
  bool after(const GeneLocation& other) const {return begin > other.end;}

  ///Assumes point is on same chromosome
  ///if  0: point inside range (including bounds)
  ///   -1: point before range
  ///    1: point after range
  int compare(unsigned long point) {
    return (begin <= point) - (end >= point);
  }

  ///if  0: genes overlap
  ///   -1: other fully before this gene
  ///    1: other fully after this gene
  int compare(const GeneLocation& other) {
    if (chromosome == other.chromosome) {
      return end < other.begin ? 1 : (begin > other.end ? -1 : 0);
    }
    return chromosome < other.chromosome ? 1 : -1;
  }

  bool in_range(const GeneLocation& other, unsigned long range) const {
    if (chromosome != other.chromosome) return false;
    if (end < other.begin) return (other.begin - end) <= range;
    if (begin > other.end) return (begin - other.end) <= range;
    return true; ///Locations overlap
  }
  
  friend ostream& operator<<(ostream& os, const GeneLocation& loc) {
    os << Utils::chr_string(loc.chromosome, false) << ':' << loc.begin << ':' << loc.end;
    return os;
  }
};

class GeneSorter {
  const vector<GeneLocation>& gene_loc;
  bool end_sort;

  struct SortObj { ///Guarantees that if A before B, A.begin <= B.begin
    GeneSorter &sorter;
    SortObj(GeneSorter &es) : sorter(es) {}
    bool operator () (const long &i, const long &j) {return sorter.gene_loc[i] < sorter.gene_loc[j];}
  };

public:
  GeneSorter(const vector<GeneLocation>& gl) : gene_loc(gl) {}
  vector<long> run() {
    vector<long> index; index.reserve(gene_loc.size());
    for (int i = 0; i < gene_loc.size(); i++) index.push_back(i);
    run(index); return index;
  }
  void run(vector<long>& geneids) {sort(geneids.begin(), geneids.end(), SortObj(*this));}
};

struct ModelState {///Values in output sort order
  enum Model {MultiModel, SNPwiseModelMeanBrown, SNPwiseModelMeanImhof, SNPwiseModelTop, SNPwiseModelTopOne, PrinCompReg, VariPow, TDTModel, BlockModel, UnknownModel, _MODELCOUNT_};
  enum Analysis {FullAnalysis, InteractionAnalysis, MainAnalysis, UnknownAnalysis, _ANALYSISCOUNT_};
  
  static string model_name(Model model);
  static string analysis_name(Analysis analysis);
  
  static Model model_type(const string& name);
  static Analysis analysis_type(const string& name);  
  
  static bool is_model(const string& name) {return model_type(name) != UnknownModel;}
  static bool is_analysis(const string& name) {return analysis_type(name) != UnknownAnalysis;}
  
  static bool is_printable(Model model) {return model_name(model) != "#UNKNOWN#";}
  static bool is_printable(Analysis analysis) {return analysis_name(analysis) != "#UNKNOWN#";}  
};

struct ModelOutputInfo {
  ModelOutputInfo(ModelState::Model model, ModelState::Analysis analysis=ModelState::MainAnalysis, int index=0, int partition=0) : analysis(analysis), model(model), partition(partition), index(index) {} 
  ModelOutputInfo* copy() const {ModelOutputInfo* out = new ModelOutputInfo(model, analysis); out->partition = partition; out->index = index; return out;}
  ModelOutputInfo* copy(int new_index) const {ModelOutputInfo* out = this->copy(); out->index = new_index; return out;}
  ModelOutputInfo* copy(int new_index, int new_partition) const {ModelOutputInfo* out = this->copy(new_index); out->partition = new_partition; return out;}

  bool operator==(const ModelOutputInfo& other) {return (analysis == other.analysis && model == other.model && partition == other.partition && index == other.index);}
  bool operator!=(const ModelOutputInfo& other) {return !(*this == other);}

  ModelState::Analysis analysis;
  ModelState::Model model;
  int partition;
  int index;
  
  int hash() const {return analysis + ModelState::_ANALYSISCOUNT_*model + ModelState::_ANALYSISCOUNT_*ModelState::_MODELCOUNT_*partition;} 
};

struct ModelInfoBlock {
  static ModelInfoBlock basic_block() {ModelInfoBlock out; out.info.push_back(new ModelOutputInfo(ModelState::UnknownModel)); return out;}
  
  enum State {HasInteraction, HasPartitions, HasMultiModel, IsUnbalanced, ModelCount, PrintableModelCount, AnalysisCount, PartitionCount, MaxIndex, _STATECOUNT_};
  vector<ModelOutputInfo*> info;
  vector<string> partition_names;
  vector<int> state;
  int no_repeat;
  
  ModelInfoBlock() : no_repeat(1) {partition_names.push_back("#UNDEFINED#"); detect_state();}
  ~ModelInfoBlock() {clear();}

  ModelInfoBlock(const ModelInfoBlock& other) {operator=(other);}  
  ModelInfoBlock& operator=(const ModelInfoBlock& other);
  ModelInfoBlock copy() {return *this;}

  bool operator==(const ModelInfoBlock& other);
  bool operator!=(const ModelInfoBlock& other) {return !(*this == other);}   

  set<ModelState::Model> get_models(); 
  set<ModelState::Analysis> get_analyses();
  
  vector<int> main_index();
  vector<pair<int,string> > main_index(const string& prefix);
  vector<int> partition_index(bool skip_root);
  
  void clear() {for (int i = 0; i < info.size(); i++) delete info[i]; info.clear();}
  void detect_state();

  int size() {return info.size() * no_repeat;}
  bool is_consistent();  
  
  string get_partition(int index) {return (index < partition_names.size()) ? partition_names[index] : "#UNDEFINED#";}
  
  ModelInfoBlock& align(const ModelInfoBlock& other);  
  ModelInfoBlock& filter(ModelState::Analysis analysis, bool keep=true);
  ModelInfoBlock& filter(ModelState::Model model, bool keep=true);  
  ModelInfoBlock& truncate(int max_index);
  ModelInfoBlock& update_partitions(const vector<string>& partitions);  
  ModelInfoBlock& strip_partitions();  
  ModelInfoBlock& set_model(ModelState::Model model);
  ModelInfoBlock& set_partition(int partition);  
  ModelInfoBlock& prune_models();
  ModelInfoBlock& shift_index(int add);
  ModelInfoBlock& append(ModelInfoBlock& other);
  ModelInfoBlock& sort(bool reindex=false);
  ModelInfoBlock& sort(vector<int> partition_rank, bool reindex=false);
  
private:
  string make_name(const string& prefix, int index, bool print_full=false);

  struct InfoSort { 
    InfoSort(vector<int>& partition_rank) : partition_rank(partition_rank) {}
    vector<int>& partition_rank;
    bool operator () (ModelOutputInfo* i, ModelOutputInfo* j);
  };
};


#endif /**GENEUTILS_H*/
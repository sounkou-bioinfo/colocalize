/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SETMODEL_H
#define SETMODEL_H

#include <set>

#include "setstats.h"
#include "setoutput.h"

class SetModel {
protected:
  SetVariables& variable_data;

  string file_prefix;
  string file_suffix;
  
  double alpha;
  int abbr_length;
  bool abbr_file;

public:
  const string error_tag;

  SetModel(SetVariables& vars, const string& prefix, const string& suffix, const string& error_tag) : variable_data(vars), file_prefix(prefix), file_suffix(suffix), alpha(0.05), abbr_length(30), abbr_file(false), error_tag(error_tag) {}
  virtual ~SetModel() {}

  virtual void validate();
  virtual void print(bool add_lead) = 0;
  virtual void run() = 0;  
  
  void set_suffix(const string& suffix) {file_suffix = suffix;}
  void set_alpha(double value) {alpha = value;}
  void set_abbreviate(int length, bool to_file) {abbr_length = length; abbr_file = to_file;}
};

class SelfContainedModel : public SetModel {
  set<int> analyse;
public:
  SelfContainedModel(SetVariables& vars, const string& prefix, const string& error_tag) : SetModel(vars, prefix, "gsa.self", error_tag) {analyse = variable_data.variables.var_list(SetStatsUtils::vc_GeneSet);}

  void validate();
  void print(bool add_lead);
  void run();
};


class RegressionConfig {
public:
  class BaseInput;
  class InputIterator;

  enum AnalysisMode {am_Mono, am_Joint, am_JointPairs, am_Interaction, am_InteractionAll, am_InteractionEach, am_InteractionPairs}; 
  enum ConditionType {ct_Internal, ct_ConditionForce, ct_ConditionOnly, ct_ConditionShow, ct_ConditionInteraction, ct_ModeDependent, ct_COUNT}; 
  enum DirectionType {dt_GeneSet, dt_GeneCovar, dt_InteractionSS, dt_InteractionSC, dt_COUNT};

private:
  class MonoInput; class JointInput; class JointPairsInput; class InteractionInput; 
  class InteractionBlocksInput; class InteractionAllInput; class InteractionEachInput; class InteractionPairsInput; 

  enum VariableUse {vu_UseAll, vu_SetsOnly, vu_CovarOnly, vu_List} variable_filter;
  enum ParseFailure {pf_NotFound, pf_NotValid, pf_Ambiguous};
  AnalysisMode analysis_mode; 

  SetVariables& variable_data;
  string error_tag;

  set<int> condition[ct_COUNT];
  set<int> variables; 
  set<int> variables_auxiliary;
  vector<set<int> > variable_blocks;
  
  vector<short> direction;
  int size_ss, size_sc; double prop_ss;
  bool residualize, allow_collinear;
  
  BaseInput* model_input;
  
  void init();
  
  set<int> filter_set(set<int>& target, set<int>& filter);
  bool parse_variables(vector<string>& names, set<int>& target, const string& label, bool do_error=true);
  void process_specifications(vector<string>& list);

  void add_conditional(int var_id, ConditionType type) {condition[type].insert(var_id);}
  void add_conditional(set<int> var_ids, ConditionType type) {condition[type].insert(var_ids.begin(), var_ids.end());}

public:
  RegressionConfig(SetVariables& vars, const string& error_tag) : variable_filter(vu_UseAll), analysis_mode(am_Mono), variable_data(vars), error_tag(error_tag), residualize(false), allow_collinear(false), model_input(0) {init();}
  ~RegressionConfig();

  void set_direction(short dir, DirectionType type) {direction[type] = dir;}
  void set_interaction(int size_sc, int size_ss, double prop_ss) {this->size_sc = size_sc; this->size_ss = size_ss; this->prop_ss = prop_ss;}

  void set_internal(string mode, vector<string> filter=vector<string>());
  void set_conditional(vector<string> names, ConditionType type);

  void set_filter(string mode);
  void set_filter(vector<string> list);
  void set_mode(AnalysisMode mode);
  void set_models(vector<string> list, AnalysisMode mode);

  void set_residualize(bool value) {residualize = value;}
  void set_collinear(bool value) {allow_collinear = value;}
  
  void validate();
  void print();
  
  void make_input();
  BaseInput& get_input();
  bool get_residualize() {return residualize;}   
  bool get_collinear() {return allow_collinear;}   

  bool has_interaction() {return analysis_mode == am_Interaction || analysis_mode == am_InteractionAll || analysis_mode == am_InteractionPairs || analysis_mode == am_InteractionEach || !condition[ct_ConditionInteraction].empty();}

  static string parse_direction(int dir, bool verbose=true);
};


class RegressionModel : public SetModel {
  RegressionConfig config; 
  bool print_genes;

  void output_results(SetResultsData& results, bool print_genes);
  
  bool check_df(int df, int eff_size);
  void check_rsq(double rsq);
  void set_gene_info(RegressionComputation& computer);
  void create_terms(vector<SetStatsUtils::VariableInfo*> vars, vector<string>& target);

public:
  RegressionModel(SetVariables& vars, const string& prefix, const string& error_tag, bool print_genes=false) : SetModel(vars, prefix, "gsa", error_tag), config(vars, error_tag), print_genes(print_genes) {}

  RegressionConfig& get_config() {return config;}

  void validate() {config.validate();}
  void print(bool add_lead);
  void run();  
};


class RegressionConfig::InputIterator {
public:
  virtual ~InputIterator() {}

  virtual bool empty() = 0;
  virtual void next() = 0;
  virtual void reset() = 0;
  
  virtual int no_vars(bool main_only) = 0;
  virtual set<int>* get_variables() = 0;
  
  virtual void print() = 0;
};  


class RegressionConfig::BaseInput {
protected: 
  RegressionConfig& config;
  int no_cond; 
public: 
  BaseInput(RegressionConfig& config) : config(config) {no_cond = config.condition[ct_ConditionShow].size() + config.condition[ct_ConditionInteraction].size() + config.condition[ct_ModeDependent].size();}
  virtual ~BaseInput() {}
  
  virtual void next() = 0;
  virtual InputIterator* get_block() = 0;

  virtual int no_tests() = 0;
  virtual int no_models() = 0;
  virtual int no_analyses() {return no_models();} /// main analyses only, for multiple testing correction
  virtual pair<int,int> var_count(bool main_only) = 0; /// range(min, max) of main (+shown covar) variables

  short get_direction(SetStatsUtils::VariableInfo& var); 
  short get_direction(DirectionType type) {return config.direction[type];}
  set<int>& get_conditional(ConditionType type) {return config.condition[type];}
};

class RegressionConfig::MonoInput : public BaseInput {
  class MonoIterator : public InputIterator {
    set<int> &analyse, wrapper;
    set<int>::iterator current, last; 
    int no_covar;
  public:
    MonoIterator(set<int>& vars, int ncov) : analyse(vars), current(analyse.begin()), last(analyse.end()), no_covar(ncov) {}

    bool empty() {return current == last;}
    void next() {if (current != last) ++current;}
    void reset() {current = analyse.begin();}
    
    int no_vars(bool main_only) {return !main_only ? no_covar + 1 : 1;}
    set<int>* get_variables() {if (current != last) {wrapper.clear(); wrapper.insert(*current); return &wrapper;} else return 0;}
    
    void print();
  } data; 
  
  bool complete;   

public:
  MonoInput(RegressionConfig& config) : BaseInput(config), data(config.variables, no_cond), complete(false) {}

  void next() {complete = true;}
  InputIterator* get_block() {return !complete ? &data : 0;}

  int no_tests() {return config.variables.size() * (no_cond + 1);}
  int no_models() {return config.variables.size();}
  pair<int,int> var_count(bool main_only);
};

class RegressionConfig::JointInput : public BaseInput {
protected:
  class JointIterator : public InputIterator {
    vector<set<int> >& blocks;
    vector<int>& analyse;
    int nvar, ncov, current;
  public:
    JointIterator(vector<set<int> >& blocks, vector<int>& use, int nvar, int ncov) : blocks(blocks), analyse(use), nvar(nvar), ncov(ncov), current(0) {}

    bool empty() {return current >= analyse.size();}
    void next() {if (current < analyse.size()) current++;}
    void reset() {current = 0;}
    
    int no_vars(bool main_only) {return !main_only ? ncov + nvar : nvar;}
    set<int>* get_variables() {return current < analyse.size() ? &blocks[analyse[current]] : 0;}
    
    void print(); 
  }; 

  vector<set<int> >& variable_blocks;
  map<int,vector<int> > index;
  map<int,vector<int> >::iterator current; 
  JointIterator* iterator;
  
  pair<int,int> range;
  int ntest;

  void process_blocks();
  void set_iterator();

  JointInput(RegressionConfig& config, vector<set<int> >& blocks) : BaseInput(config), variable_blocks(blocks), iterator(0), ntest(0) {}
public:
  JointInput(RegressionConfig& config) : BaseInput(config), variable_blocks(config.variable_blocks), iterator(0), ntest(0) {process_blocks();}
  virtual ~JointInput() {delete iterator;}

  void next() {if (current != index.end()) {++current; set_iterator();}}
  InputIterator* get_block() {return iterator;}

  int no_tests() {return ntest + variable_blocks.size() * no_cond;}
  int no_models() {return variable_blocks.size();}
  int no_analyses() {return ntest;}
  pair<int,int> var_count(bool main_only);
};

class RegressionConfig::JointPairsInput : public JointInput {
  vector<set<int> > pair_blocks;  
  void make_pairs();
public:
  JointPairsInput(RegressionConfig& config) : JointInput(config, pair_blocks) {make_pairs();}  
};

class RegressionConfig::InteractionInput : public BaseInput {
protected:
  struct PairSS {
    int id; SetStatsUtils::SetInfo *set1, *set2;
    PairSS(SetStatsUtils::SetInfo* set1, SetStatsUtils::SetInfo* set2, int id=-1) : id(id), set1(set1), set2(set2)  {}
  };
  struct PairSC {
    int id; SetStatsUtils::SetInfo* set; SetStatsUtils::CovarInfo* cov;
    PairSC(SetStatsUtils::SetInfo* set, SetStatsUtils::CovarInfo* cov, int id=-1) : id(id), set(set), cov(cov) {}
  };
 
 class InteractionIterator : public InputIterator {
    vector<set<int> > blocks;
    string log_message;
    int nvar, ncov, current;
  public:
    InteractionIterator(int nvar, int ncov, const string& msg) : log_message(msg), nvar(nvar), ncov(ncov), current(0) {}

    void set_blocks(vector<set<int> >& input) {blocks.swap(input);}

    bool empty() {return current >= blocks.size();}
    void next() {if (current < blocks.size()) current++;}
    void reset() {current = 0;}
    
    int no_vars(bool main_only) {return !main_only ? ncov + nvar : nvar;}
    set<int>* get_variables() {return current < blocks.size() ? &blocks[current] : 0;}
    
    void print() {if (!log_message.empty()) _LOG << "\t" << log_message << endl;} 
  }; 

  GeneSetData* set_info; 
  vector<GeneCovarData<float>*> covar_info; 
  
  VariableIndex* current;
  InteractionIterator* iterator;
  bool logged[2]; int counts[2];
  bool external_use, verbose, linked, delete_info;

  int ss_part; double ss_prop;
  int ss_min, ss_max, sc_min, sc_max;

  int process_ss(vector<PairSS>& blocks);
  int process_sc(vector<PairSC>& blocks);

  GeneCovarData<float>* add_covar_data(NumericID<int>*& main1, NumericID<int>*& main2);
  void create_covar_data(GeneCovarData<float>* info);

  void process_blocks();
  virtual void sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks) = 0;
  
  void make_iterator(VariableIndex* input, SetStatsUtils::VariableType type);
  virtual void set_iterator();

  InteractionInput(RegressionConfig& config, bool external=false, bool linked=false);
public:
  virtual ~InteractionInput() {delete iterator; if (delete_info) {delete set_info; Utils::delete_vector(covar_info);}}

  virtual void next() {if (current != 0) set_iterator();}
  virtual InputIterator* get_block() {return iterator;}

  virtual int no_tests() {return (counts[0] + counts[1]) * (no_cond + 3);}
  virtual int no_models() {return counts[0] + counts[1];}
  virtual pair<int,int> var_count(bool main_only);
};

class RegressionConfig::InteractionBlocksInput : public InteractionInput {
  vector<set<int> >& blocks;

  void sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks);
  void extract(vector<int>& target);
public:
  InteractionBlocksInput(RegressionConfig& config) : InteractionInput(config), blocks(config.variable_blocks) {process_blocks();}  
  InteractionBlocksInput(RegressionConfig& config, vector<set<int> >& input, vector<int>& interaction_ids) : InteractionInput(config, true), blocks(input) {
    process_blocks();
    set_iterator();
    extract(interaction_ids);
  }    
};

class RegressionConfig::InteractionAllInput : public InteractionInput {
  vector<vector<set<int> > > blocks;
  int interaction_count, analysis_count;
  int current_block;

  void sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks);
  void strip_storage(VariableIndex* info, vector<int>& used);
  void extract();
  
  void set_iterator();
public:
  InteractionAllInput(RegressionConfig& config);
  
  void next(); 
  InputIterator* get_block() {return iterator;}
  
  int no_tests() {return analysis_count * (no_cond + interaction_count + 1);}
  int no_models() {return analysis_count;}
  int no_analyses() {return analysis_count * interaction_count;}  
  pair<int,int> var_count(bool main_only);
};

class RegressionConfig::InteractionEachInput : public InteractionInput {
  void sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks);
public:
  InteractionEachInput(RegressionConfig& config) : InteractionInput(config) {process_blocks();}  
};

class RegressionConfig::InteractionPairsInput : public InteractionInput {
  void sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks);
public:
  InteractionPairsInput(RegressionConfig& config) : InteractionInput(config) {process_blocks();}  
};


#endif /** SETMODEL_H */

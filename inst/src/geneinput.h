/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEINPUT_H
#define GENEINPUT_H

#include <map>
#include <stack>
#include <deque>                                                    
#include <utility>
#include <memory>
#include <boost/unordered_map.hpp>

#include "baseinput.h"
#include "input.h"
#include "_global.h"


class GeneInputFile {
protected:
  struct InputColumn {
    virtual ~InputColumn(){}
    virtual void process_value(TextInput& fin, int gid) = 0;
  };

  struct InputUniColumn : public InputColumn {
    InputVariable* variable; string name; int index; bool allow_missing;
    InputUniColumn(InputVariable& var, string name, int index, bool allow_missing=false) : variable(&var), 
                    name(name), index(index), allow_missing(allow_missing) {variable->init();}
    void process_value(TextInput& fin, int gid);
  };

  template<typename VAR>
  struct InputMultiColumn : public InputColumn {
    BlockMultiVariable<VAR>* variable; bool allow_missing;
    vector<string> names; vector<int> indices;
    InputMultiColumn(BlockMultiVariable<VAR>& var, bool allow_missing=false) : variable(&var), allow_missing(allow_missing) {}

    void process_value(TextInput& fin, int gid) {
      for (int i = 0; i < indices.size(); i++) {
        variable->add(gid, i, fin[indices[i]]); 
        if ((allow_missing && !variable->is_valid_or_missing(true)) || (!allow_missing && !variable->is_valid(true))) fin.line_error(names[i] + " " + variable->get_msg(true), indices[i]+1);      
      }
    }

    void init() {
      if (variable->get_width() != names.size()) variable->set_width(names.size());
      if (!variable->initialised() && !names.empty()) variable->init();
    }
    void add(string name, int index) {names.push_back(name); indices.push_back(index);}
    int size() {return names.size();}
  };

  string filename;
  GeneData& gene_info;
  vector<InputColumn*> columns;  
  
  void add_column(InputVariable& var, const string& name, int index, bool allow_missing=false);
  
public:
  GeneInputFile(const string& filename) : filename(filename), gene_info(*(new GeneData())) {DataFrame::set_owner(&gene_info, this);}
  virtual ~GeneInputFile() {
    DataFrame::delete_owned(&gene_info, this);
    for (int i = 0; i < columns.size(); i++) delete columns[i];
  }

  GeneData& get_data() {return gene_info;}  
  GeneData* eject_data() {DataFrame::unset_owner(&gene_info, this); return &gene_info;}
};

class GeneInputRawFile : public GeneInputFile {
  void load(double rescale=1, set<string> covar_list=set<string>(), bool list_exclude=true);
public:
  GeneInputRawFile(const string& filename, double rescale_corrs=1) : GeneInputFile(filename) {load(rescale_corrs);}
  GeneInputRawFile(const string& filename, string covar_mode, vector<string> filter_list, double rescale_corrs=1) : GeneInputFile(filename) {
    if (covar_mode == "none") load(rescale_corrs, set<string>(), false);
    else if (covar_mode == "all") load(rescale_corrs, set<string>(), true);
    else {
      set<string> covar;
      for (int i = 0; i < filter_list.size(); i++) covar.insert(Utils::uppercase(filter_list[i]));
      load(rescale_corrs, covar, covar_mode == "exclude");
    }
  }
  
  void overwrite_stat(const string& filename);
};


class GeneInputOutFile : public GeneInputFile {
  enum PvaluePart {pp_Partition, pp_Analysis, pp_Model, pp_Unknown};
      
  map<string,int> column_index;
  bool internal_input;

  void load();
  
  bool has_column(const string& name) {return column_index.find(name) != column_index.end();}
  void map_column(InputVariable& var, const string& name, bool allow_missing=false);  

  void process_pval(vector<string>& pcol);
  bool check_consistency(vector<string>& parts, PvaluePart type);

  template<typename VAR>
  void multi_column(BlockMultiVariable<VAR>& var, const string& name1, const string& name2, bool allow_missing=false) {
    if (!has_column(name1) || !has_column(name2)) return;
    vector<string> names; names.push_back(name1); names.push_back(name2); 
    multi_column(var, names, allow_missing);
  }

  template<typename VAR>
  void multi_column(BlockMultiVariable<VAR>& var,  vector<string>& names, bool allow_missing=false) {
    InputMultiColumn<VAR>* col = new InputMultiColumn<VAR>(var, allow_missing);
    for (int i = 0; i < names.size(); i++) {if (has_column(names[i])) col->add(names[i], column_index[names[i]]);}
    if (col->size() > 0) {col->init(); columns.push_back(col);}
    else delete col;
  }

public:
  GeneInputOutFile(const string& filename, bool internal_input=false) : GeneInputFile(filename), internal_input(internal_input) {load();}
};

namespace GeneInput {
  enum Mode {Merge, Sum, Mean, CeilingMean, FloorMean, _WEIGHTEDMODES_, WeightedZstat, WeightedPval};
  static bool is_weighted(Mode mode) {return mode >= _WEIGHTEDMODES_;}

  template<typename D>
  class Aggregator {
  protected:
    Buffer<int>& access_map;
    vector<D*> data; D value;
    int size; int eff_size; 

  public:
    Aggregator(Buffer<int>& access_map) : access_map(access_map), size(access_map.ncol()), eff_size(0) {}
    virtual ~Aggregator() {}

    virtual bool set_data(vector<D*>& input);    
    virtual bool set_data(vector<BaseBuffer<D>*>& input);
    
    bool ready() {return eff_size > 0;}
    D& get_value() {return value;}
    
    virtual bool compute(int row) = 0;
  };
  
  template<typename D> 
  class MergeAggregator : public Aggregator<D> {public:
    MergeAggregator(Buffer<int>& access_map) : Aggregator<D>(access_map) {}   
    bool compute(int row);
  };
  
  template<typename D> 
  class SumAggregator : public Aggregator<D> {public:
    SumAggregator(Buffer<int>& access_map) : Aggregator<D>(access_map) {}   
    bool compute(int row);
  };

  template<typename D> 
  class MeanAggregator : public Aggregator<D> {
    short round_type;
    double integer_round(double value);    
  public:
    MeanAggregator(Buffer<int>& access_map, short round_type=0) : Aggregator<D>(access_map), round_type(round_type) {}   
    bool compute(int row);
  };
  
  template<typename D> 
  class WeightedAggregator : public Aggregator<D> {
  protected:
    double min_z; double max_z;
    vector<double*> weight_base; vector<double*> weights;
    Buffer<double>* corr_base; Buffer<double> corrs;
    
    D truncate(D value) {return value > min_z ? (value < max_z ? value : max_z) : min_z;}
  public:
    WeightedAggregator(Buffer<int>& access_map, pair<double,double> truncate) : Aggregator<D>(access_map), corr_base(0) {min_z = truncate.first; max_z = truncate.second;}   

    using Aggregator<D>::set_data;
    bool set_data(vector<D*>& input);            
    bool set_weights(vector<double*>& input);    
    bool set_corrs(Buffer<double>* corrs);
  };
  
  template<typename D> 
  class WeightedZstatAggregator : public WeightedAggregator<D> {public:
    WeightedZstatAggregator(Buffer<int>& access_map, pair<double,double> truncate) : WeightedAggregator<D>(access_map, truncate) {}   
    bool compute(int row);
  };
  
  template<typename D> 
  class WeightedPvalAggregator : public WeightedAggregator<D> {
    ConvertPvalToNorm from_pval; ConvertNormToPval to_pval;
  public:
    WeightedPvalAggregator(Buffer<int>& access_map, pair<double,double> truncate) : WeightedAggregator<D>(access_map, truncate) {}   
    bool compute(int row);
  };

}

class GeneInputCombine {
public: enum Mode {RawMode, OutMode} mode;
protected:
  Settings& settings;
  GeneData& gene_info;
    
  vector<GeneData*> input_data;
  vector<string> filenames;
  Buffer<int> access_map;
  
  string msg_ident;
  bool allow_overlap;

  void load_input(const vector<string>& filenames, bool keep_mid=false);
  
  void process_pcol(bool align_core, bool allow_partitions);
  bool check_map();
  
  virtual vector<double*> get_weights() {return vector<double*>();}
  virtual Buffer<double>* get_data_corrs() {return 0;}
  virtual pair<double,double> get_truncate() {return pair<double,double>(-500,500);}

  template<typename D> GeneInput::Aggregator<D>* make_aggregator(Buffer<int>& access_map, GeneInput::Mode mode);
  template<typename D> GeneInput::WeightedAggregator<D>* make_aggregator(Buffer<int>& access_map, vector<double*> weights, Buffer<double>* corrs, GeneInput::Mode mode);
  template<typename D, typename DATA> GeneInput::Aggregator<D>* make_aggregator(Buffer<int>& access_map, DATA input, GeneInput::Mode mode);

  template<typename VAR> VAR* get_variable(const string& name, int index);

  ///first, unnamed, parameter serves only for correct type/function determination  
  template<typename T> bool check_variable(UniVariable<T>&, const string& name);  
  template<typename VAR> bool check_variable(BlockMultiVariable<VAR>&, const string& name, bool check_width=true);
  template<typename VAR> bool check_variable(JaggedMultiVariable<VAR>&, const string& name);
  template<typename T> vector<BaseBuffer<T>*> get_data(UniVariable<T>&, const string& name);
  template<typename VAR> vector<Buffer<typename VAR::VTYPE>*> get_data(BlockMultiVariable<VAR>&, const string& name, bool check_width=true);
  template<typename VAR> vector<BaseBuffer<BaseBuffer<typename VAR::VTYPE>*>*> get_data(JaggedMultiVariable<VAR>&, const string& name);
  
  template<typename T> bool process_variable(UniVariable<T>& target, const string& name, GeneInput::Mode mode);
  template<typename VAR> bool process_variable(BlockMultiVariable<VAR>& target, const string& name, GeneInput::Mode mode);  
  template<typename VAR> bool process_variable(BlockMultiVariable<VAR>& target, vector<vector<typename VAR::VTYPE*> >& data, GeneInput::Mode mode);

  template<typename VAR> bool process_aligned_variable(BlockMultiVariable<VAR>& target, const string& name, GeneInput::Mode mode, bool align_main, bool skip_root=true);
  template<typename T> bool process_weighted_variable(UniVariable<T>& target, const string& name, vector<double>& weights, GeneInput::Mode mode);  

public:
  GeneInputCombine(Settings& settings, Mode mode, const string& identifier, bool allow_overlap=true) : mode(mode), settings(settings), gene_info(*(new GeneData())), msg_ident(identifier + " gene results files"), allow_overlap(allow_overlap) {
    DataFrame::set_owner(&gene_info, this);
  }
  virtual ~GeneInputCombine() {
    DataFrame::delete_owned(&gene_info, this);
    for (int i = 0; i < input_data.size(); i++) DataFrame::set_owner(input_data[i], this);
  }
};

class GeneInputMerge : public GeneInputCombine {
  vector<string> parse_filenames(string prefix, bool require);

  template<typename VAR> bool merge_variable(VAR& target, const string& name) {return process_variable(target, name, GeneInput::Merge);}
  template<typename VAR> bool merge_aligned_variable(BlockMultiVariable<VAR>& target, const string& name, bool align_main, bool skip_root=true) {return process_aligned_variable(target, name, GeneInput::Merge, align_main, skip_root);}  
  
  void merge_corrs();
public:
  GeneInputMerge(Settings& settings, Mode mode) : GeneInputCombine(settings, mode, "merging", false) {}
  
  bool run(const string& prefix, bool require=false);
};

class GeneInputMeta : public GeneInputCombine {
  class GeneCorrelationAccess {
    BaseBuffer<BaseBuffer<double>*>& corrs;
    int max_index; double dummy_value;  
  public:
    GeneCorrelationAccess(BaseBuffer<BaseBuffer<double>*>& corrs) : corrs(corrs) {max_index = corrs.size();}
    double& get(int gene_i, int gene_j); 
    int count(int gene_i) {return (gene_i >= 0 && gene_i < max_index && corrs[gene_i]) ? corrs[gene_i]->size() : 0;}     
  };

  vector<NonNegativeNumber<double>*> weights; 
  Buffer<double> data_corrs;

  void process_corrs();

  vector<double*> get_weights();
  Buffer<double>* get_data_corrs() {return &data_corrs;}
  pair<double,double> get_truncate() {return pair<double,double>(-settings.getn("metagene_truncate_low"), settings.getn("metagene_truncate_high"));}
  void set_weights();
  void clear_weights();

  void set_nsamp_part();
public:
  GeneInputMeta(Settings& settings, Mode mode) : GeneInputCombine(settings, mode, "meta-analysing", true) {}  
  ~GeneInputMeta() {clear_weights();}

  void run();
};

#include "geneinput.tpp"

#endif /** GENEINPUT_H */
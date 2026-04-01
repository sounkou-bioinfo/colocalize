/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SETOUTPUT_H
#define SETOUTPUT_H

#include "utils_dataframe.h"
#include "setstats.h"

class SetResultsData : public ExpandingDataFrame {
public:
  static const short CoreBlock;
  static const short RegressionBlock;

  NumericID<int> model_id;
  StringVariable model_term;
  StringVariable var_name;
  EnumVariable var_type;
  PositiveNumber<int> ngenes;
  NumericVariable<double> coefficient;
  NumericVariable<double> coefficient_std;  
  NumericVariable<double> std_error;
  PvalueVariable<double> pval;
  EnumVariable warnings;

  SetResultsData(int init_size=0) : var_type("SET COVAR INTER-SS INTER-SC"), warnings("COLLINEAR", "-") {
    register_variable(model_id, "MODEL_ID", false, RegressionBlock);
    register_variable(model_term, "MODEL_TERM", false, RegressionBlock);
    register_variable(var_name, "VARIABLE", true, CoreBlock);
    register_variable(var_type, "TYPE", true, RegressionBlock);
    register_variable(ngenes, "NGENES", true, CoreBlock);
    register_variable(coefficient, "BETA", true, CoreBlock);
    register_variable(coefficient_std, "BETA_STD", true, RegressionBlock);        
    register_variable(std_error, "SE", true, RegressionBlock); 
    register_variable(pval, "P", true, CoreBlock);    
    register_variable(warnings, "WARNINGS", true, CoreBlock);    
    if (init_size > 0) init(init_size);
  }
  
  void block_init(short block) {if (block > CoreBlock) block_init(CoreBlock); ExpandingDataFrame::block_init(block);}  
}; 

class SetOutput {
public:
  enum OutputMode {Regression, SelfContained};
private:
  SetResultsData& results;
  string prefix;

  vector<pair<string,string> > parameters;  
  vector<OutputColumn*> output; 

  int abbr_length;
  bool abbr_file;
  vector<StringVariable*> abbr_vars;

  void clear_output();
  void make_columns(OutputMode mode);

  void abbreviate_names(); 
  void write_abbreviations(string filename, bool indent=true);
  
public:
  SetOutput(SetResultsData& res, string prefix) : results(res), prefix(prefix), abbr_length(0), abbr_file(false) {}
  ~SetOutput() {clear_output();}

  void set_parameter(const string& name, string value) {parameters.push_back(pair<string,string>(name, value));}
  void set_parameter(const string& name, vector<string>& values, string sep=",") {set_parameter(name, Utils::join_string(values, sep+" "));}

  template<typename T>
  void set_parameter(const string& name, T value) {set_parameter(name, Utils::num_string(value));}
  
  void set_abbreviate(int length, bool to_file) {abbr_length = length >= 20 ? length : 0; abbr_file = to_file;}
  
  void write(OutputMode mode, string suffix, bool indent=true);
};

class SetGenesOutput {
  SetVariables& variable_data;
  GeneData& gene_info;
  
  string filename; 
  FormattedOutput output_stream;

  int no_tests;
  double alpha;
  
  int current_set;

  deque<OutputColumn*> output; 

  void clear_output();
  void write_core(const string& name, int set_size, int ngenes, double pval, int model_id, bool interaction);
  
public:
  SetGenesOutput(SetVariables& vars, string prefix, string suffix="gsa.sets.genes") : variable_data(vars), gene_info(vars.get_genes()), filename(win_txt(prefix + "." + suffix + ".out")), no_tests(0), alpha(0.05), current_set(0) {}
  ~SetGenesOutput() {clear_output();}

  void set_signif(int ntest, double thresh=0.05) {no_tests = ntest; alpha = thresh;}
  int get_count() {return current_set;}

  void write_set(SetStatsUtils::VariableInfo& var, double pval, int model_id=0);
  void print_log(bool indent=true);
};



#endif /** SETOUTPUT_H */

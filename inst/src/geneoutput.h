/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEOUTPUT_H
#define GENEOUTPUT_H

#include "data_genedata.h"
#include "engineutils.h"
#include "parse.h"

namespace OutputOptions {enum Option {HideAsymptotic, HideMultiSubPval, ShowFitted, _OPTIONCOUNT_};}

class GeneOutput {
  GeneData& gene_info;

  vector<pair<string,string> > parameters;  
  vector<double> options;
  vector<pair<int,string> > pval_info;
  vector<int> npart_info;
  vector<int> partition_rank;

  string prefix;
  int version_id;
  bool indent;
 
  void set_defaults();
  void process_pval();

public:
  GeneOutput(GeneData& gd, string prefix) : gene_info(gd), prefix(prefix) {set_defaults();}
  
  void set_option(OutputOptions::Option option, double value) {options[option] = value;}
  void set_partition_order(const vector<int>& rank) {partition_rank = rank;}  
  void set_version(int version) {version_id = version;}
  void set_prefix(const string& pref) {prefix = pref;}
  void set_indent(bool value) {indent = value;}
  
  void set_parameter(const string& name, string value) {parameters.push_back(pair<string,string>(name, value));}
  void set_parameter(const string& name, vector<string>& values, string sep=",") {set_parameter(name, Utils::join_string(values, sep+" "));}
  void clear_parameters() {parameters.clear();}

  void write_outfile(string label, string suffix);  
  void write_outfile(string suffix="genes") {write_outfile("", suffix);}
  void write_rawfile(string suffix="genes");
};

#endif /** GENEOUTPUT_H */

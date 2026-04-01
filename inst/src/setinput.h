/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SETINPUT_H
#define SETINPUT_H

#include "parse.h"
#include "data_genedata.h"
#include "setstats.h"

class SetInput {
  Settings& settings;
  GeneData& gene_info;

  void filter_genes(const string& filename, const string& mode);
  void process_missing(SetVariables* var_data);
  void prep_variables(SetVariables* var_data);  
  
  GeneSetData* load_setfile(const string& dataid, const string& filename, bool by_column);
  GeneCovarData<float>* load_covarfile(const string& dataid, const string& filename);  

public:
  SetInput(Settings& settings, GeneData& info) : settings(settings), gene_info(info) {DataFrame::set_owner(&gene_info, this);}
  ~SetInput() {DataFrame::delete_owned(&gene_info, this);}

  SetVariables* load_data();
};


#endif /** SETINPUT_H */


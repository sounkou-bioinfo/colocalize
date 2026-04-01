/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SETANALYSIS_H
#define SETANALYSIS_H

#include <stack>
#include <deque>
#include <utility>

#include "geneinput.h"
#include "setinput.h"
#include "utils.h"
#include "_global.h"
#include "setstats.h"
#include "setmodel.h"

class SetAnalysis {
  const Settings& settings;
  
  SetVariables* variable_data;
  vector<SetModel*> models;
  
  void parse_models();

public:
  SetAnalysis(const Settings& settings, SetVariables* vars) : settings(settings), variable_data(vars) {}
  ~SetAnalysis() {
    Utils::delete_vector(models);
    delete variable_data;
  }

  void run();

}; 

#endif /** SETANALYSIS_H */

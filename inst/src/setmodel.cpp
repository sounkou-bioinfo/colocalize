/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "setmodel.h"
#include "geneoutput.h"

#define MIN_DF 5
#define MAX_RSQ 0.95

using namespace SetStatsUtils;

void SetModel::validate() {
  _LOG.error(error_tag) << "model has no definition" << endl;
  die();
}

void SelfContainedModel::validate() {
  if (analyse.empty()) {_LOG.error(error_tag) << "model has no definition" << endl; die();}
}

void SelfContainedModel::print(bool add_lead) {
  if (add_lead) _LOG << "Performing self-contained gene-set analysis..." << endl;
  _LOG << "\tanalysing all available gene sets (total: " << analyse.size() << ")" <<  endl;
}

void SelfContainedModel::run() {
  SetResultsData results(analyse.size()); 
  results.block_init(SetResultsData::CoreBlock);

  SelfContainedComputation computer(variable_data);

  int success = 0;
  for (set<int>::iterator it = analyse.begin(); it != analyse.end(); ++it) {
    int row = results.add_row();

    try {computer.set_variable(*it);}
    catch (const SetException& se) {
      if (!computer.error_msg.empty()) {
        _LOG << "\tWARNING: " << computer.error_msg << "; skipping variable" << endl;
        continue;
      } else throw se;
    }

    results.var_name.set(row, computer.variable->get_name());     
    results.ngenes.set(row, computer.no_genes, true);

    try {computer.run();} 
    catch (const SetException& se) {
      if (!computer.error_msg.empty()) {
        _LOG << "\tWARNING: analysis failed for gene set '" << results.var_name.get(row) << "'; " << computer.error_msg << endl;
        continue;
      } else throw se;
    } 
    results.coefficient.set(row, computer.mu, true);
    results.pval.set(row, computer.pval, true);
    success++;
  }
  if (!success) throw SetException(error_tag, "analysis failed for all variables");

  results.lock_index();
  SetOutput output(results, file_prefix);
  output.set_abbreviate(abbr_length, abbr_file);
  int N = variable_data.mean_sampsize();
  if (N > 0) output.set_parameter("MEAN_SAMPLE_SIZE", N);  
  output.write(SetOutput::SelfContained, file_suffix);  
}


RegressionConfig::~RegressionConfig() {delete model_input;}

void RegressionConfig::init() {
  direction.assign(dt_COUNT, 1);
  direction[dt_GeneCovar] = 0;
  size_ss = 100; size_sc = 25;
  prop_ss = 0.1;
}

string RegressionConfig::parse_direction(int dir, bool verbose) {
  if (dir != 0) return string(verbose ? "one-sided, " : "") + (dir > 0 ? "positive" : "negative");    
  else return "two-sided"; 
}

set<int> RegressionConfig::filter_set(set<int>& target, set<int>& filter) {
  set<int>::iterator curr_t = target.begin(), curr_f = filter.begin(); set<int> removed;
  while (curr_t != target.end()) {
    while (curr_f != filter.end() && (*curr_f < *curr_t)) ++curr_f;
    if (curr_f == filter.end()) break;
    if (*curr_t == *curr_f) {removed.insert(*curr_t); target.erase(curr_t++);} ///post-increment to keep iterator valid after erase
    else ++curr_t;
  }
  return removed;
}

bool RegressionConfig::parse_variables(vector<string>& names, set<int>& target, const string& label, bool do_error) {
  map<string,ParseFailure> failed; 
  for (int i = 0; i < names.size(); i++) {string& var = names[i];
    set<int> ids = variable_data.variables.find_var(var, vc_AnyExternal);
    if (ids.size() != 1) {
      if (var.size() > 5) {
        if (Utils::starts_with(var, ":set:") || Utils::starts_with(var, ":SET:")) ids = variable_data.variables.find_var(var.substr(5), vc_GeneSet);
        else if (Utils::starts_with(var, ":cov:") || Utils::starts_with(var, ":COV:")) ids = variable_data.variables.find_var(var.substr(5), vc_GeneCovarExternal);
      }
      if (ids.empty()) {failed[var] = pf_NotFound; continue;}
      if (ids.size() > 1) {failed[var] = pf_Ambiguous; continue;}
    }
    if (!variable_data.get_variable(*ids.begin()).is_valid()) {failed[var] = pf_NotValid; continue;}
    target.insert(*ids.begin());
  }    
  
  if (!failed.empty()) {
    string prefix = do_error ? "\t" : "\tWARNING: ";
    if (do_error) _LOG.error(error_tag) << "errors occurred when parsing variables specified for " << label << endl;
    for (map<string,ParseFailure>::iterator it = failed.begin(); it != failed.end(); ++it) {
      _LOG << prefix << "variable '" << it->first << "' ";     
      switch (it->second) {
        case pf_NotFound: _LOG << "was not found in input" << endl; break;
        case pf_NotValid: _LOG << "was removed during preprocessing of input files" << endl; break;
        case pf_Ambiguous: _LOG << "occurs in multiple input files, please specify type" << endl; break;
      }
    }
    if (do_error) die();
    return false;
  }
  return true;
}

void RegressionConfig::set_internal(string mode, vector<string> filter) {
  if (mode == "include" && filter.empty()) mode = "none";
  if (mode == "exclude" && filter.empty()) mode = "all";  
  
  vector<string> internal; bool warning = mode == "include";
  if (mode != "none") {
    internal = InternalVariable::names();
    if (mode == "include" || mode == "exclude") {
      vector<short> index(internal.size(), false);
      for (int i = 0; i < filter.size(); i++) {
        InternalVariable::Type var_id = InternalVariable::type(filter[i], false);
        if (var_id < InternalVariable::iv_COUNT) index[var_id] = true;
        else {_LOG.error(error_tag) << "string '" << filter[i] << "' is not a recognized internal variable code" << endl; die();}
      }
      
      vector<string> retain; bool select = mode == "include";
      for (int i = 0; i < index.size(); i++) {if (index[i] == select) retain.push_back(internal[i]);}
      internal.swap(retain);
    }
  }

  for (int i = 0; i < internal.size(); i++) {bool avail = true;
    set<int> ids = variable_data.variables.find_var(internal[i], vc_Internal);
    for (set<int>::iterator it = ids.begin(); it != ids.end(); ++it) {if (!variable_data.get_variable(*it).is_valid()) avail = false;}
    if (!ids.empty() && avail) add_conditional(ids, ct_Internal);
    else if (warning) _LOG << "\tWARNING: cannot condition on requested internal variable " << internal[i] << ", variable is not available" << endl; 
  }
}

void RegressionConfig::set_conditional(vector<string> names, ConditionType type) {
  if (names.empty()) return;
  if (type == ct_Internal) {_LOG.error(error_tag) << "user-input variables cannot be used as internal covariate" << endl; die();}
  if (type == ct_ModeDependent) {_LOG.error(error_tag) << "mode-dependent covariates cannot be directly set" << endl; die();}
  if (type == ct_ConditionInteraction) {
    _LOG << "\tcreating interaction variables to condition on" << endl;
    set<set<int> > block_set; string type_tag = "conditioning on interaction";
    if (names.size() % 2 != 0) {_LOG.error(error_tag) << "conditioned-on interactions must be specified in pairs (specified number is not even)" << endl; die();}
    for (int i = 0; i < names.size(); i+= 2) {
      vector<string> curr(names.begin()+i, names.begin()+i+2); set<int> ids; 
      parse_variables(curr, ids, type_tag);
      if (ids.size() == 1) {_LOG.error(error_tag) << "cannot create interaction of variable with itself" << endl; die();}
      if (variable_data.get_variable(*ids.begin()).is_covar() && variable_data.get_variable(*ids.rbegin()).is_covar()) {
        _LOG.error(error_tag) << "cannot condition on interaction between two continuous gene properties" << endl; die();
      }
      block_set.insert(ids);        
    }

    vector<set<int> > blocks(block_set.begin(), block_set.end()); vector<int> interaction_ids;    
    InteractionBlocksInput(*this, blocks, interaction_ids);

    vector<int> failed;
    for (int i = 0; i < blocks.size(); i++) {
      if (interaction_ids[i] >= 0) {
        condition[ct_ConditionInteraction].insert(blocks[i].begin(), blocks[i].end());
        condition[ct_ConditionInteraction].insert(interaction_ids[i]);
      } else failed.push_back(i);
    }
    
    if (!failed.empty()) {
      _LOG.error(error_tag) << "following interactions failed to meet size and overlap specifications" << endl;
      for (int i = 0; i < failed.size(); i++) {set<int>& ids = blocks[failed[i]];
        _LOG << "\t" << variable_data.get_variable(*ids.begin()).get_name_label() << " -X- " << variable_data.get_variable(*ids.rbegin()).get_name_label() << endl;
      }
      die();
    }
  } else parse_variables(names, condition[type], "conditional analysis");
}

void RegressionConfig::set_filter(string mode) {
  variables.clear();
  if (mode == "all") variable_filter = vu_UseAll;
  else if  (mode == "sets") variable_filter = vu_SetsOnly;
  else if  (mode == "covar") variable_filter = vu_CovarOnly;  
  else if  (mode == "list") variable_filter = vu_List;
  else {_LOG.error(error_tag) << "unknown value '" << mode << "' for variable filter" << endl; die();}
}

void RegressionConfig::set_filter(vector<string> list) {
  variables.clear(); variable_filter = vu_List;
  parse_variables(list, variables, "variable filter list", false);  
  if (variables.empty()) {_LOG.error(error_tag) << "after applying variable filter list, no variables remain for analysis" << endl; die();}
} 

void RegressionConfig::set_mode(AnalysisMode mode) {
  analysis_mode = mode;
  variable_blocks.clear();
}

void RegressionConfig::set_models(vector<string> list, AnalysisMode mode) {set_mode(mode);
  if (analysis_mode == am_InteractionEach || analysis_mode == am_InteractionAll) parse_variables(list, variables_auxiliary, "interaction analyses");  
  else if (analysis_mode == am_Joint || analysis_mode == am_Interaction) process_specifications(list);
  else {_LOG.error(error_tag) << "invalid analysis mode type for function set_models(list, mode)" << endl; die();}
}


void RegressionConfig::process_specifications(vector<string>& list) {
  string type_tag = (analysis_mode == am_Interaction) ? "interaction model" : "multi-variable model";
  if (list.empty()) {_LOG.error(error_tag) << "list of " << type_tag << " specifications is empty" << endl; die();}
  int total = 0; bool warning[2] = {false,false};
  for (int i = 0; i < list.size(); i++) {
    vector<string> names = Utils::tokenize(list[i]); if (names.empty()) continue;
    set<int> ids; total++;
        
    if (parse_variables(names, ids, type_tag, false)) {
      if (analysis_mode == am_Interaction) {
        if (ids.size() != 2) {
          if (!warning[0]) {warning[0] = true; _LOG << "\tWARNING: interaction specifications must consist of exactly two variables" << endl;} 
          continue;
        }
        if (variable_data.get_variable(*ids.begin()).is_covar() && variable_data.get_variable(*ids.rbegin()).is_covar()) {
          if (!warning[1]) {warning[1] = true; _LOG << "\tWARNING: interactions between two continuous gene properties are currently not supported" << endl;} 
          continue;
        }
      }

      if (ids.size() < names.size()) {
        _LOG << "\tWARNING: following specification contains duplicate variable names; only single copy of each variable will be used" << endl;
        _LOG << "\t" << Utils::space("WARNING: ") << "specification: " << list[i] << endl;
      }
      variable_blocks.push_back(ids); 
    }
  }
  if (variable_blocks.empty()) {_LOG.error(error_tag) << "no valid " << type_tag << " specifications in input" << endl; die();}  
  _LOG << "\tfound valid specifications for " << Utils::plural(variable_blocks.size(), type_tag) << endl;
  if (variable_blocks.size() < total) _LOG << "\tremoved " << Utils::plural(total - variable_blocks.size(), "invalid model specification") << endl; 
}

void RegressionConfig::validate() {
  set<int> overlap;
  for (int ct = ct_ConditionForce; ct < ct_COUNT; ct++) {
    overlap = filter_set(condition[ct], condition[ct_Internal]);
    if (!overlap.empty()) {_LOG.error(error_tag) << "internal variable was included in user-specified conditional variables" << endl; die();}
  }

  filter_set(condition[ct_ConditionShow], condition[ct_ConditionInteraction]); 
  overlap = filter_set(condition[ct_ConditionInteraction], condition[ct_ConditionForce]);
  if (!overlap.empty()) {_LOG.error(error_tag) << "cannot condition on interactions involving variables included in 'condition-residualize'" << endl; die();}
  overlap = filter_set(condition[ct_ConditionShow], condition[ct_ConditionForce]);
  if (!overlap.empty()) _LOG << "\tWARNING: cannot produce output for conditioned-on variables included in 'condition-residualize'" << endl; 
  for (int ct = ct_ConditionForce; ct < ct_COUNT; ct++) {if (ct != ct_ConditionOnly) filter_set(condition[ct_ConditionOnly], condition[ct]);}

  set<int> conditioned; 
  for (int ct = 0; ct < ct_COUNT; ct++) conditioned.insert(condition[ct].begin(), condition[ct].end());

  if (analysis_mode == am_Joint || analysis_mode == am_Interaction) {variables.clear();
    string type_tag = (analysis_mode == am_Interaction) ? "interaction model" : "multi-variable model";
    int write = 0; bool all_uni = true;
    for (int read = 0; read < variable_blocks.size(); read++) {
      set<int>& curr = variable_blocks[read]; bool failed = false;
      for (set<int>::iterator it = curr.begin(); !failed && it != curr.end(); ++it) {if (conditioned.find(*it) != conditioned.end()) failed = true;}    
      if (!failed) {
        if (curr.size() > 1) all_uni = false;
        if (write < read) variable_blocks[write].swap(variable_blocks[read]);
        write++;              
      } 
    }
    if (write < variable_blocks.size()) _LOG << "\tfiltering out " << type_tag << "s containing conditioned-on variables" << endl;
    variable_blocks.resize(write);
    if (variable_blocks.empty()) {_LOG.error(error_tag) << "no models available for analysis after filtering" << endl; die();}

    if (all_uni) {
      for (int i = 0; i < variable_blocks.size(); i++) variables.insert(variable_blocks[i].begin(), variable_blocks[i].end());
      variable_blocks.clear(); 
      analysis_mode = am_Mono; variable_filter = vu_List;
    }    
  } else {variable_blocks.clear(); 
    if (variable_filter == vu_List) filter_set(variables, conditioned);
    else {
      VariableClass type = variable_filter == vu_SetsOnly ? vc_GeneSet : (variable_filter == vu_CovarOnly ? vc_GeneCovarExternal : vc_AnyExternal);
      variables = variable_data.variables.var_list(type, conditioned);
    }

    if (analysis_mode == am_InteractionEach || analysis_mode == am_InteractionAll) {
      if (!variables_auxiliary.empty()) {
        set<int> overlap = filter_set(variables_auxiliary, conditioned);
        if (!overlap.empty()) {_LOG.error(error_tag) << "conditioned-on variables cannot be used for interaction analysis" << endl; die();}      
        filter_set(variables, variables_auxiliary);
      } else {_LOG.error(error_tag) << "no interactor variables have been specified" << endl; die();}
    } else variables_auxiliary.clear();

    if (variables_auxiliary.empty()) {
      if (analysis_mode == am_InteractionEach || analysis_mode == am_InteractionAll) {_LOG.error(error_tag) << "no interactor variables have been specified" << endl; die();}
    } else {
      if (variables_auxiliary.size() == 1 && analysis_mode == am_InteractionEach) analysis_mode = am_InteractionAll;
      set<int> overlap = filter_set(variables_auxiliary, conditioned);
      if (!overlap.empty()) {
        if (analysis_mode == am_InteractionEach  || analysis_mode == am_InteractionAll) {_LOG.error(error_tag) << "conditioned-on variables cannot be used for interaction analysis" << endl; die();}      
      }
      if (analysis_mode == am_InteractionAll) condition[ct_ModeDependent] = variables_auxiliary;
    }
    if (variables.empty()) {_LOG.error(error_tag) << "no variables available for analysis with the specified model settings" << endl; die();}
  }
  
  if ((analysis_mode == am_JointPairs || analysis_mode == am_InteractionPairs) && variables.size() == 1) {_LOG.error(error_tag) << "there are no pairs to analyse, only a single variable is available for analysis" << endl; die();} 
 
  make_input();
}

void RegressionConfig::print() {
  _LOG << "\ttesting direction: " << parse_direction(direction[dt_GeneSet]) << " (sets), " << parse_direction(direction[dt_GeneCovar]) << " (covar)" << endl;
  if (has_interaction()) _LOG << "\ttesting direction interactions: " << parse_direction(direction[dt_InteractionSS]) << " (set x set), " << parse_direction(direction[dt_InteractionSC]) << " (set x cov)" << endl;

  if (!condition[ct_Internal].empty()) {
    _LOG <<"\tconditioning on internal variables:" << endl;
    vector<set<int> > types(InternalVariable::iv_COUNT); 
    for (set<int>::iterator it = condition[ct_Internal].begin(); it != condition[ct_Internal].end(); ++it) {
      VariableInfo& info = variable_data.get_variable(*it); InternalVariable::Type type;
      if (info.has_type(vc_Internal) && (type = info.as_covar()->as_internal()->get_subtype()) < InternalVariable::iv_COUNT) types[type].insert(*it);
      else _LOG.error(error_tag) << "user-input variables cannot be used as internal covariate" << endl;       
    } 
    for (int i = 0; i < types.size(); i++) {
      if (!types[i].empty()) {_LOG << "\t\t";
        for (set<int>::iterator it = types[i].begin(); it != types[i].end(); ++it) _LOG << (it != types[i].begin() ? ", " : "") << variable_data.get_variable(*it).as_covar()->as_internal()->format_name();
        _LOG << endl;
      }
    }
  }

  for (int ct = ct_ConditionForce; ct < ct_COUNT; ct++) {
    if (!condition[ct].empty()) {
      switch (ct) {
        case ct_ConditionForce: _LOG << "\tconditioning on input variables (residualizing effects):" << endl; break;
        case ct_ConditionOnly: _LOG << "\tconditioning on input variables (no output):" << endl; break;
        case ct_ConditionShow: _LOG << "\tconditioning on input variables:" << endl; break;
        case ct_ConditionInteraction: _LOG << "\tconditioning on interactions (and corresponding main effects):" << endl; break;                
        default: continue;
      }
        
      for (set<int>::iterator it = condition[ct].begin(); it != condition[ct].end(); ++it) {
        VariableInfo& info = variable_data.get_variable(*it);     
        if (ct == ct_ConditionInteraction && !info.is_interaction()) continue;
        _LOG << "\t\t" << info.get_name_label() << endl;
      }
    }
  }

  if (!variables.empty() && variable_filter != vu_UseAll) {
    switch (variable_filter) {
      case vu_SetsOnly: _LOG << "\trestricting analysis to gene set variables" << endl; break;
      case vu_CovarOnly: _LOG << "\trestricting analysis to gene covariate variables" << endl; break;
      case vu_List: _LOG << "\trestricting analysis to specified list" << endl; break;
      default:;
    }
  }

  switch (analysis_mode) {
    case am_Mono: _LOG << "\tanalysing individual variables" << endl; break;
    case am_Joint: _LOG << "\tanalysing specified multi-variable models" << endl; break;
    case am_JointPairs: _LOG << "\tanalysing all pairs of available variables" << endl; break;
    case am_Interaction: _LOG << "\tanalysing specified (valid) interactions" << endl; break;
    case am_InteractionEach: _LOG << "\tanalysing (valid) interactions for all variables with each of specified variables" << endl; break;
    case am_InteractionPairs: _LOG << "\tanalysing (valid) interactions between all pairs of available variables" << endl; break;
    case am_InteractionAll: _LOG << "\tanalysing (valid) multi-interaction models for each variable with:" << endl; break;
  }
  if (analysis_mode == am_InteractionAll) {
    for (set<int>::iterator it = variables_auxiliary.begin(); it != variables_auxiliary.end(); ++it) {
      VariableInfo& info = variable_data.get_variable(*it);     
      _LOG << "\t\t" << info.get_name_label() << endl;
    }
  }
}

void RegressionConfig::make_input() {
  delete model_input; model_input = 0;
  switch (analysis_mode) {
    case am_Mono: model_input = new MonoInput(*this); break;
    case am_Joint: model_input = new JointInput(*this); break;
    case am_JointPairs: model_input = new JointPairsInput(*this); break;    
    case am_Interaction: model_input = new InteractionBlocksInput(*this); break;
    case am_InteractionAll: model_input = new InteractionAllInput(*this); break;    
    case am_InteractionEach: model_input = new InteractionEachInput(*this); break;    
    case am_InteractionPairs: model_input = new InteractionPairsInput(*this); break;
  }
}


RegressionConfig::BaseInput& RegressionConfig::get_input() {
  if (!model_input) throw SetException(error_tag, "no input is specified for regression analysis");  
  return *model_input;
}

bool RegressionModel::check_df(int df, int eff_size) {
  if (df >= MIN_DF) {
    double prop = df / double(eff_size);
    if (prop < 0.75) {string msg;
      if (prop < 0.25) msg = "WARNING: model degrees of freedom is below 25%";
      else msg = "NOTE: model degrees of freedom is below 75%";
      _LOG << "\t" << msg << "; effective data size = " << eff_size << ", df = " << df << " (" << setprecision(4) << prop*100 << "%)" << endl;
    }
    return true;
  } else {_LOG << "\tERROR: insufficient degrees of freedom to run analyses; aborting analysis of remaining models" << endl; return false;}
}

void RegressionModel::check_rsq(double rsq) {
  if (rsq <= MAX_RSQ) {
    if (rsq > 0.25) {
      if (rsq > 0.75) _LOG << "\tWARNING: conditioned-on variables explain a very large percentage of variance of genetic associations";      
      else _LOG << "\tNOTE: conditioned-on variables explain a large percentage of variance of genetic associations";      
      _LOG << " (" << setprecision(3) << rsq*100 << "%)" << endl;
    }
  } else {
    if (rsq > 1) rsq = 1; string perc = Utils::num_string(rsq*100, 3); bool total = perc == "100";
    string msg = string("conditioned-on variables explain") + (!total ? " almost" : "") + " all variance of genetic associations";
    if (!total) msg += " (" + perc + "%)";    
    throw SetException(error_tag, msg);
  }
}   

void RegressionModel::create_terms(vector<VariableInfo*> vars, vector<string>& target) {static vector<string> letters;   
  if (letters.size() < vars.size()) {letters.resize(vars.size()); for (int i = 0; i < vars.size(); i++) letters[i] = Utils::get_letters(i, true);}
  target.resize(vars.size()); vector<int> interact; int count = 0;
  for (int i = 0; i < vars.size(); i++) {
    if (vars[i]->is_interaction()) interact.push_back(i);
    else target[i] = letters[count++]; 
  }
  for (int i = 0; i < interact.size(); i++) {
    set<int> ids = vars[interact[i]]->get_interaction()->ids(); string first = "?", second = "?";
    for (int j = 0; j < vars.size(); j++) {int id = vars[j]->get_id();
      if (id == *ids.begin()) first = target[j];
      else if (id == *ids.rbegin()) second = target[j];
    }    
    target[interact[i]] = first + "*" + second;        
  }
}

void RegressionModel::set_gene_info(RegressionComputation& computer) {
  GeneData& gene_info = variable_data.get_genes();
  gene_info.create_pval();

  Buffer<double> buffer; computer.load_fitted(buffer, RegressionComputation::mt_BaseModel);
  BlockMultiVariable<ZscoreVariable<double> >* var = new BlockMultiVariable<ZscoreVariable<double> >();
  var->set_width(2); gene_info.attach_variable(var, "ZSTAT_FITTED_BASE");
  var->get_buffer().assign(buffer);

  var = new BlockMultiVariable<ZscoreVariable<double> >();  
  var->set_width(1); gene_info.attach_variable(var, "ZSTAT_OUTCOME");
  var->get_buffer().assign(computer.get_outcome().data(), gene_info.data_used(false));  
}

void RegressionModel::print(bool add_lead) {
  if (add_lead) _LOG << "Performing regression analysis..." << endl;
  config.print();
}

void RegressionModel::run() {
  RegressionConfig::BaseInput& input = config.get_input();

  int no_tests = input.no_tests(), no_analyses = input.no_analyses(), success = 0, model = 1, block_index = 0;
  bool show_model = input.var_count(false).second > 1, show_interaction = config.has_interaction();
  if (no_tests == 0) throw SetException(error_tag, "model has nothing to analyse");

  RegressionComputation computer(variable_data, config.get_residualize(), config.get_collinear());
  computer.set_conditional(input.get_conditional(RegressionConfig::ct_Internal), RegressionComputation::pt_InternalCovariate);
  computer.set_conditional(input.get_conditional(RegressionConfig::ct_ConditionForce), RegressionComputation::pt_InternalCovariate); 
  computer.set_conditional(input.get_conditional(RegressionConfig::ct_ConditionOnly), RegressionComputation::pt_HiddenCovariate);  
  computer.set_conditional(input.get_conditional(RegressionConfig::ct_ConditionShow), RegressionComputation::pt_UsedCovariate);  
  computer.set_conditional(input.get_conditional(RegressionConfig::ct_ConditionInteraction), RegressionComputation::pt_UsedCovariate);    
  computer.set_conditional(input.get_conditional(RegressionConfig::ct_ModeDependent), RegressionComputation::pt_UsedCovariate); 
  computer.init(input.var_count(true).first);

  check_rsq(computer.get_rsq(true));
  _LOG << endl;

  SetResultsData results(no_tests);
  results.block_init(SetResultsData::RegressionBlock);
  if (show_model) results.model_id.init();
  if (show_interaction) results.model_term.init();
  
  set_gene_info(computer);
  SetGenesOutput set_genes(variable_data, file_prefix, file_suffix + (show_interaction ? ".inter" : ".sets") + ".genes");
  set_genes.set_signif(no_analyses, alpha); 
  double thresh = alpha / no_analyses;
                  
  short type_code[2][2] = {
    {results.var_type.get_code("COVAR"), results.var_type.get_code("SET")},
    {results.var_type.get_code("INTER-SC"), results.var_type.get_code("INTER-SS")}
  };

  RegressionConfig::InputIterator* block = 0; vector<string> model_terms;
  for (; (block = input.get_block()); input.next(), block_index++) {
    Timer timer;
    int main_var = block->no_vars(true), total_var = block->no_vars(false), counter = 0;
    double beta[total_var], beta_std[total_var], std_error[total_var], pval[total_var]; 
    bool size_change = computer.update_nvar(main_var);

    block->print();     
    if ((block_index == 0 || size_change) &&!check_df(computer.get_df(), computer.get_size(false))) break; ///assumes that blocks are ordered by increasing df

    for (; !block->empty(); block->next(), model++) {
      int row_offset = results.add_block(total_var);

      try {computer.set_variables(block->get_variables());}
      catch (const SetException& se) {
        if (!computer.error_msg.empty()) {
          _LOG << "\tWARNING: " << computer.error_msg << "; skipping variable" << endl;
          continue;
        } else throw se;
      }

      if (show_interaction) create_terms(computer.get_variables(), model_terms);      
      for (int i_var = 0; i_var < total_var; i_var++) {
        VariableInfo& var = computer.get_variable(i_var);    
        if (show_model) results.model_id.set(row_offset+i_var, model);     
        if (show_interaction) results.model_term.set(row_offset+i_var, model_terms[i_var]);
        if (var.is_interaction()) results.var_name.set(row_offset+i_var, var.get_interaction()->get_name());     
        else results.var_name.set(row_offset+i_var, var.get_name());     
        results.var_type.set(row_offset+i_var, type_code[var.is_interaction()][var.is_set()]);
        results.ngenes.set(row_offset+i_var, computer.get_ngenes(i_var), true);
      }

      try {computer.run();}
      catch (const SetException& se) {
        if (!computer.error_msg.empty()) {
          if (main_var == 1) {
            VariableInfo& var = computer.get_variable(0, RegressionComputation::pt_MainVariable);
            _LOG << "\tWARNING: analysis failed for '" << var.get_name() << "' " << var.get_label() << "; " << computer.error_msg << endl;
          } else _LOG << "\tWARNING: analysis failed for model " << model << "; " << computer.error_msg << endl;
          continue;
        } else throw se;
      }   

      bool skip_model = false;
      for (int i_var = 0; i_var < total_var; i_var++) {
        VariableInfo& var = computer.get_variable(i_var);    
        try {   
          beta[i_var] = computer.get_coeff(i_var); 
          beta_std[i_var] = computer.get_coeff_std(i_var);
          std_error[i_var] = computer.get_se(i_var);
          pval[i_var] = computer.get_pval(i_var, input.get_direction(var)); 
        } catch (const SetException& se) {
          if (!computer.error_msg.empty()) {
            string model_tag = (total_var > 1) ? Utils::num_string_pad(" in model ", model) : "";
            _LOG << "\tWARNING: analysis failed for '" << var.get_name() << "' " << var.get_label() << model_tag << "; " << computer.error_msg << endl;
            skip_model = true; break;
          } else throw se;
        } 
      }
      if (skip_model) continue;

      for (int i_var = 0; i_var < total_var; i_var++) {
        VariableInfo& var = computer.get_variable(i_var);    
        results.coefficient.set(row_offset+i_var, beta[i_var], true);
        results.coefficient_std.set(row_offset+i_var, beta_std[i_var], true);          
        results.std_error.set(row_offset+i_var, std_error[i_var], true);
        results.pval.set(row_offset+i_var, pval[i_var], true);
        if (computer.get_collinear(i_var)) results.warnings.set(row_offset+i_var, "COLLINEAR");

        if (pval[i_var] < thresh && computer.has_type(i_var, RegressionComputation::pt_MainVariable)) {   
          if ((main_var == 1 && var.is_set()) || var.is_interaction()) { 
            try {
              set_genes.write_set(var, pval[i_var], show_model ? model : 0);
            } catch (const BaseException& be) {
              _LOG << "\tWARNING: error when writing gene analysis results for significant gene set '" << var.get_name() << "; " << be.get_msg() << endl;
              continue;
            } 
          }
        }
      }
      counter++;
      if (timer.check() && counter > 0) {_LOG.counter(true) << "\tmodels analysed: " << counter << "           "; _LOG.counter(false);}
    }
    success += counter;
    _LOG.wipe_counter();    
  }
  if (!success) throw SetException(error_tag, "analysis failed for all variables/models");

  output_results(results, print_genes || set_genes.get_count() > 0);
  set_genes.print_log();
}


void RegressionModel::output_results(SetResultsData& results, bool print_genes) {
  RegressionConfig::BaseInput& input = config.get_input();
  results.lock_index();
  
  SetOutput output(results, file_prefix);
  output.set_abbreviate(abbr_length, abbr_file);

  GeneOutput gene_output(variable_data.get_genes(), file_prefix);
  gene_output.set_option(OutputOptions::ShowFitted, true);
  
  int N = variable_data.mean_sampsize();
  if (N > 0) output.set_parameter("MEAN_SAMPLE_SIZE", N);  
  output.set_parameter("TOTAL_GENES", variable_data.no_genes());
  
  vector<string> test_dir(2); 
  test_dir[0] = RegressionConfig::parse_direction(input.get_direction(RegressionConfig::dt_GeneSet)) + " (set)";
  test_dir[1] = RegressionConfig::parse_direction(input.get_direction(RegressionConfig::dt_GeneCovar)) + " (covar)";
  if (config.has_interaction()) {
    test_dir.resize(4);
    test_dir[2] = RegressionConfig::parse_direction(input.get_direction(RegressionConfig::dt_InteractionSS)) + " (set x set)";
    test_dir[3] = RegressionConfig::parse_direction(input.get_direction(RegressionConfig::dt_InteractionSC)) + " (set x cov)";  
  }
  output.set_parameter("TEST_DIRECTION", test_dir);

  set<int>& internal = input.get_conditional(RegressionConfig::ct_Internal); vector<string> var_names, var_hidden, var_forced;
  for (set<int>::iterator it = internal.begin(); it != internal.end(); ++it) var_names.push_back(InternalVariable::format(variable_data.get_variable(*it).get_name()));  
  if (!var_names.empty()) {
    output.set_parameter("CONDITIONED_INTERNAL", var_names); 
    gene_output.set_parameter("CONDITIONED_INTERNAL", var_names); 
    var_names.clear();
  }

  for (int c = RegressionConfig::ct_ConditionForce; c < RegressionConfig::ct_COUNT; c++) {
    set<int>& curr = input.get_conditional(static_cast<RegressionConfig::ConditionType>(c));
    for (set<int>::iterator it = curr.begin(); it != curr.end(); ++it) {
      VariableInfo& info = variable_data.variables.var_info(*it);      
      if (c == RegressionConfig::ct_ConditionForce) var_forced.push_back(info.get_name_label());
      else if (c == RegressionConfig::ct_ConditionOnly) var_hidden.push_back(info.get_name_label());
      else var_names.push_back(info.get_name_label());
    }
  } 
  if (!var_forced.empty()) {output.set_parameter("CONDITIONED_RESIDUALIZED", var_forced); gene_output.set_parameter("CONDITIONED_RESIDUALIZED", var_forced);} 
  if (!var_hidden.empty()) {output.set_parameter("CONDITIONED_HIDDEN", var_hidden); gene_output.set_parameter("CONDITIONED_HIDDEN", var_hidden);} 
  if (!var_names.empty()) {output.set_parameter("CONDITIONED_VARIABLES", var_names); gene_output.set_parameter("CONDITIONED_VARIABLES", var_names);}
  
  output.write(SetOutput::Regression, file_suffix);  
  if (print_genes) gene_output.write_outfile("gene information", file_suffix + ".genes");
}

short RegressionConfig::BaseInput::get_direction(VariableInfo& var) {
  if (var.is_interaction()) return config.direction[var.is_set() ? dt_InteractionSS : dt_InteractionSC];
  else return config.direction[var.is_set() ? dt_GeneSet : dt_GeneCovar];
}

pair<int,int> RegressionConfig::MonoInput::var_count(bool main_only) {
  int nvar = !main_only ? no_cond+1 : 1;
  return pair<int,int>(nvar,nvar);
}

void RegressionConfig::MonoInput::MonoIterator::print() {
  _LOG << "\tanalysing single-variable models (number of models: " << analyse.size() << ")" << endl;
}

void RegressionConfig::JointInput::process_blocks() {
  for (int i = 0; i < variable_blocks.size(); i++) {
    int size = variable_blocks[i].size();
    index[size].push_back(i);
    ntest += size;
  }

  if (!index.empty()) {
    range = pair<int,int>(index.begin()->first, index.rbegin()->first);
    current = index.begin();
  } else current = index.end();

  set_iterator();
}

void RegressionConfig::JointInput::set_iterator() {
  delete iterator; 
  if (current != index.end()) iterator = new JointIterator(variable_blocks, current->second, current->first, no_cond);
  else iterator = 0;  
}

pair<int,int> RegressionConfig::JointInput::var_count(bool main_only) {
  if (!main_only) return pair<int,int>(range.first+no_cond, range.second+no_cond);
  else return range;
}

void RegressionConfig::JointInput::JointIterator::print() {
  string amount = (nvar > 1) ? Utils::num_name(nvar) + " main variables" : "one main variable";
  _LOG << "\tanalysing models with " << amount << " per model (number of models: " << analyse.size() << ")" << endl;
}

void RegressionConfig::JointPairsInput::make_pairs() {
  set<int>& vars = config.variables; 
  if (vars.size() > 1) {
    pair_blocks.reserve(0.5 * vars.size() * (vars.size()-1));
    for (set<int>::iterator curr = vars.begin(); curr != vars.end(); ++curr) {
      for (set<int>::iterator other = curr; other != vars.end(); ++other) {      
        if (*curr < *other) {
          set<int> ids; ids.insert(*curr); ids.insert(*other);
          pair_blocks.push_back(ids);
        }
      }
    }
  }
  process_blocks();
}

void RegressionConfig::InteractionInput::process_blocks() {
  vector<PairSS> ss_blocks; vector<PairSC> sc_blocks;
  sort_blocks(ss_blocks, sc_blocks);

  if (!ss_blocks.empty()) counts[vt_GeneSet] = process_ss(ss_blocks);
  if (!sc_blocks.empty()) counts[vt_GeneCovar] = process_sc(sc_blocks);

  int total = counts[0] + counts[1];
  if (verbose) {
    _LOG << "\tfiltering interaction pairs on size and overlap specifications and removing duplicates" << endl;
    _LOG << "\t\tset x set: at least " << ss_part << " genes overlapping, and at least " << ss_part << " genes unique to each set" << endl;
    _LOG << "\t\tset x set: for each set, between " << ss_prop*100 << "% and " << (1-ss_prop)*100 << "% of genes in overlap with other set" << endl;
    _LOG << "\t\tset x covariate: between " << sc_min << " and " << sc_max << " genes in set" << endl;
    if (!external_use) {
      if (total > 0) _LOG << "\tafter filtering, " << Utils::plural(total, "interaction pair") << " remaining" << endl;
      else {_LOG.error(config.error_tag) << "no interactions left to analyse" << endl; die();}
    }
  }

  if (!external_use) set_iterator();
}

RegressionConfig::InteractionInput::InteractionInput(RegressionConfig& config, bool external, bool linked) : BaseInput(config), set_info(0), current(0), iterator(0), external_use(external || linked), verbose(!external_use), linked(linked), delete_info(true) {
  logged[0] = false; logged[1] = false; counts[0] = 0; counts[1] = 0;
  int no_genes = config.variable_data.no_genes(); 
  ss_part = config.size_ss; ss_prop = config.prop_ss; sc_min = config.size_sc;
  ss_min = 2*ss_part; ss_max = no_genes - ss_part; sc_max = no_genes - sc_min; 
}

int RegressionConfig::InteractionInput::process_ss(vector<PairSS>& blocks) {
  set_info = new GeneSetData(); set_info->init(1000); 
  set_info->link_data(config.variable_data.get_genes());
  NumericID<int> *main1 = new NumericID<int>(), *main2 = new NumericID<int>(), *ext_id = 0;
  set_info->attach_variable(main1, "INTERACT_MAIN1");
  set_info->attach_variable(main2, "INTERACT_MAIN2");
  if (external_use) {ext_id = new NumericID<int>(); set_info->attach_variable(ext_id, "EXTERNAL_ID");}

  for (int b = 0; b < blocks.size(); b++) {
    SetInfo &set1_var = *blocks[b].set1, &set2_var = *blocks[b].set2;
    int size1 = set1_var.get_ngenes(), size2 = set2_var.get_ngenes();  
    if (size1 >= ss_min && size1 <= ss_max && size2 >= ss_min && size2 <= ss_max) {
      set<int> &data1 = *set1_var.get_data(), &data2 = *set2_var.get_data(), intersect;
      set<int>::iterator curr1 = data1.begin(), curr2 = data2.begin();
      while (curr1 != data1.end() && curr2 != data2.end()) {
        if (*curr1 < *curr2) {while (curr1 != data1.end() && *curr1 < *curr2) ++curr1;}
        else if (*curr2 < *curr1) {while (curr2 != data2.end() && *curr2 < *curr1) ++curr2;}
        else {intersect.insert(*curr1); ++curr1; ++curr2;}
      }
  
      int size12 = intersect.size(), part1 = max(ss_part, int(round(ss_prop*size1))), part2 = max(ss_part, int(round(ss_prop*size2)));
      if (size12 >= part1 && (size1 - size12) >= part1 && size12 >= part2 && (size2 - size12) >= part2) {
        string name = string(set1_var.get_name()) + " -X- " + set2_var.get_name();
        int set_id = set_info->get_iid(name), col_id = set_info->data.insert(intersect);
  
        set_info->mapped_genes.set(set_id, size12);    
        set_info->used_genes.set(set_id, size12);    
        set_info->storage_column.set(set_id, col_id);
        main1->set(set_id, set1_var.get_id());
        main2->set(set_id, set2_var.get_id()); 
        if (ext_id) ext_id->set(set_id, blocks[b].id);   
      }
    }
  }
 
  set_info->build_index();
  return set_info->data_used(false);
}

GeneCovarData<float>* RegressionConfig::InteractionInput::add_covar_data(NumericID<int>*& main1, NumericID<int>*& main2) {
  GeneCovarData<float>* add_info = new GeneCovarData<float>(); add_info->init(1000);  
  main1 = new NumericID<int>(); add_info->attach_variable(main1, "INTERACT_MAIN1");
  main2 = new NumericID<int>(); add_info->attach_variable(main2, "INTERACT_MAIN2");
  covar_info.push_back(add_info); return covar_info.back();
}

int RegressionConfig::InteractionInput::process_sc(vector<PairSC>& blocks) {
  GeneCovarData<float>* curr_info = 0; NumericID<int> *main1 = 0, *main2 = 0, *ext_id = 0;
  int column = 0, total = 0, prev_ext = -1; bool block_swap = false; 

  for (int b = 0; b < blocks.size(); b++) {
    SetInfo& set_var = *blocks[b].set; CovarInfo& cov_var = *blocks[b].cov;
    int set_size = set_var.get_ngenes();
    if (set_size >= sc_min && set_size <= sc_max) {
      if ((column++ % 1000) == 0) block_swap = true;      
      if (block_swap && (!linked || blocks[b].id < 0 || blocks[b].id != prev_ext)) {
        block_swap = false;
        curr_info = add_covar_data(main1, main2); 
        if (external_use) {
          ext_id = new NumericID<int>(); 
          curr_info->attach_variable(ext_id, "EXTERNAL_ID");
        }  
      }
      string name = string(set_var.get_name()) + " -X- " + cov_var.get_name();
      int var_id = curr_info->get_iid(name);

      curr_info->values_observed.set(var_id, set_size);  
      curr_info->storage_column.set(var_id, var_id);
      main1->set(var_id, set_var.get_id());
      main2->set(var_id, cov_var.get_id()); 
      if (ext_id) {
        ext_id->set(var_id, blocks[b].id);       
        prev_ext = blocks[b].id;
      }
    }
  }

  for (int i = 0; i < covar_info.size(); i++) {covar_info[i]->build_index(); total += covar_info[i]->data_used(false);}
  return total;
}

void RegressionConfig::InteractionInput::create_covar_data(GeneCovarData<float>* info) {
  SetVariables& var_data = config.variable_data;
  double n_scale = var_data.no_genes() - 1;
  info->link_data(var_data.get_genes(), true);
  info->block_init(GeneCovarData<float>::StatsBlock);

  for (int i = 0; i < info->data_total(); i++) {bool success = false;
    if (info->active_id(i)) {
      VariableInfo& var = var_data.get_variable(info->variable_id.get(i), true);

      if (var.has_type(vc_InteractionSC)) {
        CovarInfo &inter = *var.as_covar(), *cov_main = inter.get_interaction()->cov; 
        set<int>& set_data = *inter.get_interaction()->set->get_data(); 
        float *cov_data = cov_main->get_data(), *inter_data = inter.get_data();

        double mean = 0, var = 0; int K = set_data.size();
        for (set<int>::iterator it = set_data.begin(); it != set_data.end(); ++it) mean += cov_data[*it];        
        mean /= K;
        for (set<int>::iterator it = set_data.begin(); it != set_data.end(); ++it) {double diff = cov_data[*it] - mean; var += diff*diff;}
        var /= K - 1;                

        double sd = var > 0 ? sqrt(var) : 0;
        if (sd > 1e-5) {
          double sd = sqrt(var), scale = sqrt(n_scale/(K-1))/sd, cov_sd = cov_main->get_orig_sd();
          info->orig_mean.set(i, cov_sd*mean + cov_main->get_orig_mean());
          info->orig_sd.set(i, cov_sd*sd*sqrt((K-1)/n_scale));
          for (set<int>::iterator it = set_data.begin(); it != set_data.end(); ++it) inter_data[*it] = scale*(cov_data[*it] - mean);
          success = true;
        } else if (verbose && !external_use) _LOG << "\tWARNING: skipping interaction " << inter.get_name() << "; variance is too low" << endl;
      }
    }
    if (!success) info->variable_id.set(i, -1);
  }
}

void RegressionConfig::InteractionInput::make_iterator(VariableIndex* input, VariableType type) {
  NumericID<int>* main1 = input->typed_variable<NumericID<int> >("INTERACT_MAIN1");
  NumericID<int>* main2 = input->typed_variable<NumericID<int> >("INTERACT_MAIN2");

  vector<set<int> > blocks; blocks.reserve(input->data_used(true));
  for (int i = 0; i < input->data_total(); i++) {
    if (input->active_id(i)) {set<int> curr;
      curr.insert(input->variable_id.get(i));
      curr.insert(main1->get(i)); curr.insert(main2->get(i));                  
      if (curr.size() == 3 && curr.find(-1) == curr.end()) blocks.push_back(curr);
    }
  }

  string type_str = (type == vt_GeneSet) ? "set by set" : "set by covariate";
  string msg = !logged[type] ? string("analysing ") + type_str + " interaction models (number of models: " + Utils::num_string(counts[type]) + ")" : "";
  iterator = new InteractionIterator(3, no_cond, msg); 
  iterator->set_blocks(blocks);

  if (!iterator->empty()) logged[type] = true;
  current = input;
}

void RegressionConfig::InteractionInput::set_iterator() {
  delete iterator; iterator = 0;
  if (current) {
    if (!external_use) current->unlink_data(true); 
    current = 0;
  }

  if (set_info) {
    config.variable_data.add_interactions(set_info);
    make_iterator(set_info, vt_GeneSet); set_info = 0;
  } else if (!covar_info.empty()) {
    GeneCovarData<float>* curr_info = covar_info[0];
    config.variable_data.add_interactions(curr_info);
    create_covar_data(curr_info);
    make_iterator(curr_info, vt_GeneCovar);
    covar_info.erase(covar_info.begin());
  }

  if (iterator && iterator->empty()) set_iterator();
}

pair<int,int> RegressionConfig::InteractionInput::var_count(bool main_only) {
  int nvar = !main_only ? no_cond + 3 : 3;
  return pair<int,int>(nvar,nvar);
}

void RegressionConfig::InteractionBlocksInput::sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks) {
  for (int i = 0; i < blocks.size(); i++) {set<int>& curr = blocks[i];
    if (curr.size() == 2) {
      VariableInfo &var1 = config.variable_data.get_variable(*curr.begin(), true), &var2 = config.variable_data.get_variable(*curr.rbegin(), true);
      if (var1.is_valid() && var2.is_valid()) {
        if (var1.is_set() && var2.is_set()) ss_blocks.push_back(PairSS(var1.as_set(), var2.as_set(), i));
        else if (var1.is_set() && var2.is_covar()) sc_blocks.push_back(PairSC(var1.as_set(), var2.as_covar(), i));
        else if (var1.is_covar() && var2.is_set()) sc_blocks.push_back(PairSC(var2.as_set(), var1.as_covar(), i));
      }
    }
  }
}

void RegressionConfig::InteractionBlocksInput::extract(vector<int>& target) {
  target.assign(blocks.size(), -1);
  for (; current != 0; next()) {
    NumericID<int> &var_id = current->variable_id, *ext_id = current->typed_variable<NumericID<int> >("EXTERNAL_ID");
    if (var_id.initialised() && ext_id) {
      for (int v = 0; v < current->data_total(); v++) {
        if (current->active_id(v)) {
          int ext = ext_id->get(v);
          if (ext >= 0) target[ext] = var_id.get(v);
        }
      }
    }
  }
}


RegressionConfig::InteractionAllInput::InteractionAllInput(RegressionConfig& config) : InteractionInput(config, true, config.variables_auxiliary.size() > 1), analysis_count(0), current_block(0) {
  interaction_count = config.variables_auxiliary.size(); verbose = true; 
  process_blocks(); delete_info = false;
  extract();
 
  if (analysis_count > 0) _LOG << "\tafter filtering, " << Utils::plural(analysis_count, "multi-interaction model") << " remaining" << endl;
  else {_LOG.error(config.error_tag) << "no interaction models left to analyse" << endl; die();}
 
  set_iterator();
}  

void RegressionConfig::InteractionAllInput::sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks) {
  set<int> &primary = config.variables, &secondary = config.variables_auxiliary; ///at this stage, no overlap between primary and secondary anymore
  
  bool has_covar = false;
  for (set<int>::iterator sec = secondary.begin(); sec != secondary.end(); ++sec) {
    VariableInfo& var = config.variable_data.get_variable(*sec, true);    
    if (!var.is_valid()) return;
    if (var.is_covar()) has_covar = true;
  }

  for (set<int>::iterator pri = primary.begin(); pri != primary.end(); ++pri) {
    VariableInfo& main = config.variable_data.get_variable(*pri, true);    
    if (main.is_valid() && (!has_covar || !main.is_covar())) {
      for (set<int>::iterator sec = secondary.begin(); sec != secondary.end(); ++sec) {       
        VariableInfo& inter = config.variable_data.get_variable(*sec, true);          
        if (main.is_set() && inter.is_set()) ss_blocks.push_back(PairSS(main.as_set(), inter.as_set(), *pri));
        else if (main.is_set() && inter.is_covar()) sc_blocks.push_back(PairSC(main.as_set(), inter.as_covar(), *pri));
        else if (main.is_covar() && inter.is_set()) sc_blocks.push_back(PairSC(inter.as_set(), main.as_covar(), *pri));
      }
    }
  }
}

void RegressionConfig::InteractionAllInput::strip_storage(VariableIndex* info, vector<int>& used) { 
  NumericID<int> &storage = info->storage_column, *main_variable = info->typed_variable<NumericID<int> >("EXTERNAL_ID");
  
  if (main_variable) {
    int miss_code = storage.get_missing().second;
    for (int v = 0; v < info->data_total(); v++) {
      if (info->active_id(v) && !used[main_variable->get(v)]) storage.set(v, miss_code);
    }
  } else storage.assign_missing();
}

void RegressionConfig::InteractionAllInput::extract() {
  int range = *config.variables.rbegin() + 1;
  vector<int> index(range, -1), block_id(range, -1); vector<set<int> > block_storage;

  for (int b = -1; b < int(covar_info.size()); b++) {VariableIndex* info = 0; 
    if (b < 0) {
      if (set_info) {config.variable_data.add_interactions(set_info); info = set_info;}
      else continue;
    } else {config.variable_data.add_interactions(covar_info[b]); info = covar_info[b];}
 
    NumericID<int> &interaction = info->variable_id, *main_variable = info->typed_variable<NumericID<int> >("EXTERNAL_ID");
    if (interaction.initialised() && main_variable) {
      for (int v = 0; v < info->data_total(); v++) {
        if (info->active_id(v)) {
          int var_id = main_variable->get(v), inter_id = interaction.get(v);
          if (var_id >= 0 && inter_id >= 0) {
            if (index[var_id] >= 0) {set<int>& ids = block_storage[index[var_id]];
              if (b >= 0 && block_id[var_id] != b) {
                if (block_id[var_id] >= 0) throw SetException("processing interactions", "incomplete variable data block in function InteractionAllInput::extract()");
                else  block_id[var_id] = b;
              }
              ids.insert(inter_id);
            } else {
              set<int> ids; ids.insert(var_id); ids.insert(inter_id);
              index[var_id] = block_storage.size(); block_storage.push_back(ids);
              block_id[var_id] = b;
            }
          }
        }
      }
    }
  }

  blocks.resize(covar_info.size()+1);
  int required = interaction_count + 1; analysis_count = 0;
  for (int i = 0; i < index.size(); i++) {
    if (index[i] >= 0) {set<int>& ids = block_storage[index[i]];
      if (ids.size() == required) {
        blocks[block_id[i]+1].push_back(ids);
        analysis_count++; index[i] = true;
      } else index[i] = false;
    } else index[i] = false;
  }
  if (set_info) strip_storage(set_info, index);
  for (int i = 0; i < covar_info.size(); i++) strip_storage(covar_info[i], index);
} 

void RegressionConfig::InteractionAllInput::set_iterator() {
  delete iterator; iterator = 0;
  string msg = (current_block == 0) ? string("analysing multi-interaction models (number of models: ") + Utils::num_string(analysis_count) + ")" : "";

  while (current_block < blocks.size() && blocks[current_block].empty()) current_block++;
  if (current_block < blocks.size()) {
    iterator = new InteractionIterator(interaction_count+1, no_cond, msg);   
    iterator->set_blocks(blocks[current_block]);
    if (current_block > 0) create_covar_data(covar_info[current_block-1]);
  }
}

void RegressionConfig::InteractionAllInput::next() {
  if (current_block < blocks.size()) {
    if (current_block > 0) covar_info[current_block-1]->unlink_data(true);
    current_block++; set_iterator();
  } 
}

pair<int,int> RegressionConfig::InteractionAllInput::var_count(bool main_only) {
  int nvar = interaction_count + (!main_only ? no_cond : 0) + 1;
  return pair<int,int>(nvar,nvar);
}

void RegressionConfig::InteractionEachInput::sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks) {
  set<int> &primary = config.variables, &secondary = config.variables_auxiliary;

  vector<short> shared(max(*primary.rbegin(), *secondary.rbegin())+1, 0);
  for (set<int>::iterator sec = secondary.begin(); sec != secondary.end(); ++sec) {if (primary.find(*sec) != primary.end()) shared[*sec] = true;}

  for (set<int>::iterator sec = secondary.begin(); sec != secondary.end(); ++sec) {
    VariableInfo& inter = config.variable_data.get_variable(*sec, true);  
    if (inter.is_valid()) {
      bool is_shared = shared[*sec];
      for (set<int>::iterator pri = primary.begin(); pri != primary.end(); ++pri) {
        if (*pri != *sec && !(is_shared && shared[*pri] && *pri < *sec)) {
          VariableInfo& main = config.variable_data.get_variable(*pri, true);
          if (main.is_valid()) {
            if (main.is_set() && inter.is_set()) ss_blocks.push_back(PairSS(main.as_set(), inter.as_set()));
            else if (main.is_set() && inter.is_covar()) sc_blocks.push_back(PairSC(main.as_set(), inter.as_covar()));
            else if (main.is_covar() && inter.is_set()) sc_blocks.push_back(PairSC(inter.as_set(), main.as_covar()));
          }
        }
      }
    }
  }
}

void RegressionConfig::InteractionPairsInput::sort_blocks(vector<PairSS>& ss_blocks, vector<PairSC>& sc_blocks) {
  set<int>& vars = config.variables; 
  for (set<int>::iterator curr = vars.begin(); curr != vars.end(); ++curr) {
    VariableInfo &var1 = config.variable_data.get_variable(*curr, true);
    if (var1.is_valid()) {
      for (set<int>::iterator other = curr; other != vars.end(); ++other) {      
        VariableInfo &var2 = config.variable_data.get_variable(*other, true);
        if (*curr < *other && var2.is_valid()) {
          if (var1.is_set() && var2.is_set()) ss_blocks.push_back(PairSS(var1.as_set(), var2.as_set()));
          else if (var1.is_set() && var2.is_covar()) sc_blocks.push_back(PairSC(var1.as_set(), var2.as_covar()));
          else if (var1.is_covar() && var2.is_set()) sc_blocks.push_back(PairSC(var2.as_set(), var1.as_covar()));
        }
      }
    }
  }
}



/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "setoutput.h"

const short SetResultsData::CoreBlock = 2;
const short SetResultsData::RegressionBlock = 3;

void SetOutput::clear_output() {
  Utils::delete_vector(output);
  Utils::delete_vector(abbr_vars);
}

void SetOutput::abbreviate_names() {
  Utils::delete_vector(abbr_vars);
  if (abbr_length > 0) {
    StringVariable& names = results.var_name; set<string> abbreviate;
    for (int i = 0; i < results.data_total(); i++) {
      if (!results.active_id(i)) continue;
      if (names.get(i).size() >= abbr_length) abbreviate.insert(names.get(i));
    }

    if (!abbreviate.empty()) {    
      map<string,int> abbr_index; vector<string> abbreviations; 
      for (set<string>::iterator it = abbreviate.begin(); it != abbreviate.end(); ++it) {
        string abbr = it->substr(0, abbr_length-2) + "...";
        abbr_index[*it] = abbreviations.size();
        abbreviations.push_back(abbr);
      }

      int dup_count = 0, last = abbreviations.size() - 1;      
      for (int i = 0; i < abbreviations.size(); i++) {
        bool match = i < last && abbreviations[i] == abbreviations[i+1];
        if (match || dup_count > 0) {
          abbreviations[i] += Utils::num_string(++dup_count);
          if (!match) dup_count = 0;
        }
      }
      
      if (abbr_file) {int row = 0;
        StringVariable* full_name = new StringVariable(abbreviations.size()); abbr_vars.push_back(full_name);
        StringVariable* abbr_name = new StringVariable(abbreviations.size()); abbr_vars.push_back(abbr_name);
        for (map<string,int>::iterator it = abbr_index.begin(); it != abbr_index.end(); ++it, ++row) {
          full_name->add(row, it->first); 
          abbr_name->add(row, abbreviations[it->second]);
        }

        for (int i = 0; i < results.data_total(); i++) {
          if (!results.active_id(i)) continue;
          if (names.get(i).size() >= abbr_length) names.set(i, abbreviations[abbr_index[names.get(i)]]);
        }
      } else {
        StringVariable* full_name = new StringVariable(); abbr_vars.push_back(full_name);
        results.align_variable(*full_name, true);
        for (int i = 0; i < results.data_total(); i++) {
          if (!results.active_id(i)) continue;
          string& curr_name = names.get(i);
          full_name->set(i, curr_name);        
          if (curr_name.size() >= abbr_length) curr_name = abbreviations[abbr_index[curr_name]];
        }
        output.push_back(full_name->get_output("FULL_NAME"));
      }
    }
  }  
}

void SetOutput::write_abbreviations(string filename, bool indent) {
  if (abbr_vars.size() >= 2) {
    _LOG << (indent ? "\twriting" : "Writing") << " variable abbreviations to file " << filename << endl; 
    FormattedOutput about(filename);

    int offset = output.size();
    output.push_back(abbr_vars[0]->get_output("FULL_NAME")); output.back()->set_header(about);    
    output.push_back(abbr_vars[1]->get_output("ABBREVIATION")); output.back()->set_header(about);
  
    about.print_header();
    while (!output[offset]->is_empty()) {for (int i = offset; i < output.size(); i++) about << *output[i];}
  }  
}

void SetOutput::make_columns(OutputMode mode) {
  Utils::delete_vector(output);
  if (mode == SelfContained) {
    output.push_back(results.var_name.get_output("GENE_SET"));
    output.push_back(results.ngenes.get_output("NGENES"));
    output.push_back(results.coefficient.get_output("MU"));
    output.push_back(results.pval.get_output("P"));
  } else if (mode == Regression) {
    output.push_back(results.var_name.get_output("VARIABLE"));
    output.push_back(results.var_type.get_output("TYPE"));  
    output.push_back(results.model_id.get_output("MODEL"));
    output.push_back(results.model_term.get_output_right("TERM"));
    output.push_back(results.ngenes.get_output("NGENES"));
  
    output.push_back(results.coefficient.get_output("BETA"));
    output.push_back(results.coefficient_std.get_output("BETA_STD"));
    output.push_back(results.std_error.get_output("SE"));
    output.push_back(results.pval.get_output("P"));
    if (results.warnings.initialised() && !results.warnings.all_missing()) output.push_back(results.warnings.get_output("WARNINGS"));  
  }
}

void SetOutput::write(OutputMode mode, string suffix, bool indent) {
  string filename = win_txt(prefix + "." + suffix + ".out");
  _LOG << (indent ? "\twriting" : "Writing") << " results to file " << filename << endl; 

  make_columns(mode);
  abbreviate_names();

  FormattedOutput fout(filename);

  for (vector<OutputColumn*>::iterator it = output.begin(); it != output.end();) {if (*it == 0 || (*it)->is_empty()) {delete *it; it = output.erase(it);} else ++it;}
  for (int i = 0; i < output.size(); i++) output[i]->set_header(fout);
  
  for (int i = 0; i < parameters.size(); i++) fout.print_param(parameters[i].first, parameters[i].second);
  fout.print_header();

  while (!output[0]->is_empty()) {for (int i = 0; i < output.size(); i++) fout << *output[i];}
  
  if (abbr_file) write_abbreviations(win_txt(prefix + "." + suffix + ".abbr"), indent);
  clear_output();
}


void SetGenesOutput::clear_output() {
  for (int i = 0; i < output.size(); i++) delete output[i];
  output.clear();
}

void SetGenesOutput::write_core(const string& name, int set_size, int ngenes, double pval, int model_id, bool interaction) {
  output_stream.clear_fields();

  if (current_set == 0) {
    output_stream.open(filename);
    output_stream.print_param("ALPHA", alpha);
    output_stream.print_param("NUMBER_OF_TESTS", no_tests);    
  }
  current_set++;  
  string prefix = string("_") + (interaction ? "INTER" : "SET") + Utils::num_string(current_set) + "_  ";

  for (deque<OutputColumn*>::iterator it = output.begin(); it != output.end();) {if (*it == 0 || (*it)->is_empty()) {delete *it; it = output.erase(it);} else ++it;}
  output.push_front(new ConstantOutput<string>(prefix, set_size));
  for (int i = 0; i < output.size(); i++) output[i]->set_header(output_stream);

  output_stream.skip_line(1);
  output_stream.print_param(prefix + "VARIABLE", name); 
  if (model_id > 0) output_stream.print_param(prefix + "MODEL", model_id); 
  output_stream.print_param(prefix + "NGENES", ngenes); 
  output_stream.print_param(prefix + "P-VALUE", pval);     

  output_stream.print_header();
  while (!output[0]->is_empty()) {for (int i = 0; i < output.size(); i++) output_stream << *output[i];}
  clear_output();
}


void SetGenesOutput::write_set(SetStatsUtils::VariableInfo& var, double pval, int model_id) {
  int masked = 0, ngenes = 0; 
  SetStatsUtils::CovarInfo* inter_sc = 0; SetStatsUtils::SetInfo* inter_ss = 0;
  if (var.is_interaction()) {
    if (var.has_type(SetStatsUtils::vc_InteractionSS)) {
      inter_ss = var.as_set(); SetStatsUtils::SetSetPair& sets = *inter_ss->get_interaction();
      set<int> joint = *sets.set1->get_data(), &other = *sets.set2->get_data();
      joint.insert(other.begin(), other.end());
      masked = gene_info.set_mask(joint); ngenes = inter_ss->get_ngenes();
    } else {
      inter_sc = var.as_covar();
      masked = gene_info.set_mask(*inter_sc->get_interaction()->set->get_data());
    }
  } else if (var.is_set()) masked = gene_info.set_mask(*var.as_set()->get_data());
  else return;
  if (!masked) return;
  if (ngenes == 0) ngenes = masked;
  
  output.push_back(gene_info.name.get_output("GENE"));
  output.push_back(gene_info.location.get_output());
  output.push_back(gene_info.nsnps.get_output("NSNPS"));
  output.push_back(gene_info.nparam.get_output("NPARAM"));  
  output.push_back(gene_info.nsamp.get_output("N"));

  if (inter_ss || inter_sc) {
    if (inter_ss) {
      SetStatsUtils::SetSetPair& sets = *inter_ss->get_interaction();
      output.push_back(sets.set1->get_storage()->data.get_output(sets.set1->get_column(), "SET1"));
      output.push_back(sets.set2->get_storage()->data.get_output(sets.set2->get_column(), "SET2"));
    } else if (inter_sc) {
      output.push_back(inter_sc->get_storage()->data.get_output(inter_sc->get_column(), "COVAR"));
    }  
  }

  typedef BlockMultiVariable<ZscoreVariable<double> > ZData;
  ZData* zstat = gene_info.typed_variable<ZData>("ZSTAT_OUTCOME");  
  if (zstat) output.push_back(zstat->get_output(0, "ZSTAT"));
  output.push_back(gene_info.pval.get_output(0, "P"));

  ZData* fit = gene_info.typed_variable<ZData>("ZSTAT_FITTED_BASE");
  if (fit) {
    ColumnBlock* block = fit->get_output();
    if (block) {
      block->add_column(0, "ZFITTED_BASE");          
      block->add_column(1, "ZRESID_BASE");
      output.push_back(block);
    }
  }
 
  write_core(var.get_name_label(), masked, ngenes, pval, model_id, var.is_interaction());
  gene_info.unset_mask();
}


void SetGenesOutput::print_log(bool indent) {
  if (current_set > 0) {
    _LOG << (indent ? "\twriting" : "Writing") <<  " gene analysis results per significant result (after multiple testing correction, at alpha = " << alpha << ") to file " << filename << endl;
  }
}



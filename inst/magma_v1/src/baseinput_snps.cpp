/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "baseinput.h"
#include "input.h"

                                                            
using namespace std;

void BaseInput::filter_snps(const string& filename, const string& mode) { 
  _LOG << "Filtering SNPs listed in file " << filename << " (mode = " << mode << ")" << endl;
  int value_index = filename.size() > 4 && (filename.substr(filename.size()-4) == ".bim") ? 1 : 0;
  
  TextInput fin(filename, value_index+1);
  fin.set_error("reading SNP filter file");

  vector<short> index(snp_info.data_total(), 0); int id;
  while (fin.process_line()) {if ((id = snp_info[fin[value_index]]) >= 0) index[id] = 1;}
  snp_info.filter_ids(index, mode == "include");

  if (snp_info.data_used(false) == 0) fin.error("no SNPs remaining after filtering");
  _LOG << "\tafter filtering, " << Utils::plural(snp_info.data_used(false), "SNP") << " remaining in data" << endl;
}

void BaseInput::load_synonyms(const string& filename, bool auto_detect) { 
  _LOG << "Reading SNP synonyms from file " << filename << (auto_detect ? " (auto-detected)" : "") << endl;

  TextInput fin(filename, 0);
  fin.set_error("reading SNP synonym file"); fin.skip_empty();
  string prefix; bool add_prefix = fin.get_param("PREFIX", prefix); 

  SNPSynonyms synonyms(settings.gets("snp_synonyms_dup", "error"));  
  synonyms.init(fin.file_size() / (snp_info.id_size()+1), snp_info.id_size());
  NumericID<int>& snpid = synonyms.snpid; 

  while (fin.process_line()) {
    while (fin.read_value()) {
      if (add_prefix) fin.value = prefix + fin.value;
      int iid = synonyms.get_iid(fin.value);
      snpid.set(iid, snp_info[fin.value]);
    }
    synonyms.process_synonym();
  }
  vector<vector<char*> > duplicates = synonyms.process(snp_info);
  int alias_count = synonyms.get_count("aliases"), target_count = synonyms.get_count("targets");

  _LOG << "\tread " << Utils::plural(alias_count, "mapped synonym") << " from file, mapping to " << Utils::plural(target_count, "SNP") << " in the data" << endl;
  if (!duplicates.empty()) {int dup_count = duplicates.size(); bool is_error = synonyms.duplicate_mode == SNPSynonyms::DupError;
    string obj, skip = Utils::space(is_error ? "ERROR - " : "WARNING: ");
    _SLOG.set_block("SNP filtering") << "# Following sets of synonymous SNP IDs were found in file " << filename << endl;  
    if (!is_error) {
      _LOG << "\tWARNING: detected " << Utils::plural(dup_count, "synonymous SNP pair") << " in the data" << endl << "\t" << skip; _SLOG.set_block("SNP synonym file") << "# Following sets of synonymous SNP IDs were found in file " << filename << endl;
      if (synonyms.duplicate_mode == SNPSynonyms::DupDropAll) {_LOG << "removed all involved synonymous SNPs in data from analysis"; _SLOG << "# All were removed from analysis";}
      else if (synonyms.duplicate_mode == SNPSynonyms::DupDropDup) {_LOG << "mapped all IDs for each synonym entry to first ID occurring in data, including other IDs present in data"; _SLOG << "# All IDs per line mapped to first ID, all but first SNP was removed from data";}   
      else if (synonyms.duplicate_mode == SNPSynonyms::DupSkipLine) {_LOG << "skipped all synonym entries involved, synonymous SNPs are kept in analysis"; _SLOG << "# All were retained in analysis, the corresponding synonym definitions were skipped";}   
      else if (synonyms.duplicate_mode == SNPSynonyms::DupSkipDup) {_LOG << "mapped all IDs for each synonym entry to first ID occurring in data, except for other IDs present in data"; _SLOG << "# All were retained in analysis, other IDs in same synonym definition were mapped to first listed ID";}   
      _LOG << endl << "\t" << skip << "writing list of detected synonyms in data to supplementary log file" << endl; _SLOG << endl; 
    } else _LOG.error("reading SNP synonym file") << "detected " << Utils::plural(dup_count, "synonymous SNP pair") << " in the data" << endl << skip << "writing list of detected synonyms in data to supplementary log file" << endl;

    for (int i = 0; i < duplicates.size(); i++) {vector<char*>& curr_dup = duplicates[i];
      for (int j = 0; j < curr_dup.size(); j++) _SLOG << (j ? "\t" : "") << curr_dup[j];
      _SLOG << endl;
    }
    if (is_error) die(); 
  } 
}

void BaseInput::load_pvalfile() { 
  string filename = settings.gets("pval_file");
  _LOG << "Reading SNP p-values from file " << filename << "... " << endl;
  
  TextInput fin(filename);
  fin.set_error("reading p-value file");
  long tot_var = fin.read_header();
  string* rsid = fin.set_var("snp", settings.gets("pval_col_snp", fin.has_header ? "SNP" : "1"));
  string* pval_raw = fin.set_var("p", settings.gets("pval_col_p", fin.has_header ? "P" : "2"));
  
  _LOG << "\tdetected " << tot_var << " variables in file" << endl;
  _LOG << "\tusing variable: " << fin.name("snp") << " (SNP id)" << endl;
  _LOG << "\tusing variable: " << fin.name("p") << " (p-value)" << endl;

  int min_n = settings.geti("min_sample_size");
  string* n_col = has_snpN ? fin.set_var("N", settings.gets("pval_col_N")) : 0;
  if (has_snpN) _LOG << "\tusing variable: " << fin.name("N") << " (sample size; discarding SNPs with N < " << min_n << ")" << endl;

  enum DupMode {DupError, DupDrop, DupFirst, DupLast} duplicate; string dubstr = settings.gets("pval_duplicate_mode", "error");
  if (dubstr == "drop") duplicate = DupDrop; else if (dubstr == "first") duplicate = DupFirst; else if (dubstr == "last") duplicate = DupLast; else duplicate = DupError;  
  bool has_synonyms = settings["snp_synonyms"]; 

  PvalueVariable<double>& pval = snp_info.pval; pval.init();
  pval.set_bounds(settings.getn("pval_truncate_low", 0), 1 - settings.getn("pval_truncate_high", 0)); 
  PositiveNumber<int>& snp_n = snp_info.nsamp; snp_n.set_bound(min_n);
  if (has_snpN) {snp_n.init(); snp_n.set_rounding(true);}
  int id; set<string> duplicates; map<int,set<char*> > dup_synonyms;

  while (fin.process_line()) {
    if ((id = snp_info[*rsid]) < 0) continue;

    if (fin.read_value()) fin.line_error("too many values");
    if (pval.get(id) != -1) {
      duplicates.insert(*rsid);
      if (duplicate == DupError) fin.line_error("duplicate SNP " + string((has_synonyms ? "(same ID or synonym)" : "ID")));        
      if (duplicate == DupFirst) continue;
      if (duplicate == DupDrop) {snp_info.pval.get(id) = -2; continue;}
    } 
    
    pval.add(id, *pval_raw);
    if (pval.get_status() == InputVariable::Missing) {snp_info.drop_id(id); continue;}
    if (!pval.is_valid()) fin.line_error(pval.get_msg() + " for variable " + fin.name("p"));

    if (has_snpN) {
      snp_n.add(id, *n_col);
      if (snp_n.get_status() == InputVariable::Missing || snp_n.get_status() == InputVariable::OutOfRange) {snp_info.drop_id(id); continue;}
      if (!snp_n.is_valid()) fin.line_error(snp_n.get_msg() + " for variable " + fin.name("N"));
    }
  }
  if (!duplicates.empty() && has_synonyms) dup_synonyms = snp_info.get_synonyms(duplicates);
  pval.filter_range(0, VariableStorageAccess::Smaller);

  int total = snp_info.data_used(false);
  if (total == 0) fin.error("no valid p-values found for any SNPs in the data");
  _LOG << "\tread " << fin.line_no << " lines from file, containing valid SNP p-values for " << total << " SNPs in data (";
  _LOG << setprecision(4) << (float(100*total) / fin.line_no) << "% of lines, ";
  _LOG << setprecision(4) << (float(100*total) / snp_info.data_total()) << "% of SNPs in data)" << endl;

  if (!duplicates.empty()) {int dup_count = has_synonyms ? dup_synonyms.size() : duplicates.size();
    string obj, plural = (dup_count != 1 ? "s" : ""), skip = Utils::space("WARNING: ");
    if (has_synonyms) obj = "SNP" + plural + " (same ID" + plural + " or synonym" + plural + ")"; else obj = "SNP ID" + plural;
    _LOG << "\tWARNING: file contained " << dup_count << " " << obj << " " << "with duplications" << endl << "\t" << skip; _SLOG.set_block("SNP p-value file") << "# Following SNPs had duplications entries in file " << filename << endl; 
    if (has_synonyms) _SLOG << "# Duplicates may have been synonyms; listing all synonymous IDs in data for SNP on each line" << endl;
    if (duplicate == DupFirst) {_LOG << "using first occurrence for each"; _SLOG << "# First occurrence in file for each was used";}
    else if (duplicate == DupLast) {_LOG << "using last occurrence for each"; _SLOG << "# Last occurrence in file for each was used";}
    else if (duplicate == DupDrop) {_LOG << "dropped all occurrences of each from analysis"; _SLOG << "# All occurrences were dropped from analysis";}       
    _LOG << endl << "\t" << skip << "writing list of duplicated IDs to supplementary log file" << endl; _SLOG << endl; 
    
    if (has_synonyms) {
      for (map<int,set<char*> >::iterator it = dup_synonyms.begin(); it != dup_synonyms.end(); ++it) {set<char*>& syn = it->second; 
        for (set<char*>::iterator jt = syn.begin(); jt != syn.end(); ++jt) {if (jt != syn.begin()) _SLOG << "\t"; _SLOG << *jt;}
        _SLOG << endl;
      }
    } else for (set<string>::iterator it = duplicates.begin(); it != duplicates.end(); ++it) _SLOG << *it << endl;
  }
}

void BaseInput::load_burden_weights() {
  string filename = settings.gets("rare_weights_file");
  _LOG << "Reading SNP weights from file " << filename << "... " << endl;
  
  TextInput fin(filename);
  fin.set_error("reading weights file");
  long tot_var = fin.read_header();
  string* rsid = fin.set_var("snp", settings.gets("rare_weights_col_snp", fin.has_header ? "SNP" : "1"));
  string* weight_raw = fin.set_var("weight", settings.gets("rare_weights_col_weight", fin.has_header ? "WEIGHT" : "2"));
  
  _LOG << "\tdetected " << tot_var << " variables in file" << endl;
  _LOG << "\tusing variable: " << fin.name("snp") << " (SNP id)" << endl;
  _LOG << "\tusing variable: " << fin.name("weight") << " (weight)" << endl;

  enum DupMode {DupError, DupDrop, DupFirst, DupLast} duplicate; string dubstr = settings.gets("rare_weights_duplicate_mode", "error");
  if (dubstr == "drop") duplicate = DupDrop; else if (dubstr == "first") duplicate = DupFirst; else if (dubstr == "last") duplicate = DupLast; else duplicate = DupError;  
  bool has_synonyms = settings["snp_synonyms"]; 
  float default_weight = settings.getn("rare_weights_default");

  PositiveWeightVariable<float>& weights = snp_info.weights; weights.init();
  int id, valid_read = 0, typed_read = 0; set<string> duplicates; map<int,set<char*> > dup_synonyms;

  map<string,pair<float,int>> type_weights; set<string> missed_types; bool use_types = settings["rare_weights_type_file"];
  enum TypeMode {TypeError, TypeDrop, TypeDefault} missing_type; string typestr = settings.gets("rare_weights_type_missing", "drop");
  if (typestr == "drop") missing_type = TypeDrop; else if (typestr == "default") missing_type = TypeDefault; else missing_type = TypeError;
  if (use_types) {
    _LOG << "\treading type weights from file " << settings.gets("rare_weights_type_file") << endl;
    TextInput wt_file(settings.gets("rare_weights_type_file"), 2); float wt;
    wt_file.set_error("reading weights file"); wt_file.skip_empty();

    while (wt_file.process_line()) {
      string type = Utils::uppercase(wt_file[0]);
      if (type_weights.find(type) != type_weights.end()) wt_file.line_error(string("duplicate type entry '") + wt_file[0] + "'");
      if (!wt_file.read_num(wt, 1) || wt < 0) wt_file.line_error(string("value '") + wt_file[1] + "' is not a valid weight value");
      type_weights[type] = pair<float,int>(wt,0);
    }
  }

  while (fin.process_line()) {
    if ((id = snp_info[*rsid]) < 0) continue;

    if (fin.read_value()) fin.line_error("too many values");
    if (weights.get(id) != -1) {
      duplicates.insert(*rsid);
      if (duplicate == DupError) fin.line_error("duplicate SNP " + string((has_synonyms ? "(same ID or synonym)" : "ID")));        
      if (duplicate == DupFirst) continue;
      if (duplicate == DupDrop) {snp_info.weights.get(id) = -2; continue;}
    } 

    if (use_types) {
      string type = Utils::uppercase(*weight_raw); map<string,pair<float,int>>::iterator found = type_weights.find(type);
      if (found == type_weights.end()) {
        if (missing_type == TypeError) fin.line_error(string("unknown type '") + *weight_raw + "'");
        if (missed_types.find(*weight_raw) == missed_types.end()) {
          if (missing_type == TypeDrop) _LOG << "WARNING: unknown type '" << *weight_raw <<"'; skipping this and future instances" << endl;
          if (missing_type == TypeDefault) _LOG << "WARNING: unknown type '" << *weight_raw <<"'; setting weight to default weight (" << default_weight << ") for this and future instances" << endl;
          missed_types.insert(*weight_raw);
        }
        if (missing_type == TypeDrop) {weights.set(id, -2); continue;}
        if (missing_type == TypeDefault) weights.set(id, default_weight);
      } else {
        weights.set(id, found->second.first);
        found->second.second++;
        typed_read++;
      }
    } else {
      weights.add(id, *weight_raw);
      if (weights.get_status() == InputVariable::Missing) {snp_info.drop_id(id); continue;}
      if (!weights.is_valid()) fin.line_error(weights.get_msg() + " for variable " + fin.name("weight"));
    }
    valid_read++;
  }
  if (!duplicates.empty() && has_synonyms) dup_synonyms = snp_info.get_synonyms(duplicates);
  if (default_weight > 0) {
    for (int i = 0; i < snp_info.data_total(); i++) {if (snp_info.active_id(i) && weights.get(i) == -1) weights.set(i, default_weight);}
  }
  weights.filter_range(0, VariableStorageAccess::SmallerEqual);

  int total = snp_info.data_used(false);
  if (valid_read == 0) fin.error("no valid weights found for any SNPs in the data");
  if (use_types && typed_read == 0) fin.error("none of SNPs had type specified in type weight file");
  _LOG << "\tread " << fin.line_no << " lines from file, containing SNP weights for " << valid_read << " SNPs in data" << endl;
  _LOG << "\tafter setting defaults (default weight = " << default_weight << ") and filtering, " << total << " SNPs remaining for analysis" << endl;

  vector<string> unused;
  for (map<string,pair<float,int>>::iterator curr = type_weights.begin(); curr != type_weights.end(); ++curr) {if (curr->second.second == 0) unused.push_back(curr->first);}
  if (!unused.empty()) _LOG << "WARNING: following types were never used: " << Utils::join_string(unused, ", ") << endl;


  if (!duplicates.empty()) {int dup_count = has_synonyms ? dup_synonyms.size() : duplicates.size();
    string obj, plural = (dup_count != 1 ? "s" : ""), skip = Utils::space("WARNING: ");
    if (has_synonyms) obj = "SNP" + plural + " (same ID" + plural + " or synonym" + plural + ")"; else obj = "SNP ID" + plural;
    _LOG << "\tWARNING: file contained " << dup_count << " " << obj << " " << "with duplications" << endl << "\t" << skip; _SLOG.set_block("SNP weight file") << "# Following SNPs had duplications entries in file " << filename << endl;
    if (has_synonyms) _SLOG << "# Duplicates may have been synonyms; listing all synonymous IDs in data for SNP on each line" << endl;
    if (duplicate == DupFirst) {_LOG << "using first occurrence for each"; _SLOG << "# First occurrence in file for each was used";}
    else if (duplicate == DupLast) {_LOG << "using last occurrence for each"; _SLOG << "# Last occurrence in file for each was used";}
    else if (duplicate == DupDrop) {_LOG << "dropped all occurrences of each from analysis"; _SLOG << "# All occurrences were dropped from analysis";}       
    _LOG << endl << "\t" << skip << "writing list of duplicated IDs to supplementary log file" << endl; _SLOG << endl; 
    
    if (has_synonyms) {
      for (map<int,set<char*> >::iterator it = dup_synonyms.begin(); it != dup_synonyms.end(); ++it) {set<char*>& syn = it->second; 
        for (set<char*>::iterator jt = syn.begin(); jt != syn.end(); ++jt) {if (jt != syn.begin()) _SLOG << "\t"; _SLOG << *jt;}
        _SLOG << endl;
      }
    } else for (set<string>::iterator it = duplicates.begin(); it != duplicates.end(); ++it) _SLOG << *it << endl;
  }
}


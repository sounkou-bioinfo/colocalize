/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "setinput.h"
#include "input.h"

#define MIN_GENE_COUNT 3

SetVariables* SetInput::load_data() {
  SetVariables* vars = new SetVariables(gene_info, settings.getn("gc_correction"), settings.getn("genelevel_prune_pc"));

  if (settings["gene_filter"]) filter_genes(settings.gets("gene_filter_file", ""), settings.gets("gene_filter"));

  int valid = 0;
  if (settings["has_set_annot"] || settings["has_gene_covar"]) {
    if (settings["has_set_annot"]) {
      _LOG << "Loading gene-set annotation..." << endl;
      GeneSetData* sets = load_setfile("SETS_MAIN", settings.gets("set_annot_file"), settings["has_set_annot_col"]);
      vars->add_sets(sets); valid += sets->count_valid();
    }

    if (settings["has_gene_covar"]) {
      _LOG << "Loading gene-level covariates..." << endl;
      GeneCovarData<float>* covar = load_covarfile("COVAR_MAIN", settings.gets("gene_covar_file"));
      vars->add_covar(covar); valid += covar->count_valid();
    }
  } else {_LOG.error("loading gene-level data") << "no input data specified" << endl; die();}
  if (valid == 0) {_LOG.error("loading gene-level data") << "no valid gene sets or covariates found in input files" << endl; die();}

  bool duplicates = !vars->lock_variables();
  if (duplicates) {_LOG.error("loading gene-level data") << "duplicate variable IDs encountered in at least one input file" << endl; die();}
 
  process_missing(vars); /// data is now complete and consecutive, and presumed in subsequent code not to change size, order or filtering
  prep_variables(vars);

  return vars;
}

void SetInput::filter_genes(const string& filename, const string& mode) { 
  _LOG << "Filtering genes listed in file " << filename << " (mode = " << mode << ")" << endl;
  
  TextInput fin(filename, 1);
  fin.set_error("reading gene filter file");

  vector<short> index(gene_info.data_total(), 0); int id;
  while (fin.process_line()) {if ((id = gene_info[fin[0]]) >= 0) index[id] = 1;}
  gene_info.filter_ids(index, mode == "include");

  if (gene_info.data_used(false) == 0) fin.error("no genes remaining after filtering");
  _LOG << "\tafter filtering, " << Utils::plural(gene_info.data_used(false), "gene")  << " remaining in data" << endl;
}

// for now, ZSTAT can be assumed to be complete from .raw file load checks
void SetInput::process_missing(SetVariables* var_data) {
  if (var_data->no_covar_blocks()) {   
    string impute_type = settings.gets("gene_covar_impute_type"); 
    int impute_value = settings.getn("gene_covar_impute_value", 0);
    bool drop_genes = impute_type == "drop" || !settings["gene_covar_fill_rows"];
      
    pair<int,int> missing = var_data->process_missing(drop_genes, impute_type, impute_value);
    if (missing != pair<int,int>(0,0)) {
      _LOG << "Processing missing values..." << endl;
  
      if (missing.first > 0) {
        _LOG << "\tfound " << Utils::plural(missing.first, "gene") << " not present in all input files: "; 
        if (drop_genes) _LOG << "removing these from analysis" << endl;
        else _LOG << "treating these as missing on corresponding variables" << endl;
      }
      if (missing.second > 0) {
        string type = missing.first > 0 ? "additional gene" : "gene";
        _LOG << "\tfound " << Utils::plural(missing.second, type) << " with partially missing values";
        if (impute_type == "drop") _LOG << ": removing these from analysis";
        _LOG << endl;
      }
      if ((missing.first > 0 && !drop_genes) || (missing.second > 0 && impute_type != "drop")) {
        _LOG << "\tsetting missing values to ";
        if (impute_type == "value") _LOG << "value '" << impute_value << "'" << endl;
        else _LOG << impute_type << " of corresponding variable" << endl;
      }
      
      if ((missing.first > 0 && drop_genes) || (missing.second > 0 && impute_type == "drop")) _LOG << "\t" << Utils::plural(gene_info.data_used(false), "gene") << " remaining in analysis" << endl;
    
      int req_obs = ceil(0.5*gene_info.data_used(false));
      for (int c = 0; c < var_data->no_covar_blocks(); c++) {
        GeneCovarData<float>& vars = var_data->get_covar(c);
        if (vars.data_file != "") {
          for (int v = 0; v < vars.data_total(); v++) {
            if (vars.is_valid(v)) {
              if (vars.values_observed.get(v) < req_obs) _LOG << "\tWARNING: discarding variable '" << vars.name.get(v) << "' from file " << vars.data_file << ", more than half of genes to be used in analysis had missing values" << endl; 
              else if (vars.values_observed.get(v) < MIN_GENE_COUNT) _LOG << "\tWARNING: discarding variable '" << vars.name.get(v) << "' from file " << vars.data_file << ", fewer than " << MIN_GENE_COUNT << " of genes to be used in analysis had non-missing values" << endl;             
              else continue;
              vars.status.set(v, VariableIndex::HighMissing);
            }
          }
        }
      }
    }
  }
}

void SetInput::prep_variables(SetVariables* var_data) {
  _LOG << "Preparing variables for analysis..." << endl;
  if (gene_info.data_used(false) < MIN_GENE_COUNT) {_LOG.error("processing input") << "less than " << MIN_GENE_COUNT << " genes available for analysis" << endl; die();}

  double min_sd = settings.getn("gene_vars_min_sd"); 
  double zlow = settings.getn("gene_zstat_truncate_low"), zhigh = settings.getn("gene_zstat_truncate_high");
  var_data->prep_outcome(zlow, zhigh);
  GeneCovarData<double>& zstat = var_data->get_outcome();
  if (zstat.is_valid(0)) {
    if (zstat.orig_sd.get(0) < min_sd) {
      _LOG.error("processing input") << "variance of the gene Z-statistic variable too low" << endl;
      zstat.status.set(0, VariableIndex::LowVariance); die();              
    }
  } else {
    if (zstat.status.get(0) == VariableIndex::Monotonic) {_LOG.error("processing input") << "the gene Z-statistic variable is (almost) constant" << endl; die();}   
    else {_LOG.error("processing input") << "the gene Z-statistic variable is invalid (unknown error)" << endl; die();}
  }
  _LOG << "\ttruncating Z-scores " << zlow << " points below zero or " << zhigh << " standard deviations above the mean" << endl;

  int total_genes = gene_info.data_used(false); bool allow_mono = settings["gene_sets_allow_mono"]; 
  if (var_data->no_set_blocks()) {
    var_data->prep_sets();
    for (int s = 0; s < var_data->no_set_blocks(); s++) {
      GeneSetData& sets = var_data->get_sets(s);
      for (int v = 0; v < sets.data_total(); v++) {
        if (sets.is_valid(v)) {
          int ngenes = sets.used_genes.get(v); string msg;

          if (ngenes == 0) msg = "no genes";
          else if (ngenes >= total_genes) msg = "all genes";
          else if (!allow_mono) {
            if (ngenes == 1) msg = "only one gene";
            else if (ngenes == total_genes - 1) msg = "all but one gene";
          }
          
          if (!msg.empty()) {
            _LOG << "\tWARNING: discarding gene set '" << sets.name.get(v) << "' from file " << sets.data_file << ", variance is too low (set contains " << msg << " used in analysis)" << endl;             
            sets.status.set(v, VariableIndex::LowVariance);
          }
        }      
      }
    }
  }
      
  if (var_data->no_covar_blocks()) {
    double trunc = settings.getn("gene_covar_truncate");
    _LOG << "\ttruncating covariate values more than " << trunc << " standard deviations from the mean" << endl;
    var_data->prep_covar(trunc);

    for (int c = 0; c < var_data->no_covar_blocks(); c++) {
      GeneCovarData<float>& vars = var_data->get_covar(c);
      bool is_internal = (c == var_data->internal_block);
      for (int v = 0; v < vars.data_total(); v++) {
        if (vars.is_valid(v)) {
          if (vars.orig_sd.get(v) < min_sd) {
            vars.status.set(v, VariableIndex::LowVariance);
            if (!is_internal && vars.data_file != "") _LOG << "\tWARNING: discarding variable '" << vars.name.get(v) << "' from file " << vars.data_file << ", variance is too low" << endl;             
          }
        } else if (vars.status.get(v) == VariableIndex::Monotonic) {
          if (!is_internal && vars.data_file != "") _LOG << "\tWARNING: discarding variable '" << vars.name.get(v) << "' from file " << vars.data_file << ", the variable is (almost) constant" << endl;                     
        }
        
        if (is_internal && !vars.is_valid(v)) {
          set<int> ids = var_data->variables.find_var(vars.name.get(v), SetStatsUtils::vc_Internal); 
          for (int v2 = 0; v2 < vars.data_total(); v2++) {
            if (v2 != v && vars.is_valid(v2) && ids.find(vars.variable_id.get(v2)) != ids.end()) vars.status.set(v2, vars.status.get(v));
          }
        }     
      }
    }
  }
  
  pair<int,int> no_sets = var_data->variables.no_sets(), no_covar = var_data->variables.no_covar(true);
  if (no_sets.first + no_covar.first > 0) {
    _LOG << "\ttotal variables available for analysis: ";
    if (no_sets.second > 0) _LOG << Utils::plural(no_sets.first, "gene set");
    if (no_covar.second > 0) _LOG << (no_sets.second > 0 ? " and " : "") << Utils::plural(no_covar.first, "gene covariate");
    _LOG << endl;
  } else {_LOG.error("processing input") << "no input variables to analyse" << endl; die();}
}

GeneSetData* SetInput::load_setfile(const string& dataid, const string& filename, bool by_column) {
  _LOG << "Reading file " << filename << "... " << endl;

  TextInput fin(filename, 1);
  fin.set_error("reading gene-set annotation file");  

  GeneSetData* set_info = new GeneSetData(dataid, filename); 
  set_info->init(1000); set_info->link_data(gene_info);
  SparseBooleanVariable& sets = set_info->data; 
  vector<short> genes_used(gene_info.data_total(), 0); int tot_sets = 0;
  if (by_column) {
    int tot_var = fin.read_header();
    if (Utils::is_int(settings.gets("set_annot_col_gene", "1")) && Utils::is_int(settings.gets("set_annot_col_set", "2"))) fin.drop_header();
    string* gid = fin.set_var("geneid", settings.gets("set_annot_col_gene", "1"));
    string* sid = fin.set_var("setid", settings.gets("set_annot_col_set", "2"));

    _LOG << "\tdetected " << tot_var << " variables in file" << endl;
    _LOG << "\tusing variable: " << fin.name("geneid") << " (gene ID)" << endl;
    _LOG << "\tusing variable: " << fin.name("setid") << " (set ID)" << endl;
  
    map<string,int> set_index; vector<set<int> > set_genes; 
    string curr_name = ""; int curr_index = 0;
    while (fin.process_line()) {
      if (*sid != curr_name) {curr_name = *sid;
        if (set_index.find(curr_name) == set_index.end()) {
          set_index[curr_name] = curr_index = set_genes.size(); 
          set_genes.push_back(set<int>());
        } else curr_index = set_index[curr_name];
      }
      
      int gene_id = gene_info[*gid];   
      if (gene_id >= 0) {
        set_genes[curr_index].insert(gene_id);
        genes_used[gene_id] = true;  
      }
    }  

    for (map<string,int>::iterator it = set_index.begin(); it != set_index.end(); ++it) {
      set<int>& genes = set_genes[it->second]; 

      int set_id = set_info->get_iid(it->first);
      set_info->mapped_genes.set(set_id, genes.size());    
      set_info->used_genes.set(set_id, genes.size());          

      if (genes.empty()) {
        _LOG << "\tWARNING: gene set " << it->first << " contains no genes defined in genotype data" << endl; 
        set_info->status.set(set_id, VariableIndex::LowVariance);
      }

      int col_id = sets.insert(genes);      
      set_info->storage_column.set(set_id, col_id);
    }
    tot_sets = set_index.size(); 
  } else {set<int> genes; 
    while (fin.process_line()) {
      if (!fin.read_value()) fin.line_error("too few values");
  
      genes.clear(); tot_sets++;
      do {
        int gene_id = gene_info[fin.value]; 
        if (gene_id >= 0) { 
          genes.insert(gene_id);
          genes_used[gene_id] = true;  
        }
      } while(fin.read_value());
         
      int set_id = set_info->get_iid(fin[0]);
      set_info->mapped_genes.set(set_id, genes.size());    
      set_info->used_genes.set(set_id, genes.size());                

      if (genes.empty()) {
        _LOG << "\tWARNING: gene set " << fin[0] << " contains no genes defined in genotype data" << endl; 
        set_info->status.set(set_id, VariableIndex::LowVariance);
      }

      int col_id = sets.insert(genes);
      set_info->storage_column.set(set_id, col_id);
    }
  } 
  
  set_info->build_index();   
  if (set_info->data_used(false) == 0) fin.error("found no gene sets containing genes defined in genotype data");

  if (!set_info->get_duplicates().empty()) {vector<char*>& duplicates = set_info->get_duplicates();
    _SLOG.set_block("Set annotation file") << "# Following gene-set names had duplications in file " << filename << endl;
    for (int i = 0; i < duplicates.size(); i++) _SLOG << duplicates[i] << endl;
    fin.error(string("file contained ") + Utils::plural(duplicates.size(), "duplicate gene-set names") + "; writing list of IDs to supplementary log file");    
  }
  
  _LOG << "\t" << Utils::plural(tot_sets, "gene-set definition") << " read from file" << endl;
  _LOG << "\tfound " << Utils::plural(set_info->count_valid(), "gene set") << " containing genes defined in genotype data (containing a total of " << MathUtils::count(genes_used) << " unique genes)" << endl;
  
  return set_info;
}


GeneCovarData<float>* SetInput::load_covarfile(const string& dataid, const string& filename) {
  _LOG << "Reading file " << filename << "... " << endl;

  TextInput fin(filename); 
  fin.set_error("reading gene covariate file");

  int col_offset = 1, tot_covar = fin.read_header(col_offset, false);
  fin.set_subset(settings.get_genecovar(true), settings.get_genecovar(false)); fin.set_map();
  
  vector<pair<int,string> > cov_names = fin.load_names();
  GeneCovarData<float>* covar_info = new GeneCovarData<float>(dataid, filename); covar_info->init(cov_names.size());
  for (int i = 0; i < cov_names.size(); i++) {
    int cov_id = covar_info->get_iid(cov_names[i].second);
    covar_info->file_column.set(cov_id, cov_names[i].first);
    covar_info->storage_column.set(cov_id, i);
  }  
  covar_info->build_index(); covar_info->link_data(gene_info);
  int no_covar = covar_info->data_total();

  if (no_covar == 0) fin.error("no variables selected");
  _LOG << "\tdetected " << tot_covar << " variables in file";
  if (no_covar == tot_covar) _LOG << " (using all)" << endl;
  else if (no_covar > 10) _LOG << " (using " << no_covar << ")" << endl;
  else _LOG << endl << "\tusing " << Utils::plural(no_covar, "variables", false) << ": " << fin.print_names() << endl;

  BlockMultiVariable<NumericVariable<float> >& covar = covar_info->data; BooleanVariable& in_file = covar_info->in_source; 
  NonNegativeNumber<int>& miss_count = covar_info->miss_count; miss_count.get_data().assign_value(no_covar);
  vector<int> non_missing(no_covar, 0); int total = 0;

  while (fin.process_line()) {int gene_id = gene_info[fin[0]]; 
    if (gene_id < 0) continue;
    if (in_file.get(gene_id)) fin.line_error("duplicate gene entry", string("ID = ") + fin[0]);
    in_file.set(gene_id, true); total++;

    int row_missing = 0;
    for (int i = 0; i < no_covar; i++) {
      covar.add(gene_id, i, fin(i));
      if (!covar.is_valid_or_missing(true)) fin.line_error(string("non-numeric value for variable ") + fin.name(i,true));
      else if (covar.is_valid(true)) non_missing[i]++;
      else row_missing++;
    }
    if (fin.read_value()) fin.line_error("too many values");
    miss_count.set(gene_id, row_missing);
  }

  float max_miss = settings.getn("gene_covar_max_miss", 0); 
  int min_obs = max(int(ceil((1-max_miss)*total)), 10);
  for (int i = 0; i < no_covar; i++) {
    covar_info->complete_column.set(i, non_missing[i] == total);
    covar_info->values_observed.set(i, non_missing[i]);
    if (non_missing[i] < min_obs) {
      covar_info->status.set(i, VariableIndex::HighMissing);
      _LOG << "\tWARNING: discarding variable '" << covar_info->name.get(i) << "', missing value rate exceeds missingness threshold (" << max_miss << ")" << endl;  
    }
  }

  _LOG << "\tfound " << Utils::plural(covar_info->count_valid(), "valid gene covariate") << ", for " << total << " genes defined in genotype data" << endl;
  return covar_info;
}


/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "parse.h"
#include <iostream>
#include <sstream>
#include <dirent.h>
#include "utils.h"
#include "input.h"

#ifndef VERNO
  #define VERNO 1.10
#endif

#ifndef VERTYPE
  #define VERTYPE custom
#endif

using namespace std;

const char* Settings::Value::type_names[6] = {"unknown", "integer", "floating point", "string", "list", "set"};

template<> long long Settings::typed_value<long long>(const Parameter& par) const {return get_value(par).geti();}
template<> double Settings::typed_value<double>(const Parameter& par) const {return get_value(par).getn();}  
template<> string Settings::typed_value<string>(const Parameter& par) const {return get_value(par).gets();}  
template<> vector<string> Settings::typed_value<vector<string> >(const Parameter& par) const {return get_value(par).getvs();}    
template<> vector<double> Settings::typed_value<vector<double> >(const Parameter& par) const {return get_value(par).getvn();}    

static string VERTYPE_HACK() {
  string str = MACRO_STRING(VERTYPE);  
  if (str == "lin") str = "linux";
  else if (str.size() >= 3 && str.substr(0,4) == "lin/") str = "linux/" + str.substr(4);
  return str;
}

int Settings::get_version() {return Utils::pad_num(Utils::num_char(MACRO_STRING(VERNO)), 100);}

void Settings::load_args() {
  if (!is_flag(raw_args[1])) parse_error("first argument is not a flag");
  for (int i = 1; i < no_args; i++) {
    if (is_flag(raw_args[i])) {
      arguments.push_back(Flag(string(raw_args[i] + 2)));
      if (flag_set(arguments.back().name())) parse_error("duplicate flag --" + arguments.back().name());
      used_flags.insert(arguments.back().name());
    } else {
      arguments.back().add_value(string(raw_args[i]));
    }
  }
  no_args = arguments.size();
}

void Settings::init_logfile() {
  for (int i = 0; i < no_args; i++) {
    if (arguments[i].name() == "out") {
      parse(arguments[i]);
      break;
    }
  }

  _LOG.set_logfile(gets("out_prefix") + ".log", true);
}

void Settings::start_output() {
  _LOG.set_logfile(gets("out_prefix").append(".log"), false);
  _SLOG.set_logfile(gets("out_prefix") + ".log.suppl");
}

string Settings::get_prefix() const {
  ostringstream out_pref;
  out_pref << gets("out_prefix_base");
  if ((*this)["batch"]) out_pref << ".batch" << geti("batch_index") << "_" << geti("batch");
  else if ((*this)["batch_chr"]) out_pref << ".batch" << Utils::chr_string(geti("batch_chr"),0) << "_chr";
  return out_pref.str();
}

void Settings::process_flags() {
  check_flag_consistency();
  set_order_rest();
  for (int i = 0; i < flag_order.size(); i++) {
    Flag& curr = arguments[flag_order[i]];
    if (curr.name("out")) continue; ///already parsed through init_logfile()
    parse(curr);
  }
  prologue();
  arguments.clear();
}

void Settings::set_order(const string& flag) {
  for (int i = 0; i < arguments.size(); i++) {
    if (flag == "" || arguments[i].name(flag)) {
      bool add = true;
      for (int j = 0; j < flag_order.size(); j++) if (flag_order[j] == i) {add = false; break;}
      if (add) flag_order.push_back(i);
    }
  }
}

void Settings::set_defaults() {
  set_value("VERSION", string(MACRO_STRING(VERNO)) + "/" + VERTYPE_HACK());
  set_value("VER_ID", get_version());

///IN/OUT
  set_marker("out_prefix");
  set_value("out_prefix_base", "magma");

///ANNOTATION
  set_value("window_up", 0);
  set_value("window_down", 0);

///PHENOTYPE / COVARIATES
  set_value("min_sample_size", 50);
  set_value("pheno_min_prop", 0.02);
  set_value("pheno_min_count", 25); 
  set_value("pheno_min_cap", 500);
  set_value("pheno_min_resvar", 0.01);   
  set_marker("use_residuals");
  set_marker("chrX_sex_covar"); 

///SNPS
  set_marker("geno_normalize");
  set_marker("rare_autoburden");
  set_value("rare_multi_max_count", 25);
  set_value("prune_cutoff", 0.95);
  set_value("prune_value", 0.1);
  set_value("prune_prop", 1);  
  set_value("prune_varimax_cutoff", 0.99);
  set_value("prune_varimax_value", 0.1);                     
  set_value("snp_max_miss", 0.05);
  set_value("snp_max_miss_max", 0.25);
  set_value("snp_min_mac", 0);
  set_value("snp_min_maf", 0);
  set_value("snp_diff_pval", 0.000001);
  set_value("snp_partition_max", 10);

///GENE ANALYSIS
  set_value("gene_model", "pcreg");
  set_value("gene_chunked_fraction", 0.5);
  set_value("gene_chunked_size", 1000);  
  set_value("corr_range", 5000000);
  set_value("corr_min_snps", 50);
  set_value("corr_max_step", 5);
  set_value("block_corr_min_snps", 50);
  set_value("block_corr_max_step", 5);
  set_value("min_perm", 1000);
  set_value("max_perm", 1000000);
  set_value("adap_count", 10);
  set_value("metagene_truncate_low", 5);
  set_value("metagene_truncate_high", 25);
  set_value("varimax_iter", 100);
  set_value("gene_top_corr_weight", 0.5);
  set_value("gene_top_corr_p1", 0.5);
  set_value("gene_top_corr_p2", 2);
  set_value("gene_top_corr_p3", 1.25);
  set_value("gene_multi_cperm", 250);
  set_value("gene_gpd_min_sim", 1000);  
  set_value("gene_gpd_nsim", 50000);
  set_value("gene_gpd_tail_min", 20000);  
  set_value("gene_gpd_adap_count", 10);  
  set_value("gene_gpd_prune", 0.999);  
  set_value("gene_gpd_tail", 0.005);    
  set_value("gene_gpd_tail_margin", 0.5);      
  set_value("gene_gpd_tail_fit", 0.025);
  set_value("gene_imhof_eigen_thresh", 1e-4);
  set_value("gene_imhof_epsilon_factor", 1e5);
  set_value("gene_imhof_perm_max", 1e9);
  set_value("gene_imhof_adap_thresh", 1000);  
  set_value("gene_imhof_pval_cutoff", 1e-25);
  set_value("gene_linked_outlier_thresh", 5); 
  set_value("gene_model_bigdata", "auto");   
  set_value("gene_model_bigdata_minsize_auto", 25000);  

///SET ANALYSIS
  set_value("gc_correction", 1.1);
  set_value("gene_covar_max_miss", 0.05);
  set_value("gene_covar_impute_type", "median");      
  set_value("gene_covar_truncate", 5);  
  set_value("gene_zstat_truncate_low", 3);
  set_value("gene_zstat_truncate_high", 6);
  set_value("gene_vars_min_sd", 1e-10);  
  set_value("genelevel_direction_sets", 1);    
  set_value("genelevel_direction_covar", 0);
  set_value("genelevel_direction_interaction-ss", 1); 
  set_value("genelevel_direction_interaction-sc", 0); 
  set_value("genelevel_abbreviate", 30);  
  set_value("genelevel_prune_pc", 0.1);
  set_value("genelevel_interaction_sc-size", 100);
  set_value("genelevel_interaction_ss-size", 25);
  set_value("genelevel_interaction_ss-prop", 0.1);
  set_marker("genelevel_residualize");

///MISC
  set_value("rng_seed", get_timestamp());

///flag order constraints  
  set_order("batch");
  set_order("bfile");  
  set_order("data");
  set_order("gene-model");
  set_order("pheno");
  set_order("covar");  
}

void Settings::check_flag_consistency() {
  bool snp_data = flag_set("bfile") || flag_set("col-data") || flag_set("row-data");
  bool gene_data = (snp_data && !flag_set("genes_only") && !flag_set("batch")) || flag_set("gene-results");
  bool set_data = flag_set("set-annot") || flag_set("gene-covar");
  bool analysis = snp_data || gene_data;
  bool annotation = flag_set("annotate");

  ///ANNOTATION
  flag_requires_set("annotate snp-loc gene-loc");
  if (analysis) flag_invalid("annotate", "when running analysis");

  ///SNP DATA
  flag_conflict_set("bfile col-data row-data");
  flag_conflict("bfile col-data row-data gene-annot", "annotate");
  flag_conflict("pheno covar burden", "pval");
  if (!snp_data) flag_invalid("gene-annot genes-only gene-settings covar pheno batch snp-wise gene-model pval rare burden partition", "except when loading SNP data");  
  flag_requires("bfile col-data row-data", "gene-annot");
  flag_requires("burden-weights", "burden");

  ///GENE DATA
  if (!gene_data) flag_invalid("set-annot gene-covar", "unless raw gene results are available (from file or analysis)");
  if (snp_data) flag_invalid("gene-results merge", "when loading SNP data");
  flag_conflict_set("gene-results annotate merge");
  
  ///SET DATA
  if (set_data) flag_invalid("genes-only batch", "when flags --set-annot or --gene-covar are set");
  else flag_invalid("settings model", "except when flags --set-annot or --gene-covar are set");
  flag_requires("set-perm", "set-annot");
 
  ///OTHER
  if (analysis) flag_invalid("meta merge", "when running analysis");
  if (annotation) flag_invalid("meta merge", "when performing annotation");
}

void Settings::parse(Flag& arg) {
  typedef pair<int,int> range;
  typedef pair<double,double> drange;
  
/// --out _FLAG_
  if (arg.name("out", 1)) {
    pair<string,string> path = Utils::split_path(arg[0]);
    if (path.second == "") parse_error("--out flag specifies a directory, should be a file prefix");
    if (path.first != "") check_dir(path.first, arg.name());
    set_value("out_prefix_base", arg[0]);

/// --annotate _FLAG_
  } else if (arg.name("annotate") && arg.set_modifiers(0)) {
    set_value("do_annotate", true);
    if (arg.check_modifier("window", range(1,2))) {
      double window_up = convert_double(arg["window"][0], arg.name(), "window");
      double window_down = arg["window"].size() > 1 ? convert_double(arg["window"][1], arg.name(), "window") : window_up;
      check_range_min(window_up, 0, arg.name(), "window", true);    
      check_range_min(window_down, 0, arg.name(), "window", true);          
      set_value("window_up", long(1000*window_up));
      set_value("window_down", long(1000*window_down));
    } 
    if (arg.check_modifier("ignore-strand")) set_marker("ignore_strand");
    if (arg.check_modifier("nonhuman")) set_marker("nonhuman_chr");
    if (arg.check_modifier("filter", 1)) set_value("annot_filter_file", check_file(arg["filter"][0], arg.name(), "filter"));

/// --gene-loc _FLAG_
  } else if (arg.name("gene-loc", 1)) {
    set_value("gene_loc_file", check_file(arg[0], arg.name()));

/// --snp-loc _FLAG_
  } else if (arg.name("snp-loc", 1)) {
    set_value("snp_loc_file", check_file(arg[0], arg.name()));

/// --bfile --col-data --row-data _FLAG_
  } else if (arg.name("bfile") && arg.set_modifiers(1)) {bool plink_data = arg.name("bfile"); 
    string filename = arg[0];
    if (plink_data) {
      set_value("data_type", "binary_geno");
      string chr_string = gets("batch_chr_string", "");
      if (chr_string != "") {
        filename = Utils::str_replace(filename, "#CHR#", chr_string);
        if (geti("batch_chr_filter", 0) == 1 && filename != arg[0]) set_value("batch_chr_filter", 0);
      }
    
      check_file(filename + ".bed", arg.name());
      check_file(filename + ".bim", arg.name());
      check_file(filename + ".fam", arg.name());
      
      set_value("plink_prefix", filename);
      set_value("do_analysis", "plink");
    }
        
    if (arg.check_modifier("synonyms", 1)) {
      if (arg["synonyms"][0] != "0") set_value("snp_synonyms", check_file(arg["synonyms"][0], arg.name(), "synonyms"));
    } else if (plink_data && Utils::is_file(filename + ".synonyms")) {
      set_value("snp_synonyms", filename + ".synonyms");
      set_marker("snp_synonyms_auto");
    }
    
    if (arg.check_modifier("synonym-dup", 1)) {
      if (!(*this)["snp_synonyms"]) mod_error(arg.name(), "synonym-dup", "modifier cannot be set unless a synonym has been specified or auto-detected");
      if (check_enum(arg["synonym-dup"][0], "drop drop-dup skip skip-dup error", arg.name(), "synonym-dup")) set_value("snp_synonyms_dup", arg["synonym-dup"][0]);
    } else set_value("snp_synonyms_dup", "skip");

/// --batch _FLAG_
  } else if (arg.name("batch") && arg.set_modifiers(2,0)) {
    if (arg[1] == "chr") {
      int chr = 0;
      bool human = !arg.check_modifier("nonhuman");
      if (!Utils::chr_val(arg[0], chr, human) || chr == 0) flag_error("batch", "chromosome code '" + arg[0] + "' not recognised");

      set_value("batch_chr", chr);
      set_value("batch_chr_string", arg[0]);      
      if (!human) set_marker("batch_chr_nonhuman");
      
      int chr_filter = 1;
      if (arg.check_modifier("chr-filter", range(0,1))) {
        if (arg["chr-filter"].size() > 0) {
          if (arg["chr-filter"][0] == "1") chr_filter = 2;
          else if (arg["chr-filter"][0] == "0") chr_filter = 0;        
          else mod_error("batch", "chr_filter", "unknown value " + Utils::quote(arg["chr-filter"][0]));
        } else chr_filter = 2;
      }
      set_value("batch_chr_filter", chr_filter);
    } else {
      int total = convert_long(arg[1], arg.name());
      int index = convert_long(arg[0], arg.name());
      check_range_min(total, 1, arg.name(), false);
      check_range(index, range(1,total), arg.name(), true);
      set_value("batch", total);
      set_value("batch_index", index);
      
      if (arg.check_modifier("nonhuman")) flag_error("batch", "modifier 'nonhuman' can only be used in chromosome batch mode");
      if (arg.check_modifier("chr-filter")) flag_error("chr-filter", "modifier 'chr-filter' can only be used in chromosome batch mode");      
    }

/// --merge _FLAG_
  } else if (arg.name("merge", range(1,2))) {
    pair<string,string> path = Utils::split_path(arg[0]);
    if (path.second == "") parse_error("--merge flag specifies a directory, should be a file prefix");
    if (path.first != "") check_dir(path.first, arg.name());
    set_marker("do_merge");
    set_value("merge_prefix", arg[0]);
    if (arg.size() >= 2) {
      if (arg[1] == "chr") set_marker("merge_batch_chr");
      else {
        int no_batches = convert_long(arg[1], arg.name());
        check_range_min(no_batches, 1, arg.name(), true);
        set_value("merge_batch_count", no_batches);
      }
    }
    
/// --meta _FLAG_
  } else if (arg.name("meta") && arg.set_modifiers(0,1)) {
    mod_conflict(arg, "genes raw genes-file raw-file");
    mod_conflict(arg, "unweighted", "genes-file raw-file");
    mod_one_required(arg, "use-partitions", "genes genes-file");
  
    set_marker("do_meta");
    vector<string> meta_genes;
    bool prefix = arg.check_modifier("prefix");
    if (arg.check_modifier("genes", 2, false) || arg.check_modifier("genes-file", 1)) {
      set_marker("do_meta_genes");
      if (arg.has_modifier("genes")) {
        for (int i = 0; i < arg["genes"].size(); i++) {
          string filename = arg["genes"][i];
          if (prefix) filename = win_txt(filename + ".genes.out");
          else check_suffix(filename, win_txt(".genes.out"), arg.name(), "genes");
          check_file(filename, arg.name(), "genes");
          meta_genes.push_back(filename);
        }
      } else parse_metaweights(arg["genes-file"][0], meta_genes, "genes", prefix);
    } else if (arg.check_modifier("raw", 2, false) || arg.check_modifier("raw-file", 1)) {
      set_marker("do_meta_raw");
      if (arg.has_modifier("raw")) {
        for (int i = 0; i < arg["raw"].size(); i++) {
          string filename = arg["raw"][i];
          if (prefix) filename.append(".genes.raw");
          else check_suffix(filename, ".genes.raw", arg.name(), "raw");
          check_file(filename, arg.name(), "raw");
          meta_genes.push_back(filename);
        }
      } else parse_metaweights(arg["raw-file"][0], meta_genes, "raw", prefix);
    } else flag_error("meta", "modifier 'genes', 'genes-file', 'raw' or 'raw-file' must be set");
    
    set<string> file_names;
    for (int i = 0; i < meta_genes.size(); i++) {
      if (file_names.find(meta_genes[i]) == file_names.end()) file_names.insert(meta_genes[i]);
      else flag_error("meta", "file '" + meta_genes[i] + "' is specified as input more than once");
    }
    set_value("gene_meta_files", meta_genes);

    if (arg.check_modifier("unweighted")) meta_weights.assign(meta_genes.size(), Triple<double>(1,1,1));
    if (arg.check_modifier("correlations", 1)) parse_metacorrs(arg["correlations"][0], meta_genes.size());


    if (arg.check_modifier("truncate", range(1,2))) {
      double val = convert_double(arg["truncate"][0], arg.name(), "truncate");
      if (val != 0) check_range_min(val, 3, arg.name(), "truncate", true);
      set_value("metagene_truncate_low", val);

      if (arg["truncate"].size() == 2) {
        val = convert_double(arg["truncate"][1], arg.name(), "truncate");
        if (val != 0) check_range_min(val, 3, arg.name(), "truncate", true);
      }
      set_value("metagene_truncate_high", val);
    }
    if (arg.check_modifier("use-partitions")) set_marker("metagene_partitions");

/// --covar _FLAG_
  } else if (arg.name("covar") && arg.set_modifiers(0,1)) {
    mod_conflict(arg, "include exclude include-all"); mod_conflict(arg, "interact interact-sex interact-snp");
    mod_required(arg, "include include-all exclude interact", "file");
  
    if (arg.check_modifier("file", 1)) set_value("covar_file", check_file(arg["file"][0], arg.name(), "file"));
    if (arg.check_modifier("chrX-use-sex", range(0,1))) {
      string value_str = (arg["chrX-use-sex"].size() > 0) ? arg["chrX-use-sex"][0] : "";
      bool value = check_bool(value_str, arg.name(), "chrX-use-sex"); set_value("chrX_sex_covar", value);
      if (!value && arg.has_modifier("use-sex")) flag_remark(arg.name(), string("setting 'chrX-use-sex=") + value_str + "' is overridden by modifier 'use-sex'");
    }
    if (arg.check_modifier("use-sex")) set_marker("covar_use_sex");
    if (arg.check_modifier("include-all")) set_marker("covar_include");
    if (arg.check_modifier("include", 1, false)) {
      covar_include = arg["include"].get_values();
      set_marker("covar_include");
    }
    if (arg.check_modifier("exclude", 1, false)) {
      covar_exclude = arg["exclude"].get_values();
      set_marker("covar_exclude");
    }
    if (arg.check_modifier("interact", 1)) {
      set_value("covar_interact", arg["interact"][0]);
      if (!arg.has_modifier("include-all")) covar_include.push_back(arg["interact"][0]);
    }
    if (arg.check_modifier("interact-sex")) {
      set_marker("covar_interact");
      set_marker("covar_interact_sex");
      set_marker("covar_use_sex");
    }

    if ((*this)["covar_interact"]) set_marker("use_residuals");
    if ((*this)["covar_use_sex"]) {
      set_marker("drop_missing_sex");
      set_value("chrX_sex_covar", false);
    }

    if (arg.check_modifier("include-snps", 1, false) + arg.check_modifier("interact-snp", 1)) {vector<string> snps;
      if (arg.has_modifier("include-snps")) snps = arg["include-snps"].get_values();
      if (arg.has_modifier("interact-snp")) {
        snps.push_back(arg["interact-snp"][0]);
        set_marker("covar_interact");
        set_value("covar_interact_snp", arg["interact-snp"][0]);        
      }
      set<string> uniq(snps.begin(), snps.end()); snps.assign(uniq.begin(), uniq.end());
      set_value("covar_snps", snps);
    }
    if ((*this)["covar_interact"] && flag_set("partition")) flag_error(arg.name(), "interaction models currently not available when using flag --partition");              

/// --pheno _FLAG_
  } else if (arg.name("pheno") && arg.set_modifiers(0,1)) {
    mod_required(arg, "file"); mod_conflict(arg, "use use-all");
    mod_required(arg, "permute", "use-all"); 
    if (arg.check_modifier("file", 1)) {
      set_value("pheno_file", check_file(arg["file"][0], arg.name(), "file"));
      set_value("alt_pheno", "1");
    }
    if (arg.check_modifier("use-all")) {
      set_marker("alt_pheno_all");
      set_value("alt_pheno", "1");          
    } else if (arg.check_modifier("use", 1)) set_value("alt_pheno", arg["use"][0]);

    if (arg.check_modifier("permute", 1)) {
      int perm_count = convert_long(arg["permute"][0], arg.name());
      check_range_min(perm_count, 1000, arg.name(), true);
      set_value("multi_pheno_permute", perm_count);
    }


/// --gene-annot _FLAG_
  } else if (arg.name("gene-annot", 1)) {
    set_value("gene_annot_file", check_file(arg[0], arg.name(), true));

/// --gene-results _FLAG_
  } else if (arg.name("gene-results") && arg.set_modifiers(1,0)) {
    set_value("gene_results_file", check_file(arg[0], arg.name()));
    set_marker("has_results_file");
    
    if (arg.check_modifier("overwrite", 1)) set_value("gene_results_overwrite", check_file(arg["overwrite"][0], arg.name(), "overwrite"));
    if (arg.check_modifier("rescale", 1)) {
      double scale = convert_double(arg["rescale"][0], arg.name(), "rescale");
      check_range(scale, drange(0, 1), arg.name(), "rescale", true);
      set_value("gene_corr_rescale", scale);
    }    

/// --set-annot _FLAG_
  } else if (arg.name("set-annot") && arg.set_modifiers(1,0)) {
    mod_conflict(arg, "col", "gene-col set-col"); mod_required(arg, "gene-col", "set-col", true);
  
    set_value("set_annot_file", check_file(arg[0], arg.name()));
    set_marker("has_set_annot");
    set_marker("do_set_analysis");
    set_marker("need_gene_raw");

    if (arg.check_modifier("col", 2)) {
      set_marker("has_set_annot_col");
      set_value("set_annot_col_gene", arg["col"][0]);
      set_value("set_annot_col_set", arg["col"][1]);
    }
    if (arg.check_modifier("gene-col", 1) && arg.check_modifier("set-col", 1)) {
      set_marker("has_set_annot_col");
      set_value("set_annot_col_gene", arg["gene-col"][0]);
      set_value("set_annot_col_set", arg["set-col"][0]);
    }
    if ((*this)["has_set_annot_col"] && (gets("set_annot_col_gene") == gets("set_annot_col_set"))) flag_error("set-annot", "same column specified for gene ID and gene-set ID");
    
    if (arg.check_modifier("allow-mono")) set_marker("gene_sets_allow_mono");    

/// --gene-covar _FLAG_
  } else if (arg.name("gene-covar") && arg.set_modifiers(1,0)) {
    mod_conflict(arg, "filter-read filter-skip");
  
    set_value("gene_covar_file", check_file(arg[0], arg.name()));
    set_marker("has_gene_covar");
    set_marker("do_set_analysis");
    set_marker("need_gene_raw");
    
    if (arg.check_modifier("filter-read", 1, false)) {
      gene_covar_include = arg["filter-read"].get_values();
      set_marker("gene_covar_include");
    }
    if (arg.check_modifier("filter-skip", 1, false)) {
      gene_covar_exclude = arg["filter-skip"].get_values();
      set_marker("gene_covar_exclude");
    }
   
    bool fill_rows = arg.check_modifier("missing-genes", 1) && check_enum(arg["missing-genes"][0], "drop fill", arg.name(), "missing-genes") && arg["missing-genes"][0] == "fill";
    if (fill_rows) set_marker("gene_covar_fill_rows"); 
    
    
    if (arg.check_modifier("missing-values", 1)) {
      string impute = arg["missing-values"][0];
      if (!check_enum(impute, "mean median drop")) { 
        impute = "value";
        set_value("gene_covar_impute_value", convert_double(arg["missing-values"][0], arg.name(), "missing-values"));
      }
      if (impute == "drop" && fill_rows) mod_error(arg.name(), "missing-genes", "cannot be set to \"fill\" if modifier 'missing-values' is set to \"drop\"");
      set_value("gene_covar_impute_type", impute);    
    }
    
    if (arg.check_modifier("max-miss", 1)) {
      double miss = convert_double(arg["max-miss"][0], arg.name(), "max-miss");
      check_range(miss, drange(0, 0.2), arg.name(), "max-miss", true);
      set_value("gene_covar_max_miss", miss);
    }

/// --model _FLAG_
  } else if (arg.name("model") && arg.set_modifiers(0,1)) {
    bool has_sets = flag_set("set-annot"), has_covar = flag_set("gene-covar");
    mod_conflict(arg, "direction", "direction-sets direction-covar");
    mod_conflict(arg, "direction-interaction", "direction-interaction-ss direction-interaction-sc");
    mod_conflict(arg, "interaction joint interaction-pairs joint-pairs interaction-all interaction-each");    
    mod_conflict(arg, "interaction joint", "analyse");

    if (arg.check_modifier("analyse", 1, false)) {
      if (check_enum(arg["analyse"][0], "all sets cov covar covariates list file", arg.name(), "analyse")) {
        string value = (arg["analyse"][0] == "covariates" || arg["analyse"][0] == "cov") ? "covar" : arg["analyse"][0];
        if (value == "sets" && !has_sets) mod_error(arg.name(), "analyse", "cannot use mode '" + arg["analyse"][0] + "' unless flag --set-annot is set"); 
        if (value == "covar" && !has_covar) mod_error(arg.name(), "analyse", "cannot use mode '" + arg["analyse"][0] + "' unless flag --gene-covar is set");         
                
        if (value == "list") {
          if (arg["analyse"].size() < 2) mod_error(arg.name(), "analyse", "when using mode 'list', at least two values should be specified");          
          vector<string> list = arg["analyse"].get_values(); list.erase(list.begin());
          set_value("genelevel_analysis_use_list", list); 
        } else if (value == "file") {
          if (arg["analyse"].size() != 2) mod_error(arg.name(), "analyse", "when using mode 'file', exactly two values should be specified");          
          vector<string> list = read_list(arg["analyse"][1], false, arg.name(), "analyse");        
          set_value("genelevel_analysis_use_list", list); value = "list";
        } else if (arg["analyse"].size() > 1) mod_error(arg.name(), "analyse", "when using mode '" + arg["analyse"][0] + "' no additional values should be supplied");
     
        set_value("genelevel_analysis_use", value);
      }
    }
    
    if (arg.check_modifier("joint", 1)) {
      vector<string> list = read_list(arg["joint"][0], true, arg.name(), "joint");
      set_value("genelevel_joint_list", list);      
      set_value("genelevel_analysis_mode", "joint_list");
    }

    if (arg.check_modifier("interaction", 1)) {
      vector<string> list = read_list(arg["interaction"][0], true, arg.name(), "interaction");
      set_value("genelevel_interaction_list", list);      
      set_value("genelevel_analysis_mode", "interaction_list");
    }

    string inter_suffix[2] = {"each", "all"};
    for (int i = 0; i < 2; i++) {
      string mod = string("interaction-") + inter_suffix[i], param = string("interaction_") + inter_suffix[i];
      if (arg.check_modifier(mod, 1, false)) {vector<string> list;
        if (arg[mod][0] == "file") {
          if (arg[mod].size() != 2) mod_error(arg.name(), mod, "when using mode 'file', exactly two values should be specified");          
          list = read_list(arg[mod][1], false, arg.name(), mod);
        } else list = arg[mod].get_values();
        set_value(string("genelevel_") + param, list);      
        set_value("genelevel_analysis_mode", param);
      }    
    }

    if (arg.check_modifier("joint-pairs")) set_value("genelevel_analysis_mode", "joint_pairs");
    if (arg.check_modifier("interaction-pairs")) set_value("genelevel_analysis_mode", "interaction_pairs");    
        
    if (arg.check_modifier("interaction-sc-size", 1)) {
      int size = convert_long(arg["interaction-sc-size"][0], arg.name(), "interaction-sc-size");
      check_range_min(size, 25, arg.name(), "interaction-sc-size", true);
      set_value("genelevel_interaction_sc-size", size);
    }
    
    if (arg.check_modifier("interaction-ss-size", 2)) {
      int size = convert_long(arg["interaction-ss-size"][0], arg.name(), "interaction-ss-size");
      double perc = convert_double(arg["interaction-ss-size"][1], arg.name(), "interaction-ss-size");
      if (size < 10) mod_error(arg.name(), "interaction-ss-size", "first value must be greater than or equal to 10");
      if (perc < 0.01 || perc > 0.49) mod_error(arg.name(), "interaction-ss-size", "second value must be between 0.01 and 0.49");
      set_value("genelevel_interaction_ss-size", size);
      set_value("genelevel_interaction_ss-prop", perc);
    }
  
    if (arg.check_modifier("correct", 1, false)) {
      string& mode = arg["correct"][0];
      if (mode == "all" || mode == "none") {
        if (arg["correct"].size() > 1) mod_error("model", "correct", "only one value is allowed when using mode " + Utils::quote(mode));
      } else if (mode == "include" || mode == "exclude") {
        if (arg["correct"].size() == 1) mod_error("model", "correct", "at least one variable must be specified when using mode " + Utils::quote(mode));
        vector<string> list;
        for (int i = 1; i < arg["correct"].size(); i++) {
          string var = Utils::lowercase(arg["correct"][i]);
          if (check_enum(var, "size genesize density mac sampsize n", arg.name(), "correct")) {
            if (var == "size" || var == "genesize") var = "nsnps";
            else if (var == "sampsize" || var == "n") var = "nsamp";
            else if (var == "density") var = "nparam";
            list.push_back(var);
          }
        }
        set_value("genelevel_correct_list", list);
      } else mod_error(arg.name(), "correct", "unknown mode " + Utils::quote(mode));
      set_value("genelevel_correct", mode);
    }
    
    string cond_suffix[3] = {"", "hide", "residualize"};
    for (int c = 0; c < 3; c++) {
      string mod = string("condition") + (cond_suffix[c] != "" ? "-" : "") + cond_suffix[c];
      if (arg.check_modifier(mod, 1, false)) {vector<string> list;
        if (arg[mod][0] == "file") {
          if (arg[mod].size() != 2) mod_error(arg.name(), mod, "when using mode 'file', exactly two values should be specified");          
          list = read_list(arg[mod][1], false, arg.name(), mod);
        } else list = arg[mod].get_values();
        set_value(string("genelevel_condition") + (cond_suffix[c] != "" ? "_" : "") + cond_suffix[c] + "_list", list);
        if (mod == "condition-residualize")  set_marker("genelevel_residualize");
      }
    }
    
    if (arg.check_modifier("condition-interaction", 2, false)) {
      vector<string> list;
      if (arg["condition-interaction"][0] == "file") {
        if (arg["condition-interaction"].size() != 2) mod_error(arg.name(), "condition-interaction", "when using mode 'file', exactly two values should be specified");          
        list = read_list(arg["condition-interaction"][1], false, arg.name(), "condition-interaction");
      } else list = arg["condition-interaction"].get_values();
      if (list.size() % 2 != 0) mod_error(arg.name(), "condition-interaction", "the number of specified variables should be even");
      set_value("genelevel_condition_interaction_list", list);
    }

    if (arg.check_modifier("self-contained")) {
      if (!has_sets) mod_error(arg.name(), "self-contained", "cannot be used unless flag --set-annot is set");
      set_marker("genelevel_sets_selfcontained");
    }
    
    if (arg.check_modifier("allow-collinear")) set_marker("genelevel_allow_collinear");
    
    set<string> dir_mods = arg.check_modifiers("direction direction-sets direction-covar direction-interaction direction-interaction-ss direction-interaction-sc", 1);
    if (!dir_mods.empty()) {
      for (set<string>::iterator mod = dir_mods.begin(); mod != dir_mods.end(); ++mod) {string value = arg[*mod][0];
        check_enum(value, "two-sided twosided both two positive pos greater negative neg smaller", arg.name(), *mod);
        int direction = !Utils::contains("positive pos greater", value) ? (!Utils::contains("negative neg smaller", value) ? 0 : -1) : 1;
        if (*mod == "direction") {
          set_value("genelevel_direction_sets", direction); set_value("genelevel_direction_covar", direction);
        } else if (*mod == "direction-sets") {
          if (!has_sets) mod_error(arg.name(), "direction-sets", "cannot be used unless flag --set-annot is set");
          set_value("genelevel_direction_sets", direction); 
        } else if (*mod == "direction-covar") {
          if (!has_covar) mod_error(arg.name(), "direction-covar", "cannot be used unless flag --gene-covar is set");
          set_value("genelevel_direction_covar", direction);
        } else if (*mod == "direction-interaction") {
          set_value("genelevel_direction_interaction-ss", direction); set_value("genelevel_direction_interaction-sc", direction);
        } else if (*mod == "direction-interaction-ss") {
          set_value("genelevel_direction_interaction-ss", direction);
        } else if (*mod == "direction-interaction-sc") {
          set_value("genelevel_direction_interaction-sc", direction);
        }
      }
    }

    double alpha = 0.05;
    if (arg.check_modifier("alpha", 1)) {
      alpha = convert_double(arg["alpha"][0], arg.name(), "alpha");
      check_range_min(alpha, 0, arg.name(), "alpha", false);
      check_range_max(alpha, 0.5, arg.name(), "alpha", true);
      set_value("genelevel_alpha", alpha);    
    }

/// --genes-only _FLAG_
  } else if (arg.name("genes-only", 0)) {
    set_marker("genes_only");
    if (gets("gene_model_bigdata", "") == "auto") unset_marker("gene_model_bigdata");

/// --big-data _FLAG_
  } else if (arg.name("big-data", range(0,1))) {
    int set = 1;
    if (arg.size() > 0) {
      if (arg[0] == "0" || arg[0] == "off") set = 0;
      else if (arg[0] == "partial") set = 1;      
      else if (arg[0] == "on" || arg[0] == "full") set = 2;
      else flag_error("big-data", "value '" + arg[0] + "' is not recognised");
    }
    if (set > 0) set_value("gene_model_bigdata", set == 1 ? "partial" : "full");
    else unset_marker("gene_model_bigdata");

/// --gene-settings _FLAG_
  } else if (arg.name("gene-settings") && arg.set_modifiers(0,1)) {
    mod_conflict(arg, "fixed-permp", "adap-permp"); mod_conflict(arg, "snp-include", "snp-exclude"); mod_conflict(arg, "indiv-include", "indiv-exclude");
    
    if (arg.check_modifier("min-perm", 1)) {
      long no_perm = convert_long(arg["min-perm"][0], arg.name(), "min-perm");
      check_range_min(no_perm, 1000, arg.name(), "min-perm", true);
      set_value("min_perm", no_perm);
    }
    if (arg.check_modifier("fixed-permp", range(0,1))) {
      set_marker("compute_permp");
      set_marker("fixed_perm");
      if (arg["fixed-permp"].size() > 0) {
        if (arg.has_modifier("min-perm")) flag_error("gene-settings", "modifier 'min-perm' is not allowed when specifying number of permutations through 'fixed-permp' modifier");
        long no_perm = convert_long(arg["fixed-permp"][0], arg.name(), "fixed-permp");
        check_range_min(no_perm, 1000, arg.name(), "fixed-permp", true);
        set_value("min_perm", no_perm);
      }
    }
    if (arg.check_modifier("adap-permp", range(0,2))) {
      set_marker("compute_permp");
      set_marker("adap_perm");

      if (arg["adap-permp"].size() > 0) {
        long max_perm = convert_long(arg["adap-permp"][0], arg.name(), "adap-permp");
        if (arg.has_modifier("min-perm") && max_perm <= geti("min_perm")) mod_error("gene-settings", "adap-permp", "first value must be greater than value for 'min-perm' modifier");        
        if (max_perm <= 1000) mod_error("gene-settings", "adap-permp", "first value must be greater than 1000");
        set_value("max_perm", max_perm);
      }
      if (arg["adap-permp"].size() > 1) {
        long adap_count = convert_long(arg["adap-permp"][1], arg.name(), "adap-permp");
        if (adap_count <= 0) mod_error("gene-settings", "adap-permp", "second value must be greater than 0");
        set_value("adap_count", adap_count);
      }
    }

    if (arg.check_modifier("snp-max-miss", 1)) {
      double miss = convert_double(arg["snp-max-miss"][0], arg.name(), "snp-max-miss");
      double max_max = this->getn("snp_max_miss_max", 1);
      check_range(miss, drange(0, max_max), arg.name(), "snp-max-miss", true);
      set_value("snp_max_miss", miss);
    }
    if (arg.check_modifier("snp-min-maf", 1)) {
      check_datatype("binary_geno", arg.name(), "snp-min-maf");
      double maf = convert_double(arg["snp-min-maf"][0], arg.name(), "snp-min-maf");
      check_range(maf, drange(0, 0.5), arg.name(), "snp-min-maf", true);
      set_value("snp_min_maf", maf);
    }
    if (arg.check_modifier("snp-max-maf", 1)) {
      check_datatype("binary_geno", arg.name(), "snp-max-maf");
      double maf = convert_double(arg["snp-max-maf"][0], arg.name(), "snp-max-maf");
      check_range_min(maf, 0, arg.name(), "snp-max-maf", false);
      check_range_max(maf, 0.5, arg.name(), "snp-max-maf", true);
      if (arg.has_modifier("snp-min-maf") && maf <= (*this)["snp_min_maf"]) mod_error("gene-settings", "snp-max-maf", "value must be higher than for modifier 'snp-min-maf'");
      set_value("snp_max_maf", maf);
    }
    if (arg.check_modifier("snp-min-mac", 1)) {
      check_datatype("binary_geno", arg.name(), "snp-min-mac");
      long mac = convert_double(arg["snp-min-mac"][0], arg.name(), "snp-min-mac");
      check_range_min(mac, 0, arg.name(), "snp-min-mac", true);
      set_value("snp_min_mac", mac);
    }
    if (arg.check_modifier("snp-max-mac", 1)) {
      check_datatype("binary_geno", arg.name(), "snp-max-mac");
      long mac = convert_double(arg["snp-max-mac"][0], arg.name(), "snp-max-mac");
      check_range_min(mac, 0, arg.name(), "snp-max-mac", false);
      if (arg.has_modifier("snp-min-mac") && mac <= (*this)["snp_min_mac"]) mod_error("gene-settings", "snp-max-mac", "value must be higher than for modifier 'snp-min-mac'"); 
      set_value("snp_max_mac", mac);
    }
    if (arg.check_modifier("snp-diff", 1)) {
      double pval = convert_double(arg["snp-diff"][0], arg.name(), "snp-diff");
      check_range(pval, drange(1e-20, 1), arg.name(), "snp-diff", false);
      set_value("snp_diff_pval", pval);
    }
    if (arg.check_modifier("prune", 1)) {
      double prune = convert_double(arg["prune"][0], arg.name(), "prune");
      check_range(prune, drange(0.01, 0.9999), arg.name(), "prune", true);
      set_value("prune_cutoff", prune);
    }
    if (arg.check_modifier("prune-powersum", 1)) {
      double prune = convert_double(arg["prune-powersum"][0], arg.name(), "prune");
      check_range(prune, drange(0.01, 0.9999), arg.name(), "prune", true);
      set_value("prune_varimax_cutoff", prune);
    }
    if (arg.check_modifier("prune-prop", 1)) {
      double prune = convert_double(arg["prune-prop"][0], arg.name(), "prune-prop");
      check_range(prune, drange(0.01, 1), arg.name(), "prune-prop", true);
      set_value("prune_prop", prune);
    }
    if (arg.check_modifier("prune-count", 1)) {
      long prune = convert_long(arg["prune-count"][0], arg.name(), "prune-count");
      check_range_min(prune, 1, arg.name(), "prune-count", true);
      set_value("prune_count", prune);
    }
    if (arg.check_modifier("prune-value", range(0,1))) {
      if (arg["prune-value"].size() > 0) {
        double prune = convert_double(arg["prune-value"][0], arg.name(), "prune-value");
        check_range_min(prune, 0, arg.name(), "prune-value", true);
        set_value("prune_value", prune);
        set_value("prune_varimax_value", prune);        
      } else {
        set_value("prune_value", 1);
        set_value("prune_varimax_value", 1);
      }
    }
    if (arg.check_modifier("snp-include", 1, true)) {
      set_value("snp_filter_file", check_file(arg["snp-include"][0], arg.name(), "snp-include"));
      set_value("snp_filter", "include");
    }
    if (arg.check_modifier("snp-exclude", 1, true)) {
      set_value("snp_filter_file", check_file(arg["snp-exclude"][0], arg.name(), "snp-exclude"));
      set_value("snp_filter", "exclude");
    }     
    if (arg.check_modifier("indiv-include", 1, true)) {
      set_value("indiv_filter_file", check_file(arg["indiv-include"][0], arg.name(), "indiv-include"));
      set_value("indiv_filter", "include");
    }
    if (arg.check_modifier("indiv-exclude", 1, true)) {
      set_value("indiv_filter_file", check_file(arg["indiv-exclude"][0], arg.name(), "indiv-exclude"));
      set_value("indiv_filter", "exclude");
    }


/// --settings _FLAG_
  } else if (arg.name("settings") && arg.set_modifiers(0,1)) {
    mod_conflict(arg, "gene-include", "gene-exclude");
    if (arg.check_modifier("outlier", 2)) {
      double val = convert_double(arg["outlier"][0], arg.name(), "outlier");
      if (val != 0) check_range_min(val, 3, arg.name(), "outlier", true);
      set_value("gene_zstat_truncate_low", val);
 
      val = convert_double(arg["outlier"][1], arg.name(), "outlier");
      if (val != 0) check_range_min(val, 3, arg.name(), "outlier", true);
      set_value("gene_zstat_truncate_high", val);      
    }
    
    if (arg.check_modifier("covar-outlier", 1)) {
      double truncate = convert_double(arg["covar-outlier"][0], arg.name(), "covar-outlier");
      check_range(truncate, drange(1, 100), arg.name(), "covar-outlier", true);
      set_value("gene_covar_truncate", truncate);
    }
  
    if (arg.check_modifier("gene-include",1)) {
      set_value("gene_filter_file", check_file(arg["gene-include"][0], arg.name(), "gene-include"));
      set_value("gene_filter", "include");
    }
    if (arg.check_modifier("gene-exclude",1)) {
      set_value("gene_filter_file", check_file(arg["gene-exclude"][0], arg.name(), "gene-exclude"));
      set_value("gene_filter", "exclude");
    }
    
    if (arg.check_modifier("abbreviate", range(1,2))) {
      if (arg["abbreviate"][0] == "file") set_marker("genelevel_abbreviate_file");      
      else if (arg["abbreviate"][0] != "col" && arg["abbreviate"][0] != "column") {
        int length = convert_long(arg["abbreviate"][0], arg.name(), "abbreviate");
        if (length != 0) check_range_min(length, 20, arg.name(), "abbreviate", true);
        set_value("genelevel_abbreviate", length);
      }
      if (arg["abbreviate"].size() == 2) {
        if (arg["abbreviate"][1] == "file") set_marker("genelevel_abbreviate_file");      
        else if (arg["abbreviate"][1] != "col" && arg["abbreviate"][1] != "column") mod_error(arg.name(), "abbreviate", string("unknown value ") + Utils::quote(arg["abbreviate"][1]));
      }
    }
    
    if (arg.check_modifier("gene-info")) set_marker("genelevel_print_geneinfo");

/// --burden _FLAG_
  } else if (arg.name("burden") && arg.set_modifiers(1)) {
    check_datatype("binary_geno", arg.name());
    unset_marker("rare_autoburden");

    double thresh = 0; bool freq_weighted = true;
    long multi_count = this->geti("rare_multi_max_count", 25);
    if (arg[0] == "all" || arg[0] == "all-plain") {
      thresh = 0.99; multi_count = 0;
      if (arg[0] == "all-plain") freq_weighted = false;
    } else if (arg[0] != "off") {
      thresh = convert_double(arg[0], arg.name());
      check_range_min(thresh, 0, arg.name(), true);
      if ((thresh > 0.5 && thresh < 1) || (thresh > 1 && thresh != floor(thresh))) flag_error(arg.name(), "value " + Utils::num_quote(thresh, '(') + " must be a decimal no greater than 0.5 (MAF) or a whole number no smaller than 1 (MAC)");
    }
    if (thresh <= 0) {if (arg.nmod() > 0) {flag_error("burden", "cannot use modifiers when turning off burden scoring"); return;}}
    else set_marker("do_rare");


    if (thresh < 1) set_value("rare_cutoff_maf", thresh);
    else set_value("rare_cutoff_mac", int(floor(thresh)));

    if (arg.check_modifier("rare-only")) set_marker("rare_only");
    if (arg.check_modifier("rare-sum")) {set_marker("rare_only"); multi_count = 0;}
    if (arg.check_modifier("rare-plain")) {set_marker("rare_only"); multi_count = 0; freq_weighted = false;}

    if (arg.check_modifier("freq-weighted", 1, true)) {
      if (arg["freq-weighted"][0] == "0") freq_weighted = false;
      else if (arg["freq-weighted"][0] == "1") freq_weighted = true;
      else if (arg["freq-weighted"][0] != "1") mod_error(arg.name(), "freq-weighted", "unknown value " + Utils::quote(arg["freq-weighted"][0]));
    }
    if (!freq_weighted) set_marker("rare_freq_unweighted");

    if (arg.check_modifier("multi-count", 1)) {
      multi_count = convert_long(arg["multi-count"][0], arg.name(), "multi-count");
      if (multi_count != 0) check_range_min(multi_count, 5, arg.name(), "multi-count", true);
    }
    set_value("rare_multi_max_count", multi_count);

    
/// --burden-weights _FLAG_
  } else if (arg.name("burden-weights") && arg.set_modifiers(1,0)) {
    mod_conflict(arg, "use", "snp-id weight");
    mod_required(arg, "snp-id", "weight", true);
    mod_required(arg, "type-missing", "type-weights");

    set_marker("has_rare_weights");
    set_value("rare_weights_file", check_file(arg[0], arg.name()));

    if (arg.check_modifier("use", 2)) {
      set_value("rare_weights_col_snp", arg["use"][0]);
      set_value("rare_weights_col_weight", arg["use"][1]);
    }
    if (arg.check_modifier("snp-id", 1) && arg.check_modifier("weight", 1)) {
      set_value("rare_weights_col_snp", arg["snp-id"][0]);
      set_value("rare_weights_col_weight", arg["weight"][0]);
    }    
    if ((*this)["rare_weights_col_snp"] && (gets("rare_weights_col_snp") == gets("rare_weights_col_weight"))) flag_error("burden-weights", "same column specified for SNP ID and SNP weight");

    if (arg.check_modifier("duplicate", 1)) {
      if (check_enum(arg["duplicate"][0], "first last drop error", arg.name(), "duplicate")) set_value("rare_weights_duplicate_mode", arg["duplicate"][0]);
    } else set_value("rare_weights_duplicate_mode", "drop");

    double default_weight = arg.check_modifier("default", 1) ? convert_double(arg["default"][0], arg.name(), "default") : 0.0;
    check_range_min(default_weight, 0, arg.name(), "default", true);
    set_value("rare_weights_default", default_weight);

    if (arg.check_modifier("type-weights", 1)) {
      set_value("rare_weights_type_file", check_file(arg["type-weights"][0], arg.name(), "type-weights"));
      if (arg.check_modifier("type-missing", 1)) {
        if (check_enum(arg["type-missing"][0], "default drop error", arg.name(), "duplicate")) set_value("rare_weights_type_missing", arg["type-missing"][0]);
      } else set_value("rare_weights_type_missing", "drop");
    }

/// --gene-model _FLAG_
  } else if (arg.name("gene-model") && arg.set_modifiers(0,1)) {
    mod_conflict(arg, "linreg snp-wise powersum multi tdt");

    if (arg.check_modifier("linreg")) { 
      if (flag_set("pval")) flag_error(arg.name(), "modifier 'linreg' cannot be set when using --pval"); 
      set_value("gene_model", "pcreg");
    }

    if (arg.check_modifier("snp-wise", range(0,2))) {
      set_value("gene_model", "snpwise");
      if (arg["snp-wise"].size() > 0) {
        if (arg["snp-wise"][0] == "top") {
          set_value("gene_snpwise", "top");
          double top = 1;
          if (arg["snp-wise"].size() > 1 && (!Utils::convert_num(arg["snp-wise"][1], top) || top <= 0)) mod_error("gene-model", "snp-wise", "second value for SNP-wise model type 'top' must be numeric and greater than 0");
          if (top >= 1) set_value("gene_snpwise_top_count", long(top));
          else set_value("gene_snpwise_top_fraction", top);      

          if (geti("gene_snpwise_top_count",0) != 1) {
            if (flag_set("partition")) flag_error(arg.name(), "when using flag --partition, snp-wise 'top' model is only available if set to single top SNP (ie. 'top' or 'top,1')");              

            set_value("gene_snpwise_top_xtx_max", 1000);
            set_marker("compute_permp");
            set_marker("no_asym");        
            if (!(*this)["fixed_perm"]) set_marker("adap_perm");            
          } else set_value("gene_model", "snpwise_top1");
        } else if (arg["snp-wise"][0] == "multi") {
          set_value("gene_model", "multi");
          set_value("gene_multi", "snpwise");
        } else if (arg["snp-wise"][0] == "sum" || arg["snp-wise"][0] == "mean" || arg["snp-wise"][0] == "imhof") {
          set_value("gene_snpwise", "imhof");
          if (arg["snp-wise"].size() > 1) {
            int power = 25;
            if (!Utils::convert_num(arg["snp-wise"][1], power) || power < 10 || power > 100) mod_error("gene-model", "snp-wise", "second value for SNP-wise model type 'mean' must be a whole number between 10 and 100");
            double cutoff = pow(10, -power);
            if (cutoff < 1e-100) cutoff = 1e-100;
            set_value("gene_imhof_pval_cutoff", cutoff);
          }
        } else if (arg["snp-wise"][0] == "perm-sum") { set_value("gene_snpwise", "imhof"); set_value("gene_snpwise_subtype", "empirical");
        } else if (arg["snp-wise"][0] == "brown1df") { set_value("gene_snpwise", "imhof"); set_value("gene_snpwise_subtype", "brown");
        } else if (arg["snp-wise"][0] == "brown2df") { set_value("gene_snpwise", "brown");                
        } else mod_error("gene-model", "snp-wise", "unknown test subtype " + Utils::quote(arg["snp-wise"][0]));
      } else set_value("gene_snpwise", "imhof");        
    }
    if (arg.check_modifier("multi", range(0,1))) { 
      set_value("gene_model", "multi");
      if (arg["multi"].size() > 0) {      
        string sub = arg["multi"][0]; if (sub == "snp-wise") sub = "snpwise";
        if (sub != "all" && sub != "snpwise" && sub != "mean") mod_error("gene-model", "multi", "unknown test subtype " + Utils::quote(arg["multi"][0]));
        if (flag_set("pval") && sub != "snpwise") flag_error(arg.name(), "when using --pval, modifier 'multi' can only be set to 'snp-wise'");            
        set_value("gene_multi", sub);
      } else if (flag_set("pval")) {set_value("gene_multi", "snpwise");
      } else set_value("gene_multi", "all");
    }

    if (arg.check_modifier("no-autoburden")) unset_marker("rare_autoburden");
    if (arg.check_modifier("multi-show-all")) {
      if (!(*this)["gene_multi"]) flag_error(arg.name(), "modifier 'multi-show-all' can only be set when using multi-model");            
      set_marker("gene_multi_show_subpval");
    }

/// --snp-wise _FLAG_
  } else if (arg.name("snp-wise") && arg.set_modifiers(0)) {
    parse_error("flag --snp-wise has been deprecated, use --gene-model with modifier 'snp-wise' instead");

/// --pval _FLAG_
  } else if (arg.name("pval") && arg.set_modifiers(1,0)) {
    mod_conflict(arg, "use", "snp-id pval"); mod_conflict(arg, "N", "ncol");
    mod_required(arg, "snp-id", "pval", true);

    set_marker("has_pval");
    set_marker("dummy_pheno");  
    set_value("pval_file", check_file(arg[0], arg.name()));
    set_value("pval_truncate_low", 1e-50);
    set_value("pval_truncate_high", 1e-5);
    unset_marker("rare_autoburden");      
    if (!flag_set("gene-model")) set_value("gene_model", "snpwise_mean");

    if (arg.check_modifier("use", 2)) {
      set_value("pval_col_snp", arg["use"][0]);
      set_value("pval_col_p", arg["use"][1]);
    }
    if (arg.check_modifier("snp-id", 1) && arg.check_modifier("pval", 1)) {
      set_value("pval_col_snp", arg["snp-id"][0]);
      set_value("pval_col_p", arg["pval"][0]);
    }    
    if ((*this)["pval_col_snp"] && (gets("pval_col_snp") == gets("pval_col_p"))) flag_error("pval", "same column specified for SNP ID and SNP p-value");

    if (arg.check_modifier("truncate", range(1,2))) {
      double val = convert_double(arg["truncate"][0], arg.name(), "truncate");
      if (val < 0) val = pow(double(10), val);
      check_range(val, drange(0.0, 0.5), arg.name(), "truncate", false);
      set_value("pval_truncate_low", val);
      
      if (arg["truncate"].size() == 2) {
        val = convert_double(arg["truncate"][1], arg.name(), "truncate");
        if (val < 0) val = pow(10.0, val);
        check_range(val, drange(0.0, 1.0), arg.name(), "truncate", false);        
        if (val > 0.5) val = 1 - val;
      }
      set_value("pval_truncate_high", val);      
    }
    
    if (arg.check_modifier("duplicate", 1)) {
      if (check_enum(arg["duplicate"][0], "first last drop error", arg.name(), "duplicate")) set_value("pval_duplicate_mode", arg["duplicate"][0]);
    } else set_value("pval_duplicate_mode", "drop");
    
    if (arg.check_modifier("N", range(1,3))) {
      long minN = geti("min_sample_size", 50);
      long size = convert_long(arg["N"][0], arg.name(), "N");
      check_range_min(size, 0, arg.name(), "N", false);
      set_value("pval_size_main", size);
      set_value("pval_size_chrX", size);
      set_value("pval_size_chrY", size);            

      if (arg["N"].size() > 1) {
        long size = convert_long(arg["N"][1], arg.name(), "N");
        check_range_min(size, 0, arg.name(), "N", false);
        set_value("pval_size_chrX", size);
        
        if (arg["N"].size() > 2) {
          long size = convert_long(arg["N"][2], arg.name(), "N");
          check_range_min(size, 0, arg.name(), "N", false);
          set_value("pval_size_chrY", size);
        } 
      }
      if (geti("pval_size_main") < minN || geti("pval_size_chrX") < minN || geti("pval_size_chrY") < minN) mod_error("pval", "N", "sample size(s) specified for modifier 'N' of flag --pval must be at least " + Utils::num_string(minN));
    }
    if (arg.check_modifier("ncol", 1)) {
      if ((*this)["pval_col_p"] && arg["ncol"][0] == gets("pval_col_p")) mod_error("pval", "ncol", "same column specified as for SNP p-value");
      if ((*this)["pval_col_snp"] && arg["ncol"][0] == gets("pval_col_snp")) mod_error("pval", "ncol", "same column specified as for SNP ID");
      set_value("pval_col_N", arg["ncol"][0]);
    }
    if (!arg.has_modifier("N") && !arg.has_modifier("ncol")) flag_error("pval", "either modifier 'N' or modifier 'ncol' must be set");

/// --debug _FLAG_
  } else if (arg.name("debug") && arg.set_modifiers(0,1)) {
    if (arg.check_modifier("dump-data", 1)) {
      pair<string,string> path = Utils::split_path(arg["dump-data"][0]);
      if (path.second == "") flag_error("debug", "'dump-data modifier' specifies a directory, should be a file prefix");
      if (path.first != "") check_dir(path.first, arg.name(), "dump-data");
      set_value("dump_data", arg["dump-data"][0]);    
    }
    if (arg.check_modifier("verbose")) set_marker("verbose_debug");
    if (arg.check_modifier("set-npar", 1, false)) {
      if (arg["set-npar"].size() == 1) set_marker(arg["set-npar"][0]);
      else {
        if (arg["set-npar"].size() % 2 == 1) flag_error("debug", "'set-npar modifier' either one or an even number of values");
        for (int i = 0; i < arg["set-npar"].size(); i += 2) {
          double val = convert_double(arg["set-npar"][i+1], arg.name(), "set-npar");
          if (floor(val) == val) set_value(arg["set-npar"][i], long(val));
          else set_value(arg["set-npar"][i], val);
        }
      }
    }
    if (arg.check_modifier("set-spar", 1, false)) {
      if (arg["set-spar"].size() == 1) set_marker(arg["set-spar"][0]);
      else {
        if (arg["set-spar"].size() % 2 == 1) flag_error("debug", "'set-spar modifier' either one or an even number of values");
        for (int i = 0; i < arg["set-spar"].size(); i += 2) set_value(arg["set-spar"][i], arg["set-spar"][i+1]);
      }
    }
    
/// --seed _FLAG_
  } else if (arg.name("seed", 1)) {
    long seed = convert_long(arg[0], arg.name());
    check_range_min(seed, 0, arg.name(), false);
    set_value("rng_seed", seed);
    
/// --version _FLAG_
  } else if (arg.name("version")) { 
    die(ExitType::NoError); ///Version flag handled elsewhere
  } else parse_error("unknown flag --" + arg.name());
}

void Settings::prologue() {
  long size = max(long(100), long(geti("gene_chunked_size", 1000)));
  size = long(geti("gene_chunked_size_override", size));
  
  double frac = min(0.5, getn("gene_chunked_fraction", 0.5));
  if ((*this)["covar_interact"]) frac = min(0.25, frac);
  if ((*this)["no_asym"]) {size = 0; frac = 0;}
  set_value("gene_chunked_fraction", frac);
  set_value("gene_chunked_size", size);

  if (gets("gene_model") == "pcreg") set_marker("gene_rsq");
  else if ((*this)["alt_pheno_all"]) {_LOG.error("configuration") << "can only analyse multiple phenotypes with PC regression model" << endl; die();}
  if (gets("gene_model") == "multi") {
    set_marker("use_residuals");
    if (!(*this)["gene_multi"]) set_value("gene_multi", "all");      
    string sub = gets("gene_multi");
    
    vector<string> multi_model;
    multi_model.push_back("snpwise_mean");
    if (sub == "all" || sub == "snpwise") multi_model.push_back("snpwise_top1");
    if (sub == "all" || sub == "mean") multi_model.push_back("pcreg");    
    if (multi_model.size() < 2) {_LOG.error("configuration") << "invalid multi-model configuration, at least two models must be set" << endl; die();}

    set_value("gene_model_multi", multi_model);
    if (!(*this)["covar_interact"] && !(*this)["snp_partition"]) set_marker("gene_multi_show_subpval");
  }

  if ((*this)["gene_model_bigdata"]) {
    set_value("gene_model_bigdata_minsize", 10000);
    set_value("gene_model_bigdata_maxstep", 25);  

    bool linreg = gets("gene_model") == "pcreg" || Utils::contains(getvs("gene_model_multi", true), "pcreg");
    bool partial = gets("gene_model_bigdata") == "partial" || gets("gene_model_bigdata") == "auto";
    if (linreg && !partial) {_LOG.error("configuration") << "--big-data cannot be set to 'full' when using principal components regression" << endl; die();}
    
    if (partial || linreg) {
      set_value("gene_model_bigdata_minsize_xtx", 10000000);      
      set_value("gene_model_bigdata_maxstep_xtx", 1);            
    } else {
      set_value("gene_model_bigdata_minsize_xtx", 20000);
      set_value("gene_model_bigdata_maxstep_xtx", 5);                  
    }
  }

  if ((*this)["fixed_perm"]) set_value("adap_perm", 0);
  if ((*this)["has_pval"]) {set_value("chrX_sex_covar", 0);}

  if ((*this)["has_rare_weights"] && !(*this)["do_rare"]) {_LOG.error("configuration") << "cannot use flag --burden-weights when burden scoring is not turned on" << endl; die();}
  if ((*this)["rare_autoburden"]) {
    set_marker("do_rare");
    set_value("rare_cutoff_maf", 0.01);
    set_value("rare_cutoff_mac", 100);
  }
  
  if ((*this)["dummy_pheno"]) set_marker("snp_skip_qc");
}

void Settings::debug_verbose() {
  FileOutput fout(gets("out_prefix").append(".debug"));
  for (map<Parameter,Value>::iterator iter = parameters.begin(); iter != parameters.end(); iter++) {
    fout << iter->first << "\t" << iter->second << endl;
  }
}

void Settings::print_version() {
  _LOG << "MAGMA version: v" << MACRO_STRING(VERNO) << " (" << VERTYPE_HACK() << ")" << endl;
}

void Settings::print_welcome() {
  _LOG << "Welcome to MAGMA v" << MACRO_STRING(VERNO) << " (" << VERTYPE_HACK() << ")" << endl;
  _LOG << "Using flags:\n";
  for (int i = 0; i < no_args; i++) _LOG << arguments[i];
  _LOG << "\nStart time is " << get_timestr(start_time) << "\n" << endl;
}

void Settings::print_remarks() {
  for (int i = 0; i < remarks.size(); i++) _LOG << "NOTE: " << remarks[i] << endl;
  if (!remarks.empty()) {remarks.clear(); _LOG << endl;}
}

void Settings::print_goodbye() {
  time_t end_time = time(0);
  _LOG << "\nEnd time is " << get_timestr(end_time);
  _LOG << " (elapsed: ";
  _LOG.time(difftime(end_time, start_time)) << ")" << endl;
}


bool Settings::flag_set(const string& flag, bool throw_error) const {
  bool out = used_flags.find(flag) != used_flags.end();
  if (throw_error && !out) {_LOG.error("configuration") << "flag --" << flag << " has not been set." << endl; die();}
  return out;
}

void Settings::required_flags(const string& flag, const string& required) {
  istringstream istr(required); // change to tokenize  
  string check;
  while (istr >> check) {
    if (!flag_set(check)) {_LOG.error("configuration") << "flag --" << flag << " requires that flag --" << check << " is also set." << endl; die();}
  }
}

void Settings::conflicting_flags(const string& flag, const string& conflicting) {
  istringstream istr(conflicting);  // change to tokenize  
  string check;
  while (istr >> check) {
    if (flag_set(check)) {_LOG.error("configuration") << "flag --" << flag << " and flag --" << check << " cannot be set simultaneously." << endl; die();} 
  }
}

void Settings::conflicting_modifiers(const string& flag, const string& mod1, const string& mod2) {
  _LOG.error("configuration") << "conflicting modifiers '" << mod1 << "' and '" << mod2 << "' for flag --" << flag << endl;
  die();
}

void Settings::required_modifier(const string& flag, const string& required, const string& mod) {
  _LOG.error("configuration") << "modifier '" << required << "' is required for flag --" << flag;
  if (mod != "") _LOG << " if modifier '" << mod << "' is set";
  _LOG << endl;
  die();
}

void Settings::required_one_modifier(const string& flag, vector<string>& req_vec, const string& mod) {string required;
  for (int i = 0; i < req_vec.size()-1; i++) {
    if (i > 0) required += ", ";
   required += "'" + req_vec[i] + "'";
  }
  required += " or '" + req_vec.back() + "'";
  
  _LOG.error("configuration") << "one of the modifiers " << required << " is required for flag --" << flag;
  if (mod != "") _LOG << " if modifier '" << mod << "' is set";
  _LOG << endl;
  die();
}

void Settings::flag_conflict(const string& flag,  const string& conflicting) {
  vector<string> flags = Utils::tokenize(flag);
  for (int i = 0; i < flags.size(); i++) {if (flag_set(flags[i])) conflicting_flags(flags[i], conflicting);}
}

void Settings::flag_conflict_set(const string& conflicting) {
  vector<string> flags = Utils::tokenize(conflicting);
  for (int i = 0; i < flags.size(); i++) {
    if (flag_set(flags[i])) {for (int j = 0; j < flags.size(); j++) {if (i != j) conflicting_flags(flags[i], flags[j]);}}
  }
}

void Settings::flag_invalid(const string& flag, const string& msg) {
  vector<string> flags = Utils::tokenize(flag);
  for (int i = 0; i < flags.size(); i++) {
    if (flag_set(flags[i])) {
      _LOG.error("configuration") << "invalid flag --" << flags[i] << ": cannot be set " << msg << endl;
      die();
    }
  }
}

void Settings::flag_requires(const string& flag, const string& required) {
  vector<string> flags = Utils::tokenize(flag);
  for (int i = 0; i < flags.size(); i++) {if (flag_set(flags[i])) required_flags(flags[i], required);}
}

void Settings::flag_requires_set(const string& required) {
  vector<string> flags = Utils::tokenize(required);
  for (int i = 0; i < flags.size(); i++) {
    if (flag_set(flags[i])) {for (int j = 0; j < flags.size(); j++) {if (i != j) required_flags(flags[i], flags[j]);}}
  }
}
 
void Settings::flag_requires(const string& flag, const string& required_flag, const string& msg) {
  if (!flag_set(required_flag)) {
    vector<string> flags = Utils::tokenize(flag);
    for (int i = 0; i < flags.size(); i++) {
      if (flag_set(flags[i])) {_LOG.error("configuration") << "flag --" << flags[i] << " requires that flag --" << required_flag << " is also set, " << msg << endl; die();} 
    }
  }
}

void Settings::parse_error(const string& msg) {_LOG.error("parsing arguments") << msg << endl; die();}
void Settings::flag_error(const string& flag, const string& msg) {_LOG.error("parsing arguments") << "for flag --" << flag << ", " << msg << endl; die();}
void Settings::mod_error(const string& flag, const string& mod, const string& msg) {_LOG.error("parsing arguments") << "for modifier '" << mod << "' of flag --" << flag << ", " << msg << endl; die();}

void Settings::add_remark(const string& msg) {remarks.push_back(msg);}
void Settings::flag_remark(const string& flag, const string& msg) {add_remark(string("for flag --") + flag + ", " + msg);}
void Settings::mod_remark(const string& flag, const string& mod, const string& msg) {add_remark(string("for modifier '") + mod + "' of flag --" + flag + ", " + msg);}


void Settings::mod_conflict(Flag& flag, const string& mod_str) {
  vector<string> mods = Utils::tokenize(mod_str);
  for (int i = 0; i < mods.size(); i++) {for (int j = i+1; j < mods.size(); j++) mod_conflict(flag, mods[i], mods[j]);}
}

void Settings::mod_conflict(Flag& flag, const string& mod_str1, const string& mod_str2) {
  vector<string> mods1 = Utils::tokenize(mod_str1), mods2 = Utils::tokenize(mod_str2);
  for (int i = 0; i < mods1.size(); i++) { for (int j = 0; j < mods2.size(); j++) {
    if (flag.has_modifier(mods1[i]) && flag.has_modifier(mods2[j]) && mods1[i] != mods2[j]) conflicting_modifiers(flag.name(), mods1[i], mods2[j]);
  } }
}

void Settings::mod_required(Flag& flag, const string& req_str) {
  vector<string> mods = Utils::tokenize(req_str);
  for (int i = 0; i < mods.size(); i++) {
    if (!flag.has_modifier(mods[i])) required_modifier(flag.name(), mods[i], "");
  }
}

void Settings::mod_required(Flag& flag, const string& mod_str, const string& required, bool symmetric) {
  vector<string> mods = Utils::tokenize(mod_str);
  for (int i = 0; i < mods.size(); i++) {
    if (flag.has_modifier(mods[i]) && !flag.has_modifier(required)) required_modifier(flag.name(), required, mods[i]);
    if (symmetric && flag.has_modifier(required) && !flag.has_modifier(mods[i])) required_modifier(flag.name(), mods[i], required);  
  }
}

void Settings::mod_one_required(Flag& flag, const string& mod_str, const string& req_str) {
  if (mod_str != "" && !flag.has_modifier(mod_str)) return;
  vector<string> mods = Utils::tokenize(req_str); 
  for (int i = 0; i < mods.size(); i++) if (flag.has_modifier(mods[i])) return;
  required_one_modifier(flag.name(), mods, mod_str);  
}

string Settings::get_timestr(const time_t& time) {
  char time_str[80];
  struct tm loc_time = *localtime(&time);
  strftime(time_str, sizeof(time_str), "%X, %A %d %b %Y", &loc_time);
  return string(time_str);
}

bool Settings::check_enum(const string& argval, const string& known_str) {
  vector<string> known = Utils::tokenize(known_str); bool hit = false;
  for (int i = 0; i < known.size(); i++) {if (argval == known[i]) {hit = true; break;}}
  return hit;
}

bool Settings::check_enum(const string& argval, const string& known_str, const string& flag, const string& modifier) {
  bool hit = check_enum(argval, known_str);
  if (!hit) mod_error(flag, modifier, "unknown value " + Utils::quote(argval));
  return hit;
}

bool Settings::check_bool(const string& argval, const string& flag, const string& modifier) {
  bool valid = argval.empty() || check_enum(argval, "0 1");
  if (!valid) mod_error(flag, modifier, "unknown value " + Utils::quote(argval));
  return argval != "0";
}


void Settings::check_datatype(const string& req_type, const string& flag, const string& modifier) {
  string used_type = gets("data_type", "unknown");
  if (req_type != used_type) {
    string req_label = (req_type == "binary_geno") ? "genotype" : req_type;
    if (modifier != "") flag_error(flag, string("modifier '") + modifier + "' can only be set when using " + req_label + " data");
    else parse_error(string("flag --") + flag + " can only be set when using " + req_label + " data");
  }
}  

long Settings::convert_long(const string& argval, const string& flag, const string& modifier) {
  long out;
  if (!Utils::convert_num(argval, out)) {
    _LOG.error("parsing arguments") << "value for ";
    if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
    _LOG << "--" << flag << " (" << argval << ") is not a valid integer." << endl;
    die();
  }  
  return out;
}

double Settings::convert_double(const string& argval, const string& flag, const string& modifier) {
  double out;
  if (!Utils::convert_num(argval, out)) {
    _LOG.error("parsing arguments") << "value for ";
    if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
    _LOG << "--" << flag << " (" << argval << ") is not a valid number." << endl;
    die();
  }  
  return out;
}

vector<string> Settings::read_list(const string& filename, bool by_line, const string& flag, const string& modifier) {
  check_file(filename, flag, modifier);
  TextInput fin(filename, 0); vector<string> list;  
  fin.set_error("reading file '" + filename + "'");  
  if (!by_line) {
    while (fin.process_line()) {
      while(fin.read_value()) list.push_back(fin.value);
    }
  } else {while (fin.process_line()) list.push_back(fin.curr_line);}       
  return list;
}
  
string Settings::check_file(const string& filename, const string& flag, const string& modifier, bool check_win) {
  if (!Utils::is_file(filename)) {
    if (check_win) {
      string winname = win_txt(filename);
      if (winname != filename && Utils::is_file(winname)) return winname; 
    }
    _LOG.error("parsing arguments") << "file '" << filename << "' specified in ";
    if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
    _LOG << "--" << flag << " does not exist or is not a file." << endl;
    die();
  }
  return filename;
}

const string& Settings::check_dir(const string& dirname, const string& flag, const string& modifier) {
  if (!Utils::is_dir(dirname)) {
    _LOG.error("parsing arguments") << "directory '" << dirname << "' specified in ";
    if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
    _LOG << "--" << flag << " does not exist or is not a directory." << endl;
    die();
  }
  return dirname;
}

vector<string> Settings::expand_prefix(const string& prefix, const string& suffix, const string& flag, long min_files) {
  vector<string> files;
  pair<string,string> path = Utils::split_path(prefix);
  long pre_len = path.second.size(), suf_len = suffix.size(), min_length = pre_len + suf_len;
  
  if (path.first != "") check_dir(path.first, flag);
  else path.first = ".";

  DIR *dir; struct dirent *file; string curr; 
  dir = opendir(path.first.c_str());
  while ( (file = readdir(dir)) ) {
    if (strlen(file->d_name) < min_length) continue;
    curr.assign(file->d_name);
    if (curr.compare(0, pre_len, path.second) != 0) continue;
    if (curr.compare(curr.size()-suf_len, suf_len, suffix) != 0) continue;
    if (path.first != ".") curr = string(path.first).append(curr);
    files.push_back(curr);    
  }
  closedir(dir);
  
  if (files.size() < min_files) {
    _LOG.error("parsing arguments") << "for flag --" << flag << ", found ";
    if (files.size() == 0) _LOG << "no files";
    else if (files.size() == 1) _LOG << "only 1 file";
    else _LOG << "only " << files.size() << " files";
    _LOG << " matching pattern '" << prefix << "*" << suffix << "' (at least " << min_files << " required)" << endl;
    die();
  }  
  return files;
}

void Settings::check_suffix(const string& filename, const string& suffix, const string& flag, const string& modifier) {
  if (suffix.size() == filename.size() || !Utils::ends_with(filename, suffix)) {
    _LOG.error("parsing arguments") << "file '" << filename << "' specified in ";
    if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
    _LOG << "--" << flag << " does not have a '" << suffix << "' suffix." << endl;
    die();
  }
}

void Settings::parse_metaweights(const string& filename, vector<string>& meta_genes, const string& type, bool prefix) {
  string modifier = string(type).append("-file");
  string suffix = type == "genes" ? win_txt(".genes.out") : ".genes.raw";

  check_file(filename, "meta", modifier);

  TextInput fin(filename);
  fin.set_error(string("reading file '").append(filename).append("'"));
  
  double value;
  while (fin.process_line()) {
    if (prefix) fin[0].append(suffix);
    else check_suffix(fin[0], suffix, "meta", modifier);
    check_file(fin[0], "meta", modifier);

    meta_genes.push_back(fin[0]);
    Triple<double> curr(0);
    for (int i = 0; i < 3 && fin.read_value(); i++) {
      if (!Utils::convert_num(fin.value, value) || value < 0) fin.line_error("value is negative or not a number");
      if (i == 0 && value == 0) fin.line_error("first value must be greater than zero");
      for (int j = i; j < 3; j++) curr[j] = value;
    }
    meta_weights.push_back(curr);
  }

  if (meta_genes.size() < 2) {_LOG.error("parsing arguments") << "file '" << filename << "' contains only one entry; at least two are required for meta-analysis" << endl; die();}
}

void Settings::parse_metacorrs(const string& filename, int file_count) {
  check_file(filename, "meta", "correlations");
  meta_corrs.set_size(file_count, file_count, true);

  TextInput fin(filename, 0);
  fin.set_error(string("reading file '").append(filename).append("'"));
  
  int curr_line; double value;
  for (curr_line = 0; fin.process_line(); curr_line++) {
    if (curr_line >= file_count) {
      if (fin.read_value()) fin.error("number of rows in matrix is larger than number of files specified for meta-analysis"); 
      else continue;
    }
    for (int i = 0; i <= curr_line; i++) {
      if (!fin.read_value()) {
        if (i == 0) {curr_line--; goto LOOP_END;}
        else fin.line_error("not enough values");
      }
      if (!Utils::convert_num(fin.value, value)) fin.line_error("value is not a number");
      if (i == curr_line && value == 0) fin.line_error("value on diagonal is zero");
      meta_corrs(curr_line,i) = value;
    } 
  }
  LOOP_END:;
  if (curr_line < file_count) fin.error("number of rows in matrix is smaller than number of files specified for meta-analysis"); 

  for (int i = 0; i < file_count; i++) { 
    double sd = sqrt(meta_corrs(i,i)); meta_corrs(i,i) = 1;
    for (int j = i+1; j < file_count; j++) {
      meta_corrs(j,i) /= sd * sqrt(meta_corrs(j,j));
      if (meta_corrs(j,i) > 1) {
        if (meta_corrs(j,i) > 1.01) fin.error("after conversion to correlation matrix, found value greater than 1"); 
        else meta_corrs(j,i) = 1;
      }
      meta_corrs(j,i) = meta_corrs(j,i)*meta_corrs(j,i);
      meta_corrs(i,j) = meta_corrs(j,i);      
    }
  }
}



/***************************************************
* Flag subclass                                    *
***************************************************/


bool Settings::Flag::check_size(int required, bool exact) {
  if ((exact && size()==required) || (!exact && size() >= required)) return true;
  _LOG.error("parsing arguments");
  if (required == 0) {
    if (is_modifier) _LOG << "modifier '" << flag << "' for flag --" << parent_flag << " should have no values";
    else _LOG << "flag --" << flag << " should have no values";
  } else {
    if (!exact) _LOG << "at least ";
    _LOG << required << " " << (required == 1 ? "value" : "values") << " required for ";
    if (is_modifier) _LOG << "modifier '" << flag << "' of flag --" << parent_flag;
    else _LOG << "flag --" << flag;
  }
  _LOG << " (found " << size() << ")." << endl;
  die();
  return true;
}

bool Settings::Flag::check_size(pair<int,int> required) {
  if (size() < required.first || size() > required.second) {
    _LOG.error("parsing arguments");
    if (is_modifier) _LOG << "modifier '" << flag << "' for ";
    _LOG << "flag --" << parent_flag << " should have ";
    if (size() < required.first) _LOG << "at least " << required.first << " " << (required.first == 1 ? "value" : "values");
    else _LOG << "at most " << required.second << " " << (required.second == 1 ? "value" : "values");
    _LOG << " (found " << size() << ")." << endl;
    die();
  }
  return true;
}

bool Settings::Flag::set_modifiers(int offset, int min_required) {
  if (processed) return true; 
  processed = true;

  check_size(offset, false);
  Flag mod; string mod_name; int is_pos;
  
  while (offset < size()) {
    mod_name = values[offset++];
    if (mod_name[0] != '=') {
      is_pos = mod_name.find('=');
      if (is_pos == string::npos) {    
        mod = Flag(mod_name, flag);
        if (offset >= size() || values[offset][0] != '=') {modifiers[mod_name] = mod; continue;}
        else mod_name = values[offset++];
      } else {   
        mod = Flag(mod_name.substr(0,is_pos), flag);
        mod_name.erase(0, is_pos);
      }
    }
    if ((mod_name.size() == 1 && offset == size()) || mod.name() == "" || mod.size() > 0) {_LOG.error("parsing arguments") << "incomplete modifier specified for flag --" << flag << endl; die();} 

    if (mod_name.size() == 1) mod_name = values[offset++];
    else mod_name.erase(0,1);
    
    while (mod_name[mod_name.size()-1] == ',') {
      if (offset == size()) {_LOG.error("parsing arguments") << "modifier for flag --" << flag << " ends with ',' character" << endl; die();} 
      mod_name.append(values[offset++]);
    }
    if (mod_name.find('=') != string::npos) {_LOG.error("parsing arguments") << "modifier for flag --" << flag << " contains more than one '=' sign" << endl; die();} 
    
    for (int begin = 0, end = mod_name.find(','); ; begin = end+1, end = mod_name.find(',', begin)) {
      mod.add_value(mod_name.substr(begin, end-begin));
      if (end == string::npos) break;
    }    
    if (modifiers.find(mod.name()) != modifiers.end()) {_LOG.error("parsing arguments") << "modifier '" << mod.name() << "' for flag --" << flag << " occurs more than once" << endl; die();} 
    modifiers[mod.name()] = mod;          
  }
  if (modifiers.size() < min_required) {
    _LOG.error("parsing arguments");
    _LOG << "at least ";
    _LOG << min_required << " " << (min_required == 1 ? "modifier" : "modifiers") << " required for ";
    _LOG << "flag --" << flag;
    _LOG << " (found " << size() << ")." << endl;
    die();
  }
  return true;
}

bool Settings::Flag::name(const string& check, int required, bool exact) {
  if (!name(check)) return false;
  check_size(required, exact);
  return true;
}

bool Settings::Flag::name(const string& check, pair<int,int> required) {
  if (!name(check)) return false;
  check_size(required);
  return true;
}


bool Settings::Flag::has_modifier(const string& name, bool process) {
  bool out = modifiers.find(name) != modifiers.end();
  if (out && process) modifiers[name].processed = true;
  return out;  
}

bool Settings::Flag::check_modifier(const string& name, int required, bool exact) {
  if (!has_modifier(name, true)) return false;
  modifiers[name].check_size(required, exact);
  return true;
}

bool Settings::Flag::check_modifier(const string& name, pair<int,int> required) {
  if (!has_modifier(name, true)) return false;
  modifiers[name].check_size(required);
  return true;
}

set<string> Settings::Flag::check_modifiers(const string& name_list, int required, bool exact) {
  vector<string> mods = Utils::tokenize(name_list); set<string> used;
  for (int i = 0; i < mods.size(); i++) {if (check_modifier(mods[i], required, exact)) used.insert(mods[i]);}
  return used;
}

set<string> Settings::Flag::check_modifiers(const string& name_list, pair<int,int> required) {
  vector<string> mods = Utils::tokenize(name_list); set<string> used;
  for (int i = 0; i < mods.size(); i++) {if (check_modifier(mods[i], required)) used.insert(mods[i]);}
  return used;
}


void Settings::Flag::check_unused_modifiers() {
  for (map<string,Flag>::iterator iter = modifiers.begin(); iter != modifiers.end(); ++iter) {
    if (!iter->second.processed) {_LOG.error("configuration") << "unknown modifier '" << iter->second.flag << "' for flag --" << flag << endl; die();} 
  }
}



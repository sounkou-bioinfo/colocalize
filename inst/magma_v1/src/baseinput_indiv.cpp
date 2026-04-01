/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "baseinput.h"
#include "input.h"


using namespace std;


void BaseInput::pheno_process(int nobs) {
  indiv_total = indiv_obs = nobs;
  indiv_missing.set_size(indiv_total,6); indiv_missing.assign_value(true); 
  indiv_missing.assign_value_col(false,1); indiv_missing.assign_value_col(false,5); 
}

void BaseInput::pheno_process(int nobs, vector<int>& missing) {
  pheno_process(nobs);
  if (!indiv_pheno.empty()) {
    indiv_missing.assign_value_col(false,0);
    for (int i = 0; i < missing.size(); i++) indiv_missing(i,0) = 1;
    check_binary();
  }
}

void BaseInput::pedigree_process(vector<pair<string, string> >& parents) {
  if (parents.size() != indiv_total) {_LOG.error("processing pedigree") << "pedigree information does not match number of individuals in data" << endl; die();} 
  
  indiv_pedigree.set_size(indiv_total, 2); indiv_pedigree.assign_value(-1);
  boost::unordered_map<string,long>::iterator found;  
  for (int i = 0; i < indiv_total; i++) {
    if ((found = indiv_map.find(parents[i].first)) != indiv_map.end()) indiv_pedigree(i,0) = found->second;
    if ((found = indiv_map.find(parents[i].second)) != indiv_map.end()) indiv_pedigree(i,1) = found->second;    
  }
}

int BaseInput::check_binary() {
  int nobs = 0; binary_pheno = true;
  
  for (int i = 0; i < indiv_total; i++) {
    if (!indiv_missing(i)) {
      nobs++;
      if (binary_pheno) {
        float& ph = indiv_pheno[i];
        if (!(ph == -9 || ph == 0 || ph == 1 || ph == 2)) binary_pheno = false;
      }
    }
  }
  
  if (binary_pheno) {
    _LOG << "\tphenotype is binary, values of -9 and 0 will be considered missing values" << endl;
    for (int i = 0; i < indiv_total; i++) {
      float& ph = indiv_pheno[i];
      if (!indiv_missing(i) && (ph == -9 || ph == 0)) {
        indiv_missing(i) = 1;
        nobs--;
      }
    }
    if (!indiv_pheno_multi.empty()) {
      _LOG << "\tassuming additional phenotypes are binary as well" << endl;    
      for (int c = 0; c < indiv_pheno_multi.ncol(); c++) {
        float* curr = indiv_pheno_multi[c];
        for (int i = 0; i < indiv_total; i++) {
          if (!indiv_missing(i)) {
            if (curr[i] == -9 || curr[i] == 0) {_LOG.error("processing phenotypes") << "missing values are not allowed for additional phenotypes" << endl; die();}  
            if (curr[i] !=  1 && curr[i] != 2) {_LOG.error("processing phenotypes") << "additional phenotypes are not all binary" << endl; die();}              
          }
        }
      }
    }
  }

  return nobs;
}

void BaseInput::set_dummy_pheno() {
  boost::normal_distribution<float> distn(0,1);
  boost::variate_generator<boost::mt19937, boost::normal_distribution<float> > norm_gen(_RNG, distn);

  indiv_pheno.reserve(indiv_total);
  for (long i = 0; i < indiv_total; i++) indiv_pheno.push_back(norm_gen());
  indiv_missing.assign_value_col(false,0);
  binary_pheno = false;  
}

void BaseInput::set_const_pheno(double value) {
  indiv_pheno.assign(indiv_total, value);
  indiv_missing.assign_value_col(false,0);
  binary_pheno = (value == 1) || (value == 2);
}
 
void BaseInput::filter_indiv(const string& filename, const string& mode) {
  _LOG << "Filtering individuals listed in file " << filename << " (mode = " << mode << ")...";

  TextInput fin(filename, 2);
  fin.set_error("reading individuals filter file");

  int match_count = 0, filter_count = 0; boost::unordered_map<string,long>::iterator found;
  vector<short> matched(indiv_total, false);
  while (fin.process_line()) {
    string cid = fin[0].append(" ").append(fin[1]);
    if ((found = indiv_map.find(cid)) == indiv_map.end()) continue;
    if (!matched[found->second]) match_count++;
    matched[found->second] = true;
  }

  short drop_value = (mode == "include") ? false : true;
  for (int i = 0; i < indiv_total; i++) {
    if ((matched[i] == drop_value)) {indiv_missing(i,0) = true; filter_count++;}
  }
  _LOG << " removing " << Utils::plural(filter_count, "individual") << endl;
}


void BaseInput::filter_missing(bool pval_only) {
  bool sex_covar = settings["covar_use_sex"], use_covar = (no_covar > 0) || sex_covar;
  int nobs = 0;

  for (long i = 0; i < indiv_total; i++) {
    if (sex_covar && indiv_gender[i] == 0) indiv_missing(i,1) = 1;
    if (use_covar && indiv_missing(i,1)) indiv_missing(i) = 1;
    if (!indiv_missing(i)) {
      nobs++;
      if (indiv_gender[i] != 0) {
        indiv_missing(i,2) = 0;
        if (indiv_gender[i] == 1) indiv_missing(i,3) = 0;
        else indiv_missing(i,4) = 0;
      }
    }
  }

  if (!pval_only) _LOG << "Filtering phenotype/covariate missing values... " << nobs << " individuals remaining" << endl;
  if (nobs < settings.geti("min_sample_size")) {_LOG.error("processing individuals") << (pval_only ? "reference" : "remaining") << " sample is too small (minimum = " << settings.geti("min_sample_size") << " individuals)" << endl; die();} 
  indiv_obs = nobs;
} 

void BaseInput::filter_trios() {
  if (!binary_pheno) {_LOG.error("processing pedigree") << "analysis of trio data only available for binary phenotype" << endl; die();} 

  bool check_case = pheno_type != ph_Trio;
  trio_map.set_size(indiv_total, 3); 
  int index = 0; vector<short> used(indiv_total, false);
  for (int i = 0; i < indiv_total; i++) {
    long &p1 = indiv_pedigree(i,0), &p2 = indiv_pedigree(i,1);
    if (p1 < 0 || p2 < 0) continue;
    if (indiv_missing(i) || indiv_missing(p1) || indiv_missing(p2)) continue;
    if (indiv_pheno[p1] != 1 || indiv_pheno[p2] != 1) continue;
    if (check_case && indiv_pheno[i] != 2) continue;
    if (used[i] || used[p1] || used[p2]) continue;
    
    trio_map(index,0) = i; trio_map(index,1) = p1; trio_map(index,2) = p2; 
    used[i] = used[p1] = used[p2] = true;
    index++;        
  }
  
  int ncol = indiv_missing.ncol() - 1;
  for (int i = 0; i < indiv_total; i++) {
    if (!used[i]) {
      indiv_missing(i) = 1;
      for (int c = 2; c < ncol; c++) indiv_missing(i,c) = 1;
    }
  }
  
  trio_map.shrink_rows(index);
  _LOG << "Mapped " << 3*index << " individuals to " << index << " valid trios" << endl;
  if (index < settings.geti("min_sample_size")) {_LOG.error("processing pedigree") << "sample is too small (minimum = " << settings.geti("min_sample_size") << " trios)" << endl; die();} 
  indiv_obs = 3*index;
}


 
void BaseInput::load_altpheno(const string& filename) {
  _LOG << "Reading phenotype file " << filename << "... " << endl;

  short offset = 2;
  TextInput fin(filename);
  fin.set_error("reading phenotype file");
  long tot_var = fin.read_header(offset);
  bool do_multi = settings["alt_pheno_all"];

  string pheno_name; string* pheno_raw = 0;
  _LOG << "\tdetected " << tot_var << " variables in file" << endl;

  if (!do_multi) {
    pheno_name = fin.name("pheno"); pheno_raw = fin.set_var("pheno", settings.gets("alt_pheno"));
    _LOG << "\tusing variable: " << fin.name("pheno") << endl;
  } else {
    pheno_name = fin.name(offset); pheno_raw = &fin[offset];
    _LOG << "\tusing all variables" << endl;
  }

  string cid_data; boost::unordered_map<string,long>::iterator found;
  long individ; float pheno;

  indiv_pheno.assign(indiv_total, 0);
  if (do_multi) indiv_pheno_multi.set_size(indiv_total, tot_var - 1, true);
  while (fin.process_line()) {
    string& cid = has_fid ? cid_data : fin[0];
    if (has_fid) cid = fin[0].append(" ").append(fin[1]);

    if ((found = indiv_map.find(cid)) == indiv_map.end()) continue;
    individ = found->second;

    if (*pheno_raw == "NA") continue;
    if (!fin.read_num(pheno, *pheno_raw)) fin.line_error("non-numeric value for variable " + pheno_name);
    if (fin.read_value()) fin.line_error("too many values");

    indiv_missing(individ) = 0;
    indiv_pheno[individ] = pheno;

    if (do_multi) {
      int from = offset+1, to = offset+tot_var;
      for (int i = from; i < to; i++) {
        if (fin[i] == "NA") fin.line_error("missing values are not allowed for additional phenotypes");

        if (!fin.read_num(pheno, *pheno_raw)) fin.line_error("non-numeric value for variable " + fin.name("pheno"));
        if (!fin.read_num(indiv_pheno_multi(individ,i-from), fin[i])) fin.line_error("non-numeric value for variable " + fin.name(i));
      }
    } 
  }
  int nobs = check_binary();

  _LOG << "\tread non-missing phenotype for " << nobs << " individuals in SNP data file" << endl;
}

void BaseInput::load_covar() { 
  if (settings["covar_file"]) {
    int nobs = load_covarfile(settings.gets("covar_file"));
    if (interact_type == vt_Regular) _LOG << "\tsetting variable " << covar_names[interact_index] << " as interactor" << endl;
    _LOG << "\tread non-missing covariates for " << nobs << " individuals in SNP data file" << endl;
  }
  if (settings["covar_use_sex"]) {
    _LOG << "Setting sex as covariate";
    if (interact_type == vt_Gender) _LOG << " and interactor";
    _LOG << endl;
  }
  if (settings["covar_snps"]) load_covarsnps(settings.getvs("covar_snps", true));  
}

long BaseInput::load_covarfile(const string& filename) {
  _LOG << "Reading covariate file " << filename << "... " << endl;

  short offset = 2;
  TextInput fin(filename);
  fin.set_error("reading covariate file");

  long tot_var = fin.read_header(offset);
  
  fin.parse_used(settings.get_covar(true));
  fin.set_subset(settings.get_covar(true), settings.get_covar(false));
  no_covar = fin.set_map();
  
  fin.load_names(covar_names);

  if (no_covar == 0) fin.error("no variables selected");
  _LOG << "\tdetected " << tot_var << " variables in file";
  if (no_covar == tot_var) _LOG << " (using all)" << endl;
  else {
    _LOG << endl << "\tusing variable" << (no_covar > 1 ? "s" : "") << ": " << fin.print_names() << endl;
  }

  if (interact_type == vt_Regular) interact_index = fin.var_index(settings.gets("covar_interact"), true);

  string cid_data; boost::unordered_map<string,long>::iterator found;
  long index = 0, individ;

  indiv_missing.assign_value_col(true,1);  
  covar.set_size(indiv_total, no_covar, true);
  while (fin.process_line()) {
    string& cid = has_fid ? cid_data : fin[0];
    if (has_fid) cid = fin[0].append(" ").append(fin[1]);

    if ((found = indiv_map.find(cid)) == indiv_map.end()) continue;
    individ = found->second;

    for (long i = 0; i < no_covar; i++) {
      if (fin(i) == "NA") goto LOOP_END;
      if (!fin.read_num(covar(individ,i), fin(i))) fin.line_error("non-numeric value for variable " + fin.name(i,true));
    }
    if (fin.read_value()) fin.line_error("too many values");

    indiv_missing(individ,1) = 0;
    index++;
    
    LOOP_END:;
  }
  return index;
}

void BaseInput::load_covarsnps(vector<string> snps) { 
  _LOG << "Adding SNPs as covariates" << endl; 
  if (snps.empty()) {_LOG.error("loading covariate SNPs") << "no SNPs have been specified" << endl; die();} 

  if (snps.size() <= 25) {
    _LOG << "\tusing SNPs: ";
    for (int i = 0; i < snps.size(); i++) {if (i > 0) _LOG << ", "; _LOG << snps[i];}
    _LOG << endl;
  } else _LOG << "\tusing " << snps.size() << " SNPs" << endl;

  int id; set<long> snp_ids; vector<string> keep, drop;
  for (int i = 0; i < snps.size(); i++) {
    if ((id = snp_info[snps[i]]) >= 0) {snp_ids.insert(id); keep.push_back(snps[i]);}
    else drop.push_back(snps[i]);
  }  
  string interactor = (interact_type == vt_SNP) ? settings.gets("covar_interact_snp") : "";

  if (snp_ids.empty()) {_LOG.error("loading covariate SNPs") << "none of the specified SNPs are present in the data" << endl; die();} 
  if (interactor != "" && !Utils::contains(keep, interactor)) {_LOG.error("loading covariate SNPs") << "interactor SNP is not present in the data" << endl; die();} 
  if (!drop.empty()) {
    _LOG << "\tWARNING: dropped SNPs not present data (";
    for (int i = 0; i < drop.size(); i++) {if (i > 0) _LOG << ", "; _LOG << drop[i];}
    _LOG << ")" << endl;
  }

  no_covar_snps = snp_ids.size();
  covar_snps.set_size(indiv_total, no_covar_snps);
  load_snpdata(snp_ids, covar_snps.begin(), sm_All);
  
  vector<int> drop_index; 
  float miss = misscode(), thresh = MathUtils::clamp(settings.getn("snp_max_miss")*2, 0, settings.getn("snp_max_miss_max"));
  for (int i = 0; i < no_covar_snps; i++) {int nmiss = 0; float sum = 0;
    for (float *curr = covar_snps[i], *end = covar_snps[i+1]; curr < end; curr++) {
      if (*curr == miss) nmiss++; 
      else sum += *curr;
    }
    if (nmiss < indiv_total*thresh) {
      float mean = sum / (indiv_total-nmiss);
      covar_names.push_back(keep[i]);      
      for (float *curr = covar_snps[i], *end = covar_snps[i+1]; curr < end; curr++) {
        if (*curr == miss) *curr = 0;
        else *curr -= mean;
      }             
    } else drop_index.push_back(i);    
  }
  no_covar_snps -= drop_index.size();

  if (interactor != "") {
    for (int i = 0; i < no_covar_snps; i++) {if (covar_names[no_covar+i] == interactor) {interact_index = i; break;}}
    if (interact_index < 0) {_LOG.error("loading covariate SNPs") << "missingness for interactor SNP exceeds missingness threshold" << endl; die();} 
  }

  if (drop_index.size() > 0) {  
    if (drop_index.size() >= no_covar_snps) {_LOG.error("loading covariate SNPs") << "missingness for all specified covariate SNPs exceeded missingness threshold" << endl; die();} 
    else {
      covar_snps.drop_cols(drop_index); 
      _LOG << "\tWARNING: dropped SNPs with missingness greater than missingness threshold (";      
      for (int i = 0; i < drop_index.size(); i++) {if (i > 0) _LOG << ", "; _LOG << keep[drop_index[i]];}
      _LOG << ")" << endl;   
    }
  }
} 


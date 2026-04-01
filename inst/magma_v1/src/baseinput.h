/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef BASEINPUT_H
#define BASEINPUT_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <boost/unordered_map.hpp>

#include "parse.h"
#include "output.h"
#include "utils.h"
#include "geneutils.h"
#include "statutils.h"
#include "mathutils.h"
#include "annotation.h"
#include "gene.h"
#include "data_snpdata.h"
#include "data_genedata.h"

using namespace std;
class GeneStats;

class BaseInput {
friend class GeneStats; friend class DefaultStats; friend class TransmissionStats;
public:
  SNPData snp_info;
  GeneData& gene_info;

  enum VarType {vt_Unknown, vt_Regular, vt_Gender, vt_SNP, vt_All};
  enum SubsetMode {sm_All, sm_AllKnown, sm_KnownGender, sm_MaleOnly, sm_FemaleOnly}; 
protected:
  Settings& settings;

///geno/pheno/covar data
  boost::unordered_map<string,long> indiv_map;         ///Maps FID+IID onto data file indiv. index/line number

  vector<short> indiv_gender;
  vector<float> indiv_pheno;
  Buffer<float> indiv_pheno_multi;
  
  Buffer<short> indiv_missing;  ///By column: 0 = missing overall, 1 = missing covariates; 2 = known gender only, 3 = males only, 4 = females only 
  Buffer<long> indiv_pedigree;
  Buffer<long> trio_map;
  
  Buffer<float> covar;
  Buffer<float> covar_snps;
  vector<string> covar_names;
  
  bool interact;
  VarType interact_type;
  int interact_index;

  long indiv_total;                   ///Total number of individuals in data file (including missings)
  long indiv_obs;                     ///Total number of individuals, after filtering
  long no_covar;
  long no_covar_snps;
  bool human_chr;
  bool binary_pheno;
  bool has_fid;
  bool has_ped_info;
  bool snp_drop_report;
    
  enum Pedigree {pd_None, pd_Family, pd_Trio} pedigree;
  enum PhenoType {ph_FromData, ph_FromFile, ph_Dummy, ph_Trio, ph_Constant} pheno_type;

  bool has_pval, has_weights;
  bool has_snpN;
  Buffer<double> curr_pval;
  Buffer<int> curr_snpN;

  void init();

  void pheno_process(int nobs);
  virtual void pheno_process(int nobs, vector<int>& missing);
  void pedigree_process(vector<pair<string, string> >& parents);  
  int check_binary();
  void load_altpheno(const string& filename);
  void filter_missing(bool pval_only=false);
  void filter_trios();
  void set_const_pheno(double value=1);
  void set_dummy_pheno(); 
 
  void load_covar();
  long load_covarfile(const string& filename);
  void load_covarsnps(vector<string> snps);
 
  void load_synonyms(const string& filename, bool auto_detect=false);
  void filter_snps(const string& filename, const string& mode);
  void filter_indiv(const string& filename, const string& mode);
  void load_pvalfile();
  void load_burden_weights();

  void load_genedef(const string& filename);
  void batch_split();   
  int batch_genes(vector<short>& keep, int curr_batch, int no_batches);
  int batch_all(vector<short>& keep, int curr_batch, int no_batches);  
  
  int check_duplicates(const string& filename, string header);
  
  void set_pval_buffers(LoadedGene* gene, const set<long>& snps);  

  void set_buffers(LoadedGene* gene, int nsnps, float*& data_buff, double*& pval_buff, int*& n_buff, long*& id_buff, float*& weight_buff);
  void load_snps_core(const set<long>& snps, SubsetMode filter_id, float*& data_buff, double*& pval_buff, int*& n_buff, long*& id_buff, float*& weight_buff);
  virtual unsigned long long load_snpdata(const set<long>& snps, float* buffer, SubsetMode filter_id) = 0;
  virtual long allele_count(const set<long>& snps, Buffer<int>& counts, SubsetMode filter_id) = 0; 

public:
  ///All derived classes must read/prepare their primary data files on construction
  BaseInput(Settings& s) : gene_info(*(new GeneData())), settings(s), interact_index(-1), indiv_total(0), no_covar(0), no_covar_snps(0), human_chr(true), binary_pheno(true), has_ped_info(false), snp_drop_report(false), pedigree(pd_None), pheno_type(ph_FromData) {
    DataFrame::set_owner(&gene_info, this);
                                        
    interact = settings["covar_interact"];
    interact_type = interact ? (!settings["covar_interact_sex"] ? (!settings["covar_interact_snp"] ? vt_Regular : vt_SNP) : vt_Gender) : vt_Unknown; 
    has_pval = settings["has_pval"];
    has_weights = settings["has_rare_weights"];
    has_snpN = settings["pval_col_N"];
    has_fid = !settings["data_no_fid"];

    if (settings.gets("gene_model") == "tdt" && settings.gets("gene_tdt", "") == "ignore_p=-eno") pheno_type = ph_Trio;
    else if (has_pval) pheno_type = ph_Dummy;
    else if (settings["alt_pheno"]) pheno_type = ph_FromFile;

    if (settings.gets("gene_model") == "tdt") pedigree = pd_Trio;
  }
  virtual ~BaseInput() {DataFrame::delete_owned(&gene_info, this);}
 
  GeneData* eject_data() {DataFrame::unset_owner(&gene_info, this); return &gene_info;}  
  virtual float misscode() = 0;
  short* get_filter(SubsetMode selection);   
  LoadedGene* load_snps(int gene_id, const set<long>& snps, LoadedGene* gene, SubsetMode filter_id);
  
  string cov_name(int index, VarType type) {
    if (type == vt_Gender) return "sex";
    if (type == vt_SNP && index < no_covar_snps) return covar_names[no_covar+index];
    if (type == vt_Regular && index < no_covar) return covar_names[index];
    return "_UNKNOWN_";
  }

  int ncov(VarType type=vt_All) {
    switch(type) {
      case vt_Regular: return no_covar;
      case vt_Gender: return settings["covar_use_sex"];
      case vt_SNP: return no_covar_snps;
      case vt_All: default: return no_covar+no_covar_snps+settings["covar_use_sex"];
    }
  }
  bool is_human() {return human_chr;}
  bool has_interactor() {return interact;}
  int samp_size(bool obs=true) {return obs ? indiv_obs : indiv_total;}
  int get_npheno() {return (!indiv_pheno_multi.empty() ? (indiv_pheno_multi.ncol() + 1) : 1);}
};

#include "baseinput_plinkinput.h"
 
#endif /** BASEINPUT_H */
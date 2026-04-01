/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef BASEINPUT_PLINKINPUT_H
#define BASEINPUT_PLINKINPUT_H

#include "baseinput.h"

class PlinkInput : public BaseInput {
private:
  ifstream geno_file;                 ///Handle to PLINK .bed file
  unsigned long bed_size;
  float geno_index[256][4];       
  short geno_count[256];       
  short geno_miss[256];       
  char* raw_snp_buffer;
  unsigned char* filter_buffer; 

  unsigned long long block_count;
  long filter_observed;
  SubsetMode active_filter;   
  
  using BaseInput::pheno_process;
  void pheno_process(int nobs, vector<int>& missing) {BaseInput::pheno_process(nobs, missing); block_count = (long) ceil(indiv_total/4.0);}

  void load_famfile(const string& filename);
  void load_bimfile(const string& filename);
  void prep_bedfile(const string& filename);

  unsigned long long load_snpdata(const set<long>& snps, float* buffer, SubsetMode filter_id);

  long allele_count(const set<long>& snps, Buffer<int>& counts, SubsetMode filter_id);
  void prep_filter(SubsetMode filter_id);
  
public:
  PlinkInput(Settings& s) : BaseInput(s), raw_snp_buffer(0), filter_buffer(0), block_count(0), filter_observed(0), active_filter(sm_All) {has_ped_info = true;
    _LOG << "Loading PLINK-format data..." << endl;
    
    string bed_file = settings.get_plinkfile("bed");
    geno_file.open(bed_file.c_str(), ios::in|ios::binary|ios::ate);
    bed_size = geno_file.tellg();
    
    load_famfile(settings.get_plinkfile("fam"));
    load_bimfile(settings.get_plinkfile("bim"));
    prep_bedfile(bed_file);

    init();
    _LOG << endl;
  }
  
  ~PlinkInput() {delete[] raw_snp_buffer; delete[] filter_buffer;}

  float misscode() {return 3;} 
};

#endif /** BASEINPUT_PLINKINPUT_H */



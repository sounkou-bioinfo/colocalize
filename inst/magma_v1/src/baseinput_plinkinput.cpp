/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "baseinput.h"
#include "baseinput_plinkinput.h"
#include "input.h"


void PlinkInput::load_famfile(const string& filename) {
  _LOG << "Reading file " << filename << "... "; _LOG.flush();
  TextInput fin(filename, 6);
  fin.set_error("reading .fam file", 2);

  string cid;
  long gender; double pheno;
  bool read_pheno = (pheno_type == ph_FromData);
  bool read_ped = (pedigree != pd_None);
  long index = 0;

  vector<int> miss_vec; vector<pair<string, string> > ped_vec;
  while (fin.process_line(false)) {
    if (!fin.read_num(gender, 4)) fin.line_error("gender value is not an integer");
    indiv_gender.push_back(gender == 1 || gender == 2 ? gender : 0);

    if (read_pheno) {
      if (fin[5] == "NA") {
        indiv_pheno.push_back(0);
        miss_vec.push_back(index);
      } else {
        if (!fin.read_num(pheno, 5)) fin.line_error("non-numeric phenotype value");
        indiv_pheno.push_back(pheno);
      }
    }
    
    if (read_ped) {
      pair<string, string> parents;
      if (fin[2] != "0") parents.first = fin[0] + " " + fin[2];
      if (fin[3] != "0") parents.second = fin[0] + " " + fin[3];      
      ped_vec.push_back(parents);
    }

    cid = fin[0] + " " + fin[1];
    if (indiv_map.find(cid) != indiv_map.end()) fin.line_error("duplicate FID+IID entry", string("ID = ").append(cid));
    indiv_map[cid] = index++;
  }
  if (index == 0) fin.error("no individuals in file");
  _LOG << index << " individuals read" << endl;

  pheno_process(index, miss_vec);
  if (read_ped) pedigree_process(ped_vec);
}

void PlinkInput::load_bimfile(const string& filename) {
  _LOG << "Reading file " << filename << "... "; _LOG.flush();

  TextInput fin(filename, 2);
  fin.set_error("reading .bim file", 2);

  snp_info.init(4*bed_size/indiv_total);
  while (fin.process_line(false)) snp_info.add_iid(fin[1]);
  snp_info.build_index();

  if (snp_info.data_used(false) == 0) fin.error("no SNPs in file");
  _LOG << snp_info.data_used(false) << " SNPs read" << endl;

  check_duplicates(filename, "PLINK .bim file");
}


void PlinkInput::prep_bedfile(const string& filename) {
  _LOG << "Preparing file " << filename << "... " << endl;
  unsigned long exp_bed_size = (unsigned long) ceil(indiv_total/4.0f) * snp_info.data_total() + 3; ///for SNP-major format
  geno_file.seekg(0, ios::beg);

  char buffer[3];
  geno_file.read(buffer, 3);
  if (geno_file.fail() || ((unsigned short) buffer[0] != 108 || (unsigned short) buffer[1] != 27)) {_LOG.error("reading .bed file") << "file is not a valid .bed file" << endl; die();} 
  if ((unsigned short) buffer[2] != 1) {
    if (buffer[2] == 0) _LOG.error("reading .bed file") << "file is in individual-major format (currently not supported)" << endl;
    else _LOG.error("reading .bed file") << "file-format specifier is invalid" << endl;
    die();
  }

  if (bed_size != exp_bed_size) {_LOG.error("reading .bed file") << "size of .bed file is inconsistent with number of SNPs and individuals in .bim and .fam files" << endl; die();} 

  float valIndex[] = {0.0,0.0,1.0,2.0}; valIndex[1] = misscode();
  for (int i = 0; i < 256; i++) {
    geno_count[i] = geno_miss[i] = 0;
    for (int j = 0; j < 4; j++) {
      geno_index[i][j] = valIndex[(i >> (2*j)) & 3];
      if (geno_index[i][j] == misscode()) geno_miss[i]++;
      else geno_count[i] += geno_index[i][j];
    }
  }
}

unsigned long long PlinkInput::load_snpdata(const set<long>& snps, float* buffer, SubsetMode filter_id) {
  unsigned long long read_index, write_index = 0;
  long sub_index = 0; float* sub_buffer = 0;
  if (!raw_snp_buffer) raw_snp_buffer = new char[block_count];
  short* filter = get_filter(filter_id);
  for (set<long>::const_iterator index = snps.begin(); index != snps.end(); ++index) {
    read_index = 3 + *index * block_count;
    geno_file.seekg(read_index);
    geno_file.read(raw_snp_buffer, block_count);

    for (long j = 0; j < indiv_total; j++) {
      if (j % 4 == 0) {
        sub_buffer = geno_index[(unsigned char) raw_snp_buffer[j/4]];
        sub_index = 0;
      }
      if (!filter[j]) {
        buffer[write_index++] = sub_buffer[sub_index];
      }
      sub_index++;
    }
  }

  return write_index;
}

long PlinkInput::allele_count(const set<long>& snps, Buffer<int>& counts, SubsetMode filter_id) { 
  if (!raw_snp_buffer) raw_snp_buffer = new char[block_count];

  prep_filter(filter_id); counts.set_size(snps.size(), 2, true); 
  unsigned long long read_index; int curr_snp = 0;
  for (set<long>::const_iterator index = snps.begin(); index != snps.end(); ++curr_snp, ++index) {
    int &curr_count = counts(curr_snp, 0), &curr_miss = counts(curr_snp, 1);
    read_index = 3 + *index * block_count;
    geno_file.seekg(read_index);
    geno_file.read(raw_snp_buffer, block_count);

    for (long i = 0; i < block_count; i++) {
      unsigned char curr = raw_snp_buffer[i] & filter_buffer[i];
      curr_count += geno_count[curr]; 
      curr_miss += geno_miss[curr];
    }
  }
  return filter_observed;
} 
            
void PlinkInput::prep_filter(SubsetMode filter_id) {
  if (filter_buffer && filter_id == active_filter) return;
  if (!filter_buffer) filter_buffer = new unsigned char[block_count];

  short* source = get_filter(filter_id);
  unsigned char keep[4] = {3, 12, 48, 192}; unsigned char* curr = filter_buffer;
  for (int i = 0; i < indiv_total; i++) {
    if (i && i % 4 == 0) {curr++; *curr = 0;}
    if (!source[i]) *curr |= keep[i%4];
  }
  active_filter = filter_id;
  filter_observed = indiv_total - MathUtils::count(source, indiv_total);
}



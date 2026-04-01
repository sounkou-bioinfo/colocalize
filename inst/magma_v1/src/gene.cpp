/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "gene.h"
#include "utils.h"
#include "mathutils.h"
#include <vector>
#include <set>

#define MIN_VAR 1e-8

using namespace std;
Buffer<float> Gene::temp_buffer;
SharedPtr<DelimitedOutput> LoadedGene::report_ptr; 

void Gene::check_block(int& offset, int& total) {
  if (offset >= no_snps) offset = 0;
  total = total ? min(no_snps-offset, total) : no_snps-offset;
}

BufferWindow<float> Gene::get_datablock(int& offset, int& total) {
  check_block(offset, total);
  return data.window(offset, total);
}

void Gene::set_rescale() {
  for (int i_snp = 0; i_snp < no_snps; i_snp++) base_stats(i_snp, bs_Rescale) = samp_size / float(samp_size - base_stats(i_snp, bs_Missing));
}

void Gene::snp_center() {
  if (state >= ProcessType::Centered) return;
  for (int i_snp = 0; i_snp < no_snps; i_snp++) {
    float& mean = base_stats(i_snp, bs_Mean);  
    for (float *curr = data[i_snp], *end = data[i_snp+1]; curr < end; curr++) *curr -= mean;
  }
  clear_missing();
  base_stats.assign_zero_col(bs_Mean);
  state = ProcessType::Centered;
}

void Gene::snp_normalize() {
  if (state >= ProcessType::Normalized) return;
  for (int i_snp = 0; i_snp < no_snps; i_snp++) {
    float& mean = base_stats(i_snp, bs_Mean), sd = sqrt(base_stats(i_snp, bs_Variance)); 
    if (sd > 0) {for (float *curr = data[i_snp], *end = data[i_snp+1]; curr < end; curr++) *curr = (*curr - mean) / sd;}
    else Utils::set_zero(data[i_snp], samp_size);
    base_stats(i_snp, bs_Variance) = 1;    
  }
  clear_missing();
  base_stats.assign_zero_col(bs_Mean);
  state = ProcessType::Normalized;
}

void Gene::clear_missing() {
  for (int i_snp = 0; i_snp < no_snps; i_snp++) {  
    if (base_stats(i_snp, bs_Missing) > 0) {
      float* curr = data[i_snp];
      deque<int>& miss_index = missing.traverse(i_snp);
      for (int i_miss = 0; i_miss < miss_index.size(); i_miss++) curr[miss_index[i_miss]] = 0;
    }    
  }
}

void Gene::prep_product() {
  if (state >= ProcessType::Normalized) return;
  
  if (normalize_geno) snp_normalize();
  else snp_center();
}

Buffer<float>& Gene::get_rescale_block(int& offset, int& total, bool norm) {
  check_block(offset, total); temp_buffer.set_size(total,1); 
  if (norm) for (int i = 0; i < total; i++) temp_buffer(i) = sqrt(base_stats(offset+i, bs_Rescale)/base_stats(offset+i, bs_Variance)); 
  else for (int i = 0; i < total; i++) temp_buffer(i) = base_stats(offset+i, bs_Rescale); 
  return temp_buffer;
}

Buffer<float>& LinkedGene::get_rescale_block(int& offset, int& total, bool norm) {
  check_block(offset, total); temp_buffer.set_size(total,1); 
  if (norm) for (int i = 0; i < total; i++) temp_buffer(i) = sqrt(base_stats(offset+i, bs_Rescale));
  else for (int i = 0; i < total; i++) temp_buffer(i) = base_stats(offset+i, bs_Rescale) * sqrt(base_stats(offset+i, bs_Variance)); 
  return temp_buffer;
}

float Gene::get_rescale_pair(Buffer<float>& bs_I, BitBuffer& miss_I, int i_snp, Buffer<float>& bs_J, BitBuffer& miss_J, int j_snp, bool samp_scale) {
  float scale = samp_scale ? 1/double(samp_size-1) : 1;
  int miss_i = bs_I(i_snp, bs_Missing), miss_j = bs_J(j_snp, bs_Missing);
  if (miss_i == 0 || miss_j == 0) return bs_I(i_snp, bs_Rescale) * bs_J(j_snp, bs_Rescale) * scale;

  int miss = miss_i + miss_j - miss_I.both(i_snp, j_snp, miss_J);
  return samp_size / float(samp_size-miss) * scale;
}

float LinkedGene::get_rescale_pair(Buffer<float>& bs_I, BitBuffer& miss_I, int i_snp, Buffer<float>& bs_J, BitBuffer& miss_J, int j_snp, bool samp_scale) {
  float base_scale = Gene::get_rescale_pair(bs_I, miss_I, i_snp, bs_J, miss_J, j_snp, false);
  base_scale *= sqrt(bs_I(i_snp, bs_Variance)) * sqrt(bs_J(j_snp, bs_Variance));
  if (samp_scale) base_scale /= double(samp_size-1);
  else base_scale *= samp_size_rescale;
  
  return base_scale;
}

float Gene::sum_rescale(int snp, bool norm) {
  if (norm) return base_stats(snp, bs_Rescale) / sqrt(base_stats(snp, bs_Variance)); 
  else return base_stats(snp, bs_Rescale);
}

float LinkedGene::sum_rescale(int snp, bool norm) {
  if (norm) return base_stats(snp, bs_Rescale);
  else return base_stats(snp, bs_Rescale) * sqrt(base_stats(snp, bs_Variance)); 
}

void Gene::do_rescale_pair(float& target, int i_snp, int j_snp, ProductType::Type mode, Gene* other) {bool samp_scale = mode >= ProductType::Covariance;
  Buffer<float>& bs_other = other ? other->base_stats : base_stats; BitBuffer& miss_other = other ? other->missing : missing;
  target *= get_rescale_pair(base_stats, missing, i_snp, bs_other, miss_other, j_snp, samp_scale);
  if (mode == ProductType::Correlation) {
    float scale = sqrt(base_stats(i_snp, bs_Variance)) * sqrt(bs_other(j_snp, bs_Variance));   
    target /= scale;
  }
}

void Gene::rescale_product(Buffer<float>& target, int offset_i, int total_i, int offset_j, int total_j, ProductType::Type mode, Gene* other) {
  for (int i_snp = 0; i_snp < total_i; i_snp++) {
    for (int j_snp = 0; j_snp < total_j; j_snp++) {
      do_rescale_pair(target(i_snp,j_snp), offset_i+i_snp, offset_j+j_snp, mode, other);
    }
  }
}

void Gene::rescale_product_skip(Buffer<float>& target, int offset_i, int step_i, int total_i, int offset_j, int step_j, int total_j, ProductType::Type mode, Gene* other) {
  for (int i = 0, i_snp = offset_i; i < total_i; i++, i_snp += step_i) {
    for (int j = 0, j_snp = offset_j; j < total_j; j++, j_snp += step_j) {
      do_rescale_pair(target(i,j), i_snp, j_snp, mode, other);
    }
  }
}      

float Gene::product_var_core(float* var, int snp, ProcessType::Type mode) {
  return MathUtils::in_product(data[snp], var, samp_size) * sum_rescale(snp, mode == ProcessType::Normalized); 
}

float LinkedGene::product_var_core(float* var, int snp, ProcessType::Type mode) {
  if (!allow_prod) block_prod();
  double out = 0; var += subset_offset; 
  for (float *curr = data[snp], *end = curr+samp_size; curr < end; curr++, var += subset_step) out += *curr * *var;
  return out * sum_rescale(snp, mode == ProcessType::Normalized); 
}

float Gene::product_var(float* var, int snp, ProcessType::Type mode) {prep_product(); 
  return product_var_core(var, snp, mode);
}

void Gene::product_var(Buffer<float>& target, Buffer<float>& var, int offset, int total, bool transpose, ProcessType::Type mode) {prep_product(); 
  int nvar = var.ncol(); check_block(offset, total);
  if (transpose) {
    target.set_size(nvar, total);
    for (int i = 0; i < total; i++) {for (int j = 0; j < nvar; j++) target(j,i) = product_var_core(var[j], i, mode);}
  } else{  
    target.set_size(total, nvar);
    for (int i = 0; i < total; i++) {for (int j = 0; j < nvar; j++) target(i,j) = product_var_core(var[j], i, mode);}
  }
}

void Gene::product(Buffer<float>& target, int offset, int total, ProductType::Type mode) {prep_product();
  MathUtils::matrix_prod(get_datablock(offset, total).view(), target, MathUtils::TransposeFirst);
  for (int i = 0; i < total; i++) {
    target(i,i) *= sum_rescale(offset+i); 
    for (int j = i+1; j < total; j++) {
      do_rescale_pair(target(j,i), offset+i, offset+j, ProductType::Centered);
      target(i,j) = target(j,i);
    }
  }

  if (mode == ProductType::Covariance) MathUtils::matrix_scale(target, 1/double(samp_size-1));
  else if (mode == ProductType::Correlation) MathUtils::cov_to_cor(target);
}

void LinkedGene::product(Buffer<float>& target, int offset, int total, ProductType::Type mode) {
  Gene::product(target, offset, total, ProductType::Centered);
  for (int i = 0; i < total; i++) target(i,i) = base_stats(offset+i,bs_Variance)*(parent_samp_size-1);
  
  if (mode == ProductType::Covariance) MathUtils::matrix_scale(target, 1/double(parent_samp_size-1));
  else if (mode == ProductType::Correlation) MathUtils::cov_to_cor(target);
}

void Gene::product(Buffer<float>& target, int offset_i, int total_i, int offset_j, int total_j, ProductType::Type mode) {prep_product(); 
  MathUtils::matrix_prod(get_datablock(offset_i, total_i).view(), get_datablock(offset_j, total_j).view(), target, MathUtils::TransposeFirst);
  rescale_product(target, offset_i, total_i, offset_j, total_j, mode);
}

void Gene::product(Buffer<float>& target, Gene* other, ProductType::Type mode) {product(target, other, 0, no_snps, 0, other->no_snps, mode);}
void Gene::product(Buffer<float>& target, Gene* other, int offset_i, int total_i, int offset_j, int total_j, ProductType::Type mode) {
  prep_product(); other->prep_product();
  MathUtils::matrix_prod(get_datablock(offset_i, total_i).view(), other->get_datablock(offset_j, total_j).view(), target, MathUtils::TransposeFirst);
  rescale_product(target, offset_i, total_i, offset_j, total_j, mode, other);
}

void Gene::product_skip(Buffer<float>& target, int offset_i, int step_i, int total_i, int offset_j, int step_j, int total_j, ProductType::Type mode) {
  if (step_i <= 1 && step_j <= 1) {product(target, offset_i, total_i, offset_j, total_j, mode); return;}
  prep_product(); check_block(offset_i, total_i); check_block(offset_j, total_j);

  int nrow = ceil(float(total_i) / step_i), ncol = ceil(float(total_j) / step_j);
  if (offset_i < offset_j) offset_i += (total_i-1) % step_i;
  else offset_j += (total_j-1) % step_j;  

  target.set_size(nrow, ncol);
  target.matrix().noalias() = data.matrix_colskip(offset_i, nrow, step_i).transpose() * data.matrix_colskip(offset_j, ncol, step_j);

  rescale_product_skip(target, offset_i, step_i, nrow, offset_j, step_j, ncol, mode);
}

void Gene::product_skip(Buffer<float>& target, Gene* other, int step_i, int step_j, ProductType::Type mode) {
  if (step_i <= 1 && step_j <= 1) product(target, other, 0, no_snps, 0, other->no_snps, mode);
  else product_skip(target, other, 0, step_i, no_snps, 0, step_j, other->no_snps, mode);  
}
  
void Gene::product_skip(Buffer<float>& target, Gene* other, int offset_i, int step_i, int total_i, int offset_j, int step_j, int total_j, ProductType::Type mode) {
  if (step_i <= 1 && step_j <= 1) {product(target, other, offset_i, total_i, offset_j, total_j, mode); return;}
  prep_product(); check_block(offset_i, total_i);
  other->prep_product(); other->check_block(offset_j, total_j); 

  int nrow = ceil(float(total_i) / step_i), ncol = ceil(float(total_j) / step_j);
  offset_i += (total_i-1) % step_i;
  
  target.set_size(nrow, ncol);
  target.matrix().noalias() = data.matrix_colskip(offset_i, nrow, step_i).transpose() * (other->data).matrix_colskip(offset_j, ncol, step_j);

  rescale_product_skip(target, offset_i, step_i, nrow, offset_j, step_j, ncol, mode, other);
}

void Gene::data_processed(Buffer<float>& target, int offset, int total, bool norm) {prep_product();
  Buffer<float>& rescale = get_rescale_block(offset, total, norm);
  target.set_size(samp_size, total);
  target.matrix().noalias() = data.matrix(offset, total) * rescale.matrix().asDiagonal();
}

void Gene::data_product(Buffer<float>& target, Buffer<float>& var, int offset, int total) {prep_product();
  Buffer<float>& rescale = get_rescale_block(offset, total, true);
  target.set_size(samp_size, total);
  target.matrix().noalias() = var.matrix().asDiagonal() * data.matrix(offset, total) * rescale.matrix().asDiagonal();
}

void Gene::data_projection(Buffer<float>& target, Buffer<float>& rotate, int offset, int total) {prep_product();
  Buffer<float>& rescale = get_rescale_block(offset, total);
  target.set_size(samp_size, rotate.ncol());
  target.matrix().noalias() = data.matrix(offset, total) * rescale.matrix().asDiagonal() * rotate.matrix();
}


void LinkedGene::init() {
  if (!parent) return;
  parent_state = parent->state;
  if (parent_state < ProcessType::Centered) {_LOG.error("processing gene") << "LinkedGene object cannot be linked to unprocessed genotype data" << endl; die();} 

  parent_samp_size = parent->samp_size; 
  samp_size = (parent_samp_size - subset_offset) / subset_step;
  if ( (parent_samp_size - subset_offset) % subset_step > 0 ) samp_size++;
  samp_size_rescale *= (parent_samp_size-1)/double(samp_size-1);
  
  no_snps = parent->no_snps; no_snps_orig = parent->no_snps_orig; 
  no_rare = parent->no_rare; miss_code = parent->miss_code;
  
  snp_pval.assign(parent->snp_pval); snp_N.assign(parent->snp_N); snp_type.assign(parent->snp_type); snp_weight.assign(parent->snp_weight);
  base_stats.assign(parent->base_stats);

  data.set_size(samp_size, no_snps);
  data.matrix().noalias() = parent->data.matrix_rowskip(subset_offset, samp_size, subset_step);
  
  missing.set_size(samp_size, no_snps, true);
  BitBuffer& parent_miss = parent->missing;  

  for (int i_snp = 0; i_snp < no_snps; i_snp++) {
    float sum = 0, square = 0; int count_pos = 0, count_neg = 0;
    for (float *curr = data[i_snp], *end = curr+samp_size; curr < end; curr++) {
      sum += *curr; square += *curr * *curr; 
      if (*curr > 0) count_pos++; else if (*curr < 0) count_neg++;
    }
   
    if (count_pos <= outlier_thresh || count_neg <= outlier_thresh) {
      if (count_pos <= outlier_thresh) {for (int i = 0; i < samp_size; i++) {if (data(i, i_snp) > 0) missing.set(i, i_snp);}}
      if (count_neg <= outlier_thresh) {for (int i = 0; i < samp_size; i++) {if (data(i, i_snp) < 0) missing.set(i, i_snp);}}
      if (!base_stats(i_snp, bs_Missing)) base_stats(i_snp, bs_Missing) = 1;
    }
       
    if (base_stats(i_snp, bs_Missing) > 0) {
      for (int i_self = 0, i_parent = subset_offset; i_self < samp_size; i_self++, i_parent += subset_step) {
        if (parent_miss.get(i_parent, i_snp)) {float &x = data(i_self, i_snp);
          missing.set(i_self, i_snp);
          sum -= x; square -= x*x;
        }
      }
      base_stats(i_snp, bs_Missing) = missing.count(i_snp);
    }
    int obs = samp_size - base_stats(i_snp, bs_Missing);
    double mean = sum / obs;
    double var = (square/obs - mean*mean) * samp_size / double(samp_size - 1);
    base_stats(i_snp, bs_Mean) = mean;
    base_stats(i_snp, bs_Variance) = var > MIN_VAR ? var : 0;
  }
  set_rescale();

  state = ProcessType::Raw;
  snp_normalize();

  if (parent_state < state) base_stats.assign_col(parent->base_stats[bs_Variance], bs_Variance);  
}


/// DATA LOADING

void LoadedGene::init() {
  full_qc = !settings["snp_skip_qc"];
  dummy_pheno = settings["dummy_pheno"];
  is_genotype = settings.gets("data_type") == "binary_geno"; //split into separate gene types later

  max_miss = settings.getn(full_qc ? "snp_max_miss" : "snp_max_miss_max") * samp_size;
  min_maf = settings.getn("snp_min_maf"); max_maf = settings.getn("snp_max_maf", 1);
  min_mac = settings.geti("snp_min_mac"); max_mac = settings.geti("snp_max_mac", 2*samp_size);

  if (is_genotype) {
    int lower = max(min_mac, int(min_maf*samp_size*2));
    int upper = min(max_mac, int(max_maf*samp_size*2));    
    
    if (lower >= upper) {_LOG.error("processing gene") << "minimum MAF/MAC thresholds exceed maximum MAF/MAC thresholds, given sample size" << endl; die();} 
  }

  has_rare = is_genotype && settings["do_rare"] && !settings["has_pval"];
  if (has_rare) {
    rare_maf_thresh = settings.getn("rare_cutoff_maf", 0);
    rare_mac_thresh = settings.geti("rare_cutoff_mac", 0);      

    rare_only = settings["rare_only"];
    rare_normalize = !settings["rare_freq_unweighted"];      
    rare_multi_max = settings.geti("rare_multi_max_count");
  }

  multi_partition = settings["snp_partition_overlap"];
  
  snp_report = settings["tmp_snps_used"];
  if (snp_report && !report_ptr.get()) {
    report_ptr.set(new DelimitedOutput(settings.gets("out_prefix") + ".snps.status", '\t'));
    DelimitedOutput& fout = *(report_ptr.get());
    fout.print_comment("status codes: negative means dropped, positive means retained (1 = normal SNP, 2 = part of burden score)");
    fout.print_comment("status codes: -1 = too many missing, -2 = monomorphic SNP, -3 = MAF/MAC too low, -4 = MAF/MAC too high");
    fout.print_comment("status codes: -5 = variance too low, -6 = differential missingness, -7 = missing/invalid SNP info value, -8 = common variant (rare-only mode)"); 
    fout << "GENE_NAME" << "SNP_INDEX" << "STATUS_CODE" << endl; 
  } 

  double t_val = ConvertPvalToT(samp_size-2).convert(settings.getn("snp_diff_pval")/2);
  max_diff_corr = t_val / sqrt(t_val*t_val + samp_size - 2);  
}

void LoadedGene::prep_gene() {
  dropped_snps.clear(); rare_snps.clear();
  snp_pval.set_empty(); snp_N.set_empty(); snp_type.set_empty(); snp_weight.set_empty();
}

Buffer<float>& LoadedGene::set_buffer(int nsnp) {no_snps = nsnp;
  prep_gene();

  data.set_size(samp_size, no_snps);
  base_stats.set_size(no_snps, bs_Count+bs_Extra, true);
  if (snp_report) snp_status.assign(no_snps, 1);

  return data;
}

int LoadedGene::count_unique(vector<int>& index, bool flip) {set<int> unique;
  if (snp_id.empty()) return -1;
  if (flip) {
    for (int i = 0; i < no_snps; i++) unique.insert(snp_id(i));            
    return unique.size() - count_unique(index);  
  } else {
    for (int i = 0; i < index.size(); i++) unique.insert(snp_id(index[i]));
    return unique.size();
  }
}

void LoadedGene::process_data(float* pheno) {state = ProcessType::Raw;
  missing.set_size(samp_size, no_snps, true);

  bool has_pheno = !dummy_pheno && pheno;
  bool check_diff_miss = full_qc && has_pheno;
//  bool allow_mono = !has_pheno;

  if (!has_pheno) pheno = data.begin();
  bool rare = false; long minor;

  for (int i_snp = 0; i_snp < no_snps; i_snp++) {
    double sum = 0, miss_pheno = 0, square = 0; int miss = 0; float* curr_snp = data[i_snp]; DropType drop = dt_Retained;
    for (int i_val = 0; i_val < samp_size; i_val++) {float& curr = curr_snp[i_val];
      if (curr == miss_code) {
        miss++; miss_pheno += pheno[i_val];
        missing.set(i_val,i_snp);
      }
      else {sum += curr; square += curr * curr;}
    }

    long obs = samp_size - miss; 
    if (miss > max_miss) drop = dt_Missing;


    if (is_genotype) {
      minor = round(min(sum, 2*obs - sum));
      if (minor == 0) {
        /*if (allow_mono) // PROCESS SNP
        else */ drop = dt_Mono;
      } else if (full_qc && !drop) {
          if (minor < min_mac || minor < min_maf*obs*2) drop = dt_LowFreq;
          else if (minor > max_mac  || minor > max_maf*obs*2) drop = dt_HighFreq;
      }
    } else minor = 0;

    double mean = sum / obs;
    double var = (square/obs - mean*mean) * samp_size / double(samp_size - 1);

    if (!drop) {
      if (var <= MIN_VAR) drop = dt_LowVar;
      if (check_diff_miss && miss > 0) {
        double var_miss = miss * (1 - double(miss) / samp_size) / (samp_size-1);
        double diff_corr = abs(miss_pheno) / (samp_size - 1) / sqrt(var_miss);
        if (diff_corr >= max_diff_corr) drop = dt_DiffMiss;
      } 
    }
    if (!drop) {      
      if (!snp_pval.empty()) {
        if (snp_pval(i_snp) < 0) drop = dt_InfoMiss;
        if (!snp_N.empty() && snp_N(i_snp) == 0) drop = dt_InfoMiss;      
      }
      
      if (!snp_type.empty() && !snp_type(i_snp)) drop = dt_InfoMiss;
      if (!snp_weight.empty() && snp_weight(i_snp) <= 0) drop = dt_InfoMiss;      
    }

    if (!drop && has_rare) {
      if (rare_maf_thresh && rare_mac_thresh) rare = (minor <= min(rare_mac_thresh, int(2*obs*rare_maf_thresh)));
      else if (rare_maf_thresh) rare = (minor <= int(2*obs*rare_maf_thresh));
      else if (rare_mac_thresh) rare = (minor <= rare_mac_thresh);      
      else rare = false;

      if (!rare && rare_only) drop = dt_NonRare;
    } 

    if (drop) dropped_snps.push_back(i_snp);
    else {
      if (rare) rare_snps.push_back(i_snp);

      base_stats(i_snp, bs_Missing) = miss;
      base_stats(i_snp, bs_Freq) = sum / float(2*obs);
      base_stats(i_snp, bs_MAF) = minor / float(2*obs);
      base_stats(i_snp, bs_MAC) = minor;      
      base_stats(i_snp, bs_Mean) = mean;      
      base_stats(i_snp, bs_Variance) = var;
    }
    if (snp_report) {
      if (drop) snp_status[i_snp] = -drop;
      else if (rare) snp_status[i_snp] = 2;
    }
  }
 
  if (multi_partition) {
    no_rare = count_unique(rare_snps);
    no_snps_orig = count_unique(dropped_snps, true);
  } else {
    no_rare = rare_snps.size();
    no_snps_orig = no_snps - dropped_snps.size();
  }
  
  if (snp_report) {
    DelimitedOutput& fout = *(report_ptr.get());
    for (int i = 0; i < no_snps; i++) fout << gene_name << snp_id(i) << snp_status[i] << endl;
  }
}

void LoadedGene::filter_snps() {
  if (!dropped_snps.empty()) {
    if (no_rare > 0) sort(dropped_snps.begin(), dropped_snps.end());
    int no_dropped = dropped_snps.size();
    data.drop_cols_by_index(dropped_snps);    
    missing.drop_cols_by_index(dropped_snps);    
    base_stats.drop_rows_by_index(dropped_snps);
    if (!snp_pval.empty()) snp_pval.drop_rows_by_index(dropped_snps);
    if (!snp_N.empty()) snp_N.drop_rows_by_index(dropped_snps);    
    if (!snp_type.empty()) snp_type.drop_rows_by_index(dropped_snps);        
    if (!snp_weight.empty()) snp_weight.drop_rows_by_index(dropped_snps);            
    dropped_snps.clear(); rare_snps.clear();
    no_snps -= no_dropped;
  }
    
  if (no_snps <= 0) throw GeneException("empty gene", "gene contains no valid SNPs after internal QC");
}

void LoadedGene::process_rare(vector<int>& rare_index, bool by_type, bool no_max) {
  if (!has_rare || no_rare <= 0 || rare_index.empty()) return;
  if (rare_index.size() == 1) {
    process_rare_single(rare_index[0]);
    return;
  }
  
  if (by_type && !snp_type.empty()) { ///currently unused, does not support weights
    int curr_type = snp_type(rare_index[0]), total = rare_index.size(); vector<int> curr_snps; 
    for (int i = 0; i <= total; i++) {int& id = rare_index[i];
      if (i == total || snp_type(id) != curr_type) {
        process_rare(curr_snps, false, false);
        if (i == total) break;
        curr_type = snp_type(id); 
        curr_snps.clear();      
      }
      curr_snps.push_back(id);
    }
  } else {
    if (!no_max && rare_multi_max && rare_index.size() > rare_multi_max) {
      int tot = rare_index.size(), nscore = tot / rare_multi_max + (tot % rare_multi_max > 0), base_size = tot / nscore, offset = 0;
      for (int i = 0; i < nscore; i++) {
        vector<int> curr_snps;
        int end = offset + base_size + (i < (tot % base_size)); if (end > tot) end = tot;
        for (int j = offset; j < end; j++) curr_snps.push_back(rare_index[j]);      
        process_rare(curr_snps, false, true);
        offset = end;
      }
    } else process_rare_core(rare_index);
  }
}

void LoadedGene::process_rare_core(vector<int>& rare_index) {  
  int write_index = rare_index[0]; float* write_target = data[write_index];
  bool use_weights = !snp_weight.empty(); 

  for (int i = 0; i < rare_index.size(); i++) {int& index = rare_index[i];
    float* geno = data[index];    
    bool flip = base_stats(index, bs_Freq) > 0.5;

    float weight = use_weights ? snp_weight(index) : 1;
    if (rare_normalize) weight /= sqrt(base_stats(index, bs_MAF)*(1-base_stats(index, bs_MAF)));
        
    if (i == 0) {
      for (float *curr = geno, *end = geno+samp_size; curr < end; curr++) {
        if (*curr == miss_code) *curr = 0;
        else if (flip) *curr = weight*(2 - *curr);
      }
    } else {
      for (float *in = geno, *out = write_target, *end = geno+samp_size; in < end; in++, out++) {                
        if (*in != miss_code) *out += weight*(flip ? 2 - *in : *in);
      }
      base_stats(write_index, bs_MAF) += base_stats(index, bs_MAF);
      base_stats(write_index, bs_MAC) += base_stats(index, bs_MAC);   /// this is left as the sum of MAC values across the variants
      dropped_snps.push_back(index);
    }              
  }
  base_stats(write_index, bs_MAF) /= rare_index.size();

  base_stats(write_index, bs_Missing) = 0;
  base_stats(write_index, bs_Freq) = base_stats(write_index, bs_MAF);  
  if (base_stats(write_index, bs_MAF) > 1) base_stats(write_index, bs_MAF) = 0;
  else if (base_stats(write_index, bs_MAF) > 0.5) base_stats(write_index, bs_MAF) = 1 - base_stats(write_index, bs_MAF);

  float mean, var; MathUtils::get_stats(write_target, samp_size, mean, var, false);
  base_stats(write_index, bs_Mean) = mean;
  base_stats(write_index, bs_Variance) = var;

  if (base_stats(write_index, bs_Variance) <= MIN_VAR) dropped_snps.push_back(write_index);
  else if (!normalize_geno) {
    double scale = sqrt(sqrt(rare_index.size()));
    MathUtils::vec_transform(write_target, -mean, scale/sqrt(var), samp_size, true);
    base_stats(write_index, bs_Mean) = 0; base_stats(write_index, bs_Variance) = scale*scale;
  }
}

void LoadedGene::process_rare_single(int index) {
  int impute = (base_stats(index, bs_Freq) > 0.5) ? 2 : 0; float* geno = data[index];      
  for (float *curr = geno, *end = geno+samp_size; curr < end; curr++) {if (*curr == miss_code) *curr = impute;}
  base_stats(index, bs_Missing) = 0;
}


void PlainBinaryGene::process(float* pheno) {
  process_data(pheno);
  filter_snps();
  set_rescale();
}

void ProcessedGene::process(float* pheno) {
  process_data(pheno);
  if (has_rare && no_rare > 0) process_rare();
  filter_snps();
  set_rescale();
  prep_product();
}



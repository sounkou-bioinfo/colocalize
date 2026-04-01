/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genestats_main.h"
#include "mathutils.h"
#include "statutils.h"
#include "exceptions.h"

#define PC_EIGEN_CUTOFF 0.0001

void DefaultStats::clear() {
  GeneStats::clear();
   
  xty.reset(); xtz.reset(); xtx.reset(); xtx_buffer.reset(); ytz.reset(); xty_multi.reset();
  phenotype.reset(); gender.reset(); covar.reset(); interactor.reset();  
  
  for (int i = 0; i < pc_data.size(); i++) delete pc_data[i];  
  pc_data.clear();
}

void DefaultStats::rewind() {
  GeneStats::rewind();
  has_xtx = has_xty = has_xtz = buffer_xtx = false;   
  for (int i = 0; i < pc_data.size(); i++) pc_data[i]->rewind();
}

GeneStats::BufferDataBlock* DefaultStats::eject_pcs(int pcid) {
  if (pcid >= 0) set_pcs(pcid);
  BufferDataBlock* out = new BufferDataBlock(baseinput.gene_info.location.get(gene_id), princomp->get_data());
  princomp->rewind();
  return out;
}

void DefaultStats::partition(vector<int>& blocks) {
  GeneStats::partition(blocks);
  int nblock = block_index.nrow(); 
  for (int i = 0; i < pc_data.size(); i++) pc_data[i]->set_partition(nblock); 
}

void DefaultStats::set_block(int index) {
  GeneStats::set_block(index);
  has_xtx = false;

}

/// Normalizes phenotype, normalizes and decorrelates covariates
void DefaultStats::init() { try {
  if (initialized) return;
  gs_id = instance_count++;

  messages.clear();
  if (Utils::chr_type(chromosome) > 0) {
    if (chrX_only_sex_covar) messages[ModelChange] = "adding gender as covariate";
  }
  
  set_filter(); 
  samp_size = data_size - MathUtils::sum(filter, data_size);

  if (samp_size < settings.geti("min_sample_size")) {
    string err = "after removing individuals with phenotype/covariate missing values, sample size is too small (minimum = " + Utils::num_string(settings.geti("min_sample_size")) + ")";
    throw DataException("missing data", err);
  }

  no_pheno = 1+baseinput.indiv_pheno_multi.ncol();      
  no_covar = baseinput.ncov(BaseInput::vt_Regular) + baseinput.ncov(BaseInput::vt_SNP) + sex_covar; 
  snp_covar_index = baseinput.ncov(BaseInput::vt_Regular);
  sex_covar_index = snp_covar_index + baseinput.ncov(BaseInput::vt_SNP);
  
  phenotype.set_size(samp_size, 1);
  gender.set_size(samp_size, 1);    
  covar.set_size(samp_size, no_covar);
  if (no_pheno > 1) phenotype_multi.set_size(samp_size, no_pheno-1);

  int offset = 0;
  for (int i = 0; i < data_size; i++) {
    if (filter[i]) continue;
    phenotype(offset) = baseinput.indiv_pheno[i];
    if (no_pheno > 1) {for (int j = 0; j < (no_pheno - 1); j++) phenotype_multi(offset,j) = baseinput.indiv_pheno_multi(i,j);}
    gender(offset) = baseinput.indiv_gender[i];
    for (int j = 0; j < baseinput.ncov(BaseInput::vt_Regular); j++) covar(offset,j) = baseinput.covar(i,j);
    for (int j = 0; j < baseinput.ncov(BaseInput::vt_SNP); j++) covar(offset,snp_covar_index+j) = baseinput.covar_snps(i,j);    
    if (sex_covar) covar(offset,sex_covar_index) = baseinput.indiv_gender[i];
    offset++;
  }

  float cutoff = 1 - max(settings.geti("pheno_min_count")/double(samp_size), settings.getn("pheno_min_prop")), cap = 1 - settings.geti("pheno_min_cap")/double(samp_size);
  if (cutoff < cap) cutoff = cap;
  if (MathUtils::get_mode(phenotype.begin(), samp_size, true) > cutoff) {
    string cov = no_covar > 0 ? " after correcting for covariates" : "";
    string note = cov + (baseinput.binary_pheno ? " (note: for binary phenotypes -9 and 0 are counted as missing values)" : "");
    throw DataException("phenotype error", "phenotype is (nearly) constant, (almost) all individuals have the same value" + note);
  }
  
  bool success = MathUtils::normalize(phenotype.begin(), samp_size);
  if (!success) throw DataException("phenotype error", "phenotype has a variance of zero");
  for (int i = 0; i < no_pheno - 1; i++) {if (!MathUtils::normalize(phenotype_multi[i], samp_size)) throw DataException("phenotype error", "one of additional phenotypes has a variance of zero");}
  if (no_covar > 0) process_covar();

  initialized = true;
  
} catch (const DataException& de) {throw;} 
  catch (const exception& e) {check_mem_error(e, "initializing data"); throw DataException("initializing data", "an error occurred when trying to initialize phenotype/covariate data", e.what());}
}

void DefaultStats::process_covar() {
  if (no_covar == 0) return;
  if (use_interactor) {
    use_residuals = true;
    if ( (baseinput.interact_type == BaseInput::vt_Gender && !sex_covar) || (baseinput.interact_type != BaseInput::vt_Gender && baseinput.interact_index < 0) ) throw DataException("covariate error", "cannot find interactor variable");
    if (baseinput.interact_type == BaseInput::vt_Gender) interactor.assign(covar[sex_covar_index], samp_size);
    else if (baseinput.interact_type == BaseInput::vt_SNP) interactor.assign(covar[snp_covar_index + baseinput.interact_index], samp_size);    
    else interactor.assign(covar[baseinput.interact_index], samp_size);
    bool success = MathUtils::normalize(interactor.begin(), samp_size);      
    if (!success) throw DataException("covariate error", "interactor variable has a variance of zero");        
  }

  vector<int> dropped;   
  for (int i = 0; i < no_covar; i++) {
    bool success = MathUtils::normalize(covar[i], samp_size);
    if (!success) {
      dropped.push_back(i);
      if (chrX_only_sex_covar && sex_covar_index == i) messages[DroppedCovarChrX] = "sex";
    }
  }
  if (dropped.size() == no_covar) throw DataException("covariate error", "all covariates have a variance of zero");

  if (!dropped.empty()) {
    covar.drop_cols_by_index(dropped); no_covar = covar.ncol();
    string note; 
    for (int i = 0; i < dropped.size(); i++) { 
      int index = dropped[i]; BaseInput::VarType type = BaseInput::vt_Regular;
      if (index >= sex_covar_index) {index = 0; type = BaseInput::vt_Gender;}
      else if (index >= snp_covar_index) {index -= snp_covar_index; type = BaseInput::vt_SNP;}
      if (i > 0) note += " ";
      note += baseinput.cov_name(index, type);      
    }
    
    messages[DroppedCovar] = note;
  }      

  if  (no_covar > 1) {
    Buffer<float> corr; corr.set_square(no_covar);
    for (int i = 0; i < no_covar; i++) {
      corr(i,i) = 1;
      for (int j = i+1; j < no_covar; j++) corr(j,i) = MathUtils::in_product(covar[i], covar[j], samp_size)/(samp_size-1);
    }
    no_covar = MathUtils::rotate(covar, corr, PC_EIGEN_CUTOFF);
  }

  ytz.set_size(no_covar, 1);
  for (int i = 0; i < no_covar; i++) ytz(i) = MathUtils::in_product(covar[i], phenotype[0], samp_size);

  Buffer<float> pheno_copy(phenotype);
  for (int i = 0; i < no_covar; i++) MathUtils::vec_product(covar[i], pheno_copy[0], -ytz(i) / (samp_size-1), samp_size, true);      
  if (MathUtils::get_var(pheno_copy[0], samp_size) < settings.getn("pheno_min_resvar")) throw DataException("phenotype error", "almost all phenotypic variance is explained by covariates, too little left to explain for genotype data");

  if (use_residuals) {
    phenotype.assign(pheno_copy);
    MathUtils::normalize(phenotype.begin(), samp_size);

    if (no_pheno > 1) {
      Buffer<float> ytz_buff(no_covar, 1);
      for (int p = 0; p < (no_pheno-1); p++) {
        for (int i = 0; i < no_covar; i++) ytz_buff(i) = MathUtils::in_product(covar[i], phenotype_multi[p], samp_size);
        for (int i = 0; i < no_covar; i++) MathUtils::vec_product(covar[i], phenotype_multi[p], -ytz_buff(i) / (samp_size-1), samp_size, true);      
        if (MathUtils::get_var(phenotype_multi[p], samp_size) < settings.getn("pheno_min_resvar")) throw DataException("phenotype error", "almost all phenotypic variancein additional phenotype is explained by covariates, too little left to explain for genotype data");
        MathUtils::normalize(phenotype_multi[p], samp_size);
      }
    }

    if (use_interactor) {
      covar.assign(interactor); no_covar = 1; ytz.set_size(1,1); 
      ytz(0) = MathUtils::in_product(covar[0], phenotype[0], samp_size);
    } else {covar.reset(); ytz.reset(); no_covar = 0;}
  } 
  if (no_pheno > 1 && (!use_residuals || use_interactor)) throw DataException("phenotype error", "cannot use unresidualized covariates or interactors when analyzing multiple phenotypes");
  

  if (no_covar > 0) ytz_prod = MathUtils::in_product(ytz[0], no_covar);
}

void DefaultStats::set_chromosome(int chr) {chromosome = chr;
  sex_covar = sex_covar_all; chrX_only_sex_covar = false; selection = sex_covar ? BaseInput::sm_KnownGender : BaseInput::sm_AllKnown;
  if (Utils::chr_type(chr) > 0) {
    if (chr == Utils::chrX_code || chr == Utils::chrZ_code) {
      sex_covar = sex_covar_all || sex_covar_chrX;
      chrX_only_sex_covar = sex_covar_chrX && !sex_covar_all;
      selection = sex_covar ? BaseInput::sm_KnownGender : BaseInput::sm_AllKnown;
    } else if (chr == Utils::chrY_code) selection = BaseInput::sm_MaleOnly;
    else if (chr == Utils::chrW_code) selection = BaseInput::sm_FemaleOnly;
  }
  if (pval_based && !read_N) {
    if (chr == Utils::chrX_code || chr == Utils::chrZ_code) orig_size = settings.geti("pval_size_chrX");        
    else if (chr == Utils::chrY_code || chr == Utils::chrW_code) orig_size = settings.geti("pval_size_chrY");        
    else orig_size = settings.geti("pval_size_main");
  }
}

void DefaultStats::process_snps() {snp_data->process(phenotype.begin());} /// only applied to first pheno, in case of multiple phenotypes

void DefaultStats::load_snp_info() {
  GeneStats::load_snp_info();  
  if (read_N) {
    BufferWindow<int> snpN = snp_data->get_snpN();
    orig_size = round(MathUtils::get_mean(snpN.view().begin(), no_var));
  }
}
 
PermutationGenerator* DefaultStats::perm_generator() {return new DefaultPermGen(this,samp_size);}
 
void DefaultStats::set_pcs(int index, bool init) {
  if (index > pc_data.size() || (!init && index == pc_data.size())) throw GeneException("computing principal components", "index for principal component block is out of range");
  if (index == pc_data.size()) pc_data.push_back(new PrinComp());
  princomp = pc_data[index];
  if (has_partition && !princomp->partitioned) princomp->set_partition(block_index.nrow());  
  pc_id = index;
}

void DefaultStats::compute_xty() {
  if (has_xty) return; has_xty = true;
  if (pval_based) load_pval_corr(xty, ProductType::Centered, false);
  else {
    snp_data->product_var(xty, phenotype);
    if (no_pheno > 1) snp_data->product_var(xty_multi, phenotype_multi, 0, 0, true); /// storing with SNPs as columns
  }
}

void DefaultStats::compute_xtz() {
  if (has_xtz) return; has_xtz = true; 
  snp_data->product_var(xtz, covar, 0, no_var, true);
}

float DefaultStats::compute_xtv(float* var, int snp, ProcessType::Type mode) {return snp_data->product_var(var, snp, mode);}
void DefaultStats::compute_xtv(Buffer<float>& buffer, Buffer<float>& var, bool transpose, ProcessType::Type mode) {
  snp_data->product_var(buffer, var, block_offset, block_size, transpose, mode);
}

void DefaultStats::compute_xtx(Buffer<float>& target, int offset, int total) {snp_data->product(target, offset, total);}
void DefaultStats::compute_xtx() {
  if (has_xtx) return; has_xtx = true;
  compute_xtx(xtx, block_offset, block_size);
}
  
void DefaultStats::compute_xtx_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode) {  
  snp_data->product(buffer, block_index(block_i,0), block_index(block_i,1), block_index(block_j,0), block_index(block_j,1), mode);
}

void DefaultStats::compute_xtx_block_skip(Buffer<float>& buffer, int block_i, int step_i, int block_j, int step_j, ProductType::Type mode) {
  snp_data->product_skip(buffer, block_index(block_i,0), step_i, block_index(block_i,1), block_index(block_j,0), step_j, block_index(block_j,1), mode);
}

void DefaultStats::compute_wtw_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode) {  
  BufferWindow<float> pc_i = princomp->get_pcs(block_i), pc_j = princomp->get_pcs(block_j); 
  MathUtils::matrix_prod(pc_i.view(), pc_j.view(), buffer, MathUtils::TransposeFirst);
  if (mode > ProductType::Centered) MathUtils::matrix_scale(buffer, 1/double(samp_size-1));
}

void DefaultStats::compute_wty(Buffer<float>& rotate) {
  if (princomp->computed(curr_block,1)) return; princomp->computed(curr_block,1) = true; 
  if (!princomp->computed(curr_block)) compute_pcs(rotate); 

  if (has_partition) {  
    BufferWindow<float> pcs = princomp->get_pcs(curr_block); 
    BufferWindow<float> write = princomp->get_wty(curr_block); 
    MathUtils::matrix_prod(pcs.view(), phenotype, write.view(), MathUtils::TransposeFirst);

    if (no_pheno > 1) {
      write = princomp->get_wty_multi(curr_block);  
      MathUtils::matrix_prod(phenotype_multi, pcs.view(), write.view(), MathUtils::TransposeFirst);
    }
  } else {
    MathUtils::matrix_prod(princomp->data, phenotype, princomp->wty, MathUtils::TransposeFirst);    
    if (no_pheno > 1) MathUtils::matrix_prod(phenotype_multi, princomp->data, princomp->wty_multi, MathUtils::TransposeFirst);    
  }
}

void DefaultStats::compute_wtz(Buffer<float>& rotate) {
  if (princomp->computed(curr_block,2)) return; princomp->computed(curr_block,2) = true; 
  if (!princomp->computed(curr_block)) compute_pcs(rotate); 
    
  if (has_partition) {  
    BufferWindow<float> pcs = princomp->get_pcs(curr_block); 
    BufferWindow<float> write = princomp->get_wtz(curr_block); 
    MathUtils::matrix_prod(covar, pcs.view(), write.view(), MathUtils::TransposeFirst);
  } else MathUtils::matrix_prod(covar, princomp->data, princomp->wtz, MathUtils::TransposeFirst);    
}

void DefaultStats::compute_pcs(Buffer<float>& rotate) {
  if (princomp->computed(curr_block)) return; 
  if (has_partition) {
    if (!princomp->computed(0) && MathUtils::sum(princomp->computed[0], princomp->computed.nrow()) == 0) {
      set_pc_buff(no_var);
      princomp->index.set_size(block_index.nrow(), 2, true);
      princomp->wty.set_size(no_var, 1);
      if (no_pheno > 1) princomp->wty_multi.set_size(no_pheno-1, no_var, true);
      princomp->wtz.set_size(no_covar, no_var);
    }

    int offset = (curr_block > 0) ? (princomp->index(curr_block-1,0) + princomp->index(curr_block-1,1)) : 0;
    int size = rotate.ncol();
    princomp->npcs += size;

    princomp->index(curr_block,0) = offset; princomp->index(curr_block,1) = size; 
    BufferWindow<float> block = princomp->data.window(offset, size);

    snp_data->data_projection(block.view(), rotate, block_offset, block_size);          
    MathUtils::normalize(block.view());   
  } else {
    set_pc_buff(rotate.ncol());
    snp_data->data_projection(princomp->data, rotate);    
    MathUtils::normalize(princomp->data);
    princomp->npcs = princomp->data.ncol();
  }
  princomp->computed(curr_block) = true;   
}   


BufferWindow<float> DefaultStats::get_xty() {
  if (!has_xty) compute_xty();
  return xty.window(block_offset, block_size, true);  
}

BufferWindow<float> DefaultStats::get_xty_multi() { /// stored with SNPs as columns, for windowing
  if (!has_xty) compute_xty();
  return xty_multi.window(block_offset, block_size);  
}     

BufferWindow<float> DefaultStats::get_xtz() {
  if (!has_xtz) compute_xtz();
  return xtz.window(block_offset, block_size);
}

BufferWindow<float> DefaultStats::get_wty(int pcid) {
  if (pcid >= 0) set_pcs(pcid);
  if (!princomp->computed(curr_block)) throw GeneException("processing gene data", "principal components for gene are not available");
  return get_wty(DUMMY_BUFFER);
}

BufferWindow<float> DefaultStats::get_wty(Buffer<float>& rotate) {
  if (!princomp->computed(curr_block,1)) compute_wty(rotate);
  if (has_partition) return princomp->get_wty(curr_block); 
  else return princomp->get_wty();
}

BufferWindow<float> DefaultStats::get_wty_multi(int pcid) {
  if (pcid >= 0) set_pcs(pcid);
  if (!princomp->computed(curr_block)) throw GeneException("processing gene data", "principal components for gene are not available");
  return get_wty_multi(DUMMY_BUFFER);
}

BufferWindow<float> DefaultStats::get_wty_multi(Buffer<float>& rotate) {
  if (!princomp->computed(curr_block,1)) compute_wty(rotate);
  if (has_partition) return princomp->get_wty_multi(curr_block); 
  else return princomp->get_wty_multi();
}

BufferWindow<float> DefaultStats::get_wtz(int pcid) {
  if (pcid >= 0) set_pcs(pcid);
  if (!princomp->computed(curr_block)) throw GeneException("processing gene data", "principal components for gene are not available");
  return get_wtz(DUMMY_BUFFER);
}
BufferWindow<float> DefaultStats::get_wtz(Buffer<float>& rotate) {
  if (!princomp->computed(curr_block,2)) compute_wtz(rotate);
  if (has_partition) return princomp->get_wtz(curr_block); 
  else return princomp->get_wtz();
}

BufferWindow<float> DefaultStats::get_xtx(bool invalidate) {
  if (!has_xtx) compute_xtx(); 
   if (invalidate && buffer_xtx) {
    xtx_buffer.assign(xtx);
    return xtx_buffer.window();
  }
  if (invalidate) has_xtx = false;
  return xtx.window();
}

BufferWindow<float> DefaultStats::get_capped_xtx(int max_size, bool invalidate) {
  if (max_size < block_size) {
    int offset = (block_size - max_size) / 2;
    compute_xtx(xtx_buffer, offset, max_size);
    return xtx_buffer.window();
  } else return get_xtx(invalidate);
}


BufferWindow<float> DefaultStats::get_pcs(int pcid) {
  if (pcid >= 0) set_pcs(pcid);
  if (!princomp->computed(curr_block)) throw GeneException("processing gene data", "principal components for gene are not available");
  return get_pcs(DUMMY_BUFFER);
}

BufferWindow<float> DefaultStats::get_pcs(Buffer<float>& rotate) {
  if (!princomp->computed(curr_block)) compute_pcs(rotate);
  if (has_partition) return princomp->get_pcs(curr_block);
  else return princomp->data.window();
}
void DefaultStats::prep_pcs(Buffer<float>& rotate) {if (!princomp->computed(curr_block,1)) compute_pcs(rotate);}

float DefaultStats::load_block_xtx(Buffer<float>& buffer, int block_i, int block_j, bool raw, bool skip, ProductType::Type mode) {float scale = 1;
  if (raw) {
    if (skip) {
      int step_i = block_stepsize(block_index(block_i,1)), step_j = block_stepsize(block_index(block_j,1));
      if (step_i > 1 || step_j > 1) {
        compute_xtx_block_skip(buffer, block_i, step_i, block_j, step_j, mode);
        scale = block_index(block_i,1)/float(buffer.nrow()) * block_index(block_j,1)/float(buffer.ncol());
      } else compute_xtx_block(buffer, block_i, block_j, mode);
    } else compute_xtx_block(buffer, block_i, block_j, mode);
  } else {
    if (!princomp->computed(block_i) || !princomp->computed(block_j)) throw GeneException("subdividing gene", "PCs for specified gene block are not available");
    compute_wtw_block(buffer, block_i, block_j, mode);
  } 
  return scale;
}

void DefaultStats::load_product(Buffer<float>& buffer, Buffer<float>& var, ProductType::Type mode, bool col_snp) {
  compute_xtv(buffer, var, col_snp, MathUtils::type_convert(mode, false));  
  if (mode > ProductType::Centered) {///var itself is assumed to be centered
    if (mode == ProductType::Correlation) {
      dat_buffer.set_size(block_size, 1);
      for (int i = 0; i < block_size; i++) dat_buffer(i) = sqrt(snp_data->get_variance(block_offset+i));
      for (int i_var = 0; i_var < var.ncol(); i_var++) {
        double curr = MathUtils::get_sd(var[i_var], samp_size) * (samp_size - 1);
        if (col_snp) for (int i_snp = 0; i_snp < block_size; i_snp++) buffer(i_var, i_snp) /= dat_buffer(i_snp) * curr;
        else for (int i_snp = 0; i_snp < block_size; i_snp++) buffer(i_snp, i_var) /= dat_buffer(i_snp) * curr;
      }  
    } else MathUtils::matrix_scale(dat_buffer, 1/double(samp_size-1));
  }
}

void DefaultStats::load_pairwise_product(Buffer<float>& buffer, Buffer<float>& var, ProductType::Type mode) {
  buffer.set_size(block_size,1);
  for (int i = 0; i < block_size; i++) buffer(i) = compute_xtv(var[i], block_offset+i, MathUtils::type_convert(mode));
  if (mode == ProductType::Covariance) MathUtils::matrix_scale(buffer, 1/double(samp_size-1));
}

void DefaultStats::load_pheno(Buffer<float>& buffer, bool blur) {
  buffer.assign(phenotype);
  if (blur && baseinput.binary_pheno) StatUtils::add_noise(buffer[0], samp_size, 0.005, true);
}

void DefaultStats::load_fitted(Buffer<float>& buffer, bool blur) {
  buffer.set_size(samp_size, 2, true);
  memcpy(buffer[1], phenotype[0], sizeof(float)*samp_size);
  if (no_covar > 0) {
    for(int i = 0; i < no_covar; i++) MathUtils::vec_product(covar[i], buffer[0], ytz(i) / (samp_size-1), samp_size, true);
    MathUtils::vec_product(buffer[0], buffer[1], -1, samp_size, true);
  }
  if (blur && baseinput.binary_pheno) StatUtils::add_noise(buffer[1], samp_size, 0.005, true, true);
}

void DefaultStats::load_data(Buffer<float>& buffer, bool norm) {snp_data->data_processed(buffer, block_offset, block_size, norm);}
void DefaultStats::load_interaction(Buffer<float>& buffer, ProcessType::Type process) {
  snp_data->data_product(buffer, interactor, block_offset, block_size);
  if (process == ProcessType::Centered) MathUtils::center(buffer);
  else if (process == ProcessType::Normalized) MathUtils::normalize(buffer);  
}

BufferWindow<double> DefaultStats::get_pval() {return snp_data->get_pval(block_offset, block_size);}

void DefaultStats::load_pval_corr(Buffer<float>& buffer, ProductType::Type mode, bool block) {
  int offset = block ? block_offset : 0; int size = block ? block_size : no_var;
  BufferWindow<double> wrapper = snp_data->get_pval(offset, size); Buffer<double>& pval = wrapper.view(); 
  ConvertPvalToT stat;

  int df = get_sampsize(false)-2; stat.set_param("df", df);
  buffer.set_size(size,1); 
  
  int* snpN = read_N ? (snp_data->get_snpN().view().begin() + offset) : 0;
  for (int i = 0; i < size; i++) {
    if (read_N) {
      df = snpN[i]-2;
      stat.set_param("df", df);
    }
    float t = stat.convert(pval(i)/2);
    buffer(i) = t / sqrt(df + t*t);
  }

  if (mode != ProductType::Correlation) {
    float base = (mode <= ProductType::Centered) ? samp_size - 1 : 1;
    for (int i = 0; i < size; i++) buffer(i) *= sqrt(snp_data->get_variance(offset+i)) * base;
  }
}  



void DefaultStats::write_vars(string prefix) {
  dump_count++;
  string suffix = ".gs" + Utils::num_string(dump_count) + ".genedata";
  dump(phenotype, prefix + "_pheno" + suffix);
  if (no_covar > 0) dump(covar, prefix + "_covar" + suffix);
}

void DefaultStats::write_gene(string prefix) {
  string file = prefix + "_" + baseinput.gene_info.name.get(gene_id) + ".gs" + Utils::num_string(dump_count) + ".genedata";
  load_data(dat_buffer);
  dump(dat_buffer, file);
}



void DefaultPermGen::run(BufferMode mode, int nperm) {
  prep_buffer(mode, nperm);
  Buffer<float>& buff = *(buffers[mode][Main]);
  for (int i = 0; i < nperm; i++) {
    if (no_covar > 0) generate_covar(buff[i]);
    else generate_plain(buff[i]);
  }
  if (no_covar > 0) MathUtils::matrix_prod(covariates, buff, *(buffers[mode][Covar]), MathUtils::TransposeFirst);
}

void DefaultPermGen::generate_plain(float* target) {
  long index, left = samp_size;
  perm_buffer.assign(phenotype.begin(), samp_size);

  while (left > 0) {
    index = _RNG() % left--;
    *(target++) = perm_buffer(index);
    perm_buffer(index) = perm_buffer(left);
  }   
}

void DefaultPermGen::generate_covar(float* target) {
  long index, left = samp_size;
  memcpy(target, phenotype[0], sizeof(float)*samp_size);
  perm_buffer.assign(phenotype[1], samp_size);

  float* curr = target;
  while (left > 0) {
    index = _RNG() % left--;
    *(curr++) += perm_buffer(index);
    perm_buffer(index) = perm_buffer(left);
  }

  MathUtils::normalize(target, samp_size);
}



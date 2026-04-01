/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genecorrelations.h"
#include "mathutils.h"


void GeneCorrelations::clear() {
  for (int i = 0; i < genes.size(); i++) delete genes[i];
  genes.clear();
  mat_buffer.reset();
}

void GeneCorrelations::compute(BaseBuffer<double>& corrs, CorrelationData* data) {
  if (inactive) return;
  if (!data || data->gene_id() < 0) {_LOG.error("computing gene correlations") << "invalid gene correlation object" << endl; die();}
  if (!data->check_stdev()) {discard(data); throw GeneException("computing gene correlations", "variance for gene is zero or negative");}
  
  GeneLocation loc = data->get_location(); 
  while (!genes.empty() && !loc.in_range(genes.front()->get_location(), max_range)) {
    discard(genes.front());
    genes.pop_front();
  }

  if (!genes.empty() && data->gene_id() <= genes.back()->gene_id()) {_LOG.error("computing gene correlations") << "gene correlation objects are not ordered" << endl; die();}
  int offset = !genes.empty() ? genes[0]->gene_id() : 0, range = !genes.empty() ? data->gene_id() - offset : 0;

  corrs.assign_value(0, range); 
  for (int i = 0; i < genes.size(); i++) corrs[genes[i]->gene_id()-offset] = data->get_correlation(genes[i]);
  genes.push_back(data);
}

void GeneCorrelations::discard(CorrelationData* data) {
  while (GeneStats::DataBlock* curr = data->eject_data()) stats->insert_data(curr);
  delete data;      
}


const double CorrelationUtils::powersum_coefficients[3][5] = {///coefficients for P = 1 - 4, plus scaling to correlation
  {3,  1,   0,  0,      4},
  {45, 60,  8,  0,      113},
  {35, 105, 56, 4,      200}
};                              

double CorrelationUtils::powersum_sum(double* rsq, int pow) {double out = 0;
  for (int i = 0; i < pow; i++) out += rsq[i]*powersum_coefficients[pow-2][i];
  return out/powersum_coefficients[pow-2][4];
}


double CorrelationData::get_correlation(CorrelationData* other) {
  check_readiness(other);

  if (!check_stdev() || !other->check_stdev()) return 0;
  double corr = get_covariance_core(other) / sd / other->get_stdev();
  return MathUtils::clamp01(corr);
}

void CorrelationData::check_readiness(CorrelationData* other) {
  if (!has_data() || !other->has_data()) {_LOG.error("computing gene correlations") << "gene correlation objects contains no source data" << endl; die();}
  if (!is_compatible(other)) {_LOG.error("computing gene correlations") << "encountered incompatible gene correlation data objects" << endl; die();}
}

int CorrelationUtils::top_power = 0;
double CorrelationUtils::top_corr_weight = 0.5;
double CorrelationUtils::top_corr_pow[3] = {1,1,1};

void CorrelationUtils::set_top_param(Settings& settings) {
  top_power = 4;
  top_corr_weight = settings.getn("gene_top_corr_weight");
  top_corr_pow[0] = settings.getn("gene_top_corr_p1");
  top_corr_pow[1] = settings.getn("gene_top_corr_p2");
  top_corr_pow[2] = settings.getn("gene_top_corr_p3");    
}


Buffer<float> CoreCorrelationData::corr_buffer;
Buffer<double> BlockCorrelationData::corr_buffer;
Buffer<double> MultiCorrelationData::corr_buffer;

double CorrelationData::get_covariance(CorrelationData* other) {
  check_readiness(other);
  return get_covariance_core(other);
}

bool CoreCorrelationData::set_block(int id) {
  if (check_index(id)) {
    if (!block_index.empty()) {
      curr_block = id;
      curr_offset = block_index(id,0); 
      curr_size = block_index(id,1);    
    } else return set_block(no_blocks);
  } else {
    curr_block = 0;
    curr_offset = 0;
    curr_size = data->size();
  }
  return id < no_blocks;
}

void CoreCorrelationData::set_partition() {
  if (no_blocks <= 1) block_index.set_size(0,2);
  else {
    block_index.assign(get_partition());
    if (block_index.nrow() != no_blocks) throw GeneException("loading partition", "gene partition is invalid");
  }
}

void CoreCorrelationData::compute_variance(Buffer<float>& xtx, int block, double scale_full_block) {
  if (!check_index(block)) return;
  double var = CorrelationUtils::variance_chisq(xtx, false);
  add_variance(var*scale_full_block, block);
}

int RawCorrelationData::get_step_size(int dim) {
  if (dim < 2*min_step_dim) return 1;
  return min(dim/min_step_dim, max_step_size);
}

double RawCorrelationData::prep_covariance(CorrelationData* other_base, Buffer<float>& target) {RawCorrelationData* other = convert<RawCorrelationData>(other_base);
  Gene *gi = other->data->get_gene(), *gj = data->get_gene();
  int &off_i = other->curr_offset, &off_j = curr_offset;
  int &total_i = other->curr_size, &total_j = curr_size;
  int step_i = other->get_step_size(total_i), step_j = get_step_size(total_j);

  gi->product_skip(target, gj, off_i, step_i, total_i, off_j, step_j, total_j, ProductType::Correlation);  
  return total_i/double(target.nrow()) * total_j/double(target.ncol());
}

double RawCorrelationData::get_covariance_core(CorrelationData* other_base) {
  double scale_full_block = prep_covariance(other_base, corr_buffer);
  return CorrelationUtils::covariance_chisq(corr_buffer, false)*scale_full_block;
}

void RawCorrelationData::compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j) {
  if (!check_index(block_i,block_j)) return;
  block_covar(block_j,block_i) = CorrelationUtils::covariance_chisq(corrs, false)*scale_full_block;
}

 
double PowerRawCorrelationData::get_covariance_core(CorrelationData* other_base) {PowerRawCorrelationData* other = convert<PowerRawCorrelationData>(other_base);
  double scale_full_block = prep_covariance(other_base, corr_buffer);
  double *sd_i = other->inverse_sds[other->curr_block], *sd_j = inverse_sds[curr_block];
  
  double merged_corr = CorrelationUtils::covariance_top(corr_buffer, sd_i, sd_j, false);  
  return merged_corr*scale_full_block;
}

void PowerRawCorrelationData::compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j) {
  if (!check_index(block_i,block_j)) return;
  double merged_corr = CorrelationUtils::covariance_top(corrs, inverse_sds[block_i], inverse_sds[block_j], false);
  block_covar(block_j,block_i) = merged_corr*scale_full_block;
}

void PowerRawCorrelationData::compute_variance(Buffer<float>& xtx, int block, double scale_full_block) {
  if (!check_index(block)) return;
  if (no_blocks == 1) block = 0;  

  bool valid = CorrelationUtils::variance_top_invert(xtx, inverse_sds[block], false, scale_full_block);
  int merged = valid ? 1 : -1;
  add_variance(merged, block);
}

double PCCorrelationData::get_covariance_core(CorrelationData* other_base) {PCCorrelationData* other = convert<PCCorrelationData>(other_base);
  BufferWindow<float> buff_i = other->get_buffer(), buff_j = get_buffer();
  MathUtils::matrix_prod(buff_i.view(), buff_j.view(), corr_buffer, MathUtils::TransposeFirst);  
  
  double scale = buff_i.view().nrow() - 1;
  return covariance_function(corr_buffer, 1/scale);
}

void PCCorrelationData::compute_covariance(Buffer<float>& corrs, double scale_full_block, int block_i, int block_j) {
  if (!check_index(block_i,block_j)) return;
  block_covar(block_j,block_i) = covariance_function(corrs)*scale_full_block; 
}

double PCCorrelationData::covariance_function(Buffer<float>& corrs, double scale_to_corr) {return CorrelationUtils::covariance_chisq(corrs, false, scale_to_corr);}
double PowerPCCorrelationData::covariance_function(Buffer<float>& corrs, double scale_to_corr) {
  pair<double,double> vars = CorrelationUtils::covariance_powersum(corrs, false, power, scale_to_corr);
  return (1-corr_weight)*vars.first + corr_weight*vars.second;
}

bool PowerPCCorrelationData::compatible_settings(CorrelationData* other_base) {PowerPCCorrelationData* other = convert<PowerPCCorrelationData>(other_base);  
  if (other) {
    return (corr_weight == other->corr_weight) && (power == other->power);    
  } else return false;
} 

void CompoundCorrelationData::compute_variance(Buffer<double>& corrs, bool process_blocks) {
  if (aggregate_type == Aggregator::Mean) {
    if (process_blocks) set_stdev(CorrelationUtils::variance_chisq(corrs, true)); 
    inverse_sds[0] = sd > 0 ? 1/sd : 1; 
  } else if (aggregate_type == Aggregator::Top) {bool valid;
    if (process_blocks) valid = CorrelationUtils::variance_top_invert(corrs, inverse_sds, true);
    else {valid = sd > 0; inverse_sds[0] = inverse_sds[1] = 1;}
    set_stdev(valid ? 1 : -1);
  }
}

double CompoundCorrelationData::compute_covariance(Buffer<double>& corrs, double* inv_sd_i, double* inv_sd_j, Aggregator::AggregateType type, bool to_corr) { /// Top automatically normalizes to correlation already
  switch (type) {
    case Aggregator::Mean: {  
      double out = CorrelationUtils::covariance_chisq(corrs, true);
      return to_corr ? out * *inv_sd_i * *inv_sd_j : out; }
    case Aggregator::Top: return CorrelationUtils::covariance_top(corrs, inv_sd_i, inv_sd_j, true);
    default: return 0;
  }
}

bool CompoundCorrelationData::compatible_settings(CorrelationData* other_base) {CompoundCorrelationData* other = convert<CompoundCorrelationData>(other_base);  
  if (other) {
    return (aggregate_type == other->aggregate_type);    
  } else return false;
}

void BlockCorrelationData::block_init() {
  data->set_partition();
  no_blocks = data->no_blocks;
  block_corr = &(data->block_covar);
  block_index = &(data->block_index);  
}

bool BlockCorrelationData::prep_variance() {
  block_sds.set_size(no_blocks,max(no_blocks,1)); set_stdev(-1);
  if (no_blocks <= 1) {sd = block_sds(0) = data->get_stdev(); return false;}
  Buffer<double>& cov = *block_corr; 
  for (int i = 0; i < no_blocks; i++) {
    block_sds(i) = cov(i,i) > 0 ? sqrt(cov(i,i)) : -1;
    if (block_sds(i) < 0) {sd = -1; return false;}
  }
  
  MathUtils::cov_to_cor(cov);
  return true;
}

double BlockCorrelationData::prep_covariance(Buffer<double>& buffer, BlockCorrelationData* other) { 
  int &size_i = no_blocks, &size_j = other->no_blocks;
  CoreCorrelationData *gene_i = data, *gene_j = other->data;

  if (size_i < 1 && size_j < 1) return 0;
  if (size_i == 1 && size_j == 1) return gene_i->get_covariance_core(gene_j);

  Buffer<double> &sds_i = block_sds, &sds_j = other->block_sds;
  buffer.set_size(size_i, size_j);
  for (int i = 0; i < size_i; i++) {
    gene_i->set_block(i);        
    for (int j = 0; j < size_j; j++) {
      gene_j->set_block(j);
      double cov = gene_i->get_covariance(gene_j);
      buffer(i,j) = cov / (sds_i(i)*sds_j(j));
    }
  }

  gene_i->set_block(-1); gene_j->set_block(-1);
  return -999;
}

double BlockCorrelationData::get_covariance_core(CorrelationData* other_base) {BlockCorrelationData* other = convert<BlockCorrelationData>(other_base);
  double corr = prep_covariance(corr_buffer, other);
  if (corr < -1) corr = compute_covariance(corr_buffer, inverse_sds, other->inverse_sds);
  return corr;
}

void MultiCorrelationData::set_corrs(Buffer<double>& corrs, Buffer<float>& perms, int use) {bool compute = (&corrs == &engine_corrs);
  if (perms.nrow() < 50) throw GeneException("multi-model variance", "insufficient permutations to compute model correlations");
  if (perms.ncol() != data.size()) throw GeneException("multi-model variance", "number of models does not match size of permutation matrix");
  if (corrs.ncol() != data.size() && !compute) throw GeneException("multi-model variance", "number of models does not match size of correlation matrix");
      
  perm_buffer.swap(perms); no_perm = (use <= 0 || use > perm_buffer.nrow()) ? perm_buffer.nrow() : use;

  if (!fixed_mode) {
    MathUtils::normalize(perm_buffer, no_perm);

    if (compute) {
      try {MathUtils::cov_matrix(perm_buffer, engine_corrs, true, no_perm);}
      catch (const MathException& me) {throw GeneException("multi-model variance", "model variance is zero");}
    } else engine_corrs.assign(corrs);

    compute_variance(engine_corrs);
  } else set_stdev(1);
}
  
bool MultiCorrelationData::prep_variance() {set_stdev(-1);
  if (data.size() <= 1) {sd = !data.empty() ? data[0]->get_stdev() : -1; return false;}
  return true;
}

double MultiCorrelationData::get_covariance_core(CorrelationData* other_base) {
  MultiCorrelationData* other = convert<MultiCorrelationData>(other_base);

  if (!fixed_mode) {
    int nrow = min(no_perm, other->no_perm), ncol = data.size();
    bool renorm = no_perm != other->no_perm;

    corr_buffer.set_square(ncol);
    for (int k = 0; k < ncol; k++) {
      corr_buffer(k,k) = data[k]->get_correlation(other->data[k]);
      for (int l = 0; l < ncol; l++) {if (l == k) continue;
        double& corr = corr_buffer(k,l);
        if (renorm) corr = MathUtils::get_cov(other->perm_buffer[k], perm_buffer[l], nrow, true);
        else corr = MathUtils::in_product(other->perm_buffer[k], perm_buffer[l], nrow) / (nrow-1);
      }
    }

    return compute_covariance(corr_buffer, inverse_sds, other->inverse_sds);
  } else {
    int ncol = data.size(); double corr = 0;
    for (int k = 0; k < ncol; k++) corr += data[k]->get_correlation(other->data[k]);
    return corr / ncol;
  }
}  

bool MultiCorrelationData::is_compatible(CorrelationData* other_base) {
  MultiCorrelationData* other = convert<MultiCorrelationData>(other_base);
  if (!other || data.size() != other->data.size() || aggregate_type != other->type()) return false;
  for (int i = 0; i < data.size(); i++) {
    if (!data[i]->is_compatible(other->data[i])) return false;
  }
  return true;
}

bool MultiCorrelationData::has_data() {
  for (int i = 0; i < data.size(); i++) {
    if (!data[i]->has_data()) return false;
  }
  return data.size() > 0;
}

GeneStats::DataBlock* MultiCorrelationData::eject_data() {
  GeneStats::DataBlock* out = 0;
  while (out == 0 && !data.empty()) {
    out = data.back()->eject_data();
    delete data.back();
    data.pop_back();
  }
  return out;    
}



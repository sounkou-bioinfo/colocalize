/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "permutation.h"

#include "exceptions.h"

#define MVN_PRUNE_MINVAL 0.01

const int PermutationEngine::core_max;
const int PermutationEngine::transient_max;

bool Permutation::set_data(GeneStats* gs) {
  data = gs->convert<DefaultStats>();
  data_id = data->get_id(); 

  delete generator;
  generator = data->perm_generator();
  
  return false;
}

void Permutation::check_failure(int total) {
  if (no_failed < 10 || no_failed/float(total) < 0.5) return;
  throw GeneException("permutation error", "unable to complete permutations for this gene");  
}

bool PermutationEngine::set_data(GeneStats* gs) {
  if (gs->equals(data_id)) return false; Permutation::set_data(gs);

  generator->set_buffer(core_perm, core_covar, PermutationGenerator::Core);
  generator->set_buffer(transient_perm, transient_covar, PermutationGenerator::Transient);  

  generator->run(PermutationGenerator::Core, min(max_perm, core_max));
  if (max_perm > core_max) generator->run(PermutationGenerator::Transient, min(max_perm-core_max, transient_max));

  return true;
}

int PermutationEngine::get_needed(int curr) {
  for (int i = 0; i < partitioning.size(); i++) {
    if (curr < partitioning[i]) return partitioning[i];
  }
  return max_perm;
}  

bool PermutationEngine::add_partition(int offset) {
  if (offset > max_perm || offset <= 0) return false;
  for (int i = 0; i < partitioning.size(); i++) {
    if (offset == partitioning[i]) return false;
    if (offset < partitioning[i]) {
      partitioning.insert(partitioning.begin()+i, offset); 
      break;
    }
  }
  if (offset < min_perm) min_perm = offset;
  return true;  
}

bool PermutationEngine::next() {
  no_failed = 0; block_curr = 0;
  if (curr_index < 0) curr_index = 0;
  else curr_index += curr_nperm;

  if ((curr_index >= max_perm) || (adaptive && sign_count[0] >= adap_count) || empty) {curr_nperm = 0; return false;}
  if (curr_index < core_max) curr_nperm = min(get_needed(curr_index), core_max) - curr_index;
  else {
    curr_nperm = min(get_needed(curr_index)-curr_index, transient_max);
    if (curr_index > core_max || curr_nperm > transient_perm.ncol()) generator->run(PermutationGenerator::Transient, curr_nperm);
  }

  return true;
}

Buffer<float>& PermutationEngine::get_product(Buffer<float>& target, Buffer<float>& data) {
  BufferWindow<float> perm;
  if (curr_index < core_max) perm = core_perm.window(curr_index, curr_nperm);
  else perm = transient_perm.window(0, curr_nperm);
  return MathUtils::matrix_prod(data, perm.view(), target, MathUtils::TransposeFirst);
}

BufferWindow<float> PermutationEngine::get_covar() {
  if (curr_index < core_max) return core_covar.window(curr_index, curr_nperm);
  else return transient_covar.window(0, curr_nperm);
}

PermutationBlock* PermutationEngine::get_block(bool from_start) {
  if (from_start) {rewind(); next();}

  PermutationBlock*& block = block_curr > 0 ? add_block : main_block; 
  if (block == 0) block = new PermutationBlock();

  bool use_core = curr_index < core_max;
  Buffer<float>& perm = use_core ? core_perm : transient_perm;
  int offset = use_core ? curr_index : 0;

  if (generator->has_covar()) block->set_perm(perm, (use_core ? core_covar : transient_covar), offset, curr_nperm);
  else block->set_perm(perm, offset, curr_nperm);

  block_curr++;
  return block;  
}

void PermutationEngine::insert(PermutationBlock* input) {
  input->rescale();
  if (input != main_block) main_block->merge(input);
}

void PermutationEngine::process_block(PermutationBlock* add) {if (add == 0) add = main_block;
  sign_count[0] += add->nsign();
  sign_count[1] += add->nperm();
}       

void PermutationBlock::transform(StatConverter* conv) {
  set_converter(conv);
  obs_stat = convert(obs_stat);
  curr_perm = perm_stats.begin();
  for (float* end = curr_perm + sign_count[1]; curr_perm < end; curr_perm++) *curr_perm = convert(*curr_perm);
  sign_count[0] = -1;
}

void PermutationBlock::merge(PermutationBlock* input) {
  if (input->sign_count[1] == 0) return;
  if (sign_count[1] > input->sign_count[1] || sign_count[1] == 0) sign_count[1] = input->sign_count[1];
  proc->obs(input->obs_stat, obs_stat);
  curr_perm = perm_stats.begin();  
  for (float *curr_input = input->perm_stats.begin(), *end = curr_input + sign_count[1]; curr_input < end; curr_input++) proc->perm(*curr_input, curr_perm);   
  sign_count[0] = -1; empty = true;
}   

double PermutationBlock::convert(double stat) {
  if (!to_pval) return stat;
  try {return logp.convert(to_pval->convert(stat));}
  catch (const StatException& se) {return logp.convert(0.5);}
}

PermutationBlock* PermutationBlock::copy_block() {
  PermutationBlock* out = new PermutationBlock();
  out->set_perm(permutations, covar);
  return out;    
}


void PermutationBlock::rescale() {
  double sd = MathUtils::get_sd(perm_stats[0], size(false));
  obs_stat /= sd;
  MathUtils::vec_scale(perm_stats[0], 1/sd, size(false)); 
}


Buffer<float>& PermutationBlock::get_product(Buffer<float>& target, Buffer<float>& data) {return MathUtils::matrix_prod(data, permutations.view(), target, MathUtils::TransposeFirst);}
BufferWindow<float> PermutationBlock::get_covar() {return covar;}


bool SimulationMVN::set_data(GeneStats* gs) {
  if (gs->equals(data_id)) return false; 
  data = gs->convert<DefaultStats>(); data_id = data->get_id(); 
  no_covar = load_covar ? data->ncovar() : 0; samp_size = data->get_sampsize();
  if (no_covar > 0) {
    covariates.assign(data->get_covar());
    data->load_fitted(pheno_fitted);
    resid_scale = MathUtils::in_product(pheno_fitted[1], samp_size);
    means.set_size(1,no_covar);
    for (int i = 0; i < no_covar; i++) means(i) = MathUtils::in_product(pheno_fitted[0], covariates[i], samp_size);
  }
  return true;
}

void SimulationMVN::simulate(int needed) {
  if (needed > max_dim) needed = max_dim;
  if (needed <= avail_dim || needed <= 0) return;
  
  int add = needed-avail_dim;
  if (avail_dim == 0) sim_buffer.set_size(tot_perm, needed);
  else sim_buffer.add_cols(add);
  
  StatUtils::generate_normal(sim_buffer[avail_dim], add*tot_perm);
  avail_dim = needed;
}

int SimulationMVN::set_input_covariance(Buffer<double>& covar) {
  projection.set_size(covar.nrow(), covar.ncol()); float* write = projection.begin();
  for (double *curr = covar.begin(), *end = covar.end(); curr < end; curr++) *(write++) = *curr;
  
  return set_input_covariance(projection);
}

int SimulationMVN::set_input_covariance(Buffer<float>& covar) {no_covar = 0; return invert_covariance(covar);}
int SimulationMVN::set_input_covariance(Buffer<float>& input, Buffer<float>& covar) {
  if (no_covar > 0) {int nsnps = covar.nrow(), nvar = no_covar + nsnps;
    if (input.ncol() != nsnps) throw GeneException("simulation error", "dimension of data covariance matrix does not match dimension of input data");    
    if (input.nrow() != samp_size) throw GeneException("simulation error", "dimension of input data does not match sample size");    

    projection.set_square(nvar);
    for (int i = 0; i < no_covar; i++) {
      projection(i,i) = samp_size-1;
      for (int j = i+1; j < no_covar; j++) projection(j,i) = 0;
    }
    projection.matrix_sub(no_covar, nsnps) = covar.matrix();
    projection.matrix_sub(no_covar, nsnps, 0, no_covar).noalias() = input.matrix().transpose() * covariates.matrix();

    MathUtils::matrix_scale(projection, resid_scale / (samp_size-1));
    
    means.total_cols(nvar);
    for (int i = 0; i < nsnps; i++) means(no_covar+i) = MathUtils::in_product(pheno_fitted[0], input[i], samp_size);
    
    return invert_covariance(projection);
  } else return invert_covariance(covar);
}

int SimulationMVN::invert_covariance(Buffer<float>& covar) {int nvar = covar.nrow(); 
  if (nvar <= 0) throw GeneException("simulation error", "dimension of covariance matrix is zero");    
  if (nvar != covar.ncol()) throw GeneException("simulation error", "covariance matrix is not square");      
  
  if (nvar == 1) {
    projection.set_size(1,1); projection(0) = 1;
    curr_dim = 1;
    return 1;
  }
  
  Map<MatrixXf> wrap(covar.begin(), nvar, nvar);
  SelfAdjointEigenSolver<MatrixXf> eig(wrap);
   
  float total = 0; const float* ev = eig.eigenvalues().data();
  for (int i = 0; i < nvar; i++) {total += ev[i];}
  
  float cutoff = MVN_PRUNE_MINVAL * total / nvar; 
  if (prune_perc < 1) {
    buffer.assign(eig.eigenvalues().data(), nvar);
    sort(buffer[0], buffer[1]);    
    float curr = 0, thresh = total * (1-prune_perc); 
    for (int i = 0; i < nvar; i++) {
      curr += buffer(i);
      if (curr > thresh) {
        if (i > 0 && cutoff < buffer(i-1)) cutoff = (buffer(i-1)+buffer(i))/2;
        break;
      }
    }
  }

  curr_dim = nvar; buffer.set_size(nvar, 1);
  for (int i = nvar-1; i >= 0; i--) {
    if (ev[i] <= cutoff) {buffer(i) = 0; curr_dim--;}
    else buffer(i) = sqrt(ev[i]);
  }
  
  if (curr_dim <= 0) {
    if (curr_dim == 0 && ev[nvar-1] > MVN_PRUNE_MINVAL * total / nvar) {
      curr_dim = 1;
      buffer(nvar-1) = sqrt(ev[nvar-1]);
    } else throw GeneException("running gene analysis", "unable to eigendecompose SNP covariance matrix");  
  }
  if (curr_dim > max_dim) curr_dim = max_dim;

  projection.set_size(nvar, curr_dim); int drop = nvar - curr_dim;
  projection.matrix().noalias() = eig.eigenvectors().block(0,drop,nvar,curr_dim) * Map<Matrix<float,Dynamic,1> >(buffer.begin()+drop, curr_dim, 1).asDiagonal();

  simulate(curr_dim);
  curr_nvar = nvar - no_covar;
  return curr_dim;
}

bool SimulationMVN::next() {
  if (curr_index < 0) curr_index = 0;
  else curr_index += curr_nperm;

  curr_nperm = min(tot_perm - curr_index, chunk_size);
  if (curr_nperm <= 0) {curr_nperm = 0; return false;} 
  return true;
}

Buffer<float>& SimulationMVN::get_sims(Buffer<float>& target) {
  target.set_size(curr_nvar, curr_nperm);
  if (no_covar > 0) {
    target.matrix().noalias() = projection.matrix_sub(no_covar, curr_nvar, 0, curr_dim) * sim_buffer.matrix_sub(curr_index, curr_nperm, 0, curr_dim).transpose();
    for (int i = 0; i < curr_nvar; i++) shift(target.begin()+i, means(no_covar+i), curr_nvar, curr_nperm);
  } else target.matrix().noalias() = projection.matrix() * sim_buffer.matrix_sub(curr_index, curr_nperm, 0, curr_dim).transpose();
  return target;
}

BufferWindow<float> SimulationMVN::get_covar() {
  covar_sims.set_size(no_covar, curr_nperm);
  if (no_covar > 0) {
    covar_sims.matrix().noalias() = projection.matrix_sub(0, no_covar, 0, curr_dim) * sim_buffer.matrix_sub(curr_index, curr_nperm, 0, curr_dim).transpose();
    for (int i = 0; i < no_covar; i++) shift(covar_sims.begin()+i, means(i), no_covar, curr_nperm);    
  }
  return covar_sims.window();
}
   
void MultiPhenoFWER::set_values(int no_values) {
  if (perm_buffer.length() != no_values) {
    boost::random::mt19937 rng_engine; rng_engine.seed(no_values);  
    boost::normal_distribution<float> distn(0, 1);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<float> > norm_gen(rng_engine, distn);

    perm_buffer.set_size(no_values, 1);
    for (float *curr = perm_buffer.begin(), *end = perm_buffer.end(); curr < end; curr++) *curr = norm_gen();  
  }
}  
   
void MultiPhenoFWER::generate(string outfile, int no_perm) {
  int N = data.get_sampsize(), no_pheno = data.npheno();
  set_values(N + no_perm); /// 1 extra value to simplify x_sum updates

  float df_m = 1, df_e = N - df_m - 1, scale = df_e / df_m;
  ConvertFToPval converter(df_m, df_e);
  min_p.set_size(no_perm, 1);

  float *curr = perm_buffer.begin(); float x_sum = MathUtils::sum(curr, N), xsq_sum = MathUtils::in_product(curr, N); 
  Buffer<float> &pheno_one = data.get_pheno(), &pheno_rest = data.get_pheno_multi(); Buffer<float> xty_rest; 
  for (int i = 0; i < no_perm; curr++, i++) {
    float sq_max = MathUtils::in_product(curr, pheno_one.begin(), N); sq_max *= sq_max;

    BufferWindow<float> x_vec(curr, N);    
    MathUtils::matrix_prod(pheno_rest, x_vec.view(), xty_rest, MathUtils::TransposeFirst);    
    for (int p = 0; p < (no_pheno-1); p++) {
      float sq = xty_rest(p) * xty_rest(p);
      if (sq > sq_max) sq_max = sq;
    }

    float x_var = (xsq_sum - x_sum*x_sum/N) / (N-1);
    float ssm = sq_max / (N-1) / x_var, sse = N - ssm - 1;
    min_p(i) = converter.convert(scale * (ssm/sse));

    x_sum += curr[N] - curr[0];
    xsq_sum += (curr[N]*curr[N]) - (curr[0]*curr[0]);  
  }
  
  ofstream out(outfile.c_str(), ios::trunc | ios::binary);
  out << N << endl << no_perm << endl << min_p;
} 


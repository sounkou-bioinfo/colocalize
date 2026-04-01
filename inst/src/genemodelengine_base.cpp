/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genemodelengine_base.h"

double GeneModelEngine::run(GeneStats* input, bool var, Permutation* perm) { try {
  set_data(input); needs_var = var; 
  sample_size = input->get_sampsize(); no_snps = input->nvar(); no_param = no_snps;   
  has_perm = (perm != 0); permutation = perm;
  no_covar = input->ncovar(); has_covar = no_covar > 0; has_interactor = (no_covar == 1) && input->has_interactor();

  init_gene();
  run_core();

  return pval(0);
} catch (const StatException& se) {throw GeneException("analyzing gene", "unable to convert test statistic to p-value", se.get_msg());} }

int CoreEngine::eigen_count(Buffer<float>& cov, double prune_prop) {
  SelfAdjointEigenSolver<Matrix<float,Dynamic,Dynamic> > eig(cov.matrix(), EigenvaluesOnly);
  vector<double> eval(eig.eigenvalues().data(), eig.eigenvalues().data()+cov.ncol());
  sort(eval.begin(), eval.end(), greater<double>());

  double cutoff = 0; int eff_size = 0;
  for (int i = 0; i < eval.size(); i++) cutoff += eval[i]; 
  cutoff *= prune_prop;

  for (int i = 0; i < eval.size(); i++) {
    eff_size++;
    cutoff -= eval[i];
    if (cutoff < 0) break;
  }  

  return max(eff_size, 1);
}

void CoreEngine::init_gene() {
  GeneModelEngine::init_gene();  

  int change_code = get_data()->set_current(gene_id, block_id, block_count);
  data_changed = change_code > 0;
  gene_changed = change_code > 1;
}

ModelInfoBlock CoreEngine::make_info(GeneStats* gs, ModelState::Model model) {ModelInfoBlock out; int index = 0;
  if (gs->ncovar() == 1 && gs->has_interactor()) {
    out.info.push_back(new ModelOutputInfo(model, ModelState::FullAnalysis, index++));   
    out.info.push_back(new ModelOutputInfo(model, ModelState::InteractionAnalysis, index++));
  } 
  out.info.push_back(new ModelOutputInfo(model, ModelState::MainAnalysis, index++));      
  out.detect_state(); return out;
}

void CoreEngine::set_variance(CoreCorrelationData* cd, double var) {
  cd->add_variance(var, block_id);  
  compute_covariance(cd);    
}

void CoreEngine::compute_variance(CoreCorrelationData* cd, Buffer<float>& xtx, double scale) {
  cd->compute_variance(xtx, block_id, scale*scale);  
  compute_covariance(cd);    
}

void CoreEngine::compute_covariance(CoreCorrelationData* cd) {
  for (int i = 0; i < block_id; i++) {
    double scale_full_block = prep_block_corr(mat_buffer, i, block_id);
    cd->compute_covariance(mat_buffer, scale_full_block, i, block_id);
  }
}

PrinCompEngine::~PrinCompEngine() {delete covar;}
void PrinCompEngine::clear() {CoreEngine::clear(); delete covar; covar = 0;}  

void PrinCompEngine::init_gene() {
  CoreEngine::init_gene();
  data->set_pcs(pcid, true);
  if (block_id == 0) {delete covar; covar = 0;}
}

double PrinCompEngine::prep_block_corr(Buffer<float>& target, int block_i, int block_j) {
  data->load_block_cor(target, block_i, block_j, false);
  return 1;
}          

CorrelationData* PrinCompEngine::corr_data() {
  CoreCorrelationData* out = covar; covar = 0;
  return out;    
}

void PrinCompEngine::prenormalize(Buffer<float>& mat, Buffer<float>& scale) {
  int dim = mat.nrow(); scale.set_size(dim, 1);
  for (int i = 0; i < dim; i++) {
    scale(i) = sqrt(1/mat(i,i));
    mat(i,i) = 1;
    
    for (int j = 0; j < i; j++) {
      mat(j,i) *= scale(i)*scale(j);
      mat(i,j) = mat(j,i);
    }
  }  
}

int PrinCompEngine::prune_invert(Buffer<float>& mat, Buffer<float>& buffer, bool root) {int size = mat.nrow();
  if (size == 1) {mat(0) = root ? 1/sqrt(mat(0)) : 1/mat(0); return 1;}

  Map<MatrixXf> wrap(mat.begin(), size, size);
  SelfAdjointEigenSolver<MatrixXf> eig(wrap);
 
  float total = 0; const float* ev = eig.eigenvalues().data();
  for (int i = 0; i < size; i++) {total += ev[i];}

  float cutoff = prune_min_val * total / size; int drop = prune_drop;
  if (prune_var_prop < 1 || drop > 0) {
    buffer.assign(eig.eigenvalues().data(), size);
    sort(buffer[0], buffer[1]);
    if (prune_var_prop < 1) {
      float curr = 0, thresh = total * (1-prune_var_prop); int i_eig = 0;
      for (; i_eig < size; i_eig++) {
        curr += buffer(i_eig);
        if (curr > thresh) break;
      }
      if (i_eig > 0) drop = max(drop, i_eig);
    }

    drop = min(drop, size-1);
    if (drop > 0) {
      int retain = 0;
      if (buffer(drop-1) > cutoff) {
        cutoff = buffer(drop-1);
        for (int i = drop; i < size; i++) {if (buffer(i) == cutoff) retain++; else break;}
      }
      drop = retain;
    }
  }

  int dim = size; buffer.set_size(size, 1);
  for (int i = size-1; i >= 0; i--) {
    if (ev[i] <= cutoff) {
      if (drop > 0 && ev[i] == cutoff) drop--;
      else {buffer(i) = 0; dim--; continue;}
    }
    buffer(i) = root ? 1/sqrt(ev[i]) : 1/ev[i];
  }
  if (dim <= 0) throw GeneException("running gene analysis", "unable to invert SNP covariance matrix");  

  Map<Matrix<float,Dynamic,Dynamic> > out_map(mat[0], size, size);
  if (root) {
    out_map.noalias() = eig.eigenvectors() * Map<Matrix<float,Dynamic,1> >(buffer[0], size, 1).asDiagonal();
    if (dim < size) mat.shift_cols(dim-size);
  } else out_map.noalias() = eig.eigenvectors() * Map<Matrix<float,Dynamic,1> >(buffer[0], size, 1).asDiagonal() * eig.eigenvectors().transpose();

  return dim;
}


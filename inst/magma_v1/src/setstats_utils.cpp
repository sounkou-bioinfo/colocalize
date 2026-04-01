/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "setstats.h"
#include "setstats_utils.h"
#include "exceptions.h"
#include "utils.h"

using namespace SetStatsUtils;

vector<string>& InternalVariable::names(bool printable) {
  static vector<string> names = Utils::tokenize("GENESIZE GENEDENSITY SAMPLESIZE INVERSE_MAC");
  static vector<string> labels = Utils::tokenize("gene size,gene density,sample size,inverse mac", ',');  
  return printable ? labels : names;  
}

InternalVariable::Type InternalVariable::type(const string& name, bool strip) {
  string var_name = Utils::uppercase(name);
  if (strip) {
    int pos = var_name.find("##");
    if (pos != string::npos) var_name = var_name.substr(pos+2);
  }
  
  vector<string>& internal = names();
  if (var_name == internal[iv_GeneSize] || var_name == "SIZE" || var_name == "NSNPS") return iv_GeneSize;
  if (var_name == internal[iv_GeneDensity] || var_name == "DENSITY" || var_name == "NPARAM") return iv_GeneDensity;
  if (var_name == internal[iv_SampleSize] || var_name == "N" || var_name == "SAMPSIZE" || var_name == "NSAMP") return iv_SampleSize;
  if (var_name == internal[iv_InverseMac] || var_name == "MAC") return iv_InverseMac;    
  return iv_UNKNOWN;  
} 

string InternalVariable::format(const string& var_name) {
  string func, out = name(type(var_name, true), true);
  if (Utils::starts_with(var_name, "#INTERNAL-LOG##")) func = "log";
  return !func.empty() ? func + "(" + out + ")" : out;
}

void CompoundGeneMatrix::CorrBlock::init() {
  corrs.set_size(size, size, true);
  for (double *curr = corrs.begin(), *end = corrs.end(); curr < end; curr += size+1) *curr = 1;
}

void CompoundGeneMatrix::CorrBlock::add_value(int row, int col, double value) {
  if (row >= size) throw SetException("inverting gene matrix", "row index is out of bounds");
  if (col >= row || col < 0) throw SetException("inverting gene matrix", "position is outside lower triangle");
  corrs(row, col) = (correction != 1) ? pow(value, correction) : value;
}

void CompoundGeneMatrix::CorrBlock::add_row(int row, BaseBuffer<double>& data) {
  if (row >= size) throw SetException("inverting gene matrix", "row index is out of bounds");
  if (data.size() > row) throw SetException("inverting gene matrix", "row size is too large");
  
  double* write = corrs[row - data.size()] + row;
  if (correction != 1) for (int i = 0; i < data.size(); write += size, i++) *write = pow(data[i], correction);
  else for (int i = 0; i < data.size(); write += size, i++) *write = data[i];
}

void CompoundGeneMatrix::CorrBlock::mirror() {
  if (!mirrored) {mirrored = true;
    for (int col = 1; col < corrs.ncol(); col++) {
      for (int row = 0; row < col; row++) corrs(row,col) = corrs(col,row);
    }    
  }
}

int CompoundGeneMatrix::CorrBlock::invert() {
  if (!inverted) {
    inverted = true; mirrored = true;
    return MathUtils::invert_matrix(corrs, prune); 
  } else return 0;
}

void CompoundGeneMatrix::init(GeneData& data) {
  if (!data.corrs.initialised()) throw SetException("inverting gene matrix", "no correlation data in input");
  if (!data.is_consecutive()) throw SetException("inverting gene matrix", "GeneData object is not consecutive");
  dimension = rank = data.data_used(false);
  
  BaseBuffer<BaseBuffer<double>*>& corrs = data.corrs.get_buffer();
  int from = 0; IndepBlock* indep = 0;
  for (int i = 1; i <= dimension; i++) {
    if (i == dimension || corrs[i]->is_empty()) {
      if (from == i-1) {
        if (!indep) {
          indep = new IndepBlock(from, 1);
          blocks.push_back(indep);
        } else indep->size++;
      } else {
        CorrBlock* curr = new CorrBlock(from, i - from, correction, prune);
        for (int j = from; j < i; j++) curr->add_row(j - from, *corrs[j]);
        blocks.push_back(curr); indep = 0;
      }
      from = i;
    }
  }
}

void CompoundGeneMatrix::init(GeneData& data, set<int>& filter) {
  if (!data.corrs.initialised()) throw SetException("inverting gene matrix", "no correlation data in input");
  if (!data.is_consecutive()) throw SetException("inverting gene matrix", "GeneData object is not consecutive");
  dimension = rank = filter.size();
  
  BaseBuffer<BaseBuffer<double>*>& corrs = data.corrs.get_buffer();
  vector<int> ids; vector<short> indep;
  for (set<int>::iterator it = filter.begin(); it != filter.end(); ++it) {
    bool curr = ids.empty() || (corrs[*it]->size() < (*it - ids.back()));
    ids.push_back(*it); indep.push_back(curr);
  }

  int from = 0; IndepBlock* ib = 0;
  for (int i = 1; i <= dimension; i++) {
    if (i == dimension || indep[i]) {
      if (from == i-1) {
        if (!ib) {
          ib = new IndepBlock(from, 1);
          blocks.push_back(ib);
        } else ib->size++;
      } else {
        CorrBlock* curr = new CorrBlock(from, i - from, correction, prune);
        for (int r = 1; r < i - from; r++) {int row_id = ids[from+r];
          BaseBuffer<double>& curr_corrs = *corrs[row_id];
          for (int c = 0; c < r; c++) {
            int index = ids[from+c] - (row_id - curr_corrs.size());
            if (index >= 0) curr->add_value(r, c, curr_corrs[index]);           
          }
        }
        blocks.push_back(curr); ib = 0;
      }
      from = i;
    }
  }
}


int CompoundGeneMatrix::invert_matrix(bool verbose) {
  if (!inverted) {rank = 0;
    for (int i = 0; i < blocks.size(); i++) {
      if (verbose) {_LOG.counter(true) << "\tprocessing block " << (i+1) << " of " << blocks.size(); _LOG.counter(false);}
      rank += blocks[i]->invert();
    }
    if (verbose) _LOG.wipe_counter();
    inverted = true;
  }
  
  return rank;
}  

void CompoundGeneMatrix::store_matrix(Buffer<double>& target) {
  target.set_size(dimension, dimension, true);
  for (int i = 0; i < blocks.size(); i++) blocks[i]->store(target); 
}

void CompoundGeneMatrix::product(set<int>* var, double* target) {
  set<int>::iterator curr = var->begin();
  for (int i = 0; i < blocks.size(); i++) blocks[i]->product(curr, var->end(), target);
}

double CompoundGeneMatrix::square(float* var) {double out = 0;
  for (int i = 0; i < blocks.size(); i++) out += blocks[i]->square(var);
  return out;
}

double CompoundGeneMatrix::square(double* var) {double out = 0;
  for (int i = 0; i < blocks.size(); i++) out += blocks[i]->square(var);
  return out;
}

double CompoundGeneMatrix::square(set<int>* var) {
  double out = 0; set<int>::iterator curr = var->begin();
  for (int i = 0; i < blocks.size(); i++) out += blocks[i]->square(curr, var->end());
  return out;
}



void CompoundGeneMatrix::CorrBlock::product(double*& target) {mirror();
  for (int i = 0; i < size; i++) *(target++) = MathUtils::sum(corrs[i], size);
}

void CompoundGeneMatrix::CorrBlock::product(float*& src, double*& target) {mirror();
  for (int i = 0; i < size; i++) target[i] = MathUtils::in_product(src, corrs[i], size);
  src += size; target += size;
}
void CompoundGeneMatrix::CorrBlock::product(double*& src, double*& target) {mirror();
  for (int i = 0; i < size; i++) target[i] = MathUtils::in_product(src, corrs[i], size);
  src += size; target += size;
}

void CompoundGeneMatrix::CorrBlock::product(set<int>::iterator& src, set<int>::iterator end, double*& target) {mirror();
  int to = from + size; 
  if (src != end && *src < to) {
    memcpy(target, corrs[*(src++) - from], sizeof(double)*size);
    while (src != end && *src < to) MathUtils::vec_sum(corrs[*(src++) - from], target, size); 
    target += size;
  } else {Utils::set_zero(target, size); target += size;}
}

void CompoundGeneMatrix::IndepBlock::product(double*& target) {
  for (double* end = target + size; target < end; target++) *target = 1;
}

void CompoundGeneMatrix::IndepBlock::product(float*& src, double*& target) {for (float* end = src + size; src < end; ++src, ++target) *target = *src;}

void CompoundGeneMatrix::IndepBlock::product(double*& src, double*& target) {
  memcpy(target, src, size*sizeof(double));
  src += size; target += size;
}

void CompoundGeneMatrix::IndepBlock::product(set<int>::iterator& src, set<int>::iterator end, double*& target) {
  int to = from + size; Utils::set_zero(target, size); 
  while (src != end && *src < to) {target[*src - from] = 1; ++src;}   
  target += size;
}

double CompoundGeneMatrix::CorrBlock::square(float*& var) {
  if (size > 1) {  
    Map<Matrix<float,Dynamic,1> > var_map(var, size, 1);
    Matrix<double,1,1> out = var_map.cast<double>().transpose() * Map<Matrix<double,Dynamic,Dynamic> >(corrs.begin(), size, size) * var_map.cast<double>();
    var += size; return out(0,0);
  } else {double out = var[0]*var[0]*corrs(0,0); var++; return out;}
}
double CompoundGeneMatrix::CorrBlock::square(double*& var) {double out = MathUtils::quad_form(corrs.begin(), var, size); var += size; return out;}
double CompoundGeneMatrix::CorrBlock::square(set<int>::iterator& var, set<int>::iterator end) {
  double diag = 0, off = 0; int to = from + size; set<int>::iterator orig = var;
  while (var != end && *var < to) {
    double* column = corrs[*var - from] - from; 
    for (set<int>::iterator curr = orig; curr != var; ++curr) off += column[*curr];
    diag += column[*(var++)];
  }
  return diag + 2*off;
}

double CompoundGeneMatrix::IndepBlock::square(float*& var) {double out = MathUtils::in_product(var, size); var += size; return out;}
double CompoundGeneMatrix::IndepBlock::square(double*& var) {double out = MathUtils::in_product(var, size); var += size; return out;}
double CompoundGeneMatrix::IndepBlock::square(set<int>::iterator& var, set<int>::iterator end) {
  double out = 0; int to = from + size;
  while (var != end && *var < to) {++out; ++var;}
  return out;
}


char DummyInfo::dummy = '\0';

char* SetInfo::get_name() {return index->name.get(iid);}
short SetInfo::get_status() {return index->status.get(iid);}
int SetInfo::get_column() {return index->storage_column.get(iid);}  
int SetInfo::get_ngenes() {return index->used_genes.get(iid);}  
set<int>* SetInfo::get_data() {int col = get_column(); return (col >= 0) ? index->data.get_data(col) : 0;}
double SetInfo::compute_sd(int tot_genes) {
  int no_genes = get_ngenes(); double mean = no_genes / double(tot_genes);
  return sqrt(no_genes * (1-mean) / (tot_genes - 1));
}

char* CovarInfo::get_name() {return index->name.get(iid);}
short CovarInfo::get_status() {return index->status.get(iid);}
int CovarInfo::get_column() {return index->storage_column.get(iid);}  
int CovarInfo::get_obs_genes() {return index->values_observed.initialised() ? index->values_observed.get(iid) : 0;}  
double CovarInfo::get_orig_mean() {return index->orig_mean.get(iid);}  
double CovarInfo::get_orig_sd() {return index->orig_sd.get(iid);}  
float* CovarInfo::get_data() {int col = get_column(); return (col >= 0) ? index->data.get_buffer()[col] : 0;}  

string InternalCovarInfo::format_name() {return InternalVariable::format(get_name());}

VariableInfo* SetSetPair::var(int i) {return (i == 0 ? set1 : (i == 1 ? set2 : 0));}
VariableInfo* SetCovPair::var(int i) {
  if (i == 0) return set;
  else if (i == 1) return cov;
  else return 0;
}


CollinearType BlockInverter::invert(Buffer<double>* X, Buffer<double>* A_inv, Buffer<double>* target) {
  int size_X = X->nrow(), size_A = A_inv->nrow(), size_D = size_X - size_A;
  if (X->ncol() != size_X) throw MathException("matrix inversion", "trying to invert a non-square matrix");
  if (A_inv->ncol() != size_A || size_A >= size_X) throw MathException("matrix inversion", "submatrix A is invalid");

  M.set_size(size_D, size_D); W.set_size(size_A, size_D); WM.set_size(size_A, size_D);
  SubMap B = X->matrix_sub(0, size_A, size_A, size_D), D = X->matrix_sub(size_A, size_D, size_A, size_D);

  W.matrix() = A_inv->matrix() * B;
  M.matrix() = D - B.transpose() * W.matrix();

  target->set_size(size_X, size_X);  

  CollinearType status = ct_None;    
  if (size_D > 1) {
    EigenInverter<double> eigen_M(M), eigen_D(*X, size_D, size_A); 
    double det_M = eigen_M.determinant(true); pair<double,double> det_D = eigen_D.determinant_ratio();

    if (det_D.second < ratio_thresh) status = ct_Variables;
    if (EigenInverter<double>::log_ratio(det_M, det_D.first) < ratio_thresh) status = ct_All;

    if (status == ct_None) {      
      eigen_M.set_cutoff(eigen_prune, true);
      eigen_M.inverse(*target, true, size_A);
  
      WM.matrix() = W.matrix() * target->matrix_sub(size_A, size_D, size_A, size_D);
      target->matrix_sub(0, size_A, 0, size_A) = A_inv->matrix() +  WM.matrix() * W.matrix().transpose();
      target->matrix_sub(0, size_A, size_A, size_D) = -1 * WM.matrix();
    }
  } else {
    if (M(0,0)/D(0,0) >= ratio_thresh) {
      double& M_inv = (*target)(size_A, size_A); M_inv = 1 / M(0,0);
          
      WM.matrix() = M_inv * W.matrix();
      target->matrix_sub(0, size_A, 0, size_A) = A_inv->matrix() + M_inv * W.matrix() * W.matrix().transpose();
      target->matrix_sub(0, size_A, size_A, size_D) = -1 * M_inv * W.matrix();
    } else status = ct_All;
  }
  if (status != ct_None) {
    if (allow_collinear) {
      EigenInverter<double> eigen(*X);
      eigen.set_cutoff(eigen_prune, true);
      eigen.inverse(*target);
    } else {
      if (status == ct_All) {
        if (size_D > 1) throw MathException("matrix inversion", "variables are collinear with conditioned covariates");
        else throw MathException("matrix inversion", "variable is collinear with conditioned covariates");
      } else throw MathException("matrix inversion", "variables are collinear with each other");
    }
  } else target->matrix_sub(size_A, size_D, 0, size_A) = target->matrix_sub(0, size_A, size_A, size_D).transpose(); 
  return status;
}

void GeneSetData::unlink_data(bool do_delete) {
  storage_column.assign_missing();
  if (data.get_owner()) data.get_owner()->unlink_variable(data);
  if (do_delete) data.clear_data(); 
}


int GeneVariableMap::add_variable(const string& prefix, const string& name, int block_index, int offset, VariableType type) {
  if (!index_built) {
    int iid = this->get_iid(prefix + " " + name); VariableInfo* info;
    if (type == vt_GeneSet) info = new SetInfo(variable_data.set_info[block_index], iid, offset); 
    else {
      if (block_index == variable_data.internal_block) {
        GeneCovarData<float>& internal = *variable_data.covar_info[block_index];
        if (!internal.using_variable("INTERNAL_TYPE")) throw DataException("processing variables", "internal variable type information is missing");
        short type = internal.typed_variable<NonNegativeNumber<short> >("INTERNAL_TYPE")->get(offset);
        info = new InternalCovarInfo(variable_data.covar_info[block_index], iid, offset, type);
      } else info = new ExternalCovarInfo(variable_data.covar_info[block_index], iid, offset);                                            
    }
    variable.set(iid, info);

    return iid;
  } else return -1;
}

int GeneVariableMap::set_interaction(int iid, int main1_id, int main2_id, int block_index, int offset, InteractionType type) {
  if (index_built && valid_id(iid)) {
    VariableInfo *info = 0, *main1 = &var_info(main1_id), *main2 = &var_info(main2_id);
    if (type == it_InteractionSS) {
      if (main1->is_set() && main2->is_set()) {
        info = new SetSetInteractionInfo(variable_data.set_info[block_index], iid, offset, main1->as_set(), main2->as_set());
      }
    } else if (type == it_InteractionSC) {
      if (main1->is_covar() && main2->is_set()) swap(main1, main2);
      if (main1->is_set() && main2->is_covar()) { 
        info = new SetCovInteractionInfo(variable_data.covar_info[block_index], iid, offset, main1->as_set(), main2->as_covar());
      }
    }
    if (info) {
      variable.set(iid, info);
      return iid;
    } else return 0;
  } else return -1;
}

set<int> GeneVariableMap::find_var(const string& name, VariableClass type) {set<int> out;
  if (type == vc_Internal) {
    GeneCovarData<float>& internal = *variable_data.covar_info[variable_data.internal_block]; InternalVariable::Type type_id = InternalVariable::type(name, true);
    BaseBuffer<short>& types = internal.typed_variable<NonNegativeNumber<short> >("INTERNAL_TYPE")->get_buffer();
    for (int v = 0; v < types.size(); v++) {if (types[v] == type_id) out.insert(internal.variable_id.get(v));}
  } else {
    if (type == vc_AnyExternal || type == vc_GeneSet) {
      for (int s = 0; s < variable_data.set_info.size(); s++) {
        int iid = (*variable_data.set_info[s])[name];
        if (iid >= 0) out.insert(variable_data.set_info[s]->variable_id.get(iid));
      }
    }
    if (type == vc_AnyExternal || type == vc_GeneCovar || type == vc_GeneCovarExternal) {
      for (int c = 0; c < variable_data.covar_info.size(); c++) {
        if (c != variable_data.internal_block || type == vc_GeneCovar) {
          int iid = (*variable_data.covar_info[c])[name];
          if (iid >= 0) out.insert(variable_data.covar_info[c]->variable_id.get(iid));
        }
      }
    }
  } 
  return out;
}

set<int> GeneVariableMap::var_list(VariableClass type, set<int>& filter) {set<int> out;
  for (int v = 0; v < data_total(); v++) {
    VariableInfo& info = var_info(v);
    if (info.is_valid() && info.has_type(type) && (filter.find(v) == filter.end())) out.insert(v);
  }
  return out;
}

pair<int,int> GeneVariableMap::no_sets() {
  int avail = 0, total = 0;
  for (int s = 0; s < variable_data.set_info.size(); s++) {
    total += variable_data.set_info[s]->data_total(); 
    avail += variable_data.set_info[s]->count_valid();
  }
  return pair<int,int>(avail,total);
}

pair<int,int> GeneVariableMap::no_covar(bool skip_internal) {
  int avail = 0, total = 0;
  for (int c = 0; c < variable_data.covar_info.size(); c++) {
    if (!skip_internal || c != variable_data.internal_block) {
      total += variable_data.covar_info[c]->data_total(); 
      avail += variable_data.covar_info[c]->count_valid();
    }
  }
  return pair<int,int>(avail,total);
}

pair<int,int> GeneVariableMap::no_vars(bool skip_internal) {
  pair<int,int> sets = no_sets(), covar = no_covar(skip_internal);
  return pair<int,int>(sets.first+covar.first, sets.second+covar.second);
}

/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#define INTERNAL_MIN_VALUE 0.001
#define MIN_COVAR_SD 1e-10

template<typename VAR> 
void SetVariables::set_internal(VAR& source, GeneCovarData<float>& var_info, SetStatsUtils::InternalVariable::Type var_id) {
  BlockMultiVariable<NumericVariable<float> >& target = var_info.data; 
  int main_index = var_id, log_index = main_index + SetStatsUtils::InternalVariable::iv_COUNT;
  bool absent = !source.initialised(), invert = var_id == SetStatsUtils::InternalVariable::iv_InverseMac; float value = target.get_missing().second;
  for (int gene_id = 0; gene_id < gene_info.data_total(); gene_id++) {
    if (!gene_info.active_id(gene_id) || source.is_missing(gene_id)) continue;
    if (!absent) {
      value = source.get(gene_id); 
      if (value < INTERNAL_MIN_VALUE) value = INTERNAL_MIN_VALUE;
      if (invert) value = 1/value;
    }
    target.set(gene_id, main_index, value);
    target.set(gene_id, log_index, log(value));
  }
  if (absent) {var_info.status.set(main_index, VariableIndex::NotAvailable); var_info.status.set(log_index, VariableIndex::NotAvailable);}
}

template<typename T>
void SetVariables::prep_covar(GeneCovarData<T>& covar_block, double trunc_low, double trunc_high, bool zstat) {
  Buffer<T>& data = covar_block.data.get_buffer();  
  covar_block.block_init(GeneCovarData<T>::ProcessBlock);
  int nrow = gene_info.data_used(false), ncol = covar_block.data_total();


  for (int v = 0; v < ncol; v++) {
    if (covar_block.is_valid(v)) {
      double sum = MathUtils::sum(data[v], nrow), mean = sum / nrow, sd = MathUtils::get_sd(data[v], nrow, mean);
      int trunc_count = 0; bool mono = true; 
      double low = zstat ? -trunc_low : mean - trunc_low*sd, high = mean + trunc_high*sd;
      if (nrow > 2) {
        T values[2] = {data(0,v), data(0,v)}; int counts[2] = {0,0}; 
        for (T *curr = data[v], *end = curr+nrow; curr < end; curr++) {
          if (*curr > high) {sum -= (*curr - high); *curr = high; trunc_count++;}
          else if (*curr < low) {sum += (low - *curr); *curr = low; trunc_count++;}
          if (mono) {
            if (*curr == values[0]) counts[0]++;
            else if (counts[1] == 0) {values[1] = *curr; counts[1]++;}
            else if (*curr == values[1]) counts[1]++;
            else mono = false;
            if (counts[0] > 1 && counts[1] > 1) mono = false;
          }          
        }
      }
      if (mono) {covar_block.status.set(v, VariableIndex::Monotonic); continue;}
      if (trunc_count > 0) {mean = sum / nrow; sd = MathUtils::get_sd(data[v], nrow, mean);}

      covar_block.orig_mean.set(v, mean); 
      covar_block.orig_sd.set(v, sd > MIN_COVAR_SD ? sd : 0); 
      covar_block.truncation_count.set(v, trunc_count);

      if (!zstat) {
        double inv_sd = sd > MIN_COVAR_SD ? 1 / sd : 1 / MIN_COVAR_SD;
        for (T *curr = data[v], *end = curr+nrow; curr < end; curr++) *curr = (*curr - mean) * inv_sd;
      }
    }
  }    
}



/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "utils_dataframe_variables.h"
#include "utils_dataframe_variables_genecorrelation.h"
  

void GeneCorrelationVariable::shrink(vector<int>& use) {
  if (!data.initialised() || !data.has_map()) return;
  BaseBuffer<BaseBuffer<double>*>& corrs = data.get_data(); BaseBuffer<short> used_rows(corrs.size(), true);
  for (int i = 0; i < use.size(); i++) used_rows[use[i]] = true;

  for (int i = 0; i < use.size(); i++) {int lid = use[i];
    BaseBuffer<double>* curr = corrs[lid]; 
    if (curr) {
      int write = 0, offset = lid - curr->size();
      for (int j = 0; j < curr->size(); j++) {int oid = offset + j;
        if (oid >= 0 && used_rows[oid]) {
          if (write < j) curr->get(write) = curr->get(j); 
          write++;
        }
      }
      curr->resize(write, true);
    }
  }

  JaggedMultiVariable<PosCorrelationVariable<double> >::shrink(use);
}

void GeneCorrelationVariable::reorder(vector<int>& use) {
  if (!data.initialised() || use.size() > data.get_data().size()) return;
  BaseBuffer<BaseBuffer<double>*>& old_corrs = data.get_data();

  BaseBuffer<BaseBuffer<double>*> new_corrs(use.size(), true); vector<int> index(old_corrs.size(), -1);
  for (int i = 0; i < use.size(); i++) index[use[i]] = i;

  for (int i = 0; i < index.size(); i++) {if (!old_corrs[i] || index[i] < 0) continue;
    BaseBuffer<double>& curr = *(old_corrs[i]); int offset = i - curr.size();
    for (int j = 0; j < curr.size(); j++) reorder_insert(new_corrs, index[i], index[offset + j], curr[j]);
  }

  if (!new_corrs[use.size()-1]) new_corrs[use.size()-1] = new BaseBuffer<double>();
  for (int i = use.size() - 2; i >= 0; i--) {
    if (!new_corrs[i] || new_corrs[i]->size() < (new_corrs[i+1]->size()-1)) reorder_pad(new_corrs[i], new_corrs[i+1]->size()-1);
  }  

  for (int i = 0; i < old_corrs.size(); i++) delete old_corrs[i];
  new_corrs.swap(old_corrs);
}

void GeneCorrelationVariable::reorder_insert(BaseBuffer<BaseBuffer<double>*>& corrs, int row, int col, double value) {
  if (col > row) swap(col, row); if (col < 0) return;
  BaseBuffer<double>*& buff = corrs[row]; int distance = row - col; 
  if (!buff || buff->size() < distance) {reorder_pad(buff, distance); buff->get(0) = value;}
  else buff->get(buff->size() - distance) = value;
}

void GeneCorrelationVariable::reorder_pad(BaseBuffer<double>*& buff, int required) {
  BaseBuffer<double>* new_buff = new BaseBuffer<double>(required, true); 
  if (buff) {memcpy(new_buff->data() + required - buff->size(), buff->data(), buff->size()*sizeof(double)); delete buff;}
  buff = new_buff;  
}  
 
void GeneCorrelationVariable::set_length(int iid, int length) {
  if (data.exists(iid)) {
    BaseBuffer<double>& buff = data.get_data(iid);
    if (buff.size() < length) {
      BaseBuffer<double> new_buff(length, true);
      memcpy(new_buff.data()+length-buff.size(), buff.data(), sizeof(double)*buff.size());
      buff.swap(new_buff);
    } else if (buff.size() > length) {int shift = buff.size() - length;
      memmove(buff.data(), buff.data()+shift, sizeof(double)*length);      
      buff.resize(length);
    }
  }
} 


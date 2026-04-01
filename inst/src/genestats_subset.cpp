/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genestats_subset.h"

#define MIN_VAR 1e-5
#define LINK_MARGIN 0.95

void SubsetStats::Subset::set(int samp_size, int step_init) {
  if (!step_init) step = max_step > 1 ? MathUtils::clamp(samp_size / min_size, 2, max_step) : 1; 
  else step = step_init;
  
  if (step > 1) {
    offset = ((samp_size - 1) % step) / 2;
    size = (samp_size - offset) / step;
    if ( (samp_size - offset) % step > 0 ) size++;
  } else {offset = 0; size = samp_size;}
}

void SubsetStats::LinkedSubset::set(SubsetStats::Subset& source, int samp_size) {
  set(samp_size); if (mode == Linked) delete proxy; proxy = &source; 
  
  if (parent == proxy->parent && step == proxy->step && offset == proxy->offset) mode = Cloned;
  else {
    mode = Indep;
    if (source.step > 1 && source.step < step) {
      int substep = max_step/source.step;
      while (substep > 1 && source.size / substep < min_size*LINK_MARGIN) substep--;
      if (substep > 1) {
        mode = Linked; link = &source;
        proxy = new Subset(settings); 
        proxy->min_size = min_size; proxy->max_step = max_step/source.step;
        proxy->set(source.size, substep);
      }
    }
  }
}

void SubsetStats::LinkedSubset::set_gene(LinkedGene* d) {
  if (mode == Indep) Subset::set_gene(d);
  else proxy->set_gene(d);
}    

void SubsetStats::Subset::set_data(Gene* snp_data) {
  if (!data) data = new LinkedGene(settings);
  data->set_subset(step, offset);
  parent = snp_data; loaded = false;
}
void SubsetStats::LinkedSubset::set_data(Gene* snp_data) {
  if (mode == Indep) Subset::set_data(snp_data);
  else if (mode == Linked) proxy->set_data(link->data);
}

void SubsetStats::Subset::load_data() {
  if (dynamic_cast<LoadedGene*>(parent)) data->set_parent(dynamic_cast<LoadedGene*>(parent));
  else if (dynamic_cast<LinkedGene*>(parent)) data->set_parent(dynamic_cast<LinkedGene*>(parent));  
  loaded = true;
}

Gene* SubsetStats::Subset::get_data() {
  if (step <= 1) return parent;
  if (!loaded) load_data();
  return data;
}

Gene* SubsetStats::LinkedSubset::get_data() {
  if (mode == Indep) return Subset::get_data(); 
  else return proxy->get_data();
}

Gene* SubsetStats::Subset::eject_data() {
  if (step > 1) {
    Gene* out = get_data(); data = 0; 
    return out;
  } else return 0;
}

Gene* SubsetStats::LinkedSubset::eject_data() {
  if (mode == Indep) return Subset::eject_data();
  else return proxy->eject_data();
}

void SubsetStats::clear() {
  DefaultStats::clear();
  
  wtw.reset();
  for (int i = 0; i < pc_sub.size(); i++) {delete pc_sub[i].data;}
  pc_sub.clear();
}  

void SubsetStats::init() {
  if (!check_size(data_size)) return;
  DefaultStats::init();
  if (!check_size(samp_size)) return;  

  sub_xtx.set(samp_size);
  sub_corrs.set(sub_xtx, samp_size);
}

bool SubsetStats::check_size(int size) {
  if ((auto_subset && size < auto_size) || (size < sub_xtx.min_size && size < sub_corrs.min_size)) {
    messages[NoSubset] = "";
    return false; 
  } else return true;
}

void SubsetStats::process_snps() {
  DefaultStats::process_snps();
  no_var = snp_data->nsnp(false);
  
  if (sub_corrs.is_empty()) sub_corrs.set_gene(dynamic_cast<LinkedGene*>(fetch_gene(no_var))); 
  sub_xtx.set_data(snp_data);
  sub_corrs.set_data(snp_data); ///may load from sub_xtx
}

LoadedGene* SubsetStats::set_gene(int size) {UNUSED(size);
  if (snp_data == 0) snp_data = make_gene();
  return snp_data;
}

GeneStats::GeneDataBlock* SubsetStats::eject_gene() {
  bool prev_ejected = sub_corrs.is_empty() && ejected_gene;
  Gene* gene = prev_ejected ? ejected_gene : sub_corrs.eject_data();
  ejected_gene = gene; 

  if (!gene) {_LOG.error("processing gene data") << "SubsetStats object has no data to eject" << endl; die();}

  GeneDataBlock* out = new GeneDataBlock(baseinput.gene_info.location.get(gene_id), gene, !prev_ejected);
  return out;  
}

GeneStats::BufferDataBlock* SubsetStats::eject_pcs(int pcid) {
  if (pcid >= 0) set_pcs(pcid);
  subset_pcs(); SubsetPC& pc = pc_sub[pc_id];

  pc.data->shrink(princomp->npcs);   
  BufferDataBlock* out = new BufferDataBlock(baseinput.gene_info.location.get(gene_id), pc.data);
  pc.data = 0; princomp->rewind();
  return out;
}  

void SubsetStats::compute_xtv(Buffer<float>& buffer, Buffer<float>& var, bool transpose, ProcessType::Type mode) {
  sub_xtx.get_data()->product_var(buffer, var, block_offset, block_size, transpose, mode);
}

void SubsetStats::compute_xtx(Buffer<float>& target, int offset, int total) {sub_xtx.get_data()->product(target, offset, total);}

void SubsetStats::compute_xtx_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode) {  
  sub_corrs.get_data()->product(buffer, block_index(block_i,0), block_index(block_i,1), block_index(block_j,0), block_index(block_j,1), mode);
}

void SubsetStats::compute_xtx_block_skip(Buffer<float>& buffer, int block_i, int step_i, int block_j, int step_j, ProductType::Type mode) {
  sub_corrs.get_data()->product_skip(buffer, block_index(block_i,0), step_i, block_index(block_i,1), block_index(block_j,0), step_j, block_index(block_j,1), mode);
}

void SubsetStats::compute_wtw(Buffer<float>& buffer, int block) {subset_pcs(); 
  MathUtils::matrix_prod(princomp->pc_block(*(pc_sub[pc_id].data), block).view(), buffer, MathUtils::TransposeFirst);
  MathUtils::matrix_scale(buffer, (samp_size-1)/double(sub_corrs.size-1));
  
  double thresh = MIN_VAR*(samp_size-1); int nvar = buffer.ncol();
  for (int i = 0; i < nvar; i++) {
    if (buffer(i,i) < thresh) {
      for (int j = 0; j < nvar; j++) buffer(i,j) = buffer(j,i) = 0;
      buffer(i,i) = samp_size-1;
    }
  }
}
 
void SubsetStats::compute_wtw_block(Buffer<float>& buffer, int block_i, int block_j, ProductType::Type mode) {subset_pcs();
  BufferWindow<float> pc_i = princomp->pc_block(*(pc_sub[pc_id].data), block_i), pc_j = princomp->pc_block(*(pc_sub[pc_id].data), block_j);
  MathUtils::matrix_prod(pc_i.view(), pc_j.view(), buffer, MathUtils::TransposeFirst);
  if (mode > ProductType::Centered) MathUtils::matrix_scale(buffer, 1/double(sub_corrs.size-1));
}

BufferWindow<float> SubsetStats::get_wtw() {
  compute_wtw(wtw, has_partition ? curr_block : -1);
  return wtw.window();
}

void SubsetStats::subset_pcs() {
  while (pc_sub.size() <= pc_id) pc_sub.push_back(SubsetPC());
  SubsetPC& pc = pc_sub[pc_id];
  
  if (!pc.data) pc.processed = 0;
  if (pc.processed > curr_block) return;
  for (int i = pc.processed; i <= curr_block; i++) {
    if (!princomp->computed(i)) {_LOG.error("processing gene data") << "no PCA data to process for SubsetStats object" << endl; die();}
  }  

  int ncol = has_partition ?  no_var : npcs(true);
  if (!pc.data) pc.data = fetch_buffer(sub_corrs.size, ncol);
  else if (pc.processed == 0) pc.data->set_size(sub_corrs.size, ncol);  

  if (has_partition) {
    for (int i = pc.processed; i <= curr_block; i++) {
      BufferWindow<float> pcs = princomp->get_pcs(i), target = princomp->pc_block(*pc.data, i);      
      target.view().matrix().noalias() = pcs.view().matrix_rowskip(sub_corrs.offset, sub_corrs.size, sub_corrs.step);  
      MathUtils::normalize(target.view());      
    }
  } else {
    BufferWindow<float> pcs = get_pcs();
    pc.data->matrix().noalias() = pcs.view().matrix_rowskip(sub_corrs.offset, sub_corrs.size, sub_corrs.step);  
    MathUtils::normalize(*(pc.data));
  }  
  pc.processed = curr_block+1;
}

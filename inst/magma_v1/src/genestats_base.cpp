/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "genestats_base.h"
#include "mathutils.h"
#include "statutils.h"
#include "exceptions.h"

int GeneStats::dump_count = 0;
int GeneStats::instance_count = 0;

void GeneStats::clear() {
  delete snp_data; snp_data = 0;
  for (int i = 0; i < gene_storage.size(); i++) delete gene_storage[i];
  gene_storage.clear();
  for (int i = 0; i < buffer_storage.size(); i++) delete buffer_storage[i];  
  buffer_storage.clear();

  dat_buffer.reset(); block_index.reset();
  initialized = false;
}

void GeneStats::insert_data(DataBlock* block) {
  insert_genes(block);
  insert_buffers(block);
  delete block;
}

void GeneStats::insert_genes(DataBlock* block) {Gene* curr_gene; 
  if (Utils::chr_type(chromosome) == 0 || chromosome == block->location.chromosome) {
    while ( (curr_gene = block->eject_gene()) != 0 ) {
      int curr_size = curr_gene->nmax();
      if (gene_storage.empty() || curr_size <= gene_storage.front()->nmax()) gene_storage.push_front(curr_gene);
      else if (curr_size >= gene_storage.back()->nmax()) gene_storage.push_back(curr_gene);
      else {
        for (deque<Gene*>::iterator curr = gene_storage.begin(); curr != gene_storage.end(); ++curr) {
          if (curr_size <= (*curr)->nmax()) {gene_storage.insert(curr, curr_gene); break;}
        }
      }
    }
  }
}

void GeneStats::insert_buffers(DataBlock* block) {
  Buffer<float>* curr_buff; 
  while ( (curr_buff = block->eject_buffer()) != 0 ) insert_buffer(curr_buff);
}

void GeneStats::insert_buffer(Buffer<float>* curr_buff) {
  if (!curr_buff) return;
  int curr_size = curr_buff->maxsize();
  if (buffer_storage.empty() || curr_size <= buffer_storage.front()->maxsize()) buffer_storage.push_front(curr_buff);
  else if (curr_size >= buffer_storage.back()->maxsize()) buffer_storage.push_back(curr_buff);
  else {
    for (deque<Buffer<float>*>::iterator curr = buffer_storage.begin(); curr != buffer_storage.end(); ++curr) {
      if (curr_size <= (*curr)->maxsize()) {buffer_storage.insert(curr, curr_buff); break;}
    }
  }
}

GeneStats::GeneDataBlock* GeneStats::eject_gene() {
  bool prev_ejected = !snp_data && ejected_gene;
  Gene* gene = prev_ejected ? ejected_gene : snp_data;
  ejected_gene = gene; snp_data = 0;

  GeneDataBlock* out = new GeneDataBlock(baseinput.gene_info.location.get(gene_id), gene, !prev_ejected);
  return out;  
}

void GeneStats::set_chromosome(int chr) {chromosome = chr;
  if (chr == Utils::chrY_code) selection = BaseInput::sm_MaleOnly;
  else if (chr == Utils::chrW_code) selection = BaseInput::sm_FemaleOnly;
  else if (Utils::chr_type(chr) > 0) selection = BaseInput::sm_KnownGender;
  else selection = BaseInput::sm_AllKnown;
}

Gene* GeneStats::fetch_gene(int size) {Gene* out = 0;
  if (!gene_storage.empty()) {
    if (size > gene_storage.front()->nmax() && size <= gene_storage.back()->nmax()) {
      for (int i = gene_storage.size()-1; i >= 0; i--) {
        if (size < gene_storage[i]->nmax()) continue;
        if (size > gene_storage[i]->nmax()) i++;
        out = gene_storage[i]; gene_storage.erase(gene_storage.begin()+i);
        break;
      }
    }
    if (!out) {out = gene_storage.front(); gene_storage.pop_front();}
  }
  if (out) out->set_N(samp_size);
  return out;
}

LoadedGene* GeneStats::set_gene(int size) {
  if (snp_data == 0) snp_data = dynamic_cast<LoadedGene*>(fetch_gene(size));
  if (snp_data == 0) snp_data = make_gene();

  return snp_data;
}

Buffer<float>* GeneStats::fetch_buffer(int nrow, int ncol, bool zero) {Buffer<float>* buff = 0; 
  if (!buffer_storage.empty()) {
    int required = nrow*ncol;
    if (required > buffer_storage.front()->maxsize() && required <= buffer_storage.back()->maxsize()) { 
      for (int i = buffer_storage.size()-1; i >= 0; i--) {
        if (required < buffer_storage[i]->maxsize()) continue;
        if (required > buffer_storage[i]->maxsize()) i++;
        buff = buffer_storage[i]; buffer_storage.erase(buffer_storage.begin()+i);
        break;
      }
    }           
    if (!buff) {buff = buffer_storage.front(); buffer_storage.pop_front();}
  }
  if (!buff) buff = new Buffer<float>;  
  buff->set_size(nrow, ncol, zero);
  return buff;
}

void GeneStats::set_buffer(Buffer<float>& target, int nrow, int ncol, bool zero) {
  if (!buffer_storage.empty() && target.maxsize() < nrow*ncol) {
    bool was_empty = target.empty();
    Buffer<float>* buff = fetch_buffer(nrow, ncol, zero);
    buff->swap(target);
    
    if (was_empty) delete buff;
    else insert_buffer(buff);
  } else target.set_size(nrow, ncol, zero);
}                                                                                            

void GeneStats::load_snps(int gid) {
  if (!initialized) init(); gene_id = gid;

  BaseBuffer<int>& snp_buff = baseinput.gene_info.snps.get_buffer(gene_id); 
  set<long> snps(snp_buff.data(), snp_buff.data()+snp_buff.size()); 
  snp_data = baseinput.load_snps(gene_id, snps, set_gene(snps.size()), selection);
  process_snps();

  load_snp_info();
  rewind();
}

void GeneStats::load_snp_info() {
  no_var = snp_data->nsnp(false); tot_snps = snp_data->nsnp(true); tot_rare = snp_data->nrare();
  gene_mac = MathUtils::get_mean(snp_data->get_mac().view().begin(), no_var);
}

void GeneStats::partition(vector<int>& blocks){
  int nblock = blocks.size(), offset = 0;
  block_index.set_size(nblock, 2);
  for (int i = 0; i < nblock; i++) {block_index(i,0) = offset; block_index(i,1) = blocks[i]; offset += blocks[i];}

  has_partition = true; curr_block = -1;
}

void GeneStats::set_block(int index) {
  if (index < 0 || (index == 0 && !has_partition)) {block_offset = 0; block_size = no_var; curr_block = 0;} 
  else if (!has_partition || index > block_index.nrow()) throw GeneException("subdividing gene", "trying to access unknown gene block");
  else {
    block_offset = block_index(index,0); 
    block_size = block_index(index,1);
    curr_block = index;
  }
}

BufferWindow<float> GeneStats::get_var() {return snp_data->get_variance_block(block_offset, block_size);}

int GeneStats::block_stepsize(int no_snps) {
  if (no_snps < 2*block_min_snps) return 1;
  return min(no_snps/block_min_snps, block_max_step);
} 

bool PermutationGenerator::prep_buffer(BufferMode mode, int nperm) {
  if (!buffers[mode][Main]) {_LOG.error("processing gene data") << "gene data objects are incompatible" << endl; die();}
  Buffer<float>& buff = *(buffers[mode][Main]);

  bool resize = (buff.nrow() != samp_size || buff.ncol() < nperm);
  if (resize) buff.set_size(samp_size, nperm);
  return resize;
}

/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "baseinput.h"
#include "input.h"

using namespace std;

void BaseInput::init() {
  _LOG.set_skip(true);
  if (settings["snp_synonyms"]) load_synonyms(settings.gets("snp_synonyms"), settings["snp_synonyms_auto"]);
  if (settings["snp_filter"]) filter_snps(settings.gets("snp_filter_file"), settings.gets("snp_filter"));
  if (has_pval) load_pvalfile();
  if (has_weights) load_burden_weights();  

  if (pheno_type == ph_FromFile) load_altpheno(settings.gets("pheno_file"));
  else if (pheno_type == ph_Dummy) set_dummy_pheno();
  else if (pheno_type == ph_Constant || pheno_type == ph_Trio) set_const_pheno();    

  load_covar();
  if (settings["indiv_filter"]) filter_indiv(settings.gets("indiv_filter_file"), settings.gets("indiv_filter"));
  filter_missing(has_pval);
  if (pedigree == pd_Trio) filter_trios();

  Utils::mem_purge(indiv_map);
  _LOG.set_skip(false);
  
  load_genedef(settings.gets("gene_annot_file"));
  if (settings["batch"]) batch_split();
}
 
void BaseInput::load_genedef(const string& filename) { 
  _LOG << "Loading gene annotation from file " << filename << "... " << endl;

  TextInput fin(filename, 2);
  fin.set_error("reading gene annotation file"); fin.skip_empty();

  gene_info.init(20000);
  GeneLocationVariable& loc_input = gene_info.location; loc_input.init();
  JaggedMultiVariable<NumericID<int> >& snp_input = gene_info.snps; snp_input.init();

  human_chr = !fin.has_param("nonhuman");
  loc_input.human() = human_chr; GeneLocation location; 

  int snp_id, no_defs = 0, chr_filter = settings.geti("batch_chr", 0); 
  set<long> snps; vector<short> chr_warning(7,false), used(snp_info.data_total(), false);

  while (fin.process_line()) {
    location = loc_input.convert(fin[1]);
    if (loc_input.get_status() == InputVariable::InvalidFormat) fin.line_error("improperly formatted location string", "format should be CHR:FROM:TO");
    if (!loc_input.is_valid()) {_LOG << "\tWARNING: on line " << fin.line_no << ", " << loc_input.get_msg() << "; skipping gene (ID = " << fin[0] << ")" << endl; continue;}
    if (location.chromosome < 1 || (chr_filter && location.chromosome != chr_filter)) continue;
    
    if (location.chromosome > Utils::chrX_code) {int chr_index = location.chromosome - Utils::chr_reserved;
      if (!chr_warning[chr_index]) _LOG << "\tWARNING: analysis of chromosome " << Utils::chr_string(location.chromosome,0) << " genes currently not supported; all such genes will be skipped" << endl;
      chr_warning[chr_index] = true; continue;
    }
     
    snps.clear();
    if (!fin.read_value()) fin.line_error("too few values");
    do {
      if ((snp_id = snp_info[fin.value]) < 0) continue;
      snps.insert(snp_id);
    } while(fin.read_value());

    no_defs++;
    if (snps.size() == 0) continue;
  
    int gid = gene_info.get_iid(fin[0]);
    loc_input.set(gid, location);
    snp_input.insert(gid, snps);
    for (set<long>::iterator it = snps.begin(); it != snps.end(); ++it) used[*it] = true; 
  }
  gene_info.build_index(); 

  if (!gene_info.get_duplicates().empty()) {vector<char*>& duplicates = gene_info.get_duplicates();
    _SLOG.set_block("Gene annotation file") << "# Following gene IDs had duplications in file " << filename << endl;
    for (int i = 0; i < duplicates.size(); i++) _SLOG << duplicates[i] << endl;
    fin.error(string("file contained ") + Utils::plural(duplicates.size(), "duplicate gene ID") + "; writing list of IDs to supplementary log file");    
  }

  if (gene_info.data_used(false) == 0) fin.error("found no genes containing SNPs in genotype data");
  _LOG << "\t" << no_defs << " gene definitions read from file" << endl;
  _LOG << "\tfound " << gene_info.data_used(false) << " genes containing valid SNPs in genotype data" << endl << endl;

  snp_info.filter_ids(used, true); used.clear();
  snp_info.clear_index();
}

void BaseInput::batch_split() {
  int curr_batch = settings.geti("batch_index"), no_batches = settings.geti("batch"); 

  vector<short> keep(gene_info.data_total(), 0);
  int used_batches = (settings["genes_only"]) ? batch_genes(keep, curr_batch-1, no_batches) : batch_all(keep, curr_batch-1, no_batches);

  if (used_batches < no_batches) {
    if (curr_batch <= used_batches) {
      _LOG << "WARNING: maximum number of batches for input data is " << used_batches << " (requested = " << no_batches << "); switching mode to --batch " << curr_batch << " " << used_batches << endl;
      settings.alter_value("batch", used_batches);
      _LOG.change_logfile(settings.gets("out_prefix") + ".log"); _SLOG.change_logfile(settings.gets("out_prefix") + ".log.suppl");      
    } else {_LOG.error("configuration") << "maximum number of batches for input data is " << used_batches << " (requested = " << no_batches << "), batch " << curr_batch << " is empty" << endl; die();} 
  }

  gene_info.filter_ids(keep, true);
  gene_info.shrink_storage();
}

int BaseInput::batch_genes(vector<short>& keep, int curr_batch, int no_batches) {int no_genes = gene_info.data_total();
  if (no_batches <= no_genes) {
    int base_size = no_genes / no_batches, remainder = no_genes % no_batches;
    int begin = curr_batch*base_size + min(curr_batch, remainder), end = begin + base_size + (curr_batch < remainder ? 1 : 0); 
    for (int i = begin; i < end; i++) keep[i] = 1;
  } else if (curr_batch < no_genes) keep[curr_batch] = 1;
  return no_genes;
}


int BaseInput::batch_all(vector<short>& keep, int curr_batch, int no_batches) {
  vector<pair<int,int> > blocks; deque<GeneLocation> glocs;
  long corr_range = settings.geti("corr_range");
  for (int i = 0; i < gene_info.data_total(); i++) {
    GeneLocation loc = gene_info.location.get(i);      
    while (!glocs.empty() && !loc.in_range(glocs.front(), corr_range)) glocs.pop_front();
    if (glocs.empty()) blocks.push_back(pair<int,int>(i,i)); 
    blocks.back().second++;

    glocs.push_back(loc);
  }

  if (no_batches > blocks.size()) {
    no_batches = blocks.size();
    if (curr_batch >= no_batches) return no_batches;
  }

  vector<int> use; 
  if (no_batches < blocks.size()) {
    vector<pair<int,float> > weights; 
    for (int i = 0; i < blocks.size(); i++) weights.push_back(pair<int,float>(i, blocks[i].second-blocks[i].first) );
    IndexSorter<int,float>(weights, false);

    vector<float> batch_weight(no_batches, 0); 
    for (int i = 0; i < blocks.size(); i++) {int index = 0;
      if (i >= no_batches) {for (int j = 1; j < no_batches; j++) {if (batch_weight[j] < batch_weight[index]) index = j;}}
      else index = i;
      batch_weight[index] += weights[i].second;
      if (index == curr_batch) use.push_back(weights[i].first);
    }
  } else use.push_back(curr_batch);

  for (int i = 0; i < use.size(); i++) {for (int j = blocks[use[i]].first; j < blocks[use[i]].second; j++) keep[j] = 1;}
  return blocks.size();
}


void BaseInput::set_buffers(LoadedGene* gene, int nsnps, float*& data_buff, double*& pval_buff, int*& n_buff, long*& id_buff, float*& weight_buff) {
  data_buff = gene->set_buffer(nsnps).begin(); 
  pval_buff = has_pval ? gene->get_buffer_pval() : 0;
  n_buff = has_pval && has_snpN ? gene->get_buffer_N() : 0;
  id_buff = gene->get_buffer_id(false);
  weight_buff = has_weights ? gene->get_buffer_weight() : 0;
}

LoadedGene* BaseInput::load_snps(int gene_id, const set<long>& snps, LoadedGene* gene, SubsetMode filter_id) {
  float* data_buff; double* pval_buff; int* n_buff; long* id_buff = 0; float* weight_buff;
  string* gene_name = new string(gene_info.name.get(gene_id)); 
  gene->set_id(gene_id, *gene_name);

  set_buffers(gene, snps.size(), data_buff, pval_buff, n_buff, id_buff, weight_buff);
  load_snps_core(snps, filter_id, data_buff, pval_buff, n_buff, id_buff, weight_buff);

  delete gene_name;
  return gene;
}

void BaseInput::load_snps_core(const set<long>& snps, SubsetMode filter_id, float*& data_buff, double*& pval_buff, int*& n_buff, long*& id_buff, float*& weight_buff) {
  unsigned long long written = load_snpdata(snps, data_buff, filter_id); data_buff += written;
  if (pval_buff) {
    for (set<long>::const_iterator iter = snps.begin(); iter != snps.end(); ++iter) {  
      *(pval_buff++) = snp_info.pval.get(*iter);
      if (n_buff) *(n_buff++) = snp_info.nsamp.get(*iter);
    }
  }
  if (id_buff) {
    for (set<long>::const_iterator iter = snps.begin(); iter != snps.end(); ++iter) *(id_buff++) = *iter;
  }
  if (weight_buff) {
    for (set<long>::const_iterator iter = snps.begin(); iter != snps.end(); ++iter) {*(weight_buff++) = snp_info.weights.get(*iter);}
  }
} 

short* BaseInput::get_filter(SubsetMode selection) {short index;
  switch (selection) {
    case sm_KnownGender:  index = 2; break;
    case sm_MaleOnly:     index = 3; break;
    case sm_FemaleOnly:   index = 4; break;    
    case sm_All:          index = 5; break;
    case sm_AllKnown: 
    default:              index = 0;
  }
  return indiv_missing[index];   
}

int BaseInput::check_duplicates(const string& filename, string header) {vector<char*>& duplicates = snp_info.get_duplicates();
  if (!duplicates.empty()) {
    _SLOG.set_block(header) << "# Following variable IDs had duplications in file " << filename << endl;
    for (int i = 0; i < duplicates.size(); i++) _SLOG << duplicates[i] << endl;
    _LOG << "\tWARNING: file contains " << Utils::plural(duplicates.size(), "duplicate SNP ID") << "; writing list of IDs to supplementary log file" << endl;
  }
  return duplicates.size();
}
  

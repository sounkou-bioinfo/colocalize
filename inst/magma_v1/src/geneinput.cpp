/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "geneinput.h"
#include "output.h"
#include "exceptions.h"
#include <dirent.h>
#include "utils.h"
#include "statutils.h"
#include "genestats.h"
#include "geneoutput.h"

#include <boost/math/distributions/normal.hpp>

void GeneInputFile::InputUniColumn::process_value(TextInput& fin, int gid) {
  variable->add(gid, fin[index]); 
  if ((allow_missing && !variable->is_valid_or_missing()) || (!allow_missing && !variable->is_valid())) fin.line_error(name + " " + variable->get_msg(), index+1);      
}

void GeneInputFile::add_column(InputVariable& var, const string& name, int index, bool allow_missing) {
  columns.push_back(new InputUniColumn(var, name, index, allow_missing));
}

void GeneInputOutFile::map_column(InputVariable& var, const string& name, bool allow_missing) {
  if (has_column(name)) add_column(var, name, column_index[name], allow_missing);
}

void GeneInputRawFile::load(double rescale, set<string> covar_list, bool list_exclude) {
  _LOG << "Reading file " << filename << "... " << endl;

  TextInput fin(filename);
  fin.set_error("reading gene results file");
  int ncol = fin.detect_size(), version_id = 0; fin.get_param("VERSION", version_id);

  gene_info.init(20000); gene_info.block_init(GeneData::CoreBlock);
  gene_info.location.human() = !fin.has_param("NONHUMAN");

  vector<double> corrs; corrs.reserve(1000); gene_info.corrs.init();
  int col_used = 4;
  add_column(gene_info.nsnps, "NSNPS", col_used++, (covar_list.find("NSNPS") != covar_list.end()) == list_exclude);
  add_column(gene_info.nparam, "NPARAM", col_used++, (covar_list.find("NPARAM") != covar_list.end()) == list_exclude);

  vector<string> covar = fin.get_param_values("COVAR");
  vector<string> zstat = fin.get_param_values("ZSTAT");  
  if (covar.empty() && zstat.empty() && version_id < 107) {int remainder = ncol - col_used - 1;
    if (remainder >= 1) covar.push_back("NSAMP");
    if (remainder >= 2) covar.push_back("MAC");
  }
  if (zstat.empty()) zstat.push_back("ZSTAT");
  for (int i = 0; i < covar.size(); i++, col_used++) {
    if (covar[i] == "NSAMP") add_column(gene_info.nsamp, covar[i], col_used, (covar_list.find("NSAMP") != covar_list.end()) == list_exclude);
    else if (covar[i] == "MAC") add_column(gene_info.mac, covar[i], col_used, (covar_list.find("MAC") != covar_list.end()) == list_exclude);
  }
  for (int i = 0; i < zstat.size(); i++, col_used++) { //for now only using first value
    if (i == 0) add_column(gene_info.zstat, zstat[i], col_used, true);  
  }
  if (ncol != col_used) fin.error("incorrect number of columns in input file; the file may be from a newer version of MAGMA or may have been altered");  

  while (fin.process_line()) {
    int gid = gene_info.get_iid(fin[0]);
    gene_info.location.add(gid, fin[1], fin[2], fin[3]); if (!gene_info.location.is_valid()) fin.line_error(gene_info.location.get_msg());
    for (int i = 0; i < columns.size(); i++) columns[i]->process_value(fin, gid);
        
    corrs.clear();
    while (fin.read_value()) {
      corrs.push_back(gene_info.corrs.convert(fin.value));
      if (!gene_info.corrs.is_valid(true)) fin.line_error("gene correlation " + gene_info.corrs.get_msg(true), col_used+corrs.size());      
    }
    if (rescale < 1) { for (int i = 0; i < corrs.size(); i++) {corrs[i] *= rescale;} }
    gene_info.corrs.insert(gid, corrs);
    
    int prev_size = (gid > 0) ? gene_info.corrs.get_buffer(gid-1).size() : 0;
    if (corrs.size() > prev_size + 1) fin.line_error("too many correlations (found " + Utils::num_string(corrs.size()) + ", maximum is " + Utils::num_string(prev_size+1)+ ")");
    if (gid > 0 && corrs.size() > 0 && !(gene_info.location.get(gid-1) <= gene_info.location.get(gid))) fin.error("genes are not ordered by genomic location");
  }  
  gene_info.build_index();
  gene_info.zstat.filter(); 

  if (gene_info.data_used(false) == 0) fin.error("found no valid genes in file");
  if (!gene_info.get_duplicates().empty()) fin.error("file contained duplicate gene IDs");  
  _LOG << "\t" << gene_info.data_used(false) << " genes read from file" << endl;
}

void GeneInputRawFile::overwrite_stat(const string& filename) {
  _LOG << "Overwriting Z-scores with values from file " << filename << "..." << endl;
  
  GeneInputOutFile file(filename, true); GeneData& stats = file.get_data();
  if (!stats.pval.initialised()) stats.upgrade_permp();
  if (!stats.zstat.initialised() && ! stats.create_zstat()) {_LOG.error("reading overwrite file") << "no Z-scores or p-values found in file" << endl; die();}

  gene_info.zstat.assign_missing();
  bool has_nsnps = stats.nsnps.initialised(), has_nparam = stats.nparam.initialised(), has_nsamp = stats.nsamp.initialised();
  if (has_nsnps) {
    if (!gene_info.nsnps.initialised()) gene_info.nsnps.init();
    gene_info.nsnps.assign_missing();
  }
  if (has_nparam) {
    if (!gene_info.nparam.initialised()) gene_info.nparam.init();
    gene_info.nparam.assign_missing();
  }
  if (has_nsamp) {
    if (!gene_info.nsamp.initialised()) gene_info.nsamp.init();
    gene_info.nsamp.assign_missing();  
  }
  for (int i = 0; i < stats.data_total(); i++) {
    if (stats.active_id(i)) {
      int iid = gene_info[stats.name.get(i)];
      if (iid >= 0) {
        gene_info.zstat.set(iid, stats.zstat.get(i), true);
        if (has_nsnps) gene_info.nsnps.set(iid, stats.nsnps.get(i), true);        
        if (has_nparam) gene_info.nparam.set(iid, stats.nparam.get(i), true);        
        if (has_nsamp) gene_info.nsamp.set(iid, stats.nsamp.get(i), true);                        
      }
    }
  }
  gene_info.zstat.filter();
  gene_info.shrink_storage();
  
  _LOG << "\toverwrote values for " << Utils::plural(gene_info.data_used(false), "gene") << "; all other genes have been discarded" << endl;
}
  
void GeneInputOutFile::load() {
  if (!internal_input) _LOG << "Reading file " << filename << "... " << endl;

  TextInput fin(filename);
  fin.set_error("reading gene results file");
  fin.read_header(); if (!fin.has_header) fin.error("file has no header");
  
  deque<string>& col_names = fin.get_names(); for (int i = 0; i < col_names.size(); i++) column_index[Utils::uppercase(col_names[i])] = i;
  if (column_index.find("NPARAM") == column_index.end() && column_index.find("EFF_SIZE") != column_index.end()) column_index["NPARAM"] = column_index["EFF_SIZE"];

  const char* req[] = {"GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N"}; vector<string> req_columns(req, req + sizeof(req)/sizeof(req[0]));
  for (int i = 0; i < (!internal_input ? req_columns.size() : 1); i++) {if (column_index.find(req_columns[i]) == column_index.end()) fin.error("required column '" + req_columns[i] + "' not found");}   
  int id_col = column_index["GENE"], chr_col = -1, start_col = -1, stop_col = -1; 

  gene_info.init(20000); 
  if (!internal_input) {
    gene_info.location.init(); gene_info.location.human() = !fin.has_param("NONHUMAN");
    chr_col = column_index["CHR"]; start_col = column_index["START"]; stop_col = column_index["STOP"];
  }
  map_column(gene_info.nsnps, "NSNPS");
  map_column(gene_info.nparam, "NPARAM");
  map_column(gene_info.nsamp, "N");
  map_column(gene_info.ndata, "DATASETS");  
  map_column(gene_info.zstat, "ZSTAT");    
  map_column(gene_info.nrare, "NRARE"); 
  multi_column(gene_info.rsq, "RSQ", "RSQ_ADJ");     

  vector<string> pcol; string permp = "";
  for (map<string,int>::iterator it = column_index.begin(); it != column_index.end(); ++it) {string name = it->first;
    if (name == "P" || name.substr(0,2) == "P_") pcol.push_back(name);
    else if (name.substr(0,5) == "PERMP") permp = name;
  }
  process_pval(pcol);

  int nsnps_part_col = -1, nsamp_part_col = -1;
  if (!pcol.empty()) {
    if (!gene_info.get_model_info().is_consistent()) {_LOG.error("reading gene results file") << "p-value columns are inconsistent" << endl; die();}
    if (pcol.size() > 1) {int npart = gene_info.get_model_info().state[ModelInfoBlock::PartitionCount] - 1;
      if (has_column("NSNPS_PART")) {nsnps_part_col = column_index["NSNPS_PART"]; gene_info.nsnps_part.init(npart);}
      if (has_column("NSAMP_PART")) {nsamp_part_col = column_index["NSAMP_PART"]; gene_info.nsamp_part.init(npart);}      
    }
  }
 
  bool has_perm = !permp.empty(); int permp_col = -1, nperm_col = -1;
  if (has_perm) {
     if (!has_column("NPERM")) {_LOG.error("reading gene results file") << "file contains permutation p-values but no NPERM column" << endl; die();}
     permp_col = column_index[permp]; gene_info.perm_pval.init(); 
     nperm_col = column_index["NPERM"]; gene_info.perm_count.init(); 
  }

  while (fin.process_line()) {
    int gid = gene_info.get_iid(fin[id_col]);
    if (chr_col >= 0) gene_info.location.add(gid, fin[chr_col], fin[start_col], fin[stop_col]); if (!gene_info.location.is_valid()) fin.line_error(gene_info.location.get_msg());  
    if (nsnps_part_col >= 0) {
      gene_info.nsnps_part.add(gid, fin[nsnps_part_col], ','); 
      if (!gene_info.nsnps_part.is_valid(false)) fin.line_error("for NSNPS_PART, " + gene_info.nsnps_part.get_msg(false),nsnps_part_col+1, gene_info.nsnps_part.get_msg_sub()); 
    }
    if (nsamp_part_col >= 0) {
      gene_info.nsamp_part.add(gid, fin[nsamp_part_col], ','); 
      if (!gene_info.nsamp_part.is_valid(false)) fin.line_error("for NSAMP_PART, " + gene_info.nsamp_part.get_msg(false),nsamp_part_col+1, gene_info.nsamp_part.get_msg_sub()); 
    }
    for (int i = 0; i < columns.size(); i++) columns[i]->process_value(fin, gid);
  
    if (has_perm) {  
      gene_info.perm_count.add(gid, fin[nperm_col]); if (!gene_info.perm_count.is_valid()) fin.line_error(gene_info.perm_count.get_msg());   
      gene_info.perm_pval.add(gid, fin[permp_col], gene_info.perm_count.get(gid)); if (!gene_info.perm_pval.is_valid()) fin.line_error(gene_info.perm_pval.get_msg());         
    }
  }
  gene_info.build_index();

  if (gene_info.data_used(false) == 0) fin.error("found no genes in file");
  if (!gene_info.get_duplicates().empty()) fin.error("file contained duplicate gene IDs");
  _LOG << "\t" << gene_info.data_used(false) << " genes read from file" << endl;
}  

void GeneInputOutFile::process_pval(vector<string>& pcol) {
  if (!pcol.empty()) {
    multi_column(gene_info.pval, pcol, true);
    if (pcol.size() > 1) {
      ModelInfoBlock info; vector<vector<string> > components; bool error = false;   
      for (int i = 0; i < pcol.size(); i++) {int col = 0; 
        vector<string> tokens = Utils::tokenize(Utils::uppercase(pcol[i]), '_');
        for (int j = 1; j < tokens.size(); j++) {
          if (tokens[j] == "") continue;
          if (i == 0) components.push_back(vector<string>());
          if (col > 0 && components[col-1].back() == "SNPWISE") components[col-1].back() += "_" + tokens[j];
          else {
            if (col >= components.size()) {error = true; break;}
            components[col++].push_back(tokens[j]);
          }
        }
        if (col < components.size() || col > 3) error = true;
        if (error) break;      
      } 
      
      vector<int> index(pp_Unknown, -1);
      if (!error) {
        vector<PvaluePart> pp_type(components.size(), pp_Unknown);
        for (int i = 0; i < components.size(); i++) {PvaluePart pp = pp_Unknown;
          if (check_consistency(components[i], pp_Partition)) pp = pp_Partition;
          else if (check_consistency(components[i], pp_Model)) pp = pp_Model;
          else if (check_consistency(components[i], pp_Analysis)) pp = pp_Analysis;
          if (pp == pp_Unknown || index[pp] >= 0) error = true;
          else index[pp] = i;
        }
      }
      if (error) {_LOG.error("reading gene results file") << "unable to convert p-value column names" << endl; die();} 
      int no_pval = components[0].size();

      vector<ModelState::Model> model(no_pval, ModelState::UnknownModel);
      if (index[pp_Model] >= 0) {vector<string>& names = components[index[pp_Model]];
        for (int i = 0; i < names.size(); i++) model[i] = ModelState::model_type(names[i]);
      }

      vector<ModelState::Analysis> analysis(no_pval, ModelState::MainAnalysis);
      if (index[pp_Analysis] >= 0) {vector<string>& names = components[index[pp_Analysis]];
        for (int i = 0; i < names.size(); i++) analysis[i] = ModelState::analysis_type(names[i]);
      }

      vector<int> partition(no_pval, 0); vector<string> partition_names; partition_names.push_back("#UNDEFINED#");
      if (index[pp_Partition] >= 0) {vector<string>& names = components[index[pp_Partition]];
        map<string,int> part_index; 
        for (int i = 0; i < names.size(); i++) {
          if (names[i] == "ALL-PART") continue;
          string part_name = names[i].substr(5);
          if (part_index.find(part_name) == part_index.end()) {part_index[part_name] = partition_names.size(); partition_names.push_back(part_name);}
          partition[i] = part_index[part_name];          
        }
      }

      for (int i = 0; i < no_pval; i++) info.info.push_back(new ModelOutputInfo(model[i], analysis[i], i, partition[i]));  
      info.partition_names = partition_names;      
      gene_info.set_model_info(info);   
    } else gene_info.set_model_info(ModelInfoBlock::basic_block());
  }
}

bool GeneInputOutFile::check_consistency(vector<string>& parts, PvaluePart type) {
  if (type == pp_Partition) {
    for (int i = 0; i < parts.size(); i++) {
      if (parts[i] == "ALL-PART") continue;
      if (parts[i].substr(0,5) == "PART-" && parts[i].size() > 5) continue;
      return false;
    }
  } else if (type == pp_Analysis) {
    for (int i = 0; i < parts.size(); i++) {if (!ModelState::is_analysis(parts[i])) return false;}
  } else if (type == pp_Model) {
    for (int i = 0; i < parts.size(); i++) {if (!ModelState::is_model(parts[i])) return false;}  
  } else return false;
  return true;
}


void GeneInputCombine::load_input(const vector<string>& files, bool keep_mid) {filenames = files;
  if (filenames.empty()) {_LOG.error(msg_ident) << "no input files specified" << endl; die();} 
  if (filenames.size() == 1) {_LOG.error(msg_ident) << "only a single file specified" << endl; die();} 

  gene_info.init(20000);
  for (int i = 0; i < filenames.size(); i++) {
    GeneData* curr_data = (mode == RawMode) ? GeneInputRawFile(filenames[i]).eject_data() : GeneInputOutFile(filenames[i]).eject_data();
    VariableStorageAccess::AccessIterator<char*>* iterator = curr_data->name.get_data().get_iterator();
    while (!iterator->empty()) gene_info.get_iid(iterator->next_value());

    DataFrame::set_owner(curr_data, this); input_data.push_back(curr_data);
    delete iterator;
  }
  gene_info.build_index(); gene_info.location.init(); 

  map<int,set<GeneLocation> > mismatch; int overlap = 0, is_human = 0;
  for (int i = 0; i < input_data.size(); i++) {
    GeneData* curr_data = input_data[i]; 
    if (curr_data->location.human()) is_human++;
    for (int j = 0; j < curr_data->data_total(); j++) {
      if (!curr_data->active_id(j)) continue;
      int iid = gene_info[curr_data->name.get(j)], used = 0;
      if (iid >= 0) {GeneLocation& loc = curr_data->location.get(j); used++;
        if (gene_info.location.is_missing(iid)) gene_info.location.set(iid, loc);
        else if (gene_info.location.get(iid) != loc) {
          if (mismatch.find(iid) == mismatch.end()) {mismatch[iid] = set<GeneLocation>(); mismatch[iid].insert(gene_info.location.get(iid));}       
          mismatch[iid].insert(loc); loc = gene_info.location.get(iid);
        }   
      }
      if (used > 1) overlap++;      
    }
  }
 
  if (!allow_overlap && overlap > 0) {_LOG.error(msg_ident) << "input files contain overlapping genes" << endl; die();} 
  if (is_human < input_data.size()) {
    gene_info.location.human() = false;
    if (is_human > 0) _LOG << "WARNING: chromosome coding set to nonhuman mode for some but not all input files (setting output to nonhuman mode)" << endl;
  }

  if (!mismatch.empty()) {int counts[] = {0,0,0};
    _SLOG.set_block(msg_ident) << "# Following genes had mismatched locations across input files" << endl;
    _SLOG << "# Listing all locations found for each mismatched gene" << endl;
    for (map<int,set<GeneLocation> >::iterator it = mismatch.begin(); it != mismatch.end(); ++it) {
      set<GeneLocation>& locs = it->second; int type = 0; GeneLocation ref = *(locs.begin());
      _SLOG << gene_info.name.get(it->first) << ":";
      for (set<GeneLocation>::iterator jt = locs.begin(); jt != locs.end(); ++jt) {
        _SLOG << " " << *jt; 
        if (type < 2 && ref.compare(*jt) != 0) type = (ref.chromosome != jt->chromosome) ? 2 : 1;
      }
      _SLOG << endl;
      counts[type]++;
    }       

    string skip = (counts[2] > 0) ? Utils::space("ERROR - ") : "\t" + Utils::space("WARNING: ");      
    if (counts[2] > 0) _LOG.error(msg_ident) << "found " << Utils::plural(counts[2], "gene") << " with mismatched chromosomes across files" << endl; 
    else {
      _LOG << "WARNING: found " << Utils::plural(counts[0] + counts[1], "gene") << " with mismatched genomic locations across files (for " << Utils::plural(counts[0], "gene") << ", locations did overlap)" << endl;
      _LOG << skip << "using first location value encountered in input for each gene" << endl;
    }
    _LOG << skip << "writing list of mismatches to supplementary log file" << endl; 
    if (counts[2] > 0) die();
  }
  gene_info.data_sort(gene_info.location);
 
  access_map.set_size(gene_info.data_used(false), input_data.size()); 
  for (int i = 0; i < input_data.size(); i++) input_data[i]->align_map(gene_info, access_map[i], gene_info.data_used(false), keep_mid);
  if (!check_map()) {_LOG.error(msg_ident) << "data access map is invalid" << endl; die();}  
  
  /// data_sort leaves gene_info data storage consecutive, with iid range from 0 to gene_info.data_used()
  /// align_map leaves input_data's with consecutive storage, with access map as index to storage 
  /// input_data's have MASTER_IID variable attached for reverse mapping (these remain valid only until next reorder() call on gene_info)
}

bool GeneInputCombine::check_map() {
  for (int i = 0; i < access_map.nrow(); i++) {bool active = false;
    for (int j = 0; j < access_map.ncol(); j++) {if (access_map(i,j) >= 0) {active = true; break;}}
    if (!active || !gene_info.active_id(i)) return false;
  }

  for (int j = 0; j < access_map.ncol(); j++) {
    for (int i = 1; i < access_map.nrow(); i++) {
      int &curr = access_map(i,j), &prev = access_map(i-1,j);
      if (curr >= 0 && prev >= 0 && curr != (prev+1)) return false;
    }
  }
  return true;
}

void GeneInputCombine::process_pcol(bool align_core, bool allow_partitions) {int total = input_data.size();
  if (!allow_partitions) {bool dropped = false;
    for (int i = 0; i < total; i++) {ModelInfoBlock& curr = input_data[i]->get_model_info();
      if (curr.state[ModelInfoBlock::HasPartitions]) {curr.strip_partitions(); dropped = true;}
    }
    if (dropped) _LOG << "\tWARNING: partitioned analysis results will be skipped (please consult MAGMA manual about meta-analysis of partition p-values)" << endl; 
  }

  vector<ModelInfoBlock> input_info; bool equal = true;
  input_info.push_back(input_data[0]->get_model_info());  
  for (int i = 1; i < total; i++) {
    input_info.push_back(input_data[i]->get_model_info());
    if (input_info.front() != input_info.back()) equal = false; 
  }

  if (!equal) {
    vector<int> models(ModelState::_MODELCOUNT_,0), analyses(ModelState::_ANALYSISCOUNT_,0); int n_part = 0;
    for (int i = 0; i < total; i++) {
      set<ModelState::Model> model_set = input_info[i].get_models(); set<ModelState::Analysis> analysis_set = input_info[i].get_analyses();
      for (set<ModelState::Model>::iterator it = model_set.begin(); it != model_set.end(); ++it) models[*it]++;
      for (set<ModelState::Analysis>::iterator it = analysis_set.begin(); it != analysis_set.end(); ++it) analyses[*it]++;      
      if (input_info[i].state[ModelInfoBlock::HasPartitions]) n_part++;
    }
    
    bool aligned[] = {true, true, true}; 
    for (int i = 0; i < models.size(); i++) {if (models[i] > 0 && models[i] != total) {aligned[0] = false; break;}}
    for (int i = 0; i < analyses.size(); i++) {if (analyses[i] > 0 && analyses[i] != total) {aligned[1] = false; break;}}
    if (n_part > 0 && n_part < total) aligned[2] = false;

    if (align_core) {
      if (!aligned[0]) {
        for (int i = 0; i < total; i++) input_info[i].prune_models();
        _LOG << "\tWARNING: gene analysis models differ across files; mapping compatible columns on each other and discarding any remainder" << endl;
        aligned[0] = true;
      }
      if (!aligned[1] && analyses[ModelState::MainAnalysis] == total) {
        for (int i = 0; i < total; i++) input_info[i].filter(ModelState::MainAnalysis, true);
        _LOG << "\tWARNING: interaction analysis results are not available in all files; using main effect p-values only" << endl;        
        aligned[1] = true;
      }
    }
    if (!aligned[0] || !aligned[1]) {_LOG.error(msg_ident) << "p-value columns are not compatible across files" << endl; die();}

    vector<string> partition_names = input_info[0].partition_names;
    if (aligned[2]) {
      for (int i = 1; i < total; i++) {
        vector<string> &curr_names = input_info[i].partition_names, unmatched; vector<string>::iterator prev = curr_names.begin();
        for (int p = 1; p < curr_names.size(); p++) {      
          vector<string>::iterator found = find(partition_names.begin(), partition_names.end(), curr_names[p]);
          if (found != partition_names.end()) {
            if (!unmatched.empty()) {
              partition_names.insert(found > prev ? found : partition_names.end(), unmatched.begin(), unmatched.end());  
              unmatched.clear();          
            }
            prev = found;
          } else unmatched.push_back(curr_names[p]);
        }
        if (!unmatched.empty()) partition_names.insert(partition_names.end(), unmatched.begin(), unmatched.end());  
      }
      if (!Utils::vec_equals(partition_names, input_info[0].partition_names)) {
        for (int i = 0; i < total; i++) input_info[i].update_partitions(partition_names); 
      }
    } else {
      for (int i = 0; i < total; i++) input_info[i].strip_partitions();
      _LOG << "\tWARNING: partitioned analysis results are not available in all files; discarding all partition p-values" << endl;        
    }
    
    ModelInfoBlock info = input_info[0].copy().strip_partitions(), core = info;
    info.partition_names = partition_names;
    for (int p = 1; p < partition_names.size(); p++) info.append(core.set_partition(p));
    info.sort(true); gene_info.set_model_info(info); 
    for (int i = 0; i < total; i++) input_data[i]->set_model_info(input_info[i]); 
  } else gene_info.set_model_info(input_info[0]);
}

vector<string> GeneInputMerge::parse_filenames(string prefix, bool require) {
  pair<string,string> path = Utils::split_path(prefix);
  prefix = path.second + ".batch";
  string suffix = (mode == RawMode) ? ".genes.raw" : ".genes.out";
  string dir_name = (path.first != "") ? path.first : ".";

  bool force_numeric = settings["merge_batch_count"], force_chromosome = settings["merge_batch_chr"];
  int prefix_length = prefix.size(), suffix_length = suffix.size(), min_length = prefix_length + suffix_length + 3;
  int no_batches = settings.geti("merge_batch_count", 0); bool chr_mode = force_chromosome;

  DIR *dir; struct dirent *file; string filename; map<int,string> files;
  dir = opendir(dir_name.c_str());
  while ( (file = readdir(dir)) ) {
    if (strlen(file->d_name) < min_length) continue;
    filename.assign(file->d_name);
    if (filename.compare(0, prefix_length, prefix) != 0) continue;
    if (filename.compare(filename.size() - suffix_length, suffix_length, suffix) != 0) continue;
    
    pair<int,int> batch_index; int curr_batch;
    string batch_name = filename.substr(prefix_length, filename.size() - prefix_length - suffix_length);
    if (!Utils::convert_num_pair(batch_name, '_', batch_index)) {
      if (force_numeric) continue;
      pair<string,string> bits = Utils::split_str(batch_name, '_');
      if (bits.second != "chr") continue;
      if (!Utils::chr_val(bits.first, curr_batch, false)) {_LOG.error(msg_ident) << "chromosome code '" << bits.first << "' not recognised" << endl; die();} 
      chr_mode = true;
    } else {
      if (force_chromosome) continue;
      if (no_batches > 0 && no_batches != batch_index.second) {
        if (!force_numeric) {_LOG.error(msg_ident) << "found multiple sets of batches for prefix '" << prefix << "'; please specify number of batches in --batch flag" << endl; die();} 
        else continue;
      }
      curr_batch = batch_index.first;
      no_batches = batch_index.second;        
    }
    if (no_batches > 0 && chr_mode) {_LOG.error(msg_ident) << "found both numbered batches and chromosome batches for prefix '" << prefix << "'; please specify number of batches or chromosome mode in --batch flag" << endl; die();} 
    files[curr_batch] = (dir_name != ".") ? dir_name + filename : filename;
  }
  closedir(dir);

  string msg_add = "";
  if (force_numeric) msg_add = " and specified batch count";
  if (force_chromosome) msg_add = " and chromosome type";  
  if (require && files.empty()) {_LOG.error(msg_ident) << "found no batch files with prefix '" << prefix << "'" << msg_add << " to merge" << endl; die();} 
  if (no_batches > 0 && files.size() != no_batches) {_LOG.error(msg_ident) << "number of " << suffix << " batch files (" << files.size() << ") with prefix '" << prefix << "'" << msg_add << " does not match expected number of batches (" << no_batches << ")" << endl; die();} 
  
  if (chr_mode) {int count = 0;
    for (int i = 1; i < 23; i++) {if (files.find(i) != files.end()) count++;}
    if (count < 22) _LOG << "\tWARNING: found only " << count << " " << suffix << " batch files with prefix '" << prefix << "' for chromosomes 1-22" << endl;
  }
  
  vector<string> out;
  for (map<int,string>::iterator it = files.begin(); it != files.end(); ++it) out.push_back(it->second);
  return out;
}

bool GeneInputMerge::run(const string& prefix, bool require) {
  vector<string> filenames = parse_filenames(prefix, require);
  if (!filenames.empty()){ 
    _LOG << "Merging .genes." << ((mode == RawMode) ? "raw" : "out") << " files with prefix '" << prefix << "'... " << endl;
    load_input(filenames); ///yields consecutive gene_info and input_data, iid runs up to gene_info.data_used()

    merge_variable(gene_info.nsnps, "NSNPS");
    merge_variable(gene_info.nrare, "NRARE");
    merge_variable(gene_info.nparam, "NPARAM");
    merge_variable(gene_info.nsamp, "N");
    merge_variable(gene_info.ndata, "DATASETS");    
    merge_variable(gene_info.mac, "MAC");    
    merge_variable(gene_info.rsq, "RSQ");    
    merge_variable(gene_info.zstat, "ZSTAT"); 
    bool permp = merge_variable(gene_info.perm_pval, "PERMP"); 
    if (permp) merge_variable(gene_info.perm_count, "NPERM");     

    if (mode == OutMode) {
      if (check_variable(gene_info.pval, "PVAL", false)) {
        process_pcol(false, true);
        merge_aligned_variable(gene_info.pval, "PVAL", true);        
        merge_aligned_variable(gene_info.nsnps_part, "NSNPS_PART", false);                
        merge_aligned_variable(gene_info.nsnps_part, "NSAMP_PART", false);                        
      } 
    } else merge_corrs();

    GeneOutput output(gene_info, settings.gets("out_prefix")); output.set_indent(false);
    if (mode == RawMode) output.write_rawfile();
    else output.write_outfile();
    
    _LOG << endl; return true;
  } else return false;
}

void GeneInputMerge::merge_corrs() {
  vector<BaseBuffer<BaseBuffer<double>*>*> data = get_data(gene_info.corrs, "CORRS");
  if (data.empty()) {_LOG.error(msg_ident) << "gene correlation information not available for all files" << endl; die();}  

  gene_info.corrs.init(); int curr_col = -1;
  for (int i = 0; i < gene_info.data_used(false); i++) {
    bool col_change = (i == 0) || (access_map(i,curr_col) < 0);
    if (col_change) {for (int j = 0; j < access_map.ncol(); j++) {if (access_map(i,j) >= 0) {curr_col = j; break;}}}
    
    BaseBuffer<double>* corrs = data[curr_col]->get(access_map(i,curr_col));
    if (corrs) {
      if (col_change && !corrs->is_empty()) {_LOG.error(msg_ident) << "encountered invalid gene correlation block" << endl; die();}         
      gene_info.corrs.assign(i, *corrs);
    }
  }
}

void GeneInputMeta::run() {
  _LOG << "Running fixed-effects meta-analysis..." << endl;
  
  vector<string> filenames = settings.getvs("gene_meta_files", true);
  if (filenames.empty()) {_LOG.error(msg_ident) << "no files to process" << endl; die();}
  if (filenames.size() == 1) {_LOG.error(msg_ident) << "only a single file to process" << endl; die();}  

  load_input(filenames, mode == RawMode); ///yields consecutive gene_info and input_data, iid runs up to gene_info.data_used()
  set_weights();

  process_variable(gene_info.nsnps, "NSNPS", GeneInput::CeilingMean);
  process_variable(gene_info.nrare, "NRARE", GeneInput::CeilingMean);
  process_variable(gene_info.nparam, "NPARAM", GeneInput::CeilingMean);
  process_variable(gene_info.nsamp, "N", GeneInput::Sum);

  for (int i = 0; i < input_data.size(); i++) {if (!input_data[i]->ndata.initialised()) {input_data[i]->ndata.init(); input_data[i]->ndata.get_buffer().assign_value(1);}}
  process_variable(gene_info.ndata, "DATASETS", GeneInput::Sum);    

  process_variable(gene_info.mac, "MAC", GeneInput::Mean);  
  process_variable(gene_info.rsq, "RSQ", GeneInput::Mean);    

  if (mode == OutMode) {bool has_pval = false;
    for (int i = 0; i < input_data.size(); i++) {GeneData& curr = *input_data[i];
      if (!curr.using_variable("PVAL") && curr.using_variable("PERMP")) curr.upgrade_permp();
    }
    
    if (check_variable(gene_info.pval, "PVAL", false)) {
      process_pcol(true, settings["metagene_partitions"]);
      set_nsamp_part();

      has_pval = process_aligned_variable(gene_info.pval, "PVAL", GeneInput::WeightedPval, true); 
       
      process_aligned_variable(gene_info.nsnps_part, "NSNPS_PART", GeneInput::CeilingMean, false);                
      process_aligned_variable(gene_info.nsamp_part, "NSAMP_PART", GeneInput::Sum, false);                      
    } 
    if (has_pval) gene_info.create_zstat();
    else {_LOG.error(msg_ident) << "unable to perform meta-analysis; p-values are not available or types do not match across files" << endl; die();}    
  } else {
    process_variable(gene_info.zstat, "ZSTAT", GeneInput::WeightedZstat);    
    process_corrs();
  }

  GeneOutput output(gene_info, settings.gets("out_prefix")); output.set_indent(false);
  if (mode == RawMode) output.write_rawfile();
  else output.write_outfile();
  _LOG << endl;
}

void GeneInputMeta::process_corrs() {int input_total = input_data.size();
  vector<int*> reverse_map; vector<double*> weights = get_weights();
  for (int i = 0; i < input_total; i++) {
    NumericID<int>* id_var = input_data[i]->typed_variable<NumericID<int> >("MASTER_IID");
    if (!id_var) {_LOG.error(msg_ident) << "for file '" << filenames[i] << "', could not find master IID map" << endl; die();} 
    reverse_map.push_back(id_var->get_buffer().data());
  }
  
  gene_info.corrs.init(); deque<pair<int,int> > impute; int prev_length = 0;
  for (int m_iid = gene_info.data_used(false)-1; m_iid >= 0; m_iid--) { ///m_iid is guaranteed to be valid from 0 to gene_info.data_used()
    vector<BaseBuffer<double>*> input_corrs(input_total, (BaseBuffer<double>*)0); int tot_length = 0; 
    for (int i = 0; i < input_total; i++) {int s_iid = access_map(m_iid, i);
      if (s_iid >= 0) {BaseBuffer<double>& curr = input_data[i]->corrs.get_buffer(s_iid);
        if (!curr.is_empty()) {
          int length = m_iid - reverse_map[i][s_iid-curr.size()];
          if (length > tot_length) tot_length = length;
          input_corrs[i] = &curr;          
        }        
      }
    }
    if (tot_length < prev_length-1) tot_length = prev_length - 1;
    prev_length = tot_length;

    if (tot_length > 0) {
      gene_info.corrs.set_length(m_iid, tot_length);
      double* corrs = gene_info.corrs.get_buffer(m_iid).data(); 
      vector<double> corr_weights(tot_length, 0);      
      int base_mid = m_iid - tot_length;

      for (int i = 0; i < input_total; i++) {int s_iid1 = access_map(m_iid, i);                
        if (input_corrs[i]) {BaseBuffer<double>& curr = *input_corrs[i]; 
          double wt_base = weights[i][s_iid1]; int offset = s_iid1 - curr.size();
          for (int j = 0; j < curr.size(); j++) {
            int s_iid2 = offset + j, target = reverse_map[i][s_iid2] - base_mid;
            double wt = wt_base * weights[i][s_iid2];
            corrs[target] += curr[j]*wt; corr_weights[target] += wt;
          }
        }
      }

      for (int i = tot_length-1; i >= 0; i--) {
        if (corr_weights[i] > 0) corrs[i] /= corr_weights[i];
        else impute.push_front(pair<int,int>(base_mid+i, m_iid)); 
      }
    }
  }

  if (!impute.empty()) {  
    GeneCorrelationAccess corrs(gene_info.corrs.get_buffer()); vector<double> corr_storage(impute.size(), 0);
    for (int i = 0; i < impute.size(); i++) {
      int gene_A = impute[i].first, gene_B = impute[i].second; double& target = corr_storage[i]; 
      for (int gene_C = gene_A - corrs.count(gene_A); gene_C < gene_info.data_used(false); gene_C++) {
        if (gene_C == gene_A || gene_C == gene_B) continue;        
        double cor = corrs.get(gene_A, gene_C) * corrs.get(gene_B, gene_C);      
        if (cor > target) target = cor;
      }
    }
    for (int i = 0; i < impute.size(); i++) corrs.get(impute[i].first, impute[i].second) = corr_storage[i];
  }
}

void GeneInputMeta::set_weights() {clear_weights();
  const vector<Triple<double> >& size_override = settings.get_metaweights(); 
  bool use_override = (size_override.size() == input_data.size());
  if (!use_override && !size_override.empty()) {_LOG.error(msg_ident) << "number of user-specified sample sizes does not correspond to number of input files" << endl; die();}
  
  for (int i = 0; i < input_data.size(); i++) {GeneData& data = *input_data[i];
    NonNegativeNumber<double>* wt = new NonNegativeNumber<double>();
    data.link_variable(*wt);
    if (use_override && size_override[i][0] > 0) {
      vector<double> weights(Utils::chr_type(Utils::chr_max), sqrt(size_override[i][2])); weights[0] = sqrt(size_override[i][0]);
      weights[Utils::chr_type(Utils::chrX_code)] = weights[Utils::chr_type(Utils::chrZ_code)] = sqrt(size_override[i][1]);
      for (int j = 0; j < data.data_total(); j++) {if (data.active_id(j)) wt->set(j, weights[Utils::chr_type(data.location.get(j).chromosome)]);}      
    } else if (data.nsamp.initialised()) {
      for (int j = 0; j < data.data_total(); j++) {if (data.active_id(j)) wt->set(j, sqrt(data.nsamp.get(j)));}      
    } else {_LOG.error(msg_ident) << "file '" << filenames[i] << "' contains no sample size column, and meta-analysis weights were not provided by user" << endl; die();} 
    weights.push_back(wt);
  }
  
  data_corrs.assign(settings.get_metacorrs());
  if (!data_corrs.empty()) {
    bool valid = (data_corrs.ncol() == data_corrs.nrow() && data_corrs.ncol() == input_data.size());
    if (valid) {
      for (int i = 0; i < data_corrs.ncol(); i++) {
        if (data_corrs(i,i) != 1) {valid = false; break;}
        for (int j = i+1; j < data_corrs.ncol(); j++) {
          if (data_corrs(j,i) != data_corrs(i,j) || data_corrs(j,i) > 1 || data_corrs(i,j) < 0) {valid = false; break;}
        } 
      }
    }
    if (!valid) {_LOG.error(msg_ident) << "specified data correlation matrix is invalid" << endl; die();}
  }
}

vector<double*> GeneInputMeta::get_weights() {vector<double*> out;
  for (int i = 0; i < weights.size(); i++) out.push_back(weights[i]->get_buffer().data());
  return out;
} 

void GeneInputMeta::clear_weights() {
  for (int i = 0; i < weights.size(); i++) {
    input_data[i]->unlink_variable(*weights[i]); 
    delete weights[i];
  } 
  weights.clear();  
}

void GeneInputMeta::set_nsamp_part() {
  ModelInfoBlock& info = gene_info.get_model_info();
  if (info.state[ModelInfoBlock::HasPartitions] && !info.partition_index(false).empty() && check_variable(gene_info.nsamp, "N") && check_variable(gene_info.nsnps_part, "NSNPS_PART", false)) {
    for (int i = 0; i < input_data.size(); i++) {
      Buffer<int> &nsnps = input_data[i]->nsnps_part.get_buffer(), &nsamp = input_data[i]->nsamp_part.get_buffer(); BaseBuffer<int>& nbuff = input_data[i]->nsamp.get_buffer();
      input_data[i]->nsamp_part.init(nsnps.ncol());
      for (int i = 0; i < nbuff.size(); i++) {int& N = nbuff[i];
        for (int j = 0; j < nsamp.ncol(); j++) {if (nsnps(i,j) > 0) nsamp(i,j) = N;}
      }       
    }
  }
}

double& GeneInputMeta::GeneCorrelationAccess::get(int gene_i, int gene_j) {
  if (gene_j < gene_i) swap(gene_i, gene_j);
  if (gene_i != gene_j && gene_j < max_index && corrs[gene_j]) {BaseBuffer<double>& curr = *corrs[gene_j];
    int offset = curr.size() - (gene_j - gene_i);
    if (offset >= 0) return curr[offset];
  }
  dummy_value = 0; return dummy_value;
}

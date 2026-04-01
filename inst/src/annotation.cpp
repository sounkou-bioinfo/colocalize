/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "annotation.h"

void Annotation::load_geneloc(const string& filename) {
  _LOG << "Reading gene locations from file " << filename << "... " << endl;
  chromosomes.reserve(27);
  for (int i = 0; i <= 22; i++) {
    chromosomes.push_back(Chromosome(i));
    chr_map[i] = i;
  }

  map<string, long>::const_iterator iter;

  string geneid(10, '\0');
  string chr_str(5, '\0');
  string begin_str(15, '\0');
  string end_str(15, '\0');  
  string strand_str(3, '\0');
  long window_up = settings.geti("window_up"), window_down = settings.geti("window_down");
  int chr, curr_chr = -1; Chromosome* active_chrom = 0; unsigned long begin, end, max_val = numeric_limits<unsigned long>::max();
  bool human_chr = !settings["nonhuman_chr"];
  unsigned long curr_up, curr_down;


  bool do_window = window_up > 0 || window_down > 0;
  bool read_strand = do_window && (window_up != window_down) && !settings["ignore_strand"];

  if (do_window) {
    _LOG << "\tadding window: ";
    if (window_down == window_up) {
      _LOG << window_up << "bp" << endl;
    } else {
      _LOG << window_up << "bp (before), " << window_down << "bp (after)" << endl;    
    }
  }

  string buffer(1000, '\0');
  ifstream gen(filename.c_str());
  istringstream istr, convert;

  char line_end = '\n';
  while (gen >> geneid) {
    if (gen.peek() == '\n') break;
    if (gen.peek() == '\r') {
      gen.get();
      if (gen.peek() != '\n') line_end = '\r';
      break;
    }
  }
  gen.seekg(0); gen.clear();

  long line_no = 0; no_genes = 0;
  bool strand_flip = false, line_fail;
  while (getline(gen, buffer, line_end)) {
    line_no++;
    istr.clear(); 
    istr.str(buffer);

    if (!(istr >> geneid)) continue;

    line_fail = (!(istr >> chr_str) || !(istr >> begin_str) || !(istr >> end_str));
    if (read_strand) {
      if (istr >> strand_str) {
        if (strand_str == "+") strand_flip = false;
        else if (strand_str == "-") strand_flip = true;
        else {
          if (no_genes == 0) read_strand = false;
          else {
            _LOG << "\tWARNING: on line " << line_no << ", strand code '" << strand_str << "' not recognised; skipping gene (ID = " << geneid << ")" << endl;
            continue;
          }
        }
      } else {
        if (no_genes == 0) read_strand = false;
        else line_fail = true;
      }
      if (!read_strand) {
        strand_flip = false;
        _LOG << "\tWARNING: no valid strand column found, assuming positive strand for all genes" << endl;
      }
    }
    if (line_fail) {_LOG.error("reading gene location file") << "too few values on line " << line_no << ":" << endl << "\tline: " << buffer << endl; die();} 
    
    if (!Utils::chr_val(chr_str, chr, human_chr)) {
      _LOG << "\tWARNING: on line " << line_no << ", chromosome code '" << chr_str << "' not recognised; skipping gene (ID = " << geneid << ")" << endl;
      continue;
    }
    
    if (chr == 0) continue;

    convert.clear(); convert.str(begin_str); convert >> begin;    
    if (!convert.eof() || begin_str[0] == '-') {
      _LOG << "\tWARNING: on line " << line_no << ", third value must be a positive integer; skipping gene (ID = " << geneid << ")" << endl;
      continue;
    }

    convert.clear(); convert.str(end_str); convert >> end;    
    if (!convert.eof() || end_str[0] == '-') {
      _LOG << "\tWARNING: on line " << line_no << ", fourth value must be a positive integer; skipping gene (ID = " << geneid << ")" << endl;
      continue;
    }

    if (end == max_val) {
      _LOG << "\tWARNING: on line " << line_no << ", detected integer overflow in gene boundaries; skipping gene (ID = " << geneid << ")" << endl;
      continue;
    }

    if (end <= begin) {
      _LOG << "\tWARNING: on line " << line_no << ", gene start must be smaller than gene end (from " << begin << ", to " << end << "); skipping gene (ID = " << geneid << ")" << endl;
      continue;
    }
    
    if (gene_names.find(geneid) != gene_names.end()) {_LOG.error("reading gene location file") << "duplicate gene entry on line " << line_no << " (ID = " << geneid << ")" << endl; die();} 
    if (do_window) {
      curr_up = strand_flip ? window_down : window_up;
      curr_down = strand_flip ? window_up : window_down;

      begin = begin < curr_up ? 0 : begin - curr_up;
      end = max_val - end < curr_down ? max_val : end + curr_down;
    }

    gene_names.insert(geneid);
    if (chr != curr_chr) {
      if (chr_map.find(chr) == chr_map.end()) {
        chr_map[chr] = chromosomes.size();
        chromosomes.push_back(Chromosome(chr));
      }
      curr_chr = chr;
      active_chrom = &(chromosomes[chr_map[curr_chr]]);
      if (curr_chr < Utils::chr_reserved && curr_chr > max_chr) max_chr = curr_chr;
    }

    active_chrom->add_gene(geneid, begin, end);
    no_genes++;
  }
  
  if (no_genes == 0) {_LOG.error("reading gene location file") << "found no valid genes in file" << endl; die();} 
  
  _LOG << "\t" << no_genes << " gene locations read from file" << endl;
  
  int pad = max(Utils::num_length(max_chr), 2);
  for(map<int,int>::iterator it = chr_map.begin(); it != chr_map.end(); ++it) {
    int chr_index = it->second;
    if (!chromosomes[chr_index].empty) {
      _LOG << "\tchromosome " << Utils::chr_string(it->first, pad);
      _LOG << ": " << chromosomes[chr_index].size() << " genes" << endl;
      if (!chromosomes[chr_index].empty) chromosomes[chr_index].prepare();
    }
  }
}

void Annotation::load_filter(const string& filename) {
  _LOG << "Reading SNP IDs to filter from file " << filename << "... " << endl;

  ifstream filter(filename.c_str());
  istringstream istr;

  string line(200, '\0');
  string rsid(1, '\0');
  long no_snps = 0; bool skip_col = filename.size() > 4 && (filename.substr(filename.size()-4) == ".bim");
  string tmp(filename.substr(filename.size()-4));

  while (getline(filter, line)) {
    istr.clear();
    istr.str(line);
    if (!(istr >> rsid)) continue;

    if (skip_col) {
      istr >> rsid; 
      if (istr.fail()) continue;
    }

    snp_filter.insert(rsid);
    no_snps++;
  }
  if (no_snps == 0) {_LOG.error("reading SNP filter file", 2) << "no SNP IDs in file" << endl; die();} 

  _LOG << "\t" << no_snps << " SNP IDs read" << endl;
}

void Annotation::set_filter(const map<string, long>* snp_map) {
  if (snp_filter.size() == 0) {
    _LOG << "Setting SNP ID filter using loaded genotype data... " << endl; 
    for (map<string,long>::const_iterator iter = snp_map->begin(); iter != snp_map->end(); ++iter) {
      snp_filter.insert(iter->first);
    } 
    _LOG << "\t" << snp_filter.size() << " SNP IDs set" << endl;    
  } else {
    _LOG << "Pruning SNP ID filter using loaded genotype data... " << endl; 
    set<string> pruned;
    for (map<string,long>::const_iterator iter = snp_map->begin(); iter != snp_map->end(); ++iter) {
      if (snp_filter.find(iter->first) == snp_filter.end()) continue;
      pruned.insert(iter->first);
    } 
    snp_filter.swap(pruned);
    _LOG << "\t" << snp_filter.size() << " SNP IDs remaining" << endl;    
  }
}


void Annotation::load_snploc(const string& filename) {
  _LOG << "Reading SNP locations from file " << filename << "... " << endl;

  string rsid(10, '\0');
  string chr_str(5, '\0');
  string loc_str(15, '\0');
  string back_up;
  int chr, curr_chr = -1; Chromosome* active_chrom = 0; unsigned long loc, max_val = numeric_limits<unsigned long>::max();
  bool do_filter = snp_filter.size() > 0;
  bool is_bim = filename.size() > 4 && (filename.substr(filename.size()-4) == ".bim");
  bool human_chr = !settings["nonhuman_chr"];

  string buffer(1000, '\0');
  ifstream snp(filename.c_str());
  istringstream istr, convert;

  clock_t last = clock(), curr;
  long line_no = 0, no_snps = 0, no_assig = 0;
  while (getline(snp, buffer)) {
    line_no++;
    istr.clear(); 
    istr.str(buffer);

    if (is_bim) {
      if (!(istr >> chr_str)) continue;
  
      if (!(istr >> rsid) || !(istr >> loc_str) || !(istr >> loc_str)) {_LOG.error("reading SNP location file") << "too few values on line " << line_no << ":" << "\tline: " << buffer << endl; die();} 
    } else {
      if (!(istr >> rsid)) continue;
  
      if (!(istr >> chr_str) || !(istr >> loc_str)) {_LOG.error("reading SNP location file") << "too few values on line " << line_no << ":" << "\tline: " << buffer << endl; die();} 
    } 

    if (do_filter && snp_filter.find(rsid) == snp_filter.end()) continue;
    
    if (!Utils::chr_val(chr_str, chr, human_chr)) {
      _LOG.screen() << back_up << _LOG.all(); back_up = "";
      _LOG << "\tWARNING: on line " << line_no << ", chromosome code '" << chr_str << "' not recognised; skipping SNP (ID = " << rsid << ")" << endl;
      continue;
    }
  
    if (chr != curr_chr) {
      if (chr == 0 || chr_map.find(chr) == chr_map.end()) {
        continue;
      }
      curr_chr = chr;
      active_chrom = &(chromosomes[chr_map[curr_chr]]);
    }

    convert.clear(); convert.str(loc_str); convert >> loc;    
    if (!convert.eof() || loc_str[0] == '-') {
      _LOG.screen() << back_up << _LOG.all(); back_up = "";
      _LOG << "\tWARNING: on line " << line_no << ", SNP location must be a positive integer; skipping SNP (ID = " << rsid << ")" << endl;
      continue;
    }

    if (loc == max_val) {
      _LOG.screen() << back_up << _LOG.all(); back_up = "";
      _LOG << "\tWARNING: on line " << line_no << ", detected integer overflow in SNP location; skipping SNP (ID = " << rsid << ")" << endl;
      continue;
    }

    no_snps++;      
    no_assig += active_chrom->add_snp(rsid, loc);
    
    curr = clock();
    if (curr - last >= CLOCKS_PER_SEC) {
      last = curr;
      _LOG.counter(true) << "\tSNPs mapped so far: " << no_assig;
      _LOG.counter(false);
    }
  }

  _LOG.counter(true) << "                                                                                                 ";
  _LOG.counter(false);

  if (no_snps == 0) {_LOG.error("reading SNP location file") << "found no valid SNPs in file" << endl; die();} 
  
  if (no_assig == 0) {_LOG.error("reading SNP location file") << "no SNPs mapped to a gene" << endl; die();} 

  _LOG.screen() << back_up << _LOG.all();
  _LOG << "\t" << no_snps << " SNP locations read from file" << endl;
  _LOG << "\tof those, " << no_assig << " (" << setprecision(4) << (float(100*no_assig) / no_snps) <<"%) mapped to at least one gene" << endl;
}

void Annotation::write_annot() {
  _LOG << "Writing annotation to file " << outname << endl;
  DelimitedOutput fout(outname, '\t');
  
  fout.print_param("window_up", settings.geti("window_up"));
  fout.print_param("window_down", settings.geti("window_down"));
  if (settings["nonhuman_chr"]) fout.print_param("nonhuman", 1);
  
  int no_empty, tot_empty = 0;
  int pad = max(Utils::num_length(max_chr), 2);
  for(map<int,int>::iterator it = chr_map.begin(); it != chr_map.end(); ++it) {
    int chr_index = it->second;
    if (!chromosomes[chr_index].empty) {
      no_empty = chromosomes[chr_index].write(fout);
      if (no_empty > 0) {
        _LOG << "\tfor chromosome " << Utils::chr_string(it->first, pad) << ", " << no_empty;
        if (no_empty == 1) _LOG << " gene is";
        else _LOG << " genes are";
        _LOG << " empty (out of " << chromosomes[chr_index].size() << ")" << endl;
        tot_empty += no_empty;
      }
    }
  }
  if (tot_empty > 0) {
    long no_mapped = no_genes - tot_empty;
    _LOG << "\tat least one SNP mapped to each of a total of " << no_mapped;
    if (no_mapped == 1) _LOG << " gene";
    else _LOG << " genes";
    _LOG << " (out of " << no_genes << ")" << endl;
  }
}


/*********************************************
* Chromosome helper class                    *
*********************************************/

void Chromosome::add_gene(const string& name, const unsigned long& begin, const unsigned long& end) {
  names.push_back(name);
  locations.push_back(GeneLocation(chr_index, begin, end));
  no_genes++;
  empty = false;
}

bool Chromosome::add_snp(const string& rsid, const unsigned long& loc) {
  if (chr_location.compare(loc) != 0) return false;
  if (loc < prev_snp_loc) {
    curr_range = (loc - chr_location.begin) * size_frac;
    while (curr_range > 0 && ranges[curr_range-1].compare(loc) <= 0) curr_range--;
    curr_range = 0;
  }
  prev_snp_loc = loc;
  
  ///starting at current range, move to first range NOT fully before SNP location
  int comp, id;
  while ((comp = ranges[curr_range].compare(loc)) > 0) curr_range++;

  if (comp == 0) {
    if (range_contents[curr_range] < 0) {
      id = -range_contents[curr_range] - 1;
      for (deque<int>::iterator iter = range_multi_contents[id].begin(); iter != range_multi_contents[id].end(); ++iter) {
        snps[*iter].push_back(rsid);
      }
    } else {
      snps[range_contents[curr_range]].push_back(rsid);
    }
    return true;
  }
  return false;
}


void Chromosome::prepare() {
  if (empty) {_LOG.error("processing genes") << "trying to access empty chromosome" << endl; die();} 

  vector<long> index; index.reserve(no_genes);
  for (long i = 0; i < no_genes; i++) index.push_back(i);
  
  GeneSorter gs(locations); ///Sorted by loc.begin, then loc.end
  gs.run(index);

  vector<string> re_name; re_name.reserve(no_genes);
  vector<GeneLocation> re_loc; re_loc.reserve(no_genes);
  for (long i = 0; i < no_genes; i++) {
    re_name.push_back(names[index[i]]);
    re_loc.push_back(locations[index[i]]);
  }
  names.swap(re_name);
  locations.swap(re_loc);

  ranges.reserve(no_genes);
  range_contents.reserve(no_genes);

  deque<unsigned long> bits;   ///alleen ends bewaren, begin is steeds hetzelfde
  deque<int> bit_contents;
     
  int multi_index = -1; unsigned long bits_begin = 0, bits_end, end; bool is_last;
  for (int i_gene = 0; i_gene < no_genes; i_gene++) {
    is_last = i_gene == no_genes-1;
    if (bits.empty()) {
      if (is_last || locations[i_gene].before(locations[i_gene+1])) {
        ranges.push_back(locations[i_gene]);
        range_contents.push_back(i_gene);
        continue;
      } else {
        bits.push_back(locations[i_gene].end);
        bit_contents.push_back(i_gene);
        bits_begin = locations[i_gene].begin;
      }
    } else {
      end = locations[i_gene].end;
      bits.push_back(end);
      bit_contents.push_back(i_gene);
      for (int i = bits.size()-1; i > 0; i--) {
        if (end >= bits[i-1]) break;
        bits[i] = bits[i-1]; bits[i-1] = end;
        bit_contents[i] = bit_contents[i-1]; bit_contents[i-1] = i_gene;
      }
    }  
      
    bits_end = is_last ? bits.back()+1 : locations[i_gene+1].begin;
    while (bits.size() > 0) {
      if (bits.front() >= bits_begin) {
        if (bits_end <= bits_begin) break;
        if (bits.size() == 1) range_contents.push_back(bit_contents.front());
        else {
          range_contents.push_back(multi_index--);
          range_multi_contents.push_back(bit_contents);
        }
        end = min(bits.front(), bits_end-1);
        ranges.push_back(GeneLocation(chr_index, bits_begin, end));
        bits_begin = end+1;
        if (bits.front() >= bits_end) break;
      }
      bits.pop_front();
      bit_contents.pop_front();
    }
  }

  no_ranges = ranges.size();
  chr_location = GeneLocation(chr_index, ranges.front().begin, ranges.back().end);
  size_frac = float(no_genes - 1) / (chr_location.end - chr_location.begin);

  snps.resize(no_genes);
}


int Chromosome::write(DelimitedOutput& fout) {
  int empty_genes = 0;
  for (int i = 0; i < no_genes; i++) {
    if (snps[i].empty()) {
      empty_genes++;
      continue;
    }
    fout << names[i] << locations[i];
    for (vector<string>::iterator iter = snps[i].begin(); iter != snps[i].end(); ++iter) {
      fout << *iter;
    }
    fout << endl;
  }
  return empty_genes;
}

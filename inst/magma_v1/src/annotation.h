/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef ANNOTATION_H
#define ANNOTATION_H

#include <vector>
#include <deque>
#include <set>

#include "parse.h"
#include "output.h"
#include "utils.h"
#include "geneutils.h"


class Chromosome {
  int chr_index;
  float size_frac;
  GeneLocation chr_location;
  long no_genes;
  long no_ranges;
  
  unsigned long prev_snp_loc;
  long curr_range;
     
  vector<string> names;
  vector<GeneLocation> locations;
  vector<vector<string> > snps;
  
  vector<GeneLocation> ranges;
  vector<int> range_contents;
  vector<deque<int> > range_multi_contents;

public:
  bool empty;
  
  Chromosome(int chr=0) : chr_index(chr), no_genes(0), no_ranges(0), prev_snp_loc(numeric_limits<unsigned long>::max()), empty(true) {}

  void add_gene(const string& name, const unsigned long& begin, const unsigned long& end);
  bool add_snp(const string& name, const unsigned long& loc);
  void prepare();
  int write(DelimitedOutput& fout);
  
  long size() {return no_genes;}
};


class Annotation {
  const Settings& settings;
  long no_genes;
  int max_chr;
  string outname;

  vector<Chromosome> chromosomes;
  map<int,int> chr_map;
  set<string> gene_names;
  set<string> snp_filter; 

  void load_geneloc(const string& filename);
  void set_filter(const map<string, long>* snp_map);
  void load_filter(const string& filename);  
  void load_snploc(const string& filename);
  void write_annot();

public:                   
  Annotation(const Settings& s, const map<string,long>* snp_filter=0) : settings(s), max_chr(0) {
    _LOG << "Starting annotation..." << endl;
    outname = settings.gets("out_prefix");
        
    load_geneloc(settings.gets("gene_loc_file"));
    if (settings["annot_filter_file"]) {
      if (settings.gets("annot_filter_file") != settings.gets("snp_loc_file")) {
        if (!settings["do_analysis"] || settings.gets("annot_filter_file") != settings.get_plinkfile("bim")) {
          load_filter(settings.gets("annot_filter_file"));
        }
      }   
    }
    outname = win_txt(outname + ".genes.annot");
    if (snp_filter && settings.gets("snp_loc_file") != settings.get_plinkfile("bim")) set_filter(snp_filter);
    load_snploc(settings.gets("snp_loc_file"));

    write_annot();
    _LOG << endl;
  }
};




#endif /**ANNOTATION_H*/


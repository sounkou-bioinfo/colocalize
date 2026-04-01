/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENEANALYSIS_H
#define GENEANALYSIS_H

#include <deque>
#include "genestats.h"
#include "permutation.h"
#include "genemodelengine.h"
#include "statutils.h"
#include "geneutils.h"


class GeneAnalysis {
  Settings& settings;
  BaseInput* baseinput; 

  GeneStatFactory statfactory;
  GeneStats* data;

  SNPData& snp_info;
  GeneData& gene_info;

  bool compute_perm; 
  bool compute_corrs; 
  bool dump_data; 
  string dump_pref;

  GeneModelEngine* engine; 
  vector<GeneModelEngine*> engine_bay;

  Permutation* perm_engine;
  GeneCorrelations corr_engine;
  BaseBuffer<double> corr_buffer;


  void print_model();
  
  void load_engine(); 
  ModelState::Model get_model_type(string model);
  GeneModelEngine* make_core_engine(ModelState::Model type);
  GeneModelEngine* make_multi_engine();  
  PartitionedEngine* blockwise(CoreEngine* core, Aggregator::AggregateType aggregate_type=Aggregator::Mean);  
  PrinCompEngine* princomp(PrinCompEngine* pc);

  void run_analysis();
  void analyse_gene(int geneid, int exp_pval);

  int set_chromosome();
  bool set_chromosome(int chr);
  bool set_data(int chr);

  float get_effsize(bool adj) {return max(engine->get_param_value(adj ? "eff_size_adjusted" : "eff_size", 0), 0.0);}
  float get_zstat(int index=0);
  void load_snpcount(int geneid, int expected);

public:
  GeneAnalysis(Settings& s, BaseInput* bi) : settings(s), baseinput(bi), statfactory(s, bi), data(0), snp_info(bi->snp_info), gene_info(bi->gene_info), perm_engine(0) {
    CorrelationUtils::set_top_param(settings);                                               
    compute_corrs = !settings["genes_only"];
    compute_perm = settings["compute_permp"];
    dump_data = settings["dump_data"];
    
    if (compute_corrs) corr_engine.set_range(settings.geti("corr_range"));
    if (dump_data) dump_pref = settings.gets("dump_data");
  }
  ~GeneAnalysis(){
    for (short i = 0; i < engine_bay.size(); i++) delete engine_bay[i];
    delete perm_engine;
  }
  
  void run();
};


#endif /** GENEANALYSIS_H */


/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "geneanalysis.h"
#include "geneoutput.h"

void GeneAnalysis::load_engine() {
  ModelState::Model model_type = get_model_type(settings.gets("gene_model"));
  engine_bay.assign(ModelState::_MODELCOUNT_, 0);
  if (model_type == ModelState::MultiModel) engine = make_multi_engine();
  else engine = make_core_engine(model_type);  
}

ModelState::Model GeneAnalysis::get_model_type(string model) {
  if (model == "pcreg") return ModelState::PrinCompReg;
  if (model == "snpwise") {string submodel = settings.gets("gene_snpwise");
    if (submodel == "top") {
      if (!settings["gene_snpwise_top_fraction"] && settings.geti("gene_snpwise_top_count", 1) == 1) return ModelState::SNPwiseModelTopOne;
      else return ModelState::SNPwiseModelTop;
    } else {
      if (submodel == "brown") return ModelState::SNPwiseModelMeanBrown;
      else return ModelState::SNPwiseModelMeanImhof;
    }
  }

  if (model == "snpwise_mean") return ModelState::SNPwiseModelMeanImhof;
  if (model == "snpwise_top1") return ModelState::SNPwiseModelTopOne;
  if (model == "varimax") return ModelState::VariPow;
  if (model == "tdt") return ModelState::TDTModel;
  if (model == "multi") return ModelState::MultiModel;

  return ModelState::UnknownModel;
}

GeneModelEngine* GeneAnalysis::make_core_engine(ModelState::Model type) {
  if (!engine_bay[type]) {
    switch (type) {
      case ModelState::PrinCompReg:
        engine_bay[type] = blockwise(princomp(new LinearRegression(settings, baseinput->get_npheno())));
        break;
      case ModelState::SNPwiseModelTop:
        engine_bay[type] = new SNPwiseTop(settings);
        break;
      case ModelState::SNPwiseModelTopOne:
        engine_bay[type] = blockwise(new SNPwiseTopOne(settings), Aggregator::Top);
        break;
      case ModelState::SNPwiseModelMeanBrown:
        engine_bay[type] = blockwise(new SNPwiseMeanBrown(settings));
        break;
      case ModelState::SNPwiseModelMeanImhof:
        engine_bay[type] = blockwise(new SNPwiseMeanImhof(settings));
        break;
      case ModelState::MultiModel:
        _LOG.error("initializing analysis") << "type MultiModel not allowed in function make_core_engine()" << endl; die();
        break; 
      default:
        _LOG.error("initializing analysis") << "unknown gene analysis model specified" << endl; die();
        break;
    }
  }
  return engine_bay[type];
}

GeneModelEngine* GeneAnalysis::make_multi_engine() {
  Aggregator::AggregateType multi_aggregate = (settings.gets("gene_multi_aggregate","") == "top") ? Aggregator::Top : Aggregator::Mean;
  MultiEngine* multi = new MultiEngine(settings, multi_aggregate);
  
  vector<string> parts = settings.getvs("gene_model_multi", true);
  for (int i = 0; i < parts.size(); i++) multi->add_engine(make_core_engine(get_model_type(parts[i])));        
  engine_bay[ModelState::MultiModel] = multi; return multi;
}

PartitionedEngine* GeneAnalysis::blockwise(CoreEngine* core, Aggregator::AggregateType aggregate_type) {
  return new BlockEngine(settings, core, aggregate_type);
}

PrinCompEngine* GeneAnalysis::princomp(PrinCompEngine* pc) {static int count = 0; pc->set_pcid(count++); return pc;}


void GeneAnalysis::print_model() {
  string model_name = settings.gets("gene_model"), prefix = "\tusing model: ";
  if (model_name == "multi") {
    _LOG << "\tusing multi-model containing:" << endl;    
    prefix = "\t\t";
  }
  
  vector<string> models = settings.getvs("gene_model_multi", true); models.push_back(model_name);
  if (Utils::contains(models, "pcreg")) _LOG << prefix << "principal components regression" << endl;
  if (Utils::contains(models, "varimax")) _LOG << prefix << "power-sum (P = " << settings.geti("gene_varimax",2) << ")" << endl;
  if (Utils::contains(models, "snpwise")) {
    string print = settings.gets("gene_snpwise");
    if (print == "sum" || print == "imhof") print = "mean";
    else if (print == "brown") print = "mean (Brown)";

    _LOG << prefix << "SNPwise-" << print; 
    if (settings.gets("gene_snpwise") == "top") {                 
      if (settings["gene_snpwise_top_count"]) _LOG << " (permutation, top " << settings.geti("gene_snpwise_top_count") << " SNPs in gene)";
      else if (settings["gene_snpwise_top_fraction"]) _LOG << " (permutation, top " << settings.getn("gene_snpwise_top_fraction")*100 << "% of SNPs in gene)";
    }
    _LOG << endl;
  } else {
    if (Utils::contains(models, "snpwise_mean")) _LOG << prefix << "SNPwise-mean" << endl;        
    if (Utils::contains(models, "snpwise_top1")) _LOG << prefix << "SNPwise-top" << endl;    
  }
  if (Utils::contains(models, "tdt")) _LOG << prefix << "transmission disequilibrium test" << endl;  

  if (settings["do_rare"]) {
    string auto_rare = settings["rare_autoburden"] ? " automatic" : "";
    string maf = (settings["rare_cutoff_maf"]) ? Utils::num_string_pad("MAF <= ", settings.getn("rare_cutoff_maf")) : "";
    string mac = (settings["rare_cutoff_mac"]) ? Utils::num_string_pad("MAC <= ", settings.geti("rare_cutoff_mac")) : "";

    _LOG << "\taggregating rare variants into" << auto_rare << " burden score per gene, for SNPs with ";
    if (maf == "" || mac == "") _LOG << maf << mac;
    else _LOG << "both " << maf << " and " << mac;

    string add = !settings["rare_freq_unweighted"] ? "normalized" : "";
    if (settings["has_rare_weights"]) add += string(!add.empty() ? ", " : "") + "weighted";
    int block_size = settings.geti("rare_multi_max_count"); //, 0);
    if (block_size > 0) add += string(!add.empty() ? ", " : "") + "max. " + Utils::num_string(block_size) + " variants per score";
    if (!add.empty()) _LOG << " (" << add << ")";
    _LOG << endl;
  }  

  if (settings["batch"]) _LOG << "\tanalysing batch " << settings.geti("batch_index") << " of " << settings.geti("batch") << ", containing " << Utils::plural(gene_info.data_used(false), "gene") << endl;
  else if (settings["batch_chr"]) _LOG << "\tanalysing chromosome " << Utils::chr_string(settings.geti("batch_chr"), 0) << ", containing " << Utils::plural(gene_info.data_used(false), "gene") << endl;
}

void GeneAnalysis::run() {
  _LOG << "Starting gene analysis... " << endl;
  print_model();

  if (compute_perm) {
    if (settings["adap_perm"]) perm_engine =  new PermutationEngine(settings.geti("min_perm"), settings.geti("max_perm"), settings.geti("adap_count"));
    else perm_engine = new PermutationEngine(settings.geti("min_perm"));
  }
  
  load_engine();
  run_analysis();
  
  _LOG << endl;                        
}


void GeneAnalysis::run_analysis() { try {
  int active_chromosome = set_chromosome();
  ModelInfoBlock model_info = engine->model_info(data);

  bool has_mac = settings.gets("data_type") == "binary_geno";
  bool store_rare = settings["do_rare"] && !settings["rare_autoburden"];
  bool has_rsq = settings["gene_rsq"];

  int tot_pval = model_info.size(); 

  gene_info.pval.set_width(tot_pval);
  gene_info.block_init(GeneData::GeneAnalysisBlock); 
  if (store_rare) gene_info.nrare.init();
  if (has_mac) gene_info.mac.init();
  if (compute_perm) {gene_info.perm_pval.init(); gene_info.perm_count.init();}
  if (has_rsq) gene_info.rsq.init();
  if (compute_corrs) gene_info.corrs.init();
  gene_info.set_model_info(model_info);

  if (settings["multi_pheno_permute"]) {
    DefaultStats* df_data = dynamic_cast<DefaultStats*>(data);
    if (df_data == 0) {_LOG.error("running analysis") << "invalid configuration for multi-phenotype permutation" << endl; die();}
    
    string outfile = settings.gets("out_prefix") + ".mpheno.minp";
    int no_perm =  settings.geti("multi_pheno_permute");
    _LOG << "\tgenerating " << no_perm << " permutations for multi-phenotype FWER correction, writing to file '" << outfile << "'" << endl;

    MultiPhenoFWER fwer(*df_data); 
    fwer.generate(outfile, no_perm);
  }

  Timer timer; int gene_count = 0;
  for (int i_gene = 0; i_gene < gene_info.data_total(); i_gene++) { try {
    if (!gene_info.active_id(i_gene)) continue;

    int curr_chromosome = gene_info.location.get(i_gene).chromosome;
    if (curr_chromosome != active_chromosome && !set_chromosome(curr_chromosome)) continue;

    gene_count++;

    analyse_gene(i_gene, tot_pval);
    gene_info.nsnps.set(i_gene, data->nsnps(), true);

    gene_info.nparam.set(i_gene, engine->get_nparam(), true);
    gene_info.nsamp.set(i_gene, data->get_sampsize(false), true);    
    if (store_rare) gene_info.nrare.set(i_gene, data->nrare(), true);
    if (has_mac) gene_info.mac.set(i_gene, data->get_mac(), true);

    gene_info.zstat.set(i_gene, get_zstat(), true);
    for (int i = 0; i < model_info.size(); i++) gene_info.pval.set(i_gene, i, engine->get_pval(i), true);

    if (compute_perm) {int nperm = perm_engine->nperm();
      gene_info.perm_count.set(i_gene, nperm, true);
      gene_info.perm_pval.set(i_gene, perm_engine->pval(), nperm);
    }

    if (has_rsq) {for (int i = 0; i < 2; i++) gene_info.rsq.set(i_gene, i, get_effsize(bool(i)), true);}
    
    if (compute_corrs) gene_info.corrs.assign(i_gene, corr_buffer, true);

    } catch (const GeneException& ge) {_LOG << "\tWARNING: analysis failed for gene " << gene_info.name.get(i_gene) << "; " << ge.get_msg() << endl; gene_info.drop_id(i_gene);}
    if (timer.check()) {
      _LOG.counter(true) << "\tprocessed genes: " << gene_count << " (" << setprecision(3) << round(float(1000*gene_count) / gene_info.data_used(false))/10  << "%)     "; _LOG.counter(false);
    }
  }
  _LOG.wipe_counter();

  if (gene_info.data_used(true) == 0) {_LOG.error("running gene analysis") << "analysis failed for all genes" << endl; die();} 
  
  GeneOutput output(gene_info, settings.gets("out_prefix"));
  if (settings["no_asym"]) output.set_option(OutputOptions::HideAsymptotic, true);
  if (model_info.state[ModelInfoBlock::HasMultiModel] && !settings["gene_multi_show_subpval"]) output.set_option(OutputOptions::HideMultiSubPval, true); 

  output.write_outfile();
  if (compute_corrs) output.write_rawfile();
} catch (const exception& e) {report_error(e, "running gene analysis");}}

void GeneAnalysis::analyse_gene(int geneid, int exp_pval) { try {
  data->load_snps(geneid);

  if (dump_data) data->write_gene(dump_pref);

  engine->run(data, compute_corrs, perm_engine);
  if (engine->ntest() != exp_pval) throw GeneException("analyzing gene", "unexpected number of p-values (expected " + Utils::num_string(exp_pval) + ", received " + Utils::num_string(engine->ntest()) + ")");

  if (compute_corrs) {
    CorrelationData* cd = engine->corr_data();
    cd->load_data(); cd->gene_id() = geneid;
    corr_engine.compute(corr_buffer, cd);
  }
} catch (const GeneException& de) {throw;} 
  catch (const exception& e) {check_mem_error(e, "analyzing gene"); throw GeneException("analyzing gene", "an internal error occurred when analyzing gene", e.what());}
}

float GeneAnalysis::get_zstat(int index) {
  static ConvertPvalToNorm zstat; static ConvertFromPermPval zstat_perm(StatType::Z);
  if (engine->perm_only()) {
    if (index == 0) {
      zstat_perm.set_param("N", perm_engine->nperm());
      return zstat_perm.convert(perm_engine->pval());    
    } else return -999;
  } else {
    double pval = engine->get_pval(index);
    return pval < 0 ? -999 : zstat.convert(pval);
  }
}

void GeneAnalysis::load_snpcount(int geneid, int expected) {
  EngineParam* sc = engine->get_param("snp_counts");
  if (sc && sc->size() == expected) {for (int i = 1; i < sc->size(); i++) gene_info.nsnps_part.set(geneid, i-1, sc->get(i), true);}
}

int GeneAnalysis::set_chromosome() {
  for (int i_gene = 0; i_gene < gene_info.data_total(); i_gene++) {
    if (!gene_info.active_id(i_gene)) continue;
    if (set_chromosome(gene_info.location.get(i_gene).chromosome)) return gene_info.location.get(i_gene).chromosome;
  }
  _LOG.error("running gene analysis") << "no analyzable genes found in data" << endl; die(); return 0;
}

bool GeneAnalysis::set_chromosome(int chr) {
  bool autosome = Utils::chr_type(chr) == 0;
  if (autosome || chr == Utils::chrX_code) {
    string chr_name = Utils::chr_string(chr,0);  
    try {
      bool changed = set_data(chr);
      if (changed) {
        if (!autosome && data->has_message(GeneStats::ModelChange)) _LOG << "\tanalyzing chromosome " << chr_name << ": " << data->get_message(GeneStats::ModelChange) << endl;      
        if (autosome && data->has_message(GeneStats::DroppedCovar)) {
          istringstream istr(data->get_message(GeneStats::DroppedCovar)); string cov;
          while (istr >> cov) _LOG << "\tWARNING: dropping covariate '" << cov << "'; variance is zero" << endl;        
        }      
        if (data->has_message(GeneStats::DroppedCovarChrX)) _LOG << "\tWARNING: dropping covariate 'sex' for chromosome " << chr_name << "; variance is zero" << endl;
        if (data->has_message(GeneStats::NoSubset)) _LOG << "\tWARNING: disabling 'big-data' mode, sample size is too small" << endl;
        else if (data->has_message(GeneStats::AutoSubset)) _LOG << "\tengaging automatic 'big-data' mode" << endl;
      }
      return true;
    } catch (const DataException& de) { 
      if (Utils::chr_type(chr) == 0) {_LOG.error("initializing phenotype/covariate data") << de.get_msg() << endl; die();} 
      else {
        _LOG << "\tWARNING: an error occurred when loading model for chromosome " << chr_name << ":\n\t\t" << de.get_msg() << endl;
        _LOG << "\tWARNING: skipping all chromosome " << chr_name << " genes" << endl;      
      }
    }
  } 
  gene_info.location.filter_chromosome(chr);
  return false;
}

bool GeneAnalysis::set_data(int chr) {
  data = statfactory.get_stats(chr);
  if (statfactory.changed) {
    if (compute_perm) perm_engine->set_data(data);
    if (compute_corrs) corr_engine.set_data(data);
    if (dump_data) data->write_vars(dump_pref);      
    return true;
  }
  return false;
}


/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef DATA_GENEDATA_H
#define DATA_GENEDATA_H

#include "utils_dataframe.h"
#include "geneutils.h"
#include "statutils.h"

class GeneData : public IndexedDataFrame {
public:
  static const short CoreBlock;
  static const short GeneAnalysisBlock;
  static const short GeneMetaBlock;  

  ModelInfoBlock model_info;

  GeneLocationVariable location;
  JaggedMultiVariable<NumericID<int> > snps;
  
  PositiveNumber<int> nsnps;
  BlockMultiVariable<NonNegativeNumber<int> > nsnps_part;  
  PositiveNumber<int> nparam;
  NonNegativeNumber<int> nrare;
  PositiveNumber<int> nsamp;
  BlockMultiVariable<NonNegativeNumber<int> > nsamp_part;    
  PositiveNumber<int> ndata;
  PositiveNumber<float> mac;
  
  ZscoreVariable<double> zstat;
  BlockMultiVariable<PvalueVariable<double> > pval;
  BlockMultiVariable<PosCorrelationVariable<double> > rsq;  
  
  PermPvalueVariable<double> perm_pval;
  PositiveNumber<int> perm_count;
  
  GeneCorrelationVariable corrs;

  ///new variables: add in GeneOutput write functions, GeneInput read/merge/meta functions
  GeneData() : IndexedDataFrame(true) {
    model_info = ModelInfoBlock::basic_block();
    id_index->set_param("drop_duplicates", true);
    id_index->set_param("tree_depth", 10);
    register_variable(location, "LOCATION", true, CoreBlock); 
    register_variable(snps, "SNP_IDS"); 
    register_variable(nsnps, "NSNPS", true, CoreBlock); 
    register_variable(nparam, "NPARAM", true, CoreBlock); 
    register_variable(mac, "MAC", false, CoreBlock);     
    register_variable(zstat, "ZSTAT", true, CoreBlock); 
    register_variable(corrs, "CORRS", false, CoreBlock);
    register_variable(nsnps_part, "NSNPS_PART", false, GeneAnalysisBlock); 
    register_variable(nrare, "NRARE", false, GeneAnalysisBlock);
    register_variable(nsamp, "N", true, GeneAnalysisBlock); 
    register_variable(pval, "PVAL", true, GeneAnalysisBlock);
    register_variable(rsq, "RSQ", false, GeneAnalysisBlock); rsq.set_width(2);
    register_variable(perm_pval, "PERMP",  false, GeneAnalysisBlock); 
    register_variable(perm_count, "NPERM", false, GeneAnalysisBlock);
    register_variable(ndata, "DATASETS", false, GeneMetaBlock); 
    register_variable(nsamp_part, "NSAMP_PART", false, GeneMetaBlock);
  }
  ~GeneData() {}
  using IndexedDataFrame::get_iid;  
  
  void build_index();  
  void block_init(short block) {if (block > CoreBlock) block_init(CoreBlock); IndexedDataFrame::block_init(block);}  

  using IndexedDataFrame::drop_id;
  using IndexedDataFrame::filter_ids;
  using IndexedDataFrame::make_consecutive;
  template<typename T> void data_sort(UniVariable<T>& var) {IndexedDataFrame::data_sort(var);}

  ModelInfoBlock& get_model_info() {return model_info;}
  void set_model_info(const ModelInfoBlock& info) {model_info = info.info.empty() ? ModelInfoBlock::basic_block() : info; model_info.detect_state();}
  
  bool create_zstat();
  bool create_pval();  
  bool upgrade_permp(bool overwrite=false);
};


#endif /** DATA_GENEDATA_H */

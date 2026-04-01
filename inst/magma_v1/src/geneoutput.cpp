/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "geneoutput.h"

using namespace OutputOptions;

void GeneOutput::set_defaults() {
  options.assign(_OPTIONCOUNT_, 0);
  options[HideAsymptotic] = false;
  options[HideMultiSubPval] = false;  
  options[ShowFitted] = false;  
  version_id = Settings::get_version();
  indent = true;
}

void GeneOutput::process_pval() {
  ModelInfoBlock info = gene_info.get_model_info();

  if (info.state[ModelInfoBlock::HasMultiModel] && options[HideMultiSubPval]) info.filter(ModelState::MultiModel);  
  if (info.state[ModelInfoBlock::HasPartitions] && !partition_rank.empty()) info.sort(partition_rank);
  else info.sort();

  pval_info = info.main_index("P");
  npart_info = info.partition_index(true);
}

void GeneOutput::write_outfile(string label, string suffix) {
  string filename = win_txt(prefix + "." + suffix + ".out");
  if (label.empty()) label = "gene analysis results";
  _LOG << (indent ? "\twriting" : "Writing") << " " << label << " to file " << filename << endl; 

  process_pval();
  vector<OutputColumn*> output; FormattedOutput fout(filename);
  output.push_back(gene_info.name.get_output("GENE"));
  output.push_back(gene_info.location.get_output());
  output.push_back(gene_info.nsnps.get_output("NSNPS"));

  output.push_back(gene_info.nsnps_part.merged_output(npart_info, "NSNPS_PART")); 
  output.push_back(gene_info.nsamp_part.merged_output(npart_info, "NSAMP_PART"));   
  output.push_back(gene_info.nrare.get_output("NRARE"));
  output.push_back(gene_info.nparam.get_output("NPARAM"));  
  output.push_back(gene_info.nsamp.get_output("N"));
  output.push_back(gene_info.ndata.get_output("DATASETS"));  
  if (!options[HideAsymptotic]) {
    typedef BlockMultiVariable<ZscoreVariable<double> > ZData;
    ZData* zstat = gene_info.typed_variable<ZData>("ZSTAT_OUTCOME");  
    if (zstat) output.push_back(zstat->get_output(0, "ZSTAT"));
    else {
      output.push_back(gene_info.zstat.get_output("ZSTAT"));
      ColumnBlock* block = gene_info.pval.get_output();
      if (block) output.push_back(block->add_columns(pval_info));
    }
  }
  
  if (options[ShowFitted]) {typedef BlockMultiVariable<ZscoreVariable<double> > FitData;
    FitData* fit = gene_info.typed_variable<FitData>("ZSTAT_FITTED_BASE");
    if (fit) {
      ColumnBlock* block = fit->get_output();
      if (block) {
        block->add_column(0, "ZFITTED_BASE");          
        block->add_column(1, "ZRESID_BASE");
        output.push_back(block);
      }
    }
  }
  
  output.push_back(gene_info.perm_pval.get_output("PERM" + pval_info[0].second));  
  output.push_back(gene_info.perm_count.get_output("NPERM"));

  ColumnBlock* rsq = gene_info.rsq.get_output(); 
  if (rsq) {rsq->add_column(0,"RSQ"); rsq->add_column(1,"RSQ_ADJ"); output.push_back(rsq);}

  for (vector<OutputColumn*>::iterator it = output.begin(); it != output.end();) {if (*it == 0 || (*it)->is_empty()) {delete *it; it = output.erase(it);} else ++it;}
  for (int i = 0; i < output.size(); i++) output[i]->set_header(fout);

  for (int i = 0; i < parameters.size(); i++) fout.print_param(parameters[i].first, parameters[i].second);
  if (!gene_info.location.human()) fout.print_param("nonhuman", 1);
  fout.print_header();

  while (!output[0]->is_empty()) {for (int i = 0; i < output.size(); i++) fout << *output[i];}
  for (int i = 0; i < output.size(); i++) delete output[i];
}

void GeneOutput::write_rawfile(string suffix) {
  string filename = prefix + "." + suffix + ".raw";
  _LOG << (indent ? "\twriting" : "Writing") << " intermediate output to file " << filename << endl;
  if (!gene_info.corrs.initialised()) {_LOG.error("writing .raw file") << "no gene correlation data available" << endl; die();}
  if (gene_info.unset_mask()) {_LOG.error("writing .raw file") << "cannot write gene correlations from masked GeneData object" << endl; die();}

  gene_info.make_consecutive();
  
  vector<OutputColumn*> output; DelimitedOutput fout(filename); string covar;
  output.push_back(gene_info.name.get_output());
  output.push_back(gene_info.location.get_output());
  output.push_back(gene_info.nsnps.get_output());
  output.push_back(gene_info.nparam.get_output());  
  if (gene_info.nsamp.initialised()) {Utils::str_append(covar, "NSAMP"); output.push_back(gene_info.nsamp.get_output());}
  if (gene_info.mac.initialised()) {Utils::str_append(covar, "MAC"); output.push_back(gene_info.mac.get_output());}
  output.push_back(gene_info.zstat.get_output()); 
  output.push_back(gene_info.corrs.get_output());

  for (vector<OutputColumn*>::iterator it = output.begin(); it != output.end();) {if (*it == 0 || (*it)->is_empty()) it = output.erase(it); else ++it;}
  if (version_id > 0) fout.print_param("VERSION", version_id);
  if (!covar.empty()) fout.print_param("COVAR", covar);
  if (!gene_info.location.human()) fout.print_param("NONHUMAN", 1);

  while (!output[0]->is_empty()) {
    for (int i = 0; i < output.size(); i++) fout << *output[i];
    fout << endl;
  }
  for (int i = 0; i < output.size(); i++) delete output[i];
}

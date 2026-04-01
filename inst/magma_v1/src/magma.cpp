/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "baseinput.h"
#include "geneinput.h"
#include "geneanalysis.h"
#include "setanalysis.h"
#include "annotation.h"
#include "parse.h"
#include "_global.h"

#include "magma.h"


int main(int argc, char* argv[]) {
try {
  Settings settings(argc, argv);
  _RNG.seed(settings.geti("rng_seed"));


  if (settings["do_analysis"] || settings["do_set_analysis"]) {GeneData* gd = 0;
    if (settings["do_analysis"]) {BaseInput* bi;
      if (settings.gets("do_analysis") == "plink") bi = new PlinkInput(settings);
      else {_LOG.error("configuration") << "No valid input data specified" << endl; die();}
      
      GeneAnalysis analysis(settings, bi);
      analysis.run();
      
      gd = bi->eject_data(); delete bi;
    } else if (settings["has_results_file"]) {
      GeneInputRawFile file(settings.gets("gene_results_file"), settings.gets("genelevel_correct", "all"), settings.getvs("genelevel_correct_list", true), settings.getn("gene_corr_rescale", 1));
      if (settings["gene_results_overwrite"]) file.overwrite_stat(settings.gets("gene_results_overwrite"));
      gd = file.eject_data();  
    } 
    if (settings["do_set_analysis"]) {
      if (!gd) {_LOG.error("configuration") << "No input for gene-level analysis found" << endl; die();}
      SetInput input(settings, *gd); gd = 0;
      SetVariables* vars = input.load_data();

      _LOG << endl;      
      SetAnalysis analysis(settings, vars); vars = 0;
      analysis.run();
    }
    delete gd; 
  } else if (settings["do_merge"]) {
    bool has_out = GeneInputMerge(settings, GeneInputCombine::OutMode).run(settings.gets("merge_prefix"));
    GeneInputMerge(settings, GeneInputCombine::RawMode).run(settings.gets("merge_prefix"), !has_out);
  } else if (settings["do_meta"]) GeneInputMeta(settings, settings["do_meta_raw"] ? GeneInputCombine::RawMode : GeneInputCombine::OutMode).run();
  else if (settings["do_annotate"]) Annotation annot(settings);  
  else {_LOG.error("configuration") << "no actionable flags have been set" << endl; die();}

  settings.print_goodbye();
  } catch (const BaseException& be) {
    _LOG.error() << endl << "FATAL ERROR - caught unprocessed exception of type " << be.get_exception() << "; error message:" << "\n\t(" << be.get_type() << ") " << be.get_msg() << endl; die();    
  } catch (const exception& e) {
    if (typeid(e) == typeid(bad_alloc)) _LOG.error() << endl << "FATAL ERROR - an error occurred when trying to allocate memory" << endl;
    else _LOG.error() << endl << "FATAL ERROR - caught unknown exception\n\tinternal error message: " << e.what() << endl;
    die();
  }
  return 0;
}



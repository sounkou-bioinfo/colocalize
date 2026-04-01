/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "setanalysis.h"
#include "output.h"
#include "exceptions.h"
#include "utils.h"
#include "mathutils.h"

void SetAnalysis::parse_models() { try {
  _LOG << "Parsing model specifications..." << endl;
  string error_tag = "running gene-level analysis";
  string output_prefix = settings.gets("out_prefix");

  if (settings["genelevel_sets_selfcontained"]) {
    if (variable_data->variables.no_sets().first <= 0) {_LOG.error(error_tag) << "self-contained gene-set analysis is requested, but no valid gene sets are available from input" << endl; die();}
    SelfContainedModel* model = new SelfContainedModel(*variable_data, output_prefix, "running self-contained analysis"); 
    model->validate();
    models.push_back(model);
  } 

  RegressionModel* model = new RegressionModel(*variable_data, output_prefix, "running gene-level regression", settings["genelevel_print_geneinfo"]); 
  RegressionConfig& config = model->get_config();
  config.set_direction(settings.geti("genelevel_direction_sets"), RegressionConfig::dt_GeneSet);
  config.set_direction(settings.geti("genelevel_direction_covar"), RegressionConfig::dt_GeneCovar);
  config.set_direction(settings.geti("genelevel_direction_interaction-ss"), RegressionConfig::dt_InteractionSS);
  config.set_direction(settings.geti("genelevel_direction_interaction-sc"), RegressionConfig::dt_InteractionSC);
  config.set_interaction(settings.geti("genelevel_interaction_sc-size"), settings.geti("genelevel_interaction_ss-size"), settings.getn("genelevel_interaction_ss-prop"));
  
  config.set_internal(settings.gets("genelevel_correct", "all"), settings.getvs("genelevel_correct_list", true)); 
  config.set_conditional(settings.getvs("genelevel_condition_residualize_list", true), RegressionConfig::ct_ConditionForce);
  config.set_conditional(settings.getvs("genelevel_condition_hide_list", true), RegressionConfig::ct_ConditionOnly);
  config.set_conditional(settings.getvs("genelevel_condition_list", true), RegressionConfig::ct_ConditionShow);
  config.set_conditional(settings.getvs("genelevel_condition_interaction_list", true), RegressionConfig::ct_ConditionInteraction);  
  config.set_residualize(settings["genelevel_residualize"]);
  config.set_collinear(settings["genelevel_allow_collinear"]);
  
  string mode = settings.gets("genelevel_analysis_mode", "uni");
  if (mode == "joint_list") config.set_models(settings.getvs("genelevel_joint_list", true), RegressionConfig::am_Joint);  
  else if (mode == "interaction_list") config.set_models(settings.getvs("genelevel_interaction_list", true), RegressionConfig::am_Interaction);  
  else {
    string use = settings.gets("genelevel_analysis_use", "all");
    if (use == "list") config.set_filter(settings.getvs("genelevel_analysis_use_list", true));
    else config.set_filter(use);
    if (mode == "joint_pairs") config.set_mode(RegressionConfig::am_JointPairs);
    else if (mode == "interaction_pairs") config.set_mode(RegressionConfig::am_InteractionPairs);
    else if (mode == "interaction_each") config.set_models(settings.getvs("genelevel_interaction_each", true), RegressionConfig::am_InteractionEach);
    else if (mode == "interaction_all") config.set_models(settings.getvs("genelevel_interaction_all", true), RegressionConfig::am_InteractionAll);    
  }

  model->validate();  
  models.push_back(model);
} catch (const SetException& se) {_LOG.error("running gene-level analysis") << se.get_msg() << endl; die();} }


void SetAnalysis::run() {
  parse_models();
  if (models.empty()) {_LOG.error("running gene-level analysis") << "no analysis models have been specified" << endl; die();}

  double alpha = settings.getn("genelevel_alpha", 0.05);
  int abbreviate = settings.geti("genelevel_abbreviate"); 
  bool abbr_file = settings["genelevel_abbreviate_file"];

  variable_data->init_corrs();
  for (int i = 0; i < models.size(); i++) {bool success = false;
    models[i]->set_alpha(alpha);    
    models[i]->set_abbreviate(abbreviate, abbr_file);
    models[i]->print(true);
    try {
      models[i]->run();
      success = true;
    } catch (const SetException& be) {
      _LOG.error(models[i]->error_tag) << be.get_msg() << endl;
    } catch (const BaseException& be) {      
     _LOG.error(models[i]->error_tag) << "caught unprocessed exception of type " << be.get_exception() << endl;
     _LOG << "\terror tag: " << be.get_type() << endl;
     _LOG << "\terror message: " << be.get_msg() << endl; 
    }  
    if (!success) {
      if (i+1 < models.size()) _LOG << endl << "Attempting to resume with next analysis model type in queue" << endl << endl;
      else die();
    }
  }
}


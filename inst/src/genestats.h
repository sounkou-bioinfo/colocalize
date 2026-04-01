/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GENESTATS_H
#define GENESTATS_H

#include "genestats_base.h"
#include "genestats_main.h"
#include "genestats_subset.h"

class GeneStatFactory {
  Settings& settings;
  BaseInput& baseinput;
  set<GeneStats::Notice> messages;

  enum StatType {st_Default, st_Subset, st_Transmission} stat_type;

  GeneStats* factory_bay[7];
  short active;

  bool sex_covar;

  void add_message(GeneStats* gs, GeneStats::Notice msg, bool once=true) {
    if (messages.find(msg) == messages.end()) {
      gs->set_message(msg);
      if (once) messages.insert(msg);
    }
  }

  GeneStats* make_stats(int chr) {
    GeneStats* out = 0;
    switch (stat_type) {
      case st_Subset: out = new SubsetStats(settings, baseinput); break;
      case st_Default: 
      default: out = new DefaultStats(settings, baseinput);
    }
    out->set_chromosome(chr);
    out->init();
    
    if (stat_type == st_Subset) {
      if (out->has_message(GeneStats::NoSubset)) {
        delete out; stat_type = st_Default;
        out = make_stats(chr);
        if (settings.gets("gene_model_bigdata") != "auto") add_message(out, GeneStats::NoSubset);
      } else if (settings.gets("gene_model_bigdata") == "auto") add_message(out, GeneStats::AutoSubset);
    }
    
    return out;
  }

public:
  GeneStatFactory(Settings& s, BaseInput* bi) : settings(s), baseinput(*bi), stat_type(st_Default), active(-1) {
    for (int i = 0; i < 7; i++) factory_bay[i] = 0;
    if (settings.gets("gene_model") == "tdt") stat_type = st_Transmission;
    else if (settings["gene_model_bigdata"]) {
      if (settings.gets("gene_model_bigdata") != "auto" || baseinput.samp_size() >= settings.geti("gene_model_bigdata_minsize_auto")) stat_type = st_Subset;
    }
  }

  ~GeneStatFactory() {
    for (int i = 0; i < 7; i++) delete factory_bay[i];
  }

  bool changed;

  GeneStats* get_stats(int chr, bool clear=true) {
    short chr_index = Utils::chr_type(chr);

    try {
      changed = (chr_index != active);
      if (changed) {
        if (active >= 0 && clear) factory_bay[active]->clear();
        active = chr_index;
        if (factory_bay[active] == 0) factory_bay[active] = make_stats(chr);
      }

      return factory_bay[active];
    } catch (const DataException& de) {
      if (active >= 0) {
        delete factory_bay[active];
        factory_bay[active] = 0;
        active = -1;
      }
      throw; 
    }
  }
};


#endif /** GENESTATS_H */
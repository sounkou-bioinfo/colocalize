/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "data_genedata.h"

const short GeneData::CoreBlock = 2;
const short GeneData::GeneAnalysisBlock = 3;
const short GeneData::GeneMetaBlock = 4;

void GeneData::build_index() {
  IndexedDataFrame::build_index();
  if (location.initialised()) data_sort(location); 
  else shrink_storage();
}

bool GeneData::create_zstat() {
  if (pval.initialised() && pval.get_width() > 0) {ConvertPvalToNorm to_stat;
    if (zstat.initialised()) zstat.assign_missing(); else zstat.init();
    for (int iid = 0; iid < index_max; iid++) {
      if (local_map[iid] > 0) { try {
        double p = pval.get(iid, 0);
        if (p >= 0 && p <= 1) zstat.set(iid, to_stat.convert(p), true);
      } catch (const exception& e) {}}
    }
    return true;
  } else return false;
}

bool GeneData::create_pval() {
  if ((!pval.initialised() || pval.get_width() == 0) && zstat.initialised()) {
    ConvertNormToPval to_pval;
    pval.set_width(1); pval.init(); pval.assign_missing();
    double miss_code = zstat.get_missing().second;
    for (int iid = 0; iid < index_max; iid++) {
      if (local_map[iid] > 0) { try {
        double z = zstat.get(iid);
        if (z != miss_code) pval.set(iid, 0, to_pval.convert(z), true);
      } catch (const exception& e) {}}
    }
    return true;
  } else return false;
}


bool GeneData::upgrade_permp(bool overwrite) {
  if (perm_pval.initialised() && (overwrite || !pval.initialised())) {
    if (pval.initialised()) {pval.set_width(1); pval.assign_missing();}
    else pval.init(1);

    for (int iid = 0; iid < index_max; iid++) {if (local_map[iid] > 0) pval.set(iid, 0, perm_pval.get(iid), true);}
    return true;
  } else return false;
}

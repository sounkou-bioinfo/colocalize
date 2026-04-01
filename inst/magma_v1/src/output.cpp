/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "output.h"
#include <algorithm>

SupplementaryLog& SupplementaryLog::set_block(const string& title, bool set_hold) {
  pending = false; log_buff.str(""); hold = set_hold;
  if (!empty) log_target() << endl;
  if (title == "") log_target() << "##### ##### ##### #####" << endl;
  else log_target() << "##### #####  " << title << "  ##### #####" << endl;
  return *this;  
}

void OutputColumn::set_header(FormattedOutput& fout) {fout.add_field(width > 0 ? width : 1, name);}

void OutputBuffer::process_buffer(vector<string>& target, bool by_line) {string str;
  if (!target.empty()) target.clear();
  out_target.clear(); out_target.seekg(0); 
  if (by_line) {while (out_target.good()) {getline(out_target, str); target.push_back(str);}}
  else {while (out_target >> str) target.push_back(str);}
  out_target.clear(); out_target.str("");
}

int FormattedOutput::add_field(long long width, const string& name, bool comp_width) {
  if (comp_width) width = width > 1 ? floor(log10(width))+1 : 1;
 
  field_width.push_back(max(width, (long long) name.size()) + ((no_fields) ? field_pad : 0));
  field_truncate.push_back(0);
  field_name.push_back(name);
  field_precision.push_back(0);
  field_align_right.push_back(true);
  return no_fields++;
}

int FormattedOutput::add_float_field(int precision, const string& name) {
  int id = add_field(precision+6, name);
  field_precision[id] = precision;
  return id;
}

int FormattedOutput::add_field(const vector<string>& str_vec, const string& name) {
  size_t len = 0;
  for (vector<string>::const_iterator iter = str_vec.begin(); iter != str_vec.end(); ++iter) len = max(len, iter->size());
  return add_field(len, name);
}

void FormattedOutput::clear_fields() {
  field_width.clear();
  field_truncate.clear();
  field_name.clear();
  field_precision.clear();
  field_align_right.clear();
  no_fields = 0;
}

void FormattedOutput::set_align_right(bool align_right, int fid) {
  if (fid < 0) fid = no_fields-1;
  field_align_right[fid] = align_right;
}

bool FormattedOutput::set_max_width(int width, int fid) {
  if (fid < 0) fid = no_fields-1;
  if (width < field_width[fid]) {
    field_width[fid] = width;
    field_truncate[fid] = 1;
    return true;
  }
  return false;
}


void FormattedOutput::print_header() {
  for (int i = 0; i < no_fields; i++) *this << field_name[i];
  out_target.precision(1);
}



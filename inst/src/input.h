/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef INPUT_H
#define INPUT_H

#include "utils.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <map>

#define _FILE_OFFSET_BITS 64

using namespace std;

class TextInput {
  string filename;
  ifstream in_file;
  istringstream istr;
  istringstream convert;
  
  map<string,string> params;
  
  int ncol;
  int col_offset;
  deque<string> col_values;
  deque<string> col_names;
  map<string, int> col_index;
  map<string, int> col_special;

  vector<string*> col_remap;
  vector<string*> col_renames;  
  vector<short> col_keep;
  bool use_all;
  
  string err_prefix;
  int err_skip;
  int err_count;
  bool multi_error;  
  
  pair<string,string> parse_param(const string& line);
  
public:
  bool has_header;
  long line_no;
  string curr_line;
  string value;

  TextInput(const string& filename, int nc=1) : filename(filename), ncol(nc), col_offset(0), use_all(true), err_prefix("reading file"), err_skip(1), err_count(0), multi_error(false), has_header(false), line_no(0), curr_line(200, '\0') {
    in_file.open(filename.c_str(), ifstream::in);
    col_values.assign(ncol, string(10, '\0'));
  }

  ~TextInput() {if (in_file.is_open()) in_file.close();}

  void set_error(const string& pref, int skip=1) {
    err_prefix = pref;
    err_skip = skip;
  }
  
  void rewind() {
    in_file.seekg(0);
    in_file.clear();
    line_no = 0;
    skip_empty();
  }
  
  bool read_value() {return bool(istr >> value);}
  bool read_line();
  bool process_line(bool skip=true);
  int process_tail();
  string& operator[](const long& index) {return col_values[index];}
  string& operator()(const long& index) {return *(col_remap[index]);}
  void skip_empty();

  unsigned long long file_size();

  bool has_param(const string& name, bool check_value=true) {string value_str;
    if (get_param(name, value_str)) {
      if (check_value && !value_str.empty()) {long value;
        if (Utils::convert_num(value_str, value) && !value) return false;   
        else return true;    
      } else return true;
    } else return false;
  }
  bool get_param(const string& name, string& target) {
    string name_upper = Utils::uppercase(name);
    map<string, string>::const_iterator found = params.find(name_upper);
    if (found == params.end()) {target = ""; return false;} 
    else {target = found->second; return true;}
  }
  template<typename T>
  bool get_param(const string& name, T& target) {string value_str;
    if (get_param(name, value_str)) {T value; 
      bool convert = Utils::convert_num(value_str, value);    
      if (convert) target = value;
      return convert;
    } else return false;
  }
  vector<string> get_param_values(const string& name) {
    vector<string> out; string values;
    if (get_param(name, values)) out = Utils::tokenize(values); 
    return out;
  }


  void set_size(int size) {ncol = size; col_values.assign(ncol, string(10, '\0'));}
  int detect_size();
  int read_header(int skip=0, bool run_check=true);
  void drop_header() {
    has_header=false; col_names.clear();
    for (int i = 0; i < ncol; i++) col_names.push_back("variable " + Utils::num_string(i - col_offset + 1));
    rewind();
  }
  bool has_var(const string& name) {return has_header && col_index.find(Utils::lowercase(name)) != col_index.end();}
  int get_var(const string& name) {
    if (!has_var(Utils::lowercase(name))) error("variable '" + name + "' not found");
    return col_index[Utils::lowercase(name)];
  }

  string* set_var(const string& id, const string& use, bool allow_dup=false);
  int var_index(const string& name, bool remap=false);
  int set_subset(const vector<string>& include, const vector<string>& exclude);
  vector<short> parse_used(const vector<string>& use);
  int set_map();
  deque<string>& get_names() {return col_names;}
  void load_names(vector<string>& name_vec) {for (int i = 0; i < ncol; i++) {if (col_keep[i] > 0) name_vec.push_back(col_names[i]);}}
  vector<pair<int,string> > load_names() {vector<pair<int,string> > index;
    for (int i = 0; i < ncol; i++) {if (col_keep[i] > 0) index.push_back(pair<int,string>(i, col_names[i]));}
    return index;
  }  
  void load_indices(vector<short>& index_vec, vector<string> use) {
    vector<short> subset = parse_used(use);
    for (int i = 0; i < ncol; i++) {
      if (col_keep[i] > 0) index_vec.push_back(subset[i]);
    }
  }

  string print_names() {
    string out;
    for (int i = 0; i < ncol; i++) {
      if (col_keep[i] > 0) {
        if (!out.empty()) out.append(", ");
        if (!col_names.empty()) out.append(col_names[i]);
        else out.append("variable " + Utils::num_string(i-col_offset+1));
      }
    }
    return out;
  }

  template<typename T>
  bool read_num(T& target, const string& input) {
    convert.clear(); convert.str(input); convert >> target;
    return convert.eof() && !convert.fail();
  }

  template<typename T>
  bool read_num(T& target, int index) {
    convert.clear(); convert.str(col_values[index]); convert >> target;
    return convert.eof() && !convert.fail();
  }

  const string& name(int index, bool remap=false) {
    if (remap) return *(col_renames[index]);
    else return col_names[index];
  }

  const string& name(const string& id) {
    return col_names[col_special[id]];
  }
  
  void error(const string& msg);
  void line_error(const string& msg, const string& suff="") {line_error(msg, -1, suff);}
  void line_error(const string& msg, int column, const string& suff="");
  void line_error(const string& msg, const vector<string>& suff) {line_error(msg, -1, suff);}
  void line_error(const string& msg, int column, const vector<string>& suff);
  void count_error(int observed, int expected);
  
  void process_exit(bool force_exit=false);
  void hold_errors(bool hold);
};









#endif /** INPUT_H */
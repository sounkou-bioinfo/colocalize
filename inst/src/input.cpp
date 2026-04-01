/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "input.h"
#include "output.h"
#include "_global.h"           

bool TextInput::read_line() {
  if (getline(in_file, curr_line)) {
    line_no++;
    istr.clear();
    istr.str(curr_line);
    return true;
  }
  return false;
}

bool TextInput::process_line(bool skip) {
  if (read_line()) {
    int val;
    for (val = 0; val < ncol; val++) {
      if (!(istr >> col_values[val])) break;
    }

    if (val == ncol) return true;
    if (val > 0) count_error(val, ncol);

    if (skip) return process_line(true);
    else {
      char ch;
      while (in_file.get(ch)) {
        if (!isspace(ch)) {
          _LOG.error(err_prefix, err_skip) << "line " << line_no << " is empty" << endl;
          process_exit();
        }
      }
    }
  }
  return false;
}

int TextInput::process_tail() {
  int count = ncol;
  while (read_value()) {
    if (count >= col_values.size()) col_values.push_back(value);
    else col_values[count].swap(value);
    count++;
  }  
  return count;
}

unsigned long long TextInput::file_size() {
  ifstream file(filename.c_str(), ios::in|ios::binary|ios::ate);
  return file.tellg();
}

int TextInput::detect_size() {  
  rewind(); skip_empty();
  if (!read_line()) error("file is empty");

  ncol = 0;
  while (istr >> value) ncol++;
  col_values.assign(ncol, string(10, '\0'));  
  
  rewind();
  return ncol; 
}

int TextInput::read_header(int skip, bool run_check) {
  col_offset = skip;
  skip_empty();
  read_line();
  
  has_header = true; double dump;
  while (istr >> value) {
    col_names.push_back(value);
    if (has_header && run_check && col_names.size() > col_offset && read_num(dump, value)) has_header = false;
  }
  ncol = col_names.size();
  col_values.assign(ncol, string(10, '\0'));

  if (ncol <= col_offset) error("file contains no variables");
  
  if (has_header) {
    for (long i = 0; i < ncol; i++) {
      if (has_var(col_names[i])) error("duplicate variable name '" + col_names[i] + "'");
      col_index[Utils::lowercase(col_names[i])] = i;
    }
  } else {
    drop_header();
  }

  col_keep.assign(ncol, 1);
  for (int i = 0; i < col_offset; i++) col_keep[i] = 0;
  return ncol - col_offset;
}

string* TextInput::set_var(const string& id, const string& use, bool allow_dup) {
  if (use_all) {
    col_keep.assign(ncol, 0);
    use_all = false;
  }

  long index = var_index(use);
  if (!allow_dup) {
    for (map<string,int>::iterator i = col_special.begin(); i != col_special.end(); ++i) {
      if (i->second == index) error("variable specified for id '" + id + "' already assigned (variable = " + col_names[index] + ")"); 
    }
  }
  col_keep[index] = -1;
  col_special[id] = index;
  return &(col_values[index]);
}

int TextInput::var_index(const string& name, bool remap) {
  long index;
  if (read_num(index, name)) {
    if (index < 1 || index > ncol-col_offset) error("invalid variable index (" + name + ")");
    index = index+col_offset-1;
  } else {
    if (has_header) {
      index = get_var(name);
    } else error("cannot select variables by name when no header is present");
  }
  if (remap) {
    string rename = col_names[index]; index = -1;
    for (int i = 0; i < col_renames.size(); i++) {
      if (rename == *(col_renames[i])) {
        index = i;
        break;
      }
    }
    if (index < 0) error("cannot find local index for variable " + rename);
  }
  
  return index;
}

///overrides manually set variables with set_var
int TextInput::set_subset(const vector<string>& include, const vector<string>& exclude) {
  col_keep.assign(ncol, 0);
  if (exclude.empty()) col_keep = parse_used(include); ///include all selected
  else {
    vector<short> drop = parse_used(exclude);
    if (!include.empty()) col_keep = parse_used(include); ///include all not selected (drop), plus all selected already (special case, exclude modifier + internally added)
    for (int i = col_offset; i < ncol; i++) col_keep[i] = col_keep[i] || !drop[i];
  }

  long used_col = col_offset;
  for (int i = col_offset; i < ncol; i++) used_col += col_keep[i];

  use_all = used_col == ncol;
  return used_col;
}


vector<short> TextInput::parse_used(const vector<string>& use) {
  long index, tot_var = ncol-col_offset;
  vector<short> keep(ncol, 0);
  if (use.empty()) {
    for (int i = col_offset; i < ncol; i++) keep[i] = 1;
    return keep;
  }
  
  for (int i = 0; i < use.size(); i++) {
    if (read_num(index, use[i])) {
      if (index < 1 || index > tot_var) error("invalid variable index (" + use[i] + ")");
      keep[index+col_offset-1] = 1;
    } else {
      if (has_var(use[i])) {keep[get_var(use[i])] = 1;}
      else {
        pair<string,string> bits = Utils::split_str(use[i], '-');
        if (bits.first.empty() && bits.second.empty()) {
          if (has_header) error("variable '" + use[i] + "' not found");
          else error("cannot select variables by name when no header is present");
        }

        long col1, col2;
        bool numeric = read_num(col1, bits.first);
        numeric = read_num(col2, bits.second) && numeric;
        if (numeric) {
          col1 += col_offset - 1;
          col2 += col_offset - 1;
        } else {
          if (!has_header) error("cannot select variables by name when no header is present");
          if (has_var(bits.first) && has_var(bits.second)) {
            col1 = get_var(bits.first);
            col2 = get_var(bits.second);
          } else {
            if (!has_var(bits.first)) error("variable '" + bits.first + "' not found");
            if (!has_var(bits.second)) error("variable '" + bits.second + "' not found");
          }
        }

        if (col1 > col2) {index = col1; col1 = col2; col2 = index;}
        if (numeric && (col1 < col_offset || col2 >= ncol)) error("invalid variable range (" + use[i] + ")");
        while (col1 <= col2) keep[col1++] = 1;
      }
    }
  }
  return keep;
}


///NOTE: map is (potentially) invalidated by changes in ncol
int TextInput::set_map() {
  if (col_values.size() < ncol) col_values.assign(ncol, string(10, '\0'));
  col_remap.clear(); col_renames.clear();
  for (long i = 0; i < ncol; i++) {
    if (col_keep[i] > 0) {
      col_remap.push_back(&(col_values[i]));
      col_renames.push_back(&(col_names[i]));      
    }
  }
  return col_remap.size();
}


void TextInput::skip_empty() {
  if (!in_file.good()) return;
  while (true) {
    while (!in_file.eof() && isspace(in_file.peek())) {
      if (in_file.get() == '\n') line_no++;
    }
    if (in_file.eof()) error("unexpected end of file");
    if (in_file.peek() != '#') return;
    read_line();
    pair<string,string> par = parse_param(curr_line);
    if (!par.first.empty()) params[par.first] = par.second;
  }
}

pair<string,string> TextInput::parse_param(const string& line) {
  pair<string,string> out; int pos = 0; 
  istringstream istr(line); string curr;

  if (istr >> curr && curr[0] == '#') {
    if (curr.size() == 1) {if (!(istr >> curr)) goto END;}
    else curr = curr.substr(1);

    if (!curr.empty()) {
      if ((pos = curr.find('=')) == string::npos) {     
        out.first = curr;
        if (!(istr >> curr)) goto END;
      } else {
        out.first = curr.substr(0, pos);
        curr = curr.substr(pos);
      }
     
      if (!curr.empty() && curr[0] == '=') {
        if (curr.size() == 1) {if (!(istr >> curr)) goto END;}
        else curr = curr.substr(1);   
        if (!curr.empty()) {
          out.second = curr;
          while (istr >> curr) out.second += " " + curr; 
        }
      }
    }
  }
  END: out.first = Utils::uppercase(out.first);
  return out;
}

void TextInput::error(const string& msg) {
  _LOG.error(err_prefix, err_skip) << msg << endl;
  process_exit();
}

void TextInput::line_error(const string& msg, int column, const string& suff) {
  _LOG.error(err_prefix, err_skip) << msg << " on line " << line_no;             
  if (column >= 0) _LOG << ", column " << column;
  if (suff.size() > 0) _LOG << " (" << suff << ")";
  _LOG.error() << endl << "\tline: " << curr_line << endl;
  process_exit();
}

void TextInput::line_error(const string& msg, int column, const vector<string>& suff) {
  _LOG.error(err_prefix, err_skip) << msg << " on line " << line_no;             
  if (column >= 0) _LOG << ", column " << column;
  if (suff.size() > 0) {
    _LOG << ":" << endl;
    for (int i = 0; i < suff.size(); i++) _LOG <<"\terror " << (i+1) << ": " << suff[i] << endl;
  }
  _LOG.error() << endl << "\tline: " << curr_line << endl;
  process_exit();
}
    
void TextInput::count_error(int observed, int expected) {
  _LOG.error(err_prefix, err_skip) << "too few values on line " << line_no << " (found " << observed << ", expecting at least " << expected << "):" << endl;
  _LOG.error() << "\tline: " << curr_line << endl;
  process_exit();
}

void TextInput::process_exit(bool force_exit) {
  if (multi_error && !force_exit) err_count++;
  else die();
}

void TextInput::hold_errors(bool hold) {
  if (err_count > 0) process_exit(true);
  multi_error = hold;
}


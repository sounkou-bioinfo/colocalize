/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>       

#include "utils.h"

using namespace std;

class Log {
  ofstream log_file;
  ostringstream log_buff;
  ostream* log_curr;
  string filename;
  
  bool batch_mode;
  bool to_screen;
  bool to_file;
  bool file_hold;
  bool add_skip;

public:
  Log(bool batch=false) : batch_mode(batch), to_screen(true), to_file(false), file_hold(false), add_skip(false) {}
  ~Log() {if (log_file.is_open()) log_file.close();}

  template<typename T>
  Log& operator<< (const T& data) {
    if (add_skip) {add_skip = false; *this << endl;}
    if (to_screen) cout << data;
    if (to_file) *log_curr << data;
    return *this;
  }

  Log& operator<< (ostream&(*manip)(ostream&)) { ///for endl
    if (to_screen) cout << manip;
    if (to_file) *log_curr << manip;
    reset();
    return *this;
  }
 
  void set_logfile(const string& fname, bool hold=false, bool append=false) {
    if (hold) {
      log_curr = &log_buff;
    } else {
      if (log_file.is_open() && fname != filename) log_file.close();
      if (!log_file.is_open()) {
        if (append) log_file.open(fname.c_str(), ofstream::out | ios::app);
        else log_file.open(fname.c_str(), ofstream::out);
      }
      log_curr = &log_file;
      if (file_hold) {
        log_file << log_buff.str();
        log_buff.str("");
      }
    }
    filename = fname;
    file_hold = hold;
    reset();
  }
  
  void change_logfile(const string& fname) {
    if (fname == filename) return;
    if (log_file.is_open()) {
      *this << "Changing log file from " << filename << " to " << fname << endl;
      ifstream src(filename.c_str(), ios::binary);
      ofstream dest(fname.c_str(), ios::trunc | ios::binary);
      if (!src.good() || !dest.good()) {*this << "WARNING: unable to change log file" << endl; return;}

      *this << endl;
      dest << src.rdbuf();
    }

    string old_log = filename;
    set_logfile(fname, false, true);
    remove(old_log.c_str());
  }
  
  void unhold() {if (file_hold) {set_logfile(filename, false);} reset();}
  void set_skip(bool on) {add_skip = on;}
  
  void reset() {
    to_screen = true;
    to_file = log_file.is_open() || file_hold;
  }
  
  Log& flush() {
    if (to_screen) cout.flush();
    if (to_file) log_file.flush();
    return *this;
  }
  
  ///Modifiers
  Log& operator<< (Log& log) {return log;}

  Log& counter(bool on=true) {
    if (on) {reset(); to_file=false;}
    else {*this << "\r"; cout.flush(); reset();}
    return *this;
  }
  void wipe_counter() {
    this->counter(true) << "                                                                                                                   ";
    this->counter(false);
  }
  Log& file(bool toggle=false) {if (!toggle) reset(); to_screen=false; return *this;}
  Log& screen(bool toggle=false) {if (!toggle) reset(); to_file=false; return *this;}
  Log& batch(bool toggle=false) {if (!toggle) reset(); to_screen=to_screen && batch_mode; to_file=to_file && batch_mode; return *this;}
  Log& no_batch(bool toggle=false) {if (!toggle) reset(); to_screen=to_screen && !batch_mode; to_file=to_file && !batch_mode; return *this;}
  Log& error() {unhold(); return *this;}
  Log& error(const string& msg, short newl=1) {
    unhold();
    for (short i = 0; i < newl; i++) *this << "\n";
    *this << "ERROR - " << msg << ": ";
    return *this;
  }
  Log& all() {reset(); return *this;}

  ///Formatted output
  Log& time(long sec) {
    flush();
    *this << setw(2) << setfill('0') << sec / 3600 << ":";
    *this << setw(2) << setfill('0') << sec / 60 % 60 << ":";
    *this << setw(2) << setfill('0') << sec % 60;
    return *this;
  }
};

class SupplementaryLog {
  ofstream log_file;
  ostringstream log_buff;
  string filename;

  bool empty;
  bool hold;
  bool pending;

  ostream& log_target() {
    if (hold) {
      pending = true;    
      return log_buff;
    } else {
      if (!log_file.is_open()) log_file.open(filename.c_str(), ofstream::out);
      if (pending) {
        log_file << log_buff.str(); log_buff.str("");
        pending = false;
      }
      empty = false;
      return log_file;
    }
  }
  
  void purge_logfile(const string& fname) {if (fname != filename && Utils::is_file(fname)) remove(fname.c_str());}

public:
  SupplementaryLog() : empty(true), hold(false), pending(false) {}    
  ~SupplementaryLog() {if (log_file.is_open()) log_file.close();}

  template<typename T>
  SupplementaryLog& operator<< (const T& data) {
    log_target() << data;
    return *this;
  }

  SupplementaryLog& operator<< (ostream&(*manip)(ostream&)) { ///for endl
    log_target() << manip;
    return *this;
  }

  SupplementaryLog& set_block(const string& title, bool set_hold=false); 
  void hold_off() {hold = false;}
  
  void set_logfile(const string& fname) {
    purge_logfile(fname);
    if (log_file.is_open() && fname != filename) change_logfile(fname);
    else filename = fname;
  }
  
  void change_logfile(const string& fname) {
    if (fname == filename) return;
    purge_logfile(fname);
    
    if (log_file.is_open()) {
      ifstream src(filename.c_str(), ios::binary);
      ofstream dest(fname.c_str(), ios::trunc | ios::binary);
      if (src.good() && dest.good()) {
        dest << src.rdbuf();
        remove(filename.c_str());
      }
    }
    filename = fname; 
  }
};


class OutputBuffer; class DelimitedOutput; class FormattedOutput;
class OutputStream {
  OutputBuffer* obuff; DelimitedOutput* dout; FormattedOutput* fout;
public:
  OutputStream() {unset();}
  void set(OutputBuffer* o) {unset(); obuff = o;}
  void set(DelimitedOutput* o) {unset(); dout = o;}
  void set(FormattedOutput* o) {unset(); fout = o;}

  void unset() {obuff = 0; dout = 0; fout = 0;}
  
  template<typename T>
  OutputStream& operator<< (const T& data);
};

class OutputColumn {
protected:
  string name;
  int width;
public:
  OutputColumn(const string& name="", int width=0) : name(name), width(width) {}
  virtual ~OutputColumn() {}

  virtual void print(OutputStream& out) = 0;
    
  virtual void set_header(FormattedOutput& fout); 
  virtual bool is_empty() = 0; 
  virtual void rewind() = 0;
  
  OutputColumn* set_name(const string& value) {name = value; return this;}
  OutputColumn* set_width(int value) {width = (value > 0 ? value : 0); return this;}
};


template<typename OT>
class Output {
protected:
  OT out_target;
  OutputStream out_stream;
  
  template<typename T> 
  static T get_max(const vector<T>& vec) {
    T curr_max = !vec.empty() ? vec[0] : 0;
    for (long i = 1; i < vec.size(); i++) curr_max = max(curr_max, vec[i]);
    return curr_max;
  }
  
  Output() {}
public:
  virtual ~Output() {}
 
  template<typename T>
  Output& operator<< (const T& data) {out_target << data; return *this;}
  Output& operator<< (ostream&(*manip)(ostream&) ) {out_target << manip; return *this;}
  Output& operator<< (OutputColumn& col) {col.print(out_stream); return *this;}  
  
  void print_comment(const string& text) {out_target <<  "# " << text << endl;}
  void print_param(const string& name) {out_target << "# " << name << endl;}
  template<typename T>
  void print_param(const string& name, const T& value, const string& suffix="") {out_target << "# " << name << " = " << value << suffix << endl;}  
};


class OutputBuffer : public Output<stringstream> {
  string separator;
  int curr_field;
  
public:
  OutputBuffer(string sep=" ") : separator(sep), curr_field(0) {out_stream.set(this);}

  template<typename T>
  OutputBuffer& operator<< (const T& data) {if (curr_field++) out_target << separator; out_target << data; return *this;}
  OutputBuffer& operator<< (ostream&(*manip)(ostream&) ) {curr_field = 0; out_target << manip; return *this;}  
  OutputBuffer& operator<< (OutputColumn& col) {col.print(out_stream); return *this;}    
  
  void process_buffer(vector<string>& target, bool by_line=false);
};


class FileOutput : public Output<ofstream> {
public:
  FileOutput() {}
  FileOutput(const string& filename, bool append=false) {open(filename, append);} 
  virtual ~FileOutput() {if (out_target.is_open()) out_target.close();}
  
  void open(const string& filename, bool append=false) {
    if (!out_target.is_open()) {
      if (append) out_target.open(filename.c_str(), ofstream::out | ios::app);
      else out_target.open(filename.c_str(), ofstream::out);
    }  
  }
  bool is_open() {return out_target.is_open();}
};

class DelimitedOutput : public FileOutput {
protected:
  char separator;
  int curr_field;
  
public:
  DelimitedOutput(char sep=' ') : separator(sep), curr_field(0) {out_stream.set(this);}  
  DelimitedOutput(const string& filename, char sep=' ', bool append=false) : FileOutput(filename, append), separator(sep), curr_field(0) {out_stream.set(this);}
  virtual ~DelimitedOutput(){}

  template<typename T>
  DelimitedOutput& operator<< (const T& data) {if (curr_field++) out_target << separator; out_target << data; return *this;}
  DelimitedOutput& operator<< (ostream&(*manip)(ostream&) ) {curr_field = 0; out_target << manip; return *this;}
  DelimitedOutput& operator<< (OutputColumn& col) {col.print(out_stream); return *this;}    
};

class FormattedOutput : public DelimitedOutput {
  int field_pad;
  int no_fields;

  vector<int> field_width;
  vector<short> field_truncate;
  vector<string> field_name;
  vector<int> field_precision;
  vector<short> field_align_right;

  template<typename T>
  FormattedOutput& process (const T& data) {
    if (field_precision[curr_field]) out_target << setprecision(field_precision[curr_field]);
    out_target << (field_align_right[curr_field] ? right : left);
    if (field_align_right[curr_field] || curr_field < (no_fields-1)) out_target << setw(field_width[curr_field]);
    out_target << data;
    if (++curr_field == no_fields) {
      out_target << '\n';
      curr_field = 0;
    } else out_target << separator;
    
    return *this;
  }

public:
  FormattedOutput() : DelimitedOutput(' '), field_pad(1), no_fields(0) {out_stream.set(this);}
  FormattedOutput(const string& filename, bool append=false) : DelimitedOutput(filename, ' ', append), field_pad(1), no_fields(0) {out_stream.set(this);}
  ~FormattedOutput(){}

  FormattedOutput& operator<< (const string& data) {
    if (field_truncate[curr_field] && data.size() > field_width[curr_field]) {
      string trunc = data.substr(0, field_width[curr_field] - 3).append("...");
      return process(trunc);
    }
    return process(data);
  }

  template<typename T>
  FormattedOutput& operator<< (const T& data) {return process(data);}
  FormattedOutput& operator<< (OutputColumn& col) {col.print(out_stream); return *this;}  

  int size() {return no_fields;}

  int add_field(long long width, const string& name="", bool comp_width=false);
  int add_float_field(int precision, const string& name="");

  int add_field(const vector<string>& str_vec, const string& name);
  int add_field(const vector<unsigned long>& num_vec, const string& name) {return add_field(get_max(num_vec), name, true);}
  template<typename T> ///Intended for numeric types
  int add_field(const vector<T>& num_vec, const string& name) {
    T abs_max = 0;
    for (long i = 0; i < num_vec.size(); i++) abs_max = max(abs_max, T(abs(num_vec[i])));
    return add_field(abs_max, name, true);
  }
  void clear_fields();

  void rename_field(int fid, const string& name) {field_name[fid] = name;}
  void set_align_right(bool align_right, int fid=-1);
  bool set_max_width(int width, int fid=-1);
  
  void print_header();
  void skip_line(int lines) {for (int i = 0; i < lines; i++) out_target << endl;}
  void fill_line(const string& val) {while (curr_field) (*this) << val;}
  void fill_to_col(const string& val, int col_id) {for (int i = curr_field; i < col_id; i++) (*this) << val;}
  
  int nfields() {return no_fields;}
};

template<typename T>
OutputStream& OutputStream::operator<< (const T& data) {
  if (fout) *fout << data; 
  else if (dout) *dout << data; 
  else if (obuff) *obuff << data; 
  return *this;
}


template<typename T>
static void dump(Buffer<T>& buff, const string& filename) {
  DelimitedOutput fout(filename);
  fout << buff;
}      

inline static string win_txt(const string& filename) {string out = filename;
#ifdef WINMAGMA
  if (!Utils::ends_with(filename, ".txt")) out.append(".txt");
#endif
  return out;
}

#endif /**OUTPUT_H*/

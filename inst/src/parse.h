/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef PARSE_H
#define PARSE_H

#include <vector>
#include <set>
#include <map>
#include <cstring>
#include "_global.h"
#include "exceptions.h"


class Settings {                         
  typedef string Parameter;
  
  class Flag {
    string flag;
    string parent_flag;
    bool is_modifier;
    bool processed;

    vector<string> values;  
    map<string, Flag> modifiers;  
 
  public:
    Flag(const string& fl="", const string& par="") : flag(fl), parent_flag(par), processed(false) {
      is_modifier = parent_flag != "";
    }
    ~Flag() {
      if (!is_modifier && processed) check_unused_modifiers();
    }
    
    void add_value(const string& val) {values.push_back(val);}  
    string& operator[](int index) {return values[index];}
    const vector<string>& get_values() {return values;}
    
    int size() {return values.size();}
    int nmod() {return modifiers.size();}
    const string& name() {return flag;}
    bool name(const string& check) {return flag == check;}
    bool name(const string& check, int required, bool exact=true);
    bool name(const string& check, pair<int,int> required);
    bool check_size(int required, bool exact=true);
    bool check_size(pair<int,int> required);

    bool set_modifiers(int offset, int min_required=0);
    bool has_modifier(const string& name, bool process=false);
    bool check_modifier(const string& name, int required=0, bool exact=true);
    bool check_modifier(const string& name, pair<int,int> required);
    set<string> check_modifiers(const string& name_list, int required=0, bool exact=true);
    set<string> check_modifiers(const string& name_list, pair<int,int> required);
    void check_unused_modifiers();
    Flag& operator[](const string& name) {modifiers[name].processed = true; return modifiers[name];}
    
    friend ostream& operator<<(ostream& os, const Flag& fl) {
      os << "\t--" << fl.flag;
      for (int i = 0; i < fl.values.size(); i++) os << " " << fl.values[i];
      os << endl;
      return os;
    }
  };
  
  class Value {
  public:
    enum Type {None, Integral, Numeric, Text, List, Set};
    static const char* type_names[6];
  private:
    union Bin {long long integral; double numeric; string* text; vector<string>* list; vector<double>* set;
      Bin() {}
      Bin(long long i) : integral(i) {}
      Bin(double d) : numeric(d) {}      
      Bin(string* s) : text(s) {}      
      Bin(vector<string>* v) : list(v) {}            
      Bin(vector<double>* v) : set(v) {}                  
    };

    Type type;
    Bin contents;

    void realloc() {
      if (type == Text) contents.text = new string(*(contents.text));
      if (type == List) contents.list = new vector<string>(*(contents.list));      
      if (type == Set) contents.set = new vector<double>(*(contents.set));            
    }
    void check_type(Type req) const {if (type != req) throw TypeMismatch(type_names[type], type_names[req]);}
    
  public:  
    Value()                   : type(None) {}
    template<typename T>
    Value(T val)              : type(Integral), contents((long long) val) {}    
    Value(float val)          : type(Numeric),  contents(val) {}        
    Value(double val)         : type(Numeric),  contents(val) {}    
    Value(const char* val)    : type(Text), contents(new string(val)) {}                
    Value(const string& val)  : type(Text), contents(new string(val)) {}            
    Value(const vector<string>& val) : type(List), contents(new vector<string>(val)) {}                
    Value(const vector<double>& val) : type(Set), contents(new vector<double>(val)) {}                    
    Value(const Value& other) : type(other.type), contents(other.contents) {realloc();}
    ~Value() {
      if (type == Text) delete contents.text;
      if (type == List) delete contents.list;      
      if (type == Set) delete contents.set;            
    }
    
    Value& operator= (const Value& other) {
      if (this != &other) {type = other.type; contents = other.contents; realloc();}
      return *this;
    }

    bool null() const {
      if (type == Text) return *(contents.text) == ""; 
      if (type == List) return contents.list->empty(); 
      if (type == Set) return contents.set->empty();       
      return ((type == Integral) && (contents.integral == 0)) || ((type == Numeric) && (contents.numeric == 0));
    }   

    vector<string>& getvs() const {check_type(List); return *(contents.list);}    
    vector<double>& getvn() const {check_type(Set); return *(contents.set);}        
    string    gets() const {check_type(Text); return *(contents.text);}
    long long geti() const {check_type(Integral); return contents.integral;}
    double    getn() const {
      if (type == Integral) return contents.integral;
      check_type(Numeric); return contents.numeric;
    }

    friend ostream& operator<<(ostream& os, const Value& val) {
      if (val.type == Integral) os << val.contents.integral;
      else if (val.type == Numeric) os << val.contents.numeric;
      else if (val.type == Text) os << *(val.contents.text);
      else if (val.type == List) {vector<string>& vec = *(val.contents.list);
        if (!vec.empty()) {
          os << vec[0];      
          for (int i = 1; i < vec.size(); i++) os << " " << vec[i];
        }
      } else if (val.type == Set) {vector<double>& vec = *(val.contents.set);
        if (!vec.empty()) {
          os << vec[0];      
          for (int i = 1; i < vec.size(); i++) os << " " << vec[i];
        }
      } 
      return os;
    }
  };

  int no_args;
  char** raw_args;

  vector<Flag> arguments;
  vector<int> flag_order;
  map<Parameter,Value> parameters;
  set<string> used_flags;
  vector<string> remarks;
  
  time_t start_time;


  void load_args();
  void init_logfile();
  void start_output();
  void debug_verbose();
  
  void set_defaults();
  void set_order_rest() {set_order("");}
  void set_order(const string& flag);
  void check_flag_consistency();
  void process_flags();
  void prologue();

  void print_version();  
  void print_welcome();
  void print_remarks();
  string get_timestr(const time_t& time);
  
  void parse(Flag& arg);
  bool is_flag(char* arg) {return strlen(arg) > 2 && arg[0] == '-' && arg[1] == '-';}
  bool flag_set(const string& flag, bool throw_error=false) const;  
  void required_flags(const string& flag, const string& required);
  void conflicting_flags(const string& flag, const string& conflicting);
  void conflicting_modifiers(const string& flag, const string& mod1, const string& mod2);
  void required_modifier(const string& flag, const string& required, const string& mod);
  void required_one_modifier(const string& flag, vector<string>& req_vec, const string& mod);
  void flag_conflict(const string& flag,  const string& conflicting);
  void flag_conflict_set(const string& conflicting);
  void flag_invalid(const string& flag, const string& msg);
  void flag_requires(const string& flag, const string& required);
  void flag_requires_set(const string& required);
  void flag_requires(const string& flag, const string& required_flag, const string& msg);

  void parse_error(const string& msg);
  void flag_error(const string& flag, const string& msg);
  void mod_error(const string& flag, const string& mod, const string& msg);  
  void add_remark(const string& msg);
  void flag_remark(const string& flag, const string& msg);
  void mod_remark(const string& flag, const string& mod, const string& msg);

  void mod_conflict(Flag& flag, const string& mod_str);
  void mod_conflict(Flag& flag, const string& mod1, const string& mod2);
  void mod_required(Flag& flag, const string& req_str);  
  void mod_required(Flag& flag, const string& mod_str, const string& required, bool symmetric=false);
  void mod_one_required(Flag& flag, const string& req_str) {mod_one_required(flag, "", req_str);}
  void mod_one_required(Flag& flag, const string& mod_str, const string& req_str);

  
  template<typename T, typename T2>
  bool check_range_min(T val, T2 lower, const string& flag, const string& modifier, bool incl) {
    if (val < lower || (!incl && val == lower)) {
      _LOG.error("parsing arguments") << "value for ";
      if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
      _LOG << "--" << flag << " (" << val << ") must be greater than ";
      if (incl) _LOG << "or equal to ";
      _LOG << lower << "." << endl;
      die();
    }
    return true;
  }

  template<typename T, typename T2>
  bool check_range_max(T val, T2 upper, const string& flag, const string& modifier, bool incl) {
    if (val > upper || (!incl && val == upper)) {
      _LOG.error("parsing arguments") << "value for ";
      if (modifier != "") _LOG << "modifier '" << modifier << "' for ";
      _LOG << "--" << flag << " (" << val << ") must be smaller than ";
      if (incl) _LOG << "or equal to ";
      _LOG << upper << "." << endl;
      die();
    }
    return true;
  }
  
  template<typename T, typename T2>
  bool check_range_min(T val, T2 lower, const string& flag, bool incl) {return check_range_min(val, lower, flag, "", incl);}
  
  template<typename T, typename T2>
  bool check_range_max(T val, T2 upper, const string& flag, bool incl) {return check_range_max(val, upper, flag, "", incl);}

  template<typename T, typename T2>
  bool check_range(T val, pair<T2,T2> range, const string& flag, const string& modifier, bool incl) {
    return check_range_min(val, range.first, flag, modifier, incl) && check_range_max(val, range.second, flag, modifier, incl);
  }

  template<typename T, typename T2>
  bool check_range(T val, pair<T2,T2> range, const string& flag, bool incl) {return check_range(val, range, flag, "", incl);}

  bool check_enum(const string& argval, const string& known_str);  
  bool check_enum(const string& argval, const string& known_str, const string& flag, const string& modifier);
  bool check_bool(const string& argval, const string& flag, const string& modifier); ///returns true if argval is empty string
  void check_datatype(const string& req_type, const string& flag) {check_datatype(req_type, flag, "");}
  void check_datatype(const string& req_type, const string& flag, const string& modifier);
  long convert_long(const string& argval, const string& flag, const string& modifier="");
  double convert_double(const string& argval, const string& flag, const string& modifier="");  
  vector<string> read_list(const string& filename, bool by_line, const string& flag, const string& modifier="");
  string check_file(const string& filename, const string& flag, bool check_win) {return check_file(filename, flag, "", check_win);}
  string check_file(const string& filename, const string& flag, const string& modifier="", bool check_win=false);
  const string& check_dir(const string& dirname, const string& flag, const string& modifier="");
  vector<string> expand_prefix(const string& prefix, const string& suffix, const string& flag, long min_files=1);
  void check_suffix(const string& filename, const string& suffix, const string& flag, const string& modifier="");

  template<typename T>
  void set_value(const Parameter& par, const T& val) {parameters[par] = Value(val);}
  void set_marker(const Parameter& par) {parameters[par] = Value(true);}
  void unset_marker(const Parameter& par) {parameters[par] = Value(false);}  
  
  const Value& get_value(const Parameter& par) const {
    map<Parameter,Value>::const_iterator iter = parameters.find(par);
    if (iter == parameters.end()) {_LOG.error("reading parameter") << "trying to access parameter (" << par << ") that doesn't exist" << endl; die();} 
    return iter->second;
  }

  template<typename T> T typed_value(const Parameter& par) const;

  template<typename T>
  T check_value(const Parameter& par) const {T out;
    try {out = typed_value<T>(par);} 
    catch (const TypeMismatch& tm) {_LOG.error("reading parameter") << "trying to access value for parameter '" << par << "' as incompatible type; " << tm.what() << endl; die();} 
    return out;
  }

  void parse_metaweights(const string& filename, vector<string>& meta_genes, const string& type, bool prefix=false);
  void parse_metacorrs(const string& filename, int file_count);

public:
  Settings(int argc, char** argv) : no_args(argc), raw_args(argv) {
    start_time = time(0);
    if (no_args > 1) load_args();
    else {_LOG.error() << "No arguments specified. Please consult manual for usage instructions." << endl; die(ExitType::BenignError);}
    
    if (flag_set("version")) {print_version(); die(ExitType::NoError);}
    
    set_defaults();    
    init_logfile();
    
    print_welcome();
    process_flags();
    start_output();
    if ((*this)["verbose_debug"]) debug_verbose();
    print_remarks();
  }

  void print_goodbye();

  bool isset(const Parameter& par) const {
    bool out = parameters.find(par) != parameters.end();
    return out;
  }

  bool operator[](const Parameter& par) const {return isset(par) && !get_value(par).null();}
  long long geti(const Parameter& par) const {return check_value<long long>(par);}
  double    getn(const Parameter& par) const {return check_value<double>(par);}
  string    gets(const Parameter& par) const {
    if (par == "out_prefix") return get_prefix();
    return check_value<string>(par);
  }
  vector<double> getvn(const Parameter& par) const {return check_value<vector<double> >(par);}  
  vector<string> getvs(const Parameter& par) const {return check_value<vector<string> >(par);}
  
  long long geti(const Parameter& par, long long def) const {return (*this)[par] ? geti(par) : def;}
  double    getn(const Parameter& par, double def) const {return (*this)[par] ? getn(par) : def;}
  string    gets(const Parameter& par, string& def) const {return (*this)[par] ? gets(par) : def;}
  string    gets(const Parameter& par, const char* def) const {return (*this)[par] ? gets(par) : string(def);}
  vector<double> getvn(const Parameter& par, bool def_empty) const {return (*this)[par] || !def_empty ? getvn(par) : vector<double>();}  
  vector<string> getvs(const Parameter& par, bool def_empty) const {return (*this)[par] || !def_empty ? getvs(par) : vector<string>();}

  template<typename T>
  void alter_value(const Parameter& par, const T& val) {set_value(par,val);}

  ///custom getters  
private:  
  vector<string> covar_include;
  vector<string> covar_exclude;
  
  vector<string> gene_covar_include;
  vector<string> gene_covar_exclude;

  vector<Triple<double> > meta_weights;
  Buffer<double> meta_corrs;
  
  string get_prefix() const;

public:
  string get_plinkfile(const string& suffix) const {return gets("plink_prefix").append(".").append(suffix);}
  const vector<string>& get_covar(bool include=true) const {return include ? covar_include : covar_exclude;}
  const vector<string>& get_covarused() const {return !covar_exclude.empty() ? covar_exclude : covar_include;}
  
  const vector<string>& get_genecovar(bool include=true) const {return include ? gene_covar_include : gene_covar_exclude;}  
  const vector<string>& get_genecovarused() const {return !gene_covar_exclude.empty() ? gene_covar_exclude : gene_covar_include;}  

  const vector<Triple<double> >& get_metaweights() const {return meta_weights;}
  const Buffer<double>& get_metacorrs() const {return meta_corrs;}
  
  static int get_version();  
};  

#endif /**PARSE_H*/

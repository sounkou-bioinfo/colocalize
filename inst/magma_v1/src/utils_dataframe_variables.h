/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef UTILS_DATAFRAME_VARIABLES_H
#define UTILS_DATAFRAME_VARIABLES_H

#include <string>
#include <sstream>
#include <limits>
#include <set>
#include <map>

#include "geneutils.h"                                                             
#include "buffer.h"
#include "output.h"

#include "utils_dataframe_variables_storage.h"

class DataFrame;
class Variable {
public: 
  enum MessageType {mt_Filtered};
protected:
  DataFrame* owner;
  bool is_initialised;
  BaseBuffer<int>* output_mask;
  
  void push_message(MessageType type, double value=0);
  void align_variable(Variable* var);
  
public:
  Variable() : owner(0), is_initialised(false), output_mask(0) {}
  virtual ~Variable() {}

  virtual void init() = 0;  
  virtual void clear_data() = 0;  

  DataFrame* get_owner() {return owner;}  
  void set_owner(DataFrame* df) {owner = df;}
  virtual void set_size(int size) = 0;
  virtual void set_map(int size, BaseBuffer<int>* local_map) = 0;  
  virtual void set_mask(BaseBuffer<int>* mask) {output_mask = mask;} ///assumes existing local_map and verified mask

  ///assumes no duplications in vector use
  virtual void shrink(vector<int>& use) = 0;
  virtual void reorder(vector<int>& use) = 0; 

  bool initialised() {return is_initialised;}
  virtual int data_size() = 0;
};


class InputVariable : public Variable {
public:
  enum InputStatus {Valid, Missing, NotNumeric, InvalidValue, InvalidFormat, OutOfRange, OutOfRangeLow, OutOfRangeHigh, TooFewValues, TooManyValues, NotConvertible, MultiError, CustomError};
protected:
  InputStatus status;
  string error_msg;
public:
  InputVariable() : status(Valid) {}
  virtual ~InputVariable() {}  
  
  virtual InputVariable::InputStatus get_status() {return status;}

  virtual void add(int iid, const string& input) = 0;  
  virtual bool exists(int iid) = 0;

  virtual string parse_msg(InputStatus is);
  virtual string& get_msg();
  virtual bool is_valid() {return get_status() == Valid;}
  virtual bool is_valid_or_missing() {return get_status() == Missing || is_valid();}  
};


template<typename T>
class TypedVariable : public InputVariable {
protected:
  VariableStorage<T>& data;
  bool has_missing;
  bool valid_missing;
  T missing_value;
  
  T custom_error(string msg) {this->error_msg = msg; this->status = InputVariable::CustomError; return missing_value;}

public:
  TypedVariable(VariableStorage<T>& data) : data(data), has_missing(false), valid_missing(false) {}
  virtual ~TypedVariable() {}

  virtual void init() {
    if (this->is_initialised) return;
    if (has_missing) data.set_default(missing_value); 
    data.init(); this->is_initialised = true; 
  }
  
  virtual void clear_data() {data.clear();}
  virtual void set_size(int size) {data.set_capacity(size);}
  virtual void set_map(int size, BaseBuffer<int>* local_map) {
    if (size < 0) size = local_map ? local_map->size() : 0;
    set_size(size); data.freeze(local_map);    
  }
   
  virtual void set_missing(T value, bool valid=false) {missing_value = value; has_missing = true; valid_missing = valid;}

  ///assumes no duplications in vector use  
  virtual void shrink(vector<int>& use) {data.remap(use);}
  virtual void reorder(vector<int>& use) {data.remap(use);}

  virtual void filter() {if (!has_missing) return; if (data.filter(missing_value, VariableStorageAccess::Equal) >= 0) Variable::push_message(Variable::mt_Filtered);}
  virtual void filter_range(T value, VariableStorageAccess::FilterMode filter) {if (data.filter(value, filter) >= 0) Variable::push_message(Variable::mt_Filtered);}  
  
  bool exists(int iid) {return this->is_initialised && data.exists(iid);}
  int data_size() {return this->is_initialised ? data.get_datasize() : 0;}
  pair<bool,T> get_missing() {return pair<bool,T>(has_missing, missing_value);}
};

template<typename T> class UniVariableOutput;
template<typename T> class BlockVariableOutput;
template<typename T> class JaggedVariableOutput;
template<typename T> class MultiVariable;
class ColumnBlock;

template<typename T>
class UniVariable : public TypedVariable<T> {
public: typedef T VTYPE;
protected:
  VariableUniStorage<T> data;

  template<typename U>
  bool string_to_value(const string& input, U& value) {
    static istringstream convert; 
    convert.clear(); convert.str(input); convert >> value;
    return convert.eof() && !convert.fail();
  }
  
  virtual T check_value(const T& value) {return value;}
  
public:                      
  UniVariable() : TypedVariable<T>(data) {}
  virtual ~UniVariable() {}
  
  VariableUniStorage<T>& get_data() {return data;}
  BaseBuffer<T>& get_buffer() {return data.get_data();}

  T& get(int iid) {return data.get(iid);}      
  virtual void add(int iid, const string& input) {data.set(iid, convert(input));}
  virtual void set(int iid, T value, bool do_check=false) {data.set(iid, do_check ? validate(value) : value);}

  void assign_missing() {if (this->has_missing) data.assign_value(this->missing_value);}
  bool is_missing(int iid) {
    if (!this->has_missing || this->valid_missing) return !data.exists(iid); 
    else return data.is_missing(iid);
  }
  bool all_missing() {
    if (this->has_missing) {
      VariableStorageAccess::AccessIterator<T>* iterator = this->data.get_iterator();
      while (!iterator->empty()) {if (iterator->next_value() != this->missing_value) return false;}
      return true;
    } else return false;
  }
                   
  T validate(const T& value) {this->status = InputVariable::Valid; return check_value(value);}
  virtual T convert(const string& input) = 0;  

  template<typename VAR>
  void dock_multi(MultiVariable<VAR>& var) {if (this->has_missing) var.set_missing(this->missing_value, this->valid_missing);}  
  virtual UniVariableOutput<T>* get_output(const string& name="", bool force=false) {UNUSED(name); UNUSED(force); return 0;}  
};


template<typename VAR>
class MultiVariable : public TypedVariable<typename VAR::VTYPE> {
  friend class JaggedVariableOutput<VAR>;
protected:
  typedef typename VAR::VTYPE BTYPE;  
  VariableMultiStorage<BTYPE>& data;
  std::set<InputVariable::InputStatus> sub_status;
  VAR base_variable;

  OutputColumn* merged_output_internal(OutputColumn* input, const string& name="", const string& sep=",");

  virtual void add(int iid, const string& input) {add(iid, 0, input);}
  void clear_status(bool clear_sub=false) {this->status = InputVariable::Valid; if (clear_sub) sub_status.clear();}
  void track_status() {if (base_variable.get_status() != InputVariable::Valid) sub_status.insert(base_variable.get_status());}

  using TypedVariable<BTYPE>::get_msg; using TypedVariable<BTYPE>::is_valid; using TypedVariable<BTYPE>::is_valid_or_missing;
  InputVariable::InputStatus get_status() {return (this->status == InputVariable::Valid && !sub_status.empty()) ? InputVariable::MultiError : this->status;}
public:
  MultiVariable(VariableMultiStorage<BTYPE>& data) : TypedVariable<BTYPE>(data), data(data) {base_variable.dock_multi(*this);}
  virtual ~MultiVariable() {}

  void set_missing(BTYPE value, bool valid=false) {TypedVariable<BTYPE>::set_missing(value, valid); base_variable.set_missing(value, valid);}  
  void set_delete(void (*func)(BTYPE&)) {if (this->has_missing) data.set_remove(this->missing_value, func);}
  virtual VariableMultiStorage<BTYPE>& get_data() {return data;}
  
  void set_mask(BaseBuffer<int>* mask) {TypedVariable<BTYPE>::set_mask(mask); base_variable.set_mask(mask);} ///assumes existing local_map and verified mask
  
  BTYPE& get(int iid, int col) {return data.get_elem(iid, col);}
  virtual void add(int iid, int col, const string& input) {clear_status(); data.set_elem(iid, col, base_variable.convert(input));}
  void set(int iid, int col, BTYPE value, bool do_check=false) {clear_status(); data.set_elem(iid, col, do_check ? base_variable.validate(value) : value);}
  VariableStorageAccess::MultiAccessElement<BTYPE>& get_row(int iid, bool to_set) {return data.get_row(iid, to_set);}
    
  InputVariable::InputStatus get_status(bool use_base) {return use_base ? base_variable.get_status() : this->get_status();}
  string& get_msg(bool use_base) {return use_base ? base_variable.get_msg() : this->get_msg();}
  bool is_valid(bool use_base) {return use_base ? base_variable.is_valid() : this->is_valid();}
  bool is_valid_or_missing(bool use_base) {return use_base ? base_variable.is_valid_or_missing() : this->is_valid_or_missing();}
  
  std::set<InputVariable::InputStatus> get_block_status() {return sub_status;}
  vector<string> get_msg_sub() {vector<string> out;
    if (get_status() == InputVariable::MultiError && !sub_status.empty()) {    
      for (std::set<InputVariable::InputStatus>::iterator it = sub_status.begin(); it != sub_status.end(); ++it) out.push_back(base_variable.parse_msg(*it));
    }
    return out;
  }

  BTYPE convert(const string& input) {return base_variable.convert(input);}
  
  virtual void assign(int iid, BaseBuffer<BTYPE>& input, bool do_check=false) = 0;
  
  virtual OutputColumn* get_output() = 0;
  virtual OutputColumn* merged_output(const string& name="", const string& sep=",") = 0;
};

template<typename VAR>
class BlockMultiVariable : public MultiVariable<VAR> {
  using typename MultiVariable<VAR>::BTYPE;
  VariableMultiBlockStorage<BTYPE> data;

public:
  BlockMultiVariable() : MultiVariable<VAR>(data) {}
  
  using MultiVariable<VAR>::init;
  void init(int width) {if (!this->is_initialised) {data.set_width(width); init();}}
  
  int get_width() {return data.width();}
  void set_width(int width) {data.set_width(width);}
  
  VariableMultiBlockStorage<BTYPE>& get_data() {return data;}  
  Buffer<BTYPE>& get_buffer() {return data.get_data();}

  void assign_missing() {if (this->has_missing) data.assign_value(this->missing_value);}  
  void assign(int iid, BaseBuffer<BTYPE>& input, bool do_check=false) {this->clear_status(true);
    int use = min(input.size(), long(data.width())); 
    for (int i = 0; i < use; i++) {this->set(iid, i, input[i], do_check); this->track_status();}
    if (input.size() != data.width()) this->status = (input.size() < data.width()) ? InputVariable::TooFewValues : InputVariable::TooManyValues;
  }

  using MultiVariable<VAR>::add;  
  void add(int iid, const string& input, char sep) {this->clear_status(true);
    vector<string> parts = Utils::tokenize(input, sep); int use = min(int(parts.size()), data.width());
    for (int i = 0; i < use; i++) {add(iid, i, parts[i]); this->track_status();}
    if (parts.size() != data.width()) this->status = (parts.size() < data.width()) ? InputVariable::TooFewValues : InputVariable::TooManyValues;    
  }
  
  ColumnBlock* get_output();
  UniVariableOutput<typename VAR::VTYPE>* get_output(int col, const string& name="");

  OutputColumn* merged_output(const string& name="", const string& sep=",");
  OutputColumn* merged_output(vector<int>& columns, const string& name="", const string& sep=",");
  OutputColumn* merged_output(vector<pair<int,string> >& columns, const string& name="", const string& sep=",");
};   

template<typename VAR>
class JaggedMultiVariable : public MultiVariable<VAR> {
protected:
  using typename MultiVariable<VAR>::BTYPE;
  VariableMultiJaggedStorage<BTYPE> data;

public:
  JaggedMultiVariable(bool do_purge=true) : MultiVariable<VAR>(data) {data.set_purge(do_purge);}

  VariableMultiJaggedStorage<BTYPE>& get_data() {return data;}    
  BaseBuffer<BaseBuffer<BTYPE>*>& get_buffer() {return data.get_data();}
  BaseBuffer<BTYPE>& get_buffer(int iid) {return data.get_data(iid);}  

  using MultiVariable<VAR>::add; using MultiVariable<VAR>::set;
  void add(int iid, const string& input) {data.add_elem(iid, this->base_variable.convert(input));}
  void set(int iid, BTYPE value, bool do_check=false) {data.add_elem(iid, do_check ? this->base_variable.validate(value) : value);}

  void assign_missing() {if (this->has_missing) data.remove_data();}
  void assign(int iid, BaseBuffer<BTYPE>& input, bool do_check=false) {
    if (do_check) {this->clear_status(true);
      for (int i = 0; i < input.size(); i++) {input[i] = this->base_variable.validate(input[i]); this->track_status();}
    }
    data.assign_data(iid, input);
  }

  template<typename C>
  void insert(int iid, C& input) {data.insert_data(iid, input);}
  
  template<typename C>
  void insert(int iid, C& input, bool do_check) {this->clear_status(true);
    if (do_check) {for (typename C::iterator it = input.begin(); it != input.end(); ++it) {*it = this->base_variable.validate(*it); this->track_status();}}
    data.insert_data(iid, input);
  }

  virtual JaggedVariableOutput<VAR>* get_output();
  OutputColumn* merged_output(const string& name="", const string& sep=",") {return this->merged_output_internal(get_output(), name, sep);}
};  

  

/** SPECIALISATIONS **/

template<typename T>
class PointerVariable : public UniVariable<T*> {
public: using typename UniVariable<T*>::VTYPE;
protected:
  bool ptr_delete;
  static void clear_pointer(T*& ptr) {delete ptr; ptr = 0;}
  
  using UniVariable<T*>::add;
public:
  PointerVariable(bool do_delete) : ptr_delete(false) {this->set_missing(0); if (do_delete) set_delete(true);}
  virtual ~PointerVariable() {}
  
  void set_delete(bool value) {ptr_delete = value; this->data.set_remove(this->missing_value, ptr_delete ? &clear_pointer : 0);}

  virtual T* convert(const string& input) {UNUSED(input); this->status = InputVariable::NotConvertible; return 0;}

  virtual string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::NotConvertible) return "cannot convert value into pointer"; 
    else return UniVariable<T*>::parse_msg(is);
  }

  template<typename VAR>  
  void dock_multi(MultiVariable<VAR>& var) {UniVariable<T*>::dock_multi(var); if (ptr_delete) var.set_delete(&clear_pointer);}    
};

class StringVariable : public UniVariable<string> {
public: using typename UniVariable<string>::VTYPE;
  StringVariable(int init_size=0) {
    this->set_missing("NA"); 
    if (init_size > 0) {
      this->set_size(init_size);
      this->init();
    }    
  }

  string convert(const string& input) {
    if (input == "NA") {this->status = InputVariable::Missing; return input;} 
    else return this->validate(input);
  }
  UniVariableOutput<string>* get_output(const string& name="", bool force=false);
  UniVariableOutput<string>* get_output_right(const string& name="", bool force=false);
};


class CStringVariable : public PointerVariable<char> {
public: using typename PointerVariable<char>::VTYPE;
private:
  ArchiveBuffer archive;
  int init_size;

public:
  CStringVariable() : PointerVariable<char>(false), init_size(0) {}
  
  void init() {if (is_initialised) return; PointerVariable::init(); archive.init(init_size);}

  void clear_data() {clear_data(false);}
  void clear_data(bool keep_archive) {PointerVariable::clear_data(); if (!keep_archive) archive.clear();}
  void set_size(int size) {PointerVariable::set_size(size); init_size = size;}

  using PointerVariable<char>::add;

  virtual char* convert(const string& input) {return archive.add(input);}
  
  UniVariableOutput<char*>* get_output(const string& name="", bool force=false);
  UniVariableOutput<char*>* get_output_right(const string& name="", bool force=false);
};  


template<typename T>
class NumericVariable : public UniVariable<T> {
public: using typename UniVariable<T>::VTYPE;
protected:
  bool do_round;
public:
  NumericVariable() : do_round(false) {this->set_missing(numeric_limits<T>::max());}
  virtual ~NumericVariable() {} 
  virtual T convert(const string& input) {T value;
    if (input == "NA") {this->status = InputVariable::Missing; return this->missing_value;}
    else if (!this->string_to_value(input, value)) {double dvalue;
      if (!do_round || !this->string_to_value(input, dvalue)) {this->status = InputVariable::NotNumeric; return this->missing_value;}
      else value = round(dvalue);
    }
    this->status = InputVariable::Valid; return this->check_value(value);
  }  
  void set_rounding(bool value) {do_round = value && numeric_limits<T>::is_integer;}
  virtual UniVariableOutput<T>* get_output(const string& name="", bool force=false);
};

template<typename T>
class PositiveNumber : public NumericVariable<T> {
public: using typename NumericVariable<T>::VTYPE;
protected:
  T min_value;
  
  virtual T check_value(const T& value) {
    if (value <= min_value) {
      this->status = value >= 0 ? InputVariable::OutOfRange : InputVariable::InvalidValue;
      return this->missing_value;
    } else return value;
  }

public:
  PositiveNumber() : min_value(0) {this->set_missing(0);}
  virtual ~PositiveNumber() {}
  
  void set_bound(T min) {min_value = max(min, (T) 0); this->set_missing(0, min_value == 0);}

  virtual string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::InvalidValue) return "value is negative"; 
    else if (is == InputVariable::OutOfRange) return "value is below lower bound of " + Utils::num_string(min_value); 
    else return NumericVariable<T>::parse_msg(is);
  }
};

template<typename T>
class NonNegativeNumber : public NumericVariable<T> {
public: using typename NumericVariable<T>::VTYPE;
private:
  virtual T check_value(const T& value) {
    if (value < 0) {this->status = InputVariable::InvalidValue; return this->missing_value;} 
    else return value;
  }

public:
  NonNegativeNumber() {this->set_missing(0, true);}

  virtual string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::InvalidValue) return "value is negative"; 
    else return NumericVariable<T>::parse_msg(is);
  }
};

template<typename T>
class NumericID : public NonNegativeNumber<T> {
public: using typename NonNegativeNumber<T>::VTYPE;
public:
  NumericID() {this->set_missing(-1);}
};

class BooleanVariable : public NumericVariable<bool> {
public: using typename NumericVariable<bool>::VTYPE;
public:
  BooleanVariable() {this->set_missing(false, true);}
};

template<typename T>
class BoundedNumericVariable : public NumericVariable<T> {
public: using typename NumericVariable<T>::VTYPE;
protected:
  T min_value; T min_allowed;
  T max_value; T max_allowed;
  
  virtual T check_value(const T& value) {
    if (value < min_allowed || value > max_allowed) {this->status = InputVariable::InvalidValue; return this->missing_value;}
    else if (value < min_value) {this->status = InputVariable::OutOfRangeLow; return min_value;}
    else if (value > max_value) {this->status = InputVariable::OutOfRangeHigh; return max_value;} 
    else return value;   
  }

public:
  BoundedNumericVariable() : min_value(numeric_limits<T>::min()), max_value(numeric_limits<T>::max()-1) {min_allowed = min_value; max_allowed = max_value;}
  BoundedNumericVariable(T min_value, T max_value) : min_value(min_value), max_value(max_value) {min_allowed = numeric_limits<T>::min(); max_allowed = numeric_limits<T>::max();}
  BoundedNumericVariable(T min_value, T max_value, T min_allowed, T max_allowed) : min_value(min_value), min_allowed(min_allowed),  max_value(max_value), max_allowed(max_allowed) {}
  virtual ~BoundedNumericVariable() {}

  void set_bounds(T min, T max) {
    min_value = min >= min_allowed ? min : min_allowed;
    max_value = max <= max_allowed ? max : max_allowed;
  }
  
  void change_bounds(T min, T max) {
    bool truncate = min > min_value || max < max_value; set_bounds(min, max);  
    if (truncate) {
      VariableStorageAccess::AccessIterator<T>* iterator = this->data.get_iterator();
      while (!iterator->empty()) {T& value = iterator->next_value();
        if (value < min_value) value = min_value;
        else if (value > max_value) value = max_value;
      }
    }
  }
  
  virtual string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::InvalidValue) return "value is out of bounds"; 
    else if (is == InputVariable::OutOfRangeLow) return "value is below minimum value of " + Utils::num_string(min_value); 
    else if (is == InputVariable::OutOfRangeHigh) return "value is above maximum value of " + Utils::num_string(max_value);     
    else return NumericVariable<T>::parse_msg(is);
  }
};

class EnumVariable : public NumericVariable<short> {
public: using typename NumericVariable<short>::VTYPE;
private:
  map<string,short> index;
  vector<string> labels;
  
  short convert(const string& input);
  short check_value(const short& value);

public:
  EnumVariable(string label_str="", string default_value=""); 

  using NumericVariable<short>::set;
  void set(int iid, const string& label) {data.set(iid, get_code(label));}

  vector<string>& get_labels() {return labels;}
  void set_labels(const string& label_str);
  
  short get_code(const string& label);

  UniVariableOutput<short>* get_output(const string& name="", bool force=false); 
}; 


template<typename T>
class PositiveWeightVariable : public NonNegativeNumber<T> {
public: using typename NonNegativeNumber<T>::VTYPE;
  PositiveWeightVariable() {this->set_missing(-1);}
};

template<typename T>
class CorrelationVariable : public BoundedNumericVariable<T> {
public: using typename BoundedNumericVariable<T>::VTYPE;
public: 
  CorrelationVariable() : BoundedNumericVariable<T>(-1, 1, -1.1, 1.1) {this->set_missing(0, true);}

  string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::InvalidValue || is == InputVariable::OutOfRangeLow || is == InputVariable::OutOfRangeLow) 
      return "value is not a valid correlation";
    else return BoundedNumericVariable<T>::parse_msg(is);
  }
};

template<typename T>
class PosCorrelationVariable : public BoundedNumericVariable<T> {
public: using typename BoundedNumericVariable<T>::VTYPE;
public: 
  PosCorrelationVariable() : BoundedNumericVariable<T>(0, 1, -0.1, 1.1) {this->set_missing(0, true);}

  string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::InvalidValue || is == InputVariable::OutOfRangeLow || is == InputVariable::OutOfRangeLow) return "value is not a valid (unsigned) correlation";
    else return BoundedNumericVariable<T>::parse_msg(is);
  }
};


template<typename T>
class PvalueVariable : public BoundedNumericVariable<T> {
public: using typename BoundedNumericVariable<T>::VTYPE;
public: 
  PvalueVariable() : BoundedNumericVariable<T>(numeric_limits<T>::min(), 1-numeric_limits<T>::epsilon(), 0, 1) {this->set_missing(-1);}
  virtual ~PvalueVariable() {}
  virtual T convert(const string& input) {
    if (input == "-1") {this->status = InputVariable::Missing; return this->missing_value;}
    return BoundedNumericVariable<T>::convert(input);
  }
  virtual string parse_msg(InputVariable::InputStatus is) {
    if (is == InputVariable::InvalidValue || is == InputVariable::OutOfRangeLow || is == InputVariable::OutOfRangeLow) return "value is not a valid p-value";
    else return "p-" + BoundedNumericVariable<T>::parse_msg(is);
  }
  bool is_valid() {InputVariable::InputStatus status = this->get_status(); return (status == InputVariable::Valid || status == InputVariable::OutOfRange || status == InputVariable::OutOfRangeLow || status == InputVariable::OutOfRangeHigh);}
};

template<typename T>
class PermPvalueVariable : public PvalueVariable<T> {
public: using typename PvalueVariable<T>::VTYPE;
private:
  T lower_bound; T upper_bound; 
  
  using PvalueVariable<T>::add;
  using PvalueVariable<T>::set;    
public: 
  PermPvalueVariable() : lower_bound(numeric_limits<T>::min()), upper_bound(1-numeric_limits<T>::epsilon()) {}

  void add(int iid, const string& input, int perm) {set_nperm(perm); add(iid, input);}
  void set(int iid, T value, int perm) {set_nperm(perm); set(iid, value, true);}
  
  void set_nperm(int perm) {perm = max(perm, 1000); 
    T lower = 1.0/perm, upper = 1 - 1.0/perm; 
    this->set_bounds(lower > 0 ? lower : lower_bound, upper < 1 ? upper : upper_bound);
  }
};


template<typename T> 
class ZscoreVariable : public BoundedNumericVariable<T> {
public: using typename BoundedNumericVariable<T>::VTYPE;
public:
  ZscoreVariable() : BoundedNumericVariable<T>(-100, 100, -500, 500) {this->set_missing(-999);}
};

class GeneLocationVariable : public UniVariable<GeneLocation> {
public: using typename UniVariable<GeneLocation>::VTYPE;
private:
  unsigned long max_value;
  bool is_human;

public:
  GeneLocationVariable(bool is_human=true) : max_value(numeric_limits<unsigned long>::max()), is_human(is_human) {this->set_missing(GeneLocation());}
  bool& human() {return is_human;}

  void filter_chromosome(int chr);

  using UniVariable<GeneLocation>::add;
  void add(int iid, const string& chr, const string& begin, const string& end) {data.set(iid, convert(chr, begin, end));}  
  GeneLocation convert(const string& input);
  GeneLocation convert(const string& chr, const string& begin, const string& end);  
  UniVariableOutput<GeneLocation>* get_output(const string& name="", bool force=false);
};

class SparseBooleanVariable : public Variable {
  vector<std::set<int> > data;
  BaseBuffer<int>* map;
  
  int sid_max;
  bool truncate;

  int get_sid(int iid) {int sid = map ? map->get(iid)-1 : iid; return sid < sid_max ? sid : -1;}
  void prune_storage(int max);
  void prune_storage(std::set<int>& values, int max);

public:
  SparseBooleanVariable() : map(0), sid_max(0), truncate(false) {}
  
  void init() {if (this->is_initialised) return; this->is_initialised = true;}
  void clear_data() {data.clear();}
  
  vector<std::set<int> >& get_data() {return data;}
  std::set<int>* get_data(int col) {return (col >= 0 && col < data.size()) ? &data[col] : 0;}
  
  void set_size(int size) {if (this->is_initialised) prune_storage(size); sid_max = size;}
  void set_map(int size, BaseBuffer<int>* local_map);

  ///assumes no duplications in vector use
  void shrink(vector<int>& use) {reorder(use);}
  void reorder(vector<int>& use);

  int count(int col) {return col <= data.size() ? data[col].size() : 0;}
  int data_size() {return sid_max;}
  int get_width() {return data.size();}
  void set_width(int ncol) {data.resize(ncol);}  

  bool get(int iid, int col) {return data[col].find(get_sid(iid)) != data[col].end();}
  void set(int iid, int col, bool value);

  int insert(std::set<int>& values) {prune_storage(values, sid_max); data.resize(data.size()+1); data.back().swap(values); return data.size()-1;}
  void update(int col, std::set<int>& values) {if (col < data.size()) {prune_storage(values, sid_max); data[col].swap(values);}}  

  ColumnBlock* get_output();  
  OutputColumn* get_output(int col, const string& name="");
};


/** OUTPUT **/

template<typename T>
class ConstantOutput : public OutputColumn {
  T value;
  string value_str;
  
  int length;
  int current;
  
  void init() {
    value_str = Utils::to_string(value);
    if (name == "") name = value_str;
  }
  
public:
  ConstantOutput(T value, int len, const string& name="") : OutputColumn(name), value(value), length(len), current(0) {init();}

  virtual void print(OutputStream& out) {out << value_str; current++;}
    
  virtual void set_header(FormattedOutput& fout) {fout.add_field(value_str.size(), name);} 
  virtual bool is_empty() {return current >= length;}
  virtual void rewind() {current = 0;}
};


template<typename T>
class UniVariableOutput : public OutputColumn {
protected:
  VariableStorageAccess::AccessIterator<T>* iterator;
public:
  UniVariableOutput(UniVariable<T>& variable, const string& name="", BaseBuffer<int>* mask=0) : OutputColumn(name) {
    if (variable.initialised()) iterator = variable.get_data().get_iterator(mask); 
    else iterator = new VariableStorageAccess::DummyAccessIterator<T>();
  }
  virtual ~UniVariableOutput() {delete iterator;}

  void set_data(VariableStorageAccess::AccessIterator<T>* update) {delete iterator; iterator = update;}
  virtual void print(OutputStream& out) {print_value(out, iterator->next_value());}
  virtual void print_value(OutputStream& out, T& value) {out << value;}
  bool is_empty() {return iterator->empty();}    
  void rewind() {iterator->reset();}
};

class ColumnBlock : public OutputColumn {
protected:
  vector<OutputColumn*> columns;

  virtual OutputColumn* get_column(int index, const string& name="") = 0;
public:
  virtual ~ColumnBlock() {for (int i = 0; i < columns.size(); i++) delete columns[i];}    

  void add_column(OutputColumn* add) {if (add && !add->is_empty()) columns.push_back(add);}
  void add_column(int index, const string& name="") {add_column(get_column(index, name));}
  
  ColumnBlock* add_columns(int from, int to) {for (int i = from; i < to; i++) add_column(i); return this;}
  ColumnBlock* add_columns(vector<int>& column_index) {for (int i = 0; i < column_index.size(); i++) add_column(column_index[i]); return this;}
  ColumnBlock* add_columns(vector<pair<int,string> >& column_index) {for (int i = 0; i < column_index.size(); i++) add_column(column_index[i].first, column_index[i].second); return this;}

  void print(OutputStream& out) {for (int i = 0; i < columns.size(); i++) columns[i]->print(out);}
  void set_header(FormattedOutput& fout) {for (int i = 0; i < columns.size(); i++) columns[i]->set_header(fout);}
  bool is_empty() {return columns.empty() || columns[0]->is_empty();}    
  void rewind() {for (int i = 0; i < columns.size(); i++) columns[i]->rewind();}
};

template<typename VAR>
class BlockVariableOutput : public ColumnBlock {
protected:
  VAR& base_variable;  
  OutputColumn* get_column(int index, const string& name="") {return base_variable.get_output(index, name);}  
public:
  BlockVariableOutput(VAR& base) : base_variable(base) {}
}; 

template<typename VAR>
class JaggedVariableOutput : public OutputColumn {
protected:
  typedef typename VAR::VTYPE CTYPE;
  UniVariableOutput<CTYPE>* base_output;
  VariableStorageAccess::AccessIterator<BaseBuffer<CTYPE>*>* iterator;
  
public:
  JaggedVariableOutput(JaggedMultiVariable<VAR>& variable, BaseBuffer<int>* mask=0) {
    base_output = variable.base_variable.get_output("", true);
    if (variable.initialised()) iterator = variable.get_data().get_iterator(mask); 
    else iterator = new VariableStorageAccess::DummyAccessIterator<BaseBuffer<CTYPE>*>();
  }
  ~JaggedVariableOutput() {delete base_output; delete iterator;}

  void print(OutputStream& out) {
    BaseBuffer<CTYPE>*& buffer = iterator->next_value();
    if (buffer) {for (int i = 0; i < buffer->size(); i++) base_output->print_value(out, buffer->get(i));}
  }
  bool is_empty() {return iterator->empty();}   
  void rewind() {iterator->reset();}   
};

class CompoundOutput : public OutputColumn {
  vector<string> data;
  string separator;
  int current;

  void process_data(OutputColumn* source);
public:
  CompoundOutput(OutputColumn* source, const string& name="", const string& sep=",") : OutputColumn(name), separator(sep), current(0) {process_data(source);}

  void print(OutputStream& out) {out << data[current++];}
  void set_header(FormattedOutput& fout);
  bool is_empty() {return current >= data.size();}
  void rewind() {current = 0;}
};

template<typename T>
class NumericOutput : public UniVariableOutput<T> {
  bool has_missing; T missing_value;  
public:
  NumericOutput(NumericVariable<T>& variable, const string& name="", BaseBuffer<int>* mask=0) : UniVariableOutput<T>(variable, name, mask), has_missing(false) {}
  NumericOutput(NumericVariable<T>& variable, T missing_value, const string& name="", BaseBuffer<int>* mask=0) : UniVariableOutput<T>(variable, name, mask), has_missing(true), missing_value(missing_value) {}
  
  void print_value(OutputStream& out, T& value) {
    if (has_missing && value == missing_value) out << "NA";
    else out << value;
  }
  void set_header(FormattedOutput& fout) {
    if (numeric_limits<T>::is_integer) {
      if (this->width == 0) {T max_neg = 0, max_pos = 0;
        while (!this->iterator->empty()) {T& value = this->iterator->next_value();
          if (value < max_neg) max_neg = value;        
          else if (value > max_pos) max_pos = value;
        }
        if (10*abs(max_neg) > max_pos) max_pos = 10*abs(max_neg);
        this->width = max_pos > 1 ? floor(log10(max_pos))+1 : 1;
        this->iterator->reset();
      }
      fout.add_field(this->width, this->name);    
    } else {
      if (this->width == 0) this->width = 5;
      fout.add_float_field(this->width, this->name);
    }
  }
};


class StringOutput : public UniVariableOutput<string> {
  bool align_right;
public:
  StringOutput(StringVariable& variable, bool right, const string& name="", BaseBuffer<int>* mask=0) : UniVariableOutput<string>(variable, name, mask), align_right(right) {}
  void set_header(FormattedOutput& fout);  
};

class CStringOutput : public UniVariableOutput<char*> {
  bool align_right;
public:
  CStringOutput(CStringVariable& variable, bool right, const string& name="", BaseBuffer<int>* mask=0) : UniVariableOutput<char*>(variable, name, mask), align_right(right) {}
  void set_header(FormattedOutput& fout);  
};

class GeneLocationOutput : public UniVariableOutput<GeneLocation> {
  bool formatted;
public:
  GeneLocationOutput(GeneLocationVariable& variable, BaseBuffer<int>* mask=0) : UniVariableOutput<GeneLocation>(variable, "#LOCATION", mask), formatted(false) {}
  
  void set_header(FormattedOutput& fout);
  void print_value(OutputStream& out, GeneLocation& value);
};

class EnumOutput : public UniVariableOutput<short> {
  vector<string>& labels;
public:
  EnumOutput(EnumVariable& variable, const string& name="", BaseBuffer<int>* mask=0) : UniVariableOutput<short>(variable, name, mask), labels(variable.get_labels()) {}
  
  void set_header(FormattedOutput& fout);
  void print_value(OutputStream& out, short& value) {out << labels[value];}
};

class SparseBooleanOutput : public OutputColumn {
  set<int>& values; BaseBuffer<int>* map;
  set<int>::iterator next_value;
  int curr_index; int index_max;

  void find_next() {if (map) {while (curr_index < index_max && map->get(curr_index) == 0) curr_index++;}} 
public:
  SparseBooleanOutput(set<int>& values, int size, const string& name="") : OutputColumn(name), values(values), map(0), index_max(size) {rewind();}
  SparseBooleanOutput(set<int>& values, BaseBuffer<int>* local_map, const string& name="") : OutputColumn(name), values(values), map(local_map), index_max(local_map->size()) {rewind();}

  void print(OutputStream& out);
  bool is_empty() {return curr_index >= index_max;}
  void rewind() {curr_index = 0; next_value = values.begin(); find_next();}
};


template<typename T>
UniVariableOutput<T>* NumericVariable<T>::get_output(const string& name, bool force) {
  if (this->is_initialised || force) {
    if (this->has_missing && !this->valid_missing) return new NumericOutput<T>(*this, this->missing_value, name, this->output_mask);
    else return new NumericOutput<T>(*this, name, this->output_mask);
  } else return 0;
}

template<typename VAR>
OutputColumn* MultiVariable<VAR>::merged_output_internal(OutputColumn* input, const string& name, const string& sep) {
  OutputColumn* output = input ? new CompoundOutput(input, name, sep) : 0;
  delete input; return output;
}

template<typename VAR>
ColumnBlock* BlockMultiVariable<VAR>::get_output() {return this->is_initialised ? new BlockVariableOutput<BlockMultiVariable<VAR> >(*this) : 0;}

template<typename VAR>
UniVariableOutput<typename VAR::VTYPE>* BlockMultiVariable<VAR>::get_output(int col, const string& name) {UniVariableOutput<typename VAR::VTYPE>* out = 0;
  if (this->is_initialised && col < data.width()) out = this->base_variable.get_output(name, true);
  if (out) out->set_data(col >= 0 ? data.get_iterator(col, this->output_mask) : data.get_constant_iterator());
  return out;
}

template<typename VAR>
OutputColumn* BlockMultiVariable<VAR>::merged_output(const string& name, const string& sep) {
  ColumnBlock* out = get_output();
  return out ? this->merged_output_internal(out->add_columns(0, data.width()), name, sep) : 0;
}

template<typename VAR>
OutputColumn* BlockMultiVariable<VAR>::merged_output(vector<int>& columns, const string& name, const string& sep) {
  ColumnBlock* out = get_output();
  return out ? this->merged_output_internal(out->add_columns(columns), name, sep) : 0;
}

template<typename VAR>
OutputColumn* BlockMultiVariable<VAR>::merged_output(vector<pair<int,string> >& columns, const string& name, const string& sep) {
  ColumnBlock* out = get_output();
  return out ? this->merged_output_internal(out->add_columns(columns), name, sep) : 0;
}


template<typename VAR>
JaggedVariableOutput<VAR>* JaggedMultiVariable<VAR>::get_output() {return this->is_initialised ? new JaggedVariableOutput<VAR>(*this, this->output_mask) : 0;}

#include "utils_dataframe_variables_genecorrelation.h"

#endif /** UTILS_DATAFRAME_VARIABLES_H */

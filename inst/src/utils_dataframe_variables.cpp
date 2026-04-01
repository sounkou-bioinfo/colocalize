/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "utils_dataframe_variables.h"
#include "utils_dataframe.h"
#include "utils.h"

void Variable::push_message(MessageType type, double value) {if (owner) owner->push_message(type, value);}
void Variable::align_variable(Variable* var) {if (owner) owner->align_variable(*var);}


string InputVariable::parse_msg(InputVariable::InputStatus status) {
  switch (status) { 
    case Missing: return "value is missing";
    case NotNumeric: return "value is not a number";
    case InvalidValue: return "value is invalid";
    case InvalidFormat: return "value string is improperly formatted";
    case OutOfRange: return "value is out of range";
    case OutOfRangeLow: return "value is too low";    
    case OutOfRangeHigh: return "value is too high";        
    case TooFewValues: return "not enough values in input";
    case TooManyValues: return "too many values in input";    
    case NotConvertible: return "value string cannot be converted";
    case MultiError: return "error(s) occurred in block input";
    case CustomError: return error_msg;
    default: return "";
  }
}

string& InputVariable::get_msg() {
  if (get_status() == Valid) error_msg = ""; 
  else if (get_status() != CustomError) error_msg = parse_msg(get_status());
  return error_msg;
}

EnumVariable::EnumVariable(string label_str, string default_value) {
  this->set_missing(0, default_value != ""); 
  set_labels(label_str);
  if (default_value != "") {
    labels[0] = default_value;
    index[default_value] = 0;
  }
}     

short EnumVariable::convert(const string& input) {
  if (input != "NA") {
    map<string,short>::iterator found = index.find(input);
    if (found != index.end()) {this->status = InputVariable::Valid; return found->second;}
    else {this->status = InputVariable::InvalidValue; return this->missing_value;}
  } else {this->status = InputVariable::Missing; return this->missing_value;}
}  

short EnumVariable::check_value(const short& value) {
  if (value == 0) {this->status = InputVariable::Missing; return this->missing_value;} 
  else if (value < 0 || value >= labels.size()) {this->status = InputVariable::InvalidValue; return this->missing_value;} 
  else return value;
}

void EnumVariable::set_labels(const string& label_str) {
  if (labels.size() < 2) {
    labels = Utils::tokenize(label_str);
    labels.insert(labels.begin(), "NA");

    index.clear();
    for (int i = 0; i < labels.size(); i++) index[labels[i]] = i;
  }
}

short EnumVariable::get_code(const string& label) {
  map<string,short>::iterator found = index.find(label);  
  return (found != index.end()) ? found->second : 0;
}


UniVariableOutput<string>* StringVariable::get_output(const string& name, bool force) {return (is_initialised || force) ? new StringOutput(*this, false, name, output_mask) : 0;} 
UniVariableOutput<string>* StringVariable::get_output_right(const string& name, bool force) {return (is_initialised || force) ? new StringOutput(*this, true, name, output_mask) : 0;} 
UniVariableOutput<char*>* CStringVariable::get_output(const string& name, bool force) {return (is_initialised || force) ? new CStringOutput(*this, false, name, output_mask) : 0;} 
UniVariableOutput<char*>* CStringVariable::get_output_right(const string& name, bool force) {return (is_initialised || force) ? new CStringOutput(*this, true, name, output_mask) : 0;} 
UniVariableOutput<GeneLocation>* GeneLocationVariable::get_output(const string& name, bool force) {UNUSED(name); return (is_initialised || force) ? new GeneLocationOutput(*this, output_mask) : 0;} 
UniVariableOutput<short>* EnumVariable::get_output(const string& name, bool force) {return (is_initialised || force) ? new EnumOutput(*this, name, output_mask) : 0;} 

void GeneLocationVariable::filter_chromosome(int chr) {
  if (!data.initialised() || !data.has_map()) return;
  VariableStorageAccess::AccessIterator<GeneLocation>* iterator = data.get_iterator();

  while (!iterator->empty()) {
    pair<int,GeneLocation> curr = iterator->next_element();
    if (curr.second.chromosome == chr) data.drop_iid(curr.first);
  }
  Variable::push_message(Variable::mt_Filtered);
}

GeneLocation GeneLocationVariable::convert(const string& input) {
  static istringstream split; static string parts[3]; string loc_str = input;
  for (int i = 0; i < loc_str.size(); i++) {if (loc_str[i] == ':') loc_str[i] = ' ';}

  split.clear(); split.str(loc_str); 
  for (int i = 0; i < 3; i++) split >> parts[i];
  if (!split.eof() || split.fail()) {this->status = InputVariable::InvalidFormat; return this->missing_value;}

  return convert(parts[0], parts[1], parts[2]);
} 

GeneLocation GeneLocationVariable::convert(const string& chr_s, const string& begin_s, const string& end_s) {
  int chr; unsigned long begin, end;
  if (!Utils::chr_val(chr_s, chr, is_human)) return custom_error(string("invalid genomic location; chromosome code ") + Utils::quote(chr_s) + " not recognised");
  if (begin_s[0] == '-' || end_s[0] == '-' || !string_to_value(begin_s, begin) || !string_to_value(end_s, end)) return custom_error("invalid genomic location; gene boundaries must be positive integers");
  if (begin >= max_value || end >= max_value) return custom_error("invalid genomic location; integer overflow in gene boundaries");
  if (end <= begin) return custom_error("invalid genomic location; gene stop position must be greater than start position");
      
  this->status = InputVariable::Valid;
  return GeneLocation(chr, begin, end);
}

void SparseBooleanVariable::prune_storage(int max) {
  if (max < sid_max) {
    for (int i = 0; i < data.size(); i++) prune_storage(data[i], max);
  }
}

void SparseBooleanVariable::prune_storage(std::set<int>& values, int max) {
  if (!values.empty() && *(values.rbegin()) >= max) {
    if (*(values.begin()) < max) values.erase(values.lower_bound(max), values.end());
    else values.clear(); 
  }
}

void SparseBooleanVariable::set_map(int size, BaseBuffer<int>* local_map) {
  if (size < 0) size = local_map ? local_map->size() : 0;
  prune_storage(size); set_size(size); map = local_map; 
}

void SparseBooleanVariable::reorder(vector<int>& use) {
  if (!this->is_initialised || use.size() >= sid_max) return;
  vector<int> remap(sid_max, -1);
  for (int i = 0; i < use.size(); i++) {remap[use[i]] = i;}

  for (int i = 0; i < data.size(); i++) {
    if (!data[i].empty()) {
      std::set<int> &curr = data[i], update;
      for (std::set<int>::iterator it = curr.begin(); it != curr.end(); ++it) {if (remap[*it] >= 0) update.insert(remap[*it]);}
      curr.swap(update);  
    }
  }
  sid_max = use.size();  
}

void SparseBooleanVariable::set(int iid, int col, bool value) {
  int sid = get_sid(iid);
  if (sid >= 0) {
    if (value) data[col].insert(sid); 
    else data[col].erase(sid);
  }
}

ColumnBlock* SparseBooleanVariable::get_output() {
  return this->is_initialised ? new BlockVariableOutput<SparseBooleanVariable>(*this) : 0;
}

OutputColumn* SparseBooleanVariable::get_output(int col, const string& name) {
  if (this->is_initialised && col < data.size()) {
    if (output_mask) return new SparseBooleanOutput(data[col], output_mask, name);
    else if (map) return new SparseBooleanOutput(data[col], map, name);
    else return new SparseBooleanOutput(data[col], sid_max, name);
  } else return 0;
}


void CompoundOutput::process_data(OutputColumn* source) {OutputBuffer buff(separator);
  if (source) {
    while (!source->is_empty()) buff << *source << endl;
    buff.process_buffer(data, (separator == " ") || (separator == "\t"));
  }
}

void CompoundOutput::set_header(FormattedOutput& fout) {
  if (width == 0) {width = 1;
    for (int i = 0; i < data.size(); i++) {if (data[i].size() > width) width = data[i].size();}
  } 
  fout.add_field(width, name);
}

void StringOutput::set_header(FormattedOutput& fout) {
  if (width == 0) {width = 1;
    while (!iterator->empty()) {
      int len = iterator->next_value().size();
      if (len > width) width = len;
    }
    iterator->reset();
  }    

  fout.add_field(width, name);
  if (!align_right) fout.set_align_right(false);
}

void CStringOutput::set_header(FormattedOutput& fout) {
  if (width == 0) {width = 1;
    while (!iterator->empty()) {
      char *str = iterator->next_value(), len = *(str-1);
      if (len == 127) {
        int cmp_width = strlen(str);
        if (cmp_width > width) width = cmp_width;
      } else if (len > width) width = len;
    }
    iterator->reset();
  }

  fout.add_field(width, name);
  if (!align_right) fout.set_align_right(false);
}

void GeneLocationOutput::set_header(FormattedOutput& fout) {formatted = true;
  int max_chr = 0; unsigned long max_begin = 0, max_end = 0; 
  while (!iterator->empty()) {
    GeneLocation& loc = iterator->next_value();
    if (loc.chromosome > max_chr && loc.chromosome < Utils::chr_reserved) max_chr = loc.chromosome;
    if (loc.begin > max_begin) max_begin = loc.begin;
    if (loc.end > max_end) max_end = loc.end;  
  }
  iterator->reset();
  fout.add_field(max(2, Utils::num_length(max_chr)), "CHR");    
  fout.add_field(max_begin, "START", true);    
  fout.add_field(max_end, "STOP", true);      
}

void GeneLocationOutput::print_value(OutputStream& out, GeneLocation& value) {
  out << Utils::chr_string(value.chromosome, formatted ? 2 : 0); 
  out << value.begin; out << value.end;
}

void EnumOutput::set_header(FormattedOutput& fout) {
  if (width == 0) {width = 1;
    for (int i = 0; i < labels.size(); i++) {
      if (labels[i].size() > width) width = labels[i].size();
    }
  }  
  fout.add_field(width, name);
}


void SparseBooleanOutput::print(OutputStream& out) {
  int sid = map ? map->get(curr_index) - 1 : curr_index;
  while (next_value != values.end() && *next_value < sid) ++next_value; 
  out << (next_value != values.end() && *next_value == sid);  
  curr_index++; find_next();    
}


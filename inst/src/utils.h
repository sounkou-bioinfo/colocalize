/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sstream>
#include <vector>
#include <set>
#include <iostream>
#include <deque>
#include <Eigen/Eigen>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>       

#include "exceptions.h"

#define STRING(s) #s
#define MACRO_STRING(s) STRING(s)
#define UNUSED(x) (void)(x)


using namespace std;
using namespace Eigen;

class Utils {
  template<typename T>
  static bool check_zero_bits() {
    T zero = 0, bits; 
    memset(&bits, 0, sizeof(T));  

    return bits == zero;
  }
         
public:
  static const int chr_reserved = 2000000;
  static const int chrX_code = 2000001;
  static const int chrY_code = 2000002;
  static const int chrXY_code = 2000003;
  static const int chrMT_code = 2000004;
  static const int chrW_code = 2000005;
  static const int chrZ_code = 2000006;
  static const int chr_max = 2000007;

  static short chr_type(int chr) {
    if (chr <= chr_reserved) return 0;
    return chr - chr_reserved;
  }


  static bool is_dir(const string& basename);
  static bool is_file(const string& basename);
  static pair<string,string> split_path(const string& filename);

  template<typename T>
  static bool convert_num_pair(const string& argval, char split, pair<T,T>& target) {
    pair<string,string> bits = split_str(argval, split);
    if (bits.first == "" || bits.second == "") return false;
    return convert_num(bits.first,target.first) && convert_num(bits.second,target.second);
  }

  static bool is_int(const string& argval) {
    long long dump; return convert_num(argval, dump);
  }

  static bool is_num(const string& argval) {
    double dump; return convert_num(argval, dump);
  }

  template<typename T>
  static bool convert_num(const string& argval, T& target) {
    istringstream istr(argval); istr >> target;
    return !istr.fail() && istr.eof();
  }

  static long num_char(const string& str) {string out_str;
    for (int i = 0; i < str.size(); i++) {
      if (str[i] >= '0' && str[i] <= '9') out_str.push_back(str[i]);
    }
    long out_val;
    return convert_num(out_str, out_val) ? out_val : 0;
  }

  static string to_string(string value) {return value;}

  template<typename T>
  static string to_string(T value) {ostringstream ostr; ostr << value; return ostr.str();}
  
//  return static_cast<ostringstream*>( &(ostringstream() << value) )->str();}

  template<typename T>
  static string num_string(T num) {ostringstream ostr; ostr << num; return ostr.str();}  
//  static string num_string(T num) {return static_cast<ostringstream*>( &(ostringstream() << num) )->str();}

  template<typename T>
  static string num_string(T num, unsigned short digits) {
    ostringstream ostr; ostr.precision(digits); 
    ostr << num; return ostr.str();
  }
  
  template<typename T>
  static string num_string_pad(const string& str, T num) {return str + num_string(num);}

  template<typename T>
  static string plural(T num, const string& str, bool print_num=true) {
    string out = print_num ? num_string(num) + " " : "";
    return out + str + (num == 1 ? "" : "s");
  }

  static string ordinal(long num) {string suff;
    switch (abs(num) % 10) {
      case 1: suff = "st"; break;
      case 2: suff = "nd"; break;
      case 3: suff = "rd"; break;            
      default: suff = "th";
    }
    return num_string(num) + suff;
  }
  
  static string num_name(int num);
  
  static int num_length(long num) {return num_string(num).size();}
  
  static string get_letters(int offset, bool uppercase);
  
  static long pad_num(long value, int min) {
    if (min < 10) min = pow(double(10),min);
    while (value < min) value *= 10;
    return value;
  } 
  
  static string& pad_string(string& str, int len, char pad=' ', bool left=true) {
    if (len > str.size()) {
      if (left) str = string(len-str.size(), pad) + str;
      else str += string(len-str.size(), pad);
    }
    return str;
  }
  static string& str_append(string& str, const string& add, const string& sep=" ") {
    if (!str.empty()) str += sep;
    str += add;
    return str;
  }    
  
  static string space(int len) {return string(len, ' ');}
  static string space(const string& ref) {return space(ref.size());}

  template<typename T>
  static string num_quote(T num, char quote='\'') {
    return Utils::quote(num_string(num), quote);
  }

  static string quote(const string& str, char quote='\'') {
    switch (quote) {
      case '(': return "(" + str + ")";
      case '{': return "{" + str + "}";      
      case '[': return "{" + str + "]";            
      default: return string(1,quote) + str + string(1,quote);
    }
  }

  static bool empty_line(const string& line, bool skip_comment=true);
  static pair<string,string> split_str(const string& str, char split);
  static vector<string> tokenize(const string& str);
  static vector<string> tokenize(const string& str, char sep);  

  static string join_string(vector<string>& values, string sep=" ");

  static string lowercase(const string& str);
  static vector<string> lowercase(const vector<string>& vec);
  static string uppercase(const string& str);  
  static vector<string> uppercase(const vector<string>& vec);
  

  static string str_replace(const string& source, const string& from, const string& to);
  static bool starts_with(const string& input, const string& match);
  static bool ends_with(const string& input, const string& match);

  template<typename T>
  static void set_zero(T* buffer, long no_elem) {static bool do_loop = !check_zero_bits<T>();
    if (do_loop) {
      T* end = buffer + no_elem;
      while (buffer < end) *(buffer++) = 0;
    } else memset(buffer, 0, no_elem*sizeof(T));      
  }

  template<typename T>
  static void fill_value(T* buffer, long no_elem, T value) {
    if (value != 0) {
      T* end = buffer + no_elem;
      while (buffer < end) *(buffer++) = value;
    } else set_zero(buffer, no_elem);
  }

  static void set_zero(string* buffer, long no_elem) {fill_value(buffer, no_elem, "");}
  static void fill_value(string* buffer, long no_elem, string value) {
    string* end = buffer + no_elem;
    while (buffer < end) *(buffer++) = value;
  }

  static string chr_string(long chr, int pad=2);
  static bool chr_val(const string& chr_str, int& chr, bool human=true) {
    if (!convert_num(chr_str, chr) || chr >= chr_reserved || (human && chr > 26)) chr = -1;

    if (chr < 0) {
      if (chr_str.size() == 1) {
        if (chr_str[0] == 'x' || chr_str[0] == 'X') chr = chrX_code;
        else if (chr_str[0] == 'y' || chr_str[0] == 'Y') chr = chrY_code;
        else if (!human) {
          if (chr_str[0] == 'z' || chr_str[0] == 'Z') chr = chrZ_code;
          else if (chr_str[0] == 'w' || chr_str[0] == 'W') chr = chrW_code;
        }
      } else if (chr_str.size() == 2) {
        if (chr_str == "xy" || chr_str == "XY") chr = chrXY_code;
        else if (chr_str == "mt" || chr_str == "MT") chr = chrMT_code;
      }
    } else if (human && chr >= 23) {
      switch(chr) {
        case 23: chr = chrX_code; break;
        case 24: chr = chrY_code; break;
        case 25: chr = chrXY_code; break;
        case 26: chr = chrMT_code; break;
      }
    }
    
    return (chr >= 0);
  }
 
  template<typename T>
  static void delete_vector(vector<T*>& ptrs) {
    for (int i = 0; i < ptrs.size(); i++) delete ptrs[i];
    ptrs.clear();
  }
 
  template<typename T>
  static bool vec_equals(const vector<T>& vec1, const vector<T>& vec2) {
    if (vec1.size() != vec2.size()) return false;
    for (int i = 0; i < vec1.size(); i++) {if (vec1[i] != vec2[i]) return false;}
    return true;
  } 
 
  template<typename T>
  static vector<long> indexify(vector<T> vec, int limit=0, bool flip=false) {
    vector<long> out; if (limit <= 0 || limit > vec.size()) limit = vec.size();
    for (long i = 0; i < limit; i++) if (bool(vec[i]) != flip) out.push_back(i);
    return out;
  }

  static bool contains(const vector<string>& vec, const string value) {return find(vec.begin(), vec.end(), value) != vec.end();}
  static bool contains(const string& target, const string value, char sep=' ') {vector<string> words = tokenize(target, sep); return contains(words, value);}
 
  template<typename T1, typename T2>
  static bool filter(vector<T1>& vec, const vector<T2>& select, bool keep=true) {
    if (vec.empty()) return true;
    if (vec.size() != select.size()) return false;
    
    long write = 0;
    if (keep) {
      for (long curr = 0; curr < vec.size(); curr++) {
        if (select[curr]) vec[write++] = vec[curr];
      }
    } else {
      for (long curr = 0; curr < vec.size(); curr++) {
        if (!select[curr]) vec[write++] = vec[curr];
      }
    }
   
    vec.resize(write);
    return true;
  }

  template<typename T1, typename T2>
  static long filter(T1* buff, const vector<T2>& select, bool keep=true) {
    if (select.empty()) return 0;

    T1* write = buff;
    if (keep) {    
      for (int curr = 0; curr < select.size(); curr++) {
        if (select[curr]) *(write++) = buff[curr];
      }
    } else {
      for (int curr = 0; curr < select.size(); curr++) {
        if (!select[curr]) *(write++) = buff[curr];
      }
    }    
  
    return write - buff;
  }


  ///index must be sorted
  template<typename T1, typename T2>
  static long filter_index(T1* buff, long len, const vector<T2>& index, bool keep=true) {
    if (index.empty()) return keep ? 0 : len;

    T1* write = buff; 
    if (keep) {T2 curr, last = -1;
      if (index[0] < 0) return 0;
      for (int i = 0; i < index.size(); i++) {curr = index[i];
        if (curr < last || curr >= len) break;
        if (curr == last) continue;
        *(write++) = buff[curr];
        last = curr;
      }
    } else {       
      T2 from = 0;
      for (int i = 0; i <= index.size(); i++) {
        T2 to = (i == index.size()) ? len : index[i]; 
        if (to == (from - 1)) continue;
        if (to > len || to < from) to = len;

        for (int i = from; i < to; i++) *(write++) = buff[i];
        if (to == len) break;
        from = to + 1;
      }
    }    
  
    return write - buff;
  }
  
  template<typename T>
  static void mem_purge(T& obj) {T temp; obj.swap(temp);}
  
  template<typename T>
  static long filter_value(T* data, long size, T skip, T* buffer=0) {
    if (!buffer) buffer = data;
    T* write = buffer;
    for (T* end = data+size; data < end; data++) {if (*data != skip) *(write++) = *data;}
    return write - buffer;
  }

  template<typename T, typename U>
  static void reorder(vector<T>& vec, const vector<U>& index) {
    vector<T> out; out.reserve(index.size());
    for (long i = 0; i < index.size(); i++) out.push_back(vec[index[i]]);
    vec.swap(out);  
  }

  template<typename T>
  static bool is_discrete(T* vec, long vsize, int max_values) {
    set<T> values;
    for (long i = 0; i < vsize; i++) {
      values.insert(vec[i]);
      if (values.size() > max_values) break;
    }

    return values.size() <= max_values;
  }
};

class Timer {
  clock_t start;
  clock_t last;
  long step_size;
  

public:
  Timer(double step=1) {
    start = clock();
    last = start;
    step_size = step*CLOCKS_PER_SEC;
  }

  bool check(int steps=1) {
    clock_t curr = clock();
    if (steps*(curr-last) > step_size) {last = curr; return true;}
    return false;
  }
};

template<typename T>
class Triple {
  T values[3];

public:
  Triple() {}
  Triple(const T& val) {
    for (int i = 0; i < 3; i++) values[i] = val;
  }
  Triple(const T& v1, const T& v2, const T& v3) {
    values[0] = v1; values[1] = v2; values[2] = v3;
  }

  T& operator[](const int& index) {return values[index];}
  const T& operator[](const int& index) const {return values[index];}
};

template<typename T>
class SharedPtr {
  T* ptr;

public:
  SharedPtr(T* ptr=0) : ptr(ptr) {}
  ~SharedPtr() {delete ptr;}

  T* get() {return ptr;}
  void set(T* input) {delete ptr; ptr = input;}
};


template<typename T1, typename T2>
class IndexSorter {
  vector<pair<T1,T2> >& pairs;
  vector<pair<T1,T2> > pair_proxy;

  static bool asc_sort(const pair<T1,T2> &i, const pair<T1,T2> &j) {return i.second < j.second;}
  static bool desc_sort(const pair<T1,T2> &i, const pair<T1,T2> &j) {return j.second < i.second;}

  void run(bool ascending) {
    if (ascending) sort(pairs.begin(), pairs.end(), asc_sort);
    else sort(pairs.begin(), pairs.end(), desc_sort);  
  }

public:
  IndexSorter(vector<pair<T1,T2> >& pairs, bool ascending) : pairs(pairs) {run(ascending);}
  IndexSorter(vector<T2>& values, bool ascending) : pairs(pair_proxy) {
    for (int i = 0; i < values.size(); i++) pairs.push_back(pair<T1,T2>(i,values[i]));    
    run(ascending);
  }
  
  vector<T1> get_index() {vector<T1> out; out.reserve(pairs.size());
    for (int i = 0; i < pairs.size(); i++) out.push_back(pairs[i].first);
    return out;
  }
};

template<typename T>
class GenericSorter {
  struct AscSort {T* obj; AscSort(T* obj) : obj(obj) {} bool operator () (const int &i, const int &j) {return obj[i] < obj[j];} };
  struct DescSort {T* obj; DescSort(T* obj) : obj(obj) {} bool operator () (const int &i, const int &j) {return obj[j] < obj[i];} };

public:
  vector<int> run(T* objects, int size, bool ascending=true) {
    vector<int> ids; ids.reserve(size);
    for (int i = 0; i < size; i++) ids.push_back(i);
    if (ascending) sort(ids.begin(), ids.end(), AscSort(objects));
    else sort(ids.begin(), ids.end(), DescSort(objects));
    return ids;
  }
};


template<typename T>
class BuffSorter {
  const T* objects;
  vector<long> index;
  
  struct SortObj {
    BuffSorter &sorter;
    SortObj(BuffSorter &bs) : sorter(bs) {}
    bool operator () (const long &i, const long &j) {return sorter.objects[j] < sorter.objects[i];}
  };

public:
  vector<long>& run(const T* obj, long size) {
    objects = obj;
    if (index.size() != size) index.resize(size);
    for (long i = 0; i < size; i++) index[i] = i;    
    sort(index.begin(), index.end(), SortObj(*this));
    return index;
  }
};

template<typename T>
class Persistent {
  T* content; int* count;

public:  
  Persistent(T* content, int* cnt=0) : content(content) { 
    if (cnt) {count = cnt; (*count)++;}
    else {count = new int(1);}
  }
  ~Persistent() {
    (*count)--;
    if (*count <= 0) {delete content; delete count;}
  }
  
  T* get_content() {return content;}
  Persistent<T>* copy() {return new Persistent<T>(content, count);}
};
                 

template<typename T> class BufferWindow;

template<typename T>
class Buffer {
  friend class BufferWindow<T>;

protected:
  T* content;
  long rows;
  long cols;
  long capacity;
  bool square;
  bool owner;

  Buffer<T>& operator=(const Buffer<T>& other) {assign(other);} ///Block accidental assignment
  void delete_content() {
    if (owner) delete[] content;
    content = 0;
  }

public:
  Buffer() : content(0), rows(0), cols(0), capacity(0), square(false), owner(true) {}
  Buffer(T* buff, long r, long c=1) : content(0), rows(0), cols(0), capacity(0), square(false), owner(true) {assign(buff, r, c);}
  Buffer(Buffer<T>& buff) : content(0), rows(0), cols(0), capacity(0), square(false), owner(true) {assign(buff);}
  Buffer(long r, long c=1, bool zero=true) : content(0), rows(0), cols(0), capacity(0), square(false), owner(true) {set_size(r, c, zero);}
  ~Buffer() {delete_content();}

  void set_empty() {rows = cols = 0;}
  void set_size(long r, long c=0, bool zero=false) {
    if (r*c > capacity) {
      if (!owner) throw BufferException("capacity of linked buffer exceeded");
      rows = r; cols = c; capacity = rows*cols;
      delete_content();

      content = new T[capacity];
      if (zero) Utils::set_zero(content, capacity);
    } else {
      rows = r; cols = c;
      if (zero && capacity > 0) Utils::set_zero(content, rows*cols);
    }
    square = square && (rows == cols);
  }

 void set_square(long rc=0, bool zero=false) {
    square = true;
    set_size(rc, rc, zero);
  }

  void make_owner(long c_add=0) {
    if (owner) return;
    if (c_add) square = false;
    capacity = rows*(cols+c_add);
    T* new_content = new T[capacity];
    memcpy(new_content, content, sizeof(T)*rows*cols);
    content = new_content;
    owner = true;
  }

  void total_cols(long required, bool zero=false) {
    if (required > cols) add_cols(required-cols, zero);
    else cols = max(required, long(1));
  }
  void total_cols(long required, bool init, T value) {
    if (required > cols) add_cols(required-cols, init, value);
    else cols = max(required, long(1));
  }

  void add_cols(long c_add, bool zero=false) {
    if (c_add <= 0) return;
    if (!owner) make_owner(c_add);

    square = false;
    if (cols > 0) {
      if ((c_add+cols)*rows > capacity) {
        capacity = (c_add+cols)*rows;
        T* new_content = new T[capacity];
        memcpy(new_content, content, sizeof(T)*rows*cols);
        delete_content();
        content = new_content;
      }
      if (zero) Utils::set_zero(content+rows*cols, rows*c_add);
      cols += c_add;
    } else {
      set_size(rows, c_add, zero);
    }
  }
  
  void add_cols(long c_add, bool init, T value) {
    if (c_add <= 0) return;
    if (!owner) make_owner(c_add);    

    long offset = cols;
    add_cols(c_add);
    if (init) for (T *curr = content+rows*offset, *end = curr+rows*c_add; curr < end; curr++) *curr = value;
  }

  bool append(Buffer<T>& add) {
    if (add.rows != rows) return false;
    long offset = cols;
    add_cols(add.cols);    
    memcpy(content+rows*offset, add.content, sizeof(T)*add.rows*add.cols);    
    return true;
  }

  void shift_cols(long shift, bool zero=false) {
    if (shift == 0) return;
    square = false;
    if (shift > 0) {
      if (cols > 0) {
        if (!owner) make_owner(shift);
        if ((shift+cols)*rows > capacity) {
          capacity = (shift+cols)*rows;
          T* new_content = new T[capacity];
          memcpy(new_content+shift*rows, content, sizeof(T)*rows*cols);
          delete_content();
          content = new_content;
        } else {
          memmove(content+shift*rows, content, sizeof(T)*rows*cols);
        }
        if (zero) Utils::set_zero(content, rows*shift);
        cols += shift;
      } else {
        set_size(rows, shift, zero);
      }
    } else {
      shift = -shift;
      if (shift < cols) {
        if (owner) memmove(content, content+shift*rows, sizeof(T)*rows*(cols-shift));
        else content += shift*rows;
        cols -= shift;
      } else {
        set_size(rows, 0);
      }
    }
  }

  T* expand(long required, long margin=0) {
    if (margin && required > cols) required = max(required, cols+margin);

    if (square) set_size(required, required, false);
    else set_size(rows, required, false);

    Utils::set_zero(content, cols*rows);
    return content;
  }
  
  void resize(int c) { ///resize to specified number of columns, retaining content (unless in discarded columns)
    if (c > cols) add_cols(c - cols);
    else shrink(c);  
  }
  
  ///shrinks to specified value (if smaller than current), retains content in remaining columns
  int shrink(int c) {
    if (c < cols) cols = c;
    square = false;
    return cols;
  }
  
  int shrink_rows(int r) {
    if (r < rows) {
      for (int c = 1; c < cols; c++) memmove(content+c*r, content+c*rows, sizeof(T)*r);
      rows = r; square = false;
    }
    return rows;
  }
  
  void fit_rows(int r) {fit_rows(r, 0, false);}
  void fit_rows(int r, T default_value, bool has_default=true) {
    if (r <= 0) {reset(); return;} if (r == rows) return; 

    capacity = r*cols;   
    T* new_content = new T[capacity]; int move = min((long) r, rows);
    for (int c = 0; c < cols; c++) memcpy(new_content+c*r, content+c*rows, sizeof(T)*move);
    if (r > rows && has_default) {for (int c = 0; c < cols; c++) Utils::fill_value(new_content+c*r+rows, r-rows, default_value);}
    delete[] content; content = new_content; rows = r;   
  }

  template <typename U>
  void drop_rows(const vector<U>& drop, bool flip=false, bool shrink=true) {drop_rows_by_index(Utils::indexify(drop, rows, flip), shrink);}

  template <typename U>
  void drop_rows_by_index(const vector<U>& drop, bool shrink=true) {
    if (drop.size() == 0 || empty()) return;
    if (shrink && !owner) make_owner();

    long nrow = shrink ? rows - drop.size() : rows;
    T* new_buff = shrink ? new T[nrow*cols] : content;

    long from = 0, to, offset = 0;
    for (long i = 0; i <= drop.size(); i++) {
      to = (i == drop.size()) ? rows : drop[i];
      if (from < to) {
        for (long c = 0; c < cols; c++) {
          memmove(new_buff + c*nrow + offset, content + c*rows + from, sizeof(T)*(to-from));
        }
        offset += (to-from);
      }
      from = to + 1;
    }

    rows = nrow;
    square = false;
    if (shrink) {
      delete[] content;
      content = new_buff;
      capacity = rows*cols;
    }
  }

  void drop_col(long drop) {
    vector<long> index; index.push_back(drop);
    drop_cols_by_index(index);
  }

  template <typename U>
  void drop_cols(const vector<U>& drop) {drop_cols_by_index(Utils::indexify(drop, cols, false));} 
  
  template <typename U>
  void drop_cols_by_index(const vector<U>& drop) {
    if (drop.size() == 0 || empty()) return;

    long no_dropped = drop.size();  
    for (long i = 0; i < no_dropped; i++) {
      long shift = i+1, from = drop[i]+1, to = (i+1) < no_dropped ? drop[i+1] : cols;
      for (long j = from; j < to; j++) memcpy(content+(j-shift)*rows, content+j*rows, sizeof(T)*rows);
    }

    cols -= no_dropped;
    square = false;
  }
  
  void assign_zero() {Utils::set_zero(content, cols*rows);}
  void assign_zero(long offset, long total=0) { 
    if (offset >= cols) return;
    if (offset < 0) offset = 0;
    total = (total > 0) ? min(total, cols - offset) : cols - offset;
    Utils::set_zero(content+rows*offset, total*rows);    
  }
  void assign_zero_col(long c) {assign_zero(c,1);} 

  void assign_value(T value) {assign_value(value, 0, cols);}
  void assign_value(T value, long offset, long total=0) { 
    if (empty() || offset >= cols) return;
    if (total == 0) total = cols;
    
    T* base = content+offset*rows;    
    for (long r = 0; r < rows; r++) base[r] = value;
    for (long c = offset+1; c < min(cols, total+offset); c++) memcpy(content+c*rows, base, sizeof(T)*rows);
  }
  void assign_value_col(T value, long c) {assign_value(value, c, 1);}  
  void assign_col(T* source, long c) {memcpy(content+c*rows, source, sizeof(T)*rows);}    

  void assign(BufferWindow<T> source) {assign(source.view());}
  void assign(const Buffer<T>& source) {if (this != &source) assign(source.content, source.rows, source.cols, source.square);}
  void assign(const Buffer<T>& source, long offset, long total=0) {
    if (this == &source) return;
    if (offset >= source.cols) offset = 0;
    total = total > 0 ? min(source.cols-offset, total) : source.cols-offset;      
    assign(source.content + offset*source.rows, source.rows, total);      
  }
  void assign(const T* source, long r, long c=1, bool square=false) {
    if (square && r == c) set_square(r, false);
    else set_size(r, c, false);
    memcpy(content, source, sizeof(T)*r*c);
  }

  void swap(Buffer<T>& other) {
    if (this == &other) return;
    std::swap(content, other.content);
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
    std::swap(capacity, other.capacity);
    std::swap(square, other.square);
    std::swap(owner, other.owner);
  }

  void unlink() {if (!owner) reset();}
  void reset(bool unsquare=false) {
    delete_content();
    rows = 0; cols = 0; capacity = 0;
    if (unsquare) square = false;
    owner = true;
  }

  T* begin() {return content;}
  T* end() {return content+cols*rows;}

  bool is_square() {return square;}
  bool is_owner() {return owner;}
  bool empty() {return rows == 0 || cols == 0;}
  long nrow() {return rows;}
  long ncol() {return cols;}
  long length() {return rows*cols;}
  long maxcol() {return rows > 0 ? capacity / rows : 0;}
  long maxsize() {return capacity;}
  bool check_size(long d) {return d == rows && (!square || d == cols);}
  bool check_size(long r, long c) {return r==rows && c == cols;}
  int type_size() {return sizeof(T);}

  Map<Matrix<T,Dynamic,Dynamic> > matrix() {return Map<Matrix<T,Dynamic,Dynamic> >(content, rows, cols);}
  Map<Matrix<T,Dynamic,Dynamic> > matrix(long offset, long total=0) {
    if (offset >= cols) throw BufferException("column offset exceeds buffer size");
    total = total > 0 ? min(cols-offset, total) : cols-offset;
    return Map<Matrix<T,Dynamic,Dynamic> >(content+rows*offset, rows, total);
  }

  Map<Matrix<T,Dynamic,Dynamic>, 0, OuterStride<> > matrix_colskip(long offset, long total=0, long step=1) {
    if (offset >= cols) throw BufferException("column offset exceeds buffer size");
    long max_steps = ceil((cols-offset)/float(step));
    total = total > 0 ? min(max_steps, total) : max_steps;
    return Map<Matrix<T,Dynamic,Dynamic>, 0, OuterStride<> >(content+rows*offset, rows, total, OuterStride<>(rows*step));
  }
  
  Map<Matrix<T,Dynamic,Dynamic>, 0, Stride<Dynamic,Dynamic> > matrix_rowskip(long offset, long total, long step=1) {
    if (offset >= rows) throw BufferException("row offset exceeds buffer size");
    long max_steps = ceil((rows-offset)/float(step));
    total = total > 0 ? min(max_steps, total) : max_steps;
    return Map<Matrix<T,Dynamic,Dynamic>, 0, Stride<Dynamic,Dynamic> >(content+offset, total, cols, Stride<Dynamic,Dynamic>(rows,step));
  }

  Map<Matrix<T,Dynamic,Dynamic>, 0, OuterStride<> > matrix_sub(long offset, long total) {return matrix_sub(offset, total, offset, total);}
  Map<Matrix<T,Dynamic,Dynamic>, 0, OuterStride<> > matrix_sub(long row_offset, long row_total, long col_offset, long col_total) {
    if (row_offset >= rows) throw BufferException("row offset exceeds buffer size");
    if (col_offset >= cols) throw BufferException("column offset exceeds buffer size");
    if (row_offset + row_total > rows) row_total = rows - row_offset;    
    if (col_offset + col_total > cols) col_total = cols - col_offset;

    return Map<Matrix<T,Dynamic,Dynamic>, 0, OuterStride<> >(content+row_offset+rows*col_offset, row_total, col_total, OuterStride<>(rows));
  }
  
  Map<Matrix<T,Dynamic,1> > matrix_vec(long col=0) {
    if (col >= cols) throw BufferException("column offset exceeds buffer size");
    return Map<Matrix<T,Dynamic,1> >(content+rows*col, rows, 1);
  }
  Map<Matrix<T,Dynamic,1> > matrix_vec(long col, long offset, long total=0) {
    if (col >= cols) throw BufferException("column offset exceeds buffer size");
    if (offset >= rows) throw BufferException("column offset exceeds buffer size");
    total = total > 0 ? min(rows-offset, total) : rows-offset;
    return Map<Matrix<T,Dynamic,1> >(content+rows*col+offset, total, 1);
  }
  
  BufferWindow<T> window(int offset=0, int size=0, bool is_vector=false) {return BufferWindow<T>(*this, offset, size, is_vector);};

  T* operator[](const long& c) {return content+(c<0 ? cols+c : c)*rows;}
  const T* operator[](const long& c) const {return content+(c<0 ? cols+c : c)*rows;}
  T& operator()(const long& index) {return content[index];}
  T& operator()(const long& r, const long& c) {return content[c*rows + r];}
  T& operator()(const long& r, const long& c) const {return content[c*rows + r];}

  template<typename U>
  T* operator+(const U& add) {return content + add;}

  friend ostream& operator<<(ostream& os, const Buffer& buff) {
    for (long r = 0; r < buff.rows; r++) {
      for (long c = 0; c < buff.cols; c++) {
        if (c > 0) os << '\t';
        os << buff(r,c);
      }
      os << endl;
    }
    return os;
  }
};

template<typename T>
class BufferWindow {
  Buffer<T>* parent;
  Buffer<T>* window;
  
  void link(T* source, long r, long c) {
    if (!window) window = new Buffer<T>; 
    window->content = source;
    window->rows = r; window->cols = c; window->capacity = r*c;
    window->square = false;
    window->owner = false;
  }  
  
  void copy(const BufferWindow<T>& other) {
    parent = other.parent; 
    if (other.window) {
      window = new Buffer<T>;    
      link(other.window->content, other.window->rows, other.window->cols);
    } else window = 0;
  }

public:
  BufferWindow() : parent(0), window(0) {}
  BufferWindow(Buffer<T>& parent, int offset=0, int size=0, bool is_vector=false) : parent(0), window(0) {link_buffer(parent, offset, size, is_vector);}
  BufferWindow(T* data, int nrow, int ncol=1) : parent(0), window(0) {link_pointer(data, nrow, ncol);}
  BufferWindow(const BufferWindow<T>& other) {copy(other);}
  ~BufferWindow() {if (window) delete window;}
  
  BufferWindow& operator=(const BufferWindow<T>& other) { 
    if (this != &other) {
      if (window) delete window;
      copy(other);
    }
    return *this;
  }
  
  Buffer<T>& view() {
    if (!window) throw BufferException("linked buffer has expired");
    return *window;
  }

  void link_buffer(Buffer<T>& source, int offset=0, int size=0, bool is_vector=false) {
    parent = &source;
    if (is_vector) {
      if (offset >= parent->rows) offset = 0;
      size = size > 0 ? min(parent->rows-offset, long(size)) : parent->rows-offset;      
      link(parent->content + offset, size, 1);      
    } else {
      if (offset >= parent->cols) offset = 0;
      size = size > 0 ? min(parent->cols-offset, long(size)) : parent->cols-offset;      
      link(parent->content + offset*parent->rows, parent->rows, size);      
    }
  }
  
  void link_pointer(T* source, long r, long c=1) {
    parent = 0;
    link(source, r, c);
  }
  
  void unlink() {
    if (window) {
      delete window;
      window = 0;
    }
  }
};

        
#endif /** UTILS_H*/

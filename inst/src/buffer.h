/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef BUFFER_H
#define BUFFER_H

#include "utils.h"

template<typename T>
class BaseBuffer {
protected:
  T* content;
  long capacity;
  long no_elem;

  void do_copy(T* to, T* from, int total);
  void set_size(long new_size);
 
  BaseBuffer<T>& operator=(const BaseBuffer<T>& other) {assign(other);} ///Block accidental assignment
  BaseBuffer(BaseBuffer<T>& source) : content(0) {assign(source.content, source.capacity); no_elem = source.no_elem;}
  
public: 
  BaseBuffer(long init_size=0, bool zero=false) : content(0), capacity(0), no_elem(0) {set_size(init_size); if (zero) set_zero();}
  virtual ~BaseBuffer() {clear();}

  void clear() {delete[] content; content = 0; no_elem = 0; capacity = 0;}
  void resize(long new_size, bool shrink_to_fit=false);
  void set_zero() {if (content) Utils::set_zero(content, no_elem);}
  
  void append(long add) {if (add > 0) resize(no_elem + add);}  
  void reserve(long total);

  void assign(BaseBuffer<T>& source) {assign(source.content, source.no_elem);}
  void assign(T* source, long len);
  void assign_value(T value, int amount=-1);

  void swap(BaseBuffer<T>& other);

  bool is_null() {return content == 0;}
  bool is_empty() {return content == 0 || no_elem == 0;}
  long size(bool avail=false) {return avail ? capacity : no_elem;}
  long remaining() {return capacity - no_elem;}
  T* data() {return content;}

  T get_value(const long& index) const {return content[index];}
  virtual T& get(const long& index) {return content[index];}  
  virtual T& operator[](const long& index) {return content[index];}  
    
  friend ostream& operator<<(ostream& os, const BaseBuffer& buff) {
    for (int i = 0; i < buff.no_elem; i++) {os << buff.get_value(i) << endl;}
    return os;
  }
};

template<typename T>
class ExpandingBuffer : public BaseBuffer<T> {typedef BaseBuffer<T> PARENT;
  const double expansion_factor;
  void expand(long min_elem=0);
 
public:
  ExpandingBuffer(long init_size=0, double factor=1.1) : BaseBuffer<T>(init_size), expansion_factor(factor) {}

  void expansion_reserve(long size);
  virtual T& get(const long& index);
  virtual T& operator[](const long& index) {return get(index);}
};

class ArchiveBuffer {
protected:
  static const double EXPANSION_FACTOR;
  static const int EXPANSION_MINIMUM;  

  vector<char*> content;
  char* active;

  long capacity;  
  long no_elem;
  long total_used;
  
  int elem_size;
  int min_elem;
  int max_elem;
  
  virtual void expand(long minimum=0);
  virtual char* add_internal(const char* str, int len);
  
public: 
  ArchiveBuffer() : active(0), capacity(0), no_elem(0), total_used(0), elem_size(10), min_elem(1000) {}
  virtual ~ArchiveBuffer() {clear();}

  virtual void init(long amount=0, int exp_size=0);
  virtual void clear();

  long size(bool total=true) {return total ? no_elem : elem_size;}
  bool initialised() {return active;}

  char* add(const string& str) {return add_internal(str.data(), str.length());}
  char* add(char* str, bool from_archive=false);    
};

class IndexedArchiveBuffer : public ArchiveBuffer {
  BaseBuffer<char*> index;
  
  void expand(long minimum=0) {ArchiveBuffer::expand(minimum); index.resize(max_elem);}
  char* add_internal(const char* str, int len);
  
public: 
  void clear() {ArchiveBuffer::clear(); index.clear();}

  char* operator[](const long& i) {return index[i];}
  char* get_last() {return index[no_elem-1];}  
  
  BaseBuffer<char*>& get_index() {index.resize(no_elem); return index;}
  void unload_index(BaseBuffer<char*>& target) {index.resize(no_elem); target.swap(index); index.clear();}
  void clear_index() {index.clear();}
};


///Assumes unsigned char = 8 bits
typedef unsigned char bitBlock;

class BitBuffer {
  static bool global_init;
  static bitBlock write_mask[8];
  static bitBlock count_mask[256];  

  bitBlock* content; 
  long rows;
  long cols;
  long row_blocks;
  long capacity;
  
  deque<int> traverse_index;

  void reset(bool init=false) {
    if (init && !global_init) mask_init();
    if (!init) delete[] content;
    content = 0;
    rows = 0; cols = 0;
    row_blocks = 0; 
    capacity = 0;
  }
  
  void mask_init();

  int count_internal(bitBlock* c);
  int both_internal(bitBlock* c1, bitBlock* c2);  
  int either_internal(bitBlock* c1, bitBlock* c2);    

public:
  BitBuffer() {reset(true);}
  BitBuffer(int r, int c=1, bool zero=true) {reset(true); set_size(r,c,zero);}
  ~BitBuffer() {delete[] content;}

  void set_size(int r, int c=0, bool zero=true) {
    rows = r; cols = c; row_blocks = rows / 8 + (rows % 8 > 0); 

    if (row_blocks*cols > capacity) {
      capacity = row_blocks*cols;
      delete[] content;
      content = new bitBlock[capacity];
    }
    if (zero && capacity > 0) Utils::set_zero(content, row_blocks*cols);
  }

  bool get(int r, int c) {return content[r/8 + c*row_blocks] & write_mask[r % 8];}
  deque<int>& traverse(int col);

  void set(int r, int c=0) {content[r/8 + c*row_blocks] |= write_mask[r % 8];}
  void unset(int r, int c=0) {content[r/8 + c*row_blocks] &= ~write_mask[r % 8];}  

  int count(int c) {return count_internal(content+c*row_blocks);}
  int both(int c1, int c2) {return both_internal(content+c1*row_blocks, content+c2*row_blocks);}
  int both(int cs, int co, BitBuffer& other) {return both_internal(content+cs*row_blocks, other.content+co*other.row_blocks);}
  int either(int c1, int c2) {return either_internal(content+c1*row_blocks, content+c2*row_blocks);}
  int either(int cs, int co, BitBuffer& other) {return either_internal(content+cs*row_blocks, other.content+co*other.row_blocks);}

  template <typename T>
  void drop_cols_by_index(const vector<T>& drop) {
    if (drop.size() == 0) return;

    int no_dropped = drop.size();  
    for (int i = 0; i < no_dropped; i++) {
      long shift = i+1, from = drop[i]+1, to = (i+1) < no_dropped ? drop[i+1] : cols;
      for (long j = from; j < to; j++) memcpy(content+(j-shift)*row_blocks, content+j*row_blocks, sizeof(bitBlock)*row_blocks);
    }
    cols -= no_dropped;
  }

  friend ostream& operator<<(ostream& os, const BitBuffer& buff) {
    for (int r = 0; r < buff.row_blocks; r++) {
      for (int c = 0; c < buff.cols; c++) {
        if (c > 0) os << '\t';
        os << int(buff.content[r + c*buff.row_blocks]);
      }
      os << endl;
    }
    
    return os;
  }
};

#include "buffer.tpp"



#endif /** BUFFER_H */

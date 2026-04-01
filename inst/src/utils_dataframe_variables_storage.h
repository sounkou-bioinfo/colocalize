/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef UTILS_DATAFRAME_VARIABLES_STORAGE_H
#define UTILS_DATAFRAME_VARIABLES_STORAGE_H

#include <string>

#include "buffer.h"

using namespace std;

template<typename T> class VariableStorage; 
template<typename T> class VariableUniStorage; 
template<typename T> class VariableMultiStorage; 

namespace VariableStorageAccess { 
  enum FilterMode {Smaller, SmallerEqual, Equal, GreaterEqual, Greater};
  
  template<typename T>
  struct Accessor {
    virtual ~Accessor() {}
    virtual bool exists(int iid) {UNUSED(iid); return true;}
    virtual int get_sid(int iid) {return iid;}
  };

  template<typename T>
  struct UniAccessor : public Accessor<T> {
    T* buffer; 
    virtual ~UniAccessor() {}
    virtual T& elem(int iid) = 0;
    virtual void set_buffer(BaseBuffer<T>& data) {buffer = data.data();}
  };

  template<typename T>  
  struct DummyAccess : public UniAccessor<T> {
    T dummy_value;
    T& elem(int iid) {UNUSED(iid); return dummy_value;}
    bool exists(int iid) {UNUSED(iid); return false;}    
  };
  
  template<typename T>
  struct DirectAccess : public UniAccessor<T> {
    DirectAccess(BaseBuffer<T>& data) {this->set_buffer(data);}
    T& elem(int iid) {return this->buffer[iid];}
  };

  template<typename T>  
  struct MappedAccess : public UniAccessor<T> {
    BaseBuffer<int>& map;
    void (*on_remove)(T&);
    bool has_default;
    T default_value; T dummy;

    MappedAccess(BaseBuffer<T>& data, BaseBuffer<int>* map) : map(*map), on_remove(0), has_default(false) {this->set_buffer(data);}          
    MappedAccess(BaseBuffer<T>& data, BaseBuffer<int>* map, T default_value, void (*delete_func)(T&)=0) : map(*map), on_remove(delete_func), has_default(true), default_value(default_value), dummy(default_value) {this->set_buffer(data);}

    T& elem(int iid) {int lid = map[iid];
      if (!lid) {
        if (on_remove) on_remove(dummy);
        if (has_default) dummy = default_value; 
        return dummy;
      } else return this->buffer[lid-1];
    }
    bool exists(int iid) {return map[iid];}
    int get_sid(int iid) {return map[iid] - 1;}
  };

  template<typename T>  
  struct ExpandSetter : public UniAccessor<T> {
    VariableUniStorage<T>& owner;
    UniAccessor<T>* setter;
    double expansion_factor;
    int capacity;
 
    ExpandSetter(VariableUniStorage<T>& owner, UniAccessor<T>* subsetter) : owner(owner), setter(subsetter), expansion_factor(0.5) {init();}
    ~ExpandSetter() {delete setter;}
    
    T& elem(int iid) {if (iid >= capacity) expand(iid); return setter->elem(iid);}
    bool exists(int iid) {return setter->exists(iid);}
    void set_buffer(BaseBuffer<T>& data) {update_capacity(); setter->set_buffer(data);}        
   
    void init();
    void update_capacity();    
    void expand(int min_iid);
  };

  template<typename T>  
  struct DeleteSetter : public UniAccessor<T> {
    void (*on_remove)(T&);
    UniAccessor<T>* setter;

    DeleteSetter(void (*delete_func)(T&), UniAccessor<T>* subsetter) : on_remove(delete_func), setter(subsetter) {}     
    ~DeleteSetter() {delete setter;}

    T& elem(int iid) {T& val = setter->elem(iid); on_remove(val); return val;}      
    bool exists(int iid) {return setter->exists(iid);}
    void set_buffer(BaseBuffer<T>& data) {UniAccessor<T>::set_buffer(data); setter->set_buffer(data);}        
  };

  template<typename T> 
  struct ElementIterator {
    T* current; T* end; int step_size;
    void set(T* start, int no_elem, int step=1) {current = start; end = start + no_elem*step; step_size = step;}
    bool empty() {return current >= end;}
    T& elem() {current += step_size; return *(current-step_size);}
    int size() {return (end-current) / step_size;}
  };

  template<typename T>
  struct MultiAccessElement {
    ElementIterator<T> iterator;
    virtual ~MultiAccessElement() {}

    virtual T& get(int col) = 0;
    virtual void set(int col, T value) {get(col) = value;}
    virtual int size() = 0;   

    virtual ElementIterator<T>& get_iterator(int lid=-1) = 0;
    virtual void set_row(int lid) = 0;   
    virtual int get_row() = 0; 
  };
  
  template<typename T>
  struct MultiAccessDummy : public MultiAccessElement<T> {
    T dummy_value;
    T& get(int col) {UNUSED(col); return dummy_value;}
    virtual int size() {return 0;}

    ElementIterator<T>& get_iterator(int lid=-1) {UNUSED(lid); this->iterator.set(0,0); return this->iterator;}      
    void set_row(int lid) {UNUSED(lid);} 
    int get_row() {return -1;}
  };
  
  template<typename T>
  struct MultiAccessBlock : public MultiAccessElement<T> {
    T *base, *current; int ncol; int step;
    MultiAccessBlock(Buffer<T>& data) : base(data.begin()), ncol(data.ncol()), step(data.nrow()) {set_row(0);}
    
    T& get(int col) {return *(current + col*step);}
    int size() {return ncol;}
    
    ElementIterator<T>& get_iterator(int lid=-1) {this->iterator.set(lid >= 0 ? base+lid : current, ncol, step); return this->iterator;}
    void set_row(int lid) {current = base + lid;}    
    int get_row() {return current - base;}
  };
  
  template<typename T>
  struct MultiAccessJagged : public MultiAccessElement<T> {
    BaseBuffer<BaseBuffer<T>*>& base; T* curr_data; int curr_row;
    MultiAccessJagged(BaseBuffer<BaseBuffer<T>*>& data) : base(data) {set_row(0);}

    bool valid_lid(int lid) {return lid < base.size() && base[lid];}
    
    T& get(int col) {return curr_data[col];}
    int size() {return curr_data ? base[curr_row]->size() : 0;}
    
    ElementIterator<T>& get_iterator(int lid=-1) {if (lid < 0) lid = curr_row;
      if (valid_lid(lid)) this->iterator.set(base[lid]->data(), base[lid]->size()); 
      else this->iterator.set(0,0);
      return this->iterator;
    }
    void set_row(int lid) {curr_row = lid; curr_data = valid_lid(lid) ? base[lid]->data() : 0;}    
    int get_row() {return curr_row;}
  };
  
  template<typename T>
  struct MultiAccessDelete : public MultiAccessElement<T> {
    MultiAccessElement<T>* child;
    void (*on_remove)(T&);    
    MultiAccessDelete(MultiAccessElement<T>* child, void (*delete_func)(T&)) : child(child), on_remove(delete_func) {}  
    ~MultiAccessDelete() {delete child;}
    
    T& get(int col) {return child->get(col);}
    void set(int col, T value) {T& val = child->get(col); on_remove(val); val = value;}    
    int size() {return child->size();}

    ElementIterator<T>& get_iterator(int lid=-1) {return child->get_iterator(lid);}
    void set_row(int lid) {child->set_row(lid);}
    int get_row() {return child->get_row();}
  };
  
  template<typename T>
  struct MultiAccessor : public Accessor<T> {
    MultiAccessElement<T>* access;
    MultiAccessor() : access(new MultiAccessDummy<T>()) {}
    virtual ~MultiAccessor() {delete access;}
    virtual void set_access(MultiAccessElement<T>* update, void (*delete_func)(T&)=0) {delete access; access = delete_func ? new MultiAccessDelete<T>(update, delete_func) : update;}
    
    virtual MultiAccessElement<T>& elem(int iid) = 0;
    virtual ElementIterator<T>& iterator(int iid) = 0 ;
  };
  
  template<typename T>
  struct DirectMultiAccess : public MultiAccessor<T> {
    MultiAccessElement<T>& elem(int iid) {this->access->set_row(iid); return *(this->access);}
    ElementIterator<T>& iterator(int iid) {return this->access->get_iterator(iid);}
  };
    
  template<typename T>
  struct MappedMultiAccess : public MultiAccessor<T> {
    BaseBuffer<int>& map; 
    MultiAccessElement<T>* dummy;      
    bool has_default; T default_value;    
    
    MappedMultiAccess(BaseBuffer<int>* map) : map(*map), dummy(new MultiAccessDummy<T>()), has_default(false) {}
    MappedMultiAccess(BaseBuffer<int>* map, T default_value) : map(*map), dummy(new MultiAccessDummy<T>()), has_default(true), default_value(default_value) {}
    ~MappedMultiAccess() {delete dummy;}

    void set_access(MultiAccessElement<T>* update, void (*delete_func)(T&)=0) {
      MultiAccessor<T>::set_access(update, delete_func);
      delete dummy; if (delete_func) dummy = new MultiAccessDelete<T>(new MultiAccessDummy<T>(), delete_func); else dummy = new MultiAccessDummy<T>();
    }
    MultiAccessElement<T>& elem(int iid) {int lid = map[iid];
      if (!lid) {
        if (has_default) dummy->set(0, default_value);
        return *dummy;
      } else {this->access->set_row(lid-1); return *(this->access);}
    }
    ElementIterator<T>& iterator(int iid) {int lid = map[iid];
      return lid ? this->access->get_iterator(lid-1) : dummy->get_iterator(0);    
    }
    int get_sid(int iid) {return map[iid] - 1;}
  };
  
  template<typename T>  
  struct ExpandingMultiAccess : public MultiAccessor<T> {
    VariableMultiStorage<T>& owner;
    MultiAccessor<T>* setter;
    double expansion_factor;
    int capacity;
 
    ExpandingMultiAccess(VariableMultiStorage<T>& owner, MultiAccessor<T>* subsetter) : owner(owner), setter(subsetter), expansion_factor(0.5) {update_capacity();}
    ~ExpandingMultiAccess() {delete setter;}
    
    void set_access(MultiAccessElement<T>* update, void (*delete_func)(T&)=0) {setter->set_access(update, delete_func);}
    MultiAccessElement<T>& elem(int iid) {if (iid >= capacity) expand(iid); return setter->elem(iid);}
    ElementIterator<T>& iterator(int iid) {if (iid >= capacity) expand(iid); return setter->iterator(iid);}
    
    void update_capacity();
    void expand(int min_iid);
  };

 
  template<typename T>
  class AccessIterator {
  public: 
    virtual ~AccessIterator() {}   
    virtual T& next_value() = 0;
    virtual pair<int,T> next_element() = 0;    
    virtual bool empty() = 0;    
    virtual void reset() = 0;
  };
  
  template<typename T>
  class DummyAccessIterator : public AccessIterator<T> {
    T dummy;
  public:
    DummyAccessIterator() {}
    DummyAccessIterator(T default_value) : dummy(default_value) {}
    T& next_value() {return dummy;}
    pair<int,T> next_element() {return pair<int,T>(-1,dummy);}
    bool empty() {return true;}      
    void reset() {}
  };
  
  template<typename T>
  class DirectAccessIterator : public AccessIterator<T> {
    T* source; T* current; T* end;
  public:
    DirectAccessIterator(T* data, int size) : source(data), current(data), end(data+size) {}
    T& next_value() {return *(current++);}
    pair<int,T> next_element() {int sid = (current++ - source); return pair<int,T>(sid,source[sid]);}    
    bool empty() {return current >= end;}
    void reset() {current = source;}    
  }; 
  
  template<typename T>
  class MappedAccessIterator : public AccessIterator<T> {
    int current; int end;
    BaseBuffer<int>& map;
    T* source; 
    void find_next() {while (current < end && map[current] == 0) current++;}
  public:
    MappedAccessIterator(T* data, BaseBuffer<int>* map) : current(0), end(map->size()), map(*map), source(data) {find_next();}
    T& next_value() {T& out = source[map[current++]-1]; find_next(); return out;}
    pair<int,T> next_element() {int sid = current++; return pair<int,T>(sid,source[map[sid]-1]);}        
    bool empty() {return current >= end;}
    void reset() {current = 0; find_next();}
  };
  
  template<typename T>
  class ConstantIterator : public AccessIterator<T> {
    T value; T source; int current; int end;
  public:
    ConstantIterator(T value, int size) : value(value), current(0), end(size) {}
    T& next_value() {current++; source = value; return source;}
    pair<int,T> next_element() {return pair<int,T>(current++,value);}    
    bool empty() {return current >= end;}
    void reset() {current = 0;}    
  }; 
  
  template<typename T>
  class MappedConstantIterator : public AccessIterator<T> {
    T value; T source; int current; int end;
    BaseBuffer<int>& map;
    void find_next() {while (current < end && map[current] == 0) current++;}
  public:
    MappedConstantIterator(T value, BaseBuffer<int>* map) : value(value), current(0), end(map->size()), map(*map) {}
    T& next_value() {current++; find_next(); source = value; return source;}
    pair<int,T> next_element() {int index = current++; find_next(); return pair<int,T>(index,value);}    
    bool empty() {return current >= end;}
    void reset() {current = 0; find_next();}      
  };
}


template<typename T>
class VariableStorage {
protected:
  BaseBuffer<int>* map;
  
  T default_value;
  void (*on_remove)(T&);
  
  int max_elem;
  bool is_initialised;
  bool is_fixed;
  bool has_default;

  virtual void resize(int size) = 0;
  
  virtual void set_buffers() = 0;
  virtual void set_access() = 0;
  
  static bool smaller(T& value, T& comp) {return value < comp;}
  static bool smaller_equal(T& value, T& comp) {return value <= comp;}
  static bool equal(T& value, T& comp) {return value == comp;}
  static bool greater_equal(T& value, T& comp) {return value >= comp;}  
  static bool greater(T& value, T& comp) {return value > comp;}  
  typedef bool (*FUNC_filter)(T&, T&);  
  
  virtual FUNC_filter get_filter(VariableStorageAccess::FilterMode mode) {
    if (mode == VariableStorageAccess::Smaller) return &smaller; 
    else if (mode == VariableStorageAccess::SmallerEqual) return &smaller_equal;     
    else if (mode == VariableStorageAccess::Equal) return &equal;
    else if (mode == VariableStorageAccess::GreaterEqual) return &greater_equal;
    else if (mode == VariableStorageAccess::Greater) return &greater;
    else return 0;
  }

  VariableStorageAccess::AccessIterator<T>* make_iterator(T* input) {return make_iterator(input, default_value, has_default);}
  VariableStorageAccess::AccessIterator<T>* make_iterator(T* input, BaseBuffer<int>* map_override) {return make_iterator(input, default_value, has_default, map_override);}
  VariableStorageAccess::AccessIterator<T>* make_constant_iterator() {return make_constant_iterator(has_default ? default_value : 0);}  
  
  template<typename U>
  VariableStorageAccess::AccessIterator<U>* make_iterator(U* input, U default_value, bool has_default, BaseBuffer<int>* map_override) { 
    if (!is_initialised) return has_default ? new VariableStorageAccess::DummyAccessIterator<U>(default_value) : new VariableStorageAccess::DummyAccessIterator<U>();
    else if (map_override) return new VariableStorageAccess::MappedAccessIterator<U>(input, map_override);    
    else if (map) return new VariableStorageAccess::MappedAccessIterator<U>(input, map);
    else return new VariableStorageAccess::DirectAccessIterator<U>(input, get_datasize());
  }
  
  template<typename U>
  VariableStorageAccess::AccessIterator<U>* make_constant_iterator(U value) { 
    if (map) return new VariableStorageAccess::MappedConstantIterator<U>(value, map);
    else return new VariableStorageAccess::ConstantIterator<U>(value, get_datasize());
  }

  int get_sid(int iid) {return map ? (map->get(iid) - 1) : iid;}     

public:
  VariableStorage(int expected=0) : map(0), on_remove(0), max_elem(expected), is_initialised(false), is_fixed(false), has_default(false) {}
  virtual ~VariableStorage() {}

  virtual void init() {
    if (is_initialised) return; is_initialised = true;
    resize(get_capacity()); set_access();
  }
  virtual void clear() {is_initialised = false; is_fixed = false; set_access();}

  void set_default(T value) {default_value = value; has_default = true;}
  void set_remove(T missing, void (*func)(T&)) {set_default(missing); on_remove = func; set_access();}
  void set_capacity(int value, bool force_init=false) {max_elem = value; if (is_initialised) resize(get_capacity()); else if (force_init) init();}
  void freeze(BaseBuffer<int>* local_map=0) {map = local_map; is_fixed = true; set_access();}

  virtual int filter(T value, VariableStorageAccess::FilterMode mode) = 0;

  ///assumes no duplications in vector use
  virtual bool remap(vector<int>& use) = 0;

  bool exists(int iid) {return iid >= 0 && iid < max_iid() && get()->exists(iid);} 
  virtual VariableStorageAccess::Accessor<T>* get() = 0; 
  virtual VariableStorageAccess::Accessor<T>* set() = 0;   
  
  bool initialised() {return is_initialised;}
  bool has_map() {return map != 0;}
  int get_capacity() {return max_elem;}
  virtual int get_datasize() = 0;
  
  void drop_iid(int iid) {if (is_initialised && map) map->get(iid) = 0;}
  int max_iid() {return map ? map->size() : max_elem;}
};

template<typename T>
class VariableUniStorage : public VariableStorage<T> {
  friend class VariableStorageAccess::ExpandSetter<T>; 
  BaseBuffer<T> data;
  
  VariableStorageAccess::UniAccessor<T>* getter;
  VariableStorageAccess::UniAccessor<T>* setter;

  void resize(int size) {this->max_elem = size;
    int used = data.size(), change = size - used;
    if (change < 0 && this->on_remove) block_remove(size, used);
    data.resize(size, true); 
    if (change > 0 && this->has_default) Utils::fill_value(data.data()+used, change, this->default_value);
    set_buffers();
  }
  
  void block_remove(int from=0, int to=0) {if (!this->on_remove) return; to = (to == 0 || to > data.size()) ? data.size() : to; for (int i = from; i < to; i++) this->on_remove(data[i]);}  
  
  void set_buffers() {this->getter->set_buffer(data); this->setter->set_buffer(data);}
  void set_access() {
    delete this->getter; delete this->setter;
    if (this->is_initialised) {
      if (this->map) {
        if (this->has_default) {this->getter = new VariableStorageAccess::MappedAccess<T>(data, this->map, this->default_value, this->on_remove); this->setter = new VariableStorageAccess::MappedAccess<T>(data, this->map, this->default_value, this->on_remove);}
        else {this->getter = new VariableStorageAccess::MappedAccess<T>(data, this->map); this->setter = new VariableStorageAccess::MappedAccess<T>(data, this->map);}
      } else {
        this->getter = new VariableStorageAccess::DirectAccess<T>(data); this->setter = new VariableStorageAccess::DirectAccess<T>(data);
        if (!this->is_fixed) this->setter = new VariableStorageAccess::ExpandSetter<T>(*this, this->setter);
      }
    } else {this->getter = new VariableStorageAccess::DummyAccess<T>(); this->setter = new VariableStorageAccess::DummyAccess<T>();}
    if (this->on_remove) this->setter = new VariableStorageAccess::DeleteSetter<T>(this->on_remove, this->setter);    
  }

public:
  VariableUniStorage(int expected=0) : VariableStorage<T>(expected), getter(new VariableStorageAccess::DummyAccess<T>()), setter(new VariableStorageAccess::DummyAccess<T>()) {}
  ~VariableUniStorage() {block_remove(); delete getter; delete setter;}

  void clear() {block_remove(); data.clear(); VariableStorage<T>::clear();}
  
  void assign_data(BaseBuffer<T>& input) {
    if (!this->is_initialised) this->init();
    int use = min(this->max_iid(), int(input.size()));
    for (int i = 0; i < use; i++) setter->elem(i) = input[i];
  }
  
  void assign_value(T value) {if (this->is_initialised) data.assign_value(value);}
  
  int filter(T value, VariableStorageAccess::FilterMode mode) {if (!this->is_initialised || !this->map) return -1; int count = 0, total = this->map->size();
    typename VariableStorage<T>::FUNC_filter filter_func = this->get_filter(mode); 
    for (int i = 0; i < total; i++) {
      int& lid = this->map->get(i);
      if (lid) {if (filter_func(data[lid-1], value)) lid = 0; else count++;} 
    }
    return count;
  }  

  ///assumes no duplications in vector use
  bool remap(vector<int>& use) {if (!this->is_initialised || use.size() > data.size()) return false;
    BaseBuffer<T> new_data(use.size()); T *read = data.data(), *write = new_data.data(); 
    for (int i = 0; i < use.size(); i++) *(write++) = read[use[i]]; 
    if (data.size() > use.size() && this->on_remove) { 
      for (int i = 0; i < use.size(); i++) read[use[i]] = this->default_value;        
      block_remove();
    } 
    data.swap(new_data); this->set_capacity(data.size()); 
    set_buffers(); return true;
  }
    
  VariableStorageAccess::AccessIterator<T>* get_iterator(BaseBuffer<int>* mask=0) {return this->make_iterator(data.data(), mask);} 
  VariableStorageAccess::AccessIterator<T>* get_constant_iterator() {return this->make_constant_iterator();}   
  VariableStorageAccess::UniAccessor<T>* get() {return getter;} 
  VariableStorageAccess::UniAccessor<T>* set() {return setter;}  
  T& get(int iid) {return getter->elem(iid);}
  void set(int iid, T value) {setter->elem(iid) = value;}
  bool is_missing(int iid) {return !this->exists(iid) || (this->has_default && get()->elem(iid) == this->default_value);}  

  BaseBuffer<T>& get_data() {return data;}
  int get_datasize() {return this->is_initialised ? data.size() : 0;}
};


template<typename T>
class VariableMultiStorage : public VariableStorage<T> {
  friend class VariableStorageAccess::ExpandingMultiAccess<T>; 
protected:
  VariableStorageAccess::MultiAccessor<T>* getter;
  VariableStorageAccess::MultiAccessor<T>* setter;
  
  void set_access() {
    delete getter; delete setter;
    if (this->is_initialised && this->map) {
      if (this->has_default) {getter = new VariableStorageAccess::MappedMultiAccess<T>(this->map, this->default_value); setter = new VariableStorageAccess::MappedMultiAccess<T>(this->map, this->default_value);}
      else {getter = new VariableStorageAccess::MappedMultiAccess<T>(this->map); setter = new VariableStorageAccess::MappedMultiAccess<T>(this->map);}
    } else {
      getter = new VariableStorageAccess::DirectMultiAccess<T>(); setter = new VariableStorageAccess::DirectMultiAccess<T>();
      if (this->is_initialised && !this->is_fixed) setter = new VariableStorageAccess::ExpandingMultiAccess<T>(*this, setter);
    }
    this->set_buffers();
  }
  
  virtual bool filter_row(int iid, T filter_value, typename VariableStorage<T>::FUNC_filter filter_func) {
    VariableStorageAccess::ElementIterator<T>& iter = getter->iterator(iid); bool drop = iter.empty();
    if (!drop) {while (!iter.empty()) {if (filter_func(iter.elem(), filter_value)) {drop = true; break;}}}    
    return drop;
  }
  
public:
  VariableMultiStorage(int expected=0) : VariableStorage<T>(expected), getter(new VariableStorageAccess::DirectMultiAccess<T>()), setter(new VariableStorageAccess::DirectMultiAccess<T>()) {}
  virtual ~VariableMultiStorage() {delete getter; delete setter;}

  int filter(T value, VariableStorageAccess::FilterMode mode) {if (!this->is_initialised || !this->map) return -1; int count = 0, total = this->map->size();
    typename VariableStorage<T>::FUNC_filter filter_func = this->get_filter(mode); 
    for (int i = 0; i < total; i++) {if (filter_row(i, value, filter_func)) this->map->get(i) = 0; else count++;}
    return count;
  }  
  
  VariableStorageAccess::MultiAccessor<T>* get() {return getter;} 
  VariableStorageAccess::MultiAccessor<T>* set() {return setter;}
  virtual T& get_elem(int iid, int col) {return getter->elem(iid).get(col);}
  virtual void set_elem(int iid, int col, T value) {setter->elem(iid).set(col, value);}
  VariableStorageAccess::MultiAccessElement<T>& get_row(int iid, bool to_set) {return to_set ? setter->elem(iid) : getter->elem(iid);}
  VariableStorageAccess::AccessIterator<T>* get_constant_iterator() {return this->make_constant_iterator();}       
};

template<typename T>
class VariableMultiBlockStorage : public VariableMultiStorage<T> {
  Buffer<T> data;
  int ncol;
  
  void resize(int size) {this->max_elem = size;
    if (data.empty()) {
      data.set_size(size, ncol);
      if (this->has_default) data.assign_value(this->default_value);
    } else {
      if (size < data.nrow() && this->on_remove) block_remove(size, data.nrow()); 
      data.fit_rows(size, this->default_value, this->has_default);
    }
    set_buffers();
  }
  
  void block_remove(int from_row=0, int to_row=0, int from_col=0, int to_col=0) {if (!this->on_remove) return; 
    to_row = (to_row == 0 || to_row > data.nrow()) ? data.nrow() : to_row; 
    to_col = (to_col == 0 || to_col > data.ncol()) ? data.ncol() : to_col;     
    for (int i = from_row; i < to_row; i++) {for (int j = from_col; j < to_col; j++) this->on_remove(data(i,j));}
  }    

  void set_buffers() {
    if (this->is_initialised) {
      this->getter->set_access(new VariableStorageAccess::MultiAccessBlock<T>(data), this->on_remove);
      this->setter->set_access(new VariableStorageAccess::MultiAccessBlock<T>(data), this->on_remove);
    }  
  }
  
public:
  VariableMultiBlockStorage(int expected=0) : VariableMultiStorage<T>(expected), ncol(0) {}
  ~VariableMultiBlockStorage() {block_remove();}

  void set_width(int width) {
    if (this->is_initialised && width != ncol) {
      if (width < ncol) {
        block_remove(0, 0, width);
        data.shrink(width); Buffer<T> new_data(data);      
        data.swap(new_data);
      } else data.total_cols(width, this->has_default, this->default_value);  
      set_buffers();
    }
    ncol = width;
  }
  void clear() {block_remove(); data.reset(); VariableMultiStorage<T>::clear();}

  void assign_value(T value) {if (this->is_initialised) {block_remove(); data.assign_value(value);}}

  ///assumes no duplications in vector use
  bool remap(vector<int>& use) {if (!this->is_initialised || use.size() > data.nrow()) return false;
    Buffer<T> new_data(use.size(), ncol);
    for (int i = 0; i < use.size(); i++) {int sid = use[i];
      for (int j = 0; j < ncol; j++) new_data(i,j) = data(sid,j);    
    }    
    if (data.nrow() > use.size() && this->on_remove) { 
      for (int i = 0; i < use.size(); i++) {int sid = use[i];
        for (int j = 0; j < ncol; j++) data(sid,j) = this->default_value;    
      }    
      block_remove();
    } 
    data.swap(new_data); this->set_capacity(data.nrow());
    set_buffers(); return true;
  }

  VariableStorageAccess::AccessIterator<T>* get_iterator(int col, BaseBuffer<int>* mask=0) {return this->make_iterator(data[col], mask);} 

  Buffer<T>& get_data() {return data;}  
  int width() {return ncol;}
  int get_datasize() {return this->is_initialised ? data.nrow() : 0;}
};

template<typename T>
class VariableMultiJaggedStorage : public VariableMultiStorage<T> {
  BaseBuffer<BaseBuffer<T>*> data;
  BaseBuffer<T> dummy_buffer;
  bool purge_filter;
  int increment;

  void resize(int size) {this->max_elem = size;
    int used = data.size(), change = size - used;
    if (change < 0) block_remove(size, used);
    data.resize(size, true);
    if (change > 0) Utils::set_zero(data.data()+used, change);
    set_buffers();
  }

  void row_remove(BaseBuffer<T>* row) {
    if (row && this->on_remove) {for (int i = 0; i < row->size(); i++) this->on_remove(row->get(i));}  
  }

  void block_remove(int from=0, int to=0) {
    to = (to == 0 || to > data.size()) ? data.size() : to; 
    for (int i = from; i < to; i++) {BaseBuffer<T>*& elem = data[i];
      row_remove(elem);
      delete elem; elem = 0;  
    }
  }   
  
  void set_buffers() {
    if (this->is_initialised) {
      this->getter->set_access(new VariableStorageAccess::MultiAccessJagged<T>(data), this->on_remove);
      this->setter->set_access(new VariableStorageAccess::MultiAccessJagged<T>(data), this->on_remove);
    }  
  } 
  
  bool purge_row(int iid, T filter_value, typename VariableStorage<T>::FUNC_filter filter_func) {
    VariableStorageAccess::MultiAccessElement<T>& elem = this->getter->elem(iid); int total = elem.size(), write = 0;
    T null_value = this->has_default ? this->default_value : 0;
    for (int i = 0; i < total; i++) {T& curr = elem.get(i);
      if (!filter_func(curr, filter_value)) {
        if (write < i) {elem.set(write, curr); curr = null_value;}
        write++;
      }
    }    
    if (write < total) {int sid = this->getter->get_sid(iid);
      if (write == 0) {delete data[sid]; data[sid] = 0;}
      else data[sid]->resize(write, true);
    }
    return (write == 0);
  }
  
  bool filter_row(int iid, T filter_value, typename VariableStorage<T>::FUNC_filter filter_func) {
    if (purge_filter) return purge_row(iid, filter_value, filter_func);
    else if (VariableMultiStorage<T>::filter_row(iid, filter_value, filter_func)) {
      int sid = this->getter->get_sid(iid);
      if (sid >= 0) {delete data[sid]; data[sid] = 0;}
      return true;
    } else return false;
  }
  
  BaseBuffer<T>* fetch_buffer(int iid) {
    int sid = this->setter->elem(iid).get_row();
    if (sid >= 0) {
      if (!data[sid]) data[sid] = new BaseBuffer<T>();
      return data[sid];
    } else return 0; 
  }

public:
  VariableMultiJaggedStorage(int expected=0) : VariableMultiStorage<T>(expected), purge_filter(true), increment(10) {}
  ~VariableMultiJaggedStorage() {block_remove();}

  void clear() {block_remove(); data.clear(); VariableMultiStorage<T>::clear();}
  void set_purge(bool value) {purge_filter = value;}

  ///assumes no duplications in vector use
  bool remap(vector<int>& use) {if (!this->is_initialised || use.size() > data.size()) return false;
    BaseBuffer<BaseBuffer<T>*> new_data(use.size());
    for (int i = 0; i < use.size(); i++) {int sid = use[i];
      new_data[i] = data[sid]; data[sid] = 0;
    }    
    block_remove();
    data.swap(new_data); this->set_capacity(data.size());
    set_buffers(); return true;    
  }

  void add_elem(int iid, T value) {BaseBuffer<T>* curr = fetch_buffer(iid);
    if (curr) {
      if (!curr->remaining()) curr->reserve(curr->size() + increment);
      curr->append(1); curr->get(curr->size()-1) = value;
    }
  }

  void assign_data(int iid, BaseBuffer<T>& input) {BaseBuffer<T>* curr = fetch_buffer(iid);
    if (curr) {
      if (!curr->is_empty()) row_remove(curr);
      curr->assign(input); 
    }
  }

  template<typename C>  
  void insert_data(int iid, C& input) {BaseBuffer<T>* curr = fetch_buffer(iid);
    if (curr) {
      if (!curr->is_empty()) row_remove(curr);
      curr->resize(input.size());
      T* write = curr->data();
      for (typename C::iterator it = input.begin(); it != input.end(); ++it) *(write++) = *it;
    }        
  }

  void remove_data(int iid) {if (!this->is_initialised) return;
    int sid = this->setter->elem(iid).get_row();
    if (sid >= 0) {
      BaseBuffer<T>*& elem = data[sid];
      row_remove(elem); delete elem; elem = 0;  
    }  
  }
  
  void remove_data() {block_remove();}
  
  VariableStorageAccess::AccessIterator<BaseBuffer<T>*>* get_iterator(BaseBuffer<int>* mask=0) {return this->make_iterator(data.data(), (BaseBuffer<T>*) 0, true, mask);} 

  BaseBuffer<BaseBuffer<T>*>& get_data() {return data;}  
  BaseBuffer<T>& get_data(int iid) {
    BaseBuffer<T>* buff = fetch_buffer(iid);
    if (!buff) {dummy_buffer.clear(); return dummy_buffer;}
    else return *buff;
  }
  int width(int iid) {return this->getter->elem(iid).size();}
  int get_datasize() {return this->is_initialised ? data.size() : 0;}
};


template<typename T>
void VariableStorageAccess::ExpandSetter<T>::init() {capacity = owner.get_capacity(); set_buffer(owner.data);}

template<typename T>
void VariableStorageAccess::ExpandSetter<T>::update_capacity() {capacity = owner.get_capacity();}

template<typename T>
void VariableStorageAccess::ExpandSetter<T>::expand(int min_iid) {
  int cap = capacity + max(int(capacity * expansion_factor), int(10));
  owner.resize(max(cap, min_iid+1));
  capacity = owner.get_capacity();
} 

template<typename T>
void VariableStorageAccess::ExpandingMultiAccess<T>::update_capacity() {capacity = owner.get_capacity();}

template<typename T>
void VariableStorageAccess::ExpandingMultiAccess<T>::expand(int min_iid) {
  int cap = capacity + max(int(capacity * expansion_factor), int(10));
  owner.resize(max(cap, min_iid+1));
  capacity = owner.get_capacity();
} 


#endif /** UTILS_DATAFRAME_VARIABLES_STORAGE_H */

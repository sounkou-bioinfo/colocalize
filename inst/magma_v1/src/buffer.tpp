/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

template<typename T>
void BaseBuffer<T>::do_copy(T* to, T* from, int total) {memcpy(to, from, sizeof(T)*total);}

template<>
void BaseBuffer<string>::do_copy(string* to, string* from, int total);

template<typename T>
void BaseBuffer<T>::set_size(long new_size) {
  if (new_size > 0) {
    T* new_content = new T[new_size];
    if (content) {do_copy(new_content, content, min(no_elem, new_size)); delete[] content;}
    content = new_content;
    no_elem = new_size;
    capacity = new_size;
  } else clear();
}

template<typename T>
void BaseBuffer<T>::resize(long new_size, bool shrink_to_fit) {
  if (new_size > capacity || (new_size < capacity && shrink_to_fit)) set_size(new_size);
  else no_elem = new_size;
}

template<typename T>
void BaseBuffer<T>::reserve(long total) {
  if (total > capacity) {
    int curr_size = no_elem;
    set_size(total);
    no_elem = curr_size;
  }
}

template<typename T>
void BaseBuffer<T>::assign(T* source, long len) {
  if (!content || len > capacity) {delete[] content; content = new T[len]; capacity = len;}
  do_copy(content, source, len); no_elem = len;
}   

template<typename T>
void BaseBuffer<T>::assign_value(T value, int amount) {
  if (amount >= 0) resize(amount); 
  if (no_elem > 0) Utils::fill_value(content, no_elem, value);
}

template<typename T>
void BaseBuffer<T>::swap(BaseBuffer<T>& other) {
  if (this == &other) return;
  std::swap(content, other.content);
  std::swap(capacity, other.capacity);    
  std::swap(no_elem, other.no_elem);
}

template<typename T>  
void ExpandingBuffer<T>::expand(long min_elem) {
  long new_size = max(long(PARENT::no_elem*expansion_factor), min_elem);
  PARENT::resize(new_size);
}

template<typename T>
void ExpandingBuffer<T>::expansion_reserve(long size) {
  if (size >= PARENT::capacity) {
    PARENT::set_size(max(long(PARENT::capacity * expansion_factor), size+1)); 
    PARENT::no_elem = size+1;
  }
}

template<typename T>
T& ExpandingBuffer<T>::get(const long& index) {
  if (index >= PARENT::capacity) expand(index+1); 
  else if (index >= PARENT::no_elem) PARENT::no_elem = index + 1; 
  return PARENT::content[index];
}  



/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "buffer.h"

template<>
void BaseBuffer<string>::do_copy(string* to, string* from, int total) {for (int i = 0; i < total; i++) to[i] = from[i];}


const double ArchiveBuffer::EXPANSION_FACTOR = 1.1;
const int ArchiveBuffer::EXPANSION_MINIMUM = 10;

void ArchiveBuffer::init(long amount, int exp_size) {
  if (!active) {
    if (exp_size > 0) elem_size = exp_size;
    if (amount > 0) min_elem = amount;
    expand();
  }
}

void ArchiveBuffer::clear() {
  for (int i = 0; i < content.size(); i++) delete[] content[i];
  content.clear(); active = 0; 
  capacity = 0; no_elem = 0; total_used = 0;
}

void ArchiveBuffer::expand(long minimum) {
  if (active) total_used += active - content.back();
  if (no_elem > 100 && total_used > 0) elem_size = ceil(total_used / double(no_elem));

  max_elem = max(max(int(no_elem*EXPANSION_FACTOR), int(no_elem+EXPANSION_MINIMUM)), min_elem);
  capacity = max(minimum, elem_size * (max_elem - no_elem));
 
  content.push_back(new char[capacity]);
  active = content.back();
}

char* ArchiveBuffer::add_internal(const char* str, int len) {
  if (len+2 > capacity) expand(len+1);
  *(active++) = (len < 127) ? char(len) : 127;
  memcpy(active, str, len); active[len] = '\0'; 
  char* out = active; active += len+1; capacity -= len+2; no_elem++;

  return out;
}

char* ArchiveBuffer::add(char* str, bool from_archive) {int len = 127;
  if (from_archive) len = *(str-1);
  if (len >= 127) len = strlen(str);
  return add_internal(str, len);
}

char* IndexedArchiveBuffer::add_internal(const char* str, int len) {
  if (len+2 > capacity || no_elem >= max_elem) expand(len+1);
  *(active++) = (len < 127) ? char(len) : 127;
  index[no_elem++] = active;
  memcpy(active, str, len); active[len] = '\0'; 
  active += len+1; capacity -= len+2; 
  return index[no_elem-1];
}
  
bool BitBuffer::global_init = false;
bitBlock BitBuffer::write_mask[8];
bitBlock BitBuffer::count_mask[256];

void BitBuffer::mask_init() {
  bitBlock value = 128;
  for (int i = 0; i < 8; value >>= 1, i++) write_mask[i] = value;
  for (int i = 0; i < 256; i++) count_mask[i] = (i & 1) + count_mask[i / 2];

  global_init = true;
}

deque<int>& BitBuffer::traverse(int col) {
  traverse_index.clear();  
  bitBlock* data = content+col*row_blocks;
  for (int b = 0; b < row_blocks; b++) {
    bitBlock& curr = data[b];
    if (curr) {for (int i = 0; i < 8; i++) {if (curr & write_mask[i]) traverse_index.push_back(8*b + i);} }  
  }
  while (!traverse_index.empty() && traverse_index.back() >= rows) traverse_index.pop_back();
  return traverse_index;
}

int BitBuffer::count_internal(bitBlock* c) {int total = 0;
  for (bitBlock *end = c + row_blocks; c < end; c++) total += count_mask[*c];
  return total;
}
 
int BitBuffer::both_internal(bitBlock* c1, bitBlock* c2) {int total = 0;
  for (bitBlock *end = c1 + row_blocks; c1 < end; c1++, c2++) total += count_mask[*c1 & *c2];
  return total;
}

int BitBuffer::either_internal(bitBlock* c1, bitBlock* c2) {int total = 0;
  for (bitBlock *end = c1 + row_blocks; c1 < end; c1++, c2++) total += count_mask[*c1 | *c2];
  return total;
}


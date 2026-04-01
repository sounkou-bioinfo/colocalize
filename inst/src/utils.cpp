/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "utils.h"
#include <sys/stat.h>
#include <iostream>

#define _FILE_OFFSET_BITS 64

bool Utils::is_dir(const string& basename) {
  struct stat status;
  return stat(basename.c_str(), &status) == 0 && S_ISDIR(status.st_mode);
}

bool Utils::is_file(const string& filename) {
  struct stat status;
  int code = stat(filename.c_str(), &status);
  return code == 0 && S_ISREG(status.st_mode);
}

pair<string,string> Utils::split_path(const string& filename) {
  pair<string,string> out;
  int offset = filename.find_last_of('/');
  if (offset == string::npos) {
    out.second = filename;  
  } else {
    if (offset != filename.size()-1) out.second = filename.substr(offset+1);
    while (filename[offset] == '/') offset--;
    out.first = filename.substr(0, offset+2);      
  }    
  return out;
}

pair<string,string> Utils::split_str(const string& str, char split) {
  pair<string,string> out;
  int offset = str.find(split);
  if (offset == string::npos) {
    out.first = str;
  } else {
    if (offset != 0) out.first = str.substr(0,offset);
    if (offset != str.size()-1) out.second = str.substr(offset+1);
  }
  return out;
}

string Utils::num_name(int num) {
  static vector<string> names = tokenize("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen nineteen twenty");
  return (num >= 0 && num < names.size()) ? names[num] : num_string(num);
}

vector<string> Utils::tokenize(const string& str) {
  vector<string> out; istringstream istr(str); string curr;
  while (istr >> curr) out.push_back(curr);
  return out;
}

vector<string> Utils::tokenize(const string& str, char sep) {
  vector<string> out; istringstream istr(str); string curr;
  while (istr.good()) {getline(istr, curr, sep); out.push_back(curr);}
  return out;
}

string Utils::join_string(vector<string>& values, string sep) {
  if (!values.empty()) {
    string output = values[0];
    for (int i = 1; i < values.size(); i++) output += sep + values[i];
    return output;
  } else return "";
}

string Utils::lowercase(const string& str) {
  string out(str);
  for (int i = 0; i < out.size(); i++) out[i] = tolower(out[i]);
  return out;
}

vector<string> Utils::lowercase(const vector<string>& vec) {
  vector<string> out; out.reserve(vec.size());
  for (int i = 0; i < vec.size(); i++) out.push_back(lowercase(vec[i]));
  return out;
}

string Utils::uppercase(const string& str) {
  string out(str);
  for (int i = 0; i < out.size(); i++) out[i] = toupper(out[i]);
  return out;
}

vector<string> Utils::uppercase(const vector<string>& vec) {
  vector<string> out; out.reserve(vec.size());
  for (int i = 0; i < vec.size(); i++) out.push_back(uppercase(vec[i]));
  return out;
}


string Utils::str_replace(const string& source, const string& from, const string& to) {
  string out = source; size_t found; int len = from.size();
  while ((found = out.find(from)) != string::npos) out.replace(found, len, to);
  return out;
}

bool Utils::starts_with(const string& input, const string& match) {return (match.size() <= input.size()) && (input.substr(0, match.size()) == match);}
bool Utils::ends_with(const string& input, const string& match) {return (match.size() <= input.size()) && (input.substr(input.size() - match.size()) == match);}

bool Utils::empty_line(const string& line, bool skip_comment) {
  int pos = 0;
  while (isspace(line[pos])) if (++pos == line.size()) return true;
  return skip_comment && line[pos] == '#';
}

string Utils::chr_string(long chr, int pad) {
  string out;
  if (chr >= chr_reserved) {
    switch(chr) {
      case chrX_code: out = "X"; break;
      case chrY_code: out = "Y"; break;
      case chrZ_code: out = "Z"; break;
      case chrW_code: out = "W"; break;
      case chrXY_code: out = "XY"; break;
      case chrMT_code: out = "MT"; break;
      default: return "?";
    }
  } else {
    out = num_string(chr);
  }

  return pad_string(out, pad);
}
 
string Utils::get_letters(int offset, bool uppercase) {
  string out = string(1, char((uppercase ? 'A' : 'a') + offset%26)); 
  return (offset >= 26) ? get_letters(offset/26-1, uppercase) + out : out;
}
 


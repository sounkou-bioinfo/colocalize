/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <exception>

using namespace std;

class BaseException : public exception {
  const string exception;
  const string type;
  const string msg;
  const string sub_msg;
  
public:
  BaseException(const string t, const string m, const string wh, const string e="BASE") : exception(e), type(t), msg(m), sub_msg(wh) {}
  virtual ~BaseException() throw() {}
  virtual const char* what() const throw() {return msg.c_str();}
  const string& get_exception() const {return exception;}
  const string& get_type() const {return type;}
  const string get_msg() const {
    if (sub_msg != "") return string(msg).append("\n\t\tinternal error message: ").append(sub_msg);
    return msg;
  }
};

struct BufferException : public BaseException {
  BufferException(const string m="unknown", const string wh="") : BaseException("accessing buffer", m, wh, "BUFFER") {}
};

struct MathException : public BaseException {
  MathException(const string t, const string m="unknown", const string wh="") : BaseException(t, m, wh, "MATH") {}
};

struct StatException : public BaseException {
  StatException(const string t, const string m="unknown", const string wh="") : BaseException(t, m, wh, "STAT") {}
};

struct DataException : public BaseException {
  DataException(const string t, const string m="unknown", const string wh="") : BaseException(t, m, wh, "DATA") {}
};

struct GeneException : public BaseException {
  GeneException(const string t, const string m="unknown", const string wh="") : BaseException(t, m, wh, "GENE") {}
};

struct EngineException : public BaseException {
  EngineException(const string t, const string m="unknown", const string wh="") : BaseException(t, m, wh, "ENGINE") {}
};

struct SetException : public BaseException {
  SetException(const string t, const string m="unknown", const string wh="") : BaseException(t, m, wh, "SET") {}
};

class TypeMismatch : public exception {
  const string found;
  const string expected;
  string msg;
  
public:
  TypeMismatch(const string found, const string expected) : found(found), expected(expected) {
    msg = "found value of type '" + found + "' (expected type '" + expected + "')";
  }
  virtual ~TypeMismatch() throw() {}
  virtual const char* what() const throw() {return msg.c_str();}
};

#endif /**EXCEPTION_H*/
/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GLOBAL_H
#define GLOBAL_H

#define EIGEN_NO_DEBUG

#include <typeinfo>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include "output.h"

extern boost::random::mt19937 _RNG;
extern Log _LOG;
extern SupplementaryLog _SLOG;

typedef unsigned long long timestamp_t;

inline static timestamp_t get_timestamp () {
  struct timeval now;
  gettimeofday (&now, NULL);
  return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

inline static float get_elapsed(timestamp_t start) {
  return (float) (get_timestamp() - start)/1000000;
}

inline static void mark(int id) {
  static timestamp_t start = get_timestamp();
  if (id == 0) start = get_timestamp();
  cerr << "mark_" << id << " " << get_elapsed(start) << endl;
} 

inline static void wait() {
  cout << endl << "PAUSED" << endl;
  if (getchar()) return;
}

namespace ExitType {
  enum Value {
    NoError = 0,
    BenignError = 1,
    GenericError = 2,
    MemError = 3
  };
}

inline static void die(int code=ExitType::GenericError) {
  if (code != ExitType::NoError) {
    if (code == ExitType::BenignError) _LOG << "\nExiting MAGMA. Goodbye." << endl;
    else _LOG << "\nTerminating program." << endl;
  }
  exit(code);
}

inline static void check_mem_error(const exception& e, const string& msg) {
  if (typeid(e) == typeid(bad_alloc)) {
    _LOG.error(msg, 2) << "an error occurred when trying to allocate memory" << endl;
    die(ExitType::MemError);
  }
}

inline static void report_error(const exception& e, const string& msg) {
  check_mem_error(e, msg);
  _LOG.error(msg, 2) << "an unknown error occurred" << endl;
  _LOG << "\tinternal error message: " << e.what() << endl;
  die();
}

#endif /**GLOBAL_H*/
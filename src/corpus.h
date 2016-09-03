#ifndef CORPUS_H
#define CORPUS_H

#include <string>
#include <vector>

#include "param.h"
#include "rand.h"

using namespace std;

class triplet {
 public:
  int docid;
  int wordid;
  int count;

 public:
  inline void read(const string& s);
  inline bool good();
};

#define USE_UCHAR_FOR_WID 0

class document {
 public:
  int docid;
  int len;
  int *wid;
#if USE_UCHAR_FOR_WID
  unsigned char *cnt;
#else
  int *cnt;
#endif
  int *zsmp;

 public:
  void set_words(vector<triplet> *v, int ntopic);
  double sample_zsmp(int i, float* zprb, int ntopic, float pr);
};

class corpus {
 public:
  param *prm;

  int ndoc;
  int nvoc;
  document *doc;

  // copy params
  int num_topics;

 public:
  corpus(param *prm);
  void load(char* fn, int ntopic);
  void show_stats();
  void init_zsmp();

 private:
  int s2i(const string& s);
  random_nr *rnd;

  // for ntst
 public:
  void error_exit(string msg, int lineno, string line);
  void load_tpc();
  inline float get_p_wk(int w, int k);
  inline void set_p_wk(int w, int k, float v);
 private:
  float *p_wk;			// to load TPC
};

inline float
corpus::get_p_wk(int w, int k) {
  int idx = w * num_topics + k;
  return p_wk[idx];
}

inline void
corpus::set_p_wk(int w, int k, float v) {
  int idx = w * num_topics + k;
  p_wk[idx] = v;
}

#endif

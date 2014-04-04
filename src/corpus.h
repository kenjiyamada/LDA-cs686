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
};

class corpus {
 public:
  param *prm;

  int ndoc;
  int nvoc;
  document *doc;

 public:
  corpus(param *prm);
  void load(char* fn, int ntopic);
  void show_stats();
  void init_zsmp(int ntopic);

 private:
  int s2i(const string& s);
  random_nr *rnd;
};

#endif

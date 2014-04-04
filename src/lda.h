#ifndef LDA_H
#define LDA_H

#include "param.h"
#include "corpus.h"
#include "gibbs.h"

using namespace std;

class lda {
  param *prm;
  corpus *cp;
  gibbs *gb;

 public:
  lda(param *prm);
  ~lda();

  void load();
  void run();

 private:
  int s2i(const string& s);

};

#endif

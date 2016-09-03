#ifndef GIBBS_H
#define GIBBS_H

#include <thread>

#include "param.h"
#include "corpus.h"
#include "rand.h"

using namespace std;

class idprob {
 public:
  int id;
  float prob;
};

class counts {
  int num_topics;    // local copy
 public:
  int *n_kj;
  int *n_wk;
  int *n_k;

  float *zprb;
  random_nr *rnd;

  // for ppx
  float *p_kj;
  double trn_logprob,tst_logprob;
  int trn_wcnt,tst_wcnt;

  void alloc(int,int);

  inline void inc_n_kj(int k, int cnt);
  inline void inc_n_wk(int w, int k, int cnt);
  inline void inc_n_k(int k, int cnt);

  inline int get_n_kj(int k);
  inline int get_n_wk(int w, int k);
  inline int get_n_k(int k);

  inline void set_n_kj(int k, int c);
  inline void set_n_wk(int w, int k, int c);
  inline void set_n_k(int k, int c);
};

class gibbs {
 public:
  corpus *cp;
  param *prm;

  // copy params
  bool verbose;
  float alpha,beta;
  int num_iterations, num_topics;
  int num_vocab;
  int num_procs;
  int ppx_every;
  int heldout_every;

 public:
  counts *gcp;    // global gibbs counts
  counts *lcnts;    // local gibbs counts

  vector<idprob> *p_wk;    // used for ppx

 public:
  gibbs(param*,corpus*);

  void iterate();
  void iterate_uni(int,int,int);

  void get_wprob(int tpc,vector<idprob>&);
  void get_tprob(document*,counts*,bool);
  void get_tprob(document*,counts*,vector<idprob>&,bool);
  void dump_wprob(const char*);
  void dump_tprob(const char*);
  void show_ppx_uni();
  void show_ppx(int);
  void parallel_ppx();

 public:
  void alloc_counts();

  void copy_counts();
  void collect_counts();
  void parallel_sweep();

  void init_n_wk();
  void init_n_kj(document *dp, counts* lcp, bool ppx);

  void remove_counts(document*,int,counts*);
  void add_counts(document*,int,counts*);
  double update_zsmp(document*,int,counts*);

  // for ntst
  void test_new_doc();
  void update_zsmp_ntst(document *dp, int i, counts *gcp);

 private:
  inline bool takeppx(int);
  inline bool heldout(int);
};

class pinfo {
 public:
  int pid;
  pthread_t thr;
  gibbs *gb;
  int res;
  double res5[5];
};

// inline definitions for counts

inline void
counts::inc_n_wk(int w, int k, int cnt){
  int idx = w * num_topics + k;
  n_wk[idx] += cnt;
}

inline void
counts::inc_n_k(int k, int cnt){
  n_k[k] += cnt;
}

inline void
counts::inc_n_kj(int k, int cnt) {
  n_kj[k] += cnt;
}

//

inline int
counts::get_n_wk(int w, int k){
  int idx = w * num_topics + k;
  return n_wk[idx];
}

inline int
counts::get_n_k(int k){
  return n_k[k];
}

inline int
counts::get_n_kj(int k) {
  return n_kj[k];
}

//

inline void
counts::set_n_wk(int w, int k, int c){
  int idx = w * num_topics + k;
  n_wk[idx] = c;
}

inline void
counts::set_n_k(int k, int c){
  n_k[k] = c;
}

inline void
counts::set_n_kj(int k, int c) {
  n_kj[k] = c;
}

inline bool
gibbs::takeppx(int i) {
  return(ppx_every==1 || (ppx_every>1 && i%ppx_every==1));
}

inline bool
gibbs::heldout(int i) {
  return(heldout_every>1 && i%heldout_every==1);
}

#endif

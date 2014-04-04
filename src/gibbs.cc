#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "gibbs.h"

void
counts::alloc(int ntopic, int nvoc) {
  num_topics = ntopic;		// make local copy
  n_kj = new int[ntopic];
  n_wk = new int[ntopic * nvoc];
  n_k = new int[ntopic];

  p_kj = new float[ntopic];
  zprb = new float[ntopic];
  rnd = new random_nr();
  rnd->set_seed(1234);
}

gibbs::gibbs(param *p, corpus *c) {
  prm=p;
  verbose=(p->verbose>0);

  // copy params
  alpha=p->alpha; beta=p->beta;
  num_iterations = p->num_iterations;
  num_topics = p->num_topics;
  num_procs = p->num_procs;
  ppx_every = p->ppx_every;
  heldout_every = p->heldout_every;

  cp=c;
  num_vocab = cp->nvoc;
}

void
gibbs::alloc_counts(){
  // alloc global counts
  gcp = new counts();
  gcp->alloc(num_topics,num_vocab);

  // alloc local counts
  lcnts = new counts[num_procs];
  for (int i=0; i<num_procs; i++) lcnts[i].alloc(num_topics,num_vocab);

  // for ppx
  p_wk = new vector<idprob>[num_topics];
}

void
gibbs::iterate() {
  int ntopic = num_topics;
  int niter = num_iterations;

  alloc_counts();		// both global and local counts

  cp->init_zsmp(ntopic);	// init by random
  init_n_wk();			// also init n_k (from zsmp), global only
  cerr << "start iteration:";
  for (int it=0; it<niter; it++) {
    if (verbose) cerr << "iteration " << it << ": ";
    copy_counts();
    parallel_sweep();
    collect_counts();

    //show_ppx_uni();
    show_ppx(it);
    if (verbose) cerr << "\n";
  }
  if (!verbose) cerr << "\n";
}

void
gibbs::init_n_wk() {
  // init both n_wk and n_k, both global and local
  if (verbose) cerr << "init_n_wk():";

  int ndoc = cp->ndoc;
  int nvoc = num_vocab;
  int ntopic = num_topics;

  // clear counts
  int len_wk = nvoc*ntopic;
  for (int i=0; i<len_wk; i++) gcp->n_wk[i]=0;
  for (int i=0; i<ntopic; i++) gcp->n_k[i]=0;

  // add all counts
  for (int j=0; j<ndoc; j++) {
    bool doppx = takeppx(j);
    if (verbose && j%10000==0) cerr << "[" << j/10000 << "]";
      document *dp = &(cp->doc[j]);
      int dlen=dp->len;
      for (int i=0; i<dlen; i++) {
	if (doppx && heldout(i)) continue;

	// don't forget to use cnt[]
	int wid = dp->wid[i];
	int cnt = dp->cnt[i];
	int tpc = dp->zsmp[i];

	// increment n_wk[wk] and n_k[k] by cnt, where k=tpc
	gcp->inc_n_wk(wid,tpc,cnt);
	gcp->inc_n_k(tpc,cnt);
      }
  }
  if (verbose) cerr << "done\n";
}

void
gibbs::init_n_kj(document *dp, counts *lcp, bool doppx) {
  int ntopic = num_topics;

  // clear counts
  for (int i=0; i<ntopic; i++) lcp->n_kj[i]=0;

  // add counts
  int dlen=dp->len;
  for (int i=0; i<dlen; i++) {
    if (doppx && heldout(i)) continue;
    int cnt = dp->cnt[i];
    int tpc = dp->zsmp[i];
    lcp->inc_n_kj(tpc,cnt);
  }
}

void
gibbs::remove_counts(document *dp, int i, counts *lcp) {
  int wid = dp->wid[i];
  int cnt = dp->cnt[i];
  int tpc = dp->zsmp[i];

  lcp->inc_n_wk(wid,tpc,-cnt);
  lcp->inc_n_k(tpc,-cnt);
  lcp->inc_n_kj(tpc,-cnt);
}

void
gibbs::add_counts(document *dp, int i, counts *lcp) {
  int wid = dp->wid[i];
  int cnt = dp->cnt[i];
  int tpc = dp->zsmp[i];

  lcp->inc_n_wk(wid,tpc,cnt);
  lcp->inc_n_k(tpc,cnt);
  lcp->inc_n_kj(tpc,cnt);
}

double
gibbs::update_zsmp(document *dp, int i, counts *lcp) {
  // will update global zsmp[], but will touch only for documents that are assigned
  // to the processor, so there will be no conflict.

  int wid = dp->wid[i];

  // init zprb array
  int ntopic = num_topics;
  float *zprb = lcp->zprb;
  for (int k=0; k<ntopic; k++) zprb[k]=0;

  // return value for debug
  double zret=0;
  bool zset=false;

  // obtain estimated zprb
  int nvoc = num_vocab;
  float cum=0;
  for (int k=0; k<ntopic; k++) {
    float prb = (alpha + lcp->get_n_kj(k)) * (beta + lcp->get_n_wk(wid,k)) / (nvoc * beta + lcp->get_n_k(k));
    cum += prb;
    zprb[k] = cum;
  }
  
  // sample topic
  float pr = lcp->rnd->get_rand();
  float ps = pr * cum;
  int kx = 0;

#define USE_BINARY_SEARCH 1
#if USE_BINARY_SEARCH

  int min=0;
  int max=ntopic-1;
  if (ps<zprb[min]) kx=min;
  else if (ps>zprb[max]) kx=max;
  else {
    while(max-min>1) {
      int mid=(min+max)/2;
      float pm=zprb[mid];
      if (pm>ps) max=mid;
      else min=mid;
    }
    if (ps>zprb[max] || ps<zprb[min]) cerr << "unexpected error in binary search\n";
    kx=max;
  }

#else
  for (int k=0; k<ntopic; k++) { 
    if (zprb[k]>ps) break;
    kx++;
  }
#endif

  // update zsmp
  dp->zsmp[i] = kx;

  // return debug value
  zret = kx;
  return zret;
}

//
bool
by_prob(const idprob &x, const idprob &y){
  return (x.prob > y.prob);
}

void
gibbs::get_wprob(int k, vector<idprob>& vec){
  int nvoc = num_vocab;
  vec.clear();
  float psum=0;
  for (int i=0; i<nvoc; i++) {
    idprob ip;
    ip.id = i;
    ip.prob = (beta + gcp->get_n_wk(i,k)) / (nvoc*beta + gcp->get_n_k(k));
    psum += ip.prob;
    vec.push_back(ip);
  }
  // normalize
  int vlen = vec.size();
  for (int i=0; i<vlen; i++) vec[i].prob /= psum;
}

void
gibbs::get_tprob(document *dp, counts *lcp, bool doppx){
  // set lcp->p_kj[k]

  int ntopic = num_topics;
  int dlen=dp->len;
  if (doppx && heldout_every>1) dlen *= (1-1.0/heldout_every);

  // get n_kj() for this doc
  init_n_kj(dp,lcp,doppx);

  float *p_kj = lcp->p_kj;
  float psum=0;
  for (int k=0; k<ntopic; k++) {
    float prob = (alpha + lcp->get_n_kj(k)) / (ntopic*alpha + dlen);
    p_kj[k] = prob;
    psum += prob;
  }

  // normalize
  for (int k=0; k<ntopic; k++) p_kj[k] /= psum;
}


void
gibbs::get_tprob(document *dp, counts *lcp, vector<idprob>& vec, bool doppx){
  // first set lcp->p_kj[k]
  get_tprob(dp,lcp,doppx);

  // then, copy lcp->p_kj[k] into vec
  int ntopic = num_topics;
  vec.clear();
  for (int k=0; k<ntopic; k++) {
    idprob ip;
    ip.id = k;
    ip.prob = lcp->p_kj[k];
    vec.push_back(ip);
  }
}

void
gibbs::dump_wprob(const char *fn){
  // dump "topic wid prob"
  ofstream os;
  os.open(fn);
  if (!os.good()) {
    cerr << "dump_word: cannot open [" << fn << "]\n";
    return;
  }

  vector<idprob> vec;
  int ntopic = num_topics;
  for (int k=0; k<ntopic; k++) {
    get_wprob(k,vec);
    sort(vec.begin(),vec.end(),by_prob);
    int vlen = vec.size();

    if (prm->full_dump) {
      // dump sorted
      for (int i=0; i<vlen; i++) {
	os << k+1 << " " << vec[i].id+1 << " " << vec[i].prob << "\n"; // print in 1-base
      }
    } else {
      // dump top-100, comma separated wid:prob
      string sep="";
      for (int i=0; i<100; i++) {
	os << sep << vec[i].id+1 << ":" << vec[i].prob;
	sep=",";
      }
      os << "\n";
    }
  }
  os.close();
}

void 
gibbs::dump_tprob(const char *fn){
  // dump "docid topic prob"
  ofstream os;
  os.open(fn);
  if (!os.good()) {
    cerr << "dump_word: cannot open [" << fn << "]\n";
    return;
  }

  int ndoc = cp->ndoc;

  vector<idprob> vec;
  for (int j=0; j<ndoc; j++) {
    document *dp=&(cp->doc[j]);
    int docid = dp->docid;
    get_tprob(dp,gcp,vec,false);
    int vlen = vec.size();

    if (prm->full_dump) {
      // dump sorted
      sort(vec.begin(),vec.end(),by_prob);
      for (int i=0; i<vlen; i++) {
	//os << j+1 << " " << vec[i].id+1 << " " << vec[i].prob << "\n"; // print in 1-base
	os << docid << " " << vec[i].id+1 << " " << vec[i].prob << "\n"; // print in 1-base
      }
    } else {
      // no sort, comma separated
      string sep="";
      for (int i=0; i<vlen; i++) {
	os << sep << vec[i].prob; 
	sep=",";
      }
      os << "\n";
    }
  }
  os.close();
}

void
gibbs::show_ppx(int it) {
  if (ppx_every<1) return;
  if (verbose) cerr << " ppx:";

  int ntopic = num_topics;
  int nproc = num_procs;

  // init for lcp->txx_logprob should be done in each thread
  
  // init global p(w|k)
  for (int k=0; k<ntopic; k++) get_wprob(k,p_wk[k]);

  // parallel logprob
  parallel_ppx();

  // collect logprob
  double trn_logprob = 0;
  double tst_logprob = 0;
  int trn_wcnt=0, tst_wcnt=0;

  for (int p=0; p<nproc; p++) {
    counts& cn = lcnts[p];
    trn_logprob += cn.trn_logprob;
    tst_logprob += cn.tst_logprob;
    trn_wcnt += cn.trn_wcnt;
    tst_wcnt += cn.tst_wcnt;
  }

  double trn_avelogprob,tst_avelogprob;
  double trn_ppx=0,tst_ppx=0;

  trn_avelogprob = trn_logprob/trn_wcnt;
  trn_ppx = pow(2,-trn_avelogprob);

  if (tst_wcnt>0) {
    tst_avelogprob = tst_logprob/tst_wcnt;
    tst_ppx = pow(2,-tst_avelogprob);
  }
  if (verbose) cerr << "[trn:" << trn_ppx << ", tst:" << tst_ppx << "]";
  else cerr << "[" << it << " " << trn_ppx << "]";
}

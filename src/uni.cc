#include <iostream>
#include <cmath>

#include "gibbs.h"

//
// This file contains old code for uni-processor mode
//

void
gibbs::iterate_uni(int niter, int ntopic, int nproc) {
  num_iterations = niter;
  num_topics = ntopic;
  num_procs = nproc;
  alloc_counts();		// both global and local counts

  cp->init_zsmp(ntopic);	// init by random
  init_n_wk();			// also init n_k (from zsmp), global only
  int ndoc = cp->ndoc;
  for (int it=0; it<niter; it++) {
    cerr << "iteration " << it << "\n";
    for (int j=0; j<ndoc; j++) {
      document *dp = &(cp->doc[j]);
      init_n_kj(dp,gcp,false);
      int dlen=dp->len;
      for (int i=0; i<dlen; i++) {
	remove_counts(dp,i,gcp);
	update_zsmp(dp,i,gcp);
	add_counts(dp,i,gcp);
      }
    }
  }
}

void
gibbs::show_ppx_uni(){
  int ndoc = cp->ndoc;
  int nvoc = num_vocab;
  int ntopic = num_topics;

  double trn_logprob = 0;
  double tst_logprob = 0;

  if (ppx_every<1) return;

  cerr << " ppx:";

  vector<idprob> *pr_wk = new vector<idprob>[ntopic];
  for (int k=0; k<ntopic; k++) get_wprob(k,pr_wk[k]);

  int trn_wcnt=0, tst_wcnt=0;
  for (int j=0; j<ndoc; j++) {
    bool doppx = takeppx(j);
    if (!doppx) continue;

    document *dp = &(cp->doc[j]);
    vector<idprob> pr_kj;
    get_tprob(dp,gcp,pr_kj,doppx);
    
    int dlen=dp->len;
    for (int i=0; i<dlen; i++) {
      int wid = dp->wid[i];
      int cnt = dp->cnt[i];
      
      // get sum(k)p(w|z)p(z|doc)
      double prob = 0;
      for (int k=0; k<ntopic; k++) {
	vector<idprob>& pr_w = pr_wk[k];
	prob += pr_w[wid].prob * pr_kj[k].prob;
      }
      double logprob = log2(prob);
      logprob *= cnt;
      
      if (doppx && heldout(i)) {
	tst_logprob += logprob; tst_wcnt += cnt;
      } else {
	trn_logprob += logprob; trn_wcnt += cnt;
      }
    }
  }
  double trn_avelogprob = trn_logprob/trn_wcnt;
  double tst_avelogprob = tst_logprob/tst_wcnt;

  double trn_ppx = pow(2,-trn_avelogprob);
  double tst_ppx = pow(2,-tst_avelogprob);

  cerr << "(trn:" << trn_ppx << ", tst:" << tst_ppx << ")";
  delete[] pr_wk;
}


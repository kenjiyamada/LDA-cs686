#include <iostream>
#include <cmath>
#include <pthread.h>

#include "gibbs.h"

void
gibbs::copy_counts(){
  // copy global counts to locals
  int nvoc = num_vocab;
  int ntopic = num_topics;
  int nproc = num_procs;

  // copy n_wk[]
  for (int k=0; k<ntopic; k++)
    for (int i=0; i<nvoc; i++) {
      int cnt = gcp->get_n_wk(i,k);
      for (int p=0; p<nproc; p++) lcnts[p].set_n_wk(i,k,cnt);
    }
  // copy n_k[]
  for (int k=0; k<ntopic; k++) {
    for (int p=0; p<nproc; p++) {
      int cnt = gcp->get_n_k(k);
      //lcnts[p].set_n_kj(k,cnt);  // bug???
      lcnts[p].set_n_k(k,cnt);
    }
  }

}

void
gibbs::collect_counts(){
  int nvoc = num_vocab;
  int ntopic = num_topics;
  int nproc = num_procs;

  // add local delta n_wkp into global
  int len_wk = nvoc*ntopic;
  for (int i=0; i<len_wk; i++)
    for (int p=0; p<nproc; p++)
      gcp->n_wk[i] += lcnts[p].n_wk[i];
  
  // get n_k from n_wk
  for (int k=0; k<ntopic; k++) gcp->n_k[k]=0;
  for (int i=0; i<nvoc; i++)
    for (int k=0; k<ntopic; k++)
      gcp->n_k[k] += gcp->get_n_wk(i,k);

}

inline bool
takeppx(int ppx_every, int i) {
  return(ppx_every==1 || (ppx_every>1 && i%ppx_every==1));
}

inline bool
heldout(int heldout_every, int i) {
  return(heldout_every>1 && i%heldout_every==1);
}

void *
local_sweep(void *arg){
  pinfo *pf = (pinfo*)arg;
  pf->res = pf->pid*10;

  gibbs *gb = pf->gb;
  corpus *cp = gb->cp;

  int ndoc = cp->ndoc;
  int nproc = gb->num_procs;
  int ppx_every = gb->ppx_every;
  int heldout_every = gb->heldout_every;

  int pid = pf->pid;
  counts *lcp = &(gb->lcnts[pid]);

  // for debub
  int xcnt=0;
  int xsum=0;

  int idx=pid;
  int dcnt=0;
  while(idx<ndoc) {
    bool doppx = takeppx(ppx_every,idx);
    document *dp = &(cp->doc[idx]);
    dcnt++;
    gb->init_n_kj(dp,lcp,doppx);
    int dlen = dp->len;
    for (int i=0; i<dlen; i++) {
      if (doppx && heldout(heldout_every,i)) continue;
      gb->remove_counts(dp,i,lcp);
      double xzd = gb->update_zsmp(dp,i,lcp);
      gb->add_counts(dp,i,lcp);
      //if (xcnt<5) pf->res5[xcnt++]=xzd;
      xsum += xzd;		// for debug
    }
    idx += nproc;
  }

  // make n_wkp delta (n_wkp <- n_wkp - n_wk)
  int nvoc = gb->num_vocab;
  int ntopic = gb->num_topics;
  counts *gcp = gb->gcp;
  int len_wk = nvoc*ntopic;
  for (int i=0; i<len_wk; i++)
    lcp->n_wk[i] -= gcp->n_wk[i];

  // put something in res
  //pf->res = dcnt;
  pf->res = xsum;
}

void
gibbs::parallel_sweep(){
  int nproc=num_procs;

  pinfo *parr = new pinfo[nproc];
  for (int p=0; p<nproc; p++) {
    parr[p].pid=p;
    parr[p].gb = this;
    parr[p].res=0;
    pthread_create(&(parr[p].thr), 0, local_sweep, (void*)&(parr[p]));
  }
  
  for (int p=0; p<nproc; p++)
    pthread_join(parr[p].thr,0);

  if (verbose) 
    for (int p=0; p<nproc; p++)
      cerr << "[" << parr[p].res << "]";

  /***
  for (int p=0; p<nproc; p++)
    //for (int xc=0; xc<5; xc++) cerr << "[" << parr[p].res5[xc] << "]";
    cerr << "[" << parr[p].res5[0] << "]";
  ***/

  delete[] parr;

  //show_ppx();
  //cerr << " sweep_done\n";

}

//
// ppx
//

void *
local_ppx(void *arg){
  pinfo *pf = (pinfo*)arg;
  pf->res = pf->pid*5;

  gibbs *gb = pf->gb;
  corpus *cp = gb->cp;

  int ndoc = cp->ndoc;
  int ntopic = gb->num_topics;
  int nproc = gb->num_procs;
  int ppx_every = gb->ppx_every;
  int heldout_every = gb->heldout_every;

  int pid = pf->pid;
  counts *lcp = &(gb->lcnts[pid]);

  // init 
  double trn_logprob = 0;
  double tst_logprob = 0;
  int trn_wcnt=0, tst_wcnt=0;

  int idx=pid;
  while(idx<ndoc) {
    bool doppx = takeppx(ppx_every,idx);
    if (!doppx) continue;

    document *dp = &(cp->doc[idx]);
    //vector<idprob> pr_kj;
    //gb->get_tprob(dp,lcp,pr_kj,doppx);
    gb->get_tprob(dp,lcp,doppx);
    float *p_kj = lcp->p_kj;

    int dlen=dp->len;
    for (int i=0; i<dlen; i++) {
      int wid = dp->wid[i];
      int cnt = dp->cnt[i];
      
      // get sum(k)p(w|z)p(z|doc)
      double prob = 0;
      for (int k=0; k<ntopic; k++) {
	vector<idprob>& pr_w = gb->p_wk[k];
	//prob += pr_w[wid].prob * pr_kj[k].prob;
	prob += pr_w[wid].prob * p_kj[k];
      }
      double logprob = log2(prob);
      logprob *= cnt;
      
      if (doppx && heldout(heldout_every,i)) {
	tst_logprob += logprob; tst_wcnt += cnt;
      } else {
	trn_logprob += logprob; trn_wcnt += cnt;
      }
    }
    idx += nproc;
  }

  // update lcp
  lcp->trn_logprob=trn_logprob;
  lcp->tst_logprob=tst_logprob;
  lcp->trn_wcnt=trn_wcnt;
  lcp->tst_wcnt=tst_wcnt;

}

void
gibbs::parallel_ppx() {

  int nproc=num_procs;

  pinfo *parr = new pinfo[nproc];
  for (int p=0; p<nproc; p++) {
    parr[p].pid=p;
    parr[p].gb = this;
    parr[p].res=0;
    pthread_create(&(parr[p].thr), 0, local_ppx, (void*)&(parr[p]));
  }
  
  for (int p=0; p<nproc; p++)
    pthread_join(parr[p].thr,0);

  /****
  for (int p=0; p<nproc; p++)
    cerr << "(" << parr[p].res << ")";
  ****/

  delete[] parr;
}

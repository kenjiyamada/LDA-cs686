#include <iostream>
#include <sstream>
#include <fstream>

#include "gibbs.h"

//
// Obtain topic distribution P(T|D) for new test documents
//
// Usage:
// % ./lda-2.0 -fd -ot tpc.txt -od dtc.txt mydocs.txt 100 10
// % ./lda-2.0 -fd -it tpc.txt -ot /dev/null -od dtc-new.txt newdoc.txt 10 10
//

void
corpus::load_tpc() {
  // init
  p_wk = new float[num_topics*nvoc];

  // load file
  const char *fn=prm->fn_tpc_in;
  ifstream ifs(fn);
  if (!ifs.good()) {
    cerr << "can't open TPC [" << fn << "]\n";
    exit(1);

  } else {
    cerr << "loading TPC [" << fn << "]...";
  }

  string line;
  int lineno=0;
  while(!ifs.eof()) {
    // must be in full-dump format
    // format: topic wid prob (1-base) 
    getline(ifs,line);
    if (line.length()==0) continue;
    lineno++;
    istringstream cs(line);

    int topic,wid;
    float prob;
    cs >> topic >> wid >> prob;

    if (topic<1 || topic>num_topics) error_exit("topid_id out of range",lineno,line);
    if (wid<1 || wid>nvoc) error_exit("word_id out of range",lineno,line);
    if (prob<0 || prob>1) error_exit("prob out of range",lineno,line);

    set_p_wk(wid-1,topic-1,prob);	// 1-base
  }
  cerr << "done (" << lineno << " lines)\n";

  if (lineno==0) {
    cerr << "TPC file is empty\n";
    exit(1);
  }
}

void
gibbs::test_new_doc() {
  // follows gibbs::iterate_uni()  while fixing P(W|T)

  int ntopic = num_topics;
  int niter = num_iterations;

  cp->load_tpc();

  alloc_counts();
  cp->init_zsmp();		// init by random
  init_n_wk();			// also init n_k (from zsmp)

  // gibbs sampling
  int ndoc = cp->ndoc;
  for (int it=0; it<niter; it++) {
    cerr << "[" << it << "]";
    for (int j=0; j<ndoc; j++) {
      document *dp = &(cp->doc[j]);
      init_n_kj(dp,gcp,false);
      int dlen=dp->len;
      for (int i=0; i<dlen; i++) {
	remove_counts(dp,i,gcp);
	update_zsmp_ntst(dp,i,gcp);
	add_counts(dp,i,gcp);
      }
    }
  }
  cerr << "\n";
}

void
gibbs::update_zsmp_ntst(document *dp, int i, counts *lcp) {
  int wid = dp->wid[i];

  // init zprb array
  int ntopic = num_topics;
  float *zprb = lcp->zprb;
  for (int k=0; k<ntopic; k++) zprb[k]=0;

  // obtain estimated zprb
  int nvoc = num_vocab;
  float cum=0;
  for (int k=0; k<ntopic; k++) {
    // P(T|D) = (alpha + n_kj(k)) / (ntopic*alpha + dlen)
    // though we ignore the denominator, since we don't need normalization.
    float prb = (alpha + lcp->get_n_kj(k)) * cp->get_p_wk(wid,k);
    cum += prb;
    zprb[k] = cum;
  }
  
  // sample topic
  float pr = lcp->rnd->get_rand();
  dp->sample_zsmp(i,lcp->zprb,ntopic,pr);
}


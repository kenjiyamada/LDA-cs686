#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "corpus.h"

inline void
triplet::read(const string& s) {
  docid=wordid=count=-1;
  istringstream cs(s);
  cs >> docid >> wordid >> count;
}

inline bool
triplet::good() {
  return (docid>0 && wordid>0 && count>0);
}

void
document::set_words(vector<triplet> *vec, int ntopic){
  len = vec->size();
  if (len==0) {
    cerr << "empty for docid=" << docid << "\n";
    return;
  }

  triplet tp0 = (*vec)[0];
  docid = tp0.docid;
  if (docid==0) {
    cerr << "unexpected docid==0\n";
    exit(1);
  }

  wid = new int[len];
#if USE_UCHAR_FOR_WID
  cnt = new unsigned char[len];
#else
  cnt = new int[len];
#endif
  zsmp = new int[len];

  for (int i=0; i<len; i++) {
    triplet tp = (*vec)[i];

    // 1-base to 0-base for wid
    int wordid = tp.wordid;
    if (wordid==0) {
      cerr << "unexpected wordid==0 \n";
      exit(1);
    }
    wordid--;			// make 0-base

    wid[i]=wordid;
#if USE_UCHAR_FOR_WID
    cnt[i]=(tp.count>255 ? 255 : tp.count); // must fit in 8 bits (unsigned char)
#else
    cnt[i]=tp.count;
#endif
  }
}

corpus::corpus(param *pm) {
  prm = pm;
  num_topics = pm->num_topics;	// copy params

  rnd = new random_nr();
  rnd->set_seed(1234);
}

void 
corpus::error_exit(string msg, int lineno, string line) {
  cerr << "\nerror: " << msg << " (line " << lineno << ") [" << line << "]\n";
  exit(1);
}

int 
corpus::s2i(const string& s) {
  int i;
  istringstream cs(s);
  cs >> i;
  return i;
}

void
corpus::load(char *fn, int ntopic){
  ifstream ifs(fn);
  if (!ifs.good()) {
    cerr << "can't open fn [" << fn << "]\n";
    exit(1);

  } else {
    cerr << "loading [" << fn << "]:";

    // read headers
    string line;
    getline(ifs,line);
    ndoc = s2i(line);

    getline(ifs,line);
    nvoc = s2i(line);

    getline(ifs,line);
    int nword = s2i(line);

    if (prm->verbose>0) cerr << "ndoc=" << ndoc << ", nvoc=" << nvoc << ", nword=" << nword << "\n";

    // alloc arrays
    doc = new document[ndoc];

    // load the content
    int docIdx = 0;
    vector<triplet> tpvec;
    int previd=-1;
    int lineno=0;
    while(!ifs.eof()) {
      getline(ifs,line);
      if (line.size()>0) {
	triplet tp;
	tp.read(line);
	if (tp.good()) {
	  lineno++;
	  if (prm->verbose>0 && lineno<5)
	    cerr << lineno << ": docid=" << tp.docid << ", wordid=" << tp.wordid << ", count=" << tp.count << "\n";

	  if (tp.docid < previd) {
	    cerr << "docid decremented: " << previd << " -> " << tp.docid << "\n";
	  } else if (tp.docid > previd) {
	    // need to flush
	    if (tpvec.size()>0) {
	      if (docIdx%10000==0) cerr << docIdx/10000 << ".";
	      doc[docIdx].set_words(&tpvec,ntopic);
	      docIdx++;
	      tpvec.clear();
	    }
	    previd = tp.docid;
	  }
	  tpvec.push_back(tp);
	}
      }
    } // end of while
    if (tpvec.size()>0) doc[docIdx].set_words(&tpvec,ntopic);

    if (prm->verbose>0) cerr << "tpvec.size()=" << tpvec.size() << "\n";
    cerr << "done (" << docIdx+1 << " docs, " << lineno << " lines)\n";
  }
}

void
corpus::show_stats(){
  cerr << "ndoc=" << ndoc << "\n";

  int sum = 0;
  int maxid = 0;
  int maxid_doc = 0;

  for (int i=0; i<ndoc; i++) {
    document *dp = &(doc[i]);
    int len = dp->len;
    if (len<1) continue;
    int lastwd = dp->wid[len-1];

    sum += len;
    if (lastwd > maxid) {
      maxid = lastwd;
      maxid_doc = dp->docid;
    }
  }

  cerr << "sum=" << sum << "\n";
  cerr << "maxid=" << maxid << " at doc=" << maxid_doc << "\n";
}

void
corpus::init_zsmp(){
  int ntopic = num_topics;	// using copied class var
  long xsum = 0;
  for (int j=0; j<ndoc; j++) {
    document *dp=&(doc[j]);
    int dlen = dp->len;
    for (int i=0; i<dlen; i++) {
      //double pr = (double)rand()/RAND_MAX;
      float pr = rnd->get_rand();
      int tpc = pr * ntopic;	// implicit truncation
      dp->zsmp[i] = tpc;
      xsum += tpc;
    }
  }
  if (prm->verbose>0) cerr << "init_zsmp(): " << xsum << "\n";
}

double
document::sample_zsmp(int i, float *zprb, int ntopic, float pr) {
  // originally defined within gibbs::update_zsmp().
  // note zprb has cumlative probs.

  float cum = zprb[ntopic-1];
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
  //dp->zsmp[i] = kx;
  zsmp[i] = kx;

  // return debug value
  double zret = kx;
  return zret;

}

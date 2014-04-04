#include <iostream>
#include <cstdlib>

#include "param.h"
#include "lda.h"
#include "rand.h"

using namespace std;

void
check_size() {
  int szint = sizeof(int);
  int szlong = sizeof(long);
  cerr << "sizeof(int)=" << szint << ", sizeof(long)=" << szlong << "\n";

  if (szint<4 || szlong<8) {
    cerr << "sizeof(int) or sizeof(long) too small\n";
    exit(1);
  }
}

void
check_random(random_nr *rnd) {
  int ntopic=20;
  for (int i=0; i<10; i++) {
    float pr = rnd->get_rand();
    int tpc = pr * ntopic;	// implicit truncation
    cout << "[" << pr << ":" << tpc << "]";
  }
  cout << "\n";
}

int 
main(int argc, char** argv) {
  check_size();

  param *prm = new param(argc,argv);

  if (prm->verbose>0) {
    random_nr *rn = new random_nr();
    rn->set_seed(1234);
    check_random(rn);
  }

  lda *lx = new lda(prm);
  lx->load();
  lx->run();
}


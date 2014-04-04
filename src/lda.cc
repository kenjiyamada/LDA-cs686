#include <iostream>
#include <sstream>
#include <cstdlib>

#include "lda.h"

lda::lda(param *pm) {
  cp = new corpus(pm);
  prm = pm;
}

lda::~lda() {
}

void
lda::load(){
  cp->load(prm->doc_fn,prm->num_topics);
  if (prm->verbose>0) cp->show_stats();
}

void
lda::run(){
  gb=new gibbs(prm,cp);
  gb->iterate();

  gb->dump_wprob("topics.txt");
  gb->dump_tprob("doctopic.txt");
}

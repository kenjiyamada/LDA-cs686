#include <iostream>
#include <sstream>
#include <cstdlib>

#include "param.h"

param::param(int argc, char** argv) {
  process_command(argc,argv);
}

bool
param::str_eq(const char *s1, const char *s2) {
  char c1 = *s1;
  char c2 = *s2;
  while(c1>0 && c2>0) {
    if (c1!=c2) return false;
    c1 = *++s1;
    c2 = *++s2;
  }
  return(c1==0 && c2==0);
}

int 
param::s2i(const string& s) {
  int i;
  istringstream cs(s);
  cs >> i;
  return i;
}

double
param::s2d(const string& s) {
  double d;
  istringstream cs(s);
  cs >> d;
  return d;
}

void
param::process_command(int argc, char** argv) {

  if (argc < 4) {
    cerr << "need 3 args\n";
    exit(1);
  }

  alpha=0.1; beta=0.1;		// alpha=0.5 may work better?
  verbose=0;
  full_dump=false;

  num_procs=4;
  ppx_every=1;
  heldout_every=0;

  fn_topic="topics.txt";
  fn_doctopic="doctopic.txt";

  fn_NULL = "";
  fn_tpc_in=fn_dtc_in=fn_NULL;

  int i=1;
  while(*argv[i]=='-') {
    char *arg = argv[i];
    if (str_eq(arg,"-p")) num_procs = s2i(argv[++i]);
    else if (str_eq(arg,"-a")) alpha = s2d(argv[++i]);
    else if (str_eq(arg,"-b")) beta = s2d(argv[++i]);
    else if (str_eq(arg,"-v")) verbose++;
    else if (str_eq(arg,"-fd")) full_dump=true;
    else if (str_eq(arg,"-px")) ppx_every = s2i(argv[++i]);
    else if (str_eq(arg,"-ho")) heldout_every = s2i(argv[++i]);

    else if (str_eq(arg,"-ot")) fn_topic = argv[++i];
    else if (str_eq(arg,"-od")) fn_doctopic = argv[++i];
    else if (str_eq(arg,"-it")) fn_tpc_in = argv[++i];
    else if (str_eq(arg,"-id")) fn_dtc_in = argv[++i];

    else {
      cerr << "unrecognized option: " << arg << "\n";
      exit(1);
    }
    i++;
  }
  doc_fn = argv[i++];
  num_iterations = s2i(argv[i++]);
  num_topics = s2i(argv[i++]);

  cerr << "num_iterations=" << num_iterations << ", num_topics=" << num_topics
       << ", num_procs=" << num_procs 
       << ", alpha:beta=" << alpha << ":" << beta 
       << ", ppx_every=" << ppx_every << ", heldout_every=" << heldout_every << "\n";
}

bool
param::test_new_doc() {
  return (fn_tpc_in != fn_NULL && fn_dtc_in == fn_NULL);
}

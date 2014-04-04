#ifndef PARAM_H
#define PARAM_H

#include <string>

using namespace std;

class param {
 public:
  char *doc_fn;
  int num_iterations;
  int num_topics;
  int num_procs;
  int ppx_every;
  int heldout_every;

  int verbose;
  bool full_dump;
  float alpha,beta;

  param(int argc, char** argv);
  void process_command(int argc, char** argv);

 private:
  bool str_eq(const char *s1, const char *s2);
  int s2i(const string& s);
  double s2d(const string& s);

};

#endif

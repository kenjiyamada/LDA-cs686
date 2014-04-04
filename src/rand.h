#ifndef RAND_H
#define RAND_H

using namespace std;

class random_nr {
  int seed;
 public:
  inline void set_seed(int i);
  inline float get_rand();
};
    
inline void 
random_nr::set_seed(int i) {
    seed=i;
}

inline float
random_nr::get_rand() {
  // from Numerical Recipes in C++, p283
  const int IA=16807, IM=2147483647, IQ=127773;
  const int IR=2836, MASK=123459876;
  const float AM=1.0/float(IM);
  int k;
  float ans;

  seed ^= MASK;
  k=seed/IQ;
  seed=IA*(seed-k*IQ)-IR*k;
  if (seed < 0) seed += IM;
  ans=AM*seed;
  seed ^= MASK;
  return ans;
}

#endif

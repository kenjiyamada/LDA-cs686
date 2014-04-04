It implements AD-LDA [Newman, et. al. NIPS-07], in C++ using POSIX thread.
No extra library is needed, but it needs -pthread option for compile.
The Makefile takes care of the flags.

To compile:
    % cd src
    % make clean
    % make

To run:
    % cp somewhere/docword.nytimes.txt ./
    % make run
or
    % make run-nips

Running with the NY Times data for 1,000 iterations (100 topics, 4-CPUs) took 19 hours on
a shared Amazon AWS server. The CPU info from /proc/cpuinfo says Intel Xeon CPU E5-2650 0
@ 2.00GHz, with 4 processors.

The training perplexity at each 100 iteration is the following:

    iteration   training_perplexity
  ---------------------------------
    100         3674.39
    200         3601.38
    300         3578.41
    400         3566.82
    500         3562.33
    600         3560.55
    700         3558.67
    800         3557.12
    900 	3555.33
   1000		3554.10

For the perplexity at each iteration, see LDA-nyt-ppx.log.


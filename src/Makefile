
# main executable
MAIN=lda-2.0

# object files
OBJX = uni.o multi.o corpus.o gibbs.o ntst.o lda.o param.o main.o 

# src files
SRCS = README.txt Makefile main.cc rand.h param.h param.cc lda.h lda.cc corpus.h corpus.cc gibbs.h gibbs.cc multi.cc uni.cc ntst.cc
TGZ = lda-2.0.tgz

# output files
OUTS = topics.txt doctopic.txt

CFLAGS = -O3 -pthread 

ECHO = /bin/echo

all : $(MAIN)

%.o : %.cc
	g++ $(CFLAGS) -c -o $@ $<

main.o : main.cc lda.h rand.h param.h

param.o : param.cc param.h

lda.o : lda.cc lda.h gibbs.h corpus.h rand.h param.h

corpus.o : corpus.cc corpus.h rand.h

gibbs.o : gibbs.cc gibbs.h corpus.h rand.h param.h

multi.o : multi.cc gibbs.h corpus.h rand.h param.h

uni.o : uni.cc gibbs.h corpus.h rand.h param.h

ntst.o : ntst.cc gibbs.h corpus.h param.h

$(MAIN) : $(OBJX)
	g++ $(CFLAGS) -o $(MAIN) $(OBJX)

# commands

run : run-nyt

run-nyt : $(MAIN)
	./$(MAIN) docword.nytimes.txt 100 100

run-nips : $(MAIN)
	./$(MAIN) docword.nips.txt 100 100

run-nyt5 : $(MAIN)
	./$(MAIN) docword.nytimes.txt 5 100

run-nips5 : $(MAIN)
	./$(MAIN) docword.nips.txt 5 100

run-verbose : $(MAIN)
	./$(MAIN) -v docword.nips.txt 5 100

run-long : $(MAIN)
	./$(MAIN) docword.nips.txt 50 100

run-ppx : $(MAIN)
	./$(MAIN) -p 4 -px 1 -ho 10 docword.nips.txt 50 100

run-nyt-ppx : $(MAIN)
	./$(MAIN) -p 4 -px 1 -ho 0 docword.nytimes.txt 5 100

run-nyt-test-ppx : $(MAIN)
	./$(MAIN) -p 4 -px 1 -ho 10 docword.nytimes.txt 10 100

run-nips-ppx : $(MAIN)
	./$(MAIN) -p 4 -px 1 -ho 0 docword.nips.txt 10 100

run-nips-ab : $(MAIN)
	./$(MAIN) -a 0.1 -b 0.1 -p 4 -px 1 -ho 0 docword.nips.txt 100 100

run-debug : 
	@echo "Args:" docword.nips.txt 5 20
	@gdb ./$(MAIN) 

debug :
	make CFLAGS="-g -DDEBUG"

no-opt :
	make CFLAGS=-pthread

tar :
	tar cvfz $(TGZ) $(SRCS) LDA-nyt-ppx.log

clean :
	rm -f *~ *.o $(MAIN) gmon.out $(OUTS) $(TGZ)

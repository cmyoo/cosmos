# Makefile

LANG=cpp
# program file name
PROG = cosmos

# the number of CPUs
#NCPUS = 1

# compiler
CXX = g++
#CXX = g++ -pg
#CXX = icpc

# compiler option
#CFLAGS = -Wall -O2
#CFLAGS = -Wall -O3 -fopenmp
#CFLAGS = -Wall -march=corei7-avx -mtune=corei7-avx -m64 -O3 -ffast-math -finline-functions -funroll-loops -fopenmp
CFLAGS = -Wall -march=native -mtune=native -m64 -O3 -ffast-math -finline-functions -funroll-loops -fopenmp -std=gnu++0x 
#CFLAGS = -Wall -march=native -mtune=native -m64 -O3 -ffast-math -finline-functions -funroll-loops -fopenmp -std=c++11
#CFLAGS = -O3 -ipo -openmp

# library link
# LIBS= -lm

# file dependence
DEP = .dep.inc
#OUTFILE=*.dat data.log

# suffix rule
.SUFFIXES: .cpp .o

# source file
SRC = $(PROG).cpp ../cosmos_bssn.cpp ../cosmos_initial.cpp ../cosmos_output.cpp ../cosmos_boundary.cpp ../cosmos_ahf.cpp ../cosmos_ipol.cpp ../cosmos_fluid.cpp ../cosmos_fmr.cpp 
OBJS = $(SRC:%.$(LANG)=%.o)

.PHONY: all
all: depend $(PROG)

# 1st rule
$(PROG): $(OBJS)
	$(CXX) $(CFLAGS) -o $(PROG) $(OBJS)

.$(LANG).o:
	$(CXX) $(CFLAGS) -o $@ -c  $<

# clean rule
.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) $(DEP) 
#$(OUTFILE)

# depend rule ( header file )
.PHONY: depend
depend: $(SRC)
	@- $(RM) $(DEP);
	@ for LIS in $^; \
	do ($(CXX) $(CFLAGS) -MM $$LIS \
	| sed "s/\ [_a-zA-Z0-9][_a-zA-Z0-9]*\.$(LANG)//g" >> $(DEP);); \
	done

# header file dependence
-include $(DEP)


BIN    = -L/opt/local/bin
LIB    = -L/opt/local/lib
INCL   = -I/opt/local/include -I./../src
SRCDIR = ../src/
EXEDIR = ../exe/
TSTDIR = ../tst/
HEADAS = /usr/local/lib/heasoft/heasoft-6.28/x86_64-apple-darwin19.6.0/

#CC     = /opt/local/bin/clang-mp-9.0 #/opt/local/bin/gcc-mp-10
#CXX    = /opt/local/bin/clang++-mp-9.0 #/opt/local/bin/g++-mp-10
CC     = /opt/local/bin/gcc-mp-10  #/opt/local/bin/clang-mp-9.0 #/opt/local/bin/gcc-mp-10
CXX    = /opt/local/bin/g++-mp-10 #/opt/local/bin/clang++-mp-9.0 #/opt/local/bin/g++-mp-10
FC     = /opt/local/bin/gfortran-mp-10
OPTS   = -std=c++17 -O2 -g -Wall -DUSE_APED

OBJ = aped_spectrum.o Aped.o FitsUtil.o
OBJTST = test_apecode_vs_heasoft.o Aped.o FitsUtil.o

LIBS   = $(LIB) -lstdc++ -lm -lcfitsio
HEADASINCL= -I$(HEADAS)include
HEADASLIB= -L$(HEADAS)lib -lXSFunctions -lXSUtil


.cpp.o:
	$(CXX) -c $(INCL) $(OPTS) $<

.f.o:
	$(FC) $(OPTS) -c $<

.c.o:
	$(CC) $(OPTS) -c $<

aped: $(OBJ)
	$(CXX) ${OPTS} -o $(EXEDIR)aped.exe $(OBJ) $(LIBS)

test: $(OBJTST)
	$(CXX) ${OPTS} -o $(EXEDIR)test_apecode_vs_heasoft.exe $(OBJTST) $(LIBS) $(HEADASLIB)

clean: 
	\rm *.o $(EXEDIR)*.exe

#------------------------------------------------------------------------------
# list of dependencies
#------------------------------------------------------------------------------

aped_spectrum.o: $(EXEDIR)aped_spectrum.cpp
	$(CXX) $(OPTS) $(INCL) -c $(EXEDIR)aped_spectrum.cpp

test_apecode_vs_heasoft.o: $(TSTDIR)test_apecode_vs_heasoft.cpp
	$(CXX) $(OPTS) $(HEADASINCL) $(INCL) -c $(TSTDIR)test_apecode_vs_heasoft.cpp

Aped.o: $(SRCDIR)Aped.cpp
	$(CXX) $(OPTS) $(INCL) -c $(SRCDIR)Aped.cpp

FitsUtil.o: $(SRCDIR)FitsUtil.cpp
	$(CXX) $(OPTS) $(INCL) -c $(SRCDIR)FitsUtil.cpp

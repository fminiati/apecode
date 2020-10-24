BIN    = -L/opt/local/bin
LIB    = -L/opt/local/lib
INCL   = -I/opt/local/include -I./../src -I./../exe
SRCDIR = ../src/
EXEDIR = ../exe/
TSTDIR = ../tst/
HEADAS = /usr/local/lib/heasoft/heasoft-6.28/x86_64-apple-darwin19.6.0/

#CC     = /opt/local/bin/clang-mp-10 #/opt/local/bin/gcc-mp-10
#CXX    = /opt/local/bin/clang++-mp-10 #/opt/local/bin/g++-mp-10
CC     = /opt/local/bin/gcc-mp-10  #/opt/local/bin/clang-mp-9.0 #/opt/local/bin/gcc-mp-10
CXX    = /opt/local/bin/g++-mp-10 #/opt/local/bin/clang++-mp-9.0 #/opt/local/bin/g++-mp-10
#FC     = /opt/local/bin/gfortran-mp-10
LD     = /opt/local/bin/gcc-mp-10
CXXFLAGS   = -std=c++17 -O2 -Wall -DUSE_APED # -DUSE_TIMER=5 -g

OBJ = aped_spectrum.o Aped.o
TOBJ = test_apecode_vs_heasoft.o Aped.o #XspecAped.o
KOBJ = test_kernel.o

LIBS   = $(LIB) -lstdc++ -lm
FITSLIBS = $(LIB) -lcfitsio
HEADASINCL= -I$(HEADAS)include
HEADASLIB= $(FITSLIBS) -L$(HEADAS)lib -lXSFunctions -lXSUtil # -lCCfits_2.5 -lXS


.cpp.o:
	$(CXX) -c $(INCL) $(OPTS) $<

.f.o:
	$(FC) $(CXXFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

aped: $(OBJ)
	$(CXX) ${CXXFLAGS} -o $(EXEDIR)aped.exe $(OBJ) $(LIBS) $(FITSLIBS)

test: $(TOBJ)
	$(CXX) ${CXXFLAGS} -o $(EXEDIR)test_apecode_vs_heasoft.exe $(TOBJ) $(LIBS) $(HEADASLIB)

kernel: $(KOBJ)
	$(CXX) ${CXXFLAGS} -o $(EXEDIR)test_kernel.exe $(KOBJ) $(LIBS)

clean: 
	\rm -f *.o $(EXEDIR)*.exe

#------------------------------------------------------------------------------
# list of dependencies
#------------------------------------------------------------------------------

aped_spectrum.o: $(EXEDIR)aped_spectrum.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c $(EXEDIR)aped_spectrum.cpp

test_apecode_vs_heasoft.o: $(TSTDIR)test_apecode_vs_heasoft.cpp
	$(CXX) $(CXXFLAGS) $(HEADASINCL) $(INCL) -c $(TSTDIR)test_apecode_vs_heasoft.cpp

test_kernel.o: $(TSTDIR)test_kernel.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c $(TSTDIR)test_kernel.cpp

Aped.o: $(SRCDIR)Aped.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c $(SRCDIR)Aped.cpp

#XspecAped.o: $(SRCDIR)XspecAped.cxx
#	$(CXX) $(CXXFLAGS) $(INCL) $(HEADASINCL) -c $(SRCDIR)XspecAped.cxx

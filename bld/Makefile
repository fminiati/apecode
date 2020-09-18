

BIN    = -L/opt/local/bin
LIB    = -L/opt/local/lib
INCL   = -I/opt/local/include -I./../src -I../../util/
EXEDIR = ../exe/

OBJ = aped_spectrum.o Aped.o FitsUtil.o

CC     = /opt/local/bin/gcc-mp-10
CXX    = /opt/local/bin/g++-mp-10
FC     = /opt/local/bin/gfortran-mp-10
OPTS   = -std=c++17 -O2 -g -Wall -DUSE_APED

LIBS   = $(LIB) -lm -lstdc++ -lcfitsio

.cpp.o:
	$(CXX) -c $(INCL) $(OPTS) $<

.f.o:
	$(FC) $(OPTS) -c $<

.c.o:
	$(CC) $(OPTS) -c $<

ox: $(OBJ)
	$(CXX) ${OPTS} -o $(EXEDIR)aped.ex $(OBJ) $(LIBS)

clean: 
	\rm *.o *.ex

#------------------------------------------------------------------------------
# list of dependencies
#------------------------------------------------------------------------------

aped_spectrum.o: ../exe/aped_spectrum.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c ../exe/aped_spectrum.cpp

Aped.o: ../src/Aped.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c ../src/Aped.cpp

FitsUtil.o: ../src/FitsUtil.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c ../src/FitsUtil.cpp

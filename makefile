#C++ compiler
CXX=g++

#Flags for linker
LD_FLAGS= -w2

#Flags for compiler
CXX_FLAGS=-diag-disable=remark -w2 -O2 


###########################################################


OBS=src/ldblock.o src/data.o src/correlations.o src/splitter.o src/output.o

ldblock: $(OBS) 
	$(CXX) $(LD_FLAGS) -o ldblock $(OBS)

%.o:	%.cpp %.h src/global.h
	$(CXX) $(CXX_FLAGS) -c $*.cpp -o $*.o


src/ldblock.o: src/global.h src/data.h src/correlations.h src/splitter.h src/output.h
src/correlations.o: src/data.h
src/splitter.o: src/data.h src/correlations.h
src/output.h: src/data.h src/splitter.h src/correlations.h
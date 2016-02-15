CC=icpc
EXE=SimLYM
SRC=SimLYM.cpp
FLAG=-std=c++0x
INC=../inc
INC_MKL=/usr/local/intel/composer_xe_2013.0.079/mkl/include/
#LIBS_PATH=../lib
LIBS_MKL_PATH=/usr/local/intel/composer_xe_2013.0.079/mkl/lib/intel64/
LIBS=../lib/libsimlym.a
LIBS_MKL=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

all:
	$(CC) $(FLAG) -o $(EXE) $(SRC) -I $(INC) -I $(INC_MKL) -L $(LIBS_MKL_PATH) $(LIBS) $(LIBS_MKL) 

clean:
	rm $(EXE)


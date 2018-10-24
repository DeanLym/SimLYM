CC=g++
STD_FLAG=-std=c++0x
MKL_LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential
OTHER_LIBS=-lpthread -lm
INCLUDE_PATH=../inc
SHARED_LIB=libsimlym.so
STATIC_LIB=libsimlym.a
OBJ=SimCtrl.o Interface_Pardiso.o CMatrixBlock.o CPardisoSolver.o CSchedule.o CStandardWell.o CSTDProdWell.o CSTDWInjWell.o CState.o CHistoryMatching.o

all: SimCtrl.so

SimCtrl.so: 
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c SimCtrl.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c Interface_Pardiso.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CMatrixBlock.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CPardisoSolver.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CSchedule.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CStandardWell.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CState.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CSTDProdWell.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CSTDWInjWell.cpp
	$(CC) $(STD_FLAG) -I $(INCLUDE_PATH) -fPIC -c CHistoryMatching.cpp
	
	$(CC) $(STD_FLAG) -shared -Wl,-soname,$(SHARED_LIB) -o $(SHARED_LIB) $(OBJ)
	ar -cvq $(STATIC_LIB) $(OBJ)
	
	mv $(SHARED_LIB) ../lib
	mv $(STATIC_LIB) ../lib

clean:
	rm *.o
	rm ../lib/$(SHARED_LIB)
	rm ../lib/$(STATIC_LIB)

#icpc -std=c++0x SimLYM.cpp SimCtrl.cpp Interface_Pardiso.cpp CMatrixBlock.cpp CPardisoSolver.cpp CSchedule.cpp CStandardWell.cpp CState.cpp CSTDProdWell.cpp CSTDWInjWell.cpp -o SimLYM -I ../mkl/include -I ../include -L../mkl/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm

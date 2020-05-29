LIBS := -lsimlym -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl -lm

INC1 = ../../inc/
INC2 = ../../inc/hdf5/

LINK_PATH = ../../lib/

all:
	gcc -I$(INC1) -I$(INC2) -L$(LINK_PATH) $(LIBS) -o nw_ms_bhp main.cpp

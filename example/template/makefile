LIBS := -lsimlym -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl -lm

INC1 = /data3/Astro/global/simlym/include/
INC2 = /data3/Astro/global/simlym/include/hdf5/

LINK_PATH = /data3/Astro/global/simlym/lib/

all:
	icpc -I$(INC1) -I$(INC2) -L$(LINK_PATH) $(LIBS) -o nw_ms_bhp main.cpp

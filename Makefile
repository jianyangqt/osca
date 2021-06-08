
# -----------------------------------------------------------------
#   Makefile for osca
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------
EIGEN_PATH = /usr/local/include
MKL_INCLUDE = /home/fanghl/.local/opt/mkl_v2021.2.0_d20210524/include
MKL_LIB = /home/fanghl/.local/opt/mkl_v2021.2.0_d20210524/lib/intel64


CXX = g++
CXXFLAGS= -Wall -O3 -fopenmp
CPPFLAGS=  -I$(EIGEN_PATH) -I$(MKL_INCLUDE)
LDFLAGS= -L$(MKL_LIB)
LIBS= -lz -lgomp -lmkl_core -lpthread -lmkl_gnu_thread -lmkl_gf_lp64


all: osca

osca: dcdflib.o l0_com.o l0_io.o l0_mem.o l0_stat.o \
	l1_op_geno.o l1_stat.o \
	l2_besd.o l2_bfile.o l2_efile.o l2_enet.o l2_reml.o \
	l3_efile.o l3_ewas.o l3_glmnet.o l3_gwas.o l3_smr.o l3_vqtl.o \
	l4_osc.o

	$(CXX) $(LDFLAGS) $(LIBS) *.o -o $@


dcdflib.o: CPP/dcdflib.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l0_com.o: CPP/l0_com.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l0_io.o: CPP/l0_io.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l0_mem.o: CPP/l0_mem.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<
	
l0_stat.o: CPP/l0_stat.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l1_op_geno.o: CPP/l1_op_geno.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l1_stat.o: CPP/l1_stat.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l2_besd.o: CPP/l2_besd.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l2_bfile.o: CPP/l2_bfile.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l2_efile.o: CPP/l2_efile.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l2_enet.o: CPP/l2_enet.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l2_reml.o: CPP/l2_reml.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l3_efile.o: CPP/l3_efile.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l3_ewas.o: CPP/l3_ewas.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l3_glmnet.o: CPP/l3_glmnet.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l3_gwas.o: CPP/l3_gwas.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l3_vqtl.o: CPP/l3_vqtl.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

l3_smr.o: CPP/l3_smr.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<
	
l4_osc.o: CPP/l4_osc.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<


clean:
	rm *.o osca

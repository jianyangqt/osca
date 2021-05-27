
# -----------------------------------------------------------------
#   Makefile for osca
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------
EIGEN_PATH = /usr/local/include
MKL_INCLUDE = /home/fanghl/.local/opt/mkl_v2021.2.0_d20210524/include
MKL_LIB = /home/fanghl/.local/opt/mkl_v2021.2.0_d20210524/lib/intel64


CXX = g++
CXXFLAGS= -Wall -g -O0 -fopenmp
CPPFLAGS= -I$(EIGEN_PATH) -I$(MKL_INCLUDE)
LDFLAGS= -L$(MKL_LIB)
LIBS= -lz -lgomp -lmkl_core -lpthread -lmkl_gnu_thread -lmkl_gf_lp64


all: osca

osca: dcdflib.o l0_com.o l0_io.o l0_mem.o l0_stat.o \
	l1_op_geno.o l1_stat.o \
	l2_besd.o l2_bfile.o l2_efile.o l2_enet.o l2_reml.o \
	l3_efile.o l3_ewas.o l3_glmnet.o l3_gwas.o l3_smr.o l3_vqtl.o \
	l4_osc.o

	$(CXX) -g $(LDFLAGS) $(LIBS) *.o -o $@


clean:
	rm *.o osca

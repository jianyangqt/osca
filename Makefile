
# -----------------------------------------------------------------
#   Makefile for osc 
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------

# Directory of the target
OUTPUT = osc

# Compiler
CXX = g++

# EIGEN library
EIGEN_PATH = /ibscratch/users/uqfzhan7/Lib/eigen

# Intel MKL library
#MKL_PATH = /opt/intel/mkl

# Compiler flags
CXXFLAGS = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG 
LIB += -static -lz -Wl,-lm -ldl
#LIB += -lz -Wl, -lm -ldl

HDR += l4_osc.h \
 	   l3_efile.h \
	   l2_efile.h \
            l2_reml.h \
            l2_bfile.h \
            l1_op_geno.h \
            l0_io.h \
            l0_mem.h \
            l0_stat.h \
            l0_com.h \
	   cdflib.h \
	   dcdflib.h \          
	   ipmpar.h 
          
            
SRC = l4_osc.cpp \
 	   l3_efile.cpp \
	   l2_efile.cpp \
            l2_reml.cpp \
            l2_bfile.cpp \
            l1_op_geno.cpp \
            l0_io.cpp \
            l0_mem.cpp \
            l0_stat.cpp \
            l0_com.cpp \
	   dcdflib.cpp

OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f *.o osc

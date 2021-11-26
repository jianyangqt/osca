# -----------------------------------------------------------------
#   Makefile for osca
#   Supported platforms: Unix / Linux
# -----------------------------------------------------------------

# You can reset this variable by change following lines or using 
# `make VARIABLE=VALUE `
EIGEN_PATH := /usr/include
MKL_INCLUDE := /usr/include
MKL_LIB := /usr/lib64/intel64
Rmath_INCLUDE := /usr/include
Rmath_LIB := /usr/lib64/lib

DEBUG :=

CXX = g++
ifdef DEBUG
CXXFLAGS = -Wall -g3 -Og -fopenmp
else
CXXFLAGS = -Wall -O2 -fopenmp
endif
CPPFLAGS =  -I$(EIGEN_PATH) -I$(MKL_INCLUDE) -I$(Rmath_INCLUDE)
LDFLAGS = -L$(MKL_LIB) -L$(Rmath_LIB)
LIBS = -lz -lgomp -lmkl_core -lpthread -lmkl_gnu_thread -lmkl_gf_lp64 -lRmath -lgsl -lgslcblas
objs = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

.PHONY: all
all: osca

osca: $(objs)
	$(CXX) $(CXXFLAGS) $(objs) $(LDFLAGS) $(LIBS) -o $@

osca_static: $(objs)
	$(CXX) $(CXXFLAGS) $(objs) $(LDFLAGS) -static -Wl,--start-group $(LIBS) -Wl,--end-group -o $@

$(objs): %.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -f src/*.o
	rm -f osca
	rm -f osca_static

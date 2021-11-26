# -----------------------------------------------------------------
#   Makefile for osca
#   Supported platforms: Unix / Linux
# -----------------------------------------------------------------

# You can reset this variable by change following lines or using 
# `make VARIABLE=VALUE `
EIGEN_PATH := /usr/include
MKL_INCLUDE := /usr/include
MKL_LIB := /usr/lib64/intel64
RMath_INCLUDE := /usr/include
RMath_LIB := /usr/lib64/lib

DEBUG :=

CXX = g++
ifdef DEBUG
CXXFLAGS = -Wall -g3 -Og -fopenmp
else
CXXFLAGS = -Wall -O3 -fopenmp
endif
CPPFLAGS =  -I$(EIGEN_PATH) -I$(MKL_INCLUDE) -I$(RMath_LIB)
LDFLAGS = -L$(MKL_LIB) -L$(RMath_LIB)
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

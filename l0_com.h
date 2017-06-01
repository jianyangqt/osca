//
//  l0_com.h
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l0_com__
#define __osc__l0_com__


#include <fstream>
#include <sstream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include <string>
#include <vector>
#include <iostream>
#include <zlib.h>
#include <bitset>
#include <map>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <iomanip>
#include <limits>

using namespace std;
using namespace Eigen;

#define WKSPACE_MIN_MB 64
#define WKSPACE_DEFAULT_MB 2048

//4 static
#define FNAMESIZE 4096
#define MAX_LINE_SIZE 0x60000
//4 dynamic
#define MAXLINEBUFLEN 0x4000000
#define MAXPROBENUM 0x100000 // >850k


typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;
typedef unsigned long		uintptr_t;

#define ONELU 1LLU
#define CACHELINE 64 
#define LOGPRINTF(...) sprintf(logbuf, __VA_ARGS__); logprintb();
#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - ONELU)))

const char err_fopen[] = "Error: Failed to open %s.\n";
const double FloatErr=numeric_limits<double>::epsilon();
extern int thread_num;
extern FILE* logfile;
extern char logbuf[];
extern char Tbuf[];
extern uint32_t g_debug_on;
extern uint32_t g_log_failed;
extern long long mem_left;

extern unsigned char* wkspace_base;
extern uintptr_t wkspace_left;


extern void logprintb();
extern void to_upper(char* str, int len);
extern void to_upper(string &str);
extern void to_lower(string &str);
extern bool has_suffix(const string &str, const string &suffix);
extern bool has_prefix(const string &str, const string &prefix);
extern bool has_prefix(const char* char1, const char* subchar);
extern string dtos(double value);

template <typename T>
extern inline string atos (T const& a)
{
    stringstream ss;
    ss << a;
    return(ss.str());
}
template <typename T>
extern inline T Abs(T const& a)
{
    return (T{} < a) ? a : -a;
}
extern int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode);
extern int32_t fwrite_checked(const void* buf, size_t len, FILE* outfile);
static inline int32_t fputs_checked(const char* ss, FILE* outfile) {
    fputs(ss, outfile);
    return ferror(outfile);
}
extern void split_str(char* tbuf, vector<string> &strs, uint32_t skip);
extern uint64_t split_string_skip(const string &str, vector<string> &vec_str, string separator, int num2skip);
extern int split_string(const string &str, vector<string> &vec_str, string separator=" ,\t;\n");
extern void getRank(vector<int> &a, vector<int> &b);
extern void getUnique(vector<int> &a);
extern void strarr2strvec(const char **a, int sizea, vector<string> &b);
extern void match(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC);
extern void match_only(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC);
extern void FileExist(char* filename);
extern bool FloatEqual(double lhs, double rhs);
extern int rand_seed();
extern double var(const vector<double> &x);
extern void get_bottom_indices(vector<int> &b, vector<double> &a, int num);
extern void strcpy2(char** to, string from);
extern void free2(char** to);
#endif /* defined(__osc__l0_com__) */

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
#include <random>
#ifndef __APPLE__
#include <omp.h>
#endif

using namespace std;
using namespace Eigen;

#define WKSPACE_MIN_MB 64
#define WKSPACE_DEFAULT_MB 2048

//4 static
#define FNAMESIZE 4096
#define MAX_LINE_SIZE 0x60000
//4 dynamic
#define MAXLINEBUFLEN 0x40000000
#define MAXPROBENUM 0x100000 // >850k

#define MISSING_PHENO -1e10
#define MISSING_PROFILE 1e10

//typedef unsigned int uint32_t;
//typedef unsigned long long uint64_t;
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
extern char* outfileName;
extern int remlstatus;
extern uintptr_t wkspace_left;
extern bool prt_mid_rlt;

extern bool loud;
extern bool remloasi;
extern void logprintb();
extern void TERMINATE();
extern void to_upper(char* str, int len);
extern void to_upper(string &str);
extern void to_lower(string &str);
extern bool has_suffix(const string &str, const string &suffix);
extern bool has_prefix(const string &str, const string &prefix);
extern bool has_prefix(const char* char1, const char* subchar);
extern string dtos(double value);
extern int MOMENT_APPROX;
template <typename T>
extern inline string atos (T const& a)
{
    stringstream ss;
    ss << a;
    return(ss.str());
}
template <typename T>
extern inline string atosm (T const& a)
{
    if(a==-9) return("NA");
    stringstream ss;
    ss << a;
    return(ss.str());
}
template <typename T>
extern inline T Abs(T const& a)
{
    return (T{} < a) ? a : -a;
}
template <typename T>
extern inline T sign (T const& a, T const& b)
{
    if(b>=0) return(abs(a));
    else return(-abs(a));
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
template <typename T>
extern void getRank(vector<T> &a, vector<T> &b)
{
    b.resize(a.size());
    //#pragma omp parallel for
    for (long i = a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        b[i] = count;
    }
}
template <typename T>
extern void getRank2(vector<T> &a, vector<T> &b)
{
    b.resize(a.size());
    //#pragma omp parallel for
    for (long i = a.size()-1; i >= 0; i--)
    {
        int count = 0 , ecount=0;
        double w=0;
        for (int j = 0; j < a.size(); j++) {
            if (a[j] < a[i]) count++;
            if (a[j] == a[i]) ecount++;
        }
        if(ecount>1)
        {
            for(int j=0;j<ecount;j++) w+=j;
            w/=ecount;
        }
        b[i] = count + w;
    }
}
template <typename T>
extern void getRank2R(vector<T> &a, vector<T> &b)
{
    b.resize(a.size());
    #pragma omp parallel for
    for (long i = a.size()-1; i >= 0; i--)
    {
        int count = 0 , ecount=0;
        double w=0;
        for (int j = 0; j < a.size(); j++) {
            if (a[j] < a[i]) count++;
            if (a[j] == a[i]) ecount++;
        }
        if(ecount>1)
        {
            for(int j=0;j<ecount;j++) w+=j;
            w/=ecount;
        }
        b[i] =1 + count + w;
    }
}

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

template <typename T>
extern void free2(T** to)
{
    if(*to)
    {
        delete(*to);
        *to=NULL;
    }
}
template<class T>
extern inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
extern inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

extern void removeRow(MatrixXd &matrix, unsigned int rowToRemove);
extern void removeColumn(MatrixXd &matrix, unsigned int colToRemove);
extern void removeRow(MatrixXf &matrix, unsigned int rowToRemove);
extern void removeColumn(MatrixXf &matrix, unsigned int colToRemove);
extern void inverse_V(MatrixXd &Vi, bool &determinant_zero);
extern void inverse_V2(MatrixXd &Vi, bool &determinant_zero);
extern void subMatrix_symm(MatrixXd &to,MatrixXd &from,vector<int> &idx);
extern double mean(const vector<double> &x);
extern double weight_mean(const vector<double> &x, const vector<double> &weight);
extern double sum(const vector<double> &x);
extern double cov(const vector<double> &x, const vector<double> &y);
extern double cor(vector<double> &y, vector<double> &x);
extern void standardise(vector<double> &data, bool divid_by_std);
#endif /* defined(__osc__l0_com__) */

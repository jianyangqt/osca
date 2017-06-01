//
//  l0_com.cpp
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l0_com.h"

int thread_num=1;
FILE* logfile = NULL;
uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
long long mem_left=0;
char logbuf[MAX_LINE_SIZE];
char Tbuf[MAX_LINE_SIZE];

struct IncGenerator {
    int current_;
    IncGenerator (int start) : current_(start) {}
    int operator() () { return current_++; }
};

void logstr(const char* ss) {
    if (!g_debug_on) {
        fputs(ss, logfile);
        if (ferror(logfile)) {
            printf("\nWarning: Logging failure on:\n%s\nFurther logging will not be attempted in this run.\n", ss);
            g_log_failed = 1;
        }
    } else {
        if (g_log_failed) {
            fputs(ss, stdout);
            fflush(stdout);
        } else {
            fputs(ss, logfile);
            if (ferror(logfile)) {
                printf("\nError: Debug logging failure.  Dumping to standard output:\n%s", ss);
                g_log_failed = 1;
            } else {
                fflush(logfile);
            }
        }
    }
}

void logprint(const char* ss) {
    logstr(ss);
    fputs(ss, stdout);
}

void logprintb() {
    logstr(logbuf);
    fputs(logbuf, stdout);
}

void to_upper(char* str, int len)
{
    int i=0;
    for(i=0; i<len; i++){
        if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
    }
}
void to_upper(string &str)
{
    int i=0;
    for(i=0; i<str.size(); i++){
        if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
    }
}

void to_lower(string &str)
{
    int i=0;
    for(i=0; i<str.size(); i++){
        if(str[i]>='A' && str[i]<='Z') str[i]-='A'-'a';
    }
}

bool has_suffix(const string &str, const string &suffix)
{
    return str.size() >= suffix.size() &&
    str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}
bool has_prefix(const string &str, const string &prefix)
{
    return str.size() >= prefix.size() &&
    str.compare(0, prefix.size(), prefix) == 0;
}
bool has_prefix(const char* char1, const char* subchar)
{
    bool match=1;
    while(*subchar != '\0')
    {
        if(*subchar++ != *char1++) match=0;
    }
    return match;
}
string dtos(double value)
{
    stringstream ss;
    ss << scientific<< value;
    // ss << fixed << setprecision(400) << __value;
    return(ss.str());
}

int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
    *target_ptr = fopen(fname, mode);
    if (!(*target_ptr)) {
        LOGPRINTF(err_fopen, fname);
        return -1;
    }
    return 0;
}
int32_t fwrite_checked(const void* buf, size_t len, FILE* outfile) {
    while (len > 0x7ffe0000) {
        // OS X can't perform >2GB writes
        fwrite(buf, 1, 0x7ffe0000, outfile);
        buf = &(((unsigned char*)buf)[0x7ffe0000]);
        len -= 0x7ffe0000;
    }
    fwrite(buf, 1, len, outfile);
    return ferror(outfile);
}

static inline char* skip_initial_spaces(char* sptr) {
    while ((*sptr == ' ') || (*sptr == '\t') || (*sptr == '\n')) {
        sptr++;
    }
    return sptr;
}

char* next_token(char* sptr) {
    if (!sptr) {
        return NULL;
    }
    while ((*sptr != ' ') && (*sptr != '\t') && (*sptr != '\n')) {
        if (!(*sptr)) {
            return sptr--; //in case of no delimiter befor '\0'
        }
        sptr++;
    }
    return sptr;
}

void split_str(char* tbuf, vector<string> &strs, uint32_t skip)
{
    strs.clear();
    char* start_ptr;
    char* end_ptr;
    tbuf=skip_initial_spaces(tbuf);
    end_ptr=start_ptr=tbuf;
    uint32_t strnum=0;
    while(*end_ptr !='\0')
    {
        end_ptr=next_token(end_ptr);
        string str(start_ptr,end_ptr-start_ptr);
        strnum++;
        end_ptr= skip_initial_spaces(end_ptr);
        start_ptr=end_ptr;
        if(strnum>skip) strs.push_back(str);
    }
}
uint64_t split_string_skip(const string &str, vector<string> &vec_str, string separator, int num2skip) //form head
{
    if(str.empty()) return 0;
    vec_str.clear();
    
    bool look=false;
    string str_buf;
    int count=0;
    
    for(int i=0; i<str.size(); i++)
    {
        if( separator.find(str[i])==string::npos )
        {
            if(!look) look=true;
            str_buf += str[i];
        }
        else
        {
            if(look)
            {
                look=false;
                if(++count>num2skip) vec_str.push_back(str_buf);
                str_buf.erase(str_buf.begin(), str_buf.end());
            }
        }
    }
    if(look) vec_str.push_back(str_buf);
    
    return vec_str.size();
}

int split_string(const string &str, vector<string> &vec_str, string separator)
{
    if(str.empty()) return 0;
    vec_str.clear();
    
    int i=0;
    bool look=false;
    string str_buf;
    string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
    string::size_type pos;
    
    for(i=0; i<separator.size(); i++){
        pos=symbol_pool.find(separator[i]);
        if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
    }
    
    for(i=0; i<str.size(); i++){
        if( symbol_pool.find(str[i])!=string::npos ){
            if(!look) look=true;
            str_buf += str[i];
        }
        else{
            if(look){
                look=false;
                vec_str.push_back(str_buf);
                str_buf.erase(str_buf.begin(), str_buf.end());
            }
        }
    }
    if(look) vec_str.push_back(str_buf);
    
    return (int)vec_str.size();
}



void getRank(vector<int> &a, vector<int> &b)
{
    b.resize(a.size());
    for (long i = a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        b[i] = count;
    }
}
void get_bottom_indices(vector<int> &b, vector<double> &a, int num)
{
    vector<double> c;
     long asize=a.size();
    if(asize==0)
    {
        LOGPRINTF("Error: no element found in the vector.\n");
        exit(EXIT_FAILURE);
    }
   
    if(num>=asize)
    {
        b.resize(asize);
        for(int i=0;i<asize;i++) b[i]=i;
    } else {
        b.resize(num);
        c.resize(num);
        int ptr=0; //index of c
        for (int i =0; i < asize; i++)
        {
            if(i<num) {
                b[i]=i;
                c[i]=a[i];
                if(c[ptr]<a[i]) ptr=i;
            } else {
                if(c[ptr]>a[i])
                {
                    c[ptr]=a[i];
                    b[ptr]=i;
                    for(int j=0;j<num;j++)
                        if(c[j]>c[ptr]) ptr=j;
                }
            }
        }
    }
}
void getUnique(vector<int> &a)
{
    sort(a.begin(),a.end());
    vector<int> ::iterator it=unique(a.begin(),a.end());
    a.erase(it,a.end());
}

/*********************************
 * const char *flgs[] = { "a","bb","ccc"};
 * vector<string> b;
 * strarr2strvec(flgs,sizeof(flgs),b);
 *********************************/
void strarr2strvec(const char **a, int sizea, vector<string> &b)
{
    b.insert(b.begin(),a, a+sizea/sizeof(a[0]));
}
void match(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC)
{
    int i=0;
    map<string, int> id_map;
    map<string, int>::iterator iter;
    VecC.clear();
    for(i=0; i<VecB.size(); i++) id_map.insert(pair<string,int>(VecB[i], i));
    for(i=0; i<VecA.size(); i++){
        iter=id_map.find(VecA[i]);
        if(iter==id_map.end()) VecC.push_back(-9);
        else VecC.push_back(iter->second);
    }
    
}


void match_only(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC)
{
    int i = 0;
    map<string, int> id_map;
    map<string, int>::iterator iter;
    VecC.clear();
    for (i = 0; i<VecB.size(); i++) id_map.insert(pair<string, int>(VecB[i], i));
    for (i = 0; i<VecA.size(); i++){
        iter = id_map.find(VecA[i]);
        if (iter != id_map.end()) VecC.push_back(iter->second);
    }
    
}
void FileExist(char* filename)
{
    ifstream ifile(filename);
    if(!ifile){
        LOGPRINTF("Error: can not open the file %s to read.",filename);
        exit(EXIT_FAILURE);
    }
}
bool FloatEqual(double lhs, double rhs)
{
    if (Abs(lhs - rhs) < FloatErr) return true;
    return false;
}
int rand_seed()
{
    stringstream str_strm;
    str_strm<<time(NULL);
    string seed_str=str_strm.str();
    reverse(seed_str.begin(), seed_str.end());
    seed_str.erase(seed_str.begin()+7, seed_str.end());
    return(abs(atoi(seed_str.c_str())));
}

double var(const vector<double> &x)
{
    long size = x.size();
    if(size<=1) return(0.0);
    int i=0;
    double mu=0.0, s2=0.0;
    for(i=0; i<size; i++) mu+=x[i];
    mu/=(double)size;
    for(i=0; i<size; i++) s2+=(x[i]-mu)*(x[i]-mu);
    s2/=(double)(size-1);
    return (double)s2;
}
/*
vector<uint32_t> sort_re_index(const vector<double> &x) {
    vector<uint32_t> x_index(x.size());
    
    //iota(x_index.begin(), x_index.end(), 0); //c++11
    IncGenerator g(0);
    generate(x_index.begin(), x_index.end(),g);
    
    sort(x_index.begin(), x_index.end(), [&x](uint32_t index1, uint32_t index2) {return x[index1] > x[index2];});
    
    return x_index;
}
 */
void strcpy2(char** to, string from)
{
    char* tmp=new char[from.size() + 1];
    copy(from.begin(), from.end(), tmp);
    tmp[from.size()]='\0';
    *to=tmp;
}
void free2(char** to)
{
    if(*to)
    {
        delete(*to);
        *to=NULL;
    }
}

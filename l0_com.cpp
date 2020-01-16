//
//  l0_com.cpp
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l0_com.h"

int thread_num=1;
double PFISHER = 0.05;
FILE* logfile = NULL;
char* outfileName=NULL;
int remlstatus = 0; // -1: matrix is not invertible. -2:more than half of the variance components are constrained. -3: variance is <0 or >1 . -4: not coverge. -5 constrained
uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
long long mem_left=0;
char logbuf[MAX_LINE_SIZE];
char Tbuf[MAX_LINE_SIZE];
bool prt_mid_rlt=false;
bool loud = false;
int MOMENT_APPROX = 0x0;
bool remloasi = true;
struct IncGenerator {
    int current_;
    IncGenerator (int start) : current_(start) {}
    int operator() () { return current_++; }
};

void logstr(const char* ss) {
    if (!g_debug_on) {
        fputs(ss, logfile);
        fflush(stdout);
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
void TERMINATE() {
    if(logfile) {
        fflush(logfile);
        fclose(logfile);
    }
    exit(EXIT_FAILURE);
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
    while ((*sptr == ' ') || (*sptr == '\t') || (*sptr == '\n') || (*sptr == '\r') ) {
        sptr++;
    }
    return sptr;
}

char* next_token(char* sptr) {
    if (!sptr) {
        return NULL;
    }
    while ((*sptr != ' ') && (*sptr != '\t') && (*sptr != '\n') && (*sptr != '\r')) {
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

void get_bottom_indices(vector<int> &b, vector<double> &a, int num)
{
    vector<double> c;
     long asize=a.size();
    if(asize==0)
    {
        LOGPRINTF("Error: no element found in the vector.\n");
        TERMINATE();
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
        TERMINATE();
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
void subMatrix_symm(MatrixXd &to,MatrixXd &from,vector<int> &idx)
{
    long num = idx.size();
    if(num>from.cols())
    {
        LOGPRINTF("ERROR: to extract sub-matrix. the size of index %ld is larger than the matrix dimension (%ld,%ld).\n",num, from.rows(),from.cols());
        exit(EXIT_FAILURE);
    } else if(num==0) {
        LOGPRINTF("ERROR: the size of index is 0.\n");
        exit(EXIT_FAILURE);
    } else {
        to.resize(num,num);
        for(int i=0;i<num;i++)
            for(int j=i;j<num;j++)
                to(i,j)=to(j,i)=from(idx[i],idx[j]);
    }
}
void removeRow(MatrixXd &matrix, unsigned int rowToRemove)
{
    unsigned long numRows = matrix.rows()-1;
    unsigned long numCols = matrix.cols();
    
    if( rowToRemove < numRows ) {
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
        matrix.conservativeResize(numRows,numCols);
    } else if(rowToRemove == numRows) {
        matrix.conservativeResize(numRows,numCols);
    } else {
        LOGPRINTF("ERROR: out of bound. row id %u is no less than total row numbers %ld \n",rowToRemove, matrix.rows());
        exit(EXIT_FAILURE);
    }
}

void removeColumn(MatrixXd &matrix, unsigned int colToRemove)
{
    unsigned long numRows = matrix.rows();
    unsigned long numCols = matrix.cols()-1;
    
    if( colToRemove < numCols ) {
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
        
        matrix.conservativeResize(numRows,numCols);
    } else if(colToRemove == numCols) {
        matrix.conservativeResize(numRows,numCols);
    } else {
        LOGPRINTF("ERROR: out of bound. column id %u is no less than total column numbers %ld \n",colToRemove, matrix.cols());
        exit(EXIT_FAILURE);
    }
    
}
void removeRow(MatrixXf &matrix, unsigned int rowToRemove)
{
    unsigned long numRows = matrix.rows()-1;
    unsigned long numCols = matrix.cols();
    
    if( rowToRemove < numRows ) {
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
        
        matrix.conservativeResize(numRows,numCols);
    } else if(rowToRemove == numRows) {
        matrix.conservativeResize(numRows,numCols);
    } else {
        LOGPRINTF("ERROR: out of bound. row id %u is no less than total row numbers %ld \n",rowToRemove, matrix.rows());
        exit(EXIT_FAILURE);
    }
    
}

void removeColumn(MatrixXf &matrix, unsigned int colToRemove)
{
    unsigned long numRows = matrix.rows();
    unsigned long numCols = matrix.cols()-1;
    
    if( colToRemove < numCols ) {
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
        
        matrix.conservativeResize(numRows,numCols);
    } else if(colToRemove == numCols) {
        matrix.conservativeResize(numRows,numCols);
    } else {
        LOGPRINTF("ERROR: out of bound. column id %u is no less than total column numbers %ld \n",colToRemove, matrix.cols());
        exit(EXIT_FAILURE);
    }
    
}

void inverse_V(MatrixXd &Vi, bool &determinant_zero)
{
    SelfAdjointEigenSolver<MatrixXd> eigensolver(Vi);
    VectorXd eval = eigensolver.eigenvalues();
    for(int i=0;i<eval.size();i++)
    {
        //if(abs(eval(i))<1e-6) {
         if(abs(eval(i))<1e-8) { // I observed a scenario with the min eigenval as 1.393741e-09, the matrix is still invertible.
             /***for test***/
             
             string filena=string(outfileName)+".H.mat";
             FILE* tmpfi=fopen(filena.c_str(),"w");
             if(!tmpfi)
             {
                 LOGPRINTF("error open file.\n");
                 TERMINATE();
             }
             for( int j=0;j<Vi.rows();j++)
             {
                 string str ="";
                 for(int t=0;t<Vi.cols();t++) str +=atos(Vi(j,t)) + '\t';
                 str += '\n';
                 fputs(str.c_str(),tmpfi);
             }
             fclose(tmpfi);
              
             /***end of test***/
             
            determinant_zero=true;
            eval(i)=0;
        } else {
            eval(i) = 1.0 / eval(i);
        }
    }
    Vi = eigensolver.eigenvectors() * DiagonalMatrix<double, Dynamic, Dynamic>(eval) * eigensolver.eigenvectors().transpose();
}
void inverse_V2(MatrixXd &Vi, bool &determinant_zero)
{
    JacobiSVD<MatrixXd> svd(Vi, ComputeThinU | ComputeThinV);
    VectorXd svdd=svd.singularValues();
    MatrixXd svdU=svd.matrixU();
    MatrixXd svdV=svd.matrixV();
    
    for(int i=0;i<svdd.size();i++)
    {
        //if(abs(svdd(i))<1e-6) {
            if(abs(svdd(i))<1e-10) { // I observed a scenario with the min eigenval as 1.393741e-09, the matrix is still invertible.
            determinant_zero=true;
            svdd(i)=0;
        } else {
            svdd(i) = 1.0 / svdd(i);
        }
    }
    
    Vi = svdU * DiagonalMatrix<double, Dynamic, Dynamic>(svdd) * svdV.transpose();
}

double mean(const vector<double> &x)
{
    long size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    d_buf/=(double)size;
    return (double)d_buf;
}
double weight_mean(const vector<double> &x, const vector<double> &weight)
{
    if(x.size()!=weight.size())
    {
        LOGPRINTF("ERROR: different size of x and weight.\n")
        TERMINATE();
    }
    long size = x.size();
    int i=0;
    double d_buf=0.0,e_buf=0.0;
    for(i=0; i<size; i++) {
        d_buf+=x[i]*weight[i];
        e_buf+=weight[i];
    }
    d_buf/=e_buf;
    return (double)d_buf;
}
double sum(const vector<double> &x)
{
    long size = x.size();
    double d_buf=0.0;
    for(int i=0; i<size; i++) d_buf+=x[i];
    return d_buf;
}
double cov(const vector<double> &x, const vector<double> &y)
{
    long size = x.size();
    int i=0;
    double mu1=0.0, mu2=0.0, c=0.0;
    for(i=0; i<size; i++){
        mu1+=x[i];
        mu2+=y[i];
    }
    mu1/=(double)size;
    mu2/=(double)size;
    
    for(i=0; i<size; i++) c+=(x[i]-mu1)*(y[i]-mu2);
    c/=(double)(size-1);
    return c;
}

double cor(vector<double> &y, vector<double> &x)
{
    long N = x.size();
    if (N != y.size() || N < 1) {
        LOGPRINTF("Error: The lengths of x and y do not match.\n");
        TERMINATE();
    }
    
    int i = 0;
    double d_buf = 0.0, y_mu = 0.0, x_mu = 0.0, x_var = 0.0, y_var = 0.0, cov = 0.0;
    for (i = 0; i < N; i++) {
        x_mu += x[i];
        y_mu += y[i];
    }
    x_mu /= (double) N;
    y_mu /= (double) N;
    for (i = 0; i < N; i++) {
        d_buf = (x[i] - x_mu);
        x_var += d_buf*d_buf;
        d_buf = (y[i] - y_mu);
        y_var += d_buf*d_buf;
    }
    x_var /= (double) (N - 1.0);
    y_var /= (double) (N - 1.0);
    for (i = 0; i < N; i++) cov += (x[i] - x_mu)*(y[i] - y_mu);
    cov /= (double) (N - 1);
    double a = 0.0, b = 0.0, sse = 0.0, a_se = 0.0, b_se = 0.0, r = 0.0;
    if (x_var > 0.0) b = cov / x_var;
    a = y_mu - b*x_mu;
    for (i = 0; i < N; i++) {
        d_buf = y[i] - a - b * x[i];
        sse += d_buf*d_buf;
    }
    if (x_var > 0.0) {
        a_se = sqrt((sse / (N - 2.0))*(1.0 / N + x_mu * x_mu / (x_var * (N - 1.0))));
        b_se = sqrt(sse / x_var / (N - 1.0) / (N - 2.0));
    }
    if (x_var > 0.0 && y_var > 0.0) {
        r = cov / sqrt(y_var * x_var);
    }
    
    return (r);
}
void standardise(vector<double> &data, bool divid_by_std)
{
    // missing is labelled as 1e10
    double mu=0.0, nonmiss=0;
    long n = data.size();
        for(int i=0; i<n; i++){
            if(data[i]<1e9)
            {
                mu += data[i];
                nonmiss += 1.0;
            }
        }
        if(nonmiss>0) mu /=nonmiss;
    
    for(int i=0; i<n; i++){
            if(data[i]<1e9) data[i] -= mu;
            else data[i] = 0.0;
        }
    
    if(divid_by_std)
    {
        double sd = sqrt(var(data));
        #pragma omp parallel for
        for(int i=0; i<n; i++){
           
                if(fabs(sd)>1e-30) data[i] /= sd;
                else data[i] = 0.0;
        }
    }
    
}

void coeff_mat(const vector<string> &vec, MatrixXd &coeff_mat, string errmsg1, string errmsg2) {
    vector<string> value(vec);
    stable_sort(value.begin(), value.end());
    value.erase(unique(value.begin(), value.end()), value.end());
    if (value.size() > 0.5 * vec.size()) {
        printf("%s\n",errmsg1.c_str());
        exit(EXIT_FAILURE);
    } // throw("Error: too many classes for the envronmental factor. \nPlease make sure you input a discrete variable as the environmental factor.");
    if (value.size() == 1) {
        printf("%s\n",errmsg2.c_str());
        exit(EXIT_FAILURE);
    } //throw("Error: the envronmental factor should has more than one classes.");
    
    int i = 0, j = 0, row_num = vec.size(), column_num = value.size();
    map<string, int> val_map;
    for (i = 0; i < value.size(); i++) val_map.insert(pair<string, int>(value[i], i));
    
    coeff_mat.resize(row_num, column_num);
    coeff_mat.setZero(row_num, column_num);
    map<string, int>::iterator iter;
    for (i = 0; i < row_num; i++) {
        iter = val_map.find(vec[i]);
        coeff_mat(i, iter->second) = 1.0;
    }
}

double ssy(vector<double> &y)
{
    double rs = 0.0;
    for(int i=0;i<y.size();i++)
    {
        double tmp=y[i];
        rs += tmp*tmp;
    }
    return rs;
}

double msy(vector<double> &y)
{
    double rs = 0.0;
    for(int i=0;i<y.size();i++) rs += y[i];
    return rs*rs/y.size();
}


void progress(int &cur, double &disp, int ttl)
{
    double desti=1.0*cur/(ttl-1);
    if(desti>=disp)
    {
        printf("%3.0f%%\r", 100.0*desti);
        fflush(stdout);
        if(disp==0) disp+=0.05;
        else if(disp==0.05) disp+=0.2;
        else if(disp==0.25) disp+=0.5;
        else disp+=0.25;
    }
}

double  cor(VectorXd &Y, VectorXd &X, bool centered, bool standardised)
{
    if(Y.size()!= X.size())
    {
        printf("The lenght of vectors not match.\n");
        exit(EXIT_FAILURE);
    }
    double ld=-9;
    long n=Y.size();
    if(standardised)
    {
        double xy=X.dot(Y);
        ld=xy/(n-1);
    }
    else if(centered)
    {
        double xx=X.dot(X), yy=Y.dot(Y), xy=X.dot(Y);
        ld=xy/(sqrt(xx*yy));
    }
    else
    {
        double ysum=Y.sum();
        double xsum=X.sum();
        double xx=X.dot(X), yy=Y.dot(Y), xy=X.dot(Y);
        ld=(n*xy-xsum*ysum)/(sqrt(n*xx-xsum*xsum)*sqrt(n*yy-ysum*ysum));
    }
    return ld;
}

void  cor(VectorXd &CORR , VectorXd &x, MatrixXd &X, bool centered, bool standardised)
{
    if(x.size()!= X.rows())
    {
        printf("The lenght not match.\n");
        exit(EXIT_FAILURE);
    }
    long n=x.size();
    if(standardised)
    {
        CORR = x.transpose()*X/(n-1);
    }
    else if(centered)
    {
        printf("NOT IMPLEMENTED.\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("NOT IMPLEMENTED.\n");
        exit(EXIT_FAILURE);
    }

}

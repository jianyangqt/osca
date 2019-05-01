//
//  l1_stat.cpp
//  osc
//
//  Created by Futao Zhang on 23/12/2016.
//  Copyright © 2016 Futao Zhang. All rights reserved.
//

#include "l1_stat.hpp"

void reg(VectorXd &y, VectorXd &x,vector<double> &rst)
{
    // THis function only for fast linear, missing is not allowed.
    if(y.size()!=x.size())
    {
        LOGPRINTF("ERROR: y has different size from x.\n" );
        TERMINATE();
    }
    if(y.size()<1)
    {
        LOGPRINTF("ERROR: the lenght of y <1.\n" );
        TERMINATE();
    }
    long nmiss = y.size();
    double mu_y=0.0, mu_x=0.0;
    rst.clear();
    for(int j=0; j<nmiss; j++)
    {
        double phval=y(j);
        double val=x(j);
        mu_y+=phval;
        mu_x+=val;
    }
    mu_y/=nmiss;
    mu_x/=nmiss;
    double xty=0.0, xtx=0.0, yty=0.0;
    for(int j=0; j<nmiss; j++)
    {
        double phval=y(j);
        double val=x(j);
        yty+=(phval-mu_y)*(phval-mu_y);
        xty+=(val-mu_x)*(phval-mu_y);
        xtx+=(val-mu_x)*(val-mu_x);
    }
    xty/=nmiss-1;
    xtx/=nmiss-1;
    yty/=nmiss-1;
    
    // Test statistics
    double beta = 0, se=-9, t=0, t_p=1;
    if(xtx>0)
    {
         beta = xty / xtx;
         se = sqrt(( yty/xtx - (xty*xty)/(xtx*xtx) ) / (nmiss-2));
         t = beta / se;
         t_p=pT(t,nmiss-2);
    }
    rst.push_back(beta);
    rst.push_back(se);
    rst.push_back(t_p);
}
void adjusted_reg(VectorXd &yresi, VectorXd &xresi,vector<double> &rst, int exdf)
{
    // THis function only for fast linear, missing is not allowed.
    if(yresi.size()!=xresi.size())
    {
        LOGPRINTF("ERROR: y has different size from x.\n" );
        TERMINATE();
    }
    if(yresi.size()<1)
    {
        LOGPRINTF("ERROR: the lenght of y <1.\n" );
        TERMINATE();
    }
    long n = yresi.size();
    rst.clear();
    double xty=0.0, xtx=0.0, yty=0.0;
    for(int j=0; j<n; j++)
    {
        double phval=yresi(j);
        double val=xresi(j);
        yty+=phval*phval;
        xty+=val*phval;
        xtx+=val*val;
    }
    xty/=n-1;
    xtx/=n-1;
    yty/=n-1;
    
    // Test statistics
    double beta = 0, se=-9, t=0, t_p=1;
    if(xtx>0)
    {
        beta = xty / xtx;
        se = sqrt(( yty/xtx - (xty*xty)/(xtx*xtx) ) / (n-exdf-2));
        t = beta / se;
        t_p=pT(t,n-exdf-2);
    }
    rst.push_back(beta);
    rst.push_back(se);
    rst.push_back(t_p);
}

MatrixXd reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table)
{
    long N=x.size();
    if(N!=y.size() || N<1) throw("Error: The lengths of x and y do not match.");
    
    int i=0;
    double d_buf=0.0, y_mu=0.0, x_mu=0.0, x_var=0.0, y_var=0.0, cov=0.0;
    for(i=0; i<N; i++){ x_mu+=x[i]; y_mu+=y[i]; }
    x_mu/=(double)N;
    y_mu/=(double)N;
    for(i=0; i<N; i++){
        d_buf=(x[i]-x_mu);
        x_var+=d_buf*d_buf;
        d_buf=(y[i]-y_mu);
        y_var+=d_buf*d_buf;
    }
    x_var/=(double)(N-1.0);
    y_var/=(double)(N-1.0);
    for(i=0; i<N; i++) cov+=(x[i]-x_mu)*(y[i]-y_mu);
    cov/=(double)(N-1);
    double a=0.0, b=0.0, sse=0.0, a_se=0.0, b_se=0.0, p=0.0, rsq=0.0, r=0.0;
    if(x_var>0.0) b=cov/x_var;
    a=y_mu-b*x_mu;
    for(i=0; i<N; i++){
        d_buf=y[i]-a-b*x[i];
        sse+=d_buf*d_buf;
    }
    if(x_var>0.0){
        a_se=sqrt((sse/(N-2.0))*(1.0/N + x_mu*x_mu/(x_var*(N-1.0))));
        b_se=sqrt(sse/x_var/(N-1.0)/(N-2.0));
    }
    if(x_var>0.0 && y_var>0.0){
        r=cov/sqrt(y_var*x_var);
        rsq=r*r;
    }
    double t=0.0;
    if(b_se>0.0) t=fabs(b/b_se);
    p=t_prob(N-2.0, t, true);
    rst.clear();
    rst.push_back(b); rst.push_back(b_se); rst.push_back(p); rst.push_back(rsq); rst.push_back(r);
    
    MatrixXd reg_sum(3,3);
    if(table){
        reg_sum(2,0)=rsq;
        reg_sum(1,0)=b;
        reg_sum(1,1)=b_se;
        reg_sum(1,2)=p;
        if(a_se>0.0) t=fabs(a/a_se);
        p=t_prob(N-2.0, t, true);
        reg_sum(0,0)=a;
        reg_sum(0,1)=a_se;
        reg_sum(0,2)=p;
        return(reg_sum);
    }
    return(reg_sum);
}

bool lin(VectorXd &y, MatrixXd &C, VectorXd &x, vector<double> &rst)
{
    long N=y.size();
    if(N!=x.size() || N<1) {
        LOGPRINTF("Error: The lengths of x and y do not match.\n");
        TERMINATE();
    }
    if(N!=C.rows() || N<1) {
        LOGPRINTF("Error: The row number of C and the length of y do not match.\n");
        TERMINATE();
    }
    
    MatrixXd X(C.rows(), C.cols()+1);
   
    X.block(0,0,C.rows(),C.cols())=C;
    long x_idx=X.cols()-1;
    X.col(x_idx)=x;
    MatrixXd XtX_i;
    XtX_i=X.transpose()*X;
    bool determinant_zero=false;
    inverse_V(XtX_i, determinant_zero);
    VectorXd b_hat=XtX_i*X.transpose()*y;
    VectorXd residual=(y-X*b_hat);
    residual=residual.array()*residual.array();
    double sy=sqrt(residual.sum()/(y.size()-X.cols()));
    VectorXd se=sy*XtX_i.diagonal().array().sqrt();
    VectorXd t=b_hat.array()/se.array();
    double z=t(x_idx)*t(x_idx);
    double p=pchisq(z,1);
    //double p=t_prob(Y.size()-2.0, fabs(t(x_idx)), true);
    rst.clear();
    rst.push_back(b_hat(x_idx)); rst.push_back(se(x_idx)); rst.push_back(p);
    return determinant_zero;
}
void lin2(VectorXd &y, MatrixXd &X, vector<double> &rst) // the last column of X is x, other columns are covariates
{
    long N=y.size();
    
    if(N!=X.rows() || N<1) {
        LOGPRINTF("Error: The row number of C and the length of y do not match.\n");
        TERMINATE();
    }
    
    long x_idx=X.cols()-1;
    MatrixXd XtX_i;
    XtX_i=X.transpose()*X;
    bool determinant_zero=false;
    inverse_V(XtX_i, determinant_zero);
    VectorXd b_hat=XtX_i*X.transpose()*y;
    VectorXd residual=(y-X*b_hat);
    residual=residual.array()*residual.array();
    double sy=sqrt(residual.sum()/(y.size()-X.cols()));
    VectorXd se=sy*XtX_i.diagonal().array().sqrt();
    VectorXd t=b_hat.array()/se.array();
    double z=t(x_idx)*t(x_idx);
    double p=pchisq(z,1);
    //double p=t_prob(Y.size()-2.0, fabs(t(x_idx)), true);
    rst.clear();
    rst.push_back(b_hat(x_idx)); rst.push_back(se(x_idx)); rst.push_back(p);
}

void cov(MatrixXd &mat, MatrixXd &covar)
{
    MatrixXd centered = mat.rowwise() - mat.colwise().mean();
    covar = (centered.adjoint() * centered) / double(mat.rows() - 1);
}
double var(VectorXd &x)
{
    double xmean=x.mean();
    x=x.array()-xmean;
    double variance=x.dot(x)/(x.size()-1);
    return variance;
}
double var( vector<double> &x)
{
    long N=x.size();
    double d_buf=0.0, x_mu=0.0, x_var=0.0;
    for(int i=0; i<N; i++){ x_mu+=x[i]; }
    x_mu/=(double)N;
    for(int i=0; i<N; i++){
        d_buf=(x[i]-x_mu);
        x_var+=d_buf*d_buf;
    }
    x_var/=(double)(N-1.0);
    return x_var ;
}
double mean( vector<double> &x)
{
    long N=x.size();
    double x_mu=0.0;
    for(int i=0; i<N; i++){ x_mu+=x[i]; }
    x_mu/=(double)N;
    return x_mu ;
}

int bartlett(vector<double> &y,vector<double> &x, vector<double> &rst, double freq) //norm: rutrun 1; k=1 return -1; n[i]=1 return -2; N=k return -3
{
    int N=(int)y.size();
    vector<int> cat(N);
    #pragma omp parallel for
    for(int i=0;i<x.size();i++) cat[i]=(int)x[i];
    getUnique(cat);
    int k=(int)cat.size();
    if(k==1) return -1;
    vector<int> n(k);
    vector<double> sigmai(k),catd(k);
    for(int i=0;i<k;i++) {
        n[i]=0;
        int category=cat[i];
        catd[i]=(double)cat[i];
        vector<double> tmp;
        for(int j=0;j<N;j++) if((int)x[j]==category) {n[i]++; tmp.push_back(y[j]);}
        sigmai[i]=var(tmp);
    }
    double sigmap=0.0, a=0.0, b=0.0;
    for( int i=0;i<k;i++) {
        if(n[i]==1) return -2;
        sigmap+=(n[i]-1.0)*sigmai[i];
        a+=(n[i]-1.0)*log(sigmai[i]);
        b+=1/(n[i]-1.0);
    }
    if(N==k) return -3;
    sigmap*=(1.0/(N-k));
    a=(N-k)*log(sigmap)-a;
    b=1+ (1.0/(3*(k-1)))*(b-1.0/(N-k));
    double z2=a/b;
    double p=pchisq(z2, k-1);
    
    double zest=sqrt(qchisq(p,1));
    double domin=sqrt(2*freq*(1-freq)*(N+zest*zest)); //can't be zero here
    double best=zest/domin;
    double seest=1/domin;
    reg(sigmai,catd,rst);
    if(rst[0]<0) best *= -1.0; //updated by futao.zhang 20180321
    rst.clear();
    rst.push_back(best);
    rst.push_back(seest);
    rst.push_back(p);
    rst.push_back(z2);
    rst.push_back(k-1);
    return 1;
}

int leveneTest_mean(vector<double> &y,vector<double> &x, vector<double> &rst,double freq) //norm: rutrun 1; k=1 return -1;
{
    int N=(int)y.size();
    vector<int> cat(N);
    #pragma omp parallel for
    for(int i=0;i<x.size();i++) cat[i]=(int)x[i];
    getUnique(cat);
    int k=(int)cat.size();
    if(k==1) return -1;
    vector<int> n(k);
    vector<double> sigmai(k),catd(k);
    vector<double> ymean(k);
    vector<int> idx;
    vector<double> tmp;
    double nzi_z=0.0, zij_zi=0.0;
    for(int i=0;i<k;i++) {
        idx.clear();
        tmp.clear();
        n[i]=0;
        int category=cat[i];
        catd[i]=(double)cat[i];
        for(int j=0;j<N;j++) if((int)x[j]==category) { n[i]++; idx.push_back(j); tmp.push_back(y[j]);}
        sigmai[i]=var(tmp);
        ymean[i]=mean(tmp);
        double zit=0.0;
        for(int j=0;j<idx.size();j++) {
            double absmu=abs(y[idx[j]]-ymean[i]);
            y[idx[j]]=absmu; //now y is z_ij
            zit+=absmu;
        }
        zit/=idx.size();
        ymean[i]=zit; // now ymean is z_i
        for(int j=0;j<idx.size();j++) {
            zij_zi+=(y[idx[j]]-zit)*(y[idx[j]]-zit);
        }
    }
    double z=mean(y);
    
    for(int i=0;i<k;i++)
        nzi_z+=n[i]*(ymean[i]-z)*(ymean[i]-z);
    double stat=(N-k)*nzi_z/(k-1)/zij_zi;
    int df1=k-1;
    int df2=N-k;
    double p=F_prob(df1,df2,stat);
    
    double zest=sqrt(qchisq(p,1));
    double domin=sqrt(2*freq*(1-freq)*(N+zest*zest));
    double best=zest/domin;
    double seest=1/domin;
    reg(sigmai,catd,rst);
    if(rst[0]<0) best *= -1.0; //updated by futao.zhang 20180321
    rst.clear();
    rst.push_back(best);
    rst.push_back(seest);
    rst.push_back(p);
    rst.push_back(stat);
    rst.push_back(df1);
    rst.push_back(df2);
    return 1;
}
int leveneTest_median(vector<double> &y,vector<double> &x, vector<double> &rst,double freq) //norm: rutrun 1; k=1 return -1;
{
    // a little different with the R when the individual number is even.
    //in R, the mean of the middle two is taken as the median
    //in here, the value of the bigger of the middle two is tanken as the median
    int N=(int)y.size();
    vector<int> cat(N);
    for(int i=0;i<x.size();i++) cat[i]=(int)x[i];
    getUnique(cat);
    int k=(int)cat.size();
    if(k==1) return -1;
    vector<int> n(k);
    vector<double> sigmai(k),catd(k);
    vector<double> ymedian(k);
    vector<int> idx;
    vector<double> tmp;
    double nzi_z=0.0, zij_zi=0.0;
    for(int i=0;i<k;i++) {
        idx.clear();
        tmp.clear();
        n[i]=0;
        int category=cat[i];
        catd[i]=(double)cat[i];
        for(int j=0;j<N;j++) if((int)x[j]==category) { n[i]++; idx.push_back(j); tmp.push_back(y[j]);}
        sigmai[i]=var(tmp);
        ymedian[i]=median(tmp);
        double zit=0.0;
        for(int j=0;j<idx.size();j++) {
            double absmu=abs(y[idx[j]]-ymedian[i]);
            y[idx[j]]=absmu; //now y is z_ij
            zit+=absmu;
        }
        zit/=idx.size();
        ymedian[i]=zit; // now ymedian is z_i
        for(int j=0;j<idx.size();j++) {
            zij_zi+=(y[idx[j]]-zit)*(y[idx[j]]-zit);
        }
    }
    double z=mean(y);
    
    for(int i=0;i<k;i++)
        nzi_z+=n[i]*(ymedian[i]-z)*(ymedian[i]-z);
    double stat=(N-k)*nzi_z/(k-1)/zij_zi;
    int df1=k-1;
    int df2=N-k;
    double p=F_prob(df1,df2,stat);
    
    double zest=sqrt(qchisq(p,1));
    double domin=sqrt(2*freq*(1-freq)*(N+zest*zest));
    double best=zest/domin;
    double seest=1/domin;
    reg(sigmai,catd,rst);
    if(rst[0]<0) best *= -1.0; //updated by futao.zhang 20180321
    rst.clear();
    rst.push_back(best);
    rst.push_back(seest);
    rst.push_back(p);
    rst.push_back(stat);
    rst.push_back(df1);
    rst.push_back(df2);
    return 1;
}
int flignerTest(vector<double> &y,vector<double> &x, vector<double> &rst,double freq) //norm: rutrun 1
{
    // a little different from R. Because of median and rank.
    //in R. rank(c(1,2,2,3,4,5)) should be c(1,2.5,2.5,3,4,5)
    //here, would be c(1,2,2,3,4,5)
    
    int N=(int)y.size();
    vector<int> cat(N);
    #pragma omp parallel for
    for(int i=0;i<x.size();i++) cat[i]=(int)x[i];
    getUnique(cat);
    int k=(int)cat.size();
    if(k==1) return -1;
    vector<int> n(k);
    vector<double> sigmai(k),catd(k);
    vector<double> ymedian(k);
    vector<int> idx;
    vector<double> tmp;
    for(int i=0;i<k;i++) {
        idx.clear();
        tmp.clear();
        n[i]=0;
        int category=cat[i];
        catd[i]=(double)cat[i];
        for(int j=0;j<N;j++) if((int)x[j]==category) { n[i]++; idx.push_back(j); tmp.push_back(y[j]);}
        sigmai[i]=var(tmp);
        ymedian[i]=median(tmp);
        #pragma omp parallel for
        for(int j=0;j<idx.size();j++) {
            double absmu=abs(y[idx[j]]-ymedian[i]);
            y[idx[j]]=absmu; //now y is z_ij
        }
    }
    vector<double> rankzij;
    getRank2R(y, rankzij);
    #pragma omp parallel for
    for(int i=0;i<N;i++) {
        y[i]=qnorm((1+rankzij[i]/(N+1))/2);
    }
    double z=mean(y);
    for(int i=0;i<k;i++) {
        int category=cat[i];
        tmp.clear();
        for(int j=0;j<N;j++) if((int)x[j]==category) { tmp.push_back(y[j]);}
        ymedian[i]=mean(tmp);
    }
    double stat=0.0;
    for(int i=0;i<k;i++)
    {
        stat+=n[i]*(ymedian[i]-z)*(ymedian[i]-z);
    }
    stat/=var(y);
    double p=pchisq(stat, k-1);
    
    double zest=sqrt(qchisq(p,1));
    double domin=sqrt(2*freq*(1-freq)*(N+zest*zest));
    double best=zest/domin;
    double seest=1/domin;
    reg(sigmai,catd,rst);
    if(rst[0]<0) best *= -1.0; //updated by futao.zhang 20180321
    rst.clear();
    rst.push_back(best);
    rst.push_back(seest);
    rst.push_back(p);
    rst.push_back(stat);
    rst.push_back(k-1);
    return 1;
}
void  getResidual(VectorXd &y, MatrixXd &_X)
{
    if(y.size()!=_X.rows() || y.size()<1) {
        LOGPRINTF("Error: The row number of C and the length of y do not match.\n");
        TERMINATE();
    }
        vector<double> cvec;
        vector<double> xvec;
        vector<int> NMISS;
        MatrixXd X=_X;
        double nonmiss=0.0;
        int miss=0;
        for(int j=0; j<y.size(); j++)
        {
            double val=y(j);
            if(val<1e9){
                NMISS.push_back(j);
                xvec.push_back(val);
                nonmiss+=1.0;
            } else {
                removeRow(X, j-miss);
                miss++;
            }
        }
        
        VectorXd x(xvec.size());
        for(int j=0;j<xvec.size();j++) x(j)=xvec[j];
    
        MatrixXd XtX_i;
        XtX_i=X.transpose()*X;
        bool determinant_zero=false;
        inverse_V(XtX_i, determinant_zero);
        if(determinant_zero) {
            LOGPRINTF("ERROR: Maybe there is multicollinearity in the covariates.\n");
            TERMINATE();
        }
        
        VectorXd b_hat=XtX_i*X.transpose()*x;
        VectorXd residual=(x-X*b_hat);
        for(int j=0;j<residual.size();j++) y[NMISS[j]]=residual(j);
}
void  getResidual(MatrixXd &Y, MatrixXd &COV)
{
    //Y is the geno matrix
    //COV has no missing here
    if(Y.rows()!=COV.rows() || Y.rows()<1) {
        LOGPRINTF("Error: The row number of C and Y not match.\n");
        TERMINATE();
    }
    for(int i=0;i<=Y.cols();i++)
    {
        vector<double> xvec;
        vector<int> NMISS;
        MatrixXd X=COV;
        double nonmiss=0.0;
        int miss=0;
        for(int j=0; j<Y.rows(); j++)
        {
            double bval=Y(i,j);
            if(bval<1e5){
                NMISS.push_back(j);
                xvec.push_back(bval);
                nonmiss+=1.0;
            } else {
                removeRow(X, j-miss);
                miss++;
            }
        }
        
        VectorXd x(xvec.size());
        for(int j=0;j<xvec.size();j++) x(j)=xvec[j];
        
        if(x.size()!=X.rows() || x.size()<1) {
            LOGPRINTF("Error: The row number of C and the length of y do not match.\n");
            TERMINATE();
        }
        
        MatrixXd XtX_i;
        XtX_i=X.transpose()*X;
        bool determinant_zero=false;
        inverse_V(XtX_i, determinant_zero);
        if(determinant_zero) {
            LOGPRINTF("ERROR: Maybe there is multicollinearity in the covariates.\n");
            TERMINATE();
        }
        
        VectorXd b_hat=XtX_i*X.transpose()*x;
        VectorXd residual=(x-X*b_hat);
        for(int j=0;j<residual.size();j++) Y(i,NMISS[j])=residual(j);
    }
}
void  fast_getResidual(MatrixXd &Y, MatrixXd &X_XtXi_Xt)
{
    //Y must have no missing
    if(Y.rows()!=X_XtXi_Xt.rows() || Y.rows()<1) {
        LOGPRINTF("Error: The row number of C and Y not match.\n");
        TERMINATE();
    }
    Y=Y-X_XtXi_Xt*Y;
//    for(int i=0;i<Y.cols();i++)
//    {
//        VectorXd x=Y.col(i);
//        VectorXd residual=(x-X_XtXi_Xt*x);
//        Y.col(i)=residual;
//    }
}
void  fast_getResidual(MatrixXd &Y, MatrixXd &X,MatrixXd &XtXi)
{
    //Y must have no missing
    if(Y.rows()!=X.rows() || Y.rows()<1) {
        LOGPRINTF("Error: The row number of C and Y not match.\n");
        TERMINATE();
    }
   MatrixXd XtXiXt=XtXi*X.transpose();
        for(int i=0;i<Y.cols();i++)
        {
            VectorXd x=Y.col(i);
            VectorXd b_hat=XtXiXt*x;
            Y.col(i)=(x-X*b_hat);
        }
}
void mlma_cal_stat(VectorXd &_Y, VectorXd &_x, MatrixXd &_Vi, double &beta, double &se, double &pval)
{
    if(_Y.size()!=_x.size())
    {
        LOGPRINTF("Error: The size of Y and the size of X do not match.\n");
        TERMINATE();
    }
    if(_Y.size()!=_Vi.rows())
    {
        LOGPRINTF("Error: The size of Y and the rows of Vi do not match.\n");
        TERMINATE();
    }
    double Xt_Vi_X=0.0, chisq=0.0;
    uint64_t n=_Y.size();
    VectorXd X(n);
    VectorXd Vi_X(n);
    double nonmiss=0.0;
    double mu=0.0;
    //calc the mean
    for(int j = 0; j < n; j++) {
        double val=_x(j);
        if(val<1e9){
            mu+=val;
            nonmiss+=1.0;
        }
    }
    mu/=nonmiss;
    for(int j = 0; j < n; j++) {
        double val=_x(j);
        if(val<1e9){
            X(j) = val-mu;
        } else {
            X(j)=0.0;
        }
    }
    Vi_X=_Vi*X;
    Xt_Vi_X=Vi_X.dot(X);
    se=1.0/Xt_Vi_X;
    beta=se*(Vi_X.dot(_Y));
    if(se>1.0e-30){
        se=sqrt(se);
        chisq=beta/se;
        pval=pchisq(chisq*chisq, 1);
    }
}
bool comput_inverse_logdet_LU(MatrixXd &Vi, double &logdet)
{
    long n = Vi.cols();
    
    FullPivLU<MatrixXd> lu(Vi);
    if (!lu.isInvertible()) return false;
    VectorXd u = lu.matrixLU().diagonal();
    logdet = 0.0;
    for (int i = 0; i < n; i++) logdet += log(fabs(u[i]));
    Vi = lu.inverse();
    return true;
}

void mlma_cal_stat_covar(VectorXd &_Y,VectorXd &_x,  MatrixXd &_Vi,  MatrixXd &_X, double &beta, double &se, double &pval)
{
    
    long _X_c = _X.cols();
    unsigned long col_num=_X_c+1;
    double chisq=0.0;
    uint64_t n=_Y.size();

    MatrixXd X(n,col_num);
    MatrixXd Vi_X(n,col_num);
    MatrixXd Xt_Vi_X(col_num,col_num);
    VectorXd Xt_Vi_y(col_num);
    VectorXd b_vec(col_num);
    
    X.block(0, 0, n, _X_c)=_X;
    
        double nonmiss=0.0;
        double mu=0.0;
        //calc the mean
        for(int j = 0; j < n; j++) {
            double val=_x(j);
            if(val<1e9){
                mu+=val;
                nonmiss+=1.0;
            }
        }
        mu/=nonmiss;
        for(int j = 0; j < n; j++) {
            double val=_x(j);
            if(val<1e9){
                X(j,_X_c) = val-mu;
            } else {
                X(j,_X_c)=0.0;
            }
        }
        Vi_X=_Vi*X;
        Xt_Vi_X=X.transpose()*Vi_X;
        double logdt=0.0;
        if(!comput_inverse_logdet_LU( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
        Xt_Vi_y=Vi_X.transpose()*_Y;
        b_vec=Xt_Vi_X*Xt_Vi_y;
        
        se=Xt_Vi_X(_X_c,_X_c);
        beta=b_vec[_X_c];
        if(se>1.0e-30){
            se=sqrt(se);
            chisq=beta/se;
            pval=pchisq(chisq*chisq, 1);
        }
}
void  LR(VectorXd &_Y,VectorXd &_x,  MatrixXd &_X, double &beta, double &se, double &pval)
{
    long _n=_Y.size();
    vector<double> yvec;
    vector<double> cvec;
    vector<double> xvec;
    MatrixXd X=_X;
    double nonmiss=0.0;
    int miss=0;
    for(int j=0; j<_n; j++)
    {
        double phval=_Y(j);
        double val=_x(j);
        if(val<1e9){
            yvec.push_back(phval);
            xvec.push_back(val);
            nonmiss+=1.0;
        } else {
            removeRow(X, j-miss);
            miss++;
        }
    }
    VectorXd y(yvec.size()),x(xvec.size());
    for(int j=0;j<xvec.size();j++)
    {
        y(j)=yvec[j];
        x(j)=xvec[j];
    }
    vector<double> rst;
    bool notInvertible=lin(y, X, x, rst);
    if(notInvertible) {
        LOGPRINTF("ERROR: the X’X matrix is not invertible likely because of the multicollinearity in the covariates.\n");
        TERMINATE();
    }
    beta=rst[0];
    se=rst[1];
    pval=rst[2];
}

double SQR(double a)
{
    return a*a;
}


double pythag(const double a, const double b)
{
    double absa,absb;
    
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

bool svdcmp(MatrixXd & a,
            vector<double> & w,
            vector<vector<double> > &v)
{
    bool flag;
    int i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    double volatile temp;
    
    int m=a.rows();
    if (m==0){
        printf("Internal problem in SVD function (no observations left?)");
        exit(EXIT_FAILURE);
    }
    int n=a.cols();
    
    vector<double> rv1(n);
    g=scale=anorm=0.0;
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += fabs(a(k,i));
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    a(k,i) /= scale;
                    s += a(k,i)*a(k,i);
                }
                f=a(i,i);
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a(i,i)=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
                    f=s/h;
                    for (k=i;k<m;k++) a(k,j) += f*a(k,i);
                }
                for (k=i;k<m;k++) a(k,i) *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += fabs(a(i,k));
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    a(i,k) /= scale;
                    s += a(i,k)*a(i,k);
                }
                f=a(i,l-1);
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a(i,l-1)=f-g;
                for (k=l-1;k<n;k++) rv1[k]=a(i,k)/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += a(j,k)*a(i,k);
                    for (k=l-1;k<n;k++) a(j,k) += s*rv1[k];
                }
                for (k=l-1;k<n;k++) a(i,k) *= scale;
            }
        }
        anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++)
                    v[j][i]=(a(i,j)/a(i,l))/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += a(i,k)*v[k][j];
                    for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=min(m,n)-1;i>=0;i--) {
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) a(i,j)=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += a(k,i)*a(k,j);
                f=(s/a(i,i))*g;
                for (k=i;k<m;k++) a(k,j) += f*a(k,i);
            }
            for (j=i;j<m;j++) a(j,i) *= g;
        } else for (j=i;j<m;j++) a(j,i)=0.0;
        ++a(i,i);
    }
    for (k=n-1;k>=0;k--) {
        for (its=0;its<30;its++) {
            flag=true;
            for (l=k;l>=0;l--) {
                nm=l-1;
                temp=fabs(rv1[l])+anorm;
                if (temp == anorm) {
                    flag=false;
                    break;
                }
                temp=fabs(w[nm])+anorm;
                if (temp == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    temp = fabs(f)+anorm;
                    if (temp == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=a(j,nm);
                        z=a(j,i);
                        a(j,nm)=y*c+z*s;
                        a(j,i)=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 29) 
                return false; // cannot converge: multi-collinearity?
            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=a(jj,j);
                    z=a(jj,i);
                    a(jj,j)=y*c+z*s;
                    a(jj,i)=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    return true;
}

bool svdcmp(vector<vector<double> > & a,
            vector<double> & w,
            vector<vector<double> > &v)
{
    bool flag;
    int i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    double volatile temp;
    
    int m=a.size();
    if (m==0) printf("Internal problem in SVD function (no observations left?)");
    int n=a[0].size();
    
    vector<double> rv1(n);
    g=scale=anorm=0.0;
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += fabs(a[k][i]);
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l-1];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l-1]=f-g;
                for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
                    for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l-1;k<n;k++) a[i][k] *= scale;
            }
        }
        anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++)
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=MIN(m,n)-1;i>=0;i--) {
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) a[i][j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<m;j++) a[j][i] *= g;
        } else for (j=i;j<m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n-1;k>=0;k--) {
        for (its=0;its<30;its++) {
            flag=true;
            for (l=k;l>=0;l--) {
                nm=l-1;
                temp=fabs(rv1[l])+anorm;
                if (temp == anorm) {
                    flag=false;
                    break;
                }
                temp=fabs(w[nm])+anorm;
                if (temp == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    temp = fabs(f)+anorm;
                    if (temp == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 29)
                return false; // cannot converge: multi-collinearity?
            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    return true;
}
vector< vector<double> > svd_inverse(vector< vector<double> > & u , bool & flag )
{
    
    const double eps = 1e-24;
    
    if (u.size() == 0)
        printf("Internal problem: matrix with no rows (inverse function)");
    if (u.size() != u[0].size() )
        printf("Internal problem: Cannot invert non-square matrix");
    int n = u.size();
    
    vector<double> w(n,0);
    
    vector<vector<double> > v(n);
    for (int i=0; i<n; i++)
        v[i].resize(n,0);
    
    flag = svdcmp(u,w,v);
    
    // Look for singular values
    double wmax = 0;
    for (int i=0; i<n; i++)
        wmax = w[i] > wmax ? w[i] : wmax;
    double wmin = wmax * eps;
    for (int i=0; i<n; i++)
    {
        w[i] = w[i] < wmin ? 0 : 1/w[i];
    }
    
    
    // u w t(v)
    
    // row U * 1/w
    
    // results matrix
    vector<vector<double> > r(n);
    for (int i=0; i<n; i++)
    {
        r[i].resize(n,0);
        for (int j=0; j<n; j++)
            u[i][j] = u[i][j] * w[j];
    }
    
    // [nxn].[t(v)]
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            for (int k=0; k<n; k++)
                r[i][j] += u[i][k] * v[j][k];
    
    return r;
}


MatrixXd svd_inverse(MatrixXd & u , bool & flag )
{
    
    const double eps = 1e-24;
    long n = u.rows();
    
    vector<double> w(n,0);
    vector<vector<double> > v(n);
    for (int i=0; i<n; i++)
        v[i].resize(n,0);
    flag = svdcmp(u,w,v);
    
    // Look for singular values
    double wmax = 0;
    for (int i=0; i<n; i++)
        wmax = w[i] > wmax ? w[i] : wmax;
    double wmin = wmax * eps;
    for (int i=0; i<n; i++)
    {
        w[i] = w[i] < wmin ? 0 : 1/w[i];
    }
    
    
    // u w t(v)
    
    // row U * 1/w
    
    // results matrix
    MatrixXd r=MatrixXd::Zero(n,n);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
            u(i,j) = u(i,j) * w[j];
    }
    
    // [nxn].[t(v)]
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            for (int k=0; k<n; k++)
                r(i,j) += u(i,k) * v[j][k];
    
    return r;
}

void LogisticReg(VectorXi &_Y,  MatrixXd &_X, double &beta, double &se, double &pval)
{
    long _n=_Y.size();
    long m=_X.cols();
    VectorXd coef=VectorXd::Zero(m), p(_n),V(_n);
    MatrixXd S(m,m);
    
    // Newton-Raphson to fit logistic model
    
    bool converge = false;
    int it = 0;
    
    while ( ! converge && it < 20 )
    {
        
        // Determine p and V
        for (int i=0; i<_n; i++)
        {
            double t = 0;
            for (int j=0; j<m; j++)
                t += coef[j] * _X(i,j);
            p(i) = 1/(1+exp(-t));
            V(i) = p(i) * (1-p(i));
        }
        
        MatrixXd T(m,m);
        
        for (int j=0; j<m; j++)
            for (int k=j; k<m; k++)
            {
                double sum = 0;
                for (int i=0; i<_n; i++)
                    sum += _X(i,j) * V[i] * _X(i,k) ;
                T(j,k) = T(k,j) = sum;
            }
        /****/
//        string filename=string(outfileName)+".mat";
//        FILE* tmpfile=fopen(filename.c_str(),"w");
//        if(!tmpfile)
//        {
//            LOGPRINTF("error open file.\n");
//            TERMINATE();
//        }
//        for (int j=0; j<m; j++)
//        {
//            string str="";
//            for (int k=0; k<m; k++)
//            {
//                
//                str+=atos(T(j,k))+"\t";
//            }
//            str += '\n';
//            fputs(str.c_str(),tmpfile);
//        }
//        fclose(tmpfile);
        /****/
        /*
        bool detzero =false;
        inverse_V(T, detzero);
        if(detzero) {
            printf("ERROR: Maybe there is multicollinearity in the covariates.\n");
            exit(1);
        }
        */
        bool flag = true;
        T = svd_inverse(T,flag);
        if(!flag)
        {
            printf("Warning: fitted probabilities numerically 0 or 1 occurred.\n");
            //exit(1);
            beta=0;
            se=-9;
            pval=-9;
            return;
        }
        
        MatrixXd T2 = T*_X.transpose();
        
        vector<double> t3(_n);
        for (int i=0; i<_n; i++)
            t3[i] = _Y(i) - p(i);
        
        vector<double>  ncoef(m);
        for (int j=0; j<m; j++)
        {
            ncoef[j]=0;
            for (int i=0; i<_n; i++)
                ncoef[j] += T2(j,i) * t3[i];
        }
        
        
        double delta = 0;
        for (int j=0; j<m; j++)
        {
            delta += abs(ncoef[j]);
            coef[j] += ncoef[j];
        }
        
        if ( delta < 1e-6 )
            converge = true;
        
        it++;
    }
    
    MatrixXd Xt = _X.transpose()*V.asDiagonal();
    S=Xt*_X;
    /*
    bool flag = false;
    inverse_V(S, flag);
    if(flag) {
        printf("ERROR: Maybe there is multicollinearity in the covariates.\n");
        exit(1);
    }
    */
    bool flag = true;
    S = svd_inverse(S,flag);
    if(!flag)
    {
        printf("Warning: fitted probabilities numerically 0 or 1 occurred.\n");
        //exit(1);
        beta=0;
        se=-9;
        pval=-9;
        return;
    }
    
    beta=coef(m-1);
    se=sqrt(S(m-1,m-1));
    double z=beta/se;
    pval=pchisq(z*z, 1);
    
}

void mlm_stat_covar(VectorXd &_Y,  MatrixXd &_Vi,  MatrixXd &_X, double &beta, double &se, double &pval)
{
    //_X includes covariates and the probe of interest
    long n=_X.rows(), _X_c=_X.cols()-1;
    if(_Y.size()!=n || _Vi.cols()!=n || _Vi.rows()!=n)
    {
        LOGPRINTF("ERROE: ");
        TERMINATE();
        
    }
    MatrixXd Vi_X=_Vi*_X;
    MatrixXd Xt_Vi_X=_X.transpose()*Vi_X;
    double logdt=0.0;
    if(!comput_inverse_logdet_LU( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
    VectorXd Xt_Vi_y=Vi_X.transpose()*_Y;
    VectorXd b_vec=Xt_Vi_X*Xt_Vi_y;
    
    se=Xt_Vi_X(_X_c,_X_c);
    beta=b_vec[_X_c];
    if(se>1.0e-30){
        se=sqrt(se);
        double chisq=beta/se;
        pval=pchisq(chisq*chisq, 1);
    }
}
void mlm_stat(VectorXd &_y,MatrixXd &_Vi, VectorXd &_x, double &beta, double &se, double &pval)
{
    if(_y.size()!=_x.size() || _y.size()!=_Vi.rows() || _y.size()!=_Vi.cols())
    {
        LOGPRINTF("ERROR: ");
        TERMINATE();
    }
    VectorXd Vi_X=_Vi*_x;
    double Xt_Vi_X=Vi_X.dot(_x);
    se=1.0/Xt_Vi_X;
    beta=se*(Vi_X.dot(_y));
    if(se>1.0e-30){
        se=sqrt(se);
        double chisq=beta/se;
        pval=pchisq(chisq*chisq, 1);
    }
}


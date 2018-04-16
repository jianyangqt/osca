//
//  l1_stat.cpp
//  osc
//
//  Created by Futao Zhang on 23/12/2016.
//  Copyright © 2016 Futao Zhang. All rights reserved.
//

#include "l1_stat.hpp"

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
    getRank(y, rankzij);
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
    return 1;
}


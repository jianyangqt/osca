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

void lin(VectorXd &y, MatrixXd &C, VectorXd &x, vector<double> &rst)
{
    long N=y.size();
    if(N!=x.size() || N<1) {
        LOGPRINTF("Error: The lengths of x and y do not match.\n");
        exit(EXIT_FAILURE);
    }
    if(N!=C.rows() || N<1) {
        LOGPRINTF("Error: The row number of C and the length of y do not match.\n");
        exit(EXIT_FAILURE);
    }
    
    MatrixXd X(C.rows(), C.cols()+1);
   
    X.block(0,0,C.rows(),C.cols())=C;
    long x_idx=X.cols()-1;
    X.col(x_idx)=x;
    MatrixXd XtX_i;
    XtX_i=(X.transpose()*X).inverse();
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



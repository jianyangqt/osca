//
//  l1_stat.hpp
//  osc
//
//  Created by Futao Zhang on 23/12/2016.
//  Copyright © 2016 Futao Zhang. All rights reserved.
//

#ifndef l1_stat_hpp
#define l1_stat_hpp

#include "l0_stat.h"
#endif /* l1_stat_hpp */

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

void reg(VectorXd &y, VectorXd &x,vector<double> &rst);
MatrixXd reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table = false);
bool lin(VectorXd &y, MatrixXd &C, VectorXd &x, vector<double> &rst);
void lin2(VectorXd &y, MatrixXd &X, vector<double> &rst);
double var(VectorXd &x);
void cov(MatrixXd &mat, MatrixXd &covar);
int bartlett(vector<double> &y,vector<double> &x, vector<double> &rst, double freq);
int leveneTest_mean(vector<double> &y,vector<double> &x, vector<double> &rst, double freq);
int leveneTest_median(vector<double> &y,vector<double> &x, vector<double> &rst, double freq) ;
int flignerTest(vector<double> &y,vector<double> &x, vector<double> &rst, double freq);
void  getResidual(VectorXd &y, MatrixXd &_X);
bool comput_inverse_logdet_LU(MatrixXd &Vi, double &logdet);
void mlma_cal_stat(VectorXd &_Y, VectorXd &_x, MatrixXd &_Vi, double &beta, double &se, double &pval);
void mlma_cal_stat_covar(VectorXd &_Y,VectorXd &_x,  MatrixXd &_Vi,  MatrixXd &_X, double &beta, double &se, double &pval);
void LR(VectorXd &_Y,VectorXd &_x,  MatrixXd &_X, double &beta, double &se, double &pval);
void  fast_getResidual(MatrixXd &Y, MatrixXd &X_XtXi_Xt);
void  fast_getResidual(MatrixXd &Y, MatrixXd &X,MatrixXd &XtXi);
void adjusted_reg(VectorXd &yresi, VectorXd &xresi,vector<double> &rst, int exdf);
void LogisticReg(VectorXi &_Y,  MatrixXd &_X, double &beta, double &se, double &pval);
void mlm_stat_covar(VectorXd &_Y,  MatrixXd &_Vi,  MatrixXd &_X, double &beta, double &se, double &pval);
void mlm_stat(VectorXd &_y,MatrixXd &_Vi, VectorXd &_x, double &beta, double &se, double &pval);
MatrixXd svd_inverse(MatrixXd & u , bool & flag );
vector< vector<double> > svd_inverse(vector< vector<double> > & u , bool & flag );

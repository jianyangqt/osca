//
//  l3_mlma.cpp
//  osc
//
//  Created by Futao Zhang on 13/04/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l2_reml.h"

namespace EFILE {

void init_varcomp(eInfo* einfo, vector<double> &reml_priors_var, vector<double> &reml_priors, VectorXd &varcmp) {
        int i = 0, pos = 0;
        double d_buf = 0.0;
        
        varcmp = VectorXd::Zero(einfo->_r_indx.size());
        if (einfo->_bivar_reml) {
            if (!reml_priors_var.empty()) {
                for (i = 0; i < einfo->_r_indx.size(); i++) varcmp[i] = reml_priors_var[i];
            }
            else if (!reml_priors.empty()) {
                for (i = 0, d_buf = 0; i < einfo->_bivar_pos[0].size() - 1; i++) {
                    pos = einfo->_bivar_pos[0][i];
                    varcmp[pos] = reml_priors[pos] * einfo->_y_Ssq;
                    d_buf += reml_priors[pos];
                }
                if (d_buf > 1.0) throw ("\nError: --reml-priors. The sum of all prior values for trait 1 should not exceed 1.0.");
                varcmp[einfo->_bivar_pos[0][einfo->_bivar_pos[0].size() - 1]] = (1.0 - d_buf) * einfo->_y_Ssq;
                for (i = 0, d_buf = 0; i < einfo->_bivar_pos[1].size() - 1; i++) {
                    pos = einfo->_bivar_pos[1][i];
                    varcmp[pos] = reml_priors[pos] * einfo->_y_Ssq;
                    d_buf += reml_priors[pos];
                }
                if (d_buf > 1.0) throw ("\nError: --reml-priors. The sum of all prior values for trait 2 should not exceed 1.0.");
                varcmp[einfo->_bivar_pos[1][einfo->_bivar_pos[1].size() - 1]] = (1.0 - d_buf) * einfo->_y2_Ssq;
                for (i = 0; i < einfo->_bivar_pos[2].size(); i++) varcmp[einfo->_bivar_pos[2][i]] = reml_priors[einfo->_bivar_pos[2][i]] * sqrt(einfo->_y_Ssq * einfo->_y2_Ssq);
            }
            else {
                for (i = 0; i < einfo->_bivar_pos[0].size(); i++) varcmp[einfo->_bivar_pos[0][i]] = einfo->_y_Ssq / einfo->_bivar_pos[0].size();
                for (i = 0; i < einfo->_bivar_pos[1].size(); i++) varcmp[einfo->_bivar_pos[1][i]] = einfo->_y2_Ssq / einfo->_bivar_pos[1].size();
                for (i = 0; i < einfo->_bivar_pos[2].size(); i++) varcmp[einfo->_bivar_pos[2][i]] = 0.5 * sqrt(varcmp[einfo->_bivar_pos[0][i]] * varcmp[einfo->_bivar_pos[1][i]]);
            }
            
            return;
        }
        
        if (!reml_priors_var.empty()) {
            for (i = 0; i < einfo->_r_indx.size() - 1; i++) varcmp[i] = reml_priors_var[i];
            if (reml_priors_var.size() < einfo->_r_indx.size()) varcmp[einfo->_r_indx.size() - 1] = einfo->_y_Ssq - varcmp.sum();
            else varcmp[einfo->_r_indx.size() - 1] = reml_priors_var[einfo->_r_indx.size() - 1];
        }
        else if (!reml_priors.empty()) {
            for (i = 0, d_buf = 0; i < einfo->_r_indx.size() - 1; i++) {
                varcmp[i] = reml_priors[i] * einfo->_y_Ssq;
                d_buf += reml_priors[i];
            }
            if (d_buf > 1.0) throw ("\nError: --reml-priors. The sum of all prior values should not exceed 1.0.");
            varcmp[einfo->_r_indx.size() - 1] = (1.0 - d_buf) * einfo->_y_Ssq;
        }
        else varcmp.setConstant(einfo->_y_Ssq / (einfo->_r_indx.size()));
    }
    


void update_A(eInfo* einfo,VectorXd &prev_varcmp) {
    int i = 0;
    double g1 = 0.0, g2 = 0.0;
    for (i = 0; i < einfo->_fixed_rg_val.size(); i++) {
        g1 = prev_varcmp[einfo->_bivar_pos[0][i]];
        g2 = prev_varcmp[einfo->_bivar_pos[1][i]];
        einfo->_Asp[einfo->_bivar_pos_prev[0][i]] = einfo->_Asp_prev[einfo->_bivar_pos_prev[0][i]] + einfo->_Asp_prev[einfo->_bivar_pos_prev[2][i]]*(0.5 * einfo->_fixed_rg_val[i] * sqrt(g2 / g1));
        einfo->_Asp[einfo->_bivar_pos_prev[1][i]] = einfo->_Asp_prev[einfo->_bivar_pos_prev[1][i]] + einfo->_Asp_prev[einfo->_bivar_pos_prev[2][i]]*(0.5 * einfo->_fixed_rg_val[i] * sqrt(g1 / g2));
    }
    
}


bool comput_inverse_logdet_LDLT(MatrixXd &Vi, double &logdet) {
    int i = 0, n = Vi.cols();
    LDLT<MatrixXd> ldlt(Vi);
    VectorXd d = ldlt.vectorD();
    
    if (d.minCoeff() < 0) return false;
    else {
        logdet = 0.0;
        for (i = 0; i < n; i++) logdet += log(d[i]);
        Vi.setIdentity();
        ldlt.solveInPlace(Vi);
    }
    return true;
}


bool bending_eigenval(VectorXd &eval) {
    int j = 0;
    double eval_m = eval.mean();
    if (eval.minCoeff() > 0.0) return false;
    double S = 0.0, P = 0.0;
    for (j = 0; j < eval.size(); j++) {
        if (eval[j] >= 0) continue;
        S += eval[j];
        P = -eval[j];
    }
    double W = S * S * 100.0 + 1;
    for (j = 0; j < eval.size(); j++) {
        if (eval[j] >= 0) continue;
        eval[j] = P * (S - eval[j])*(S - eval[j]) / W;
    }
    eval *= eval_m / eval.mean();
    return true;
}


void bend_V(MatrixXd &Vi)
{
    SelfAdjointEigenSolver<MatrixXd> eigensolver(Vi);
    VectorXd eval = eigensolver.eigenvalues();
    bending_eigenval(eval);
    eval.array() = 1.0 / eval.array();
    Vi = eigensolver.eigenvectors() * DiagonalMatrix<double, Dynamic, Dynamic>(eval) * eigensolver.eigenvectors().transpose();
}


bool calcu_Vi(eInfo* einfo, MatrixXd &Vi, VectorXd &prev_varcmp, double &logdet, int &iter,vector<MatrixXd> &_A)
{
    int i = 0;
    string errmsg = "\nError: the V (variance-covariance) matrix is not invertible.";
    int _n=(int)einfo->_eii_include.size();
    Vi = MatrixXd::Zero(_n, _n);
   
    if (einfo->_r_indx.size() == 1) {
        Vi.diagonal() = VectorXd::Constant(_n, 1.0 / prev_varcmp[0]);
        logdet = _n * log(prev_varcmp[0]);
    }
     else {
        for (i = 0; i < einfo->_r_indx.size(); i++) Vi += (_A[einfo->_r_indx[i]]) * prev_varcmp[i];
        
        if (einfo->_V_inv_mtd == 0) {
            if (!comput_inverse_logdet_LDLT(Vi, logdet)) {
                if(einfo->_reml_force_inv) {
                    LOGPRINTF("Warning: the variance-covaraince matrix V is non-positive definite.\n" );
                    einfo->_V_inv_mtd = 1;
                    LOGPRINTF("\nSwitching to the \"bending\" approach to invert V. This method hasn't been tested. The results might not be reliable!\n" );
                }
                /*else {
                 if(_reml_no_converge){
                 cout<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<endl;
                 double d_buf = Vi.diagonal().mean() * 0.01;
                 for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
                 if(!comput_inverse_logdet_LDLT_mkl(Vi, logdet)) return false;
                 }
                 else throw("Error: the variance-covaraince matrix V is not positive definite.");
                 }*/
            }
        }
        if (einfo->_V_inv_mtd == 1) bend_V(Vi);
         /*if (_V_inv_mtd == 2) {
          if(!_reml_force_converge){
          cout << "Switching from Cholesky to LU decomposition approach. The results might not be reliable!" << endl;
          if (!comput_inverse_logdet_LU_mkl(Vi, logdet)){
          if(_reml_no_converge){
          cout<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<endl;
          double d_buf = Vi.diagonal().mean() * 0.01;
          for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
          if(!comput_inverse_logdet_LDLT_mkl(Vi, logdet)) return false;
          }
          else throw ("Error: the variance-covaraince matrix V is still not invertible using LU decomposition.");
          }
          }
          else{
          cout<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<endl;
          double d_buf = Vi.diagonal().mean() * 0.01;
          for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
          comput_inverse_logdet_LU_mkl(Vi, logdet);
          }
          }*/
    }
    
    return true;
}


bool comput_inverse_logdet_LU(MatrixXd &Vi, double &logdet)
{
    int n = Vi.cols();
 
    FullPivLU<MatrixXd> lu(Vi);
    if (!lu.isInvertible()) return false;
    VectorXd u = lu.matrixLU().diagonal();
    logdet = 0.0;
    for (int i = 0; i < n; i++) logdet += log(fabs(u[i]));
    Vi = lu.inverse();
    return true;
}


void calcu_Vi_bivar(eInfo* einfo, MatrixXd &Vi, VectorXd &prev_varcmp, double &logdet, int &iter)
{
    int i = 0, n = Vi.cols();
    double d_buf = 0.0;
    logdet = 0.0;
    string errmsg = "\nError: the V (variance-covariance) matrix is not invertible.";
    int _n=(int)einfo->_eii_include.size();
    Vi = MatrixXd::Zero(_n, _n);
    for (i = 0; i < einfo->_r_indx.size(); i++) Vi += (einfo->_Asp[einfo->_r_indx[i]]) * prev_varcmp[i];
    
    if (einfo->_V_inv_mtd == 0) {
        if (!comput_inverse_logdet_LDLT(Vi, logdet)) {
            //cout<<"Note: the variance-covaraince matrix V is non-positive definite. Switching to Cholesky to LU decomposition approach."<<endl;
            einfo->_V_inv_mtd = 1;
        }
    }
    if (einfo->_V_inv_mtd == 1) {
        if (!comput_inverse_logdet_LU(Vi, logdet)) {
            LOGPRINTF("Error: the variance-covaraince matrix V is not invertible.\n");
            remlstatus=-1;
        }
    }
}


bool calcu_Vi_within_family(eInfo* einfo, MatrixXd &Vi, VectorXd &prev_varcmp, double &logdet, int &iter)
{
    int i=0;
    double logdet_buf=0.0;
    string errmsg="\nError: the V (variance-covariance) matrix is not invertible.";
    int _n=(int)einfo->_eii_include.size();
    Vi=MatrixXd::Zero(_n, _n);
    if(einfo->_r_indx.size()==1){
        Vi.diagonal()=VectorXd::Constant(_n, 1.0/prev_varcmp[0]);
        logdet=_n*log(prev_varcmp[0]);
    }
    else{
        for(i=0; i<einfo->_r_indx.size(); i++) Vi+=(einfo->_Asp[einfo->_r_indx[i]])*prev_varcmp[i];
        int prev_pnt=0, subn=0;
        logdet=0.0;
        for(i=0; i<einfo->_fam_brk_pnt.size()+1; i++){
            if(i==einfo->_fam_brk_pnt.size()) subn=_n-prev_pnt;
            else subn=einfo->_fam_brk_pnt[i]-prev_pnt+1;
            MatrixXd subVi=Vi.block(prev_pnt, prev_pnt, subn, subn);
            if(!comput_inverse_logdet_LDLT(subVi, logdet_buf)) {
                LOGPRINTF("Error: the sub-matrix of V for the %d-th family is not invertible.\n",i+1);
                remlstatus=-1;
            }
            logdet+=logdet_buf;
            //logdet+=comput_inverse_logdet_LU(subVi, errmsg.str());
            Vi.block(prev_pnt, prev_pnt, subn, subn)=subVi;
            
            if(i<einfo->_fam_brk_pnt.size()) prev_pnt=einfo->_fam_brk_pnt[i]+1;
        }
    }
    
    return true;
}


bool inverse_H(MatrixXd &H)
{
    double d_buf = 0.0;
    if (!comput_inverse_logdet_LDLT(H, d_buf)) return false;
    
    else return true;
}


// input P, calculate PA and Hi
void calcu_Hi(eInfo* einfo, MatrixXd &P, MatrixXd &Hi, vector<MatrixXd> &_A)
{
    int i = 0, j = 0, k = 0, l = 0;
    double d_buf = 0.0;
    int _n=(int)einfo->_eii_include.size();
    // Calculate PA
    vector<MatrixXd> PA(einfo->_r_indx.size());
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        (PA[i]).resize(_n, _n);
        if (einfo->_bivar_reml || einfo->_within_family) (PA[i]) = P * (einfo->_Asp[einfo->_r_indx[i]]);
        else (PA[i]) = P * (_A[einfo->_r_indx[i]]);
    }
    
    // Calculate Hi
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        for (j = 0; j <= i; j++) {
            d_buf = 0.0;
            for (k = 0; k < _n; k++) {
                for (l = 0; l < _n; l++) d_buf += (PA[i])(k, l)*(PA[j])(l, k);
            }
            Hi(i, j) = Hi(j, i) = d_buf;
        }
    }
    
    if (!inverse_H(Hi)){
        if(einfo->_reml_force_converge){
            cout << "Warning: the information matrix is not invertible." << endl;
            einfo->_reml_AI_not_invertible = true;
        }
        else {LOGPRINTF("Error: the information matrix is not invertible."); TERMINATE();}
    }
}


double calcu_P(MatrixXd &Vi, MatrixXd &Vi_X, MatrixXd &Xt_Vi_X_i, MatrixXd &P, MatrixXd &_X)
{
    Vi_X = Vi*_X;
    Xt_Vi_X_i = _X.transpose() * Vi_X;
    double logdet_Xt_Vi_X = 0.0;
    if(!comput_inverse_logdet_LU(Xt_Vi_X_i, logdet_Xt_Vi_X))
    {
        LOGPRINTF("\nError: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).");
        remlstatus=-1;
        return logdet_Xt_Vi_X;
    }
    P = Vi - Vi_X * Xt_Vi_X_i * Vi_X.transpose();
    return logdet_Xt_Vi_X;
}



// input P, calculate tr(PA)
void calcu_tr_PA(eInfo* einfo, MatrixXd &P, VectorXd &tr_PA,vector<MatrixXd> &_A) {
    int i = 0, k = 0, l = 0;
    double d_buf = 0.0;
    int _n=(int)einfo->_eii_include.size();
    
    // Calculate trace(PA)
    tr_PA.resize(einfo->_r_indx.size());
   
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        
        if (einfo->_bivar_reml || einfo->_within_family) tr_PA(i) = (P * (einfo->_Asp[einfo->_r_indx[i]])).matrix().diagonal().sum(); //modified 20160418
        else {
            d_buf = 0.0;
            for (k = 0; k < _n; k++) {
                for (l = 0; l < _n; l++) d_buf += P(k, l)*(_A[einfo->_r_indx[i]])(k, l);
            }
            tr_PA(i) = d_buf;
        }
    }
    
}



// use Fisher-scoring to estimate variance component
// input P, calculate PA, H, R and varcmp


void em_reml(eInfo* einfo, MatrixXd &P, VectorXd &Py, VectorXd &prev_varcmp, VectorXd &varcmp,VectorXd &_y,vector<MatrixXd> &_A)
{
    int i = 0;
    int _n=(int)einfo->_eii_include.size();
    // Calculate trace(PA)
    VectorXd tr_PA;
    calcu_tr_PA(einfo, P, tr_PA, _A);
    
    // Calculate R
    Py = P*_y;
    VectorXd R(einfo->_r_indx.size());
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        if (einfo->_bivar_reml || einfo->_within_family) R(i) = (Py.transpose()*(einfo->_Asp[einfo->_r_indx[i]]) * Py)(0, 0);
        else R(i) = (Py.transpose()*(_A[einfo->_r_indx[i]]) * Py)(0, 0);
    }
    
    // Calculate variance component
    for (i = 0; i < einfo->_r_indx.size(); i++) varcmp(i) = (prev_varcmp(i) * _n - prev_varcmp(i) * prev_varcmp(i) * tr_PA(i) + prev_varcmp(i) * prev_varcmp(i) * R(i)) / _n;
    
    // added by Jian Yang Dec 2014
    //varcmp = (varcmp.array() - prev_varcmp.array())*2 + prev_varcmp.array();
}


void reml_equation(eInfo* einfo, MatrixXd &P, MatrixXd &Hi, VectorXd &Py, VectorXd &varcmp,VectorXd &_y,vector<MatrixXd> &_A)
{
    // Calculate Hi
    int _n=(int)einfo->_eii_include.size();
    calcu_Hi(einfo, P, Hi, _A);
    if(einfo->_reml_AI_not_invertible) return;
    
    // Calculate R
    Py = P*_y;
    VectorXd R(einfo->_r_indx.size());
    for (int i = 0; i < einfo->_r_indx.size(); i++) {
        if (einfo->_bivar_reml || einfo->_within_family) R(i) = (Py.transpose()*(einfo->_Asp[einfo->_r_indx[i]]) * Py)(0, 0);
        else R(i) = (Py.transpose()*(_A[einfo->_r_indx[i]]) * Py)(0, 0);
    }
    
    // Calculate variance component
    varcmp = Hi*R;
    Hi = 2 * Hi; // for calculation of SE
}


void ai_reml(eInfo* einfo, MatrixXd &P, MatrixXd &Hi, VectorXd &Py, VectorXd &prev_varcmp, VectorXd &varcmp, double dlogL,VectorXd &_y,vector<MatrixXd> &_A)
{
    int i = 0, j = 0;
    
    Py = P*_y;
    int _n=(int)einfo->_eii_include.size();
    VectorXd cvec(_n);
    MatrixXd APy(_n, einfo->_r_indx.size());
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        if (einfo->_bivar_reml || einfo->_within_family) (APy.col(i)) = (einfo->_Asp[einfo->_r_indx[i]]) * Py;
        else (APy.col(i)) = (_A[einfo->_r_indx[i]]) * Py;
    }
    
    // Calculate Hi
    VectorXd R(einfo->_r_indx.size());
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        R(i) = (Py.transpose()*(APy.col(i)))(0, 0);
        cvec = P * (APy.col(i));
        Hi(i, i) = ((APy.col(i)).transpose() * cvec)(0, 0);
        for (j = 0; j < i; j++) Hi(j, i) = Hi(i, j) = ((APy.col(j)).transpose() * cvec)(0, 0);
    }
    Hi = 0.5 * Hi;
    
    // Calcualte tr(PA) and dL
    VectorXd tr_PA;
    calcu_tr_PA(einfo, P, tr_PA, _A);
    R = -0.5 * (tr_PA - R);
    
    // Calculate variance component
    if (!inverse_H(Hi)){
        if(einfo->_reml_force_converge){
            LOGPRINTF( "Warning: the information matrix is not invertible.\n");
            einfo->_reml_AI_not_invertible = true;
            return;
        }
        else { LOGPRINTF("Error: the information matrix is not invertible."); remlstatus=-1;}
    }
    
    VectorXd delta(einfo->_r_indx.size());
    delta = Hi*R;
    if (dlogL > 1.0) varcmp = prev_varcmp + 0.316 * delta;
    else varcmp = prev_varcmp + delta;
}


int constrain_varcmp(eInfo* einfo, VectorXd &varcmp) {
    int pos = 0;
    double delta = 0.0, constr_scale = 1e-6;
    int i = 0, num = 0;
    vector<int> constrain(einfo->_r_indx.size());
    
    if (einfo->_bivar_reml) {
        for (i = 0, num = 0; i < einfo->_bivar_pos[0].size(); i++) {
            pos = einfo->_bivar_pos[0][i];
            if (varcmp[pos] < 0) {
                delta += einfo->_y_Ssq * constr_scale - varcmp[pos];
                varcmp[pos] = einfo->_y_Ssq * constr_scale;
                constrain[i] = 1;
                num++;
            }
        }
        delta /= (einfo->_bivar_pos[0].size() - num);
        for (i = 0; i < einfo->_bivar_pos[0].size(); i++) {
            pos = einfo->_bivar_pos[0][i];
            if (constrain[pos] < 1 && varcmp[pos] > delta) varcmp[pos] -= delta;
        }
        
        for (i = 0, num = 0; i < einfo->_bivar_pos[1].size(); i++) {
            pos = einfo->_bivar_pos[1][i];
            if (varcmp[pos] < 0) {
                delta += einfo->_y_Ssq * constr_scale - varcmp[pos];
                varcmp[pos] = einfo->_y_Ssq * constr_scale;
                constrain[i] = 1;
                num++;
            }
        }
        delta /= (einfo->_bivar_pos[1].size() - num);
        for (i = 0; i < einfo->_bivar_pos[1].size(); i++) {
            pos = einfo->_bivar_pos[1][i];
            if (constrain[pos] < 1 && varcmp[pos] > delta) varcmp[pos] -= delta;
        }
        
        for (i = 0, num = 0; i < constrain.size(); i++) {
            if (constrain[i] == 1) num++;
            
        }
        return num;
    }
    
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        if (varcmp[i] < 0) {
            delta += einfo->_y_Ssq * constr_scale - varcmp[i];
            varcmp[i] = einfo->_y_Ssq * constr_scale;
            constrain[i] = 1;
            num++;
        }
    }
    delta /= (einfo->_r_indx.size() - num);
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        if (constrain[i] < 1 && varcmp[i] > delta) varcmp[i] -= delta;
    }
    
    return num;
}


void constrain_rg(eInfo* einfo, VectorXd &varcmp) {
    static int count = 0;
    int v_pos = 0, c_pos = 0, d = einfo->_bivar_pos[0].size() - einfo->_bivar_pos[2].size();
    if (einfo->_ignore_Ce) d--;
    MatrixXd G(2, 2);
    
    for (c_pos = einfo->_bivar_pos[2].size() - 1; c_pos >= 0; c_pos--) {
        v_pos = c_pos + d;
        G(0, 0) = varcmp[einfo->_bivar_pos[0][v_pos]];
        G(1, 1) = varcmp[einfo->_bivar_pos[1][v_pos]];
        G(0, 1) = G(1, 0) = varcmp[einfo->_bivar_pos[2][c_pos]];
        
        SelfAdjointEigenSolver<MatrixXd> eigensolver(G);
        VectorXd eval = eigensolver.eigenvalues();
        if (eval.minCoeff() <= 0.0) {
            if (count == 0) {
                cout << "Note: to constrain the correlation being from -1 to 1, a genetic (or residual) variance-covariance matrix is bended to be positive definite. In this case, the SE is unreliable." << endl;
                count++;
            }
            bending_eigenval(eval);
            G = eigensolver.eigenvectors() * DiagonalMatrix<double, Dynamic, Dynamic>(eval) * eigensolver.eigenvectors().transpose();
            varcmp[einfo->_bivar_pos[0][v_pos]] = G(0, 0);
            varcmp[einfo->_bivar_pos[1][v_pos]] = G(1, 1);
            varcmp[einfo->_bivar_pos[2][c_pos]] = G(0, 1);
        }
    }
}


double reml_iteration(eInfo* einfo, MatrixXd &Vi_X, MatrixXd &Xt_Vi_X_i, MatrixXd &Hi, VectorXd &Py, VectorXd &varcmp, bool prior_var_flag, bool no_constrain, bool reml_bivar_fix_rg,MatrixXd &_Vi,vector<MatrixXd> &_A, MatrixXd &_X,VectorXd &_y)
{
 
    char *mtd_str[3] = {"AI-REML", "Fisher-scoring REML", "EM-REML"};
    int i = 0, constrain_num = 0, iter = 0, reml_mtd_tmp = einfo->_reml_mtd;
    double logdet = 0.0, logdet_Xt_Vi_X = 0.0, prev_lgL = -1e20, lgL = -1e20, dlogL = 1000.0;
    VectorXd prev_prev_varcmp(varcmp), prev_varcmp(varcmp), varcomp_init(varcmp);
    bool converged_flag = false;
    for (iter = 0; iter < einfo->_reml_max_iter; iter++) {
        if (reml_bivar_fix_rg) update_A(einfo, prev_varcmp);
        if (iter == 0) {
            prev_varcmp = varcomp_init;
            if (prior_var_flag){
                if(einfo->_reml_fixed_var) cout << "Variance components are fixed at: " << varcmp.transpose() << endl;
                else cout << "Prior values of variance components: " << varcmp.transpose() << endl;
            }
            else {
                einfo->_reml_mtd = 2;
                LOGPRINTF("Calculating prior values of variance components by EM-REML ...\n");
            }
        }
        if (iter == 1) {
            einfo->_reml_mtd = reml_mtd_tmp;
            LOGPRINTF("Running %s algorithm ...\nIter.\tlogL\t",mtd_str[einfo->_reml_mtd]);
            for (i = 0; i < einfo->_r_indx.size(); i++) {LOGPRINTF("%s\t",einfo->_var_name[einfo->_r_indx[i]].c_str())};
            LOGPRINTF("\n");
        }
        if (einfo->_bivar_reml) calcu_Vi_bivar(einfo, _Vi, prev_varcmp, logdet, iter); // Calculate Vi, bivariate analysis
        else if (einfo->_within_family) calcu_Vi_within_family(einfo, _Vi, prev_varcmp, logdet, iter); // within-family REML
        else {
            if (!calcu_Vi(einfo, _Vi, prev_varcmp, logdet, iter, _A )){ // Calculate Vi
                LOGPRINTF("Warning: V matrix is not positive-definite.\n");
                varcmp = prev_prev_varcmp;
                if(!calcu_Vi(einfo, _Vi, varcmp, logdet, iter,_A)) {LOGPRINTF("Error: V matrix is not positive-definite."); remlstatus=-1;}
                calcu_Hi(einfo, einfo->_P, Hi,_A);
                Hi = 2 * Hi;
                break;
            }
        }
        if(remlstatus!=0) return lgL;
        logdet_Xt_Vi_X = calcu_P(_Vi, Vi_X, Xt_Vi_X_i, einfo->_P,_X); // Calculate P
        if(remlstatus!=0) return lgL;
        
        if (einfo->_reml_mtd == 0) ai_reml(einfo, einfo->_P, Hi, Py, prev_varcmp, varcmp, dlogL, _y,_A );
        else if (einfo->_reml_mtd == 1) reml_equation(einfo, einfo->_P, Hi, Py, varcmp, _y, _A);
        else if (einfo->_reml_mtd == 2) em_reml(einfo, einfo->_P, Py, prev_varcmp, varcmp, _y,_A);
        lgL = -0.5 * (logdet_Xt_Vi_X + logdet + (_y.transpose() * Py)(0, 0));
        if(remlstatus!=0) return lgL;
        if(einfo->_reml_force_converge && einfo->_reml_AI_not_invertible) break;
       
        
        // output log
        if (!no_constrain) constrain_num = constrain_varcmp(einfo, varcmp);
        if (einfo->_bivar_reml && !einfo->_bivar_no_constrain) constrain_rg(einfo, varcmp);
        if (iter > 0) {
            //cout << iter << "\t" << setiosflags(ios::fixed) << setprecision(2) << lgL << "\t";
            LOGPRINTF("%d\t%.2f\t",iter,lgL);
            for (i = 0; i < einfo->_r_indx.size(); i++)
            {
                LOGPRINTF("%.5f\t",varcmp[i]);
            }//cout << setprecision(5) << varcmp[i] << "\t";
            
            if (constrain_num > 0) {
                LOGPRINTF("(%d component(s) constrained)\n",constrain_num);
            } else {
                LOGPRINTF("\n");
            }
        } else {
            if (!prior_var_flag) //cout << "Updated prior values: " << varcmp.transpose() << endl;
            {
                LOGPRINTF("Updated prior values: ");
                for(int k=0;k<varcmp.size();k++) {
                    LOGPRINTF("%f\t",varcmp(k));
                }
                LOGPRINTF("\n");
            }
            LOGPRINTF("logL: %.3f\n", lgL);
            //if(_reml_max_iter==1) cout<<"logL: "<<lgL<<endl;
        }
        if(einfo->_reml_fixed_var){
            varcmp = prev_varcmp;
            break;
        }
        if (constrain_num * 2 > einfo->_r_indx.size())
        {
            LOGPRINTF("Error: analysis stopped because more than half of the variance components are constrained. The result would be unreliable.\n Please have a try to add the option --reml-no-constrain.\n");
            remlstatus = -2;
        }
        // added by Jian Yang on 22 Oct 2014
        //if (constrain_num == _r_indx.size()) throw ("Error: analysis stopped because all variance components are constrained. You may have a try of adding the option --reml-no-constrain.");
        
        if((einfo->_reml_force_converge || einfo->_reml_no_converge) && prev_lgL > lgL){
            varcmp = prev_varcmp;
            calcu_Hi(einfo, einfo->_P, Hi,  _A);
            Hi = 2 * Hi;
            break;
        }
        
        // convergence
        dlogL = lgL - prev_lgL;
        if ((varcmp - prev_varcmp).squaredNorm() / varcmp.squaredNorm() < 1e-8 && (fabs(dlogL) < 1e-4 || (fabs(dlogL) < 1e-2 && dlogL < 0)))
        {
            converged_flag = true;
            for(int k=0;k<varcmp.size();k++)
                if(varcmp(k)<1e-6) remlstatus = -3;
            double vp=0;
            for(int k=0;k<varcmp.size();k++) vp+=varcmp(k);
            if(abs(vp)<1e-6) remlstatus = -3;
            else {
                for(int k=0;k<(varcmp.size()-1);k++) // the last one is v(e)
                {
                    if(varcmp(k)/vp > 0.999999) remlstatus = -3;
                }
            }
             if (constrain_num > 0) remlstatus = -5;
            if ( einfo->_reml_mtd == 2) {
                 calcu_Hi(einfo, einfo->_P, Hi, _A);
                Hi = 2 * Hi;
            } // for calculation of SE
            break;
        }
        prev_prev_varcmp = prev_varcmp;
        prev_varcmp = varcmp;
        prev_lgL = lgL;
    }
    
    if(einfo->_reml_fixed_var) {
        LOGPRINTF("Warning: the model is evaluated at fixed variance components. The likelihood might not be maximised.\n");
    }
    else {
        if(converged_flag) {
            LOGPRINTF( "Log-likelihood ratio converged.\n" );
        }
        else {
            if(einfo->_reml_force_converge || einfo->_reml_no_converge) {
                LOGPRINTF( "Warning: Log-likelihood not converged. Results are not reliable.\n" );
            }
            else if(iter == einfo->_reml_max_iter){
                if (einfo->_reml_max_iter > 1) {
                    LOGPRINTF("Error: Log-likelihood not converged (stop after %d iteractions). \nYou can specify the option --reml-maxit to allow for more iterations.\n",einfo->_reml_max_iter);
                    remlstatus = -4;
                }
            }
        }
    }
    return lgL;
}


void calcu_Vp(eInfo* einfo, double &Vp, double &Vp2, double &VarVp, double &VarVp2, VectorXd &varcmp, MatrixXd &Hi) {
    int i = 0, j = 0;
    Vp = 0.0;
    VarVp = 0.0;
    Vp2 = 0.0;
    VarVp2 = 0.0;
    if (einfo->_bivar_reml) {
        for (i = 0; i < einfo->_bivar_pos[0].size(); i++) {
            Vp += varcmp[einfo->_bivar_pos[0][i]];
            for (j = 0; j < einfo->_bivar_pos[0].size(); j++) VarVp += Hi(einfo->_bivar_pos[0][i], einfo->_bivar_pos[0][j]);
        }
        for (i = 0; i < einfo->_bivar_pos[1].size(); i++) {
            Vp2 += varcmp[einfo->_bivar_pos[1][i]];
            for (j = 0; j < einfo->_bivar_pos[1].size(); j++) VarVp2 += Hi(einfo->_bivar_pos[1][i], einfo->_bivar_pos[1][j]);
        }
        return;
    }
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        Vp += varcmp[i];
        for (j = 0; j < einfo->_r_indx.size(); j++) VarVp += Hi(i, j);
    }
}


void calcu_hsq(eInfo* einfo, int i, double Vp, double Vp2, double VarVp, double VarVp2, double &hsq, double &var_hsq, VectorXd &varcmp, MatrixXd &Hi) {
    int j = 0;
    double V1 = varcmp[i], VarV1 = Hi(i, i), Cov12 = 0.0;
    
    if (einfo->_bivar_reml) {
        vector<int>::iterator iter;
        iter = find(einfo->_bivar_pos[0].begin(), einfo->_bivar_pos[0].end(), i);
        if (iter != einfo->_bivar_pos[0].end()) {
            for (j = 0; j < einfo->_bivar_pos[0].size(); j++) {
                Cov12 += Hi(*iter, einfo->_bivar_pos[0][j]);
            }
            hsq = V1 / Vp;
            var_hsq = (V1 / Vp)*(V1 / Vp)*(VarV1 / (V1 * V1) + VarVp / (Vp * Vp)-(2 * Cov12) / (V1 * Vp));
            return;
        }
        iter = find(einfo->_bivar_pos[1].begin(), einfo->_bivar_pos[1].end(), i);
        if (iter != einfo->_bivar_pos[1].end()) {
            for (j = 0; j < einfo->_bivar_pos[1].size(); j++) {
                Cov12 += Hi(*iter, einfo->_bivar_pos[1][j]);
            }
            hsq = V1 / Vp2;
            var_hsq = (V1 / Vp2)*(V1 / Vp2)*(VarV1 / (V1 * V1) + VarVp2 / (Vp2 * Vp2)-(2 * Cov12) / (V1 * Vp2));
            return;
        }
        hsq = var_hsq = -2;
        return;
    }
    
    for (j = 0; j < einfo->_r_indx.size(); j++) {
        Cov12 += Hi(i, j);
    }
    hsq = V1 / Vp;
    var_hsq = (V1 / Vp)*(V1 / Vp)*(VarV1 / (V1 * V1) + VarVp / (Vp * Vp)-(2 * Cov12) / (V1 * Vp));
}


double lgL_reduce_mdl(eInfo* einfo,bool no_constrain,int _X_c,MatrixXd &_Vi,vector<MatrixXd> &_A, MatrixXd &_X,VectorXd &_y) {
    if (einfo->_r_indx.size() - 1 == 0) return 0;
    int _n=(int)einfo->_eii_include.size();
    bool multi_comp = (einfo->_r_indx.size() - einfo->_r_indx_drop.size() > 1);
    LOGPRINTF( "\nCalculating the logLikelihood for the reduced model ...\n(variance component %s",(multi_comp ? "s " : " "));
    for (int i = 0; i < einfo->_r_indx.size() - 1; i++) {
        if (find(einfo->_r_indx_drop.begin(), einfo->_r_indx_drop.end(), einfo->_r_indx[i]) == einfo->_r_indx_drop.end()) cout << einfo->_r_indx[i] + 1 << " ";
    }
    LOGPRINTF( "%s dropped from the model)\n",(multi_comp ? "are" : "is") );
    vector<int> vi_buf(einfo->_r_indx);
    einfo->_r_indx = einfo->_r_indx_drop;
    MatrixXd Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(einfo->_r_indx.size(), einfo->_r_indx.size());
    VectorXd Py(_n);
    VectorXd varcmp;
    vector<double> reml_priors_var, reml_priors;
    init_varcomp(einfo, reml_priors_var, reml_priors, varcmp );
    double lgL = reml_iteration(einfo, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, false, no_constrain,false, _Vi,  _A,  _X,_y );
    einfo->_r_indx = vi_buf;
    return lgL;
}


double lgL_fix_rg(eInfo* einfo, VectorXd &prev_varcmp, bool no_constrain,int _X_c,MatrixXd &_Vi,vector<MatrixXd> &_A, MatrixXd &_X,VectorXd &_y) {
    int i = 0, j = 0;
    int _n=(int)einfo->_eii_include.size();
    if (einfo->_fixed_rg_val.size() > einfo->_bivar_pos[0].size() - 1) {
        vector<double> rg_val_buf(einfo->_fixed_rg_val);
        einfo->_fixed_rg_val.clear();
        for (i = 0; i < einfo->_bivar_pos[0].size() - 1; i++) einfo->_fixed_rg_val.push_back(rg_val_buf[i]);
    }
    
    LOGPRINTF("\nCalculating the logLikelihood for the model with the genetic correlation %s being fixed at ",(einfo->_fixed_rg_val.size() > 1 ? "s" : ""));
    for (int i = 0; i < einfo->_fixed_rg_val.size() - 1; i++) {
        LOGPRINTF("%f\t", einfo->_fixed_rg_val[i] );
    }
    LOGPRINTF("%f\n",einfo->_fixed_rg_val[einfo->_fixed_rg_val.size() - 1]);
    
    vector<int> r_indx_buf(einfo->_r_indx);
    einfo->_bivar_pos_prev = einfo->_bivar_pos;
    for (i = einfo->_fixed_rg_val.size() - 1; i >= 0; i--) {
        for (j = i + 1; j < einfo->_bivar_pos_prev[0].size(); j++) einfo->_bivar_pos[0][j]--;
        for (j = i + 1; j < einfo->_bivar_pos_prev[1].size(); j++) einfo->_bivar_pos[1][j]--;
        for (j = i + 1; j < einfo->_bivar_pos_prev[2].size(); j++) einfo->_bivar_pos[2][j]--;
        einfo->_r_indx.erase(einfo->_r_indx.begin() + einfo->_bivar_pos[2][i]);
        einfo->_bivar_pos[2].erase(einfo->_bivar_pos[2].begin() + i);
    }
    
    einfo->_Asp_prev = einfo->_Asp;
    MatrixXd Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(einfo->_r_indx.size(), einfo->_r_indx.size());
    VectorXd Py(_n);
    VectorXd varcmp(einfo->_r_indx.size());
    for (i = 0; i < einfo->_r_indx.size(); i++) {
        varcmp[i] = fabs(prev_varcmp[einfo->_r_indx[i]]);
        if (varcmp[i] < 1.0e-30) varcmp[i] = 0.1;
    }
    double lgL = reml_iteration(einfo, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, false, no_constrain, true,_Vi,_A, _X,_y);
    einfo->_r_indx = r_indx_buf;
    einfo->_bivar_pos = einfo->_bivar_pos_prev;
    
    return lgL;
}


void eigenVector2Vector(VectorXd &x, vector<double> &y) {
    y.resize(x.size());
    for (int i = 0; i < x.size(); i++) y[i] = x[i];
}


void calcu_sum_hsq(eInfo* einfo, double Vp, double VarVp, double &sum_hsq, double &var_sum_hsq, VectorXd &varcmp, MatrixXd &Hi) {
    int i = 0, j = 0;
    double V1 = 0.0, VarV1 = 0.0, Cov12 = 0.0;
    for(i = 0; i < einfo->_r_indx.size()-1; i++) {
        V1 += varcmp[i];
        for(j = 0; j < einfo->_r_indx.size()-1; j++) VarV1 += Hi(i, j);
        for(j = 0; j < einfo->_r_indx.size(); j++) Cov12 += Hi(i, j);
    }
    sum_hsq = V1/Vp;
    var_sum_hsq = (V1/Vp)*(V1/Vp)*(VarV1/(V1*V1)+VarVp/(Vp*Vp)-(2*Cov12)/(V1*Vp));
}


double transform_hsq_L(double P, double K, double hsq) {
    double t = qnorm(1.0 - K);
    double z = dnorm(t);
    double C = (K * (1 - K) / (z * z))*(K * (1 - K) / (P * (1 - P)));
    return (hsq * C);
}


void calcu_rg(VectorXd &varcmp, MatrixXd &Hi, VectorXd &rg, VectorXd &rg_var, vector<string> &rg_name,vector< vector<int> > &_bivar_pos) {
    int i = 0, j = 0;
    double V1 = 0, V2 = 0, C = 0, VarV1 = 0, VarV2 = 0, VarC = 0.0, CovV1V2 = 0.0, CovV1C = 0.0, CovV2C = 0.0;
    
    rg = VectorXd::Zero(_bivar_pos[0].size() - 1);
    rg_var = rg;
    for (i = 0; i < _bivar_pos[0].size() - 1; i++) {
        V1 = varcmp[_bivar_pos[0][i]];
        V2 = varcmp[_bivar_pos[1][i]];
        C = varcmp[_bivar_pos[2][i]];
        VarV1 = Hi(_bivar_pos[0][i], _bivar_pos[0][i]);
        VarV2 = Hi(_bivar_pos[1][i], _bivar_pos[1][i]);
        VarC = Hi(_bivar_pos[2][i], _bivar_pos[2][i]);
        CovV1V2 = Hi(_bivar_pos[0][i], _bivar_pos[1][i]);
        CovV1C = Hi(_bivar_pos[0][i], _bivar_pos[2][i]);
        CovV2C = Hi(_bivar_pos[1][i], _bivar_pos[2][i]);
        if (V1 * V2 > 0) {
            rg[i] = sqrt(V1 * V2);
            if (rg[i] > 0) rg[i] = C / rg[i];
            rg_var[i] = rg[i] * rg[i]*(VarV1 / (4 * V1 * V1) + VarV2 / (4 * V2 * V2) + VarC / (C * C) + CovV1V2 / (2 * V1 * V2) - CovV1C / (V1 * C) - CovV2C / (V2 * C));
        }
        
        if (_bivar_pos[0].size() == 2) rg_name.push_back("rG");
        else {
            stringstream strstrm;
            strstrm << "rG" << i + 1;
            rg_name.push_back(strstrm.str());
        }
    }
}


void reml( eInfo* einfo, bool pred_rand_eff, bool est_fix_eff, vector<double> &reml_priors, vector<double> &reml_priors_var, double prevalence, double prevalence2, bool no_constrain, bool no_lrt, bool mlmassoc, int _X_c,MatrixXd &_X, VectorXd &_y,vector<MatrixXd> &_A, MatrixXd &_Vi, string _out)
{
    int i = 0, j = 0;
    int _n=(int)einfo->_eii_include.size();
    // Initialize variance component
    // 0: AI; 1: Fisher-scoring; 2: EM
    stringstream errmsg;
    VectorXd y_tmp = _y.array() - _y.mean();
    if (!einfo->_bivar_reml) {
        einfo->_y_Ssq = y_tmp.squaredNorm() / (_n - 1.0);
        if (!(fabs(einfo->_y_Ssq) < 1e30)) {LOGPRINTF ("Error: the phenotypic variance is infinite. Please check the missing data in your phenotype file. Missing values should be represented by \"NA\" or \"-9\"."); TERMINATE();}
    }
    bool reml_priors_flag = !reml_priors.empty(), reml_priors_var_flag = !reml_priors_var.empty();
    if (reml_priors_flag && reml_priors.size() < einfo->_r_indx.size() - 1) {
        LOGPRINTF("Error: in option --reml-priors. There are %ld variance components. At least %ld prior values should be specified.\n", einfo->_r_indx.size(),einfo->_r_indx.size() - 1 );
        TERMINATE();
    }
    if (reml_priors_var_flag && reml_priors_var.size() < einfo->_r_indx.size() - 1) {
         LOGPRINTF("Error: in option --reml-priors-var. There are %ld variance components. At least %ld prior values should be specified.\n", einfo->_r_indx.size(),einfo->_r_indx.size() - 1 );
        TERMINATE();
    }
    
    LOGPRINTF("\nPerforming %s REML analysis ... (Note: may take hours depending on sample size).\n",(einfo->_bivar_reml ? "bivariate" : ""));
    if (_n < 10) {LOGPRINTF ("Error: sample size is too small.");exit(EXIT_FAILURE);}
    LOGPRINTF("%d observations, %d fixed effect(s), and %ld variance component(s)(including residual variance).\n",_n,_X_c,einfo->_r_indx.size());
    MatrixXd Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(einfo->_r_indx.size(), einfo->_r_indx.size());
    VectorXd Py(_n), varcmp;
    init_varcomp(einfo,reml_priors_var, reml_priors, varcmp);
    double lgL = reml_iteration(einfo, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, reml_priors_var_flag | reml_priors_flag, no_constrain,false,_Vi, _A,   _X,_y );
    if(remlstatus!=0 && remlstatus!=-3 && remlstatus!=-5) return;
    
    MatrixXd u;
    if (pred_rand_eff) {
        u.resize(_n, einfo->_r_indx.size());
        for (i = 0; i < einfo->_r_indx.size(); i++) {
            if (einfo->_bivar_reml || einfo->_within_family)(u.col(i)) = (((einfo->_Asp[einfo->_r_indx[i]]) * Py) * varcmp[i]);
            else (u.col(i)) = (((_A[einfo->_r_indx[i]]) * Py) * varcmp[i]);
        }
    }
    if (est_fix_eff) einfo->_b = Xt_Vi_X_i * (Vi_X.transpose() * _y);
    
    // calculate Hsq and SE
    double Vp = 0.0, Vp2 = 0.0, VarVp = 0.0, VarVp2 = 0.0, Vp_f = 0.0, VarVp_f = 0.0;
    vector<double> Hsq(einfo->_r_indx.size() - 1), VarHsq(einfo->_r_indx.size() - 1);
    calcu_Vp(einfo, Vp, Vp2, VarVp, VarVp2, varcmp, Hi);
    for (i = 0; i < Hsq.size(); i++) calcu_hsq(einfo, i, Vp, Vp2, VarVp, VarVp2, Hsq[i], VarHsq[i], varcmp, Hi);
   
    // calculate the logL for a reduce model
    double lgL_rdu_mdl = 0.0, LRT = 0.0;
    if (!no_lrt) {
        lgL_rdu_mdl = lgL_reduce_mdl(einfo, no_constrain, _X_c,_Vi, _A,  _X,_y);
        if(remlstatus!=0 && remlstatus!=-3) return;
        LRT = 2.0 * (lgL - lgL_rdu_mdl);
        if (LRT < 0.0) LRT = 0.0;
    }
    
    // calcuate the logL given a rG in a bivariate analysis
    double lgL_fixed_rg = 0.0;
    if (einfo->_bivar_reml && !einfo->_fixed_rg_val.empty()) {
        lgL_fixed_rg = lgL_fix_rg(einfo, varcmp, no_constrain, _X_c, _Vi, _A,  _X,_y);
        if(remlstatus!=0 && remlstatus!=-3) return;
        LRT = 2.0 * (lgL - lgL_fixed_rg);
        if (LRT < 0.0) LRT = 0.0;
    }
    
    if (mlmassoc) {
        if(remlstatus==0) eigenVector2Vector(varcmp, einfo->_varcmp);
        return;
    }
    // output results
    double sum_hsq = 0.0, var_sum_hsq = 0.0;
    if (!einfo->_bivar_reml && einfo->_r_indx.size() > 2) calcu_sum_hsq(einfo, Vp, VarVp, sum_hsq, var_sum_hsq, varcmp, Hi);
    LOGPRINTF("\nSummary result of REML analysis:\n");
    cout << "Source\tVariance\tSE" << setiosflags(ios::fixed) << setprecision(6) << endl;
    for (i = 0; i < einfo->_r_indx.size(); i++) cout << einfo->_var_name[i] << "\t" << varcmp[i] << "\t" << sqrt(Hi(i, i)) << endl;
    if (einfo->_bivar_reml) {
        cout << "Vp_tr1\t" << Vp << "\t" << sqrt(VarVp) << endl;
        cout << "Vp_tr2\t" << Vp2 << "\t" << sqrt(VarVp2) << endl;
        for (i = 0, j = 0; i < einfo->_bivar_pos[0].size() - 1; i++, j += 2) {
            cout << einfo->_hsq_name[j] << "\t" << Hsq[einfo->_bivar_pos[0][i]] << "\t" << sqrt(VarHsq[einfo->_bivar_pos[0][i]]) << endl;
            cout << einfo->_hsq_name[j + 1] << "\t" << Hsq[einfo->_bivar_pos[1][i]] << "\t" << sqrt(VarHsq[einfo->_bivar_pos[1][i]]) << endl;
        }
    } else {
        cout << "Vp\t" << Vp << "\t" << sqrt(VarVp) << endl;
        for (i = 0; i < Hsq.size(); i++) cout << einfo->_hsq_name[i] << "\t" << Hsq[i] << "\t" << sqrt(VarHsq[i]) << endl;
        if (einfo->_r_indx.size() > 2) cout << "\nSum of V(G)/Vp\t" << sum_hsq << "\t" << sqrt(var_sum_hsq) << endl;
    }
    if ((einfo->_flag_CC && prevalence>-1) || (einfo->_flag_CC2 && prevalence2>-1)) {
        cout << "The estimate of variance explained on the observed scale is transformed to that on the underlying scale:" << endl;
        if (einfo->_bivar_reml) {
            if (einfo->_flag_CC) cout << "Proportion of cases in the sample = " << einfo->_ncase << " for trait #1; User-specified disease prevalence = " << prevalence << " for trait #1" << endl;
            if (einfo->_flag_CC2) cout << "Proportion of cases in the sample = " << einfo->_ncase2 << " for trait #2; User-specified disease prevalence = " << prevalence2 << " for trait #2" << endl;
            for (i = 0, j = 0; i < einfo->_bivar_pos[0].size() - 1; i++, j += 2) {
                if (einfo->_flag_CC) cout << einfo->_hsq_name[j] << "_L\t" << transform_hsq_L(einfo->_ncase, prevalence, Hsq[einfo->_bivar_pos[0][i]]) << "\t" << transform_hsq_L(einfo->_ncase, prevalence, sqrt(VarHsq[einfo->_bivar_pos[0][i]])) << endl;
                if (einfo->_flag_CC2) cout << einfo->_hsq_name[j + 1] << "_L\t" << transform_hsq_L(einfo->_ncase2, prevalence2, Hsq[einfo->_bivar_pos[1][i]]) << "\t" << transform_hsq_L(einfo->_ncase2, prevalence2, sqrt(VarHsq[einfo->_bivar_pos[1][i]])) << endl;
            }
        } else {
            cout << "(Proportion of cases in the sample = " << einfo->_ncase << "; User-specified disease prevalence = " << prevalence << ")" << endl;
            for (i = 0; i < Hsq.size(); i++) cout << einfo->_hsq_name[i] << "_L\t" << transform_hsq_L(einfo->_ncase, prevalence, Hsq[i]) << "\t" << transform_hsq_L(einfo->_ncase, prevalence, sqrt(VarHsq[i])) << endl;
            if (einfo->_r_indx.size() > 2)  cout << "\nSum of V(G)_L/Vp\t" << transform_hsq_L(einfo->_ncase, prevalence, sum_hsq) << "\t" << transform_hsq_L(einfo->_ncase, prevalence, sqrt(var_sum_hsq)) << endl;
        }
    }
    // output genetic correlation
    VectorXd rg, rg_var;
    vector<string> rg_name;
    if (einfo->_bivar_reml) {
        calcu_rg(varcmp, Hi, rg, rg_var, rg_name,einfo->_bivar_pos);
        for (i = 0; i < rg_name.size(); i++) {
            cout << rg_name[i] << "\t" << rg[i] << "\t" << sqrt(rg_var[i]) << endl;
        }
    }
    if(!einfo->_reml_force_converge || !einfo->_reml_AI_not_invertible){
        cout << "\nSampling variance/covariance of the estimates of variance components:" << endl;
        for (i = 0; i < einfo->_r_indx.size(); i++) {
            for (j = 0; j < einfo->_r_indx.size(); j++) printf("%e\t",Hi(i, j));
            cout << endl;
        }
    }
    if (est_fix_eff) {
        LOGPRINTF("Estimate of fixed effect %s:\n",(_X_c > 1 ? "s" : ""))
        LOGPRINTF("\nSource\tEstimate\tSE\n");
        for (i = 0; i < _X_c; i++) {
            if (i == 0){ LOGPRINTF("mean\t");}
            else { LOGPRINTF( "X_%d\t", i+1);}
            LOGPRINTF("%f\t%f\n", einfo->_b[i] , sqrt(Xt_Vi_X_i(i, i)));
        }
    }
    
    // save summary result into a file
    string reml_rst_file = _out + ".hsq";
    ofstream o_reml(reml_rst_file.c_str());
    if (!o_reml) throw ("Error: can not open the file [" + reml_rst_file + "] to write.");
    o_reml << "Source\tVariance\tSE" << setiosflags(ios::fixed) << setprecision(6) << endl;
    for (i = 0; i < einfo->_r_indx.size(); i++) o_reml << einfo->_var_name[i] << "\t" << varcmp[i] << "\t" << sqrt(Hi(i, i)) << endl;
    if (einfo->_bivar_reml) {
        o_reml << "Vp_tr1\t" << Vp << "\t" << sqrt(VarVp) << endl;
        o_reml << "Vp_tr2\t" << Vp2 << "\t" << sqrt(VarVp2) << endl;
        for (i = 0, j = 0; i < einfo->_bivar_pos[0].size() - 1; i++, j += 2) {
            o_reml << einfo->_hsq_name[j] << "\t" << Hsq[einfo->_bivar_pos[0][i]] << "\t" << sqrt(VarHsq[einfo->_bivar_pos[0][i]]) << endl;
            o_reml << einfo->_hsq_name[j + 1] << "\t" << Hsq[einfo->_bivar_pos[1][i]] << "\t" << sqrt(VarHsq[einfo->_bivar_pos[1][i]]) << endl;
        }
    } else {
        o_reml << "Vp\t" << Vp << "\t" << sqrt(VarVp) << endl;
        for (i = 0; i < Hsq.size(); i++) o_reml << einfo->_hsq_name[i] << "\t" << Hsq[i] << "\t" << sqrt(VarHsq[i]) << endl;
        if (einfo->_r_indx.size() > 2) o_reml << "\nSum of V(G)/Vp\t" << sum_hsq << "\t" << sqrt(var_sum_hsq) << endl;
    }
    if (einfo->_flag_CC && prevalence>-1) {
        o_reml << "The estimate of variance explained on the observed scale is transformed to that on the underlying scale:" << endl;
        if (einfo->_bivar_reml) {
            o_reml << "(Proportion of cases in the sample = " << einfo->_ncase << "; User-specified disease prevalence = " << prevalence << " for disease 1 and = " << prevalence2 << " for disease 2)" << endl;
            for (i = 0, j = 0; i < einfo->_bivar_pos[0].size() - 1; i++, j += 2) {
                o_reml << einfo->_hsq_name[j] << "_L\t" << transform_hsq_L(einfo->_ncase, prevalence, Hsq[einfo->_bivar_pos[0][i]]) << "\t" << transform_hsq_L(einfo->_ncase, prevalence, sqrt(VarHsq[einfo->_bivar_pos[0][i]])) << endl;
                o_reml << einfo->_hsq_name[j + 1] << "_L\t" << transform_hsq_L(einfo->_ncase2, prevalence2, Hsq[einfo->_bivar_pos[1][i]]) << "\t" << transform_hsq_L(einfo->_ncase2, prevalence2, sqrt(VarHsq[einfo->_bivar_pos[1][i]])) << endl;
            }
        } else {
            o_reml << "(Proportion of cases in the sample = " << einfo->_ncase << "; User-specified disease prevalence = " << prevalence << ")" << endl;
            for (i = 0; i < Hsq.size(); i++) o_reml << einfo->_hsq_name[i] << "_L\t" << transform_hsq_L(einfo->_ncase, prevalence, Hsq[i]) << "\t" << transform_hsq_L(einfo->_ncase, prevalence, sqrt(VarHsq[i])) << endl;
            if (einfo->_r_indx.size() > 2)  o_reml << "\nSum of V(G)_L/Vp\t" << transform_hsq_L(einfo->_ncase, prevalence, sum_hsq) << "\t" << transform_hsq_L(einfo->_ncase, prevalence, sqrt(var_sum_hsq)) << endl;
        }
    }
    if (einfo->_bivar_reml) {
        for (i = 0; i < rg_name.size(); i++) o_reml << rg_name[i] << "\t" << rg[i] << "\t" << sqrt(rg_var[i]) << endl;
    }
    o_reml << "logL\t" << setprecision(3) << lgL << endl;
    if (!no_lrt && einfo->_r_indx.size() - 1 > 0) {
        o_reml << "logL0\t" << setprecision(3) << lgL_rdu_mdl << endl;
        o_reml << "LRT\t" << setprecision(3) << LRT << endl;
        o_reml << "df\t" << setprecision(1) << einfo->_r_indx.size() - einfo->_r_indx_drop.size() << endl;
        o_reml << "Pval\t" << dtos( 0.5 * chi_prob(einfo->_r_indx.size() - einfo->_r_indx_drop.size(), LRT)) << endl;
    }
    if (einfo->_bivar_reml && !einfo->_fixed_rg_val.empty()) {
        o_reml << "logL0\t" << setprecision(3) << lgL_fixed_rg << " (when rG fixed at ";
        for (i = 0; i < einfo->_fixed_rg_val.size() - 1; i++) o_reml << einfo->_fixed_rg_val[i] << "\t";
        o_reml << einfo->_fixed_rg_val[einfo->_fixed_rg_val.size() - 1] << ")" << endl;
        o_reml << "LRT\t" << setprecision(3) << LRT << endl;
        o_reml << "df\t" << setprecision(1) << einfo->_fixed_rg_val.size() << endl;
        o_reml << "Pval\t" << dtos( 0.5 * chi_prob(einfo->_fixed_rg_val.size(), LRT) ) << " (one-tailed test)" << endl;
    }
    o_reml << "n\t" << _n << endl;
    if (est_fix_eff) {
        o_reml << "\nFix_eff\tSE" << endl;
        for (i = 0; i < _X_c; i++) o_reml << setprecision(6) << einfo->_b[i] << "\t" << sqrt(Xt_Vi_X_i(i, i)) << endl;
        o_reml.close();
    }
    cout << "\nSummary result of REML analysis has been saved in the file [" + reml_rst_file + "]." << endl;
    
    // save random effect to a file
    if (pred_rand_eff) {
        string rand_eff_file = _out + ".indi.blp";
        ofstream o_rand_eff(rand_eff_file.c_str());
        for (i = 0; i < einfo->_eii_include.size(); i++)
        {
            o_rand_eff<<einfo->_eii_fid[einfo->_eii_include[i]]<<'\t'<<einfo->_eii_iid[einfo->_eii_include[i]]<<'\t';
            for (j = 0; j < einfo->_r_indx.size(); j++) o_rand_eff << setprecision(6) << Py[i] * varcmp[j] << "\t" << u(i, j) << "\t";
            o_rand_eff << endl;
        }
        cout << "\nBLUP of the genetic effects for " << einfo->_eii_include.size() << " individuals has been saved in the file [" + rand_eff_file + "]." << endl;
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

int construct_X(eInfo* einfo, vector<MatrixXd> &E_float, MatrixXd &qE_float, MatrixXd &_X) {
    
    int n=(int)einfo->_eii_include.size();
    
    int i = 0, j = 0;
    map<string, int>::iterator iter;
    stringstream errmsg;
    
   int  _X_c = 1;
    // quantitative covariates
    MatrixXd X_q;
    if (einfo->_eii_qcov_num>0) {
        X_q.resize(n, einfo->_eii_qcov_num);
        for (i = 0; i < n; i++) {
            for (j = 0; j < einfo->_eii_qcov_num; j++) X_q(i, j) = einfo->_eii_qcov[j*einfo->_eii_num+einfo->_eii_include[i]];
        }
        LOGPRINTF("%d quantitative variable(s) included as covariate(s).\n",einfo->_eii_qcov_num);
        _X_c += einfo->_eii_qcov_num;
    }

    // discrete covariates
    vector<MatrixXd> X_d;
    if (einfo->_eii_cov_num>0) {
        vector< vector<string> > covar_tmp(einfo->_eii_cov_num);
        for (i = 0; i < einfo->_eii_cov_num; i++) covar_tmp[i].resize(n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < einfo->_eii_cov_num; j++) covar_tmp[j][i] = einfo->_eii_cov[j*einfo->_eii_num+einfo->_eii_include[i]];
        }
        LOGPRINTF("%d discrete variable(s) included as covariate(s).\n",einfo->_eii_cov_num);
        X_d.resize(einfo->_eii_cov_num);
        for (i = 0; i < einfo->_eii_cov_num; i++) {
            stringstream errmsg;
            errmsg << "Error: too many classes for the " << i + 1 << "th discrete variable. \nPlease use the --qcovar if it is a quantitative covariate.";
            string errmsg1 = errmsg.str();
            errmsg.str("");
            errmsg << "Error: the " << i + 1 << "th discrete variable has only one class.";
            string errmsg2 = errmsg.str();
            coeff_mat(covar_tmp[i], X_d[i], errmsg1, errmsg2);
            _X_c += (X_d[i]).cols() - 1;
        }
    }
    // E factor
    _X_c += qE_float.cols();
    for (i = 0; i < E_float.size(); i++) _X_c += (E_float[i]).cols() - 1;
    
    // Construct _X
    int col = 0;
    _X.resize(n, _X_c);
    _X.block(0, col, n, 1) = MatrixXd::Ones(n, 1);
    col++;
    if (einfo->_eii_qcov_num>0) {
        _X.block(0, col, n, X_q.cols()) = X_q;
        col += X_q.cols();
    }
    for (i = 0; i < X_d.size(); i++) {
        _X.block(0, col, n, (X_d[i]).cols() - 1) = (X_d[i]).block(0, 1, n, (X_d[i]).cols() - 1);
        col += (X_d[i]).cols() - 1;
    }
    if (qE_float.cols() > 0) {
        _X.block(0, col, n, qE_float.cols()) = qE_float;
        col += qE_float.cols();
    }
    for (i = 0; i < E_float.size(); i++) {
        _X.block(0, col, n, (E_float[i]).cols() - 1) = (E_float[i]).block(0, 1, n, (E_float[i]).cols() - 1);
        col += (E_float[i]).cols() - 1;
    }
    return _X_c;
}

void detect_family(eInfo* einfo,  vector<MatrixXd> &_A)
    {
        cout<<"Detecting sub-matrix for each family from the ORM ..."<<endl;
        int _n=(int)einfo->_eii_include.size();
        int i=0, j=0, k=0, l=0, prev_pnt=0;
        double d_buf1=0.0, d_buf2=0.0;
        einfo->_fam_brk_pnt.clear();
        for(i=0; i<_n-1; i++){
            d_buf1=_A[0].row(i).tail(_n-i-1).sum();
            d_buf2=_A[0].col(i).tail(_n-i-1).sum();
            if(FloatEqual(d_buf1, 0.0) && FloatEqual(d_buf2, 0.0)) einfo->_fam_brk_pnt.push_back(i);
        }
        
        einfo->_Asp.resize(einfo->_r_indx.size());
        for(i=0; i<einfo->_r_indx.size(); i++) (einfo->_Asp[i]).resize(_n, _n);
        
        int pos=0;
        for(l=0; l<einfo->_r_indx.size()-1; l++){
            pos=einfo->_r_indx[l];
            prev_pnt=0;
            for(k=0; k<einfo->_fam_brk_pnt.size(); k++){
                for(j=prev_pnt; j<=einfo->_fam_brk_pnt[k]; j++){
                    (einfo->_Asp[pos]).startVec(j);
                    for(i=prev_pnt; i<=einfo->_fam_brk_pnt[k]; i++) (einfo->_Asp[pos]).insertBack(i,j)=(_A[pos])(i,j);
                }
                prev_pnt=einfo->_fam_brk_pnt[k]+1;
            }
            for(j=prev_pnt; j<_n; j++){
                (einfo->_Asp[pos]).startVec(j);
                for(i=prev_pnt; i<_n; i++) (einfo->_Asp[pos]).insertBack(i,j)=(_A[pos])(i,j);
            }
            (einfo->_Asp[pos]).finalize();
        }
        pos=einfo->_r_indx[einfo->_r_indx.size()-1];
        for(i=0; i<_n; i++){
            (einfo->_Asp[pos]).startVec(i);
            (einfo->_Asp[pos]).insertBack(i,i)=1.0;
        }
        (einfo->_Asp[pos]).finalize();
        cout<<"There are "<<einfo->_fam_brk_pnt.size()+1<<" sub-matrices detected."<<endl;
        
        // release momery
        _A.clear();
    }

void mlma_calcu_stat_covar(VectorXd &_Y, double *predictor, unsigned long n, unsigned long m, int _X_c,  MatrixXd &_Vi,  MatrixXd &_X, VectorXd &beta, VectorXd &se, VectorXd &pval)
    {
        
        unsigned long col_num=_X_c+1;
        double chisq=0.0;
        
        MatrixXd X(n,col_num);
        MatrixXd Vi_X(n,col_num);
        MatrixXd Xt_Vi_X(col_num,col_num);
        VectorXd Xt_Vi_y(col_num);
        VectorXd b_vec(col_num);
        
        X.block(0, 0, n, _X_c)=_X;
        
        beta.resize(m);
        se=VectorXd::Zero(m);
        pval=VectorXd::Constant(m,2);
        LOGPRINTF("\nRunning association tests for %ld probes ...\n",m);
        
        
        for(int i = 0; i < m; i++ ){
            
            for(int j = 0; j < n; j++) X(j,_X_c) = predictor[j*m+i];
            Vi_X=_Vi*_X;
            Xt_Vi_X=X.transpose()*Vi_X;
            double logdt=0.0;
            if(!comput_inverse_logdet_LU( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
            Xt_Vi_y=Vi_X*_Y;
            b_vec=Xt_Vi_X*Xt_Vi_y;
            
            se[i]=Xt_Vi_X(_X_c,_X_c);
            beta[i]=b_vec[_X_c];
            if(se[i]>1.0e-30){
                se[i]=sqrt(se[i]);
                chisq=beta[i]/se[i];
                pval[i]=pchisq(chisq*chisq, 1);
            }
        }
        
    }

void mlma_calcu_stat_covar(VectorXd &_Y, eInfo* einfo, int _X_c,  MatrixXd &_Vi,  MatrixXd &_X, VectorXd &beta, VectorXd &se, VectorXd &pval)
    {
        
        unsigned long col_num=_X_c+1;
        double chisq=0.0;
        uint64_t n=einfo->_eii_include.size();
        uint64_t m=einfo->_epi_include.size();
        
        MatrixXd X(n,col_num);
        MatrixXd Vi_X(n,col_num);
        MatrixXd Xt_Vi_X(col_num,col_num);
        VectorXd Xt_Vi_y(col_num);
        VectorXd b_vec(col_num);

        X.block(0, 0, n, _X_c)=_X;
        
        beta.resize(m);
        se=VectorXd::Zero(m);
        pval=VectorXd::Constant(m,2);
        LOGPRINTF("\nRunning association tests for %llu probes ...\n",m);
     
       
        for(int i = 0; i < m; i++ ){
            double nonmiss=0.0;
            double mu=0.0;
            //calc the mean
            for(int j = 0; j < n; j++) {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            for(int j = 0; j < n; j++) {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    X(j,_X_c) = val-mu;
                } else {
                   X(j,_X_c)=0.0;
                }
            }
            Vi_X=_Vi*_X;
            Xt_Vi_X=X.transpose()*Vi_X;
            double logdt=0.0;
            if(!comput_inverse_logdet_LU( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
            Xt_Vi_y=Vi_X*_Y;
            b_vec=Xt_Vi_X*Xt_Vi_y;
          
            se[i]=Xt_Vi_X(_X_c,_X_c);
            beta[i]=b_vec[_X_c];
            if(se[i]>1.0e-30){
                se[i]=sqrt(se[i]);
                chisq=beta[i]/se[i];
                pval[i]=pchisq(chisq*chisq, 1);
            }
        }
    }
    
    
void mlma_calcu_stat(VectorXd &_Y, eInfo* einfo, MatrixXd &_Vi, VectorXd &beta, VectorXd &se, VectorXd &pval)
    {
  
        double Xt_Vi_X=0.0, chisq=0.0;
        uint64_t n=einfo->_eii_include.size();
        uint64_t m=einfo->_epi_include.size();
        VectorXd X(n);
        VectorXd Vi_X(n);
        
        beta.resize(m);
        se=VectorXd::Zero(m);
        pval=VectorXd::Constant(m,2);
        LOGPRINTF("\nRunning association tests for %llu probes ...\n",m);
        
        for(int i = 0; i < m; i++){
            double nonmiss=0.0;
            double mu=0.0;
            //calc the mean
            for(int j = 0; j < n; j++) {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            for(int j = 0; j < n; j++) {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    X(j) = val-mu;
                } else {
                    X(j)=0.0;
                }
            }
            Vi_X=_Vi*X;
            Xt_Vi_X=Vi_X.dot(X);
            se[i]=1.0/Xt_Vi_X;
            beta[i]=se[i]*(Vi_X.dot(_Y));
            if(se[i]>1.0e-30){
                se[i]=sqrt(se[i]);
                chisq=beta[i]/se[i];
                pval[i]=pchisq(chisq*chisq, 1);
            }
        }
    }
    void mlma_calcu_stat(VectorXd &_Y, double *predictor, unsigned long n, unsigned long m, MatrixXd &_Vi, VectorXd &beta, VectorXd &se, VectorXd &pval)
    {
        
        double Xt_Vi_X=0.0, chisq=0.0;
       
        VectorXd X(n);
        VectorXd Vi_X(n);
        
        beta.resize(m);
        se=VectorXd::Zero(m);
        pval=VectorXd::Constant(m,2);
        LOGPRINTF("\nRunning association tests for %ld probes ...\n",m);
        
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++) X(j) = predictor[j*m+i];
            Vi_X=_Vi*X;
            
            Xt_Vi_X=Vi_X.dot(X);
            se[i]=1.0/Xt_Vi_X;
            beta[i]=se[i]*(Vi_X.dot(_Y));
            if(se[i]>1.0e-30){
                se[i]=sqrt(se[i]);
                chisq=beta[i]/se[i];
                pval[i]=pchisq(chisq*chisq, 1);
            }
        }

    }
    
    void read_indi_blup(eInfo* einfo, string blup_indi_file) {
        vector< vector<string> > g_buf;
        ifstream i_indi_blup(blup_indi_file.c_str());
        if (!i_indi_blup) throw ("Error: can not open the file [" + blup_indi_file + "] to read.");
        string str_buf, id_buf;
        vector<string> id, vs_buf;
        int i = 0, j = 0, k = 0, col_num = 0;
        while (i_indi_blup) {
            i_indi_blup >> str_buf;
            if (i_indi_blup.eof()) break;
            id_buf = str_buf + ":";
            i_indi_blup >> str_buf;
            id_buf += str_buf;
            getline(i_indi_blup, str_buf);
            col_num = split_string(str_buf, vs_buf, " \t\n");
            if (col_num < 1) continue;
            id.push_back(id_buf);
            g_buf.push_back(vs_buf);
        }
        i_indi_blup.close();
        
        update_map_kp(id, einfo->_eii_map, einfo->_eii_include);
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for (i = 0; i < einfo->_eii_include.size(); i++) uni_id_map.insert(pair<string, int>(einfo->_eii_fid[einfo->_eii_include[i]] + ":" + einfo->_eii_iid[einfo->_eii_include[i]], i));
        einfo->_varcmp_Py.setZero(einfo->_eii_include.size(), col_num / 2);
        for (i = 0; i < id.size(); i++) {
            iter = uni_id_map.find(id[i]);
            if (iter == uni_id_map.end()) continue;
            for (j = 0, k = 0; j < col_num; j += 2, k++) einfo->_varcmp_Py(iter->second, k) = atof(g_buf[i][j].c_str());
        }
        cout << "BLUP solution to the total genetic effects for " << einfo->_eii_include.size() << " individuals have been read from [" + blup_indi_file + "]." << endl;
    }
    
    // blup estimate of probe effect
    void output_blup_snp(char* outFileName,  eInfo* einfo, MatrixXd &b_probe) {
        string o_b_snp_file = string(outFileName) + ".probe.blp";
        ofstream o_b_snp(o_b_snp_file.c_str());
        if (!o_b_snp) throw ("Error: can not open the file " + o_b_snp_file + " to write.");
        long col_num = b_probe.cols();
        cout << "Writing BLUP solutions of probe effects for " << einfo->_epi_include.size() << " probes to [" + o_b_snp_file + "]." << endl;
        for (int i = 0; i < einfo->_epi_include.size(); i++) {
            o_b_snp << einfo->_epi_prb[einfo->_epi_include[i]] << "\t" ;
            for (int j = 0; j < 1; j++) o_b_snp << b_probe(i, j) << "\t"; // not output the residual effect
            o_b_snp << endl;
        }
        o_b_snp.close();
        cout << "BLUP solutions of probe effects for " << einfo->_epi_include.size() << " probes have been saved in the file [" + o_b_snp_file + "]." << endl;
    }
    
    
    void blup_probe_geno( eInfo* einfo, char* outFileName, char* blup_indi_file) {
        
        read_indi_blup(einfo, blup_indi_file);
        
        if (einfo->_mu.empty() || einfo->_var.empty()) cal_var_mean(einfo, true, true);
        
        long col_num = einfo->_varcmp_Py.cols();
        double fcount = 0.0;
        
        // Calcuate A matrix
        cout << "Calculating the BLUP solution to probe effects ..." << endl;
        MatrixXd b_probe = MatrixXd::Zero(einfo->_epi_include.size(), col_num); // variance of each probe
        
        einfo->_var.resize(einfo->_epi_num);
        for (int j = 0; j < einfo->_epi_include.size(); j++) {
            if (fabs(einfo->_var[einfo->_epi_include[j]] ) < 1.0e-50) einfo->_var[einfo->_epi_include[j]]  = 0.0;
            else einfo->_var[einfo->_epi_include[j]]  = 1.0 / einfo->_var[einfo->_epi_include[j]] ;
        }
        
        for (int k = 0; k < einfo->_epi_include.size(); k++) {
            fcount = 0.0;
            for (int i = 0; i < einfo->_eii_include.size(); i++) {
                double x=einfo->_val[einfo->_epi_include[k]*einfo->_eii_num+einfo->_eii_include[i]];
                if(x<1e9) {
                    x = (x - einfo->_mu[einfo->_epi_include[k]]);
                    for (int j = 0; j < col_num; j++) b_probe(k, j) += x * einfo->_varcmp_Py(i, j);
                    fcount += 1.0;
                }
            }
            for (int j = 0; j < col_num; j++) b_probe(k, j) = (b_probe(k, j) * einfo->_var[einfo->_epi_include[k]] / fcount)*((double) einfo->_eii_include.size() / (double) einfo->_epi_include.size());
        }
        output_blup_snp(outFileName, einfo,b_probe);
    }


}

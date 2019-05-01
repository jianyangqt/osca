//
//  l3_mlma.h
//  osc
//
//  Created by Futao Zhang on 13/04/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l3_mlma__
#define __osc__l3_mlma__

#include "l2_efile.h"
#include "l0_stat.h"
namespace EFILE {
    void detect_family(eInfo* einfo, vector<MatrixXd> &_A);
    bool calcu_Vi(MatrixXd &Vi, VectorXd &prev_varcmp, double &logdet,MatrixXd &_U, VectorXd &_S);
    void reml( eInfo* einfo, bool pred_rand_eff, bool est_fix_eff, vector<double> &reml_priors, vector<double> &reml_priors_var, double prevalence, double prevalence2, bool no_constrain, bool no_lrt, bool mlmassoc, int _X_c,MatrixXd &_X, VectorXd &_y,vector<MatrixXd> &_A, MatrixXd &_Vi, string _out);
    void reml( bool pred_rand_eff, bool est_fix_eff, vector<double> &reml_priors, vector<double> &reml_priors_var, bool no_constrain, int _X_c,MatrixXd &_X, VectorXd &_y,vector<MatrixXd> &_A, MatrixXd &_U, VectorXd &_eval, MatrixXd &_Vi, int reml_mtd, int reml_max_iter,VectorXd &_b,VectorXd &_se);
    void mlma_calcu_stat_covar(VectorXd &_Y, double *predictor, unsigned long n, unsigned long m, int _X_c,  MatrixXd &_Vi,  MatrixXd &_X, VectorXd &beta, VectorXd &se, VectorXd &pval);
    void mlma_calcu_stat_covar(VectorXd &_Y, eInfo* einfo, int _X_c,  MatrixXd &_Vi,  MatrixXd &_X, VectorXd &beta, VectorXd &se, VectorXd &pval);
    void mlma_calcu_stat(VectorXd &_Y, eInfo* einfo, MatrixXd &_Vi, VectorXd &beta, VectorXd &se, VectorXd &pval);
    void mlma_calcu_stat(VectorXd &_Y, double *predictor, unsigned long n, unsigned long m, MatrixXd &_Vi, VectorXd &beta, VectorXd &se, VectorXd &pval);
    void blup_probe_geno( eInfo* einfo, char* outFileName, char* blup_indi_file);
    void eigenVector2Vector(VectorXd &x, vector<double> &y);
}
#endif /* defined(__osc__l3_mlma__) */

//
//  l3_gwas.hpp
//  osc
//
//  Created by Futao Zhang on 8/5/19.
//  Copyright © 2019 Futao Zhang. All rights reserved.
//

#ifndef l3_gwas_hpp
#define l3_gwas_hpp

#include "l2_reml.h"
#include "l2_bfile.h"

namespace BFILE{
    void moment_gwas(char* outFileName, char* bFileName,double maf,char* snplstName,char* snplst2exclde,int chr,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool nopreadj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_num, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_flag, bool approximate_stepwise,int grm_alg, double swthresh, double swfdr, bool swlogit,bool swforward,double swrsq, char* grm_file,bool grm_bin_flag);
    void mlm_cal_stat(VectorXd &_Y, bInfo* bdata, MatrixXd &_Vi, VectorXd &beta, VectorXd &se, VectorXd &pval);
    void mlm_cal_stat_covar(VectorXd &_Y,bInfo* bdata,  MatrixXd &_Vi,  MatrixXd &_COV, VectorXd &beta, VectorXd &se, VectorXd &pval);
}

#endif /* l3_gwas_hpp */

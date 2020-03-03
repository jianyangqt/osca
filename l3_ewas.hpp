//
//  l3_ewas.hpp
//  osc
//
//  Created by Futao Zhang on 16/04/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#ifndef l3_ewas_hpp
#define l3_ewas_hpp

#include "l3_efile.h"
namespace EFILE {
    void moment(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool nopreadj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_num, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_flag, bool approximate_stepwise,int erm_alg, double swthresh, double swfdr, bool swlogit, bool stepforwardonly, bool Baptiste, double sw_rsq, double thresh_pcc, char* orm_file,bool orm_bin_flag, bool force_moment);
    void moment_exact(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_stepwise,int erm_alg, double swthresh, int tsk_ttl, int tsk_id, double swfdr, bool stepforwardonly, double sw_rsq);
}

#endif /* l3_ewas_hpp */

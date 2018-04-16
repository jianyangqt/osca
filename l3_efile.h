//
//  l3_efile.h
//  osc
//
//  Created by Futao Zhang on 11/04/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l3_efile__
#define __osc__l3_efile__
#include "l2_efile.h"
#include "l2_reml.h"
#include "l2_bfile.h"

namespace EFILE {
    void merge_beed(char* outfileName, char* befileFlistName, char* problstName, char* problst2exclde,char* genelistName,  int chr, char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename,char* probe2rm,char* indilstName,char* indilst2remove,bool beta2m,bool m2beta);
    
    void make_beed(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta,double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,bool adjprb, char* covfileName,char* qcovfileName, bool enveff, char* effprblstfname,char* efffname);
    void make_efile(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,bool t_flag, bool impute_mean_flag);
    void make_erm(char* outFileName, char* efileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,bool transposed, int efileType,bool no_fid_flag,int valueType);
    
    void mlma(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* subtract_erm_file, bool m_erm_flag, bool within_family,char* priors,char* priors_var,bool no_constrain, int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_no_adj_covar, double percentage_out, double lambda_wind,bool fastlinear);
    void mlma_loco(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, int MaxIter, char* priors, char* priors_var, bool no_constrain, bool no_adj_covar,int reml_mtd,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, int autosome_num, double percentage_out, double lambda_wind, bool fastlinear);
    void pca(char* outFileName, char* grm_file, char* indilstName, char* indilst2remove, double grm_cutoff, bool erm_cutoff_2sides, bool merge_grm_flag, int out_pc_num);
    
    void diff(char* befileName1,char* befileName2);
    void getRefactor(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* covfileName,char* qcovfileName, int celltype_num, int dmr_num,int out_pc_num);
    void  assoc(char* outfileName, char* befileName, char* problstName, char* problst2exclde, char* genelistName,  int chr,char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename, char* probe2exclde, char* indilstName, char* indilst2remove,char* phenofileName,char* mpheno,  char* covfileName, char* qcovfileName, double std_thresh ,double upperBeta,double lowerBeta,bool est_eff_n=false);
    void  linear(char* outfileName, char* befileName, char* problstName, char* problst2exclde, char* genelistName,  int chr,char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename, char* probe2exclde, char* indilstName, char* indilst2remove,char* phenofileName,char* mpheno,  char* covfileName, char* qcovfileName, double std_thresh ,double upperBeta,double lowerBeta,int tsk_ttl,int tsk_id,bool fastlinear);
    
    void EWAS_simu(char* outFileName, char* befileName, int simu_num, char* sigCpG_file, int case_num, int control_num, double hsq, double K, int seed, int eff_mod,bool simu_residual=false);
    
    void fit_reml(char* outFileName, char* phenofileName,char* mpheno,bool erm_bin_flag,bool grm_bin_flag, int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* grm_file,bool m_erm_flag, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_no_adj_covar, bool pred_rand_eff, bool est_fix_eff,bool no_lrt,double prevalence, bool mlmassoc,vector<int> drop,char* indilstName, char* indilst2remove, char* sex_file, double grm_cutoff, bool erm_cutoff_2sides, double adj_grm_fac, int dosage_compen,bool prt_residiual);
    
    void extract_inden_probes(char* outFileName, char* befileName, int extr_num, double ldr,  int seed, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool rand_eff=false);
    void  testQAssoc(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, bool est_eff_n=false);
    void  testLinear(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, MatrixXd &COV_plus);
    void  testLinear_fast(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, MatrixXd &COV_plus);
    void getPrbVarianceMean(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool varflag,bool meanflag);
    void blup_probe(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,double std_thresh,double upperBeta,double lowerBeta, char* blup_indi_file);
    void scoreIndividuals(char* outFileName,char* befileName, char* betaFileName, int col_prb, int col_score, bool hasHeader, char* phenofileName, char* mpheno,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, double std_thresh, bool impute_mean_flag);
    void  adjprobe(eInfo* einfo);
    void  addEnvEff(eInfo* einfo, char* prblstfname,char* efffname);
    double get_lambda(vector<ASSOCRLT> &assoc_rlts);
    void write_assoc_rlt(vector<ASSOCRLT> &assoc_rlts,char* outfileName);
}
#endif /* defined(__osc__l3_efile__) */

//
//  l3_vqtl.hpp
//  osc
//
//  Created by Futao Zhang on 29/08/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#ifndef l3_vqtl_hpp
#define l3_vqtl_hpp

#include "l2_efile.h"
#include "l2_bfile.h"
#include "l2_reml.h"
#include "l2_besd.hpp"
using namespace BFILE;
using namespace EFILE;
using namespace SMR;
namespace VQTL {
    void load_vqtl_workspace(eInfo* einfo,bInfo* bdata, char* efileName, \
        char* befileName, char* phenofileName, char* bFileName, bool transposed, \
        int efileType, char* problstName, char * problst2exclde, char * genelistName, \
        int chr, char * prbname, char* fromprbname, char * toprbname, int prbWind, \
        int fromprbkb, int toprbkb, bool prbwindFlag, char * genename, \
        char * probe2exclde, char * indilstName, char * indilst2remove, \
        bool no_fid_flag, int valueType, bool beta2m, bool m2beta, double std_thresh, \
        double upperBeta, double lowerBeta, char * dpvalfName, double dp_thresh, \
        double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, \
        int autosome_num, char * snplstName, char * snplst2exclde, int tsk_ttl, \
        int tsk_id, char * covfileName, char * qcovfileName, char * grm_file=NULL, \
        int xqtlNO = 0, double zeroratio = 0.2, eInfo * eCov = NULL, \
        char * covbodfileName = NULL, char * covefileName = NULL, \
        bool transopse_ecov = false);

    void V_QTL(char * outFileName, char * efileName, char * befileName, \
        char * phenofileName, char * bFileName, bool transposed, int efileType, \
        char * problstName, char * problst2exclde, char * genelistName, int chr, \
        char * prbname, char * fromprbname, char * toprbname, int prbWind, \
        int fromprbkb, int toprbkb, bool prbwindFlag, char * genename, \
        char * probe2exclde, char * indilstName, char * indilst2remove, \
        bool no_fid_flag, int valueType, bool beta2m, bool m2beta, \
        double std_thresh, double upperBeta, double lowerBeta, char * dpvalfName, \
        double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, \
        double mssratio_prob, int autosome_num, double maf, char * snplstName, \
        char * snplst2exclde, int tsk_ttl, int tsk_id, int vqtl_mtd, \
        char * covfileName, char * qcovfileName, bool tosmrflag, bool cis_flag, \
        int cis_itvl);

    void eQTL(char * outFileName,  char * efileName, char * befileName, char * bFileName, \
        bool transposed, int efileType, char * problstName, char * problst2exclde, \
        char * genelistName, int chr, char * prbname, char * fromprbname, char * toprbname, \
        int prbWind, int fromprbkb, int toprbkb, bool prbwindFlag, char * genename, \
        char * probe2exclde, char * indilstName, char * indilst2remove, bool no_fid_flag, \
        int valueType, bool beta2m, bool m2beta, double std_thresh, double upperBeta, \
        double lowerBeta, char * dpvalfName, double dp_thresh, double prb_thresh, \
        double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num, \
        double maf, char * snplstName, char * snplst2exclde, int tsk_ttl, int tsk_id, \
        char * covfileName, char * qcovfileName, bool tosmrflag, bool nofastlinear, \
        bool cis_flag, int cis_itvl);

    void eQTL_MLM(char * outFileName, char * efileName, char * befileName, \
        char * bFileName, bool transposed, int efileType, char * problstName, \
        char * problst2exclde, char * genelistName, int chr, char * prbname, \
        char * fromprbname, char * toprbname, int prbWind, int fromprbkb, \
        int toprbkb, bool prbwindFlag, char * genename, char * probe2exclde, \
        char * indilstName, char * indilst2remove, bool no_fid_flag, \
        int valueType, bool beta2m, bool m2beta, double std_thresh, double upperBeta, \
        double lowerBeta, char * dpvalfName, double dp_thresh, double prb_thresh, \
        double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num, \
        double maf, char * snplstName, char * snplst2exclde, int tsk_ttl, int tsk_id, \
        char * covfileName, char * qcovfileName, bool tosmrflag, bool cis_flag, \
        int cis_itvl, char * grm_file, bool grm_bin_flag, bool no_constrain, \
        int reml_mtd, int MaxIter, bool nopreadj_covar);

    void sQTL(char* outFileName, char * efileName, char * befileName, char * bFileName, \
        bool transposed, int efileType, char * problstName, char * problst2exclde, \
        char * genelistName, int chr, char * prbname, char * fromprbname, char * toprbname, \
        int prbWind, int fromprbkb, int toprbkb, bool prbwindFlag, char * genename, \
        char * probe2exclde, char * indilstName, char * indilst2remove, \
        bool no_fid_flag, int valueType, bool beta2m, bool m2beta, double std_thresh, \
        double upperBeta, double lowerBeta, char * dpvalfName, double dp_thresh, \
        double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, \
        int autosome_num, double maf, char * snplstName, char * snplst2exclde, \
        int tsk_ttl, int tsk_id, char * covfileName, char * qcovfileName, bool tosmrflag, \
        bool nofastlinear, bool cis_flag, int cis_itvl, double zeroratio, double call, \
        char * annofileName, char * covbodfileName, char* covefileName, bool transopse_ecov, bool use_top_p);

    void ssQTL(char* outFileName, char* beqtlFileName, char* problstName,char* problst2exclde, \
        char * genelistName, int chr, int prbchr, char * prbname, char * fromprbname, \
        char * toprbname, int prbWind, int fromprbkb, int toprbkb, bool prbwindFlag, \
        char * genename, char * probe2exclde, int autosome_num, double maf, \
        char * snplstName, char * snplst2exclde, int snpchr, char * snprs, \
        char * fromsnprs, char * tosnprs, int snpWind, int fromsnpkb, int tosnpkb, \
        bool snpwindFlag, char * snprs2exclde, int tsk_ttl, int tsk_id, \
        bool tosmrflag, bool nofastlinear, bool cis_flag,int cis_itvl, \
        char * annofileName, double pmecs, int nmecs, bool use_top_p);

}
#endif /* l3_vqtl_hpp */

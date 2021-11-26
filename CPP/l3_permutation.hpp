#ifndef l3_permutation_hpp
#define l3_permutation_hpp

#include "l2_efile.h"
#include "l2_bfile.h"
#include "l2_reml.h"
#include "l2_besd.hpp"
#include "l3_vqtl.hpp"
#include <algorithm>
#include <ctime>
#include <cstdlib>
using namespace BFILE;
using namespace EFILE;
using namespace SMR;
using namespace VQTL;

typedef struct OUTDT {
    char probe_id[1024];
    uint32_t probe_pos;
    char gene_name[1024];
    char chrom;
    char orientation;
    uint32_t snp_contained;
    char top_snp_id[1024];
    char top_snp_chrom;
    uint32_t top_snp_pos;
    double p_nominal;
    double beta_ml1;
    double beta_ml2;
    double pemp;
    double pbml;
    struct OUTDT * next;
} output_data;

namespace PERMU 
{
    void 
    permu_sqtl (char* outFileName, char * efileName, char * befileName, char * bFileName,
        bool transposed, int efileType, char * problstName, char * problst2exclde,
        char * genelistName, int chr, char * prbname, char * fromprbname, char * toprbname,
        int prbWind, int fromprbkb, int toprbkb, bool prbwindFlag, char * genename,
        char * probe2exclde, char * indilstName, char * indilst2remove,
        bool no_fid_flag, int valueType, bool beta2m, bool m2beta, double std_thresh,
        double upperBeta, double lowerBeta, char * dpvalfName, double dp_thresh,
        double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,
        int autosome_num, double maf, char * snplstName, char * snplst2exclde,
        int tsk_ttl, int tsk_id, char * covfileName, char * qcovfileName,
        bool nofastlinear, bool cis_flag, int cis_itvl, double zeroratio, double call,
        char * annofileName, char * covbodfileName, char* covefileName, bool transopse_ecov, 
        bool use_top_p, bool trans_flag, int trans_itvl, uint32_t permute_times);

}
#endif

//
//  l4_osc.h
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l4_osc__
#define __osc__l4_osc__

#include "l3_vqtl.hpp"
#include "l3_smr.hpp"
#include "l3_ewas.hpp"
#include "l3_efile.h"
#include "l3_gwas.hpp"

void option(int option_num, char* option_str[]);

static inline bool not_in_flags(vector<string> &flags, string str)
{
    return find(flags.begin(),flags.end(),str) == flags.end();
}

static inline void FLAGS_VALID_CK(int option_num, char* option_str[])
{
    const char *flgs[] = {"--efile","--tefile","--make-bod","--befile",
        "--gene-expression","--methylation","--methylation-beta","--methylation-m",
        "--out","--make-efile", "--probe","--probe-rm","--pheno","--mpheno",
        "--make-orm","--make-orm-bin","--make-orm-gz","--orm-alg","--moa",
        "--mlma-loco","--covar","--qcovar","--orm","--orm-bin","--keep",
        "--extract-probe","--exclude-probe","--no-fid","--genes","--gene",
        "--remove","--m2beta","--beta2m", "--befile-flist","--diff","--refactor",
        "--celltype-num","--dmr-num","--thread-num","--autosome-num","--update-opi",
        "--assoc","--sd-min" ,"--linear","--pca","--simu-qt","--simu-cc",
        "--simu-causal-loci","--simu-rsq","--simu-k","--simu-seed","--simu-eff-mod",
        "--reml","--reml-priors","--reml-priors-var","--reml-pred-rand",
        "--reml-est-fix","--reml-no-lrt","--prevalence","--probes-independent",
        "--ld-rsq", "--lxpo","--upper-beta","--lower-beta","--detection-pval-file",
        "--dpval-thresh","--ratio-probe","--ratio-sample","--dpval-mth","--chr",
        "--eff-n","--get-variance","--make-tefile","--get-mean","--missing-ratio-probe",
        "--blup-probe","--score","--impute-mean","--score-has-header","--from-probe",
        "--to-probe","--probe-wind","--from-probe-kb","--to-probe-kb","--orm-cutoff",
        "--multi-orm","--reml-maxit","--reml-alg","--orm-cutoff","--reml-no-constrain",
        "--task-num","--task-id","--vqtl","--bfile","--maf","--extract-snp","--exclude-snp",
        "--vqtl-mtd","--adj-probe","--grm","--grm-bin","--output-residual",
        "--simu-residual","--eqtl","--sqtl","--query","--beqtl-summary","--probe-chr",
        "--snp-chr","--snp","--from-snp","--to-snp","--snp-wind","--snp-rm","--cis-wind",
        "--besd-flist","--meta","--mecs","--pmecs","--gwas-flist","--mecs-mth","--besd-shrink",
        "--make-besd-dense","--to-smr","--make-besd","--add-eff","--eff-file","--eff-probe",
        "--cor-mat","--pcc-z","--prt-mid","--std-probe","--freq-file","--var-file",
        "--lambda-range","--fast-linear","--orm-cutoff-2sides","--force-mlm","--bin-num",
        "--bin-mth","--npcs","--fixed-pc","--slct-dom-pc","--no-preadj-covar",
        "--clustering-mth","--pthresh","--moment","--moment-wind","--moment-num",
        "--moment-alt-pcs","--simu-causal-loci2","--simu-rsq2","--moment-cor",
        "--num-rand-comp","--rint-probe","--moment-r2","--feature-slct-mtd",
        "--moment-percent","--moment-prior","--approximate-num","--stepwise-slct",
        "--no-prior-var","--logistic","--missing-ratio-indi","--no-fast-linear",
        "--cis","--mlm","--ewas-flist","--pairwise-common","--all-common",
        "--moa-exact","--ewas-summary","--gc","--pstep","--simu-reverse","--moment-exact",
        "--loud","--zero-ratio-probe","--call","--bed","--stepwise-fdr","--stepwise-rsq",
        "--moment2-beta","--tpm","--stepwise-logistic","--stepwise-forward","--covar-bod",
        "--covar-efile","--covar-tefile","--nmecs","--reverse-assoc","--fdr","--cor-r2",
        "--make-bld","--save-r2","--r2-thresh","--moment-force"
    };

    vector<string> flags(flgs, flgs + sizeof(flgs)/sizeof(flgs[0]));

    if(option_num<2)
    {
        LOGPRINTF("Flags include:\n");
        int cur_mark=0;
        for(int i=0;i<flags.size();i++)
        {
            int tmp=i>>2;
            if(tmp>cur_mark)
            {
                cout<<endl;
                cur_mark=tmp;
            }
            LOGPRINTF("%s,",flags[i].c_str());
        }
        LOGPRINTF("\n");
        TERMINATE();
    }
    for(int i=0;i<option_num;i++)
    {
        if(has_prefix(option_str[i],"--"))
            if(not_in_flags(flags, option_str[i]))
            {
                LOGPRINTF("%s: Invalid option\n",option_str[i]);
                TERMINATE();
            }
    }

}

static inline void FLAG_VALID_CK(string str, char* flag)
{
    if(flag==NULL || has_prefix(flag, "--"))
    {
        LOGPRINTF("Please verify the flag %s: \n",str.c_str());
        TERMINATE();

    }
}

#endif /* defined(__osc__l4_osc__) */

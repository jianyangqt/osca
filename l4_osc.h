//
//  l4_osc.h
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l4_osc__
#define __osc__l4_osc__


#include "l3_efile.h"
#include "l3_vqtl.hpp"

void option(int option_num, char* option_str[]);
static inline bool not_in_flags(vector<string> &flags, string str)
{
    return find(flags.begin(),flags.end(),str) == flags.end();
}

static inline void FLAGS_VALID_CK(int option_num, char* option_str[])
{
    const char *flgs[] = {"--efile","--tefile","--make-bod","--befile","--gene-expression","--methylation","--methylation-beta","--methylation-m","--out","--make-efile", "--probe","--probe-rm","--pheno","--mpheno","--make-orm","--make-orm-bin","--make-orm-gz","--orm-alg","--mlma","--mlma-loco","--covar","--qcovar","--orm","--orm-bin","--keep","--extract-probe","--exclude-probe","--no-fid","--genes","--gene","--remove","--m2beta","--beta2m", "--befile-flist","--diff","--refactor","--celltype-num","--dmr-num","--thread-num","--autosome-num","--update-opi","--assoc","--std" ,"--linear","--pca","--simu-qt","--simu-cc","--simu-causal-loci","--simu-hsq","--simu-k","--simu-seed","--simu-eff-mod","--reml","--reml-pred-rand","--reml-est-fix","--reml-no-lrt","--prevalence","--probes-independent","--ld-rsq", "--lxpo","--upper-beta","--lower-beta","--detection-pval-file","--dpval-thresh","--ratio-probe","--ratio-sample","--dpval-mth","--chr","--eff-n","--get-variance","--make-tefile","--get-mean","--missing-ratio-probe","--blup-probe","--score","--impute-mean","--score-has-header","--from-probe","--to-probe","--probe-wind","--from-probe-kb","--to-probe-kb","--orm-cutoff","--merge-orm","--reml-maxit","--reml-alg","--orm-cutoff","--reml-no-constrain","--task-total","--task-id","--vqtl","--bfile","--maf","--extract-snp","--exclude-snp","--vqtl-mtd","--adj-probe","--grm","--grm-bin","--output-residual","--simu-residual"};
    
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

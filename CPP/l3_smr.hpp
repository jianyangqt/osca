//
//  l3_smr.hpp
//  osc
//
//  Created by Futao Zhang on 20/11/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#ifndef l3_smr_hpp
#define l3_smr_hpp

#include "l2_besd.hpp"

namespace SMR {

    typedef struct{
        int* ptr; //vector of int. using in combining besd files and meta analysis to save the probe index of the according besd file.
        char* probeId;
        char* genename;
        char* esdpath;
        char* bfilepath;
        int probechr;
        int gd;
        int bp;
        char orien;
    } smr_probeinfo;

    typedef struct{
        int* rstr; //vector of int. using in combining besd files and meta analysis to save the rs index of the according besd file.
        bool* revs; //vector of bool. using in combining besd files and meta analysis to save if the sign of beta should be changed the according besd file.
        char* snprs;
        char* a1;
        char* a2;
        int snpchr;
        int gd;
        int bp;
        float beta;
        float se;
        float freq;
        float estn;
    } smr_snpinfo;

    typedef struct{
        char* snprs;
        char* a1;
        char* a2;
        float* beta;
        float* se;
        float freq;
        float estn;
        int bp;
    } gwasinfo;

    void query_besd(char* outFileName,char* beqtlFileName, char* snplstName, char* snplst2exclde, char* problstName,char* problst2exclde, char* genelistName, double plookup, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl, char* probe2exclde, char* snprs2exclde);
     void make_besd(char* outFileName,char* beqtlFileName, char* snplstName, char* snplst2exclde, char* problstName,char* problst2exclde, char* genelistName, double plookup, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl, char* probe2exclde, char* snprs2exclde, bool save_dense_flag,bool tosmrflag, bool besd_shrink_flag,bool stdprb, char* frqFName,char* varFName);
    void meta(char* besdlistFileName, char* outFileName, int meta_mth, double pthresh, bool cis_flag, int cis_itvl,int nmecs,char* problstName, char* problst2exclde, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* snplstName, char* snplst2exclde,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb, bool smpwindFlag, char* probe2excld, char* snprs2excld, bool trans_meta);
    void meta_gwas(char* gwaslistFileName, char* ewaslistFileName, char* outFileName, int meta_mth, double pthresh, int mecs_mth, char* corMatFName, char* snplstName,char* problstName, bool zflag, bool out_comm_flag, bool all_comm_flag);
    void gc_ewas(char* outFileName, char* ewasFileName);
}

#endif /* l3_smr_hpp */

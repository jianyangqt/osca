//
//  l2_besd.hpp
//  osc
//
//  Created by Futao Zhang on 10/11/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#ifndef l2_besd_hpp
#define l2_besd_hpp

#include "l2_efile.h"
#include "l2_bfile.h"
#include "l2_reml.h"

#define SMR_DENSE_1 0 // 0x00000000 + floats  : beta values (followed by se values) for each probe across all the snps are adjacent.
#define SMR_SPARSE_3F 0x40400000 // 0x40400000: uint32_t + uint64_t + uint64_ts + uint32_ts + floats
#define SMR_DENSE_3 5  // RESERVEDUNITS*ints + floats (indicator+samplesize+snpnumber+probenumber+ 12*-9s + values) [SMR default and OSCA default]
#define SMR_SPARSE_3 3 // RESERVEDUNITS*ints + uint64_t + uint64_ts + uint32_ts + floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s +valnumber+cols+rowids+betases) [SMR default]
#define RESERVEDUNITS 16

#define OSCA_SPARSE_1 1 // 0x00000001: RESERVEDUNITS*ints + uint64_t  + uint64_ts + uint32_ts + floats: value number + (half uint64_ts and half uint32_ts of SMR_SPARSE_3) [OSCA default]
#define OSCA_DENSE_1 4 // 0x00000004: RESERVEDUNITS*ints  + floats  :  <beta, se> for each SNP across all the probes are adjacent.


using namespace BFILE;
using namespace EFILE;
namespace SMR {
    typedef struct{
        vector<int> _esi_chr;
        vector<string> _esi_rs;
        vector<int> _esi_gd;
        vector<int> _esi_bp;
        vector<string> _esi_allele1;
        vector<string> _esi_allele2;
        vector<int> _esi_include; // initialized in the readesi
        map<string,int> _snp_name_map;
        vector<float> _esi_freq;
        
        vector<int> _epi_chr;
        vector<string> _epi_prbID;
        vector<int> _epi_gd;
        vector<int> _epi_bp;
        vector<string> _epi_gene;
        vector<char> _epi_orien;
        vector<int> _include; // initialized in the readepi
        map<string,int> _probe_name_map;
        vector<double> _epi_var;
        vector<int> _epi_start; // if no probe sequence region input, its size should be 0. for the probe not in probe sequence file, the value should be set as -9, no technical eQTL would be removed from this probe.
        vector<int> _epi_end;
        
        //for sparse
        vector<uint64_t> _cols; // different content from SMR
        vector<uint32_t> _rowid; // different content from SMR
        vector<float> _val; // the same content as SMR
        // for dense
        vector< vector<float> > _bxz; // first dimension is probe, second is snp
        vector< vector<float> > _sexz;
        
        uint64_t _probNum;
        uint64_t _snpNum;
        uint64_t _valNum;
        int _sampleNum;
        
    } eqtlInfo;

    typedef struct{
        long snpNum;
        vector<string> snpName;
        vector<int> snpBp;
        vector<string> allele_1;
        vector<string> allele_2;
        vector<double> freq;
        vector<double> byz;
        vector<double> seyz;
        vector<double> pvalue;
        vector<int> splSize;
        vector<int> _include;
        map<string,int> _snp_name_map;
    } gwasData;

    void write_smr_epi(char* outFileName, eInfo* einfo);
    void write_smr_esi(char* outFileName, bInfo* binfo);
    void read_smr_epifile(eqtlInfo* eqtlinfo, char* epiFileName);
    void read_smr_esifile(eqtlInfo* eqtlinfo, char* esiFileName);
    void read_smr_besdfile(eqtlInfo* eqtlinfo, char* besdFileName);
    void update_epi(eqtlInfo* eqtlinfo);
    void update_esi(eqtlInfo* eqtlinfo);
    void smr_esi_man(eqtlInfo* eqtlinfo,char* snplstName, char* snplst2exclde, int chr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, int cis_itvl,const char* prbname,char* snprs2exclde);
    void smr_epi_man(eqtlInfo* eqtlinfo,char* problstName,char* problst2exclde,char* genelistName, int chr, int probchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde);
    void get_BesdHeaders(char* besdFileName, vector<int> &headers);
    void read_gwas_data(gwasData* gdata, char* gwasFileName);
    void read_ewas_data(gwasData* gdata, char* ewasFileName);
    void write_smr_esi(char* outFileName, eqtlInfo* eqtlinfo);
    void write_smr_epi(char* outFileName, eqtlInfo* eqtlinfo);
    void write_s2s_besd(char* outFileName, eqtlInfo* eqtlinfo, bool tosmrflag=false);
    void get_shrink_null(eqtlInfo* eqtlinfo,vector<string> &nullprbs, vector<string> &nullsnps);
    void extract_prb_dense(FILE* fptr,  uint64_t pid, uint64_t epinum,uint64_t esinum, vector<float> &betases);
    void write_d2d_besd(char* outFileName, eqtlInfo* eqtlinfo, char* inputname, bool stdprb);
    void extract_prb_sparse(FILE* fptr, uint64_t pid, uint64_t probnum,vector<uint32_t> &row_ids, vector<float> &betases);
    void read_probevarfile(eqtlInfo* eqtlinfo, char* vpFileName);
}

#endif /* l2_besd_hpp */

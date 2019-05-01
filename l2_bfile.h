//
//  l2_bfile.h
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l2_bfile__
#define __osc__l2_bfile__

#include "l1_op_geno.h"
#include "l2_efile.h"

namespace BFILE{

    typedef struct{
        // bim file
        int _autosome_num;
        vector<int> _chr;
        vector<string> _snp_name;
        map<string, int> _snp_name_map;
        vector<double> _genet_dst;
        vector<int> _bp;
        vector<string> _allele1;
        vector<string> _allele2;
        vector<string> _ref_A; // reference allele
        vector<string> _other_A; // the other allele
        uint64_t _snp_num;
        vector<double> _rc_rate;
        vector<int> _include; // initialized in the read_bimfile()
        VectorXf _maf;
        
        // fam file
        vector<string> _fid;
        vector<string> _pid;
        map<string, int> _id_map;
        vector<string> _fa_id;
        vector<string> _mo_id;
        vector<int> _sex;
        vector<double> _pheno;
        uint64_t _indi_num;
        vector<int> _keep; // initialized in the read_famfile()
        MatrixXf _varcmp_Py; // BLUP solution to the total genetic effects of individuals
        
        // bed file
        vector< vector<bool> > _snp_1;
        vector< vector<bool> > _snp_2;
        
        // imputed data
        bool _dosage_flag;
        vector< vector<float> > _geno_dose;
        vector<double> _impRsq;
        
        // genotypes
        MatrixXf _geno;
        
        vector<double> _mu;
        
        MatrixXf _grm_N;
        MatrixXd _grm;
        
    } bInfo;
    
    void keep_indi(bInfo* bdata,string indi_list_file);
    void remove_indi(bInfo* bdata, string indi_list_file);
    void extract_snp(bInfo* bdata,string snplistfile);
    void exclude_snp(bInfo* bdata,string snplistfile);
    void calcu_mu(bInfo* bdata, bool ssq_flag=false);
    bool make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool mu=false);
    bool make_XMat(bInfo* bdata, int start, int slide_wind, MatrixXd &X, bool mu=false);
    bool make_XMat_subset(bInfo* bdata, int start, int slide_wind, MatrixXd &X, bool standardise=false);
    bool make_XMat_subset(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool standardise=false);
    void read_famfile(bInfo* bdata, string famfile);
    void read_bimfile(bInfo* bdata,string bimfile);
    void read_bedfile(bInfo* bdata, string bedfile);
    void filter_snp_maf(bInfo* bdata,double maf);
    void make_grm(bInfo* bdata, int grm_mtd,  bool diag_f3_flag=false);
    void read_grm(bInfo* bdata,bool grm_bin_flag, string grm_file, vector<string> &grm_id, bool read_id_only, bool dont_read_N);
    int read_grm_id(bInfo* bdata,string grm_file, vector<string> &grm_id);
}

#endif /* defined(__osc__l2_bfile__) */

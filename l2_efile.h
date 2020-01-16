//
//  l2_efile.h
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l2_efile__
#define __osc__l2_efile__

#include "l1_op_geno.h"
#include "l1_stat.hpp"

#define GENEEXPRESSION 0
#define METHYLATION 1
#define OTHERDATA 2
#define BETAVALUE 0
#define MVALUE 1
#define VALUE 2
#define TPM 0

namespace EFILE {
    typedef struct{
        
        uint32_t _eType;
        // epi file
        uint32_t autosome_num;
        uint64_t _epi_num;
        vector<int> _epi_chr;
        vector<string> _epi_prb;
        vector<int> _epi_gd;
        vector<int> _epi_bp;
        vector<string> _epi_gene;
        vector<char> _epi_orien;
        vector<int> _epi_include;
        map<string,int> _epi_map;
        
        // eii file
        uint64_t _eii_num;
        vector<string> _eii_fid;
        vector<string> _eii_iid;
        vector<string> _eii_fa_id;
        vector<string> _eii_mo_id;
        vector<int> _eii_sex;
        vector<double> _eii_pheno;
        uint32_t _eii_pheno_num;
        vector<string> _eii_cov; //cov major, values belong to the same covariate are adjacent
        uint32_t _eii_cov_num;
        vector<double> _eii_qcov;
        uint32_t _eii_qcov_num;
        vector<int> _eii_include;
        map<string, int> _eii_map;
        
        // eed file
        uint32_t _valType; // 0 beta, 1 m, 2 other
        vector< double > _val; //probe major
        vector< double > _mu; //probe mean
        vector< double > _var; //probe variance
        
        MatrixXf _grm_N;
        MatrixXd _grm;
        double * _profile;
        double * _grm_ptr;
        
        //
        int _reml_mtd;
        int _reml_max_iter;
        int _V_inv_mtd;
        bool _reml_force_inv;
        bool _reml_force_converge;
        bool _reml_no_converge;
        bool _reml_AI_not_invertible;
        MatrixXd _P;
        vector<int> _r_indx;
        vector<int> _r_indx_drop;
        vector<string> _var_name;
        vector<double> _varcmp;
        vector<string> _hsq_name;
        
        bool _within_family;
        vector<int> _fam_brk_pnt;
        vector< SparseMatrix<double> > _Asp;
        vector< SparseMatrix<double> > _Asp_prev;
        double _y_Ssq;
        vector<double> _fixed_rg_val;
        bool _reml_fixed_var;
        VectorXd _b;
        VectorXd _se;
        
        //reserved
        bool _bivar_reml;
        bool _bivar_no_constrain;
        bool _ignore_Ce;
        double _y2_Ssq;
        vector< vector<int> > _bivar_pos;
        vector< vector<int> > _bivar_pos_prev;
        
        double _ncase;
        double _ncase2;
        bool _flag_CC;
        bool _flag_CC2;
        
        MatrixXd _varcmp_Py; // BLUP solution to the total genetic effects of individuals
        
    } eInfo;
    
    typedef struct{
        int* itr;
        char* fid;
        char* iid;
        char* fa_id;
        char* mo_id;
        int sex;
        double pheno;
    } indiinfolst;
    
    typedef struct{
        int* ptr;
        char* probeId;
        char* genename;
        
        int probechr;
        int gd;
        int bp;
        char orien;
    } probeinfolst;
    
    typedef struct{
        char* PROBE;
        char* GENE;
        double BETA;
        double SE;
        double T;
        double R2;
        double PVAL;
        double NMISS;
        int BP;
        int CHR;
        char OREN;
    } ASSOCRLT;
    

    string getFileType(uint32_t tid);
    string getValType(uint32_t fid, uint32_t vid);
    void  read_efile(char* eFileName, eInfo* einfo,uint32_t filetype, bool no_fid_flag, int valueType);
    void  read_pheno2(char* eFileName, eInfo* einfo,int colid);
    void  read_efile_t(char* eFileName, eInfo* einfo,uint32_t filetype, bool no_fid_flag,int valueType);
    void read_eii(char* eiiFileName,eInfo* einfo);
    void read_epi(char* epiFileName,eInfo* einfo);
    void update_epifile(char* befileName,char* s_epiFileName);
    void read_beed(char* beedFileName,eInfo* einfo);
    void write_eii(char* outFileName, eInfo* einfo);
    void write_epi(char* outFileName, eInfo* einfo);
    void write_beed(char* outFileName, eInfo* einfo);
    void init_einfo(eInfo* einfo);
    void update_eii(eInfo* einfo);
    void update_epi(eInfo* einfo);
    void keep_indi(eInfo* einfo,string indi_list_file);
    void remove_indi(eInfo* einfo, string indi_list_file);
    void extract_probe_by_chr(eInfo* einfo, int prbchr);
    void extract_probe(eInfo* einfo,string problstName);
    void extract_probe_by_gene(eInfo* einfo, string genelistName);
    void extract_probe(eInfo* einfo, string prbname, int prbWind);
    void extract_single_probe(eInfo* einfo, string prbname);
    void extract_probe(eInfo* einfo, string fromprbname, string toprbname);
    void extract_probe(eInfo* einfo, int fromprbkb, int toprbkb, int chr);
    void extract_probe_by_single_gene(eInfo* einfo, string genename);
    void extract_probe(eInfo* einfo, int tsk_ttl, int tsk_id);
    void exclude_probe(eInfo* einfo, string problstName);
    void exclude_single_probe(eInfo* einfo, string prbname);
    void epi_man(eInfo* einfo,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde);
    void eii_man(eInfo* einfo,char* indilstName,char* indilst2remove);
    void write_efile(char* outFileName, eInfo* einfo,bool impute_mean_flag);
    void write_tefile(char* outFileName, eInfo* einfo,bool impute_mean_flag);
    void read_phen(eInfo* einfo, string phen_file, char* mpheno, bool mvFlg=false);
    void read_cc(eInfo* einfo, string phen_file, char* mpheno, bool mvFlg=false);
    void make_erm(eInfo* einfo, int erm_mtd=0, bool output_bin=true, char* outFileName=NULL, bool output_profile=false, bool have_stand = false);
    void make_erm(eInfo* einfo, MatrixXd &VZ, int erm_mtd=0, bool output_bin=true, char* outFileName=NULL, bool output_profile=false, bool have_stand = false);
    void output_grm(eInfo* einfo, string _out, bool output_grm_bin);
    void read_cov(eInfo* einfo, string cov_file, bool qcovFlg);
    void read_grm_gz(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only);
    void read_grm_bin(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only, bool dont_read_N);
    void read_grm(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only, bool dont_read_N, bool grm_bin_flag);
    void merge_grm(eInfo* einfo, char* merge_grm_file);
    void rm_cor_indi(eInfo* einfo, double grm_cutoff,bool erm_cutoff_2sides);
    void update_sex(eInfo* einfo, char* sex_file);
    void adj_grm(eInfo* einfo, double adj_grm_fac);
    void dc(eInfo* einfo, int dosage_compen);
    bool moment_eligibility_ck(eInfo* einfo);
    void std_probe(eInfo* einfo, vector< vector<bool> > &X_bool, bool divid_by_std, MatrixXd &_probe_data,bool output_profile=false);
    void std_probe_in_place( bool divid_by_std, MatrixXd &_probe_data);
    void beta_2_m(eInfo* einfo);
    void m_2_beta(eInfo* einfo);
    void std_probe_filtering( eInfo* einfo, double std_thresh);
    void filtering_constitutive_probes( eInfo* einfo, double upper_beta_thresh,double lower_beta_thresh);
    void filtering_with_detpval(eInfo* einfo, char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh,int mth, bool no_fid_flag);
    void filtering_with_missingratio(eInfo* einfo,double missratioprobe);
    void filtering_indi_missingratio(eInfo* einfo,double missratioindi);
    void cal_var_mean(eInfo* einfo, bool mean_flag, bool var_flag);
    void load_workspace(eInfo* einfo,char* efileName, char* befileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num);
    void  adjprobe(eInfo* einfo);
    void  fast_adjprobe(eInfo* einfo);
    void  fast_adjprobe(eInfo* einfo, MatrixXd &X,  MatrixXd &XtXi);
    int construct_X(eInfo* einfo, vector<MatrixXd> &E_float, MatrixXd &qE_float, MatrixXd &_X);
    void free_indilist(vector<indiinfolst> &a);
    void free_probelist(vector<probeinfolst> &a);
    void free_assoclist(vector<ASSOCRLT> &a);
    bool check_case_control(double &ncase,  VectorXd &y) ;
    void update_startend(int length, int tsk_ttl, int tsk_id, int &start, int &end);
    void extract_sqtl_probe(eInfo* einfo,int tsk_ttl,int tsk_id);
    void filtering_with_zeroratio(eInfo* einfo,double zeroratioprobe);
}
#endif /* defined(__osc__l2_efile__) */

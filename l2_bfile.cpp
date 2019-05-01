//
//  l2_bfile.cpp
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l2_bfile.h"
namespace BFILE{
   
    void read_famfile(bInfo* bdata, string famfile) {
        bdata->_autosome_num = 22;
        ifstream Fam(famfile.c_str());
        if (!Fam)
        {
            sprintf(logbuf, "Error: can not open the file %s to read.\n", famfile.c_str());
            logprintb();
           TERMINATE();

        }
        sprintf(logbuf, "Reading PLINK FAM file from %s .\n", famfile.c_str());
        logprintb();
        string str_buf;
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        while (Fam) {
            Fam >> str_buf;
            if (Fam.eof()) break;
            bdata->_fid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_pid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_fa_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_mo_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_sex.push_back(atoi(str_buf.c_str()));
            Fam >> str_buf;
            bdata->_pheno.push_back(atoi(str_buf.c_str()));
        }
        Fam.clear();
        Fam.close();
        bdata->_indi_num = bdata->_fid.size();
        sprintf(logbuf, "%llu individuals to be included from %s .\n", bdata->_indi_num, famfile.c_str());
        logprintb();
        
        // Initialize _keep
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        uint64_t size = 0;
        for (int i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
            if (size == bdata->_id_map.size())
            {
                sprintf(logbuf, "Error: Duplicate individual ID found: \" %s \t %s \".\n", bdata->_fid[i].c_str(), bdata->_pid[i].c_str());
                logprintb();
               TERMINATE();
            }
            size = bdata->_id_map.size();
        }
    }
    void update_bim(bInfo* bdata,vector<int> &rsnp) {
        
        //update bim information
        vector<int> chr_buf, bp_buf;
        vector<string> a1_buf, a2_buf, ref_A_buf, other_A_buf;
        vector<string> snp_name_buf;
        vector<double> genet_dst_buf, impRsq_buf;
        for (int i = 0; i < bdata->_snp_num; i++) {
            if (!rsnp[i]) continue;
            chr_buf.push_back(bdata->_chr[i]);
            snp_name_buf.push_back(bdata->_snp_name[i]);
            genet_dst_buf.push_back(bdata->_genet_dst[i]);
            bp_buf.push_back(bdata->_bp[i]);
            a1_buf.push_back(bdata->_allele1[i]);
            a2_buf.push_back(bdata->_allele2[i]);
            ref_A_buf.push_back(bdata->_ref_A[i]);
            other_A_buf.push_back(bdata->_other_A[i]);
            if(bdata->_impRsq.size()>0) impRsq_buf.push_back(bdata->_impRsq[i]);
        }
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        bdata->_ref_A.clear();
        bdata->_other_A.clear();
        bdata->_impRsq.clear();
        bdata->_chr = chr_buf;
        bdata->_snp_name = snp_name_buf;
        bdata->_genet_dst = genet_dst_buf;
        bdata->_bp = bp_buf;
        bdata->_allele1 = a1_buf;
        bdata->_allele2 = a2_buf;
        bdata->_ref_A = ref_A_buf;
        bdata->_other_A = other_A_buf;
        bdata->_impRsq=impRsq_buf;
        bdata->_snp_num = bdata->_chr.size();
        bdata->_include.clear();
        bdata-> _include.resize(bdata->_snp_num);
        bdata->_snp_name_map.clear();
        
        for (int i = 0; i < bdata->_snp_num; i++) {
            bdata->_include[i] = i;
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }    
    void update_fam(bInfo* bdata,vector<int> &rindi) {
        //update fam information
        int i = 0;
        vector<string> fid_buf, pid_buf, fa_id_buf, mo_id_buf;
        vector<int> sex_buf;
        vector<double> pheno_buf;
        for (i = 0; i < bdata->_indi_num; i++) {
            if (!rindi[i]) continue;
            fid_buf.push_back(bdata->_fid[i]);
            pid_buf.push_back(bdata->_pid[i]);
            fa_id_buf.push_back(bdata->_fa_id[i]);
            mo_id_buf.push_back(bdata->_mo_id[i]);
            sex_buf.push_back(bdata->_sex[i]);
            pheno_buf.push_back(bdata->_pheno[i]);
        }
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        bdata->_fid = fid_buf;
        bdata->_pid = pid_buf;
        bdata->_fa_id = fa_id_buf;
        bdata->_mo_id = mo_id_buf;
        bdata->_sex = sex_buf;
        bdata->_pheno = pheno_buf;
        
        bdata->_indi_num = bdata->_fid.size();
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        for (i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
        }
    }
    void read_bimfile(bInfo* bdata,string bimfile) {
        // Read bim file: recombination rate is defined between SNP i and SNP i-1
        int ibuf = 0;
        string cbuf = "0";
        double dbuf = 0.0;
        string str_buf;
        ifstream Bim(bimfile.c_str());
        if (!Bim)
        {
            sprintf(logbuf, "Error: can not open the file %s to read.\n", bimfile.c_str());
            logprintb();
           TERMINATE();
        }
        sprintf(logbuf, "Reading PLINK BIM file from %s .\n", bimfile.c_str());
        logprintb();
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        while (Bim) {
            Bim >> ibuf;
            if (Bim.eof()) break;
            bdata->_chr.push_back(ibuf);
            Bim >> str_buf;
            bdata->_snp_name.push_back(str_buf);
            Bim >> dbuf;
            bdata->_genet_dst.push_back(dbuf);
            Bim >> ibuf;
            bdata->_bp.push_back(ibuf);
            Bim >> cbuf;
            to_upper(cbuf);
            bdata->_allele1.push_back(cbuf.c_str());
            Bim >> cbuf;
            to_upper(cbuf);
            bdata->_allele2.push_back(cbuf.c_str());
        }
        Bim.close();
        bdata->_snp_num = bdata->_chr.size();
        bdata->_ref_A = bdata->_allele1;
        bdata->_other_A = bdata->_allele2;
        sprintf(logbuf, "%llu SNPs to be included from %s .\n", bdata->_snp_num, bimfile.c_str());
        logprintb();
        
        // Initialize _include
        bdata->_include.clear();
        bdata->_include.resize( bdata->_snp_num);
        bdata->_snp_name_map.clear();
        uint64_t size = 0;
        for (int i = 0; i <  bdata->_snp_num; i++) {
            bdata->_include[i] = i;
            bdata->_snp_name_map.insert(pair<string, int>( bdata->_snp_name[i], i));
            if (size ==  bdata->_snp_name_map.size())
            {
                sprintf(logbuf, "Error: Duplicated SNP IDs found: %s .\n", bdata->_snp_name[i].c_str());
                logprintb();
               TERMINATE();
            }
            size =  bdata->_snp_name_map.size();
        }
    }
    // some code are adopted from PLINK with modifications
    void read_bedfile(bInfo* bdata, string bedfile) {
        int i = 0, j = 0, k = 0;
        
        // Flag for reading individuals and SNPs
        vector<int> rindi, rsnp;
        //get_rindi
        rindi.clear();
        rindi.resize(bdata->_indi_num);
        for (int i = 0; i < bdata->_indi_num; i++) {
            if (bdata->_id_map.find(bdata->_fid[i] + ":" + bdata->_pid[i]) != bdata->_id_map.end()) rindi[i] = 1;
            else rindi[i] = 0;
        }
        //get_rsnp
        rsnp.clear();
        rsnp.resize(bdata->_snp_num);
        for (int i = 0; i < bdata->_snp_num; i++) {
            if (bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()) rsnp[i] = 1;
            else rsnp[i] = 0;
        }
        
        if (bdata->_include.size() == 0)
        {
            sprintf(logbuf, "Error: No SNP is retained for analysis.\n");
            logprintb();
           TERMINATE();
        }
        if (bdata->_keep.size() == 0)
        {
            sprintf(logbuf, "Error: No individual is retained for analysis.\n");
            logprintb();
           TERMINATE();
        }
        
        // Read bed file
        char ch[1];
        bitset<8> b;
        bdata->_snp_1.resize(bdata->_include.size());
        bdata->_snp_2.resize(bdata->_include.size());
        for (i = 0; i < bdata->_include.size(); i++) {
            bdata->_snp_1[i].reserve(bdata->_keep.size());
            bdata->_snp_2[i].reserve(bdata->_keep.size());
        }
        fstream BIT(bedfile.c_str(), ios::in | ios::binary);
        if (!BIT) {
            sprintf(logbuf, "Error: can not open the file %s to read.\n", bedfile.c_str());
            logprintb();
           TERMINATE();
        }
        sprintf(logbuf, "Reading PLINK BED file from  %s in SNP-major format ...\n", bedfile.c_str());
        logprintb();
        for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
        int snp_indx = 0, indi_indx = 0;
        for (j = 0, snp_indx = 0; j < bdata->_snp_num; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
            if (!rsnp[j]) {
                for (i = 0; i < bdata->_indi_num; i += 4) BIT.read(ch, 1);
                continue;
            }
            for (i = 0, indi_indx = 0; i < bdata->_indi_num;) {
                BIT.read(ch, 1);
                if (!BIT)
                {
                    sprintf(logbuf, "Error: problem with the BED file ... has the FAM/BIM file been changed?\n" );
                    logprintb();
                   TERMINATE();
                }
                b = ch[0];
                k = 0;
                while (k < 7 && i < bdata->_indi_num) { // change code: 11 for AA; 00 for BB;
                    if (!rindi[i]) k += 2;
                    else {
                        bdata->_snp_2[snp_indx][indi_indx] = (!b[k++]);
                        bdata->_snp_1[snp_indx][indi_indx] = (!b[k++]);
                        indi_indx++;
                    }
                    i++;
                }
            }
            if (snp_indx == bdata->_include.size()) break;
            snp_indx++;
        }
        BIT.clear();
        BIT.close();
        sprintf(logbuf, "Genotype data for %ld  individuals and %ld SNPs to be included from %s .\n",bdata->_keep.size(), bdata->_include.size(),bedfile.c_str() );
        logprintb();
        update_fam(bdata, rindi);
        update_bim(bdata, rsnp);
    }
    void read_indi_list(string indi_list_file, vector<string> &indi_list) {
        ifstream i_indi_list(indi_list_file.c_str());
        if(!i_indi_list) {
            sprintf(logbuf, "Error: can not open the file %s to read.\n", indi_list_file.c_str());
            logprintb();
           TERMINATE();
        }
        string str_buf, id_buf;
        indi_list.clear();
        while(i_indi_list){
            i_indi_list>>str_buf;
            if(i_indi_list.eof()) break;
            id_buf=str_buf+":";
            i_indi_list>>str_buf;
            id_buf+=str_buf;
            indi_list.push_back(id_buf);
            getline(i_indi_list, str_buf);
        }
        i_indi_list.close();
    }
    void keep_indi(bInfo* bdata,string indi_list_file) {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        update_map_kp(indi_list, bdata->_id_map, bdata->_keep);
        sprintf(logbuf, "%ld individuals are kept from %s.\n", bdata->_keep.size(), indi_list_file.c_str());
        logprintb();
    }
    void remove_indi(bInfo* bdata, string indi_list_file) {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        long prev_size = bdata->_keep.size();
        update_map_rm(indi_list, bdata->_id_map, bdata->_keep);
        sprintf(logbuf, "%ld individuals are removed from %s and there are %ld individuals remaining.\n", prev_size - bdata->_keep.size(), indi_list_file.c_str(),bdata->_keep.size());
        logprintb();
    }
    void extract_snp(bInfo* bdata,string snplistfile) {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        update_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        sprintf(logbuf, "%ld SNPs are extracted from %s.\n", bdata->_include.size(), snplistfile.c_str());
        logprintb();
    }
    void exclude_snp(bInfo* bdata,string snplistfile) {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        long prev_size = bdata->_include.size();
        update_map_rm(snplist, bdata->_snp_name_map, bdata->_include);
        sprintf(logbuf, "%ld SNPs are excluded from  %s MB  and there are %ld SNPs remaining.\n", prev_size - bdata->_include.size(), snplistfile.c_str(), bdata->_include.size());
        logprintb();
        
    }
    void mu_func(bInfo* bdata, int j, vector<double> &fac) {
        int i = 0;
        bdata->_dosage_flag = 0;
        double fcount = 0.0, f_buf = 0.0;
        if (bdata->_dosage_flag) {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                    bdata->_mu[bdata->_include[j]] += fac[i] * bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    fcount += fac[i];
                }
            }
        } else {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                    f_buf = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    if (bdata->_allele2[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) f_buf = 2.0 - f_buf;
                    bdata->_mu[bdata->_include[j]] += fac[i] * f_buf;
                    fcount += fac[i];
                }
            }
        }
        
        if (fcount > 0.0)bdata->_mu[bdata->_include[j]] /= fcount;
    }
    void calcu_mu(bInfo* bdata, bool ssq_flag) {
        int i = 0, j = 0;
        
        vector<double> auto_fac(bdata->_keep.size()), xfac(bdata->_keep.size()), fac(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++)
        {
            auto_fac[i] = 1.0;
            if (bdata->_sex[bdata->_keep[i]] == 1) xfac[i] = 0.5;
            else if (bdata->_sex[bdata->_keep[i]] == 2) xfac[i] = 1.0;
            fac[i] = 0.5;
        }
        sprintf(logbuf, "Calculating allele frequencies ...\n");
        logprintb();
        bdata->_mu.clear();
        bdata->_mu.resize(bdata->_snp_num);
        
        #pragma omp parallel for
        for (j = 0; j < bdata->_include.size(); j++)
        {
            if (bdata->_chr[bdata->_include[j]]<(bdata->_autosome_num + 1)) mu_func(bdata, j, auto_fac);
            else if (bdata->_chr[bdata->_include[j]] == (bdata->_autosome_num + 1)) mu_func(bdata,j, xfac);
            else mu_func(bdata, j, fac);
        }
    }
    void filter_snp_maf(bInfo* bdata,double maf) {
        if(bdata->_mu.empty()) calcu_mu(bdata);
        sprintf(logbuf, "Pruning SNPs with MAF > %lf ...\n", maf);
        logprintb();
        map<string, int> id_map_buf(bdata->_snp_name_map);
        map<string, int>::iterator iter, end=id_map_buf.end();
        long prev_size=bdata->_include.size();
        double fbuf=0.0;
        bdata->_include.clear();
        bdata->_snp_name_map.clear();
        for(iter=id_map_buf.begin(); iter!=end; iter++){
            fbuf=bdata->_mu[iter->second]*0.5;
            if(fbuf<=maf || (1.0-fbuf)<=maf) continue;
            bdata->_snp_name_map.insert(*iter);
            bdata->_include.push_back(iter->second);
        }
        if(bdata->_include.size()==0)
        {
            sprintf(logbuf, "Error: No SNP is retained for analysis.\n");
            logprintb();
           TERMINATE();
        }
        else{
            stable_sort(bdata->_include.begin(), bdata->_include.end());
            sprintf(logbuf, "After pruning SNPs with MAF > %lf, there are %ld SNPs (%ld SNPs with MAF < %lf).\n",maf, bdata->_include.size(), prev_size-bdata->_include.size(),maf);
            logprintb();
        }
    }
    bool make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool mu) {
        // Eigen is column-major by default. here row of X is individual, column of X is SNP.
        if (mu && bdata->_mu.empty()) calcu_mu(bdata);
        bool have_mis = false;
        long snpNum=snpids.size();
        X.resize(bdata->_keep.size(),snpNum);
        #pragma omp parallel for
        for (int i = 0; i < snpNum ; i++)
        {
            uint32_t snpid=snpids[i];
            for (int j = 0; j < bdata->_keep.size() ; j++)
            {
                
                if (!bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] || bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]])
                {
                    if (bdata->_allele1[bdata->_include[snpid]] == bdata->_ref_A[bdata->_include[snpid]]) X(j,i)= bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]];
                    else X(j,i)= 2.0 - (bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]]);
                } else {
                    if(mu) X(j,i) = bdata->_mu[bdata->_include[snpid]];
                    else X(j,i)=1e6;
                    have_mis = true;
                }
            }
        }
        return have_mis;
    }
    bool make_XMat(bInfo* bdata, int start, int slide_wind, MatrixXd &X, bool mu)
    {
        if (mu && bdata->_mu.empty()) calcu_mu(bdata);
        
        bool have_mis = false;
        long  snpNum = 0;
        if(bdata->_include.size()-start>slide_wind) snpNum=slide_wind;
        else snpNum=bdata->_include.size()-start;
        X.resize(bdata->_keep.size(), snpNum);
        #pragma omp parallel for
        for (int i = 0; i < snpNum ; i++)
        {
            int idx=start+i;
            for (int j = 0; j < bdata->_keep.size() ; j++)
            {
                
                if (!bdata->_snp_1[bdata->_include[idx]][bdata->_keep[j]] || bdata->_snp_2[bdata->_include[idx]][bdata->_keep[j]])
                {
                    if (bdata->_allele1[bdata->_include[idx]] == bdata->_ref_A[bdata->_include[idx]]) X(j,i)= bdata->_snp_1[bdata->_include[idx]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[idx]][bdata->_keep[j]];
                    else X(j,i)= 2.0 - (bdata->_snp_1[bdata->_include[idx]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[idx]][bdata->_keep[j]]);
                } else {
                    if(mu) X(j,i) = bdata->_mu[bdata->_include[idx]];
                    else X(j,i)=1e6;
                    have_mis = true;
                }
            }
        }
        //cout<<"["<<start<<","<<(start+snpNum-1)<<"]"<<":"<<(have_mis?"havemiss":"nomiss")<<":"<<X.cols()<<":"<<X.rows()<<endl;
        return have_mis;
    }
    bool make_XMat_subset(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool standardise) {
        // Eigen is column-major by default. here row of X is individual, column of X is SNP.
        if (bdata->_mu.empty()) calcu_mu(bdata);
        bool have_mis = false;
        long snpNum=snpids.size();
        X.resize(bdata->_keep.size(),snpNum);
        #pragma omp parallel for
        for (int i = 0; i < snpNum ; i++)
        {
            uint32_t snpid=snpids[i];
            for (int j = 0; j < bdata->_keep.size() ; j++)
            {
                
                if (!bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] || bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]])
                {
                    if (bdata->_allele1[bdata->_include[snpid]] == bdata->_ref_A[bdata->_include[snpid]]) X(j,i)= bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]];
                    else X(j,i)= 2.0 - (bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]]);
                    X(j,i) -= bdata->_mu[bdata->_include[snpid]];
                } else {
                     X(j,i)=0;
                    have_mis = true;
                }
            }
        }
        if(standardise)
        {
            vector<double> sd_SNP(snpNum);
            for (int j = 0; j < snpNum; j++){
                uint32_t snpid=snpids[j];
                sd_SNP[j] = bdata->_mu[bdata->_include[snpid]]*(1.0 - 0.5 * bdata->_mu[bdata->_include[snpid]]);
                if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
                else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
            }
            for (int j = 0; j < snpNum; j++) X.col(j) = X.col(j).array() * sd_SNP[j];
        }
        return have_mis;
    }
    bool make_XMat_subset(bInfo* bdata, int start, int slide_wind, MatrixXd &X, bool standardise)
    {
        if ( bdata->_mu.empty()) calcu_mu(bdata);
        
        bool have_mis = false;
        long  snpNum = 0;
        if(bdata->_include.size()-start>slide_wind) snpNum=slide_wind;
        else snpNum=bdata->_include.size()-start;
        X.resize(bdata->_keep.size(), snpNum);
        #pragma omp parallel for
        for (int i = 0; i < snpNum ; i++)
        {
            int idx=start+i;
            for (int j = 0; j < bdata->_keep.size() ; j++)
            {
                
                if (!bdata->_snp_1[bdata->_include[idx]][bdata->_keep[j]] || bdata->_snp_2[bdata->_include[idx]][bdata->_keep[j]])
                {
                    if (bdata->_allele1[bdata->_include[idx]] == bdata->_ref_A[bdata->_include[idx]]) X(j,i)= bdata->_snp_1[bdata->_include[idx]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[idx]][bdata->_keep[j]];
                    else X(j,i)= 2.0 - (bdata->_snp_1[bdata->_include[idx]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[idx]][bdata->_keep[j]]);
                    X(j,i) -= bdata->_mu[bdata->_include[idx]];
                } else {
                    X(j,i)=0;
                    have_mis = true;
                }
            }
        }
        if(standardise)
        {
            vector<double> sd_SNP(snpNum);
            for (int j = 0; j < snpNum; j++){
                int idx=start+j;
                sd_SNP[j] = bdata->_mu[bdata->_include[idx]]*(1.0 - 0.5 * bdata->_mu[bdata->_include[idx]]);
                if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
                else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
            }
            for (int j = 0; j < snpNum; j++) X.col(j) = X.col(j).array() * sd_SNP[j];
        }
        return have_mis;
    }
    void check_autosome(bInfo* bdata) {
        for (int i = 0; i <bdata->_include.size(); i++) {
            if (bdata->_chr[bdata->_include[i]] > bdata->_autosome_num) throw ("Error: this option is for the autosomal SNPs only. Please check the option --autosome.");
        }
    }
    bool make_XMat(bInfo* bdata, MatrixXf &X)
    {
        if (bdata->_mu.empty()) calcu_mu(bdata);
        
        cout << "Recoding genotypes (individual major mode) ..." << endl;
        bool have_mis = false;
        unsigned long n = bdata->_keep.size(), m = bdata->_include.size();
        
        X.resize(0,0);
        X.resize(n, m);
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < bdata->_include.size(); j++) {
                    if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                        if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) X(i,j) = bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]];
                        else X(i,j) = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    }
                    else {
                        X(i,j) = 1e6;
                        have_mis = true;
                    }
                }
            
        }
        return have_mis;
    }
    void std_XMat(bInfo* bdata, MatrixXf &X, VectorXd &sd_SNP, bool miss_with_mu, bool divid_by_std)
    {
        if (bdata->_mu.empty()) calcu_mu(bdata);
        
        unsigned long n = bdata->_keep.size(), m = bdata->_include.size();
        sd_SNP.resize(m);
        
            for (int j = 0; j < m; j++) sd_SNP[j] = bdata->_mu[bdata->_include[j]]*(1.0 - 0.5 * bdata->_mu[bdata->_include[j]]);
        
        if (divid_by_std) {
            for (int j = 0; j < m; j++) {
                if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
                else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
            }
        }
        
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (X(i,j) < 1e5) {
                    X(i,j) -= bdata->_mu[bdata->_include[j]];
                    if (divid_by_std) X(i,j) *= sd_SNP[j];
                }
                else if (miss_with_mu) X(i,j) = 0.0;
            }
        }
    }
    void make_grm(bInfo* bdata, int grm_mtd,  bool diag_f3_flag)
    {
        bool have_mis = false;
        check_autosome(bdata);
        unsigned long  n = bdata->_keep.size(), m = bdata->_include.size();
        have_mis = make_XMat(bdata,bdata->_geno);
        VectorXd sd_SNP;
        if (grm_mtd == 0)
            std_XMat(bdata, bdata->_geno, sd_SNP, false, true);
        else
            std_XMat(bdata, bdata->_geno, sd_SNP, false, false);
         cout << "\nCalculating the genetic relationship matrix (GRM) ... " << endl;
        
        // count the number of missing genotypes
        vector< vector<int> > miss_pos;
        vector< vector<bool> > X_bool;
        if(have_mis){
            miss_pos.resize(n);
            X_bool.resize(n);
#pragma omp parallel for
            for (int i = 0; i < n; i++) {
                X_bool[i].resize(m);
                for (int j = 0; j < m; j++) {
                    if (bdata->_geno(i,j) < 1e5) X_bool[i][j] = true;
                    else {
                        bdata->_geno(i,j) = 0.0;
                        miss_pos[i].push_back(j);
                        X_bool[i][j] = false;
                    }
                }
            }
        }
        
        // Calculate A_N matrix
        if(have_mis){
            bdata->_grm_N.resize(n, n);
#pragma omp parallel for
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= i; j++) {
                    int miss_j = 0;
                    for (int k = 0; k < miss_pos[j].size(); k++) miss_j += (int) X_bool[i][miss_pos[j][k]];
                    bdata->_grm_N(i,j) = m - miss_pos[i].size() - miss_j;
                    bdata->_grm_N(j,i) = bdata->_grm_N(i,j);
                }
            }
        }
        else bdata->_grm_N = MatrixXf::Constant(n,n,m);
        
        // Calcuate WW'
#ifdef SINGLE_PRECISION
        bdata->_grm = _geno * _geno.transpose();
#else
        bdata->_grm = (bdata->_geno * bdata->_geno.transpose()).cast<double>();
#endif
        
        // Calculate A matrix
        if (grm_mtd == 1) bdata->_grm_N = bdata->_grm_N.array() * sd_SNP.mean();
        
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                if(bdata->_grm_N(i,j) > 0) bdata->_grm(i,j) /= bdata->_grm_N(i,j);
                else bdata->_grm(i,j) = 0.0;
                bdata->_grm(j,i) = bdata->_grm(i,j);
            }
        }
        
        
        // re-calcuate the diagonals (Fhat3+1)
        if (diag_f3_flag) {
#pragma omp parallel for
            for (int i = 0; i < n; i++) {
                bdata->_grm(i,i) = 0.0;
                double non_missing = 0.0;
                for (int j = 0; j < m; j++) {
                    if (bdata->_geno(i,j) < 1e5){
                        bdata->_grm(i,i) += bdata->_geno(i,j)*(bdata->_geno(i,j)+(bdata->_mu[bdata->_include[j]] - 1.0) * sd_SNP[j]);
                        non_missing += 1.0;
                    }
                }
                bdata->_grm(i,i) /= non_missing;
            }
        }
        
        
        if ( grm_mtd == 0) {
            for (int j = 0; j < m; j++) {
                if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
                else sd_SNP[j] = 1.0 / sd_SNP[j];
            }
#pragma omp parallel for
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    if (bdata->_geno(i,j) < 1e5) bdata->_geno(i,j) *= sd_SNP[j];
                    else bdata->_geno(i,j) = 0.0;
                }
            }
        }
    }
    void init_keep(bInfo* bdata) {
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        long size = 0;
        for (int i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
            if (size == bdata->_id_map.size()) throw ("Error: Duplicate individual ID found: \"" + bdata->_fid[i] + "\t" +bdata-> _pid[i] + "\".");
            size = bdata->_id_map.size();
        }
    }
    int read_grm_id(bInfo* bdata,string grm_file, vector<string> &grm_id)
    {
        // read GRM IDs
        string grm_id_file = grm_file + ".grm.id";
        ifstream i_grm_id(grm_id_file.c_str());
        if (!i_grm_id) throw ("Error: can not open the file [" + grm_id_file + "] to read.");
        string str_buf, id_buf;
        vector<string> fid, pid;
        grm_id.clear();
        while (i_grm_id) {
            i_grm_id >> str_buf;
            if (i_grm_id.eof()) break;
            fid.push_back(str_buf);
            id_buf = str_buf + ":";
            i_grm_id >> str_buf;
            pid.push_back(str_buf);
            id_buf += str_buf;
            grm_id.push_back(id_buf);
            getline(i_grm_id, str_buf);
        }
        i_grm_id.close();
        int n = (int)grm_id.size();
        
        if (bdata->_id_map.empty()) {
            bdata->_fid = fid;
            bdata->_pid = pid;
            bdata->_indi_num = bdata->_fid.size();
            bdata->_sex.resize(bdata->_fid.size());
            init_keep(bdata);
        }else{
            update_map_kp(grm_id, bdata->_id_map, bdata->_keep);
        }
        
        return (n);
    }
    
    void read_grm_gz(bInfo* bdata,string grm_file, vector<string> &grm_id, bool read_id_only) {
        int n = read_grm_id(bdata, grm_file, grm_id);
        
        if (read_id_only) return;
        
        string grm_gzfile = grm_file + ".grm.gz", str_buf;
        const int MAX_LINE_LENGTH = 1000;
        char buf[MAX_LINE_LENGTH];
        gzFile zinf= gzopen(grm_gzfile.c_str(), "rb");
        
        if (!(zinf)) {
            LOGPRINTF("ERROR: Couldn't open file %s\n", grm_gzfile.c_str());
            TERMINATE();
        }
        
        int indx1 = 0, indx2 = 0, nline = 0;
        double grm_buf = 0.0, grm_N_buf;
        string errmsg = "Error: failed to read [" + grm_gzfile + "]. The format of the GRM file has been changed?\nError occurs in line:\n";
        cout << "Reading the GRM from [" + grm_gzfile + "]." << endl;
        bdata->_grm.resize(n, n);
        bdata->_grm_N.resize(n, n);
        while(!gzeof(zinf)) {
            gzgets(zinf, buf, MAX_LINE_SIZE);
            if(buf[0]!='\0')
            {
            stringstream ss(buf);
            if (!(ss >> indx1)) throw (errmsg + buf);
            if (!(ss >> indx2)) throw (errmsg + buf);
            if (!(ss >> grm_N_buf)) throw (errmsg + buf);
            if (!(ss >> grm_buf)) throw (errmsg + buf);
            if (indx1 < indx2 || indx1 > n || indx2 > n) throw (errmsg + buf);
            if (grm_N_buf == 0) cout << "Warning: " << buf << endl;
            bdata->_grm_N(indx1 - 1, indx2 - 1) = bdata->_grm_N(indx2 - 1, indx1 - 1) = grm_N_buf;
            bdata->_grm(indx1 - 1, indx2 - 1) = bdata->_grm(indx2 - 1, indx1 - 1) = grm_buf;
            nline++;
            if (ss >> str_buf) throw (errmsg + buf);
            }
        }
         gzclose(zinf);
       
        cout << "GRM for " << n << " individuals are included from [" + grm_gzfile + "]." << endl;
    }
    
    void read_grm_bin(bInfo* bdata,string grm_file, vector<string> &grm_id, bool read_id_only, bool dont_read_N)
    {
        int i = 0, j = 0, n = read_grm_id(bdata,grm_file, grm_id);
        
        if (read_id_only) return;
        
        string grm_binfile = grm_file + ".grm.bin";
        ifstream A_bin(grm_binfile.c_str(), ios::in | ios::binary);
        if (!A_bin.is_open()) throw ("Error: can not open the file [" + grm_binfile + "] to read.");
        bdata->_grm.resize(n, n);
        cout << "Reading the GRM from [" + grm_binfile + "]." << endl;
        int size = sizeof (float);
        float f_buf = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                if (!(A_bin.read((char*) &f_buf, size))) throw ("Error: the size of the [" + grm_binfile + "] file is incomplete?");
                bdata->_grm(j, i) = bdata->_grm(i, j) = f_buf;
            }
        }
        A_bin.close();
        
        if(!dont_read_N){
            string grm_Nfile = grm_file + ".grm.N.bin";
            ifstream N_bin(grm_Nfile.c_str(), ios::in | ios::binary);
            if (!N_bin.is_open()) throw ("Error: can not open the file [" + grm_Nfile + "] to read.");
            bdata->_grm_N.resize(n, n);
            cout << "Reading the number of SNPs for the GRM from [" + grm_Nfile + "]." << endl;
            size = sizeof (float);
            f_buf = 0.0;
            for (i = 0; i < n; i++) {
                for (j = 0; j <= i; j++) {
                    if (!(N_bin.read((char*) &f_buf, size))) throw ("Error: the size of the [" + grm_Nfile + "] file is incomplete?");
                    bdata->_grm_N(j, i) = bdata->_grm_N(i, j) = f_buf;
                }
            }
            N_bin.close();
        }
        
        cout << "GRM for " << n << " individuals are included from [" + grm_binfile + "]." << endl;
    }
    
    void read_grm(bInfo* bdata,bool grm_bin_flag, string grm_file, vector<string> &grm_id, bool read_id_only, bool dont_read_N)
    {
        if (grm_bin_flag) read_grm_bin(bdata,grm_file, grm_id, read_id_only, dont_read_N);
        else read_grm_gz(bdata,grm_file, grm_id, read_id_only);
    }
}

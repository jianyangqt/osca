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
        vector<char> a1_buf, a2_buf, ref_A_buf, other_A_buf;
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
            bdata->_allele1.push_back(cbuf.c_str()[0]);
            Bim >> cbuf;
            to_upper(cbuf);
            bdata->_allele2.push_back(cbuf.c_str()[0]);
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
    void make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool minus_2p) {
        // Eigen is column-major by default. here row of X is individual, column of X is SNP.
        if (minus_2p && bdata->_mu.empty()) calcu_mu(bdata);
        uint64_t snpNum=snpids.size();
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
                } else X(j,i) = bdata->_mu[bdata->_include[snpid]];
                if (minus_2p) X(j,i) -= bdata->_mu[bdata->_include[snpid]];
                
            }
        }
    }
    bool make_XMat(bInfo* bdata, int start, int slide_wind, MatrixXf &X, bool mu)
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
    
}

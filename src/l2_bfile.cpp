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
                sprintf(logbuf, "Error: Duplicated SNP IDs found: %s .(Untill now, sex chromsomes is not supported)\n", bdata->_snp_name[i].c_str());
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
    void extract_snp(bInfo* bdata, int chr)
    {
        vector<string> snplist;
        for(int i = 0; i < bdata->_include.size(); i++){
            int j = bdata->_include[i];
            if(bdata->_chr[j] == chr) snplist.push_back(bdata->_snp_name[j]);
        }
        if(snplist.empty()) throw ("Error: on SNP found in this region.");
        update_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout << bdata->_include.size() << " SNPs are extracted from chromosome "<<chr<< "."<<endl;
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
        double fcount = 0.0, f_buf = 0.0, call = 0.0;
        if (bdata->_dosage_flag) {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                    bdata->_mu[bdata->_include[j]] += fac[i] * bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    fcount += fac[i];
                    call += 1;
                }
            }
        } else {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                    f_buf = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    if (bdata->_allele2[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) f_buf = 2.0 - f_buf;
                    bdata->_mu[bdata->_include[j]] += fac[i] * f_buf;
                    fcount += fac[i];
                    call += 1;
                }
            }
        }
        
        if (fcount > 0.0) bdata->_mu[bdata->_include[j]] /= fcount;
        bdata->_mr[bdata->_include[j]] = call/bdata->_keep.size();
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
        //sprintf(logbuf, "Calculating allele frequencies ...\n");
        //logprintb();
        bdata->_mu.clear();
        bdata->_mr.clear();
        bdata->_mu.resize(bdata->_snp_num);
        bdata->_mr.resize(bdata->_snp_num);
        
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
        //sprintf(logbuf, "Pruning SNPs with MAF > %lf ...\n", maf);
        //logprintb();
        map<string, int> id_map_buf(bdata->_snp_name_map);
        map<string, int>::iterator iter, end=id_map_buf.end();
        long prev_size=bdata->_include.size();
        double fbuf=0.0;
        bdata->_include.clear();
        bdata->_snp_name_map.clear();
        for(iter=id_map_buf.begin(); iter!=end; iter++){
            fbuf=bdata->_mu[iter->second]*0.5;
            if(fbuf<maf || (1.0-fbuf)<maf) continue;
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
            //sprintf(logbuf, "After pruning SNPs with MAF >= %lf, there are %ld SNPs (%ld SNPs with MAF =< %lf).\n",maf, bdata->_include.size(), prev_size-bdata->_include.size(),maf);
            //logprintb();
            if(prev_size > bdata->_include.size()) {LOGPRINTF("%ld SNPs are excluded due to MAF < %4.2f. \n", prev_size-bdata->_include.size(),maf);}
        }
    }
    void filter_snp_call(bInfo* bdata,double call)
    {
        if(bdata->_mr.empty()) calcu_mu(bdata);
        //sprintf(logbuf, "Pruning SNPs with MAF > %lf ...\n", maf);
        //logprintb();
        map<string, int> id_map_buf(bdata->_snp_name_map);
        map<string, int>::iterator iter, end=id_map_buf.end();
        long prev_size=bdata->_include.size();
        double fbuf=0.0;
        bdata->_include.clear();
        bdata->_snp_name_map.clear();
        for(iter=id_map_buf.begin(); iter!=end; iter++){
            fbuf=bdata->_mr[iter->second];
            if(fbuf<call) continue;
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
            //sprintf(logbuf, "After pruning SNPs with MAF >= %lf, there are %ld SNPs (%ld SNPs with MAF =< %lf).\n",maf, bdata->_include.size(), prev_size-bdata->_include.size(),maf);
            //logprintb();
            if(prev_size > bdata->_include.size()) {LOGPRINTF("%ld SNPs are excluded due to call < %4.2f. \n", prev_size-bdata->_include.size(),call);}
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
    bool make_Xvec_subset(bInfo* bdata,uint32_t &includeid, VectorXd &vec, bool standardise)
    {
        if (bdata->_mu.empty()) calcu_mu(bdata);
        bool have_mis = false;
        vec.resize(bdata->_keep.size());
        int snpid =bdata->_include[includeid];
        for (int j = 0; j < bdata->_keep.size() ; j++)
        {
                
            if (!bdata->_snp_1[snpid][bdata->_keep[j]] || bdata->_snp_2[snpid][bdata->_keep[j]])
            {
                if (bdata->_allele1[snpid] == bdata->_ref_A[snpid]) vec(j)= bdata->_snp_1[snpid][bdata->_keep[j]] + bdata->_snp_2[snpid][bdata->_keep[j]];
                else vec(j)= 2.0 - (bdata->_snp_1[snpid][bdata->_keep[j]] + bdata->_snp_2[snpid][bdata->_keep[j]]);
                vec(j) -= bdata->_mu[snpid];
            } else {
                vec(j)=0;
                have_mis = true;
            }
        }
        if(standardise)
        {
            
                double sd_SNP = bdata->_mu[snpid]*(1.0 - 0.5 * bdata->_mu[snpid]);
                if (fabs(sd_SNP) < 1.0e-50) sd_SNP = 0.0;
                else sd_SNP = sqrt(1.0 / sd_SNP);
            
             vec = vec * sd_SNP;
        }
        return have_mis;
    }
    bool make_Xvec_subset(bInfo* bdata,int &snpid, VectorXd &vec, bool standardise)
    {
        if (bdata->_mu.empty()) calcu_mu(bdata);
        bool have_mis = false;
        vec.resize(bdata->_keep.size());

        for (int j = 0; j < bdata->_keep.size() ; j++)
        {
            
            if (!bdata->_snp_1[snpid][bdata->_keep[j]] || bdata->_snp_2[snpid][bdata->_keep[j]])
            {
                if (bdata->_allele1[snpid] == bdata->_ref_A[snpid]) vec(j)= bdata->_snp_1[snpid][bdata->_keep[j]] + bdata->_snp_2[snpid][bdata->_keep[j]];
                else vec(j)= 2.0 - (bdata->_snp_1[snpid][bdata->_keep[j]] + bdata->_snp_2[snpid][bdata->_keep[j]]);
                vec(j) -= bdata->_mu[snpid];
            } else {
                vec(j)=0;
                have_mis = true;
            }
        }
        if(standardise)
        {
            
            double sd_SNP = bdata->_mu[snpid]*(1.0 - 0.5 * bdata->_mu[snpid]);
            if (fabs(sd_SNP) < 1.0e-50) sd_SNP = 0.0;
            else sd_SNP = sqrt(1.0 / sd_SNP);
            
            vec = vec * sd_SNP;
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
    bool make_XMat(bInfo* bdata, MatrixXd &X)
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
    void std_XMat(bInfo* bdata, MatrixXd &X, VectorXd &sd_SNP, bool miss_with_mu, bool divid_by_std)
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
    void std_XMat_(bInfo* bdata, MatrixXd &X,  VectorXd &sd_SNP, bool miss_with_mu, bool divid_by_std)
    {
        unsigned long n = bdata->_keep.size(), m = bdata->_include.size();
        sd_SNP.resize(m);
        for(int i=0;i<m;i++)
        {
            double nonmiss=0.0, mu=0.0,sd=0.0;
            for(int j=0; j<n; j++)
            {
                double val=X(j,i);
                if(val<1e5){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            if(nonmiss>0)
            {
                mu/=nonmiss;
                for(int j=0; j<n; j++)
                {
                    double val=X(j,i);
                    if(val<1e5)
                        X(j,i)=val-mu;
                    else if(miss_with_mu)
                        X(j,i) = 0.0;
                }
            }
            if(divid_by_std)
            {
                if(nonmiss>1)
                {
                    for(int j=0; j<n; j++)
                    {
                        double val=X(j,i);
                        if(val<1e5) sd+=val*val;
                    }
                    sd=sqrt(sd/(nonmiss-1.0));
                    if(sd>1e-30) sd_SNP[i]=1/sd;
                    else sd_SNP[i]= 0.0;
                    if(sd>1e-30)
                    {
                        for(int j=0; j<n; j++)
                        {
                            double val=X(j,i);
                            if(val<1e5) X(j,i)=val/sd;
                        }
                    }
                }
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
            //std_XMat_(bdata, bdata->_geno, sd_SNP, false, true); // for test, with EWAS
        else
            std_XMat(bdata, bdata->_geno, sd_SNP, false, false);
         cout << "Calculating the genetic relationship matrix (GRM) ... " << endl;
        
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
        bdata->_grm = (bdata->_geno * bdata->_geno.transpose());
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
        
        /*
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
         */
    }
    void make_grm(bInfo* bdata, MatrixXd &VZ, int grm_mtd,  bool diag_f3_flag)
    {
        bool have_mis = false;
        check_autosome(bdata);
        unsigned long  n = bdata->_keep.size(), m = bdata->_include.size();
        have_mis = make_XMat(bdata,bdata->_geno);
        VectorXd sd_SNP;
        if (grm_mtd == 0)
            std_XMat(bdata, bdata->_geno, sd_SNP, false, true);
        //std_XMat_(bdata, bdata->_geno, sd_SNP, false, true); // for test, with EWAS
        else
            std_XMat(bdata, bdata->_geno, sd_SNP, false, false);
        cout << "Calculating the genetic relationship matrix (GRM) ... " << endl;
        
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
        bdata->_grm = (bdata->_geno * bdata->_geno.transpose());
        VZ = bdata->_grm;
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
    
    void read_phen(bInfo* bdata, string phen_file, char* mpheno, bool mvFlg)
    {
        
        if(bdata->_keep.size()==0 )
        {
            LOGPRINTF("Error: no individual found in genotype data.\n");
            TERMINATE();
        }
        // Read phenotype data
        ifstream in_phen(phen_file.c_str());
        if (!in_phen)
        {
            LOGPRINTF("Error: can not open the file %s to read.\n",phen_file.c_str());
            TERMINATE();
        }
        
        vector<string>  vs_buf, indi_list;
        string str_buf, fid_buf, iid_buf;
        vector<int> pheno_ids;
        
        LOGPRINTF("Reading phenotypes from %s ...\n",phen_file.c_str());
        getline(in_phen, str_buf);
        
        int phen_num = split_string(str_buf, vs_buf) - 2;
        if (phen_num <= 0)
        {
            LOGPRINTF("Error: no phenotype data is found.\n");
            TERMINATE();
        }
        if(mpheno==NULL)
        {
            if(mvFlg) for(int i=0;i<phen_num;i++) pheno_ids.push_back(i);
            else pheno_ids.push_back(0);
        }else {
            int tmpint=split_string(mpheno, vs_buf, ",");
            map<int,int> tmpidsmap;
            long mapsize=0;
            for(int i=0;i<tmpint;i++)
            {
                int curid=atoi(vs_buf[i].c_str())-1;
                if(curid<0)
                {
                    LOGPRINTF("Error:  Trait number specified by --mpheno should be >=1.\n");
                    TERMINATE();
                }
                if(curid>=phen_num)
                {
                    LOGPRINTF("Error:  can not find the %dth trait in the file %s.\n",curid+1,phen_file.c_str());
                    TERMINATE();
                }
                tmpidsmap.insert(pair<int,int>(curid,i));
                if(tmpidsmap.size()==mapsize)
                {
                    LOGPRINTF("Warning: duplicate phenotype id specified by --mpheno.\n");
                } else {
                    
                    pheno_ids.push_back(curid);
                    mapsize=tmpidsmap.size();
                }
            }
            if(pheno_ids.size()>1 && !mvFlg)
            {
                LOGPRINTF("Error: multiple traits are specified by --mpheno in unitrait analysis.\n");
                TERMINATE();
            }
        }
        
        bdata->_pheno_num=(int)pheno_ids.size();
        if(!sudoAllocReserved( (pheno_ids.size()-1)*bdata->_indi_num*sizeof(double), "phenotype")) {TERMINATE();}
        bdata->_pheno.resize(pheno_ids.size()*bdata->_indi_num);
        for(int i=0;i<pheno_ids.size()*bdata->_indi_num;i++) bdata->_pheno[i]=MISSING_PHENO;
        in_phen.seekg(ios::beg);
        int line = 1;
        map<string,int> missindi;
        map<string, int>::iterator iter;
        while (in_phen) {
            line++;
            in_phen >> fid_buf;
            if (in_phen.eof()) break;
            in_phen >> iid_buf;
            getline(in_phen, str_buf);
            if (split_string(str_buf, vs_buf) != phen_num) {
                LOGPRINTF("Error: %ld phenotype values are missing in line #%d in the file %s.\n",vs_buf.size() - phen_num , line, phen_file.c_str());
                TERMINATE();
            }
            iter = bdata->_id_map.find(fid_buf + ":" + iid_buf);
            if (iter != bdata->_id_map.end())
            {
                int idx = iter->second;
                for(int i=0;i<pheno_ids.size();i++)
                {
                    //if (vs_buf[pheno_ids[i]] != "-9" && vs_buf[pheno_ids[i]] != "NA")
                    if (vs_buf[pheno_ids[i]] != "NA" && vs_buf[pheno_ids[i]] != "na")
                    {
                        bdata->_pheno[i*bdata->_indi_num+idx]=atof(vs_buf[pheno_ids[i]].c_str());
                    }
                }
            }
        }
        in_phen.close();
        
        for(int i=0;i<pheno_ids.size();i++)
            for(int j=0;j<bdata->_indi_num;j++)
                if(bdata->_pheno[i*bdata->_indi_num+j]==MISSING_PHENO) missindi.insert(pair<string,int>(bdata->_fid[j] + ":" + bdata->_pid[j],j));
        
        for (iter = missindi.begin(); iter != missindi.end(); iter++) indi_list.push_back(iter->first);
        update_map_rm(indi_list, bdata->_id_map, bdata->_keep);
        if(bdata->_keep.size()==0)
        {
            LOGPRINTF("Error: no individual is included.\n");
            TERMINATE();
        }
        LOGPRINTF("Non-missing phenotypes of %ld individuals are included from %s. \n",bdata->_keep.size(), phen_file.c_str());
    }
    
    void read_cov(bInfo* bdata, string cov_file, bool qcovFlg)
    {
        
        if(bdata->_keep.size()==0 )
        {
            LOGPRINTF("Error: no individual found in genotype data.\n");
            TERMINATE();
        }
        // Read cov data
        ifstream in_cov(cov_file.c_str());
        if (!in_cov)
        {
            LOGPRINTF("Error: can not open the file %s to read.\n",cov_file.c_str());
            TERMINATE();
        }
        
        vector<string>  vs_buf, indi_list;
        string str_buf, fid_buf, iid_buf, covtypestr;
        
        if(qcovFlg) covtypestr="quantitative";
        else covtypestr="discrete";
        
        LOGPRINTF("Reading %s covariate(s) from %s ...\n",covtypestr.c_str(),cov_file.c_str());
        getline(in_cov, str_buf);
        
        int cov_num = split_string(str_buf, vs_buf) - 2;
        if (cov_num <= 0)
        {
            LOGPRINTF("Error: no covariate data is found.\n");
            TERMINATE();
        }
        
        if(qcovFlg)
        {
            bdata->_qcov_num=cov_num;
            if(!sudoAllocReserved( cov_num*bdata->_indi_num*sizeof(double), covtypestr)) TERMINATE();
            bdata->_qcov.resize(cov_num*bdata->_indi_num);
            LOGPRINTF("There are %d %s covariate(s) in the file %s. \n",cov_num,covtypestr.c_str(), cov_file.c_str());
            for(int i=0;i<cov_num*bdata->_indi_num;i++) bdata->_qcov[i]=-9;
            in_cov.seekg(ios::beg);
            int line = 1;
            map<string,int> missindi;
            map<string, int>::iterator iter;
            while (in_cov) {
                line++;
                in_cov >> fid_buf;
                if (in_cov.eof()) break;
                in_cov >> iid_buf;
                getline(in_cov, str_buf);
                if (split_string(str_buf, vs_buf) != cov_num) {
                    LOGPRINTF("Error: %ld %s covariate values are missing in line #%d in the file %s.\n",vs_buf.size() - cov_num , covtypestr.c_str(),line, cov_file.c_str());
                    TERMINATE();
                }
                iter = bdata->_id_map.find(fid_buf + ":" + iid_buf);
                if (iter != bdata->_id_map.end())
                {
                    int idx = iter->second;
                    for(int i=0;i<cov_num;i++)
                    {
                        if (vs_buf[i] != "-9" && vs_buf[i] != "NA")
                        {
                            bdata->_qcov[i*bdata->_indi_num+idx]=atof(vs_buf[i].c_str());
                        }
                    }
                }
            }
            for(int i=0;i<cov_num;i++)
                for(int j=0;j<bdata->_indi_num;j++)
                    if(bdata->_qcov[i*bdata->_indi_num+j]==-9) missindi.insert(pair<string,int>(bdata->_fid[j] + ":" + bdata->_pid[j],j));
            for (iter = missindi.begin(); iter != missindi.end(); iter++) indi_list.push_back(iter->first);
            update_map_rm(indi_list, bdata->_id_map, bdata->_keep);
        }
        else
        {
            bdata->_cov_num=cov_num;
            if(!sudoAllocReserved( cov_num*bdata->_indi_num*sizeof(double), covtypestr)) TERMINATE();
            bdata->_cov.resize(cov_num*bdata->_indi_num);
            LOGPRINTF("There are %d %s covariate(s) in the file %s. \n",cov_num, covtypestr.c_str(),cov_file.c_str());
            for(int i=0;i<cov_num*bdata->_indi_num;i++) bdata->_cov[i]="-9";
            in_cov.seekg(ios::beg);
            int line = 1;
            map<string,int> missindi;
            map<string, int>::iterator iter;
            while (in_cov) {
                line++;
                in_cov >> fid_buf;
                if (in_cov.eof()) break;
                in_cov >> iid_buf;
                getline(in_cov, str_buf);
                if (split_string(str_buf, vs_buf) != cov_num) {
                    LOGPRINTF("Error: %ld %s covariate values are missing in line #%d in the file %s.\n",vs_buf.size() - cov_num , covtypestr.c_str(),line, cov_file.c_str());
                    TERMINATE();
                }
                iter = bdata->_id_map.find(fid_buf + ":" + iid_buf);
                if (iter != bdata->_id_map.end())
                {
                    int idx = iter->second;
                    for(int i=0;i<cov_num;i++)
                    {
                        if (vs_buf[i] != "-9" && vs_buf[i] != "NA")
                        {
                            bdata->_cov[i*bdata->_indi_num+idx]=vs_buf[i];
                        }
                    }
                }
            }
            for(int i=0;i<cov_num;i++)
                for(int j=0;j<bdata->_indi_num;j++)
                    if(bdata->_cov[i*bdata->_indi_num+j]=="-9") missindi.insert(pair<string,int>(bdata->_fid[j] + ":" + bdata->_pid[j],j));
            for (iter = missindi.begin(); iter != missindi.end(); iter++) indi_list.push_back(iter->first);
            update_map_rm(indi_list, bdata->_id_map, bdata->_keep);
            
        }
        in_cov.close();
        if(bdata->_keep.size()==0)
        {
            LOGPRINTF("Error: no individual is included.\n");
            TERMINATE();
        }
        LOGPRINTF("Non-missing %s covariates of %ld individuals are included from %s. \n",covtypestr.c_str(),bdata->_keep.size(), cov_file.c_str());
    }
    
    int construct_X(bInfo* bdata, vector<MatrixXd> &E_float, MatrixXd &qE_float, MatrixXd &_X) {
        
        int n=(int)bdata->_keep.size();
        
        int i = 0, j = 0;
        map<string, int>::iterator iter;
        stringstream errmsg;
        
        int  _X_c = 1;
        // quantitative covariates
        MatrixXd X_q;
        if (bdata->_qcov_num>0) {
            X_q.resize(n, bdata->_qcov_num);
            for (i = 0; i < n; i++) {
                for (j = 0; j < bdata->_qcov_num; j++) X_q(i, j) = bdata->_qcov[j*bdata->_indi_num+bdata->_keep[i]];
            }
            LOGPRINTF("%d quantitative variable(s) included as covariate(s).\n",bdata->_qcov_num);
            _X_c += bdata->_qcov_num;
        }
        
        // discrete covariates
        vector<MatrixXd> X_d;
        if (bdata->_cov_num>0) {
            vector< vector<string> > covar_tmp(bdata->_cov_num);
            for (i = 0; i < bdata->_cov_num; i++) covar_tmp[i].resize(n);
            for (i = 0; i < n; i++) {
                for (j = 0; j < bdata->_cov_num; j++) covar_tmp[j][i] = bdata->_cov[j*bdata->_indi_num+bdata->_keep[i]];
            }
            LOGPRINTF("%d discrete variable(s) included as covariate(s).\n",bdata->_cov_num);
            X_d.resize(bdata->_cov_num);
            for (i = 0; i < bdata->_cov_num; i++) {
                stringstream errmsg;
                errmsg << "Error: too many classes for the " << i + 1 << "th discrete variable. \nPlease use the --qcovar if it is a quantitative covariate.";
                string errmsg1 = errmsg.str();
                errmsg.str("");
                errmsg << "Error: the " << i + 1 << "th discrete variable has only one class.";
                string errmsg2 = errmsg.str();
                coeff_mat(covar_tmp[i], X_d[i], errmsg1, errmsg2);
                _X_c += (X_d[i]).cols() - 1;
            }
        }
        // E factor
        _X_c += qE_float.cols();
        for (i = 0; i < E_float.size(); i++) _X_c += (E_float[i]).cols() - 1;
        
        // Construct _X
        int col = 0;
        _X.resize(n, _X_c);
        _X.block(0, col, n, 1) = MatrixXd::Ones(n, 1);
        col++;
        if (bdata->_qcov_num>0) {
            _X.block(0, col, n, X_q.cols()) = X_q;
            col += X_q.cols();
        }
        for (i = 0; i < X_d.size(); i++) {
            _X.block(0, col, n, (X_d[i]).cols() - 1) = (X_d[i]).block(0, 1, n, (X_d[i]).cols() - 1);
            col += (X_d[i]).cols() - 1;
        }
        if (qE_float.cols() > 0) {
            _X.block(0, col, n, qE_float.cols()) = qE_float;
            col += qE_float.cols();
        }
        for (i = 0; i < E_float.size(); i++) {
            _X.block(0, col, n, (E_float[i]).cols() - 1) = (E_float[i]).block(0, 1, n, (E_float[i]).cols() - 1);
            col += (E_float[i]).cols() - 1;
        }
        return _X_c;
    }
    int comp_assoc(const void *a,const void *b){ return ((*(ASSOCRLT *)a).PVAL < (*(ASSOCRLT *)b).PVAL)?1:-1; } //decend
    void  testLinear(vector<ASSOCRLT> &assoc_rlts,char* outfileName, bInfo* bdata, MatrixXd &COV_plus)
    {
        FILE* assoc=NULL;
        long write_count=0;
        string outstr="";
        string assocfile="";
        if(outfileName!=NULL){
            assocfile = string(outfileName)+".linear";
            assoc = fopen(assocfile.c_str(), "w");
            if (!(assoc)) {
                LOGPRINTF("ERROR: open error %s\n", assocfile.c_str());
                TERMINATE();
            }
            outstr="Chr\tSNP\tbp\ta1\ta2\tb\tse\tp\tNMISS\n";
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
        }
        assoc_rlts.clear();
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(bdata, E_float, qE_float,_X);
        if(COV_plus.cols()>0)
        {
            LOGPRINTF("Adding additional covariates...\n");
            if(COV_plus.rows()!=_X.rows())
            {
                LOGPRINTF("ERROR: the row number of the additional covariates matrix %ld is not identical with the individual number %ld.\n",COV_plus.rows(),_X.rows());
                TERMINATE();
            }
            MatrixXd _X_PC(_X.rows(), _X.cols()+COV_plus.cols());
            _X_PC << _X, COV_plus;
            _X=_X_PC;
            _X_c=(int)_X.cols();
            //LOGPRINTF("%ld covariates are attached.\n",COV_plus.cols());
        }
        LOGPRINTF("\nPerforming linear regression analysis...\n");
        double cr=0;
        for(int i=0;i<bdata->_include.size();i++)
        {
            double desti=1.0*i/(bdata->_include.size()-1);
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            vector<double> yvec;
            vector<double> cvec;
            vector<double> xvec;
            MatrixXd X=_X;
            VectorXd _x;
            uint32_t snpid = i;
            bool has_missing = make_Xvec_subset(bdata, snpid, _x,false);
            double nonmiss=0.0;
            int miss=0;
            int chr=bdata->_chr[bdata->_include[i]];
            string rsid=bdata->_snp_name[bdata->_include[i]];
            int BP=bdata->_bp[bdata->_include[i]];
            string a1 = bdata->_allele1[bdata->_include[i]];
            string a2 = bdata->_allele2[bdata->_include[i]];
            double freq = bdata->_mu[bdata->_include[i]]/2;
            for(int j=0; j<bdata->_keep.size(); j++)
            {
                double phval=bdata->_pheno[bdata->_keep[j]];
                if(abs(phval+MISSING_PHENO)>1e-8) // no phval missing here, cos individuals with missing phenotype were removed when read phenotype file.
                {
                    double val=_x(j);
                    if(val!=0){
                        yvec.push_back(phval);
                        xvec.push_back(val);
                        nonmiss+=1.0;
                    } else {
                        removeRow(X, j-miss);
                        miss++;
                    }
                }
            }
            
            VectorXd y(yvec.size()),x(xvec.size());
            for(int j=0;j<xvec.size();j++)
            {
                y(j)=yvec[j];
                x(j)=xvec[j];
            }
            vector<double> rst;
            bool notInvertible=lin(y, X, x, rst);
            if(notInvertible) {
                LOGPRINTF("ERROR: The matrix for probe %s is not invertible. Maybe there are multicollinearity in the covariates.\n",rsid.c_str());
                TERMINATE();
            }
            ASSOCRLT currlt;
            if(assoc) {
                string outstr = atos(chr) + '\t' + rsid + '\t' + atos(BP)  + '\t' +atos(a1) +'\t'  +atos(a2) +'\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            }
            
            currlt.BETA=rst[0];
            currlt.SE=rst[1];
            currlt.PVAL=rst[2];
            currlt.CHR=chr;
            strcpy2(&currlt.SNP, rsid);
            strcpy2(&currlt.a1, a1);
            strcpy2(&currlt.a2, a2);
            currlt.freq = freq;
            currlt.BP=BP;
            currlt.NMISS=nonmiss;
            assoc_rlts.push_back(currlt);
            
        }
        if(assoc){
            LOGPRINTF("Results of %ld SNPs have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            LOGPRINTF("Results of %ld SNPs have been returned.\n",assoc_rlts.size());
        }
    }
    void  testLogit(vector<ASSOCRLT> &assoc_rlts,char* outfileName, bInfo* bdata, MatrixXd &COV_plus)
    {
        FILE* assoc=NULL;
        long write_count=0;
        string outstr="";
        string assocfile="";
        if(outfileName!=NULL){
            assocfile = string(outfileName)+".logistic";
            assoc = fopen(assocfile.c_str(), "w");
            if (!(assoc)) {
                LOGPRINTF("ERROR: open error %s\n", assocfile.c_str());
                TERMINATE();
            }
            outstr="Chr\tSNP\tbp\ta1\ta2\tb\tse\tp\tNMISS\n";
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
        }
        
        vector<double> bvec(bdata->_include.size()), sevec(bdata->_include.size()), nmvec(bdata->_include.size());
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(bdata, E_float, qE_float,_X);
        if(COV_plus.cols()>0)
        {
            LOGPRINTF("Adding additional covariates...\n");
            if(COV_plus.rows()!=_X.rows())
            {
                LOGPRINTF("ERROR: the row number of the additional covariates matrix %ld is not identical with the individual number %ld.\n",COV_plus.rows(),_X.rows());
                TERMINATE();
            }
            MatrixXd _X_PC(_X.rows(), _X.cols()+COV_plus.cols());
            _X_PC << _X, COV_plus;
            _X=_X_PC;
            _X_c=(int)_X.cols();
            //LOGPRINTF("%ld covariates are attached.\n",COV_plus.cols());
        }
        
        LOGPRINTF("\nPerforming logistic regression analysis...\n");
        #pragma omp parallel for
        for(int i=0;i<bdata->_include.size();i++)
        {
            printf("%3.0f%%\r", 100.0*i/bdata->_include.size());
            fflush(stdout);
            
            vector<int> yvec;
            vector<double> cvec;
            vector<double> xvec;
            MatrixXd X=_X;
            VectorXd _x;
            uint32_t snpid = i;
            bool has_missing = make_Xvec_subset(bdata, snpid, _x,false);
            double nonmiss=0.0;
            int miss=0;
            int chr=bdata->_chr[bdata->_include[i]];
            string rsid=bdata->_snp_name[bdata->_include[i]];
            int BP=bdata->_bp[bdata->_include[i]];
            for(int j=0; j<bdata->_keep.size(); j++)
            {
                int phval=bdata->_pheno[bdata->_keep[j]];
                if(abs(phval+MISSING_PHENO)>1e-8) // no phval missing here, cos individuals with missing phenotype were removed when read phenotype file.
                {
                    double val=_x(j);
                    if(val!=0){
                        yvec.push_back(phval);
                        xvec.push_back(val);
                        nonmiss+=1.0;
                    } else {
                        removeRow(X, j-miss);
                        miss++;
                    }
                }
            }
            
            VectorXi y(yvec.size());
            VectorXd x(xvec.size());
            for(int j=0;j<xvec.size();j++)
            {
                y(j)=yvec[j];
                x(j)=xvec[j];
            }
            MatrixXd LOGITX(X.rows(), X.cols()+1);
            LOGITX.block(0,0,X.rows(),X.cols())=X;
            long x_idx=LOGITX.cols()-1;
            LOGITX.col(x_idx)=x;
            double beta, se, pval;
            LogisticReg(y, LOGITX, beta, se, pval);
            bvec[i]=beta;
            sevec[i]=se;
            nmvec[i]=nonmiss;
        }
        
        if(assoc){
            for(int i=0;i<bdata->_include.size();i++)
            {
                int idx=bdata->_include[i];
                double z=bvec[i]/sevec[i];
                double pval=pchisq(z*z, 1);
                string outstr = atos(bdata->_chr[idx]) + '\t' + bdata->_snp_name[idx] + '\t' + atos(bdata->_bp[idx]) + '\t' + atos(bdata->_allele1[idx]) + '\t' + atos(bdata->_allele2[idx]) + '\t' + (sevec[i]<0?"NA":atos(exp(bvec[i]))) + '\t' + (sevec[i]<0?"NA":atos(sevec[i])) + '\t' + (sevec[i]<0?"NA":dtos(pval)) + '\t' + atos(nmvec[i]) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            }
            LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            assoc_rlts.clear();
            assoc_rlts.resize(bdata->_include.size());
            for(int i=0;i<bdata->_include.size();i++)
            {
                int chr=bdata->_chr[bdata->_include[i]];
                string prbid=bdata->_snp_name[bdata->_include[i]];
                string a1=bdata->_allele1[bdata->_include[i]];
                int BP=bdata->_bp[bdata->_include[i]];
                string a2=bdata->_allele2[bdata->_include[i]];
                double freq = bdata->_mu[bdata->_include[i]]/2;
                double beta=bvec[i], se=sevec[i];
                double z=beta/se;
                double pval=pchisq(z*z, 1);
                double OR=exp(beta);
                assoc_rlts[i].BETA=OR;
                assoc_rlts[i].SE=se;
                assoc_rlts[i].CHR=chr;
                assoc_rlts[i].BP=BP;
                assoc_rlts[i].NMISS=nmvec[i];
                assoc_rlts[i].PVAL=pval;
                strcpy2(&assoc_rlts[i].SNP, prbid);
                strcpy2(&assoc_rlts[i].a1, a1);
                strcpy2(&assoc_rlts[i].a2, a1);
                assoc_rlts[i].freq = freq;
            }
            
            LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
        }
        
        
    }
    void backward_elimn(bInfo* bdata, vector<int> &slct, MatrixXd &snp_profile, double p_thresh,double fdr_thresh, vector<int> &dump, double &rsq )
    {
        int nindi=(int)bdata->_keep.size();
        while(slct.size()>1)
        {
            int X_c=1+(int)slct.size();;
            MatrixXd X(nindi, X_c);
            X.block(0, 0, nindi, 1) = MatrixXd::Ones(nindi, 1);
            for(int i=0;i<slct.size();i++) X.col(i+1)=snp_profile.col(slct[i]);
            MatrixXd XtX_i;
            XtX_i=X.transpose()*X;
            bool determinant_zero=false;
            inverse_V(XtX_i, determinant_zero);
            
            VectorXd y(bdata->_keep.size());
            #pragma omp parallel for
            for(int j=0; j<bdata->_keep.size(); j++)
            {
                y(j)=bdata->_pheno[bdata->_keep[j]];
            }
            VectorXd b_hat=XtX_i*X.transpose()*y;
            double sst=y.transpose()*y;
            double ssr=y.transpose()*X*b_hat;
            rsq=ssr/sst;
            VectorXd residual=(y-X*b_hat);
            residual=residual.array()*residual.array();
            double sy=sqrt(residual.sum()/(y.size()-X.cols()));
            VectorXd se=sy*XtX_i.diagonal().array().sqrt();
            VectorXd t=b_hat.array()/se.array();
            vector<double> p(t.size()-1);
            for(int j=1;j<t.size();j++) p[j-1]=pchisq( t(j)*t(j),1 );
            int m = (int)(max_element(p.begin(), p.end()) - p.begin());
            if(p_thresh<0)
            {
                vector<double> jointq;
                getQval(p,jointq);
                if(jointq[m]>fdr_thresh)
                {
                    dump.push_back(slct[m]);
                    slct.erase(slct.begin()+m);
                } else break;
            }
            else
            {
                if(p[m]>p_thresh)
                {
                    dump.push_back(slct[m]);
                    slct.erase(slct.begin()+m);
                } else break;
            }
        }
        
    }
    bool forward_slct(bInfo* bdata, vector<int> &slct,vector<int> &remain, MatrixXd &snp_profile, double p_thresh, double fdr_thresh, vector<int> &dump , double rsqthresh)
    {
        vector<double> rlt;
        int nindi=(int)bdata->_keep.size(), X_c=1+(int)slct.size();;
        MatrixXd X(nindi, X_c);
        X.block(0, 0, nindi, 1) = MatrixXd::Ones(nindi, 1);
        for(int i=0;i<slct.size();i++) X.col(i+1)=snp_profile.col(slct[i]);
        
        MatrixXd XtX_i;
        XtX_i=X.transpose()*X;
        bool determinant_zero=false;
        inverse_V(XtX_i, determinant_zero);
        if(determinant_zero)
        {
            dump.push_back(slct[slct.size()-1]);
            slct.erase(slct.end()-1);
            slct.push_back(remain[0]);
            remain.erase(remain.begin());
            return false;
        }
        MatrixXd XtXiXt=XtX_i*X.transpose();
        
        VectorXd y(bdata->_keep.size());
#pragma omp parallel for
        for(int j=0; j<bdata->_keep.size(); j++)
        {
            y(j)=bdata->_pheno[bdata->_keep[j]];
        }
        VectorXd b_hat=XtXiXt*y;
        double sst=y.transpose()*y;
        double ssr=y.transpose()*X*b_hat;
        double rsq=ssr/sst;
        rsq=1-(1-rsq)*(nindi-1)/(nindi-X_c); // adjusted R-squared
        if(rsq > rsqthresh)
        {
            if(loud)
            {
                printf("Adjusted R-squared: %f\n",rsq);
            }
            return false;
        }
        VectorXd yresi=(y-X*b_hat);// if no covariate yresi=y-mean(y)
        vector<double> condp;
        for(int i=0;i<remain.size();i++)
        {
            VectorXd x=snp_profile.col(remain[i]);
            VectorXd b_hat=XtXiXt*x;
            VectorXd residual=(x-X*b_hat);
            vector<double> rst;
            adjusted_reg(yresi, residual, rst,X_c-1);
            condp.push_back(rst[2]);
        }
        int m = (int)(min_element(condp.begin(), condp.end()) - condp.begin());
        if(p_thresh<0)
        {
            vector<double> condq;
            getQval(condp,condq);
            if(condq[m]>fdr_thresh)
            {
                for(int j=0;j<remain.size();j++) dump.push_back(remain[j]);
                remain.clear();
                return false;
            }
            else {
                slct.push_back(remain[m]);
                remain.erase(remain.begin()+m);
                return true;
            }
        }
        else
        {
            if(condp[m]>p_thresh)
            {
                for(int j=0;j<remain.size();j++) dump.push_back(remain[j]);
                remain.clear();
                return false;
            }
            else {
                slct.push_back(remain[m]);
                remain.erase(remain.begin()+m);
                return true;
            }
        }

    }
   
    void stepwise_slct(bInfo* bdata,vector<int> &sig, vector<int> &slctid, vector<int> &rmid, vector<ASSOCRLT> &assoc_rlts, double p_cutoff,double fdr_cutoff, bool updatesig, bool swforward, double swrsq)
    {
        vector<double> p_buf;
        map<string, int>::iterator iter;
        MatrixXd sig_profile(bdata->_keep.size(),sig.size());
        int p=1,m=0;
        for(int i=0;i<sig.size();i++)
        {
            if(p>assoc_rlts[sig[i]].PVAL)
            {
                p=assoc_rlts[sig[i]].PVAL;
                m=i;
            }
            iter=bdata->_snp_name_map.find(assoc_rlts[sig[i]].SNP);
            if(iter!=bdata->_snp_name_map.end())
            {
                vector<int> MISS;
                double mu= 0.0, sum=0.0;
                int idx=iter->second;
                VectorXd vec;
                make_Xvec_subset(bdata, idx, vec, false);
                for(int j=0; j<bdata->_keep.size(); j++)
                {
                    double val=vec(j);
                    if(val != 0)
                    {
                        sig_profile(j,i)=val;
                        sum+=val;
                    } else {
                        MISS.push_back(j);
                    }
                    
                    double ratio=1.0*MISS.size()/bdata->_keep.size();
                    if(ratio>0.1)
                    {
                        //printf("WARNING: %f %% of the values in snp %s are missing.\n",ratio*100, bdata->_snp_name[idx].c_str());
                    }
                    if(MISS.size()>0)
                    {
                        if(bdata->_keep.size() > MISS.size()) mu=sum/(bdata->_keep.size()-MISS.size());
                        for(int k=0;k<MISS.size();k++)
                        {
                            sig_profile(MISS[k],i)=mu;
                        }
                    }
                }
            }
            else
            {
                LOGPRINTF("ERROR: probe %s can't be found.\n",assoc_rlts[sig[i]].SNP);
                TERMINATE();
            }
        }
        
        vector<int> slctidx, remain;
        if(p_cutoff >0)
        {
            if (assoc_rlts[sig[m]].PVAL >= p_cutoff) return;
        }
        else
        {
            vector<double> ptmp(sig.size()),qtmp;
            for(int i=0;i<sig.size();i++) ptmp[i]=assoc_rlts[sig[i]].PVAL;
            getQval(ptmp, qtmp);
            if(qtmp[m] >= fdr_cutoff) return;
        }
        
        slctidx.push_back(m);
        for (int i = 1; i < sig.size(); i++) {
            remain.push_back(i);
        }
        while (!remain.empty()) {
            if (forward_slct(bdata,slctidx, remain, sig_profile,p_cutoff,fdr_cutoff, rmid, swrsq)) {
                if(loud){
                    LOGPRINTF("Forward: Selected %ld, remained %ld, removed %ld \n",slctidx.size(),remain.size(),rmid.size());
                }
                if(!swforward) {
                    double rsq = 0.0;
                    backward_elimn(bdata,slctidx, sig_profile,p_cutoff,fdr_cutoff,rmid, rsq);
                    if(rsq > swrsq)
                    {
                        if(loud){
                            LOGPRINTF("Multiple R-squared: %f\n",rsq);
                        }
                        for(int j=0;j<remain.size();j++) rmid.push_back(remain[j]);
                        remain.clear();
                    }
                }
                if(!swforward && loud){
                    LOGPRINTF("Backward: Selected %ld, remained %ld, removed %ld \n",slctidx.size(),remain.size(),rmid.size());
                }
            } else {
                for(int j=0;j<remain.size();j++) rmid.push_back(remain[j]);
                remain.clear();
                if(loud){
                    LOGPRINTF("Forward: Selected %ld, remained %ld, removed %ld \n",slctidx.size(),remain.size(),rmid.size());
                }
            }
        }
        stable_sort(slctidx.begin(), slctidx.end());
        stable_sort(rmid.begin(), rmid.end());
        remain.resize(rmid.size());
        for(int i=0;i<rmid.size();i++) remain[i]=sig[rmid[i]];
        rmid.swap(remain);
        remain.resize(slctidx.size());
        slctid.resize(slctidx.size());
        for(int i=0;i<slctidx.size();i++)
        {
            remain[i]=sig[slctidx[i]];
            slctid[i]=sig[slctidx[i]];
        }
        if(updatesig)
        {
            sig.swap(remain);
            LOGPRINTF("%ld SNPs are selected from stepwise selection.\n", slctidx.size());
        }
    }
    void detect_family(bInfo* bdata,  vector<MatrixXd> &_A)
    {
        cout<<"Detecting sub-matrix for each family from the GRM ..."<<endl;
        int _n=(int)bdata->_keep.size();
        int i=0, j=0, k=0, l=0, prev_pnt=0;
        double d_buf1=0.0, d_buf2=0.0;
        bdata->_fam_brk_pnt.clear();
        for(i=0; i<_n-1; i++){
            d_buf1=_A[0].row(i).tail(_n-i-1).sum();
            d_buf2=_A[0].col(i).tail(_n-i-1).sum();
            if(FloatEqual(d_buf1, 0.0) && FloatEqual(d_buf2, 0.0)) bdata->_fam_brk_pnt.push_back(i);
        }
        
        bdata->_Asp.resize(bdata->_r_indx.size());
        for(i=0; i<bdata->_r_indx.size(); i++) (bdata->_Asp[i]).resize(_n, _n);
        
        int pos=0;
        for(l=0; l<bdata->_r_indx.size()-1; l++){
            pos=bdata->_r_indx[l];
            prev_pnt=0;
            for(k=0; k<bdata->_fam_brk_pnt.size(); k++){
                for(j=prev_pnt; j<=bdata->_fam_brk_pnt[k]; j++){
                    (bdata->_Asp[pos]).startVec(j);
                    for(i=prev_pnt; i<=bdata->_fam_brk_pnt[k]; i++) (bdata->_Asp[pos]).insertBack(i,j)=(_A[pos])(i,j);
                }
                prev_pnt=bdata->_fam_brk_pnt[k]+1;
            }
            for(j=prev_pnt; j<_n; j++){
                (bdata->_Asp[pos]).startVec(j);
                for(i=prev_pnt; i<_n; i++) (bdata->_Asp[pos]).insertBack(i,j)=(_A[pos])(i,j);
            }
            (bdata->_Asp[pos]).finalize();
        }
        pos=bdata->_r_indx[bdata->_r_indx.size()-1];
        for(i=0; i<_n; i++){
            (bdata->_Asp[pos]).startVec(i);
            (bdata->_Asp[pos]).insertBack(i,i)=1.0;
        }
        (bdata->_Asp[pos]).finalize();
        cout<<"There are "<<bdata->_fam_brk_pnt.size()+1<<" sub-matrices detected."<<endl;
        
        // release momery
        _A.clear();
    }
    void init_binfo(bInfo* bdata)
    {
        bdata->_snp_num=0;
        bdata->_indi_num=0;
        bdata->_cov_num=0;
        bdata->_qcov_num=0;
        bdata->_pheno_num=0;
        bdata->_mu.clear();

        
        bdata->_autosome_num=22;
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_snp_name_map.clear();
        bdata->_bp.clear();
        bdata->_include.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        
        
        
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        bdata->_cov.clear();
        bdata->_qcov.clear();
        bdata->_keep.clear();
        bdata->_id_map.clear();
        
        
        //
        bdata->_reml_mtd = 0;
        bdata->_reml_max_iter = 100;
        bdata->_V_inv_mtd = 0;
        bdata->_reml_force_inv = false;
        bdata->_reml_force_converge = false;
        bdata->_reml_no_converge = false;
        bdata->_reml_AI_not_invertible = false;
        
        bdata->_r_indx.clear();
        bdata->_r_indx_drop.clear();
        bdata->_var_name.clear();
        bdata->_varcmp.clear();
        bdata->_hsq_name.clear();
        
        bdata->_within_family = false;
        bdata->_fam_brk_pnt.clear();
        bdata->_Asp.clear();
        bdata->_Asp_prev.clear();
        bdata->_y_Ssq = 0;
        bdata->_fixed_rg_val.clear();
        bdata->_reml_fixed_var = 0;
        
        bdata->_bivar_reml = false;
        bdata->_bivar_no_constrain =  false;
        bdata->_ignore_Ce = false;
        bdata->_y2_Ssq = 0.0;
        bdata->_bivar_pos.clear();
        bdata->_bivar_pos_prev.clear();
        
        bdata->_ncase = 0.0;
        bdata->_ncase2 = 0.0;
        bdata->_flag_CC = false;
        bdata->_flag_CC2 = false;
        
    }
    void free_assoclist(vector<ASSOCRLT> &a)
    {
        for(int i=0;i<a.size();i++)
        {
            if(a[i].SNP!=NULL) free2(&a[i].SNP);
        }
    }
}

//
//  l3_vqtl.cpp
//  osc
//
//  Created by Futao Zhang on 29/08/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#include "l3_vqtl.hpp"

namespace VQTL {
    void indi_check(bInfo* bdata,eInfo* einfo)
    {
        map<string, int>::iterator iter;
        vector<string> indi_list;
        vector<int> bkeep,ekeep;
        map<string,int> bmap,emap;
        
        for(int i=0;i<bdata->_keep.size();i++)
        {
            int idx=bdata->_keep[i];
            string iid=bdata->_fid[idx]+":"+bdata->_pid[idx];
            iter = einfo->_eii_map.find(iid);
            if (iter != einfo->_eii_map.end()) {
                indi_list.push_back(iid);
                bkeep.push_back(idx);
                ekeep.push_back(iter->second);
                bmap.insert(pair<string,int>(iid,idx));
                emap.insert(pair<string,int>(iid,iter->second));
            }
        }
        
        if(indi_list.size()==0){
            LOGPRINTF("No individual in common.\n");
            TERMINATE();
        }
        bdata->_keep.swap(bkeep);
        einfo->_eii_include.swap(ekeep);
        bdata->_id_map.swap(bmap);
        einfo->_eii_map.swap(emap);
        
        LOGPRINTF("%ld individuals in common are kept.\n",  einfo->_eii_include.size());
    }
    void indi_check2(bInfo* bdata,eInfo* einfo)
    {
        map<string, int>::iterator iter;
        vector<string> indi_list;
        vector<int> bkeep,ekeep;
        map<string,int> bmap,emap;
        
        for(int i=0;i<einfo->_eii_include.size();i++)
        {
            int idx=einfo->_eii_include[i];
            string iid=einfo->_eii_fid[idx]+":"+einfo->_eii_iid[idx];
            iter = bdata->_id_map.find(iid);
            if (iter != bdata->_id_map.end()) {
                indi_list.push_back(iid);
                ekeep.push_back(idx);
                bkeep.push_back(iter->second);
                emap.insert(pair<string,int>(iid,idx));
                bmap.insert(pair<string,int>(iid,iter->second));
            }
        }
        if(indi_list.size()==0){
            LOGPRINTF("No individual in common.\n");
            TERMINATE();
        }

        bdata->_keep.swap(bkeep);
        einfo->_eii_include.swap(ekeep);
        bdata->_id_map.swap(bmap);
        einfo->_eii_map.swap(emap);
        
        LOGPRINTF("%ld individuals in common are kept.\n",  einfo->_eii_include.size());
    }

    void load_vqtl_workspace(eInfo* einfo,bInfo* bdata, char* efileName, char* befileName,  char* bFileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id, char* covfileName, char* qcovfileName)
    {
        if(befileName==NULL && efileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile.\n");
            TERMINATE();
        }
        if(bFileName==NULL )
        {
            LOGPRINTF("Error: please input the Plink file by the option --bfile.\n");
            TERMINATE();
        }

        if(efileName!=NULL)
        {
            if(transposed) read_efile_t(efileName,einfo,efileType,no_fid_flag,valueType);
            else read_efile(efileName,einfo,efileType,no_fid_flag,valueType);
            eii_man(einfo,indilstName,indilst2remove);
            read_famfile(bdata, string(bFileName)+".fam");
            if(covfileName != NULL) read_cov(einfo, covfileName, false);
            if(qcovfileName != NULL) read_cov(einfo, qcovfileName, true);
            //commn individuals
            indi_check(bdata,einfo);
            epi_man(einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            if(tsk_ttl>1) extract_probe(einfo,  tsk_ttl,  tsk_id);
            if(dpvalfName!=NULL) filtering_with_detpval(einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(einfo, upperBeta, lowerBeta);
            if(mssratio_prob<1) filtering_with_missingratio(einfo, mssratio_prob);
            read_bimfile(bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(bdata, snplst2exclde);
            read_bedfile(bdata, string(bFileName)+".bed");
            
        }else{
            char inputname[FNAMESIZE];
            memcpy(inputname,befileName,strlen(befileName)+1);
            char* suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,einfo);
            eii_man(einfo,indilstName,indilst2remove);
            if(covfileName != NULL) read_cov(einfo, covfileName, false);
            if(qcovfileName != NULL) read_cov(einfo, qcovfileName, true);
            read_famfile(bdata, string(bFileName)+".fam");
            indi_check(bdata,einfo);
            memcpy(suffix,".opi",5);
            read_epi(inputname,einfo);
            epi_man(einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            if(tsk_ttl>1) extract_probe(einfo,  tsk_ttl,  tsk_id);
            memcpy(suffix,".bod",5);
            clock_t begin_time = clock();
            read_beed(inputname,einfo);
            //LOGPRINTF("read_bod: %f ms.\n",float( clock () - begin_time ) /  1000);
            if(dpvalfName!=NULL) filtering_with_detpval(einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(einfo, upperBeta, lowerBeta);
            if(mssratio_prob<1) filtering_with_missingratio(einfo, mssratio_prob);
            read_bimfile(bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(bdata, snplst2exclde);
             begin_time = clock();
            read_bedfile(bdata, string(bFileName)+".bed");
            //LOGPRINTF("read_bedfile: %f ms.\n",float( clock () - begin_time ) /  1000);
        }
        //test alignment.
        if(einfo->_eii_include.size() != bdata->_keep.size())
        {
            LOGPRINTF("Error: Failed in extracting samples in common. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
            TERMINATE();
        }
        for(int i=0;i<einfo->_eii_include.size();i++)
            if(einfo->_eii_fid[einfo->_eii_include[i]]!=bdata->_fid[bdata->_keep[i]]) {
                LOGPRINTF("Error: Failed in alignment. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                TERMINATE();
            }
        if(einfo->_eType == METHYLATION)
        {
            if(beta2m && m2beta){
                //no chance to enter here
                LOGPRINTF("Error: --m2beta should not be with --beta2m.\n");
                TERMINATE();
            }
            if(beta2m && einfo->_valType==BETAVALUE) beta_2_m(einfo);
            if(m2beta && einfo->_valType==MVALUE) m_2_beta(einfo);
        }
        einfo->autosome_num = autosome_num;
        if(bdata->_mu.empty()) calcu_mu(bdata);
    }

    
     void V_QTL(char* outFileName,  char* efileName, char* befileName,  char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,int vqtl_mtd,char* covfileName, char* qcovfileName, bool besd_snp_major,bool cis_flag, int cis_itvl)
    {
        
        setNbThreads(thread_num);
        eInfo einfo;
        bInfo bdata;
        init_einfo(&einfo);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, bFileName,transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num,snplstName,snplst2exclde,tsk_ttl,tsk_id,covfileName,qcovfileName); // using _keep and _eii_include, the individuals are aligned.
        if (maf > 0) filter_snp_maf(&bdata, maf);
        
        char outputname[FNAMESIZE];
        outputname[0]='\0';
        if(tsk_ttl>1) {
            if(outFileName!=NULL) { // outFileName could be the specified name or osca
                string tmp=  string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outFileName=outputname;
            }
        }
        
        if(vqtl_mtd==0) { LOGPRINTF("\nPerforming vQTL analysis with Bartlett’s test...\n");}
        else if (vqtl_mtd==1) {LOGPRINTF("\nPerforming vQTL analysis with Levene’s test (mean)...\n");}
        else if (vqtl_mtd==2) {LOGPRINTF("\nPerforming vQTL analysis with Levene’s test (median)...\n");}
        else  {LOGPRINTF("\nPerforming vQTL analysis with Fligner-Killeen test...\n");}
        
        LOGPRINTF("\nThe results would be saved in BESD format ...\n");
        if(outFileName!=NULL){
            write_smr_esi(outFileName, &bdata);
            write_smr_epi(outFileName, &einfo);
        }
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=SMR_DENSE_3;
        if(besd_snp_major && !cis_flag) filetype=OSCA_DENSE_1;
        else if(!besd_snp_major && cis_flag) filetype=OSCA_SPARSE_1;
        else if(besd_snp_major && cis_flag) {
            LOGPRINTF("ERROR: --besd-snp-besd and --cis-wind are mutual exclusive.\n");
            TERMINATE();
        } else {
            filetype=SMR_DENSE_3;
        }
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=(int)bdata._keep.size();
        ten_ints[2]=(int)bdata._include.size();
        ten_ints[3]=(int)einfo._epi_include.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
            TERMINATE();
        }
        int slide_wind=10000;
        int loops=ceil(1.0*bdata._include.size()/slide_wind);
        int warn1=0, warn2=0, warn3=0;
        clock_t begin_time = clock();
        if(filetype==SMR_DENSE_3)
        {
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                printf("%3.0f%%\r", 100.0*jj/einfo._epi_include.size());
                fflush(stdout);
                string prbid=einfo._epi_prb[einfo._epi_include[jj]];
                string gene=einfo._epi_gene[einfo._epi_include[jj]];

                vector<float> beta,se;
                for(int ii=0;ii<loops;ii++)
                {
                    MatrixXf _X;
                    int snpstart=ii*slide_wind;
                    make_XMat(&bdata,snpstart,slide_wind, _X);
                    
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        string snprs=bdata._snp_name[bdata._include[snpstart+kk]];
                        double snpfreq=bdata._mu[bdata._include[snpstart+kk]]/2;
                        if(snpfreq==0 || snpfreq==1) {
                            LOGPRINTF("ERROR: MAF found 0 or 1 with SNP %s.\n",snprs.c_str());
                            TERMINATE();
                        }
                        double nonmiss=0.0;
                        vector<double> yvec,bvec;
                        for(int ll=0; ll<einfo._eii_include.size(); ll++)
                        {
                            if(einfo._eii_fid[einfo._eii_include[ll]]!=bdata._fid[bdata._keep[ll]]) printf("ERROR: Alignment failed. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                            double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]]; // einfo._eii_include[ll] should equals ll. because of update_eii all the time after reading beed file.
                            double bval=_X(ll,kk);
                            if(val<1e9 && bval<1e5){
                                yvec.push_back(val);
                                bvec.push_back(bval);
                                nonmiss+=1.0;
                            }
                        }
                        vector<double> rst;
                        if(vqtl_mtd==0) {
                            int flag=bartlett(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                beta.push_back(-9);
                                se.push_back(-9);
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else if(flag==-2) {
                                beta.push_back(-9);
                                se.push_back(-9);
                                if(warn2<10){
                                    LOGPRINTF("WARNING: One of the allele counts is 1 of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn2++;
                                } else if (warn2==10) {
                                    LOGPRINTF("WARNING: More than 10 SNPs one of whose allele counts is 1. \n");
                                    warn2++;
                                }
                                
                            } else if(flag==-3) {
                                beta.push_back(-9);
                                se.push_back(-9);
                                if(warn3<10){
                                    LOGPRINTF("WARNING: The number of allele categories equals the sample size of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn3++;
                                } else if(warn3==10) {
                                    LOGPRINTF("WARNING: More than 10 SNPs whose allele number equals the sample size. \n");
                                    warn3++;
                                }
                            } else {
                                beta.push_back(rst[0]);
                                se.push_back(rst[1]);
                            }
                        } else if(vqtl_mtd==1)
                        {
                            int flag=leveneTest_mean(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                beta.push_back(-9);
                                se.push_back(-9);
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else {
                                beta.push_back(rst[0]);
                                se.push_back(rst[1]);
                            }
                        } else if(vqtl_mtd==2)
                        {
                            int flag=leveneTest_median(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                beta.push_back(-9);
                                se.push_back(-9);
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else {
                                beta.push_back(rst[0]);
                                se.push_back(rst[1]);
                            }
                            
                        } else if(vqtl_mtd==3)
                        {
                            int flag=flignerTest(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                beta.push_back(-9);
                                se.push_back(-9);
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else {
                                beta.push_back(rst[0]);
                                se.push_back(rst[1]);
                            }
                            
                        }
                        
                    }
                }
                if (fwrite_checked(&beta[0],beta.size()*sizeof(float), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                if (fwrite_checked(&se[0], se.size()*sizeof(float), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
            }
        }
        else if(filetype==OSCA_DENSE_1)
        {
            vector<float> betases;
            betases.resize(2*einfo._epi_include.size());
            for(int ii=0;ii<loops;ii++)
            {
                MatrixXf _X;
                int snpstart=ii*slide_wind;
                make_XMat(&bdata,snpstart,slide_wind, _X);
                
                for(int kk=0;kk<_X.cols();kk++)
                {
                    printf("%3.0f%%\r", 100.0*(ii*slide_wind+kk)/bdata._include.size());
                    fflush(stdout);
                    string snprs=bdata._snp_name[bdata._include[snpstart+kk]];
                    double snpfreq=bdata._mu[bdata._include[snpstart+kk]]/2;
                    if(snpfreq==0 || snpfreq==1) {
                        LOGPRINTF("ERROR: MAF found 0 or 1 with SNP %s.\n",snprs.c_str());
                        TERMINATE();
                    }
                    #pragma omp parallel for shared(warn1, warn2, warn3)
                    for(int jj=0;jj<einfo._epi_include.size();jj++)
                    {
                        string prbid=einfo._epi_prb[einfo._epi_include[jj]];
                        string gene=einfo._epi_gene[einfo._epi_include[jj]];

                        double nonmiss=0.0;
                        vector<double> yvec,bvec;
                        for(int ll=0; ll<einfo._eii_include.size(); ll++)
                        {
                            if(einfo._eii_fid[einfo._eii_include[ll]]!=bdata._fid[bdata._keep[ll]]) {
                                LOGPRINTF("ERROR: Alignment failed. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                                TERMINATE();
                            }
                            double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]]; // einfo._eii_include[ll] should equals ll. because of update_eii all the time after reading beed file.
                            double bval=_X(ll,kk);
                            if(val<1e9 && bval<1e5){
                                yvec.push_back(val);
                                bvec.push_back(bval);
                                nonmiss+=1.0;
                            }
                        }
                        vector<double> rst;
                        if(vqtl_mtd==0) {
                            int flag=bartlett(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else if(flag==-2) {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                if(warn2<10){
                                    LOGPRINTF("WARNING: One of the allele counts is 1 of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn2++;
                                } else if (warn2==10) {
                                    LOGPRINTF("WARNING: More than 10 SNPs one of whose allele counts is 1. \n");
                                    warn2++;
                                }
                                
                            } else if(flag==-3) {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                if(warn3<10){
                                    LOGPRINTF("WARNING: The number of allele categories equals the sample size of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn3++;
                                } else if(warn3==10) {
                                    LOGPRINTF("WARNING: More than 10 SNPs whose allele number equals the sample size. \n");
                                    warn3++;
                                }
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                            }
                        } else if(vqtl_mtd==1)
                        {
                            int flag=leveneTest_mean(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                            }
                        } else if(vqtl_mtd==2)
                        {
                            int flag=leveneTest_median(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                            }
                            
                        } else if(vqtl_mtd==3)
                        {
                            int flag=flignerTest(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                if(warn1<10){
                                    LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    warn1++;
                                } else if(warn1==10){
                                    LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                    warn1++;
                                }
                                
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                            }
                        }
                    }
                    if (fwrite_checked(&betases[0],betases.size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                }
                
            }
        } else if (filetype == OSCA_SPARSE_1)
        {
            vector<uint64_t> cols;
            vector<uint32_t> rowid;
            vector<float> vals;
            cols.resize(einfo._epi_include.size()+1);
            cols[0]=0;
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                printf("%3.0f%%\r", 100.0*jj/einfo._epi_include.size());
                fflush(stdout);
                
                MatrixXf _X;
                vector<float> se;
                string prbid=einfo._epi_prb[einfo._epi_include[jj]];
                string gene=einfo._epi_gene[einfo._epi_include[jj]];
                int prbchr=einfo._epi_chr[einfo._epi_include[jj]];
                int prbbp=einfo._epi_bp[einfo._epi_include[jj]];
                int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
                int cisend=prbbp+cis_itvl*1000;
                vector<uint32_t> snpids; //to save _include id not _include value
                for(int kk=0;kk<bdata._include.size();kk++)
                {
                    int snpchr=bdata._chr[bdata._include[kk]];
                    int snpbp=bdata._bp[bdata._include[kk]];
                    if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend) snpids.push_back(kk);
                }
                if(snpids.size()>0) {
                    LOGPRINTF("%ld SNPs are included for the cis-region of probe %s [%d:%d,%d:%d].\n",snpids.size(),prbid.c_str(),prbchr,cisstart,prbchr,cisend);
                } else {
                    LOGPRINTF("No SNP is included for the cis-region of probe %s [%d:%d,%d:%d].\n",prbid.c_str(),prbchr,cisstart,prbchr,cisend);
                    cols[jj+1]=cols[jj];
                    continue;
                }
                
                make_XMat(&bdata,snpids, _X);
                for(int kk=0;kk<_X.cols();kk++) //_X.cols() ==snpids.size()
                {
                    uint32_t snpid=snpids[kk];
                    string snprs=bdata._snp_name[bdata._include[snpid]];
                    double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                    if(snpfreq==0 || snpfreq==1) {
                        LOGPRINTF("ERROR: MAF found 0 or 1 with SNP %s.\n",snprs.c_str());
                        TERMINATE();
                    }
                    double nonmiss=0.0;
                    vector<double> yvec,bvec;
                    for(int ll=0; ll<einfo._eii_include.size(); ll++)
                    {
                        if(einfo._eii_fid[einfo._eii_include[ll]]!=bdata._fid[bdata._keep[ll]]) printf("ERROR: Alignment failed. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                        double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]]; // einfo._eii_include[ll] should equals ll. because of update_eii all the time after reading beed file.
                        double bval=_X(ll,kk);
                        if(val<1e9 && bval<1e5){
                            yvec.push_back(val);
                            bvec.push_back(bval);
                            nonmiss+=1.0;
                        }
                    }
                    vector<double> rst;
                    if(vqtl_mtd==0) {
                        int flag=bartlett(yvec,bvec, rst,snpfreq);
                        if(flag==-1)
                        {
                            if(warn1<10){
                                LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                warn1++;
                            } else if(warn1==10){
                                LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                warn1++;
                            }
                            
                        } else if(flag==-2) {
                            if(warn2<10){
                                LOGPRINTF("WARNING: One of the allele counts is 1 of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                warn2++;
                            } else if (warn2==10) {
                                LOGPRINTF("WARNING: More than 10 SNPs one of whose allele counts is 1. \n");
                                warn2++;
                            }
                            
                        } else if(flag==-3) {
                            if(warn3<10){
                                LOGPRINTF("WARNING: The number of allele categories equals the sample size of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                warn3++;
                            } else if(warn3==10) {
                                LOGPRINTF("WARNING: More than 10 SNPs whose allele number equals the sample size. \n");
                                warn3++;
                            }
                        } else {
                            rowid.push_back(snpid);
                            vals.push_back(rst[0]);
                            se.push_back(rst[1]);
                        }
                    } else if(vqtl_mtd==1)
                    {
                        int flag=leveneTest_mean(yvec,bvec, rst,snpfreq);
                        if(flag==-1)
                        {
                            if(warn1<10){
                                LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                warn1++;
                            } else if(warn1==10){
                                LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                warn1++;
                            }
                            
                        } else {
                            rowid.push_back(snpid);
                            vals.push_back(rst[0]);
                            se.push_back(rst[1]);
                        }
                    } else if(vqtl_mtd==2)
                    {
                        int flag=leveneTest_median(yvec,bvec, rst,snpfreq);
                        if(flag==-1)
                        {
                            if(warn1<10){
                                LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                warn1++;
                            } else if(warn1==10){
                                LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                warn1++;
                            }
                            
                        } else {
                            rowid.push_back(snpid);
                            vals.push_back(rst[0]);
                            se.push_back(rst[1]);
                        }
                        
                    } else if(vqtl_mtd==3)
                    {
                        int flag=flignerTest(yvec,bvec, rst,snpfreq);
                        if(flag==-1)
                        {
                            if(warn1<10){
                                LOGPRINTF("WARNING: No variation with SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                warn1++;
                            } else if(warn1==10){
                                LOGPRINTF("WARNING: More than 10 SNPs with no variation. \n");
                                warn1++;
                            }
                            
                        } else {
                            rowid.push_back(snpid);
                            vals.push_back(rst[0]);
                            se.push_back(rst[1]);
                        }
                    }
                    
                    
                    
                }
                for(int kk=0;kk<se.size();kk++)
                {
                    vals.push_back(se[kk]);
                }
                uint64_t real_num=se.size()*2;
                cols[jj+1]=real_num+cols[jj];
            }
            uint64_t valNum=vals.size();
            if (fwrite_checked(&valNum,sizeof(uint64_t), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if (fwrite_checked(&cols[0],cols.size()*sizeof(uint64_t), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if (fwrite_checked(&rowid[0],rowid.size()*sizeof(uint32_t), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if (fwrite_checked(&vals[0],vals.size()*sizeof(float), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
           
        }
        LOGPRINTF("vQTL tests cost: %f ms.\n",float( clock () - begin_time ) /  1000);
        fclose(besd);
        if (filetype == SMR_SPARSE_3) {
             LOGPRINTF("vQTL results in the cis-regions of %ld probes have been saved in sparse binary file %s.\n",einfo._epi_include.size(), besdName.c_str());
        } else {
             LOGPRINTF("vQTL results of %ld probes and %ld SNPs have been saved in file %s.\n",einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
        }
       
    }
    
    void eQTL(char* outFileName,  char* efileName, char* befileName,  char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,char* covfileName, char* qcovfileName, bool besd_snp_major)
    {
        
        setNbThreads(thread_num);
        eInfo einfo;
        bInfo bdata;
        init_einfo(&einfo);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, bFileName,transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num,snplstName,snplst2exclde,tsk_ttl,tsk_id,covfileName,qcovfileName); // using _keep and _eii_include, the individuals are aligned.
        if (maf > 0) filter_snp_maf(&bdata, maf);
        
        char outputname[FNAMESIZE];
        outputname[0]='\0';
        if(tsk_ttl>1) {
            if(outFileName!=NULL) {
                string tmp=  string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outFileName=outputname;
            }
        }
        
        LOGPRINTF("\nPerforming eQTL analysis ...\n");
        LOGPRINTF("\nThe results would be saved in dense BESD format ...\n");
        if(outFileName!=NULL){
            write_smr_esi(outFileName, &bdata);
            write_smr_epi(outFileName, &einfo);
        }
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=SMR_DENSE_3;
        if(besd_snp_major) filetype=OSCA_DENSE_1;
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=(int)bdata._keep.size();
        ten_ints[2]=(int)bdata._include.size();
        ten_ints[3]=(int)einfo._epi_include.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
            TERMINATE();
        }
       
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _Cov;
        int _X_c;
        _X_c=construct_X(&einfo, E_float, qE_float,_Cov);
        
        int slide_wind=10000;
        int loops=ceil(1.0*bdata._include.size()/slide_wind);
        //clock_t begin_time = clock();
        if(!besd_snp_major) {
            
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                printf("%3.0f%%\r", 100.0*jj/einfo._epi_include.size());
                fflush(stdout);
                
                vector<float> beta,se;
                for(int ii=0;ii<loops;ii++)
                {
                    MatrixXf _X;
                    int snpstart=ii*slide_wind;
                    make_XMat(&bdata,snpstart,slide_wind, _X);
                    
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        double nonmiss=0.0;
                        int miss=0;
                        
                        MatrixXd X(_Cov.rows(), _Cov.cols()+1);
                        X.block(0,0,_Cov.rows(),_Cov.cols())=_Cov;
                        long x_idx=X.cols()-1;
                        VectorXd y(einfo._eii_include.size());
                        
                        for(int ll=0; ll<einfo._eii_include.size(); ll++)
                        {
                            if(einfo._eii_fid[einfo._eii_include[ll]]!=bdata._fid[bdata._keep[ll]]) printf("ERROR: Alignment failed. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                            double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]]; // einfo._eii_include[ll] should equals ll. because of update_eii all the time after reading beed file.
                            double bval=_X(ll,kk);
                            if(val<1e9 && bval<1e5){
                                y(ll-miss)=val;
                                X(ll-miss,x_idx)=bval;
                                nonmiss+=1.0;
                            } else {
                                removeRow(X, ll-miss);
                                miss++;
                            }
                        }
                        if(miss>0) y.conservativeResize(einfo._eii_include.size()-miss);
                        vector<double> rst;
                        lin2(y, X, rst);
                        if(rst.size()!=3) {
                            LOGPRINTF("ERROR: bugs in function of linear regression found. Please report this.\n");
                            TERMINATE();
                        }
                        beta.push_back(rst[0]);
                        se.push_back(rst[1]);
                    }
                }
                if (fwrite_checked(&beta[0],beta.size()*sizeof(float), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                if (fwrite_checked(&se[0], se.size()*sizeof(float), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
            }
        } else {
            for(int ii=0;ii<loops;ii++) // one SNP one recoding
            {
                MatrixXf _X;
                int snpstart=ii*slide_wind;
                make_XMat(&bdata,snpstart,slide_wind, _X);
                for(int kk=0;kk<_X.cols();kk++)
                {
                    printf("%3.0f%%\r", 100.0*(ii*slide_wind+kk)/bdata._include.size());
                    fflush(stdout);
                    vector<float> betases;
                    for(int jj=0;jj<einfo._epi_include.size();jj++)
                    {
                        double nonmiss=0.0;
                        int miss=0;
                        MatrixXd X(_Cov.rows(), _Cov.cols()+1);
                        X.block(0,0,_Cov.rows(),_Cov.cols())=_Cov;
                        long x_idx=X.cols()-1;
                        VectorXd y(einfo._eii_include.size());

                        for(int ll=0; ll<einfo._eii_include.size(); ll++)
                        {
                            if(einfo._eii_fid[einfo._eii_include[ll]]!=bdata._fid[bdata._keep[ll]]) printf("ERROR: Alignment failed. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                            double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]]; // einfo._eii_include[ll] should equals ll. because of update_eii all the time after reading beed file.
                            double bval=_X(ll,kk);
                            if(val<1e9 && bval<1e5){
                                y(ll-miss)=val;
                                X(ll-miss,x_idx)=bval;
                                nonmiss+=1.0;
                            } else {
                                removeRow(X, ll-miss);
                                miss++;
                            }
                        }
                        if(miss>0) y.conservativeResize(einfo._eii_include.size()-miss);
                        
                        vector<double> rst;
                        lin2(y, X, rst);
                        if(rst.size()!=3) {
                            LOGPRINTF("ERROR: bugs in function of linear regression found. Please report this.\n");
                            TERMINATE();
                        }
                        betases.push_back(rst[0]);
                        betases.push_back(rst[1]);
                    }
                    if (fwrite_checked(&betases[0],betases.size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                }
                
            }
        }
        
        fclose(besd);
        LOGPRINTF("eQTL summary statistics for %ld probes and %ld SNPs are saved in file %s.\n", einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
       // LOGPRINTF("%ld probes cost: %f ms.\n",einfo._epi_include.size(),float( clock () - begin_time ) /  1000);
    }
}

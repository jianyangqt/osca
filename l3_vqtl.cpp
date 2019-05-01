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

    void load_vqtl_workspace(eInfo* einfo,bInfo* bdata, char* efileName, char* befileName, char* phenofileName, char* bFileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id, char* covfileName, char* qcovfileName,char* grm_file)
    {
        
        if(befileName==NULL && efileName==NULL && phenofileName==NULL)
        {
            LOGPRINTF("Error: please input the phenotype / Gene expression / Methylation data by the option --pheno, --efile or --befile.\n");
            TERMINATE();
        }
        if((befileName!=NULL || efileName!=NULL) && phenofileName!=NULL)
        {
            LOGPRINTF("Error: please do not input the phenotype or Gene expression / Methylation data together.\n");
            TERMINATE();
        }
        if(bFileName==NULL )
        {
            LOGPRINTF("Error: please input the Plink file by the option --bfile.\n");
            TERMINATE();
        }

        if(efileName!=NULL || phenofileName!=NULL)
        {
            vector<string> grm_id;
            if(phenofileName!=NULL)
            {
                read_pheno2(phenofileName, einfo,0);
            }
            else
            {
                if(transposed) read_efile_t(efileName,einfo,efileType,no_fid_flag,valueType);
                else read_efile(efileName,einfo,efileType,no_fid_flag,valueType);
            }
            eii_man(einfo,indilstName,indilst2remove);
            read_famfile(bdata, string(bFileName)+".fam");
            if(grm_file!=NULL) read_grm_id(bdata,grm_file, grm_id);
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
            
        }else if(befileName!=NULL){
            vector<string> grm_id;
            char inputname[FNAMESIZE];
            memcpy(inputname,befileName,strlen(befileName)+1);
            char* suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,einfo);
            eii_man(einfo,indilstName,indilst2remove);
            if(covfileName != NULL) read_cov(einfo, covfileName, false);
            if(qcovfileName != NULL) read_cov(einfo, qcovfileName, true);
            read_famfile(bdata, string(bFileName)+".fam");
            if(grm_file!=NULL) read_grm_id(bdata,grm_file, grm_id);
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

    
     void V_QTL(char* outFileName,  char* efileName, char* befileName,  char* phenofileName, char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,int vqtl_mtd,char* covfileName, char* qcovfileName, bool tosmrflag, bool cis_flag, int cis_itvl)
    {
        
        setNbThreads(thread_num);
        eInfo einfo;
        bInfo bdata;
        init_einfo(&einfo);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, phenofileName,bFileName,transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num,snplstName,snplst2exclde,tsk_ttl,tsk_id,covfileName,qcovfileName); // using _keep and _eii_include, the individuals are aligned.
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
        
        uint32_t filetype=OSCA_DENSE_1;
        if(tosmrflag && !cis_flag) filetype=SMR_DENSE_3;
        else if(!tosmrflag && cis_flag) filetype=OSCA_SPARSE_1;
        else if(tosmrflag && cis_flag) filetype=SMR_SPARSE_3;
        
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=(int)bdata._keep.size();
        ten_ints[2]=(int)bdata._include.size();
        ten_ints[3]=(int)einfo._epi_include.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        FILE* novar=NULL;
        string novarFName=string(outFileName)+".novar.list"; // like all are AA, or all are AC, or all are CC
        if(fopen_checked(&novar, novarFName.c_str(),"w")) TERMINATE();
        FILE* singleton=NULL;
        string singletonFName=string(outFileName)+".singleton.list"; // like only one individula is AA
        if(fopen_checked(&singleton, singletonFName.c_str(),"w")) TERMINATE();
        
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        string outstr="";
        if(phenofileName!=NULL)
        {
            besdName=string(outFileName)+".vqtl";
            if(fopen_checked(&besd, besdName.c_str(),"w")) TERMINATE();
            outstr="Chr\tSNP\tbp";
            if(vqtl_mtd==0 || vqtl_mtd==3 ) outstr+="\tstatistic\tdf\tbeta\tse\tP\tNMISS\n";
            else outstr+="\tF-statistic\tdf1\tdf2\tbeta\tse\tP\tNMISS\n";
            if(fputs_checked(outstr.c_str(),besd))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                TERMINATE();
            }
            
        } else {
            LOGPRINTF("\nThe results would be saved in BESD format ...\n");
            if(outFileName!=NULL){
                write_smr_esi(outFileName, &bdata);
                write_smr_epi(outFileName, &einfo);
            }
            
            if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
            if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
        }
        
        int slide_wind=10000;
        int loops=ceil(1.0*bdata._include.size()/slide_wind);
       //clock_t begin_time = clock();
        long write_count=0;
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
                    MatrixXd _X;
                    int snpstart=ii*slide_wind;
                    make_XMat(&bdata,snpstart,slide_wind, _X);
                    
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        string snprs=bdata._snp_name[bdata._include[snpstart+kk]];
                        double snpfreq=bdata._mu[bdata._include[snpstart+kk]]/2;
                        if(snpfreq==0 || snpfreq==1)
                        {
                            string tmpstr=snprs+'\n';
                            if(fputs_checked(tmpstr.c_str(),novar))
                            {
                                LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                TERMINATE();
                            }
                            beta.push_back(-9);
                            se.push_back(-9);
                            continue;
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
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
                                }
                                
                            } else if(flag==-2) {
                                beta.push_back(-9);
                                se.push_back(-9);
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),singleton))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", singletonFName.c_str());
                                    TERMINATE();
                                }
                            } else if(flag==-3) {
                                beta.push_back(-9);
                                se.push_back(-9);
                                printf("WARNING: The number of allele categories equals the sample size of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str()); // would hardly happen
                                
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
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
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
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
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
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
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
                MatrixXd _X;
                int snpstart=ii*slide_wind;
                make_XMat(&bdata,snpstart,slide_wind, _X);
                
                for(int kk=0;kk<_X.cols();kk++)
                {
                    printf("%3.0f%%\r", 100.0*(ii*slide_wind+kk)/bdata._include.size());
                    fflush(stdout);
                    string snprs=bdata._snp_name[bdata._include[snpstart+kk]];
                    int snpchr=bdata._chr[bdata._include[snpstart+kk]];
                    int snpbp=bdata._bp[bdata._include[snpstart+kk]];
                    double snpfreq=bdata._mu[bdata._include[snpstart+kk]]/2;
                    if(snpfreq==0 || snpfreq==1)
                    {
                        string tmpstr=snprs+'\n';
                        if(fputs_checked(tmpstr.c_str(),novar))
                        {
                            LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                            TERMINATE();
                        }
                        if(phenofileName==NULL)
                        {
                            for(int jj=0;jj<betases.size();jj++) betases[jj]=-9;
                            if (fwrite_checked(&betases[0],betases.size()*sizeof(float), besd))
                            {
                                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                TERMINATE();
                            }
                        }
                        continue;
                    }

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
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
                                }
                                
                            } else if(flag==-2) {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),singleton))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", singletonFName.c_str());
                                    TERMINATE();
                                }
                                
                            } else if(flag==-3) {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                               printf("WARNING: The number of allele categories equals the sample size of SNP %s. vQTL of %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                                    
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                                if(phenofileName!=NULL && besd) {
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }

                            }
                        } else if(vqtl_mtd==1)
                        {
                            int flag=leveneTest_mean(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
                                }
                                
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                                if(phenofileName!=NULL && besd) {
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[5]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + atos(rst[2]) + '\t'  + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }

                            }
                        } else if(vqtl_mtd==2)
                        {
                            int flag=leveneTest_median(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
                                }
                                
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                                if(phenofileName!=NULL && besd) {
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[5]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + atos(rst[2]) + '\t'  + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }
                            }
                            
                        } else if(vqtl_mtd==3)
                        {
                            int flag=flignerTest(yvec,bvec, rst,snpfreq);
                            if(flag==-1)
                            {
                                betases[2*jj]=-9;
                                betases[2*jj+1]=-9;
                                string tmpstr=snprs+'\n';
                                if(fputs_checked(tmpstr.c_str(),novar))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                    TERMINATE();
                                }
                                
                            } else {
                                betases[2*jj]=rst[0];
                                betases[2*jj+1]=rst[1];
                                if(phenofileName!=NULL && besd) {
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }

                            }
                        }
                    }
                    if(phenofileName==NULL)
                    {
                        if (fwrite_checked(&betases[0],betases.size()*sizeof(float), besd))
                        {
                            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                            TERMINATE();
                        }
                    }
                }
                
            }
        }
        else if (filetype == OSCA_SPARSE_1)
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
                
                MatrixXd _X;
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
                        string tmpstr=snprs+'\n';
                        if(fputs_checked(tmpstr.c_str(),novar))
                        {
                            LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                            TERMINATE();
                        }
                        continue;
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
                            string tmpstr=snprs+'\n';
                            if(fputs_checked(tmpstr.c_str(),novar))
                            {
                                LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                TERMINATE();
                            }
                            
                        } else if(flag==-2) {
                            string tmpstr=snprs+'\n';
                            if(fputs_checked(tmpstr.c_str(),singleton))
                            {
                                LOGPRINTF("ERROR: in writing file %s .\n", singletonFName.c_str());
                                TERMINATE();
                            }
                            
                        } else if(flag==-3) {
                                printf("WARNING: The number of allele categories equals the sample size of SNP %s. vQTL of probe %s and SNP %s is skipped.\n", snprs.c_str(),prbid.c_str(),snprs.c_str());
                            
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
                            string tmpstr=snprs+'\n';
                            if(fputs_checked(tmpstr.c_str(),novar))
                            {
                                LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                TERMINATE();
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
                            string tmpstr=snprs+'\n';
                            if(fputs_checked(tmpstr.c_str(),novar))
                            {
                                LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                TERMINATE();
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
                            string tmpstr=snprs+'\n';
                            if(fputs_checked(tmpstr.c_str(),novar))
                            {
                                LOGPRINTF("ERROR: in writing file %s .\n", novarFName.c_str());
                                TERMINATE();
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
        else if(filetype==SMR_SPARSE_3)
        {
            LOGPRINTF("Would be released soon!\n")
        }
        //LOGPRINTF("vQTL tests cost: %f ms.\n",float( clock () - begin_time ) /  1000);
        fclose(besd);
        fclose(novar);
        fclose(singleton);
        if (tosmrflag) {
             LOGPRINTF("vQTL results in the cis-regions of %ld probes have been saved in sparse binary file %s.\n",einfo._epi_include.size(), besdName.c_str());
        } else {
            if(phenofileName!=NULL){
                LOGPRINTF("vQTL results of %ld SNPs have been saved in file %s.\n",write_count, besdName.c_str());
            } else {
                LOGPRINTF("vQTL results of %ld probes and %ld SNPs have been saved in file %s.\n",einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
            }
            
        }
       
    }
    int cis_eQTL_num(eInfo* einfo, bInfo* bdata, int cis_itvl, vector<int> &cis_num)
    {
        int maxr=0;
        cis_num.resize(einfo->_epi_include.size());
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int prbchr=einfo->_epi_chr[einfo->_epi_include[i]];
            int prbbp=einfo->_epi_bp[einfo->_epi_include[i]];
            int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
            int cisend=prbbp+cis_itvl*1000;
            int r=0, r1=0;
            for(int kk=0;kk<bdata->_include.size();kk++)
            {
                int snpchr=bdata->_chr[bdata->_include[kk]];
                int snpbp=bdata->_bp[bdata->_include[kk]];
                if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend) r++;
                if(snpchr==prbchr && snpbp>=cisstart && snpbp<=prbbp) r1++;
            }
            cis_num[i]=r;
            if(r>maxr) maxr=r;
        }
        return maxr;
    }
    long est_cis_eQTL_num(eInfo* einfo, bInfo* bdata, int cis_itvl)
    {
        srand((int)time(0));
        long idx = 0, ncount =0, testnum=einfo->_epi_include.size()>10?10:einfo->_epi_include.size();
        while(idx++ < testnum) {
            int r = rand() % einfo->_epi_include.size();
            int prbchr=einfo->_epi_chr[einfo->_epi_include[r]];
            int prbbp=einfo->_epi_bp[einfo->_epi_include[r]];
            int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
            int cisend=prbbp+cis_itvl*1000;
            r=0;
            for(int kk=0;kk<bdata->_include.size();kk++)
            {
                int snpchr=bdata->_chr[bdata->_include[kk]];
                int snpbp=bdata->_bp[bdata->_include[kk]];
                if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend) r++;
            }
            if(r>ncount) ncount=r;
        }
        return ncount;
    }
    void eQTL(char* outFileName,  char* efileName, char* befileName,  char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,char* covfileName, char* qcovfileName, bool tosmrflag, bool nofastlinear,bool cis_flag,int cis_itvl)
    {
        //default: OSCA format
        setNbThreads(thread_num);
        eInfo einfo;
        bInfo bdata;
        char* phenofileName=NULL;
        init_einfo(&einfo);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, phenofileName, bFileName,transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num,snplstName,snplst2exclde,tsk_ttl,tsk_id,covfileName,qcovfileName); // using _keep and _eii_include, the individuals are aligned.
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
        //LOGPRINTF("\nThe results would be saved in dense BESD format ...\n");
        if(outFileName!=NULL){
            write_smr_esi(outFileName, &bdata);
            write_smr_epi(outFileName, &einfo);
        }
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=OSCA_DENSE_1;
        if(tosmrflag && !cis_flag) filetype=SMR_DENSE_3;
        else if(!tosmrflag && cis_flag) filetype=OSCA_SPARSE_1;
        else if(tosmrflag && cis_flag) filetype=SMR_SPARSE_3;
        
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
        MatrixXd XtX_i;
        XtX_i=_Cov.transpose()*_Cov;
        bool determinant_zero=false;
        inverse_V(XtX_i, determinant_zero);
        if(determinant_zero)
        {
            LOGPRINTF("The matrix is not invertible. please check the correlation of the covariates.\n");
            TERMINATE();
        }
        if(!nofastlinear) fast_adjprobe(&einfo,_Cov, XtX_i);
        
        int slide_wind=10000;
        int loops=ceil(1.0*bdata._include.size()/slide_wind);
        
        vector<double> rst;
        if(filetype==SMR_DENSE_3)
        {
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                printf("%3.0f%%\r", 100.0*jj/einfo._epi_include.size());
                fflush(stdout);
                
                VectorXd y(einfo._eii_include.size());
                if(!nofastlinear) for(int ll=0; ll<einfo._eii_include.size(); ll++) y(ll)=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]];
                vector<float> beta,se;
                for(int ii=0;ii<loops;ii++)
                {
                    MatrixXd _X;
                    int snpstart=ii*slide_wind;
                    if(!nofastlinear)
                    {
                        make_XMat(&bdata,snpstart,slide_wind, _X,true);
                        //clock_t begin_time = clock();
                        fast_getResidual(_X, _Cov, XtX_i);
                       //LOGPRINTF(" cost: %f ms.\n",float( clock () - begin_time ) /  1000);
                        for(int kk=0;kk<_X.cols();kk++)
                        {
                            VectorXd x=_X.col(kk);
                            rst.clear();
                            adjusted_reg(y,x,rst,_X_c-1);
                            beta.push_back(rst[0]);
                            se.push_back(rst[1]);
                        }
                    }
                    else
                    {
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
                            rst.clear();
                            lin2(y, X, rst);
                            if(rst.size()!=3) {
                                LOGPRINTF("ERROR: bugs in function of linear regression found. Please report this.\n");
                                TERMINATE();
                            }
                            beta.push_back(rst[0]);
                            se.push_back(rst[1]);
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
            bool warned = false;
            for(int ii=0;ii<loops;ii++) // one SNP one record
            {
                MatrixXd _X;
                int snpstart=ii*slide_wind;
                if(!nofastlinear)
                {
                    make_XMat(&bdata,snpstart,slide_wind, _X,true);
                    fast_getResidual(_X, _Cov, XtX_i);
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        printf("%3.0f%%\r", 100.0*(ii*slide_wind+kk)/bdata._include.size());
                        fflush(stdout);
                         int snpid=snpstart+kk;
                        string snprs=bdata._snp_name[bdata._include[snpid]];
                        double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                        if(snpfreq==0 || snpfreq==1) {
                            if(!warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                        }
                        vector<float> betases;
                        VectorXd x=_X.col(kk);
                        for(int jj=0;jj<einfo._epi_include.size();jj++)
                        {
                           VectorXd y(einfo._eii_include.size());
                           for(int ll=0; ll<einfo._eii_include.size(); ll++) y(ll)=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]];
                           rst.clear();
                           adjusted_reg(y,x,rst,_X_c-1);
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
                else
                {
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
                            rst.clear();
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
        }
        else if(filetype==OSCA_SPARSE_1)
        {
            /*
            bool warned = false;
            vector<uint64_t> cols;
            vector<uint32_t> rowid;
            vector<float> vals;
            long ncount=est_cis_eQTL_num(&einfo, &bdata, cis_itvl);
            vals.reserve(ncount*einfo._epi_include.size()*2);
            rowid.reserve(ncount*einfo._epi_include.size());
            cols.resize(einfo._epi_include.size()+1);
            cols[0]=0;
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                printf("%3.0f%%\r", 100.0*jj/einfo._epi_include.size());
                fflush(stdout);
                
                MatrixXd _X;
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
                    //LOGPRINTF("%ld SNPs are included for the cis-region of probe %s [%d:%d,%d:%d].\n",snpids.size(),prbid.c_str(),prbchr,cisstart,prbchr,cisend);
                } else {
                    //LOGPRINTF("No SNP is included for the cis-region of probe %s [%d:%d,%d:%d].\n",prbid.c_str(),prbchr,cisstart,prbchr,cisend);
                    cols[jj+1]=cols[jj];
                    continue;
                }
               
                if(!nofastlinear)
                {
                    VectorXd y(einfo._eii_include.size());
                    for(int ll=0; ll<einfo._eii_include.size(); ll++) y(ll)=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]];
                    make_XMat(&bdata,snpids, _X,true);
                    fast_getResidual(_X, _Cov, XtX_i);
                    se.reserve(_X.cols());
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        uint32_t snpid=snpids[kk];
                        double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                        if(snpfreq==0 || snpfreq==1) {
                            if(!warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                            continue;
                        }
                        VectorXd x=_X.col(kk);
                        rst.clear();
                        adjusted_reg(y,x,rst,_X_c-1);
                        rowid.push_back(snpid);
                        vals.push_back(rst[0]);
                        se.push_back(rst[1]);
                    }
                }
                else
                {
                    make_XMat(&bdata,snpids, _X);
                    for(int kk=0;kk<_X.cols();kk++) //_X.cols() ==snpids.size()
                    {
                        uint32_t snpid=snpids[kk];
                        string snprs=bdata._snp_name[bdata._include[snpid]];
                        double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                        if(snpfreq==0 || snpfreq==1) {
                            if(!warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                            continue;
                        }
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
                        rst.clear();
                        lin2(y, X, rst);
                        if(rst.size()!=3) {
                            LOGPRINTF("ERROR: bugs in function of linear regression found. Please report this.\n");
                            TERMINATE();
                        }
                        rowid.push_back(snpid);
                        vals.push_back(rst[0]);
                        se.push_back(rst[1]);
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
             */
            
            vector< vector<uint32_t> > rowids(einfo._epi_include.size());
            vector< vector<float> > betas(einfo._epi_include.size());
            vector< vector<float> > ses(einfo._epi_include.size());
            vector<int> cis_num;
            cis_eQTL_num(&einfo,&bdata,cis_itvl,cis_num);
            for(int ii=0;ii<einfo._epi_include.size();ii++)
            {
                rowids[ii].reserve(cis_num[ii]);
                betas[ii].reserve(cis_num[ii]);
                ses[ii].reserve(cis_num[ii]);
                cis_num[ii]=0;
            }
            bool warned = false;
            VectorXd x, y(einfo._eii_include.size());
            MatrixXd _X;
            for(int ii=0;ii<loops;ii++) // one SNP one record
            {
                int snpstart=ii*slide_wind;
                if(!nofastlinear)
                {
                    make_XMat(&bdata,snpstart,slide_wind, _X,true);
                    fast_getResidual(_X, _Cov, XtX_i);
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        printf("%3.0f%%\r", 100.0*(ii*slide_wind+kk)/bdata._include.size());
                        fflush(stdout);
                        
                        int snpid=snpstart+kk;
                        int snpchr=bdata._chr[bdata._include[snpid]];
                        int snpbp=bdata._bp[bdata._include[snpid]];
                        double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                        if(snpfreq==0 || snpfreq==1) {
                            if(!warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                            continue;
                        }
                        x=_X.col(kk);
                        for(int jj=0;jj<einfo._epi_include.size();jj++)
                        {
                            int prbchr=einfo._epi_chr[einfo._epi_include[jj]];
                            int prbbp=einfo._epi_bp[einfo._epi_include[jj]];
                            int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
                            int cisend=prbbp+cis_itvl*1000;
                            if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend)
                            {
                                for(int ll=0; ll<einfo._eii_include.size(); ll++) y(ll)=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]];
                                rst.clear();
                                adjusted_reg(y,x,rst,_X_c-1);
                                rowids[jj].push_back(snpid);
                                betas[jj].push_back(rst[0]);
                                ses[jj].push_back(rst[1]);
                            }
                        }
                    }
                }
                else
                {
                    make_XMat(&bdata,snpstart,slide_wind, _X);
                    for(int kk=0;kk<_X.cols();kk++)
                    {
                        printf("%3.0f%%\r", 100.0*(ii*slide_wind+kk)/bdata._include.size());
                        fflush(stdout);
                        int snpid=snpstart+kk;
                        string snprs=bdata._snp_name[bdata._include[snpid]];
                        int snpchr=bdata._chr[bdata._include[snpid]];
                        int snpbp=bdata._bp[bdata._include[snpid]];
                        double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                        if(snpfreq==0 || snpfreq==1) {
                            if(!warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                            continue;
                        }
                        for(int jj=0;jj<einfo._epi_include.size();jj++)
                        {
                            int prbchr=einfo._epi_chr[einfo._epi_include[jj]];
                            int prbbp=einfo._epi_bp[einfo._epi_include[jj]];
                            int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
                            int cisend=prbbp+cis_itvl*1000;
                            if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend)
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
                                rst.clear();
                                lin2(y, X, rst);
                                if(rst.size()!=3) {
                                    LOGPRINTF("ERROR: bugs in function of linear regression found. Please report this.\n");
                                    TERMINATE();
                                }
                                rowids[jj].push_back(snpid);
                                betas[jj].push_back(rst[0]);
                                ses[jj].push_back(rst[1]);
                            }
                        }
                    }
                }
            }
            uint64_t valNum=0;
            vector<uint64_t> cols;
            cols.resize(einfo._epi_include.size()+1);
            cols[0]=0;
            for(int jj=0;jj<betas.size();jj++)
            {
                long co=2*betas[jj].size();
                valNum+=co;
                cols[jj+1]=cols[jj]+co;
            }
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
            for(int jj=0;jj<rowids.size();jj++)
            {
                if(rowids[jj].size())
                    if (fwrite_checked(&rowids[jj][0],rowids[jj].size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
            }
            for(int jj=0;jj<betas.size();jj++)
            {
                if(betas[jj].size())
                {
                    if (fwrite_checked(&betas[jj][0],betas[jj].size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&ses[jj][0],ses[jj].size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                }
            }
        }
        else if(filetype==SMR_SPARSE_3)
        {
            LOGPRINTF("Would be released soon!\n")
        }
        fclose(besd);
        LOGPRINTF("eQTL summary statistics for %ld probes and %ld SNPs are saved in file %s.\n", einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
       // LOGPRINTF("%ld probes cost: %f ms.\n",einfo._epi_include.size(),float( clock () - begin_time ) /  1000);
    }
    
    void eQTL_MLM(char* outFileName,  char* efileName, char* befileName,  char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,char* covfileName, char* qcovfileName, bool tosmrflag,bool cis_flag,int cis_itvl,char* grm_file,bool grm_bin_flag,bool no_constrain,int reml_mtd,int MaxIter,bool mlma_adj_covar)
    {
        setNbThreads(thread_num);
        eInfo einfo;
        bInfo bdata;
        char* phenofileName=NULL;
        init_einfo(&einfo);
        vector<string> grm_id, grm_files;
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, phenofileName,bFileName,transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num,snplstName,snplst2exclde,tsk_ttl,tsk_id,covfileName,qcovfileName, grm_file);
        if (maf > 0) filter_snp_maf(&bdata, maf);
        
        einfo._r_indx.clear();
        vector<MatrixXd> _A;
        vector<int> kp;
        vector<string> uni_id;
        map<string, int>::iterator iter;
        for(int i=0; i<bdata._keep.size(); i++)
            uni_id.push_back(bdata._fid[bdata._keep[i]]+":"+bdata._pid[bdata._keep[i]]);
        
        int _n=(int)bdata._keep.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        if(grm_file!=NULL)
        {
            grm_files.push_back(grm_file);
            read_grm(&bdata,grm_bin_flag, grm_file, grm_id, false, true);
            for(int i=0; i < grm_files.size() + 1; i++) einfo._r_indx.push_back(i);
            _A.resize(einfo._r_indx.size());
                match(uni_id, grm_id, kp);
                (_A[0]).resize(_n, _n);
                #pragma omp parallel for
                for(int i=0; i<_n; i++){
                    for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=bdata._grm(kp[i],kp[j]);
                }
                bdata._grm.resize(0,0);
        }
        else
        {
            grm_files.push_back("NA");
            make_grm(&bdata,0);
            for(int i=0; i < grm_files.size() + 1; i++) einfo._r_indx.push_back(i);
            _A.resize(einfo._r_indx.size());
            for(int i=0; i<einfo._eii_include.size(); i++) grm_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
            match(uni_id, grm_id, kp);
            (_A[0]).resize(_n, _n);
            #pragma omp parallel for
            for(int i=0; i<_n; i++){
                for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=bdata._grm(kp[i],kp[j]);
            }
        }
        _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
        SelfAdjointEigenSolver<MatrixXd> eigensolver(_A[0]);
        VectorXd eval = eigensolver.eigenvalues();
        MatrixXd U=eigensolver.eigenvectors();
        char outputname[FNAMESIZE];
        outputname[0]='\0';
        if(tsk_ttl>1) {
            if(outFileName!=NULL) {
                string tmp=  string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outFileName=outputname;
            }
        }
        LOGPRINTF("\nPerforming eQTL analysis with Mixed Linear Model ...\n");
        if(outFileName!=NULL){
            write_smr_esi(outFileName, &bdata);
            write_smr_epi(outFileName, &einfo);
        }
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=OSCA_DENSE_1;
        if(tosmrflag && !cis_flag) filetype=SMR_DENSE_3;
        else if(!tosmrflag && cis_flag) filetype=OSCA_SPARSE_1;
        else if(tosmrflag && cis_flag) filetype=SMR_SPARSE_3;
        
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
        
        //construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _Cov;
        int _X_c;
        _X_c=construct_X(&einfo, E_float, qE_float,_Cov);
        MatrixXd XtX_i;
        XtX_i=_Cov.transpose()*_Cov;
        bool determinant_zero=false;
        inverse_V(XtX_i, determinant_zero);
        if(determinant_zero)
        {
            LOGPRINTF("The matrix is not invertible. please check the correlation of the covariates.\n");
            TERMINATE();
        }
        //if(mlma_adj_covar) fast_adjprobe(&einfo,_Cov, XtX_i);
        
        int slide_wind=10000;
        vector<double> rst;
        if(filetype==SMR_DENSE_3)
        {
            LOGPRINTF("Would be released soon!\n")
        }
        else if(filetype==OSCA_DENSE_1)
        {
            LOGPRINTF("Would be released soon!\n")
        }
        else if(filetype==OSCA_SPARSE_1)
        {
            int omp_num_threads=1;
        #ifndef __APPLE__
            omp_num_threads=omp_get_max_threads();
        #endif
            vector< vector<uint32_t> > rowids(einfo._epi_include.size());
            vector< vector<float> > betas(einfo._epi_include.size());
            vector< vector<float> > ses(einfo._epi_include.size());
            vector<int> cis_num;
            int maxr=cis_eQTL_num(&einfo,&bdata,cis_itvl,cis_num);
            int maxr4=maxr*3;
            slide_wind=slide_wind>maxr4?slide_wind:maxr4;
            for(int ii=0;ii<einfo._epi_include.size();ii++)
            {
                rowids[ii].resize(cis_num[ii]);
                betas[ii].resize(cis_num[ii]);
                ses[ii].resize(cis_num[ii]);
                cis_num[ii]=0;
            }
            bool warned = false;
            MatrixXd Covx(_n,_X_c+1);
            Covx.block(0,0,_n,_X_c)=_Cov;
            if(omp_num_threads>2)
            {
                #pragma omp parallel for
                for(int jj=0;jj<einfo._epi_include.size();jj++)
                {
                    printf("%3.0f%%\r", 100.0*jj/einfo._epi_include.size());
                    fflush(stdout);
                    MatrixXd _X;
                    string prbid=einfo._epi_prb[einfo._epi_include[jj]];
                    int prbchr=einfo._epi_chr[einfo._epi_include[jj]];
                    int prbbp=einfo._epi_bp[einfo._epi_include[jj]];
                    int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
                    int cisend=prbbp+cis_itvl*1000;
                    vector<uint32_t> snpids;
                    for(int kk=0;kk<bdata._include.size();kk++)
                    {
                        int snpchr=bdata._chr[bdata._include[kk]];
                        int snpbp=bdata._bp[bdata._include[kk]];
                        if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend) snpids.push_back(kk);
                    }
                    if(snpids.size()==0) continue;
                    VectorXd y(einfo._eii_include.size());
                    for(int ll=0; ll<einfo._eii_include.size(); ll++) y(ll)=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[ll]];
                    MatrixXd _Vi;
                    VectorXd _b,_se;
                    vector<double> reml_priors,reml_priors_var;
                    mute=true;
                    reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,_Cov, y,_A, U, eval, _Vi,  reml_mtd,  MaxIter,_b,_se);
                    if(remlstatus==0 || remlstatus==-5 || (remlstatus==-3))
                    {
                        make_XMat_subset(&bdata,snpids, _X,false);
                        int realcisc=0;
                        bool missing=false;
                        for(int kk=0;kk<_X.cols();kk++)
                        {
                            uint32_t snpid=snpids[kk];
                            double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                            if(snpfreq==0 || snpfreq==1) {
                                if(!mute && !warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                                missing=true;
                                continue;
                            }
                            VectorXd y_buf=y;
                            Covx.col(_X_c)=_X.col(kk);
                            VectorXd x_buf=_X.col(kk);
                            double beta, se, pval;
                            if(mlma_adj_covar) mlm_stat(y_buf, _Vi, x_buf, beta, se, pval);
                            else mlm_stat_covar(y_buf, _Vi, Covx, beta, se, pval);
                            rowids[jj][realcisc]=snpid;
                            betas[jj][realcisc]=beta;
                            ses[jj][realcisc++]=se;
                        }
                        if(missing)
                        {
                            rowids[jj].resize(realcisc);
                            betas[jj].resize(realcisc);
                            ses[jj].resize(realcisc);
                        }
                    }
                    else
                    {
                        if(!mute) printf("WARNING: OREML failed in probe %s.\n ", prbid.c_str());
                        rowids[jj].clear();
                        betas[jj].clear();
                        ses[jj].clear();
                    }
                }
            }
            else
            {
                int recodc=0;
                vector<int> remain(einfo._epi_include), start, end;
                vector< vector<int> > idx;
                while(recodc<bdata._include.size())
                {
                    int snpstart=recodc-maxr>0?recodc-maxr:0;
                    int chrrecodc=bdata._chr[bdata._include[recodc]];
                    int chrstart=bdata._chr[bdata._include[snpstart]];
                    if(chrrecodc!=chrstart) snpstart=recodc;
                    chrstart=bdata._chr[bdata._include[snpstart]];
                    int snpend=(snpstart+slide_wind)>(bdata._include.size()-1)?((int)bdata._include.size()-1):(snpstart+slide_wind);
                    int chrend=bdata._chr[bdata._include[snpend]];
                    while(chrend!=chrstart) chrend=bdata._chr[bdata._include[--snpend]];
                    int bpstart=bdata._bp[bdata._include[snpstart]];
                    int bpend=bdata._bp[bdata._include[snpend]];
                    vector<int> slct;
                    for(int i=0;i<remain.size();i++)
                    {
                        if(remain[i]>=0)
                        {
                            int prbchr=einfo._epi_chr[remain[i]];
                            int prbbp=einfo._epi_bp[remain[i]];
                            int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
                            if(cisstart<bpstart) cisstart=bpstart;
                            int cisend=prbbp+cis_itvl*1000;
                            if(snpend-snpstart<slide_wind) cisend=prbbp;
                            if(prbchr==chrstart && cisstart>=bpstart && cisend<=bpend)
                            {
                                    slct.push_back(remain[i]);
                                    remain[i]=-9;
                            }
                        }
                    }
                    if(slct.size()>0)
                    {
                        start.push_back(snpstart);
                        end.push_back(snpend);
                        idx.push_back(slct);
                    }
                    recodc=snpend+1;
                }
                #pragma omp parallel for
                for(int ii=0;ii<start.size();ii++)
                {
                    printf("%3.0f%%\r", 100.0*ii/start.size());
                    fflush(stdout);
                    int snpstart=start[ii];
                    int window=end[ii]-snpstart+1;
                    if(window<=0){LOGPRINTF("ERROR: please report to futao.zhang@imb.uq.edu.au.\n");TERMINATE();};
                    MatrixXd _X;
                    VectorXd y(einfo._eii_include.size());
                    make_XMat_subset(&bdata,snpstart,window, _X,false);
                    for(int i=0;i<idx[ii].size();i++)
                    {
                        string prbid=einfo._epi_prb[idx[ii][i]];
                        int prbchr=einfo._epi_chr[idx[ii][i]];
                        int prbbp=einfo._epi_bp[idx[ii][i]];
                        int cisstart=(prbbp-cis_itvl*1000>0)?(prbbp-cis_itvl*1000):0;
                        int cisend=prbbp+cis_itvl*1000;
                        for(int ll=0; ll<einfo._eii_include.size(); ll++) y(ll)=einfo._val[idx[ii][i]*einfo._eii_num+einfo._eii_include[ll]];
                        MatrixXd _Vi;
                        VectorXd _b,_se;
                        vector<double> reml_priors,reml_priors_var;
                        mute=true;
                        reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,_Cov, y,_A, U, eval, _Vi,  reml_mtd,  MaxIter,_b,_se);
                        if(remlstatus==0 || remlstatus==-5 || (remlstatus==-3))
                        {
                            VectorXd y_buf=y;
                            int realcisc=0;
                            bool missing=false;
                            if(mlma_adj_covar) y_buf=y.array()-(_Cov*_b).array(); // adjust phenotype for covariates
                            for(int kk=0;kk<_X.cols();kk++)
                            {
                                int snpid=snpstart+kk;
                                int snpchr=bdata._chr[bdata._include[snpid]];
                                int snpbp=bdata._bp[bdata._include[snpid]];
                                double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                                double beta, se, pval;
                                if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend)
                                {
                                    if(snpfreq==0 || snpfreq==1) {
                                        if(!mute && !warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                                        missing=true;
                                        continue;
                                    }
                                    Covx.col(_X_c)=_X.col(kk);
                                    VectorXd x_buf=_X.col(kk);
                                    if(mlma_adj_covar) mlm_stat(y_buf, _Vi, x_buf, beta, se, pval);
                                    else mlm_stat_covar(y_buf, _Vi, Covx, beta, se, pval);
                                    rowids[idx[ii][i]][realcisc]=snpid;
                                    betas[idx[ii][i]][realcisc]=beta;
                                    ses[idx[ii][i]][realcisc++]=se;
                                }
                            }
                            if(missing)
                            {
                                rowids[idx[ii][i]].resize(realcisc);
                                betas[idx[ii][i]].resize(realcisc);
                                ses[idx[ii][i]].resize(realcisc);
                            }
                        }
                        else
                        {
                            if(!mute) printf("WARNING: OREML failed in probe %s.\n ", prbid.c_str());
                            rowids[idx[ii][i]].clear();
                            betas[idx[ii][i]].clear();
                            ses[idx[ii][i]].clear();
                        }
                    }
                }
                
            }
            uint64_t valNum=0;
            vector<uint64_t> cols;
            cols.resize(einfo._epi_include.size()+1);
            cols[0]=0;
            for(int jj=0;jj<betas.size();jj++)
            {
                long co=2*betas[jj].size();
                valNum+=co;
                cols[jj+1]=cols[jj]+co;
            }
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
            for(int jj=0;jj<rowids.size();jj++)
            {
                if(rowids[jj].size())
                    if (fwrite_checked(&rowids[jj][0],rowids[jj].size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
            }
            for(int jj=0;jj<betas.size();jj++)
            {
                if(betas[jj].size())
                {
                    if (fwrite_checked(&betas[jj][0],betas[jj].size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&ses[jj][0],ses[jj].size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                }
            }
        }
        else if(filetype==SMR_SPARSE_3)
        {
            LOGPRINTF("Would be released soon!\n")
        }
        fclose(besd);
        LOGPRINTF("cis-eQTL summary statistics for %ld probes are saved in file %s.\n", einfo._epi_include.size(), besdName.c_str());
        // LOGPRINTF("%ld probes cost: %f ms.\n",einfo._epi_include.size(),float( clock () - begin_time ) /  1000);
    }
}

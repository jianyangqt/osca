//
//  l3_vqtl.cpp
//  osc
//
//  Created by Futao Zhang on 29/08/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#include "l3_vqtl.hpp"
#include <time.h>
#include <fstream>
#include <sstream>

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
    void indi_check(bInfo* bdata,eInfo* einfo, eInfo* eCov)
    {
        map<string, int>::iterator iter0, iter1;
        vector<string> indi_list;
        vector<int> bkeep,ekeep,ckeep;
        map<string,int> bmap,emap,cmap;

        for(int i=0;i<bdata->_keep.size();i++)
        {
            int idx=bdata->_keep[i];
            string iid=bdata->_fid[idx]+":"+bdata->_pid[idx];
            iter0 = einfo->_eii_map.find(iid);
            iter1 = eCov->_eii_map.find(iid);
            if (iter0 != einfo->_eii_map.end() && iter1 != eCov->_eii_map.end()) {
                indi_list.push_back(iid);
                bkeep.push_back(idx);
                ekeep.push_back(iter0->second);
                ckeep.push_back(iter1->second);
                bmap.insert(pair<string,int>(iid,idx));
                emap.insert(pair<string,int>(iid,iter0->second));
                cmap.insert(pair<string,int>(iid,iter1->second));
            }
        }

        if(indi_list.size()==0){
            LOGPRINTF("No individual in common.\n");
            TERMINATE();
        }
        bdata->_keep.swap(bkeep);
        einfo->_eii_include.swap(ekeep);
        eCov->_eii_include.swap(ckeep);
        bdata->_id_map.swap(bmap);
        einfo->_eii_map.swap(emap);
        eCov->_eii_map.swap(cmap);

        LOGPRINTF("%ld individuals in common are kept.\n",  einfo->_eii_include.size());
    }
    void indi_check(eInfo* edata,eInfo* einfo)
    {
        map<string, int>::iterator iter;
        vector<string> indi_list;
        vector<int> bkeep,ekeep;
        map<string,int> bmap,emap;

        for(int i=0;i<edata->_eii_include.size();i++)
        {
            int idx=edata->_eii_include[i];
            string iid=edata->_eii_fid[idx]+":"+edata->_eii_iid[idx];
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
        edata->_eii_include.swap(bkeep);
        einfo->_eii_include.swap(ekeep);
        edata->_eii_map.swap(bmap);
        einfo->_eii_map.swap(emap);

        LOGPRINTF("%ld individuals in common are kept.\n",  einfo->_eii_include.size());
    }
    void indi_check(eInfo* edata,eInfo* einfo, eInfo* eCov)
    {
        map<string, int>::iterator iter0, iter1;
        vector<string> indi_list;
        vector<int> bkeep,ekeep,ckeep;
        map<string,int> bmap,emap,cmap;

        for(int i=0;i<edata->_eii_include.size();i++)
        {
            int idx=edata->_eii_include[i];
            string iid=edata->_eii_fid[idx]+":"+edata->_eii_iid[idx];
            iter0 = einfo->_eii_map.find(iid);
            iter1 = eCov->_eii_map.find(iid);
            if (iter0 != einfo->_eii_map.end() && iter1 != eCov->_eii_map.end()) {
                indi_list.push_back(iid);
                bkeep.push_back(idx);
                ekeep.push_back(iter0->second);
                ckeep.push_back(iter1->second);
                bmap.insert(pair<string,int>(iid,idx));
                emap.insert(pair<string,int>(iid,iter0->second));
                cmap.insert(pair<string,int>(iid,iter1->second));
            }
        }

        if(indi_list.size()==0){
            LOGPRINTF("No individual in common.\n");
            TERMINATE();
        }
        edata->_eii_include.swap(bkeep);
        einfo->_eii_include.swap(ekeep);
        eCov->_eii_include.swap(ckeep);
        edata->_eii_map.swap(bmap);
        einfo->_eii_map.swap(emap);
        eCov->_eii_map.swap(cmap);

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

    void load_vqtl_workspace(eInfo* einfo, bInfo* bdata, char* efileName, \
        char* befileName, char* phenofileName, char* bFileName, bool transposed, \
        int efileType,char* problstName,char* problst2exclde,char* genelistName, \
        int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind, \
        int fromprbkb, int toprbkb,bool prbwindFlag, char* genename, \
        char* probe2exclde,char* indilstName,char* indilst2remove, \
        bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh, \
        double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, \
        double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, \
        int autosome_num,char* snplstName,char* snplst2exclde,int tsk_ttl, \
        int tsk_id, char* covfileName, char* qcovfileName,char* grm_file, \
        int xqtlNO, double zeroratio, eInfo* eCov, char* covbodfileName, \
        char* covefileName, bool transopse_ecov)
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
        char bcname[FNAMESIZE];
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
            if(covbodfileName!=NULL ||  covefileName!=NULL)
            {
                if(covbodfileName!=NULL)
                {
                    memcpy(bcname,covbodfileName,strlen(covbodfileName)+1);
                    char* suffix=bcname+strlen(covbodfileName);
                    memcpy(suffix,".oii",5);
                    read_eii(bcname,eCov);
                } else {
                    if(transopse_ecov) read_efile_t(covefileName,eCov,efileType,no_fid_flag,valueType);
                    else read_efile(covefileName,eCov,efileType,no_fid_flag,valueType);
                }
            }
            //commn individuals
            if(covbodfileName!=NULL || covefileName!=NULL) indi_check(bdata,einfo,eCov);
            else indi_check(bdata,einfo);
            if(covbodfileName!=NULL)
            {
                char* suffix=bcname+strlen(covbodfileName);
                memcpy(suffix,".opi",5);
                read_epi(bcname,eCov);
                memcpy(suffix,".bod",5);
                read_beed(bcname,eCov);
            }
            epi_man(einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            if(xqtlNO == 3) extract_sqtl_probe(einfo,tsk_ttl,  tsk_id);
            else if(tsk_ttl>1) extract_probe(einfo,  tsk_ttl,  tsk_id);
            if(dpvalfName!=NULL) filtering_with_detpval(einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(einfo, upperBeta, lowerBeta);
            if(xqtlNO == 3 && einfo->_valType ==0) filtering_with_zeroratio(einfo, zeroratio);
            if(mssratio_prob<1) filtering_with_missingratio(einfo, mssratio_prob);
            read_bimfile(bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(bdata, snplst2exclde);
            read_bedfile(bdata, string(bFileName)+".bed");

        //here
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
            if(covbodfileName!=NULL || covefileName!=NULL)
            {
                if(covbodfileName!=NULL)
                {
                    memcpy(bcname,covbodfileName,strlen(covbodfileName)+1);
                    char* suffix=bcname+strlen(covbodfileName);
                    memcpy(suffix,".oii",5);
                    read_eii(bcname,eCov);
                } else {
                    if(transopse_ecov) read_efile_t(covefileName,eCov,efileType,no_fid_flag,valueType);
                    else read_efile(covefileName,eCov,efileType,no_fid_flag,valueType);
                }
            }
            //commn individuals
            if(covbodfileName!=NULL || covefileName!=NULL)
                indi_check(bdata,einfo,eCov);
            else
                indi_check(bdata,einfo);
            if(covbodfileName!=NULL)
            {
                char* suffix=bcname+strlen(covbodfileName);
                memcpy(suffix,".opi",5);
                read_epi(bcname,eCov);
                memcpy(suffix,".bod",5);
                read_beed(bcname,eCov);
            }
            memcpy(suffix,".opi",5);
            read_epi(inputname,einfo);
            epi_man(einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            if(xqtlNO == 3)
                extract_sqtl_probe(einfo,tsk_ttl,  tsk_id);
            else if(tsk_ttl>1)
                extract_probe(einfo,  tsk_ttl,  tsk_id);
            memcpy(suffix,".bod",5);
            clock_t begin_time = clock();
            read_beed(inputname,einfo);
            //LOGPRINTF("read_bod: %f ms.\n",float( clock () - begin_time ) /  1000);
            if(dpvalfName!=NULL) filtering_with_detpval(einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(einfo, upperBeta, lowerBeta);
            if(mssratio_prob<1) filtering_with_missingratio(einfo, mssratio_prob);
            if(xqtlNO == 3 && einfo->_valType ==0) filtering_with_zeroratio(einfo, zeroratio);
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


     void V_QTL(char* outFileName,  char* efileName, char* befileName, char* phenofileName, char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,int vqtl_mtd,char* covfileName, char* qcovfileName, bool tosmrflag,bool cis_flag, int cis_itvl)
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
        FILE* ma=NULL;
        string besdName=string(outFileName)+".besd";
        string maName = string(outFileName)+".ma";
        string outstr="";
        if(phenofileName!=NULL)
        {
            besdName=string(outFileName)+".vqtl";
            if(fopen_checked(&besd, besdName.c_str(),"w")) TERMINATE();
            outstr="Chr\tSNP\tA1\tA2\tfreq\tbp";
            if(vqtl_mtd==0 || vqtl_mtd==3 ) outstr+="\tstatistic\tdf\tbeta\tse\tP\tNMISS\n";
            else outstr+="\tF-statistic\tdf1\tdf2\tbeta\tse\tP\tNMISS\n";
            if(fputs_checked(outstr.c_str(),besd))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if(fopen_checked(&ma, maName.c_str(),"w")) TERMINATE();
            outstr="SNP\tA1\tA2\tfreq\tbeta\tse\tP\tNMISS\n";
            if(fputs_checked(outstr.c_str(),ma))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", maName.c_str());
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
                    string snpa1 = bdata._allele1[bdata._include[snpstart+kk]];
                    string snpa2 = bdata._allele2[bdata._include[snpstart+kk]];
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
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t'+ atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }
                                if(phenofileName!=NULL && ma) {
                                    string outstr = snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),ma))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", maName.c_str());
                                        TERMINATE();
                                    }
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
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[5]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + atos(rst[2]) + '\t'  + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }
                                if(phenofileName!=NULL && ma) {
                                    string outstr = snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),ma))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", maName.c_str());
                                        TERMINATE();
                                    }
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
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[5]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + atos(rst[2]) + '\t'  + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }
                                if(phenofileName!=NULL && ma) {
                                    string outstr = snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),ma))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", maName.c_str());
                                        TERMINATE();
                                    }
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
                                    string outstr = atos(snpchr) + '\t' + snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(snpbp) + '\t' + atos(rst[3]) + '\t' + atos(rst[4]) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),besd))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", besdName.c_str());
                                        TERMINATE();
                                    }
                                    write_count++;
                                }
                                if(phenofileName!=NULL && ma) {
                                    string outstr = snprs + '\t' +snpa1 + '\t' +snpa2+ '\t' + atos(snpfreq) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                    if(fputs_checked(outstr.c_str(),ma))
                                    {
                                        LOGPRINTF("ERROR: in writing file %s .\n", maName.c_str());
                                        TERMINATE();
                                    }
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
        if(ma!=NULL) fclose(ma);
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


    int cis_eQTL_num_2(eInfo* einfo, eqtlInfo* eqtlinfo, int cis_itvl, vector< vector<int>> &tranids, vector< vector<int> > &snpids, vector<int> &cis_num)
    {
        //bp as start and gd as end
        int maxr=0;
        cis_num.resize(einfo->_epi_include.size());
        snpids.resize(einfo->_epi_include.size());
        vector<int>::iterator it;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int prbchr=einfo->_epi_chr[einfo->_epi_include[i]];
            int prbstart=einfo->_epi_bp[einfo->_epi_include[i]];
            int prbend=einfo->_epi_gd[einfo->_epi_include[i]];
            int cisstart=(prbstart-cis_itvl*1000>0)?(prbstart-cis_itvl*1000):0;
            int cisend=prbend+cis_itvl*1000;
            vector< vector<int>> tmpsnpid(tranids[i].size());
            for(int kk=0;kk<tranids[i].size();kk++)
            {
                int prbid=tranids[i][kk];
                long end=eqtlinfo->_cols[prbid+1];
                long start=eqtlinfo->_cols[prbid];
                uint64_t num=(end-start)>>1;
                uint64_t rowpos=start>>1;
                for(int j=0;j<num;j++)
                {
                    int esi_id = eqtlinfo->_rowid[rowpos+j];
                    int snpbp = eqtlinfo->_esi_bp[esi_id];
                    int snpchr = eqtlinfo->_esi_chr[esi_id];
                    if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend)
                        tmpsnpid[kk].push_back(esi_id);
                }
                stable_sort(tmpsnpid[kk].begin(),tmpsnpid[kk].end());
            }
            //when cis, snp is intersection
            vector<int> intersection(tmpsnpid[0].size());
            for(int kk=1;kk<tranids[i].size();kk++)
            {
                it=set_intersection(tmpsnpid[kk].begin(),tmpsnpid[kk].end(), tmpsnpid[0].begin(), tmpsnpid[0].end(), intersection.begin());
                intersection.resize(it-intersection.begin());
                intersection.swap(tmpsnpid[0]);
            }
            snpids[i].swap(tmpsnpid[0]);
            cis_num[i]=(int)snpids[i].size();
            if(cis_num[i]>maxr) maxr=cis_num[i];
        }
        return maxr;
    }

    int trans_eQTL_num(eInfo* einfo, eqtlInfo* eqtlinfo, int trans_itvl, vector< vector<int>> &tranids, vector< vector<int> > &snpids, vector<int> &trans_num)
    {
        //bp as start and gd as end
        int maxr=0;
        trans_num.resize(einfo->_epi_include.size());
        snpids.resize(einfo->_epi_include.size());
        vector<int>::iterator it;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int prbchr=einfo->_epi_chr[einfo->_epi_include[i]];
            int prbstart=einfo->_epi_bp[einfo->_epi_include[i]];
            int prbend=einfo->_epi_gd[einfo->_epi_include[i]];
            int trans_start=(prbstart-trans_itvl*1000>0)?(prbstart-trans_itvl*1000):0;
            int trans_end=prbend+trans_itvl*1000;
            vector< vector<int>> tmpsnpid(tranids[i].size());
            for(int kk=0;kk<tranids[i].size();kk++)
            {
                int prbid=tranids[i][kk];
                long end=eqtlinfo->_cols[prbid+1];
                long start=eqtlinfo->_cols[prbid];
                uint64_t num=(end-start)>>1;
                uint64_t rowpos=start>>1;
                for(int j=0;j<num;j++)
                {
                    int esi_id = eqtlinfo->_rowid[rowpos+j];
                    int snpbp = eqtlinfo->_esi_bp[esi_id];
                    int snpchr = eqtlinfo->_esi_chr[esi_id];
                    if(snpchr != prbchr || snpbp <= trans_start || snpbp >= trans_end)
                        tmpsnpid[kk].push_back(esi_id);
                }
                stable_sort(tmpsnpid[kk].begin(),tmpsnpid[kk].end());
            }
            //when trans, snp is union.
            map<int, int> union_map;
            map<int, int>::iterator mit;
            vector<int> union_my;
            for(int kk = 1; kk < tranids[i].size(); kk++)
            {
                for (int ll = 0; ll < tmpsnpid[kk].size(); ll++) {
                    union_map.insert(pair<int, int>(tmpsnpid[kk][ll], 1));
                }
            }
            for (mit = union_map.begin(); mit != union_map.end(); mit++) {
                union_my.push_back(mit -> first);
            }
            snpids[i].swap(union_my);
            trans_num[i]=(int)snpids[i].size();
            if(trans_num[i]>maxr) maxr=trans_num[i];
        }
        return maxr;
    }


    void make_bs(eqtlInfo* eqtlinfo, vector<int> &prids, vector<int> &snpids, MatrixXd &eqtlb, MatrixXd &eqtls, vector<float> &freqs)
    {
        freqs.resize(snpids.size());
        map<int ,int > snp_map;
        for(int i=0;i<snpids.size();i++) snp_map.insert(pair<int,int>(snpids[i],i));
        map<int ,int >::iterator it;

        string filename=string(outfileName)+"."+atos(eqtlinfo->_epi_gene[prids[0]])+".txt";
        FILE* tmpfile= NULL;
        if(loud)
        {
            tmpfile=fopen(filename.c_str(),"w");
            if(!tmpfile)
            {
                printf("ERROR: open file %s.\n",filename.c_str());
                exit(EXIT_FAILURE);
            }
        }

            for(int kk=0;kk<prids.size();kk++)
            {
                int prbid=prids[kk];
                int prbchr= eqtlinfo->_epi_chr[prbid];
                string prbID = eqtlinfo->_epi_prbID[prbid];
                string gene = eqtlinfo->_epi_gene[prbid];
                int bp = eqtlinfo->_epi_bp[prbid];

                long end=eqtlinfo->_cols[prbid+1];
                long start=eqtlinfo->_cols[prbid];
                uint64_t num=(end-start)>>1;
                uint64_t rowpos=start>>1;
                for(int j=0;j<num;j++)
                {
                    int esi_id = eqtlinfo->_rowid[rowpos+j];
                    string snprs = eqtlinfo->_esi_rs[esi_id];
                    float snpfrq = eqtlinfo->_esi_freq[esi_id];
                    if(snpfrq==-9)
                    {
                        LOGPRINTF("ERROR: one or more frequencies are missing, please update the .esi file.\n");
                        exit(EXIT_FAILURE);
                    }
                    it = snp_map.find(esi_id);
                    if(it==snp_map.end()) continue;
                    int idx=it->second;
                    double beta=eqtlinfo->_val[start+j];
                    double se=eqtlinfo->_val[start+j+num];
                    eqtlb(idx,kk)=beta;
                    eqtls(idx,kk)=se;
                    if(kk==0) freqs[idx]=snpfrq;
                    if(loud)
                    {
                        string str=prbID + '\t'+ snprs +'\t'+atos(snpfrq)+'\t'+atos(beta)+'\t'+atos(se)+'\n';
                        fputs(str.c_str(),tmpfile);
                    }
                }
            }
        if(loud) fclose(tmpfile);
    }


    int cis_eQTL_num_2(eInfo* einfo, bInfo* bdata, int cis_itvl, vector< vector<uint32_t> > &snpids, vector<int> &cis_num)
    {
        //bp as start and gd as end
        int maxr=0;
        cis_num.resize(einfo->_epi_include.size());
        snpids.resize(einfo->_epi_include.size());
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int prbchr=einfo->_epi_chr[einfo->_epi_include[i]];
            int prbstart=einfo->_epi_bp[einfo->_epi_include[i]];
            int prbend=einfo->_epi_gd[einfo->_epi_include[i]];
            int cisstart=(prbstart-cis_itvl*1000>0)?(prbstart-cis_itvl*1000):0;
            int cisend=prbend+cis_itvl*1000;
            int r=0;
            for(int kk=0;kk<bdata->_include.size();kk++)
            {
                int snpchr=bdata->_chr[bdata->_include[kk]];
                int snpbp=bdata->_bp[bdata->_include[kk]];
                if(snpchr==prbchr && snpbp>=cisstart && snpbp<=cisend)
                {
                    snpids[i].push_back(kk);
                    r++;
                }
            }
            cis_num[i]=r;
            if(r>maxr) maxr=r;
        }
        return maxr;
    }

    int trans_eQTL_num(eInfo* einfo, bInfo* bdata, int trans_itvl, vector< vector<uint32_t> > &snpids, vector<int> &trans_num)
    {
        //bp as start and gd as end
        int maxr=0;
        trans_num.resize(einfo -> _epi_include.size());
        snpids.resize(einfo -> _epi_include.size());
        for(int i = 0; i < einfo -> _epi_include.size(); i++) {
            int prbchr = einfo -> _epi_chr[einfo -> _epi_include[i]];
            int prbstart = einfo -> _epi_bp[einfo -> _epi_include[i]];
            int prbend = einfo -> _epi_gd[einfo -> _epi_include[i]];
            int trans_left = (prbstart - trans_itvl * 1000 > 0)? (prbstart - trans_itvl * 1000): 0;
            int trans_right = prbend + trans_itvl * 1000;
            int r = 0;
            for(int kk = 0; kk < bdata -> _include.size(); kk++) {
                int snpchr = bdata -> _chr[bdata -> _include[kk]];
                int snpbp = bdata -> _bp[bdata -> _include[kk]];
                if(snpchr != prbchr || snpbp <= trans_left || snpbp >= trans_right) {
                    snpids[i].push_back(kk);
                    r++;
                }
            }
            trans_num[i] = r;
            if(r > maxr) {
                maxr = r;
            }
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
            LOGPRINTF("Enter OSCA eQTL SMR_SPARSE_3\n");
            int probe_num = einfo._epi_include.size();
            vector < vector <uint32_t> > rowids;
            vector < vector <float> > betas, ses;
            vector < uint64_t > beta_se_offset;
            beta_se_offset.push_back(0);
            uint64_t beta_se_offset_head = 0;

            uint32_t snp_num = bdata._include.size();
            uint64_t num_beta_se = 0;
            vector <uint32_t> snp_idx_per_probe;
            uint32_t probe_chr = 0;
            int64_t probe_bp = 0;
            uint32_t probe_bp_start = 0;
            uint32_t probe_bp_end = 0;
            uint32_t snp_chr = 0;
            uint32_t snp_pos = 0;
            uint32_t snp_ok_num = 0;
            for (int jj = 0; jj < probe_num; jj++)
            {
                printf("%3.0f%%\r", 100.0 * jj / einfo._epi_include.size());
                fflush(stdout);

                probe_chr = einfo._epi_chr[einfo._epi_include[jj]];
                probe_bp = (int64_t)einfo._epi_bp[einfo._epi_include[jj]];
                probe_bp_start = (probe_bp - cis_itvl * 1000) > 0? probe_bp - cis_itvl * 1000: 1;
                probe_bp_end = probe_bp + cis_itvl * 1000;
                
                snp_idx_per_probe.clear();
                for (uint32_t snp_i = 0; snp_i < snp_num; snp_i++) {
                    snp_chr = bdata._chr[bdata._include[snp_i]];
                    snp_pos = bdata._bp[bdata._include[snp_i]];
                    if (snp_chr == probe_chr && snp_pos >= probe_bp_start && snp_pos <= probe_bp_end) {
                        snp_idx_per_probe.push_back(snp_i);
                    }

                }
                snp_ok_num = snp_idx_per_probe.size();
                num_beta_se += snp_ok_num * 2;
                rowids.push_back(snp_idx_per_probe);
                beta_se_offset_head += snp_ok_num;
                beta_se_offset.push_back(beta_se_offset_head);
                beta_se_offset_head += snp_ok_num;
                beta_se_offset.push_back(beta_se_offset_head);

                VectorXd y(einfo._eii_include.size());
                if (!nofastlinear)
                    for (int ll = 0; ll < einfo._eii_include.size(); ll++)
                        y(ll) = einfo._val[einfo._epi_include[jj] * einfo._eii_num + einfo._eii_include[ll]];
                vector<float> beta, se;

                
                MatrixXd _X;
                if (!nofastlinear)
                {
                    make_XMat(&bdata, snp_idx_per_probe, _X, true);
                    // clock_t begin_time = clock();
                    fast_getResidual(_X, _Cov, XtX_i);
                    // LOGPRINTF(" cost: %f ms.\n",float( clock () - begin_time ) /  1000);
                    for (int kk = 0; kk < _X.cols(); kk++)
                    {
                        VectorXd x = _X.col(kk);
                        rst.clear();
                        adjusted_reg(y, x, rst, _X_c - 1);
                        beta.push_back(rst[0]);
                        se.push_back(rst[1]);
                    }
                    betas.push_back(beta);
                    ses.push_back(se);
                }
                else
                {
                    make_XMat(&bdata, snp_idx_per_probe, _X);
                    for (int kk = 0; kk < _X.cols(); kk++)
                    {
                        double nonmiss = 0.0;
                        int miss = 0;
                        MatrixXd X(_Cov.rows(), _Cov.cols() + 1);
                        X.block(0, 0, _Cov.rows(), _Cov.cols()) = _Cov;
                        long x_idx = X.cols() - 1;
                        VectorXd y(einfo._eii_include.size());

                        for (int ll = 0; ll < einfo._eii_include.size(); ll++)
                        {
                            if (einfo._eii_fid[einfo._eii_include[ll]] != bdata._fid[bdata._keep[ll]])
                                printf("ERROR: Alignment failed. Please report this bug to Futao < futao.zhang@imb.uq.edu.au >. Thanks.\n");
                            double val = einfo._val[einfo._epi_include[jj] * einfo._eii_num + einfo._eii_include[ll]]; // einfo._eii_include[ll] should equals ll. because of update_eii all the time after reading beed file.
                            double bval = _X(ll, kk);
                            if (val < 1e9 && bval < 1e5)
                            {
                                y(ll - miss) = val;
                                X(ll - miss, x_idx) = bval;
                                nonmiss += 1.0;
                            }
                            else
                            {
                                removeRow(X, ll - miss);
                                miss++;
                            }
                        }
                        if (miss > 0)
                            y.conservativeResize(einfo._eii_include.size() - miss);
                        rst.clear();
                        lin2(y, X, rst);
                        if (rst.size() != 3)
                        {
                            LOGPRINTF("ERROR: bugs in function of linear regression found. Please report this.\n");
                            TERMINATE();
                        }
                        beta.push_back(rst[0]);
                        se.push_back(rst[1]);
                    }
                    betas.push_back(beta);
                    ses.push_back(se);
                }
            }

            if (fwrite_checked(&num_beta_se, 1 * sizeof(uint64_t), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }

            if (fwrite_checked(&beta_se_offset[0], beta_se_offset.size() * sizeof(uint64_t), besd))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }

            for (uint32_t probe_cnt = 0; probe_cnt < probe_num; probe_cnt++)
            {
                if (fwrite_checked(&rowids[probe_cnt][0], sizeof(uint32_t) * rowids[probe_cnt].size(), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                if (fwrite_checked(&rowids[probe_cnt][0], sizeof(uint32_t) * rowids[probe_cnt].size(), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
            }

            for (uint32_t probe_cnt = 0; probe_cnt < probe_num; probe_cnt++)
            {
                if (fwrite_checked(&betas[probe_cnt][0], sizeof(float) * betas[probe_cnt].size(), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                if (fwrite_checked(&ses[probe_cnt][0], sizeof(float) * ses[probe_cnt].size(), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
            }
        }
        fclose(besd);
        LOGPRINTF("eQTL summary statistics for %ld probes and %ld SNPs are saved in file %s.\n", einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
       // LOGPRINTF("%ld probes cost: %f ms.\n",einfo._epi_include.size(),float( clock () - begin_time ) /  1000);
    }

    void eQTL_MLM(char* outFileName,  char* efileName, char* befileName,  char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,char* covfileName, char* qcovfileName, bool tosmrflag,bool cis_flag,int cis_itvl,char* grm_file,bool grm_bin_flag,bool no_constrain,int reml_mtd,int MaxIter,bool nopreadj_covar)
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
            string filena=string(outFileName)+".failed.probe.list";
            FILE* tmpfi=fopen(filena.c_str(),"w");
            if(!tmpfi)
            {
                LOGPRINTF("error open file.\n");
                TERMINATE();
            }

            //if(omp_num_threads>2)
            //{
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
                    VectorXd y_buf=y;
                    MatrixXd Covx(_n,_X_c+1);
                    Covx.block(0,0,_n,_X_c)=_Cov;
                    reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,_Cov, y,_A, U, eval, _Vi,  reml_mtd,  MaxIter,_b,_se);
                    if(remlstatus==0 || remlstatus==-5 || (remlstatus==-3))
                    {
                        if(!nopreadj_covar) y_buf=y.array()-(_Cov*_b).array();
                        make_XMat_subset(&bdata,snpids, _X,false);
                        int realcisc=0;
                        bool missing=false;
                        for(int kk=0;kk<_X.cols();kk++)
                        {
                            uint32_t snpid=snpids[kk];
                            double snpfreq=bdata._mu[bdata._include[snpid]]/2;
                            if(snpfreq==0 || snpfreq==1) {
                                if(loud && !warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                                missing=true;
                                continue;
                            }
                            Covx.col(_X_c)=_X.col(kk);
                            VectorXd x_buf=_X.col(kk);
                            double beta, se, pval;
                            if(!nopreadj_covar) mlm_stat(y_buf, _Vi, x_buf, beta, se, pval);
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
                        if(loud) printf("WARNING: OREML failed in probe %s.\n ", prbid.c_str());
                        rowids[jj].clear();
                        betas[jj].clear();
                        ses[jj].clear();
                        string str = prbid +'\n';
                        #pragma omp critical
                        fputs(str.c_str(),tmpfi);

                    }
                }
            fclose(tmpfi);
            /* }
            else
            {
                // irfahan found this part can exclude some probes from analysis.
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
                    MatrixXd Covx(_n,_X_c+1);
                    Covx.block(0,0,_n,_X_c)=_Cov;
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
                        reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,_Cov, y,_A, U, eval, _Vi,  reml_mtd,  MaxIter,_b,_se);
                        if(remlstatus==0 || remlstatus==-5 || (remlstatus==-3))
                        {
                            VectorXd y_buf=y;
                            int realcisc=0;
                            bool missing=false;
                            if(!no_preadj_covar) y_buf=y.array()-(_Cov*_b).array(); // adjust phenotype for covariates
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
                                        if(loud && !warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s). These results would be labeled as missing.\n"); warned=1;}
                                        missing=true;
                                        continue;
                                    }
                                    Covx.col(_X_c)=_X.col(kk);
                                    VectorXd x_buf=_X.col(kk);
                                    if(!no_preadj_covar) mlm_stat(y_buf, _Vi, x_buf, beta, se, pval);
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
                            if(loud) printf("WARNING: OREML failed in probe %s.\n ", prbid.c_str());
                            rowids[idx[ii][i]].clear();
                            betas[idx[ii][i]].clear();
                            ses[idx[ii][i]].clear();
                        }
                    }
                }

            }
            */
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

    void read_annofile(eqtlInfo* eqtlinfo, char* bedFileName)
    {
        FILE* epifile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum=6;
        if(fopen_checked(&epifile, bedFileName,"r")) TERMINATE();
        LOGPRINTF("Reading annotation information from %s ...\n", bedFileName);
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        eqtlinfo->_snpNum = 0 ;
        eqtlinfo->_probNum = 0 ;
        eqtlinfo->_valNum = 0 ;
        eqtlinfo->_sampleNum = 0 ;
        bool chrwarning=false;
        bool genewarning=false;
        bool orienwarning=false;
        while(fgets(Tbuf, MAX_LINE_SIZE, epifile))
        {
            split_str(Tbuf,strlist,0);
            if(Tbuf[0]=='\0') {
                LOGPRINTF("ERROR: Line %u is blank.\n", line_idx);
                TERMINATE();
            }
            if(strlist.size()>colnum)
            {
                //LOGPRINTF("WARNING: Line %u has more than %d items. The first %d columns would be used. \n", line_idx,colnum,colnum);
            }
            if(strlist.size()<colnum)
            {
                LOGPRINTF("ERROR: Line %u has less than %d items.\n", line_idx,colnum);
                TERMINATE();
            }
            to_upper(strlist[3]);
            eqtlinfo->_probe_name_map.insert(pair<string,int>(strlist[3], line_idx));
            if(eqtlinfo->_probe_name_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate probe : %s.\n", strlist[1].c_str());
                TERMINATE();
            }
            if(strlist[0]=="X" || strlist[0]=="x") eqtlinfo->_epi_chr.push_back(23);
            else if(strlist[0]=="Y" || strlist[0]=="y") eqtlinfo->_epi_chr.push_back(24);
            else if(strlist[0]=="NA" || strlist[0]=="na"){
                eqtlinfo->_epi_chr.push_back(-9);
                if(!chrwarning) {
                    LOGPRINTF("WARNING: At least one probe chr is missing.\n");
                    chrwarning=true;
                }
            } else if (atoi(strlist[0].c_str())==0 ) {
                //LOGPRINTF("ERROR: unrecongized chromosome found:\n");
                //LOGPRINTF("WARNING: unrecongized chromosome found. This chromosome is set to 0:\n");
                //LOGPRINTF("%s\n",Tbuf);
                //TERMINATE();
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            } else if ( atoi(strlist[0].c_str())>24 || atoi(strlist[0].c_str())<0) {
                //LOGPRINTF("WARNING: abmormal chromosome found:\n");
                //LOGPRINTF("%s\n",Tbuf);
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            } else eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));

            if(strlist[1]=="NA" || strlist[1]=="na" || strlist[2]=="NA" || strlist[2]=="na") {
                LOGPRINTF("ERROR: NA start / end position found:\n");
                LOGPRINTF("%s\n",Tbuf);
                TERMINATE();
            }
            int startp = atoi(strlist[1].c_str());
            int endp = atoi(strlist[2].c_str());
            int pos = startp + (endp-startp)/2;
            eqtlinfo->_epi_bp.push_back(startp);
            eqtlinfo->_epi_prbID.push_back(strlist[3]);
            eqtlinfo->_epi_gd.push_back(endp);

            if(strlist[5]=="NA" || strlist[5]=="na") {
                if(!genewarning) {
                    LOGPRINTF("WARNING: at least one gene id is missing.\n");
                    genewarning=true;
                }
            }
            eqtlinfo->_epi_gene.push_back(strlist[5].c_str());
            if(strlist[4]=="NA") {
                eqtlinfo->_epi_orien.push_back('*');
                if(!orienwarning) {
                    LOGPRINTF("WARNING: At least one gene strand is missing.\n");
                    orienwarning=true;
                }
            } else eqtlinfo->_epi_orien.push_back(strlist[4][0]);
            eqtlinfo->_include.push_back(line_idx);
            line_idx++;
        }
        eqtlinfo->_probNum =line_idx;
        fclose(epifile);
        LOGPRINTF("%lu probes to be included from  %s .\n", eqtlinfo->_probNum, bedFileName);
    }

    void gene_check(eInfo* sqtlinfo,vector< vector<int>> &tranids, char* annofileName,eqtlInfo* eqtlinfo)
    {
        eqtlInfo tmpinfo;
        if(annofileName != NULL) read_annofile(&tmpinfo, annofileName);

        int ids = 0;
        map<int, int> echr_map;
        map<int, int>::iterator iter1;
        for(int i=0;i<eqtlinfo->_include.size();i++)
            echr_map.insert(pair<int,int>(eqtlinfo->_epi_chr[eqtlinfo->_include[i]],ids++));

        map<string, int> gene_count_map;
        map<string, int>::iterator iter;
        vector<string> gene;
        vector< vector< int> > idx;
        ids = 0;
        for(int i=0; i<eqtlinfo->_include.size();i++)
        {
            string gn = eqtlinfo->_epi_gene[eqtlinfo->_include[i]];
            to_upper(gn);
            iter = gene_count_map.find(gn);
            if(iter == gene_count_map.end())
            {
                gene_count_map.insert(pair<string, int>(gn, ids++));
                gene.push_back(gn);
                vector<int> tmpidx;
                tmpidx.push_back(eqtlinfo->_include[i]);
                idx.push_back(tmpidx);
            }
            else
            {
                int curid = iter->second;
                idx[curid].push_back(eqtlinfo->_include[i]);
            }
        }
        vector<int> inids;
        int tmpc=0;
        for(int i=0;i<idx.size();i++)
            if(idx[i].size()>1)
            {
                iter = tmpinfo._probe_name_map.find(gene[i]);
                if(iter!=tmpinfo._probe_name_map.end())
                {
                    int tid=iter->second;
                    int tchr = tmpinfo._epi_chr[tid];
                    iter1 = echr_map.find(tchr);
                    if(iter1 != echr_map.end())
                    {
                        sqtlinfo->_epi_prb.push_back(tmpinfo._epi_prbID[tid]);
                        sqtlinfo->_epi_chr.push_back(tchr);
                        sqtlinfo->_epi_gd.push_back(tmpinfo._epi_gd[tid]);
                        sqtlinfo->_epi_gene.push_back(tmpinfo._epi_gene[tid]);
                        sqtlinfo->_epi_orien.push_back(tmpinfo._epi_orien[tid]);
                        sqtlinfo->_epi_bp.push_back(tmpinfo._epi_bp[tid]);
                        sqtlinfo->_epi_include.push_back(tmpc++);
                        inids.push_back(i);
                    }
                }
            }

        for(int i=0;i<inids.size();i++)
        {
            vector<int> tmp;
            tmp.swap(idx[inids[i]]);
            tranids.push_back(tmp);
        }
        LOGPRINTF("%ld genes are retained after gene check with annaotation information and eQTL summary statistics.\n",tranids.size());
    }

    void gene_check(eInfo* sqtlinfo,vector< vector<int>> &tranids, char* annofileName, bInfo* bdata,eInfo* einfo, bool cis_flag, bool trans_flag)
    {
        eqtlInfo tmpinfo;
        if (annofileName) {
            read_annofile(&tmpinfo, annofileName);
        } else {
            LOGPRINTF("a annotation file is needed.\n");
        }

        int ids = 0;
        map<int, int> bchr_map;
        map<int, int>::iterator iter1;
        for(int i=0;i<bdata->_include.size();i++)
            bchr_map.insert(pair<int,int>(bdata->_chr[bdata->_include[i]],ids++));

        map<string, int> gene_count_map;
        map<string, int>::iterator iter;
        vector<string> gene;
        vector< vector< int> > idx;
        ids = 0;
        for(int i=0; i<einfo->_epi_include.size();i++)
        {
            string gn = einfo->_epi_gene[einfo->_epi_include[i]];
            to_upper(gn);
            iter = gene_count_map.find(gn);
            if(iter == gene_count_map.end())
            {
                gene_count_map.insert(pair<string, int>(gn, ids++));
                gene.push_back(gn);
                vector<int> tmpidx;
                tmpidx.push_back(einfo->_epi_include[i]);
                idx.push_back(tmpidx);
            }
            else
            {
                int curid = iter->second;
                idx[curid].push_back(einfo->_epi_include[i]);
            }
        }
        vector<int> inids;
        int tmpc=0;
        for(int i = 0; i < idx.size(); i++)
            if(idx[i].size()>1)
            {
                iter = tmpinfo._probe_name_map.find(gene[i]);
                if(iter!=tmpinfo._probe_name_map.end())
                {
                    int tid=iter->second;
                    int tchr = tmpinfo._epi_chr[tid];
                    iter1 = bchr_map.find(tchr);
                    if(cis_flag && !trans_flag) {
                        if (iter1 != bchr_map.end()) {
                            //printf("inter cis mode\n");
                            sqtlinfo->_epi_prb.push_back(tmpinfo._epi_prbID[tid]);
                            sqtlinfo->_epi_chr.push_back(tchr);
                            sqtlinfo->_epi_gd.push_back(tmpinfo._epi_gd[tid]);
                            sqtlinfo->_epi_gene.push_back(tmpinfo._epi_gene[tid]);
                            sqtlinfo->_epi_orien.push_back(tmpinfo._epi_orien[tid]);
                            sqtlinfo->_epi_bp.push_back(tmpinfo._epi_bp[tid]);
                            sqtlinfo->_epi_include.push_back(tmpc++);
                            inids.push_back(i);
                        }
                    } else if (trans_flag && !cis_flag) {
                        //printf("inter trans mode\n");
                        sqtlinfo->_epi_prb.push_back(tmpinfo._epi_prbID[tid]);
                        sqtlinfo->_epi_chr.push_back(tchr);
                        sqtlinfo->_epi_gd.push_back(tmpinfo._epi_gd[tid]);
                        sqtlinfo->_epi_gene.push_back(tmpinfo._epi_gene[tid]);
                        sqtlinfo->_epi_orien.push_back(tmpinfo._epi_orien[tid]);
                        sqtlinfo->_epi_bp.push_back(tmpinfo._epi_bp[tid]);
                        sqtlinfo->_epi_include.push_back(tmpc++);
                        inids.push_back(i);

                    } else {
                        fprintf(stderr, "cis or trans must choosed and choose one\n");
                    }
                }
            }

        for(int i=0;i<inids.size();i++)
        {
            vector<int> tmp;
            tmp.swap(idx[inids[i]]);
            tranids.push_back(tmp);
        }
        LOGPRINTF("%ld genes are retained after gene check with annaotation information, genotype information, and gene expression information.\n",tranids.size());
    }

    void rankr(VectorXd &a, VectorXd &b)
    {
        b.resize(a.size());
        #pragma omp parallel for
        for (long i = a.size()-1; i >= 0; i--)
        {
            int count = 0 , ecount=0;
            double w=0;
            for (int j = 0; j < a.size(); j++) {
                if (a[j] < a[i]) count++;
                if (a[j] == a[i]) ecount++;
            }
            if(ecount>1)
            {
                for(int j=0;j<ecount;j++) w+=j;
                w/=ecount;
            }
            b[i] =1 + count + w;
        }
    }


static double
get_var_mean(vector <double>& se, MatrixXd& cor_null)
{
    double part1 = 0;
    double part2 = 0;
    double dt_out = 0;
    int i = 0, j = 0;
    int vector_len = se.size();
    for (i = 0; i < vector_len; i++) {
        part1 += se[i] * se[i];
    }
    for (i = 0; i < (vector_len - 1); i++) {
        for (j = (i + 1); j < vector_len; j++) {
            part2 += 2 * se[i] * se[j] * cor_null(i, j);
        }
    }

    dt_out = (part1 + part2) / (vector_len * vector_len);
    return dt_out;
}


static double
get_cov_beta_mean(int target, vector <double>& se, MatrixXd& cor_null)
{
    double dt_out = 0;
    int vector_len = se.size();
    int i = 0;
    for (i = 0; i < vector_len; i++) {
        dt_out += se[target] * se[i] * cor_null(target, i);
    }

    dt_out = dt_out / vector_len;
    return dt_out;
}


/*
    for every snp of a probe, there are severl transctripts(isoform),
    print snprs(snp name), probe name,  isoform name, beta of this isoform, se of
    this isoform.
    isoform name, beta, se will repeat number of isoform times.
    --Benjamin Fang

    This founction was abandoned by fanghl. 20220212
*/
#include <stdlib.h>
static void
output_beta_se(const double pdev, string snprs, string prbid, vector< double > beta,
    vector < double > se, vector< int > tranids_prb, eInfo & einfo,
    FILE * beta_se_of_trans, const double p_cutoff)
{

    struct ORDER_LIST {
        double beta;
        double se;
        int tranids_prb;
        struct ORDER_LIST * prev;
        struct ORDER_LIST * next;
    };
    if (pdev >= p_cutoff)
        return;

    int len = beta.size();
    if (len != se.size() || len != tranids_prb.size()) {
        fprintf(stderr, "error, length of beta se and tranids should equal, snp %s passed\n", snprs.c_str());
        return;
    }

    int i = 0, j = 0;
    struct ORDER_LIST * head = NULL, * new_ele = NULL, * ptr = NULL, * tmp = NULL;
    for (i = 0; i < len; i++) {
        new_ele = (struct ORDER_LIST *)malloc(sizeof(struct ORDER_LIST));
        new_ele -> beta = beta[i];
        new_ele -> se = se[i];
        new_ele -> tranids_prb = tranids_prb[i];
        new_ele -> prev = NULL;
        new_ele -> next = NULL;
        if (head) {
            ptr -> next = new_ele;
            new_ele -> prev = ptr;
            ptr = new_ele;
        } else {
            head = ptr = new_ele;
        }
    }
    struct ORDER_LIST * order_head = NULL;
    order_head = head;
    head = head -> next;
    head -> prev = NULL;
    order_head -> next = NULL;
    char finder_marker = 0;
    while(head) {
        //printf("%lf\n", head -> beta);
        finder_marker = 0;
        tmp = head -> next;
        head -> prev = NULL;
        head -> next = NULL;
        ptr = order_head;
        while (true) {
            if (head -> beta >= ptr -> beta) {
                if (ptr -> prev) {
                    ptr -> prev -> next = head;
                    head -> prev = ptr -> prev;
                } else {
                    order_head = head;
                }
                head -> next = ptr;
                ptr -> prev = head;
                finder_marker = 1;
                break;
            }
            if (ptr -> next)
                ptr = ptr -> next;
            else
                break;
        }
        if (!finder_marker) {
            ptr -> next = head;
            head -> prev = ptr;
        }
        head = tmp;
    }
    fprintf(beta_se_of_trans, "%s\t%s", snprs.c_str(), prbid.c_str());
    if (len <= 10) {
        while(order_head) {
            fprintf(beta_se_of_trans, "\t%s\t%lf\t%lf", (einfo._epi_prb[order_head -> tranids_prb]).c_str(), order_head -> beta, order_head -> se);
            order_head = order_head -> next;
        }
        fprintf(beta_se_of_trans, "\n");
    } else {
        i = 0;
        j = len - 5;
        while(order_head) {
            if (i < 5 || i >= j) {
                fprintf(beta_se_of_trans, "\t%s\t%lf\t%lf", (einfo._epi_prb[order_head -> tranids_prb]).c_str(), order_head -> beta, order_head -> se);
            }
            order_head = order_head -> next;
        }
        fprintf(beta_se_of_trans, "\n");

    }

    return;
}


    static void
    convert_chrom_num_to_string(int chr_num, char * chr_string)
    {
        // This function convert a number to string, 23 will be convert to X,
        // 24 --> Y, 25 --> Mt, other except 1-22, will be convert to NA

        if (chr_num == 23) {
            strcpy(chr_string, "X");
        } else if (chr_num == 24) {
            strcpy(chr_string, "Y");
        } else if (chr_num == 25) {
            strcpy(chr_string, "Mt");
        } else if (chr_num < 1 || chr_num > 25) {
            fprintf(stderr, "Chromsome %d not recgnized, will convert it to NA", chr_num);
            strcpy(chr_string, "NA");
        } else {
            sprintf(chr_string, "%d", chr_num);
        }

        return;
    }


    void sQTL(char* outFileName, char* efileName, char* befileName,
        char* bFileName, bool transposed, int efileType, char* problstName,
        char* problst2exclde, char* genelistName, int chr, char* prbname,
        char* fromprbname, char* toprbname, int prbWind, int fromprbkb,
        int toprbkb, bool prbwindFlag, char* genename, char* probe2exclde,
        char* indilstName, char* indilst2remove, bool no_fid_flag, int valueType,
        bool beta2m, bool m2beta, double std_thresh, double upperBeta,
        double lowerBeta, char* dpvalfName, double dp_thresh, double prb_thresh,
        double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num,
        double maf, char* snplstName, char* snplst2exclde, int tsk_ttl,
        int tsk_id, char* covfileName, char* qcovfileName, bool tosmrflag,
        bool nofastlinear, bool cis_flag, int cis_itvl, double zeroratio, double call,
        char* annofileName, char* covbodfileName, char* covefileName,
        bool transopse_ecov, bool use_top_p, bool trans_flag, int trans_itvl, bool no_isoform_eqtl, bool not_use_top)
    {
        LOGPRINTF("\n\033[0;32mEnter sQTL analysis...\033[0m\n");
        if (cis_flag && trans_flag) {
            fprintf(stderr, "The cis_flag and trans_flag can not be true at same time.");
            TERMINATE();
        }
        //by default, will set cis_flag = true
        if (!cis_flag && !trans_flag) {
            cis_flag = true;
        }
        if (!use_top_p && !not_use_top) {
            use_top_p = true;
        } else if (use_top_p && not_use_top) {
            not_use_top = false;
        }
        setNbThreads(thread_num);
        LOGPRINTF("Using %d thread(s) to conduct analysis ...\n", thread_num);

        eInfo einfo;
        bInfo bdata;
        eInfo eCov;
        char * phenofileName = NULL;
        int xqtlNO = 3;
        init_einfo(&einfo);
        init_einfo(&eCov);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, phenofileName,
            bFileName,transposed, efileType, problstName, problst2exclde, genelistName,
            chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb,
            prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag,
            valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName,
            dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob, autosome_num,
            snplstName, snplst2exclde, tsk_ttl, tsk_id, covfileName, qcovfileName,
            NULL, xqtlNO, zeroratio, &eCov, covbodfileName, covefileName,
            transopse_ecov); // using _keep and _eii_include, the individuals are aligned.
        if(maf > 0)
            filter_snp_maf(&bdata, maf);
        if(call > 0)
            filter_snp_call(&bdata, call);
        char outputname[FNAMESIZE];
        outputname[0] = '\0';
        if(tsk_ttl > 1) {
            if(outFileName != NULL) {
                string tmp=  string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outFileName=outputname;
            }
        }

        eInfo sqtlinfo;
        vector< vector<int>> tranids;  //smaller einfo._epi_include
        gene_check( &sqtlinfo,tranids, annofileName, &bdata, &einfo, cis_flag, trans_flag);

        LOGPRINTF("\nPerforming sQTL analysis ...\n");
        if(outFileName!=NULL){
            write_smr_esi(outFileName, &bdata);
            write_smr_epi(outFileName, &sqtlinfo,true);
        }
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=OSCA_SPARSE_1;
        if(tosmrflag) filetype=SMR_SPARSE_3;

        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=(int)bdata._keep.size();
        ten_ints[2]=(int)bdata._include.size();
        ten_ints[3]=(int)sqtlinfo._epi_include.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
            TERMINATE();
        }
        if(covfileName!=NULL || qcovfileName!=NULL) adjprobe(&einfo);

        vector< vector<uint32_t> > rowids(sqtlinfo._epi_include.size());
        vector< vector<float> > betas(sqtlinfo._epi_include.size());
        vector< vector<float> > ses(sqtlinfo._epi_include.size());
        vector< vector<uint32_t> > snpids; //to save _include id not _include value
        if (cis_flag) {
            vector<int> cis_num;
            cis_eQTL_num_2(&sqtlinfo, &bdata, cis_itvl, snpids, cis_num);

            for(int ii=0;ii<sqtlinfo._epi_include.size();ii++)
            {
                rowids[ii].resize(cis_num[ii]);
                betas[ii].resize(cis_num[ii]);
                ses[ii].resize(cis_num[ii]);
            }
        } else if (trans_flag) {
            vector <int> trans_num;
            trans_eQTL_num(&sqtlinfo, &bdata, trans_itvl, snpids, trans_num);
            for (int i = 0; i < sqtlinfo._epi_include.size(); i++) {
                rowids[i].resize(trans_num[i]);
                betas[i].resize(trans_num[i]);
                ses[i].resize(trans_num[i]);

            }
        }

        /*
        FILE * f_out = fopen("dtobs", "w");
        for (int i = 0; i < rowids.size(); i++) {
            fprintf(f_out, "%d\n", rowids[i].size());
        }
        fclose(f_out);
        exit(0);
        */

        bool warned = false;
        int nindi = (int)einfo._eii_include.size();
        double cr=0.0;
        FILE * fmid_data = NULL;
        if (!no_isoform_eqtl) {
            char snp_data_file[1024] = "";
            if (strlen(outFileName) > 1000) {
                fprintf(stderr, "file name too long, using snp_data.txt instead\n");
                strcpy(snp_data_file, "isoform_eQTL_effect.txt");
            } else {
                strcpy(snp_data_file, outFileName);
                strcat(snp_data_file, "_isoform_eQTL_effect.txt");
            }
            
            fmid_data = fopen(snp_data_file, "w");
            if (!fmid_data) {
                fprintf(stderr, "open snp_data.txt failed.\n");
                exit(1);
            }
            fprintf(fmid_data,
                "#SNP\tChr\tBP\tA1\tA2\tFreq\t"
                "Probe\tProbe_chr\tgene\tOritation\tb\tse\tp\n"
                );
        }
        
        int probe_num_ok = sqtlinfo._epi_include.size();        
        int probe_counter = 0;
        #pragma omp parallel for private(cr)
        for(int jj = 0; jj < sqtlinfo._epi_include.size(); jj++)
        {
            
            string prbid = sqtlinfo._epi_prb[sqtlinfo._epi_include[jj]];
            string gene_name = sqtlinfo._epi_gene[sqtlinfo._epi_include[jj]];
            int prb_chr = sqtlinfo._epi_chr[sqtlinfo._epi_include[jj]];
            char chrom_probe[8];
            convert_chrom_num_to_string(prb_chr, chrom_probe);
            char prb_ori = sqtlinfo._epi_orien[sqtlinfo._epi_include[jj]];
            int numTrans = (int)tranids[jj].size();
            vector <int> transcripts_idx = tranids[jj];
            uint32_t snp_num = snpids[jj].size();
            vector <uint32_t> snp_idx = snpids[jj];

            #pragma omp critical
            {
                LOGPRINTF("\n\033[0;32m>\033[0mProcessing gene:%s, probe_id:%s (%d/%d).\n",
                          gene_name.c_str(), prbid.c_str(), ++probe_counter, probe_num_ok);
                LOGPRINTF("    This gene contains %d isoforms and %d SNPs.\n",
                          numTrans, snp_num);
            }

            double desti = 1.0 * jj / (sqtlinfo._epi_include.size() - 1);
            if(desti >= cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr == 0)
                    cr += 0.05;
                else if(cr == 0.05)
                    cr += 0.2;
                else if(cr == 0.25)
                    cr += 0.5;
                else
                    cr += 0.25;
            }
            
            if(snp_num == 0){
                printf("gene %s was passed, because no snp contained.\n", gene_name.c_str());
                continue;
            }

            MatrixXd _X;
            make_XMat(&bdata, snpids[jj], _X);

            vector<double> tpm(numTrans * nindi);
            VectorXd overall;
            vector<int> missidx;
            map<int,int> missidx_map;
            map<int, int>::iterator iter;
            map<string, int>::iterator citer;
            vector< vector<double>> trpv;

            trpv.resize(numTrans);
            for(int kk = 0; kk < numTrans; kk++)
                trpv[kk].resize(einfo._eii_include.size());
            for( int kk = 0; kk < numTrans; kk++)
            {
                for(int ll = 0; ll < einfo._eii_include.size(); ll++)
                    trpv[kk][ll] = einfo._val[tranids[jj][kk] * einfo._eii_num + \
                        einfo._eii_include[ll]];
            }

            vector< vector<double> > cor_null;
            cor_null.resize(numTrans);
            for(int kk = 0; kk < numTrans; kk++)
                cor_null[kk].resize(numTrans);

            for(int kk = 0; kk < numTrans; kk++)
            {
                for(int ll = kk + 1; ll < numTrans; ll++)
                {
                    vector<double> y, x;
                    for(int mm = 0; mm < einfo._eii_include.size(); mm++)
                    {
                        if(trpv[kk][mm] < 1e9 && trpv[ll][mm] < 1e9)
                        {
                            y.push_back(trpv[kk][mm]);
                            x.push_back(trpv[ll][mm]);
                        }
                    }
                    cor_null[kk][ll] = cor_null[ll][kk] = cor(y,x);
                }
            }

            //modified by fanghl
            //remove row and columns which contain value is 1
            
            vector < int > need_remove;
            int i = 0, j = 0, k = 0, l = 0;
            vector < int >::iterator it;

            for (i = 0; i < cor_null.size(); i++){
                for (j = 1 + i; j < cor_null[i].size(); j++){
                    if (abs(cor_null[i][j] - 1) < 1e-15){
                        need_remove.push_back(i);
                        break;
                    }
                }
            }

            vector <int> transcripts_idx_tmp;
            if (need_remove.size() > 0) {
                for (int i = 0; i < transcripts_idx.size(); i++) {
                    it = find(need_remove.begin(), need_remove.end(), i);
                    if (it == need_remove.end()) {
                        transcripts_idx_tmp.push_back(transcripts_idx[i]);
                    }
                }
                transcripts_idx = transcripts_idx_tmp;
            }

            numTrans -= need_remove.size();
            vector < vector <double> > cor_null_clean;
            vector < vector < double > > trpv_clean;
            vector < double > tmp;
            if (need_remove.size() > 0 && numTrans > 1){
                cor_null_clean.resize(numTrans);
                for (i = 0; i < numTrans; i++){
                    cor_null_clean[i].resize(numTrans);
                }

                k = 0;
                for (i = 0; i < cor_null.size(); i++){
                    l = 0;
                    it = find(need_remove.begin(), need_remove.end(), i);
                    if (it != need_remove.end()) {
                        k++;
                    } else {
                        for (j = 0; j < cor_null[i].size(); j++) {
                            it = find(need_remove.begin(), need_remove.end(), j);
                            if (it != need_remove.end()) {
                                l++;
                            } else {
                                cor_null_clean[i - k][j - l] = cor_null[i][j];
                            }
                        }
                    }
                }

                for(int kk = 0; kk < numTrans; kk++)
                    cor_null_clean[kk][kk] = 1.0;

                for (i = 0; i < trpv.size(); i++) {
                    tmp.clear();
                    it = find(need_remove.begin(), need_remove.end(), i);
                    if (it == need_remove.end()) {
                        for (j = 0; j < trpv[i].size(); j++){
                            tmp.push_back(trpv[i][j]);
                        }
                        trpv_clean.push_back(tmp);
                    }
                }
            } else {
                cor_null_clean = cor_null;
                for(int kk = 0; kk < numTrans; kk++)
                    cor_null_clean[kk][kk] = 1.0;
                trpv_clean = trpv;
            }

            //_X.cols() ==snpids[jj].size()
            for(int kk = 0; kk < _X.cols(); kk ++) {
                uint32_t snpid = snp_idx[kk];
                string snprs = bdata._snp_name[bdata._include[snpid]];
                char snp_chr = bdata._chr[bdata._include[snpid]];
                uint32_t snp_pos = bdata._bp[bdata._include[snpid]];
                string snp_a1 = bdata._allele1[bdata._include[snpid]];
                string snp_a2 = bdata._allele2[bdata._include[snpid]];
                double snpfreq = bdata._mu[bdata._include[snpid]] / 2;
                if (snpfreq == 0 || snpfreq == 1) {
                    if(!warned) {LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s).\n"); warned=1;}
                    rowids[jj][kk] = snpid;
                    betas[jj][kk] = 0;
                    ses[jj][kk] = 1;
                    continue;
                }
                vector<double> beta(numTrans), se(numTrans);
                for(int ll = 0; ll < numTrans; ll++)
                {
                    vector<double> y, x, rst;
                    for(int mm=0; mm < einfo._eii_include.size(); mm++)
                    {
                        double bval=_X(mm,kk), tval=trpv_clean[ll][mm];
                        if(bval < 1e5 && tval < 1e9)
                        {
                            y.push_back(tval);
                            x.push_back(bval);
                        }
                    }
                    reg(y, x, rst);

                    beta[ll]=rst[0];
                    se[ll]=rst[1];
                }

                int varnum = numTrans * (numTrans - 1) / 2;
                int i = 0, j = 0, k = 0;
                VectorXd d(varnum), vardev(varnum), chisq_dev(varnum);
                MatrixXd vdev(varnum, varnum);
                MatrixXd corr_dev(varnum, varnum);
                double beta_mean = 0;
                double var_mean = 0;
                double tmp1 = 0;
                double tmp2 = 0;
                if (use_top_p) {
                    d.resize(numTrans);
                    vardev.resize(numTrans);
                    chisq_dev.resize(numTrans);
                    vdev.resize(numTrans, numTrans);
                    corr_dev.resize(numTrans, numTrans);

                    for (i = 0; i < numTrans; i++) {
                        beta_mean += beta[i];
                    }
                    beta_mean = beta_mean / numTrans;

                    tmp1 = 0;
                    tmp2 = 0;
                    for (i = 0; i < numTrans; i++) {
                        tmp1 += se[i] * se[i];
                    }
                    for (i = 0; i < numTrans - 1; i++) {
                        for (j = i + 1; j < numTrans; j++) {
                            tmp2 += 2 * se[i] * se[j] * cor_null_clean[i][j];
                        }
                    }
                    var_mean = (tmp1 + tmp2) / (numTrans * numTrans);

                    for (i = 0; i < numTrans; i++) {
                        tmp1 = 0;
                        d[i] = beta[i] - beta_mean;
                        for (j = 0; j < numTrans; j++) {
                            tmp1 += se[i] * se[j] * cor_null_clean[i][j];
                        }
                        vardev[i] = se[i] * se[i] + var_mean - 2 * (tmp1 / numTrans);
                    }

                    for (i = 0; i < numTrans; i++) {
                        chisq_dev[i] = d[i] * d[i] / vardev[i];
                    }

                    for (i = 0; i < numTrans; i++) {
                        for (j = 0; j < numTrans; j++) {
                            tmp1 = 0;
                            tmp2 = 0;
                            for (k = 0; k < numTrans; k++){
                                tmp1 += se[i] * se[k] * cor_null_clean[i][k];
                            }
                            tmp1 = tmp1 / numTrans;

                            for (k = 0; k < numTrans; k++) {
                                tmp2 += se[j] * se[k] * cor_null_clean[i][k];
                            }
                            tmp2 = tmp2 / numTrans;

                            vdev(i, j) = se[i] * se[j] * cor_null_clean[i][j] - \
                                tmp1 - tmp2 + var_mean;
                        }
                    }

                    for (i = 0; i < numTrans; i++) {
                        for (j = i; j < numTrans; j++) {
                            corr_dev(i, j) = corr_dev(j, i) = \
                                vdev(i, j) / sqrt(vdev(i, i) * vdev(j, j));
                        }
                    }

                } else {
                    k = 0;
                    for(int m1 = 0; m1 < numTrans - 1; m1++) {
                        for(int m2 = m1 + 1; m2 < numTrans; m2++) {
                            d[k] = beta[m1] - beta[m2];
                            vardev[k] = se[m1] * se[m1] + se[m2] * se[m2] - \
                                2 * cor_null_clean[m1][m2] * se[m1] * se[m2];
                            k++;
                        }
                    }

                    for(int m1 = 0; m1 < varnum; m1++)
                        chisq_dev[m1] = d[m1] * d[m1] / vardev[m1];

                    int mi = 0, mj =0;
                    for( int m1 = 0; m1 < numTrans - 1; m1++) {
                        for( int m2 = m1 + 1; m2 < numTrans; m2++) {
                            mj = 0;
                            for(int m3 = 0; m3 < numTrans - 1; m3++){
                                for(int m4 = m3 + 1; m4 < numTrans; m4++) {
                                    vdev(mi,mj) = se[m1] * se[m3] * cor_null_clean[m1][m3] - \
                                        se[m1] * se[m4] * cor_null_clean[m1][m4] - \
                                        se[m2] * se[m3] * cor_null_clean[m2][m3] + \
                                        se[m2] * se[m4] * cor_null_clean[m2][m4];
                                    mj++;
                                }
                            }
                            mi++;
                        }
                    }

                    for( int m1 = 0; m1 < varnum; m1++) {
                        for( int m2 = m1; m2 < varnum; m2++){
                            corr_dev(m1,m2) = corr_dev(m2,m1) = vdev(m1, m2) / sqrt(vdev(m1, m1) * \
                                vdev(m2, m2));
                        }
                    }
                }

                VectorXd lambda;
                #pragma omp critical
                {
                    SelfAdjointEigenSolver<MatrixXd> es;
                    es.compute(corr_dev, EigenvaluesOnly);
                    lambda=es.eigenvalues();

                    double z = 0.0;

                    double sumChisq_dev = chisq_dev.sum();
                    double pdev = 0.0;
                    pdev = pchisqsum(sumChisq_dev,lambda);
                    z = sqrt(qchisq(pdev,1));

                    double beta_hat = z / sqrt(2 * snpfreq * (1 - snpfreq) * (nindi + z * z));
                    double se_hat = 1 / sqrt(2 * snpfreq * (1 - snpfreq) * (nindi + z * z));

                    rowids[jj][kk] = snpid;
                    betas[jj][kk] = beta_hat;
                    ses[jj][kk] = se_hat;

                    if (fmid_data) {
                        char chrom_snp[8] = "";
                        double snp_pval = 0;
                        double snp_p_tmp = 0;
                        convert_chrom_num_to_string(snp_chr, chrom_snp);

                        fprintf(fmid_data,
                                "%s\t%s\t%u\t%s\t%s\t%le\t"
                                "%s\t%s\t%s\t%c\t%le\t%le\t%le",
                                snprs.c_str(), chrom_snp, snp_pos, snp_a1.c_str(),
                                snp_a2.c_str(), snpfreq,
                                prbid.c_str(), chrom_probe, gene_name.c_str(), prb_ori,
                                beta_hat, se_hat, pdev
                                );
                        for (int i = 0; i < numTrans; i++) {
                            snp_p_tmp = beta[i] / se[i];
                            snp_p_tmp *= snp_p_tmp;
                            snp_pval = pchisq(snp_p_tmp, 1);
                            fprintf(fmid_data, 
                                "\t%s\t%le\t%le\t%le",
                                einfo._epi_prb[transcripts_idx[i]].c_str(), beta[i], se[i], snp_pval);
                        }
                        fprintf(fmid_data, "\n");
                    }
                }

            }
        }
        if (fmid_data) {
            fclose(fmid_data);
        }

        if(tosmrflag)
        {
            uint64_t valNum=0;
            vector<uint64_t> cols;
            cols.resize(2*sqtlinfo._epi_include.size()+1);
            cols[0]=0;
            for(int jj=0;jj<betas.size();jj++)
            {
                long co=betas[jj].size();
                valNum+=co;
                cols[2*jj+1]=cols[2*jj]+co;
                valNum+=co;
                cols[2*(jj+1)]=cols[2*jj+1]+co;
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
                {
                    if (fwrite_checked(&rowids[jj][0],rowids[jj].size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&rowids[jj][0],rowids[jj].size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
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
        else
        {
            uint64_t valNum=0;
            vector<uint64_t> cols;
            cols.resize(sqtlinfo._epi_include.size()+1);
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


        fclose(besd);
        LOGPRINTF("sQTL summary statistics for %ld probes and %ld SNPs are saved in file %s.\n", einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
    }


    bool pcc(MatrixXd& PCC, MatrixXd& buffer_beta, MatrixXd& buffer_se,
        double pmecs, int nmecs)
    {
        long snpnum = buffer_beta.rows();
        long cohortnum = buffer_beta.cols();
        //pearson correlation with pairwise.complete.obs
        double zmecs=qchisq(pmecs,1);
        vector<double> beta1,beta2;
        vector<int> pairsNoCor1,pairsNoCor2;
        double sumcor=0.0;
        int pairHasCorNUm=0;
        for( int i=0;i<cohortnum;i++)
            for(int j=i+1;j<cohortnum;j++)
            {

                beta1.clear();
                beta2.clear();
                for(int k=0;k<snpnum;k++)
                {
                    double sei=buffer_se(k,i);
                    double betai=buffer_beta(k,i);
                    double sej=buffer_se(k,j);
                    double betaj=buffer_beta(k,j);
                    if(abs(sei+9)>1e-6 && abs(sej+9)>1e-6) {
                        double zi=betai/sei;
                        double zj=betaj/sej;
                        zi*=zi;
                        zj*=zj;
                        if(zi < zmecs && zj < zmecs)
                        {

                            beta1.push_back(betai);
                            beta2.push_back(betaj);
                        }
                    }
                }
                if(beta1.size()<nmecs) {
                    //LOGPRINTF("WARNING: %ld SNP in common between cohort %i and cohort %d (cohort number stats form 0).\n",beta1.size(),i,j);
                    //LOGPRINTF("The correlation value of cohort %d and cohort %d would be imputed with the mean of all the correlation values excluding the diagnoal.\n",i,j);
                    pairsNoCor1.push_back(i);
                    pairsNoCor2.push_back(j);
                    PCC(i,j)=PCC(j,i)=0;
                } else {
                    double corrtmp=cor(beta1,beta2);
                    sumcor +=corrtmp;
                    PCC(i,j)=PCC(j,i)=corrtmp;
                    pairHasCorNUm++;
                }
            }
        if(pairsNoCor1.size()>0) {
            //LOGPRINTF("WARNING: %ld cohort pairs didn't get enough common SNPs to calcualte the correlation.\n",pairsNoCor1.size());
            if(pairHasCorNUm==0) {
                //LOGPRINTF("ERROR: Every pair of cohort has not enough common SNPs to calcualte the correlation.\n");
                //TERMINATE();
                return false;
            }
            double corMean=sumcor/pairHasCorNUm;
            //LOGPRINTF("WARNING: These missing correlation values are imputed with the mean %f.\n",corMean);
            for(int i=0;i<pairsNoCor1.size();i++)
            {
                int p1=pairsNoCor1[i];
                int p2=pairsNoCor2[i];
                PCC(p1,p2)=PCC(p2,p1)=corMean;
            }
        }
        for( int i=0;i<cohortnum;i++) PCC(i,i)=1;
        return true;
    }


    void ssQTL(char * outFileName, char * beqtlFileName, char * problstName, \
        char * problst2exclde, char * genelistName, int chr, int prbchr, \
        char * prbname, char * fromprbname, char * toprbname, int prbWind, \
        int fromprbkb, int toprbkb, bool prbwindFlag, char * genename, \
        char * probe2exclde, int autosome_num, double maf, char * snplstName, \
        char * snplst2exclde, int snpchr, char * snprs, char * fromsnprs, \
        char * tosnprs, int snpWind, int fromsnpkb, int tosnpkb, bool snpwindFlag, \
        char * snprs2exclde, int tsk_ttl, int tsk_id, bool tosmrflag, \
        bool nofastlinear, bool cis_flag, int cis_itvl,  char * annofileName, \
        double pmecs, int nmecs, bool use_top_p, bool trans_flag, int trans_itvl, bool no_isoform_eqtl, bool not_use_top)
    {

        setNbThreads(thread_num);
        LOGPRINTF("Using %d thread(s) to conduct analysis ...\n", thread_num);
        if (cis_flag && trans_flag) {
            fprintf(stderr, "The cis_flag and trans_flag can not be true at same time.");
            TERMINATE();
        }
        if (!cis_flag && !trans_flag) {
            cis_flag = true;
        }
        if (!use_top_p && !not_use_top) {
            use_top_p = true;
        } else if (use_top_p && not_use_top) {
            not_use_top = false;
        }
        eqtlInfo eqtlinfo;
        LOGPRINTF("\nReading eQTL summary data...\n");
        char inputname[FNAMESIZE];
        memcpy(inputname,beqtlFileName,strlen(beqtlFileName)+1);
        char* suffix=inputname+strlen(beqtlFileName);
        memcpy(suffix,".epi",5);
        read_smr_epifile(&eqtlinfo, inputname);
        smr_epi_man(&eqtlinfo, problstName, problst2exclde, genelistName,  chr,
            prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,
            toprbkb, prbwindFlag,  genename, probe2exclde);
        extract_sqtl_probe(&eqtlinfo,tsk_ttl,  tsk_id);
        memcpy(suffix,".esi",5);
        read_smr_esifile(&eqtlinfo, inputname);
        smr_esi_man(&eqtlinfo, snplstName, snplst2exclde,chr, snpchr,  snprs,
            fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag,
            cis_flag,  cis_itvl, prbname,snprs2exclde); //cis will work when prbname and cis_flag is ture.
        if(eqtlinfo._include.size()==0)
        {
            LOGPRINTF("Error: no probe included.\n");
            TERMINATE();
        }
        if(eqtlinfo._esi_include.size()==0)
        {
            LOGPRINTF("Error: no SNP included.\n");
            TERMINATE();
        }
        memcpy(suffix,".besd",6);
        vector<int> headers;
        get_BesdHeaders(inputname, headers);
        int indicator = headers[0];
        if(indicator==SMR_DENSE_1 || indicator==SMR_DENSE_3 || indicator==OSCA_DENSE_1)
        {
            LOGPRINTF("Error: Summary data based sQTL analysis does not support dense BESD format.\n");
            TERMINATE();
        }
        eqtlinfo._sampleNum=headers[1];
        if(eqtlinfo._sampleNum==-9)
        {
            LOGPRINTF("Error: no sample size found. Please upadte the BESD file using --add-n.\n");
            TERMINATE();
        }
        char outputname[FNAMESIZE];
        outputname[0]='\0';
        if(tsk_ttl>1) {
            if(outFileName!=NULL) {
                string tmp=  string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outFileName=outputname;
            }
        }

        LOGPRINTF("Loading the file %s into memory...\n",beqtlFileName);
        read_smr_besdfile(&eqtlinfo, inputname);
        if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
        {
            LOGPRINTF("No data included from %s under current condition.\n",beqtlFileName);
            TERMINATE();
        }
        eInfo sqtlinfo;
        vector< vector<int>> tranids;  //smaller einfo._epi_include
        /*different from sQTL gene_check. the chromosome map is generated by eqtlinfo, and probe's chr of eqtinfo must in this chr map, even is call back from annotation data again.
        */
        gene_check( &sqtlinfo,tranids, annofileName, &eqtlinfo);


        LOGPRINTF("\nPerforming sQTL analysis ...\n");
        if(outFileName!=NULL){
            write_smr_esi(outFileName, &eqtlinfo);
            write_smr_epi(outFileName, &sqtlinfo,true);
        }
        FILE* besd=NULL;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=OSCA_SPARSE_1;
        if(tosmrflag) filetype=SMR_SPARSE_3;
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=eqtlinfo._sampleNum;
        ten_ints[2]=(int)eqtlinfo._esi_include.size();
        ten_ints[3]=(int)sqtlinfo._epi_include.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
            TERMINATE();
        }

        vector< vector<uint32_t> > rowids(sqtlinfo._epi_include.size());
        vector< vector<float> > betas(sqtlinfo._epi_include.size());
        vector< vector<float> > ses(sqtlinfo._epi_include.size());
        vector< vector<int> > snpids;
        vector<int> snp_num;

        if (cis_flag){
            cis_eQTL_num_2(&sqtlinfo,&eqtlinfo,cis_itvl,tranids, snpids, snp_num);
            for(int ii=0;ii<sqtlinfo._epi_include.size();ii++)
            {
                rowids[ii].resize(snp_num[ii]);
                betas[ii].resize(snp_num[ii]);
                ses[ii].resize(snp_num[ii]);
            }
        } else if (trans_flag) {
            trans_eQTL_num(&sqtlinfo, &eqtlinfo, trans_itvl, tranids, snpids, snp_num);
            for (int ii = 0; ii < sqtlinfo._epi_include.size(); ii++) {
                rowids[ii].resize(snp_num[ii]);
                betas[ii].resize(snp_num[ii]);
                ses[ii].resize(snp_num[ii]);
            }

        }
        bool warned = false;
        int nindi = eqtlinfo._sampleNum;
        double cr=0.0;
        unsigned int i = 0, j = 0, k = 0, l = 0;
        vector < int > need_remove;
        vector < int >::iterator it;
        FILE * fmid_data = NULL;
        if (!no_isoform_eqtl) {
            char snp_data_out_file[1024] = "";
            if (strlen(outFileName) < 1000) {
                strcpy(snp_data_out_file, outFileName);
                strcat(snp_data_out_file, "_isoform_eQTL_effect.txt");
            } else {
                fprintf(stderr, "outFileName too long, using snp_data.txt instead.\n");
                strcpy(snp_data_out_file, "isoform_eQTL_effect.txt");
            }
            fmid_data = fopen(snp_data_out_file, "w");
            if (!fmid_data) {
                fprintf(stderr, "open snp_data.txt failed.\n");
                exit(1);
            }
            fprintf(fmid_data, "#SNP\tchr_snp\tBP_snp\tA1\tA2\tFreq\t"
                "Probe\tProbe_chr\tgene\torientation\tbeta\tse\tp_val\n"
            );
        } 

        #pragma omp parallel for private(cr)
        for (int jj = 0; jj < sqtlinfo._epi_include.size(); jj++) {
            string gene_name = sqtlinfo._epi_gene[sqtlinfo._epi_include[jj]];
            string prb_name = sqtlinfo._epi_prb[sqtlinfo._epi_include[jj]];
            int prb_chr = sqtlinfo._epi_chr[sqtlinfo._epi_include[jj]];
            char prb_ori = sqtlinfo._epi_orien[sqtlinfo._epi_include[jj]];

            LOGPRINTF(">%s \n", gene_name.c_str());

            double desti = 1.0 * jj / (sqtlinfo._epi_include.size() - 1);
            if (desti >= cr) {
                printf("%3.0f%%\r", 100.0 * desti);
                fflush(stdout);
                if(cr == 0)
                    cr+=0.05;
                else if(cr == 0.05)
                    cr+=0.2;
                else if(cr == 0.25)
                    cr+=0.5;
                else
                    cr+=0.25;
            }
            string prbid = sqtlinfo._epi_prb[sqtlinfo._epi_include[jj]];
            if(snpids[jj].size()==0){
                continue;
            }

            int numTrans = (int)tranids[jj].size();
            vector <int> iso_idx = tranids[jj];
            int numSNP = snp_num[jj];
            if (numSNP != snpids[jj].size()) {
                fprintf(stderr, "snp_num not eqtal to snpids.size\n");
                TERMINATE();
            }

            MatrixXd eqtlb(numSNP, numTrans), eqtls(numSNP, numTrans);
            for (int i = 0; i < numSNP; i++) {
                for (int j = 0; j < numTrans; j++) {
                    eqtlb(i, j) = -9;
                    eqtls(i, j) = -9;
                }
            }
            vector < float > eqtlfreq(numSNP);
            make_bs(&eqtlinfo, tranids[jj], snpids[jj], eqtlb, eqtls, eqtlfreq);

            LOGPRINTF("numTrans %d\n", numTrans);

            // estimate probe correlaton
            MatrixXd cor_null(numTrans, numTrans);
            pcc(cor_null,  eqtlb, eqtls, pmecs, nmecs);


            //filter cor_null, eqtls eqtlb
            need_remove.clear();
            for (i = 0; i < cor_null.rows(); i++) {
                for (j = i + 1; j < cor_null.cols(); j++){
                    if (abs(cor_null(i, j) - 1) < 1e-15) {
                        need_remove.push_back(i);
                        break;
                    }
                }
            }
            vector <int> iso_idx_tmp;
            if (need_remove.size() > 0) {
                for (int i = 0; i < iso_idx.size(); i++)
                {
                    it = find(need_remove.begin(), need_remove.end(), i);
                    if (it == need_remove.end())
                    {
                        iso_idx_tmp.push_back(iso_idx[i]);
                    }
                }
                iso_idx = iso_idx_tmp;
            }
            

            numTrans -= need_remove.size();
            LOGPRINTF("numTrans after filter: %u\n", numTrans);
            MatrixXd cor_null_clean(numTrans, numTrans);
            MatrixXd eqtlb_clean(numSNP, numTrans);
            MatrixXd eqtls_clean(numSNP, numTrans);
            k = 0;
            if (need_remove.size() > 0 && numTrans > 1) {
                l = 0;
                for (i = 0; i < cor_null.rows(); i++){
                    it = find(need_remove.begin(), need_remove.end(), i);
                    if (it != need_remove.end()){
                        k++;
                    } else {
                        for (j = 0; j < cor_null.cols(); j++){
                            it = find(need_remove.begin(), need_remove.end(), j);
                            if (it != need_remove.end()){
                                l++;
                            } else{
                                cor_null_clean(i - k, j - l) = cor_null(i, j);
                            }
                        }
                    }
                }

                for (i = 0; i < eqtlb.rows(); i++) {
                    k = 0;
                    for (j = 0; j < eqtlb.cols(); j++) {
                        it = find(need_remove.begin(), need_remove.end(), j);
                        if (it != need_remove.end()) {
                            k++;
                        } else {
                            eqtlb_clean(i, j - k) = eqtlb(i, j);
                        }
                    }
                }

                for (i = 0; i < eqtls.rows(); i++) {
                    k = 0;
                    for (j = 0; j < eqtls.cols(); j++) {
                        it = find(need_remove.begin(), need_remove.end(), j);
                        if (it != need_remove.end()) {
                            k++;
                        } else {
                            eqtls(i, j - k) = eqtls(i, j);
                        }
                    }
                }

                cor_null = cor_null_clean;
                eqtlb = eqtlb_clean;
                eqtls = eqtls_clean;

            } else {
                numTrans += need_remove.size();
            }


            if(loud)
            {
                string filename=string(outfileName)+"."+atos(prbid)+".corr.txt";
                FILE* tmpfile =fopen(filename.c_str(),"w");
                if(!tmpfile)
                {
                    printf("ERROR: open file %s.\n",filename.c_str());
                    exit(EXIT_FAILURE);
                }
                for(int kk=0;kk<cor_null.rows();kk++)
                {
                    string str="";
                    for(int jj=0;jj<cor_null.cols();jj++)
                    {
                        str+= atos(cor_null(kk,jj))+'\t';
                    }
                    str+='\n';
                    fputs(str.c_str(),tmpfile);
                }
                fclose(tmpfile);
            }

            LOGPRINTF("numSNP: %u\n",  numSNP);
            for(int kk = 0; kk < numSNP; kk++)
            {
                uint32_t snpid = snpids[jj][kk];
                string snprs = eqtlinfo._esi_rs[snpid];
                double snpfreq = eqtlinfo._esi_freq[snpid];
                int snp_chr = eqtlinfo._esi_chr[snpid];
                uint64_t snp_bp = eqtlinfo._esi_bp[snpid];
                string a1 = eqtlinfo._esi_allele1[snpid];
                string a2 = eqtlinfo._esi_allele2[snpid];

                if(snpfreq == 0 || snpfreq == 1) {
                    if (!warned) {
                        LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s).\n");
                        warned=1;
                    }
                    rowids[jj][kk] = snpid;
                    betas[jj][kk] = 0;
                    ses[jj][kk] = 1;
                    continue;
                }

                vector<double> beta, se;
                vector <int> beta_se_rm;
                double beta_tmp, se_tmp;
                for(int ll = 0; ll < numTrans; ll++)
                {
                    beta_tmp = eqtlb(kk, ll);
                    se_tmp = eqtls(kk, ll);
                    if (abs(se_tmp + 9) < 1e-15){
                        beta_se_rm.push_back(ll);
                    } else {
                        beta.push_back(beta_tmp);
                        se.push_back(se_tmp);
                    }
                }
                if (beta_se_rm.size() > 0){
                    iso_idx_tmp.clear();
                    for (int i = 0; i < iso_idx.size(); i++) {
                        it = find(beta_se_rm.begin(), beta_se_rm.end(), i);
                        if (it == beta_se_rm.end()) {
                            iso_idx_tmp.push_back(iso_idx[i]);
                        }
                    }
                    iso_idx = iso_idx_tmp;
                }
                int numTrans_snp = 0;
                numTrans_snp = numTrans - beta_se_rm.size();
                MatrixXd cor_null_snp(numTrans_snp, numTrans_snp);
                if (beta_se_rm.size() > 0) {
                    k = 0;
                    for (int i = 0; i < cor_null.rows(); i++) {
                        l = 0;
                        it = find(beta_se_rm.begin(), beta_se_rm.end(), i);
                        if (it != beta_se_rm.end()) {
                            k++;
                        } else {
                            for (int j = 0; j < cor_null.cols(); j++) {
                                it = find(beta_se_rm.begin(), beta_se_rm.end(), j);
                                if (it != beta_se_rm.end()) {
                                    l++;
                                } else {
                                    cor_null_snp(i - k, j - l) = cor_null(i, j);
                                }
                            }
                        }
                    }
                } else {
                    cor_null_snp = cor_null;
                }
                int varnum = numTrans_snp * (numTrans_snp - 1) / 2;
                VectorXd d(varnum), vardev(varnum), chisq_dev(varnum);
                int i = 0, j = 0, k = 0;
                MatrixXd vdev (varnum, varnum);
                MatrixXd corr_dev (varnum, varnum);
                double beta_mean = 0;
                double var_mean = 0;

                if (use_top_p) {
                    d.resize(numTrans_snp);
                    vardev.resize(numTrans_snp);
                    chisq_dev.resize(numTrans_snp);
                    vdev.resize(numTrans_snp, numTrans_snp);
                    corr_dev.resize(numTrans_snp, numTrans_snp);
                    beta_mean = 0;
                    for (i = 0; i < numTrans_snp; i++) {
                        beta_mean += beta[i];
                    }
                    beta_mean = beta_mean / numTrans_snp;
                    var_mean = 0;
                    var_mean = get_var_mean(se, cor_null_snp);
                    for (i = 0; i < numTrans_snp; i++) {
                        d[i] = beta[i] - beta_mean;
                        vardev[i] = se[i] * se[i] + var_mean - \
                            2 * get_cov_beta_mean(i, se, cor_null_snp);
                    }

                    for (i = 0; i < numTrans_snp; i++) {
                        chisq_dev[i] = d[i] * d[i] / vardev[i];
                    }

                    for (i = 0; i < numTrans_snp; i++) {
                        for (j = 0; j < numTrans_snp; j++) {
                            vdev(i, j) = se[i] * se[j] * cor_null_snp(i, j) - \
                                get_cov_beta_mean(i, se, cor_null_snp) - \
                                get_cov_beta_mean(j, se, cor_null_snp) + \
                                var_mean;
                        }
                    }

                    for (i = 0; i < numTrans_snp; i++) {
                        for (j = i; j < numTrans_snp; j++) {
                            corr_dev(i, j) = corr_dev(j, i) = \
                                vdev(i, j) / sqrt(vdev(i, i) * vdev(j, j));
                        }
                    }


                } else {
                    k = 0;
                    for(int m1 = 0; m1 < numTrans_snp - 1; m1++) {
                        for(int m2 = m1 + 1; m2 < numTrans_snp; m2++){
                            d[k]=beta[m1]-beta[m2];
                            vardev[k] = se[m1] * se[m1] + se[m2] * se[m2] - 2 * \
                                cor_null_snp(m1, m2) * se[m1] * se[m2];
                            k++;
                        }
                    }
                    for(int m1=0;m1<varnum;m1++)
                        chisq_dev[m1] = d[m1]*d[m1]/vardev[m1];

                    int mi = 0, mj =0;
                    for( int m1=0;m1< numTrans_snp-1;m1++) {
                        for( int m2=m1+1;m2<numTrans_snp;m2++) {
                            mj=0;
                            for(int m3=0; m3< numTrans_snp-1; m3++){
                                for(int m4=m3+1;m4<numTrans_snp;m4++) {
                                    vdev(mi,mj) = se[m1] * se[m3] * cor_null_snp(m1,m3) - \
                                        se[m1] * se[m4] * cor_null_snp(m1,m4) - \
                                        se[m2] * se[m3] * cor_null_snp(m2,m3) + \
                                        se[m2] * se[m4] * cor_null_snp(m2,m4);
                                    mj++;
                                }
                            }
                            mi++;
                        }
                    }

                    for( int m1=0;m1<varnum;m1++) {
                        for( int m2=m1;m2<varnum;m2++){
                            corr_dev(m1,m2) = corr_dev(m2,m1) = vdev(m1,m2)/sqrt(vdev(m1,m1)*vdev(m2,m2));
                        }
                    }
                }

                VectorXd lambda;
#pragma omp critical
                {

                    SelfAdjointEigenSolver<MatrixXd> es(corr_dev, EigenvaluesOnly);
                    lambda=es.eigenvalues();
                
                    double sumChisq_dev=chisq_dev.sum();
                    double pdev= 0.0;
                    pdev=pchisqsum(sumChisq_dev,lambda);
                    double z=0.0;

                    z=sqrt(qchisq(pdev,1));
                    double beta_hat=z/sqrt(2*snpfreq*(1-snpfreq)*(nindi+z*z));
                    double se_hat=1/sqrt(2*snpfreq*(1-snpfreq)*(nindi+z*z));

                    rowids[jj][kk]=snpid;
                    betas[jj][kk]=beta_hat;
                    ses[jj][kk]=se_hat;

                    if (fmid_data) {
                        fprintf(fmid_data, 
                            "%s\t%d\t%lu\t%s\t%s\t%le\t"       //snp, chr_snp, bp_snp, a1, a2, Freq
                            "%s\t%d\t%s\t%c\t%le\t%le\t%le",   //probe, probe_chr, gene, orien, beta, se, pval
                            snprs.c_str(), snp_chr, snp_bp, a1.c_str(), a2.c_str(), snpfreq,
                            prb_name.c_str(), prb_chr, gene_name.c_str(), prb_ori, beta_hat, se_hat, pdev
                        );
                        for (int i = 0; i < iso_idx.size(); i++)
                        {
                            double snp_t = beta[i] / se[i];
                            snp_t = snp_t * snp_t;
                            double snp_p = pchisq(snp_t, 1);
                            fprintf(fmid_data,
                                "\t%s\t%le\t%le\t%le",
                                eqtlinfo._epi_prbID[iso_idx[i]].c_str(),
                                beta[i], se[i], snp_p
                            );
                        }
                        fprintf(fmid_data, "\n");

                        if (numTrans_snp != iso_idx.size())
                        {
                            fprintf(stderr, "Waring, isoform num should equle to numTrans_snp, %s\t%s\t%s\t%u\t%u\n",
                                    snprs.c_str(), prbid.c_str(),
                                    prb_name.c_str(), numTrans, numTrans_snp, iso_idx.size());
                        }
                    }
                    
                }
            
            }
        }
        if (fmid_data) {
            fclose(fmid_data);
        }

        if (tosmrflag) {
            uint64_t valNum=0;
            vector<uint64_t> cols;
            cols.resize(2*sqtlinfo._epi_include.size()+1);
            cols[0]=0;
            for(int jj=0;jj<betas.size();jj++)
            {
                long co=betas[jj].size();
                valNum+=co;
                cols[2*jj+1]=cols[2*jj]+co;
                valNum+=co;
                cols[2*(jj+1)]=cols[2*jj+1]+co;
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
                {
                    if (fwrite_checked(&rowids[jj][0],rowids[jj].size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&rowids[jj][0],rowids[jj].size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
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
        else {
            uint64_t valNum=0;
            vector<uint64_t> cols;
            cols.resize(sqtlinfo._epi_include.size()+1);
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


        fclose(besd);
        //LOGPRINTF("sQTL summary statistics for %ld probes and %ld SNPs are saved in file %s.\n", einfo._epi_include.size(),bdata._include.size(), besdName.c_str());
    }
}

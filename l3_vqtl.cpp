//
//  l3_vqtl.cpp
//  osc
//
//  Created by Futao Zhang on 29/08/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#include "l3_vqtl.hpp"
using namespace BFILE;
using namespace EFILE;
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

    void load_vqtl_workspace(eInfo* einfo,bInfo* bdata, char* efileName, char* befileName,  char* bFileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id)
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
            read_famfile(bdata, string(bFileName)+".fam");
            indi_check(bdata,einfo);
            memcpy(suffix,".opi",5);
            read_epi(inputname,einfo);
            epi_man(einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            if(tsk_ttl>1) extract_probe(einfo,  tsk_ttl,  tsk_id);
            memcpy(suffix,".bod",5);
            clock_t begin_time = clock();
            read_beed(inputname,einfo);
            LOGPRINTF("read_beed: %f ms.\n",float( clock () - begin_time ) /  1000);
            if(dpvalfName!=NULL) filtering_with_detpval(einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(einfo, upperBeta, lowerBeta);
            if(mssratio_prob<1) filtering_with_missingratio(einfo, mssratio_prob);
            read_bimfile(bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(bdata, snplst2exclde);
             begin_time = clock();
            read_bedfile(bdata, string(bFileName)+".bed");
            LOGPRINTF("read_bedfile: %f ms.\n",float( clock () - begin_time ) /  1000);
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
        
    }

    
     void V_QTL(char* outFileName,  char* efileName, char* befileName,  char* bFileName, bool transposed,  int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob,int autosome_num, double maf,char* snplstName,char* snplst2exclde,int tsk_ttl,int tsk_id,int vqtl_mtd)
    {
        
        setNbThreads(thread_num);
        eInfo einfo;
        bInfo bdata;
        init_einfo(&einfo);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, bFileName,transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num,snplstName,snplst2exclde,tsk_ttl,tsk_id); // using _keep and _eii_include, the individuals are aligned.
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
        FILE* assoc=NULL;
        long write_count=0;
        string outstr="";
        string assocfile="";
        if(outFileName!=NULL){
            assocfile = string(outFileName)+".vqtl";
            assoc = fopen(assocfile.c_str(), "w");
            if (!(assoc)) {
                LOGPRINTF("ERROR: open error %s\n", assocfile.c_str());
                TERMINATE();
            }
            outstr="ProbeChr\tProbeID\tProbe_bp\tSNPChr\tSNP\tSNP_bp";
            if(vqtl_mtd==0 || vqtl_mtd==3 ) outstr+="\tStatistic\tdf\tP\tNMISS\n";
            else outstr+="\tStatistic\tdf1\tdf2\tP\tNMISS\n";
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
        } 
        
        if(vqtl_mtd==0) { LOGPRINTF("\nPerforming vQTL analysis with Bartlett’s test...\n");}
        else if (vqtl_mtd==1) {LOGPRINTF("\nPerforming vQTL analysis with Levene’s test (mean)...\n");}
        else if (vqtl_mtd==2) {LOGPRINTF("\nPerforming vQTL analysis with Levene’s test (median)...\n");}
        else  {LOGPRINTF("\nPerforming vQTL analysis with Fligner-Killeen test...\n");}
        
        int slide_wind=10000;
        int loops=ceil(1.0*bdata._include.size()/slide_wind);
        int warn1=0, warn2=0, warn3=0;
        clock_t begin_time = clock();
        for(int ii=0;ii<loops;ii++)
        {
            MatrixXf _X;
            int snpstart=ii*slide_wind;
            make_XMat(&bdata,snpstart,slide_wind, _X);
            
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                printf("%3.0f%%\r", 100.0*(ii*einfo._epi_include.size()+jj)/(loops*einfo._epi_include.size()));
                fflush(stdout);
                
                int chr=einfo._epi_chr[einfo._epi_include[jj]];
                string prbid=einfo._epi_prb[einfo._epi_include[jj]];
                string gene=einfo._epi_gene[einfo._epi_include[jj]];
                int BP=einfo._epi_bp[einfo._epi_include[jj]];
                char oren=einfo._epi_orien[einfo._epi_include[jj]];
                for(int kk=0;kk<_X.cols();kk++)
                {
                    int snpchr=bdata._chr[bdata._include[snpstart+kk]];
                    string snprs=bdata._snp_name[bdata._include[snpstart+kk]];
                    int snpbp=bdata._bp[bdata._include[snpstart+kk]];
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
                        
                        int flag=bartlett(yvec,bvec, rst);
                        
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
                            if(assoc) {
                                string outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                if(fputs_checked(outstr.c_str(),assoc))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                                    TERMINATE();
                                }
                                write_count++;
                            }
                        }
                    } else if(vqtl_mtd==1)
                    {
                        int flag=leveneTest_mean(yvec,bvec, rst);
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
                            if(assoc) {
                                string outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + atos(rst[2]) + '\t' + dtos(rst[3]) + '\t' + atos(nonmiss) +'\n';
                                if(fputs_checked(outstr.c_str(),assoc))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                                    TERMINATE();
                                }
                                write_count++;
                            }
                        }
                    } else if(vqtl_mtd==2)
                    {
                        int flag=leveneTest_median(yvec,bvec, rst);
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
                            if(assoc) {
                                string outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + atos(rst[2]) + '\t' + dtos(rst[3]) + '\t' + atos(nonmiss) +'\n';
                                if(fputs_checked(outstr.c_str(),assoc))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                                    TERMINATE();
                                }
                                write_count++;
                            }
                        }

                    } else if(vqtl_mtd==3)
                    {
                        int flag=flignerTest(yvec,bvec, rst);
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
                            if(assoc) {
                                string outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(snpchr) + '\t' + snprs + '\t' + atos(snpbp) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                                if(fputs_checked(outstr.c_str(),assoc))
                                {
                                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                                    TERMINATE();
                                }
                                write_count++;
                            }
                        }

                    }
                   
                }
            }

        }
        LOGPRINTF("%ld tests cost: %f ms.\n",write_count,float( clock () - begin_time ) /  1000);
        if(assoc){
            LOGPRINTF("Results of %ld vQTLs have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        }
    }
}

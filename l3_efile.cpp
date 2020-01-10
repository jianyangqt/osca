//
//  l3_efile.cpp
//  osc
//
//  Created by Futao Zhang on 11/04/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l3_efile.h"
namespace EFILE {
    int comp_eii(const void *a,const void *b){ return ( ( (*(indiinfolst *)a).fid > (*(indiinfolst *)b).fid ) || ( ((*(indiinfolst *)a).fid==(*(indiinfolst *)b).fid) && ( (*(indiinfolst *)a).iid>(*(indiinfolst *)b).iid ) ) )?1:-1; }
    int comp_epi(const void *a,const void *b){ return ((*(probeinfolst *)a).probeId>(*(probeinfolst *)b).probeId)?1:-1; }
    int comp_epi2(const void *a,const void *b){ return (((*(probeinfolst *)a).probechr>(*(probeinfolst *)b).probechr) || ( ((*(probeinfolst *)a).probechr ==(*(probeinfolst *)b).probechr) && ((*(probeinfolst *)a).bp > (*(probeinfolst *)b).bp) ))?1:-1; }
    
    void read_beflist(vector<string> &befNames, char* befileFlistName)
    {
        FILE* fptr=NULL;
        vector<string> strlist;
        map<string, int> befile_map;
        long mapsize=0;
        if(fopen_checked(&fptr, befileFlistName,"r")) TERMINATE();
        LOGPRINTF("Reading file names from %s ...\n", befileFlistName);
        while(fgets(Tbuf, MAX_LINE_SIZE, fptr))
        {
            split_str(Tbuf,strlist,0);
            if(strlist.size()==0)
            {
                LOGPRINTF("Blank row found and skipped!\n");
                continue;
            } else if( strlist.size() >1) {
                LOGPRINTF("Column number is not correct with this row\n %s\n",Tbuf);
                TERMINATE();
            } else {
                befile_map.insert(pair<string,int>(strlist[0],mapsize));
                if(mapsize<befile_map.size())
                {
                    befNames.push_back(strlist[0]);
                    mapsize=befile_map.size();
                }
                else
                {
                    LOGPRINTF("WARNING: duplicate summary file name %s found and skipped.\n",strlist[0].c_str());
                }
                
            }
        }
        LOGPRINTF("%ld file names are included from %s.\n", befNames.size(), befileFlistName);
        fclose(fptr);
    }
    void combine_eii_outer(vector<indiinfolst> &indiinfo, vector<string> &befNames, vector<uint64_t> &nindi)
    {
        long counter = 0;
        indiinfo.clear();
        map<string, int> in_map;
        map<string, int>::iterator iter;
        nindi.clear();
        char inputname[FNAMESIZE];
        for (int i = 0; i < befNames.size(); i++)
        {
            eInfo etmp;
            memcpy(inputname,befNames[i].c_str(),befNames[i].length()+1);
            char* suffix=inputname+befNames[i].length();
            memcpy(suffix,".oii",5);
            read_eii(inputname,&etmp);
            nindi.push_back(etmp._eii_num);
            for (int j = 0; j<etmp._eii_num; j++)
            {
                string fiidstr=etmp._eii_fid[j]+":"+atos(etmp._eii_iid[j]);
                iter=in_map.find(fiidstr.c_str());
                if(iter==in_map.end())
                {
                    in_map.insert(pair<string, int>(fiidstr.c_str(), indiinfo.size()));//the second is snpinfo id
                    
                    indiinfolst iinfotmp;
                    counter = in_map.size();
                    strcpy2(&iinfotmp.fid, etmp._eii_fid[j]);
                    strcpy2(&iinfotmp.iid, etmp._eii_iid[j]);
                    strcpy2(&iinfotmp.fa_id, etmp._eii_fa_id[j]);
                    strcpy2(&iinfotmp.mo_id, etmp._eii_mo_id[j]);
                    iinfotmp.sex=etmp._eii_sex[j];
                    iinfotmp.pheno=etmp._eii_pheno[j];
                    iinfotmp.itr=new int[befNames.size()];
                    for(int k=0;k<befNames.size();k++){
                        if(i==k){
                            iinfotmp.itr[k]=j;
                        } else {
                            iinfotmp.itr[k]=-9;
                        }
                    }
                    indiinfo.push_back(iinfotmp);
                } else {
                        indiinfo[iter->second].itr[i]=j;
                }
            }
        }
        if(in_map.size()  != indiinfo.size()){
            LOGPRINTF("ERROR: bugs found. please report.\n") ;
            TERMINATE();
        }
        LOGPRINTF("Total %ld individuals to be included from %ld eii files.\n",indiinfo.size(),befNames.size());
        
    }
    void combine_epi(vector<probeinfolst> &probeinfo, vector<string> &befNames, vector<uint64_t> &nprb)
    {
        long counter = 0;
        map<string, int> prb_map;
        map<string, int> prbbp_map;
        map<string, int>::iterator iter;
        nprb.clear();
        char inputname[FNAMESIZE];
        for (int i = 0; i < befNames.size(); i++)
        {
            eInfo etmp;
            memcpy(inputname,befNames[i].c_str(),befNames[i].length()+1);
            char* suffix=inputname+befNames[i].length();
            memcpy(suffix,".opi",5);
            read_epi(inputname, &etmp);
            nprb.push_back(etmp._epi_num);
            for (int j = 0; j<etmp._epi_num; j++)
            {
                string crsbpstr=etmp._epi_prb[j]+":"+atos(etmp._epi_bp[j]);
                prb_map.insert(pair<string, int>(etmp._epi_prb[j].c_str(), counter));
                prbbp_map.insert(pair<string, int>(crsbpstr.c_str(), counter));
                if(prb_map.size() != prbbp_map.size())
                {
                    LOGPRINTF("ERROR: inconsistent position for the probe %s  in different .epi files. Please check.\n", etmp._epi_prb[j].c_str()) ;
                    TERMINATE();
                }
                
                if (counter < prb_map.size())
                {
                    probeinfolst probinfotmp;
                    counter=prb_map.size();
                    probinfotmp.probechr=etmp._epi_chr[j];
                    strcpy2(&probinfotmp.probeId, etmp._epi_prb[j]);
                    probinfotmp.bp=etmp._epi_bp[j];
                    probinfotmp.gd=etmp._epi_gd[j];
                    strcpy2(&probinfotmp.genename, etmp._epi_gene[j]);
                    probinfotmp.orien=etmp._epi_orien[j];
                    probinfotmp.ptr=new int[befNames.size()];
                    for(int k=0;k<befNames.size();k++){
                        if(i==k){
                            probinfotmp.ptr[k]=j;
                        } else {
                            probinfotmp.ptr[k]=-9;
                        }
                    }
                    probeinfo.push_back(probinfotmp);
                    
                } else {
                    iter=prb_map.find(etmp._epi_prb[j]);
                    if(iter!=prb_map.end())
                    {
                        probeinfo[iter->second].ptr[i]=j; //probeinfo with prb_map
                    }
                    else
                    {
                        LOGPRINTF("ERROR: This would never happen. please help to report this bug.\n") ;
                        TERMINATE();
                    }
                }
            }
        }
        LOGPRINTF("Total %ld probes to be included from %ld opi files.\n",probeinfo.size(),befNames.size());
    }
    void write_epi(char* outFileName, vector<probeinfolst> &probeinfo)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".opi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        for(int i=0;i<probeinfo.size();i++)
        {
            string chrstr;
            if(probeinfo[i].probechr==23) chrstr="X";
            else if(probeinfo[i].probechr==24) chrstr="Y";
            else chrstr=atosm(probeinfo[i].probechr);

            string str=chrstr+'\t'+probeinfo[i].probeId+'\t'+atosm(probeinfo[i].bp)+'\t'+probeinfo[i].genename+'\t'+(probeinfo[i].orien=='*'?"NA":atos(probeinfo[i].orien))+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld probes have been saved in the file %s .\n", probeinfo.size(), epiName.c_str());
    }
    void write_eii(char* outFileName, vector<indiinfolst> &indiinfo)
    {
        FILE* efile=NULL;
        string eiiName=string(outFileName)+".oii";
        if(fopen_checked(&efile, eiiName.c_str(),"w")) TERMINATE();
        for(int i=0;i<indiinfo.size();i++)
        {
            string str=string(indiinfo[i].fid)+'\t'+string(indiinfo[i].iid)+'\t'+string(indiinfo[i].fa_id)+'\t'+string(indiinfo[i].mo_id)+'\t'+atosm(indiinfo[i].sex)+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                
                LOGPRINTF("ERROR: in writing file %s .\n", eiiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld individuals have been saved in the file %s .\n", indiinfo.size(), eiiName.c_str());
        
    }

    void extract_prb(FILE* fptr,  uint64_t pid, uint64_t epinum,uint64_t eiinum, vector<double> &profval)
    {
        if(epinum==0 || eiinum==0) {
            LOGPRINTF("ERROR: .opi file or .oii file is empty. please check.\n");
            TERMINATE();
        }
        if(pid>epinum) {
            LOGPRINTF("ERROR: probe index %llu is larger than the totoal probe number %llu.\n", pid, epinum);
            TERMINATE();
        }
        int infoLen=sizeof(uint32_t)*3;
        profval.resize(eiinum);
        fseek(fptr,((pid*eiinum)<<3)+infoLen, SEEK_SET);
        fread(&profval[0], sizeof(double),eiinum,fptr);
    }

    void merge_beed(char* outfileName, char* befileFlistName, char* problstName, char* problst2exclde,char* genelistName,  int chr, char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename,char* probe2rm,char* indilstName,char* indilst2remove,bool beta2m,bool m2beta)
    {
        vector<string> befNames;
        vector<indiinfolst> indiinfo;
        vector<probeinfolst> probeinfo;
        vector<uint64_t> nprb,nindi;
        vector< vector<int> > lookup;
        read_beflist(befNames, befileFlistName);
        if(befNames.size()==0){
            LOGPRINTF("No file names included from %s ...\n", befileFlistName);
            TERMINATE();
        }
        combine_eii_outer(indiinfo, befNames, nindi);
        if(indiinfo.size()==0)
        {
            LOGPRINTF("ERROR: No individuals to be included!\n");
            TERMINATE();
        }
        indiinfolst* eiiptr=&indiinfo[0];
        qsort(eiiptr,indiinfo.size(),sizeof(indiinfolst),comp_eii);
    
        combine_epi(probeinfo, befNames, nprb);
        if(probeinfo.size()==0)
        {
            LOGPRINTF("ERROR: No probes to be included!\n");
            TERMINATE();
        }
        int cotmp=0;
        for(int i=0;i<probeinfo.size();i++)
            if(probeinfo[i].probechr<0 || probeinfo[i].bp<0) cotmp++;
        if(cotmp*1.0/probeinfo.size() < 0.5) cotmp=0;
        probeinfolst* epiptr=&probeinfo[0];
        if(cotmp) qsort(epiptr,probeinfo.size(),sizeof(probeinfolst),comp_epi);
        else qsort(epiptr,probeinfo.size(),sizeof(probeinfolst),comp_epi2);
        
        LOGPRINTF("\nGenerating opi file...\n");
        write_epi(outfileName, probeinfo);
        LOGPRINTF("\nGenerating oii file...\n");
        write_eii(outfileName, indiinfo);
        LOGPRINTF("\nGenerating bod file...\n");
        
        lookup.resize(befNames.size());
        for(int i=0;i<befNames.size();i++) lookup[i].resize(nindi[i]);
        for(int i=0;i<befNames.size();i++)
            for(int j=0;j<lookup[i].size();j++)
                lookup[i][j]=-9;
        for(int i=0;i<indiinfo.size();i++)
        {
            for(int j=0;j<befNames.size();j++)
            {
                int tmpval=indiinfo[i].itr[j];
                if(tmpval>=0) {
                    if(tmpval >=nindi[j])
                    {
                        LOGPRINTF("ERROR: bug found in snpinfo. Please report.\n");
                        TERMINATE();
                    }
                    lookup[j][tmpval]=i;
                }
            }
        }

        uint32_t indicator=0;
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * befNames.size());
        for(int i=0;i<befNames.size();i++) {
            string bodFileName=befNames[i]+".bod";
            if(fopen_checked(&fptrs[i],bodFileName.c_str(),"rb")) {
                LOGPRINTF("ERROR: in opening file %s .\n", bodFileName.c_str());
                TERMINATE();
            }
            uint32_t tmpindr = readuint32(fptrs[i]);
            uint32_t tmpi= readuint32(fptrs[i]);
            uint32_t tmpp= readuint32(fptrs[i]);
            if(i==0) indicator=tmpindr;
            else {
                if(indicator!=tmpindr) {
                    LOGPRINTF("ERROR: the file %s has different file type from the other files.\n",befNames[i].c_str());
                    TERMINATE();
                }
            }
            if(tmpi != nindi[i])
            {
                LOGPRINTF("ERROR: in %s : the individual number in .bod file is not consistent with .oii file.\n",befNames[i].c_str());
                TERMINATE();
            }
            if(tmpp != nprb[i])
            {
                LOGPRINTF("ERROR: in %s : the probe number in .bod file is not consistent with .opi file.\n",befNames[i].c_str());
                TERMINATE();
            }
        }
        uint32_t mni = (uint32_t)indiinfo.size();
        uint32_t mnp = (uint32_t)probeinfo.size();
        FILE * bodfptr;
        string bodName=string(outfileName)+".bod";
        if(fopen_checked(&bodfptr, bodName.c_str(),"wb")) { TERMINATE();}
        
        if (fwrite_checked(&indicator, sizeof(uint32_t), bodfptr))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", bodName.c_str());
            TERMINATE();
        }
        if (fwrite_checked(&mni, sizeof(uint32_t), bodfptr))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", bodName.c_str());
            TERMINATE();
        }
        
        if (fwrite_checked(&mnp, sizeof(uint32_t), bodfptr))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", bodName.c_str());
            TERMINATE();
        }
        vector<double> profval, mprofval(indiinfo.size());
        for(int i=0;i<probeinfo.size();i++)
        {
            //printf("%3.0f%%\r", 100.0*i/(probeinfo.size()));
            fflush(stdout);
            LOGPRINTF("processing with probe %s...\n",probeinfo[i].probeId);
            for(int j=0;j<mprofval.size();j++) mprofval[j]=1e10;
            for(int j=0;j<befNames.size();j++)
            {
                int pid=probeinfo[i].ptr[j];
                profval.clear();
                if(pid>=0)
                {
                    extract_prb(fptrs[j], (uint64_t)pid, nprb[j], nindi[j], profval);
                    if(loud) printf("%llu individuals  of probe %s extracted from the file %s.\n",nindi[j],probeinfo[i].probeId,befNames[j].c_str());
                    for(int k=0;k<nindi[j];k++)
                    {
                        int idx=lookup[j][k];
                        if(idx>=0)
                        {
                            if(mprofval[idx]<1e9)
                            {
                                LOGPRINTF("ERROR: Duplicate individual ( %s : %s) of probe %s in different .bod files found.\n", indiinfo[idx].fid, indiinfo[idx].iid,probeinfo[i].probeId);
                                TERMINATE();
                                
                            } else {
                                mprofval[idx]=profval[k];
                            }
                        }
                    }
                } else {
                    if(loud) printf("probe %s is not in the file %s.\n",probeinfo[i].probeId,befNames[j].c_str());
                }
            }
            if(fwrite_checked(&mprofval[0],sizeof(double)*mprofval.size(), bodfptr)){
                LOGPRINTF("ERROR: in writing binary file %s .\n", bodName.c_str());
                TERMINATE();
            }
        }
        
        fclose(bodfptr);
        for(int i=0;i<befNames.size();i++)
        {
            fclose(fptrs[i]);
        }
        free_probelist(probeinfo);
        free_indilist(indiinfo);
        LOGPRINTF("\nThe infomation of %ld probes and %ld individuals has been in binary file %s.\n",probeinfo.size(),indiinfo.size(),bodName.c_str());

    }
    void make_beed(char* outFileName, char* efileName, char* befileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, double missing_ratio_indi, bool adjprb, char* covfileName,char* qcovfileName, bool enveff, char* effprblstfname,char* efffname, bool stdprb, bool rint)
    {
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL && efileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile.\n");
            TERMINATE();
        }
        if(adjprb && covfileName==NULL && qcovfileName==NULL)
        {
            LOGPRINTF("Error: due to --adj-probe is specified. please input the discrete / Continuous covariates by the option --covar or --qcovar.\n");
            TERMINATE();
        }
        if(efileName!=NULL)
        {
            if(transposed) read_efile_t(efileName,&einfo,efileType,no_fid_flag,valueType);
            else read_efile(efileName,&einfo,efileType,no_fid_flag,valueType);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            eii_man(&einfo,indilstName,indilst2remove);
            if(adjprb && covfileName != NULL) read_cov(&einfo, covfileName, false);
            if(adjprb && qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
            if(dpvalfName!=NULL) filtering_with_detpval(&einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
            if(mssratio_prob<1) filtering_with_missingratio(&einfo, mssratio_prob);
            if(missing_ratio_indi<1) filtering_indi_missingratio(&einfo, missing_ratio_indi);
        }else{
            char inputname[FNAMESIZE];
            memcpy(inputname,befileName,strlen(befileName)+1);
            char* suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,&einfo);
            eii_man(&einfo,indilstName,indilst2remove);
            memcpy(suffix,".opi",5);
            read_epi(inputname,&einfo);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            if(adjprb && covfileName != NULL) read_cov(&einfo, covfileName, false);
            if(adjprb && qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
            memcpy(suffix,".bod",5);
            read_beed(inputname,&einfo);
            if(dpvalfName!=NULL) filtering_with_detpval(&einfo, dpvalfName, dp_thresh, prb_thresh, spl_thresh,filter_mth, no_fid_flag);
            if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
            if(mssratio_prob<1) filtering_with_missingratio(&einfo, mssratio_prob);
            if(missing_ratio_indi<1) filtering_indi_missingratio(&einfo, missing_ratio_indi);
        }
        if(einfo._eType == METHYLATION)
        {
            if(beta2m && m2beta){
                //no chance to enter here
                LOGPRINTF("Error: --m2beta should not be with --beta2m.\n");
                TERMINATE();
            }
            if(beta2m && einfo._valType==BETAVALUE) beta_2_m(&einfo);
            if(m2beta && einfo._valType==MVALUE) m_2_beta(&einfo);
        }
        if(stdprb) stdprobe(&einfo);
        if(rint) rintprobe(&einfo);
        if(adjprb) adjprobe(&einfo);
        if(enveff) addEnvEff(&einfo, effprblstfname,efffname);
        string typestr=getFileType(einfo._eType);
        LOGPRINTF("Saving %s data ...\n",typestr.c_str());
        write_eii(outFileName, &einfo);
        write_epi(outFileName, &einfo);
        write_beed(outFileName, &einfo);
    }
    void make_efile(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,bool t_flag, bool impute_mean_flag)
    {
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL && efileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile.\n");
            TERMINATE();
        }
        if(efileName!=NULL)
        {
            if(transposed) read_efile_t(efileName,&einfo,efileType,no_fid_flag,valueType);
            else read_efile(efileName,&einfo,efileType,no_fid_flag,valueType);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            eii_man(&einfo,indilstName,indilst2remove);
            if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        }else{
            char inputname[FNAMESIZE];
            memcpy(inputname,befileName,strlen(befileName)+1);
            char* suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,&einfo);
            eii_man(&einfo,indilstName,indilst2remove);
            memcpy(suffix,".opi",5);
            read_epi(inputname,&einfo);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            memcpy(suffix,".bod",5);
            read_beed(inputname,&einfo);
            if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        }
        
        if(einfo._eType == METHYLATION)
        {
            if(beta2m && m2beta){
                //no chance to enter here
                LOGPRINTF("Error: --m2beta should not be with --beta2m.\n");
                TERMINATE();
            }
            if(beta2m && einfo._valType==BETAVALUE) beta_2_m(&einfo);
            if(m2beta && einfo._valType==MVALUE) m_2_beta(&einfo);
        }
        
        if(t_flag) write_tefile(outFileName,&einfo,impute_mean_flag);
        else write_efile(outFileName,&einfo,impute_mean_flag);
        
    }
    
     void make_erm(char* outFileName, char* efileName, char* befileName, char* erm_file,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,bool transposed, int efileType,bool no_fid_flag,int valueType, double grm_cutoff, bool erm_cutoff_2sides, bool merge_grm_flag)
    {
        eInfo einfo;
        init_einfo(&einfo);
        
        if(befileName==NULL && efileName==NULL && erm_file==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile or input ORM by --orm.\n");
            TERMINATE();
        }
        if(erm_file!=NULL)
        {
            manipulate_orm(&einfo, erm_file, indilstName, indilst2remove, NULL, grm_cutoff,erm_cutoff_2sides, -2.0, -2, merge_grm_flag, false);
            output_grm(&einfo, outFileName ,erm_bin_flag);
        }
        else
        {
            if(efileName!=NULL)
            {
                if(transposed) read_efile_t(efileName,&einfo,efileType,no_fid_flag,valueType);
                else read_efile(efileName,&einfo,efileType,no_fid_flag,valueType);
                epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
                eii_man(&einfo,indilstName,indilst2remove);
                if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
                if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
            }else{
                char inputname[FNAMESIZE];
                memcpy(inputname,befileName,strlen(befileName)+1);
                char* suffix=inputname+strlen(befileName);
                memcpy(suffix,".oii",5);
                read_eii(inputname,&einfo);
                eii_man(&einfo,indilstName,indilst2remove);
                memcpy(suffix,".opi",5);
                read_epi(inputname,&einfo);
                epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
                if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
                memcpy(suffix,".bod",5);
                read_beed(inputname,&einfo);
                if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
                if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
            }
            if(einfo._eType == METHYLATION)
            {
                if(beta2m && m2beta){
                    //no chance to enter here
                    LOGPRINTF("Error: --m2beta should not be with --beta2m.\n");
                    TERMINATE();
                }
                if(beta2m && einfo._valType==BETAVALUE) beta_2_m(&einfo);
                if(m2beta && einfo._valType==MVALUE) m_2_beta(&einfo);
            }
            make_erm(&einfo, erm_alg, erm_bin_flag,  outFileName, false);
        }
        

    }
    int comp_assoc(const void *a,const void *b){ return ((*(ASSOCRLT *)a).PVAL < (*(ASSOCRLT *)b).PVAL)?1:-1; } //decend
    int ascend_assoc(const void *a,const void *b){ return ((*(ASSOCRLT *)a).PVAL > (*(ASSOCRLT *)b).PVAL)?1:-1; } //ascend
    bool  testMOA(vector<ASSOCRLT> &assoc_rlts, eInfo* einfo, int erm_alg)
    {
        assoc_rlts.clear();
        vector<string> grm_id;
        vector<string> erm_files;
        erm_files.push_back("NA");
        make_erm( einfo,erm_alg);
        for(int i=0; i<einfo->_eii_include.size(); i++) grm_id.push_back(einfo->_eii_fid[einfo->_eii_include[i]]+":"+einfo->_eii_iid[einfo->_eii_include[i]]);
        vector<string> uni_id;
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for(int i=0; i<einfo->_eii_include.size(); i++){
            uni_id.push_back(einfo->_eii_fid[einfo->_eii_include[i]]+":"+einfo->_eii_iid[einfo->_eii_include[i]]);
            uni_id_map.insert(pair<string,int>(einfo->_eii_fid[einfo->_eii_include[i]]+":"+einfo->_eii_iid[einfo->_eii_include[i]], i));
        }
        int _n=(int)einfo->_eii_include.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        einfo->_r_indx.clear();
        vector<MatrixXd> _A;
        vector<int> kp;
        for(int i=0; i < erm_files.size() + 1; i++) einfo->_r_indx.push_back(i);
        _A.resize(einfo->_r_indx.size());
        match(uni_id, grm_id, kp);
        (_A[0]).resize(_n, _n);
#pragma omp parallel for
        for(int i=0; i<_n; i++){
            for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo->_grm(kp[i],kp[j]);
        }
        
        _A[einfo->_r_indx.size()-1]=MatrixXd::Identity(_n, _n);
        
        VectorXd _y; // only for univariate now, using _eii_pheno_num later for multiple variates.
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=einfo->_eii_pheno[einfo->_eii_include[i]];
        }
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(einfo, E_float, qE_float,_X);
        
        // names of variance component
        for (int i = 0; i < erm_files.size(); i++) {
            stringstream strstrm;
            if (erm_files.size() == 1) strstrm << "";
            else strstrm << i + 1;
            einfo->_var_name.push_back("V(O" + strstrm.str() + ")");
            einfo->_hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        einfo->_var_name.push_back("V(e)");
        
        // run REML algorithm
        MatrixXd _Vi;
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        reml(einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, false, true, true, _X_c, _X,_y,_A,_Vi,outfileName);
        if(remlstatus==0 || remlstatus==-5 || remlstatus==-3 )
        {
            einfo->_P.resize(0,0);
            _A.clear();
            unsigned long m=einfo->_epi_include.size();
            
            VectorXd y_buf=_y;
            
                y_buf=_y.array()-(_X*einfo->_b).array(); // adjust phenotype for covariates
                if(einfo->_eii_qcov_num>0 || einfo->_eii_cov_num>0) adjprobe(einfo);
            
            
            VectorXd beta, se, pval;
            
            mlma_calcu_stat(y_buf, einfo, _Vi, beta, se, pval);
            for(int i=0; i<m; i++)
            {
                int j=einfo->_epi_include[i];
                string chrstr;
                if(einfo->_epi_chr[j]==23) chrstr="X";
                else if(einfo->_epi_chr[j]==24) chrstr="Y";
                else chrstr=atosm(einfo->_epi_chr[j]);
                
                ASSOCRLT currlt;
                currlt.BETA=beta[i];
                currlt.SE=se[i];
                if(pval[i]>1.5) pval[i]=1;
                currlt.PVAL=pval[i];
                currlt.CHR=einfo->_epi_chr[j];
                strcpy2(&currlt.GENE, einfo->_epi_gene[j]);
                strcpy2(&currlt.PROBE, einfo->_epi_prb[j]);
                currlt.BP=einfo->_epi_bp[j];
                currlt.OREN=einfo->_epi_orien[j];
                currlt.NMISS=einfo->_eii_include[j];
                assoc_rlts.push_back(currlt);
               
            }
            
        } else {
            return false;
        }
        return true;
    }
    void mlma(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* subtract_erm_file, bool m_erm_flag, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool nopreadj_covar, double percentage_out, double lambda_wind, bool fastlinear, bool force_mlm, bool stdprb)
    {
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        if(priors!=NULL){
            int tmpint=split_string(priors, vs_buf, ",");
            for(int i=0;i<tmpint;i++)
            {
                double tmpd=atof(vs_buf[i].c_str());
                if(tmpd<0 || tmpd>1)
                {
                    LOGPRINTF("Error:  --reml-priors. Prior values of variance explained should be between 0 and 1.\n");
                    TERMINATE();
                }
                reml_priors.push_back(tmpd);
            }
        }
        if(priors_var!=NULL){
            int tmpint=split_string(priors_var, vs_buf, ",");
            for(int i=0;i<tmpint;i++)
            {
                double tmpd=atof(vs_buf[i].c_str());
                if(tmpd<0 )
                {
                    LOGPRINTF("Error:  --reml-priors-var. Prior values of variance components should be over 0.\n");
                    TERMINATE();
                }
                reml_priors_var.push_back(tmpd);
            }
        }
        
        eInfo einfo;
        init_einfo(&einfo);
        einfo._reml_mtd=reml_mtd;
        einfo._reml_max_iter=MaxIter;
        einfo._reml_fixed_var=reml_fixed_var_flag;
        einfo._reml_force_inv=reml_force_inv_fac_flag;
        einfo._reml_force_converge=reml_force_converge_flag;
        einfo._reml_no_converge=reml_no_converge_flag;
       
        
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("Error: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }

        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo); // eii and epi are updated in it.
        if(stdprb) stdprobe(&einfo);
        
        vector<string> grm_id;
        vector<string> erm_files;
        if(percentage_out > 0)
        {
            if(erm_file!=NULL) {
                LOGPRINTF("WARNING: ORM file(s) would be disabled at the presence of --lxpo %f.\n",percentage_out);
            }
            LOGPRINTF("Performing linear regression test ...\n");
            vector<ASSOCRLT> assoc_rlts;
            testQAssoc(assoc_rlts,NULL, &einfo);
            ASSOCRLT* sortptr=&assoc_rlts[0];
            qsort(sortptr,assoc_rlts.size(),sizeof(ASSOCRLT),comp_assoc);
            long numlo=ceil(einfo._epi_include.size()*percentage_out);
            long numslct=einfo._epi_include.size()-numlo;
            vector<string> erm_prbs(numslct);
            for(int i=0;i<numslct;i++) erm_prbs[i]=assoc_rlts[i].PROBE;
            LOGPRINTF("%ld probes are included to make relationship matrix.\n",numslct);
            free_assoclist(assoc_rlts);
            vector<int> include_o(einfo._epi_include);
            map<string, int> snp_name_map_o(einfo._epi_map);
            update_map_kp(erm_prbs, einfo._epi_map, einfo._epi_include);
            erm_files.push_back("NA"); //important
            make_erm(&einfo,erm_alg);
            for(int i=0; i<einfo._eii_include.size(); i++) grm_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
            einfo._epi_include=include_o;
            einfo._epi_map=snp_name_map_o;
            
        } else {
            if(subtract_erm_file){
                erm_files.push_back(erm_file);
                erm_files.push_back(subtract_erm_file);
                for (int i = 0; i < erm_files.size(); i++) {
                    read_grm(&einfo, erm_files[i], grm_id, false, true, true,erm_bin_flag);  // maybe rm grm_id later
                }
            }
            else{
                if(erm_file!=NULL && !m_erm_flag){
                    erm_files.push_back(erm_file);
                    read_grm(&einfo,erm_file, grm_id, true, false, true,erm_bin_flag);
                }
                else if (erm_file!=NULL && m_erm_flag) {
                    read_msglist(erm_file, erm_files,"ORM file names");
                    for (int i = 0; i < erm_files.size(); i++) {
                        read_grm(&einfo,erm_files[i], grm_id, false, true, true,erm_bin_flag);
                    }
                }
                else{
                    erm_files.push_back("NA");
                    make_erm( &einfo,erm_alg);
                    for(int i=0; i<einfo._eii_include.size(); i++) grm_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
                }
            }
        }
        
       
        vector<string> uni_id;
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for(int i=0; i<einfo._eii_include.size(); i++){
            uni_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
            uni_id_map.insert(pair<string,int>(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]], i));
        }
        int _n=(int)einfo._eii_include.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        einfo._r_indx.clear();
        vector<MatrixXd> _A;
        vector<int> kp;

        if(percentage_out > 0)
        {
            for(int i=0; i < erm_files.size() + 1; i++) einfo._r_indx.push_back(i);
            _A.resize(einfo._r_indx.size());
            match(uni_id, grm_id, kp);
            (_A[0]).resize(_n, _n);
#pragma omp parallel for
            for(int i=0; i<_n; i++){
                for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(kp[i],kp[j]);
            }
            //einfo._grm.resize(0,0);

        } else {
            if (subtract_erm_file) {
                for(int i=0; i < 2; i++) einfo._r_indx.push_back(i);
                _A.resize(einfo._r_indx.size());
                
                LOGPRINTF("\nReading the primary ORM from %s ...\n",erm_files[1].c_str() );
                read_grm(&einfo, erm_files[1], grm_id, true, false, false,erm_bin_flag);
                
                match(uni_id, grm_id, kp);
                (_A[0]).resize(_n, _n);
                MatrixXf A_N_buf(_n, _n);
#pragma omp parallel for
                for (int j = 0; j < _n; j++) {
                    for (int k = 0; k <= j; k++) {
                        if (kp[j] >= kp[k]){
                            (_A[0])(k, j) = (_A[0])(j, k) = einfo._grm(kp[j], kp[k]);
                            A_N_buf(k, j) = A_N_buf(j, k) = einfo._grm_N(kp[j], kp[k]);
                        }
                        else{
                            (_A[0])(k, j) = (_A[0])(j, k) = einfo._grm(kp[k], kp[j]);
                            A_N_buf(k, j) = A_N_buf(j, k) = einfo._grm_N(kp[k], kp[j]);
                        }
                    }
                }
                
                LOGPRINTF("\nReading the secondary ORM from %s ...\n",erm_files[0].c_str());
                read_grm(&einfo, erm_files[0], grm_id, true, false, false,erm_bin_flag);
                LOGPRINTF("\nSubtracting %s from %s ...\n",erm_files[1].c_str(),erm_files[0].c_str());
                match(uni_id, grm_id, kp);
#pragma omp parallel for
                for (int j = 0; j < _n; j++) {
                    for (int k = 0; k <= j; k++) {
                        if (kp[j] >= kp[k]) (_A[0])(k, j) = (_A[0])(j, k) = ((_A[0])(j, k) * A_N_buf(j, k)  - einfo._grm(kp[j], kp[k]) * einfo._grm_N(kp[j], kp[k])) / (A_N_buf(j, k) - einfo._grm_N(kp[j], kp[k]));
                        else (_A[0])(k, j) = (_A[0])(j, k) = ((_A[0])(j, k) * A_N_buf(j, k) - einfo._grm(kp[k], kp[j]) * einfo._grm_N(kp[k], kp[j])) / (A_N_buf(j, k) - einfo._grm_N(kp[k], kp[j]));
                    }
                }
                einfo._grm.resize(0,0);
                einfo._grm_N.resize(0,0);
            }
            else {
                for(int i=0; i < erm_files.size() + 1; i++) einfo._r_indx.push_back(i);
                _A.resize(einfo._r_indx.size());
                if(erm_file!=NULL && !m_erm_flag){
                    match(uni_id, grm_id, kp);
                    (_A[0]).resize(_n, _n);
#pragma omp parallel for
                    for(int i=0; i<_n; i++){
                        for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(kp[i],kp[j]);
                    }
                    einfo._grm.resize(0,0);
                }
                else if(erm_file!=NULL && m_erm_flag){
                    LOGPRINTF("There are %ld ORM file names specified in the file %s.\n",erm_files.size(),erm_file);
                    for (int i = 0; i < erm_files.size(); i++) {
                        LOGPRINTF("Reading the ORM from the %dth file ...\n",i + 1);
                        read_grm(&einfo, erm_files[i], grm_id, true, false, true,erm_bin_flag);
                        match(uni_id, grm_id, kp);
                        (_A[i]).resize(_n, _n);
#pragma omp parallel for
                        for (int j = 0; j < _n; j++) {
                            for (int k = 0; k <= j; k++) {
                                if (kp[j] >= kp[k]) (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[j], kp[k]);
                                else (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[k], kp[j]);
                            }
                        }
                    }
                }
                else{
                    match(uni_id, grm_id, kp);
                    (_A[0]).resize(_n, _n);
#pragma omp parallel for
                    for(int i=0; i<_n; i++){
                        for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(kp[i],kp[j]);
                    }
                    //einfo._grm.resize(0,0);
                }
            }
        }
       _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
        
        VectorXd _y; // only for univariate now, using _eii_pheno_num later for multiple variates.
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=einfo._eii_pheno[einfo._eii_include[i]];
        }
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(&einfo, E_float, qE_float,_X); //_X has the dimension of [einfo._eii_include.size(),_X_c]

        // names of variance component
        for (int i = 0; i < erm_files.size(); i++) {
            stringstream strstrm;
            if (erm_files.size() == 1) strstrm << "";
            else strstrm << i + 1;
            einfo._var_name.push_back("V(O" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        einfo._within_family=within_family;
        if(within_family) detect_family(&einfo, _A);
        
        // run REML algorithm
        LOGPRINTF("\nPerforming MOA analyses %s ...\n",(subtract_erm_file?"":" (including the candidate probes)"));
        MatrixXd _Vi;
        loud = true;
        reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
        /******/
        /*
        string filename=string(outfileName)+".cor.mat";
        FILE* tmpfile=fopen(filename.c_str(),"w");
        if(!tmpfile)
        {
            LOGPRINTF("error open file.\n");
            TERMINATE();
        }
        for(int t=0;t<_Vi.rows();t++)
        {
            string str="";
            for(int k=0;k<_Vi.cols();k++)
            {
                str +=atos(_Vi(t,k)) + '\t';
            }
            str += '\n';
            fputs(str.c_str(),tmpfile);
        }
        fclose(tmpfile);
         */
        /*****/
        if(remlstatus == 0 || remlstatus == -3 || remlstatus == -5 )
        {
            einfo._P.resize(0,0);
            _A.clear();
            unsigned long m=einfo._epi_include.size();
            
            VectorXd y_buf=_y;
            if(!nopreadj_covar){
                y_buf=_y.array()-(_X*einfo._b).array(); // adjust phenotype for covariates
                if(einfo._eii_qcov_num>0 || einfo._eii_cov_num>0) adjprobe(&einfo);
            }
            
            VectorXd beta, se, pval;
            if(!nopreadj_covar) mlma_calcu_stat(y_buf, &einfo, _Vi, beta, se, pval);
            else mlma_calcu_stat_covar(y_buf, &einfo, _X_c, _Vi, _X, beta, se, pval);
            
            string filename=string(outFileName)+".moa";
            LOGPRINTF("\nSaving the association analysis results of %ld probes to %s ...\n",m,filename.c_str());
            ofstream ofile(filename.c_str());
            if(!ofile) throw("Can not open the file ["+filename+"] to write.");
            ofile<<"Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp"<<endl;
            for(int i=0; i<m; i++){
                int j=einfo._epi_include[i];
                string chrstr;
                if(einfo._epi_chr[j]==23) chrstr="X";
                else if(einfo._epi_chr[j]==24) chrstr="Y";
                else chrstr=atosm(einfo._epi_chr[j]);
                
                ofile<<chrstr<<"\t"<<einfo._epi_prb[j]<<"\t"<<einfo._epi_bp[j]<<"\t"<<einfo._epi_gene[j]<<"\t"<<(einfo._epi_orien[j]=='*'?"NA":atos(einfo._epi_orien[j]))<<"\t";
                if(pval[i]>1.5) ofile<<"NA\tNA\tNA"<<endl;
                else ofile<<beta[i]<<"\t"<<se[i]<<"\t"<<pval[i]<<endl;
            }
            ofile.close();

        }
        else
        {
            if(remlstatus==-1) {
                LOGPRINTF("ERROR: The matrix is not invertible.\n");
            } else if(remlstatus==-2) {
                LOGPRINTF("ERROR: More than half of the variance components are constrained.\n");
            }else if(remlstatus==-4) {
                LOGPRINTF("ERROR:Log-likelihood not converged.\n");
            }
            TERMINATE();
        }
    }
    
    void output_simu_par(char* outFileName, vector<string> &qtl_name, vector<double> &qtl_eff)
    {
        int i = 0;
        string out_parfile = string(outFileName) + ".par";
        ofstream out_par(out_parfile.c_str());
        if (!out_par)
        {
            LOGPRINTF("Error: can not open par file %s to write!\n",out_parfile.c_str());
            TERMINATE();
        }
        out_par << "Causal\tEffect" << endl;
        for (i = 0; i < qtl_eff.size(); i++) out_par << qtl_name[i] << "\t" << qtl_eff[i] << endl;
        out_par.close();
        LOGPRINTF("Simulated Causal effect(s) have been saved in %s.\n",  out_parfile.c_str() );
    }
    void output_simu_par(char* outFileName, vector< vector<string> > &qtl_name, vector< vector<double> > &qtl_eff)
    {
        string out_parfile = string(outFileName) + ".par";
        ofstream out_par(out_parfile.c_str());
        if (!out_par) {
            LOGPRINTF("Error: can not open par file %s to write!\n",out_parfile.c_str());
            TERMINATE();
        }
        out_par << "Causal\tEffect" << endl;
        
        if(qtl_name.size()!=qtl_eff.size())
        {
            LOGPRINTF("ERROR: different dimension of causal probe and causal effect. please report.\n");
            TERMINATE();
        }
        for(int i=0;i<qtl_name.size();i++)
            for(int j=0;j<qtl_name[i].size();j++) out_par << qtl_name[i][j] << "\t" << qtl_eff[i][j] << endl;
        out_par.close();
        cout << "Simulated Causal effect(s) have been saved in [" + out_parfile + "]." << endl;
    }

    void save_phenfile(eInfo* einfo,char* outFileName, vector< vector<double> > &y)
    {
        string phenfile = string(outFileName) + ".phen";
        ofstream phen(phenfile.c_str());
        if (!phen) throw ("Error: can not open the file [" + phenfile + "] to write.");
        int i = 0, j = 0;
        for (i = 0; i < einfo->_eii_include.size(); i++) {
            phen << einfo->_eii_fid[einfo->_eii_include[i]] << " " << einfo->_eii_iid[einfo->_eii_include[i]] << " ";
            for (j = 0; j < y.size(); j++) phen << y[j][i] << " ";
            phen << endl;
        }
        phen.close();
    }
    
    void diff(char* befileName1,char* befileName2)
    {
        if(befileName1==NULL || befileName2==NULL)
        {
            LOGPRINTF("please input gene expression / methylation data for diff by the flag --befile.\n");
            TERMINATE();
        }
        map<string,int>::iterator iter;
        eInfo edata1;
        eInfo edata2;
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName1,strlen(befileName1)+1);
        char* suffix=inputname+strlen(befileName1);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&edata1);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&edata1);
        memcpy(suffix,".bod",5);
        FILE *fptr1=fopen(inputname, "rb");
        if(!fptr1)
        {
            LOGPRINTF("ERROR: Couldn't open file %s\n", inputname);
            TERMINATE();
        }
        memcpy(inputname,befileName2,strlen(befileName2)+1);
        suffix=inputname+strlen(befileName2);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&edata2);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&edata2);
        memcpy(suffix,".bod",5);
        FILE *fptr2=fopen(inputname, "rb");
        if(!fptr2)
        {
            LOGPRINTF("ERROR: Couldn't open file %s\n", inputname);
            TERMINATE();
        }
        for(int i=0;i<edata1._eii_num;i++)
        {
            //.oii is sorted by FID:IID
            string tmpstr1=edata1._eii_fid[i]+":"+edata1._eii_iid[i];
            string tmpstr2=edata2._eii_fid[i]+":"+edata2._eii_iid[i];
            if(tmpstr1.compare(tmpstr2.c_str())) {
                LOGPRINTF("ERROR: no consistence in .oii files.\n");
                TERMINATE();
            }
        }
        for(int i=0;i<edata1._epi_num;i++)
        {
            // .opi is sorted by PID or <chr,pos>
            if(edata2._epi_map.find(edata1._epi_prb[i].c_str())==edata2._epi_map.end()) {
                LOGPRINTF("ERROR: can't find probe %s in file %s.\n",edata1._epi_prb[i].c_str(),befileName2);
                TERMINATE();
            }
        }
        char* readBuf1;
        string msg="Reading Buffer";
        if(!allocReserved(&readBuf1, edata1._eii_num*sizeof(double),msg)) TERMINATE();
        char* readBuf2;
        if(!allocReserved(&readBuf2, edata2._eii_num*sizeof(double),msg)) TERMINATE();
        
        
        uint32_t tmpuint1=0, tmpuint2=0;
        for(int i=0;i<3;i++)
        {
            if(fread(&tmpuint1, sizeof(uint32_t),1, fptr1)!=1)
            {
                LOGPRINTF("ERROR: File %s read failed!\n", befileName1);
                TERMINATE();
            }
            if(fread(&tmpuint2, sizeof(uint32_t),1, fptr2)!=1)
            {
                LOGPRINTF("ERROR: File %s read failed!\n", befileName2);
                TERMINATE();
            }
            if(tmpuint1!=tmpuint2)
            {
                if(i==0)
                {
                    LOGPRINTF("ERROR: no consistence in file type.\n");
                    TERMINATE();
                } else if(i==1) {
                    LOGPRINTF("ERROR: no consistence in individual number.\n");
                    TERMINATE();
                } else {
                    LOGPRINTF("ERROR: no consistence in probe number.\n");
                    TERMINATE();
                }
            }
        }
        uint64_t rsize=edata1._eii_num*sizeof(double);
        for(int i=0;i<edata1._epi_num;i++)
        {
            memset(readBuf1,0, rsize);
            memset(readBuf2,0, rsize);
            if(fread(readBuf1, 1,rsize, fptr1)!=rsize)
            {
                LOGPRINTF("ERROR: File %s read failed!\n", befileName1);
                TERMINATE();
            }
            iter=edata2._epi_map.find(edata1._epi_prb[i].c_str());
            if(iter!=edata2._epi_map.end())
            {
                
                uint64_t readpos=3*sizeof(uint32_t)+iter->second*edata2._eii_num*sizeof(double);
                fseek( fptr2, readpos, SEEK_SET );
                if(fread(readBuf2, 1,rsize, fptr2)!=rsize)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", befileName2);
                    TERMINATE();
                }
            }
            double* ptrtmp1=(double *)readBuf1;
            double* ptrtmp2=(double *)readBuf2;
            
            for(int j=0;j<edata1._eii_num;j++)
            {
                if(abs(ptrtmp1[j]-ptrtmp2[j])>1e-6)
                {
                    LOGPRINTF("ERROR: no consistence in values of probe %s and individual %s.\n",edata1._epi_prb[i].c_str(),(edata1._eii_fid[j]+":"+edata1._eii_iid[j]).c_str());
                    TERMINATE();
                }
            }
        }
        deallocReserved(&readBuf1, edata1._eii_num*sizeof(double));
        deallocReserved(&readBuf2, edata2._eii_num*sizeof(double));
        LOGPRINTF("PASSED: the files identify with each other.\n");
    }
    
    void getRefactor(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* covfileName,char* qcovfileName, int celltype_num, int dmr_num, int out_pc_num)
    {
        setNbThreads(thread_num);
        bool rawpaperpcs=false;
        eInfo einfo;
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(celltype_num<=0)
        {
            LOGPRINTF("Error: please input the cell type number by the option --celltype-num.\n");
            TERMINATE();
        }
        if(dmr_num<=0)
        {
            LOGPRINTF("Error: please input the DMR number by the option --dmr-num.\n");
            TERMINATE();
        }
        if(out_pc_num<=0)
        {
            out_pc_num=celltype_num;
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        
        //Adjust the data for the covariates
        //add later
        //end of adjustment
        
        
        bool divid_by_std = true;
        if(divid_by_std){
            LOGPRINTF("Centering and standardizing probes...\n");
        } else {
            LOGPRINTF("Centering probes...\n");
        }
        vector< vector<bool> > X_bool;
        MatrixXd _probe_data;
        std_probe(&einfo, X_bool, divid_by_std, _probe_data);
        
        LOGPRINTF("Generating low rank approximation matrix...\n");
        JacobiSVD<MatrixXd> svd(_probe_data, ComputeThinU | ComputeThinV);
        MatrixXd approx;
        approx=(svd.matrixU()*svd.singularValues().asDiagonal()).leftCols(celltype_num)*svd.matrixV().leftCols(celltype_num).transpose();
        if(divid_by_std){
            LOGPRINTF("Centering and standardizing probes of rank approximation matrix...\n");
        } else {
            LOGPRINTF("Centering probes of rank approximation matrix...\n");
        }
        std_probe_in_place( divid_by_std, approx);
        
        LOGPRINTF("Caculating the distance between the raw matrix and the low rank approximation matrix...\n");
        approx=_probe_data-approx;
        approx=approx.array()*approx.array();
        VectorXd dis=approx.colwise().sum();
        vector<double> distance(dis.size());
        for(int i=0;i<dis.size();i++) distance[i]=dis(i);
        vector<int> indics;
        get_bottom_indices(indics, distance, dmr_num);
        stable_sort(indics.begin(), indics.end());
        
        LOGPRINTF("Getting the sparse matrix %ld by %ld...\n",_probe_data.rows(),indics.size());
        approx.resize(_probe_data.rows(),indics.size());
        for(int i=0;i<indics.size();i++)
            approx.col(i)=_probe_data.col(indics[i]);
        if(!rawpaperpcs)
        {
            LOGPRINTF("Computing the Refactor PCs...\n");
            /***method of GCTA****/
            _probe_data=approx*approx.transpose();
            
            // JacobiSVD has the same effect as SelfAdjointEigenSolver
            //JacobiSVD<MatrixXd> svd3(_probe_data, ComputeThinU | ComputeThinV);
            //cout<<svd3.matrixU().block(0,0,5,5)<<endl;
            
            SelfAdjointEigenSolver<MatrixXd> eigensolver(_probe_data);
            MatrixXd evec = (eigensolver.eigenvectors());
            VectorXd eval = eigensolver.eigenvalues();
            //cout<<evec.block(0,evec.cols()-5,5,5)<<endl;
            
            FILE* efile=NULL;
            string evalName=string(outFileName)+".eigenval";
            if(fopen_checked(&efile, evalName.c_str(),"w")) TERMINATE();
            for (long i = einfo._eii_num - 1; i >= 0; i--) {
                string str=atos(eval(i))+'\n';
                if(fputs_checked(str.c_str(),efile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", evalName.c_str());
                    TERMINATE();
                }
                
            }
            fclose(efile);
            LOGPRINTF("Eigenvalues of %llu  individuals have been saved in %s.\n",einfo._eii_num ,evalName.c_str());
            
            evalName = string(outFileName) + ".eigenvec";
            if(fopen_checked(&efile, evalName.c_str(),"w")) TERMINATE();
            
            if (out_pc_num > einfo._eii_num) out_pc_num = einfo._eii_num;
            for (int i = 0; i < einfo._eii_num; i++) {
                string str=einfo._eii_fid[i]+'\t'+einfo._eii_iid[i];
                
                for (long j = einfo._eii_num - 1; j >= (einfo._eii_num - out_pc_num); j--) str += '\t'+atos(evec(i, j));
                str+='\n';
                if(fputs_checked(str.c_str(),efile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", evalName.c_str());
                    TERMINATE();
                }
            }
            fclose(efile);
            LOGPRINTF("The first %d eigenvectors of %llu individuals have been saved in %s.\n",out_pc_num ,einfo._eii_num, evalName.c_str());
            /***/

        } else {
            /**method of Elior Rahmani's paper**/
            
            JacobiSVD<MatrixXd> svd2(approx, ComputeThinU | ComputeThinV);
            
            //cout<<svd2.matrixU().block(0,0,5,5)<<endl;  // has some difference in signs from svd() in R
            //cout<<svd2.singularValues()<<endl; // same as svd()$d in R
            
            MatrixXd pcs=(svd2.matrixU()*svd2.singularValues().asDiagonal()).leftCols(celltype_num);
            //cout<<pcs.block(0,0,5,5)<<endl;
             FILE* efile=NULL;
            string evalName = string(outFileName) + ".eigenvec";
            if(fopen_checked(&efile, evalName.c_str(),"w")) TERMINATE();
            
            if (out_pc_num > einfo._eii_num) out_pc_num = einfo._eii_num;
            for (int i = 0; i < einfo._eii_num; i++) {
                string str=einfo._eii_fid[i]+'\t'+einfo._eii_iid[i];
                
                for (int j =0; j< out_pc_num; j++) str += '\t'+atos(pcs(i, j));
                str+='\n';
                if(fputs_checked(str.c_str(),efile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", evalName.c_str());
                    TERMINATE();
                }
            }
            fclose(efile);
            LOGPRINTF("The first %d eigenvectors of %llu individuals have been saved in %s.\n",out_pc_num ,einfo._eii_num, evalName.c_str());
            
            /***/
        }
       
    }
    
    void mlma_loco(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, int MaxIter, char* priors, char* priors_var, bool no_constrain, bool no_adj_covar,int reml_mtd,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, int autosome_num, double percentage_out, double lambda_wind,bool fastlinear)
    {
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        if(priors!=NULL){
            int tmpint=split_string(priors, vs_buf, ",");
            for(int i=0;i<tmpint;i++)
            {
                double tmpd=atof(vs_buf[i].c_str());
                if(tmpd<0 || tmpd>1)
                {
                    LOGPRINTF("Error:  --reml-priors. Prior values of variance explained should be between 0 and 1.\n");
                    TERMINATE();
                }
                reml_priors.push_back(tmpd);
            }
        }
        if(priors_var!=NULL){
            int tmpint=split_string(priors_var, vs_buf, ",");
            for(int i=0;i<tmpint;i++)
            {
                double tmpd=atof(vs_buf[i].c_str());
                if(tmpd<0 )
                {
                    LOGPRINTF("Error:  --reml-priors-var. Prior values of variance components should be over 0.\n");
                    TERMINATE();
                }
                reml_priors_var.push_back(tmpd);
            }
        }

        eInfo einfo;
        init_einfo(&einfo);
        einfo._reml_mtd=reml_mtd;
        einfo._reml_max_iter=MaxIter;
        einfo._reml_fixed_var=reml_fixed_var_flag;
        einfo._reml_force_inv=reml_force_inv_fac_flag;
        einfo._reml_force_converge=reml_force_converge_flag;
        einfo._reml_no_converge=reml_no_converge_flag;

        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("Error: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }

        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        if (covfileName == NULL && qcovfileName == NULL) no_adj_covar=false;
        
        int _n=(int)einfo._eii_include.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        
        vector<int> include_raw(einfo._epi_include);
        map<string, int> snp_name_map_raw(einfo._epi_map);
        if(percentage_out > 0)
        {
            LOGPRINTF("Performing linear regression test ...\n");
            vector<ASSOCRLT> assoc_rlts;
            testQAssoc(assoc_rlts,NULL, &einfo);
            ASSOCRLT* sortptr=&assoc_rlts[0];
            qsort(sortptr,assoc_rlts.size(),sizeof(ASSOCRLT),comp_assoc);
            long numlo=ceil(einfo._epi_include.size()*percentage_out);
            long numslct=einfo._epi_include.size()-numlo;
            vector<string> erm_prbs(numslct);
            for(int i=0;i<numslct;i++) erm_prbs[i]=assoc_rlts[i].PROBE;
            LOGPRINTF("%ld probes are included to make relationship matrix.\n",numslct);
            update_map_kp(erm_prbs, einfo._epi_map, einfo._epi_include);
            free_assoclist(assoc_rlts);
        }
        
        vector<string> grm_id;
        vector<int> chrs, vi_buf(einfo._epi_chr);
        stable_sort(vi_buf.begin(), vi_buf.end());
        vi_buf.erase(unique(vi_buf.begin(), vi_buf.end()), vi_buf.end());
        if(vi_buf.size()<2) {
            LOGPRINTF("Error: There is only one chromosome. The MLM leave-on-chromosome-out (LOCO) analysis requires at least two chromosomes.\n");
            TERMINATE();
        }
        
        for(int i=0; i<vi_buf.size(); i++){
            if(vi_buf[i]<=autosome_num && vi_buf[i]>0) chrs.push_back(vi_buf[i]); // remove
            else {
                if(vi_buf[i]>autosome_num) {
                    string chrstr;
                    if(vi_buf[i]==23) chrstr="X";
                    else if(vi_buf[i]==24) chrstr="Y";
                    else chrstr=atos(vi_buf[i]);
                    LOGPRINTF("\nWARNING: Chromosome %s is out of autosomes (%d). This chromosome would be excluded from this analysis.\n",chrstr.c_str(),autosome_num);
                    LOGPRINTF("Please use --autosome-num to update the autosome number.\n");
                } else if(vi_buf[i] == -9) {
                    LOGPRINTF("\nWARNING: Missing chromosome (NA) found. This chromosome would be excluded from this analysis.\n");
                } else if(vi_buf[i] == 0) {
                    LOGPRINTF("\nWARNING: Unrecognized chromosome (Chr1,M,*,etc.) found. This chromosome would be excluded from this analysis.\n");
                }
            }
        }
        vector<int> include_o(einfo._epi_include);
        map<string, int> snp_name_map_o(einfo._epi_map);
        vector<float> m_chrs_f(chrs.size());
        vector<double *> grm_chrs(chrs.size());
        vector<double *> geno_chrs(chrs.size());
        vector< vector<int> > icld_chrs(chrs.size());
        LOGPRINTF("\nCalculating the omics relationship matrix for each of the %ld chromosomes ... \n",chrs.size());
        for(int c1=0; c1<chrs.size(); c1++)
        {
            cout<<"Chr "<<chrs[c1]<<":"<<endl;
            extract_probe_by_chr(&einfo,chrs[c1]);
            LOGPRINTF("%ld probes on chromosome %d are included to make ORM.\n",einfo._epi_include.size(), chrs[c1]);
            make_erm(&einfo, 0, true, NULL, true); // output the profile
            /**test**/
            //for(int dd=0;dd<20;dd++) printf("%f,",einfo._profile[dd*einfo._epi_include.size()+0]);
            //printf("\n");
            //for(int dd=0;dd<15;dd++) printf("%f,",einfo._grm_ptr[dd*einfo._eii_include.size()+0]);
            //printf("\n");
            /***/
            if(percentage_out > 0){
                m_chrs_f[c1]=(float)einfo._epi_include.size(); // used to calculate matrix A
                einfo._epi_include=include_raw;
                einfo._epi_map=snp_name_map_raw;
                delete einfo._profile;
                einfo._profile=NULL;
                vector< vector<bool> > X_bool;
                MatrixXd prbdata;
                extract_probe_by_chr(&einfo,chrs[c1]);
                std_probe(&einfo, X_bool, false, prbdata, true); // to out the profile. in the profile miss value is 0.0
                icld_chrs[c1]=einfo._epi_include; //for output
                geno_chrs[c1]=einfo._profile;
                LOGPRINTF("%ld probes on chromosome %d are ready to test.\n",einfo._epi_include.size(), chrs[c1]);
            } else {
                m_chrs_f[c1]=(float)einfo._epi_include.size();
                icld_chrs[c1]=einfo._epi_include;
                geno_chrs[c1]=einfo._profile; // raw nonmissing profile minus the mean. the miss was set to 0. dimension:  einfo._eii_include.size() * probe number current chromosome
            }
            grm_chrs[c1]=einfo._grm_ptr;
            einfo._epi_include=include_o;
            einfo._epi_map=snp_name_map_o;
            einfo._profile=NULL;
            einfo._grm_ptr=NULL;
        }
        for(int i=0; i<einfo._eii_include.size(); i++) grm_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
        
        vector<string> uni_id;
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for(int i=0; i<einfo._eii_include.size(); i++){
            uni_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
            uni_id_map.insert(pair<string,int>(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]], i));
        }
        
        // construct model terms
        VectorXd _y; // only for univariate now, using _eii_pheno_num later for multiple variates.
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=einfo._eii_pheno[einfo._eii_include[i]];
        }
        

        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(&einfo, E_float, qE_float,_X);
        
        
        // names of variance component
        einfo._var_name.push_back("V(O)");
        einfo._hsq_name.push_back("V(O)/Vp");
        einfo._var_name.push_back("V(e)");
        
        // MLM association
        LOGPRINTF("\nPerforming MLM association analyses (leave-one-chromosome-out) ...\n");
        
        vector<MatrixXd> _A;
        vector<int> kp;
        match(uni_id, grm_id, kp); // oh, man, uni_id is the same as grm_id
        einfo._r_indx.resize(2);
        for(int i=0; i<2; i++) einfo._r_indx[i]=i;
        _A.resize(einfo._r_indx.size());
        _A[1]=MatrixXd::Identity(_n, _n);
        
        VectorXd y_buf=_y;
        float *y=new float[_n];
        vector<VectorXd> beta(chrs.size()), se(chrs.size()), pval(chrs.size());
        for(int c1=0; c1<chrs.size(); c1++)
        {
            if(percentage_out > 0) {
                einfo._epi_include=include_raw;
                einfo._epi_map=snp_name_map_raw;
            }
            LOGPRINTF("\n-----------------------------------\n#Chr %d:\n",chrs[c1]);
            extract_probe_by_chr(&einfo,chrs[c1]);
            _A[0]=MatrixXd::Zero(_n, _n);
            double d_buf=0;
            for(int c2=0; c2<chrs.size(); c2++)
            {
                if(chrs[c1]==chrs[c2]) continue;
#pragma omp parallel for
                for(int i=0; i<_n; i++)
                {
                    for(int j=0; j<=i; j++)
                    {
                        (_A[0])(i,j)+=(grm_chrs[c2])[kp[i]*_n+kp[j]]*m_chrs_f[c2];
                    }
                }
                d_buf+=m_chrs_f[c2];
            }
            
#pragma omp parallel for
            for(int i=0; i<_n; i++){
                for(int j=0; j<=i; j++){
                    (_A[0])(i,j)/=d_buf;
                    (_A[0])(j,i)=(_A[0])(i,j);
                }
            }
            
            // run REML algorithm
            MatrixXd _Vi;
            reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
            if(remlstatus!=0) break;
            /**test**/
            /*
            cout<<"_X_c: "<<_X_c<<endl;
            cout<<"-----X----"<<endl;
            for(int k=0;k<15;k++)
            {
                for(int l=0;l<_X_c;l++) cout<<_X(k,l)<<",";
                cout<<endl;
            }
            cout<<"-----y----"<<endl;
            for(int k=0;k<15;k++) cout<<_y(k)<<",";
            cout<<endl;
            cout<<"-----Vi----"<<endl;
            MatrixXd cao;
            cao=_Vi.block(0,0,15,15);
            cout<<cao<<endl;
             */
            /****/
            if(!no_adj_covar){
                y_buf=_y.array()-(_X*einfo._b).array(); // adjust phenotype for covariates
                if(einfo._eii_qcov_num>0 || einfo._eii_cov_num>0) adjprobe(&einfo);
            }
           
            for(int i=0; i<_n; i++) y[i]=y_buf[i];
            reml_priors.clear();
            reml_priors_var=einfo._varcmp;
            einfo._P.resize(0,0);
            _A[0].resize(0,0);
            unsigned long n=einfo._eii_include.size(), m=einfo._epi_include.size();
            mlma_calcu_stat(y_buf, (geno_chrs[c1]), n, m, _Vi, beta[c1], se[c1], pval[c1]);
        
            einfo._epi_include=include_o;
            einfo._epi_map=snp_name_map_o;
            LOGPRINTF("-----------------------------------\n");
        }
        
        delete[] y;
        for(int c1=0; c1<chrs.size(); c1++)
        {
            delete[] (grm_chrs[c1]);
            delete[] (geno_chrs[c1]);
        }
        if(remlstatus==0)
        {
            string filename=string(outFileName)+".loco.moa";
            LOGPRINTF("\nSaving the results of the mixed linear model association analyses of %ld probes to %s ...\n",include_raw.size(),filename.c_str());
            ofstream ofile(filename.c_str());
            if(!ofile) throw("Can not open the file ["+filename+"] to write.");
            ofile<<"Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp"<<endl;
            for(int c1=0; c1<chrs.size(); c1++){
                for(int i=0; i<icld_chrs[c1].size(); i++){
                    int j=icld_chrs[c1][i];
                    string strand=((einfo._epi_orien[j]=='*')?"NA":atos(einfo._epi_orien[j]));
                    string chrstr;
                    if(einfo._epi_chr[j]==23) chrstr="X";
                    else if(einfo._epi_chr[j]==24) chrstr="Y";
                    else chrstr=atosm(einfo._epi_chr[j]);
                    ofile<<chrstr<<"\t"<<einfo._epi_prb[j]<<"\t"<<einfo._epi_bp[j]<<"\t"<<einfo._epi_gene[j]<<"\t"<<strand<<"\t";
                    if(pval[c1][i]>1.5) ofile<<"NA\tNA\tNA\tNA"<<endl;
                    else ofile<<beta[c1][i]<<"\t"<<se[c1][i]<<"\t"<<pval[c1][i]<<endl;
                }
            }
            ofile.close();

        }
        else
        {
            einfo._epi_include=include_o;
            einfo._epi_map=snp_name_map_o;
            if(remlstatus==-1) {
                LOGPRINTF("\nThe matrix is not invertible, ");
            } else if(remlstatus==-2) {
                LOGPRINTF("\nMore than half of the variance components are constrained, ");
            } else if(remlstatus==-3) {
                LOGPRINTF("\nVariance component going to 0 or 1, ");
            }else if(remlstatus==-4) {
                LOGPRINTF("\nLog-likelihood not converged, ");
            }
            LOGPRINTF("standard PCA based linear regression actived \n");
            double upperlambda=1+lambda_wind;
            double lowerlambda=(1-lambda_wind)>0?(1-lambda_wind):0;
            LOGPRINTF("Iteratively identify the number of PCs to achieve lambda between the region [%f,%f].\n ", lowerlambda,upperlambda);
            LOGPRINTF("\nPerforming principal component analysis ...\n");
            make_erm( &einfo);
            SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
            MatrixXd evec = (eigensolver.eigenvectors());
            VectorXd eval = eigensolver.eigenvalues();
            long maxpc2test = evec.cols()/2;
            long pc2test = (maxpc2test>2)?2:maxpc2test;
            long regionL = pc2test, regionR=pc2test;
            char outputname[FNAMESIZE];
            outputname[0]='\0';
            double outlambda=1e6;
            vector<ASSOCRLT> out_rlts;
            while(pc2test <= maxpc2test)
            {
                LOGPRINTF("\nPerforming association analysis with %ld PCs ...\n",pc2test);
                if(outFileName!=NULL) {
                    string tmp=  string(outFileName)+"_PC"+atos(pc2test);
                    strcpy(outputname,tmp.c_str());
                }
                MatrixXd _X_PC = evec.block(0,evec.cols()-pc2test, evec.rows(),pc2test);
                /****/
                if(_X_PC.cols()>0)
                {
                    LOGPRINTF("Saving the %ld PCs...\n",pc2test);
                    string filename=string(outfileName)+"."+atos(pc2test)+"PC.eigenvec";
                    FILE* tmpfile=fopen(filename.c_str(),"w");
                    if(!tmpfile)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<einfo._eii_include.size();t++)
                    {
                        string str=einfo._eii_fid[einfo._eii_include[t]]+'\t'+einfo._eii_iid[einfo._eii_include[t]];
                        for(long k=(_X_PC.cols()-1);k>=0;k--)
                        {
                            str +='\t'+atos(_X_PC(t,k));
                        }
                        str += '\n';
                        fputs(str.c_str(),tmpfile);
                    }
                    fclose(tmpfile);
                    LOGPRINTF("These PCs are saved in the file %s.\n",filename.c_str());
                }
                /****/
                vector<ASSOCRLT> assoc_rlts;
                if(fastlinear) testLinear_fast(assoc_rlts,outputname, &einfo,_X_PC);
                else testLinear(assoc_rlts,outputname, &einfo,_X_PC);
                double lambda=get_lambda(assoc_rlts);
                LOGPRINTF("The genomic inflation factor lambda with %ld PCs is %f.\n",pc2test, lambda);
                if(abs(1-lambda) < abs(1-outlambda)) {
                    out_rlts = assoc_rlts;
                    outlambda = lambda;
                }
                if(lambda>upperlambda)
                {
                    regionL=pc2test;
                    if(regionR==pc2test)
                    {
                        pc2test*=2;
                        regionR=pc2test;
                    } else {
                        pc2test=(regionR+regionL)/2;
                    }
                    if(pc2test==regionL || pc2test>maxpc2test)
                    {
                        write_assoc_rlt(out_rlts, outFileName);
                        break;
                    }
                }
                else if(lambda<lowerlambda)
                {
                    regionR=pc2test;
                    if(regionL==pc2test)
                    {
                        pc2test/=2;
                        regionL=pc2test;
                    } else {
                        pc2test=(regionR+regionL)/2;
                    }
                    if(pc2test==regionR)
                    {
                        write_assoc_rlt(out_rlts, outFileName);
                        break;
                    }
                    
                }
                else
                {
                    write_assoc_rlt(out_rlts, outFileName);
                    break;
                }
            }
        }
    }
    void manipulate_orm(eInfo* einfo, char* grm_file, char* indilstName, char* indilst2remove, char* sex_file, double grm_cutoff, bool erm_cutoff_2sides, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool dont_read_N)
    {
        int i = 0, j = 0;
        
        vector<string> grm_id;
        if (merge_grm_flag) merge_grm(einfo, grm_file);
        else read_grm(einfo,grm_file, grm_id, true, false, dont_read_N,true);
        
        if (indilstName!=NULL) keep_indi(einfo, indilstName);
        if (indilst2remove!=NULL) remove_indi(einfo, indilst2remove);
        if (grm_cutoff>-1.0) rm_cor_indi(einfo, grm_cutoff, erm_cutoff_2sides);
        if (sex_file!=NULL) update_sex(einfo,sex_file);
        if (adj_grm_fac>-1.0) adj_grm(einfo, adj_grm_fac);
        if (dosage_compen>-1) dc(einfo, dosage_compen);
        if (grm_cutoff>-1.0 || indilstName!=NULL || indilst2remove!=NULL) {
            MatrixXd grm_buf(einfo->_grm);
            einfo->_grm.resize(einfo->_eii_include.size(), einfo->_eii_include.size());
            for (i = 0; i < einfo->_eii_include.size(); i++) {
                for (j = 0; j <= i; j++) einfo->_grm(i, j) = grm_buf(einfo->_eii_include[i], einfo->_eii_include[j]);
            }
            grm_buf.resize(0,0);
            MatrixXf grm_N_buf = einfo->_grm_N;
            einfo->_grm_N.resize(einfo->_eii_include.size(), einfo->_eii_include.size());
            for (i = 0; i < einfo->_eii_include.size(); i++) {
                for (j = 0; j <= i; j++) einfo->_grm_N(i, j) = grm_N_buf(einfo->_eii_include[i], einfo->_eii_include[j]);
            }
        }
    }
    void pca(char* outFileName, char* efileName, char* befileName, char* erm_file,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,bool transposed, int efileType,bool no_fid_flag,int valueType, double grm_cutoff, bool erm_cutoff_2sides, bool merge_grm_flag, int out_pc_num)
    {
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL && efileName==NULL && erm_file==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile or input ORM by --orm.\n");
            TERMINATE();
        }
        if(erm_file!=NULL)
        {
        manipulate_orm(&einfo, erm_file, indilstName, indilst2remove, NULL, grm_cutoff,erm_cutoff_2sides, -2.0, -2, merge_grm_flag, true);
        }
        else
        {
            if(efileName!=NULL)
            {
                if(transposed) read_efile_t(efileName,&einfo,efileType,no_fid_flag,valueType);
                else read_efile(efileName,&einfo,efileType,no_fid_flag,valueType);
                epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
                eii_man(&einfo,indilstName,indilst2remove);
                if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
                if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
            }else{
                char inputname[FNAMESIZE];
                memcpy(inputname,befileName,strlen(befileName)+1);
                char* suffix=inputname+strlen(befileName);
                memcpy(suffix,".oii",5);
                read_eii(inputname,&einfo);
                eii_man(&einfo,indilstName,indilst2remove);
                memcpy(suffix,".opi",5);
                read_epi(inputname,&einfo);
                epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
                if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
                memcpy(suffix,".bod",5);
                read_beed(inputname,&einfo);
                if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
                if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
            }
            if(einfo._eType == METHYLATION)
            {
                if(beta2m && m2beta){
                    //no chance to enter here
                    LOGPRINTF("Error: --m2beta should not be with --beta2m.\n");
                    TERMINATE();
                }
                if(beta2m && einfo._valType==BETAVALUE) beta_2_m(&einfo);
                if(m2beta && einfo._valType==MVALUE) m_2_beta(&einfo);
            }
            make_erm(&einfo, erm_alg, erm_bin_flag,  outFileName, false);
        }
        einfo._grm_N.resize(0, 0);
        int i = 0, j = 0, n = (int)einfo._eii_include.size();
        LOGPRINTF("\nPerforming principal component analysis ...\n");
        
        SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
        MatrixXd evec = (eigensolver.eigenvectors());
        VectorXd eval = eigensolver.eigenvalues();
        
        string eval_file = string(outFileName) + ".eigenval";
        ofstream o_eval(eval_file.c_str());
        if (!o_eval) {
            LOGPRINTF("Error: can not open the file %s to read.",eval_file.c_str());
            TERMINATE();
        }
        for (i = n - 1; i >= 0; i--) o_eval << eval(i) << endl;
        o_eval.close();
        LOGPRINTF("%d Eigenvalues of  individuals have been saved in %s.\n",n,eval_file.c_str());
        string evec_file = string(outFileName) + ".eigenvec";
        ofstream o_evec(evec_file.c_str());
        if (!o_evec) {
            LOGPRINTF ("Error: can not open the file %s to read.\n",evec_file.c_str());
            TERMINATE();
        }
        if (out_pc_num > n) out_pc_num = n;
        for (i = 0; i < n; i++) {
            o_evec << einfo._eii_fid[einfo._eii_include[i]] << " " << einfo._eii_iid[einfo._eii_include[i]];
            for (j = n - 1; j >= (n - out_pc_num); j--) o_evec << " " << evec(i, j);
            o_evec << endl;
        }
        o_evec.close();
        LOGPRINTF("The first %d eigenvectors of %d individuals have been saved in %s.\n" ,out_pc_num,n,evec_file.c_str());
    }
    void  testQAssoc(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, bool est_eff_n)
    {
        FILE* assoc=NULL;
        long write_count=0;
        string outstr="";
        string assocfile="";
        if(outfileName!=NULL){
            assocfile = string(outfileName)+".qassoc";
            assoc = fopen(assocfile.c_str(), "w");
            if (!(assoc)) {
                LOGPRINTF("ERROR: open error %s\n", assocfile.c_str());
                TERMINATE();
            }
            if (est_eff_n) outstr="probeChr\tProbeID\tProbe_bp\tBETA\tSE\tR2\tT\tP\tNMISS\tEFFN\n";
            else outstr="probeChr\tProbeID\tProbe_bp\tBETA\tSE\tR2\tT\tP\tNMISS\n";
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
        } else {
            assoc_rlts.clear();
        }

        //LOGPRINTF("Performing association analysis by linear regression...\n");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            //printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            //fflush(stdout);
            double mu_y=0.0, mu_x=0.0;
            double nonmiss=0.0;
            int chr=einfo->_epi_chr[einfo->_epi_include[i]];
            string prbid=einfo->_epi_prb[einfo->_epi_include[i]];
            string gene=einfo->_epi_gene[einfo->_epi_include[i]];
            int BP=einfo->_epi_bp[einfo->_epi_include[i]];
            char oren=einfo->_epi_orien[einfo->_epi_include[i]];
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double phval=einfo->_eii_pheno[einfo->_eii_include[j]];
                if(abs(phval+MISSING_PHENO)>1e-8)
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9){
                        mu_y+=phval;
                        mu_x+=val;
                        nonmiss+=1.0;
                    }
                }
            }
            mu_y/=nonmiss;
            mu_x/=nonmiss;
            double xty=0.0, xtx=0.0, yty=0.0;
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double phval=einfo->_eii_pheno[einfo->_eii_include[j]];
                if(abs(phval+MISSING_PHENO)>1e-8)
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9){
                        yty+=(phval-mu_y)*(phval-mu_y);
                        xty+=(val-mu_x)*(phval-mu_y);
                        xtx+=(val-mu_x)*(val-mu_x);
                    }
                }
            }
            xty/=nonmiss-1;
            xtx/=nonmiss-1;
            yty/=nonmiss-1;
            
            // Test statistics
            double beta = xty / xtx;
            double vbeta = ( yty/xtx - (xty*xty)/(xtx*xtx) ) / (nonmiss-2);
            double t = beta / sqrt(vbeta);
            double t_p=pT(t,nonmiss-2);
            double r2 = (xty * xty ) / ( yty * xtx );
            double estn=yty/(xtx*vbeta);
           
            ASSOCRLT currlt;
            if(assoc) {
                if (est_eff_n) outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(beta) + '\t' + atos(sqrt(vbeta)) + '\t' + atos(r2) + '\t' + atos(t) + '\t' + dtos(t_p) + '\t' + atos(nonmiss) +'\t'+ atos(estn) +'\n';
                else  outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(beta) + '\t' + atos(sqrt(vbeta)) + '\t' + atos(r2) + '\t' + atos(t) + '\t' + dtos(t_p) + '\t' + atos(nonmiss) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            } else {
                currlt.BETA=beta;
                currlt.SE=sqrt(vbeta);
                currlt.T=t;
                currlt.R2=r2;
                currlt.PVAL=t_p;
                currlt.CHR=chr;
                strcpy2(&currlt.GENE, gene);
                strcpy2(&currlt.PROBE, prbid);
                currlt.BP=BP;
                currlt.OREN=oren;
                currlt.NMISS=nonmiss;
                assoc_rlts.push_back(currlt);
            }
        }
        if(assoc){
            //LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            //LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
        }


    }
    void  assoc(char* outfileName, char* befileName, char* problstName, char* problst2exclde, char* genelistName,  int chr,char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename, char* probe2exclde, char* indilstName, char* indilst2remove,char* phenofileName,char* mpheno,  char* covfileName, char* qcovfileName, double std_thresh ,double upperBeta,double lowerBeta, bool est_eff_n )
    {
        eInfo einfo;
        init_einfo(&einfo);
        
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("Error: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
         if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        
        vector<ASSOCRLT> assoc_rlts;
        testQAssoc(assoc_rlts,outfileName, &einfo,est_eff_n);
        if(outfileName==NULL) free_assoclist(assoc_rlts);
       
    }
    void  stdprobe(eInfo* einfo)
    {
        if(loud) {LOGPRINTF("\nStandardizing probes...\n");}
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            //printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            //fflush(stdout);
            
            double nonmiss=0.0, mu=0.0,sd=0.0;
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            if(nonmiss>0)
            {
                mu/=nonmiss;
                for(int j=0; j<einfo->_eii_include.size(); j++)
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9)
                        einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]]=val-mu;
                }
            }
            
            if(nonmiss>1)
            {
                for(int j=0; j<einfo->_eii_include.size(); j++)
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9) sd+=val*val;
                }
                sd=sqrt(sd/(nonmiss-1.0));
                if(sd>1e-30)
                {
                    for(int j=0; j<einfo->_eii_include.size(); j++)
                    {
                        double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                        if(val<1e9) einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]]=val/sd;
                    }
                }
                
            }
            
        }
        if(loud) {LOGPRINTF("%ld probes have been standardized.\n",einfo->_epi_include.size());}
    }
    void rintprobe(eInfo* einfo)
    {
        LOGPRINTF("\nRank-Based Inverse Normal Transformation for each probe...\n");
        double c=0.5;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            fflush(stdout);
            
            vector<double> xvec, rxvec;
            vector<int> NMISS;
            double nonmiss=0.0;
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    NMISS.push_back(einfo->_eii_include[j]);
                    xvec.push_back(val);
                    nonmiss+=1.0;
                }
            }
            getRank2R(xvec,rxvec);
            #pragma omp parallel for
            for(int j=0;j<rxvec.size();j++)
                xvec[j]=qnorm((rxvec[j]-c)/(nonmiss-2*c+1));
            for(int j=0;j<xvec.size();j++) einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+NMISS[j]]=xvec[j];
        }
        LOGPRINTF("%ld probes have been normalised.\n",einfo->_epi_include.size());
    }
    
    void  addEnvEff(eInfo* einfo, char* prblstfname,char* efffname)
    {
        if(efffname !=NULL)
        {
            read_phen(einfo, efffname,NULL,false);
            if(prblstfname==NULL)
            {
                LOGPRINTF("\nAdjusting all the probes by adding environmental effect from file %s ...\n",efffname);
                for(int i=0;i<einfo->_epi_include.size();i++)
                {
                    for(int j=0; j<einfo->_eii_include.size(); j++)
                    {
                        double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                        if(val<1e9) einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]]=val+einfo->_eii_pheno[einfo->_eii_include[j]];
                    }
                }
                 LOGPRINTF("%ld probes have been adjusted.\n",einfo->_epi_include.size());
            }
            else
            {
                LOGPRINTF("\nAdjusting the specified probes from file %s by adding environmental effect from file %s ...\n",prblstfname,efffname);
                vector<string> problist;
                string msg="probes";
                read_msglist(prblstfname, problist,msg);
                vector<int> keep;
                map<string, int>::iterator iter;
                for(int i=0; i<problist.size(); i++) {
                    iter=einfo->_epi_map.find(problist[i]);
                    if(iter != einfo->_epi_map.end())
                    {
                         keep.push_back(iter->second);
                    }
                }
                stable_sort(keep.begin(), keep.end());
                LOGPRINTF("%ld probes are included.\n",keep.size());
                for(int i=0;i<keep.size();i++)
                {
                    for(int j=0; j<einfo->_eii_include.size(); j++)
                    {
                        double val=einfo->_val[keep[i]*einfo->_eii_num+einfo->_eii_include[j]];
                        if(val<1e9) einfo->_val[keep[i]*einfo->_eii_num+einfo->_eii_include[j]]=val+einfo->_eii_pheno[einfo->_eii_include[j]];
                    }
                }
                LOGPRINTF("%ld probes have been adjusted.\n",keep.size());
            }
        }
        else
        {
            LOGPRINTF("\nAdjusting probes by adding normal distributed environmental effect with the mean as 0 and the sd as the one of the corresponding probe...\n");
            for(int i=0;i<einfo->_epi_include.size();i++)
            {
                printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
                fflush(stdout);
                
                int Seed = -rand_seed();
                default_random_engine generator(Seed);
                
                vector<double> xvec;
                vector<int> NMISS;
                double nonmiss=0.0;
                for(int j=0; j<einfo->_eii_include.size(); j++)
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9){
                        NMISS.push_back(einfo->_eii_include[j]);
                        xvec.push_back(val);
                        nonmiss+=1.0;
                    }
                }
                double para=var(xvec);
                normal_distribution<double> distribution(0.0,sqrt(para));
                vector<double> enveff;
                for(int j=0;j<xvec.size();j++) enveff.push_back(distribution(generator));
                
                for(int j=0;j<xvec.size();j++) einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+NMISS[j]]=xvec[j]+enveff[j];
            }
            LOGPRINTF("%ld probes have been adjusted.\n",einfo->_epi_include.size());

        }
        
    }

    void write_assoc_rlt(vector<ASSOCRLT> &assoc_rlts,char* outfileName)
    {
        if(outfileName==NULL)
        {
            LOGPRINTF("ERROR: the output file name is null.\n");
            TERMINATE();
        }
        FILE* assoc=NULL;
        long write_count=0;
        string outstr="";
        string assocfile = string(outfileName)+".linear";
        assoc = fopen(assocfile.c_str(), "w");
        if (!(assoc)) {
            LOGPRINTF("ERROR: open error %s\n", assocfile.c_str());
            TERMINATE();
        }
        outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\tNMISS\n";
        if(fputs_checked(outstr.c_str(),assoc))
        {
            LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
            TERMINATE();
        }
        for(int i=0;i<assoc_rlts.size();i++)
        {
            string outstr = atos(assoc_rlts[i].CHR) + '\t' + assoc_rlts[i].PROBE + '\t' + atos(assoc_rlts[i].BP) + '\t' + atos(assoc_rlts[i].GENE) + '\t' + atos(assoc_rlts[i].OREN) + '\t' + atos(assoc_rlts[i].BETA) + '\t' + atos(assoc_rlts[i].SE) + '\t' + dtos(assoc_rlts[i].PVAL) + '\t' + atos(assoc_rlts[i].NMISS) +'\n';
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
            write_count++;
        }
        LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
        fclose(assoc);
    }
    double get_lambda(vector<ASSOCRLT> &assoc_rlts)
    {
        vector<double> chis;
        chis.resize(assoc_rlts.size());
        for(int i=0;i<assoc_rlts.size();i++)
        {
            double z=assoc_rlts[i].BETA/assoc_rlts[i].SE;
            chis[i]=z*z;
        }
        double lambda=median(chis)/0.455;
        return lambda;
    }
    void  testLinear_fast(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, MatrixXd &COV_plus, bool reverse, bool fdrflag)
    {
        LOGPRINTF("NOTE: the missing value will be imputed by the mean of each probe. This might lead to some differences in result if there are a substantial number of missing values in the data.\n");
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
            outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\tNMISS\n";
            if(fdrflag)  outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\tq\tNMISS\n";
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
        _X_c=construct_X(einfo, E_float, qE_float,_X);
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
        vector<double> pval, qval;
        LOGPRINTF("\nAdjusting the phenotype...\n");
        /*****/
//        LOGPRINTF("Saving the _X ...\n");
//        string filename=string(outfileName)+"._X.mat";
//        FILE* tmpfile=fopen(filename.c_str(),"w");
//        if(!tmpfile)
//        {
//            LOGPRINTF("error open file.\n");
//            TERMINATE();
//        }
//        for(int t=0;t<_X.rows();t++)
//        {
//            string str="";
//            for(int k=0;k<_X.cols();k++)
//            {
//                str +=atos(_X(t,k)) + '\t';
//            }
//            str += '\n';
//            fputs(str.c_str(),tmpfile);
//        }
//        fclose(tmpfile);
//        LOGPRINTF("the _X matrix is saved in the file %s.\n",filename.c_str());
        /*****/
        MatrixXd XtX_i;
        XtX_i=_X.transpose()*_X;
        /*****/
//        LOGPRINTF("Saving the XtX ...\n");
//        filename=string(outfileName)+"._XtX.mat";
//        tmpfile=fopen(filename.c_str(),"w");
//        if(!tmpfile)
//        {
//            LOGPRINTF("error open file.\n");
//            TERMINATE();
//        }
//        for(int t=0;t<XtX_i.rows();t++)
//        {
//            string str="";
//            for(int k=0;k<XtX_i.cols();k++)
//            {
//                str +=atos(XtX_i(t,k)) + '\t';
//            }
//            str += '\n';
//            fputs(str.c_str(),tmpfile);
//        }
//        fclose(tmpfile);
//        LOGPRINTF("the XtX matrix is saved in the file %s.\n",filename.c_str());
        /*****/
        bool determinant_zero=false;
        inverse_V(XtX_i, determinant_zero);
        /*****/
//        LOGPRINTF("Saving the XtXi ...\n");
//        filename=string(outfileName)+"._XtXi.mat";
//        tmpfile=fopen(filename.c_str(),"w");
//        if(!tmpfile)
//        {
//            LOGPRINTF("error open file.\n");
//            TERMINATE();
//        }
//        for(int t=0;t<XtX_i.rows();t++)
//        {
//            string str="";
//            for(int k=0;k<XtX_i.cols();k++)
//            {
//                str +=atos(XtX_i(t,k)) + '\t';
//            }
//            str += '\n';
//            fputs(str.c_str(),tmpfile);
//        }
//        fclose(tmpfile);
//        LOGPRINTF("the XtXi matrix is saved in the file %s.\n",filename.c_str());
        /*****/
        if(determinant_zero)
        {
            LOGPRINTF("ERROR: The matrix is not invertible. please check the correlation of the covariates.\n");
            TERMINATE();
        }
        MatrixXd XtXiXt=XtX_i*_X.transpose();
        
        VectorXd y(einfo->_eii_include.size());
        #pragma omp parallel for
        for(int j=0; j<einfo->_eii_include.size(); j++)
        {
            y(j)=einfo->_eii_pheno[einfo->_eii_include[j]]; // no phval missing here, cos individuals with missing phenotype were removed when read phenotype file.
        }
        VectorXd b_hat=XtXiXt*y;
        VectorXd yresi=(y-_X*b_hat);// if no covariate yresi=y-mean(y)
        
        LOGPRINTF("\nPerforming assocaiton analysis...\n");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            fflush(stdout);
            
            int chr=einfo->_epi_chr[einfo->_epi_include[i]];
            string prbid=einfo->_epi_prb[einfo->_epi_include[i]];
            string gene=einfo->_epi_gene[einfo->_epi_include[i]];
            int BP=einfo->_epi_bp[einfo->_epi_include[i]];
            char oren=einfo->_epi_orien[einfo->_epi_include[i]];
            
            vector<int> MISS;
            VectorXd x(einfo->_eii_include.size());
            double mu= 0.0, sum=0.0;
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9)
                {
                    x(j)=val;
                    sum+=val;
                } else {
                    MISS.push_back(j);
                }
            }
            double ratio=1.0*MISS.size()/einfo->_eii_include.size();
            if(ratio>0.1)
            {
                printf("WARNING: %f %% of the values in probe %s are missing.\n",ratio*100, einfo->_epi_prb[einfo->_epi_include[i]].c_str());
            }
            if(MISS.size()>0)
            {
                if(einfo->_eii_include.size() > MISS.size()) mu=sum/(einfo->_eii_include.size()-MISS.size());
                for(int k=0;k<MISS.size();k++)
                {
                    x(MISS[k])=mu;
                }
            }
            
            if(x.size()!=_X.rows() || x.size()<1) {
                LOGPRINTF("Error: The row number of C and the length of y not match.\n");
                TERMINATE();
            }
            VectorXd b_hat=XtXiXt*x;
            VectorXd residual=(x-_X*b_hat);
            vector<double> rst;
            long nmiss=einfo->_eii_include.size()-MISS.size();
            if(reverse) adjusted_reg(residual, yresi, rst,_X_c-1);
            else adjusted_reg(yresi, residual, rst,_X_c-1);
            
            if(outfileName!=NULL && !fdrflag)
            {
                outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nmiss) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            }
            
          ASSOCRLT currlt;
            currlt.BETA=rst[0];
            currlt.SE=rst[1];
            currlt.T=-9;
            currlt.R2=-9;
            currlt.PVAL=rst[2];
            if(fdrflag) pval.push_back(rst[2]);
            currlt.CHR=chr;
            strcpy2(&currlt.GENE, gene);
            strcpy2(&currlt.PROBE, prbid);
            currlt.BP=BP;
            currlt.OREN=oren;
            currlt.NMISS=nmiss;
            assoc_rlts.push_back(currlt);
        }
        if(assoc && fdrflag) {
            getQval(pval,qval);
            for(int i=0;i<assoc_rlts.size();i++)
            {
                string outstr = atos(assoc_rlts[i].CHR) + '\t' + assoc_rlts[i].PROBE + '\t' + atos(assoc_rlts[i].BP) + '\t' + atos(assoc_rlts[i].GENE) + '\t' + atos(assoc_rlts[i].OREN) + '\t' + atos(assoc_rlts[i].BETA) + '\t' + atos(assoc_rlts[i].SE) + '\t' + dtos(assoc_rlts[i].PVAL) + '\t' + dtos(qval[i]) + '\t' + atos(assoc_rlts[i].NMISS) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            }
        }
        if(assoc){
            LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
        }
        
        
    }
    void  testLinear(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, MatrixXd &COV_plus, bool reverse, bool fdrflag)
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
            outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\tNMISS\n";
            if(fdrflag)  outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\tq\tNMISS\n";
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
        _X_c=construct_X(einfo, E_float, qE_float,_X);
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
        vector<double> pval,qval;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            double desti=1.0*i/(einfo->_epi_include.size()-1);
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }

            vector<double> yvec, cvec, xvec;
            MatrixXd X=_X;
            double nonmiss=0.0;
            int miss=0;
            int chr=einfo->_epi_chr[einfo->_epi_include[i]];
            string prbid=einfo->_epi_prb[einfo->_epi_include[i]];
            string gene=einfo->_epi_gene[einfo->_epi_include[i]];
            int BP=einfo->_epi_bp[einfo->_epi_include[i]];
            char oren=einfo->_epi_orien[einfo->_epi_include[i]];
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double phval=einfo->_eii_pheno[einfo->_eii_include[j]];
                if(abs(phval+MISSING_PHENO)>1e-8) // no phval missing here, cos individuals with missing phenotype were removed when read phenotype file.
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9){
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
            bool notInvertible= false;
            if(reverse) notInvertible=lin(x, X, y, rst);
            else notInvertible=lin(y, X, x, rst);
            if(notInvertible) {
                LOGPRINTF("ERROR: The matrix for probe %s is not invertible. Maybe there are multicollinearity in the covariates.\n",prbid.c_str());
                TERMINATE();
            }
            ASSOCRLT currlt;
            if(assoc && !fdrflag) {
                string outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
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
            if(fdrflag) pval.push_back(rst[2]);
            currlt.CHR=chr;
            strcpy2(&currlt.GENE, gene);
            strcpy2(&currlt.PROBE, prbid);
            currlt.BP=BP;
            currlt.OREN=oren;
            currlt.NMISS=nonmiss;
            assoc_rlts.push_back(currlt);
        }
        if(assoc && fdrflag) {
            getQval(pval,qval);
            for(int i=0;i<assoc_rlts.size();i++)
            {
                string outstr = atos(assoc_rlts[i].CHR) + '\t' + assoc_rlts[i].PROBE + '\t' + atos(assoc_rlts[i].BP) + '\t' + atos(assoc_rlts[i].GENE) + '\t' + atos(assoc_rlts[i].OREN) + '\t' + atos(assoc_rlts[i].BETA) + '\t' + atos(assoc_rlts[i].SE) + '\t' + dtos(assoc_rlts[i].PVAL) + '\t' + dtos(qval[i]) + '\t' + atos(assoc_rlts[i].NMISS) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            }
        }
        if(assoc){
            LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
        }
    }

    void  linear(char* outfileName, char* befileName, char* problstName, char* problst2exclde, char* genelistName,  int chr,char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename, char* probe2exclde, char* indilstName, char* indilst2remove,char* phenofileName,char* mpheno,  char* covfileName, char* qcovfileName, double std_thresh ,double upperBeta,double lowerBeta ,int tsk_ttl,int tsk_id,bool fastlinear, bool stdprb, bool reverse, bool fdrflag)
    {
        eInfo einfo;
        init_einfo(&einfo);
      
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("Error: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        char outputname[FNAMESIZE];
        outputname[0]='\0';
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(tsk_ttl>1) {
            extract_probe(&einfo,  tsk_ttl,  tsk_id);
            if(outfileName!=NULL) {
                string tmp=  string(outfileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outfileName=outputname;
            }
        }
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
        if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        if(stdprb) stdprobe(&einfo);
        
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        if(fastlinear) testLinear_fast(assoc_rlts,outfileName, &einfo,COV_plus, reverse, fdrflag);
        else testLinear(assoc_rlts,outfileName, &einfo,COV_plus,reverse,fdrflag);
    }
    void  testLogit(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo, MatrixXd &COV_plus)
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
            outstr="Chr\tProbe\tbp\tGene\tOrientation\tOR\tse\tp\tNMISS\n";
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
        }
        
        vector<double> bvec(einfo->_epi_include.size()), sevec(einfo->_epi_include.size()), nmvec(einfo->_epi_include.size());
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(einfo, E_float, qE_float,_X);
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
        double cr=0.0;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            double desti=1.0*i/(einfo->_epi_include.size()-1);
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            
            vector<int> yvec;
            vector<double> cvec;
            vector<double> xvec;
            MatrixXd X=_X;
            double nonmiss=0.0;
            int miss=0;
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                int phval=einfo->_eii_pheno[einfo->_eii_include[j]];
                if(abs(phval+MISSING_PHENO)>1e-8) // no phval missing here, cos individuals with missing phenotype were removed when read phenotype file.
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9){
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
            for(int i=0;i<einfo->_epi_include.size();i++)
            {
                int idx=einfo->_epi_include[i];
                double z=bvec[i]/sevec[i];
                double pval=pchisq(z*z, 1);
                string outstr = atos(einfo->_epi_chr[idx]) + '\t' + einfo->_epi_prb[idx] + '\t' + atos(einfo->_epi_bp[idx]) + '\t' + atos(einfo->_epi_gene[idx]) + '\t' + atos(einfo->_epi_orien[idx]) + '\t' + (sevec[i]<0?"NA":atos(exp(bvec[i]))) + '\t' + (sevec[i]<0?"NA":atos(sevec[i])) + '\t' + (sevec[i]<0?"NA":dtos(pval)) + '\t' + atos(nmvec[i]) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            }
        }
        assoc_rlts.clear();
        assoc_rlts.resize(einfo->_epi_include.size());
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int chr=einfo->_epi_chr[einfo->_epi_include[i]];
            string prbid=einfo->_epi_prb[einfo->_epi_include[i]];
            string gene=einfo->_epi_gene[einfo->_epi_include[i]];
            int BP=einfo->_epi_bp[einfo->_epi_include[i]];
            char oren=einfo->_epi_orien[einfo->_epi_include[i]];
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
            assoc_rlts[i].OREN=oren;
            strcpy2(&assoc_rlts[i].GENE, gene);
            strcpy2(&assoc_rlts[i].PROBE, prbid);
        }
        if(assoc)
        {
            LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
        }
    }
    void  logistic(char* outfileName, char* befileName, char* problstName, char* problst2exclde, char* genelistName,  int chr,char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename, char* probe2exclde, char* indilstName, char* indilst2remove,char* phenofileName,char* mpheno,  char* covfileName, char* qcovfileName, double std_thresh ,double upperBeta,double lowerBeta ,int tsk_ttl,int tsk_id)
    {
        eInfo einfo;
        init_einfo(&einfo);
        
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("Error: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        char outputname[FNAMESIZE];
        outputname[0]='\0';
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(tsk_ttl>1) {
            extract_probe(&einfo,  tsk_ttl,  tsk_id);
            if(outfileName!=NULL) {
                string tmp=  string(outfileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id);
                strcpy(outputname,tmp.c_str());
                outfileName=outputname;
            }
        }
        if(phenofileName !=NULL) read_cc(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
        if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        testLogit(assoc_rlts,outfileName, &einfo,COV_plus);
        
    }
    long read_QTL_file(eInfo* einfo,string qtl_file, vector<string> &qtl_name, vector<double> &qtl_eff, vector<int> &have_eff)
    {
        qtl_name.clear();
        qtl_eff.clear();
        have_eff.clear();
        
        ifstream i_qtl(qtl_file.c_str());
        if (!i_qtl) throw ("Error: can not open the file [" + qtl_file + "] to read.");
        string qtl_buf, str_buf;
        double qtl_eff_buf = 0.0;
        cout << "Reading a list of probs (as causal variants) from [" + qtl_file + "]." << endl;
        map<string, int>::iterator iter, End = einfo->_epi_map.end();
        vector<string> vs_buf;
        vector<int> confirm(einfo->_epi_num);
        while (i_qtl) {
            i_qtl >> qtl_buf;
            if (i_qtl.eof()) break;
            iter = einfo->_epi_map.find(qtl_buf);
            if (getline(i_qtl, str_buf) && split_string(str_buf, vs_buf, " \t\n") > 0) {
                have_eff.push_back(1);
                qtl_eff_buf = atof(vs_buf[0].c_str());
                if (fabs(qtl_eff_buf) > 1e5) throw ("Error: invalid effect size specified for the causal variant [" + str_buf + "].");
            } else {
                have_eff.push_back(0);
                qtl_eff_buf = 0.0;
            }
            if (iter != End) {
                qtl_name.push_back(qtl_buf);
                qtl_eff.push_back(qtl_eff_buf);
            }
        }
        vector<string> qtl_name_buf=qtl_name;
        stable_sort(qtl_name_buf.begin(), qtl_name_buf.end());
        qtl_name_buf.erase(unique(qtl_name_buf.begin(), qtl_name_buf.end()), qtl_name_buf.end());
        i_qtl.close();
        
        if(qtl_name_buf.size() < qtl_name.size()) throw("Error: there are duplicated probe IDs.");
        cout << qtl_name.size() << " Probes (as causal variants) to be included from [" + qtl_file + "]." << endl;
        return (qtl_name.size());
    }
    
    
    
    void EWAS_simu(char* outFileName, char* befileName, int simu_num, char* sigCpG_file, int case_num, int control_num, double hsq, double K, int seed, int eff_mod, bool simu_residual)
    {
       
        bool cc_flag = false;
        if (case_num > 0 || control_num > 0) cc_flag = true;
        
        cout << "Simulation parameters:" << endl;
        cout << "Number of simulation replicate(s) = " << simu_num << " (Default = 1)" << endl;
        cout << "R-squared " << (cc_flag ? "of liability = " : " = ") << hsq << " (Default = 0.1)" << endl;
        if (cc_flag) {
            cout << "Disease prevalence = " << K << " (Default = 0.1)" << endl;
            cout << "Number of cases = " << case_num << endl;
            cout << "Number of controls = " << control_num << endl;
        }
        cout << endl;
        
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        // Read sig CpG sites
   
        vector<string> qtl_name;
        vector<int> have_eff;
        vector<double> qtl_eff;
        long qtl_num = read_QTL_file(&einfo,sigCpG_file, qtl_name, qtl_eff, have_eff);
        update_map_kp(qtl_name, einfo._epi_map, einfo._epi_include);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
       
        // Generate QTL effects
        int Seed = -rand_seed();
        if (hsq > 0.0) {
            int num_gener_qtl_eff = 0;
            for (int i = 0; i < qtl_num; i++) {
                if (have_eff[i] == 0) {
                    qtl_eff[i] = gasdev(Seed);
                    num_gener_qtl_eff++;
                }
            }
            if (qtl_num - num_gener_qtl_eff > 0) cout << qtl_num - num_gener_qtl_eff << " user-specified QTL effects." << endl;
            if (num_gener_qtl_eff > 0) cout << num_gener_qtl_eff << " unspecified QTL effects are generated from standard normal distribution." << endl;
            
            vector<string> vs_buf(qtl_num);
            for (int i = 0; i < qtl_num; i++) vs_buf[i] = einfo._epi_prb[einfo._epi_include[i]];
            vector<int> indx;
            match(vs_buf, qtl_name, indx);
            qtl_name = vs_buf;
            vector<double> qtl_eff_buf(qtl_eff);
            for (int i = 0; i < qtl_num; i++){
                qtl_eff[i] = qtl_eff_buf[indx[i]];
            }
           //now qtl_name,qtl_eff have the same order as einfo._epi_prb
        } else {
            qtl_eff.clear();
            qtl_eff.resize(qtl_num, 0.0);
        }
     
        
    
        vector< vector<bool> > X_bool;
        MatrixXd X;
        if(eff_mod == 0){ //standardise
            LOGPRINTF("Simulating phenotype with standardised probes...\n");
            LOGPRINTF("Standardising probes...\n");
            std_probe(&einfo, X_bool, true, X);
        } else {
            LOGPRINTF("Simulating phenotype with raw probes...\n");
            uint64_t n = einfo._eii_include.size(), m = einfo._epi_include.size();
            if(n!=einfo._eii_num)
            {
                LOGPRINTF("Inconsistent individual number when standardising probes. Please report this bug!\n");
                TERMINATE();
            }
            Map<MatrixXd> X_map(&einfo._val[0],n,m);  // n should equal with einfo->_eii_num. otherwise it is wrong
            X.resize(n,m);
            bool missingwarnflag=false;
            #pragma omp parallel for
            for(int j=0; j<m; j++){
                double mu=0.0, nonmiss=0.0;
                for(int i=0; i<n; i++)
                {
                    double val=X_map(einfo._eii_include[i], einfo._epi_include[j]);
                    X(i,j)=val;
                    if(val<1e9) {
                        mu+=val;
                        nonmiss+=1.0;
                    }
                }
                if(nonmiss>0) mu/=nonmiss;
                if(nonmiss<n) {
                    if(!missingwarnflag) LOGPRINTF("WARNING: Missing value found.\n");
                    for(int i=0; i<n; i++) if(X(i,j)>=1e9){
                             X(i,j)=mu; //missing value with mean
                             LOGPRINTF("probe %s of individual %s was inputed with mean %f.\n",einfo._epi_prb[einfo._epi_include[j]].c_str(),einfo._eii_iid[einfo._eii_include[i]].c_str(),mu);
                    }
                }
            }

        }
        
        // Calculate Ve and threhold
        double var_g = 0.0, var_e = 1.0;
        vector<double> g(einfo._eii_include.size());
        if (hsq > 0.0) {
            for (int i = 0; i < einfo._eii_include.size(); i++) {
                for (int j = 0; j < qtl_num; j++) g[i] += X(i,j) * qtl_eff[j];  //sigma(xb)
            }
            var_g = var(g);
            var_e = var_g * (1.0 / hsq - 1.0);
        }
        double sd_e = sqrt(var_e);
        // Output par file
        output_simu_par(outFileName, qtl_name, qtl_eff);
        
        // Output phenotype file
        cout << "Simulating EWAS based on the real profile data with " << simu_num << " replicate(s) ..." << endl;
        vector< vector<double> > y(simu_num);
        int case_num_buf = 0, control_num_buf = 0;
        for (int i = 0; i < simu_num; i++) {
            y[i].resize(einfo._eii_include.size());
            for (int j = 0; j < einfo._eii_include.size(); j++) {
                if(!simu_residual){
                    if (hsq < 1.0) y[i][j] = g[j] + sd_e * gasdev(Seed);
                    else y[i][j] = g[j];
                } else {
                    if (hsq < 1.0) y[i][j] =  sd_e * gasdev(Seed);
                    else y[i][j] = 0;
                }
                
            }
            if (cc_flag) {
                case_num_buf = 0;
                control_num_buf = 0;
                vector<double> y_buf(y[i]);
                stable_sort(y_buf.begin(), y_buf.end());
                int n = (int) (einfo._eii_include.size() * (1.0 - K));
                double Th = 0.5 * (y_buf[n] + y_buf[n - 1]);
                for (int j = 0; j < einfo._eii_include.size(); j++) {
                    if (y[i][j] > Th) {
                        if (case_num_buf < case_num) {
                            y[i][j] = 2;
                            case_num_buf++;
                        } else y[i][j] = -9;
                    } else {
                        if (control_num_buf < control_num) {
                            y[i][j] = 1;
                            control_num_buf++;
                        } else y[i][j] = -9;
                    }
                }
            }
        }
        
    
            save_phenfile(&einfo,outFileName, y);
            if (cc_flag) cout << "Simulated " << case_num_buf << " cases and " << control_num << " controls have been saved in [" + string(outFileName) + ".phen" + "]." << endl;
            else cout << "Simulated phenotypes of " << einfo._eii_include.size() << " individuals have been saved in [" + string(outFileName) + ".phen" + "]." << endl;
            
    }
    
    void EWAS_simu2(char* outFileName, char* befileName, int simu_num, char* sigCpG_file, char* sigConfouder_file, int case_num, int control_num, double hsq1, double hsq2, double K, int seed, int eff_mod, bool simu_residual)
    {
        
        bool cc_flag = false;
        if (case_num > 0 || control_num > 0) cc_flag = true;
        if(hsq1+hsq2 > 1)
        {
            LOGPRINTF("Error: the sum of R2 > 1.\n");
            TERMINATE();
        }
        LOGPRINTF("Simulation parameters:\n");
        LOGPRINTF("Number of simulation replicate(s) = %d (Default = 1) \n", simu_num);
        LOGPRINTF("R-squear of the causal is %f (Default = 0.1) \n",hsq1);
        LOGPRINTF("R-squear of the confounder is %f (Default = 0.0) \n",hsq2);
        
        if (cc_flag) {
            cout << "Disease prevalence = " << K << " (Default = 0.1)" << endl;
            cout << "Number of cases = " << case_num << endl;
            cout << "Number of controls = " << control_num << endl;
        }
        cout << endl;
        
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        // Read sig CpG sites
        vector< vector<double> > g;
        vector< vector<string> > qtl_name;
        vector< vector<int> > have_eff;
        vector< vector<double> > qtl_eff;
        vector< MatrixXd > X;
        vector<long> qtl_num;
        vector<double> hsq;
        g.resize(2);
        X.resize(2);
        g[0].resize(einfo._eii_include.size());
        g[1].resize(einfo._eii_include.size());
        qtl_name.resize(2);
        have_eff.resize(2);
        qtl_eff.resize(2);
        qtl_num.resize(2);
        hsq.push_back(hsq1);
        hsq.push_back(hsq2);
        qtl_num[0] = read_QTL_file(&einfo,sigCpG_file, qtl_name[0], qtl_eff[0], have_eff[0]);
        qtl_num[1] = read_QTL_file(&einfo,sigConfouder_file, qtl_name[1], qtl_eff[1], have_eff[1]);
        
        vector<string> ttlnames;
        ttlnames.reserve( qtl_name[0].size() + qtl_name[1].size() );
        ttlnames.insert( ttlnames.end(), qtl_name[0].begin(), qtl_name[0].end() );
        ttlnames.insert( ttlnames.end(), qtl_name[1].begin(), qtl_name[1].end() );
        update_map_kp(ttlnames, einfo._epi_map, einfo._epi_include);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        
        vector<int> include_o(einfo._epi_include);
        map<string, int> snp_name_map_o(einfo._epi_map);
        int Seed = -rand_seed();
        for(int cidx=0;cidx<2;cidx++)
        {
            // Generate QTL effects
            update_map_kp(qtl_name[cidx], einfo._epi_map, einfo._epi_include);
            
            if (hsq[cidx] > 0.0) {
                int num_gener_qtl_eff = 0;
                for (int i = 0; i < qtl_num[cidx]; i++) {
                    if (have_eff[cidx][i] == 0) {
                        qtl_eff[cidx][i] = gasdev(Seed);
                        num_gener_qtl_eff++;
                    }
                }
                if (qtl_num[cidx] - num_gener_qtl_eff > 0) cout << qtl_num[cidx] - num_gener_qtl_eff << " user-specified QTL effects." << endl;
                if (num_gener_qtl_eff > 0) cout << num_gener_qtl_eff << " unspecified QTL effects are generated from standard normal distribution." << endl;
                
                vector<string> vs_buf(qtl_num[cidx]);
                for (int i = 0; i < qtl_num[cidx]; i++) vs_buf[i] = einfo._epi_prb[einfo._epi_include[i]];
                vector<int> indx;
                match(vs_buf, qtl_name[cidx], indx);
                qtl_name[cidx] = vs_buf;
                vector<double> qtl_eff_buf(qtl_eff[cidx]);
                for (int i = 0; i < qtl_num[cidx]; i++){
                    qtl_eff[cidx][i] = qtl_eff_buf[indx[i]];
                }
                //now qtl_name,qtl_eff have the same order as einfo._epi_prb
            } else {
                qtl_eff[cidx].clear();
                qtl_eff[cidx].resize(qtl_num[cidx], 0.0);
            }

            vector< vector<bool> > X_bool;
            if(eff_mod == 0){ //standardise
                LOGPRINTF("Simulating phenotype with standardised probes...\n");
                LOGPRINTF("Standardising probes...\n");
                std_probe(&einfo, X_bool, true, X[cidx]);
            }
            else
            {
                LOGPRINTF("Simulating phenotype with raw probes...\n");
                uint64_t n = einfo._eii_include.size(), m = einfo._epi_include.size();
                if(n!=einfo._eii_num)
                {
                    LOGPRINTF("Inconsistent individual number when standardising probes. Please report this bug!\n");
                    TERMINATE();
                }
                Map<MatrixXd> X_map(&einfo._val[0],n,m);  // n should equal with einfo->_eii_num. otherwise it is wrong
                X[cidx].resize(n,m);
                bool missingwarnflag=false;
#pragma omp parallel for
                for(int j=0; j<m; j++){
                    double mu=0.0, nonmiss=0.0;
                    for(int i=0; i<n; i++)
                    {
                        double val=X_map(einfo._eii_include[i], einfo._epi_include[j]);
                        X[cidx](i,j)=val;
                        if(val<1e9) {
                            mu+=val;
                            nonmiss+=1.0;
                        }
                    }
                    if(nonmiss>0) mu/=nonmiss;
                    if(nonmiss<n) {
                        if(!missingwarnflag) LOGPRINTF("WARNING: Missing value found.\n");
                        for(int i=0; i<n; i++) if(X[cidx](i,j)>=1e9){
                            X[cidx](i,j)=mu; //missing value with mean
                            LOGPRINTF("probe %s of individual %s was inputed with mean %f.\n",einfo._epi_prb[einfo._epi_include[j]].c_str(),einfo._eii_iid[einfo._eii_include[i]].c_str(),mu);
                        }
                    }
                }
            }
            if (hsq[cidx] > 0.0)
            {
                for (int i = 0; i < einfo._eii_include.size(); i++)
                {
                    for (int j = 0; j < qtl_num[cidx]; j++) g[cidx][i] += X[cidx](i,j) * qtl_eff[cidx][j];  //sigma(xb)
                }
            }
            
            einfo._epi_include=include_o;
            einfo._epi_map=snp_name_map_o;
        }
        //adj var_g
        double var_g0=var(g[0]);
        double var_g1=var(g[1]);
        double var_e=1, wsquared=1;
        if(hsq1>0 && hsq2>0)
        {
            double w=var_g0*hsq[1]/(var_g1*hsq[0]);
            wsquared=sqrt(w);
            var_e = var_g0*(1-hsq[0]-hsq[1])/hsq[0];
        }
        else if(hsq1>0)
        {
            var_e = var_g0*(1-hsq[0])/hsq[0];
        }
        else if(hsq2>0)
        {
            var_e = var_g1*(1-hsq[1])/hsq[1];
        }
        
        //vector<double> adjg;
        //for(int kk=0;kk<g[1].size();kk++)
        //    adjg.push_back(wsquared*g[1][kk]);
        //cout<<var(adjg)<<":"<<var_g0<<endl;
        
               
        double sd_e = sqrt(var_e);
        // Output par file
        
        output_simu_par(outFileName, qtl_name, qtl_eff);
        
        // Output phenotype file
        cout << "Simulating EWAS based on the real profile data with " << simu_num << " replicate(s) ..." << endl;
        vector< vector<double> > y(simu_num);
        int case_num_buf = 0, control_num_buf = 0;
        for (int i = 0; i < simu_num; i++)
        {
            y[i].resize(einfo._eii_include.size());
            for (int j = 0; j < einfo._eii_include.size(); j++)
            {
                if(!simu_residual){
                    y[i][j] = g[0][j] +wsquared*g[1][j] + sd_e * gasdev(Seed);
                    
                } else {
                    y[i][j] =  sd_e * gasdev(Seed);
                }
                
            }
            if (cc_flag)
            {
                case_num_buf = 0;
                control_num_buf = 0;
                vector<double> y_buf(y[i]);
                stable_sort(y_buf.begin(), y_buf.end());
                int n = (int) (einfo._eii_include.size() * (1.0 - K));
                double Th = 0.5 * (y_buf[n] + y_buf[n - 1]);
                for (int j = 0; j < einfo._eii_include.size(); j++) {
                    if (y[i][j] > Th) {
                        if (case_num_buf < case_num) {
                            y[i][j] = 2;
                            case_num_buf++;
                        } else y[i][j] = -9;
                    } else {
                        if (control_num_buf < control_num) {
                            y[i][j] = 1;
                            control_num_buf++;
                        } else y[i][j] = -9;
                    }
                }
            }
        }
        
        
        save_phenfile(&einfo,outFileName, y);
        if (cc_flag) cout << "Simulated " << case_num_buf << " cases and " << control_num << " controls have been saved in [" + string(outFileName) + ".phen" + "]." << endl;
        else cout << "Simulated phenotypes of " << einfo._eii_include.size() << " individuals have been saved in [" + string(outFileName) + ".phen" + "]." << endl;
        
    }
    
    void reverse_causal_simu(char* outFileName, char* befileName, char* phenofileName, char* mpheno, char* sigCpG_file, double rsq, int eff_mod)
    {
        
       if(mpheno!=NULL)
       {
           LOGPRINTF("Error: --mphen isn't supported.\n");
           TERMINATE();
       }
        cout << "Simulation parameters:" << endl;
        cout << "R-squared "  << rsq << " (Default = 0.1)" << endl;
        cout << endl;
        
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        read_phen(&einfo, phenofileName, mpheno,false);
        
        // Read sig CpG sites
        vector<string> qtl_name;
        vector<double> reverse_eff; vector<int> have_eff;
        read_QTL_file(&einfo,sigCpG_file, qtl_name, reverse_eff, have_eff);
        vector<int> include_o(einfo._epi_include);
        map<string, int> map_o(einfo._epi_map);
        update_map_kp(qtl_name, einfo._epi_map, einfo._epi_include);
        uint64_t n = einfo._eii_include.size(), m = einfo._epi_include.size();
        double var_y = var(einfo._eii_pheno);
        for(int i=0; i<m; i++)
        {
            vector<double> probe_data(n);
            for(int j=0; j<n; j++) probe_data[j]=einfo._val[einfo._epi_include[i]*einfo._eii_num+einfo._eii_include[j]];
            if(eff_mod == 0){
                standardise(probe_data, true);
            } else {
                // scan missing
                double mu=0.0, nonmiss=0;
                long n = probe_data.size();
                for(int i=0; i<n; i++){
                    if(probe_data[i]<1e9)
                    {
                        mu += probe_data[i];
                        nonmiss += 1.0;
                    }
                }
                if(nonmiss>0) mu /=nonmiss;
                for(int i=0; i<n; i++)
                    if(probe_data[i]>=1e9) probe_data[i] = mu;
            }
            double var_x = var(probe_data);
            double b = var_x * rsq / (var_y*(1-rsq));
            reverse_eff[i] = b;
            for(int j=0; j<n; j++) einfo._val[einfo._epi_include[i]*einfo._eii_num+einfo._eii_include[j]] = b*einfo._eii_pheno[j] + probe_data[j];
        }
        
        output_simu_par(outFileName, qtl_name, reverse_eff);
        
        einfo._epi_include=include_o;
        einfo._epi_map=map_o;
        write_eii(outFileName, &einfo);
        write_epi(outFileName, &einfo);
        write_beed(outFileName, &einfo);
        cout << "Simulated probes along with the rest probes of " << einfo._eii_include.size() << " individuals have been saved in [" <<  string(outFileName) <<"]." << endl;
    }
    
    void drop_comp(eInfo* einfo, vector<int> &drop) {
        int i = 0;
        stringstream errmsg;
        einfo->_r_indx_drop = einfo->_r_indx;
        stable_sort(drop.begin(), drop.end());
        drop.erase(unique(drop.begin(), drop.end()), drop.end());
        for (i = drop.size() - 1; i >= 0; i--) {
            if (drop[i] < 1 || drop[i] > einfo->_r_indx.size() - 1) {
                errmsg << "Error: there " << (einfo->_r_indx.size() > 2 ? "are" : "is") << " only " << einfo->_r_indx.size() - 1 << " genetic variance component in the model. You can't drop the " << drop[i] << "-th component.";
                throw (errmsg.str());
            }
            einfo->_r_indx_drop.erase(einfo->_r_indx_drop.begin() + drop[i] - 1);
        }
        if (einfo->_r_indx.size() == einfo->_r_indx_drop.size()) throw ("Error: no component has been dropped from the model. Please check the --reml-lrt option.");
    }
    
    
    void fit_reml(char* outFileName, char* befileName, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde, char* phenofileName,char* mpheno,bool erm_bin_flag, bool grm_bin_flag,int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* grm_file, bool m_erm_flag, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool pred_rand_eff, bool est_fix_eff,bool no_lrt,double prevalence, bool mlmassoc,vector<int> drop, char* indilstName, char* indilst2remove, char* sex_file, double grm_cutoff, bool erm_cutoff_2sides, double adj_grm_fac, int dosage_compen,bool prt_residiual)
    {
        loud = true;
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        if(priors!=NULL){
            int tmpint=split_string(priors, vs_buf, ",");
            for(int i=0;i<tmpint;i++)
            {
                double tmpd=atof(vs_buf[i].c_str());
                if(tmpd<0 || tmpd>1)
                {
                    LOGPRINTF("Error:  --reml-priors. Prior values of variance explained should be between 0 and 1.\n");
                    TERMINATE();
                }
                reml_priors.push_back(tmpd);
            }
        }
        if(priors_var!=NULL){
            int tmpint=split_string(priors_var, vs_buf, ",");
            for(int i=0;i<tmpint;i++)
            {
                double tmpd=atof(vs_buf[i].c_str());
                if(tmpd<0 )
                {
                    LOGPRINTF("Error:  --reml-priors-var. Prior values of variance components should be over 0.\n");
                    TERMINATE();
                }
                reml_priors_var.push_back(tmpd);
            }
        }
        if(prt_residiual) est_fix_eff=true;
        
        eInfo einfo;
        init_einfo(&einfo);
        einfo._reml_mtd=reml_mtd;
        einfo._reml_max_iter=MaxIter;
        einfo._reml_fixed_var=reml_fixed_var_flag;
        einfo._reml_force_inv=reml_force_inv_fac_flag;
        einfo._reml_force_converge=reml_force_converge_flag;
        einfo._reml_no_converge=reml_no_converge_flag;
        if(prevalence>-1) einfo._flag_CC=true;
        
        vector<string> grm_id;
        vector<string> erm_files;
        char inputname[FNAMESIZE];
        char* suffix = NULL;
        if(befileName != NULL)
        {
            memcpy(inputname,befileName,strlen(befileName)+1);
            suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,&einfo);
            eii_man(&einfo,indilstName,indilst2remove);
            memcpy(suffix,".opi",5);
            read_epi(inputname,&einfo);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            erm_files.push_back(befileName);
        } else {
            if (m_erm_flag || grm_file!=NULL) {
                read_msglist(erm_file, erm_files,"ORM file names");
                if(grm_file!=NULL) {
                    erm_files.push_back(grm_file);
                }
                for (int i = 0; i < erm_files.size(); i++) {
                    read_grm(&einfo,erm_files[i], grm_id, false, true, true,erm_bin_flag);
                }
            } else if(erm_file!=NULL){
                erm_files.push_back(erm_file);
                read_grm(&einfo,erm_file, grm_id, true, false, true,erm_bin_flag);
            }
        }
        
        
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        else {
            LOGPRINTF("ERROR: no phenoytpe file inputed.\n");
            TERMINATE();
        }
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        if (indilstName!=NULL) keep_indi(&einfo, indilstName);
        if (indilst2remove!=NULL) remove_indi(&einfo, indilst2remove);
        if(befileName != NULL)
        {
            memcpy(suffix,".bod",5);
            read_beed(inputname,&einfo);
            make_erm(&einfo, erm_alg);
            for(int i=0; i<einfo._eii_include.size(); i++)
                grm_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
        } else {
            if(erm_file!=NULL) {
                if (grm_cutoff>-1.0) rm_cor_indi(&einfo,grm_cutoff,erm_cutoff_2sides);
                if (sex_file!=NULL) update_sex(&einfo,sex_file);
                if (adj_grm_fac>-1.0) adj_grm(&einfo,adj_grm_fac);
                if (dosage_compen>-1) dc(&einfo,dosage_compen);
            }
        }
        
        
        vector<string> uni_id;
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for(int i=0; i<einfo._eii_include.size(); i++){
            uni_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
            uni_id_map.insert(pair<string,int>(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]], i));
        }
        
        int _n=(int)einfo._eii_include.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        
        einfo._r_indx.clear();
        vector<MatrixXd> _A;
        vector<int> kp;
      
             if(m_erm_flag){
                cout << "There are " << erm_files.size() << " ORM file names specified in the file [" + string(erm_file) + "]." << endl;
                for(int i=0; i < (1+0+0)*erm_files.size() + 1; i++) einfo._r_indx.push_back(i);
                if (!no_lrt) drop_comp(&einfo,drop);
                _A.resize(einfo._r_indx.size());
                for (int i = 0; i < erm_files.size(); i++) {
                    cout << "Reading the ORM from the " << i + 1 << "th file ..." << endl;
                    read_grm(&einfo, erm_files[i], grm_id, true, false, true,erm_bin_flag);
                    match(uni_id, grm_id, kp);
                    (_A[i]).resize(_n, _n);
                    for (int j = 0; j < _n; j++) {
                        for (int k = 0; k <= j; k++) {
                            if (kp[j] >= kp[k]) (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[j], kp[k]);
                            else (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[k], kp[j]);
                        }
                    }
                }
            } else if(erm_file!=NULL || befileName != NULL){
                for(int i=0; i < 1 + 1; i++) einfo._r_indx.push_back(i);
                if (!no_lrt) drop_comp(&einfo,drop);
                _A.resize(einfo._r_indx.size());
                match(uni_id, grm_id, kp);
                (_A[0]).resize(_n, _n);
                for(int i=0; i<_n; i++){
                    for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(kp[i],kp[j]);
                }
                einfo._grm.resize(0,0);
            }
            else{
                einfo._r_indx.push_back(0);
                 _A.resize(einfo._r_indx.size());
            }
        
        _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
        
        VectorXd _y; // only for univariate now, using _eii_pheno_num later for multiple variates.
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=einfo._eii_pheno[einfo._eii_include[i]];
        }
        if(prevalence>-1) check_case_control(einfo._ncase,  _y);
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(&einfo, E_float, qE_float,_X);
        
        // names of variance component
        for (int i = 0; i < erm_files.size(); i++) {
            stringstream strstrm;
            if (erm_files.size() == 1) strstrm << "";
            else strstrm << i + 1;
            einfo._var_name.push_back("V(O" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        einfo._within_family=within_family;
        if(within_family) detect_family(&einfo, _A);
        
        // run REML algorithm
        MatrixXd _Vi;
        reml(&einfo, pred_rand_eff, est_fix_eff, reml_priors, reml_priors_var, prevalence, -2.0, no_constrain, no_lrt, mlmassoc, _X_c, _X,_y,_A,_Vi,outFileName);
        
        
        if(prt_residiual) {
            VectorXd y_buf=_y.array()-(_X*einfo._b).array(); // adjust phenotype for covariates
            
            if(y_buf.size()!=einfo._eii_include.size()) {
                LOGPRINTF("ERROR: the individual number of adjusted phenotype %ld is inconsistent with the individual numebr seleced %ld.\nplease report this bug.\n",y_buf.size(),einfo._eii_include.size());
                TERMINATE();
            }
            string adjphenfilename=string(outFileName)+".adj.phen";
            LOGPRINTF("\nSaving adjusted non-missing phenotype to %s ...\n",adjphenfilename.c_str());
            ofstream phfile(adjphenfilename.c_str());
            if(!phfile) throw("Can not open the file ["+adjphenfilename+"] to write.");
            for(int i=0; i<y_buf.size(); i++){
                phfile<<einfo._eii_fid[einfo._eii_include[i]]<<"\t"<<einfo._eii_iid[einfo._eii_include[i]]<<"\t"<<y_buf(i)<<endl;
            }
            phfile.close();
            
        }

    }
    
    void extract_inden_probes(char* outFileName, char* befileName, int extr_num, double ldr,  int seed,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool rand_eff)
    {
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        
     
        vector<string> slct_prbs;
        map<string,int> slct_prbs_map;
        map<string,int>::iterator iter;
        int map_size=0;
        vector<int> slct_idx;
        int Seed = -rand_seed();
        srand(Seed);
        vector<double> x, y, rst;
        long repcount=0;
        while(slct_prbs.size()<extr_num) {
            if(repcount>ceil(1.5*einfo._eii_include.size()))
            {
                LOGPRINTF("ERROR: total probe number is too small or independent probes to extract are too many. \n");
                LOGPRINTF("Please loosen the input correlation R-squer to have another try. \n");
                TERMINATE();
            }
            int idx=rand()%einfo._epi_include.size();
            if(slct_prbs.size()==0) {
                slct_idx.push_back(idx);
                slct_prbs.push_back(einfo._epi_prb[einfo._epi_include[idx]]);
                slct_prbs_map.insert(pair<string,int>(einfo._epi_prb[einfo._epi_include[idx]],map_size));
                map_size=(int)slct_prbs_map.size();
            } else {
                bool found=true;
                for(int j=0;j<slct_prbs.size();j++)
                {
                    x.clear();
                    y.clear();
                    for(int k=0;k<einfo._eii_include.size();k++) {
                        double valy=einfo._val[einfo._epi_include[idx]*einfo._eii_num+einfo._eii_include[k]];
                        double valx=einfo._val[einfo._epi_include[slct_idx[j]]*einfo._eii_num+einfo._eii_include[k]];
                        if(valy<1e9 && valx<1e9){
                            y.push_back(valy);
                            x.push_back(valx);
                        }
                    }
                    reg(y, x, rst);
                    if(rst[3]>ldr) found=false;
                }
                if(found){
                    slct_prbs_map.insert(pair<string,int>(einfo._epi_prb[einfo._epi_include[idx]],map_size));
                    if(slct_prbs_map.size()>map_size)
                    {
                        slct_idx.push_back(idx);
                        slct_prbs.push_back(einfo._epi_prb[einfo._epi_include[idx]]);
                        map_size=(int)slct_prbs_map.size();
                    } else {
                        LOGPRINTF("OOPS: The same probe was picked up. Resample again. \n");
                    }
                    
                }
            }
            repcount++;
        }
        vector<double> qtl_eff;
        if(rand_eff)
        {
            Seed = -rand_seed();
            for (int i = 0; i < slct_prbs.size(); i++) qtl_eff.push_back(gasdev(Seed));
            
        }
        // Output par file
        
        FILE* ofile=NULL;
        ofile = fopen(outFileName, "w");
        if (!(ofile)) {
            LOGPRINTF("Open error %s\n", outFileName);
            TERMINATE();
        }
        string outstr="";
        for(int i=0;i<slct_prbs.size();i++)
        {
            if(rand_eff) outstr=slct_prbs[i]+'\t'+atos(qtl_eff[i])+'\n';
            else outstr=slct_prbs[i]+'\n';
            if(fputs_checked(outstr.c_str(),ofile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                TERMINATE();
            }
        }
        fclose(ofile);
        LOGPRINTF("%ld probes have been saved in file %s.\n",slct_prbs.size(),outFileName);
    }
    void getPrbVarianceMean(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool varflag,bool meanflag)
    {
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL && efileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile.\n");
            TERMINATE();
        }
        if(efileName!=NULL)
        {
            if(transposed) read_efile_t(efileName,&einfo,efileType,no_fid_flag,valueType);
            else read_efile(efileName,&einfo,efileType,no_fid_flag,valueType);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            eii_man(&einfo,indilstName,indilst2remove);
        }else{
            char inputname[FNAMESIZE];
            memcpy(inputname,befileName,strlen(befileName)+1);
            char* suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,&einfo);
            eii_man(&einfo,indilstName,indilst2remove);
            memcpy(suffix,".opi",5);
            read_epi(inputname,&einfo);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            memcpy(suffix,".bod",5);
            read_beed(inputname,&einfo);
        }
        
        LOGPRINTF("Calculating variance of probes ...\n");
        FILE* ofile=NULL;
        FILE* mfile=NULL;
        char outputname[FNAMESIZE];
        memcpy(outputname,outFileName,strlen(outFileName)+1);
        char* suffix=outputname+strlen(outFileName);
        if(varflag){
            memcpy(suffix,".var.txt",9);
            ofile = fopen(outputname, "w");
            if (!(ofile)) {
                LOGPRINTF("Open error %s\n", outputname);
                TERMINATE();
            }
        }
        if(meanflag){
            memcpy(suffix,".mean.txt",10);
            mfile = fopen(outputname, "w");
            if (!(mfile)) {
                LOGPRINTF("Open error %s\n", outputname);
                TERMINATE();
            }
        }
        string outstr="";
        long n = einfo._eii_include.size(), m = einfo._epi_include.size();
        uint64_t ttl_n = einfo._eii_num;
        for(int j=0; j<m; j++){
            double mu=0.0;
            double nonmiss=0.0;
            for(int i=0; i<n; i++){
                double val=einfo._val[einfo._epi_include[j]*ttl_n+einfo._eii_include[i]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            
            double sd =0.0;
            for(int i=0; i<n; i++){
                double val=einfo._val[einfo._epi_include[j]*ttl_n+einfo._eii_include[i]];
                if(val<1e9) sd += (val-mu)*(val-mu);
            }
            sd/=(nonmiss - 1.0);
            if(varflag){
                 outstr=einfo._epi_prb[einfo._epi_include[j]]+'\t'+atos(sd)+'\n';
                 if(fputs_checked(outstr.c_str(),ofile))
                 {
                     LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                     TERMINATE();
                 }
             }
            if(meanflag)
            {
                outstr=einfo._epi_prb[einfo._epi_include[j]]+'\t'+atos(mu)+'\n';
                if(fputs_checked(outstr.c_str(),mfile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                    TERMINATE();
                }

            }
        }
        if(ofile) fclose(ofile);
        if(mfile) fclose(mfile);
        if(varflag) LOGPRINTF("Variances of %ld probes have been save in file %s.\n", einfo._epi_include.size(), (string(outFileName)+".var.txt").c_str());
        if(meanflag) LOGPRINTF("Mean values of %ld probes have been save in file %s.\n", einfo._epi_include.size(), (string(outFileName)+".mean.txt").c_str());

    }

    void blup_probe(char* outFileName, char* efileName, char* befileName,bool transposed, int efileType, char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,double std_thresh,double upperBeta,double lowerBeta, char* blup_indi_file)
    {
        eInfo einfo;
        init_einfo(&einfo);
        if(befileName==NULL && efileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --efile or --befile.\n");
            TERMINATE();
        }
        if(efileName!=NULL)
        {
            if(transposed) read_efile_t(efileName,&einfo,efileType,no_fid_flag,valueType);
            else read_efile(efileName,&einfo,efileType,no_fid_flag,valueType);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            eii_man(&einfo,indilstName,indilst2remove);
            if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        }else{
            char inputname[FNAMESIZE];
            memcpy(inputname,befileName,strlen(befileName)+1);
            char* suffix=inputname+strlen(befileName);
            memcpy(suffix,".oii",5);
            read_eii(inputname,&einfo);
            eii_man(&einfo,indilstName,indilst2remove);
            memcpy(suffix,".opi",5);
            read_epi(inputname,&einfo);
            epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
            memcpy(suffix,".bod",5);
            read_beed(inputname,&einfo);
            if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
            if(upperBeta<1 ||lowerBeta>0) filtering_constitutive_probes(&einfo, upperBeta, lowerBeta);
        }
        
        blup_probe_geno(&einfo, outFileName,blup_indi_file);

    }
    
    /* betaFileName can be the output file of an association analysis or BLUP probe.
     *
     */
    void read_eff(char* betaFileName,int col_prb, int col_score, bool hasHeader, eInfo* einfo, vector<double> &scores)
    {
        int cpp_pos_prb=col_prb-1;
        int cpp_pos_score=col_score-1;
        scores.resize(einfo->_epi_prb.size()); // according to _epi_prb.size(), not _epi_include.size()
        if(einfo->_epi_include.size()==0 )
        {
            LOGPRINTF("Error: no probe found in gene expression / methylation data.\n");
            TERMINATE();
        }
        FILE* betafile=NULL;
        map<string, int>::iterator iter;
        vector<string> strlist;
        vector<string> scprbs;
        int line_idx = 0;
        int probe_match_count=0;
        int missing_score_count=0;
        if(fopen_checked(&betafile, betaFileName,"r")) {
             LOGPRINTF("Error: can not open the file %s to read.\n",betaFileName);
            TERMINATE();
        }
        LOGPRINTF("Reading effect size information from %s ...\n", betaFileName);
        long colnum=0;
        if(hasHeader) {
            if(fgets(Tbuf, MAX_LINE_SIZE, betafile)){
                split_str(Tbuf,strlist,0);
                if(strlist.size()<2)
                {
                    LOGPRINTF("ERROR: Total number of columns is less than 2.\n");
                    TERMINATE();
                }
                if(strlist.size()<col_prb)
                {
                    LOGPRINTF("ERROR: Total number of columns is less than specified probe column %d.\n", col_prb);
                    TERMINATE();
                }
                if(strlist.size()<col_score)
                {
                    LOGPRINTF("ERROR: Total number of columns is less than specified score column %d.\n", col_score);
                    TERMINATE();
                }
                colnum=strlist.size();
            } else {
                LOGPRINTF("ERROR: File %s read failed.\n", betaFileName);
                TERMINATE();
            }
            
        }
        memset(Tbuf,0,sizeof(char)*(MAX_LINE_SIZE));
        while(fgets(Tbuf, MAX_LINE_SIZE, betafile))
        {
            split_str(Tbuf,strlist,0);
            if( line_idx==0 && !hasHeader){
                // if no header and the first row
                if(strlist.size()<2)
                {
                    LOGPRINTF("ERROR: Total number of columns is less than 2.\n");
                    TERMINATE();
                }
                if(strlist.size()<col_prb)
                {
                    LOGPRINTF("ERROR: Total number of columns is less than specified probe column %d.\n", col_prb);
                    TERMINATE();
                }
                if(strlist.size()<col_score)
                {
                    LOGPRINTF("ERROR: Total number of columns is less than specified score column %d.\n", col_score);
                    TERMINATE();
                }
                colnum=strlist.size();
            } else {
                if(strlist.size()<colnum)
                {
                    LOGPRINTF("ERROR: Line %u doesn't have the same number of items (%ld) as the header.\n", line_idx, colnum);
                    TERMINATE();
                }
            }
            string scostr=strlist[cpp_pos_score];
            to_upper(scostr);
            if(scostr=="NA") {
                missing_score_count++;
                continue;
            }
            iter = einfo->_epi_map.find(strlist[cpp_pos_prb].c_str());
            if (iter != einfo->_epi_map.end()){
                probe_match_count++;
                int curpIdx=iter->second;
                scprbs.push_back(strlist[cpp_pos_prb].c_str());
                scores[curpIdx]=atof(scostr.c_str());
            }
            line_idx++;
            memset(Tbuf,0,sizeof(char)*(MAX_LINE_SIZE));
        }
        string typestr=getFileType(einfo->_eType);
        LOGPRINTF("%u probes have been read from the file %s.\n", line_idx, betaFileName);
        LOGPRINTF("among %u probes %d probes have missing value and were excluded.\n", line_idx, missing_score_count);
        LOGPRINTF("%d probes were matched with %ld probes of the %s data.\n", probe_match_count, einfo->_epi_include.size(),typestr.c_str());
        LOGPRINTF("%d probes were unmatched and excluded.\n", line_idx-missing_score_count-probe_match_count);
        update_map_kp(scprbs, einfo->_epi_map, einfo->_epi_include);
        if(einfo->_epi_include.size()<einfo->_epi_num)
        {
            vector<double> doubletmp;
            for (int i = 0; i < einfo->_epi_include.size(); i++) doubletmp.push_back(scores[einfo->_epi_include[i]]);
            scores.clear();
            scores.swap(doubletmp);
        }
    }
    void calculateProfile(eInfo* einfo,vector<double> scores, vector<double> &profile,vector<int> &cnt,bool impute_mean_flag)
    {
        profile.resize(einfo->_eii_include.size());
        cnt.resize(einfo->_eii_include.size());
        if (einfo->_mu.empty()) cal_var_mean(einfo, true, false);
        for(int i=0;i<einfo->_eii_include.size();i++)
        {
            double score=0.0;
            int count=0;
            for(int j=0;j<einfo->_epi_include.size();j++)
            {
                double val=einfo->_val[einfo->_epi_include[j]*einfo->_eii_num+einfo->_eii_include[i]];
                if(val<1e9){
                    score+=val*scores[einfo->_epi_include[j]];
                    count++;
                } else {
                    if(impute_mean_flag) {
                        score+=einfo->_mu[einfo->_epi_include[j]]*scores[einfo->_epi_include[j]];
                        count++;
                    }
                }
            }
            score/=(double)count;
            profile[i]=score;
            cnt[i]=count;
        }
        
    }
    void scoreIndividuals(char* outFileName,char* befileName, char* betaFileName, int col_prb, int col_score, bool hasHeader, char* phenofileName, char* mpheno,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, double std_thresh, bool impute_mean_flag)
    {
        eInfo einfo;
        vector<double> scores;
        init_einfo(&einfo);
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(betaFileName==NULL)
        {
            LOGPRINTF("Error: please input the probe effect data by the option --score.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        read_eff(betaFileName,col_prb, col_score, hasHeader, &einfo,scores);
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo);
        if(std_thresh>0) std_probe_filtering( &einfo, std_thresh);
        
        vector<double> profile;
        vector<int> cnt;
        calculateProfile(&einfo,scores, profile,cnt,impute_mean_flag);
        
        FILE* ofile=NULL;
        string outFN=string(outFileName)+".profile";
        if(fopen_checked(&ofile, outFN.c_str(),"w")) TERMINATE();
        string str ="FID\tIID\tPHENO\tCNT\tSCORE\n";
        if(fputs_checked(str.c_str(),ofile))
        {
            LOGPRINTF("ERROR: in writing file %s .\n", outFN.c_str());
            TERMINATE();
        }
        for(int i=0;i<einfo._eii_include.size();i++)
        {
            str=atos(einfo._eii_fid[einfo._eii_include[i]])+'\t'+einfo._eii_iid[einfo._eii_include[i]]+'\t'+((phenofileName==NULL)?"-9":atos(einfo._eii_pheno[einfo._eii_include[i]]))+'\t'+atos(cnt[i])+'\t'+dtos(profile[i])+'\n';
            if(fputs_checked(str.c_str(),ofile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", outFN.c_str());
                TERMINATE();
            }
        }
        fclose(ofile);
        LOGPRINTF("%ld profiles have been saved in the file %s .\n", einfo._eii_include.size(), outFN.c_str());

    }
    
    void calcu_r_rsq(char* outFileName, char* efileName, char* befileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, int wind_size, int autosome_num, bool rflag, bool r2flag)
    {
        eInfo einfo;
        init_einfo(&einfo);
        load_workspace(&einfo, efileName, befileName, transposed, efileType, problstName, problst2exclde, genelistName, chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag, valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName, dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob,autosome_num);
        
        for (int i = 0; i < einfo._epi_include.size(); i++) {
            if (einfo._epi_chr[einfo._epi_include[i]] > einfo.autosome_num) {
                LOGPRINTF("Error: this option is for the autosomal probes only. Please check the option --autosome.\n");
                TERMINATE();
            }
        }
        write_epi(outFileName, &einfo);
        string bldname=string(outFileName)+".bcor";
        FILE* outfile=fopen(bldname.c_str(), "wb");
        if (!(outfile)) {
            printf("Error: Failed to open file %s.\n",bldname.c_str());
            exit(EXIT_FAILURE);
        }

        
    }
    void moa(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* subtract_erm_file, bool m_erm_flag,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool force_mlm,int tsk_ttl,int tsk_id)
    {
        eInfo einfo;
        init_einfo(&einfo);
        einfo._reml_mtd=reml_mtd;
        einfo._reml_max_iter=MaxIter;
        einfo._reml_fixed_var=reml_fixed_var_flag;
        einfo._reml_force_inv=reml_force_inv_fac_flag;
        einfo._reml_force_converge=reml_force_converge_flag;
        einfo._reml_no_converge=reml_no_converge_flag;
        einfo._within_family=false;
        
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("Error: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }
        
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".oii",5);
        read_eii(inputname,&einfo);
        eii_man(&einfo,indilstName,indilst2remove);
        memcpy(suffix,".opi",5);
        read_epi(inputname,&einfo);
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo); // eii and epi are updated in it.
        
        vector<string> grm_id;
        vector<string> erm_files;
            if(subtract_erm_file){
                erm_files.push_back(erm_file);
                erm_files.push_back(subtract_erm_file);
                for (int i = 0; i < erm_files.size(); i++) {
                    read_grm(&einfo, erm_files[i], grm_id, false, true, true,erm_bin_flag);  // maybe rm grm_id later
                }
            }
            else{
                if(erm_file!=NULL && !m_erm_flag){
                    erm_files.push_back(erm_file);
                    read_grm(&einfo,erm_file, grm_id, true, false, true,erm_bin_flag);
                }
                else if (erm_file!=NULL && m_erm_flag) {
                    read_msglist(erm_file, erm_files,"ORM file names");
                    for (int i = 0; i < erm_files.size(); i++) {
                        read_grm(&einfo,erm_files[i], grm_id, false, true, true,erm_bin_flag);
                    }
                }
                else{
                    erm_files.push_back("NA");
                    make_erm( &einfo,erm_alg);
                    for(int i=0; i<einfo._eii_include.size(); i++) grm_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
                }
            }
        vector<string> uni_id;
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for(int i=0; i<einfo._eii_include.size(); i++){
            uni_id.push_back(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]]);
            uni_id_map.insert(pair<string,int>(einfo._eii_fid[einfo._eii_include[i]]+":"+einfo._eii_iid[einfo._eii_include[i]], i));
        }
        int _n=(int)einfo._eii_include.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        einfo._r_indx.clear();
        vector<MatrixXd> _A;
        vector<int> kp;
      
            if (subtract_erm_file) {
                for(int i=0; i < 2; i++) einfo._r_indx.push_back(i);
                _A.resize(einfo._r_indx.size());
                
                LOGPRINTF("\nReading the primary ORM from %s ...\n",erm_files[1].c_str() );
                read_grm(&einfo, erm_files[1], grm_id, true, false, false,erm_bin_flag);
                
                match(uni_id, grm_id, kp);
                (_A[0]).resize(_n, _n);
                MatrixXf A_N_buf(_n, _n);
                #pragma omp parallel for
                for (int j = 0; j < _n; j++) {
                    for (int k = 0; k <= j; k++) {
                        if (kp[j] >= kp[k]){
                            (_A[0])(k, j) = (_A[0])(j, k) = einfo._grm(kp[j], kp[k]);
                            A_N_buf(k, j) = A_N_buf(j, k) = einfo._grm_N(kp[j], kp[k]);
                        }
                        else{
                            (_A[0])(k, j) = (_A[0])(j, k) = einfo._grm(kp[k], kp[j]);
                            A_N_buf(k, j) = A_N_buf(j, k) = einfo._grm_N(kp[k], kp[j]);
                        }
                    }
                }
                
                LOGPRINTF("\nReading the secondary ORM from %s ...\n",erm_files[0].c_str());
                read_grm(&einfo, erm_files[0], grm_id, true, false, false,erm_bin_flag);
                LOGPRINTF("\nSubtracting %s from %s ...\n",erm_files[1].c_str(),erm_files[0].c_str());
                match(uni_id, grm_id, kp);
                #pragma omp parallel for
                for (int j = 0; j < _n; j++) {
                    for (int k = 0; k <= j; k++) {
                        if (kp[j] >= kp[k]) (_A[0])(k, j) = (_A[0])(j, k) = ((_A[0])(j, k) * A_N_buf(j, k)  - einfo._grm(kp[j], kp[k]) * einfo._grm_N(kp[j], kp[k])) / (A_N_buf(j, k) - einfo._grm_N(kp[j], kp[k]));
                        else (_A[0])(k, j) = (_A[0])(j, k) = ((_A[0])(j, k) * A_N_buf(j, k) - einfo._grm(kp[k], kp[j]) * einfo._grm_N(kp[k], kp[j])) / (A_N_buf(j, k) - einfo._grm_N(kp[k], kp[j]));
                    }
                }
                einfo._grm.resize(0,0);
                einfo._grm_N.resize(0,0);
            }
            else {
                for(int i=0; i < erm_files.size() + 1; i++) einfo._r_indx.push_back(i);
                _A.resize(einfo._r_indx.size());
                if(erm_file!=NULL && !m_erm_flag){
                    match(uni_id, grm_id, kp);
                    (_A[0]).resize(_n, _n);
                    #pragma omp parallel for
                    for(int i=0; i<_n; i++){
                        for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(kp[i],kp[j]);
                    }
                    einfo._grm.resize(0,0);
                }
                else if(erm_file!=NULL && m_erm_flag){
                    LOGPRINTF("There are %ld ORM file names specified in the file %s.\n",erm_files.size(),erm_file);
                    for (int i = 0; i < erm_files.size(); i++) {
                        LOGPRINTF("Reading the ORM from the %dth file ...\n",i + 1);
                        read_grm(&einfo, erm_files[i], grm_id, true, false, true,erm_bin_flag);
                        match(uni_id, grm_id, kp);
                        (_A[i]).resize(_n, _n);
                        #pragma omp parallel for
                        for (int j = 0; j < _n; j++) {
                            for (int k = 0; k <= j; k++) {
                                if (kp[j] >= kp[k]) (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[j], kp[k]);
                                else (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[k], kp[j]);
                            }
                        }
                    }
                }
                else{
                    match(uni_id, grm_id, kp);
                    (_A[0]).resize(_n, _n);
                    #pragma omp parallel for
                    for(int i=0; i<_n; i++){
                        for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(kp[i],kp[j]);
                    }
                    //einfo._grm.resize(0,0);
                }
            }
        
        _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
        SelfAdjointEigenSolver<MatrixXd> eigensolver(_A[0]);
        VectorXd eval = eigensolver.eigenvalues();
        MatrixXd U=eigensolver.eigenvectors();
        
        VectorXd _y; // only for univariate now, using _eii_pheno_num later for multiple variates.
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=einfo._eii_pheno[einfo._eii_include[i]];
        }
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(&einfo, E_float, qE_float,_X); //_X has the dimension of [einfo._eii_include.size(),_X_c]
        
        // names of variance component
        for (int i = 0; i < erm_files.size(); i++) {
            stringstream strstrm;
            if (erm_files.size() == 1) strstrm << "";
            else strstrm << i + 1;
            einfo._var_name.push_back("V(O" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        if(tsk_ttl>1) extract_probe(&einfo,  tsk_ttl,  tsk_id);
        vector<double> betas(einfo._epi_include.size()),ses(einfo._epi_include.size());
        //#pragma omp parallel
        //{
            MatrixXd XX(_X.rows(), _X.cols()+1);
            XX.block(0,0,_X.rows(),_X.cols())=_X;
            int x_idx=_X_c++;
            double cr=0;
            #pragma omp parallel for private(remlstatus, cr)
            for(int jj=0;jj<einfo._epi_include.size();jj++)
            {
                double desti=1.0*jj/(einfo._epi_include.size()-1);
                if(desti>=cr)
                {
                    printf("%3.0f%%\r", 100.0*desti);
                    fflush(stdout);
                    if(cr==0) cr+=0.05;
                    else if(cr==0.05) cr+=0.2;
                    else if(cr==0.25) cr+=0.5;
                    else cr+=0.25;
                }
                
                string prbid=einfo._epi_prb[einfo._epi_include[jj]];
                vector<double> reml_priors,reml_priors_var;
                MatrixXd X=XX;
                double nonmiss=0.0;
                double mu=0.0, sd=0.0;
                long n=einfo._eii_include.size();
                VectorXd x(n);
                //calc the mean
                for(int j = 0; j < n; j++) {
                    double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[j]];
                    if(val<1e9){
                        mu+=val;
                        nonmiss+=1.0;
                    }
                }
                if(nonmiss>1)
                {
                    mu/=nonmiss;
                    for(int j = 0; j < n; j++) {
                        double val=einfo._val[einfo._epi_include[jj]*einfo._eii_num+einfo._eii_include[j]];
                        if(val<1e9){
                            double centval=val-mu;
                            x(j) = centval;
                            sd+=centval*centval;
                        } else {
                            x(j)=0.0;
                        }
                    }
                    sd=sqrt(sd/(nonmiss-1.0));
                    if(sd>1e-30)
                    {
                        for(int j=0; j<n; j++) x(j)/=sd;
                    }
                    else
                    {
                        LOGPRINTF("ERROR: no variance of probe %s.\n",prbid.c_str());
                        //do something here;
                        betas[jj]=0;
                        ses[jj]=1;
                        continue;
                    }
                }
                else
                {
                    LOGPRINTF("ERROR: too many missing values in probe %s.\n",prbid.c_str());
                    TERMINATE();
                }
                X.col(x_idx)=x;
            
                remlstatus=0; //reset reml status
                MatrixXd _Vi;
                VectorXd _b,_se;
                reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,X, _y,_A, U, eval, _Vi,  reml_mtd,  MaxIter,_b,_se);
                //reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, X,_y,_A,_Vi,outFileName);
                if(remlstatus==0 || remlstatus==-5 || remlstatus==-3)
                {
                    betas[jj]=_b[_X_c-1];
                    ses[jj]=sqrt(_se(_X_c-1));
                   
                }
                else
                {
                    if(remlstatus==-1) {
                        LOGPRINTF("ERROR: The matrix is not invertible.\n");
                    } else if(remlstatus==-2) {
                        LOGPRINTF("ERROR: More than half of the variance components are constrained.\n");
                    }else if(remlstatus==-4) {
                        LOGPRINTF("ERROR:Log-likelihood not converged.\n");
                    }
                    TERMINATE();
                    //probes which failed in REML
                    //betas[jj]=0;
                    //ses[jj]=1;
                }
            }
       // }
        
            string filename=string(outFileName)+".moa";
            if(tsk_ttl>1) filename=string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id)+".moa";
            LOGPRINTF("\nSaving the association analysis results of %ld probes to %s ...\n",einfo._epi_include.size(),filename.c_str());
            ofstream ofile(filename.c_str());
            if(!ofile) throw("Can not open the file ["+filename+"] to write.");
            ofile<<"Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp"<<endl;
            for(int i=0; i<einfo._epi_include.size(); i++){
                int j=einfo._epi_include[i];
                string chrstr;
                if(einfo._epi_chr[j]==23) chrstr="X";
                else if(einfo._epi_chr[j]==24) chrstr="Y";
                else chrstr=atosm(einfo._epi_chr[j]);
                double z2=betas[i]/ses[i];
                z2*=z2;
                double pval=pchisq(z2, 1);
                ofile<<chrstr<<"\t"<<einfo._epi_prb[j]<<"\t"<<einfo._epi_bp[j]<<"\t"<<einfo._epi_gene[j]<<"\t"<<(einfo._epi_orien[j]=='*'?"NA":atos(einfo._epi_orien[j]))<<"\t";
                if(pval>1.5) ofile<<"NA\tNA\tNA"<<endl;
                else ofile<<betas[i]<<"\t"<<ses[i]<<"\t"<<pval<<endl;
            }
            ofile.close();
            
        }
    
    
}

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
    void combine_eii(vector<indiinfolst> &indiinfo, vector<string> &befNames, bool &indi_uni)
    {
        long counter = 0;
        map<string, int> indi_map;
        long f2r=befNames.size();
        char inputname[FNAMESIZE];
        bool diffidsflg=false;
        for (int i = 0; i < f2r; i++)
        {
            eInfo etmp;
            string eiifile = befNames[i]+".oii";
            strncpy(inputname, eiifile.c_str(), sizeof(inputname));
            inputname[sizeof(inputname) - 1] = 0;
            read_eii(inputname, &etmp);

            for (int j = 0; j<etmp._eii_num; j++)
            {
                string fiidstr=etmp._eii_fid[j]+":"+atos(etmp._eii_iid[j]);
                indi_map.insert(pair<string, int>(fiidstr.c_str(), counter));
             
                if (counter < indi_map.size())
                {
                    indiinfolst iinfotmp;
                    counter = indi_map.size();
                    strcpy2(&iinfotmp.fid, etmp._eii_fid[j]);
                    strcpy2(&iinfotmp.iid, etmp._eii_iid[j]);
                    strcpy2(&iinfotmp.fa_id, etmp._eii_fa_id[j]);
                    strcpy2(&iinfotmp.mo_id, etmp._eii_mo_id[j]);
                    iinfotmp.sex=etmp._eii_sex[j];
                    iinfotmp.pheno=etmp._eii_pheno[j];
                    indiinfo.push_back(iinfotmp);
                    if(i && !diffidsflg) diffidsflg=true;
                }
            }
        }
        if(!diffidsflg) indi_uni=true;
        LOGPRINTF("Total %ld individuals to be included from %ld eii files.\n",indiinfo.size(),f2r);
    }
    void combine_epi(vector<probeinfolst> &probeinfo, vector<string> &befNames, bool &prb_uni) {
        
        long counter = 0;
        map<string, int> prb_map;
        map<string, int>::iterator iter;
        long f2r=befNames.size();
        char inputname[FNAMESIZE];
        bool diffidsflg=false;
        for (int i = 0; i < f2r; i++)
        {
            eInfo etmp;
            string eiifile = befNames[i]+".opi";
            strncpy(inputname, eiifile.c_str(), sizeof(inputname));
            inputname[sizeof(inputname) - 1] = 0;
            read_epi(inputname, &etmp);
            
            for (int j = 0; j<etmp._epi_num; j++)
            {
                string pidstr=etmp._epi_prb[j];
                prb_map.insert(pair<string, int>(pidstr.c_str(), counter));
                
                if (counter < prb_map.size())
                {
                    probeinfolst pinfotmp;
                    counter = prb_map.size();
                    strcpy2(&pinfotmp.probeId, etmp._epi_prb[j]);
                    pinfotmp.probechr=etmp._epi_chr[j];
                    pinfotmp.bp=etmp._epi_bp[j];
                    pinfotmp.gd=etmp._epi_gd[j];
                    strcpy2(&pinfotmp.genename, etmp._epi_gene[j]);
                    pinfotmp.orien=etmp._epi_orien[j];
                    pinfotmp.bepath.reserve(f2r);
                    pinfotmp.bepath.push_back(befNames[i]);
                    probeinfo.push_back(pinfotmp);
                    if(i && !diffidsflg) diffidsflg=true;
                } else {
                    iter=prb_map.find(pidstr);
                    if(iter!=prb_map.end())
                    {
                        probeinfo[iter->second].bepath.push_back(befNames[i]);
                    }
                }
            }
        }
        if(!diffidsflg) prb_uni=true;
        LOGPRINTF("Total %ld probes to be included from %ld epi files.\n",probeinfo.size(),f2r);
    }
    void write_epi(char* outFileName, vector<probeinfolst> &probeinfo)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".opi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        for(int i=0;i<probeinfo.size();i++)
        {
            string str=atosm(probeinfo[i].probechr)+'\t'+probeinfo[i].probeId+'\t'+atosm(probeinfo[i].bp)+'\t'+probeinfo[i].genename+'\t'+(probeinfo[i].orien=='*'?"NA":atos(probeinfo[i].orien))+'\n';
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

    void retrive_summary(string beedFileName, uint32_t &indicator,uint32_t &inum,uint32_t &pnum)
    {
        FILE* beedfptr=NULL;
        if(fopen_checked(&beedfptr,beedFileName.c_str(),"rb")) TERMINATE();
        if(fread(&indicator, sizeof(uint32_t),1, beedfptr)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", beedFileName.c_str());
            TERMINATE();
        }
        if(fread(&inum, sizeof(uint32_t),1, beedfptr)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", beedFileName.c_str());
            TERMINATE();
        }
        if(fread(&pnum, sizeof(uint32_t),1, beedfptr)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", beedFileName.c_str());
            TERMINATE();
        }
        fclose(beedfptr);
    }
    void save_beed_indi_uni(char* outFileName, vector<string> &befNames, vector<indiinfolst> &indiinfo, vector<probeinfolst> &probeinfo)
    {
        string beedfile=string(outFileName)+string(".bod");
        FILE * beedfptr;
        char inputname[FNAMESIZE];
        
        beedfptr = fopen (beedfile.c_str(), "wb");
        if (!(beedfptr)) {
            LOGPRINTF("ERROR: Failed to open file %s.\n",beedfile.c_str());
            TERMINATE();
        }

        //check the summary of bod files
        uint32_t indicator=0,inum=(uint32_t)indiinfo.size(), pnum=0;
        for(int i=0;i<befNames.size();i++)
        {
            
            string tmpname=befNames[i]+".bod";
            uint32_t tmpindr,tmpi, tmpp;
            retrive_summary( tmpname,  tmpindr, tmpi, tmpp);
            if(tmpi!=inum) {
                LOGPRINTF("ERROR: in %s : .bod file is not consistent with .oii file.\n",befNames[i].c_str());
                TERMINATE();
            } else {
                if(i==0)
                {
                    indicator=tmpindr;
                    pnum=tmpp;
                } else {
                    if(indicator!=tmpindr) {
                        LOGPRINTF("ERROR: the file %s has different file tyoe from the other files.\n",befNames[i].c_str());
                        TERMINATE();
                    } else {
                        pnum+=tmpp;
                    }
                }
            }
            
        }
        if(pnum!=probeinfo.size())
        {
            LOGPRINTF("ERROR: .bod file is not consistent with .opi file found.\n");
            TERMINATE();
        }
        //end of check.
        uint32_t valType=indicator & 0xF;
        uint32_t eType=(indicator & 0xFFF0)>>8;
        string valtypestr=getValType(valType);
        string eTypestr=getFileType(eType);
        
        if (fwrite_checked(&indicator, sizeof(uint32_t), beedfptr))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", beedfile.c_str());
            TERMINATE();
        }
        if (fwrite_checked(&inum, sizeof(uint32_t), beedfptr))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", beedfile.c_str());
            TERMINATE();
        }

        if (fwrite_checked(&pnum, sizeof(uint32_t), beedfptr))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", beedfile.c_str());
            TERMINATE();
        }

        eInfo einfo;
        string eiifile=string(probeinfo[0].bepath[0])+".oii";
        strncpy(inputname, eiifile.c_str(), sizeof(inputname));
        inputname[sizeof(inputname) - 1] = 0;
        read_eii(inputname , &einfo);
        
        char* writeBuf;
        uint64_t wBufsize=allocReserved(&writeBuf,inum*pnum*sizeof(double));
        double* dptr=(double *)writeBuf;
        map<string, int>::iterator iter;
        long sizeperprb=inum*sizeof(double);
        long prbperloop=wBufsize/sizeperprb;
        int loops=ceil(1.0*pnum/prbperloop);
        long writepnum=0;
        for(int j=0;j<loops;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/loops);
            fflush(stdout);
            memset(writeBuf,0,sizeof(char)*wBufsize);
            uint64_t numprbcurloop=prbperloop;
            if(j==loops-1) numprbcurloop=pnum-j*numprbcurloop;
            writepnum+=numprbcurloop;
            map<string,int> fcurloop;
            long fnum=0;
            vector<string> fpaths;
            vector< vector<int> > pids;
            for(int k=0;k<numprbcurloop;k++)
            {
                long curPrid=j*prbperloop+k;
                if(curPrid>=probeinfo.size()) break;
                string prbname=probeinfo[curPrid].probeId;
                if(probeinfo[curPrid].bepath.size()>1)
                {
                    LOGPRINTF("ERROR: duplicate probe %s found in one more beed files.\n",probeinfo[curPrid].probeId);
                    for(int j=0;j<probeinfo[curPrid].bepath.size();j++) LOGPRINTF("%s\n",probeinfo[curPrid].bepath[j].c_str());
                    TERMINATE();
                }
                iter=fcurloop.find(probeinfo[curPrid].bepath[0]);
                if(iter!=fcurloop.end())
                {
                    pids[iter->second].push_back(k);
                    
                } else {
                    fcurloop.insert(pair<string,int>(probeinfo[curPrid].bepath[0],fnum));
                    
                    fpaths.push_back(probeinfo[curPrid].bepath[0]);
                    vector<int> pidtmp;
                    pidtmp.reserve(numprbcurloop);
                    pidtmp.push_back(k);
                    
                    pids.push_back(pidtmp);
                    fnum++;
                }
            }
            for(int l=0;l<fpaths.size();l++) {
                
                string epifile=fpaths[l]+".opi";
                strncpy(inputname, epifile.c_str(), sizeof(inputname));
                inputname[sizeof(inputname) - 1] = 0;
                read_epi(inputname , &einfo);
                epifile=fpaths[l]+".bod";
                strncpy(inputname, epifile.c_str(), sizeof(inputname));
                inputname[sizeof(inputname) - 1] = 0;
                FILE *fptr=fopen(inputname, "rb");
                if(!fptr)
                {
                    LOGPRINTF( "ERROR: Couldn't open file %s\n", inputname);
                    TERMINATE();
                }
                for(int m=0;m<pids[l].size();m++)
                {
                    long curPrid=j*prbperloop+pids[l][m];
                    string curprb=probeinfo[curPrid].probeId;
                    iter=einfo._epi_map.find(curprb);
                    if(iter!=einfo._epi_map.end())
                    {
                        uint64_t rawpid=iter->second;
                        fseek(fptr,3*sizeof(uint32_t)+rawpid*einfo._eii_num*sizeof(double), SEEK_SET);
                        double* wptr=dptr+pids[l][m]*einfo._eii_num; //here einfo._eii_num == indiinfo.size()
                        if(fread(wptr, sizeof(double),einfo._eii_num,fptr)!=einfo._eii_num)
                        {
                            LOGPRINTF("ERROR: File %s read failed!\n", inputname);
                            TERMINATE();
                        }

                    }
                }
                fclose(fptr);
            }
            if (fwrite_checked(writeBuf, numprbcurloop*einfo._eii_num*sizeof(double), beedfptr))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedfile.c_str());
                TERMINATE();
            }
        }
        if(writepnum!=pnum)
        {
            LOGPRINTF("ERROR: wrote probe number %ld is not consistent with the total probe number %d. please report this bug.\n",writepnum, pnum);
            TERMINATE();
        }
        deallocReserved(&writeBuf, wBufsize);
        fclose(beedfptr);
        LOGPRINTF("%s of %s for %d probes and %d individuals have been saved in the binary file %s.\n",valtypestr.c_str(),eTypestr.c_str(),pnum,inum,beedfile.c_str());
    }

    void merge_beed(char* outfileName, char* befileFlistName, char* problstName, char* problst2exclde,char* genelistName,  int chr, char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename,char* probe2rm,char* indilstName,char* indilst2remove,bool beta2m,bool m2beta)
    {
        vector<string> befNames;
        vector<indiinfolst> indiinfo;
        vector<probeinfolst> probeinfo;
        read_beflist(befNames, befileFlistName);
        if(befNames.size()==0){
            LOGPRINTF("No file names included from %s ...\n", befileFlistName);
            TERMINATE();
        }
        bool indi_uni=false;
        combine_eii(indiinfo, befNames, indi_uni);
        if(indiinfo.size()==0)
        {
            LOGPRINTF("ERROR: No individuals to be included!\n");
            TERMINATE();
        }
        if(!indi_uni)
        {
            //individuals in each sub-eii file should be in consistent order.
            indiinfolst* eiiptr=&indiinfo[0];
            qsort(eiiptr,indiinfo.size(),sizeof(indiinfolst),comp_eii);
        }
        
        bool prb_uni=false;
        combine_epi(probeinfo, befNames, prb_uni);
        if(probeinfo.size()==0)
        {
            LOGPRINTF("ERROR: No probes to be included!\n");
            TERMINATE();
        }
        if(!prb_uni)
        {
            //probes in each sub-epi file should be in consistent order.
            probeinfolst* epiptr=&probeinfo[0];
            qsort(epiptr,probeinfo.size(),sizeof(probeinfolst),comp_epi);
        }
        
        
        
        printf("\nGenerating epi file...\n");
        write_epi(outfileName, probeinfo);
        printf("\nGenerating eii file...\n");
        write_eii(outfileName, indiinfo);
        printf("\nGenerating beed file...\n");
        map<string, int> eii_map;
        if(!indi_uni)
            for(int i=0;i<indiinfo.size();i++) eii_map.insert(pair<string,int>(string(indiinfo[i].fid)+":"+string(indiinfo[i].iid),i));
        
        if(indi_uni){
            
            save_beed_indi_uni(outfileName, befNames, indiinfo, probeinfo);
        }
        else if(prb_uni) {
            LOGPRINTF("Release soon!\n");
            TERMINATE();
        }
        else {
            LOGPRINTF("Release soon!\n");
            TERMINATE();
        }

        free_probelist(probeinfo);
        free_indilist(indiinfo);
        

    }
    void make_beed(char* outFileName, char* efileName, char* befileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, bool adjprb, char* covfileName,char* qcovfileName)
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
            LOGPRINTF("Error: due to --adj-probe is specified. please input the Gene discrete / Continuous covariates by the option --covar or --qcovar.\n");
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
        if(adjprb) adjprobe(&einfo);
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
    
     void make_erm(char* outFileName, char* efileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,bool transposed, int efileType,bool no_fid_flag,int valueType)
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
    int comp_assoc(const void *a,const void *b){ return ((*(ASSOCRLT *)a).PVAL < (*(ASSOCRLT *)b).PVAL)?1:-1; } //decend
    
    void mlma(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno,bool erm_bin_flag, int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* subtract_erm_file, bool m_erm_flag, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_no_adj_covar, double percentage_out)
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
        
        if (covfileName == NULL && qcovfileName == NULL) mlma_no_adj_covar=false;
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo); // eii and epi are updated in it.
        
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
            make_erm(&einfo);
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
                    make_erm( &einfo);
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
            LOGPRINTF("ERROR: no individual is in common in the input files.\n");
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
            einfo._grm.resize(0,0);

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
                    einfo._grm.resize(0,0);
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
            einfo._var_name.push_back("V(G" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(G" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        einfo._within_family=within_family;
        if(within_family) detect_family(&einfo, _A);
        
        // run REML algorithm
        LOGPRINTF("\nPerforming MLM association analyses %s ...\n",(subtract_erm_file?"":" (including the candidate probes)"));
        MatrixXd _Vi;
        reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
        einfo._P.resize(0,0);
        _A.clear();
        unsigned long m=einfo._epi_include.size();
       
        VectorXd y_buf=_y;
        if(!mlma_no_adj_covar) y_buf=_y.array()-(_X*einfo._b).array(); // adjust phenotype for covariates
       
        VectorXd beta, se, pval;
        if(mlma_no_adj_covar) mlma_calcu_stat_covar(y_buf, &einfo, _X_c, _Vi, _X, beta, se, pval);
        else mlma_calcu_stat(y_buf, &einfo, _Vi, beta, se, pval);
        
        
        string filename=string(outFileName)+".mlma";
        LOGPRINTF("\nSaving the results of the mixed linear model association analyses of %ld probes to %s ...\n",m,filename.c_str());
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
    
    void output_simu_par(char* outFileName, vector<string> &qtl_name, vector<double> &qtl_eff)
    {
        int i = 0;
        string out_parfile = string(outFileName) + ".par";
        ofstream out_par(out_parfile.c_str());
        if (!out_par) throw ("Error: can not open par file [" + out_parfile + "] to write!");
        out_par << "QTL\tEffect" << endl;
        for (i = 0; i < qtl_eff.size(); i++) out_par << qtl_name[i] << "\t" << qtl_eff[i] << endl;
        out_par.close();
        cout << "Simulated QTL effect(s) have been saved in [" + out_parfile + "]." << endl;
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
            for (int i = einfo._eii_num - 1; i >= 0; i--) {
                string str=atos(eval(i))+'\n';
                if(fputs_checked(str.c_str(),efile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", evalName.c_str());
                    TERMINATE();
                }
                
            }
            fclose(efile);
            LOGPRINTF("Eigenvalues of %d  individuals have been saved in %s.\n",einfo._eii_num ,evalName.c_str());
            
            evalName = string(outFileName) + ".eigenvec";
            if(fopen_checked(&efile, evalName.c_str(),"w")) TERMINATE();
            
            if (out_pc_num > einfo._eii_num) out_pc_num = einfo._eii_num;
            for (int i = 0; i < einfo._eii_num; i++) {
                string str=einfo._eii_fid[i]+'\t'+einfo._eii_iid[i];
                
                for (int j = einfo._eii_num - 1; j >= (einfo._eii_num - out_pc_num); j--) str += '\t'+atos(evec(i, j));
                str+='\n';
                if(fputs_checked(str.c_str(),efile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", evalName.c_str());
                    TERMINATE();
                }
            }
            fclose(efile);
            LOGPRINTF("The first %d eigenvectors of %d individuals have been saved in %s.\n",out_pc_num ,einfo._eii_num, evalName.c_str());
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
            LOGPRINTF("The first %d eigenvectors of %d individuals have been saved in %s.\n",out_pc_num ,einfo._eii_num, evalName.c_str());
            
            /***/
        }
       
    }
    
    void mlma_loco(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, int MaxIter, char* priors, char* priors_var, bool no_constrain, bool no_adj_covar,int reml_mtd,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, int autosome_num, double percentage_out)
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
            LOGPRINTF("ERROR: no individual is in common in the input files.\n");
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
        einfo._var_name.push_back("V(G)");
        einfo._hsq_name.push_back("V(G)/Vp");
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
            if(!no_adj_covar) y_buf=_y.array()-(_X*einfo._b).array(); // adjust phenotype for covariates
           
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
        for(int c1=0; c1<chrs.size(); c1++){
            delete[] (grm_chrs[c1]);
            delete[] (geno_chrs[c1]);
        }
        
        string filename=string(outFileName)+".loco.mlma";
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
    void manipulate_grm(eInfo* einfo, char* grm_file, char* indilstName, char* indilst2remove, char* sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool dont_read_N)
    {
        int i = 0, j = 0;
        
        vector<string> grm_id;
        if (merge_grm_flag) merge_grm(einfo, grm_file);
        else read_grm(einfo,grm_file, grm_id, true, false, dont_read_N,true);
        
        if (indilstName!=NULL) keep_indi(einfo, indilstName);
        if (indilst2remove!=NULL) remove_indi(einfo, indilst2remove);
        if (grm_cutoff>-1.0) rm_cor_indi(einfo, grm_cutoff);
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
    void pca(char* outFileName, char* grm_file, char* indilstName, char* indilst2remove, double grm_cutoff, bool merge_grm_flag, int out_pc_num)
    {
        eInfo einfo;
        manipulate_grm(&einfo, grm_file, indilstName, indilst2remove, NULL, grm_cutoff, -2.0, -2, merge_grm_flag, true);
        einfo._grm_N.resize(0, 0);
        int i = 0, j = 0, n = einfo._eii_include.size();
        cout << "\nPerforming principal component analysis ..." << endl;
        
        SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
        MatrixXd evec = (eigensolver.eigenvectors());
        VectorXd eval = eigensolver.eigenvalues();
        
        string eval_file = string(outFileName) + ".eigenval";
        ofstream o_eval(eval_file.c_str());
        if (!o_eval) throw ("Error: can not open the file [" + eval_file + "] to read.");
        for (i = n - 1; i >= 0; i--) o_eval << eval(i) << endl;
        o_eval.close();
        cout << "Eigenvalues of " << n << " individuals have been saved in [" + eval_file + "]." << endl;
        string evec_file = string(outFileName) + ".eigenvec";
        ofstream o_evec(evec_file.c_str());
        if (!o_evec) throw ("Error: can not open the file [" + evec_file + "] to read.");
        if (out_pc_num > n) out_pc_num = n;
        for (i = 0; i < n; i++) {
            o_evec << einfo._eii_fid[einfo._eii_include[i]] << " " << einfo._eii_iid[einfo._eii_include[i]];
            for (j = n - 1; j >= (n - out_pc_num); j--) o_evec << " " << evec(i, j);
            o_evec << endl;
        }
        o_evec.close();
        cout << "The first " << out_pc_num << " eigenvectors of " << n << " individuals have been saved in [" + evec_file + "]." << endl;
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

        LOGPRINTF("\nPerforming association analysis...\n");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            fflush(stdout);
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
                if(abs(phval+9)>1e-8)
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
                if(abs(phval+9)>1e-8)
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
            LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
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
    void  adjprobe(eInfo* einfo)
    {
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(einfo, E_float, qE_float,_X);
        LOGPRINTF("\nAdjusting probes...\n");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            fflush(stdout);
            
            vector<double> cvec;
            vector<double> xvec;
            vector<int> NMISS;
            MatrixXd X=_X;
            double nonmiss=0.0;
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                if(val<1e9){
                    NMISS.push_back(einfo->_eii_include[j]);
                    xvec.push_back(val);
                    nonmiss+=1.0;
                } else {
                    removeRow(X, j);
                }
            }
            
            VectorXd x(xvec.size());
            for(int j=0;j<xvec.size();j++) x(j)=xvec[j];
            
            if(x.size()!=X.rows() || x.size()<1) {
                LOGPRINTF("Error: The row number of C and the length of y do not match.\n");
                TERMINATE();
            }
            
            MatrixXd XtX_i;
            //XtX_i=(X.transpose()*X).inverse(); //DO NOT USE IT, IT WOULD GIVE VERY VERY WRONG RESULT WHEN THE MATRIX IS NOT INVERTIBLE
            XtX_i=X.transpose()*X;
            bool determinant_zero=false;
            inverse_V(XtX_i, determinant_zero);
            VectorXd b_hat=XtX_i*X.transpose()*x;
            VectorXd residual=(x-X*b_hat);
            for(int j=0;j<residual.size();j++) einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+NMISS[j]]=residual(j);
        }
         LOGPRINTF("%ld probes have been adjusted.\n",einfo->_epi_include.size());
    }

    void  testLinear(vector<ASSOCRLT> &assoc_rlts,char* outfileName, eInfo* einfo)
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
            outstr="probeChr\tProbeID\tProbe_bp\tBETA\tSE\tP\tNMISS\n";
            if(fputs_checked(outstr.c_str(),assoc))
            {
                LOGPRINTF("ERROR: error in writing file %s .\n", assocfile.c_str());
                TERMINATE();
            }
        } else {
            assoc_rlts.clear();
        }
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(einfo, E_float, qE_float,_X);
        LOGPRINTF("\nPerforming linear regression analysis...\n");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            printf("%3.0f%%\r", 100.0*i/einfo->_epi_include.size());
            fflush(stdout);

            vector<double> yvec;
            vector<double> cvec;
            vector<double> xvec;
            MatrixXd X=_X;
            double nonmiss=0.0;
            int chr=einfo->_epi_chr[einfo->_epi_include[i]];
            string prbid=einfo->_epi_prb[einfo->_epi_include[i]];
            string gene=einfo->_epi_gene[einfo->_epi_include[i]];
            int BP=einfo->_epi_bp[einfo->_epi_include[i]];
            char oren=einfo->_epi_orien[einfo->_epi_include[i]];
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                double phval=einfo->_eii_pheno[einfo->_eii_include[j]];
                if(abs(phval+9)>1e-8) // no phval missing here, cos individuals with missing phenotype were removed when read phenotype file.
                {
                    double val=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9){
                        yvec.push_back(phval);
                        xvec.push_back(val);
                        nonmiss+=1.0;
                    } else {
                        removeRow(X, j);
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
            lin(y, X, x, rst);
            
            ASSOCRLT currlt;
            if(assoc) {
                string outstr = atos(chr) + '\t' + prbid + '\t' + atos(BP) + '\t' + atos(rst[0]) + '\t' + atos(rst[1]) + '\t' + dtos(rst[2]) + '\t' + atos(nonmiss) +'\n';
                if(fputs_checked(outstr.c_str(),assoc))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", assocfile.c_str());
                    TERMINATE();
                }
                write_count++;
            } else {
                currlt.BETA=rst[0];
                currlt.SE=rst[1];
                currlt.PVAL=rst[2];
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
            LOGPRINTF("Results of %ld probes have been saved in file %s.\n",write_count,assocfile.c_str());
            fclose(assoc);
        } else {
            LOGPRINTF("Results of %ld probes have been returned.\n",assoc_rlts.size());
        }
        
        
    }

    void  linear(char* outfileName, char* befileName, char* problstName, char* problst2exclde, char* genelistName,  int chr,char* prbname,  char* fromprbname,  char* toprbname, int prbWind, int fromprbkb,  int toprbkb, bool prbwindFlag,  char* genename, char* probe2exclde, char* indilstName, char* indilst2remove,char* phenofileName,char* mpheno,  char* covfileName, char* qcovfileName, double std_thresh ,double upperBeta,double lowerBeta ,int tsk_ttl,int tsk_id)
    {
        eInfo einfo;
        init_einfo(&einfo);
      
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
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
        
        vector<ASSOCRLT> assoc_rlts;
        testLinear(assoc_rlts,outfileName, &einfo);
        
        
    }
    /*
    void LogisticfitLM()
    {
        vector<double> coef;
        
        coef.resize(np);
        sizeMatrix(S,np,np);
        
        if (np==0 || nind==0 || ! all_valid )
            return;
        
        if (par::verbose)
        {
            for (int i=0; i<nind; i++)
            {
                cout << i << "\t"
                << Y[i] << "\t";
                for (int j=0; j<np; j++)
                    cout << X[i][j] << "\t";
                cout << "\n";
            }
        }
        
        ///////////////////////////////////////
        // Newton-Raphson to fit logistic model
        
        bool converge = false;
        int it = 0;
        
        while ( ! converge && it < 20 )
        {
            
            // Determine p and V
            for (int i=0; i<nind; i++)
            {
                double t = 0;
                for (int j=0; j<np; j++)
                    t += coef[j] * X[i][j];
                p[i] = 1/(1+exp(-t));
                V[i] = p[i] * (1-p[i]);
            }
            
            // Update coefficients
            // b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p )
            
            matrix_t T;
            sizeMatrix(T,np,np);
            
            for (int j=0; j<np; j++)
                for (int k=j; k<np; k++)
                {
                    double sum = 0;
                    for (int i=0; i<nind; i++)
                        sum += X[i][j] * V[i] * X[i][k] ;
                    T[j][k] = T[k][j] = sum;
                }
            
            bool flag = true;
            T = svd_inverse(T,flag);
            if ( ! flag )
            {
                all_valid = false;
                return;
            }
            
            matrix_t T2;
            // Resize and set elements to 0
            sizeMatrix(T2,np,nind);
            
            // note implicit transpose of X
            for (int i=0; i<np; i++)
                for (int j=0; j<nind; j++)
                    for (int k=0; k<np; k++)
                        T2[i][j] += T[i][k] * X[j][k];
            
            vector_t t3(nind);
            for (int i=0; i<nind; i++)
                t3[i] = Y[i] - p[i];
            
            vector_t ncoef(np);
            for (int j=0; j<np; j++)
                for (int i=0; i<nind; i++)
                    ncoef[j] += T2[j][i] * t3[i];
            
            // Update coefficients, and check for
            // convergence
            double delta = 0;
            for (int j=0; j<np; j++) 	
            {
                delta += abs(ncoef[j]);
                coef[j] += ncoef[j];
            }
            
            if ( delta < 1e-6 )
                converge = true;
            
            // Next iteration
            it++;
        }
        
        
        /////////////////////////////////////////
        // Obtain covariance matrix of estimates
        
        // S <- solve( t(X) %*% V %*% X )    
        
        // Transpose X and multiple by diagonal V
        matrix_t Xt;
        sizeMatrix(Xt, np, nind);
        for (int i=0; i<nind; i++)
            for (int j=0; j<np; j++) 
                Xt[j][i] = X[i][j] * V[i];
        
        multMatrix(Xt,X,S);  
        bool flag = true;
        S = svd_inverse(S,flag);     
        if ( ! flag ) 
        {
            all_valid = false;
            return;
        }
        if ( cluster ) 
            HuberWhite();
        
        if (par::verbose)
        {
            cout << "beta\n";
            display(coef);
            cout << "Sigma\n";
            display(S);
            cout << "\n";
        }
    }
    */
    int read_QTL_file(eInfo* einfo,string qtl_file, vector<string> &qtl_name, vector<double> &qtl_eff, vector<int> &have_eff)
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
        int icount = 0;
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
        cout << "Heritability " << (cc_flag ? "of liability = " : " = ") << hsq << " (Default = 0.1)" << endl;
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
        int qtl_num = read_QTL_file(&einfo,sigCpG_file, qtl_name, qtl_eff, have_eff);
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
    
    
    void fit_reml(char* outFileName, char* phenofileName,char* mpheno,bool erm_bin_flag, bool grm_bin_flag,int erm_alg, char* covfileName,char* qcovfileName, char* erm_file, char* grm_file, bool m_erm_flag, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_no_adj_covar, bool pred_rand_eff, bool est_fix_eff,bool no_lrt,double prevalence, bool mlmassoc,vector<int> drop, char* indilstName, char* indilst2remove, char* sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen,bool prt_residiual)
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
        if(prt_residiual) est_fix_eff=true;
        
        eInfo einfo;
        init_einfo(&einfo);
        einfo._reml_mtd=reml_mtd;
        einfo._reml_max_iter=MaxIter;
        einfo._reml_fixed_var=reml_fixed_var_flag;
        einfo._reml_force_inv=reml_force_inv_fac_flag;
        einfo._reml_force_converge=reml_force_converge_flag;
        einfo._reml_no_converge=reml_no_converge_flag;

        
        vector<string> grm_id;
        vector<string> erm_files;
        
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
        
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        else {
            LOGPRINTF("ERROR: no phenoytpe file inputed.\n");
            TERMINATE();
        }
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        if (covfileName == NULL && qcovfileName == NULL) mlma_no_adj_covar=false;
        if (indilstName!=NULL) keep_indi(&einfo, indilstName);
        if (indilst2remove!=NULL) remove_indi(&einfo, indilst2remove);
        if(erm_file!=NULL) {
            if (grm_cutoff>-1.0) rm_cor_indi(&einfo,grm_cutoff);
            if (sex_file!=NULL) update_sex(&einfo,sex_file);
            if (adj_grm_fac>-1.0) adj_grm(&einfo,adj_grm_fac);
            if (dosage_compen>-1) dc(&einfo,dosage_compen);
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
            LOGPRINTF("ERROR: no individual is in common in the input files.\n");
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
#pragma omp parallel for
                    for (int j = 0; j < _n; j++) {
                        for (int k = 0; k <= j; k++) {
                            if (kp[j] >= kp[k]) (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[j], kp[k]);
                            else (_A[i])(k, j) = (_A[i])(j, k) = einfo._grm(kp[k], kp[j]);
                        }
                    }
                }
            } else if(erm_file!=NULL){
                for(int i=0; i < 1 + 1; i++) einfo._r_indx.push_back(i);
                if (!no_lrt) drop_comp(&einfo,drop);
                _A.resize(einfo._r_indx.size());
                match(uni_id, grm_id, kp);
                (_A[0]).resize(_n, _n);
#pragma omp parallel for
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
            einfo._var_name.push_back("V(G" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(G" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        einfo._within_family=within_family;
        if(within_family) detect_family(&einfo, _A);
        
        // run REML algorithm
        MatrixXd _Vi;
        reml(&einfo, pred_rand_eff, est_fix_eff, reml_priors, reml_priors_var, prevalence, -2.0, no_constrain, no_lrt, mlmassoc, _X_c, _X,_y,_A,_Vi,outFileName);
        
        
        if(prt_residiual) {
            VectorXd y_buf=_y;
            if(!mlma_no_adj_covar) y_buf=_y.array()-(_X*einfo._b).array(); // adjust phenotype for covariates
            
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
    void get_ld_blk_pnt(eInfo* einfo, vector<int> &brk_pnt1, vector<int> &brk_pnt2, vector<int> &brk_pnt3, int wind_bp)
    {
        unsigned long i = 0, j = 0, k = 0, m = einfo->_epi_include.size();
        
        brk_pnt1.clear();
        brk_pnt1.push_back(0);
        bool chr_start = true;
        for (i = 1, j = 0; i < m; i++) {
            if (i == (m - 1)){
                if(chr_start
                   || ((einfo->_epi_bp[einfo->_epi_include[i]] - einfo->_epi_bp[einfo->_epi_include[brk_pnt1[j]]] > 0.5 * wind_bp))) brk_pnt1.push_back(m - 1);
                else brk_pnt1[j - 1] = brk_pnt1[j] = m - 1;
            }
            else if (einfo->_epi_chr[einfo->_epi_include[i]] != einfo->_epi_chr[einfo->_epi_include[brk_pnt1[j]]] || einfo->_epi_bp[einfo->_epi_include[i]] - einfo->_epi_bp[einfo->_epi_include[i-1]] > 1e6) {
                if(chr_start
                   || ((einfo->_epi_bp[einfo->_epi_include[i-1]] - einfo->_epi_bp[einfo->_epi_include[brk_pnt1[j]]] > 0.5 * wind_bp))){
                       brk_pnt1.push_back(i - 1);
                       j++;
                       brk_pnt1.push_back(i);
                       j++;
                   }
                else{
                    brk_pnt1[j - 1] = i - 1;
                    brk_pnt1[j] = i;
                }
                chr_start = true;
            }
            else if ((einfo->_epi_bp[einfo->_epi_include[i]] - einfo->_epi_bp[einfo->_epi_include[brk_pnt1[j]]] > wind_bp)) {
                chr_start = false;
                brk_pnt1.push_back(i - 1);
                j++;
                brk_pnt1.push_back(i);
                j++;
            }
        }
        stable_sort(brk_pnt1.begin(), brk_pnt1.end());
        brk_pnt1.erase(unique(brk_pnt1.begin(), brk_pnt1.end()), brk_pnt1.end());
        
        brk_pnt2.clear();
        brk_pnt3.clear();
        for (i = 1; i < brk_pnt1.size() && brk_pnt1.size() > 2; i++) {
            if ((einfo->_epi_chr[einfo->_epi_include[brk_pnt1[i - 1]]] == einfo->_epi_chr[einfo->_epi_include[brk_pnt1[i]]]) && (brk_pnt1[i] - brk_pnt1[i - 1] > 1)) {
                int i_buf = (brk_pnt1[i - 1] + brk_pnt1[i]) / 2;
                brk_pnt2.push_back(i_buf);
                brk_pnt2.push_back(i_buf + 1);
                brk_pnt3.push_back(brk_pnt1[i]);
                brk_pnt3.push_back(brk_pnt1[i]);
            }
        }
    }
    /*
    void calcu_ld_blk(eInfo* einfo, vector<int> &brk_pnt, vector<int> &brk_pnt3, VectorXd &mean_rsq, VectorXd &snp_num, VectorXd &max_rsq, bool second, double rsq_cutoff)
    {
        int i = 0, j = 0, k = 0, s1 = 0, s2 = 0, n = einfo->_eii_include.size(), m = einfo->_epi_include.size(), size = 0, size_limit = 10000;
        
        for (i = 0; i < brk_pnt.size() - 1; i++) {
            if (einfo->_epi_chr[einfo->_epi_include[brk_pnt[i]]] != einfo->_epi_chr[einfo->_epi_include[brk_pnt[i + 1]]]) continue;
            size = brk_pnt[i + 1] - brk_pnt[i] + 1;
            if (size < 3) continue;
            if (second) {
                s1 = brk_pnt3[i] - brk_pnt[i];
                s2 = s1 + 1;
            }
            else {
                s1 = 0;
                s2 = size - 1;
            }
            
            VectorXd rsq_size(size), mean_rsq_sub(size), max_rsq_sub = VectorXd::Constant(size, -1.0);
            
            // make genotype matrix
            vector<int> snp_indx(size);
            for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++) snp_indx[k] = j;
            MatrixXf X_sub;
             make_XMat_subset(X_sub, snp_indx, true);
            VectorXd ssx_sqrt_i_sub(size);
            for (j = 0; j < size; j++){
                ssx_sqrt_i_sub[j] = X_sub.col(j).squaredNorm();
                if (ssx_sqrt_i_sub[j] < 1.0e-30) ssx_sqrt_i_sub[j] = 0.0;
                else ssx_sqrt_i_sub[j] = 1.0 / sqrt(ssx_sqrt_i_sub[j]);
            }
            
            if (size > size_limit) calcu_ld_blk_split(size, size_limit, X_sub, ssx_sqrt_i_sub, rsq_cutoff, rsq_size, mean_rsq_sub, max_rsq_sub, s1, s2, second);
            else {
                MatrixXf rsq_sub = X_sub.transpose() * X_sub;
#pragma omp parallel for private(k)
                for (j = 0; j < size; j++) {
                    rsq_size[j] = 0.0;
                    mean_rsq_sub[j] = 0.0;
                    for (k = 0; k < size; k++) {
                        if (second) {
                            if (j <= s1 && k <= s1) continue;
                            if (j >= s2 && k >= s2) continue;
                        }
                        if (k == j) continue;
                        rsq_sub(j,k) *= (ssx_sqrt_i_sub[j] * ssx_sqrt_i_sub[k]);
                        rsq_sub(j,k) = rsq_sub(j,k) * rsq_sub(j,k);
                        if (rsq_sub(j,k) >= rsq_cutoff) {
                            if(_ldscore_adj) mean_rsq_sub[j] += rsq_sub(j,k) - (1.0 - rsq_sub(j,k)) / (n - 2.0);
                            else mean_rsq_sub[j] += rsq_sub(j,k);
                            rsq_size[j] += 1.0;
                        }
                        if (rsq_sub(j,k) > max_rsq_sub[j]) max_rsq_sub[j] = rsq_sub(j,k);
                    }
                    if (rsq_size[j] > 0.0) mean_rsq_sub[j] /= rsq_size[j];
                }
            }
            
            for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
                if (second) {
                    if (rsq_size[j] > 0.0) {
                        mean_rsq[k] = (mean_rsq[k] * snp_num[k] + mean_rsq_sub[j] * rsq_size[j]) / (snp_num[k] + rsq_size[j]);
                        snp_num[k] = (snp_num[k] + rsq_size[j]);
                        if(max_rsq[k] < max_rsq_sub[j]) max_rsq[k] = max_rsq_sub[j];
                    }
                }
                else {
                    mean_rsq[k] = mean_rsq_sub[j];
                    snp_num[k] = rsq_size[j];
                    max_rsq[k] = max_rsq_sub[j];
                }
            }
        }
    }
    */
    
    void calcu_mean_rsq(char* outFileName, char* efileName, char* befileName, bool transposed, int efileType,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, bool no_fid_flag,int valueType,bool beta2m,bool m2beta, double std_thresh,double upperBeta,double lowerBeta,char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh, int filter_mth, double mssratio_prob, int wind_size, double rsq_cutoff, int autosome_num, bool _ldscore_adj)
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
        int i = 0, m = einfo._epi_include.size();
        
        LOGPRINTF("\nCalculating correlation score for probes (block size of %d Kb with an overlap of %d Kb between blocks); correlation rsq threshold = %f ) ... \n", wind_size / 1000,wind_size/2000, rsq_cutoff );
        if(_ldscore_adj) cout << "correlation rsq will be adjusted for chance correlation, i.e. rsq_adj = rsq - (1 - rsq) / (n -2)." << endl;
        vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
        get_ld_blk_pnt(&einfo, brk_pnt1, brk_pnt2, brk_pnt3, wind_size);
        /*
        VectorXd mean_rsq = VectorXd::Zero(m), snp_num = VectorXd::Zero(m), max_rsq = VectorXd::Zero(m);
        calcu_ld_blk(brk_pnt1, brk_pnt3, mean_rsq, snp_num, max_rsq, false, rsq_cutoff, dominance_flag);
        if (brk_pnt2.size() > 1) calcu_ld_blk(brk_pnt2, brk_pnt3, mean_rsq, snp_num, max_rsq, true, rsq_cutoff, dominance_flag);
        
        string mrsq_file = "";
        if(dominance_flag) mrsq_file = _out + ".d.score.ld";
        else mrsq_file = _out + ".score.ld";
        ofstream o_mrsq(mrsq_file.data());
        o_mrsq<<"SNP chr bp MAF mean_rsq snp_num max_rsq ldscore"<<endl;
        double ldscore = 0.0;
        for (i = 0; i < m; i++){
            o_mrsq << _snp_name[_include[i]] << " " << _chr[_include[i]] << " " << _bp[_include[i]] << " ";
            double MAF = 0.5 * _mu[_include[i]];
            if(MAF > 0.5) MAF = 1.0 - MAF;
            ldscore = 1.0 + mean_rsq[i] * snp_num[i];
            o_mrsq << MAF << " " << mean_rsq[i] << " " << snp_num[i] << " " << max_rsq[i] << " " << ldscore << "\n";
        }
        o_mrsq << endl;
        cout << "LD score for " << m << " SNPs have been saved in the file [" + mrsq_file + "]." << endl;
         */
    }
    
    
}

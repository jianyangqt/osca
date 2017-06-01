//
//  l2_efile.cpp
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l2_efile.h"

namespace EFILE {
    
    // read eii, epi, pheno,beed
    
   string getFileType(uint32_t tid)
    {
        return (tid==0?"Gene expression":"Methylation");
    }
    string getValType(uint32_t tid)
    {
        if(tid==0) return("Beta-values");
        else if(tid==1) return("M-values");
        else return("Expression values");
    }
    void  read_efile(char* eFileName, eInfo* einfo, uint32_t filetype, bool no_fid_flag,int valueType){
        
        FILE* efile=NULL;
        einfo->_eType=filetype;
        einfo->_valType=valueType;
        string typestr=getFileType(einfo->_eType);
        vector<string> strlist;
        uint32_t line_idx = 0;
        char* readBuf;
        string msg="Reading Buffer";
        if(!allocReserved(&readBuf, MAXLINEBUFLEN,msg)) exit(EXIT_FAILURE);
        if(fopen_checked(&efile, eFileName,"r")) exit(EXIT_FAILURE);
        LOGPRINTF("Reading %s data from %s ...\n", typestr.c_str(),eFileName);
        if(fgets(readBuf, MAXLINEBUFLEN, efile))
        {
            if(no_fid_flag) split_str(readBuf,strlist,1);
            else split_str(readBuf,strlist,2);
            strlist.swap(einfo->_epi_prb);
            if(!sudoAllocReserved( einfo->_epi_prb.size()*(sizeof(int)*5+sizeof(string)*2+sizeof(char)), "probe infomation")) exit(EXIT_FAILURE);
            einfo->_epi_include.resize(einfo->_epi_prb.size());
            einfo->_epi_chr.resize(einfo->_epi_prb.size());
            einfo->_epi_bp.resize(einfo->_epi_prb.size());
            einfo->_epi_gd.resize(einfo->_epi_prb.size());
            einfo->_epi_gene.resize(einfo->_epi_prb.size());
            einfo->_epi_orien.resize(einfo->_epi_prb.size());
            for(int i=0;i<einfo->_epi_prb.size();i++)
            {
                einfo->_epi_include[i]=i;
                einfo->_epi_chr[i]=-9;
                einfo->_epi_bp[i]=-9;
                einfo->_epi_gd[i]=-9;
                einfo->_epi_gene[i]="NA";
                einfo->_epi_orien[i]='*';
                einfo->_epi_map.insert(pair<string,int>(einfo->_epi_prb[i],i));
                if(einfo->_epi_map.size()==i)
                {
                    LOGPRINTF("ERROR: Duplicate probe %s .\n", einfo->_epi_prb[i].c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }
        while(fgets(readBuf, MAXLINEBUFLEN, efile)) line_idx++;
        einfo->_eii_num=line_idx;
        einfo->_epi_num=(uint32_t)einfo->_epi_prb.size();
        einfo->_eii_pheno_num=1;
        msg="Expression Values";
        uint64_t sizeneed=einfo->_epi_num*einfo->_eii_num*sizeof(double)+einfo->_eii_num*(sizeof(string)*5+sizeof(int)*3+sizeof(double));
        if(!sudoAllocReserved( sizeneed, msg)) exit(EXIT_FAILURE);
        einfo->_val.resize(einfo->_epi_num*einfo->_eii_num);
        einfo->_eii_fid.resize(einfo->_eii_num);
        einfo->_eii_iid.resize(einfo->_eii_num);
        einfo->_eii_include.resize(einfo->_eii_num);
        einfo->_eii_fa_id.resize(einfo->_eii_num);
        einfo->_eii_mo_id.resize(einfo->_eii_num);
        einfo->_eii_sex.resize(einfo->_eii_num);
        einfo->_eii_pheno.resize(einfo->_eii_num);
        
        rewind(efile);
        fgets(readBuf, MAXLINEBUFLEN, efile); //head
        line_idx=0;
        memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
        while(fgets(readBuf, MAXLINEBUFLEN, efile))
        {
            split_str(readBuf,strlist,0);
            if(no_fid_flag){
                if(strlist.size()!=1+einfo->_epi_num)
                {
                    LOGPRINTF("ERROR: Inconsistent element number found in line %u with IID %s .\n", line_idx+1, strlist[0].c_str());
                    exit(EXIT_FAILURE);
                }
                einfo->_eii_fid[line_idx]=strlist[0].c_str();
                einfo->_eii_iid[line_idx]=strlist[0].c_str();
            }else {
                if(strlist.size()!=2+einfo->_epi_num)
                {
                    LOGPRINTF("ERROR: Inconsistent element number found in line %u with FID: %s and IID %s .\n", line_idx+1, strlist[0].c_str(),strlist[1].c_str());
                    exit(EXIT_FAILURE);
                }
                einfo->_eii_fid[line_idx]=strlist[0].c_str();
                einfo->_eii_iid[line_idx]=strlist[1].c_str();
            }
            einfo->_eii_include[line_idx]=line_idx;
            einfo->_eii_fa_id[line_idx]="0";
            einfo->_eii_mo_id[line_idx]="0";
            einfo->_eii_sex[line_idx]=-9;
            einfo->_eii_pheno[line_idx]=-9;
            einfo->_eii_map.insert(pair<string,int>(strlist[0]+":"+strlist[1],line_idx));
            if(einfo->_eii_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate individual with FID: %s and IID %s.\n", strlist[0].c_str(),strlist[1].c_str());
                exit(EXIT_FAILURE);
            }
            int valStartPos=2;
            if(no_fid_flag) valStartPos=1;
            for(int i=valStartPos;i<strlist.size();i++)
            {
                long pid=i-valStartPos;
                if(strlist[i]=="na" ||strlist[i]=="NA" || strlist[i]=="." || strlist[i]=="-" || strlist[i]=="-9") einfo->_val[pid*einfo->_eii_num+line_idx]=1e10; // -9 as miss value should be questioned. 1e10 presents as miss value also should be questioned.
                else {
                    double valtmp=atof(strlist[i].c_str());
                    if(!valueType && valtmp<0) {
                        LOGPRINTF("ERROR: Beta value of probe %s should be positive.\n", einfo->_epi_prb[line_idx].c_str());
                        exit(EXIT_FAILURE);
                    } else {
                        einfo->_val[pid*einfo->_eii_num+line_idx]=valtmp;
                    }
                }
            }
            line_idx++;
            memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
        }
        fclose(efile);
        deallocReserved(&readBuf, MAXLINEBUFLEN);
        LOGPRINTF("%s data for %u probes of %u individuals have been included from the file %s.\n",typestr.c_str(), einfo->_epi_num, einfo->_eii_num, eFileName);
    }
    
    void  read_efile_t(char* eFileName, eInfo* einfo, uint32_t filetype, bool no_fid_flag,int valueType){
        
        FILE* efile=NULL;
         einfo->_eType=filetype;
        einfo->_valType=valueType;
        vector<string> strlist;
        uint32_t line_idx = 0;
        string typestr=getFileType(einfo->_eType);
        char* readBuf;
        string msg="Reading Buffer";
        if(!allocReserved(&readBuf, MAXLINEBUFLEN,msg)) exit(EXIT_FAILURE);
        if(fopen_checked(&efile, eFileName,"r")) exit(EXIT_FAILURE);
        LOGPRINTF("Reading %s data from transposed file %s ...\n", typestr.c_str(), eFileName);
        if(no_fid_flag)
        {
            if(fgets(readBuf, MAXLINEBUFLEN, efile))
            {
                split_str(readBuf,strlist,1);
                strlist.swap(einfo->_eii_iid);
                einfo->_eii_fid=einfo->_eii_iid;
            }
        } else {
            if(fgets(readBuf, MAXLINEBUFLEN, efile))
            {
                split_str(readBuf,strlist,1);
                strlist.swap(einfo->_eii_fid);
            }
            memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
            if(fgets(readBuf, MAXLINEBUFLEN, efile))
            {
                split_str(readBuf,strlist,1);
                strlist.swap(einfo->_eii_iid);
            }
        }
      
        if(!sudoAllocReserved( einfo->_eii_iid.size()*(sizeof(int)*3+sizeof(string)*3+sizeof(double)), "individual information")) exit(EXIT_FAILURE);
        einfo->_eii_include.resize(einfo->_eii_iid.size());
        einfo->_eii_fa_id.resize(einfo->_eii_iid.size());
        einfo->_eii_mo_id.resize(einfo->_eii_iid.size());
        einfo->_eii_sex.resize(einfo->_eii_iid.size());
        einfo->_eii_pheno.resize(einfo->_eii_iid.size());
        for(int i=0;i<einfo->_eii_iid.size();i++)
        {
            einfo->_eii_include[i]=i;
            einfo->_eii_fa_id[i]="0";
            einfo->_eii_mo_id[i]="0";
            einfo->_eii_sex[i]=-9;
            einfo->_eii_pheno[i]=-9;
            einfo->_eii_map.insert(pair<string,int>(einfo->_eii_fid[i]+":"+einfo->_eii_iid[i],i));
            if(einfo->_eii_map.size()==i)
            {
                LOGPRINTF("ERROR: Duplicate individual with FID: %s and IID %s.\n", einfo->_eii_fid[i].c_str(),einfo->_eii_iid[i].c_str());
                exit(EXIT_FAILURE);
            }
        }
        while(fgets(readBuf, MAXLINEBUFLEN, efile))  line_idx++;
        einfo->_eii_num=(uint32_t)einfo->_eii_iid.size();
        einfo->_epi_num=line_idx;
        einfo->_eii_pheno_num=1;
        msg="Expression Values";
        uint64_t sizeneed=einfo->_epi_num*einfo->_eii_num*sizeof(double)+line_idx*(sizeof(string)*3+sizeof(int)*5+sizeof(char));
        if(!sudoAllocReserved( sizeneed, msg)) exit(EXIT_FAILURE);
        einfo->_val.resize(einfo->_epi_num*einfo->_eii_num);
        einfo->_epi_prb.resize(line_idx);
        einfo->_epi_include.resize(line_idx);
        einfo->_epi_chr.resize(line_idx);
        einfo->_epi_bp.resize(line_idx);
        einfo->_epi_gd.resize(line_idx);
        einfo->_epi_gene.resize(line_idx);
        einfo->_epi_orien.resize(line_idx);

        rewind(efile);
        if(!no_fid_flag) fgets(readBuf, MAXLINEBUFLEN, efile); //fid
        fgets(readBuf, MAXLINEBUFLEN, efile); //iid
        memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
        line_idx=0;
        while(fgets(readBuf, MAXLINEBUFLEN, efile))
        {
            split_str(readBuf,strlist,0);
            if(strlist.size()!=1+einfo->_eii_num)
            {
                LOGPRINTF("ERROR: Inconsistent element number found in line %u with probe ID %s .\n", line_idx+1, strlist[0].c_str());
                exit(EXIT_FAILURE);
            }
            einfo->_epi_prb[line_idx]=strlist[0].c_str();
            einfo->_epi_include[line_idx]=line_idx;
            einfo->_epi_chr[line_idx]=-9;
            einfo->_epi_bp[line_idx]=-9;
            einfo->_epi_gd[line_idx]=-9;
            einfo->_epi_gene[line_idx]="NA";
            einfo->_epi_orien[line_idx]='*';
            einfo->_epi_map.insert(pair<string,int>(strlist[0],line_idx));
            if(einfo->_epi_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate probe %s .\n", einfo->_epi_prb[line_idx].c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=1;i<strlist.size();i++)
            {
                long iid=i-1;
                if(strlist[i]=="NA" || strlist[i]=="." || strlist[i]=="-" || strlist[i]=="-9") einfo->_val[line_idx*einfo->_eii_num+iid]=1e10;
                else {
                    double valtmp=atof(strlist[i].c_str());
                    if(!valueType && valtmp<0) {
                        LOGPRINTF("ERROR: Beta value of probe %s should be positive.\n", einfo->_epi_prb[line_idx].c_str());
                        exit(EXIT_FAILURE);
                    } else {
                        einfo->_val[line_idx*einfo->_eii_num+iid]=valtmp;
                    }
                }
                
            }
            line_idx++;
            memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
        }
        fclose(efile);
        deallocReserved(&readBuf, MAXLINEBUFLEN);
        LOGPRINTF("%s data for %u probes of %u individuals have been included from the file %s.\n", typestr.c_str(),  einfo->_epi_num, einfo->_eii_num, eFileName);
    }
    void beta_2_m(eInfo* einfo)
    {
        if(einfo->_valType==BETAVALUE)
        {
            LOGPRINTF("Transform beta-values to m-values...\n");
            for(long i=0;i<einfo->_val.size();i++) {
                double beta=einfo->_val[i];
                einfo->_val[i]=log2(beta/(1-beta));
            }
            einfo->_valType=MVALUE;
            LOGPRINTF("Transformation of %ld beta-values to m-values has finished.\n",einfo->_val.size());
        } else {
            LOGPRINTF("WARNING: This data are not beta-values. Transformation would not be performed.\n");
        }
        
    }
    void m_2_beta(eInfo* einfo)
    {
        if(einfo->_valType==MVALUE) {
            LOGPRINTF("Transform m-values to beta-values...\n");
            for(long i=0;i<einfo->_val.size();i++) {
                double m=einfo->_val[i];
                double z=pow(2, m);
                einfo->_val[i]=z/(z+1);
            }
            einfo->_valType=BETAVALUE;
            LOGPRINTF("Transformation of %ld m-values to beta-values has finished.\n",einfo->_val.size());
        } else {
            LOGPRINTF("WARNING: This data are not m-values. Transformation would not be performed.\n");
        }
        
    }
    void write_efile(char* outFileName, eInfo* einfo,bool impute_mean_flag)
    {
        if (impute_mean_flag && einfo->_mu.empty()) cal_var_mean(einfo, true, false);
        FILE* outfile = NULL;
        if (fopen_checked(&outfile, outFileName, "w")) exit(EXIT_FAILURE);
        char* writeBuf;
        string msg="Writing Buffer";
        if(!allocReserved(&writeBuf, MAXLINEBUFLEN,msg)) exit(EXIT_FAILURE);
        string typestr=getFileType(einfo->_eType);
        LOGPRINTF("Writing %s data to %s ...\n", typestr.c_str(),outFileName);
        char* wptr=writeBuf;
        char delim=' ';
        wptr +=sprintf(wptr,"FID IID");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            wptr +=sprintf(wptr,"%c%s",delim, einfo->_epi_prb[einfo->_epi_include[i]].c_str());
        }
        wptr +=sprintf(wptr,"%c",'\n');
        if(fwrite_checked(writeBuf,wptr-writeBuf,outfile))
        {
            LOGPRINTF("ERROR: in writing text file %s .\n", outFileName);
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<einfo->_eii_include.size();i++)
        {
            wptr=writeBuf;
            wptr +=sprintf(wptr,"%s %s",einfo->_eii_fid[einfo->_eii_include[i]].c_str(),einfo->_eii_iid[einfo->_eii_include[i]].c_str());
            for( int j=0;j<einfo->_epi_include.size();j++)
            {
                if(einfo->_val[einfo->_epi_include[j]*einfo->_eii_num+einfo->_eii_include[i]]>1e9) {
                    if(impute_mean_flag) wptr +=sprintf(wptr,"%c%.8lf",delim, einfo->_mu[einfo->_epi_include[j]]);
                    else wptr +=sprintf(wptr,"%c%s",delim, "NA");
                }
                else wptr +=sprintf(wptr,"%c%.8lf",delim, einfo->_val[einfo->_epi_include[j]*einfo->_eii_num+einfo->_eii_include[i]]);
            }
            wptr +=sprintf(wptr,"%c",'\n');
            if(fwrite_checked(writeBuf,wptr-writeBuf,outfile))
            {
                LOGPRINTF("ERROR: in writing text file %s .\n", outFileName);
                exit(EXIT_FAILURE);
            }
        }
        deallocReserved(&writeBuf, MAXLINEBUFLEN);
        fclose(outfile);
        string valtypestr=getValType(einfo->_valType);
        LOGPRINTF("%s of %s data for %ld probes of %ld individuals have been save in the file %s.\n",valtypestr.c_str(), typestr.c_str(), einfo->_epi_include.size(), einfo->_eii_include.size(), outFileName);
    }
    void write_tefile(char* outFileName, eInfo* einfo,bool impute_mean_flag)
    {
        if (impute_mean_flag && einfo->_mu.empty()) cal_var_mean(einfo, true, false);
        FILE* outfile = NULL;
        if (fopen_checked(&outfile, outFileName, "w")) exit(EXIT_FAILURE);
        char* writeBuf;
        string msg="Writing Buffer";
        if(!allocReserved(&writeBuf, MAXLINEBUFLEN,msg)) exit(EXIT_FAILURE);
        string typestr=getFileType(einfo->_eType);
        LOGPRINTF("Writing %s data to %s ...\n", typestr.c_str(),outFileName);
        char* wptr=writeBuf;
        char delim='\t';
        wptr +=sprintf(wptr,"ID");
        for(int i=0;i<einfo->_eii_include.size();i++)
        {
            wptr +=sprintf(wptr,"%c%s",delim, einfo->_eii_iid[einfo->_eii_include[i]].c_str());
        }
        wptr +=sprintf(wptr,"%c",'\n');
        if(fwrite_checked(writeBuf,wptr-writeBuf,outfile))
        {
            LOGPRINTF("ERROR: in writing text file %s .\n", outFileName);
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            wptr=writeBuf;
            wptr +=sprintf(wptr,"%s",einfo->_epi_prb[einfo->_epi_include[i]].c_str());
            for( int j=0;j<einfo->_eii_include.size();j++)
            {
                if(einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]]>1e9){
                   if(impute_mean_flag) wptr +=sprintf(wptr,"%c%.8lf",delim, einfo->_mu[einfo->_epi_include[i]]);
                   else wptr +=sprintf(wptr,"%c%s",delim, "NA");
                } else {
                    wptr +=sprintf(wptr,"%c%.8lf",delim, einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]]);
                }
            }
            wptr +=sprintf(wptr,"%c",'\n');
            if(fwrite_checked(writeBuf,wptr-writeBuf,outfile))
            {
                LOGPRINTF("ERROR: in writing text file %s .\n", outFileName);
                exit(EXIT_FAILURE);
            }
        }
        deallocReserved(&writeBuf, MAXLINEBUFLEN);
        fclose(outfile);
        string valtypestr=getValType(einfo->_valType);
        LOGPRINTF("%s of %s data for %ld probes of %ld individuals have been save in the file %s.\n",valtypestr.c_str(), typestr.c_str(), einfo->_epi_include.size(), einfo->_eii_include.size(), outFileName);
    }

    void write_eii(char* outFileName, eInfo* einfo)
    {
        // not save 23 to X and 24 for Y
        FILE* efile=NULL;
        string eiiName=string(outFileName)+".eii";
        if(fopen_checked(&efile, eiiName.c_str(),"w")) exit(EXIT_FAILURE);
        for(int i=0;i<einfo->_eii_include.size();i++)
        {
            string str=einfo->_eii_fid[einfo->_eii_include[i]]+'\t'+einfo->_eii_iid[einfo->_eii_include[i]]+'\t'+einfo->_eii_fa_id[einfo->_eii_include[i]]+'\t'+einfo->_eii_mo_id[einfo->_eii_include[i]]+'\t'+atos(einfo->_eii_sex[einfo->_eii_include[i]])+'\t'+atos(einfo->_eii_pheno[einfo->_eii_include[i]])+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                
                LOGPRINTF("ERROR: in writing file %s .\n", eiiName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        fclose(efile);
        LOGPRINTF("%ld individuals have been saved in the file %s .\n", einfo->_eii_include.size(), eiiName.c_str());

    }
    void write_epi(char* outFileName, eInfo* einfo)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".epi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) exit(EXIT_FAILURE);
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            string str=atos(einfo->_epi_chr[einfo->_epi_include[i]])+'\t'+einfo->_epi_prb[einfo->_epi_include[i]]+'\t'+atos(einfo->_epi_gd[einfo->_epi_include[i]])+'\t'+atos(einfo->_epi_bp[einfo->_epi_include[i]])+'\t'+einfo->_epi_gene[einfo->_epi_include[i]]+'\t'+atos(einfo->_epi_orien[einfo->_epi_include[i]])+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        fclose(efile);
        LOGPRINTF("%ld probes have been saved in the file %s .\n", einfo->_epi_include.size(), epiName.c_str());
    }
    void write_beed(char* outFileName, eInfo* einfo)
    {
        FILE* efile=NULL;
        string beedName=string(outFileName)+".beed";
        if(fopen_checked(&efile, beedName.c_str(),"wb")) exit(EXIT_FAILURE);
        uint32_t indicator=(einfo->_eType << 8)| einfo->_valType;
        
        if(einfo->_eii_include.size()==einfo->_eii_num && einfo->_epi_include.size()==einfo->_epi_num )
        {
            if (fwrite_checked(&indicator, sizeof(uint32_t), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
            if (fwrite_checked(&einfo->_eii_num, sizeof(uint32_t), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
            if (fwrite_checked(&einfo->_epi_num, sizeof(uint32_t), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
            if (fwrite_checked(&einfo->_val[0], einfo->_val.size()*sizeof(double), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
        }else{
            
            if (fwrite_checked(&indicator, sizeof(uint32_t), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
            uint32_t inum=(uint32_t)einfo->_eii_include.size();
            if (fwrite_checked(&inum, sizeof(uint32_t), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
            uint32_t pnum=(uint32_t)einfo->_epi_include.size();
            if (fwrite_checked(&pnum, sizeof(uint32_t), efile))
            {
                LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                exit(EXIT_FAILURE);
            }
            char* writeBuf;
            uint64_t wBufsize=allocReserved(&writeBuf,einfo->_epi_include.size()*einfo->_eii_include.size()*sizeof(double));
            double* dptr=(double *)writeBuf;
            for(int i=0;i<einfo->_epi_include.size();i++)
             {
                 for(int j=0;j<einfo->_eii_include.size();j++)
                 {
                     *dptr++=einfo->_val[einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]];
                     if((char*)dptr-writeBuf==wBufsize)
                     {
                         if (fwrite_checked(writeBuf, wBufsize, efile))
                         {
                             LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                             exit(EXIT_FAILURE);
                         }
                         dptr=(double*)writeBuf;
                     }                     
                 }
             }
            if((char*)dptr-writeBuf>0)
            {
                if (fwrite_checked(writeBuf, (char*)dptr-writeBuf, efile))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", beedName.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            deallocReserved(&writeBuf, wBufsize);
        }
        fclose(efile);
        string typestr=getFileType(einfo->_eType);
        string valtypestr=getValType(einfo->_valType);
        LOGPRINTF("%s of %s data for %lu probes of %lu individuals have been saved in the file %s.\n", valtypestr.c_str(), typestr.c_str(), einfo->_epi_include.size(), einfo->_eii_include.size(), beedName.c_str());
    }
    void init_einfo(eInfo* einfo)
    {
        einfo->_epi_num=0;
        einfo->_eii_num=0;
        einfo->_eii_cov_num=0;
        einfo->_eii_qcov_num=0;
        einfo->_eii_pheno_num=0;
        einfo->_val.clear();
        einfo->_mu.clear();
        einfo->_var.clear();
    }
    void read_eii(char* eiiFileName,eInfo* einfo)
    {
        FILE* eiifile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        if(fopen_checked(&eiifile, eiiFileName,"r")) exit(EXIT_FAILURE);
        LOGPRINTF("Reading individual information from %s ...\n", eiiFileName);
        einfo->_eii_fid.clear();
        einfo->_eii_iid.clear();
        einfo->_eii_fa_id.clear();
        einfo->_eii_mo_id.clear();
        einfo->_eii_sex.clear();
        einfo->_eii_pheno.clear();
        einfo->_eii_include.clear();
        einfo->_eii_map.clear();
        while(fgets(Tbuf, MAX_LINE_SIZE, eiifile))
        {
            split_str(Tbuf,strlist,0);
            if(strlist.size()>6)
            {
                LOGPRINTF("WARNING: Line %u has more than 6 items.\n", line_idx);
            } else if(strlist.size()<6)
            {
               LOGPRINTF("ERROR: Line %u has less than 6 items.\n", line_idx);
               exit(EXIT_FAILURE);
            }
            einfo->_eii_map.insert(pair<string,int>(strlist[0]+":"+strlist[1],line_idx));
            if(einfo->_eii_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate individual with FID: %s and IID %s.\n", strlist[0].c_str(),strlist[1].c_str());
                exit(EXIT_FAILURE);
            }
            einfo->_eii_fid.push_back(strlist[0]);
            einfo->_eii_iid.push_back(strlist[1]);
            einfo->_eii_fa_id.push_back(strlist[2]);
            einfo->_eii_mo_id.push_back(strlist[3]);
            einfo->_eii_sex.push_back(atoi(strlist[4].c_str()));
            einfo->_eii_pheno.push_back(atof(strlist[5].c_str()));
            einfo->_eii_include.push_back(line_idx);
            line_idx++;
        }
        einfo->_eii_num=line_idx;
        einfo->_eii_pheno_num=1;
        fclose(eiifile);
        LOGPRINTF("%u individuals to be included from  %s .\n", einfo->_eii_num, eiiFileName);
    }
    void read_epi(char* epiFileName,eInfo* einfo)
    {
        FILE* epifile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        if(fopen_checked(&epifile, epiFileName,"r")) exit(EXIT_FAILURE);
        LOGPRINTF("Reading probe information from %s ...\n", epiFileName);
        einfo->_epi_chr.clear();
        einfo->_epi_prb.clear();
        einfo->_epi_gd.clear();
        einfo->_epi_bp.clear();
        einfo->_epi_gene.clear();
        einfo->_epi_orien.clear();
        einfo->_epi_include.clear();
        einfo->_epi_map.clear();
        while(fgets(Tbuf, MAX_LINE_SIZE, epifile))
        {
            split_str(Tbuf,strlist,0);
            if(strlist.size()>6)
            {
                LOGPRINTF("WARNING: Line %u has more than 6 items.\n", line_idx);
            }
            if(strlist.size()<6)
            {
                LOGPRINTF("ERROR: Line %u has less than 6 items.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            einfo->_epi_map.insert(pair<string,int>(strlist[1],line_idx));
            if(einfo->_epi_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate probe : %s.\n", strlist[1].c_str());
                exit(EXIT_FAILURE);
            }
            if(strlist[0]=="X" || strlist[0]=="x") einfo->_epi_chr.push_back(23);
            else if(strlist[0]=="Y" || strlist[0]=="y") einfo->_epi_chr.push_back(24);
            else einfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            einfo->_epi_prb.push_back(strlist[1]);
            einfo->_epi_gd.push_back(atoi(strlist[2].c_str()));
            einfo->_epi_bp.push_back(atoi(strlist[3].c_str()));
            einfo->_epi_gene.push_back(strlist[4].c_str());
            einfo->_epi_orien.push_back(strlist[5][0]);
            einfo->_epi_include.push_back(line_idx);
            line_idx++;
        }
        einfo->_epi_num=line_idx;
        fclose(epifile);
        LOGPRINTF("%u probes to be included from  %s .\n", einfo->_epi_num, epiFileName);
    }
    void read_beed(char* beedFileName,eInfo* einfo)
    {
        FILE* beedfile=NULL;
        if(fopen_checked(&beedfile,beedFileName,"rb")) exit(EXIT_FAILURE);
        uint32_t indicator;
        if(fread(&indicator, sizeof(uint32_t),1, beedfile)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", beedFileName);
            exit(EXIT_FAILURE);
        }
        einfo->_valType=indicator & 0xF;
        einfo->_eType=(indicator & 0xFFF0)>>8;
        uint32_t inum;
        if(fread(&inum, sizeof(uint32_t),1, beedfile)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", beedFileName);
            exit(EXIT_FAILURE);
        }
        uint32_t pnum;
        if(fread(&pnum, sizeof(uint32_t),1, beedfile)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", beedFileName);
            exit(EXIT_FAILURE);
        }
        
        uint64_t cur_pos = ftell( beedfile );
        fseek( beedfile, 0L, SEEK_END );
        uint64_t size_file = ftell( beedfile );
        fseek( beedfile, cur_pos, SEEK_SET );
        string typestr=getFileType(einfo->_eType);
        uint64_t val_ttl_file=einfo->_epi_num*einfo->_eii_num;
        if(inum!=einfo->_eii_num || pnum!=einfo->_epi_num)
        {
            LOGPRINTF("ERROR: .eii or .epi file is not consistent with .beed file %s!\n", beedFileName);
            exit(EXIT_FAILURE);
        }
        if(size_file!=val_ttl_file*sizeof(double)+3*sizeof(uint32_t))
        {
            LOGPRINTF("ERROR: File %s is broken!\n", beedFileName);
            exit(EXIT_FAILURE);
        }
        uint64_t val_ttl_icld=einfo->_epi_include.size()*einfo->_eii_include.size();
        if(!sudoAllocReserved( val_ttl_icld*sizeof(double), typestr)) exit(EXIT_FAILURE);
        einfo->_val.resize(val_ttl_icld);
        if(einfo->_eii_include.size()==einfo->_eii_num && einfo->_epi_include.size()==einfo->_epi_num )
        {
            if(fread(&einfo->_val[0], sizeof(double),val_ttl_icld, beedfile)!=val_ttl_icld)
            {
                LOGPRINTF("ERROR: File %s read failed!\n", beedFileName);
                exit(EXIT_FAILURE);
            }
        }
        else if(einfo->_epi_include.size()<einfo->_epi_num/2)
        {
            char* readBuf;
            string msg="Reading Buffer";
            if(!allocReserved(&readBuf, einfo->_eii_num*sizeof(double),msg)) exit(EXIT_FAILURE);
            for(int i=0; i<einfo->_epi_include.size();i++)
            {
                memset(readBuf,0,sizeof(double)*(einfo->_eii_num));
                int prbid=einfo->_epi_include[i];
                uint64_t readpos=3*sizeof(uint32_t)+prbid*einfo->_eii_num*sizeof(double);
                fseek( beedfile, readpos, SEEK_SET );
                if(fread(readBuf, sizeof(double),einfo->_eii_num, beedfile)!=einfo->_eii_num)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", beedFileName);
                    exit(EXIT_FAILURE);
                }
                for(int j=0;j<einfo->_eii_include.size();j++)
                {
                    double* ptr=(double *)readBuf;
                    einfo->_val[i*einfo->_eii_include.size()+j]=*(ptr+einfo->_eii_include[j]);
                }
                
            }
            deallocReserved(&readBuf, einfo->_eii_num*sizeof(double));
        }
        else
        {
            char* readBuf;
            uint64_t rBufsize=allocReserved(&readBuf,einfo->_epi_num*einfo->_eii_num*sizeof(double));
            double* dptr=(double *)readBuf;
            
            if(rBufsize==einfo->_epi_num*einfo->_eii_num*sizeof(double))
            {
                if(fread(readBuf, 1,rBufsize, beedfile)!=rBufsize)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", beedFileName);
                    exit(EXIT_FAILURE);
                }
                for(int i=0; i<einfo->_epi_include.size();i++)
                    for(int j=0;j<einfo->_eii_include.size();j++)
                        einfo->_val[i*einfo->_eii_include.size()+j]=*(dptr+einfo->_epi_include[i]*einfo->_eii_num+einfo->_eii_include[j]);
            }
            else
            {
                uint64_t prb_per_time=rBufsize/(einfo->_eii_num*sizeof(double));
                uint64_t prb_lower_boundary=0, prb_upper_boundary=0;
                int i=0;
                while (!feof(beedfile)) {
                    memset(readBuf,0,sizeof(char)*rBufsize);
                    uint64_t nSize = fread(readBuf, sizeof(double), prb_per_time*einfo->_eii_num, beedfile);
                    uint64_t nPrb=nSize/einfo->_eii_num;
                    prb_upper_boundary+=nPrb;
                    for(; i<einfo->_epi_include.size();i++)
                    {
                        if(einfo->_epi_include[i]>=prb_upper_boundary) break;
                        for(int j=0;j<einfo->_eii_include.size();j++)
                            einfo->_val[i*einfo->_eii_include.size()+j]=*(dptr+(einfo->_epi_include[i]-prb_lower_boundary)*einfo->_eii_num+einfo->_eii_include[j]);
                    }
                    prb_lower_boundary=prb_upper_boundary;                    
                }
            }
             deallocReserved(&readBuf, rBufsize);
            
        }
        fclose(beedfile);
        if(einfo->_eii_include.size()<einfo->_eii_num) update_eii(einfo);
        if(einfo->_epi_include.size()<einfo->_epi_num) update_epi(einfo);
        string valtypestr=getValType(einfo->_valType);
        LOGPRINTF("%s of %s data for %u probes of %u individuals have been included from the file %s.\n",valtypestr.c_str(), typestr.c_str(),  einfo->_epi_num, einfo->_eii_num, beedFileName);
    }
    void update_eii(eInfo* einfo)
    {
        vector<int> inttmp;
        vector<string> strtmp;
        vector<double> dtmp;
        for (int i = 0; i < einfo->_eii_include.size(); i++) strtmp.push_back(einfo->_eii_fid[einfo->_eii_include[i]]);
        einfo->_eii_fid.clear();
        einfo->_eii_fid.swap(strtmp);
        for (int i = 0; i < einfo->_eii_include.size(); i++) strtmp.push_back(einfo->_eii_iid[einfo->_eii_include[i]]);
        einfo->_eii_iid.clear();
        einfo->_eii_iid.swap(strtmp);
        for (int i = 0; i < einfo->_eii_include.size(); i++) strtmp.push_back(einfo->_eii_fa_id[einfo->_eii_include[i]]);
        einfo->_eii_fa_id.clear();
        einfo->_eii_fa_id.swap(strtmp);
        for (int i = 0; i < einfo->_eii_include.size(); i++) strtmp.push_back(einfo->_eii_mo_id[einfo->_eii_include[i]]);
        einfo->_eii_mo_id.clear();
        einfo->_eii_mo_id.swap(strtmp);
        for (int i = 0; i < einfo->_eii_include.size(); i++) inttmp.push_back(einfo->_eii_sex[einfo->_eii_include[i]]);
        einfo->_eii_sex.clear();
        einfo->_eii_sex.swap(inttmp);
        for(int ii=0;ii<einfo->_eii_pheno_num;ii++)
            for (int i = 0; i < einfo->_eii_include.size(); i++) dtmp.push_back(einfo->_eii_pheno[ii*einfo->_eii_num+einfo->_eii_include[i]]);
        einfo->_eii_pheno.clear();
        einfo->_eii_pheno.swap(dtmp);
        for(int ii=0;ii<einfo->_eii_cov_num;ii++)
            for (int i = 0; i < einfo->_eii_include.size(); i++) strtmp.push_back(einfo->_eii_cov[ii*einfo->_eii_num+einfo->_eii_include[i]]);
        einfo->_eii_cov.clear();
        einfo->_eii_cov.swap(strtmp);
        for(int ii=0;ii<einfo->_eii_qcov_num;ii++)
            for (int i = 0; i < einfo->_eii_include.size(); i++) dtmp.push_back(einfo->_eii_qcov[ii*einfo->_eii_num+einfo->_eii_include[i]]);
        einfo->_eii_qcov.clear();
        einfo->_eii_qcov.swap(dtmp);
        einfo->_eii_num =(uint32_t)einfo->_eii_include.size();
        
        einfo->_eii_include.clear();
        einfo->_eii_map.clear();
        for (int i = 0; i < einfo->_eii_num; i++)
        {
            einfo->_eii_include.push_back(i);
            einfo->_eii_map.insert(pair<string,int>( einfo->_eii_fid[i]+":"+einfo->_eii_iid[i],i));
        }
        

    }
    void update_epi(eInfo* einfo)
    {
        vector<int> inttmp;
        vector<string> strtmp;
        vector<char> chartmp;
        for (int i = 0; i < einfo->_epi_include.size(); i++) inttmp.push_back(einfo->_epi_chr[einfo->_epi_include[i]]);
        einfo->_epi_chr.clear();
        einfo->_epi_chr.swap(inttmp);
        for (int i = 0; i < einfo->_epi_include.size(); i++) strtmp.push_back(einfo->_epi_prb[einfo->_epi_include[i]]);
        einfo->_epi_prb.clear();
        einfo->_epi_prb.swap(strtmp);
        for (int i = 0; i < einfo->_epi_include.size(); i++) inttmp.push_back(einfo->_epi_gd[einfo->_epi_include[i]]);
        einfo->_epi_gd.clear();
        einfo->_epi_gd.swap(inttmp);
        for (int i = 0; i < einfo->_epi_include.size(); i++) inttmp.push_back(einfo->_epi_bp[einfo->_epi_include[i]]);
        einfo->_epi_bp.clear();
        einfo->_epi_bp.swap(inttmp);
        for (int i = 0; i < einfo->_epi_include.size(); i++) strtmp.push_back(einfo->_epi_gene[einfo->_epi_include[i]]);
        einfo->_epi_gene.clear();
        einfo->_epi_gene.swap(strtmp);
        for (int i = 0; i < einfo->_epi_include.size(); i++) chartmp.push_back(einfo->_epi_orien[einfo->_epi_include[i]]);
        einfo->_epi_orien.clear();
        einfo->_epi_orien.swap(chartmp);
        einfo->_epi_num =(uint32_t)einfo->_epi_include.size();
        
        einfo->_epi_include.clear();
        einfo->_epi_map.clear();
        for (int i = 0; i < einfo->_epi_num; i++)
        {
            einfo->_epi_include.push_back(i);
            einfo->_epi_map.insert(pair<string,int>( einfo->_epi_prb[i],i));
        }

    }
    void update_epifile(char* befileName,char* s_epiFileName)
    {
        eInfo einfo;
        if(befileName==NULL)
        {
            LOGPRINTF("Error: please input the Gene expression / Methylation data by the option --befile.\n");
            exit(EXIT_FAILURE);
        }
        if(s_epiFileName==NULL)
        {
            LOGPRINTF("Error: please input the probe information file by the option --update-epi.\n");
            exit(EXIT_FAILURE);
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,befileName,strlen(befileName)+1);
        char* suffix=inputname+strlen(befileName);
        memcpy(suffix,".epi",5);
        read_epi(inputname,&einfo);
        memcpy(suffix,".bak",5);
        write_epi(inputname, &einfo);
        
        FILE* epifile=NULL;
        map<string, int>::iterator iter;
        vector<string> strlist;
        uint32_t line_idx = 0;
        if(fopen_checked(&epifile, s_epiFileName,"r")) exit(EXIT_FAILURE);
        LOGPRINTF("Reading probe information from %s ...\n", s_epiFileName);
        while(fgets(Tbuf, MAX_LINE_SIZE, epifile))
        {
            split_str(Tbuf,strlist,0);
            if(strlist.size()>6)
            {
                LOGPRINTF("WARNING: Line %u has more than 6 items.\n", line_idx);
            }
            if(strlist.size()<6)
            {
                LOGPRINTF("ERROR: Line %u has less than 6 items.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            iter=einfo._epi_map.find(strlist[1]);
            if(iter!=einfo._epi_map.end())
            {
                int probeidx=iter->second;
                if(strlist[0]=="X" || strlist[0]=="x") einfo._epi_chr[probeidx]=23;
                else if(strlist[0]=="Y" || strlist[0]=="y") einfo._epi_chr[probeidx]=24;
                else einfo._epi_chr[probeidx]=atoi(strlist[0].c_str());
                einfo._epi_gd[probeidx]=atoi(strlist[2].c_str());
                einfo._epi_bp[probeidx]=atoi(strlist[3].c_str());
                einfo._epi_gene[probeidx]=strlist[4].c_str();
                einfo._epi_orien[probeidx]=strlist[5][0];
            }
       
            line_idx++;
        }
        fclose(epifile);
        memcpy(inputname,befileName,strlen(befileName)+1);
        write_epi(inputname, &einfo);
    }
    void keep_indi(eInfo* einfo,string indi_list_file)
    {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        update_map_kp(indi_list, einfo->_eii_map, einfo->_eii_include);
        LOGPRINTF("%ld individuals are kept from %s.\n",  einfo->_eii_include.size(), indi_list_file.c_str());
    }
    void remove_indi(eInfo* einfo, string indi_list_file)
    {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        long prev_size = einfo->_eii_include.size();
        update_map_rm(indi_list, einfo->_eii_map, einfo->_eii_include);
        LOGPRINTF("%ld individuals are removed from %s and there are %ld individuals remaining.\n",  prev_size - einfo->_eii_include.size(), indi_list_file.c_str(),einfo->_eii_include.size());
    }
    void extract_probe_by_chr(eInfo* einfo, int prbchr)
    {
        vector<string> prblst;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int tmpint=einfo->_epi_include[i];
            if(einfo->_epi_chr[tmpint]==prbchr ) prblst.push_back(einfo->_epi_prb[tmpint]);
        }
        update_map_kp(prblst, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from chromosome %d.\n",  einfo->_epi_include.size(), prbchr);
    }
    void extract_probe(eInfo* einfo,string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        update_map_kp(problist, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from %s.\n",  einfo->_epi_include.size(), problstName.c_str());
    }
    void extract_probe_by_gene(eInfo* einfo, string genelistName)
    {
        vector<string> genelist;
        string msg="genes";
        read_msglist(genelistName, genelist,msg);
        
        vector<string> prblst;
        for(int i=0;i<genelist.size();i++)
        {
            string tmpname1=genelist[i];
            for(int j=0;j<einfo->_epi_include.size();j++)
            {
                string tmpname2=einfo->_epi_gene[einfo->_epi_include[j]];
                if(tmpname1==tmpname2)  prblst.push_back(einfo->_epi_prb[einfo->_epi_include[j]]);
                else
                {
                    vector<string> substrs;
                    uint64_t tmpnum=split_string_skip(tmpname2,substrs,",;",0);
                    if(tmpnum>1)
                        for(int k=0;k<tmpnum;k++)
                            if(tmpname1==substrs[k])
                            {
                                prblst.push_back(einfo->_epi_prb[einfo->_epi_include[j]]);
                                break;
                            }
                }
            }
        }
        update_map_kp(prblst, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from %s.\n",  einfo->_epi_include.size(), genelistName.c_str());
    }
    void extract_probe(eInfo* einfo, string prbname, int prbWind)
    {
        long idx=find(einfo->_epi_prb.begin(), einfo->_epi_prb.end(), prbname)-einfo->_epi_prb.begin();
        if(idx==einfo->_epi_prb.size())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  prbname.c_str());
            exit(EXIT_FAILURE);
        }
        int prbbp=einfo->_epi_bp[idx];
        int prbchr=einfo->_epi_chr[idx];
        int upbound=prbbp+prbWind*1000;
        int tmpint=prbbp-prbWind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<string> prblst;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            tmpint=einfo->_epi_include[i];
            if(einfo->_epi_chr[tmpint]==prbchr && einfo->_epi_bp[tmpint]>=lowbound && einfo->_epi_bp[tmpint]<=upbound) prblst.push_back(einfo->_epi_prb[einfo->_epi_include[i]]);
        }
        update_map_kp(prblst, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from the region %d Kb around %s.\n",  einfo->_epi_include.size(), prbWind, prbname.c_str());
    }
    void extract_single_probe(eInfo* einfo, string prbname)
    {
        long idx=find(einfo->_epi_prb.begin(), einfo->_epi_prb.end(), prbname)-einfo->_epi_prb.begin();
        if(idx==einfo->_epi_prb.size())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  prbname.c_str());
            exit(EXIT_FAILURE);
        }
        einfo->_epi_include.clear();
        einfo->_epi_map.clear();
        einfo->_epi_include.push_back((int)idx);
        einfo->_epi_map.insert(pair<string,int>(prbname,0));
        cout << prbname << " is extracted. " << endl;
    }
    void extract_probe(eInfo* einfo, string fromprbname, string toprbname)
    {
        long fromidx=find(einfo->_epi_prb.begin(), einfo->_epi_prb.end(), fromprbname)-einfo->_epi_prb.begin();
        if(fromidx==einfo->_epi_prb.size())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  fromprbname.c_str());
            exit(EXIT_FAILURE);
        }
        int fromprbbp=einfo->_epi_bp[fromidx];
        int prbchr=einfo->_epi_chr[fromidx];
        
        long toidx=find(einfo->_epi_prb.begin(), einfo->_epi_prb.end(), toprbname)-einfo->_epi_prb.begin();
        if(toidx==einfo->_epi_prb.size())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  toprbname.c_str());
            exit(EXIT_FAILURE);
        }
        int toprbbp=einfo->_epi_bp[toidx];
        int toprbchr=einfo->_epi_chr[toidx];
        if(toprbchr != prbchr)
        {
            LOGPRINTF("ERROR: probe %s and probe %s are not from the same chromosome.\n", fromprbname.c_str(), toprbname.c_str());
            exit(EXIT_FAILURE);
        }
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<string> prblst;
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            int tmpint=einfo->_epi_include[i];
            if(einfo->_epi_chr[tmpint]==prbchr && einfo->_epi_bp[tmpint]>=fromprbbp && einfo->_epi_bp[tmpint]<=toprbbp) prblst.push_back(einfo->_epi_prb[tmpint]);
        }
        update_map_kp(prblst, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from probe %s to probe %s.\n",  einfo->_epi_include.size(), fromprbname.c_str(), toprbname.c_str());
    }
    void extract_probe(eInfo* einfo, int fromprbkb, int toprbkb, int chr)
    {
        int fromprbbp=fromprbkb*1000;
        int toprbbp=toprbkb*1000;
        
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<string> prblst;
        for(int i=0;einfo->_epi_include.size();i++)
        {
            int tmpint=einfo->_epi_include[i];
            if(einfo->_epi_chr[tmpint]==chr && einfo->_epi_bp[tmpint]>=fromprbbp && einfo->_epi_bp[tmpint]<=toprbbp) prblst.push_back(einfo->_epi_prb[tmpint]);
        }
        update_map_kp(prblst, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from probe bp %d to probe bp %d.\n",  einfo->_epi_include.size(), fromprbkb, toprbkb);
    }
    void extract_probe_by_single_gene(eInfo* einfo, string genename)
    {
        vector<string> prblst;
        for(int j=0;j<einfo->_epi_include.size();j++)
        {
            string tmpname2=einfo->_epi_gene[einfo->_epi_include[j]];
            if(genename==tmpname2)  prblst.push_back(einfo->_epi_prb[einfo->_epi_include[j]]);
            else
            {
                vector<string> substrs;
                uint64_t tmpnum=split_string_skip(tmpname2,substrs,",;",0);
                if(tmpnum>1)
                    for(int k=0;k<tmpnum;k++)
                        if(genename==substrs[k])
                        {
                            prblst.push_back(einfo->_epi_prb[einfo->_epi_include[j]]);
                            break;
                        }
            }
        }
        update_map_kp(prblst, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are extracted from gene %s.\n",  einfo->_epi_include.size(), genename.c_str());
    }
    void exclude_probe(eInfo* einfo, string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        long pre_num=einfo->_epi_include.size();
        update_map_rm(problist, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("%ld probes are excluded from %s and there are %ld probe remaining.\n",  pre_num - einfo->_epi_include.size(), problstName.c_str(),einfo->_epi_include.size());
    }
    void exclude_single_probe(eInfo* einfo, string prbname)
    {
        vector<string> problist;
        problist.push_back(prbname);
        update_map_rm(problist, einfo->_epi_map, einfo->_epi_include);
        LOGPRINTF("Probe %s are excluded and there are %ld probe remaining.\n", prbname.c_str(), einfo->_epi_include.size());
    }
    void epi_man(eInfo* einfo,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde)
    {
        string logstr;
        int flags4prb=0;
        if(problstName != NULL) flags4prb++;
        if(prbname != NULL) flags4prb++;
        if(fromprbname!=NULL) flags4prb++;
        if(fromprbkb>=0) flags4prb++;
        if(genename != NULL) flags4prb++;
        if(flags4prb>1)
        {
             LOGPRINTF("WARNING: Flags for probes in this section are mutual exclusive. The priority order (from high to low) is: --extract-probe, --gene-list, --probe-wind, --probe, --from(to)--probe, --from(to)-probe-kb, --gene.\n");
        }
        if(chr>0)
        {
            extract_probe_by_chr(einfo, chr);
        }
        
        if(problstName != NULL || genelistName != NULL)
        {
            if(problstName != NULL) extract_probe(einfo, problstName);
            if(genelistName != NULL) extract_probe_by_gene(einfo, genelistName);
        }
        else if(prbwindFlag)
        {
            if(prbname==NULL)
            {
                LOGPRINTF("ERROR: Please identify the probe name by --probe when using --probe-wind.\n");
                exit(EXIT_FAILURE);
            }
            extract_probe(einfo, prbname, prbWind);
        }
        else if(prbname!=NULL)
        {
            extract_single_probe(einfo, prbname);
        }
        else if(fromprbname!=NULL)
        {
            if(toprbname==NULL)
            {
                LOGPRINTF("ERROR: Please identify the probe name by --to-probe.\n");
                exit(EXIT_FAILURE);
            }
            extract_probe(einfo, fromprbname, toprbname);
        }
        else if(fromprbkb>=0)
        {
            if(fromprbkb>=0 && chr==0 ) {
                LOGPRINTF("ERROR: Please identify the chromosome by --chr.\n");
                exit(EXIT_FAILURE);
            }
            if(toprbkb<0)
            {
                LOGPRINTF("ERROR: probe BP can't be negative.\n");
                exit(EXIT_FAILURE);
            }
            extract_probe(einfo, fromprbkb, toprbkb,chr);
        }
        else if(genename!=NULL)
        {
            extract_probe_by_single_gene(einfo, genename);
        }
        
        if(problst2exclde!=NULL)
        {
            exclude_probe(einfo, problst2exclde);
        }
        else if(probe2exclde!=NULL)
        {
            exclude_single_probe(einfo,probe2exclde);
        }
    }
     void eii_man(eInfo* einfo,char* indilstName,char* indilst2remove)
    {
        if(indilstName!=NULL) keep_indi(einfo,indilstName);
        if(indilst2remove!=NULL) remove_indi(einfo,indilst2remove);
    }
    
    /*Standardise probes and save the standardised probes in an eigen matrix. missing values in the standardised matrix are saved as 0.0 that means mean*/
    void std_probe(eInfo* einfo, vector< vector<bool> > &X_bool, bool divid_by_std, MatrixXd &_probe_data, bool output_profile)
    {
        Map<MatrixXd> X(&einfo->_val[0],einfo->_eii_num,einfo->_epi_num);
        
        uint64_t n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        _probe_data.resize(n,m);
       
        
#pragma omp parallel for
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++) _probe_data(i,j)=X(einfo->_eii_include[i], einfo->_epi_include[j]);
        }
        VectorXd mu(m), nonmiss(m);
        
        X_bool.resize(n);
        for(int i=0; i<n; i++) X_bool[i].resize(m);
        for(int j=0; j<m; j++){
            mu(j)=0.0;
            nonmiss(j)=0.0;
            for(int i=0; i<n; i++){
                if(_probe_data(i,j)<1e9){
                    mu(j)+=_probe_data(i,j);
                    nonmiss(j)+=1.0;
                    X_bool[i][j] = true;
                }
                else X_bool[i][j] = false;
            }
           if(nonmiss(j)>0) mu(j)/=nonmiss(j);
        }
        for(int j=0; j<m; j++)
            if(nonmiss(j)<n)
            {
                LOGPRINTF("WARNING: Missing value found.\n");
                break;
            }
#pragma omp parallel for
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                if(_probe_data(i,j)<1e9) _probe_data(i,j) -= mu(j);
                else _probe_data(i,j) = 0.0;
            }
        }
        if(output_profile)
        {
            einfo->_profile=new double[n*m]; // alloc memory
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    einfo->_profile[i*m+j]=_probe_data(i,j);
                }
            }
            
        }
       
        if(divid_by_std){
            VectorXd sd(m);
#pragma omp parallel for
            for(int j=0; j<m; j++){
                if(nonmiss(j)>1) sd(j) = sqrt((_probe_data.col(j).dot(_probe_data.col(j))) / (nonmiss(j) - 1.0));
                else sd(j)=0.0;
            }
#pragma omp parallel for
            for(int i=0; i<n; i++){
                for(int j=0; j<m; j++) {
                    if(fabs(sd(j))>1e-30) _probe_data(i,j) /= sd(j);
                    else _probe_data(i,j) = 0.0;
                }
            }
        }
        
    }
    /* standardisze eigen matrix */
    void std_probe_in_place( bool divid_by_std, MatrixXd &_probe_data)
    {
        
        uint64_t n = _probe_data.rows(), m = _probe_data.cols();
        for(int j=0; j<m; j++){
            double mu=0.0;
            double nonmiss=0.0;
            for(int i=0; i<n; i++){
                if(_probe_data(i,j)<1e9){
                    mu+=_probe_data(i,j);
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            
            #pragma omp parallel for
            for(int i=0; i<n; i++){
                if(_probe_data(i,j)<1e9) _probe_data(i,j) -= mu;
                else _probe_data(i,j) = 0.0;
            }
            
            if(divid_by_std){
                double sd = sqrt((_probe_data.col(j).dot(_probe_data.col(j))) / (nonmiss - 1.0));
                for(int i=0; i<n; i++){
                    _probe_data(i,j) /= sd;
                }
            }
        }
    }
    /* QC step: remove the probes which standard deviation is less than the threshold*/
    void std_probe_filtering( eInfo* einfo, double std_thresh)
    {
        LOGPRINTF("Excluding sites with low variance (std < %f)...\n",std_thresh);
        long n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        uint64_t ttl_n = einfo->_eii_num;
        vector<int> keep;
        for(int j=0; j<m; j++){
            double mu=0.0;
            double nonmiss=0.0;
            for(int i=0; i<n; i++){
                double val=einfo->_val[einfo->_epi_include[j]*ttl_n+einfo->_eii_include[i]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            
            double sd =0.0;
            for(int i=0; i<n; i++){
                double val=einfo->_val[einfo->_epi_include[j]*ttl_n+einfo->_eii_include[i]];
                if(val<1e9) sd += (val-mu)*(val-mu);
            }
            sd=sqrt(sd/(nonmiss - 1.0));
            if(sd>=std_thresh) keep.push_back(einfo->_epi_include[j]);
        }
        einfo->_epi_map.clear();
        for(int i=0;i<keep.size();i++)
            einfo->_epi_map.insert(pair<string,int>(einfo->_epi_prb[keep[i]],keep[i]));
        einfo->_epi_include.swap(keep);
        
        LOGPRINTF("%ld of %ld sites are Excluded with low variance (std < %f) and %ld sites left.\n",m-einfo->_epi_include.size(),m,std_thresh,einfo->_epi_include.size());
    }
    /* QC step: remove the probes whose mean value is larger than upper beta value or smaller than lower beta value*/
    void filtering_constitutive_probes( eInfo* einfo, double upper_beta_thresh,double lower_beta_thresh)
    {
        if(einfo->_valType==0) {
            LOGPRINTF("Excluding constitutively methylated/unmethylated sites with upper beta threshold %f and lower beta threshold %f ...\n",upper_beta_thresh, lower_beta_thresh);
        }
        else if(einfo->_valType==1) {
            LOGPRINTF("The data input are of methylation m-values.\n Transfering beta thresholds to m thresholds...\n");
            upper_beta_thresh=log2(upper_beta_thresh/(1-upper_beta_thresh));
            lower_beta_thresh=log2(lower_beta_thresh/(1-lower_beta_thresh));
            LOGPRINTF("Excluding constitutively methylated/unmethylated sites with upper m threshold %f and lower m threshold %f ...\n",upper_beta_thresh, lower_beta_thresh);
            
        } else {
             LOGPRINTF("ERROR: No methylation data found.\n");
            exit(EXIT_FAILURE);
        }
        long n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        uint64_t ttl_n = einfo->_eii_num;
        vector<int> keep;
        for(int j=0; j<m; j++){
            bool rmflag=false;
            double mu=0.0;
            double nonmiss=0.0;
            for(int i=0; i<n; i++){
                double val=einfo->_val[einfo->_epi_include[j]*ttl_n+einfo->_eii_include[i]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            if(mu>=upper_beta_thresh || mu<=lower_beta_thresh) rmflag=true;
            
            if(!rmflag) keep.push_back(einfo->_epi_include[j]);
        }
        einfo->_epi_map.clear();
        for(int i=0;i<keep.size();i++)
            einfo->_epi_map.insert(pair<string,int>(einfo->_epi_prb[keep[i]],keep[i]));
        einfo->_epi_include.swap(keep);
        
        LOGPRINTF("%ld of %ld sites are Excluded with upper beta threshold %f and lower beta threshold %f and %ld sites left.\n",m-einfo->_epi_include.size(),m,upper_beta_thresh,lower_beta_thresh,einfo->_epi_include.size());
    }
    /* QC step: remove the probe which has one or more detection pvlaues voilate the threshold*/
    void filtering_with_detpval(eInfo* einfo, char* dpvalfName, double dp_thresh, double prb_thresh, double spl_thresh,int mth, bool no_fid_flag)
    {
        // mth ==0 : not to remove individuals. to remove the probes that has one or more detection pvlaues voilate the threshold ( default as 0.05)
        // mth ==1 : drop samples that failed in 1% (default) of probes. drop probes that are unsuccessfully measured in 1% (default) samples.
        if(einfo->_eType!=METHYLATION)
        {
            LOGPRINTF("ERROR: No methylation data found.\n");
            exit(EXIT_FAILURE);
        }
        FILE* pfile=NULL;  //each row of this file is for a probe. each column is for an individual.
   
        vector<string> strlist;
        uint32_t line_idx = 0;
        char* readBuf;
        string msg="Reading Buffer";
        if(!allocReserved(&readBuf, MAXLINEBUFLEN,msg)) exit(EXIT_FAILURE);
        if(fopen_checked(&pfile, dpvalfName,"r")) exit(EXIT_FAILURE);
        
        LOGPRINTF("Reading detection p-values from %s ...\n",dpvalfName);
        vector<string> fid;
        vector<string> iid;
        if(no_fid_flag)
        {
            if(fgets(readBuf, MAXLINEBUFLEN, pfile))
            {
                split_str(readBuf,strlist,1);
                fid.swap(strlist);
                iid=fid;
               
            }
        } else {
            if(fgets(readBuf, MAXLINEBUFLEN, pfile))
            {
                split_str(readBuf,strlist,1);
                fid.swap(strlist);
            }
            memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
            if(fgets(readBuf, MAXLINEBUFLEN, pfile))
            {
                split_str(readBuf,strlist,1);
                iid.swap(strlist);
            }
        }
        vector<string> pfiid;
        map<string,int> eiimap;
        for(int i=0;i<fid.size();i++) {
            pfiid.push_back(fid[i]+":"+iid[i]);
            eiimap.insert(pair<string,int>(fid[i]+":"+iid[i],i));
            if(eiimap.size()==i)
            {
                LOGPRINTF("ERROR: Duplicate individual with FID: %s and IID %s.\n", fid[i].c_str(),iid[i].c_str());
                exit(EXIT_FAILURE);
            }

        }
        vector<string> bfiid;
        for(int i=0;i<einfo->_eii_include.size();i++) bfiid.push_back(einfo->_eii_fid[einfo->_eii_include[i]]+":"+einfo->_eii_fid[einfo->_eii_include[i]]);
        vector<int> idx;
        match_only(bfiid, pfiid, idx);
        if(idx.size()==0)
        {
            LOGPRINTF("ERROR: No individuals in common between detection p-value file and beed file.\n");
            exit(EXIT_FAILURE);
        }
        if(idx.size()<0.2*bfiid.size() && idx.size()<0.2*pfiid.size())
        {
            LOGPRINTF("WARNING: Low overlap of individuals between detection p-value file and beed file.\n");
        }
        LOGPRINTF("Total %ld individuals are matched between methylation data (%ld individuals) and detection p-value file (%ld individuals).\n",idx.size(),bfiid.size(),pfiid.size());
        if(idx.size() != bfiid.size()) { LOGPRINTF("Only detection p-values of common individuals are taken into account.\n")};
        memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
        line_idx=0;
        vector<string> pid_rm;
        map<string, int>::iterator iter;
        vector<int> spl_count(idx.size(),0), prb_count(einfo->_epi_num,0); // It is hard to determine the idx of _eip_include when _epi_inlcude.size()<_epi_num. so _epi_num is used here.
        int probe_match_count=0;
        while(fgets(readBuf, MAXLINEBUFLEN, pfile))
        {
            int curpIdx=-9;
            split_str(readBuf,strlist,0);
            if(strlist.size()!=1+pfiid.size())
            {
                LOGPRINTF("ERROR: Inconsistent element number found in line %u with probe ID %s .\n", line_idx+1, strlist[0].c_str());
                exit(EXIT_FAILURE);
            }
            iter = einfo->_epi_map.find(strlist[0].c_str());
            if (iter != einfo->_epi_map.end()){
                probe_match_count++;
                curpIdx=iter->second;
                bool hit=false;
                for(int i=0;i<idx.size();i++)
                {
                    if(strlist[idx[i]]=="NA" || strlist[idx[i]]=="." || strlist[idx[i]]=="-" || strlist[idx[i]]=="-9") {
                        LOGPRINTF("ERROR: Detection p-value is missing in line %u with probe ID %s .\n", line_idx+1, strlist[0].c_str());
                        exit(EXIT_FAILURE);
                    } else {
                        double valtmp=atof(strlist[idx[i]].c_str());
                        if( valtmp<0 || valtmp >1 ) {
                            LOGPRINTF("ERROR: Detection p-value value should be between 0 and 1 in line %u.\n", line_idx+1 );
                            exit(EXIT_FAILURE);
                        } else {
                            if(valtmp>dp_thresh) {
                                if(!hit){
                                    pid_rm.push_back(einfo->_epi_prb[curpIdx]);
                                    hit=true;
                                }
                                prb_count[curpIdx]++;
                                spl_count[i]++; // use pfiid[idx[]] to get fiid
                            }
                        }
                    }
                }
            }
            line_idx++;
            memset(readBuf,0,sizeof(char)*(MAXLINEBUFLEN));
        }
        LOGPRINTF("Total %u probes have been read from file %s and %d probes were matched with .epi file.\n", line_idx,  dpvalfName,probe_match_count);
        fclose(pfile);
        deallocReserved(&readBuf, MAXLINEBUFLEN);
        if(mth==0)
        {
            if(pid_rm.size()>0) update_map_rm(pid_rm, einfo->_epi_map, einfo->_epi_include);
            LOGPRINTF("%ld probes have been excluded and %ld probes remained.\n", pid_rm.size(),  einfo->_epi_include.size());
        } else if( mth==1) {
            int maxnum=ceil(prb_thresh*einfo->_eii_include.size());
            vector<string> pid_rm;
            for(int i=0;i<prb_count.size();i++){
                if(prb_count[i]>maxnum) pid_rm.push_back(einfo->_epi_prb[i]);
            }
            update_map_rm(pid_rm, einfo->_epi_map, einfo->_epi_include);
            
            maxnum=ceil(spl_thresh*einfo->_epi_include.size());
            vector<string> iid_rm;
            for(int i=0;i<spl_count.size();i++){
                if(spl_count[i]>maxnum) iid_rm.push_back(pfiid[idx[i]]);
            }
             update_map_rm(iid_rm, einfo->_eii_map, einfo->_eii_include);
            LOGPRINTF("%ld probes have been excluded and %ld probes remained. \n", pid_rm.size(),  einfo->_epi_include.size());
            LOGPRINTF("%ld individuals have been excluded and %ld individuals remained. \n", iid_rm.size(),  einfo->_eii_include.size());
        }
        
    }
    void std_probe_ind(eInfo* einfo, vector< vector<bool> > &X_bool, bool divid_by_std, MatrixXd &_probe_data, bool output_profile)
    {
        int i = 0, j = 0;
        uint64_t n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        Map<MatrixXd> X(&einfo->_val[0],n,m);
        _probe_data.resize(n, m);
        
#pragma omp parallel for private(j)
        for(i=0; i<n; i++){
            for(j=0; j<m; j++) _probe_data(i,j)=X(einfo->_eii_include[i], einfo->_epi_include[j]);
        }
        X.resize(0,0);
        
        VectorXd mu(n), nonmiss(n);
        
        X_bool.resize(n);
        for(i=0; i<n; i++) X_bool[i].resize(m);
#pragma omp parallel for private(j)
        for(i=0; i<n; i++){
            mu(i)=0.0;
            nonmiss(i)=0.0;
            for(j=0; j<m; j++){
                if(_probe_data(i,j)<1e9){
                    mu(i)+=_probe_data(i,j);
                    nonmiss(i)+=1.0;
                    X_bool[i][j] = true;
                }
                else X_bool[i][j] = false;
            }
            mu(i)/=nonmiss(i);
        }
        
#pragma omp parallel for private(j)
        for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                if(_probe_data(i,j)<1e9) _probe_data(i,j) -= mu(i);
                else _probe_data(i,j) = 0.0;
            }
        }
        if(output_profile)
        {
            einfo->_profile=new double[n*m]; // alloc memory
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    einfo->_profile[i*m+j]=_probe_data(i,j);
                }
            }
            
        }
        if(divid_by_std){
            VectorXd sd(n);
#pragma omp parallel for
            for(i=0; i<n; i++){
                sd(i) = sqrt((_probe_data.row(i).dot(_probe_data.row(i))) / (nonmiss(i) - 1.0));
            }
            
#pragma omp parallel for private(j)
            for(i=0; i<n; i++){
                for(j=0; j<m; j++) _probe_data(i,j) /= sd(i);
            }
        }
    }
    /* QC step: remove the probe with missing rate over a threshold*/
    void filtering_with_missingratio(eInfo* einfo,double missratioprobe)
    {
        LOGPRINTF("Excluding probes with missing rate over %f ...\n",missratioprobe);
        long n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        uint64_t ttl_n = einfo->_eii_num;
        vector<int> keep;
        vector<string> rmprb;
        for(int j=0; j<m; j++){
            double nmiss=0.0;
            for(int i=0; i<n; i++){
                double val=einfo->_val[einfo->_epi_include[j]*ttl_n+einfo->_eii_include[i]];
                if(val>=1e9){
                    nmiss+=1.0;
                }
            }
            
            double mrate=nmiss/n;
            if(mrate<=missratioprobe) keep.push_back(einfo->_epi_include[j]);
            else rmprb.push_back(einfo->_epi_prb[einfo->_epi_include[j]]);
        }
        if(rmprb.size()>0)
        {
            char prb_oname[FNAMESIZE];
            memcpy(prb_oname,"rm.probe.list",14);
            write_msglist(prb_oname,rmprb);
            LOGPRINTF("%ld probes to remove have been saved in file %s.\n",rmprb.size(),prb_oname);
        }
        einfo->_epi_map.clear();
        for(int i=0;i<keep.size();i++)
            einfo->_epi_map.insert(pair<string,int>(einfo->_epi_prb[keep[i]],keep[i]));
        einfo->_epi_include.swap(keep);
        
        LOGPRINTF("%ld of %ld sites are Excluded with a missing value threshold (%f) and %ld sites left.\n",m-einfo->_epi_include.size(),m,missratioprobe,einfo->_epi_include.size());
    }
    /***********************************
    * this file is with no header
     mpheno is a list of comma-delimited trait numbers. e.g. "1,3".
     The input trait number should be >=1.
     mvFlg indicates mulitple variates analysis or not.
     if mpheno is null and mvFlg is false, the first trait is to be read.
     if mpheno is null and mvFlg is true, all traits are to be read.
     either included phenotype of the individual is missing, this 
     individual is deemed as missning individual.
     */
    void read_phen(eInfo* einfo, string phen_file, char* mpheno, bool mvFlg) {
        
        if(einfo->_eii_include.size()==0 )
        {
            LOGPRINTF("Error: no individual found in gene expression / methylation data.\n");
            exit(EXIT_FAILURE);
        }
        // Read phenotype data
        ifstream in_phen(phen_file.c_str());
        if (!in_phen)
        {
            LOGPRINTF("Error: can not open the file %s to read.\n",phen_file.c_str());
            exit(EXIT_FAILURE);
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
            exit(EXIT_FAILURE);
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
                    exit(EXIT_FAILURE);
                }
                if(curid>=phen_num)
                {
                    LOGPRINTF("Error:  can not find the %dth trait in the file %s.\n",curid+1,phen_file.c_str());
                    exit(EXIT_FAILURE);
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
                exit(EXIT_FAILURE);
            }
        }
        
        einfo->_eii_pheno_num=(int)pheno_ids.size();
        if(!sudoAllocReserved( (pheno_ids.size()-1)*einfo->_eii_num*sizeof(double), "phenotype")) exit(EXIT_FAILURE);
        einfo->_eii_pheno.resize(pheno_ids.size()*einfo->_eii_num);
        LOGPRINTF("There are %d traits in the file %s and %ld trait(s) specified. \n",phen_num, phen_file.c_str(),pheno_ids.size());
        for(int i=0;i<pheno_ids.size()*einfo->_eii_num;i++) einfo->_eii_pheno[i]=-9;
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
                 exit(EXIT_FAILURE);
            }
            iter = einfo->_eii_map.find(fid_buf + ":" + iid_buf);
            if (iter != einfo->_eii_map.end())
            {
                int idx = iter->second;
                for(int i=0;i<pheno_ids.size();i++)
                {
                    if (vs_buf[pheno_ids[i]] != "-9" && vs_buf[pheno_ids[i]] != "NA")
                    {
                        einfo->_eii_pheno[i*einfo->_eii_num+idx]=atof(vs_buf[pheno_ids[i]].c_str());
                    }
                }
            }
        }
        in_phen.close();
        
        for(int i=0;i<pheno_ids.size();i++)
            for(int j=0;j<einfo->_eii_num;j++)
                if(einfo->_eii_pheno[i*einfo->_eii_num+j]==-9) missindi.insert(pair<string,int>(einfo->_eii_fid[j] + ":" + einfo->_eii_iid[j],j));
        
        for (iter = missindi.begin(); iter != missindi.end(); iter++) indi_list.push_back(iter->first);        
        update_map_rm(indi_list, einfo->_eii_map, einfo->_eii_include);
        if(einfo->_eii_include.size()==0)
        {
            LOGPRINTF("Error: no individual is included.\n");
            exit(EXIT_FAILURE);
        }
        LOGPRINTF("Non-missing phenotypes of %ld individuals are included from %s. \n",einfo->_eii_include.size(), phen_file.c_str());
    }
    
    void read_cov(eInfo* einfo, string cov_file, bool qcovFlg) {
        
        if(einfo->_eii_include.size()==0 )
        {
            LOGPRINTF("Error: no individual found in gene expression / methylation data.\n");
            exit(EXIT_FAILURE);
        }
        // Read cov data
        ifstream in_cov(cov_file.c_str());
        if (!in_cov)
        {
            LOGPRINTF("Error: can not open the file %s to read.\n",cov_file.c_str());
            exit(EXIT_FAILURE);
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
            exit(EXIT_FAILURE);
        }
        
        if(qcovFlg)
        {
            einfo->_eii_qcov_num=cov_num;
            if(!sudoAllocReserved( cov_num*einfo->_eii_num*sizeof(double), covtypestr)) exit(EXIT_FAILURE);
            einfo->_eii_qcov.resize(cov_num*einfo->_eii_num);
            LOGPRINTF("There are %d %s covariate(s) in the file %s. \n",cov_num,covtypestr.c_str(), cov_file.c_str());
            for(int i=0;i<cov_num*einfo->_eii_num;i++) einfo->_eii_qcov[i]=-9;
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
                    exit(EXIT_FAILURE);
                }
                iter = einfo->_eii_map.find(fid_buf + ":" + iid_buf);
                if (iter != einfo->_eii_map.end())
                {
                    int idx = iter->second;
                    for(int i=0;i<cov_num;i++)
                    {
                        if (vs_buf[i] != "-9" && vs_buf[i] != "NA")
                        {
                            einfo->_eii_qcov[i*einfo->_eii_num+idx]=atof(vs_buf[i].c_str());
                        }
                    }
                }
            }
            for(int i=0;i<cov_num;i++)
                for(int j=0;j<einfo->_eii_num;j++)
                    if(einfo->_eii_qcov[i*einfo->_eii_num+j]==-9) missindi.insert(pair<string,int>(einfo->_eii_fid[j] + ":" + einfo->_eii_iid[j],j));
            for (iter = missindi.begin(); iter != missindi.end(); iter++) indi_list.push_back(iter->first);
            update_map_rm(indi_list, einfo->_eii_map, einfo->_eii_include);
        }
        else
        {
            einfo->_eii_cov_num=cov_num;
            if(!sudoAllocReserved( cov_num*einfo->_eii_num*sizeof(double), covtypestr)) exit(EXIT_FAILURE);
            einfo->_eii_cov.resize(cov_num*einfo->_eii_num);
            LOGPRINTF("There are %d %s covariate(s) in the file %s. \n",cov_num, covtypestr.c_str(),cov_file.c_str());
            for(int i=0;i<cov_num*einfo->_eii_num;i++) einfo->_eii_cov[i]="-9";
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
                    exit(EXIT_FAILURE);
                }
                iter = einfo->_eii_map.find(fid_buf + ":" + iid_buf);
                if (iter != einfo->_eii_map.end())
                {
                    int idx = iter->second;
                    for(int i=0;i<cov_num;i++)
                    {
                        if (vs_buf[i] != "-9" && vs_buf[i] != "NA")
                        {
                            einfo->_eii_cov[i*einfo->_eii_num+idx]=vs_buf[i];
                        }
                    }
                }
            }
            for(int i=0;i<cov_num;i++)
                for(int j=0;j<einfo->_eii_num;j++)
                    if(einfo->_eii_cov[i*einfo->_eii_num+j]=="-9") missindi.insert(pair<string,int>(einfo->_eii_fid[j] + ":" + einfo->_eii_iid[j],j));
            for (iter = missindi.begin(); iter != missindi.end(); iter++) indi_list.push_back(iter->first);
            update_map_rm(indi_list, einfo->_eii_map, einfo->_eii_include);

        }
        in_cov.close();
        if(einfo->_eii_include.size()==0)
        {
            LOGPRINTF("Error: no individual is included.\n");
            exit(EXIT_FAILURE);
        }
        LOGPRINTF("Non-missing %s covariates of %ld individuals are included from %s. \n",covtypestr.c_str(),einfo->_eii_include.size(), cov_file.c_str());
    }
    
    void output_grm(eInfo* einfo, string _out, bool output_grm_bin)
    {
        int i = 0, j = 0;
        string grm_file;
        if (output_grm_bin) {
            // Save matrix A in binary file
            grm_file = _out + ".grm.bin";
            fstream A_Bin(grm_file.c_str(), ios::out | ios::binary);
            if (!A_Bin) throw ("Error: can not open the file [" + grm_file + "] to write.");
            float f_buf = 0.0;
            int size = sizeof (float);
            for (i = 0; i < einfo->_eii_include.size(); i++) {
                for (j = 0; j <= i; j++) {
                    f_buf = (float) ( einfo->_grm(i, j));
                    A_Bin.write((char*) &f_buf, size);
                }
            }
            A_Bin.close();
            cout << "ERM of " << einfo->_eii_include.size() << " individuals has been saved in the file [" + grm_file + "] (in binary format)." << endl;
            
            string grm_N_file = _out + ".grm.N.bin";
            fstream N_Bin(grm_N_file.c_str(), ios::out | ios::binary);
            if (!N_Bin) throw ("Error: can not open the file [" + grm_N_file + "] to write.");
            f_buf = 0.0;
            size = sizeof (int);
            for (i = 0; i < einfo->_eii_include.size(); i++) {
                for (j = 0; j <= i; j++) {
                    f_buf = (float) (einfo->_grm_N(i, j));
                    N_Bin.write((char*) &f_buf, size);
                }
            }
            N_Bin.close();
            cout << "Number of SNPs to calcuate the genetic relationship between each pair of individuals has been saved in the file [" + grm_N_file + "] (in binary format)." << endl;
        }
        else {
            // Save A matrix in txt format
            grm_file = _out + ".grm.gz";
            gzFile zoutf = gzopen(grm_file.c_str(),"wb");
            if (!(zoutf)) {
                LOGPRINTF("Error: Couldn't open file %s to write.\n", grm_file.c_str());
                exit (EXIT_FAILURE);
            }
            cout << "Saving the genetic relationship matrix to the file [" + grm_file + "] (in compressed text format)." << endl;
            for (i = 0; i < einfo->_eii_include.size(); i++) {
                 string tmpstr="";
                if (einfo->_grm_N.rows() > 0){
                    for (j = 0; j <= i; j++) {
                        tmpstr=atos(i + 1)+'\t'+atos(j + 1)+'\t'+dtos(einfo->_grm_N(i, j))+'\t'+ dtos(einfo->_grm(i, j))+'\n';
                        if(gzputs(zoutf, tmpstr.c_str()) == -1) cout<<"error"<<endl;
                    }
                }
                else{
                    for (j = 0; j <= i; j++) {
                        tmpstr=atos(i + 1)+'\t'+atos(j + 1)+'\t'+atos(0)+'\t' + dtos(einfo->_grm(i, j))+'\n';
                        if(gzputs(zoutf, tmpstr.c_str()) == -1) cout<<"error"<<endl;
                    }
                }
                
            }
            gzclose(zoutf);
            cout << "The genetic relationship matrix has been saved in the file [" + grm_file + "] (in compressed text format)." << endl;
        }
        
        string famfile = _out + ".grm.id";
        ofstream Fam(famfile.c_str());
        if (!Fam) throw ("Error: can not open the file [" + famfile + "] to write.");
        for (i = 0; i < einfo->_eii_include.size(); i++) Fam << einfo->_eii_fid[einfo->_eii_include[i]] + "\t" + einfo->_eii_iid[einfo->_eii_include[i]] << endl;
        Fam.close();
        cout << "IDs for the GRM file [" + grm_file + "] have been saved in the file [" + famfile + "]." << endl;
    }

    
    void make_erm(eInfo* einfo, int erm_mtd, bool output_bin, char* outFileName, bool output_profile)
    {
        uint64_t n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        
        MatrixXd _probe_data;
        string typestr=getFileType(einfo->_eType);
        LOGPRINTF("Recoding %s data...\n",typestr.c_str());
      
        bool divid_by_std = false;
        vector< vector<bool> > X_bool;
        if(erm_mtd < 2){
            if(erm_mtd == 0) divid_by_std = true;
            else if(erm_mtd == 1) divid_by_std = false;
            std_probe(einfo, X_bool, divid_by_std, _probe_data, output_profile);
        }
        else std_probe_ind(einfo, X_bool, false, _probe_data, output_profile);
        /*
        FILE* tmpfile=fopen("osc.txt","w");
        if(!tmpfile)
        {
            printf("error open file.\n");
            exit(EXIT_FAILURE);
        }
        printf("\n**************\n");
        for(int i=0;i<_probe_data.rows();i++) {
            string str=atos(_probe_data(i,0));
            for(int j=1;j<_probe_data.cols();j++)
            {
                str+='\t'+atos(_probe_data(i,j));
            }
            str+='\n';
            fputs(str.c_str(),tmpfile);
        }
        fclose(tmpfile);
        */
        
        cout << "\nCalculating expression relationship matrix (ERM) ... " << endl;
        // count the number of missing genotypes
        vector< vector<int> > miss_pos(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (X_bool[i][j] == false) miss_pos[i].push_back(j);
            }
        }
        
        // Calculate A_N matrix
        einfo->_grm_N.resize(n, n);
        if(erm_mtd == 0){
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= i; j++) {
                    int miss_j = 0;
                    for (int k = 0; k < miss_pos[j].size(); k++) miss_j += (int) X_bool[i][miss_pos[j][k]];
                    einfo->_grm_N(i,j) = m - miss_pos[i].size() - miss_j;
                }
            }
        }
        else if (erm_mtd > 0){
            if (erm_mtd == 1){
                VectorXd nonmiss(m);
                for (int j = 0; j < m; j++) {
                    nonmiss[j] = 0.0;
                    for (int i = 0; i < n; i++) {
                        if (X_bool[i][j] == true) nonmiss[j] += 1.0;
                    }
                }
                VectorXd var(m);
                for(int j=0; j<m; j++){
                    if((nonmiss(j) - 1.0) > 0.0) var(j) = (_probe_data.col(j).dot(_probe_data.col(j))) / (nonmiss(j) - 1.0);
                    else var(j) = 0.0;
                    
                }
                double sum_var = var.sum();
                for (int i = 0; i < n; i++) {
                    double i_miss_sum_var = 0.0;
                    for (int k = 0; k < miss_pos[i].size(); k++) i_miss_sum_var += var[miss_pos[i][k]];
                    for (int j = 0; j <= i; j++) {
                        double j_miss_sum_var = 0.0;
                        for (int k = 0; k < miss_pos[j].size(); k++){
                            if (X_bool[i][miss_pos[j][k]] == true) j_miss_sum_var += var[miss_pos[j][k]];
                        }
                        einfo->_grm_N(i,j) = sum_var - i_miss_sum_var - j_miss_sum_var;
                    }
                }
            }
            else{
                VectorXd ssx(n);
                for(int i=0; i<n; i++) ssx(i) = _probe_data.row(i).dot(_probe_data.row(i));
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j <= i; j++) {
                        double ssx_i = ssx(i);
                        double ssx_j = ssx(j);
                        for (int k = 0; k < miss_pos[j].size(); k++){
                            int l = miss_pos[j][k];
                            if(X_bool[i][l] == true) ssx_i -= _probe_data(i,l) * _probe_data(i,l);
                        }
                        for (int k = 0; k < miss_pos[i].size(); k++){
                            int l = miss_pos[i][k];
                            if(X_bool[j][l] == true) ssx_j -= _probe_data(j,l) * _probe_data(j,l);
                        }
                        einfo->_grm_N(i,j) = sqrt(ssx_i*ssx_j);
                    }
                }
            }
        }
        // Calculate A matrix
        einfo->_grm = _probe_data * _probe_data.transpose();
    
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                if (einfo->_grm_N(i,j) > 0.0) einfo->_grm(i,j) /= einfo->_grm_N(i,j);
                else einfo->_grm(i,j) = 0.0;
            }
        }
         einfo->_grm = einfo->_grm.array() / einfo->_grm.diagonal().mean();
        // Output A_N and A
        if(outFileName!=NULL) output_grm(einfo, outFileName ,output_bin);
       
        if(output_profile)
        {
            einfo->_grm_ptr=new double[n*n]; // alloc memory
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    einfo->_grm_ptr[i*n+j]=einfo->_grm(i,j);
                }
            }
            
        }
    }
    
    void read_grm(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only, bool dont_read_N, bool grm_bin_flag)
    {
        if (grm_bin_flag) read_grm_bin(einfo, grm_file, grm_id, out_id_log, read_id_only, dont_read_N);
        else read_grm_gz(einfo, grm_file, grm_id, out_id_log, read_id_only);
    }
    
    int read_grm_id(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log)
    {
        // read GRM IDs
        string grm_id_file = grm_file + ".grm.id";
        if (out_id_log) cout << "Reading IDs of the ERM from [" + grm_id_file + "]." << endl;
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
        int n = grm_id.size();
        if (out_id_log) cout << n << " IDs read from [" + grm_id_file + "]." << endl;
        
        if (einfo->_eii_map.empty()) {
            einfo->_eii_fid = fid;
            einfo->_eii_iid = pid;
            einfo->_eii_num =(int)fid.size();
            einfo->_eii_sex.resize(fid.size());
            einfo->_eii_include.clear();
            einfo->_eii_include.resize(einfo->_eii_num);
            int i = 0, size = 0;
            for (i = 0; i < einfo->_eii_num; i++) {
                einfo->_eii_include[i] = i;
                einfo->_eii_map.insert(pair<string, int>(einfo->_eii_fid[i] + ":" + einfo->_eii_iid[i], i));
                if (size ==  einfo->_eii_map.size()) throw ("Error: Duplicate individual ID found: \"" + einfo->_eii_fid[i] + "\t" + einfo->_eii_iid[i] + "\".");
                size = einfo->_eii_map.size();
            }
        }else{
            update_map_kp(grm_id, einfo->_eii_map, einfo->_eii_include);
        }
        
        return (n);
    }

    void read_grm_gz(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only) {
        int n = read_grm_id(einfo,grm_file, grm_id, out_id_log);
        if (read_id_only) return;
        
        string grm_gzfile = grm_file + ".grm.gz", str_buf;
        char buf[MAX_LINE_SIZE];
        gzFile zinf= gzopen(grm_gzfile.c_str(), "rb");
      
        if (!(zinf)) {
            LOGPRINTF("ERROR: Couldn't open file %s\n", grm_gzfile.c_str());
            exit (EXIT_FAILURE);
        }
        
        int indx1 = 0, indx2 = 0, nline = 0;
        double grm_buf = 0.0, grm_N_buf;
        string errmsg = "Error: failed to read [" + grm_gzfile + "]. The format of the ERM file has been changed?\nError occurs in line:\n";
        cout << "Reading the GRM from [" + grm_gzfile + "]." << endl;
        einfo->_grm.resize(n, n);
        einfo->_grm_N.resize(n, n);
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
                einfo->_grm_N(indx1 - 1, indx2 - 1) = einfo->_grm_N(indx2 - 1, indx1 - 1) = grm_N_buf;
                einfo->_grm(indx1 - 1, indx2 - 1) = einfo->_grm(indx2 - 1, indx1 - 1) = grm_buf;
                nline++;
                if (ss >> str_buf) throw (errmsg + buf);
            }
        }
        gzclose(zinf);
        if (!einfo->_within_family && nline != (int) (n * (n + 1)*0.5)){
            stringstream errmsg_tmp;
            errmsg_tmp << "Error: there are " << nline << " lines in the [" << grm_gzfile << "] file. The expected number of lines is " << (int) (n * (n + 1)*0.5) << "." << endl;
            throw(errmsg_tmp.str());
            //throw ("Error: incorrect number of lines in the grm file. *.grm.gz file and *.grm.id file are mismatched?");
        }
        cout << "ERM for " << n << " individuals are included from [" + grm_gzfile + "]." << endl;
    }
    
    void read_grm_bin(eInfo* einfo, string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only, bool dont_read_N)
    {
        int i = 0, j = 0, n = read_grm_id(einfo,grm_file, grm_id, out_id_log);
        
        if (read_id_only) return;
        
        string grm_binfile = grm_file + ".grm.bin";
        ifstream A_bin(grm_binfile.c_str(), ios::in | ios::binary);
        if (!A_bin.is_open()) throw ("Error: can not open the file [" + grm_binfile + "] to read.");
        einfo->_grm.resize(n, n);
        cout << "Reading the GRM from [" + grm_binfile + "]." << endl;
        int size = sizeof (float);
        float f_buf = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                if (!(A_bin.read((char*) &f_buf, size))) throw ("Error: the size of the [" + grm_binfile + "] file is incomplete?");
                einfo->_grm(j, i) = einfo->_grm(i, j) = f_buf;
            }
        }
        A_bin.close();
        
        if(!dont_read_N){
            string grm_Nfile = grm_file + ".grm.N.bin";
            ifstream N_bin(grm_Nfile.c_str(), ios::in | ios::binary);
            if (!N_bin.is_open()) throw ("Error: can not open the file [" + grm_Nfile + "] to read.");
            einfo->_grm_N.resize(n, n);
            cout << "Reading the number of SNPs for the GRM from [" + grm_Nfile + "]." << endl;
            size = sizeof (float);
            f_buf = 0.0;
            for (i = 0; i < n; i++) {
                for (j = 0; j <= i; j++) {
                    if (!(N_bin.read((char*) &f_buf, size))) throw ("Error: the size of the [" + grm_Nfile + "] file is incomplete?");
                    einfo->_grm_N(j, i) = einfo->_grm_N(i, j) = f_buf;
                }
            }
            N_bin.close();
        }
        
        cout << "ERM for " << n << " individuals are included from [" + grm_binfile + "]." << endl;
    }
    void merge_grm(eInfo* einfo, char* merge_grm_file) {
        vector<string> grm_files, grm_id;
        read_msglist(merge_grm_file, grm_files, "ERM file");
        
        
        for (int f = 0; f < grm_files.size(); f++) {
            read_grm(einfo, grm_files[f], grm_id, false, true, false,true);
        }
        
        vector<string> uni_id;
        map<string, int> uni_id_map;
        map<string, int>::iterator iter;
        for(int i=0; i<einfo->_eii_include.size(); i++){
            uni_id.push_back(einfo->_eii_fid[einfo->_eii_include[i]]+":"+einfo->_eii_iid[einfo->_eii_include[i]]);
            uni_id_map.insert(pair<string,int>(einfo->_eii_fid[einfo->_eii_include[i]]+":"+einfo->_eii_iid[einfo->_eii_include[i]], i));
        }
        
        int _n=(int)einfo->_eii_include.size();
        if(_n<1) {
            LOGPRINTF("ERROR: no individual is in common in the input files.\n");
            exit(EXIT_FAILURE);
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        
        vector<int> kp;
        MatrixXd grm = MatrixXd::Zero(_n, _n);
        MatrixXd grm_N = MatrixXd::Zero(_n, _n);
        for (int f = 0; f < grm_files.size(); f++) {
            cout << "Reading the ERM from the " << f + 1 << "th file ..." << endl;
            read_grm(einfo, grm_files[f], grm_id,true,false,false,true);
            match(uni_id, grm_id, kp);
            for (int i = 0; i < _n; i++) {
                for (int j = 0; j <= i; j++) {
                    if (kp[i] >= kp[j]) {
                        grm(i, j) += einfo->_grm(kp[i], kp[j]) * einfo->_grm_N(kp[i], kp[j]);
                        grm_N(i, j) += einfo->_grm_N(kp[i], kp[j]);
                    } else {
                        grm(i, j) += einfo->_grm(kp[j], kp[i]) * einfo->_grm_N(kp[j], kp[i]);
                        grm_N(i, j) += einfo->_grm_N(kp[j], kp[i]);
                    }
                }
            }
        }
        for (int i = 0; i < _n; i++) {
            for (int j = 0; j <= i; j++) {
                if (grm_N(i, j) == 0) einfo->_grm(i, j) = 0;
                else einfo->_grm(i, j) = grm(i, j) / grm_N(i, j);
                einfo->_grm_N(i, j) = grm_N(i, j);
            }
        }
        grm.resize(0, 0);
        grm_N.resize(0, 0);
        cout << "\n" << grm_files.size() << " ERMs have been merged together." << endl;
    }
    
    void rm_cor_indi(eInfo* einfo, double grm_cutoff) {
        cout << "Pruning the ERM with a cutoff of " << grm_cutoff << " ..." << endl;
        
        int i = 0, j = 0, i_buf = 0;
        
        // identify the positions where you see a value > than the threshold
        vector<int> rm_grm_ID1, rm_grm_ID2;
        for (i = 0; i < einfo->_eii_include.size(); i++) {
            for (j = 0; j < i; j++) {
                if (einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[j]) > grm_cutoff) {
                    rm_grm_ID1.push_back(einfo->_eii_include[i]);
                    rm_grm_ID2.push_back(einfo->_eii_include[j]);
                }
            }
        }
        
        // count the number of appearance of each "position" in the vector, which involves a few steps
        vector<int> rm_uni_ID(rm_grm_ID1);
        rm_uni_ID.insert(rm_uni_ID.end(), rm_grm_ID2.begin(), rm_grm_ID2.end());
        stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
        rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
        map<int, int> rm_uni_ID_count;
        for (i = 0; i < rm_uni_ID.size(); i++) {
            i_buf = count(rm_grm_ID1.begin(), rm_grm_ID1.end(), rm_uni_ID[i]) + count(rm_grm_ID2.begin(), rm_grm_ID2.end(), rm_uni_ID[i]);
            rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
        }
        
        // swapping
        map<int, int>::iterator iter1, iter2;
        for (i = 0; i < rm_grm_ID1.size(); i++) {
            iter1 = rm_uni_ID_count.find(rm_grm_ID1[i]);
            iter2 = rm_uni_ID_count.find(rm_grm_ID2[i]);
            if (iter1->second < iter2->second) {
                i_buf = rm_grm_ID1[i];
                rm_grm_ID1[i] = rm_grm_ID2[i];
                rm_grm_ID2[i] = i_buf;
            }
        }
        
        stable_sort(rm_grm_ID1.begin(), rm_grm_ID1.end());
        rm_grm_ID1.erase(unique(rm_grm_ID1.begin(), rm_grm_ID1.end()), rm_grm_ID1.end());
        vector<string> removed_ID;
        for (i = 0; i < rm_grm_ID1.size(); i++) removed_ID.push_back(einfo->_eii_fid[rm_grm_ID1[i]] + ":" + einfo->_eii_iid[rm_grm_ID1[i]]);
        
        // update _keep and _id_map
        update_map_rm(removed_ID, einfo->_eii_map, einfo->_eii_include);
        
        cout << "After pruning the ERM, there are " << einfo->_eii_include.size() << " individuals (" << removed_ID.size() << " individuals removed)." << endl;
    }
    
    void update_sex(eInfo* einfo, char* sex_file) {
        ifstream isex(sex_file);
        if (!isex) throw ("Error: can not open the file [" + string(sex_file) + "] to read.");
        int sex_buf = 0, icount = 0;
        string str_buf, fid, pid;
        cout << "Reading sex information from [" + string(sex_file) + "]." << endl;
        map<string, int>::iterator iter, End = einfo->_eii_map.end();
        einfo->_eii_sex.clear();
        einfo->_eii_sex.resize(einfo->_eii_num);
        vector<int> confirm(einfo->_eii_num);
        while (isex) {
            isex >> fid;
            if (isex.eof()) break;
            isex >> pid;
            isex >> str_buf;
            if (str_buf != "1" && str_buf != "2" && str_buf != "M" && str_buf != "F") throw ("Error: unrecognized sex code: \"" + fid + " " + pid + " " + str_buf + "\" in [" + sex_file + "].");
            iter = einfo->_eii_map.find(fid + ":" + pid);
            if (iter != End) {
                if (str_buf == "M" || str_buf == "1") einfo->_eii_sex[iter->second] = 1;
                else if (str_buf == "F" || str_buf == "2") einfo->_eii_sex[iter->second] = 2;
                confirm[iter->second] = 1;
                icount++;
            }
            getline(isex, str_buf);
        }
        isex.close();
        
        for (int i = 0; i < einfo->_eii_include.size(); i++) {
            if (confirm[einfo->_eii_include[i]] != 1) throw ("Error: sex information for all of the included individuals should be updated.");
        }
        cout << "Sex information for " << icount << " individuals are update from [" + string(sex_file) + "]." << endl;
    }
    
    void adj_grm(eInfo* einfo, double adj_grm_fac) {
        cout << "Adjusting the ERM for sampling errors ..." << endl;
        int i = 0, j = 0, n = einfo->_eii_include.size();
        double off_mean = 0.0, diag_mean = 0.0, off_var = 0.0, diag_var = 0.0, d_buf = 0.0;
        for (i = 0; i < n; i++) {
            diag_mean += einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[i]);
            for (j = 0; j < i; j++) off_mean += einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[j]);
        }
        diag_mean /= n;
        off_mean /= 0.5 * n * (n - 1.0);
        for (i = 0; i < n; i++) {
            d_buf = einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[i]) - diag_mean;
            diag_var += d_buf*d_buf;
            for (j = 0; j < i; j++) {
                d_buf = einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[j]) - off_mean;
                off_var += d_buf*d_buf;
            }
        }
        diag_var /= n - 1.0;
        off_var /= 0.5 * n * (n - 1.0) - 1.0;
        for (i = 0; i < einfo->_eii_include.size(); i++) {
            d_buf = 1.0 - (adj_grm_fac + 1.0 / einfo->_grm_N(einfo->_eii_include[i], einfo->_eii_include[i])) / diag_var;
            if (einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[i]) > 0) einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[i]) = 1.0 + d_buf * (einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[i]) - 1.0);
            for (j = 0; j < i; j++) {
                if (einfo->_grm_N(einfo->_eii_include[i], einfo->_eii_include[j]) > 0) einfo->_grm(einfo->_eii_include[i], einfo->_eii_include[j]) *= 1.0 - (adj_grm_fac + 1.0 / einfo->_grm_N(einfo->_eii_include[i], einfo->_eii_include[j])) / off_var;
            }
        }
    }
    
    
    void dc(eInfo* einfo, int dosage_compen) {
        cout << "Parameterizing the ERM under the assumption of ";
        if (dosage_compen == 1) cout << "full dosage compensation ..." << endl;
        else if (dosage_compen == 0) cout << "no dosage compensation ..." << endl;
        
        int i = 0, j = 0, i_buf = 0;
        double c1 = 1.0, c2 = 1.0;
        if (dosage_compen == 1) {
            c1 = 2.0;
            c2 = sqrt(2.0);
        }// full dosage compensation
        else if (dosage_compen == 0) {
            c1 = 0.5;
            c2 = sqrt(0.5);
        } // on dosage compensation
        for (i = 0; i < einfo->_eii_include.size(); i++) {
            for (j = 0; j <= i; j++) {
                i_buf = einfo->_eii_sex[einfo->_eii_include[i]] * einfo->_eii_sex[einfo->_eii_include[j]];
                if (i_buf == 1) einfo->_grm(i, j) *= c1;
                else if (i_buf == 2) einfo->_grm(i, j) *= c2;
            }
        }
    }
    
    void cal_var_mean(eInfo* einfo, bool mean_flag, bool var_flag)
    {
        long n = einfo->_eii_include.size(), m = einfo->_epi_num;
        uint64_t ttl_n = einfo->_eii_num;
        if(mean_flag) einfo->_mu.resize(m);
        if(var_flag) einfo->_var.resize(m);
        for(int j=0; j<m; j++){
            double mu=0.0;
            double nonmiss=0.0;
            for(int i=0; i<n; i++){
                double val=einfo->_val[j*ttl_n+einfo->_eii_include[i]];
                if(val<1e9){
                    mu+=val;
                    nonmiss+=1.0;
                }
            }
            mu/=nonmiss;
            if(mean_flag)  einfo->_mu[j]=mu;
            if(var_flag)
            {
                double var =0.0;
                for(int i=0; i<n; i++){
                    double val=einfo->_val[einfo->_epi_include[j]*ttl_n+einfo->_eii_include[i]];
                    if(val<1e9) var += (val-mu)*(val-mu);
                }
                var/=(nonmiss - 1.0);
                einfo->_var[j]=var;
            }
        }
    }

    void free_indilist(vector<indiinfolst> &a){
        for(int i=0;i<a.size();i++)
        {
            if(a[i].fid) free2(&a[i].fid);
            if(a[i].iid) free2(&a[i].iid);
            if(a[i].fa_id) free2(&a[i].fa_id);
            if(a[i].mo_id) free2(&a[i].mo_id);
        }
    }
    void free_probelist(vector<probeinfolst> &a){
        for(int i=0;i<a.size();i++)
        {
            if(a[i].probeId) free2(&a[i].probeId);
            if(a[i].genename) free2(&a[i].genename);
        }
    }
    void free_assoclist(vector<ASSOCRLT> &a){
        for(int i=0;i<a.size();i++)
        {
            if(a[i].PROBE) free2(&a[i].PROBE);
            if(a[i].GENE) free2(&a[i].GENE);
        }
    }


}


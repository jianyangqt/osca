//
//  l2_besd.cpp
//  osc
//
//  Created by Futao Zhang on 10/11/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#include "l2_besd.hpp"

namespace SMR {
    
    void read_smr_epifile(eqtlInfo* eqtlinfo, char* epiFileName)
    {
        FILE* epifile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum=6;
        if(fopen_checked(&epifile, epiFileName,"r")) TERMINATE();
        LOGPRINTF("Reading eQTL probe information from %s ...\n", epiFileName);
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
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
                LOGPRINTF("WARNING: Line %u has more than %d items. The first %d columns would be used. \n", line_idx,colnum,colnum);
            }
            if(strlist.size()<colnum)
            {
                LOGPRINTF("ERROR: Line %u has less than %d items.\n", line_idx,colnum);
                TERMINATE();
            }
            eqtlinfo->_probe_name_map.insert(pair<string,int>(strlist[1],line_idx));
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
                LOGPRINTF("WARNING: unrecongized chromosome found. This chromosome is set to 0:\n");
                LOGPRINTF("%s\n",Tbuf);
                //TERMINATE();
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            } else if ( atoi(strlist[0].c_str())>24 || atoi(strlist[0].c_str())<0) {
                LOGPRINTF("WARNING: abmormal chromosome found:\n");
                LOGPRINTF("%s\n",Tbuf);
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            } else eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            
            if(strlist[1]=="NA" || strlist[1]=="na") {
                LOGPRINTF("ERROR: NA probe ID found:\n");
                LOGPRINTF("%s\n",Tbuf);
                TERMINATE();
            }
            eqtlinfo->_epi_prbID.push_back(strlist[1]);
            eqtlinfo->_epi_gd.push_back(atoi(strlist[2].c_str()));
            if(strlist[3]=="NA" || strlist[3]=="na") eqtlinfo->_epi_bp.push_back(-9);
            else eqtlinfo->_epi_bp.push_back(atoi(strlist[3].c_str()));
            if(strlist[4]=="NA" || strlist[4]=="na") {
                if(!genewarning) {
                    LOGPRINTF("WARNING: at least one gene id is missing.\n");
                    genewarning=true;
                }
            }
            eqtlinfo->_epi_gene.push_back(strlist[4].c_str());
            if(strlist[5]=="NA") {
                eqtlinfo->_epi_orien.push_back('*');
                if(!orienwarning) {
                    LOGPRINTF("WARNING: At least one gene strand is missing.\n");
                    orienwarning=true;
                }
            } else eqtlinfo->_epi_orien.push_back(strlist[5][0]);
            eqtlinfo->_include.push_back(line_idx);
            line_idx++;
        }
        eqtlinfo->_probNum =line_idx;
        fclose(epifile);
        LOGPRINTF("%llu probes to be included from  %s .\n", eqtlinfo->_probNum, epiFileName);
    }
    
    void read_smr_esifile(eqtlInfo* eqtlinfo, char* esiFileName)
    {
        FILE* esifile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum=7;
        if(fopen_checked(&esifile, esiFileName,"r")) TERMINATE();
        LOGPRINTF("Reading eQTL SNP information from %s ...\n", esiFileName);
        eqtlinfo->_esi_chr.clear();
        eqtlinfo->_esi_rs.clear();
        eqtlinfo->_esi_gd.clear();
        eqtlinfo->_esi_bp.clear();
        eqtlinfo->_esi_allele1.clear();
        eqtlinfo->_esi_allele2.clear();
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        eqtlinfo->_esi_freq.clear();

        bool chrwarning=false;
        bool feqwarning=false, allele1warning=false, allele2warning=false;
        bool orienwarning=false;
        while(fgets(Tbuf, MAX_LINE_SIZE, esifile))
        {
            split_str(Tbuf,strlist,0);
            if(Tbuf[0]=='\0') {
                LOGPRINTF("ERROR: Line %u is blank.\n", line_idx);
                TERMINATE();
            }
            if(strlist.size()<colnum-1)
            {
                LOGPRINTF("ERROR: Line %u has less than %d items.\n", line_idx,colnum-1);
                TERMINATE();
            } else if(strlist.size()==colnum-1) {
                if(!feqwarning) {
                    LOGPRINTF("WARNING: Maybe this is an old .esi file which doesn't contain frequency information. \n");
                    feqwarning=true;
                }
            } else if(strlist.size()>colnum) {
                LOGPRINTF("WARNING: Line %u has more than %d items. The first %d columns would be used. \n", line_idx,colnum,colnum);
            }
            
            eqtlinfo->_snp_name_map.insert(pair<string,int>(strlist[1],line_idx));
            if(eqtlinfo->_snp_name_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate SNP : %s.\n", strlist[1].c_str());
                TERMINATE();
            }
            if(strlist[0]=="X" || strlist[0]=="x") eqtlinfo->_esi_chr.push_back(23);
            else if(strlist[0]=="Y" || strlist[0]=="y") eqtlinfo->_esi_chr.push_back(24);
            else if(strlist[0]=="NA" || strlist[0]=="na"){
                eqtlinfo->_esi_chr.push_back(-9);
                if(!chrwarning) {
                    LOGPRINTF("WARNING: At least one SNP chr is missing.\n");
                    chrwarning=true;
                }
            } else if (atoi(strlist[0].c_str())==0 ) {
                LOGPRINTF("WARNING: unrecongized chromosome found. This chromosome is set to 0:\n");
                LOGPRINTF("%s\n",Tbuf);
                eqtlinfo->_esi_chr.push_back(atoi(strlist[0].c_str()));
            } else if ( atoi(strlist[0].c_str())>24 || atoi(strlist[0].c_str())<0) {
                LOGPRINTF("WARNING: abmormal chromosome found:\n");
                LOGPRINTF("%s\n",Tbuf);
                eqtlinfo->_esi_chr.push_back(atoi(strlist[0].c_str()));
            } else eqtlinfo->_esi_chr.push_back(atoi(strlist[0].c_str()));
            
            if(strlist[1]=="NA" || strlist[1]=="na") {
                LOGPRINTF("ERROR: NA SNP ID found:\n");
                LOGPRINTF("%s\n",Tbuf);
                TERMINATE();
            }
            eqtlinfo->_esi_rs.push_back(strlist[1]);
            eqtlinfo->_esi_gd.push_back(atoi(strlist[2].c_str()));
            if(strlist[3]=="NA" || strlist[3]=="na") eqtlinfo->_esi_bp.push_back(-9);
            else eqtlinfo->_esi_bp.push_back(atoi(strlist[3].c_str()));
            if(strlist[4]=="NA" || strlist[4]=="na") {
                if(!allele1warning) {
                    LOGPRINTF("WARNING: At least one reference allele is missing.\n");
                    allele1warning=true;
                }
            }
            to_upper(strlist[4]);
            eqtlinfo->_esi_allele1.push_back(strlist[4].c_str());
            if(strlist[5]=="NA" || strlist[5]=="na") {
                if(!allele2warning) {
                    LOGPRINTF("WARNING: At least one alternative allele is missing.\n");
                    allele2warning=true;
                }
            }
            to_upper(strlist[5]);
            eqtlinfo->_esi_allele2.push_back(strlist[5].c_str());
            if(strlist.size()==colnum)
            {
                if(strlist[6]=="NA" || strlist[6]=="na"){
                    if(!orienwarning){
                        LOGPRINTF("WARNING: frequency is \"NA\" in one or more rows.\n");
                        orienwarning=true;
                    }
                    eqtlinfo->_esi_freq.push_back(-9);
                } else {
                    eqtlinfo->_esi_freq.push_back(atof(strlist[6].c_str()));
                }
            } else {
                eqtlinfo->_esi_freq.push_back(-9);
            }
            eqtlinfo->_esi_include.push_back(line_idx);
            line_idx++;
        }
        eqtlinfo->_snpNum =line_idx;
        fclose(esifile);
        LOGPRINTF("%llu SNPs to be included from  %s .\n", eqtlinfo->_snpNum, esiFileName);
    }
    void read_indicators(FILE** besd, eqtlInfo* eqtlinfo)
    {
        int length=(RESERVEDUNITS-1)*sizeof(int);
        char* indicators=new char[length];
        fread(indicators, sizeof(int),(RESERVEDUNITS-1), *besd);
        int* tmp=(int *)indicators;
        int ss=*tmp++;
        if(ss!=-9)
        {
            eqtlinfo->_sampleNum=ss;
            LOGPRINTF("The sample size is %d.\n",ss);
        }
        if(*tmp++!=eqtlinfo->_snpNum)
        {
            LOGPRINTF("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file.\n");
            TERMINATE();
        }
        if(*tmp++!=eqtlinfo->_probNum)
        {
            LOGPRINTF("ERROR: The probes in your .epi file are not in consistency with the one in .besd file.\n");
            TERMINATE();
        }
        delete[] indicators;

    }
    void read_smr_besdfile(eqtlInfo* eqtlinfo, char* besdFileName)
    {
        if (eqtlinfo->_include.size() == 0) {
            LOGPRINTF("Error: No probe is retained for analysis.\n");
            TERMINATE();
        }
        if (eqtlinfo->_esi_include.size() == 0) {
            LOGPRINTF("Error: No SNP is retained for analysis.\n");
            TERMINATE();
        }
        
        eqtlinfo->_cols.clear();
        eqtlinfo->_rowid.clear();
        eqtlinfo->_val.clear();
        eqtlinfo->_valNum = 0;
        eqtlinfo->_bxz.clear();
        eqtlinfo->_sexz.clear();
        eqtlinfo->_sampleNum=-9;
        
        FILE* besd=NULL;
        if(fopen_checked(&besd,besdFileName,"rb")) {
            TERMINATE();
        }
        LOGPRINTF("Reading eQTL summary data from %s. \n",besdFileName);
        uint32_t indicator;
        if(fread(&indicator, sizeof(uint32_t),1, besd)!=1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
            TERMINATE();
        }
        if(indicator == SMR_SPARSE_3F || indicator == SMR_SPARSE_3)
        {
            char* buffer;
            uint64_t colNum=(eqtlinfo->_probNum<<1)+1;
            uint64_t valNum;
            
            uint64_t cur_pos = ftell( besd );
            fseek( besd, 0L, SEEK_END );
            uint64_t size_file = ftell( besd );
            fseek( besd, cur_pos, SEEK_SET );
            //if( indicator == SMR_SPARSE_3) read_indicators(&besd, eqtlinfo);
            if( indicator == SMR_SPARSE_3)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), besd);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    eqtlinfo->_sampleNum=ss;
                    LOGPRINTF("The sample size is %d.\n",ss);
                }
                if(*tmp++!=eqtlinfo->_snpNum)
                {
                    LOGPRINTF("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                    TERMINATE();
                }
                if(*tmp++!=eqtlinfo->_probNum)
                {
                    LOGPRINTF("ERROR: The probes in your .epi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                    TERMINATE();
                }
                delete[] indicators;
            }
            if(fread(&valNum, sizeof(uint64_t),1, besd)!=1)
            {
                LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                TERMINATE();
            }
            int descriptive=1;
            if(indicator == SMR_SPARSE_3) descriptive=RESERVEDUNITS;
            if( size_file - (descriptive*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {
                LOGPRINTF("ERROR: File %s is broken!\n", besdFileName);
                TERMINATE();
            }
            
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {
                
                buffer = (char*) malloc (sizeof(uint64_t)*(colNum));
                if (buffer == NULL) {
                    LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                    TERMINATE();
                }
                if(fread(buffer, sizeof(uint64_t),colNum, besd)!=colNum)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                    TERMINATE();
                }

                uint64_t* ptr=(uint64_t *)buffer;
                if(ptr[colNum-1] != valNum)
                {
                    LOGPRINTF("ERROR: File %s is broken!\n", besdFileName);
                    TERMINATE();
                }
                eqtlinfo->_cols.resize(eqtlinfo->_include.size()+1);
                eqtlinfo->_cols[0]=*ptr;
                
                map<int, int > _incld_id_map;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) {
                        LOGPRINTF("Error: Duplicated SNP IDs found: %s.\n",eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]].c_str());
                        TERMINATE();
                    }
                    size = _incld_id_map.size();
                }
                
                uint64_t rowSTART=descriptive*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                uint64_t valSTART=descriptive*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    uint32_t pid=eqtlinfo->_include[i];
                    uint64_t pos=*(ptr+(pid<<1)); //BETA START
                    uint64_t pos1=*(ptr+(pid<<1)+1); //SE START
                    uint64_t num=pos1-pos;
                    uint64_t real_num=0;
                    if(num==0) {
                        eqtlinfo->_cols[i+1]=eqtlinfo->_cols[i];
                        //LOGPRINTF("WARNING: Probe %s with no eQTL found.\n",eqtlinfo->_epi_prbID[pid].c_str());
                        continue;
                        
                    }
                    char* row_char_ptr;
                    row_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(uint32_t));
                    if (row_char_ptr == NULL) {
                        LOGPRINTF("Memory error.\n");
                        TERMINATE();
                    }
                    char* val_char_ptr;
                    val_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(float));
                    if (val_char_ptr == NULL) {
                        LOGPRINTF("Memory error.\n");
                        TERMINATE();
                    }
                    fseek( besd, rowSTART+pos*sizeof(uint32_t), SEEK_SET );
                    if(fread(row_char_ptr, sizeof(uint32_t),2*num, besd)!=2*num)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                        TERMINATE();
                    }
                    uint32_t* row_ptr=(uint32_t *)row_char_ptr;
                    fseek( besd, valSTART+pos*sizeof(float), SEEK_SET );
                    if(fread(val_char_ptr, sizeof(float),2*num, besd)!=2*num)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                        TERMINATE();
                    }
                    float* val_ptr=(float*)val_char_ptr;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=*(row_ptr+j);
                        
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;
                            
                            if(j<num) eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+j));
                            real_num++;
                        }
                    }
                    eqtlinfo->_cols[i+1]=real_num+eqtlinfo->_cols[i];
                    free(row_char_ptr);
                    free(val_char_ptr);
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
                
                LOGPRINTF("eQTL summary data of %ld Probes to be included from %s.\n",eqtlinfo->_include.size(), besdFileName);
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                uint64_t size2read=sizeof(char)*(size_file)-descriptive*sizeof(int) - sizeof(uint64_t);
                buffer = (char*) malloc (size2read);
                if (buffer == NULL) {
                    LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                    TERMINATE();
                }
                if(fread(buffer, sizeof(char),size2read, besd)!=size2read)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                    TERMINATE();
                }
                
                uint64_t* ptr=(uint64_t *)buffer;
                
                eqtlinfo->_cols.resize(eqtlinfo->_probNum+1);
                eqtlinfo->_rowid.resize(valNum>>1);
                eqtlinfo->_val.resize(valNum);
                
                for(int i=0;i<colNum;i++) {
                    uint64_t tmp= *ptr++;
                    if((i & 1) == 0) eqtlinfo->_cols[i>>1]= tmp;
                }
                uint32_t* ptr4B=(uint32_t *)ptr;
                long count=0;
                for(int i=0;i<eqtlinfo->_probNum;i++) {
                    long num=eqtlinfo->_cols[i+1]-eqtlinfo->_cols[i];
                    long num2=num>>1;
                    for(int j=0;j<num;j++) {
                        uint32_t tmp=*ptr4B++;
                        if(j<num2) eqtlinfo->_rowid[count++]=tmp;
                    }
                }
                if(count != valNum>>1) {
                    LOGPRINTF("ERROR: bugs in reading sparse 3F. help to report pls.\n");
                    TERMINATE();
                }
                float* val_ptr=(float*)ptr4B;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*val_ptr++;
                eqtlinfo->_valNum = valNum;
                LOGPRINTF("eQTL summary data of %llu Probes to be included from %s.\n",eqtlinfo->_probNum,besdFileName);
            }
            // terminate
            free (buffer);
        }
        else if(indicator == SMR_DENSE_1 || indicator == SMR_DENSE_3 )
        {
            if(indicator==SMR_DENSE_3)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), besd);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    eqtlinfo->_sampleNum=ss;
                    LOGPRINTF("The sample size is %d.\n",ss);
                }
                if(*tmp++!=eqtlinfo->_snpNum)
                {
                    LOGPRINTF("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                    TERMINATE();
                }
                if(*tmp++!=eqtlinfo->_probNum)
                {
                    LOGPRINTF("ERROR: The probes in your .epi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                    TERMINATE();
                }
                delete[] indicators;
            }
            int infoLen=sizeof(uint32_t);
            if(indicator==SMR_DENSE_3) infoLen=RESERVEDUNITS*sizeof(int);
            
            uint64_t cur_pos = ftell( besd );
            fseek( besd, 0L, SEEK_END );
            uint64_t size_file = ftell( besd );
            fseek( besd, cur_pos, SEEK_SET );
            if(size_file!= eqtlinfo->_probNum*eqtlinfo->_snpNum*sizeof(float)*2 + infoLen) {
                LOGPRINTF("ERROR: File %s is broken!\n", besdFileName);
                TERMINATE();
            }
            uint64_t mem2use=eqtlinfo->_include.size()*eqtlinfo->_esi_include.size()*2*sizeof(float);
            if(mem2use>0x200000000){
                LOGPRINTF("WARNING: %llu GB should be allocated for your besd file.\n",mem2use>>30);
            }

            eqtlinfo->_bxz.resize(eqtlinfo->_include.size());
            eqtlinfo->_sexz.resize(eqtlinfo->_include.size());
            for(unsigned int i=0;i<eqtlinfo->_include.size();i++)
            {
                eqtlinfo->_bxz[i].resize(eqtlinfo->_esi_include.size());
                eqtlinfo->_sexz[i].resize(eqtlinfo->_esi_include.size());
            }
            char* buffer;
            buffer = (char*) malloc (sizeof(char)*eqtlinfo->_snpNum<<3);
            if (buffer == NULL) {
                LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                TERMINATE();
            }
            float* ft;
            float* se_ptr;
            if (eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)  //means with the parameter --extract-probe. This also can read all the probes, but currently I don't think it is good for too many I/Os.
            {
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    unsigned long pid=eqtlinfo->_include[i];
                    uint64_t readpos=((pid*eqtlinfo->_snpNum)<<3)+infoLen;
                    
                    fseek( besd, readpos, SEEK_SET );
                    if(fread(buffer, sizeof(char),eqtlinfo->_snpNum<<3, besd)!=eqtlinfo->_snpNum<<3)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                        TERMINATE();
                    }
                    ft=(float *)buffer;
                    for (int j = 0; j<eqtlinfo->_esi_include.size(); j++) eqtlinfo->_bxz[i][j] = *(ft + eqtlinfo->_esi_include[j]);
                    se_ptr = ft + eqtlinfo->_snpNum;
                    for (int j = 0; j<eqtlinfo->_esi_include.size(); j++) eqtlinfo->_sexz[i][j] = *(se_ptr + eqtlinfo->_esi_include[j]);
                }
                LOGPRINTF("eQTL summary-level statistics of %ld Probes and %ld SNPs to be included from %s.\n" , eqtlinfo->_include.size(), eqtlinfo->_esi_include.size(),besdFileName);
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                char* buff;
                uint64_t buffszie=0x40000000;
                buff = (char*) malloc (sizeof(char)*buffszie);
                if (buff == NULL) {
                    LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                    TERMINATE();
                }
                memset(buff,0,sizeof(char)*buffszie);
                
                uint64_t perbeta=(eqtlinfo->_snpNum<<2);
                uint64_t probonce=sizeof(char)*buffszie/perbeta;  //should be even number
                probonce>>=1;
                probonce<<=1;
                uint64_t readsize=perbeta*probonce;
                uint64_t probcount=0;
                while(!feof(besd))
                {
                    uint64_t Bread=fread(buff, sizeof(char),readsize, besd);
                    
                    char* rptr=buff;
                    while(Bread)
                    {
                        memcpy(&eqtlinfo->_bxz[probcount][0],rptr,perbeta);
                        rptr+=perbeta;
                        memcpy(&eqtlinfo->_sexz[probcount++][0],rptr,perbeta);
                        rptr+=perbeta;
                        Bread-=(perbeta<<1);
                    }
                    printf("Redinging... %3.0f%%\r", 100.0*probcount/eqtlinfo->_probNum);
                    fflush(stdout);
                }
                LOGPRINTF("\neQTL summary data of %llu Probes and %llu SNPs to be included from %s.\n",eqtlinfo->_probNum,eqtlinfo->_snpNum,besdFileName );
                free(buff);
                
            }
            free(buffer);

        }
        else if(indicator == OSCA_DENSE_1)
        {
          
            int length=(RESERVEDUNITS-1)*sizeof(int);
            char* indicators=new char[length];
            fread(indicators, sizeof(int),(RESERVEDUNITS-1), besd);
            int* tmp=(int *)indicators;
            int ss=*tmp++;
            if(ss!=-9)
            {
                eqtlinfo->_sampleNum=ss;
                LOGPRINTF("The sample size is %d.\n",ss);
            }
            if(*tmp++!=eqtlinfo->_snpNum)
            {
                LOGPRINTF("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                TERMINATE();
            }
            if(*tmp++!=eqtlinfo->_probNum)
            {
                LOGPRINTF("ERROR: The probes in your .epi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                TERMINATE();
            }
            delete[] indicators;
            int infoLen=RESERVEDUNITS*sizeof(int);
            
            uint64_t cur_pos = ftell( besd );
            fseek( besd, 0L, SEEK_END );
            uint64_t size_file = ftell( besd );
            fseek( besd, cur_pos, SEEK_SET );
            if(size_file!= eqtlinfo->_probNum*eqtlinfo->_snpNum*sizeof(float)*2 + infoLen) {
                LOGPRINTF("ERROR: File %s is broken!\n", besdFileName);
                TERMINATE();
            }
            uint64_t mem2use=eqtlinfo->_include.size()*eqtlinfo->_esi_include.size()*2*sizeof(float);
            if(mem2use>0x200000000){
                LOGPRINTF("WARNING: %llu GB should be allocated for your besd file.\n",mem2use>>30);
            }
            
            eqtlinfo->_bxz.resize(eqtlinfo->_include.size());
            eqtlinfo->_sexz.resize(eqtlinfo->_include.size());
            for(unsigned int i=0;i<eqtlinfo->_include.size();i++)
            {
                eqtlinfo->_bxz[i].resize(eqtlinfo->_esi_include.size());
                eqtlinfo->_sexz[i].resize(eqtlinfo->_esi_include.size());
            }
            char* buffer;
            buffer = (char*) malloc (sizeof(char)*eqtlinfo->_probNum<<3);
            if (buffer == NULL) {
                LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                TERMINATE();
            }
            
            float* ft;
            for(int j=0;j<eqtlinfo->_esi_include.size();j++)
            {
                uint64_t sid=eqtlinfo->_esi_include[j];
                uint64_t readpos=((sid*eqtlinfo->_probNum)<<3)+infoLen;
                fseek( besd, readpos, SEEK_SET );
                if(fread(buffer, sizeof(char),eqtlinfo->_probNum<<3, besd)!=eqtlinfo->_probNum<<3)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                    TERMINATE();
                }
                ft=(float *)buffer;
                for (int i = 0; i<eqtlinfo->_include.size(); i++)
                {
                    int pos=eqtlinfo->_include[i];
                    eqtlinfo->_bxz[i][j] = *(ft + 2*pos);
                    eqtlinfo->_sexz[i][j] = *(ft + 2*pos +1);
                }
            }
            LOGPRINTF("eQTL summary-level statistics of %ld Probes and %ld SNPs to be included from %s.\n",eqtlinfo->_include.size(),eqtlinfo->_esi_include.size(),besdFileName);
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
            if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            
            free(buffer);
        }
        else if(indicator == OSCA_SPARSE_1)
        {
            int length=(RESERVEDUNITS-1)*sizeof(int);
            char* indicators=new char[length];
            fread(indicators, sizeof(int),(RESERVEDUNITS-1), besd);
            int* tmp=(int *)indicators;
            int ss=*tmp++;
            if(ss!=-9)
            {
                eqtlinfo->_sampleNum=ss;
                LOGPRINTF("The sample size is %d.\n",ss);
            }
            if(*tmp++!=eqtlinfo->_snpNum)
            {
                LOGPRINTF("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                TERMINATE();
            }
            if(*tmp++!=eqtlinfo->_probNum)
            {
                LOGPRINTF("ERROR: The probes in your .epi file are not in consistency with the one in .besd file %s.\n", besdFileName);
                TERMINATE();
            }
            delete[] indicators;
            int infoLen=RESERVEDUNITS*sizeof(int);
            
            char* buffer;
            uint64_t colNum=eqtlinfo->_probNum+1;
            uint64_t valNum;
            
            uint64_t cur_pos = ftell( besd );
            fseek( besd, 0L, SEEK_END );
            uint64_t size_file = ftell( besd );
            fseek( besd, cur_pos, SEEK_SET );
            
            if(fread(&valNum, sizeof(uint64_t),1, besd)!=1)
            {
                LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                TERMINATE();
            }
            
            if( size_file - (infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t)/2 + valNum*sizeof(float)) != 0) {
                LOGPRINTF("ERROR: File %s is broken!\n", besdFileName);
                TERMINATE();
            }
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {
                
                buffer = (char*) malloc (sizeof(uint64_t)*(colNum));
                if (buffer == NULL) {
                    LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                    TERMINATE();
                }
                if(fread(buffer, sizeof(uint64_t),colNum, besd)!=colNum)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                    TERMINATE();
                }
                
                uint64_t* ptr=(uint64_t *)buffer;
                
                eqtlinfo->_cols.resize(eqtlinfo->_include.size()+1);
                eqtlinfo->_cols[0]=*ptr;
                
                map<int, int > _incld_id_map;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) {
                        LOGPRINTF("Error: Duplicated SNP IDs found: %s.\n",eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]].c_str());
                        TERMINATE();
                    }
                    size = _incld_id_map.size();
                }
                
                uint64_t rowSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                uint64_t valSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t)/2;
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    uint32_t pid=eqtlinfo->_include[i];
                    uint64_t pos=*(ptr+pid); //BETA START
                    uint64_t pos1=*(ptr+pid+1); //next BETA START
                    uint64_t num=(pos1-pos)/2; //num must be even
                    uint64_t colpos=pos>>1;
                    uint64_t real_num=0;
                    if(num==0) {
                        eqtlinfo->_cols[i+1]=eqtlinfo->_cols[i];
                        LOGPRINTF("WARNING: Probe %s with no eQTL found.\n",eqtlinfo->_epi_prbID[pid].c_str());
                        continue;
                        
                    }
                    char* row_char_ptr;
                    row_char_ptr = (char*) malloc (sizeof(char)*num*sizeof(uint32_t));
                    if (row_char_ptr == NULL) {
                        LOGPRINTF("Memory error.\n");
                        TERMINATE();
                    }
                    char* val_char_ptr;
                    val_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(float));
                    if (val_char_ptr == NULL) {
                        LOGPRINTF("Memory error.\n");
                        TERMINATE();
                    }
                    
                    fseek( besd, rowSTART+colpos*sizeof(uint32_t), SEEK_SET );
                    if(fread(row_char_ptr, sizeof(uint32_t),num, besd)!=num)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                        TERMINATE();
                    }
                    
                    fseek( besd, valSTART+pos*sizeof(float), SEEK_SET );
                    if(fread(val_char_ptr, sizeof(float),2*num, besd)!=2*num)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                        TERMINATE();
                    }
                    uint32_t* row_ptr=(uint32_t *)row_char_ptr;
                    float* val_ptr=(float *)val_char_ptr;
                    vector<float> setmp;
                    for(int j=0;j<num;j++)
                    {
                        uint32_t rid=*(row_ptr+j);
                        
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;
                            
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+j));
                            setmp.push_back(*(val_ptr+num+j));
                            real_num++;
                        }
                    }
                    for(int j=0;j<setmp.size();j++) eqtlinfo->_val.push_back(setmp[j]);
                    eqtlinfo->_cols[i+1]=(real_num<<1)+eqtlinfo->_cols[i];
                    free(row_ptr);
                    free(val_ptr);
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
                
                LOGPRINTF("eQTL summary data of %ld Probes to be included from %s.\n",eqtlinfo->_include.size(), besdFileName);
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                uint64_t size2read=sizeof(char)*(size_file)- infoLen - sizeof(uint64_t);
                buffer = (char*) malloc (size2read);
                if (buffer == NULL) {
                    LOGPRINTF("ERROR: memory allocation failed to read %s.\n", besdFileName);
                    TERMINATE();
                }
                if(fread(buffer, sizeof(char),size2read, besd)!=size2read)
                {
                    LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
                    TERMINATE();
                }
                
                uint64_t* ptr=(uint64_t *)buffer;
                
                eqtlinfo->_cols.resize(colNum);
                eqtlinfo->_rowid.resize(valNum>>1);
                eqtlinfo->_val.resize(valNum);
                
                for(int i=0;i<colNum;i++) eqtlinfo->_cols[i]=*ptr++;
                uint32_t* ptr4B=(uint32_t *)ptr;
                for(int i=0;i<valNum/2;i++) eqtlinfo->_rowid[i]=*ptr4B++;
                float* val_ptr=(float*)ptr4B;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*val_ptr++;
                eqtlinfo->_valNum = valNum;
                LOGPRINTF("eQTL summary data of %llu Probes to be included from %s.\n",eqtlinfo->_probNum,besdFileName);
            }
            free (buffer);
        }
        else if(indicator == 0x40000000 || indicator == 0x3f800000)
        {
            LOGPRINTF("The file is an old BESD format, please use SMR to deal with it.\n");
            TERMINATE();
        } else {
            LOGPRINTF("The file is broken.\n");
            TERMINATE();
        }
        fclose(besd);
    }
    void get_BesdHeaders(char* besdFileName, vector<int> &headers)//
    {
        headers.resize(RESERVEDUNITS);
        FILE* besd=NULL;
        if(fopen_checked(&besd,besdFileName,"rb")) {
            TERMINATE();
        }
        LOGPRINTF("Reading eQTL summary data from %s. \n",besdFileName);
        if(fread(&headers[0], sizeof(int),RESERVEDUNITS, besd)<1)
        {
            LOGPRINTF("ERROR: File %s read failed!\n", besdFileName);
            TERMINATE();
        }
        fclose(besd);
    }
    void update_epi(eqtlInfo* eqtlinfo)
    {
        eqtlinfo->_probNum = eqtlinfo->_include.size();
        
        vector<int> chr_buf, gd_buf,bp_buf, start, end;
        vector<string> prbID_buf, gene_buf;
        vector<char> orien_buf;
        for (int i = 0; i < eqtlinfo->_probNum; i++)
        {
            chr_buf.push_back(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]);
            gd_buf.push_back(eqtlinfo->_epi_gd[eqtlinfo->_include[i]]);
            bp_buf.push_back(eqtlinfo->_epi_bp[eqtlinfo->_include[i]]);
            prbID_buf.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
            gene_buf.push_back(eqtlinfo->_epi_gene[eqtlinfo->_include[i]]);
            orien_buf.push_back(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]);
            if(eqtlinfo->_epi_start.size()>0) start.push_back(eqtlinfo->_epi_start[eqtlinfo->_include[i]]);
            if(eqtlinfo->_epi_end.size()>0) start.push_back(eqtlinfo->_epi_end[eqtlinfo->_include[i]]);
        }
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
        eqtlinfo->_epi_start.clear();
        eqtlinfo->_epi_end.clear();
        eqtlinfo->_epi_chr.swap(chr_buf);
        eqtlinfo->_epi_gd.swap(gd_buf);
        eqtlinfo->_epi_bp.swap(bp_buf);
        eqtlinfo->_epi_prbID.swap(prbID_buf);
        eqtlinfo->_epi_gene.swap(gene_buf);
        eqtlinfo->_epi_orien.swap(orien_buf);
        eqtlinfo->_epi_start.swap(start);
        eqtlinfo->_epi_end.swap(end);
        
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        for (int i = 0; i < eqtlinfo->_probNum; i++)
        {
            eqtlinfo->_include.push_back(i);
            eqtlinfo->_probe_name_map.insert(pair<string,int>(eqtlinfo->_epi_prbID[i],i));
        }
    }
    
    void update_esi(eqtlInfo* eqtlinfo)
    {
        eqtlinfo->_snpNum = eqtlinfo->_esi_include.size();
        
        vector<int> chr_buf, gd_buf, bp_buf;
        vector<string> rs_buf;
        vector<string> allele1_buf, allele2_buf;
        vector<float> freq_buf;
        for (int i = 0; i < eqtlinfo->_snpNum; i++)
        {
            chr_buf.push_back(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]);
            gd_buf.push_back(eqtlinfo->_esi_gd[eqtlinfo->_esi_include[i]]);
            bp_buf.push_back(eqtlinfo->_esi_bp[eqtlinfo->_esi_include[i]]);
            rs_buf.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
            allele1_buf.push_back(eqtlinfo->_esi_allele1[eqtlinfo->_esi_include[i]]);
            allele2_buf.push_back(eqtlinfo->_esi_allele2[eqtlinfo->_esi_include[i]]);
            freq_buf.push_back(eqtlinfo->_esi_freq[eqtlinfo->_esi_include[i]]);
        }
        eqtlinfo->_esi_chr.clear();
        eqtlinfo->_esi_gd.clear();
        eqtlinfo->_esi_bp.clear();
        eqtlinfo->_esi_rs.clear();
        eqtlinfo->_esi_allele1.clear();
        eqtlinfo->_esi_allele2.clear();
        eqtlinfo->_esi_freq.clear();
        eqtlinfo->_esi_chr.swap(chr_buf);
        eqtlinfo->_esi_gd.swap(gd_buf);
        eqtlinfo->_esi_bp.swap(bp_buf);
        eqtlinfo->_esi_rs.swap(rs_buf);
        eqtlinfo->_esi_allele1.swap(allele1_buf);
        eqtlinfo->_esi_allele2.swap(allele2_buf);
        eqtlinfo->_esi_freq.swap(freq_buf);
        
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        for (int i = 0; i < eqtlinfo->_snpNum; i++)
        {
            eqtlinfo->_esi_include.push_back(i);
            eqtlinfo->_snp_name_map.insert(pair<string, int>(eqtlinfo->_esi_rs[i], i));
        }
    }


    void write_smr_epi(char* outFileName, eInfo* einfo, bool bpAsstart)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".epi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        LOGPRINTF("Saving probe information ...\n");
        for(int i=0;i<einfo->_epi_include.size();i++)
        {
            string chrstr;
            if(einfo->_epi_chr[einfo->_epi_include[i]]==23) chrstr="X";
            else if(einfo->_epi_chr[einfo->_epi_include[i]]==24) chrstr="Y";
            else chrstr=atosm(einfo->_epi_chr[einfo->_epi_include[i]]);
            int curbp = einfo->_epi_bp[einfo->_epi_include[i]];
            if(bpAsstart) curbp = einfo->_epi_bp[einfo->_epi_include[i]] + (einfo->_epi_gd[einfo->_epi_include[i]] -einfo->_epi_bp[einfo->_epi_include[i]])/2;
            string str=chrstr+'\t'+einfo->_epi_prb[einfo->_epi_include[i]]+'\t'+atos(0)+'\t'+atosm(curbp)+'\t'+einfo->_epi_gene[einfo->_epi_include[i]]+'\t'+(einfo->_epi_orien[einfo->_epi_include[i]]=='*'?"NA":atos(einfo->_epi_orien[einfo->_epi_include[i]]))+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld probes have been saved in the file %s .\n", einfo->_epi_include.size(), epiName.c_str());
    }
    
    void write_smr_esi(char* outFileName, bInfo* binfo)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".esi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        LOGPRINTF("Saving SNP information ...\n");
        for(int i=0;i<binfo->_include.size();i++)
        {
            string chrstr;
            if(binfo->_chr[binfo->_include[i]]==23) chrstr="X";
            else if(binfo->_chr[binfo->_include[i]]==24) chrstr="Y";
            else chrstr=atosm(binfo->_chr[binfo->_include[i]]);
            string freqstr;
            if(binfo->_mu.empty()) freqstr="NA";
            else freqstr=atos(binfo->_mu[binfo->_include[i]] /2 );
            string str=chrstr+'\t'+binfo->_snp_name[binfo->_include[i]]+'\t'+atos(0)+'\t'+atosm(binfo->_bp[binfo->_include[i]])+'\t'+binfo->_allele1[binfo->_include[i]]+'\t'+binfo->_allele2[binfo->_include[i]]+'\t'+freqstr+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld SNPs have been saved in the file %s .\n", binfo->_include.size(), epiName.c_str());
    }
    void write_smr_esi(char* outFileName, eqtlInfo* eqtlinfo)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".esi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        LOGPRINTF("Saving SNP information ...\n");
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            string chrstr;
            if(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]==23) chrstr="X";
            else if(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]==24) chrstr="Y";
            else chrstr=atosm(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]);
            
            string str=chrstr+'\t'+eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]+'\t'+atos(0)+'\t'+atosm(eqtlinfo->_esi_bp[eqtlinfo->_esi_include[i]])+'\t'+eqtlinfo->_esi_allele1[eqtlinfo->_esi_include[i]]+'\t'+eqtlinfo->_esi_allele2[eqtlinfo->_esi_include[i]]+'\t'+atosm(eqtlinfo->_esi_freq[eqtlinfo->_esi_include[i]])+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld SNPs have been saved in the file %s .\n", eqtlinfo->_esi_include.size(), epiName.c_str());
    }
    void write_smr_epi(char* outFileName, eqtlInfo* eqtlinfo)
    {
        FILE* efile=NULL;
        string epiName=string(outFileName)+".epi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        LOGPRINTF("Saving probe information ...\n");
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            string chrstr;
            if(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]==23) chrstr="X";
            else if(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]==24) chrstr="Y";
            else chrstr=atosm(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]);
            
            string str=chrstr+'\t'+eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]+'\t'+atos(0)+'\t'+atosm(eqtlinfo->_epi_bp[eqtlinfo->_include[i]])+'\t'+eqtlinfo->_epi_gene[eqtlinfo->_include[i]]+'\t'+(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]=='*'?"NA":atos(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]))+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld probes have been saved in the file %s .\n", eqtlinfo->_include.size(), epiName.c_str());
    }
    void write_s2s_besd(char* outFileName, eqtlInfo* eqtlinfo, bool tosmrflag)
    {
        FILE * besd;
        string besdName = string(outFileName)+".besd";
        if(fopen_checked(&besd, besdName.c_str(),"wb")) {
            LOGPRINTF("ERROR: in opening file %s to write .\n", besdName.c_str());
            TERMINATE();
        }
        if(eqtlinfo->_valNum!=0)
        {
            if(eqtlinfo->_snpNum==eqtlinfo->_esi_include.size() && eqtlinfo->_probNum==eqtlinfo->_include.size())
            {
                if(tosmrflag)
                {
                    uint64_t valNum=eqtlinfo->_val.size();
                    vector<int> ten_ints(RESERVEDUNITS);
                    ten_ints[0]=SMR_SPARSE_3;
                    ten_ints[1]=eqtlinfo->_sampleNum;
                    ten_ints[2]=(int)eqtlinfo->_esi_include.size();
                    ten_ints[3]=(int)eqtlinfo->_include.size();
                    for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;

                    vector<uint64_t> cols;
                    vector<uint32_t> rowid;
                    cols.resize(2*eqtlinfo->_probNum+1);
                    rowid.resize(valNum);
                    cols[0]=0;
                    for(int i=1;i<eqtlinfo->_cols.size();i++)
                    {
                        long end=eqtlinfo->_cols[i];
                        long start=eqtlinfo->_cols[i-1];
                        long num=end-start;
                        long betanum=num>>1;
                        cols[2*i-1]=start+betanum;
                        cols[2*i]=end;
                        long rstart=start>>1;
                        memcpy(&rowid[start],&eqtlinfo->_rowid[rstart],betanum*sizeof(uint32_t));
                        memcpy(&rowid[start+betanum],&eqtlinfo->_rowid[rstart],betanum*sizeof(uint32_t));
                    }

                    if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
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
                    if (fwrite_checked(&rowid[0],rowid.size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&eqtlinfo->_val[0],eqtlinfo->_val.size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }

                }
                else
                {
                    uint64_t valNum=eqtlinfo->_val.size();
                    vector<int> ten_ints(RESERVEDUNITS);
                    ten_ints[0]=OSCA_SPARSE_1;
                    ten_ints[1]=eqtlinfo->_sampleNum;
                    ten_ints[2]=(int)eqtlinfo->_esi_include.size();
                    ten_ints[3]=(int)eqtlinfo->_include.size();
                    for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                    
                    if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&valNum,sizeof(uint64_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&eqtlinfo->_cols[0],eqtlinfo->_cols.size()*sizeof(uint64_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&eqtlinfo->_rowid[0],eqtlinfo->_rowid.size()*sizeof(uint32_t), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                    if (fwrite_checked(&eqtlinfo->_val[0],eqtlinfo->_val.size()*sizeof(float), besd))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }

                }
            }
            else
            {
                //not often here. it would be used when _include is changed after reading the besd. for example, matrix shrinking
                vector<uint64_t> cols;
                vector<uint32_t> rowid;
                vector<float> val;
                map<int, int > _incld_id_map;
                map<int, int>::iterator iter;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) {
                        LOGPRINTF("Error: Duplicated SNP IDs found: %s.\n",eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]].c_str());
                        TERMINATE();
                    }
                    size = _incld_id_map.size();
                }
                vector<float> sestmp;
                vector<uint32_t> ridtmp;
                if(tosmrflag)
                {
                    cols.resize((eqtlinfo->_probNum<<1)+1);
                    cols[0]=0;
                    for(int i=0;i<eqtlinfo->_include.size();i++)
                    {
                        int proid=eqtlinfo->_include[i];
                        uint64_t pos=eqtlinfo->_cols[proid];
                        uint64_t pos1=eqtlinfo->_cols[proid+1];
                        uint64_t num=(pos1-pos)>>1;
                        uint64_t colpos=pos>>1;
                        uint64_t real_num=0;
                        sestmp.clear();
                        for(int j=0;j<num;j++)
                        {
                            float beta=eqtlinfo->_val[pos+j];
                            float se=eqtlinfo->_val[pos+j+num];
                            int rid=eqtlinfo->_rowid[colpos+j];
                            iter=_incld_id_map.find(rid);
                            if(iter!=_incld_id_map.end())
                            {
                                int sid=iter->second;
                                rowid.push_back(sid);
                                val.push_back(beta);
                                sestmp.push_back(se);
                                ridtmp.push_back(sid);
                                real_num++;
                            }
                        }
                        for(int j=0;j<sestmp.size();j++) {
                            val.push_back(sestmp[j]);
                            rowid.push_back(ridtmp[j]);
                        }
                        cols[(i<<1)+1]=real_num+cols[i<<1];
                        cols[i+1<<1]=(real_num<<1)+cols[i<<1];
                    }
                }
                else
                {
                    cols.resize(eqtlinfo->_include.size()+1);
                    cols[0]=0;
                    for(int i=0;i<eqtlinfo->_include.size();i++)
                    {
                        int proid=eqtlinfo->_include[i];
                        uint64_t pos=eqtlinfo->_cols[proid];
                        uint64_t pos1=eqtlinfo->_cols[proid+1];
                        uint64_t num=(pos1-pos)>>1;
                        uint64_t colpos=pos>>1;
                        uint64_t real_num=0;
                        sestmp.clear();
                        for(int j=0;j<num;j++)
                        {
                            float beta=eqtlinfo->_val[pos+j];
                            float se=eqtlinfo->_val[pos+j+num];
                            int rid=eqtlinfo->_rowid[colpos+j];
                            iter=_incld_id_map.find(rid);
                            if(iter!=_incld_id_map.end())
                            {
                                int sid=iter->second;
                                rowid.push_back(sid);
                                val.push_back(beta);
                                sestmp.push_back(se);
                                real_num++;
                            }
                        }
                        for(int j=0;j<sestmp.size();j++) val.push_back(sestmp[j]);
                        cols[i+1]=(real_num<<1)+cols[i];
                    }
                }
                uint32_t fileType=OSCA_SPARSE_1;
                if(tosmrflag) fileType=SMR_SPARSE_3;
                vector<int> ten_ints(RESERVEDUNITS);
                ten_ints[0]=fileType;
                ten_ints[1]=eqtlinfo->_sampleNum;
                ten_ints[2]=(int)eqtlinfo->_esi_include.size();
                ten_ints[3]=(int)eqtlinfo->_include.size();
                for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                uint64_t valNum=val.size();
                if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
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
                if (fwrite_checked(&rowid[0],rowid.size()*sizeof(uint32_t), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                if (fwrite_checked(&val[0],val.size()*sizeof(float), besd))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }

            }
        }
        else
        {
            LOGPRINTF("ERROR: No information to write.\n");
            TERMINATE();
        }
        
        fclose (besd);
        LOGPRINTF("eQTL summary statistics have been saved in the file %s .\n", besdName.c_str());
    }
    void extract_prb_sparse(FILE* fptr, uint64_t pid, uint64_t probnum,vector<uint32_t> &row_ids, vector<float> &betases)
    {
        if(pid>probnum) {
            LOGPRINTF("ERROR: probe index %llu is larger than the totoal probe number %llu.\n", pid, probnum);
            TERMINATE();
        }
        row_ids.clear();
        betases.clear();
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        if(indicator==SMR_SPARSE_3F || indicator==SMR_SPARSE_3)
        {
            int infoLen=sizeof(uint32_t);
            if(indicator==SMR_SPARSE_3)
            {
                infoLen=RESERVEDUNITS*sizeof(int);
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), fptr);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    LOGPRINTF("The sample size is %d.\n",ss);
                }
                delete[] indicators;
            }
            
            uint64_t colNum=(probnum<<1)+1;
            uint64_t valnum=readuint64(fptr);
            fseek(fptr,(pid<<1)*sizeof(uint64_t),SEEK_CUR);
            uint64_t betaStart=readuint64(fptr);
            uint64_t seStart=readuint64(fptr);
            long num=seStart-betaStart;
            if(num>0)
            {
                row_ids.resize(num);
                betases.resize(2*num);
                uint64_t rowSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                uint64_t valSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valnum*sizeof(uint32_t);
                fseek(fptr, rowSTART+betaStart*sizeof(uint32_t), SEEK_SET);
                fread(&row_ids[0], sizeof(uint32_t),num,fptr);
                fseek(fptr,valSTART+betaStart*sizeof(float),SEEK_SET);
                fread(&betases[0],sizeof(float), 2*num,fptr);
                
            }
        }
        else if(indicator==OSCA_SPARSE_1) {
            
            int infoLen=RESERVEDUNITS*sizeof(int);
            int length=(RESERVEDUNITS-1)*sizeof(int);
            char* indicators=new char[length];
            fread(indicators, sizeof(int),(RESERVEDUNITS-1), fptr);
            int* tmp=(int *)indicators;
            int ss=*tmp++;
            if(ss!=-9)
            {
                LOGPRINTF("The sample size is %d.\n",ss);
            }
            delete[] indicators;
            
            uint64_t colNum=probnum+1;
            int64_t valnum=readuint64(fptr);
            fseek(fptr,pid*sizeof(uint64_t),SEEK_CUR);
            uint64_t betaStart=readuint64(fptr);
            uint64_t betanextStart=readuint64(fptr);
            long num=(betanextStart-betaStart)/2;
            if(num>0)
            {
                row_ids.resize(num);
                betases.resize(2*num);
                uint64_t rowSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                uint64_t valSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valnum*sizeof(uint32_t)/2;
                fseek(fptr, rowSTART+betaStart*sizeof(uint32_t)/2, SEEK_SET);
                fread(&row_ids[0], sizeof(uint32_t),num,fptr);
                fseek(fptr,valSTART+betaStart*sizeof(float),SEEK_SET);
                fread(&betases[0],sizeof(float), 2*num,fptr);
            }
        }
    }
    void extract_prb_dense(FILE* fptr,  uint64_t pid, uint64_t epinum,uint64_t esinum, vector<float> &betases)
    {
        if(epinum==0 || esinum==0) {
            LOGPRINTF("ERROR: .epi file or .esi file is empty. please check.\n");
            TERMINATE();
        }
        if(pid>epinum) {
            LOGPRINTF("ERROR: probe index %llu is larger than the totoal probe number %llu.\n", pid, epinum);
            TERMINATE();
        }
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        if(indicator==SMR_DENSE_1 || indicator==SMR_DENSE_3) {
            if(indicator==SMR_DENSE_3)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), fptr);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    LOGPRINTF("The sample size is %d.\n",ss);
                }
                delete[] indicators;
            }
            int infoLen=sizeof(uint32_t);
            if(indicator==SMR_DENSE_3) infoLen=RESERVEDUNITS*sizeof(int);
            
            betases.resize(2*esinum);
            fseek(fptr,((pid*esinum)<<3)+infoLen, SEEK_SET);
            fread(&betases[0], sizeof(float),2*esinum,fptr);
        }
        else if(indicator==OSCA_DENSE_1)
        {
            int length=(RESERVEDUNITS-1)*sizeof(int);
            char* indicators=new char[length];
            fread(indicators, sizeof(int),(RESERVEDUNITS-1), fptr);
            int* tmp=(int *)indicators;
            int ss=*tmp++;
            if(ss!=-9)
            {
                LOGPRINTF("The sample size is %d.\n",ss);
            }
            delete[] indicators;
            int infoLen=RESERVEDUNITS*sizeof(int);
            betases.resize(2*esinum);
            for(int i=0;i<esinum;i++)
            {
                fseek(fptr,infoLen+2*(i*epinum+pid)*sizeof(float), SEEK_SET);
                betases[i]=readfloat(fptr);
                betases[i+esinum]=readfloat(fptr);
            }
        }
    }
    void write_d2d_besd(char* outFileName, eqtlInfo* eqtlinfo, char* inputname, bool stdprb)
    {
        FILE * besdout;
        string besdName = string(outFileName)+".besd";
        if(fopen_checked(&besdout, besdName.c_str(),"wb")) {
            LOGPRINTF("ERROR: in opening file %s to write .\n", besdName.c_str());
            TERMINATE();
        }
        FILE * besdin;
        if(fopen_checked(&besdin, inputname,"rb")) {
            LOGPRINTF("ERROR: in opening file %s to read .\n", inputname);
            TERMINATE();
        }
        uint32_t filetype = SMR_DENSE_3;
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=eqtlinfo->_sampleNum;
        ten_ints[2]=(int)eqtlinfo->_esi_include.size();
        ten_ints[3]=(int)eqtlinfo->_include.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;

        if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), besdout))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
            TERMINATE();
        }

        if(eqtlinfo->_include.size()<(eqtlinfo->_probNum>>1))
        {
            vector<float> betases;
            vector<float> betasestmp;
            long realnum=eqtlinfo->_esi_include.size();
            betasestmp.resize(2*realnum);
            LOGPRINTF("The number of probes to save (%ld) is less than half size of the total probes (%llu).\n",eqtlinfo->_include.size(),eqtlinfo->_probNum);
            LOGPRINTF("Probe-by-probe strategy is exploited.\n");
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                betases.clear();
                extract_prb_dense(besdin,  eqtlinfo->_include[i], eqtlinfo->_probNum,eqtlinfo->_snpNum, betases);
                if(stdprb)
                {
                    double prbvar_sqrt=sqrt(eqtlinfo->_epi_var[eqtlinfo->_include[i]]);
                    for(int j=0;j<betases.size();j++) if(abs(betases[j]+9)>1e-6) betases[j]/=prbvar_sqrt; //would be updated. in case that beta is -9 but se not
                    
                }
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
                {
                    for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                    {
                        betasestmp[j]=betases[eqtlinfo->_esi_include[j]];
                        betasestmp[j+realnum]=betases[eqtlinfo->_esi_include[j]+eqtlinfo->_snpNum];
                    }
                    
                   
                    if (fwrite_checked(&betasestmp[0],betasestmp.size()*sizeof(float), besdout))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                }
                else
                {
                    if (fwrite_checked(&betases[0],betases.size()*sizeof(float), besdout))
                    {
                        LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                        TERMINATE();
                    }
                }
            }
        }
        else
        {
            LOGPRINTF("The number of probes to save (%ld) is more than half size of the total probes (%llu).\n",eqtlinfo->_include.size(),eqtlinfo->_probNum);
            vector<float> betase4w;
            long realnum=eqtlinfo->_esi_include.size();
            betase4w.resize(2*realnum);
            fseek(besdin,0L,SEEK_SET);
            uint32_t indicator=readuint32(besdin);
            if(indicator==SMR_DENSE_1 || indicator==SMR_DENSE_3)
            {
                if(indicator==SMR_DENSE_3)
                {
                    int length=(RESERVEDUNITS-1)*sizeof(int);
                    char* indicators=new char[length];
                    fread(indicators, sizeof(int),(RESERVEDUNITS-1), besdin);
                    int* tmp=(int *)indicators;
                    int ss=*tmp++;
                    if(ss!=-9)
                    {
                        LOGPRINTF("The sample size is %d.\n",ss);
                    }
                    delete[] indicators;
                }
                char* readBuf;
                long sizebesd=2*eqtlinfo->_probNum*eqtlinfo->_snpNum*sizeof(float);
                uint64_t rBufsize=allocReserved(&readBuf,sizebesd);
                float* fptr=(float *)readBuf;
                if(rBufsize==sizebesd) // read the whole
                {
                    LOGPRINTF("Reading the whole file %s into the memory.\n",inputname);
                    if(fread(readBuf, 1,rBufsize, besdin)!=rBufsize)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", inputname);
                        TERMINATE();
                    }
                    
                    for(int i=0; i<eqtlinfo->_include.size();i++)
                    {
                        if(stdprb)
                        {
                            int pos=eqtlinfo->_include[i];
                            //printf("Standardizing probe %s ...",eqtlinfo->_epi_prbID[pos].c_str());
                            double prbvar_sqrt=sqrt(eqtlinfo->_epi_var[pos]);
                            for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                            {
                                *(fptr+2*pos*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j])/=prbvar_sqrt;
                                double tmpval=*(fptr+(2*pos+1)*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j]);
                                if(abs(tmpval+9)>1e-6) *(fptr+(2*pos+1)*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j])/=prbvar_sqrt;
                            }
                        }

                        if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
                        {
                            for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                            {
                                int pos=eqtlinfo->_include[i];
                                betase4w[j]=*(fptr+2*pos*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j]);
                                betase4w[j+realnum]=*(fptr+(2*pos+1)*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j]);
                            }
                            if (fwrite_checked(&betase4w[0],betase4w.size()*sizeof(float), besdout))
                            {
                                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                TERMINATE();
                            }
                        }
                        else
                        {
                            if (fwrite_checked(fptr+2*eqtlinfo->_include[i]*eqtlinfo->_snpNum,2*eqtlinfo->_snpNum*sizeof(float), besdout))
                            {
                                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                TERMINATE();
                            }

                        }
                    }
                }
                else
                {
                    LOGPRINTF("Reading the file %s in batches ...\n",inputname);
                    uint64_t prb_per_time=rBufsize/(2*eqtlinfo->_snpNum*sizeof(float));
                    uint64_t prb_lower_boundary=0, prb_upper_boundary=0;
                    int i=0;
                    while (!feof(besdin)) {
                        uint64_t nSize = fread(readBuf, sizeof(float), prb_per_time*2*eqtlinfo->_snpNum, besdin);
                        uint64_t nPrb=nSize/(2*eqtlinfo->_snpNum);
                        prb_upper_boundary+=nPrb;
                        for(; i<eqtlinfo->_include.size();i++)
                        {
                            int pid=eqtlinfo->_include[i];
                            if(pid>=prb_upper_boundary) break;
                            long idcur=pid-prb_lower_boundary;
                            if(stdprb)
                            {
                                //printf("Standardizing probe %s ...",eqtlinfo->_epi_prbID[pid].c_str());
                                double prbvar_sqrt=sqrt(eqtlinfo->_epi_var[pid]);
                                for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                                {
                                    *(fptr+2*idcur*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j])/=prbvar_sqrt;
                                    double tmpval=*(fptr+(2*idcur+1)*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j]);
                                    if(abs(tmpval+9)>1e-6) *(fptr+(2*idcur+1)*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j])/=prbvar_sqrt;
                                }
                            }

                            if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
                            {
                                for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                                {
                                    betase4w[j]=*(fptr+2*idcur*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j]);
                                    betase4w[j+realnum]=*(fptr+(2*idcur+1)*eqtlinfo->_snpNum+eqtlinfo->_esi_include[j]);
                                }
                                if (fwrite_checked(&betase4w[0],betase4w.size()*sizeof(float), besdout))
                                {
                                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                    TERMINATE();
                                }
                            }
                            else
                            {
                                if (fwrite_checked(fptr+2*idcur*eqtlinfo->_snpNum,2*eqtlinfo->_snpNum*sizeof(float), besdout))
                                {
                                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                    TERMINATE();
                                }
                            }
                        }
                        prb_lower_boundary=prb_upper_boundary;
                    }
                }
                deallocReserved(&readBuf, rBufsize);
            }
            else if(indicator==OSCA_DENSE_1)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), besdin);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    LOGPRINTF("The sample size is %d.\n",ss);
                }
                delete[] indicators;

                char* IOBuf;
                long sizebesd=2*eqtlinfo->_probNum*eqtlinfo->_snpNum*sizeof(float);
                uint64_t rBufsize=allocReserved(&IOBuf,sizebesd);
                float* floatptr=(float *)IOBuf;
                if(rBufsize==sizebesd) // read the whole
                {
                    LOGPRINTF("Reading the whole file %s into the memory.\n",inputname);
                    if(fread(IOBuf, 1,rBufsize, besdin)!=rBufsize)
                    {
                        LOGPRINTF("ERROR: File %s read failed!\n", inputname);
                        TERMINATE();
                    }
                    
                    for(int i=0; i<eqtlinfo->_include.size();i++)
                    {
                        if(stdprb)
                        {
                            int pid=eqtlinfo->_include[i];
                            //printf("Standardizing probe %s ...",eqtlinfo->_epi_prbID[pid].c_str());
                            double prbvar_sqrt=sqrt(eqtlinfo->_epi_var[pid]);
                            for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                            {
                                int pos=eqtlinfo->_esi_include[j];
                                *(floatptr+2*pos*eqtlinfo->_probNum+2*eqtlinfo->_include[i])/=prbvar_sqrt;
                                double tmpval=*(floatptr+2*pos*eqtlinfo->_probNum+2*eqtlinfo->_include[i]+1);
                                if(abs(tmpval+9)>1e-6) *(floatptr+2*pos*eqtlinfo->_probNum+2*eqtlinfo->_include[i]+1)/=prbvar_sqrt;
                            }
                        }

                        for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                        {
                            int pos=eqtlinfo->_esi_include[j];
                            betase4w[j]=*(floatptr+2*pos*eqtlinfo->_probNum+2*eqtlinfo->_include[i]);
                            betase4w[j+realnum]=*(floatptr+2*pos*eqtlinfo->_probNum+2*eqtlinfo->_include[i]+1);
                        }
                        if (fwrite_checked(&betase4w[0],betase4w.size()*sizeof(float), besdout))
                        {
                            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                            TERMINATE();
                        }
                    }
                     deallocReserved(&IOBuf, rBufsize);
                }
                else
                {
                    LOGPRINTF("Transforming the data to probe major ...\n");
                    LOGPRINTF("Probe-by-probe strategy is exploited.\n");
                    deallocReserved(&IOBuf, rBufsize);
                    vector<float> betases;
                    vector<float> betasestmp;
                    long realnum=eqtlinfo->_esi_include.size();
                    betasestmp.resize(2*realnum);
                    for(int i=0;i<eqtlinfo->_include.size();i++)
                    {
                        printf("%3.0f%%\r", 100.0*i/eqtlinfo->_include.size());
                        fflush(stdout);
                        if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
                        {
                            for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                            {
                                fseek(besdin,RESERVEDUNITS*sizeof(int)+2*(eqtlinfo->_esi_include[j]*eqtlinfo->_probNum+eqtlinfo->_include[i])*sizeof(float), SEEK_SET);
                                betasestmp[j]=readfloat(besdin);
                                betasestmp[j+realnum]=readfloat(besdin);
                            }
                            if(stdprb)
                            {
                                int pid=eqtlinfo->_include[i];
                                //printf("Standardizing probe %s ...",eqtlinfo->_epi_prbID[pid].c_str());
                                double prbvar_sqrt=sqrt(eqtlinfo->_epi_var[pid]);
                                for(int j=0;j<betasestmp.size();j++) if(abs(betasestmp[j]+9)>1e-6) betasestmp[j]/=prbvar_sqrt; // would be updated in case that beta is -9 but se not
                            }

                            if (fwrite_checked(&betasestmp[0],betasestmp.size()*sizeof(float), besdout))
                            {
                                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                TERMINATE();
                            }
                        }
                        else
                        {
                            betases.clear();
                            extract_prb_dense(besdin,  eqtlinfo->_include[i], eqtlinfo->_probNum,eqtlinfo->_snpNum, betases);
                            if(stdprb)
                            {
                                int pid=eqtlinfo->_include[i];
                                double prbvar_sqrt=sqrt(eqtlinfo->_epi_var[pid]);
                                for(int j=0;j<betases.size();j++) if(abs(betasestmp[j]+9)>1e-6) betases[j]/=prbvar_sqrt; // would be updated in case that beta is -9 but se not
                            }
                            if (fwrite_checked(&betases[0],betases.size()*sizeof(float), besdout))
                            {
                                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                                TERMINATE();
                            }
                        }
                    }

                }
               
            }
        }
        fclose (besdout);
        fclose (besdin);
        LOGPRINTF("eQTL summary statistics have been saved in the file %s .\n", besdName.c_str());
    }

    void extract_epi_by_chr(eqtlInfo* eqtlinfo, int prbchr)
    {
        vector<string> prblst;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr ) prblst.push_back(eqtlinfo->_epi_prbID[tmpint]);
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from chromosome %d.\n",  eqtlinfo->_include.size(), prbchr);
    }
    void extract_smr_probe(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        update_map_kp(problist, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from %s.\n",  eqtlinfo->_include.size(), problstName.c_str());
    }
    void extract_smr_probe_by_gene(eqtlInfo* eqtlinfo, string genelistName)
    {
        vector<string> genelist;
        string msg="genes";
        read_msglist(genelistName, genelist,msg);
        
        vector<string> prblst;
        for(int i=0;i<genelist.size();i++)
        {
            string tmpname1=genelist[i];
            for(int j=0;j<eqtlinfo->_include.size();j++)
            {
                string tmpname2=eqtlinfo->_epi_gene[eqtlinfo->_include[j]];
                if(tmpname1==tmpname2)  prblst.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[j]]);
                else
                {
                    vector<string> substrs;
                    uint64_t tmpnum=split_string_skip(tmpname2,substrs,",;",0);
                    if(tmpnum>1)
                        for(int k=0;k<tmpnum;k++)
                            if(tmpname1==substrs[k])
                            {
                                prblst.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[j]]);
                                break;
                            }
                }
            }
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from %s.\n",  eqtlinfo->_include.size(), genelistName.c_str());
    }
    void extract_smr_probe(eqtlInfo* eqtlinfo, string prbname, int prbWind)
    {
        map<string, int>::iterator iter;
        iter=eqtlinfo->_probe_name_map.find(prbname);
        if(iter==eqtlinfo->_probe_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  prbname.c_str());
            TERMINATE();
        }
        long idx=iter->second;
        int prbbp=eqtlinfo->_epi_bp[idx];
        int prbchr=eqtlinfo->_epi_chr[idx];
        int upbound=prbbp+prbWind*1000;
        int tmpint=prbbp-prbWind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<string> prblst;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr && eqtlinfo->_epi_bp[tmpint]>=lowbound && eqtlinfo->_epi_bp[tmpint]<=upbound) prblst.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from the region %d Kb around %s.\n",  eqtlinfo->_include.size(), prbWind, prbname.c_str());
    }
    void extract_smr_single_probe(eqtlInfo* eqtlinfo, string prbname)
    {
        map<string, int>::iterator iter;
        iter=eqtlinfo->_probe_name_map.find(prbname);
        if(iter==eqtlinfo->_probe_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  prbname.c_str());
            TERMINATE();
        }
        long idx=iter->second;
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        eqtlinfo->_include.push_back((int)idx);
        eqtlinfo->_probe_name_map.insert(pair<string,int>(prbname,(int)idx));
        LOGPRINTF( "%s is extracted.\n",prbname.c_str());
    }
    void extract_smr_probe(eqtlInfo* eqtlinfo, string fromprbname, string toprbname)
    {
        map<string, int>::iterator iter;
        iter=eqtlinfo->_probe_name_map.find(fromprbname);
        if(iter==eqtlinfo->_probe_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  fromprbname.c_str());
            TERMINATE();
        }
        int fromprbbp=eqtlinfo->_epi_bp[iter->second];
        int prbchr=eqtlinfo->_epi_chr[iter->second];
        
        iter=eqtlinfo->_probe_name_map.find(toprbname);
        if(iter==eqtlinfo->_probe_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  toprbname.c_str());
            TERMINATE();
        }
        int toprbbp=eqtlinfo->_epi_bp[iter->second];
        int toprbchr=eqtlinfo->_epi_chr[iter->second];
        if(toprbchr != prbchr)
        {
            LOGPRINTF("ERROR: probe %s and probe %s are not from the same chromosome.\n", fromprbname.c_str(), toprbname.c_str());
            TERMINATE();
        }
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<string> prblst;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr && eqtlinfo->_epi_bp[tmpint]>=fromprbbp && eqtlinfo->_epi_bp[tmpint]<=toprbbp) prblst.push_back(eqtlinfo->_epi_prbID[tmpint]);
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from probe %s to probe %s.\n",  eqtlinfo->_include.size(), fromprbname.c_str(), toprbname.c_str());
    }

    void extract_smr_probe(eqtlInfo* eqtlinfo, int fromprbkb, int toprbkb, int chr)
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
        for(int i=0;eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==chr && eqtlinfo->_epi_bp[tmpint]>=fromprbbp && eqtlinfo->_epi_bp[tmpint]<=toprbbp) prblst.push_back(eqtlinfo->_epi_prbID[tmpint]);
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from probe bp %d to probe bp %d.\n",  eqtlinfo->_include.size(), fromprbkb, toprbkb);
    }
    void extract_smr_probe_by_single_gene(eqtlInfo* eqtlinfo, string genename)
    {
        vector<string> prblst;
        for(int j=0;j<eqtlinfo->_include.size();j++)
        {
            string tmpname2=eqtlinfo->_epi_gene[eqtlinfo->_include[j]];
            if(genename==tmpname2)  prblst.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[j]]);
            else
            {
                vector<string> substrs;
                uint64_t tmpnum=split_string_skip(tmpname2,substrs,",;",0);
                if(tmpnum>1)
                    for(int k=0;k<tmpnum;k++)
                        if(genename==substrs[k])
                        {
                            prblst.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[j]]);
                            break;
                        }
            }
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are extracted from gene %s.\n",  eqtlinfo->_include.size(), genename.c_str());
    }
    void exclude_smr_probe(eqtlInfo* eqtlinfo, string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        long pre_num=eqtlinfo->_include.size();
        update_map_rm(problist, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld probes are excluded from %s and there are %ld probe remaining.\n",  pre_num - eqtlinfo->_include.size(), problstName.c_str(),eqtlinfo->_include.size());
    }
    void exclude_smr_single_probe(eqtlInfo* eqtlinfo, string prbname)
    {
        vector<string> problist;
        problist.push_back(prbname);
        update_map_rm(problist, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("Probe %s are excluded and there are %ld probe remaining.\n", prbname.c_str(), eqtlinfo->_include.size());
    }

    void smr_epi_man(eqtlInfo* eqtlinfo,char* problstName,char* problst2exclde,char* genelistName, int chr, int probchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde)
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
        if(probchr>0)
        {
            extract_epi_by_chr(eqtlinfo, probchr);
        } else if (chr>0) {
            extract_epi_by_chr(eqtlinfo, chr);
        }
        
        if(problstName != NULL || genelistName != NULL)
        {
            if(problstName != NULL) extract_smr_probe(eqtlinfo, problstName);
            if(genelistName != NULL) extract_smr_probe_by_gene(eqtlinfo, genelistName);
        }
        else if(prbwindFlag)
        {
            if(prbname==NULL)
            {
                LOGPRINTF("ERROR: Please identify the probe name by --probe when using --probe-wind.\n");
                TERMINATE();
            }
            extract_smr_probe(eqtlinfo, prbname, prbWind);
        }
        else if(prbname!=NULL)
        {
            extract_smr_single_probe(eqtlinfo, prbname);
        }
        else if(fromprbname!=NULL)
        {
            if(toprbname==NULL)
            {
                LOGPRINTF("ERROR: Please identify the probe name by --to-probe.\n");
                TERMINATE();
            }
            extract_smr_probe(eqtlinfo, fromprbname, toprbname);
        }
        else if(fromprbkb>=0)
        {
            if(fromprbkb>=0 && chr==0 ) {
                LOGPRINTF("ERROR: Please identify the chromosome by --chr.\n");
                TERMINATE();
            }
            if(toprbkb<0)
            {
                LOGPRINTF("ERROR: probe BP can't be negative.\n");
                TERMINATE();
            }
            extract_smr_probe(eqtlinfo, fromprbkb, toprbkb,chr);
        }
        else if(genename!=NULL)
        {
            extract_smr_probe_by_single_gene(eqtlinfo, genename);
        }
        
        if(problst2exclde!=NULL)
        {
            exclude_smr_probe(eqtlinfo, problst2exclde);
        }
        else if(probe2exclde!=NULL)
        {
            exclude_smr_single_probe(eqtlinfo,probe2exclde);
        }
    }
    void extract_smr_esi_by_chr(eqtlInfo* eqtlinfo, int snpchr)
    {
        vector<int> newIcld;
        eqtlinfo->_snp_name_map.clear();
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==snpchr ) {
                newIcld.push_back(tmpint);
                eqtlinfo->_snp_name_map.insert(pair<string,int>(eqtlinfo->_esi_rs[tmpint],tmpint));
            }
        }
        eqtlinfo->_esi_include.swap(newIcld);
        LOGPRINTF("%ld SNPs are extracted from chromosome %d. \n",eqtlinfo->_esi_include.size(),snpchr);
    }
    void extract_smr_eqtl_snp(eqtlInfo* eqtlinfo, string snporprb, int Wind, string msg)
    {
        string logstr;
        
        int bp=-9;
        int chr=-9;
        map<string, int>::iterator iter;
        if(msg=="SNP")
        {
            iter=eqtlinfo->_snp_name_map.find(snporprb);
            if(iter==eqtlinfo->_snp_name_map.end())
            {
                LOGPRINTF("ERROR: Can't find %s %s.\n", msg.c_str(), snporprb.c_str());
                TERMINATE();
            }
            bp=eqtlinfo->_esi_bp[iter->second];
            chr=eqtlinfo->_esi_chr[iter->second];
            if(chr < 0) {
                LOGPRINTF("ERROR: Missing Chromosome found of %s %s.\n", msg.c_str(), snporprb.c_str());
                TERMINATE();
            }
            if(bp < 0) {
                LOGPRINTF("ERROR: Missing BP found of %s %s.\n", msg.c_str(), snporprb.c_str());
                TERMINATE();
            }
        }
        else if(msg=="probe")
        {
            iter=eqtlinfo->_probe_name_map.find(snporprb);
            if(iter==eqtlinfo->_probe_name_map.end())
            {
                LOGPRINTF("ERROR: Can't find %s %s.\n", msg.c_str(), snporprb.c_str());
                TERMINATE();
            }
            bp=eqtlinfo->_epi_bp[iter->second];
            chr=eqtlinfo->_epi_chr[iter->second];
            if(chr < 0) {
                LOGPRINTF("ERROR: Missing Chromosome found of %s %s.\n", msg.c_str(), snporprb.c_str());
                TERMINATE();
            }
            if(bp < 0) {
                LOGPRINTF("ERROR: Missing BP found of %s %s.\n", msg.c_str(), snporprb.c_str());
                TERMINATE();
            }
        }
        
        int upbound=bp+Wind*1000;
        int tmpint=bp-Wind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<int> newIcld;
        eqtlinfo->_snp_name_map.clear();
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==chr && eqtlinfo->_esi_bp[tmpint]>=lowbound && eqtlinfo->_esi_bp[tmpint]<=upbound) {
                newIcld.push_back(tmpint);
                eqtlinfo->_snp_name_map.insert(pair<string,int>(eqtlinfo->_esi_rs[tmpint],tmpint));
            }
        }
        eqtlinfo->_esi_include.swap(newIcld);
        LOGPRINTF("%ld SNPs are extracted from the region: %d Kb around %s %s.\n", eqtlinfo->_esi_include.size() ,Wind,msg.c_str(), snporprb.c_str());
    }
    void extract_smr_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
    {
        vector<string> snplist;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        update_map_kp(snplist, eqtlinfo->_snp_name_map, eqtlinfo->_esi_include);
        LOGPRINTF("%ld SNPs are extracted from %s.\n", eqtlinfo->_esi_include.size() ,snplstName.c_str() );
    }
    void extract_smr_eqtl_single_snp(eqtlInfo* eqtlinfo, string snprs)
    {
        map<string, int>::iterator iter;
        iter=eqtlinfo->_snp_name_map.find(snprs);
        if(iter==eqtlinfo->_snp_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find SNP %s.\n", snprs.c_str());
            TERMINATE();
        }
        long idx=iter->second;
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include.push_back((int)idx);
        eqtlinfo->_snp_name_map.clear();
        eqtlinfo->_snp_name_map.insert(pair<string,int>(snprs,idx));
        cout << snprs << " is extracted. " << endl;
    }
    void extract_smr_eqtl_snp(eqtlInfo* eqtlinfo, string fromsnprs, string tosnprs)
    {
        map<string, int>::iterator iter;
        iter=eqtlinfo->_snp_name_map.find(fromsnprs);
        if(iter==eqtlinfo->_snp_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  fromsnprs.c_str());
            TERMINATE();
        }
        int fromsnpbp=eqtlinfo->_esi_bp[iter->second];
        int snpchr=eqtlinfo->_esi_chr[iter->second];
        
        iter=eqtlinfo->_snp_name_map.find(tosnprs);
        if(iter==eqtlinfo->_snp_name_map.end())
        {
            LOGPRINTF("ERROR: Can't find probe %s.\n",  tosnprs.c_str());
            TERMINATE();
        }
        int tosnpbp=eqtlinfo->_esi_bp[iter->second];
        int tosnpchr=eqtlinfo->_esi_chr[iter->second];
        if(tosnpchr != snpchr)
        {
            LOGPRINTF("ERROR: SNP %s and SNP %s are not from the same chromosome.\n", fromsnprs.c_str(), tosnprs.c_str());
            TERMINATE();
        }
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
        
        vector<string> snplst;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==snpchr && eqtlinfo->_esi_bp[tmpint]>=fromsnpbp && eqtlinfo->_esi_bp[tmpint]<=tosnpbp) snplst.push_back(eqtlinfo->_esi_rs[tmpint]);
        }
        update_map_kp(snplst, eqtlinfo->_snp_name_map, eqtlinfo->_esi_include);
        LOGPRINTF("%ld SNPs are extracted from SNP %s to SNP %s.\n",  eqtlinfo->_esi_include.size(), fromsnprs.c_str(), tosnprs.c_str());
    }
    void extract_smr_eqtl_snp(eqtlInfo* eqtlinfo, int chr, int fromsnpkb, int tosnpkb)
    {
        int fromsnpbp=fromsnpkb*1000;
        int tosnpbp=tosnpkb*1000;
        
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
        
        vector<int> newIcld;
         eqtlinfo->_snp_name_map.clear();
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if( eqtlinfo->_esi_chr[tmpint]==chr &&eqtlinfo->_esi_bp[tmpint]>=fromsnpbp && eqtlinfo->_esi_bp[tmpint]<=tosnpbp) {
                newIcld.push_back(tmpint);
                eqtlinfo->_snp_name_map.insert(pair<string,int>(eqtlinfo->_esi_rs[tmpint],tmpint));
            }
        }
        eqtlinfo->_esi_include.swap(newIcld);
        LOGPRINTF("%ld SNPs are extracted from SNP BP:  %dKb to SNP BP: %dKb on chromosome %d.\n",eqtlinfo->_esi_include.size(),fromsnpkb,tosnpkb, chr);
    }
    void exclude_smr_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
    {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplstName, snplist,msg);
        long pre_num=eqtlinfo->_esi_include.size();
        update_map_rm(snplist, eqtlinfo->_snp_name_map, eqtlinfo->_esi_include);
        LOGPRINTF("%ld SNPs are excluded from %s and there are %ld SNPs remaining.\n",  pre_num - eqtlinfo->_esi_include.size(), snplstName.c_str(),eqtlinfo->_esi_include.size());
    }
    void exclude_smr_single_snp(eqtlInfo* eqtlinfo, string snprs2exclde)
    {
        vector<string> snplist;
        snplist.push_back(snprs2exclde);
        update_map_rm(snplist, eqtlinfo->_snp_name_map, eqtlinfo->_esi_include);
        LOGPRINTF("SNP %s are excluded and there are %ld probe remaining.\n", snprs2exclde.c_str(), eqtlinfo->_esi_include.size());
    }

    void smr_esi_man(eqtlInfo* eqtlinfo,char* snplstName, char* snplst2exclde, int chr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, int cis_itvl,const char* prbname,char* snprs2exclde)
    {
        string logstr;
        int flags4snp=0;
        if(snplstName != NULL) flags4snp++;
        if(snprs != NULL) flags4snp++;
        if(fromsnprs!=NULL) flags4snp++;
        if(fromsnpkb>=0) flags4snp++;
        if(flags4snp>1)
        {
            LOGPRINTF("WARNING: Flags for SNPs in this section are mutual exclusive. The priority order (from high to low) is: --extract-snp, --snp-wind, --snp, --from(to)--snp, --from(to)-snp-kb.\n");
        }
        if(snpchr>0)
        {
            extract_smr_esi_by_chr(eqtlinfo, snpchr);
        }
        else if(chr>0)
        {
            extract_smr_esi_by_chr(eqtlinfo, chr);
        }
        
        if(prbname!=NULL && cis_flag)
        {
            extract_smr_eqtl_snp(eqtlinfo, prbname, cis_itvl, "probe"); // extract cis eQTLs
        }
        else if (snplstName != NULL) extract_smr_eqtl_snp(eqtlinfo, snplstName);
        else if (snpwindFlag)
        {
            if(snprs==NULL)
            {
                LOGPRINTF("ERROR: please specify the SNP name by --snp when using --snp-wind.\n");
                TERMINATE();
            }
            extract_smr_eqtl_snp(eqtlinfo, snprs, snpWind, "SNP");
        }
        else if(snprs!=NULL)
        {
            extract_smr_eqtl_single_snp(eqtlinfo, snprs);
        }
        else if(fromsnprs!=NULL)
        {
            if(tosnprs==NULL)
            {
                LOGPRINTF("ERROR: please specify the SNP name by --to-snp.\n");
                TERMINATE();
            }
            extract_smr_eqtl_snp(eqtlinfo, fromsnprs, tosnprs);
        }
        else if(fromsnpkb>=0)
        {
            
            if(fromsnpkb>=0 && chr==0 && snpchr==0) {
                LOGPRINTF("ERROR: please specify the chromosome by --snp-chr or --chr.\n");
                TERMINATE();
            }
            
            if(tosnpkb<0)
            {
                LOGPRINTF("ERROR: SNP BP can't be negative.\n");
                TERMINATE();
            }
            if(snpchr!=0) extract_smr_eqtl_snp(eqtlinfo, snpchr, fromsnpkb, tosnpkb);
            else if (chr!=0) extract_smr_eqtl_snp(eqtlinfo, chr, fromsnpkb, tosnpkb);
        }
        if(snplst2exclde!=NULL) exclude_smr_eqtl_snp(eqtlinfo,snplst2exclde);
        if(snprs2exclde!=NULL) exclude_smr_single_snp(eqtlinfo,snprs2exclde);
        
    }
    
    void read_gwas_data(gwasData* gdata, char* gwasFileName)
    {
        bool warnnullfreq=false;
        FILE* gwasFile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum=8;
        if(fopen_checked(&gwasFile, gwasFileName,"r")) { TERMINATE(); }
        
        LOGPRINTF("Reading GWAS summary data from %s.\n",gwasFileName);
        gdata->_include.clear();
        gdata->snpName.clear();
        gdata->snpBp.clear();
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        gdata->_snp_name_map.clear();
        if(fgets(Tbuf, MAX_LINE_SIZE, gwasFile)) // the header
        {
            if(Tbuf[0]=='\0')
            {
                LOGPRINTF("ERROR: the first row of the file %s is empty.\n",gwasFileName);
                TERMINATE();
            }
            split_string(Tbuf, strlist, " \t\n");
            to_upper(strlist[0]);
            if(strlist[0]!="SNP") {
                LOGPRINTF("ERROR: %s should have headers that start with \"snp\".\n", gwasFileName);
                TERMINATE();
            }
        }
        while(fgets(Tbuf, MAX_LINE_SIZE, gwasFile))
        {
            if(Tbuf[0]=='\0') {
                LOGPRINTF("ERROR: Line %u is blank.\n", line_idx+2);
                TERMINATE();
            }
             split_str(Tbuf,strlist,0);
            if(strlist.size()!=colnum)
            {
                LOGPRINTF("ERROR: Line %u has %d items.\n", line_idx+2,colnum);
                TERMINATE();
            }
            if(strlist[0]=="NA" || strlist[0]=="na"){
                LOGPRINTF("ERROR: SNP name is \'NA\' in row %d.\n", line_idx+2);
                TERMINATE();
            }
            gdata->_snp_name_map.insert(pair<string,int>(strlist[0],line_idx));
            if(gdata->_snp_name_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate SNP : %s.\n", strlist[0].c_str());
                TERMINATE();
            }
            gdata->snpName.push_back(strlist[0]);
            if(strlist[1]=="NA" || strlist[1]=="na"){
                LOGPRINTF("ERROR: allele1 is \'NA\' in row %d.\n", line_idx+2);
                TERMINATE();
            }
            to_upper(strlist[1]);
            gdata->allele_1.push_back(strlist[1]);
            
            if(strlist[2]=="NA" || strlist[2]=="na"){
                LOGPRINTF("ERROR: allele2 is \'NA\' in row %d.\n", line_idx+2);
                TERMINATE();
            }
            to_upper(strlist[2]);
            gdata->allele_2.push_back(strlist[2]);
            
            if(strlist[3]=="NA" || strlist[3]=="na"){
                if(!warnnullfreq){
                    warnnullfreq=true;
                    LOGPRINTF("WARNING: frequency is \'NA\' in one or more rows.\n");
                }
                gdata->freq.push_back(-9);
            } else {
                gdata->freq.push_back(atof(strlist[3].c_str()));
            }
            
            if(strlist[4]=="NA" || strlist[4]=="na"){
                LOGPRINTF("WARNING: effect size is \'NA\' in row %d.\n", line_idx+2);
                gdata->byz.push_back(0);
            } else {
                gdata->byz.push_back(atof(strlist[4].c_str()));
            }
            if(strlist[5]=="NA" || strlist[5]=="na"){
                LOGPRINTF("WARNING: standard error is \'NA\' in row %d.\n", line_idx+2);
                gdata->seyz.push_back(-9);
            } else {
                gdata->seyz.push_back(atof(strlist[5].c_str()));
            }
            
            gdata->pvalue.push_back(atof(strlist[6].c_str()));
            
            if(strlist[7]=="NA" || strlist[7]=="na"){
                gdata->splSize.push_back(-9);
            } else {
                gdata->splSize.push_back(atoi(strlist[7].c_str()));
            }
            gdata->_include.push_back(line_idx);
            line_idx++;
        }
        gdata->snpNum=gdata->_include.size();
        LOGPRINTF("GWAS summary data of %ld SNPs to be included from %s.\n" ,gdata->snpNum ,gwasFileName);
        fclose(gwasFile);
    }
    
    void read_ewas_data(gwasData* gdata, char* ewasFileName)
    {
        FILE* ewasFile=NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum=8;
        if(fopen_checked(&ewasFile, ewasFileName,"r")) { TERMINATE(); }
        
        LOGPRINTF("Reading EWAS summary data from %s.\n",ewasFileName);
        gdata->_include.clear();
        gdata->snpName.clear();
        gdata->snpBp.clear();
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        gdata->_snp_name_map.clear();
        if(fgets(Tbuf, MAX_LINE_SIZE, ewasFile)) // the header
        {
            if(Tbuf[0]=='\0')
            {
                LOGPRINTF("ERROR: the first row of the file %s is empty.\n",ewasFileName);
                TERMINATE();
            }
        }
        while(fgets(Tbuf, MAX_LINE_SIZE, ewasFile))
        {
            if(Tbuf[0]=='\0') {
                LOGPRINTF("ERROR: Line %u is blank.\n", line_idx+2);
                TERMINATE();
            }
            split_str(Tbuf,strlist,0);
            if(strlist.size()!=colnum && strlist.size()!=colnum+1)
            {
                LOGPRINTF("ERROR: Line %u has %d items.\n", line_idx+2,colnum);
                TERMINATE();
            }
            if(strlist[1]=="NA" || strlist[1]=="na"){
                LOGPRINTF("ERROR: probe name is \'NA\' in row %d.\n", line_idx+2);
                TERMINATE();
            }
            gdata->_snp_name_map.insert(pair<string,int>(strlist[1],line_idx));
            if(gdata->_snp_name_map.size()==line_idx)
            {
                LOGPRINTF("ERROR: Duplicate probe : %s.\n", strlist[1].c_str());
                TERMINATE();
            }
            gdata->snpName.push_back(strlist[1]);
            int chr=0;
            if(strlist[0]=="NA" || strlist[0]=="na") chr=-9;
            else if (strlist[0]=="X" || strlist[0]=="x") chr=23;
            else if (strlist[0]=="Y" || strlist[0]=="y") chr=24;
            else chr=atoi(strlist[0].c_str());
            gdata->freq.push_back(chr);
            gdata->snpBp.push_back(atoi(strlist[2].c_str()));
            gdata->allele_1.push_back(strlist[3]);
            gdata->allele_2.push_back(strlist[4]);
            if(strlist[5]=="NA" || strlist[5]=="na"){
                LOGPRINTF("WARNING: effect size is \'NA\' in row %d.\n", line_idx+2);
                gdata->byz.push_back(0);
            } else {
                gdata->byz.push_back(atof(strlist[5].c_str()));
            }
            if(strlist[6]=="NA" || strlist[6]=="na"){
                LOGPRINTF("WARNING: standard error is \'NA\' in row %d.\n", line_idx+2);
                gdata->seyz.push_back(-9);
            } else {
                gdata->seyz.push_back(atof(strlist[6].c_str()));
            }
            gdata->pvalue.push_back(atof(strlist[7].c_str()));
            if(strlist.size()==colnum+1)
            {
                if(strlist[8]=="NA" || strlist[8]=="na"){
                    gdata->splSize.push_back(-9);
                } else {
                    gdata->splSize.push_back(atoi(strlist[8].c_str()));
                }
            }
            else gdata->splSize.push_back(-9);
            
            gdata->_include.push_back(line_idx);
            line_idx++;
        }
        gdata->snpNum=gdata->_include.size();
        LOGPRINTF("EWAS summary data of %ld probes to be included from %s.\n" ,gdata->snpNum ,ewasFileName);
        fclose(ewasFile);
    }
    
    void get_shrink_null(eqtlInfo* eqtlinfo,vector<string> &nullprbs, vector<string> &nullsnps)
    {
        //working on sparse besd
        if(eqtlinfo->_val.size()==0)
        {
            LOGPRINTF("ERROR: No information obtained.\n");
            TERMINATE();
        }
        if(eqtlinfo->_include.size()!=eqtlinfo->_probNum || eqtlinfo->_esi_include.size()!=eqtlinfo->_snpNum)
        {
            LOGPRINTF("ERROR: Shrink function needs to update. please report to Futao Zhang (futao.zhang@imb.uq.edu.au) for further update.\n");
            TERMINATE();
        }

        nullprbs.clear();
        nullsnps.clear();
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int proid=eqtlinfo->_include[i];
            string prb=eqtlinfo->_epi_prbID[proid];
            uint64_t pos=eqtlinfo->_cols[proid];
            uint64_t pos1=eqtlinfo->_cols[proid+1];
            if(pos1==pos) {
                nullprbs.push_back(prb);
                continue;
            }
        }
        vector<bool> miss;
        miss.resize(eqtlinfo->_snpNum);
        for(int i=0;i<eqtlinfo->_rowid.size();i++)
        {
            int rid=eqtlinfo->_rowid[i];
            if(rid<eqtlinfo->_snpNum) {
                miss[rid]=true;
            }
            else
            {
                LOGPRINTF("ERROR: row id >= eqtlinfo._snpNum. please report.\n");
                TERMINATE();
            }
        }
        for(int i=0;i<miss.size();i++)
        {
            if(!miss[i])
            {
                string snprs=eqtlinfo->_esi_rs[i];
                nullsnps.push_back(snprs);
            }
        }
        
    }
    
    
    void read_probevarfile(eqtlInfo* eqtlinfo, char* vpFileName)
    {
        LOGPRINTF("Reading variance information from %s.\n",vpFileName);
        FILE* flptr=NULL;
        if(fopen_checked(&flptr, vpFileName,"r")) { TERMINATE(); }
        
        int lineNum(0), hit(0);
        vector<string> strlist;
        map<string, int>::iterator iter;
        eqtlinfo->_epi_var.resize(eqtlinfo->_epi_prbID.size());
        while(fgets(Tbuf, MAX_LINE_SIZE, flptr))
        {
            if(Tbuf[0]=='\0') {
                LOGPRINTF("ERROR: Line %u is blank.\n", lineNum+2);
                TERMINATE();
            }
            split_str(Tbuf,strlist,0);
            if(strlist.size()!=2)
            {
                LOGPRINTF("ERROR: Line %u has %d items.\n", lineNum+2,2);
                TERMINATE();
            }
            iter=eqtlinfo->_probe_name_map.find(strlist[0]);
            if(iter!=eqtlinfo->_probe_name_map.end())
            {
                eqtlinfo->_epi_var[iter->second]=atof(strlist[1].c_str());
                hit++;
            }
            lineNum++;
        }
        fclose(flptr);
        if(hit!=eqtlinfo->_include.size())
        {
            LOGPRINTF("ERROR: %d of %ld probes have the variance.\n",hit,eqtlinfo->_include.size());
            TERMINATE();
        }
        LOGPRINTF("%d probes variance information to be included from %s.\n",hit,vpFileName);
    }
    void extract_sqtl_probe(eqtlInfo* eqtlinfo,int tsk_ttl,int tsk_id)
    {
        map<string, int> gene_count_map;
        map<string, int>::iterator iter;
        vector<string> gene;
        vector< vector< string> > idx;
        int ids = 0;
        for(int i=0; i<eqtlinfo->_include.size();i++)
        {
            string gn = eqtlinfo->_epi_gene[eqtlinfo->_include[i]];
            to_upper(gn);
            if( gn != "NA" )
            {
                iter = gene_count_map.find(gn);
                if(iter == gene_count_map.end())
                {
                    gene_count_map.insert(pair<string, int>(gn, ids++));
                    gene.push_back(gn);
                    vector<string> tmpidx;
                    tmpidx.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
                    idx.push_back(tmpidx);
                }
                else
                {
                    int curid = iter->second;
                    idx[curid].push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
                }
            }
            else
            {
                if(loud) {
                    LOGPRINTF("NA gene name found with the probe %s. \n", eqtlinfo->_epi_prbID[eqtlinfo->_include[i]].c_str());
                }
            }
        }
        vector<int> inids;
        for(int i=0;i<idx.size();i++)
            if(idx[i].size()>1) inids.push_back(i);
        
        LOGPRINTF("%ld genes with more than 1 transcript are extracted.\n", inids.size());
        int unicount=(int)inids.size();
        if(tsk_ttl>1)
        {
            if(tsk_ttl>unicount) {
                LOGPRINTF("WARNING: Total task number %d is larger than the gene number %d.\nThe task number shrinks to %d.\n", tsk_ttl, unicount,unicount);
                tsk_ttl=unicount;
                if(tsk_id>tsk_ttl) {
                    LOGPRINTF("Task id %d is larger than shrinked total task number %d. Task terminated.\n", tsk_id,tsk_ttl);
                    TERMINATE();
                }
            }
            long num_per_tsk=ceill(unicount/(1.0*tsk_ttl));
            
            long fromidx=(tsk_id-1)*num_per_tsk;
            if(fromidx>=unicount) {
                LOGPRINTF("The whole analysis can be accomplished by the tasks ahead. This task %d would be no more needed.\n", tsk_id);
                TERMINATE();
            }
            long toidx=tsk_id*num_per_tsk-1;
            if(toidx >= unicount) toidx=unicount-1;
            LOGPRINTF("The %ldth gene to the %ldth gene of total %d genes are extracted from this sub-task.\n", fromidx+1,toidx+1,unicount);
            vector<int> subinids;
            for(long i=fromidx;i<=toidx;i++)
                subinids.push_back(inids[i]);
            inids.swap(subinids);
        }
        vector<string> prblst;
        for(int i=0;i<inids.size();i++)
        {
            if(idx[inids[i]].size()<=1)
            {
                LOGPRINTF("Error: here is a bug. please report to futao.zhang@imb.uq.edu.au\n");
            }
            for(int j=0;j<idx[inids[i]].size();j++)
            {
                prblst.push_back(idx[inids[i]][j]);
            }
        }
        update_map_kp(prblst, eqtlinfo->_probe_name_map, eqtlinfo->_include);
        LOGPRINTF("%ld genes are extracted from task id %d of total task number %d (total %ld probes).\n",  inids.size(), tsk_id, tsk_ttl, eqtlinfo->_include.size());
    }
    
    bool pcc(MatrixXd &PCC, float* buffer_beta,float* buffer_se,long snpnum, long cohortnum, double pmecs, int nmecs)
    {
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
                    double sei=buffer_se[i*snpnum+k];
                    double betai=buffer_beta[i*snpnum+k];
                    double sej=buffer_se[j*snpnum+k];
                    double betaj=buffer_beta[j*snpnum+k];
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
}

//
//  l3_smr.cpp
//  osc
//
//  Created by Futao Zhang on 20/11/2017.
//  Copyright © 2017 Futao Zhang. All rights reserved.
//

#include "l3_smr.hpp"
namespace SMR {

    int comp_epi(const void *a,const void *b){ return (((*(smr_probeinfo *)a).probechr>(*(smr_probeinfo *)b).probechr) || ( ((*(smr_probeinfo *)a).probechr ==(*(smr_probeinfo *)b).probechr) && ((*(smr_probeinfo *)a).bp > (*(smr_probeinfo *)b).bp) ))?1:-1; }
    int comp_esi(const void *a,const void *b){ return (((*(smr_snpinfo *)a).snpchr >(*(smr_snpinfo *)b).snpchr) || ( ((*(smr_snpinfo *)a).snpchr ==(*(smr_snpinfo *)b).snpchr) && ((*(smr_snpinfo *)a).bp > (*(smr_snpinfo *)b).bp) ))?1:-1; }

    void query_besd(char* outFileName,char* beqtlFileName, char* snplstName, char* snplst2exclde, char* problstName,char* problst2exclde, char* genelistName, double plookup, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl, char* probe2exclde, char* snprs2exclde)
    {
        string logstr;
        bool readwhole=false;
        int flag4chr=0, miss=0;
        if(chr > 0) flag4chr++;
        if(prbchr > 0 || snpchr > 0) flag4chr++;
        if(flag4chr==2)
        {
            chr=-9;
            LOGPRINTF("WARNING: --chr is not surpposed to use together with --probe-chr or --snp-chr. --chr will be disabled.\n");
        }

        eqtlInfo eqtlinfo;
        LOGPRINTF("\nReading eQTL summary data...\n");
        if(beqtlFileName == NULL)
         {
            LOGPRINTF("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
            TERMINATE();
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,beqtlFileName,strlen(beqtlFileName)+1);
        char* suffix=inputname+strlen(beqtlFileName);
        memcpy(suffix,".epi",5);

        read_smr_epifile(&eqtlinfo, inputname);
        smr_epi_man(&eqtlinfo, problstName, problst2exclde, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2exclde);
        memcpy(suffix,".esi",5);
        read_smr_esifile(&eqtlinfo, inputname);
        smr_esi_man(&eqtlinfo, snplstName, snplst2exclde,chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname,snprs2exclde);
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
        get_BesdHeaders(inputname,headers);
        int indicator = headers[0];
        if(indicator==SMR_DENSE_1 || indicator==SMR_DENSE_3 || indicator == OSCA_DENSE_1 )
        {
            if(((eqtlinfo._probNum*eqtlinfo._snpNum*sizeof(float)) >> 30) <1) readwhole=true;
        }
        if(indicator==SMR_SPARSE_3 || indicator==SMR_SPARSE_3F || indicator== OSCA_SPARSE_1) readwhole=true;
        if(readwhole)
        {
            LOGPRINTF("Loading the whole file %s into memory...\n",beqtlFileName);
            read_smr_besdfile(&eqtlinfo, inputname);
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                LOGPRINTF("No data included from %s under current condition.\n",beqtlFileName);
                TERMINATE();
            }

            vector<int> out_esi_id;
            vector<int> out_epi_id;
            vector<float> out_beta;
            vector<float> out_se;
            vector<double> out_pval;
            if(eqtlinfo._valNum==0)
            {
                for(int i=0;i<eqtlinfo._probNum;i++)
                {
                    for(int j=0;j<eqtlinfo._snpNum;j++)
                    {
                        double beta=eqtlinfo._bxz[i][j];
                        double se=eqtlinfo._sexz[i][j];
                        if(abs(se+9)<1e-6) {
                            miss++;
                            continue;
                        }
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=plookup)
                        {
                            out_esi_id.push_back(j);
                            out_epi_id.push_back(i);
                            out_beta.push_back(beta);
                            out_se.push_back(se);
                            out_pval.push_back(pxz);
                        }
                    }
                }
            }
            else
            {
                if(eqtlinfo._val.size()==0)
                {
                    LOGPRINTF("Error: No data extracted from the input, please check.\n");
                    TERMINATE();
                }

                for(int i=0;i<eqtlinfo._include.size();i++)
                {
                    int proid=eqtlinfo._include[i];
                    uint64_t pos=eqtlinfo._cols[proid];
                    uint64_t pos1=eqtlinfo._cols[proid+1];
                    uint64_t num=(pos1-pos)>>1;
                    uint64_t colpos=pos>>1;
                    for(int j=0;j<num;j++)
                    {
                        double beta=eqtlinfo._val[pos+j];
                        double se=eqtlinfo._val[pos+j+num];
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=plookup)
                        {
                            out_esi_id.push_back(eqtlinfo._rowid[colpos+j]);
                            out_epi_id.push_back(proid);
                            out_beta.push_back(beta);
                            out_se.push_back(se);
                            out_pval.push_back(pxz);
                        }
                    }
                }
            }

            ofstream smr(outFileName);
            if (!smr) {
                LOGPRINTF("Error: can not open the file %s to save!\n",outFileName);
                TERMINATE();
            }

            smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' <<"Freq"<<'\t'<< "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"Orientation"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';

            for (int i = 0;i <out_esi_id.size(); i++) {
                smr<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_bp[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<((eqtlinfo._esi_freq[out_esi_id[i]]+9>1e-6)?atos(eqtlinfo._esi_freq[out_esi_id[i]]):"NA")<<'\t'<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_chr[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_bp[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_gene[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_orien[out_epi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
            }

            smr.close();
            if(miss) {
                LOGPRINTF("%d missing values found.\n",miss);
            }
            LOGPRINTF("Extracted results of %ld SNPs have been saved in the file %s.\n",out_esi_id.size(),outFileName);

        }
        else
        {
            LOGPRINTF("Loading the file %s probe by probe ...\n",beqtlFileName);
            LOGPRINTF("Only the Non-missing SNPs would be writen to the output file.\n");
            long wcount=0;
            FILE* besd=NULL;
            if(fopen_checked(&besd,inputname,"rb")) {
                TERMINATE();
            }
            ofstream smr(outFileName);
            if (!smr) {
                LOGPRINTF("Error: can not open the file %s to save!\n",outFileName);
                TERMINATE();
            }
            smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' <<"Freq"<<'\t'<< "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"Orientation"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
            vector<float> betases;
            for(int i=0;i<eqtlinfo._include.size();i++)
            {
                int epid=eqtlinfo._include[i];
                betases.clear();
                extract_prb_dense(besd,  eqtlinfo._include[i], eqtlinfo._probNum,eqtlinfo._snpNum, betases);
                for(int j=0;j<eqtlinfo._esi_include.size();j++)
                {
                    int esid=eqtlinfo._esi_include[j];
                    double beta=betases[esid];
                    double se=betases[esid+eqtlinfo._snpNum];
                    if(abs(se+9)>1e-6)
                    {
                        double z=beta/se;
                        double p=pchisq(z*z, 1);
                        smr<<eqtlinfo._esi_rs[esid]<<'\t'<<eqtlinfo._esi_chr[esid]<<'\t'<<eqtlinfo._esi_bp[esid]<<'\t'<<eqtlinfo._esi_allele1[esid]<<'\t'<<eqtlinfo._esi_allele2[esid]<<'\t'<<((eqtlinfo._esi_freq[esid]+9>1e-6)?atos(eqtlinfo._esi_freq[esid]):"NA")<<'\t'<<eqtlinfo._epi_prbID[epid]<<'\t'<<eqtlinfo._epi_chr[epid]<<'\t'<<eqtlinfo._epi_bp[epid]<<'\t'<<eqtlinfo._epi_gene[epid]<<'\t'<<eqtlinfo._epi_orien[epid]<<'\t'<<beta<<'\t'<<se<<'\t'<<p<< '\n';
                        wcount++;
                    } else {
                        miss++;
                    }

                }
            }
            smr.close();
            if(miss) {
                LOGPRINTF("%d missing values found.\n",miss);
            }
            LOGPRINTF("Extracted results of %ld SNPs have been saved in the file %s.\n",wcount,outFileName);
        }

    }
    void check_besds_format( vector<string> &besds, vector<int> &format, vector<int> &smpsize) {

        char inputname[FNAMESIZE];
        format.clear();
        smpsize.clear();
        vector<int> headers;
        for(int i=0;i<besds.size();i++)
        {
            string tmpstr=besds[i]+".besd";
            memcpy(inputname,tmpstr.c_str(),tmpstr.length()+1);
            get_BesdHeaders(inputname, headers);
            format.push_back(headers[0]);
            if(headers[0]==SMR_DENSE_1 || headers[0]==SMR_SPARSE_3F) smpsize.push_back(-9);
            else smpsize.push_back(headers[1]);
        }
    }
    void combine_epi(vector<smr_probeinfo> &probeinfo, vector<string> &besds, vector<uint64_t> &nprb,char* problstName, char* problst2exclde, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename, char* probe2exclde)
    {
        long counter = 0;
        map<string, int> prb_map;
        map<string, int> prbbp_map;
        map<string, int>::iterator iter;
        nprb.clear();
        char inputname[FNAMESIZE];
        for (int i = 0; i < besds.size(); i++)
        {
            eqtlInfo etmp;
            memcpy(inputname,besds[i].c_str(),besds[i].length()+1);
            char* suffix=inputname+besds[i].length();
            memcpy(suffix,".epi",5);
            read_smr_epifile(&etmp, inputname);
            smr_epi_man(&etmp, problstName, problst2exclde, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2exclde);
            nprb.push_back(etmp._probNum);
            for (int j = 0; j<etmp._probNum; j++)
            {
                string crsbpstr=etmp._epi_prbID[j]+":"+atos(etmp._epi_bp[j]);
                prb_map.insert(pair<string, int>(etmp._epi_prbID[j].c_str(), counter));
                prbbp_map.insert(pair<string, int>(crsbpstr.c_str(), counter));
                if(prb_map.size() != prbbp_map.size())
                {
                    LOGPRINTF("ERROR: inconsistent position for the probe %s  in different .epi files. Please check.\n", etmp._epi_prbID[j].c_str()) ;
                    TERMINATE();
                }

                if (counter < prb_map.size())
                {
                    smr_probeinfo probinfotmp;
                    counter=prb_map.size();
                    probinfotmp.probechr=etmp._epi_chr[j];
                    strcpy2(&probinfotmp.probeId, etmp._epi_prbID[j]);
                    probinfotmp.bp=etmp._epi_bp[j];
                    probinfotmp.gd=etmp._epi_gd[j];
                    strcpy2(&probinfotmp.genename, etmp._epi_gene[j]);
                    probinfotmp.orien=etmp._epi_orien[j];
                    probinfotmp.bfilepath=NULL;
                    probinfotmp.esdpath=NULL;
                    probinfotmp.ptr=new int[besds.size()];
                    for(int k=0;k<besds.size();k++){
                        if(i==k){
                            probinfotmp.ptr[k]=j;
                        } else {
                            probinfotmp.ptr[k]=-9;
                        }
                    }
                    probeinfo.push_back(probinfotmp);

                } else {
                    iter=prb_map.find(etmp._epi_prbID[j]);
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
        LOGPRINTF("Total %ld probes to be included from %ld epi files.\n",probeinfo.size(),besds.size());
    }
    void combine_esi(vector<smr_snpinfo> &snpinfo, vector<string> &besds, vector<uint64_t> &nsnp,char* snplstName, char* snplst2exclde, int chr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb, bool smpwindFlag, char* snprs2excld)
    {
        long ex_counter = 0;
        snpinfo.clear();
        vector<string> ex_snp;
        map<string, int> in_map;
        map<string, int> ex_map;
        map<string, int>::iterator iter;
        nsnp.clear();
        char inputname[FNAMESIZE];
        LOGPRINTF("\nPerforming Allele checking. This step could be a little long....\n");
        for (int i = 0; i < besds.size(); i++)
        {
            eqtlInfo etmp;
            memcpy(inputname,besds[i].c_str(),besds[i].length()+1);
            char* suffix=inputname+besds[i].length();
            memcpy(suffix,".esi",5);
            read_smr_esifile(&etmp, inputname);
            smr_esi_man(&etmp, snplstName, snplst2exclde, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, smpwindFlag, false,  0, NULL, snprs2excld);
            nsnp.push_back(etmp._snpNum);
            for (int j = 0; j<etmp._esi_include.size(); j++)
            {
                if(ex_map.size()>0) {
                    iter=ex_map.find(etmp._esi_rs[j]);
                    if(iter!=ex_map.end()) continue;
                }

                iter=in_map.find(etmp._esi_rs[j]);
                if(iter==in_map.end())
                {
                    in_map.insert(pair<string, int>(etmp._esi_rs[j].c_str(), snpinfo.size()));//the second is snpinfo id

                    smr_snpinfo snpinfotmp;
                    snpinfotmp.snpchr=etmp._esi_chr[j];
                    strcpy2(&snpinfotmp.snprs, etmp._esi_rs[j]);
                    snpinfotmp.bp=etmp._esi_bp[j];
                    snpinfotmp.gd=etmp._esi_gd[j];
                    strcpy2(&snpinfotmp.a1, etmp._esi_allele1[j]);
                    strcpy2(&snpinfotmp.a2, etmp._esi_allele2[j]);
                    snpinfotmp.freq=etmp._esi_freq[j];

                    snpinfotmp.rstr=new int[besds.size()];
                    snpinfotmp.revs=new bool[besds.size()];
                    for(int k=0;k<besds.size();k++){
                        if(i==k){
                            snpinfotmp.rstr[k]=j;
                            snpinfotmp.revs[k]=false;
                        } else {
                            snpinfotmp.rstr[k]=-9;
                            snpinfotmp.revs[k]=false;
                        }
                    }
                    snpinfo.push_back(snpinfotmp);
                } else {

                    if((snpinfo[iter->second].snpchr != etmp._esi_chr[j]) ||(snpinfo[iter->second].bp != etmp._esi_bp[j]))
                    {
                        // SNP in one esi file has BP1 but in another esi file has BP2
                        LOGPRINTF("ERROR: inconsistent chromosome or position for the SNP %s in different .epi files. Please check.\n", etmp._esi_rs[j].c_str()) ;
                        TERMINATE();
                    }
                    string a1=etmp._esi_allele1[j];
                    string a2=etmp._esi_allele2[j];
                    string a3=snpinfo[iter->second].a1;
                    string a4=snpinfo[iter->second].a2;
                    if(a1==a3 && a2==a4) {
                        snpinfo[iter->second].rstr[i]=j;
                    }
                    else if(a1==a4 && a2==a3 ){
                        snpinfo[iter->second].rstr[i]=j;
                        snpinfo[iter->second].revs[i]=true;
                    }
                    else {
                        //allele check failed. the SNP would be removed
                        ex_snp.push_back(etmp._esi_rs[j]);
                        ex_map.insert(pair<string,int>(etmp._esi_rs[j],ex_counter));
                        ex_counter=ex_map.size();
                        in_map.erase(iter->first);
                    }
                }
            }
        }
        //do not use multiple .erase() on vector, shuffle costs O(n*m), n is the vector size, m is the size of elememts to remove
        if(in_map.size() + ex_map.size() != snpinfo.size()){
            LOGPRINTF("ERROR: bugs found in allele check. please report.\n") ;
            TERMINATE();
        }
        long ttl_snp_common=snpinfo.size();
        FILE* failfptr=NULL;
        string failName=string(outfileName)+".failed.snp.list";
        if(fopen_checked(&failfptr, failName.c_str(),"w"))
        {
            LOGPRINTF("ERROR: failed in open file %s.\n",failName.c_str()) ;
            TERMINATE();
        }
        vector<smr_snpinfo> snpinfo_adj;
        snpinfo_adj.resize(in_map.size());
        int ids=0;
        for(int i=0; i<snpinfo.size();i++)
        {
            iter=in_map.find(snpinfo[i].snprs);
            if(iter!=in_map.end()) {
                snpinfo_adj[ids++] = snpinfo[i];
            } else {
                string snpstr=string(snpinfo[i].snprs) + '\n';
                if(fputs_checked(snpstr.c_str(),failfptr))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", failName.c_str());
                    TERMINATE();
                }
                //free the space of SNPs failed in allele check.
                if(snpinfo[i].a1) free2(&snpinfo[i].a1);
                if(snpinfo[i].a2) free2(&snpinfo[i].a2);
                if(snpinfo[i].snprs) free2(&snpinfo[i].snprs);
                if(snpinfo[i].rstr) free2(&snpinfo[i].rstr);
                if(snpinfo[i].revs) free2(&snpinfo[i].revs);
            }
        }
        snpinfo.swap(snpinfo_adj);
        snpinfo_adj.clear();
        fclose(failfptr);

        LOGPRINTF("Total %ld SNPs to be included from %ld esi files. %ld SNPs failed in allele check. %ld SNPs included in analysis.\n",ttl_snp_common,besds.size(),ex_snp.size(),in_map.size());
        LOGPRINTF("%ld SNPs that failed in allele check were saved in file %s.\n",ex_snp.size(),failName.c_str());
    }
    void extract_gwas_snp(gwasData* gdata, string snplstName)
    {
        vector<string> snplist;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        update_map_kp(snplist, gdata->_snp_name_map, gdata->_include);
        LOGPRINTF("%ld SNPs are extracted from %s.\n",gdata->_include.size(),snplstName.c_str());
    }
    void update_gwas(gwasData* gdata){

        bool hasBP=false;
        if(gdata->snpBp.size()>0) hasBP=true;
        vector<int> snpBp;
        if(hasBP) snpBp.resize(gdata->_include.size());
        vector<string> snpName(gdata->_include.size());
        vector<string> allele_1(gdata->_include.size());
        vector<string> allele_2(gdata->_include.size());
        vector<double> freq(gdata->_include.size());
        vector<double> byz(gdata->_include.size());
        vector<double> seyz(gdata->_include.size());
        vector<double> pvalue(gdata->_include.size());
        vector<int> splSize(gdata->_include.size());

        gdata->snpNum=gdata->_include.size();
        for(int i=0;i<gdata->_include.size();i++ )
        {
            snpName[i]=gdata->snpName[gdata->_include[i]];
            allele_1[i]=gdata->allele_1[gdata->_include[i]];
            allele_2[i]=gdata->allele_2[gdata->_include[i]];
            freq[i]=gdata->freq[gdata->_include[i]];
            byz[i]=gdata->byz[gdata->_include[i]];
            seyz[i]=gdata->seyz[gdata->_include[i]];
            pvalue[i]=gdata->pvalue[gdata->_include[i]];
            splSize[i]=gdata->splSize[gdata->_include[i]];
            if(hasBP) snpBp[i]=gdata->snpBp[gdata->_include[i]];
        }

        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();

        gdata->snpName.swap(snpName);
        gdata->allele_1.swap(allele_1);
        gdata->allele_2.swap(allele_2);
        gdata->freq.swap(freq);
        gdata->byz.swap(byz);
        gdata->seyz.swap(seyz);
        gdata->pvalue.swap(pvalue);
        gdata->splSize.swap(splSize);
        if(hasBP) gdata->snpBp=snpBp;

        gdata->_snp_name_map.clear();
        for(int i=0;i<gdata->snpNum;i++) {
            gdata->_include[i]=i;
            gdata->_snp_name_map.insert(pair<string, int>(gdata->snpName[i], i));
        }
        LOGPRINTF("The GWAS summary data is updated.\n");
    }
    void combine_gwas(vector<gwasinfo> &snpinfo, vector<string> &gwass, char* snplstName)
    {
        long ex_counter = 0;
        snpinfo.clear();
        vector<string> ex_snp;
        map<string, int> in_map;
        map<string, int> ex_map;
        map<string, int>::iterator iter;
        char inputname[FNAMESIZE];
        LOGPRINTF("\nPerforming Allele checking. This step could be a little long....\n");
        for (int i = 0; i < gwass.size(); i++)
        {
            gwasData gtmp;
            memcpy(inputname,gwass[i].c_str(),gwass[i].length()+1);
            read_gwas_data(&gtmp, inputname);
            if(snplstName!=NULL) {
                extract_gwas_snp(&gtmp, snplstName);
                update_gwas(&gtmp);
            }
            if(gtmp.snpNum==0) {
                LOGPRINTF("WARNING: No SNP included from the file %s.\n",inputname);
                continue;
            }
            for (int j = 0; j<gtmp.snpNum; j++)
            {
                if(ex_map.size()>0) {
                    iter=ex_map.find(gtmp.snpName[j]);
                    if(iter!=ex_map.end()) continue;
                }

                iter=in_map.find(gtmp.snpName[j]);
                if(iter==in_map.end())
                {
                    in_map.insert(pair<string, int>(gtmp.snpName[j].c_str(), snpinfo.size()));//the second is snpinfo id

                    gwasinfo snpinfotmp;
                    strcpy2(&snpinfotmp.snprs, gtmp.snpName[j]);
                    strcpy2(&snpinfotmp.a1, gtmp.allele_1[j]);
                    strcpy2(&snpinfotmp.a2, gtmp.allele_2[j]);
                    snpinfotmp.freq=gtmp.freq[j];
                    snpinfotmp.estn=gtmp.splSize[j];

                    snpinfotmp.beta =new float[gwass.size()];
                    snpinfotmp.se=new float[gwass.size()];
                    for(int k=0;k<gwass.size();k++){
                        if(i==k){
                            snpinfotmp.beta[k]=gtmp.byz[j];
                            snpinfotmp.se[k]=gtmp.seyz[j];
                        } else {
                            snpinfotmp.beta[k]=0;
                            snpinfotmp.se[k]=-9;
                        }
                    }
                    snpinfo.push_back(snpinfotmp);
                } else {
                    string a1=gtmp.allele_1[j];
                    string a2=gtmp.allele_2[j];
                    string a3=snpinfo[iter->second].a1;
                    string a4=snpinfo[iter->second].a2;
                    if(a1==a3 && a2==a4) {
                        snpinfo[iter->second].beta[i]=gtmp.byz[j];
                        snpinfo[iter->second].se[i]=gtmp.seyz[j];
                    }
                    else if(a1==a4 && a2==a3 ){
                        snpinfo[iter->second].beta[i]=-1.0*gtmp.byz[j];
                        snpinfo[iter->second].se[i]=gtmp.seyz[j];
                    }
                    else {
                        //allele check failed. the SNP would be removed
                        ex_snp.push_back(gtmp.snpName[j]);
                        ex_map.insert(pair<string,int>(gtmp.snpName[j],ex_counter));
                        ex_counter=ex_map.size();
                        in_map.erase(iter->first);
                    }
                    if(abs(gtmp.seyz[j]+9)>1e-6 && gtmp.splSize[j] != -9) snpinfo[iter->second].estn += gtmp.splSize[j];
                }
            }
        }
        //do not use multiple .erase() on vector, shuffle costs O(n*m), n is the vector size, m is the size of elememts to remove
        if(in_map.size() + ex_map.size() != snpinfo.size()){
            LOGPRINTF("ERROR: bugs found in allele check. please report.\n") ;
            TERMINATE();
        }
        long ttl_snp_common=snpinfo.size();
        FILE* failfptr=NULL;
        string failName=string(outfileName)+".failed.snp.list";
        if(fopen_checked(&failfptr, failName.c_str(),"w"))
        {
            LOGPRINTF("ERROR: failed in open file %s.\n",failName.c_str()) ;
            TERMINATE();
        }
        vector<gwasinfo> snpinfo_adj;
        snpinfo_adj.resize(in_map.size());
        int ids=0;
        for(int i=0; i<snpinfo.size();i++)
        {
            iter=in_map.find(snpinfo[i].snprs);
            if(iter!=in_map.end()) {
                snpinfo_adj[ids++] = snpinfo[i];
            } else {
                string snpstr=string(snpinfo[i].snprs) + '\n';
                if(fputs_checked(snpstr.c_str(),failfptr))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", failName.c_str());
                    TERMINATE();
                }
                //free the space of SNPs failed in allele check.
                if(snpinfo[i].a1) free2(&snpinfo[i].a1);
                if(snpinfo[i].a2) free2(&snpinfo[i].a2);
                if(snpinfo[i].snprs) free2(&snpinfo[i].snprs);
                if(snpinfo[i].beta) free2(&snpinfo[i].beta);
                if(snpinfo[i].se) free2(&snpinfo[i].se);
            }
        }
        snpinfo.swap(snpinfo_adj);
        snpinfo_adj.clear();
        fclose(failfptr);

        LOGPRINTF("Total %ld SNPs to be included from %ld GWAS summary files. %ld SNPs failed in allele check. %ld SNPs included in analysis.\n",ttl_snp_common,gwass.size(),ex_snp.size(),in_map.size());
        LOGPRINTF("%ld SNPs that failed in allele check were saved in file %s.\n",ex_snp.size(),failName.c_str());
    }
    void combine_ewas(vector<gwasinfo> &snpinfo, vector<string> &ewass, char* problstName)
    {
        snpinfo.clear();
        map<string, int> in_map;
        map<string, int>::iterator iter;
        char inputname[FNAMESIZE];
        for (int i = 0; i < ewass.size(); i++)
        {
            gwasData gtmp;
            memcpy(inputname,ewass[i].c_str(),ewass[i].length()+1);
            read_ewas_data(&gtmp, inputname);
            if(problstName!=NULL) {
                extract_gwas_snp(&gtmp, problstName);
                update_gwas(&gtmp);
            }
            if(gtmp.snpNum==0) {
                LOGPRINTF("WARNING: No probe included from the file %s.\n",inputname);
                continue;
            }
            for (int j = 0; j<gtmp.snpNum; j++)
            {
                iter=in_map.find(gtmp.snpName[j]);
                if(iter==in_map.end())
                {
                    in_map.insert(pair<string, int>(gtmp.snpName[j].c_str(), snpinfo.size()));//the second is snpinfo id

                    gwasinfo snpinfotmp;
                    strcpy2(&snpinfotmp.snprs, gtmp.snpName[j]);
                    strcpy2(&snpinfotmp.a1, gtmp.allele_1[j]);
                    strcpy2(&snpinfotmp.a2, gtmp.allele_2[j]);
                    snpinfotmp.freq=gtmp.freq[j];
                    snpinfotmp.estn=gtmp.splSize[j];
                    snpinfotmp.bp=gtmp.snpBp[j];

                    snpinfotmp.beta =new float[ewass.size()];
                    snpinfotmp.se=new float[ewass.size()];
                    for(int k=0;k<ewass.size();k++){
                        if(i==k){
                            snpinfotmp.beta[k]=gtmp.byz[j];
                            snpinfotmp.se[k]=gtmp.seyz[j];
                        } else {
                            snpinfotmp.beta[k]=0;
                            snpinfotmp.se[k]=-9;
                        }
                    }
                    snpinfo.push_back(snpinfotmp);
                } else {
                        snpinfo[iter->second].beta[i]=gtmp.byz[j];
                        snpinfo[iter->second].se[i]=gtmp.seyz[j];
                    if(abs(gtmp.seyz[j]+9)>1e-6 && gtmp.splSize[j] != -9) snpinfo[iter->second].estn += gtmp.splSize[j];
                }
            }
        }
        LOGPRINTF("Total %ld probes to be included from %ld EWAS summary files. %ld probes included in analysis.\n",snpinfo.size(),ewass.size(),in_map.size());
    }
    //Practical aspects of imputation-driven meta-analysis of genome-wide association studies. de Bakker PI1, Ferreira MA, Jia X, Neale BM, Raychaudhuri S, Voight BF.Hum Mol Genet. 2008 Oct 15;17(R2):R122-8. doi: 10.1093/hmg/ddn288.
    void meta_per_prob(float* buffer_beta,float* buffer_se, long snpnum, long cohortnum)
    {
        #pragma omp parallel for
        for(int j=0;j<snpnum;j++)
        {
            double numerator=0.0;
            double deno=0.0;
            int nmiss=0;
            for(int k=0;k<cohortnum;k++)
            {
                double se=buffer_se[k*snpnum+j];
                double beta=buffer_beta[k*snpnum+j];
                buffer_se[k*snpnum+j]=-9;
                if(abs(se+9)>1e-6){
                    double tmp2=se*se;
                    deno+=1/tmp2;
                    numerator+=beta/tmp2;
                    nmiss++;
                }
            }
            if(nmiss>0)
            {
                buffer_beta[j]=numerator/deno;
                buffer_se[j]=1/sqrt(deno);
            }
        }
    }
    void meta_gwas_fun(vector<gwasinfo> &snpinfo, long cohortnum, vector<int> &pairwise_comm, vector<int> &all_comm, bool outcom, bool allcom)
    {
        string filena=string(outfileName)+".betas";
        FILE* tmpfi= NULL;
        if(loud)
        {
            tmpfi = fopen(filena.c_str(),"w");
            if(!tmpfi)
            {
                LOGPRINTF("error open file.\n");
                TERMINATE();
            }
        }

        pairwise_comm.clear();
        all_comm.clear();
        long snpnum=snpinfo.size();
        #pragma omp parallel for
        for(int j=0;j<snpnum;j++)
        {
            double numerator=0.0;
            double deno=0.0;
            int nmiss=0;
            string str = atos(snpinfo[j].snprs) + '\t';
            for(int k=0;k<cohortnum;k++)
            {
                double se=snpinfo[j].se[k];
                double beta=snpinfo[j].beta[k];
                if(abs(se+9)>1e-6){
                    double tmp2=se*se;
                    deno+=1/tmp2;
                    numerator+=beta/tmp2;
                    nmiss++;
                    if(loud) str +=atos(beta) + '\t';
                }
            }
            if(nmiss>0)
            {
                *snpinfo[j].beta=numerator/deno;
                *snpinfo[j].se=1/sqrt(deno);
                if(loud) str += atos(numerator/deno) + '\n';
            }
            if(nmiss>1) pairwise_comm.push_back(j);
            if(nmiss==cohortnum) all_comm.push_back(j);
            if(loud) {
                if(outcom)
                {
                    if(nmiss>1) fputs(str.c_str(),tmpfi);
                }
                else if(allcom)
                {
                    if(nmiss==cohortnum) fputs(str.c_str(),tmpfi);
                }
                else
                {
                    fputs(str.c_str(),tmpfi);
                }
            }
        }

        if(loud) {
            if(tmpfi != NULL) fclose(tmpfi);
            LOGPRINTF("Beta values are saved in the file %s.\n",filena.c_str());
        }
    }

    void est_cor(MatrixXd &Cor,vector<gwasinfo> &snpinfo, long cohortnum, double pmecs)
    {
        // the off-diagonal can be larger than 1. DO NOT USE anymore in eQTL. no p-value threshold is needed
        long snpnum=snpinfo.size();
        //double zmecs=qchisq(pmecs,1);
        vector<double> se1sq,se2sq, beta1,beta2, beta1_2;
        double ve1, ve2;
        for( int i=0;i<cohortnum;i++)
            for(int j=i+1;j<cohortnum;j++)
            {
                se1sq.clear();
                se2sq.clear();
                beta1.clear();
                beta2.clear();
                beta1_2.clear();
                for(int k=0;k<snpnum;k++)
                {
                    double sei=snpinfo[k].se[i];
                    double betai=snpinfo[k].beta[i];
                    double sej=snpinfo[k].se[j];
                    double betaj=snpinfo[k].beta[j];
                    if(abs(sei+9)>1e-6 && abs(sej+9)>1e-6) {

                        se1sq.push_back(sei*sei);
                        se2sq.push_back(sej*sej);
                        beta1.push_back(betai);
                        beta2.push_back(betaj);
                        beta1_2.push_back(betai-betaj);

                        /*
                         //with p-value threshold
                        double zi=betai/sei;
                        double zj=betaj/sej;
                        zi*=zi;
                        zj*=zj;
                        if(zi < zmecs && zj < zmecs)
                        {
                            se1sq.push_back(sei*sei);
                            se2sq.push_back(sej*sej);
                            beta1.push_back(betai);
                            beta2.push_back(betaj);
                            beta1_2.push_back(betai-betaj);
                        }
                         */
                    }
                }
                ve1=mean(se1sq);
                ve2=mean(se2sq);
                Cor(i,j)=Cor(j,i)=(ve1+ve2-cov(beta1_2,beta1)+cov(beta1_2,beta2))/(2*sqrt(ve1*ve2));
            }
        for( int i=0;i<cohortnum;i++) Cor(i,i)=1;
    }


    void est_cor(MatrixXd &V, float* buffer_beta,float* buffer_se,long snpnum, long cohortnum, double pmecs)
    {
        // the off-diagonal can be larger than 1. DO NOT USE anymore.
        double zmecs=qchisq(pmecs,1);
        vector<double> se1sq,se2sq, beta1,beta2, beta1_2;
        double ve1, ve2;
        for( int i=0;i<cohortnum;i++)
            for(int j=i+1;j<cohortnum;j++)
            {
                se1sq.clear();
                se2sq.clear();
                beta1.clear();
                beta2.clear();
                beta1_2.clear();
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
                            se1sq.push_back(sei*sei);
                            se2sq.push_back(sej*sej);
                            beta1.push_back(betai);
                            beta2.push_back(betaj);
                            beta1_2.push_back(betai-betaj);
                        }
                    }
                }
                ve1=mean(se1sq);
                ve2=mean(se2sq);
                V(i,j)=V(j,i)=(ve1+ve2-cov(beta1_2,beta1)+cov(beta1_2,beta2))/(2*sqrt(ve1*ve2));
            }
        for( int i=0;i<cohortnum;i++) V(i,i)=1;
    }

    bool mecs_per_prob(float* buffer_beta,float* buffer_se, long snpnum, long cohortnum,double pmecs,vector<int> &noninvertible, vector<int> &negativedeno, int nmecs)
    {

        MatrixXd Corr(cohortnum,cohortnum);
        //LOGPRINTF("Estimate the cohort correlation using the beta values of pairwised common SNPs.\n");
        //LOGPRINTF("We exclude the significant common SNPs with a p-value threshold %e to estimate the correlation matrix.\n",pmecs);
        bool enoughsnp = pcc(Corr,buffer_beta,buffer_se,snpnum,cohortnum,pmecs, nmecs);
        if(!enoughsnp) return false;
        //cout<<Corr<<endl;
        #pragma omp parallel for
        for(int j=0;j<snpnum;j++)
        {
            vector<double> ses, betas;
            vector<int> keep;
            //MatrixXd Corr_work = Corr;
            MatrixXd Corr_work;
            int nmiss=0, miss=0;
            for(int k=0;k<cohortnum;k++)
            {
                double se=buffer_se[k*snpnum+j];
                double beta=buffer_beta[k*snpnum+j];
                buffer_se[k*snpnum+j]=-9;
                if(abs(se+9)>1e-6){
                    ses.push_back(se);
                    betas.push_back(beta);
                    keep.push_back(k);
                    nmiss++;
                } else {
                    //removeRow(Corr_work, k-miss);
                    //removeColumn(Corr_work, k-miss);
                    miss++;
                }
            }

            if(nmiss==1)
            {
                buffer_beta[j]=betas[0];
                buffer_se[j]=ses[0];
            }
            else if(nmiss>1)
            {
                if(nmiss<cohortnum) subMatrix_symm(Corr_work, Corr, keep);
                else if(nmiss==cohortnum) Corr_work = Corr;
                else {
                    LOGPRINTF("Can't happen. I can guarantee!\n");
                }
                VectorXd sev(ses.size());
                for(int k=0;k<ses.size();k++) sev(k)=ses[k];
                MatrixXd W=sev*sev.transpose();
                W=W.array()*Corr_work.array();
                bool determinant_zero=false;
                inverse_V(W,determinant_zero);
                if(determinant_zero) noninvertible.push_back(j);
                double deno=W.sum();
                if(deno<=0) {
                    negativedeno.push_back(j);
                } else {
                    VectorXd colsum=W.colwise().sum();
                    double numerator=0.0;
                    for(int k=0;k<betas.size();k++) numerator+=colsum(k)*betas[k];
                    buffer_beta[j]=numerator/deno;
                    buffer_se[j]=1/sqrt(deno);
                }

            }
        }
        return true;
    }

    void pcc(MatrixXd &PCC, vector<gwasinfo> &snpinfo, long cohortnum, double pmecs, bool zflag)
    {
        //pearson correlation with pairwise.complete.obs
        if(zflag) {
            LOGPRINTF("Using z score to calcualte the pearson correlation.\n");
        }
        else {
             LOGPRINTF("Using beta to calcualte the pearson correlation.\n");
        }
        double zmecs=qchisq(pmecs,1);
        long snpnum = snpinfo.size();
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
                    double sei=snpinfo[k].se[i];
                    double betai=snpinfo[k].beta[i];
                    double sej=snpinfo[k].se[j];
                    double betaj=snpinfo[k].beta[j];
                    if(abs(sei+9)>1e-6 && abs(sej+9)>1e-6) {
                        double zi=betai/sei;
                        double zj=betaj/sej;
                        double zi2=zi*zi;
                        double zj2=zj*zj;
                        if(zi2 < zmecs && zj2 < zmecs)
                        {
                            if(zflag) {
                                beta1.push_back(zi);
                                beta2.push_back(zj);
                            } else {
                                beta1.push_back(betai);
                                beta2.push_back(betaj);
                            }
                        }
                    }
                }
                if(beta1.size()<1) {
                    LOGPRINTF("WARNING: %ld SNPs in common between cohort %i and cohort %d (cohort number starts form 0).\n",beta1.size(),i,j);
                    LOGPRINTF("The correlation value of cohort %d and cohort %d would be imputed with the mean of all the correlation values excluding the diagnoal.\n",i,j);
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
            LOGPRINTF("WARNING: %ld cohort pairs didn't get enough common SNPs to calcualte the correlation.\n",pairsNoCor1.size());
            if(pairHasCorNUm==0) {
                LOGPRINTF("ERROR: Every pair of cohort has not enough common SNPs to calcualte the correlation.\n");
                TERMINATE();
            }
            double corMean=sumcor/pairHasCorNUm;
            LOGPRINTF("WARNING: These missing correlation values are imputed with the mean %f.\n",corMean);
            for(int i=0;i<pairsNoCor1.size();i++)
            {
                int p1=pairsNoCor1[i];
                int p2=pairsNoCor2[i];
                PCC(p1,p2)=PCC(p2,p1)=corMean;
            }
        }
        for( int i=0;i<cohortnum;i++) PCC(i,i)=1;
    }

    void mecs_gwas_fun(vector<gwasinfo> &snpinfo, long cohortnum,double pmecs,vector<int> &noninvertible, vector<int> &negativedeno, int mecs_mth,char* corMatFName, bool zflag)
    {

        MatrixXd Corr(cohortnum,cohortnum);
        long snpnum = snpinfo.size();
        if(corMatFName!=NULL)
        {
            LOGPRINTF("Reading the correlation matrix information from %s ...\n", corMatFName);
            FILE* corfile=NULL;
            uint32_t line_idx = 0;
            vector<string> strlist;
            if(fopen_checked(&corfile, corMatFName,"r")) {
                LOGPRINTF("ERROR: can't open file %s to read. \n", corMatFName);
                TERMINATE();
            }
            while(fgets(Tbuf, MAX_LINE_SIZE, corfile))
            {
                if(line_idx>=cohortnum)
                {
                    LOGPRINTF("ERROR: the file %s has more rows than the cohort number %ld.\n", corMatFName,cohortnum);
                    TERMINATE();
                }
                split_str(Tbuf,strlist,0);
                if(strlist.size()!=cohortnum)
                {
                    LOGPRINTF("ERROR: the number of elemnets of Line %u is not equal the cohort number %ld.\n", line_idx,cohortnum);
                    TERMINATE();
                }
                for(int j=0;j<strlist.size();j++) Corr(line_idx,j)=atof(strlist[j].c_str());
                line_idx++;
            }
            LOGPRINTF(" The correlation matrix has been included from the file %s.\n",corMatFName)
        }
        else
        {
            LOGPRINTF("Estimate the cohort correlation using the beta values of pairwised common SNPs.\n")
            if(mecs_mth==1) {
                LOGPRINTF("The correlation would be estimated.\n")
                est_cor(Corr,snpinfo,cohortnum,pmecs);
            } else {
                LOGPRINTF("The correlation would be computed using Pearson Correlaton Coefficient method.\n")
                LOGPRINTF("We exclude the significant common SNPs with a p-value threshold %e to estimate the correlation matrix.\n",pmecs);
                pcc(Corr,snpinfo,cohortnum,pmecs,zflag);
            }
        }

        LOGPRINTF("Saving the correlation matrix...\n");
        string filename=string(outfileName)+".cor.mat";
        FILE* tmpfile=fopen(filename.c_str(),"w");
        if(!tmpfile)
        {
            LOGPRINTF("error open file.\n");
            TERMINATE();
        }
        for(int t=0;t<Corr.rows();t++)
        {
            string str="";
            for(int k=0;k<Corr.cols();k++)
            {
                str +=atos(Corr(t,k)) + '\t';
            }
            str += '\n';
            fputs(str.c_str(),tmpfile);
        }
        fclose(tmpfile);
        LOGPRINTF("These correlation matrix is saved in the file %s.\n",filename.c_str());

        //int outScount=0;
        #pragma omp parallel for
        for(int j=0;j<snpnum;j++)
        {
            vector<double> ses, betas;
            vector<int> keep;
            MatrixXd Corr_work;
            int nmiss=0, miss=0;
            for(int k=0;k<cohortnum;k++)
            {
                double se=snpinfo[j].se[k];
                double beta=snpinfo[j].beta[k];
                snpinfo[j].se[k]=-9;
                if(abs(se+9)>1e-6){
                    ses.push_back(se);
                    betas.push_back(beta);
                    keep.push_back(k);
                    nmiss++;
                } else {
                    miss++;
                }
            }

            if(nmiss==1)
            {
                *snpinfo[j].beta=betas[0];
                *snpinfo[j].se=ses[0];
            }
            else if(nmiss>1)
            {
                if(prt_mid_rlt)
                {
                    string filena=string(outfileName)+"."+string(snpinfo[j].snprs)+".keep";
                    FILE* tmpfi=fopen(filena.c_str(),"w");
                    if(!tmpfi)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<keep.size();t++)
                    {

                        string str =atos(keep[t]) + '\n';
                        fputs(str.c_str(),tmpfi);
                    }
                    fclose(tmpfi);
                    LOGPRINTF("These keep IDs are saved in the file %s.\n",filena.c_str());
                }
                if(nmiss<cohortnum) subMatrix_symm(Corr_work, Corr, keep);
                else if(nmiss==cohortnum) Corr_work = Corr;
                else {
                    LOGPRINTF("Can't happen. I can guarantee!\n")
                }

                if(prt_mid_rlt){
                    LOGPRINTF("Saving the Cor matrix of SNP %s...\n",snpinfo[j].snprs);
                    string filenam=string(outfileName)+"."+string(snpinfo[j].snprs)+".cor.mat";
                    FILE* tmpfil=fopen(filenam.c_str(),"w");
                    if(!tmpfil)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<Corr_work.rows();t++)
                    {
                        string str="";
                        for(int k=0;k<Corr_work.cols();k++)
                        {
                            str +=atos(Corr_work(t,k)) + '\t';
                        }
                        str += '\n';
                        fputs(str.c_str(),tmpfil);
                    }
                    fclose(tmpfil);
                    LOGPRINTF("The cor matrix is saved in the file %s.\n",filename.c_str());
                }

                VectorXd sev(ses.size());
                for(int k=0;k<ses.size();k++) sev(k)=ses[k];
                MatrixXd W=sev*sev.transpose();
                W=W.array()*Corr_work.array();

                if(prt_mid_rlt)
                {
                    LOGPRINTF("Saving the S matrixn of SNP %s...\n",snpinfo[j].snprs);
                    string filename=string(outfileName)+"."+string(snpinfo[j].snprs)+".S.mat";
                    FILE* tmpfile=fopen(filename.c_str(),"w");
                    if(!tmpfile)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<W.rows();t++)
                    {
                        string str="";
                        for(int k=0;k<W.cols();k++)
                        {
                            str +=atos(W(t,k)) + '\t';
                        }
                        str += '\n';
                        fputs(str.c_str(),tmpfile);
                    }
                    fclose(tmpfile);
                    LOGPRINTF("These S matrix is saved in the file %s.\n",filename.c_str());
                }
                /*
                 //disable this becasue the determinant is toooooo stringent.
                double determinant = W.determinant();
                if(abs(determinant)<1e-6)
                {
                    LOGPRINTF("WARING: matrix S of %s SNP is not invertible because the determinant is %e.\n",snpinfo[j].snprs, determinant)
                    noninvertible.push_back(j);
                    *snpinfo[j].beta=-9;
                    *snpinfo[j].se=-9;
                }
                else
                {
                    bool determinant_zero=false;
                    inverse_V(W,determinant_zero);
                    if(determinant_zero)
                    {
                        LOGPRINTF("WARING: matrix S of %s SNP is not invertible because at least one of the eigenvalues is 0.\n",snpinfo[j].snprs)
                        noninvertible.push_back(j);
                        *snpinfo[j].beta=-9;
                        *snpinfo[j].se=-9;
                    }
                    else
                    {
                        double deno=W.sum();
                        if(deno<=0) {
                            LOGPRINTF("WARING: the sum of the inverse matrix of %s SNP is negative.\n",snpinfo[j].snprs)
                            negativedeno.push_back(j);
                            *snpinfo[j].beta=-9;
                            *snpinfo[j].se=-9;
                        }
                        else
                        {
                            VectorXd colsum=W.colwise().sum();
                            double numerator=0.0;
                            for(int k=0;k<betas.size();k++) numerator+=colsum(k)*betas[k];
                            *snpinfo[j].beta=numerator/deno;
                            *snpinfo[j].se=1/sqrt(deno);
                        }
                    }
                }
                 */
                bool determinant_zero=false;
                inverse_V(W,determinant_zero);
                if(prt_mid_rlt)
                {
                    LOGPRINTF("Saving the S inverse matrix of SNP %s...\n",snpinfo[j].snprs);
                    filename=string(outfileName)+"."+string(snpinfo[j].snprs)+".S.inv.mat";
                    tmpfile=fopen(filename.c_str(),"w");
                    if(!tmpfile)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<W.rows();t++)
                    {
                        string str="";
                        for(int k=0;k<W.cols();k++)
                        {
                            str +=atos(W(t,k)) + '\t';
                        }
                        str += '\n';
                        fputs(str.c_str(),tmpfile);
                    }
                    fclose(tmpfile);
                    LOGPRINTF("These inverse matrix is saved in the file %s.\n",filename.c_str());
                    LOGPRINTF("Saving the  betas...\n");
                    filename=string(outfileName)+"."+string(snpinfo[j].snprs)+".beta";
                    tmpfile=fopen(filename.c_str(),"w");
                    if(!tmpfile)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<betas.size();t++)
                    {

                        string str =atos(betas[t]) + '\n';
                        fputs(str.c_str(),tmpfile);
                    }
                    fclose(tmpfile);
                    LOGPRINTF("These betas is saved in the file %s.\n",filename.c_str());
                    filename=string(outfileName)+"."+string(snpinfo[j].snprs)+".se";
                    tmpfile=fopen(filename.c_str(),"w");
                    if(!tmpfile)
                    {
                        LOGPRINTF("error open file.\n");
                        TERMINATE();
                    }
                    for(int t=0;t<ses.size();t++)
                    {

                        string str =atos(ses[t]) + '\n';
                        fputs(str.c_str(),tmpfile);
                    }
                    fclose(tmpfile);
                    LOGPRINTF("These ses is saved in the file %s.\n",filename.c_str());

                }
                if(determinant_zero)
                {
                    LOGPRINTF("WARING: matrix S of %s SNP is not invertible because at least one of the eigenvalues is 0.\n",snpinfo[j].snprs)
                    noninvertible.push_back(j);
                    *snpinfo[j].beta=-9;
                    *snpinfo[j].se=-9;
                }
                else
                {
                    double deno=W.sum();
                    if(deno<=0) {
                        LOGPRINTF("WARING: the sum of the inverse matrix of SNP %s is negative (%f).\n",snpinfo[j].snprs,deno);
                        negativedeno.push_back(j);
                        *snpinfo[j].beta=-9;
                        *snpinfo[j].se=-9;
                    }
                    else
                    {
                        VectorXd colsum=W.colwise().sum();
                        double numerator=0.0;
                        for(int k=0;k<betas.size();k++) numerator+=colsum(k)*betas[k];
                        *snpinfo[j].beta=numerator/deno;
                        *snpinfo[j].se=1/sqrt(deno);
                    }
                }
            }
        }
    }

    // save to DENSE1 or SPARSE (SMR format)
    void meta(char* besdlistFileName, char* outFileName, int meta_mth, double pthresh, bool cis_flag, int cis_itvl,int nmecs,char* problstName, char* problst2exclde, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* snplstName, char* snplst2exclde,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb, bool smpwindFlag, char* probe2excld, char* snprs2excld, bool trans_meta)
    {
        if(meta_mth && !trans_meta) {
            cis_flag=true; // for later update. !!!!
            LOGPRINTF("NOTE: --mecs is applied. Only the information in the cis-region would be used.\n");
        }
        string analysisType="";
        if(meta_mth) analysisType="MeCS";
        else analysisType="Meta";
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;

        read_msglist(besdlistFileName, besds,"eQTL summary file names");
        if(besds.size()<=1) {
            LOGPRINTF("Less than 2 BESD files list in %s.\n",besdlistFileName);
            TERMINATE();
        }
        LOGPRINTF("%ld eQTL summary file names are included.\n",besds.size());

        LOGPRINTF("Checking the BESD format...\n");
        vector<int> format, smpsize;
        check_besds_format(besds, format, smpsize);
        int label=-1;
        for(int i=0;i<format.size();i++)
        {
            if(format[i]==SMR_DENSE_1 || format[i]==SMR_DENSE_3 || format[i]==OSCA_DENSE_1) {
                if(label==-1) {
                    label=0;
                } else if(label==1) {
                    label=2;
                    break;
                }

            } else if (format[i]==SMR_SPARSE_3 || format[i]==SMR_SPARSE_3F || format[i]==OSCA_SPARSE_1) {
                if(label==-1) {
                    label=1;
                } else if(label==0) {
                    label=2;
                    break;
                }
            } else {
                LOGPRINTF("Some BESDs are of old sparse format. please use SMR to re-make it.\n");
                TERMINATE();
            }
        }
        if(label==0) {
            LOGPRINTF("All the BESDs are of dense format.\n");
            LOGPRINTF("To save the memory, the strategy is processing (reading information, meta analysis, writing result) probe by probe.\n");
        } else if(label==1) {
            LOGPRINTF("All the BESDs are of sparse format.\n");
        } else {
            LOGPRINTF("Some of the BESDs are of sparse format and others are of dense format.\n");
            LOGPRINTF("The result would be saved in dense format.\n");
        }
        if(meta_mth) label=1; // for later update. !!!!

        combine_epi(probeinfo, besds,nprb, problstName, problst2exclde, genelistName,  chr, prbchr, prbname, fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2excld);
        combine_esi(snpinfo, besds,nsnp,snplstName, snplst2exclde,  chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb,  smpwindFlag,snprs2excld);

        if(probeinfo.size()==0)
        {
            LOGPRINTF("ERROR: No probe to be included!\n");
            TERMINATE();
        }
        smr_probeinfo* epiptr=&probeinfo[0];
        qsort(epiptr,probeinfo.size(),sizeof(smr_probeinfo),comp_epi);
        if(snpinfo.size()==0)
        {
            LOGPRINTF("ERROR: No SNP to be included!\n");
            TERMINATE();
        }
        smr_snpinfo* esiptr=&snpinfo[0];
        qsort(esiptr,snpinfo.size(),sizeof(smr_snpinfo),comp_esi);

        long besdNum=besds.size();
        long metaPrbNum=probeinfo.size();
        long metaSNPnum=snpinfo.size();

        LOGPRINTF("\nGenerating epi file...\n");
        FILE* efile=NULL;
        string epiName=string(outFileName)+".epi";
        if(fopen_checked(&efile, epiName.c_str(),"w")) TERMINATE();
        for(int i=0;i<probeinfo.size();i++)
        {
            string chrstr;
            if(probeinfo[i].probechr==23) chrstr="X";
            else if(probeinfo[i].probechr==24) chrstr="Y";
            else chrstr=atosm(probeinfo[i].probechr);

            string str=chrstr+'\t'+probeinfo[i].probeId+'\t'+atos(0)+'\t'+atosm(probeinfo[i].bp)+'\t'+probeinfo[i].genename+'\t'+(probeinfo[i].orien=='*'?"NA":atos(probeinfo[i].orien))+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", epiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld probes have been saved in the file %s .\n", probeinfo.size(), epiName.c_str());


        LOGPRINTF("\nGenerating esi file...\n");
        string esiName=string(outFileName)+".esi";
        if(fopen_checked(&efile, esiName.c_str(),"w")) TERMINATE();
        for(int i=0;i<snpinfo.size();i++)
        {
            string chrstr;
            if(snpinfo[i].snpchr==23) chrstr="X";
            else if(snpinfo[i].snpchr==24) chrstr="Y";
            else chrstr=atosm(snpinfo[i].snpchr);
            string str=chrstr+'\t'+snpinfo[i].snprs+'\t'+atos(0)+'\t'+atosm(snpinfo[i].bp)+'\t'+snpinfo[i].a1+'\t'+snpinfo[i].a2+'\t'+(abs(snpinfo[i].freq+9)>1e-6?atos(snpinfo[i].freq):"NA")+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", esiName.c_str());
                TERMINATE();
            }
        }
        fclose(efile);
        LOGPRINTF("%ld SNPs have been saved in the file %s .\n", snpinfo.size(), esiName.c_str());

        lookup.resize(besdNum);
        for(int i=0;i<besdNum;i++) lookup[i].resize(nsnp[i]);
        for(int i=0;i<besdNum;i++)
            for(int j=0;j<lookup[i].size();j++)
                lookup[i][j]=-9;
        for(int i=0;i<metaSNPnum;i++)
        {
            for(int j=0;j<besdNum;j++)
            {
                int tmpval=snpinfo[i].rstr[j];
                if(tmpval>=0) {
                    if(tmpval >=nsnp[j])
                    {
                        LOGPRINTF("ERROR: bug found in snpinfo. Please report.\n");
                        TERMINATE();
                    }
                    lookup[j][tmpval]=i;
                }
            }
        }

        LOGPRINTF("\nPerforming %s analysis (the result will be saved in BESD format)....\n",analysisType.c_str());
        if(meta_mth){
            LOGPRINTF("Estimate the cohort correlation using the beta values of pair-wised common SNPs.\n");
            LOGPRINTF("We exclude the significant common SNPs with a p-value threshold %6.2e to estimate the correlation matrix.\n",pthresh);
        }
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            if(fopen_checked(&fptrs[i],besdFileName.c_str(),"rb")) {
                LOGPRINTF("ERROR: in opening file %s .\n", besdFileName.c_str());
                TERMINATE();
            }
        }
        vector<uint64_t> cols;
        vector<uint32_t> rowids;
        vector<float> val;
        string besdName=string(outFileName)+".besd";
        if(fopen_checked(&efile, besdName.c_str(),"wb")) TERMINATE();
        uint32_t filetype=SMR_DENSE_3;
        if(label==1)
        {
            filetype=SMR_SPARSE_3;
            cols.resize((metaPrbNum<<1)+1);
            cols[0]=0;
        }
        int ttl_smp_size=0;
        for(int i=0;i<smpsize.size();i++)
        {
            if(smpsize[i]!=-9)
            {
                ttl_smp_size+=smpsize[i];
            } else {
                ttl_smp_size=-9;
                LOGPRINTF("WARNING: At lease one BESD file doesn't contain the sample size information. Please ues --add-n to update the sample size after the BESD file is generated.\n");
                break;
            }
        }
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=ttl_smp_size;
        ten_ints[2]=(int)snpinfo.size();
        ten_ints[3]=(int)probeinfo.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        if (fwrite_checked(&ten_ints[0],RESERVEDUNITS*sizeof(int), efile))
        {
            LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
            TERMINATE();
        }


        float* buffer_beta=(float *)malloc(sizeof(float)*besdNum*metaSNPnum);
        if (buffer_beta == NULL) {
            LOGPRINTF("Memory buffer for beta values error.\n");
            TERMINATE();
        } //probe major

        float* buffer_se=(float *)malloc(sizeof(float)*besdNum*metaSNPnum);
        if (buffer_se == NULL) {
            LOGPRINTF("Memory buffer for SEs error.\n");
            TERMINATE();
        }
        for(int i=0;i<besdNum*metaSNPnum;i++) buffer_se[i]=-9;
        vector<string> noninvtb_prbs;
        vector<string> nega_prbs;
        vector<string> snpdeficent;
        double cr=0.0;
        for(int i=0;i<metaPrbNum;i++)
        {
            progress(i, cr, (int)metaPrbNum);
            //LOGPRINTF("processing with probe %s...\n",probeinfo[i].probeId);
            vector<float> betases;
            vector<uint32_t> row_ids;
            int probebp=probeinfo[i].bp;
            int probechr=probeinfo[i].probechr;
            long cohortnum=0;
            for(int j=0;j<besds.size();j++)
            {
                int pid=probeinfo[i].ptr[j];
                betases.clear();
                row_ids.clear();
                if(pid>=0)
                {
                    if(format[j]==SMR_SPARSE_3F || format[j]==SMR_SPARSE_3 || format[j]==OSCA_SPARSE_1 )
                    {
                        extract_prb_sparse(fptrs[j], (uint64_t)pid, nprb[j],row_ids, betases);
                        long num=row_ids.size();
                        if(num==0) {
                            //LOGPRINTF("WARNING: empty probe %s found in the sparse file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                            continue;
                        }
                        long alignnum=0;
                        if(meta_mth && cis_flag) {
                            //LOGPRINTF("Extract the cis-region of probe %s in the sparse file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        } else {
                            //LOGPRINTF("Extract the eQTLs of probe %s in the sparse file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        }
                        for(int k=0;k<num;k++)
                        {
                            //align each cohort to the buffer. we don't have -9 of se from sparse.
                            int idx=lookup[j][row_ids[k]];
                            if(idx>=0)
                            {
                                int snpchr=snpinfo[idx].snpchr;
                                int snpbp=snpinfo[idx].bp;
                                if(meta_mth && cis_flag) {
                                    //extract cis if mecs
                                    int startpos=(probebp-cis_itvl*1000)>0?(probebp-cis_itvl*1000):0;
                                    int endpos=probebp+cis_itvl*1000;
                                    if(snpchr==probechr && snpbp >= startpos && snpbp <= endpos) {
                                        if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                        else buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                        buffer_se[cohortnum*metaSNPnum+idx]=betases[k+num];
                                        alignnum++;
                                    }
                                } else {
                                    if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                    else buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                    buffer_se[cohortnum*metaSNPnum+idx]=betases[k+num];
                                    alignnum++;
                                }
                            }
                        }
                        /*
                        if(alignnum==0) {
                            LOGPRINTF("No eQTLs extracted from probe %s in the sparse file %s for %s analysis.\n",probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                            continue;
                        } else {
                            LOGPRINTF("%ld eQTLs extracted from probe %s in the sparse file %s for %s analysis.\n",alignnum,probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                        }
                         */
                    }
                    else if(format[j]==SMR_DENSE_1 || format[j]==SMR_DENSE_3 ||  format[j]==OSCA_DENSE_1 )
                    {
                        extract_prb_dense(fptrs[j], (uint64_t)pid, nprb[j], nsnp[j], betases);
                       // LOGPRINTF("%llu eQTLs extracted from probe %s in the dense file %s.\n",nsnp[j],probeinfo[i].probeId,besds[j].c_str());
                        long alignnum=0;
                        if(meta_mth && cis_flag) {
                            LOGPRINTF("Extract the cis-region of probe %s in the file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        } else {
                            LOGPRINTF("Extract the eQTLs of probe %s in the file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        }
                        for(int k=0;k<nsnp[j];k++)
                        {
                            int idx=lookup[j][k];
                            if(idx>=0)
                            {
                                int snpchr=snpinfo[idx].snpchr;
                                int snpbp=snpinfo[idx].bp;
                                if(meta_mth && cis_flag)
                                {
                                    //extract cis if mecs
                                    int startpos=(probebp-cis_itvl*1000)>0?(probebp-cis_itvl*1000):0;
                                    int endpos=probebp+cis_itvl*1000;
                                    if(snpchr==probechr && snpbp >= startpos && snpbp <= endpos) {
                                        if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                        else buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                        buffer_se[cohortnum*metaSNPnum+idx]=betases[k+nsnp[j]];
                                        alignnum++;
                                    }
                                } else
                                {
                                    if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                    else  buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                    buffer_se[cohortnum*metaSNPnum+idx]=betases[k+nsnp[j]];
                                    alignnum++;
                                }
                            }
                        }
                        /*
                        if(alignnum==0) {
                            LOGPRINTF("no eQTLs extracted from probe %s in the dense file %s for %s analysis.\n",probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                            continue;
                        } else {
                            LOGPRINTF("%ld eQTLs extracted from probe %s in the dense file %s for %s analysis.\n",alignnum,probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                        }
                         */
                    }
                    cohortnum++;
                }
                else {
                    //LOGPRINTF("probe %s is not in the file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                }
            }

            //****test**
//            string filename=string(probeinfo[i].probeId)+".txt";
//            FILE* tmpfile=fopen(filename.c_str(),"w");
//            if(!tmpfile)
//            {
//                printf("error open file.\n");
//                exit(EXIT_FAILURE);
//            }
//            for(int t=0;t<metaSNPnum;t++)
//            {
//                string str=snpinfo[t].snprs;
//                for(int tt=0;tt<cohortnum;tt++)
//                {
//                    if( abs(buffer_se[tt*metaSNPnum+t]+9)<1e-6) {
//                        str+="\tNA\tNA";
//                    } else {
//                        str+='\t'+atos(buffer_beta[tt*metaSNPnum+t])+'\t'+atos(buffer_se[tt*metaSNPnum+t]);
//                    }
//
//                }
//                str+='\n';
//                fputs(str.c_str(),tmpfile);
//            }
//
//            fclose(tmpfile);
            //***end test**


            if(cohortnum==0) {
                //LOGPRINTF("No information of probe %s is included from any cohort for %s analysis.\n\n",probeinfo[i].probeId,analysisType.c_str());
            } else if(cohortnum==1) {
                //LOGPRINTF("The information of probe %s is included from %ld / %ld cohorts for %s analysis.\n",probeinfo[i].probeId,cohortnum,besds.size(),analysisType.c_str());
                //LOGPRINTF("The information of the probe %s would be saved in the result.\n\n",probeinfo[i].probeId);

            } else if(cohortnum>1) {
                //LOGPRINTF("The information of probe %s is included from %ld / %ld cohorts for %s analysis.\n",probeinfo[i].probeId,cohortnum,besds.size(),analysisType.c_str());
                if(meta_mth){
                    //LOGPRINTF("Performing %s analysis of probe %s...\n",analysisType.c_str(), probeinfo[i].probeId);
                    vector<int> noninvertible, negativedeno;
                    bool mecsflag = mecs_per_prob( buffer_beta, buffer_se, metaSNPnum, cohortnum, pthresh,noninvertible,negativedeno, nmecs);
                    if(!mecsflag) snpdeficent.push_back(probeinfo[i].probeId);
                    else {
                        if(noninvertible.size()>0) {
                            //LOGPRINTF("%ld SNPs of probe %s have non-invertible S matrix.\n",noninvertible.size(), probeinfo[i].probeId);
                            noninvtb_prbs.push_back(probeinfo[i].probeId);
                        }
                        if(negativedeno.size()>0) {
                            //LOGPRINTF("%ld SNPs of probe %s have negative 1'inv(S)1 .\n",negativedeno.size(), probeinfo[i].probeId);
                            nega_prbs.push_back(probeinfo[i].probeId);
                        }
                    }

                    //LOGPRINTF("end of %s analysis of probe %s.\n\n",analysisType.c_str(), probeinfo[i].probeId);
                }
                else {
                    //LOGPRINTF("Performing %s analysis of probe %s...\n",analysisType.c_str(), probeinfo[i].probeId);
                    meta_per_prob(buffer_beta,buffer_se, metaSNPnum, cohortnum);
                    //LOGPRINTF("End of %s analysis of probe %s.\n\n",analysisType.c_str(), probeinfo[i].probeId);
                }
            }

            if(label!=1){
                if (fwrite_checked(buffer_beta,metaSNPnum*sizeof(float), efile))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                if (fwrite_checked(buffer_se, metaSNPnum*sizeof(float), efile))
                {
                    LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                    TERMINATE();
                }
                for(int k=0;k<metaSNPnum;k++) buffer_se[k]=-9;
            } else {
                vector<float> tmpse;
                vector<uint32_t> tmprid;
                for(int k=0;k<metaSNPnum;k++)
                {
                    if(abs(buffer_se[k]+9)>1e-6)
                    {
                        val.push_back(buffer_beta[k]);
                        rowids.push_back(k);
                        tmpse.push_back(buffer_se[k]);
                        tmprid.push_back(k);
                        buffer_se[k]=-9;
                    }
                }
                for(int k=0;k<tmpse.size();k++)
                {
                    val.push_back(tmpse[k]);
                    rowids.push_back(tmprid[k]);
                }
                uint64_t real_num=tmpse.size();
                cols[(i<<1)+1]=real_num+cols[i<<1];
                cols[i+1<<1]=(real_num<<1)+cols[i<<1];
            }
        }

        if(label==1)
        {
            uint64_t valNum=val.size();
            if(fwrite_checked(&valNum,sizeof(uint64_t), efile)) {
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if(fwrite_checked(&cols[0],sizeof(uint64_t)*cols.size(), efile)){
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if(fwrite_checked(&rowids[0],sizeof(uint32_t)*rowids.size(), efile)){
                LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                TERMINATE();
            }
            if(fwrite_checked(&val[0],sizeof(float)*val.size(), efile)){
                   LOGPRINTF("ERROR: in writing binary file %s .\n", besdName.c_str());
                   TERMINATE();
               }
        }
        if(snpdeficent.size()>0)
        {
            string filename=string(outFileName)+".nmecs"+atos(nmecs)+".probe.list";
            FILE* tmpfile=fopen(filename.c_str(),"w");
            if(!tmpfile)
            {
                printf("ERROR: open file %s.\n",filename.c_str());
                exit(EXIT_FAILURE);
            }
            for(int t=0;t<snpdeficent.size();t++)
            {
                string str=snpdeficent[t]+'\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
        }
        if(noninvtb_prbs.size()>0)
        {
            //LOGPRINTF("\nWARNING: %ld probes have at least one eQTL whose S is non-invertible.\n",noninvtb_prbs.size());
            string filename=string(outFileName)+".non-invertible.probe.list";
            FILE* tmpfile=fopen(filename.c_str(),"w");
            if(!tmpfile)
            {
                LOGPRINTF("error open file.\n");
                TERMINATE();
            }
            for(int t=0;t<noninvtb_prbs.size();t++)
            {
                string str=noninvtb_prbs[t]+'\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
            //LOGPRINTF("These probes are saved in file %s.\n",filename.c_str());
        }

        if(nega_prbs.size()>0)
        {
            //LOGPRINTF("\nWARNING: %ld probes have at least one eQTL whose 1'inv(S)1 is negative .\n",nega_prbs.size());
            //LOGPRINTF("WARNING: That means we can't get SE by doing square-root.\n");
            //LOGPRINTF("WARNING: In term of such case we set effect size as 0 and SE as missing (-9).\n");
            string filename=string(outFileName)+".negative.probe.list";
            FILE* tmpfile=fopen(filename.c_str(),"w");
            if(!tmpfile)
            {
                printf("error open file.\n");
                exit(EXIT_FAILURE);
            }
            for(int t=0;t<nega_prbs.size();t++)
            {
                string str=nega_prbs[t]+'\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
            //LOGPRINTF("These probes are saved in file %s.\n",filename.c_str());
        }


        for(int i=0;i<probeinfo.size();i++)
        {
            if(probeinfo[i].genename) free2(&probeinfo[i].genename);
            if(probeinfo[i].probeId) free2(&probeinfo[i].probeId);
            if(probeinfo[i].ptr) free2(&probeinfo[i].ptr);
            if(probeinfo[i].bfilepath) free2(&probeinfo[i].bfilepath);
            if(probeinfo[i].esdpath) free2(&probeinfo[i].esdpath);
        }
        for(int i=0;i<snpinfo.size();i++)
        {
            if(snpinfo[i].a1) free2(&snpinfo[i].a1);
            if(snpinfo[i].a2) free2(&snpinfo[i].a2);
            if(snpinfo[i].snprs) free2(&snpinfo[i].snprs);
            if(snpinfo[i].rstr) free2(&snpinfo[i].rstr);
            if(snpinfo[i].revs) free2(&snpinfo[i].revs);
        }
        for(int i=0;i<besds.size();i++)
        {
            fclose(fptrs[i]);
        }
        fclose(efile);
        free(fptrs);
        free(buffer_se);
        free(buffer_beta);

        LOGPRINTF("\nThe eQTL infomation of %ld probes and %ld SNPs has been in binary file %s.\n",metaPrbNum,metaSNPnum,besdName.c_str());
    }

    void meta_gwas(char* gwaslistFileName, char* ewaslistFileName, char* outFileName, int meta_mth, double pthresh, int mecs_mth, char* corMatFName, char* snplstName, char* problstName, bool zflag,bool out_comm_flag, bool all_comm_flag)
    {
        if(corMatFName!=NULL)
        {
            LOGPRINTF("--cor-mat is active. it would disbale --mecs-mth, --pcc-z and --pmecs.\n")
        }
        if(snplstName!=NULL)
        {
            if(corMatFName==NULL)
            {
                LOGPRINTF("ERROR: --extract-snp should be with --cor-mat.\n");
                TERMINATE();
            }
        }
        if(gwaslistFileName!=NULL && ewaslistFileName!=NULL)
        {
            LOGPRINTF("ERROR: please specify only one statistic summary file list.\n");
            TERMINATE();
        }
        string analysisType="";
        if(meta_mth) analysisType="MeCS";
        else analysisType="Meta";
        vector<string> gwass;
        vector<gwasinfo> snpinfo;

        if(gwaslistFileName!=NULL) read_msglist(gwaslistFileName, gwass,"GWAS summary file names");
        else read_msglist(ewaslistFileName, gwass,"EWAS summary file names");
        if(gwass.size()<=1) {
            LOGPRINTF("Less than 2 GWAS/EWAS summary files list in %s.\n",gwaslistFileName);
            TERMINATE();
        }
        LOGPRINTF("%ld GWAS/EWAS summary file names are included.\n",gwass.size());

        if(gwaslistFileName!=NULL) combine_gwas(snpinfo, gwass,snplstName);
        else combine_ewas(snpinfo, gwass,problstName);

        long cohortNum=gwass.size();
        long metaSNPnum=snpinfo.size();
        if(metaSNPnum==0) {
            LOGPRINTF("ERROR: No SNPs/probes included.\n");
            TERMINATE();
        }

        LOGPRINTF("\nPerforming %s analysis and save the result in text file....\n",analysisType.c_str());

        FILE* rltfile=NULL;
        if(fopen_checked(&rltfile, outFileName,"w"))
        {
            LOGPRINTF("ERROR: open result %s file error.\n",outFileName);
            TERMINATE();
        }
        string str="SNP\tA1\tA2\tfreq\tb\tse\tp\tn\n";
        if(ewaslistFileName!=NULL) str="Chr\tProbe\tBP\tGene\tOrien\tb\tse\tp\tn\n";
        if(fputs_checked(str.c_str(),rltfile))
        {

            LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
            TERMINATE();
        }

        vector<int> itscid,allcomm;
        if(meta_mth){
            vector<int> noninvertible, negativedeno;
            LOGPRINTF("NOTE: MeCS could be sensitive and give a biased result when QC is not pre-performed. e.g. extravagant se such as 1e15...\n");
            mecs_gwas_fun( snpinfo, cohortNum, pthresh,noninvertible,negativedeno,mecs_mth,corMatFName,zflag);
            if(noninvertible.size()>0) {
                LOGPRINTF("%ld SNPs have non-invertible S matrix.\n",noninvertible.size());
            }
            if(negativedeno.size()>0) {
                LOGPRINTF("%ld SNPs have negative 1'inv(S)1 .\n",negativedeno.size());
            }
            LOGPRINTF("End of %s analysis.\n",analysisType.c_str());
        }
        else {
            meta_gwas_fun(snpinfo, cohortNum,itscid, allcomm, out_comm_flag, all_comm_flag);
            LOGPRINTF("End of %s analysis.\n",analysisType.c_str());
        }
        LOGPRINTF("Saving %s results...\n",analysisType.c_str());
        if(out_comm_flag)
        {
            metaSNPnum=itscid.size();
            for(int ii=0;ii<itscid.size();ii++)
            {
                int i=itscid[ii];
                double pval= -9;
                if(abs(*snpinfo[i].se+9)>1e-6) {
                    double z= *snpinfo[i].beta / *snpinfo[i].se;
                    pval=pchisq(z*z, 1);
                }

                if(gwaslistFileName!=NULL) str = string(snpinfo[i].snprs) + '\t' + snpinfo[i].a1 +'\t' + snpinfo[i].a2 +'\t' + atosm(snpinfo[i].freq) +'\t' + atosm(*snpinfo[i].beta) +'\t' + atosm(*snpinfo[i].se) +'\t' + ((pval==-9)?"NA":dtos(pval)) +'\t' + atosm((snpinfo[i].estn)) + '\n';
                else
                {
                    string chr=atosm(round(snpinfo[i].freq));
                    if(snpinfo[i].freq==-9) chr="NA";
                    else if(snpinfo[i].freq==23) chr="X";
                    else if(snpinfo[i].freq==24) chr="Y";
                    str = chr + '\t' + snpinfo[i].snprs +'\t' + atosm(snpinfo[i].bp) +'\t'+ snpinfo[i].a1 +'\t' + snpinfo[i].a2 +'\t' + atosm(*snpinfo[i].beta) +'\t' + atosm(*snpinfo[i].se) +'\t' + ((pval==-9)?"NA":dtos(pval)) + '\n';
                }
                if(fputs_checked(str.c_str(),rltfile))
                {

                    LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                    TERMINATE();
                }
            }
        }
        else if(all_comm_flag)
        {
            metaSNPnum=allcomm.size();
            for(int ii=0;ii<allcomm.size();ii++)
            {
                int i=allcomm[ii];
                double pval= -9;
                if(abs(*snpinfo[i].se+9)>1e-6) {
                    double z= *snpinfo[i].beta / *snpinfo[i].se;
                    pval=pchisq(z*z, 1);
                }

                if(gwaslistFileName!=NULL) str = string(snpinfo[i].snprs) + '\t' + snpinfo[i].a1 +'\t' + snpinfo[i].a2 +'\t' + atosm(snpinfo[i].freq) +'\t' + atosm(*snpinfo[i].beta) +'\t' + atosm(*snpinfo[i].se) +'\t' + ((pval==-9)?"NA":dtos(pval)) +'\t' + atosm((snpinfo[i].estn)) + '\n';
                else
                {
                    string chr=atosm(round(snpinfo[i].freq));
                    if(snpinfo[i].freq==-9) chr="NA";
                    else if(snpinfo[i].freq==23) chr="X";
                    else if(snpinfo[i].freq==24) chr="Y";
                    str = chr + '\t' + snpinfo[i].snprs +'\t' + atosm(snpinfo[i].bp) +'\t'+ snpinfo[i].a1 +'\t' + snpinfo[i].a2 +'\t' + atosm(*snpinfo[i].beta) +'\t' + atosm(*snpinfo[i].se) +'\t' + ((pval==-9)?"NA":dtos(pval)) + '\n';
                }
                if(fputs_checked(str.c_str(),rltfile))
                {

                    LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                    TERMINATE();
                }
            }
        }
        else
        {
            for(int i=0;i<metaSNPnum;i++)
            {
                double pval= -9;
                if(abs(*snpinfo[i].se+9)>1e-6) {
                    double z= *snpinfo[i].beta / *snpinfo[i].se;
                    pval=pchisq(z*z, 1);
                }

                if(gwaslistFileName!=NULL) str = string(snpinfo[i].snprs) + '\t' + snpinfo[i].a1 +'\t' + snpinfo[i].a2 +'\t' + atosm(snpinfo[i].freq) +'\t' + atosm(*snpinfo[i].beta) +'\t' + atosm(*snpinfo[i].se) +'\t' + ((pval==-9)?"NA":dtos(pval)) +'\t' + atosm((snpinfo[i].estn)) + '\n';
                else
                {
                    string chr=atosm(round(snpinfo[i].freq));
                    if(snpinfo[i].freq==-9) chr="NA";
                    else if(snpinfo[i].freq==23) chr="X";
                    else if(snpinfo[i].freq==24) chr="Y";
                    str = chr + '\t' + snpinfo[i].snprs +'\t' + atosm(snpinfo[i].bp) +'\t'+ snpinfo[i].a1 +'\t' + snpinfo[i].a2 +'\t' + atosm(*snpinfo[i].beta) +'\t' + atosm(*snpinfo[i].se) +'\t' + ((pval==-9)?"NA":dtos(pval)) + '\n';
                }
                if(fputs_checked(str.c_str(),rltfile))
                {

                    LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                    TERMINATE();
                }
            }
        }

        for(int i=0;i<snpinfo.size();i++)
        {
            if(snpinfo[i].a1) free2(&snpinfo[i].a1);
            if(snpinfo[i].a2) free2(&snpinfo[i].a2);
            if(snpinfo[i].snprs) free2(&snpinfo[i].snprs);
            if(snpinfo[i].beta) free2(&snpinfo[i].beta);
            if(snpinfo[i].se) free2(&snpinfo[i].se);
        }

        fclose(rltfile);

        LOGPRINTF("The GWAS summary statistics of %ld SNPs have been saved in text file %s.\n",metaSNPnum,outFileName);
    }

    void make_besd(char* outFileName,char* beqtlFileName, char* snplstName, char* snplst2exclde, char* problstName,char* problst2exclde, char* genelistName, double plookup, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl, char* probe2exclde, char* snprs2exclde,bool save_dense_flag, bool tosmrflag,bool besd_shrink_flag, bool stdprb, char* frqFName,char* varFName)
    {
        string logstr;
        if(beqtlFileName == NULL) {
            LOGPRINTF("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
            TERMINATE();
        }
        int flag4chr=0;
        if(chr > 0) flag4chr++;
        if(prbchr > 0 || snpchr > 0) flag4chr++;
        if(flag4chr==2)
        {
            chr=-9;
            LOGPRINTF("WARNING: --chr is not surpposed to use together with --probe-chr or --snp-chr. --chr will be disabled.\n");
        }

        eqtlInfo eqtlinfo;
        LOGPRINTF("\nReading eQTL summary data...\n");
        char inputname[FNAMESIZE];
        memcpy(inputname,beqtlFileName,strlen(beqtlFileName)+1);
        char* suffix=inputname+strlen(beqtlFileName);
        memcpy(suffix,".epi",5);
        read_smr_epifile(&eqtlinfo, inputname);
        smr_epi_man(&eqtlinfo, problstName, problst2exclde, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2exclde);
        memcpy(suffix,".esi",5);
        read_smr_esifile(&eqtlinfo, inputname);
        smr_esi_man(&eqtlinfo, snplstName, snplst2exclde,chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname,snprs2exclde);
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
        if(stdprb)
        {
            if(frqFName==NULL && varFName==NULL)
            {
                LOGPRINTF("Error: please input the fequency data or the variance data for the standardisation by the flag --freq-file or --var-file.\n");
                TERMINATE();
            }
        }

        if(stdprb)
        {
            if(varFName!=NULL) read_probevarfile(&eqtlinfo,varFName);
            else {
                LOGPRINTF("Error: please input the variance data for the standardisation by the flag --var-file .\n");
                LOGPRINTF("WARNING: The software tool doesn't support --freq-file for the time being.\n");
                TERMINATE();
            }
        }

        memcpy(suffix,".besd",6);
        vector<int> headers;
        get_BesdHeaders(inputname, headers);
        int indicator = headers[0];
        if(indicator==SMR_DENSE_1 || indicator==SMR_SPARSE_3F) eqtlinfo._sampleNum=-9;
        else eqtlinfo._sampleNum=headers[1];

        //if(!save_dense_flag) {
        //    LOGPRINTF("Release soon!\n");
        //    TERMINATE();

       // }
       // else
       // {
       //     LOGPRINTF("Save dense flag is actived.\n");
            if(indicator==SMR_SPARSE_3 || indicator==SMR_SPARSE_3F || indicator==OSCA_SPARSE_1)
            {
                LOGPRINTF("file %s is in sparse BESD format.\n",inputname);
                read_smr_besdfile(&eqtlinfo, inputname);
                if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
                {
                    LOGPRINTF("No data included from %s under current condition.\n",beqtlFileName);
                    TERMINATE();
                }
                if(stdprb)
                {
                    LOGPRINTF("Standarding probes...\n"); // I don't think there is missing value in sparse format.
                    for(int i=0;i<eqtlinfo._include.size();i++)
                    {
                        int proid=eqtlinfo._include[i];
                        double prbvar_sqrt=sqrt(eqtlinfo._epi_var[proid]);
                        uint64_t pos=eqtlinfo._cols[proid];
                        uint64_t pos1=eqtlinfo._cols[proid+1];
                        uint64_t num=(pos1-pos)>>1;
                        for(int j=0;j<num;j++)
                        {
                            double beta=eqtlinfo._val[pos+j];
                            double se=eqtlinfo._val[pos+j+num];
                            eqtlinfo._val[pos+j]=beta/prbvar_sqrt;
                            eqtlinfo._val[pos+j+num]=se/prbvar_sqrt;
                        }
                    }
                }
                if(besd_shrink_flag)
                {
                    vector<string> nullprbs;
                    vector<string> nullsnps;
                     LOGPRINTF("Checking null probes and null SNPs ...\n");
                    get_shrink_null(&eqtlinfo,nullprbs,nullsnps);
                    if(nullprbs.size()>0)
                    {
                        LOGPRINTF("%ld probes have no SNPs .\n",nullprbs.size());
                        string filename=string(outFileName)+".null.probe.list";
                        FILE* tmpfile=fopen(filename.c_str(),"w");
                        if(!tmpfile)
                        {
                            LOGPRINTF("error open file.\n");
                            TERMINATE();
                        }
                        for(int t=0;t<nullprbs.size();t++)
                        {
                            string str=nullprbs[t]+'\n';
                            fputs(str.c_str(),tmpfile);
                        }
                        fclose(tmpfile);
                        LOGPRINTF("These probes are saved in file %s.\n",filename.c_str());
                        update_map_rm(nullprbs, eqtlinfo._probe_name_map, eqtlinfo._include);
                    } else {
                        LOGPRINTF("All probes have SNPs.\n");
                    }
                    if(nullsnps.size()>0)
                    {
                        LOGPRINTF("%ld SNPs have no value across all the probes .\n",nullsnps.size());
                        string filename=string(outFileName)+".null.snp.list";
                        FILE* tmpfile=fopen(filename.c_str(),"w");
                        if(!tmpfile)
                        {
                            LOGPRINTF("error open file.\n");
                            TERMINATE();
                        }
                        for(int t=0;t<nullsnps.size();t++)
                        {
                            string str=nullsnps[t]+'\n';
                            fputs(str.c_str(),tmpfile);
                        }
                        fclose(tmpfile);
                        LOGPRINTF("These SNPs are saved in file %s.\n",filename.c_str());
                        update_map_rm(nullsnps, eqtlinfo._snp_name_map, eqtlinfo._esi_include);
                    } else {
                         LOGPRINTF("No null SNP found.\n");
                    }
                }
                write_smr_esi(outFileName, &eqtlinfo);
                write_smr_epi(outFileName, &eqtlinfo);
                write_s2s_besd(outFileName, &eqtlinfo,tosmrflag);

            }
            else if(indicator==SMR_DENSE_1 || indicator == SMR_DENSE_3 || indicator == OSCA_DENSE_1)
            {
                LOGPRINTF("file %s is in dense BESD format.\n",inputname);
                if(besd_shrink_flag){
                    LOGPRINTF("WARNING: --besd-shrink doesn't work on dense BESD file.\n");
                }
                if(!stdprb && indicator==SMR_DENSE_3 && eqtlinfo._snpNum==eqtlinfo._esi_include.size() && eqtlinfo._probNum==eqtlinfo._include.size())
                {
                    LOGPRINTF("The output would be identical with the input. You can make a copy with termial command.\n");
                    TERMINATE();
                }
                write_smr_esi(outFileName, &eqtlinfo);
                write_smr_epi(outFileName, &eqtlinfo);
                write_d2d_besd(outFileName, &eqtlinfo,inputname, stdprb);

            }
        //}

        //LOGPRINTF("Extracted results of %ld SNPs have been saved in the file %s.\n",out_esi_id.size(),outFileName);

    }

    void gc_ewas(char* outFileName, char* ewasFileName)
    {
        gwasData gtmp;
        read_ewas_data(&gtmp, ewasFileName);
        vector<double> z2;
        for(int i=0;i<gtmp.byz.size();i++)
        {
            if(gtmp.seyz[i]>0)
            {
                double z=gtmp.byz[i]/gtmp.seyz[i];
                z2.push_back(z*z);
            }
        }
        double lambda=median(z2)/0.455;
        LOGPRINTF("The genomic inflation factor is %7.2f.\n", lambda);
        double lambdasq=sqrt(lambda);
        for(int i=0;i<gtmp.byz.size();i++)
        {
            if(gtmp.seyz[i]>0)
            {
                 gtmp.byz[i]/=lambdasq;
                double chi=gtmp.byz[i]/gtmp.seyz[i];
                chi*=chi;
                gtmp.pvalue[i]=pchisq(chi,1);
            }
            else
            {
                gtmp.seyz[i]=-9;
                gtmp.pvalue[i]=-9;
            }

        }
        FILE* rltfile=NULL;
        if(fopen_checked(&rltfile, outFileName,"w"))
        {
            LOGPRINTF("ERROR: open result %s file error.\n",outFileName);
            TERMINATE();
        }
        string str="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\n";
        if(fputs_checked(str.c_str(),rltfile))
        {

            LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
            TERMINATE();
        }
        for(int i=0;i<gtmp.byz.size();i++)
        {
            string chr=atosm(gtmp.freq[i]);
            if(gtmp.freq[i]==-9) chr="NA";
            else if(gtmp.freq[i]==23) chr="X";
            else if(gtmp.freq[i]==24) chr="Y";
            str = chr + '\t' + gtmp.snpName[i] +'\t' + atosm(gtmp.snpBp[i]) +'\t'+ gtmp.allele_1[i] +'\t' + gtmp.allele_2[i] +'\t' + ((gtmp.seyz[i]<0)?"NA":atosm(gtmp.byz[i])) +'\t' + ((gtmp.seyz[i]<0)?"NA":atosm(gtmp.seyz[i]))  +'\t' + ((gtmp.seyz[i]<0)?"NA":dtos(gtmp.pvalue[i]))  + '\n';
            if(fputs_checked(str.c_str(),rltfile))
            {

                LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
                TERMINATE();
            }
        }
        fclose(rltfile);
        LOGPRINTF("GC adjusted results have been saved in file %s.\n", outFileName);
    }

}

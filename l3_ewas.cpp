//
//  l3_ewas.cpp
//  osc
//
//  Created by Futao Zhang on 16/04/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#include "l3_ewas.hpp"
namespace EFILE
{
     void backward_elimn(eInfo* einfo, vector<int> &slct, MatrixXd &prob_profile, double p_thresh, vector<int> &dump )
    {
        int nindi=(int)einfo->_eii_include.size();
        while(slct.size()>2)
        {
            int X_c=1+(int)slct.size();;
            MatrixXd X(nindi, X_c);
            X.block(0, 0, nindi, 1) = MatrixXd::Ones(nindi, 1);
            for(int i=0;i<slct.size();i++) X.col(i+1)=prob_profile.col(slct[i]);
            MatrixXd XtX_i;
            XtX_i=X.transpose()*X;
            bool determinant_zero=false;
            inverse_V(XtX_i, determinant_zero);
            
            VectorXd y(einfo->_eii_include.size());
            #pragma omp parallel for
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                y(j)=einfo->_eii_pheno[einfo->_eii_include[j]];
            }
            VectorXd b_hat=XtX_i*X.transpose()*y;
            VectorXd residual=(y-X*b_hat);
            residual=residual.array()*residual.array();
            double sy=sqrt(residual.sum()/(y.size()-X.cols()));
            VectorXd se=sy*XtX_i.diagonal().array().sqrt();
            VectorXd t=b_hat.array()/se.array();
            vector<double> p(t.size()-1);
            for(int j=1;j<t.size();j++) p[j-1]=pchisq( t(j)*t(j),1 );
            int m = (int)(max_element(p.begin(), p.end()) - p.begin());
            if(p[m]>p_thresh)
            {
                dump.push_back(slct[m]);
                slct.erase(slct.begin()+m);
            } else break;
        }
        
    }
    bool forward_slct(eInfo* einfo, vector<int> &slct,vector<int> &remain, MatrixXd &prob_profile, double p_thresh, vector<int> &dump )
    {
        vector<double> rlt;
        int nindi=(int)einfo->_eii_include.size(), X_c=1+(int)slct.size();;
        MatrixXd X(nindi, X_c);
        X.block(0, 0, nindi, 1) = MatrixXd::Ones(nindi, 1);
        for(int i=0;i<slct.size();i++) X.col(i+1)=prob_profile.col(slct[i]);
        
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
        
        VectorXd y(einfo->_eii_include.size());
        #pragma omp parallel for
        for(int j=0; j<einfo->_eii_include.size(); j++)
        {
            y(j)=einfo->_eii_pheno[einfo->_eii_include[j]];
        }
        VectorXd b_hat=XtXiXt*y;
        VectorXd yresi=(y-X*b_hat);// if no covariate yresi=y-mean(y)
        vector<double> condp;
        for(int i=0;i<remain.size();i++)
        {
            VectorXd x=prob_profile.col(remain[i]);
            VectorXd b_hat=XtXiXt*x;
            VectorXd residual=(x-X*b_hat);
            vector<double> rst;
            adjusted_reg(yresi, residual, rst,X_c-1);
            condp.push_back(rst[2]);
        }
        int m = (int)(min_element(condp.begin(), condp.end()) - condp.begin());
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
    void MakeX(eInfo* einfo,vector<string> &prbs, MatrixXd &X)
    {
        map<string, int>::iterator iter;
        X.resize(einfo->_eii_include.size(),prbs.size());
        for(int i=0;i<prbs.size();i++)
        {
            iter=einfo->_epi_map.find(prbs[i]);
            if(iter!=einfo->_epi_map.end())
            {
                vector<int> MISS;
                double mu= 0.0, sum=0.0;
                int idx=iter->second;
                for(int j=0; j<einfo->_eii_include.size(); j++)
                {
                    double val=einfo->_val[idx*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9)
                    {
                        X(j,i)=val;
                        sum+=val;
                    } else {
                        MISS.push_back(j);
                    }
                    
                    double ratio=1.0*MISS.size()/einfo->_eii_include.size();
                    if(ratio>0.1)
                    {
                        printf("WARNING: %f %% of the values in probe %s are missing.\n",ratio*100, einfo->_epi_prb[idx].c_str());
                    }
                    if(MISS.size()>0)
                    {
                        if(einfo->_eii_include.size() > MISS.size()) mu=sum/(einfo->_eii_include.size()-MISS.size());
                        for(int k=0;k<MISS.size();k++)
                        {
                            X(MISS[k],i)=mu;
                        }
                    }
                }
            } else {
                LOGPRINTF("ERROR: probe %s can't be found.\n",prbs[i].c_str());
                TERMINATE();
            }
        }

    }
    void stepwise_slct(eInfo* einfo,vector<int> &sig, vector<int> &slctid, vector<int> &rmid, vector<ASSOCRLT> &assoc_rlts, double p_cutoff, bool updatesig)
    {
        vector<double> p_buf;
        map<string, int>::iterator iter;
        MatrixXd sig_profile(einfo->_eii_include.size(),sig.size());
        int p=1,m=0;
        for(int i=0;i<sig.size();i++)
        {
            if(p>assoc_rlts[sig[i]].PVAL)
            {
                p=assoc_rlts[sig[i]].PVAL;
                m=i;
            }
            iter=einfo->_epi_map.find(assoc_rlts[sig[i]].PROBE);
            if(iter!=einfo->_epi_map.end())
            {
                vector<int> MISS;
                double mu= 0.0, sum=0.0;
                int idx=iter->second;
                for(int j=0; j<einfo->_eii_include.size(); j++)
                {
                    double val=einfo->_val[idx*einfo->_eii_num+einfo->_eii_include[j]];
                    if(val<1e9)
                    {
                        sig_profile(j,i)=val;
                        sum+=val;
                    } else {
                        MISS.push_back(j);
                    }
                    
                    double ratio=1.0*MISS.size()/einfo->_eii_include.size();
                    if(ratio>0.1)
                    {
                        printf("WARNING: %f %% of the values in probe %s are missing.\n",ratio*100, einfo->_epi_prb[idx].c_str());
                    }
                    if(MISS.size()>0)
                    {
                        if(einfo->_eii_include.size() > MISS.size()) mu=sum/(einfo->_eii_include.size()-MISS.size());
                        for(int k=0;k<MISS.size();k++)
                        {
                            sig_profile(MISS[k],i)=mu;
                        }
                    }
                }
            }
            else
            {
                LOGPRINTF("ERROR: probe %s can't be found.\n",assoc_rlts[sig[i]].PROBE);
                TERMINATE();
            }
        }
        
        vector<int> slctidx, remain;
        if (assoc_rlts[sig[m]].PVAL >= p_cutoff) return;
        slctidx.push_back(m);
        for (int i = 1; i < sig.size(); i++) {
             remain.push_back(i);
        }
        while (!remain.empty()) {
            if (forward_slct(einfo,slctidx, remain, sig_profile,p_cutoff,rmid)) {
                backward_elimn(einfo,slctidx, sig_profile,p_cutoff,rmid);
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
            LOGPRINTF("%ld probes are selected.\n", slctidx.size());
        }
    }
    
    //backgound control
    void moment2(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_no_adj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_num, int expect_pcs, double r2thresh,int erm_alg)
    {
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        expect_wind>>=1;
        int nrandcomp =1;
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
        double secondthresh=1e-5;
        
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
            LOGPRINTF("ERROR: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("ERROR: please input the phenotype data by the option --pheno.\n");
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
        
        map<string, int>::iterator iter;
        
        int _n=(int)einfo._eii_include.size();
        if(_n<1)
        {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        VectorXd _y;
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=einfo._eii_pheno[einfo._eii_include[i]];
        }
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _Xin;
        int _Xin_c;
        _Xin_c=construct_X(&einfo, E_float, qE_float,_Xin); //_X has the dimension of [einfo._eii_include.size(),_X_c]
        MatrixXd _X_PC;
        MatrixXd evec;
        
        vector< vector<int> > mapids;
        vector< vector<string> > erm_prbs;
        vector<MatrixXd> AZERO, _A, A0;
        mapids.resize(2); // one for significant one for insignificant
        erm_prbs.resize(2);
        AZERO.resize(1);
        LOGPRINTF("Performing feature selection.\n");
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        testLinear(assoc_rlts, NULL, &einfo,COV_plus);
        //testQAssoc(assoc_rlts,NULL, &einfo);
        ASSOCRLT* sortptr=&assoc_rlts[0];
        qsort(sortptr,assoc_rlts.size(),sizeof(ASSOCRLT),comp_assoc);
        if(r2thresh<0 && bcthresh<0) bcthresh=0.05/assoc_rlts.size();
        else if(r2thresh > 0) {
            double chi=einfo._eii_include.size()* r2thresh / (1 - r2thresh) + 1;
            bcthresh = pchisq(chi, 1);
        }
        
        for(int k=(int)(assoc_rlts.size()-1);k>=0;k--)
        {
            if(assoc_rlts[k].PVAL<=bcthresh) {
                mapids[0].push_back(k); //sig
                erm_prbs[0].push_back(assoc_rlts[k].PROBE);
            }
            else {
                mapids[1].push_back(k);
                erm_prbs[1].push_back(assoc_rlts[k].PROBE);
            }
        }
        if( mapids[0].size()>0)
        {
            vector<int> slctids, rmsig;
            stepwise_slct(&einfo,mapids[0], slctids, rmsig, assoc_rlts,  bcthresh, true);
            erm_prbs[0].clear();
            for(int j=0;j<mapids[0].size();j++)
            {
                erm_prbs[0].push_back(assoc_rlts[mapids[0][j]].PROBE);
            }
            for(int j=0;j<rmsig.size();j++)
            {
                mapids[1].push_back(rmsig[j]);
                erm_prbs[1].push_back(assoc_rlts[rmsig[j]].PROBE);
            }
            string tmpfname=string(outFileName)+".AS.txt";
            
            FILE* tmpfile=fopen(tmpfname.c_str(),"w");
            if(!tmpfile)
            {
                LOGPRINTF("error open file %s.\n",tmpfname.c_str());
                TERMINATE();
            }
            for(int i=0;i< erm_prbs[0].size();i++) {
                string str=atos( erm_prbs[0][i]);
                str+='\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);

        }
        vector<int> include_o(einfo._epi_include);
        map<string, int> snp_name_map_o(einfo._epi_map);
        //random
        A0.resize(2);
        
            if(erm_prbs[1].size()>0)
            {
                update_map_kp(erm_prbs[1], einfo._epi_map, einfo._epi_include);
                make_erm(&einfo,erm_alg);
                (AZERO[0]).resize(_n, _n);
                (A0[0]).resize(_n, _n);
                #pragma omp parallel for
                for(int k=0; k<_n; k++)
                {
                    for(int l=0; l<=k; l++){
                        (AZERO[0])(l,k)=(AZERO[0])(k,l)=einfo._grm(k,l)*einfo._grm_N(k,l);
                        (A0[0])(l,k)=(A0[0])(k,l)=einfo._grm(k,l);
                    }
                }
                einfo._epi_include=include_o;
                einfo._epi_map=snp_name_map_o;
            }
        
        A0[1]=MatrixXd::Identity(_n, _n);
        
        LOGPRINTF("\nPerforming the MOMENT analysis ...\n");
        LOGPRINTF("For each probe, the analysis will exclude probes in %d Kb region of centered at the probe to be tested.\n",expect_wind);
        int exwind=expect_wind*1000;
        
        string filename=string(outFileName)+".mlma";
        FILE* ofile = fopen(filename.c_str(), "w");
        if (!(ofile)) {
            LOGPRINTF("ERROR: open error %s\n", filename.c_str());
            TERMINATE();
        }
        string outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\n";
        if(fputs_checked(outstr.c_str(),ofile))
        {
            LOGPRINTF("ERROR: error in writing file %s .\n", filename.c_str());
            TERMINATE();
        }
        
        long lrcount=0, wcount=0;
        int outid=0;
        long num2exp=expect_num>mapids[outid].size()?mapids[outid].size():expect_num;
        double cr=0;
        for(int i=0;i<num2exp;i++)
        {
            double desti=1.0*wcount/(einfo._epi_include.size()-1);
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            
            
            int id=mapids[outid][i];
            int targetChr=assoc_rlts[id].CHR;
            int targetBP=assoc_rlts[id].BP;
            string targetPrb=assoc_rlts[id].PROBE;
            string gene=assoc_rlts[id].GENE;
            char oren=assoc_rlts[id].OREN;
            
            iter=einfo._epi_map.find(assoc_rlts[id].PROBE);
            VectorXd x_buf(einfo._eii_include.size());
            if(iter!=einfo._epi_map.end())
            {
                int curidx=iter->second;
                for(int j=0; j<einfo._eii_include.size(); j++)
                    x_buf(j)=einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
            } else {
                LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                TERMINATE();
            }
            
            vector<string> fxprb,rndexprb;
            int norms=1,excount=0;
            for(int jj=0;jj<mapids.size();jj++)
            {
                for(int j=0;j<mapids[jj].size();j++)
                {
                    int tid=mapids[jj][j];
                    if(targetChr == assoc_rlts[tid].CHR && abs(targetBP-assoc_rlts[tid].BP)<=exwind) {
                        if(jj) rndexprb.push_back(assoc_rlts[tid].PROBE);
                        excount++;
                    } else {
                        if(!jj) fxprb.push_back(assoc_rlts[tid].PROBE);
                    }
                }
                
            }
            
            einfo._r_indx.clear();
            for(int j=0; j < norms + 1; j++) einfo._r_indx.push_back(j);
            _A.resize(einfo._r_indx.size());
            int aid=0;
            
                long N=mapids[1].size()-rndexprb.size();
                if(N)
                {
                    if(rndexprb.size()>0)
                    {
                        update_map_kp(rndexprb, einfo._epi_map, einfo._epi_include);
                        make_erm(&einfo,erm_alg);
                        (_A[aid]).resize(_n, _n);
                        #pragma omp parallel for
                        for(int k=0; k<_n; k++)
                        {
                            for(int l=0; l<=k; l++) (_A[aid])(l,k)=(_A[aid])(k,l)=(AZERO[aid](k,l)-einfo._grm(k,l)*einfo._grm_N(k,l))/N;
                        }
                        einfo._epi_include=include_o;
                        einfo._epi_map=snp_name_map_o;
                    } else {
                        _A[aid]=AZERO[aid]/N;
                    }
                }
            
            _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
            
            
            MatrixXd _X(_Xin.rows(), _Xin.cols()+fxprb.size());
            if(fxprb.size()>0)
            {
                MatrixXd Xprb2;
                MakeX(&einfo,fxprb,Xprb2);
                _X << _Xin, Xprb2;
            } else {
                _X  = _Xin;
            }
            int _X_c=(int)_X.cols();
            
            // names of variance component
            einfo._var_name.clear();
            einfo._hsq_name.clear();
            for (int j = 0; j < norms; j++) {
                stringstream strstrm;
                if (norms == 1) strstrm << "";
                else strstrm << j + 1;
                einfo._var_name.push_back("V(O" + strstrm.str() + ")");
                einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
            }
            einfo._var_name.push_back("V(e)");
            
            einfo._within_family=within_family;
            if(within_family) detect_family(&einfo, _A);
            MatrixXd _Vi;
            remlstatus=0; //reset reml status
            if(norms) reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
            double beta, se, pval;
            if(norms && (remlstatus==0 || remlstatus==-5 || (remlstatus==-3 && force_mlm)))
            {
                einfo._P.resize(0,0);
                _A.clear();
                VectorXd y_buf=_y;
                if(!mlma_no_adj_covar)
                {
                    y_buf=_y.array()-(_X*einfo._b).array();
                    if(_X_c>1) getResidual(x_buf, _X);
                }
                if(mlma_no_adj_covar) mlma_cal_stat_covar(y_buf, x_buf, _Vi, _X, beta, se, pval);
                else mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
            }
            else
            {
                if(remlstatus==-1) {
                    LOGPRINTF("\nThe matrix is not invertible, ");
                } else if(remlstatus==-2) {
                    LOGPRINTF("\nMore than half of the variance components are constrained, ");
                } else if(remlstatus==-3) {
                    LOGPRINTF("\nVariance component going to 0 or 1, ");
                }else if(remlstatus==-4) {
                    LOGPRINTF("\nLog-likelihood not converged, ");
                }
                if(i==lrcount && lrcount>=10)
                {
                    expect_num=i;
                    break;
                }
                LOGPRINTF("\nPerforming association analysis with %d PCs ...\n",expect_pcs);
                VectorXd y_buf=_y;
                if(_X_PC.cols()==0)
                {
                    make_erm( &einfo,erm_alg); //use the whole
                    SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
                    evec = (eigensolver.eigenvectors());
                    VectorXd eval = eigensolver.eigenvalues();
                    //if(einfo._epi_include.size()<=einfo._eii_include.size())
                    if(einfo._epi_include.size()<=expect_pcs)
                    {
                        LOGPRINTF("WARNING: probe number is too few for PCs...\n");
                        expect_pcs = (int)einfo._epi_include.size();
                    }
                    MatrixXd _PCs = evec.block(0,evec.cols()-expect_pcs, evec.rows(),expect_pcs);
                    _X_PC.resize(_X.rows(), _X.cols()+expect_pcs);
                    _X_PC << _X, _PCs;
                }
                
                LR(y_buf, x_buf, _X_PC, beta, se, pval);
                lrcount++;
            }
            string chrstr;
            if(targetChr==23) chrstr="X";
            else if(targetChr==24) chrstr="Y";
            else chrstr=atosm(targetChr);
            outstr = chrstr + '\t' + targetPrb + '\t' + atosm(targetBP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
            if(fputs_checked(outstr.c_str(),ofile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                TERMINATE();
            }
            wcount++;
        }
        vector<int> leftids;
        if(mapids[outid].size()>expect_num)
            for(int i=expect_num;i<mapids[outid].size();i++)
                leftids.push_back(mapids[outid][i]);
        for(int outid=1;outid<mapids.size();outid++)
        {
            for(int i=0;i<mapids[outid].size();i++)
                leftids.push_back(mapids[outid][i]);
        }
        
        einfo._r_indx.clear();
        for(int j=0; j < 1 + 1; j++) einfo._r_indx.push_back(j);
        // names of variance component
        einfo._var_name.clear();
        einfo._hsq_name.clear();
        for (int j = 0; j < 1; j++) {
            stringstream strstrm;
            if (1 == 1) strstrm << "";
            else strstrm << j + 1;
            einfo._var_name.push_back("V(O" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        einfo._within_family=within_family;
        if(within_family) detect_family(&einfo, _A);
        MatrixXd _Vi;
        //probes as fixed
        MatrixXd _X(_Xin.rows(), _Xin.cols()+erm_prbs[0].size());
        if(erm_prbs[0].size()>0)
        {
            MatrixXd Xprb;
            MakeX(&einfo,erm_prbs[0],Xprb);
            _X << _Xin, Xprb;
        } else {
            _X  = _Xin;
        }
        int _X_c=(int)_X.cols();
        remlstatus=0; //reset reml status
        reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,A0,_Vi,outFileName);
        
        if( (remlstatus==0 || remlstatus==-5 || (remlstatus==-3 && force_mlm)))
        {
            double beta, se, pval;
            VectorXd y_buf=_y;
            if(!mlma_no_adj_covar)
                y_buf=_y.array()-(_X*einfo._b).array();
            for(int i=0;i<leftids.size();i++)
            {
                double desti=1.0*wcount/(einfo._epi_include.size()-1);
                if(desti>=cr)
                {
                    printf("%3.0f%%\r", 100.0*desti);
                    fflush(stdout);
                    if(cr==0) cr+=0.05;
                    else if(cr==0.05) cr+=0.2;
                    else if(cr==0.25) cr+=0.5;
                    else cr+=0.25;
                }
                
                int id=leftids[i];
                int targetChr=assoc_rlts[id].CHR;
                int targetBP=assoc_rlts[id].BP;
                string targetPrb=assoc_rlts[id].PROBE;
                string gene=assoc_rlts[id].GENE;
                char oren=assoc_rlts[id].OREN;
                
                iter=einfo._epi_map.find(assoc_rlts[id].PROBE);
                VectorXd x_buf(einfo._eii_include.size());
                if(iter!=einfo._epi_map.end())
                {
                    int curidx=iter->second;
                    for(int j=0; j<einfo._eii_include.size(); j++)
                        x_buf(j)=einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                } else {
                    LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                    TERMINATE();
                }
                
                if(!mlma_no_adj_covar)
                    if(_X_c>1) getResidual(x_buf, _X);
                
                if(mlma_no_adj_covar) mlma_cal_stat_covar(y_buf, x_buf, _Vi, _X, beta, se, pval);
                else mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                
                string chrstr;
                if(targetChr==23) chrstr="X";
                else if(targetChr==24) chrstr="Y";
                else chrstr=atosm(targetChr);
                outstr = chrstr + '\t' + targetPrb + '\t' + atosm(targetBP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
                if(fputs_checked(outstr.c_str(),ofile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                    TERMINATE();
                }
                wcount++;
            }
            
        }
        else
        {
            erm_prbs[outid].clear();
            for(int i=0;i<leftids.size();i++)
            {
                int id=leftids[i];
                erm_prbs[outid].push_back(assoc_rlts[id].PROBE);
            }
            update_map_kp(erm_prbs[outid], einfo._epi_map, einfo._epi_include);
            long maxpc2test = (evec.cols()>128)?64:(evec.cols()>>1);
            long pc2test = (maxpc2test>2)?2:maxpc2test;
            long regionL = pc2test, regionR=pc2test;
            vector<ASSOCRLT> out_rlts, otemp;
            double outlambda=1e6;
            while(pc2test <= maxpc2test)
            {
                LOGPRINTF("\nPerforming association analysis with %ld PCs ...\n",pc2test);
                double upperlambda=1+lambda_wind;
                double lowerlambda=(1-lambda_wind)>0?(1-lambda_wind):0;
                MatrixXd _X_PC = evec.block(0,evec.cols()-pc2test, evec.rows(),pc2test);
                if(fastlinear) testLinear_fast(otemp,NULL, &einfo,_X_PC);
                else testLinear(otemp,NULL, &einfo,_X_PC);
                double lambda=get_lambda(otemp);
                LOGPRINTF("The lambda for %ld PCs is %f.\n",pc2test, lambda);
                double curdis=abs(1-lambda);
                double predis= abs(1-outlambda);
                if( curdis < predis)
                {
                    out_rlts = otemp;
                    outlambda = lambda;
                    LOGPRINTF("The lambda for %ld PCs is most close to 1.\n",pc2test);
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
                    if(pc2test==regionL || pc2test>maxpc2test) break;
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
                    if(pc2test==regionR) break;
                    
                }
                else break;
            }
            for(int k=0;k<out_rlts.size();k++)
            {
                string chrstr;
                if(out_rlts[k].CHR==23) chrstr="X";
                else if(out_rlts[k].CHR==24) chrstr="Y";
                else chrstr=atosm(out_rlts[k].CHR);
                outstr = chrstr + '\t' + out_rlts[k].PROBE + '\t' + atosm(out_rlts[k].BP) + '\t' + atos(out_rlts[k].GENE) + '\t' + atos(out_rlts[k].OREN) + '\t' + atos(out_rlts[k].BETA) + '\t' + atos(out_rlts[k].SE) + '\t' + dtos(out_rlts[k].PVAL)  +'\n';
                if(fputs_checked(outstr.c_str(),ofile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                    TERMINATE();
                }
                wcount++;
            }
            free_assoclist(otemp);
        }
        LOGPRINTF("Results of %ld probes have been saved in file %s.\n",wcount,filename.c_str());
        fclose(ofile);
        free_assoclist(assoc_rlts);
    }

    void moment(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_preadj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_num, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_flag, bool approximate_stepwise,int erm_alg, double swthresh)
    {
        //bcthresh: default as bonfferoni
        //r2thresh: default as -9, disabled.
        //percent: default as -9, disabled.
        //slctmtd: feature selection with linear or mlma. default as linear.
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        expect_wind>>=1;
        
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
        double secondthresh=1e-5;
        
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
            LOGPRINTF("ERROR: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("ERROR: please input the phenotype data by the option --pheno.\n");
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
        if(!moment_eligibility_ck(&einfo))
        {
            LOGPRINTF("ERROR: the .opi file contains NA chromosome or NA probe position. Plsease annotate it using flag --update-opi.\n");
            TERMINATE();
        }
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo); // eii and epi are updated in it.
        stdprobe(&einfo);
        if(einfo._eType==0 && expect_wind==50) expect_wind=100;
        map<string, int>::iterator iter;
    
        int _n=(int)einfo._eii_include.size();
        if(_n<1)
        {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        //LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        VectorXd _y;
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
        
         /*
        //LOGPRINTF("Performing principal component analysis ...\n");
        make_erm( &einfo,erm_alg); //use the whole
        SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
        MatrixXd evec = (eigensolver.eigenvectors());
        VectorXd eval = eigensolver.eigenvalues();
        if(einfo._epi_include.size()<=expect_pcs)
        {
            LOGPRINTF("WARNING: probe number is too few for PCs...\n");
            expect_pcs = (int)einfo._epi_include.size();
        }
         
        MatrixXd _PCs = evec.block(0,evec.cols()-expect_pcs, evec.rows(),expect_pcs);
        MatrixXd _X_PC(_X.rows(), _X.cols()+expect_pcs);
        _X_PC << _X, _PCs;
        */
        
        vector< vector<int> > mapids;
        vector< vector<string> > erm_prbs;
        vector<MatrixXd> AZERO, _A, A0;
        mapids.resize(nrandcomp); // one for significant one for insignificant
        erm_prbs.resize(nrandcomp);
        AZERO.resize(nrandcomp);
        
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        if(slctmtd)
        {
            bool flagtmp=testMOA(assoc_rlts, &einfo, erm_alg);
            if(!flagtmp) testLinear(assoc_rlts, NULL, &einfo,COV_plus);
        }
        else
        {
            testLinear(assoc_rlts, NULL, &einfo,COV_plus);
        }
        
        //testQAssoc(assoc_rlts,NULL, &einfo);
        ASSOCRLT* sortptr=&assoc_rlts[0];
        qsort(sortptr,assoc_rlts.size(),sizeof(ASSOCRLT),comp_assoc);
        
            if(nrandcomp==2)
            {
                if(percent>0)
                {
                    int number=ceil(assoc_rlts.size()*percent);
                    for(int k=(int)(assoc_rlts.size()-1);k>=(int)(assoc_rlts.size()-number);k--)
                    {
                        mapids[0].push_back(k); //sig
                        erm_prbs[0].push_back(assoc_rlts[k].PROBE);
                    }
                    for(int k=(int)(assoc_rlts.size()-number-1);k>=0;k--)
                    {
                        mapids[1].push_back(k);
                        erm_prbs[1].push_back(assoc_rlts[k].PROBE);
                    }
                    LOGPRINTF("%ld probes incldued in the first component and %ld in the second component.\n",mapids[0].size(),mapids[1].size());
                }
                else
                {
                    if(r2thresh<0 && bcthresh<0) bcthresh=0.05/assoc_rlts.size();
                    else if(r2thresh > 0) {
                        double chi=einfo._eii_include.size()* r2thresh / (1 - r2thresh) + 1;
                        bcthresh = pchisq(chi, 1);
                    }
                    for(int k=(int)(assoc_rlts.size()-1);k>=0;k--)
                    {
                        if(assoc_rlts[k].PVAL<=bcthresh) {
                            mapids[0].push_back(k); //sig
                            erm_prbs[0].push_back(assoc_rlts[k].PROBE);
                        }
                        else {
                            mapids[1].push_back(k);
                            erm_prbs[1].push_back(assoc_rlts[k].PROBE);
                        }
                    }
                    LOGPRINTF("%ld probes incldued in the first component and %ld in the second component.\n",mapids[0].size(),mapids[1].size());
                }
                
            }
            else if(nrandcomp==3)
            {
                // for test only
                double bonthresh=0.05/assoc_rlts.size();
                
                for(int k=(int)(assoc_rlts.size()-1);k>=0;k--)
                {
                    if(assoc_rlts[k].PVAL<=bonthresh)
                    {
                        mapids[0].push_back(k);
                        erm_prbs[0].push_back(assoc_rlts[k].PROBE);
                    }
                    else if(assoc_rlts[k].PVAL<=secondthresh)
                    {
                        mapids[1].push_back(k);
                        erm_prbs[1].push_back(assoc_rlts[k].PROBE);
                    }
                    else
                    {
                        mapids[2].push_back(k);
                        erm_prbs[2].push_back(assoc_rlts[k].PROBE);
                    }
                }
                LOGPRINTF("%ld probes included in the 1st random component.\n",mapids[0].size());
                LOGPRINTF("%ld probes included in the 2nd random component.\n",mapids[1].size());
                LOGPRINTF("%ld probes included in the 3rd random component.\n",mapids[2].size());
            }
        
        if( approximate_stepwise && mapids[0].size()>0)
        {
            // stepwise for the components
            vector<int> slctids, rmsig;
            if(swthresh<0) swthresh = bcthresh;
            stepwise_slct(&einfo,mapids[0], slctids, rmsig, assoc_rlts,  swthresh, true);
            erm_prbs[0].clear();
            for(int j=0;j<mapids[0].size();j++)
            {
                erm_prbs[0].push_back(assoc_rlts[mapids[0][j]].PROBE);
            }
            for(int j=0;j<rmsig.size();j++)
            {
                mapids[1].push_back(rmsig[j]);
                erm_prbs[1].push_back(assoc_rlts[rmsig[j]].PROBE);
            }
            string tmpfname=string(outFileName)+".AS.txt";
            
            FILE* tmpfile=fopen(tmpfname.c_str(),"w");
            if(!tmpfile)
            {
                LOGPRINTF("error open file %s.\n",tmpfname.c_str());
                TERMINATE();
            }
            for(int i=0;i< erm_prbs[0].size();i++) {
                string str=atos( erm_prbs[0][i]);
                str+='\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
            
        }
        
        vector<int> include_o(einfo._epi_include);
        map<string, int> snp_name_map_o(einfo._epi_map);
        int rawORMc=0;
        for(int i=0;i<mapids.size();i++) if(mapids[i].size()) rawORMc++;
        A0.resize(rawORMc+1);
        rawORMc=0;
        for(int j=0;j<erm_prbs.size();j++)
        {
            if(erm_prbs[j].size()>0)
            {
                update_map_kp(erm_prbs[j], einfo._epi_map, einfo._epi_include);
                make_erm(&einfo,erm_alg);
                (AZERO[j]).resize(_n, _n);
                (A0[rawORMc]).resize(_n, _n);
                #pragma omp parallel for
                for(int k=0; k<_n; k++)
                {
                    for(int l=0; l<=k; l++){
                        (AZERO[j])(l,k)=(AZERO[j])(k,l)=einfo._grm(k,l)*einfo._grm_N(k,l);
                         (A0[rawORMc])(l,k)=(A0[rawORMc])(k,l)=einfo._grm(k,l);
                    }
                }
                einfo._epi_include=include_o;
                einfo._epi_map=snp_name_map_o;
                rawORMc++;
            }
        }
        A0[rawORMc]=MatrixXd::Identity(_n, _n);
        
        LOGPRINTF("\nPerforming the MOMENT analysis ...\n");
        LOGPRINTF("For each probe, the analysis will exclude probes in %d Kb region of centered at the probe to be tested.\n",expect_wind);
        int exwind=expect_wind*1000;
       
        string filename=string(outFileName)+".mlma";
        FILE* ofile = fopen(filename.c_str(), "w");
        if (!(ofile)) {
            LOGPRINTF("ERROR: open error %s\n", filename.c_str());
            TERMINATE();
        }
        string outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\n";
        if(fputs_checked(outstr.c_str(),ofile))
        {
            LOGPRINTF("ERROR: error in writing file %s .\n", filename.c_str());
            TERMINATE();
        }
        
        long lrcount=0, wcount=0;
        int outid=0;
        /*
        if(approximate_flag)
        {
            // bonfferoni for componets and setwise for excluding the targets
                vector<int> slcted, rmsig;
                if(r2thresh<0 && bcthresh<0) bcthresh=0.05/assoc_rlts.size();
                else if(r2thresh > 0) {
                    double chi=einfo._eii_include.size()* r2thresh / (1 - r2thresh) + 1;
                    bcthresh = pchisq(chi, 1);
                }
            if(swthresh<0) swthresh = bcthresh;
                stepwise_slct(&einfo,mapids[0],slcted, rmsig, assoc_rlts,  swthresh,false);
                for(int j=0; j<slcted.size();j++) mapids[0][j]=slcted[j];
                for(int j=0; j<rmsig.size();j++) mapids[0][j+slcted.size()]=rmsig[j];
                expect_num=(int)slcted.size();
        }
         */
        if(expect_num<0) expect_num = (int)mapids[outid].size();
        long num2exp=expect_num>mapids[outid].size()?mapids[outid].size():expect_num;
        
        MatrixXd XX(_X.rows(), _X.cols()+1);
        XX.block(0,0,_X.rows(),_X.cols())=_X;
        int x_idx=_X_c;
        
        double cr=0;
        for(int i=0;i<num2exp;i++)
        {
            double desti=1.0*wcount/(einfo._epi_include.size()-1);
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }

                
                int id=mapids[outid][i];
                int targetChr=assoc_rlts[id].CHR;
                int targetBP=assoc_rlts[id].BP;
                string targetPrb=assoc_rlts[id].PROBE;
                string gene=assoc_rlts[id].GENE;
                char oren=assoc_rlts[id].OREN;
                MatrixXd X=XX;
                iter=einfo._epi_map.find(assoc_rlts[id].PROBE);
                VectorXd x_buf(einfo._eii_include.size());
                if(iter!=einfo._epi_map.end())
                {
                    int curidx=iter->second;
                    for(int j=0; j<einfo._eii_include.size(); j++)
                        x_buf(j)=einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                } else {
                    LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                    TERMINATE();
                }
            X.col(x_idx)=x_buf;
            for(int jj=0;jj<nrandcomp;jj++) erm_prbs[jj].clear();
            int norms=0,excount=0;
            for(int jj=0;jj<nrandcomp;jj++)
            {
                for(int j=0;j<mapids[jj].size();j++)
                {
                    int tid=mapids[jj][j];
                    if(targetChr == assoc_rlts[tid].CHR && abs(targetBP-assoc_rlts[tid].BP)<=exwind) {
                        erm_prbs[jj].push_back(assoc_rlts[tid].PROBE);
                        excount++;
                    }
                }
            }
            //LOGPRINTF("%d probe(s) excluded from calculating the ORM(s) for the probe %s.\n", excount,assoc_rlts[id].PROBE);
            for(int jj=0;jj<nrandcomp;jj++) if(mapids[jj].size()-erm_prbs[jj].size()>0) norms++;
            
                einfo._r_indx.clear();
                for(int j=0; j < norms + 1; j++) einfo._r_indx.push_back(j);
                _A.resize(einfo._r_indx.size());
                int aid=0;
                for(int j=0; j < erm_prbs.size(); j++)
                {
                    long N=mapids[j].size()-erm_prbs[j].size();
                    if(N)
                    {
                        if(erm_prbs[j].size()>0)
                        {
                            update_map_kp(erm_prbs[j], einfo._epi_map, einfo._epi_include);
                            make_erm(&einfo,erm_alg);
                            (_A[aid]).resize(_n, _n);
                            #pragma omp parallel for
                            for(int k=0; k<_n; k++)
                            {
                                for(int l=0; l<=k; l++) (_A[aid])(l,k)=(_A[aid])(k,l)=(AZERO[j](k,l)-einfo._grm(k,l)*einfo._grm_N(k,l))/N;
                            }
                            einfo._epi_include=include_o;
                            einfo._epi_map=snp_name_map_o;
                        } else {
                            _A[aid]=AZERO[j]/N;
                        }
                        aid++;
                    }
                }
                _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
                
                // names of variance component
                einfo._var_name.clear();
                einfo._hsq_name.clear();
                for (int j = 0; j < norms; j++) {
                    stringstream strstrm;
                    if (norms == 1) strstrm << "";
                    else strstrm << j + 1;
                    einfo._var_name.push_back("V(O" + strstrm.str() + ")");
                    einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
                }
                einfo._var_name.push_back("V(e)");
                
                einfo._within_family=within_family;
                if(within_family) detect_family(&einfo, _A);
                MatrixXd _Vi;
                remlstatus=0; //reset reml status
                //if(norms) reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c+1, X,_y,_A,_Vi,outFileName); //exact model for test
                if(norms) reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName); // null model
                double beta=0, se=-9, pval=1;
                if(norms && (remlstatus==0  || remlstatus==-5 || remlstatus==-3 ))
                {
                    einfo._P.resize(0,0);
                    _A.clear();
                    VectorXd y_buf=_y;
                    if(mlma_preadj_covar)
                    {
                        y_buf=_y.array()-(_X*einfo._b).array();
                        if(_X_c>1) getResidual(x_buf, _X);
                    }
                    if(mlma_preadj_covar) mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                    else mlma_cal_stat_covar(y_buf, x_buf, _Vi, _X, beta, se, pval);
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
                string chrstr;
                if(targetChr==23) chrstr="X";
                else if(targetChr==24) chrstr="Y";
                else chrstr=atosm(targetChr);
                outstr = chrstr + '\t' + targetPrb + '\t' + atosm(targetBP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
                if(fputs_checked(outstr.c_str(),ofile))
                {
                    LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                    TERMINATE();
                }
                wcount++;
        }
        vector<int> leftids;
        if(mapids[outid].size()>expect_num)
            for(int i=expect_num;i<mapids[outid].size();i++)
                leftids.push_back(mapids[outid][i]);
        for(int outid=1;outid<mapids.size();outid++)
        {
            for(int i=0;i<mapids[outid].size();i++)
                leftids.push_back(mapids[outid][i]);
        }
        
        einfo._r_indx.clear();
        for(int j=0; j < rawORMc + 1; j++) einfo._r_indx.push_back(j);
        // names of variance component
        einfo._var_name.clear();
        einfo._hsq_name.clear();
        for (int j = 0; j < rawORMc; j++) {
            stringstream strstrm;
            if (rawORMc == 1) strstrm << "";
            else strstrm << j + 1;
            einfo._var_name.push_back("V(O" + strstrm.str() + ")");
            einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        einfo._var_name.push_back("V(e)");
        
        einfo._within_family=within_family;
        if(within_family) detect_family(&einfo, _A);
        MatrixXd _Vi;
        remlstatus=0; //reset reml status
        remloasi = false;
        reml_priors_var.clear();
        reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,A0,_Vi,outFileName);
        
        if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ))
        {
            double beta, se, pval;
            VectorXd y_buf=_y;
            if(mlma_preadj_covar)
                y_buf=_y.array()-(_X*einfo._b).array();
            for(int i=0;i<leftids.size();i++)
            {
                double desti=1.0*wcount/(einfo._epi_include.size()-1);
                if(desti>=cr)
                {
                    printf("%3.0f%%\r", 100.0*desti);
                    fflush(stdout);
                    if(cr==0) cr+=0.05;
                    else if(cr==0.05) cr+=0.2;
                    else if(cr==0.25) cr+=0.5;
                    else cr+=0.25;
                }

                int id=leftids[i];
                int targetChr=assoc_rlts[id].CHR;
                int targetBP=assoc_rlts[id].BP;
                string targetPrb=assoc_rlts[id].PROBE;
                string gene=assoc_rlts[id].GENE;
                char oren=assoc_rlts[id].OREN;
                    
                iter=einfo._epi_map.find(assoc_rlts[id].PROBE);
                VectorXd x_buf(einfo._eii_include.size());
                if(iter!=einfo._epi_map.end())
                {
                    int curidx=iter->second;
                    for(int j=0; j<einfo._eii_include.size(); j++)
                        x_buf(j)=einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                    } else {
                        LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                        TERMINATE();
                    }
                    
                    if(mlma_preadj_covar)
                        if(_X_c>1) getResidual(x_buf, _X);
                    
                    if(mlma_preadj_covar) mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                    else mlma_cal_stat_covar(y_buf, x_buf, _Vi, _X, beta, se, pval);
                    
                    string chrstr;
                    if(targetChr==23) chrstr="X";
                    else if(targetChr==24) chrstr="Y";
                    else chrstr=atosm(targetChr);
                    outstr = chrstr + '\t' + targetPrb + '\t' + atosm(targetBP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
                    if(fputs_checked(outstr.c_str(),ofile))
                    {
                        LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                        TERMINATE();
                    }
                    wcount++;
                }
            
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
        LOGPRINTF("Results of %ld probes have been saved in file %s.\n",wcount,filename.c_str());
        fclose(ofile);
        free_assoclist(assoc_rlts);
    }
    
    void moment_exact(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool mlma_preadj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_stepwise,int erm_alg, double swthresh, int tsk_ttl, int tsk_id)
    {
        //bcthresh: default as bonfferoni
        //r2thresh: default as -9, disabled.
        //percent: default as -9, disabled.
        //slctmtd: feature selection with linear or mlma. default as linear.
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        expect_wind>>=1;
        
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
        double secondthresh=1e-5;
        
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
            LOGPRINTF("ERROR: please input the Gene expression / Methylation data by the option --befile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("ERROR: please input the phenotype data by the option --pheno.\n");
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
        if(!moment_eligibility_ck(&einfo))
        {
            LOGPRINTF("ERROR: the .opi file contains NA chromosome or NA probe position. Plsease annotate it using flag --update-opi.\n");
            TERMINATE();
        }
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo); // eii and epi are updated in it.
        if(einfo._eType==0 && expect_wind==50) expect_wind=100;
        map<string, int>::iterator iter;
        
        int _n=(int)einfo._eii_include.size();
        if(_n<1)
        {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        LOGPRINTF("%d individuals are in common in the input files.\n",_n);
        VectorXd _y;
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
        
        /*
        //LOGPRINTF("Performing principal component analysis ...\n");
        make_erm( &einfo,erm_alg); //use the whole
     
        SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
        MatrixXd evec = (eigensolver.eigenvectors());
        VectorXd eval = eigensolver.eigenvalues();
        //if(einfo._epi_include.size()<=einfo._eii_include.size())
        if(einfo._epi_include.size()<=expect_pcs)
        {
            LOGPRINTF("WARNING: probe number is too few for PCs...\n");
            expect_pcs = (int)einfo._epi_include.size();
        }
        MatrixXd _PCs = evec.block(0,evec.cols()-expect_pcs, evec.rows(),expect_pcs);
        MatrixXd _X_PC(_X.rows(), _X.cols()+expect_pcs);
        _X_PC << _X, _PCs;
        */
        
        vector< vector<int> > mapids;
        vector< vector<string> > erm_prbs;
        vector<MatrixXd> AZERO, _A, A0;
        mapids.resize(nrandcomp); // one for significant one for insignificant
        erm_prbs.resize(nrandcomp);
        AZERO.resize(nrandcomp);
        
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        if(slctmtd)
        {
            bool flagtmp=testMOA(assoc_rlts, &einfo, erm_alg);
            if(!flagtmp) testLinear(assoc_rlts, NULL, &einfo,COV_plus);
        }
        else
        {
            testLinear(assoc_rlts, NULL, &einfo,COV_plus);
        }
        
        //testQAssoc(assoc_rlts,NULL, &einfo);
        ASSOCRLT* sortptr=&assoc_rlts[0];
        qsort(sortptr,assoc_rlts.size(),sizeof(ASSOCRLT),ascend_assoc);
        
        if(nrandcomp==2)
        {
            if(percent>0)
            {
                int number=ceil(assoc_rlts.size()*percent);
                for(int k=(int)(assoc_rlts.size()-1);k>=(int)(assoc_rlts.size()-number);k--)
                {
                    mapids[0].push_back(k); //sig
                    erm_prbs[0].push_back(assoc_rlts[k].PROBE);
                }
                for(int k=(int)(assoc_rlts.size()-number-1);k>=0;k--)
                {
                    mapids[1].push_back(k);
                    erm_prbs[1].push_back(assoc_rlts[k].PROBE);
                }
                LOGPRINTF("%ld probes incldued in the first component and %ld in the second component.\n",mapids[0].size(),mapids[1].size());
            }
            else
            {
                if(r2thresh<0 && bcthresh<0) bcthresh=0.05/assoc_rlts.size();
                else if(r2thresh > 0) {
                    double chi=einfo._eii_include.size()* r2thresh / (1 - r2thresh) + 1;
                    bcthresh = pchisq(chi, 1);
                }
                for(int k=(int)(assoc_rlts.size()-1);k>=0;k--)
                {
                    if(assoc_rlts[k].PVAL<=bcthresh) {
                        mapids[0].push_back(k); //sig
                        erm_prbs[0].push_back(assoc_rlts[k].PROBE);
                    }
                    else {
                        mapids[1].push_back(k);
                        erm_prbs[1].push_back(assoc_rlts[k].PROBE);
                    }
                }
                LOGPRINTF("%ld probes incldued in the first component and %ld in the second component.\n",mapids[0].size(),mapids[1].size());
            }
            
        }
        else if(nrandcomp==3)
        {
            double bonthresh=0.05/assoc_rlts.size();
            
            for(int k=(int)(assoc_rlts.size()-1);k>=0;k--)
            {
                if(assoc_rlts[k].PVAL<=bonthresh)
                {
                    mapids[0].push_back(k);
                    erm_prbs[0].push_back(assoc_rlts[k].PROBE);
                }
                else if(assoc_rlts[k].PVAL<=secondthresh)
                {
                    mapids[1].push_back(k);
                    erm_prbs[1].push_back(assoc_rlts[k].PROBE);
                }
                else
                {
                    mapids[2].push_back(k);
                    erm_prbs[2].push_back(assoc_rlts[k].PROBE);
                }
            }
            LOGPRINTF("%ld probes included in the 1st random component.\n",mapids[0].size());
            LOGPRINTF("%ld probes included in the 2nd random component.\n",mapids[1].size());
            LOGPRINTF("%ld probes included in the 3rd random component.\n",mapids[2].size());
        }
        
        if( approximate_stepwise && mapids[0].size()>0)
        {
            // stepwise for the components
            vector<int> slctids, rmsig;
            if(swthresh<0) swthresh = bcthresh;
            stepwise_slct(&einfo,mapids[0], slctids, rmsig, assoc_rlts,  swthresh, true);
            erm_prbs[0].clear();
            for(int j=0;j<mapids[0].size();j++)
            {
                erm_prbs[0].push_back(assoc_rlts[mapids[0][j]].PROBE);
            }
            for(int j=0;j<rmsig.size();j++)
            {
                mapids[1].push_back(rmsig[j]);
                erm_prbs[1].push_back(assoc_rlts[rmsig[j]].PROBE);
            }
            string tmpfname=string(outFileName)+".AS.txt";
            
            FILE* tmpfile=fopen(tmpfname.c_str(),"w");
            if(!tmpfile)
            {
                LOGPRINTF("error open file %s.\n",tmpfname.c_str());
                TERMINATE();
            }
            for(int i=0;i< erm_prbs[0].size();i++) {
                string str=atos( erm_prbs[0][i]);
                str+='\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
        }
        
        vector<int> include_o(einfo._epi_include);
        map<string, int> snp_name_map_o(einfo._epi_map);
        int rawORMc=0;
        for(int i=0;i<mapids.size();i++) if(mapids[i].size()) rawORMc++;
        A0.resize(rawORMc+1);
        rawORMc=0;
        for(int j=0;j<erm_prbs.size();j++)
        {
            if(erm_prbs[j].size()>0)
            {
                update_map_kp(erm_prbs[j], einfo._epi_map, einfo._epi_include);
                make_erm(&einfo,erm_alg);
                (AZERO[j]).resize(_n, _n);
                (A0[rawORMc]).resize(_n, _n);
                #pragma omp parallel for
                for(int k=0; k<_n; k++)
                {
                    for(int l=0; l<=k; l++){
                        (AZERO[j])(l,k)=(AZERO[j])(k,l)=einfo._grm(k,l)*einfo._grm_N(k,l);
                        (A0[rawORMc])(l,k)=(A0[rawORMc])(k,l)=einfo._grm(k,l);
                    }
                }
                einfo._epi_include=include_o;
                einfo._epi_map=snp_name_map_o;
                rawORMc++;
            }
        }
        A0[rawORMc]=MatrixXd::Identity(_n, _n);
        
        LOGPRINTF("\nPerforming the MOMENT analysis ...\n");
        LOGPRINTF("For each probe, the analysis will exclude probes in %d Kb region of centered at the probe to be tested.\n",expect_wind);
        int exwind=expect_wind*1000;
        
        string filename=string(outFileName)+".mlma";
        if(tsk_ttl>1) filename=string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id)+".mlma";
        FILE* ofile = fopen(filename.c_str(), "w");
        if (!(ofile)) {
            LOGPRINTF("ERROR: open error %s\n", filename.c_str());
            TERMINATE();
        }
        string outstr="Chr\tProbe\tbp\tGene\tOrientation\tb\tse\tp\n";
        if(fputs_checked(outstr.c_str(),ofile))
        {
            LOGPRINTF("ERROR: error in writing file %s .\n", filename.c_str());
            TERMINATE();
        }
        
        long wcount=0;
        
        int start=0,end=(int)assoc_rlts.size()-1;
        if(tsk_ttl>1) {
            update_startend((int)assoc_rlts.size(),  tsk_ttl,  tsk_id,  start, end);
        }
        int nrun=end-start+1;
        MatrixXd XX(_X.rows(), _X.cols()+1);
        XX.block(0,0,_X.rows(),_X.cols())=_X;
        int x_idx=_X_c++;
        
        double cr=0;
        #pragma omp parallel for private(remlstatus)
        for(int id=start;id<=end;id++)
        {
            double desti=1.0*wcount/nrun;
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            
            int targetChr=assoc_rlts[id].CHR;
            int targetBP=assoc_rlts[id].BP;
            string targetPrb=assoc_rlts[id].PROBE;
            string gene=assoc_rlts[id].GENE;
            char oren=assoc_rlts[id].OREN;
            double beta=0, se=1;
            MatrixXd X=XX;
            double nonmiss=0.0;
            double mu=0.0, sd=0.0;
            long n=einfo._eii_include.size();
            iter=einfo._epi_map.find(assoc_rlts[id].PROBE);
            VectorXd x(n);
            if(iter!=einfo._epi_map.end())
            {
                int curidx=iter->second;
                for(int j=0; j<einfo._eii_include.size(); j++) {
                    double val =einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                    if(val<1e9){
                        mu+=val;
                        nonmiss+=1.0;
                    }
                }
                if(nonmiss>1)
                {
                    mu/=nonmiss;
                    for(int j = 0; j < n; j++) {
                        double val =einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
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
                        //LOGPRINTF("ERROR: no variance of probe %s.\n",targetPrb.c_str());
                        //do something here;
                        double pval=1;
                        string chrstr;
                        if(targetChr==23) chrstr="X";
                        else if(targetChr==24) chrstr="Y";
                        else chrstr=atosm(targetChr);
                        outstr = chrstr + '\t' + targetPrb + '\t' + atosm(targetBP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
                        if(fputs_checked(outstr.c_str(),ofile))
                        {
                            LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                            TERMINATE();
                        }
                        wcount++;
                        continue;
                    }
                }
                else
                {
                    LOGPRINTF("ERROR: too many missing values in probe %s.\n",targetPrb.c_str());
                    TERMINATE();
                }
            }
            else
            {
                LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                TERMINATE();
            }
            X.col(x_idx)=x;
            
            for(int jj=0;jj<nrandcomp;jj++) erm_prbs[jj].clear();
            int norms=0,excount=0;
            for(int jj=0;jj<nrandcomp;jj++)
            {
                for(int j=0;j<mapids[jj].size();j++)
                {
                    int tid=mapids[jj][j];
                    if(targetChr == assoc_rlts[tid].CHR && abs(targetBP-assoc_rlts[tid].BP)<=exwind) {
                        erm_prbs[jj].push_back(assoc_rlts[tid].PROBE);
                        excount++;
                    }
                }
            }
            //LOGPRINTF("%d probe(s) excluded from calculating the ORM(s) for the probe %s.\n", excount,assoc_rlts[id].PROBE);
            for(int jj=0;jj<nrandcomp;jj++) if(mapids[jj].size()-erm_prbs[jj].size()>0) norms++;
            
            einfo._r_indx.clear();
            for(int j=0; j < norms + 1; j++) einfo._r_indx.push_back(j);
            _A.resize(einfo._r_indx.size());
            int aid=0;
            for(int j=0; j < erm_prbs.size(); j++)
            {
                long N=mapids[j].size()-erm_prbs[j].size();
                if(N)
                {
                    if(erm_prbs[j].size()>0)
                    {
                        update_map_kp(erm_prbs[j], einfo._epi_map, einfo._epi_include);
                        make_erm(&einfo,erm_alg);
                        (_A[aid]).resize(_n, _n);
                        #pragma omp parallel for
                        for(int k=0; k<_n; k++)
                        {
                            for(int l=0; l<=k; l++) (_A[aid])(l,k)=(_A[aid])(k,l)=(AZERO[j](k,l)-einfo._grm(k,l)*einfo._grm_N(k,l))/N;
                        }
                        einfo._epi_include=include_o;
                        einfo._epi_map=snp_name_map_o;
                    } else {
                        _A[aid]=AZERO[j]/N;
                    }
                    aid++;
                }
            }
            _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
            
            // names of variance component
            einfo._var_name.clear();
            einfo._hsq_name.clear();
            for (int j = 0; j < norms; j++) {
                stringstream strstrm;
                if (norms == 1) strstrm << "";
                else strstrm << j + 1;
                einfo._var_name.push_back("V(O" + strstrm.str() + ")");
                einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
            }
            einfo._var_name.push_back("V(e)");
            
            einfo._within_family=within_family;
            if(within_family) detect_family(&einfo, _A);
            MatrixXd _Vi;
            remlstatus=0; //reset reml status
            if(norms) reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, X,_y,_A,_Vi,outFileName);
            if(remlstatus==0 || remlstatus==-5 || remlstatus==-3)
            {
                beta=einfo._b[_X_c-1];
                se=sqrt(einfo._se(_X_c-1));
                
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
            double chisq=beta/se;
            double pval=pchisq(chisq*chisq, 1);
            string chrstr;
            if(targetChr==23) chrstr="X";
            else if(targetChr==24) chrstr="Y";
            else chrstr=atosm(targetChr);
            outstr = chrstr + '\t' + targetPrb + '\t' + atosm(targetBP) + '\t' + atos(gene) + '\t' + atos(oren) + '\t' + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
            if(fputs_checked(outstr.c_str(),ofile))
            {
                LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                TERMINATE();
            }
            wcount++;
        }
        
        LOGPRINTF("Results of %ld probes have been saved in file %s.\n",wcount,filename.c_str());
        fclose(ofile);
        free_assoclist(assoc_rlts);
    }
}

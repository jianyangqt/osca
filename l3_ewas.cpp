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
     void backward_elimn(eInfo* einfo, vector<int> &slct, MatrixXd &prob_profile, double p_thresh, double fdr_thresh, vector<int> &dump, double &rsq )
    {
        int nindi=(int)einfo->_eii_include.size();
        int nCOV=1;
        while(slct.size()>1)
        {
            int X_c=nCOV+(int)slct.size();;
            MatrixXd X(nindi, X_c);
            X.block(0, 0, nindi, nCOV) = MatrixXd::Ones(nindi, 1);
            for(int i=0;i<slct.size();i++) X.col(i+nCOV)=prob_profile.col(slct[i]);
            MatrixXd XtX_i;
            XtX_i=X.transpose()*X;
            bool determinant_zero=false;
            inverse_V(XtX_i, determinant_zero);
            if(determinant_zero)
            {
                LOGPRINTF("Error: please input a strigent threshold for stepwise selection.\n");
                TERMINATE();
            }
            VectorXd y(einfo->_eii_include.size());
            #pragma omp parallel for
            for(int j=0; j<einfo->_eii_include.size(); j++)
            {
                y(j)=einfo->_eii_pheno[einfo->_eii_include[j]];
            }
            VectorXd b_hat=XtX_i*X.transpose()*y;
            double sst=y.transpose()*y;
            double ssr=y.transpose()*X*b_hat;
            rsq=ssr/sst;
            rsq=1-(1-rsq)*(nindi-1)/(nindi-X_c); // adjusted R-squared
            if(loud)
            {
                printf("Adjusted R-squared: %f explained by %ld selected probes.\n",rsq, slct.size());
            }
            VectorXd residual=(y-X*b_hat);
            residual=residual.array()*residual.array();
            double sy=sqrt(residual.sum()/(y.size()-X.cols()));
            VectorXd se=sy*XtX_i.diagonal().array().sqrt();
            VectorXd t=b_hat.array()/se.array();
            vector<double> p(t.size()-1);
            for(int j=1;j<t.size();j++) p[j-1]=pchisq( t(j)*t(j),1 );
            int m = (int)(max_element(p.begin(), p.end()) - p.begin());
            if(p_thresh<0)
            {
                vector<double> jointq;
                getQval(p,jointq);
                if(jointq[m]>fdr_thresh)
                {
                    dump.push_back(slct[m]);
                    slct.erase(slct.begin()+m);
                } else break;
            }
            else
            {
                if(p[m]>p_thresh)
                {
                    dump.push_back(slct[m]);
                    slct.erase(slct.begin()+m);
                } else break;
            }
            
        }
        
    }
    bool forward_slct(eInfo* einfo, vector<int> &slct,vector<int> &remain, MatrixXd &prob_profile, double p_thresh, double fdr_thresh, vector<int> &dump , double rsqthresh)
    {
        vector<double> rlt;
        int nCOV=1;
        int nindi=(int)einfo->_eii_include.size(), X_c=nCOV+(int)slct.size();;
        MatrixXd X(nindi, X_c);
        X.block(0, 0, nindi, nCOV) = MatrixXd::Ones(nindi, 1);
        for(int i=0;i<slct.size();i++) X.col(i+nCOV)=prob_profile.col(slct[i]);
        
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
        double sst=y.transpose()*y;
        double ssr=y.transpose()*X*b_hat;
        double rsq=ssr/sst;
        rsq=1-(1-rsq)*(nindi-1)/(nindi-X_c); // adjusted R-squared
        if(loud)
        {
            printf("Adjusted R-squared: %f explained by %ld selected probes.\n",rsq, slct.size());
        }
        if(rsq > rsqthresh) return false;
        
        VectorXd yresi=(y-X*b_hat);
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
        if(p_thresh<0)
        {
            vector<double> condq;
            getQval(condp,condq);
            if(condq[m]>fdr_thresh)
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
        else
        {
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
    void stepwise_slct(eInfo* einfo, vector<int> &sig, vector<int> &slctid, vector<int> &rmid, vector<ASSOCRLT> &assoc_rlts, double p_cutoff, double fdr_cutoff, bool updatesig, bool stepforwardonly, double rsqthresh)
    {
        // as long as p_cutoff>0, fdr_cutoff is disabled.
        vector<double> p_buf;
        map<string, int>::iterator iter;
        MatrixXd sig_profile(einfo->_eii_include.size(),sig.size());
        double p=1;int m=0;
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
        if (p_cutoff > 0)
        {
            if(assoc_rlts[sig[m]].PVAL >= p_cutoff) return;
        }
        else
        {
            vector<double> ptmp(sig.size()),qtmp;
            for(int i=0;i<sig.size();i++) ptmp[i]=assoc_rlts[sig[i]].PVAL;
            getQval(ptmp, qtmp);
            if(qtmp[m] >= fdr_cutoff) return;
        }
        slctidx.push_back(m);// m must be 0 because assor_rlts is already sorted.
        for (int i = 1; i < sig.size(); i++) {
             remain.push_back(i);
        }
        while (!remain.empty()) {
            if (forward_slct(einfo,slctidx, remain, sig_profile,p_cutoff,fdr_cutoff, rmid,rsqthresh)) {
                if(loud){
                    LOGPRINTF("Forward: Selected %ld, remained %ld, removed %ld \n",slctidx.size(),remain.size(),rmid.size());
                }
                if(!stepforwardonly) {
                    double rsq = 0.0;
                    backward_elimn(einfo,slctidx, sig_profile,p_cutoff,fdr_cutoff, rmid, rsq);
                    if(rsq > rsqthresh)
                    {
                        if(loud){
                            LOGPRINTF("Multiple R-squared: %f\n",rsq);
                        }
                        for(int j=0;j<remain.size();j++) rmid.push_back(remain[j]);
                        remain.clear();
                    }
                }
                if(!stepforwardonly && loud){
                    LOGPRINTF("Backward: Selected %ld, remained %ld, removed %ld \n",slctidx.size(),remain.size(),rmid.size());
                }
            } else {
                for(int j=0;j<remain.size();j++) rmid.push_back(remain[j]);
                remain.clear();
                if(loud){
                    LOGPRINTF("Forward: Selected %ld, remained %ld, removed %ld \n",slctidx.size(),remain.size(),rmid.size());
                }
                
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
            LOGPRINTF("%ld probes are selected from stepwise selection.\n", slctidx.size());
        }
    }
    void cal_cor(eInfo* einfo, int target, VectorXd &pcc )
    {
        uint64_t n = einfo->_eii_include.size(), m = einfo->_epi_include.size();
        VectorXd xt(n),yt(n);
        pcc.resize(m);
        for(int j=0; j<n; j++)
        {
            double val = einfo->_val[target*einfo->_eii_num+einfo->_eii_include[j]];
            if(val>1e9)
                val = 0;
            xt(j) = val;
        }
        #pragma omp parallel for
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++)
            {
                double val = einfo->_val[i*einfo->_eii_num+einfo->_eii_include[j]];
                if(val>1e9)
                    val = 0;
                yt(j) = val;
            }
            pcc[i] = cor(xt,yt,true,true);
        }
    }
    void moment(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool nopreadj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_num, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_flag, bool approximate_stepwise,int erm_alg, double swthresh, double swfdr, bool swlogit,bool stepforwardonly, bool Baptiste, double sw_rsq, double thresh_pcc, char* orm_file,bool orm_bin_flag, bool force_moment)
    {
       
        //bcthresh: default as bonfferoni
        //r2thresh: default as -9, disabled.
        //percent: default as -9, disabled.
        //slctmtd: feature selection with linear or mlma. default as linear.
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        vector<string> orm_files;
        vector<string> orm_id;
        vector< vector<int> > mapids;
        vector< vector<string> > erm_prbs;
        vector<MatrixXd> AZERO, _A, A0;
        vector<MatrixXf> AN;
        MatrixXd _Vi;
        mapids.resize(nrandcomp); // one for significant one for insignificant
        erm_prbs.resize(nrandcomp);
        AZERO.resize(nrandcomp);
        AN.resize(nrandcomp);
        
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
        if(!Baptiste && !moment_eligibility_ck(&einfo))
        {
            LOGPRINTF("NOTE: the .opi file contains NA chromosome or NA probe position. Window-based exclusion switches to correlation-based exclusion.\n");
            Baptiste = true;
        }
        epi_man(&einfo,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,probe2exclde);
        if(phenofileName !=NULL) read_phen(&einfo, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&einfo, covfileName, false);
        if(qcovfileName != NULL) read_cov(&einfo, qcovfileName, true);
        if(orm_file!=NULL){
            orm_files.push_back(orm_file);
            read_grm(&einfo,orm_file, orm_id, true, false, false,orm_bin_flag);
            AN[nrandcomp-1] = einfo._grm_N;
            AZERO[nrandcomp-1] = einfo._grm.array()*einfo._grm_N.cast<double>().array();
        }
        memcpy(suffix,".bod",5);
        read_beed(inputname,&einfo); // eii and epi are updated in it.
        stdprobe(&einfo);
        if(einfo._eType==0 && expect_wind==50) expect_wind=100;
        map<string, int>::iterator iter, iter2;
    
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
        
        
        //LOGPRINTF("Performing principal component analysis ...\n");
        //make_erm( &einfo,erm_alg); //use the whole
        //SelfAdjointEigenSolver<MatrixXd> eigensolver(einfo._grm.cast<double>());
        //MatrixXd evec = (eigensolver.eigenvectors());
        //VectorXd eval = eigensolver.eigenvalues();
        //if(einfo._epi_include.size()<=expect_pcs)
        //{
       //     LOGPRINTF("WARNING: probe number is too few for PCs...\n");
       //     expect_pcs = (int)einfo._epi_include.size();
       // }
         
        //MatrixXd _PCs = evec.block(0,evec.cols()-expect_pcs, evec.rows(),expect_pcs);
        //MatrixXd _X_PC(_X.rows(), _X.cols()+expect_pcs);
        //_X_PC << _X, _PCs;
        
        
        map<double, long> levels;
        long cursize=0,flag=0;
        
        bool qtrait=false;
        for(int i=0;i<_n;i++)
        {
            levels.insert(pair<double,int>(_y(i),levels.size()));
            
            if(levels.size()>2)
            {
                qtrait=true;
                break;
            }
            if(levels.size()>cursize) {
                flag+=(int)_y(i);
                cursize = levels.size();
            }
        }
        if(swlogit && !qtrait)
        {
            for(int i=0; i<_n; i++){
                if(flag==2 && _y(i)==2)
                {
                    einfo._eii_pheno[einfo._eii_include[i]]= 1;
                    _y(i)=1;
                }
                else if(flag==3)
                {
                    einfo._eii_pheno[einfo._eii_include[i]]-=1;
                    _y(i)=-1;
                }
            }
        }
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        if(slctmtd)
        {
            bool flagtmp=testMOA(assoc_rlts, &einfo, erm_alg);
            if(!flagtmp) {
                if(swlogit && !qtrait) testLogit(assoc_rlts, NULL, &einfo,COV_plus);
                else testLinear(assoc_rlts,NULL, &einfo,COV_plus);
            }
        }
        else
        {
            if(swlogit && !qtrait) testLogit(assoc_rlts, NULL, &einfo,COV_plus);
            else testLinear(assoc_rlts,NULL, &einfo,COV_plus);
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
                        double chi=einfo._eii_include.size()* r2thresh / (1 - r2thresh) ;
                        bcthresh = pchisq(chi, 1);
                        if(loud)
                        {
                            LOGPRINTF("Using  Rsq = %6.2f (equivalent p-value threshod = %6.2e) as the criterion for the feature selection.\n",r2thresh,bcthresh);
                        }
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
        double RSQVO =1;
        for(int i=0; i < 2; i++) einfo._r_indx.push_back(i);
        if(orm_file==NULL)
        {
            make_erm(&einfo,AZERO[nrandcomp-1], erm_alg, true, NULL, false, true);
            AN[nrandcomp-1]=einfo._grm_N;
        }
        
        if( approximate_stepwise && mapids[0].size()>0)
        {
            _A.resize(einfo._r_indx.size());
            (_A[0]).resize(_n, _n);
            #pragma omp parallel for
            for(int i=0; i<_n; i++)
            {
                for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=einfo._grm(i,j);
            }
            _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
            einfo._var_name.clear();
            einfo._hsq_name.clear();
            for (int i = 0; i < 1; i++) {
                stringstream strstrm;
                strstrm << "";
                einfo._var_name.push_back("V(O" + strstrm.str() + ")");
                einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
            }
            einfo._var_name.push_back("V(e)");
            reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
            _A.clear();
            einfo._r_indx.clear();
            reml_priors_var.clear();
            if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ))
            {
                RSQVO=einfo._varcmp[0]/(einfo._varcmp[0] + einfo._varcmp[1]);
                einfo._varcmp_Py = _Vi;
            }
            else RSQVO = 1;
            sw_rsq *= RSQVO;
         
            // stepwise for the components
                vector<int> slctids, rmsig;
                if(swfdr<0 && swthresh<0) swthresh = PFISHER/mapids[0].size();
                if(loud)
                {
                    if(swthresh>0) {LOGPRINTF("Using p-value threshod = %6.2e as the criterion for the stepwise selection.\n",swthresh);}
                    else {LOGPRINTF("Using FDR threshod = %6.2e as the criterion for the stepwise selection.\n",swfdr);}
                }
                if(_X_c>1)
                {
                    VectorXd y_buf =_y;
                    getResidual(y_buf, _X);
                    for(int i=0; i<_n; i++){
                        einfo._eii_pheno[einfo._eii_include[i]]=y_buf[i];
                    }
                    //if(loud) {LOGPRINTF("Phenotype reset.\n");}
                } else {
                    // center y for RSQ calculating
                    double ymean=_y.mean();
                    VectorXd y_buf =_y.array() - ymean;
                    for(int i=0; i<_n; i++){
                        einfo._eii_pheno[einfo._eii_include[i]]=y_buf[i];
                    }
                }
                stepwise_slct(&einfo, mapids[0], slctids, rmsig, assoc_rlts,  swthresh,swfdr, true, stepforwardonly,sw_rsq);
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
                
                for(int i=0; i<_n; i++){
                    einfo._eii_pheno[einfo._eii_include[i]]=_y[i];
                    //if(loud) {LOGPRINTF("Phenotype restored.\n");}
                }
        }
        
        vector<int> include_o(einfo._epi_include);
        map<string, int> snp_name_map_o(einfo._epi_map);
        map<string,int> componet1_map;
        for(int i=0;i< erm_prbs[0].size();i++) componet1_map.insert(pair<string,int>(erm_prbs[0][i],mapids[0][i]));
        int rawORMc=0;
        for(int i=0;i<mapids.size();i++) if(mapids[i].size()) rawORMc++;
        A0.resize(rawORMc+1);
        rawORMc=0;
        for(int j=0;j<erm_prbs.size()-1;j++)
        {
            if(erm_prbs[j].size()>0)
            {
                update_map_kp(erm_prbs[j], einfo._epi_map, einfo._epi_include);
                make_erm(&einfo,AZERO[j], erm_alg, true, NULL, false, true);
                AN[j] = einfo._grm_N;
                A0[rawORMc] = einfo._grm;
                
                //make_erm(&einfo, erm_alg);
                //(AZERO[j]).resize(_n, _n);
                //(A0[rawORMc]).resize(_n, _n);
                //#pragma omp parallel for
                //for(int k=0; k<_n; k++)
                //{
                 //   for(int l=0; l<=k; l++){
                  //      (AZERO[j])(l,k)=(AZERO[j])(k,l)=einfo._grm(k,l)*einfo._grm_N(k,l);
                  //       (A0[rawORMc])(l,k)=(A0[rawORMc])(k,l)=einfo._grm(k,l);
                  //  }
                //}
                 
                einfo._epi_include=include_o;
                einfo._epi_map=snp_name_map_o;
                rawORMc++;
            }
        }
        long j =erm_prbs.size()-1;
        if(erm_prbs[j].size()>0)
        {
            for(int ii=0;ii<erm_prbs.size()-1;ii++)
            {
                if(erm_prbs[ii].size()>0)
                {
                    AN[j] -= AN[ii];
                    AZERO[j] -= AZERO[ii];
                }
            }
            A0[rawORMc].resize(_n,_n);
            #pragma omp parallel for
            for (int ii = 0; ii < _n; ii++) {
                for (int jj = 0; jj <= ii; jj++) {
                    if(AN[j](ii,jj) > 0) A0[rawORMc](ii,jj) = AZERO[j](ii,jj) / AN[j](ii,jj);
                    else A0[rawORMc](ii,jj) = 0.0;
                    A0[rawORMc](jj,ii) = A0[rawORMc](ii,jj);
                }
            }
            rawORMc++;
            
        }
        A0[rawORMc]=MatrixXd::Identity(_n, _n);
        /*
        if(loud)
        {
            LOGPRINTF("Saving the A matrix ...\n");
            for(int ii=0;ii<A0.size()-1;ii++)
            {
                string filenam=string(outfileName)+".A"+atos(ii)+".mat";
                FILE* tmpfil=fopen(filenam.c_str(),"w");
                if(!tmpfil)
                {
                    LOGPRINTF("error open file.\n");
                    TERMINATE();
                }
                for(int t=0;t<A0[ii].rows();t++)
                {
                    string str="";
                    for(int k=0;k<A0[ii].cols();k++)
                    {
                        str +=atos(A0[ii](t,k)) + '\t';
                    }
                    str += '\n';
                    fputs(str.c_str(),tmpfil);
                }
                fclose(tmpfil);
                LOGPRINTF("The A matrix is saved in the file %s.\n",filenam.c_str());
            }
            if(A0.size()>1)
            {
                for(int tt=0;tt<2;tt++)
                {
                    for(int i=0; i < 2; i++) einfo._r_indx.push_back(i);
                    _A.resize(einfo._r_indx.size());
                    _A[0]=A0[tt];
                    _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
                    einfo._var_name.clear();
                    einfo._hsq_name.clear();
                    for (int i = 0; i < 1; i++) {
                        stringstream strstrm;
                        strstrm << "";
                        einfo._var_name.push_back("V(O" + strstrm.str() + ")");
                        einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
                    }
                    einfo._var_name.push_back("V(e)");
                    string tmpn = string(outFileName) + ".A" + atos(tt);
                    reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,tmpn);
                    _A.clear();
                    einfo._r_indx.clear();
                    reml_priors_var.clear();
                }
            }
        }
        */
        LOGPRINTF("\nPerforming the MOMENT analysis ...\n");
        LOGPRINTF("For each probe, the analysis will exclude probes in %d Kb region of centered at the probe to be tested.\n",expect_wind);
        int exwind=expect_wind*1000;
       
        string filename=string(outFileName)+".moment";
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
        int outid=0;
        
        if(approximate_flag && mapids[0].size()>0)
        {
            
            if(loud) {LOGPRINTF("Proiritize probes in the first componet.\n");}
            vector<int> slcted, rmsig;
            if(swfdr<0 && swthresh<0) swthresh = PFISHER/mapids[0].size();
                stepwise_slct(&einfo, mapids[0],slcted, rmsig, assoc_rlts, swthresh, swfdr, false, stepforwardonly,sw_rsq);
                for(int j=0; j<slcted.size();j++) mapids[0][j]=slcted[j];
                for(int j=0; j<rmsig.size();j++) mapids[0][j+slcted.size()]=rmsig[j];
                expect_num=(int)slcted.size();
        }
        
        if(expect_num<0) expect_num = (int)mapids[outid].size();
        long num2exp=expect_num>mapids[outid].size()?mapids[outid].size():expect_num;
        
        MatrixXd XX(_X.rows(), _X.cols()+1);
        XX.block(0,0,_X.rows(),_X.cols())=_X;
        int x_idx=_X_c;
        double cr=0;
        vector<int> freml;
        FILE* batist = NULL;
        if(Baptiste && loud)
        {
            string filenam=string(outfileName)+".excluded.txt";
             batist=fopen(filenam.c_str(),"w");
            if(!batist)
            {
                LOGPRINTF("error open file.\n");
                TERMINATE();
            }
        }
        
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
                    {
                        double val = einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                        x_buf(j) = val;
                        if(val>1e9)
                            val = 0;
                        X(j,x_idx) = val;
                    }
                } else {
                    LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                    TERMINATE();
                }
            
            for(int jj=0;jj<nrandcomp;jj++) erm_prbs[jj].clear();
            int norms=0,excount=0;
            if(Baptiste)
            {
                //calculate correlations
                VectorXd pcc;
                int curidx=iter->second;
                cal_cor(&einfo,curidx,pcc);
                for(int i=0;i<pcc.size();i++)
                {
                    double pccr2 = pcc[i]*pcc[i];
                    if(pccr2>thresh_pcc)
                    {
                        string prob=einfo._epi_prb[einfo._epi_include[i]];
                        int probchr=einfo._epi_chr[einfo._epi_include[i]];
                        int probbp=einfo._epi_bp[einfo._epi_include[i]];
                        iter2=componet1_map.find(prob);
                        if(iter2!=componet1_map.end()) erm_prbs[0].push_back(prob);
                        else erm_prbs[1].push_back(prob);
                        excount++;
                        
                        if(loud)
                        {
                            string str = targetPrb + '\t' + atos(targetChr) + '\t'+ atos(targetBP) + '\t' + prob + '\t' + atos(probchr) + '\t'+ atos(probbp) + '\t' + atos(pcc[i]) + '\n';
                            fputs(str.c_str(),batist);
                        }
                    }
                }
            } else {
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
                            make_erm(&einfo,_A[aid], erm_alg, true, NULL, false, true);
                            einfo._grm_N = AN[j] - einfo._grm_N;
                            _A[aid] = AZERO[j] - _A[aid];
                            #pragma omp parallel for
                            for (int ii = 0; ii < _n; ii++) {
                                for (int jj = 0; jj <= ii; jj++) {
                                    if(einfo._grm_N(ii,jj) > 0) _A[aid](ii,jj) /= einfo._grm_N(ii,jj);
                                    else _A[aid](ii,jj) = 0.0;
                                    _A[aid](jj,ii) = _A[aid](ii,jj);
                                }
                            }
                            einfo._epi_include=include_o;
                            einfo._epi_map=snp_name_map_o;
                            
                        } else {
                            _A[aid].resize(_n, _n);
                            #pragma omp parallel for
                            for (int ii = 0; ii < _n; ii++) {
                                for (int jj = 0; jj <= ii; jj++) {
                                    if (AN[j](ii,jj) > 0.0) _A[aid](ii,jj) = AZERO[j](ii,jj)/AN[j](ii,jj);
                                    else _A[aid](ii,jj) = 0.0;
                                    _A[aid](jj,ii) = _A[aid](ii,jj);
                                }
                            }
                        }
                        //_A[aid] /= _A[aid].diagonal().mean();// disabled by Futao on 4 Jun, 2019. because of the update in make_erm().
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
                    double rox=0;
                    for( int ii=0;ii<einfo._varcmp.size()-1;ii++)
                        rox +=einfo._varcmp[ii];
                    rox = rox / (rox + einfo._varcmp[einfo._varcmp.size()-1]); // if rox is too small, it means REML mihgt fail in two componets.
                    if(!approximate_stepwise || force_moment || rox/RSQVO >= 0.05)
                    {
                        einfo._P.resize(0,0);
                        _A.clear();
                        VectorXd y_buf=_y;
                        if(!nopreadj_covar)
                        {
                            y_buf=_y.array()-(_X*einfo._b).array();
                            if(_X_c>1) getResidual(x_buf, _X);
                        }
                        if(!nopreadj_covar) mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
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
                    else
                    {
                        if(loud) printf("%s needs recalcualtion.\n", targetPrb.c_str());
                        freml.push_back(id);
                    }
                    
                }
                else
                {
                    freml.push_back(id);
                }
        }
        if(loud && batist!=NULL) fclose(batist);
        if(loud)
        {
            LOGPRINTF("%ld probes in the first componet failed in REML.\n", freml.size());
        }
        
        long leftnum = mapids[outid].size()-num2exp+freml.size();
        for(int otid=1;otid<mapids.size();otid++) leftnum += mapids[otid].size();
        vector<int> leftids(leftnum);
        long count = 0;
        for(int i=0;i<freml.size();i++) {
            leftids[count++] = freml[i];
        }
        if(mapids[outid].size()>num2exp)
            for(long i=num2exp;i<mapids[outid].size();i++)
            {
                int id=mapids[outid][i];
                leftids[count] = id;
                count++;
            }
        for(int otid=1;otid<mapids.size();otid++)
        {
            for(int i=0;i<mapids[otid].size();i++)
            {
                int id=mapids[otid][i];
                leftids[count] = id;
                count++;
            }
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
        remlstatus=0; //reset reml status
        remloasi = false;
        reml_priors_var.clear();
        loud = true;
        reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,A0,_Vi,outFileName);
        bool rsqck=false;
        if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ))
        {
            double rox=0;
            for( int ii=0;ii<einfo._varcmp.size()-1;ii++)
                rox +=einfo._varcmp[ii];
            rox = rox / (rox + einfo._varcmp[einfo._varcmp.size()-1]);
            //if(loud) printf("%f:%f\n",rox,RSQVO);
            if(!approximate_stepwise || force_moment || rox/RSQVO >= 0.05 || mapids[0].size()==0) rsqck=true;
        }
        if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ) && rsqck)
        {
            double beta, se, pval;
            MatrixXd Vi_C;
            MatrixXd A;
            VectorXd Ct_vi_y;
            
            VectorXd y_buf=_y;
            if(!nopreadj_covar)
                y_buf=_y.array()-(_X*einfo._b).array();
            else {
                Vi_C = _Vi*_X; // n by m
                A = _X.transpose()*Vi_C; // Ct_Vi_C : m by m
                bool determinant_zero=false;
                inverse_V(A, determinant_zero);
                if(determinant_zero)
                {
                    LOGPRINTF("The matrix is not invertible.\n");
                    TERMINATE();
                }
                Ct_vi_y = _X.transpose()*_Vi*y_buf;
            }
            #ifndef __APPLE__
            int n = (int)y_buf.size();
            float* y_mkl = new float[n];
            float* x_mkl = new float[n];
            float* Vi_mkl = new float[n*n];
            #pragma omp parallel for
            for (int i = 0; i < n; i++) y_mkl[i] = y_buf[i];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Vi_mkl[i * n + j] = _Vi(i, j);
                }
            }
            #endif
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
                    for(int j=0; j<einfo._eii_include.size(); j++) x_buf(j) = einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                } else {
                        LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                        TERMINATE();
                }
                    
                    if(!nopreadj_covar)
                        if(_X_c>1)
                            getResidual(x_buf, _X);
                #ifndef __APPLE__
                for(int j=0; j<einfo._eii_include.size(); j++) x_mkl[j]= x_buf(j);
                #endif
                if(!nopreadj_covar) {
                    #ifndef __APPLE__
                    mlma_cal_stat_mkl(y_mkl, x_mkl, Vi_mkl,n, beta, se, pval);
                    #else
                    mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                    #endif
                }
                else mlma_cal_stat_covar(y_buf, x_buf, _Vi, Vi_C, A,Ct_vi_y, beta, se, pval);
                    
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
            #ifndef __APPLE__
            delete[] y_mkl;
            delete[] x_mkl;
            delete[] Vi_mkl;
             #endif
            }
        else
        {
            
           LOGPRINTF("Degrade to MOA\n");
           if(!approximate_stepwise && mapids[0].size()>0)
           {
               einfo._r_indx.clear();
               for(int i=0; i < 2; i++) einfo._r_indx.push_back(i);
              
               _A.resize(einfo._r_indx.size());
               (_A[0]).resize(_n, _n);
               for(int ii=1;ii<AN.size();ii++)
               {
                   AZERO[0] += AZERO[ii];
                   AN[0] += AN[ii];
               }
                #pragma omp parallel for
               for (int ii = 0; ii < _n; ii++) {
                   for (int jj = 0; jj <= ii; jj++) {
                       if(AN[0](ii,jj) > 0) _A[0](ii,jj) = AZERO[0](ii,jj) / AN[0](ii,jj);
                       else _A[0](ii,jj) = 0.0;
                       _A[0](jj,ii) = _A[0](ii,jj);
                   }
               }
               _A[einfo._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
               einfo._var_name.clear();
               einfo._hsq_name.clear();
               for (int i = 0; i < 1; i++) {
                   stringstream strstrm;
                   strstrm << "";
                   einfo._var_name.push_back("V(O" + strstrm.str() + ")");
                   einfo._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
               }
               einfo._var_name.push_back("V(e)");
               reml_priors_var.clear();
               remlstatus=0;
               reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
               _A.clear();
               if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 )) einfo._varcmp_Py = _Vi;
           }
            if(einfo._varcmp_Py.size()>0)
            {
                
                _Vi = einfo._varcmp_Py;
                MatrixXd Vi_C;
                MatrixXd A;
                VectorXd Ct_vi_y;
                VectorXd y_buf=_y;
                if(!nopreadj_covar)
                {
                    y_buf=_y.array()-(_X*einfo._b).array();
                } else {
                    Vi_C = _Vi*_X; // n by m
                    A = _X.transpose()*Vi_C; // Ct_Vi_C : m by m
                    bool determinant_zero=false;
                    inverse_V(A, determinant_zero);
                    if(determinant_zero)
                    {
                        LOGPRINTF("The matrix is not invertible.\n");
                        TERMINATE();
                    }
                    Ct_vi_y = _X.transpose()*_Vi*y_buf;
                }
                #ifndef __APPLE__
                int n = (int)y_buf.size();
                float* y_mkl = new float[n];
                float* x_mkl = new float[n];
                float* Vi_mkl = new float[n*n];
                #pragma omp parallel for
                for (int i = 0; i < n; i++) y_mkl[i] = y_buf[i];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        Vi_mkl[i * n + j] = _Vi(i, j);
                    }
                }
                #endif
                VectorXd x_buf(einfo._eii_include.size());
                for(int i=0;i<leftids.size();i++)
                {
                    int id=leftids[i];
                    int targetChr=assoc_rlts[id].CHR;
                    int targetBP=assoc_rlts[id].BP;
                    string targetPrb=assoc_rlts[id].PROBE;
                    string gene=assoc_rlts[id].GENE;
                    char oren=assoc_rlts[id].OREN;
                    
                    iter=einfo._epi_map.find(assoc_rlts[id].PROBE);
                    int curidx = -9;
                    if(iter!=einfo._epi_map.end())
                    {
                         curidx=iter->second;
                    } else {
                        LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                        TERMINATE();
                    }
                    
                    for(int j=0; j<einfo._eii_include.size(); j++) x_buf(j)=einfo._val[curidx*einfo._eii_num+einfo._eii_include[j]];
                    if(!nopreadj_covar)
                        if(_X_c>1)
                            getResidual(x_buf, _X);
                    #ifndef __APPLE__
                    for(int j=0; j<einfo._eii_include.size(); j++) x_mkl[j]= x_buf(j);
                    #endif
                    double beta,se,pval;
                    if(!nopreadj_covar) {
                        #ifndef __APPLE__
                        mlma_cal_stat_mkl(y_mkl, x_mkl, Vi_mkl,n, beta, se, pval);
                        #else
                        mlma_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                        #endif
                    }
                    else mlma_cal_stat_covar(y_buf, x_buf, _Vi, Vi_C, A,Ct_vi_y, beta, se, pval);
                    
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
                #ifndef __APPLE__
                delete[] y_mkl;
                delete[] x_mkl;
                delete[] Vi_mkl;
                #endif
                
            } else {
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
        LOGPRINTF("Results of %ld probes have been saved in file %s.\n",wcount,filename.c_str());
        fclose(ofile);
        free_assoclist(assoc_rlts);
       
    }
    
    void moment_exact(char* outFileName, char* befileName,char* problstName,char* problst2exclde,char* genelistName, int chr,char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,char* probe2exclde,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_stepwise,int erm_alg, double swthresh, int tsk_ttl, int tsk_id, double swfdr, bool stepforwardonly, double sw_rsq)
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
        vector<MatrixXf> AN;
        mapids.resize(nrandcomp); // one for significant one for insignificant
        erm_prbs.resize(nrandcomp);
        AZERO.resize(nrandcomp);
        AN.resize(nrandcomp);
        
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
                    double chi=einfo._eii_include.size()* r2thresh / (1 - r2thresh) ;
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
            if(swfdr<0 && swthresh<0) swthresh = PFISHER/mapids[0].size();
            stepwise_slct(&einfo, mapids[0], slctids, rmsig, assoc_rlts,  swthresh, swfdr, true, stepforwardonly,sw_rsq);
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
                make_erm(&einfo,AZERO[j], erm_alg);
                AN[j] = einfo._grm_N;
                A0[rawORMc] = einfo._grm;
                einfo._epi_include=include_o;
                einfo._epi_map=snp_name_map_o;
                rawORMc++;
            }
        }
        A0[rawORMc]=MatrixXd::Identity(_n, _n);
        
        LOGPRINTF("\nPerforming the MOMENT analysis ...\n");
        LOGPRINTF("For each probe, the analysis will exclude probes in %d Kb region of centered at the probe to be tested.\n",expect_wind);
        int exwind=expect_wind*1000;
        
        string filename=string(outFileName)+".moment";
        if(tsk_ttl>1) filename=string(outFileName)+"_"+atos(tsk_ttl)+"_"+atos(tsk_id)+".moment";
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
        //#pragma omp parallel for private(remlstatus)
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
                        make_erm(&einfo,_A[aid], erm_alg);
                        einfo._grm_N = AN[j] - einfo._grm_N;
                        _A[aid] = AZERO[j] - _A[aid];
                        _A[aid] = _A[aid].array()/einfo._grm_N.cast<double>().array();
                        einfo._epi_include=include_o;
                        einfo._epi_map=snp_name_map_o;
                    } else {
                        _A[aid]=AZERO[j].array()/AN[j].cast<double>().array();
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

//
//  l3_gwas.cpp
//  osc
//
//  Created by Futao Zhang on 8/5/19.
//  Copyright © 2019 Futao Zhang. All rights reserved.
//

#include "l3_gwas.hpp"

namespace BFILE{
    void mlm_cal_stat(VectorXd &_Y, VectorXd &_x, MatrixXd &_Vi, double &beta, double &se, double &pval)
    {
        if(_Y.size()!=_x.size())
        {
            LOGPRINTF("Error: The size of Y and the size of X do not match.\n");
            TERMINATE();
        }
        if(_Y.size()!=_Vi.rows())
        {
            LOGPRINTF("Error: The size of Y and the rows of Vi do not match.\n");
            TERMINATE();
        }
        double Xt_Vi_X=0.0, chisq=0.0;
        uint64_t n=_Y.size();
        VectorXd Vi_X(n);
        Vi_X=_Vi*_x;
        Xt_Vi_X=Vi_X.dot(_x);
        se=1.0/Xt_Vi_X;
        beta=se*(Vi_X.dot(_Y));
        if(se>1.0e-30){
            se=sqrt(se);
            chisq=beta/se;
            pval=pchisq(chisq*chisq, 1);
        }
    }
    
    void mlm_cal_stat_covar(VectorXd &_Y,VectorXd &_x,  MatrixXd &_Vi,  MatrixXd &_X, double &beta, double &se, double &pval)
    {
        if(_Y.size()!=_x.size())
        {
            LOGPRINTF("Error: The size of Y and the size of X do not match.\n");
            TERMINATE();
        }
        if(_Y.size()!=_Vi.rows())
        {
            LOGPRINTF("Error: The size of Y and the rows of Vi do not match.\n");
            TERMINATE();
        }
        long _X_c = _X.cols();
        unsigned long col_num=_X_c+1;
        double chisq=0.0;
        uint64_t n=_Y.size();
        
        MatrixXd X(n,col_num);
        MatrixXd Vi_X(n,col_num);
        MatrixXd Xt_Vi_X(col_num,col_num);
        VectorXd Xt_Vi_y(col_num);
        VectorXd b_vec(col_num);
        
        X.block(0, 0, n, _X_c)=_X;
        X.col(_X_c) = _x;
        
        Vi_X=_Vi*X;
        Xt_Vi_X=X.transpose()*Vi_X;
        double logdt=0.0;
         #ifndef __APPLE__
        if(!comput_inverse_logdet_LU_mkl( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
        #else
        if(!comput_inverse_logdet_LU( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
        #endif
        Xt_Vi_y=Vi_X.transpose()*_Y;
        b_vec=Xt_Vi_X*Xt_Vi_y;
        
        se=Xt_Vi_X(_X_c,_X_c);
        beta=b_vec[_X_c];
        if(se>1.0e-30){
            se=sqrt(se);
            chisq=beta/se;
            pval=pchisq(chisq*chisq, 1);
        }
    }
    void mlm_cal_stat(VectorXd &_Y, bInfo* bdata, MatrixXd &_Vi, VectorXd &beta, VectorXd &se, VectorXd &pval)
    {
        int max_block_size = 10000;
        int m = (int)bdata->_include.size();
        int n = (int)bdata->_keep.size();
        VectorXd X, Vi_X;
        if(_Y.size()!=n)
        {
            LOGPRINTF("Error: The size of Y and the size of X do not match.\n");
            TERMINATE();
        }
        if(_Y.size()!=_Vi.rows())
        {
            LOGPRINTF("Error: The size of Y and the rows of Vi do not match.\n");
            TERMINATE();
        }
        beta.resize(m);
        se=VectorXd::Zero(m);
        pval=VectorXd::Constant(m,2);
        cout<<"\nRunning association tests for "<<m<<" SNPs ..."<<endl;
        int new_start = 0, block_size = 0, block_col = 0, k = 0, l = 0;
        MatrixXd X_block;
        vector<uint32_t> indx;
         double Xt_Vi_X;
        #ifndef __APPLE__
        float * Vi_mkl = new float[n*n];
        float *X_mkl=new float[n];
        float *Vi_X_mkl=new float[n];
        float *y = new float[n];
        for(int i=0; i<n; i++) y[i] = _Y[i];
        #pragma omp parallel for
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++) Vi_mkl[i*n+j]=_Vi(i,j);
        }
        #endif
        for(int i = 0; i < m; i++, block_col++)
        {
            // get a block of SNPs
            if(i == new_start){
                block_col = 0;
                new_start = i + max_block_size;
                block_size = max_block_size;
                if(new_start > m) block_size = m - i;
                indx.resize(block_size);
                for(k = i, l = 0; l < block_size; k++, l++) indx[l] = k;
                make_XMat_subset(bdata, indx,X_block,  false);
            }
            
            #ifndef __APPLE__
            for(int j = 0; j < n; j++) X_mkl[j] = X_block(j, block_col);
            cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, Vi_mkl, n, X_mkl, 1, 0.0, Vi_X_mkl, 1);
            Xt_Vi_X=cblas_sdot(n, X_mkl, 1, Vi_X_mkl, 1);
            se[i]=1.0/Xt_Vi_X;
            beta[i]=se[i]*cblas_sdot(n, y, 1, Vi_X_mkl, 1);
            if(se[i]>1.0e-30){
                se[i]=sqrt(se[i]);
                double chisq=beta[i]/se[i];
                pval[i]=pchisq(chisq*chisq, 1);
            }
            #else
            X = X_block.col(block_col);
            Vi_X=_Vi*X;
            Xt_Vi_X=Vi_X.dot(X);
            se(i)=1.0/Xt_Vi_X;
            beta(i)=se(i)*(Vi_X.dot(_Y));
            if(se(i)>1.0e-30){
                se(i)=sqrt(se(i));
                double chisq=beta(i)/se(i);
                pval(i)=pchisq(chisq*chisq, 1);
            }
            #endif
        }
        #ifndef __APPLE__
        delete[] X_mkl;
        delete[] Vi_X_mkl;
        delete[] Vi_mkl;
        delete[] y;
        #endif
    }
    
    void mlm_cal_stat_covar(VectorXd &_Y,bInfo* bdata,  MatrixXd &_Vi,  MatrixXd &_COV, VectorXd &beta, VectorXd &se, VectorXd &pval)
    {
        
        int max_block_size = 10000;
        long _X_c = _COV.cols();
        unsigned long col_num=_X_c+1;
        double chisq=0.0;
        int m = (int)bdata->_include.size();
        int n = (int)bdata->_keep.size();
        if(_Y.size()!=n)
        {
            LOGPRINTF("Error: The size of Y and the size of X do not match.\n");
            TERMINATE();
        }
        if(_Y.size()!=_Vi.rows())
        {
            LOGPRINTF("Error: The size of Y and the rows of Vi do not match.\n");
            TERMINATE();
        }
        
        MatrixXd X(n,col_num);
        MatrixXd Vi_X(n,col_num);
        MatrixXd Xt_Vi_X(col_num,col_num);
        VectorXd Xt_Vi_y(col_num);
        VectorXd b_vec(col_num);
        X.block(0, 0, n, _X_c)=_COV;
        
        beta.resize(m);
        se=VectorXd::Zero(m);
        pval=VectorXd::Constant(m,2);
        
        cout<<"\nRunning association tests for "<<m<<" SNPs ..."<<endl;
        int new_start = 0, block_size = 0, block_col = 0, k = 0, l = 0;
        MatrixXd X_block;
        vector<uint32_t> indx;
        for(int i = 0; i < m; i++, block_col++)
        {
            // get a block of SNPs
            if(i == new_start)
            {
                block_col = 0;
                new_start = i + max_block_size;
                block_size = max_block_size;
                if(new_start > m) block_size = m - i;
                indx.resize(block_size);
                for(k = i, l = 0; l < block_size; k++, l++) indx[l] = k;
                make_XMat_subset(bdata, indx,X_block,  false);
            }
            X.col(_X_c)=X_block.col(block_col);
            Vi_X=_Vi*X;
            Xt_Vi_X=X.transpose()*Vi_X;
            double logdt=0.0;
             #ifndef __APPLE__
            if(!comput_inverse_logdet_LU_mkl( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
            #else
             if(!comput_inverse_logdet_LU( Xt_Vi_X, logdt)) throw("Error: Xt_Vi_X is not invertable.");
            #endif
            Xt_Vi_y=Vi_X.transpose()*_Y;
            b_vec=Xt_Vi_X*Xt_Vi_y;
            se(i)=Xt_Vi_X(_X_c,_X_c);
            beta(i)=b_vec[_X_c];
            if(se(i)>1.0e-30){
                se(i)=sqrt(se(i));
                chisq=beta(i)/se(i);
                pval(i)=pchisq(chisq*chisq, 1);
            }
        }
    }
    
    void cast( eInfo* einfo, bInfo* bdata)
    {
        einfo->_eii_include = bdata->_keep;
        einfo->_r_indx = bdata->_r_indx;
        einfo->_reml_mtd = bdata->_reml_mtd;
        einfo->_reml_max_iter = bdata->_reml_max_iter;
        einfo->_reml_fixed_var = bdata->_reml_fixed_var;
        einfo->_reml_force_inv = bdata->_reml_force_inv;
        einfo->_bivar_reml = false;
        einfo->_y_Ssq =bdata->_y_Ssq;
        einfo->_Asp = bdata->_Asp;
        einfo->_fixed_rg_val = bdata->_fixed_rg_val;
        einfo->_var_name = bdata->_var_name;
        einfo->_within_family = bdata->_within_family;
        einfo->_V_inv_mtd = bdata->_V_inv_mtd;
        einfo->_P = bdata->_P;
        einfo->_reml_force_converge = bdata->_reml_force_converge;
        einfo->_reml_AI_not_invertible = bdata->_reml_AI_not_invertible;
        einfo->_varcmp = bdata->_varcmp;
        einfo->_hsq_name = bdata->_hsq_name;
        einfo->_b = bdata->_b;
        einfo->_se = bdata->_se;
    }
    void moment_gwas(char* outFileName, char* bFileName,double maf,char* snplstName,char* snplst2exclde,int chr,char* indilstName,char* indilst2remove, char* phenofileName,char* mpheno, char* covfileName,char* qcovfileName, bool within_family,char* priors,char* priors_var, bool no_constrain,int reml_mtd,int MaxIter,bool reml_fixed_var_flag,bool reml_force_inv_fac_flag, bool reml_force_converge_flag, bool  reml_no_converge_flag, bool nopreadj_covar,  double lambda_wind, bool fastlinear, bool force_mlm,double bcthresh,int expect_wind, int expect_num, int expect_pcs, int nrandcomp,int slctmtd, double r2thresh, double percent, bool approximate_flag, bool approximate_stepwise,int grm_alg, double swthresh, double swfdr, bool swlogit, bool swforward, double swrsq, char* grm_file,bool grm_bin_flag)
    {
        setNbThreads(thread_num);
        if(bFileName==NULL)
        {
            LOGPRINTF("ERROR: please input the genotype data by the option --bfile.\n");
            TERMINATE();
        }
        if(phenofileName ==NULL)
        {
            LOGPRINTF("ERROR: please input the phenotype data by the option --pheno.\n");
            TERMINATE();
        }
        
        bInfo bdata;
        init_binfo(&bdata);
        bdata._reml_mtd=reml_mtd;
        bdata._reml_max_iter=MaxIter;
        bdata._reml_fixed_var=reml_fixed_var_flag;
        bdata._reml_force_inv=reml_force_inv_fac_flag;
        bdata._reml_force_converge=reml_force_converge_flag;
        bdata._reml_no_converge=reml_no_converge_flag;
        
        read_famfile(&bdata, string(bFileName)+".fam");
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(chr>0) extract_snp(&bdata, chr);
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        
        
        vector<double> reml_priors;
        vector<double> reml_priors_var;
        vector<string> vs_buf;
        vector<string> grm_id;
        vector<string> grm_files;
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
        
        bool succ=1;
        for(int i=0;i<bdata._include.size();i++)
        {
            if(bdata._chr[bdata._include[i]]==-9 ||  bdata._bp[bdata._include[i]]==-9)
            {
                succ=0;
                LOGPRINTF("ERROR: NA chromosome or NA snp position found.\n");
                TERMINATE();
            }
        }
      
        if(phenofileName !=NULL) read_phen(&bdata, phenofileName, mpheno,false);
        if(covfileName != NULL) read_cov(&bdata, covfileName, false);
        if(qcovfileName != NULL) read_cov(&bdata, qcovfileName, true);
        if(grm_file!=NULL){
            grm_files.push_back(grm_file);
            read_grm(&bdata, grm_bin_flag,  grm_file, grm_id, false,  false);
            AN[nrandcomp-1] = bdata._grm_N;
            AZERO[nrandcomp-1] = bdata._grm.array()*bdata._grm_N.cast<double>().array();
        }
        read_bedfile(&bdata, string(bFileName)+".bed");
        if(bdata._mu.empty()) calcu_mu(&bdata);
        if (maf >=0 ) filter_snp_maf(&bdata, maf);
        
        map<string, int>::iterator iter;
        
        int _n=(int)bdata._keep.size();
        if(_n<1)
        {
            LOGPRINTF("ERROR: no individual is in common among the input files.\n");
            TERMINATE();
        }
        VectorXd _y;
        _y.setZero(_n);
        for(int i=0; i<_n; i++){
            _y(i)=bdata._pheno[bdata._keep[i]];
        }
        
        // construct X matrix
        vector<MatrixXd> E_float;
        MatrixXd qE_float;
        MatrixXd _X;
        int _X_c;
        _X_c=construct_X(&bdata, E_float, qE_float,_X); //_X has the dimension of [einfo._eii_include.size(),_X_c]
        
        
        
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
                    bdata._pheno[ bdata._keep[i]]= 1;
                    _y(i)=1;
                }
                else if(flag==3)
                {
                     bdata._pheno[ bdata._keep[i]]-=1;
                    _y(i)=-1;
                }
            }
        }
        vector<ASSOCRLT> assoc_rlts;
        MatrixXd COV_plus;
        if(swlogit && !qtrait) testLogit(assoc_rlts, NULL, &bdata,COV_plus);
        else testLinear(assoc_rlts, NULL, &bdata,COV_plus);
        
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
                    erm_prbs[0].push_back(assoc_rlts[k].SNP);
                }
                for(int k=(int)(assoc_rlts.size()-number-1);k>=0;k--)
                {
                    mapids[1].push_back(k);
                    erm_prbs[1].push_back(assoc_rlts[k].SNP);
                }
                LOGPRINTF("%ld probes incldued in the first component and %ld in the second component.\n",mapids[0].size(),mapids[1].size());
            }
            else
            {
                if(r2thresh<0 && bcthresh<0) bcthresh=0.05/assoc_rlts.size();
                else if(r2thresh > 0) {
                    double chi=bdata._keep.size()* r2thresh / (1 - r2thresh) ;
                    bcthresh = pchisq(chi, 1);
                }
                for(int k=(int)(assoc_rlts.size()-1);k>=0;k--)
                {
                    if(assoc_rlts[k].PVAL<=bcthresh) {
                        mapids[0].push_back(k); //sig
                        erm_prbs[0].push_back(assoc_rlts[k].SNP);
                    }
                    else {
                        mapids[1].push_back(k);
                        erm_prbs[1].push_back(assoc_rlts[k].SNP);
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
                    erm_prbs[0].push_back(assoc_rlts[k].SNP);
                }
                else if(assoc_rlts[k].PVAL<=secondthresh)
                {
                    mapids[1].push_back(k);
                    erm_prbs[1].push_back(assoc_rlts[k].SNP);
                }
                else
                {
                    mapids[2].push_back(k);
                    erm_prbs[2].push_back(assoc_rlts[k].SNP);
                }
            }
            LOGPRINTF("%ld probes included in the 1st random component.\n",mapids[0].size());
            LOGPRINTF("%ld probes included in the 2nd random component.\n",mapids[1].size());
            LOGPRINTF("%ld probes included in the 3rd random component.\n",mapids[2].size());
        }
        
         double RSQVO =-9;
        VectorXd mlma_b,mlma_se;
        vector<double> mlma_varcmp;
        int rawORMc=0;
        vector<int> include_o(bdata._include);
        map<string, int> snp_name_map_o(bdata._snp_name_map);
        
        if( approximate_stepwise && mapids[0].size()>0)
        {
            for(int i=0; i < 2; i++) bdata._r_indx.push_back(i);
            if(grm_file==NULL) {
                make_grm( &bdata,AZERO[nrandcomp-1], grm_alg);
                AN[nrandcomp-1]=bdata._grm_N;
            }
            bdata._geno.resize(0,0);
                _A.resize(bdata._r_indx.size());
                (_A[0]).resize(_n, _n);
                #pragma omp parallel for
                for(int i=0; i<_n; i++)
                {
                    for(int j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=bdata._grm(i,j);
                }
                _A[bdata._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
                bdata._var_name.clear();
                bdata._hsq_name.clear();
                for (int i = 0; i < 1; i++) {
                    stringstream strstrm;
                    strstrm << "";
                    bdata._var_name.push_back("V(O" + strstrm.str() + ")");
                    bdata._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
                }
            
                bdata._var_name.push_back("V(e)");
                remlstatus=0;
                reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,_X, _y,_A, _Vi,  reml_mtd,  MaxIter,mlma_b,mlma_se, mlma_varcmp);
                _A.clear();
                bdata._r_indx.clear();
            if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ))
            {
                RSQVO=mlma_varcmp[0]/(mlma_varcmp[0] + mlma_varcmp[1]);
                bdata._varcmp_Py = _Vi;
            } else RSQVO=1;
                swrsq *= RSQVO;
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
                    bdata._pheno[bdata._keep[i]]=y_buf[i];
                }
                //if(loud) {LOGPRINTF("Phenotype reset.\n");}
            } else
            {
                double ymean=_y.mean();
                VectorXd y_buf =_y.array() - ymean;
                for(int i=0; i<_n; i++){
                    bdata._pheno[bdata._keep[i]]=y_buf[i];
                }
            }
            stepwise_slct(&bdata,mapids[0], slctids, rmsig, assoc_rlts,  swthresh,swfdr, true,swforward, swrsq);
            erm_prbs[0].clear();
            for(int j=0;j<mapids[0].size();j++)
            {
                erm_prbs[0].push_back(assoc_rlts[mapids[0][j]].SNP);
            }
            for(int j=0;j<rmsig.size();j++)
            {
                mapids[1].push_back(rmsig[j]);
                erm_prbs[1].push_back(assoc_rlts[rmsig[j]].SNP);
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
                    bdata._pheno[bdata._keep[i]]=_y[i];
                }
                //if(loud) {LOGPRINTF("Phenotype restored.\n");}
            
        }
        
        for(int i=0;i<mapids.size();i++) if(mapids[i].size()) rawORMc++;
        A0.resize(rawORMc+1);
        rawORMc=0;
        for(int j=0;j<erm_prbs.size()-1;j++)
        {
            if(erm_prbs[j].size()>0)
            {
                update_map_kp(erm_prbs[j], bdata._snp_name_map, bdata._include);
                make_grm(&bdata,AZERO[j],grm_alg);
                bdata._geno.resize(0,0);
                AN[j] = bdata._grm_N;
                A0[rawORMc] = bdata._grm;
                bdata._include=include_o;
                bdata._snp_name_map=snp_name_map_o;
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
        
        LOGPRINTF("\nPerforming the MOMENT analysis ...\n");
        LOGPRINTF("For each probe, the analysis will exclude probes in %d Kb region of centered at the probe to be tested.\n",expect_wind);
        int exwind=expect_wind*1000;
        
        string filename=string(outFileName)+".moment";
        FILE* ofile = fopen(filename.c_str(), "w");
        if (!(ofile)) {
            LOGPRINTF("ERROR: open error %s\n", filename.c_str());
            TERMINATE();
        }
        string outstr="Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\tp\n";
        if(fputs_checked(outstr.c_str(),ofile))
        {
            LOGPRINTF("ERROR: error in writing file %s .\n", filename.c_str());
            TERMINATE();
        }
        
        long wcount=0;
        int outid=0;
        if(expect_num<0) expect_num = (int)mapids[outid].size();
        long num2exp=expect_num>mapids[outid].size()?mapids[outid].size():expect_num;
        
        double cr=0;
        vector<string> freml;
        
        for(int i=0;i<num2exp;i++)
        {
            double desti=1.0*wcount/(bdata._include.size()-1);
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
            string target=assoc_rlts[id].SNP;
            string targetA1=assoc_rlts[id].a1;
            string targetA2=assoc_rlts[id].a2;
            double targetFreq=assoc_rlts[id].freq;
        
            iter=bdata._snp_name_map.find(assoc_rlts[id].SNP);
            VectorXd x_buf(bdata._keep.size());
            if(iter!=bdata._snp_name_map.end())
            {
                int curidx=iter->second;
                make_Xvec_subset(&bdata, curidx, x_buf, false);
            } else {
                LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                TERMINATE();
            }
          
            for(int jj=0;jj<nrandcomp;jj++) erm_prbs[jj].clear();
            int norms=0,excount=0;
            for(int jj=0;jj<nrandcomp;jj++)
            {
                for(int j=0;j<mapids[jj].size();j++)
                {
                    int tid=mapids[jj][j];
                    if(targetChr == assoc_rlts[tid].CHR && abs(targetBP-assoc_rlts[tid].BP)<=exwind) {
                        erm_prbs[jj].push_back(assoc_rlts[tid].SNP);
                        excount++;
                    }
                }
            }
            
            //LOGPRINTF("%d probe(s) excluded from calculating the ORM(s) for the probe %s.\n", excount,assoc_rlts[id].PROBE);
            for(int jj=0;jj<nrandcomp;jj++) if(mapids[jj].size()-erm_prbs[jj].size()>0) norms++;
            
            bdata._r_indx.clear();
            for(int j=0; j < norms + 1; j++) bdata._r_indx.push_back(j);
            _A.resize(bdata._r_indx.size());
            int aid=0;
            for(int j=0; j < erm_prbs.size(); j++)
            {
                long N=mapids[j].size()-erm_prbs[j].size();
                if(N)
                {
                    if(erm_prbs[j].size()>0)
                    {
                        update_map_kp(erm_prbs[j], bdata._snp_name_map, bdata._include);
                        make_grm(&bdata,_A[aid],grm_alg);
                        bdata._geno.resize(0,0);
                        bdata._grm_N = AN[j] - bdata._grm_N;
                        _A[aid] = AZERO[j] - _A[aid];
                        //_A[aid] = _A[aid].array()/bdata._grm_N.cast<double>().array(); // divided by 0
                        #pragma omp parallel for
                        for (int ii = 0; ii < _n; ii++) {
                            for (int jj = 0; jj <= ii; jj++) {
                                if(bdata._grm_N(ii,jj) > 0) _A[aid](ii,jj) /= bdata._grm_N(ii,jj);
                                else _A[aid](ii,jj) = 0.0;
                                _A[aid](jj,ii) = _A[aid](ii,jj);
                            }
                        }
                        bdata._include=include_o;
                        bdata._snp_name_map=snp_name_map_o;
                    } else {
                        //_A[aid]=AZERO[j].array()/AN[j].cast<double>().array();
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
                    aid++;
                }
            }
            _A[bdata._r_indx.size()-1]=MatrixXd::Identity(_n, _n);
            
            // names of variance component
            bdata._var_name.clear();
            bdata._hsq_name.clear();
            for (int j = 0; j < norms; j++) {
                stringstream strstrm;
                if (norms == 1) strstrm << "";
                else strstrm << j + 1;
                bdata._var_name.push_back("V(O" + strstrm.str() + ")");
                bdata._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
            }
            bdata._var_name.push_back("V(e)");
            
            bdata._within_family=within_family;
            if(within_family) detect_family(&bdata, _A);
           
            VectorXd _b,_se;
            vector<double> _varcmp;
            remlstatus=0; //reset reml status
            //cast(&einfo, &bdata);
            //if(norms) reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,_A,_Vi,outFileName);
            reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c,_X, _y,_A, _Vi,  reml_mtd,  MaxIter,_b,_se,_varcmp);
            double beta=0, se=-9, pval=1;
            if(norms && (remlstatus==0  || remlstatus==-5 || remlstatus==-3 ))
            {
                double rox=0;
                for( int ii=0;ii<_varcmp.size()-1;ii++)
                rox +=_varcmp[ii];
                rox = rox / (rox + _varcmp[_varcmp.size()-1]);
                if(rox/RSQVO >= 0.1)
                {
                    bdata._P.resize(0,0);
                    _A.clear();
                    VectorXd y_buf=_y;
                    if(!nopreadj_covar)
                    y_buf=_y.array()-(_X*_b).array();
                    
                    if(!nopreadj_covar) mlm_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                    else mlm_cal_stat_covar(y_buf, x_buf, _Vi, _X, beta, se, pval);
                    
                    string chrstr;
                    if(targetChr==23) chrstr="X";
                    else if(targetChr==24) chrstr="Y";
                    else chrstr=atosm(targetChr);
                    outstr = chrstr + '\t' + target + '\t' + atosm(targetBP)  + '\t'+ targetA1  + '\t' +targetA2  + '\t' +atosm(targetFreq)  + '\t'  + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
                    if(fputs_checked(outstr.c_str(),ofile))
                    {
                        LOGPRINTF("ERROR: in writing file %s .\n", filename.c_str());
                        TERMINATE();
                    }
                    wcount++;
                } else {
                    if(loud) printf("%s needs recalcualtion.\n", target.c_str());
                    freml.push_back(target);
                }
            }
            else
            {
                freml.push_back(target);
            }
        }
        
        if(loud)
        {
            LOGPRINTF("%ld probes in the first componet failed in REML.\n", freml.size());
        }
        long leftnum = mapids[outid].size()-num2exp+freml.size();
        for(int outid=1;outid<mapids.size();outid++) leftnum += mapids[outid].size();
        vector<string> leftrss(leftnum);
        long count = 0;
        for(int i=0;i<freml.size();i++) {
            leftrss[count++] = freml[i];
        }
        if(mapids[outid].size()>num2exp)
            for(long i=num2exp;i<mapids[outid].size();i++)
            {
                int id=mapids[outid][i];
                leftrss[count] = assoc_rlts[id].SNP;
                count++;
            }
        for(int outid=1;outid<mapids.size();outid++)
        {
            for(int i=0;i<mapids[outid].size();i++)
            {
                int id=mapids[outid][i];
                leftrss[count] = assoc_rlts[id].SNP;
                count++;
            }
        }
        
        bdata._r_indx.clear();
        for(int j=0; j < rawORMc + 1; j++) bdata._r_indx.push_back(j);
        // names of variance component
        bdata._var_name.clear();
        bdata._hsq_name.clear();
        for (int j = 0; j < rawORMc; j++) {
            stringstream strstrm;
            if (rawORMc == 1) strstrm << "";
            else strstrm << j + 1;
            bdata._var_name.push_back("V(O" + strstrm.str() + ")");
            bdata._hsq_name.push_back("V(O" + strstrm.str() + ")/Vp");
        }
        bdata._var_name.push_back("V(e)");
        
        bdata._within_family=within_family;
        if(within_family) detect_family(&bdata, _A);
        VectorXd _b,_se;
        remlstatus=0; //reset reml status
        remloasi = false;
        reml_priors_var.clear();
        vector<double> _varcmp;
        loud = true;
        //cast(&einfo, &bdata);
        //reml(&einfo, false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true, _X_c, _X,_y,A0,_Vi,outFileName);
        reml( false, true, reml_priors, reml_priors_var,  no_constrain,  _X_c, _X, _y,A0, _Vi,  reml_mtd,  MaxIter,_b,_se,_varcmp);
        bool rsqck=false;
        if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ))
        {
            double rox=0;
            for( int ii=0;ii<_varcmp.size()-1;ii++)
            rox +=_varcmp[ii];
            rox = rox / (rox + _varcmp[_varcmp.size()-1]);
            if(loud) printf("%f:%f\n",rox,RSQVO);
            if(rox/RSQVO >= 0.1) rsqck=true;
        }
        if( (remlstatus==0 || remlstatus==-3  || remlstatus==-5 ) && rsqck)
        {
            VectorXd beta, se, pval;
            VectorXd y_buf=_y;
            if(!nopreadj_covar)
                y_buf=_y.array()-(_X*_b).array();
            update_map_kp(leftrss, bdata._snp_name_map, bdata._include);
            if(!nopreadj_covar) mlm_cal_stat(y_buf, &bdata, _Vi, beta, se, pval);
            else mlm_cal_stat_covar(y_buf, &bdata, _Vi, _X, beta, se, pval);
            
            for(int i=0;i<bdata._include.size();i++)
            {
                
                int id=bdata._include[i];
                int targetChr=bdata._chr[id];
                int targetBP=bdata._bp[id];
                string target=bdata._snp_name[id];
                string targetA1=bdata._allele1[id];
                string targetA2=bdata._allele2[id];
                double targetFreq=bdata._mu[id]/2;
                string chrstr;
                if(targetChr==23) chrstr="X";
                else if(targetChr==24) chrstr="Y";
                else chrstr=atosm(targetChr);
                outstr = chrstr + '\t' + target + '\t' + atosm(targetBP)  + '\t' + targetA1  + '\t' +targetA2  + '\t' +atosm(targetFreq)  + '\t' +atos(beta(i)) + '\t' + atos(se(i)) + '\t' + dtos(pval(i))  +'\n';
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
            
            if(loud) printf("Degrade to MLMA\n");
            if(bdata._varcmp_Py.size()>0)
            {
                _Vi = bdata._varcmp_Py;
                VectorXd y_buf=_y;
                if(!nopreadj_covar)
                    y_buf=_y.array()-(_X*mlma_b).array();
                
                VectorXd x_buf(bdata._keep.size());
                for(int i=0;i<leftrss.size();i++)
                {
                    string target = leftrss[i];
                    int curidx = -9;
                    iter=bdata._snp_name_map.find(leftrss[i]);
                    if(iter!=bdata._snp_name_map.end())
                    {
                        curidx=iter->second;
                        make_Xvec_subset(&bdata, curidx, x_buf, false);
                    } else {
                        LOGPRINTF("ERROR: bugs found in MOMENT(). please report.\n");
                        TERMINATE();
                    }
                    int targetChr=bdata._chr[curidx];
                    int targetBP=bdata._bp[curidx];
                    string targetA1=bdata._allele1[curidx];
                    string targetA2=bdata._allele2[curidx];
                    double targetFreq=bdata._mu[curidx]/2;
                    
                    double beta,se,pval;
                    if(!nopreadj_covar) mlm_cal_stat(y_buf, x_buf, _Vi, beta, se, pval);
                    else mlm_cal_stat_covar(y_buf, x_buf, _Vi, _X, beta, se, pval);
                    
                    string chrstr;
                    if(targetChr==23) chrstr="X";
                    else if(targetChr==24) chrstr="Y";
                    else chrstr=atosm(targetChr);
                    outstr = chrstr + '\t' + target + '\t' + atosm(targetBP)  + '\t'+ targetA1  + '\t' +targetA2  + '\t' +atosm(targetFreq)  + '\t'  + atos(beta) + '\t' + atos(se) + '\t' + dtos(pval)  +'\n';
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
      
        }
             
        LOGPRINTF("Results of %ld probes have been saved in file %s.\n",wcount,filename.c_str());
        fclose(ofile);
        free_assoclist(assoc_rlts);
 
    }
    
}

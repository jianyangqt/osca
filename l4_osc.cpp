//
//  l4_osc.cpp
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l4_osc.h"

using namespace EFILE;

int main(int argc, char * argv[])
{
    /*
     MatrixXd A = MatrixXd::Random(3,3);
     cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
     EigenSolver<MatrixXd> es(A);
     cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
     cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
     complex<double> lambda = es.eigenvalues()[0];
     cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
     VectorXcd v = es.eigenvectors().col(0);
     cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
     cout << "... and A * v = " << endl << A.cast<complex<double> >() * v << endl << endl;
     MatrixXcd D = es.eigenvalues().asDiagonal();
     MatrixXcd V = es.eigenvectors();
     cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
     */
    
    /*
    MatrixXf m = MatrixXf::Random(3,8);
    cout << "Here is the matrix m:" << endl << m << endl;
    JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
    Vector3f rhs(1, 0, 0);
    cout << "Now consider this rhs vector:" << endl << rhs << endl;
    cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;

    MatrixXf ff=svd.singularValues().asDiagonal();
    cout<<ff<<endl;
    */
    cout << "*******************************************************************" << endl;
    cout << "* Omics-data-based Complex Trait Analysis (OCTA)" << endl;
    cout << "* version 0.30" << endl;
    cout << "* (C) 2016 Futao Zhang, Zhihong Zhu and Jian Yang" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;
    

    char tmpname[4]="osc";
    char* outname=tmpname;
    string logfname=string(outname)+".log";
    logfile=fopen(logfname.c_str(),"w");
    if (!logfile) {
        printf("Error: Failed to open %s.  Try ", logfname.c_str());
        exit(1);
    }
    printf("Logging to %s.\n", logfname.c_str());
    
    
    /**************************/
    /*
    MatrixXf m1 = MatrixXf::Random(2,6);
    srand((unsigned int) time(0));
    MatrixXf m2 = MatrixXf::Random(2,6);
    cout<<m1<<endl;
    cout<<m2<<endl;
    m1=m1-m2;
    cout<<m1<<endl;
    m1=m1.array()*m1.array();
    cout<<m1<<endl;
    VectorXf n1=m1.colwise().sum();
    cout<<n1<<endl;
    
    vector<double> a;
    for(int i=0;i<n1.size();i++) a.push_back(n1[i]);
    vector<int> b;
    get_bottom_indices(b,a,10);
   */
    /***************************/
    
    
    FLAGS_VALID_CK(argc, argv);
    
    long int time_used = 0, start = time(NULL);
    time_t curr = time(0);
    LOGPRINTF("Analysis started: %s \n", ctime(&curr));
    
    unsigned char* wkspace_ua;
    uint64_t mb=2048;
    wkspace_ua = (unsigned char*)malloc(mb * 1048576 * sizeof(char));
    memset(wkspace_ua, 0, mb * 1048576 * sizeof(char));
    uint64_t llxx=getMemSize_Plink();
    mem_left=getAllocMB_Plink(llxx);
    if (llxx) {
        sprintf(logbuf, "%llu MB RAM detected; reserving %lld MB for main workspace.\n", llxx, mem_left);
        logprintb();
    }
    
    LOGPRINTF("\nOptions:\n");
    option(argc, argv);
    
   
    curr = time(0);
    time_used = time(NULL) - start;
    LOGPRINTF("\nAnalysis finished: %s \nComputational time: %ld:%ld:%ld\n", ctime(&curr),time_used / 3600,(time_used % 3600) / 60,time_used % 60);
    
    if(logfile) fclose(logfile);
    return 0;

}

void option(int option_num, char* option_str[])
{
    
    char* efileName=NULL;
    char* befileName=NULL;
    char* befileName2=NULL;
    char* outfileName=NULL;
    char* indilstName=NULL;
    char* snplstName=NULL;
    char* indilst2remove=NULL;
    char* snplst2exclde=NULL;
    char* problstName=NULL;
    char* problst2exclde=NULL;
    char* genelistName=NULL;
    char* prbname=NULL;
    char* probe2rm=NULL;
    char* fromprbname=NULL;
    char* toprbname=NULL;
    char* genename=NULL;
    int fromprbkb=-9;
    int toprbkb=-9;
    char* phenofileName=NULL;
    char* mpheno=NULL;
    char* befileFlstName=NULL;
    char* epifname=NULL;
    
    int prbWind=1000; //Kb
   

    int chr=-9;
    double maf=0.0;
    double std_thresh = 0.0;

    bool no_fid_flag=false;
    bool make_beed_flag=false;
    bool make_efile_flag=false;
    bool prbwindFlag=false;
    bool merge_beed_flag=false;

    
    int efileType=GENEEXPRESSION;
    int valueType=VALUE;
    bool transposedin=false;
    bool transposedout=false;
    
    // for ERM
    bool make_erm_flag=false;
    bool erm_bin_flag=true;
    int erm_alg=0;
    char* erm_file=NULL; // can be an erm file name or an erm filelist file name
    char* subtract_erm_file=NULL;
    bool m_erm_flag=false;
    char* priors= NULL;
    char* priors_var= NULL;
    bool reml_fixed_var_flag=false;
    double erm_cutoff = -2.0;
    //for MLMA
    bool mlma_flag=false;
    bool mlma_loco_flag=false;
    char* qcovfileName=NULL;
    char* covfileName=NULL;
    bool within_family = false;
    int reml_mtd=0;
    bool no_constrain=false;
    int  MaxIter = 100;
    bool reml_force_inv_fac_flag = false;
    bool reml_force_converge_flag = false;
    bool  reml_no_converge_flag = false;
    bool mlma_no_adj_covar = false;
    
    bool ignore_Ce=false; // bivar
    bool bivar_no_constrain = false;
    
    bool m2beta = false;
    bool beta2m = false;
    bool diffflag =false;
    bool update_epi_file_flag = false;
    bool pca_flag = false;
    bool refacotr_flag = false;
    int out_pc_num = 20;
    int celltype_num = -9;
    int dmr_num = -9;
    int autosome_num = 22;
    //for assoc
    bool assoc_flag = false;
    bool linear_flag = false;
    bool simu_qt_flag = false, simu_cc = false;
    int simu_rep = 1, simu_case_num = 0, simu_control_num = 0, simu_eff_mod = 1;  // 0 standarise the probes. 1 use raw probe profiles.
    char* simu_causal = NULL;
    double simu_h2 = 0.1, simu_K = 0.1, simu_seed = -rand_seed();
    
    bool reml_flag = false, pred_rand_eff = false, est_fix_eff = false, no_lrt = false;
    double prevalence = -2.0, prevalence2 = -2.0;
    vector<int> reml_drop;
    reml_drop.push_back(1);

    int inde_num=0;
    bool ext_inde_flag=false;
    double ldrsq=1.0;
    
    double percentage_out = 0.0; // 0: leave top 0% of probes of uncorrected association out from the ERM calculation.  1: leave top 100% of probes out from the ERM calculation.
    double upperBeta=1;
    double lowerBeta=0;
    
    char* dpvalfName=NULL;
    double thresh_det_pval=0.05;
    double thresh_prpt_prb=0.01;
    double thresh_prpt_spl=0.01;
    int filter_det_pval_mth=0;
    
    double missing_ratio_prob=1.0; // missing value rate threshold to QC probes. 1.0 means filtering out no probe. 0.0 means filtering out probes even with 1 missing value.
    //
    bool estn_flag=false;
    
    bool getvariance_flag = false;
    bool getmean_flag =false;
    
    //
    bool blup_probe_flag = false;
    char* blup_indi_file =NULL;
    //
    bool score_flag = false;
    char* score_file = NULL;
    int col_prb = -9;
    int col_score = -9;
    bool impute_mean_flag = false;
    bool hasHeader = false;
    
    for(int i=0;i<option_num;i++)
    {
        if(0==strcmp(option_str[i],"--efile")){
            efileName=option_str[++i];
            FLAG_VALID_CK("--efile", efileName);
            FileExist(efileName);
            LOGPRINTF("--efile %s\n",efileName);
        }
        else if(0==strcmp(option_str[i],"--tefile")){
            efileName=option_str[++i];
            transposedin=true;
            FLAG_VALID_CK("--tefile", efileName);
            FileExist(efileName);
            LOGPRINTF("--tefile %s\n",efileName);
        }
        else if(0==strcmp(option_str[i],"--gene-expression")){
            efileType=GENEEXPRESSION;
            LOGPRINTF("--gene-expression\n");
        }
        else if(0==strcmp(option_str[i],"--methylation")){
            efileType=METHYLATION;
            valueType=BETAVALUE;
            LOGPRINTF("--methylation\n");
        }
        else if(0==strcmp(option_str[i],"--methylation-beta")){
            efileType=METHYLATION;
            valueType=BETAVALUE;
            LOGPRINTF("--methylation-beta\n");
        }
        else if(0==strcmp(option_str[i],"--methylation-m")){
            efileType=METHYLATION;
            valueType=MVALUE;
            LOGPRINTF("--methylation-m\n");
        }
        else if(0==strcmp(option_str[i],"--make-beed")){
            make_beed_flag=true;
            LOGPRINTF("--make-beed\n");
        }
        else if(0==strcmp(option_str[i],"--make-efile")){
            make_efile_flag=true;
            LOGPRINTF("--make-efile\n");
        }
        else if(0==strcmp(option_str[i],"--make-tefile")){
            make_efile_flag=true;
            transposedout=true;
            LOGPRINTF("--make-tefile\n");
        }
        else if (0 == strcmp(option_str[i], "--out")){
            outfileName =option_str[++i];;
            if(outfileName !=NULL && has_prefix(outfileName, "--"))
            {
                outfileName=NULL;
                i--;
            }
            LOGPRINTF("--out %s\n", outfileName);
        }
        else if(0==strcmp(option_str[i],"--befile")){
            if(befileName==NULL) {
                befileName=option_str[++i];
                FLAG_VALID_CK("--befile", befileName);
                LOGPRINTF("--befile %s\n",befileName);
            } else {
                befileName2=option_str[++i];
                FLAG_VALID_CK("--befile", befileName2);
                LOGPRINTF("--befile %s\n",befileName2);
            }
           
        }
        else if(0==strcmp(option_str[i],"--befile-flist")){
            befileFlstName=option_str[++i];
            merge_beed_flag=true;
            FLAG_VALID_CK("--befile-flist", befileFlstName);
            LOGPRINTF("--befile-flist %s\n",befileFlstName);
        }
        else if(strcmp(option_str[i],"--keep")==0){
            indilstName=option_str[++i];
            FLAG_VALID_CK("--keep", indilstName);           
            FileExist(indilstName);
            LOGPRINTF("--keep %s\n", indilstName);
        }
        else if(strcmp(option_str[i],"--remove")==0){
            indilst2remove=option_str[++i];
            FLAG_VALID_CK("--remove", indilst2remove);
            FileExist(indilst2remove);
            LOGPRINTF("--remove %s\n",indilst2remove);
        }
        else if(strcmp(option_str[i],"--extract-snp")==0){
            snplstName=option_str[++i];
            FLAG_VALID_CK("--extract-snp", snplstName);
            FileExist(snplstName);
            LOGPRINTF("--extract-snp %s\n",snplstName);
        }
        else if(strcmp(option_str[i],"--extract-probe")==0){
            problstName=option_str[++i];
            FLAG_VALID_CK("--extract-probe", problstName);
            FileExist(problstName);
            LOGPRINTF("--extract-probe %s\n",problstName);
        }
        else if(strcmp(option_str[i],"--exclude-probe")==0){
            problst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-probe", problst2exclde);
            FileExist(problst2exclde);
            LOGPRINTF("--exclude-probe %s\n",problst2exclde);
        }
        else if(strcmp(option_str[i],"--maf")==0){
            maf=atof(option_str[++i]);
            LOGPRINTF("--maf %lf\n",maf);
            if(maf<0 || maf>0.5)
            {
                LOGPRINTF("Error: --maf should be within the range from 0 to 0.5.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--genes")==0){
            genelistName=option_str[++i];
            FLAG_VALID_CK("--genes", genelistName);
            LOGPRINTF("--genes %s\n",genelistName);
            FileExist(genelistName);
        }
        else if(0 == strcmp(option_str[i],"--thread-num")){
            thread_num=atoi(option_str[++i]);
            LOGPRINTF("--thread-num %d\n",thread_num);
        }
        else if(strcmp(option_str[i],"--chr")==0){
            char* tmpstr=option_str[++i];
            if(strncmp(tmpstr,"X",1)==0) chr=23;
            else if(strncmp(tmpstr,"Y",1)==0) chr=24;
            else chr=atoi(tmpstr);
            FLAG_VALID_CK("--chr", tmpstr);
            LOGPRINTF("--chr %s\n",tmpstr);
        }
        else if (0 == strcmp(option_str[i], "--probe")){
            prbname = option_str[++i];
            FLAG_VALID_CK("--probe", prbname);
            LOGPRINTF("--probe %s\n", prbname);
        }
        else if (0 == strcmp(option_str[i], "--from-probe")){
            fromprbname = option_str[++i];
            FLAG_VALID_CK("--from-probe", fromprbname);
            LOGPRINTF("--from-probe %s\n", fromprbname);
        }
        else if (0 == strcmp(option_str[i], "--to-probe")){
            toprbname = option_str[++i];
            FLAG_VALID_CK("--to-probe", toprbname);
            LOGPRINTF("--to-probe %s\n", toprbname);
        }
        else if(strcmp(option_str[i],"--probe-wind")==0){
            prbwindFlag=true;
            char* tmpstr=option_str[++i];
            if(tmpstr==NULL || has_prefix(tmpstr, "--")) i--;
            else prbWind=atoi(tmpstr);
            LOGPRINTF("--probe-wind %d Kb\n", prbWind);
            if(prbWind<0 )
            {
                LOGPRINTF ("Error: --probe-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--gene")){
            genename = option_str[++i];
            FLAG_VALID_CK("--gene", genename);
            LOGPRINTF("--gene %s\n", genename);
        }
        else if(strcmp(option_str[i],"--from-probe-kb")==0){
            fromprbkb=atoi(option_str[++i]);
            LOGPRINTF("--from-probe-kb %d Kb\n", fromprbkb);
            if(fromprbkb<0 )
            {
                LOGPRINTF ("Error: --from-probe-kb should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--to-probe-kb")==0){
            toprbkb=atoi(option_str[++i]);
            LOGPRINTF("--to-probe-kb %d Kb\n", toprbkb);
            if(toprbkb<0 )
            {
                LOGPRINTF ("Error: --to-probe-kb should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--probe-rm")==0){
            probe2rm=option_str[++i];
            FLAG_VALID_CK("--probe-rm", probe2rm);
            LOGPRINTF("--probe-rm %s\n", probe2rm);
        }
        else if(0==strcmp(option_str[i],"--pheno")){
            phenofileName=option_str[++i];
            FLAG_VALID_CK("--pheno", phenofileName);
            FileExist(phenofileName);
            LOGPRINTF("--pheno %s\n",phenofileName);
        }
        else if(0==strcmp(option_str[i],"--mpheno")){
            mpheno=option_str[++i];
            FLAG_VALID_CK("--mpheno", mpheno);
            LOGPRINTF("--mpheno %s\n",mpheno);
        }
        else if(0==strcmp(option_str[i],"--make-erm") || 0==strcmp(option_str[i],"--make-erm-bin")){
            make_erm_flag=true;
            erm_bin_flag=true;
            LOGPRINTF("--make-erm\n");
        }
        else if(0==strcmp(option_str[i],"--make-erm-gz")){
            make_erm_flag=true;
            erm_bin_flag=false;
            LOGPRINTF("--make-erm-gz\n");
        }
        else if(0==strcmp(option_str[i],"--erm-alg")){
            erm_alg=atoi(option_str[++i])-1;
            if(erm_alg<0 || erm_alg>2)
            {
                LOGPRINTF("Error: --erm-alg should be 1, 2 or 3.\n");
                exit (EXIT_FAILURE);
            }
            LOGPRINTF("--erm-alg %d\n",erm_alg);
        }
        else if(0==strcmp(option_str[i],"--erm-cutoff")){
            erm_cutoff=atof(option_str[++i]);
            if(erm_cutoff>=-1 || erm_cutoff<=2)
            {
                LOGPRINTF("--erm_cutoff %f\n",erm_cutoff);
            } else erm_cutoff=-2;
            
        }
        else if(0==strcmp(option_str[i],"--mlma")){
            reml_flag = false;
            mlma_flag=true;
            LOGPRINTF("--mlma\n");
        }
        else if(0==strcmp(option_str[i],"--mlma-loco")){
            mlma_loco_flag=true;
            reml_flag = false;
            LOGPRINTF("--mlma_loco\n");
        }
        else if(0==strcmp(option_str[i],"--no-fid")){
            no_fid_flag=true;
            LOGPRINTF("--no-fid\n");
        }
        else if (0 == strcmp(option_str[i], "--covar")){
            covfileName = option_str[++i];
            FLAG_VALID_CK("--covar", covfileName);
            LOGPRINTF("--covar %s\n", covfileName);
        }
        else if (0 == strcmp(option_str[i], "--qcovar")){
            qcovfileName = option_str[++i];
            FLAG_VALID_CK("--qcovar", qcovfileName);
            LOGPRINTF("--qcovar %s\n", qcovfileName);
        }
        else if(0==strcmp(option_str[i],"--erm") || 0==strcmp(option_str[i],"--erm-bin")){
            erm_file= option_str[++i];
            erm_bin_flag=true;
            FLAG_VALID_CK("--erm", erm_file);
            LOGPRINTF("--erm %s\n",erm_file);
        }
        else if(0==strcmp(option_str[i],"--subtract-erm")){
            subtract_erm_file= option_str[++i];
            erm_bin_flag=true;
            FLAG_VALID_CK("--subtract-erm", subtract_erm_file);
            LOGPRINTF("--subtract-erm %s\n",subtract_erm_file);
        }
        else if(0==strcmp(option_str[i],"--merge-erm")){
            erm_file= option_str[++i];
            m_erm_flag=true;
            FLAG_VALID_CK("--merge-erm", erm_file);
            LOGPRINTF("--merge-erm %s\n",erm_file);
        }
        else if (0==strcmp(option_str[i], "--reml-wfam") ) {
            within_family = true;
            LOGPRINTF("--reml-wfam \n");
        }
        else if (0==strcmp(option_str[i],"--reml-priors") ) {
            priors=option_str[++i];
            FLAG_VALID_CK("--reml-priors", priors);
            LOGPRINTF("--reml-priors %s\n",priors);
            
        } else if (0==strcmp(option_str[i], "--reml-priors-var")  || 0==strcmp(option_str[i], "--reml-fixed-var")) {
            string s_buf = option_str[i];
            if(s_buf == "--reml-fixed-var") reml_fixed_var_flag = true;
            priors_var=option_str[++i];
            FLAG_VALID_CK(s_buf, priors_var);
            LOGPRINTF("%s %s\n",s_buf.c_str(),priors_var);
        } else if (0==strcmp(option_str[i], "--reml-alg") ) {
            reml_mtd = atoi(option_str[++i]);
            if (reml_mtd < 0 || reml_mtd > 2)
            {
                LOGPRINTF("Error: --reml-alg should be 0, 1 or 2.\n");
                exit(EXIT_FAILURE);
            }
            LOGPRINTF("--reml-alg %d\n",reml_mtd);
          
        } else if (0==strcmp(option_str[i], "--reml-no-constrain")) {
            no_constrain = true;
            cout << "--reml-no-constrain" << endl;
        } else if (0==strcmp(option_str[i], "--reml-maxit") ) {
            MaxIter = atoi(option_str[++i]);
            if (MaxIter < 1 || MaxIter > 10000)
            {
                LOGPRINTF("Error: --reml-maxit should be within the range from 1 to 10000.\n");
                exit(EXIT_FAILURE);
            }
            LOGPRINTF("--reml-maxit %d\n",MaxIter);
        } else if (0 == strcmp(option_str[i], "--reml-bendV") ) {
            reml_force_inv_fac_flag = true;
            LOGPRINTF("--reml-bendV \n");
        } else if (0 == strcmp(option_str[i], "--reml-force-converge")) {
            reml_force_converge_flag = true;
             LOGPRINTF("--reml-force-converge \n");
        } else if (0==strcmp(option_str[i], "--reml-allow-no-converge") ) {
            reml_no_converge_flag = true;
            LOGPRINTF("--reml-allow-no-converge \n");
        }  else if (0==strcmp(option_str[i], "--reml-bivar-nocove") ) {
            ignore_Ce = true;
            LOGPRINTF("--reml-bivar-nocove \n");
        } else if (0==strcmp(option_str[i], "--reml-bivar-no-constrain")) {
            bivar_no_constrain = true;
            LOGPRINTF("--reml-bivar-no-constrain \n");            
        } else if (0==strcmp(option_str[i], "--mlma-no-adj-covar")) {
            mlma_no_adj_covar = true;
            LOGPRINTF("--mlma-no-adj-covar \n");
        } else if (0==strcmp(option_str[i], "--m2beta")) {
            m2beta = true;
            beta2m = false;
            LOGPRINTF("--m2beta \n");
        } else if (0==strcmp(option_str[i], "--beta2m")) {
            if(m2beta) {
                LOGPRINTF("Error: --m2beta should not be with --beta2m.\n");
                exit(EXIT_FAILURE);
            }
            beta2m = true;
            m2beta = false;
            LOGPRINTF("--beta2m \n");
        } else if(strcmp(option_str[i],"--diff")==0){
            diffflag=true;
            cout<<"--diff "<<endl;
        }  else if(0==strcmp(option_str[i],"--refactor")){
            refacotr_flag=true;
            out_pc_num=atoi(option_str[++i]);
            if (out_pc_num <= 0 )
            {
                LOGPRINTF("Error: --refactor should over 0.\n");
                exit(EXIT_FAILURE);
            }
            LOGPRINTF("--refactor %d\n",out_pc_num);
        } else if(0==strcmp(option_str[i],"--celltype-num")){
            celltype_num=atoi(option_str[++i]);
            if (celltype_num <= 0 )
            {
                LOGPRINTF("Error: --celltype-num should over 0.\n");
                exit(EXIT_FAILURE);
            }
            LOGPRINTF("--celltype-num %d\n",celltype_num);
        } else if(0==strcmp(option_str[i],"--dmr-num")){
            dmr_num=atoi(option_str[++i]);
            if (dmr_num <= 0 )
            {
                LOGPRINTF("Error: --dmr-num should over 0.\n");
                exit(EXIT_FAILURE);
            }
            LOGPRINTF("--dmr-num %d\n", dmr_num);
        } else if(0==strcmp(option_str[i],"--autosome-num")){
            autosome_num=atoi(option_str[++i]);
            if (autosome_num <= 0 || autosome_num > 100)
            {
                LOGPRINTF("Error: invalid number specified after the option --autosome-num.\n");
                exit(EXIT_FAILURE);
            }
            LOGPRINTF("--autosome-num %d\n", autosome_num);
        } else if(strcmp(option_str[i],"--update-epi")==0){
            update_epi_file_flag=true;
            epifname=option_str[++i];
            FLAG_VALID_CK("--update-epi", epifname);
            FileExist(epifname);
            LOGPRINTF("--update-epi %s\n",epifname);
        } else if(0==strcmp(option_str[i],"--pca")){
            pca_flag=true;
            char* tmpstr=NULL;
            tmpstr=option_str[++i];
            if(tmpstr !=NULL && has_prefix(tmpstr, "--")) i--;
            else if (tmpstr !=NULL){
                out_pc_num = atoi(tmpstr);
                if (out_pc_num <= 0 )
                {
                    LOGPRINTF("Error: the value to be specified after --pca should be positive.\n");
                    exit(EXIT_FAILURE);
                }

            }
            LOGPRINTF("--pca %d\n",out_pc_num);
        } else if(0==strcmp(option_str[i],"--assoc")){
            assoc_flag=true;
            LOGPRINTF("--assoc \n");
        }
        else if(0==strcmp(option_str[i],"--linear")){
            linear_flag=true;
            LOGPRINTF("--linear \n");
        } else if(strcmp(option_str[i],"--std")==0){
            std_thresh=atof(option_str[++i]);
            LOGPRINTF("--std %lf\n",std_thresh);
            if(std_thresh<0 || std_thresh>0.5)
            {
                LOGPRINTF("Error: --std should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        } else if (strcmp(option_str[i], "--simu-qt") == 0) {
            simu_qt_flag = true;
            cout << "--simu-qt" << endl;
        } else if (strcmp(option_str[i], "--simu-cc") == 0) {
            simu_cc = true;
            simu_case_num = atoi(option_str[++i]);
            simu_control_num = atoi(option_str[++i]);
            cout << "--simu-cc " << simu_case_num << " " << simu_control_num << endl;
            if (simu_case_num < 10) throw ("Error: --simu-cc, Invalid number of cases. Minimun number 10.");
            if (simu_control_num < 10) throw ("Error: --simu-cc, Invalid number of controls. Minimum number 10.");
        } else if (strcmp(option_str[i], "--simu-causal-loci") == 0) {
            simu_causal = option_str[++i];
            FLAG_VALID_CK("--simu-causal-loci",simu_causal);
            FileExist(simu_causal);
            LOGPRINTF("--simu-causal-loci %s\n",simu_causal);
        } else if (strcmp(option_str[i], "--simu-hsq") == 0) { //heritability
            simu_h2 = atof(option_str[++i]);
            cout << "--simu-hsq " << simu_h2 << endl;
            if (simu_h2 > 1.0 || simu_h2 < 0.0) throw ("Error: --simu-h2 should be within the range from 0 to 1.");
        } else if (strcmp(option_str[i], "--simu-k") == 0) {  // disease prevalence
            simu_K = atof(option_str[++i]);
            cout << "--simu-k " << simu_K << endl;
            if (simu_K > 0.5 || simu_K < 0.0001) throw ("Error: --simu-K should be within the range from 0.0001 to 0.5.");
        } else if (strcmp(option_str[i], "--simu-seed") == 0) {
            simu_seed = atof(option_str[++i]);
            cout << "--simu-seed " << simu_seed << endl;
            if (simu_seed <= 100) throw ("Error: --simu-seed should be >100.");
        } else if (strcmp(option_str[i], "--simu-eff-mod") == 0) {
            simu_eff_mod = atoi(option_str[++i]);
            cout << "--simu-eff-mod " << simu_eff_mod << endl;
            if (simu_eff_mod != 0 && simu_eff_mod !=1) throw ("Error: --simu-eff-mod should be 0 or 1.");
        } else if (strcmp(option_str[i], "--reml") == 0) {
            reml_flag = true;
            cout << "--reml" << endl;
            if (m_erm_flag) no_lrt = true;
        } else if (strcmp(option_str[i], "--reml-pred-rand") == 0) {
            pred_rand_eff = true;
            cout << "--reml-pred-rand" << endl;
        } else if (strcmp(option_str[i], "--reml-est-fix") == 0) {
            est_fix_eff = true;
            cout << "--reml-est-fix" << endl;
        } else if (strcmp(option_str[i], "--reml-no-lrt") == 0) {
            no_lrt = true;
            cout << "--reml-no-lrt" << endl;
        } else if (strcmp(option_str[i], "--prevalence") == 0) {
            prevalence = atof(option_str[++i]);
            cout << "--prevalence " << prevalence << endl;
            if (prevalence <= 0 || prevalence >= 1) throw ("\nError: --prevalence should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--probes-independent") == 0) {
            inde_num = atoi(option_str[++i]);
            ext_inde_flag=true;
            cout << "--probes-independent " << inde_num << endl;
            if (inde_num <= 0 ) throw ("\nError: --probes-independent should be over 0.\n");
        }
        else if (strcmp(option_str[i], "--ld-rsq") == 0) {
            ldrsq = atof(option_str[++i]);
          
            cout << "--ld-rsq " << ldrsq << endl;
            if (ldrsq <= 0 || ldrsq >= 1) throw ("\nError: --ld-rsq should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--lxpo") == 0) {
            percentage_out = atof(option_str[++i]);
            cout << "--lxpo " << percentage_out << endl;
            if (percentage_out < 0 || percentage_out > 1) {
                LOGPRINTF("\nError: --lxpo should be between 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (strcmp(option_str[i], "--upper-beta") == 0) {
            upperBeta = atof(option_str[++i]);
            cout << "--upper-beta " << upperBeta << endl;
            if (upperBeta < 0 || upperBeta > 1) throw ("\nError: --upper-beta should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--lower-beta") == 0) {
            lowerBeta = atof(option_str[++i]);
            cout << "--lower-beta " << lowerBeta << endl;
            if (lowerBeta < 0 || lowerBeta > 1) throw ("\nError: --lower-beta should be between 0 to 1.\n");
        }
        else if(0==strcmp(option_str[i],"--detection-pval-file")){
            dpvalfName=option_str[++i];
            FLAG_VALID_CK("--detection-pval-file", dpvalfName);
            FileExist(dpvalfName);
            LOGPRINTF("--detection-pval-file %s\n",dpvalfName);
        }
        else if (strcmp(option_str[i], "--dpval-thresh") == 0) {
            thresh_det_pval = atof(option_str[++i]);
            cout << "--dpval-thresh " << thresh_det_pval << endl;
            if (thresh_det_pval < 0 || thresh_det_pval > 1) throw ("\nError: --dpval-thresh should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--ratio-probe") == 0) {
            thresh_prpt_prb = atof(option_str[++i]);
            cout << "--ratio-probe " << thresh_prpt_prb << endl;
            if (thresh_prpt_prb < 0 || thresh_prpt_prb > 1) throw ("\nError: --ratio-probe should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--ratio-sample") == 0) {
            thresh_prpt_spl = atof(option_str[++i]);
            cout << "--ratio-sample " << thresh_prpt_spl << endl;
            if (thresh_prpt_spl < 0 || thresh_prpt_spl > 1) throw ("\nError: --ratio-sample should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--dpval-mth") == 0) {
            filter_det_pval_mth = atoi(option_str[++i]);
            cout << "--dpval-mth " << filter_det_pval_mth << endl;
            if (filter_det_pval_mth !=0 && filter_det_pval_mth != 1) throw ("\nError: --dpval-mth should be 0 or 1.\n");
        }
        else if (strcmp(option_str[i], "--eff-n") == 0) {
            estn_flag = true;
            cout << "--eff-n " << endl;
        }
        else if (strcmp(option_str[i], "--get-variance") == 0) {
            getvariance_flag = true;
            cout << "--get-variance " << endl;
        }
        else if (strcmp(option_str[i], "--get-mean") == 0) {
            getmean_flag = true;
            cout << "--get-mean " << endl;
            
        }
        else if (strcmp(option_str[i], "--missing-ratio-probe") == 0) {
            missing_ratio_prob = atof(option_str[++i]);
            cout << "--missing-ratio-probe " << missing_ratio_prob << endl;
            if (missing_ratio_prob < 0 || missing_ratio_prob > 1) throw ("\nError: --missing-ratio-probe should be between 0 to 1.\n");
        }
        else if (strcmp(option_str[i], "--blup-probe") == 0) {
            blup_probe_flag = true;
            blup_indi_file = option_str[++i];
            FileExist(blup_indi_file);
            cout << "--blup-probe " << blup_indi_file << endl;
        }
        else if (strcmp(option_str[i], "--score") == 0) {
            score_flag = true;
            score_file = option_str[++i];
            FileExist(score_file);
            
            char* strtmp =option_str[++i];
            if(strtmp==NULL || has_prefix(strtmp, "--"))
            {
                col_prb=1;
                col_score=2;
                i--;
                cout << "--score " << score_file << endl;
            } else {
                col_prb=atoi(strtmp);
                strtmp =option_str[++i];
                if(strtmp==NULL || has_prefix(strtmp, "--"))
                {
                    col_score=col_prb+1;
                    i--;
                    cout << "--score " << score_file <<" "<< col_prb <<endl;
                } else {
                    col_score=atoi(strtmp);
                    cout << "--score " << score_file <<" "<< col_prb<<" "<< col_score << endl;
                }
            }
        }
        else if (strcmp(option_str[i], "--impute-mean") == 0) {
            impute_mean_flag = true;
            cout << "--impute-mean " << endl;
        }
        else if (strcmp(option_str[i], "--score-has-header") == 0) {
            hasHeader = true;
            cout << "--score-has-header " << endl;
        }


    }
    
#ifndef __APPLE__
#if defined _WIN64 || defined _WIN32
    omp_set_num_threads(thread_num);
#else
    stringstream ss;
    ss << thread_num;
    setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
    omp_set_num_threads(thread_num);
    
#endif
#endif

    char tmpch[5]="orca";
    if(outfileName == NULL) outfileName=tmpch;
    
    if(merge_beed_flag) merge_beed( outfileName, befileFlstName, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, beta2m, m2beta);
    else if(make_beed_flag) make_beed( outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType, beta2m, m2beta, std_thresh,upperBeta,lowerBeta, dpvalfName,  thresh_det_pval,  thresh_prpt_prb,  thresh_prpt_spl,  filter_det_pval_mth,missing_ratio_prob);
    else if(make_efile_flag) make_efile(outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType, beta2m, m2beta,std_thresh,upperBeta,lowerBeta,transposedout,impute_mean_flag);
    else if(make_erm_flag) make_erm(outfileName, efileName,befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno,erm_bin_flag, erm_alg, beta2m, m2beta,std_thresh,upperBeta,lowerBeta, transposedin, efileType, no_fid_flag, valueType);
    else if(mlma_flag) mlma(outfileName, befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno, erm_bin_flag,  erm_alg,  covfileName, qcovfileName,  erm_file,  subtract_erm_file,  m_erm_flag,within_family,priors,priors_var,no_constrain,reml_mtd,MaxIter,reml_fixed_var_flag,reml_force_inv_fac_flag,reml_force_converge_flag, reml_no_converge_flag, mlma_no_adj_covar,percentage_out);
    else if(mlma_loco_flag) mlma_loco(outfileName, befileName,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname,prbWind,fromprbkb, toprbkb, prbwindFlag, genename, probe2rm, indilstName, indilst2remove, phenofileName, mpheno, covfileName, qcovfileName, MaxIter, priors, priors_var, no_constrain, mlma_no_adj_covar, reml_mtd, reml_fixed_var_flag, reml_force_inv_fac_flag, reml_force_converge_flag,  reml_no_converge_flag, autosome_num,percentage_out);
    else if(pca_flag) pca(outfileName, erm_file, indilstName, indilst2remove,  erm_cutoff,  m_erm_flag,  out_pc_num);
    else if(diffflag) diff(befileName,befileName2);
    else if(update_epi_file_flag) update_epifile(befileName, epifname);
    else if(refacotr_flag) getRefactor(outfileName, befileName,problstName, problst2exclde, genelistName,  chr, prbname, fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm,indilstName,indilst2remove, covfileName, qcovfileName, celltype_num,dmr_num,out_pc_num);
    else if(assoc_flag) assoc(outfileName, befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno,  covfileName, qcovfileName,std_thresh ,upperBeta,lowerBeta,estn_flag);
    else if(linear_flag) linear(outfileName, befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno,  covfileName, qcovfileName,std_thresh ,upperBeta,lowerBeta);
    else if(simu_qt_flag || simu_cc) EWAS_simu(outfileName, befileName, simu_rep, simu_causal,  simu_case_num,  simu_control_num,  simu_h2,  simu_K,  simu_seed, simu_eff_mod);
    else if(reml_flag) fit_reml(outfileName, phenofileName, mpheno, erm_bin_flag, erm_alg, covfileName, qcovfileName, erm_file, m_erm_flag, within_family, priors, priors_var, no_constrain, reml_mtd, MaxIter, reml_fixed_var_flag, reml_force_inv_fac_flag, reml_force_converge_flag, reml_no_converge_flag, mlma_no_adj_covar, pred_rand_eff, est_fix_eff, no_lrt, prevalence, mlma_flag,reml_drop,indilstName,indilst2remove,  NULL, erm_cutoff, -2.0,-2);
    else if (ext_inde_flag)  extract_inden_probes(outfileName, befileName, inde_num, ldrsq,   simu_seed, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove);
    else if(getvariance_flag || getmean_flag) getPrbVarianceMean(outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType,getvariance_flag,getmean_flag);
    else if(blup_probe_flag) blup_probe(outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType, std_thresh,upperBeta,lowerBeta,blup_indi_file);
    else if(score_flag) scoreIndividuals(outfileName,befileName, score_file, col_prb, col_score, hasHeader, phenofileName,  mpheno, problstName, problst2exclde,genelistName, chr, prbname,  fromprbname, toprbname, prbWind,fromprbkb, toprbkb,prbwindFlag, genename, probe2rm, indilstName, indilst2remove,  std_thresh, impute_mean_flag);
   /* unsigned char* wkspace_ua = NULL;
    
    unsigned char* wkspace;
    unsigned char* wkspace_base;
    uintptr_t wkspace_left;
    
    wkspace_ua = (unsigned char*)malloc((malloc_size_mb<<20) * sizeof(char));
   
    while (!wkspace_ua) {
        malloc_size_mb = (malloc_size_mb * 3) / 4;
        if (malloc_size_mb < WKSPACE_MIN_MB) {
            malloc_size_mb = WKSPACE_MIN_MB;
        }
        wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
        if (wkspace_ua) {
            LOGPRINTF("Allocated %llu MB successfully, after larger attempt(s) failed.\n", malloc_size_mb);
        } else if (malloc_size_mb == WKSPACE_MIN_MB) {
            exit(1);
        }
    }    
    
    // force 64-byte align to make cache line sensitivity work
    wkspace = (unsigned char*)CACHEALIGN((uintptr_t)wkspace_ua);
    wkspace_base = wkspace;
    wkspace_left = ((malloc_size_mb<<20) - (uintptr_t)(wkspace - wkspace_ua)) & (~(CACHELINE - ONELU));

     cout<<wkspace_left<<endl;
    */
    
    

   

}

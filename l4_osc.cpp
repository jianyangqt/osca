//
//  l4_osc.cpp
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l4_osc.h"

using namespace EFILE;
using namespace VQTL;

int main(int argc, char * argv[])
{
    /*
    FILE* fptr=fopen("test.rs10501886.S.mat","r");
    vector<string> strlist;
    vector<double> val;
    while(fgets(Tbuf, MAX_LINE_SIZE, fptr))
    {
        split_str(Tbuf,strlist,0);
        for(int i=0;i<strlist.size();i++) val.push_back(atof(strlist[i].c_str()));;
    }
    fclose(fptr);
    long dim=sqrt(val.size());
    Map<MatrixXd> raw(&val[0],dim,dim);
    cout<<raw.cols()<<":"<<raw.rows()<<endl;
     cout<<raw.determinant()<<endl;

    //inverse using eigen
    MatrixXd iraw=raw;
    bool determinant_zero = false;
    SelfAdjointEigenSolver<MatrixXd> eigensolver(iraw);
    VectorXd eval = eigensolver.eigenvalues();
    cout<<eval<<endl;
    for(int i=0;i<eval.size();i++)
    {
        if(abs(eval(i))<1e-3) {
            determinant_zero=true;
            eval(i)=0;
        } else {
            eval(i) = 1.0 / eval(i);
        }
    }
    iraw = eigensolver.eigenvectors() * DiagonalMatrix<double, Dynamic, Dynamic>(eval) * eigensolver.eigenvectors().transpose();
    MatrixXd Iden=raw*iraw;
    string filename="test.rs10501886.eigen.inv.mat";
    FILE* tmpfile=fopen(filename.c_str(),"w");
    if(!tmpfile)
    {
        LOGPRINTF("error open file.\n");
        TERMINATE();
    }
    for(int t=0;t<iraw.rows();t++)
    {
        string str="";
        for(int k=0;k<iraw.cols();k++)
        {
            str +=atos(iraw(t,k)) + '\t';
        }
        str += '\n';
        fputs(str.c_str(),tmpfile);
    }
    fclose(tmpfile);
    
    //inverse using SVD
    iraw=raw;
    JacobiSVD<MatrixXd> svd(iraw, ComputeThinU | ComputeThinV);
    VectorXd sin=svd.singularValues();
    cout<<sin<<endl;
    for(int i=0;i<sin.size();i++) {
        if(abs(sin[i])<1e-6) {
            sin[i] =0;
        } else {
            sin[i] = 1/sin[i];
        }
    }
    MatrixXd dd=svd.matrixU()*sin.asDiagonal()*svd.matrixV().transpose();
    Iden=raw*dd;
    filename="test.rs10501886.svd.inv.mat";
    tmpfile=fopen(filename.c_str(),"w");
    if(!tmpfile)
    {
        LOGPRINTF("error open file.\n");
        TERMINATE();
    }
    for(int t=0;t<dd.rows();t++)
    {
        string str="";
        for(int k=0;k<dd.cols();k++)
        {
            str +=atos(dd(t,k)) + '\t';
        }
        str += '\n';
        fputs(str.c_str(),tmpfile);
    }
    fclose(tmpfile);
     */

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
    cout << "* OmicS-data-based Complex trait Analysis (OSCA)" << endl;
    cout << "* version 0.392" << endl;
    cout << "* (C) 2016 Futao Zhang, Zhihong Zhu and Jian Yang" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;
    
    char tmpname[5]="osca";
    char* outname=tmpname;
    string logfname=string(outname)+".log";
    logfile=fopen(logfname.c_str(),"w");
    if (!logfile) {
        printf("Error: Failed to open log file %s.\n", logfname.c_str());
        exit(EXIT_FAILURE);
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
    
    // for ORM
    bool make_erm_flag=false;
    bool erm_bin_flag=true, grm_bin_flag=true;
    int erm_alg=0; // 0 for standardizing probes, 1 for centering probes, 2 for standardizing individuals, 3 for iteratively standardizing probes and individuals
    char* erm_file=NULL; // can be an erm file name or an erm filelist file name
    char* grm_file=NULL;
    char* subtract_erm_file=NULL;
    bool m_erm_flag=false;
    char* priors= NULL;
    char* priors_var= NULL;
    bool reml_fixed_var_flag=false;
    double erm_cutoff = -2.0;
    bool erm_cutoff_2sides = false;
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
    int simu_rep = 1, simu_case_num = 0, simu_control_num = 0, simu_eff_mod = 0;  // 0 standarise the probes. 1 use raw probe profiles.
    char* simu_causal = NULL;
    double simu_h2 = 0.1, simu_K = 0.1, simu_seed = -rand_seed();
    
    bool reml_flag = false, pred_rand_eff = false, est_fix_eff = false, no_lrt = false;
    double prevalence = -2.0, prevalence2 = -2.0;
    vector<int> reml_drop;
    reml_drop.push_back(1);

    int inde_num=0;
    bool ext_inde_flag=false;
    double ldrsq=1.0;
    
    double percentage_out = 0.0; // 0: leave top 0% of probes of uncorrected association out from the ORM calculation.  1: leave top 100% of probes out from the ORM calculation.
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
    //
    int tsk_ttl = 1;
    int tsk_id = 1;
    bool vqtl=false;
    int vqtl_mtd=0; // 0 for Bartlett’s test; 1 for Levene’s test with mean; 2 for Levene’s test with median; 3 for Fligner-Killeen test
    char* bFileName = NULL;
    bool adjprb=false;
    
    //
    bool eqtl=false;
    
    bool prt_residiual=false;
    bool simu_residual_only=false;
    
    bool besd_snp_major = false;
    
    //for smr
    bool queryBesd = false;
    double pQueryBesd = 5e-8;
    char* beqtlFileName = NULL;
    char* beqtlFileName2 = NULL;
    int prbchr = -9;
    int snpchr = -9;
    char* snprs = NULL;
    char* fromsnprs = NULL;
    char* tosnprs = NULL;
    int snpWind=50;   //Kb
    bool snpwindFlag = false;
    int fromsnpkb = -9;
    int tosnpkb = -9;
    char* snp2rm=NULL;
    bool cis_flag=false;
    int cis_itvl=2000;
    char* besdflstName = NULL;
    bool combineFlg =false;
    bool metaflg = false;
    bool smr_flag= false;
    char* gwasflstName = NULL;
    
    
    int meta_mtd = 0; //0 for traditional meta, 1 for MeCS
    double pmecs = 0.01;
    int mecs_mth = 0; // 1 for est_cor, 0 for pcc
    char* corMatFName = NULL;
    bool zflag = false; // estimate pcc with z
    
    bool make_besd_flag = false;
    bool save_dense_flag = false;
    bool to_smr_flag = false;
    bool besd_shrink_flag = false; //remove null probes and null SNPs
    
    
    
    //private use
    bool enveff = false;
    char* efffname = NULL;
    char* effproblstName = NULL;
    
    bool stdprb =false;
    char* freqFName = NULL;
    char* varFName =NULL;
    
    double lambda_wind = 0.01;
    bool fastlinear =false;
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
        else if(0==strcmp(option_str[i],"--make-bod")){
            make_beed_flag=true;
            LOGPRINTF("--make-bod\n");
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
        else if(strcmp(option_str[i],"--exclude-snp")==0){
            snplst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-snp", snplst2exclde);
            LOGPRINTF("--exclude-snp %s\n",snplst2exclde);
            FileExist(snplst2exclde);
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
                TERMINATE();
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
            if(strncmp(tmpstr,"X",1)==0 || strncmp(tmpstr,"x",1)==0) chr=23;
            else if(strncmp(tmpstr,"Y",1)==0 || strncmp(tmpstr,"y",1)==0) chr=24;
            else chr=atoi(tmpstr);
            FLAG_VALID_CK("--chr", tmpstr);
            if(chr<=0)
            {
                LOGPRINTF("Error: --chr should be over 0.\n");
                TERMINATE();
            }
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
                TERMINATE();
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
                TERMINATE();
            }
        }
        else if(strcmp(option_str[i],"--to-probe-kb")==0){
            toprbkb=atoi(option_str[++i]);
            LOGPRINTF("--to-probe-kb %d Kb\n", toprbkb);
            if(toprbkb<0 )
            {
                LOGPRINTF ("Error: --to-probe-kb should be over 0.\n");
                TERMINATE();
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
        else if(0==strcmp(option_str[i],"--make-orm") || 0==strcmp(option_str[i],"--make-orm-bin")){
            make_erm_flag=true;
            erm_bin_flag=true;
            LOGPRINTF("--make-orm\n");
        }
        else if(0==strcmp(option_str[i],"--make-orm-gz")){
            make_erm_flag=true;
            erm_bin_flag=false;
            LOGPRINTF("--make-orm-gz\n");
        }
        else if(0==strcmp(option_str[i],"--orm-alg")){
            erm_alg=atoi(option_str[++i])-1;
            if(erm_alg<0 || erm_alg>2)
            {
                LOGPRINTF("Error: --orm-alg should be 1, 2 or 3.\n");
                TERMINATE();
            }
            LOGPRINTF("--orm-alg %d\n",erm_alg+1);
        }
        else if(0==strcmp(option_str[i],"--orm-cutoff")){
            erm_cutoff=atof(option_str[++i]);
            if(erm_cutoff>=-1 || erm_cutoff<=2)
            {
                LOGPRINTF("--orm-cutoff %f\n",erm_cutoff);
            } else erm_cutoff=-2;
        }
        else if(0==strcmp(option_str[i],"--orm-cutoff-2sides")){
            erm_cutoff_2sides=true;
            LOGPRINTF("--orm-cutoff-2sides\n");
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
        else if(0==strcmp(option_str[i],"--orm") || 0==strcmp(option_str[i],"--orm-bin")){
            erm_file= option_str[++i];
            erm_bin_flag=true;
            FLAG_VALID_CK("--orm", erm_file);
            LOGPRINTF("--orm %s\n",erm_file);
        }
        else if(0==strcmp(option_str[i],"--grm") || 0==strcmp(option_str[i],"--grm-bin")){
            grm_file= option_str[++i];
            grm_bin_flag=true;
            FLAG_VALID_CK("--grm", grm_file);
            LOGPRINTF("--grm %s\n",grm_file);
        }
        else if(0==strcmp(option_str[i],"--subtract-orm")){
            subtract_erm_file= option_str[++i];
            erm_bin_flag=true;
            FLAG_VALID_CK("--subtract-orm", subtract_erm_file);
            LOGPRINTF("--subtract-orm %s\n",subtract_erm_file);
        }
        else if(0==strcmp(option_str[i],"--merge-orm")){
            erm_file= option_str[++i];
            m_erm_flag=true;
            FLAG_VALID_CK("--merge-orm", erm_file);
            LOGPRINTF("--merge-orm %s\n",erm_file);
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
                TERMINATE();
            }
            LOGPRINTF("--reml-alg %d\n",reml_mtd);
          
        } else if (0==strcmp(option_str[i], "--reml-no-constrain")) {
            no_constrain = true;
            LOGPRINTF("--reml-no-constrain\n");
        } else if (0==strcmp(option_str[i], "--reml-maxit") ) {
            MaxIter = atoi(option_str[++i]);
            if (MaxIter < 1 || MaxIter > 100000)
            {
                LOGPRINTF("Error: --reml-maxit should be within the range from 1 to 100000.\n");
                TERMINATE();
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
                TERMINATE();
            }
            beta2m = true;
            m2beta = false;
            LOGPRINTF("--beta2m \n");
        } else if(strcmp(option_str[i],"--diff")==0){
            diffflag=true;
            LOGPRINTF("--diff \n");
        }  else if(0==strcmp(option_str[i],"--refactor")){
            refacotr_flag=true;
            out_pc_num=atoi(option_str[++i]);
            if (out_pc_num <= 0 )
            {
                LOGPRINTF("Error: --refactor should over 0.\n");
                TERMINATE();
            }
            LOGPRINTF("--refactor %d\n",out_pc_num);
        } else if(0==strcmp(option_str[i],"--celltype-num")){
            celltype_num=atoi(option_str[++i]);
            if (celltype_num <= 0 )
            {
                LOGPRINTF("Error: --celltype-num should over 0.\n");
                TERMINATE();
            }
            LOGPRINTF("--celltype-num %d\n",celltype_num);
        } else if(0==strcmp(option_str[i],"--dmr-num")){
            dmr_num=atoi(option_str[++i]);
            if (dmr_num <= 0 )
            {
                LOGPRINTF("Error: --dmr-num should over 0.\n");
                TERMINATE();
            }
            LOGPRINTF("--dmr-num %d\n", dmr_num);
        } else if(0==strcmp(option_str[i],"--autosome-num")){
            autosome_num=atoi(option_str[++i]);
            if (autosome_num <= 0 || autosome_num > 100)
            {
                LOGPRINTF("Error: invalid number specified after the option --autosome-num.\n");
                TERMINATE();
            }
            LOGPRINTF("--autosome-num %d\n", autosome_num);
        } else if(strcmp(option_str[i],"--update-opi")==0){
            update_epi_file_flag=true;
            epifname=option_str[++i];
            FLAG_VALID_CK("--update-opi", epifname);
            FileExist(epifname);
            LOGPRINTF("--update-opi %s\n",epifname);
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
                    TERMINATE();
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
        } else if(strcmp(option_str[i],"--sd-min")==0){
            std_thresh=atof(option_str[++i]);
            LOGPRINTF("--sd-min %lf\n",std_thresh);
            if(std_thresh<0 || std_thresh>1)
            {
                LOGPRINTF("Error: --sd-min should be within the range from 0 to 1.\n");
                TERMINATE();
            }
        } else if (strcmp(option_str[i], "--simu-qt") == 0) {
            simu_qt_flag = true;
            LOGPRINTF("--simu-qt\n");
        } else if (strcmp(option_str[i], "--simu-cc") == 0) {
            simu_cc = true;
            simu_case_num = atoi(option_str[++i]);
            simu_control_num = atoi(option_str[++i]);
            LOGPRINTF( "--simu-cc %d\t%d\n",simu_case_num , simu_control_num);
            if (simu_case_num < 10) {
                LOGPRINTF("Error: --simu-cc, Invalid number of cases. Minimun number 10.\n");
                TERMINATE();
            }
            if (simu_control_num < 10) {
                LOGPRINTF("Error: --simu-cc, Invalid number of controls. Minimum number 10.\n");
                TERMINATE();
            }
        } else if (strcmp(option_str[i], "--simu-causal-loci") == 0) {
            simu_causal = option_str[++i];
            FLAG_VALID_CK("--simu-causal-loci",simu_causal);
            FileExist(simu_causal);
            LOGPRINTF("--simu-causal-loci %s\n",simu_causal);
        } else if (strcmp(option_str[i], "--simu-hsq") == 0) { //heritability
            simu_h2 = atof(option_str[++i]);
            LOGPRINTF("--simu-hsq %f\n", simu_h2);
            if (simu_h2 > 1.0 || simu_h2 < 0.0) {
                LOGPRINTF("Error: --simu-h2 should be within the range from 0 to 1.\n");
                TERMINATE();
            }
        } else if (strcmp(option_str[i], "--simu-k") == 0) {  // disease prevalence
            simu_K = atof(option_str[++i]);
            LOGPRINTF("--simu-k %f\n", simu_K);
            if (simu_K > 0.5 || simu_K < 0.0001) {
                LOGPRINTF("Error: --simu-K should be within the range from 0.0001 to 0.5.\n");
                TERMINATE();
            }
        } else if (strcmp(option_str[i], "--simu-seed") == 0) {
            simu_seed = atof(option_str[++i]);
            LOGPRINTF("--simu-seed %f\n", simu_seed);
            if (simu_seed <= 100) {
                LOGPRINTF("Error: --simu-seed should be >100.\n");
                TERMINATE();
            }
        } else if (strcmp(option_str[i], "--simu-eff-mod") == 0) {
            simu_eff_mod = atoi(option_str[++i]);
            LOGPRINTF("--simu-eff-mod %d\n", simu_eff_mod);
            if (simu_eff_mod != 0 && simu_eff_mod !=1) {
                LOGPRINTF("Error: --simu-eff-mod should be 0 or 1.\n");
                TERMINATE();
            }
        } else if (strcmp(option_str[i], "--reml") == 0) {
            reml_flag = true;
            LOGPRINTF( "--reml\n");
            if (m_erm_flag) no_lrt = true;
        } else if (strcmp(option_str[i], "--reml-pred-rand") == 0) {
            pred_rand_eff = true;
            LOGPRINTF("--reml-pred-rand\n");
        } else if (strcmp(option_str[i], "--reml-est-fix") == 0) {
            est_fix_eff = true;
            LOGPRINTF("--reml-est-fix\n");
        } else if (strcmp(option_str[i], "--reml-no-lrt") == 0) {
            no_lrt = true;
            LOGPRINTF("--reml-no-lrt\n");
        } else if (strcmp(option_str[i], "--prevalence") == 0) {
            prevalence = atof(option_str[++i]);
            LOGPRINTF("--prevalence %f\n", prevalence);
            if (prevalence <= 0 || prevalence >= 1) {
                LOGPRINTF("\nError: --prevalence should be between 0 to 1.\n")
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--probes-independent") == 0) {
            inde_num = atoi(option_str[++i]);
            ext_inde_flag=true;
            LOGPRINTF("--probes-independent %d\n", inde_num);
            if (inde_num <= 0 ) {
                LOGPRINTF("\nError: --probes-independent should be over 0.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--ld-rsq") == 0) {
            ldrsq = atof(option_str[++i]);
          
            LOGPRINTF( "--ld-rsq %f\n", ldrsq);
            if (ldrsq <= 0 || ldrsq >= 1) {
                LOGPRINTF("\nError: --ld-rsq should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--lxpo") == 0) {
            percentage_out = atof(option_str[++i]);
            if (percentage_out < 0 || percentage_out > 100) {
                LOGPRINTF("\nError: --lxpo should be between 0 to 100.\n");
                TERMINATE();
            }
            LOGPRINTF("--lxpo %f%%\n", percentage_out);
            percentage_out*=0.01;
        }
        else if (strcmp(option_str[i], "--upper-beta") == 0) {
            upperBeta = atof(option_str[++i]);
            LOGPRINTF("--upper-beta %f\n", upperBeta);
            if (upperBeta < 0 || upperBeta > 1) {
                LOGPRINTF("\nError: --upper-beta should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--lower-beta") == 0) {
            lowerBeta = atof(option_str[++i]);
            LOGPRINTF("--lower-beta %f\n", lowerBeta);
            if (lowerBeta < 0 || lowerBeta > 1) {
                LOGPRINTF("\nError: --lower-beta should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if(0==strcmp(option_str[i],"--detection-pval-file")){
            dpvalfName=option_str[++i];
            FLAG_VALID_CK("--detection-pval-file", dpvalfName);
            FileExist(dpvalfName);
            LOGPRINTF("--detection-pval-file %s\n",dpvalfName);
        }
        else if (strcmp(option_str[i], "--dpval-thresh") == 0) {
            thresh_det_pval = atof(option_str[++i]);
            LOGPRINTF("--dpval-thresh %f\n", thresh_det_pval );
            if (thresh_det_pval < 0 || thresh_det_pval > 1) {
                LOGPRINTF("\nError: --dpval-thresh should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--ratio-probe") == 0) {
            thresh_prpt_prb = atof(option_str[++i]);
            LOGPRINTF("--ratio-probe %f\n", thresh_prpt_prb);
            if (thresh_prpt_prb < 0 || thresh_prpt_prb > 1) {
                LOGPRINTF("\nError: --ratio-probe should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--ratio-sample") == 0) {
            thresh_prpt_spl = atof(option_str[++i]);
            LOGPRINTF("--ratio-sample %f\n",  thresh_prpt_spl);
            if (thresh_prpt_spl < 0 || thresh_prpt_spl > 1) {
                LOGPRINTF("\nError: --ratio-sample should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--dpval-mth") == 0) {
            filter_det_pval_mth = atoi(option_str[++i]);
            LOGPRINTF( "--dpval-mth %d\n",  filter_det_pval_mth);
            if (filter_det_pval_mth !=0 && filter_det_pval_mth != 1) {
                LOGPRINTF("\nError: --dpval-mth should be 0 or 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--eff-n") == 0) {
            estn_flag = true;
            LOGPRINTF("--eff-n \n");
        }
        else if (strcmp(option_str[i], "--get-variance") == 0) {
            getvariance_flag = true;
            LOGPRINTF("--get-variance \n");
        }
        else if (strcmp(option_str[i], "--get-mean") == 0) {
            getmean_flag = true;
            LOGPRINTF("--get-mean \n");
            
        }
        else if (strcmp(option_str[i], "--missing-ratio-probe") == 0) {
            missing_ratio_prob = atof(option_str[++i]);
            LOGPRINTF( "--missing-ratio-probe %f\n", missing_ratio_prob );
            if (missing_ratio_prob < 0 || missing_ratio_prob > 1) {
                LOGPRINTF("\nError: --missing-ratio-probe should be between 0 to 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--blup-probe") == 0) {
            blup_probe_flag = true;
            blup_indi_file = option_str[++i];
            FileExist(blup_indi_file);
            LOGPRINTF( "--blup-probe %s\n", blup_indi_file);
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
                LOGPRINTF("--score %s\n", score_file);
            } else {
                col_prb=atoi(strtmp);
                strtmp =option_str[++i];
                if(strtmp==NULL || has_prefix(strtmp, "--"))
                {
                    col_score=col_prb+1;
                    i--;
                    LOGPRINTF( "--score %s\t%d\n", score_file ,col_prb);
                } else {
                    col_score=atoi(strtmp);
                    LOGPRINTF("--score %s\t%d\t%d\n",  score_file , col_prb, col_score);
                }
            }
        }
        else if (strcmp(option_str[i], "--impute-mean") == 0) {
            impute_mean_flag = true;
            LOGPRINTF("--impute-mean \n");
        }
        else if (strcmp(option_str[i], "--score-has-header") == 0) {
            hasHeader = true;
            LOGPRINTF( "--score-has-header \n");
        }
        else if (strcmp(option_str[i], "--task-total") == 0) {
            tsk_ttl = atoi(option_str[++i]);
            LOGPRINTF( "--task-total %d\n", tsk_ttl);
            if(tsk_ttl<1 )
            {
                LOGPRINTF ("Error: --task-total should be over 1.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--task-id") == 0) {
            tsk_id = atoi(option_str[++i]);
            LOGPRINTF( "--task-id %d\n", tsk_id);
            if(tsk_id<1 )
            {
                LOGPRINTF ("Error: --task-id should be over 1.\n");
                TERMINATE();
            }
            if(tsk_id > tsk_ttl )
            {
                LOGPRINTF ("Error: --task-id should not be larger than --task-total.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--vqtl") == 0) {
            vqtl = true;
            LOGPRINTF( "--vqtl \n");
        }
        else if(strcmp(option_str[i],"--bfile")==0){
            bFileName=option_str[++i];
            FLAG_VALID_CK("--bfile", bFileName);
            LOGPRINTF("--bfile %s\n",bFileName);
        }
        else if(0==strcmp(option_str[i],"--vqtl-mtd")){
            vqtl_mtd=atoi(option_str[++i]);
            if(vqtl_mtd<0 || vqtl_mtd>4)
            {
                LOGPRINTF("Error: --vqtl-mtd should be 0, 1, 2 or 3.\n");
                TERMINATE();
            }
            LOGPRINTF("--vqtl-mtd %d\n",vqtl_mtd);
        }
        else if (strcmp(option_str[i], "--adj-probe") == 0) {
            adjprb = true;
            LOGPRINTF( "--adj-probe \n");
        }
        else if (strcmp(option_str[i], "--output-residual") == 0) {
            prt_residiual = true;
            LOGPRINTF( "--output-residual \n");
        }
        else if (strcmp(option_str[i], "--simu-residual") == 0) {
            simu_residual_only = true;
            LOGPRINTF( "--simu-residua \n");
        }
        else if (strcmp(option_str[i], "--eqtl") == 0) {
            eqtl = true;
            LOGPRINTF( "--eqtl \n");
        }
        else if (strcmp(option_str[i], "--besd-snp-major") == 0) {
            besd_snp_major = true;
            LOGPRINTF( "--besd-snp-major \n");
        }
        else if(strcmp(option_str[i], "--query")==0) {
            queryBesd=true;
            if(i+1==option_num || has_prefix(option_str[i+1],"--"))  pQueryBesd = 5.0e-8;
            else pQueryBesd = atof (option_str[++i]);
            LOGPRINTF("--query %10.2e\n", pQueryBesd);
            if(pQueryBesd<0 || pQueryBesd>1)
            {
                LOGPRINTF("Error: --query should be within the range from 0 to 1.\n");
                TERMINATE();
            }
        }
        else if(strcmp(option_str[i],"--beqtl-summary")==0){
            if(beqtlFileName==NULL)
            {
                beqtlFileName=option_str[++i];
                FLAG_VALID_CK("--beqtl-summary", beqtlFileName);
                LOGPRINTF("--beqtl-summary %s\n",beqtlFileName);
            }else{
                beqtlFileName2=option_str[++i];
                FLAG_VALID_CK("--beqtl-summary", beqtlFileName2);
                LOGPRINTF("--beqtl-summary %s\n",beqtlFileName2);
            }
        }
        else if(strcmp(option_str[i],"--probe-chr")==0){
            char* tmpstr=option_str[++i];
            if(strncmp(tmpstr,"X",1)==0 || strncmp(tmpstr,"x",1)==0) prbchr=23;
            else if(strncmp(tmpstr,"Y",1)==0 || strncmp(tmpstr,"y",1)==0) prbchr=24;
            else prbchr=atoi(tmpstr);
            FLAG_VALID_CK("--probe-chr", tmpstr);
            if(chr<=0)
            {
                LOGPRINTF("Error: --probe-chr should be over 0.\n");
                TERMINATE();
            }
            LOGPRINTF("--probe-chr %s\n",tmpstr);
        }
        else if(strcmp(option_str[i],"--snp-chr")==0){
            char* tmpstr=option_str[++i];
            if(strncmp(tmpstr,"X",1)==0 || strncmp(tmpstr,"x",1)==0) snpchr=23;
            else if(strncmp(tmpstr,"Y",1)==0 || strncmp(tmpstr,"y",1)==0) snpchr=24;
            else snpchr=atoi(tmpstr);
            FLAG_VALID_CK("--snp-chr", tmpstr);
            if(chr<=0)
            {
                LOGPRINTF("Error: --snp-chr should be over 0.\n");
                TERMINATE();
            }
            LOGPRINTF("--snp-chr %s\n",tmpstr);
        }
        else if (strcmp(option_str[i], "--snp")==0){
            snprs = option_str[++i];
            FLAG_VALID_CK("--snp", snprs);
            LOGPRINTF("--snp %s\n", snprs);
        }
        else if (0 == strcmp(option_str[i], "--from-snp")){
            fromsnprs = option_str[++i];
            FLAG_VALID_CK("--from-snp", fromsnprs);
            printf("--from-snp %s\n", fromsnprs);
        }
        else if (0 == strcmp(option_str[i], "--to-snp")){
            tosnprs = option_str[++i];
            FLAG_VALID_CK("--to-snp", tosnprs);
            printf("--to-snp %s\n", tosnprs);
        }
        else if(strcmp(option_str[i],"--snp-wind")==0){
            
            snpwindFlag=true;
            char* tmpstr=option_str[++i];
            if(tmpstr==NULL || has_prefix(tmpstr, "--")) i--;
            else snpWind=atoi(tmpstr);
            LOGPRINTF("--snp-wind %d Kb\n", snpWind);
            if(snpWind<0 )
            {
                LOGPRINTF ("Error: --snp-wind should be over 0.\n");
                TERMINATE();
            }
        }
        else if(strcmp(option_str[i],"--from-snp-kb")==0){
            fromsnpkb=atoi(option_str[++i]);
            printf("--from-snp-kb %d Kb\n", fromsnpkb);
            if(fromsnpkb<0 )
            {
                LOGPRINTF("Error: --from-snp-kb should be over 0.\n");
                TERMINATE();
            }
        }
        else if(strcmp(option_str[i],"--to-snp-kb")==0){
            tosnpkb=atoi(option_str[++i]);
            printf("--to-snp-kb %d Kb\n", tosnpkb);
            if(tosnpkb<0 )
            {
                LOGPRINTF("Error: --to-snp-kb should be over 0.\n");
                TERMINATE();
            }
        }
        else if(strcmp(option_str[i],"--snp-rm")==0){
            snp2rm=option_str[++i];
            FLAG_VALID_CK("--snp-rm", snp2rm);
            LOGPRINTF("--snp-rm %s\n", snp2rm);
        }
        else if(strcmp(option_str[i],"--cis-wind")==0){
            cis_flag=true;
            cis_itvl=atoi(option_str[++i]);
            LOGPRINTF("--cis-wind %d Kb\n", cis_itvl);
            if(cis_itvl<=0 )
            {
                LOGPRINTF("Error: --cis-wind should be over 0.\n");
                TERMINATE();
            }
        }
        else if (strcmp(option_str[i], "--besd-flist")==0){
            combineFlg=true;
            besdflstName = option_str[++i];
            FLAG_VALID_CK("--besd-flist", besdflstName);
            LOGPRINTF("--besd-flist %s\n", besdflstName);
            FileExist(besdflstName);
        }
        else if (strcmp(option_str[i], "--gwas-flist")==0){
            gwasflstName = option_str[++i];
            FLAG_VALID_CK("--gwas-flist", gwasflstName);
            LOGPRINTF("--gwas-flist %s\n", gwasflstName);
            FileExist(gwasflstName);
        }
        else if(strcmp(option_str[i],"--meta")==0){
            metaflg=true;
            combineFlg=false;
            smr_flag=false;
            LOGPRINTF("--meta\n");
        }
        else if(0==strcmp(option_str[i],"--mecs")){
            metaflg=true;
            meta_mtd=1;
            LOGPRINTF("--mecs\n");
        }
        else if(0==strcmp(option_str[i],"--pmecs")){
            pmecs=atof(option_str[++i]);
            if (pmecs < 0 || pmecs > 1) {
                LOGPRINTF("\nError: --pmecs should be between 0 to 1.\n");
                TERMINATE();
            }
            LOGPRINTF("--pmecs %0.2e\n", pmecs);
        }
        else if(0==strcmp(option_str[i],"--mecs-mth")){
            mecs_mth=atoi(option_str[++i]);
            LOGPRINTF( "--mecs-mth %d\n",  mecs_mth);
            if (mecs_mth !=0 && mecs_mth != 1) {
                LOGPRINTF("\nError: --mecs-mth should be 0 or 1.\n");
                TERMINATE();
            }
        }
        else if(0==strcmp(option_str[i],"--make-besd")){
            make_besd_flag=true;
            LOGPRINTF("--make-besd\n");
        }
        else if(0==strcmp(option_str[i],"--to-smr")){
            make_besd_flag=true;
            to_smr_flag=true;
            LOGPRINTF("--to-smr\n");
        }
        else if (0 == strcmp(option_str[i], "--make-besd-dense")){
            make_besd_flag=true;
            save_dense_flag=true;
            LOGPRINTF("--make-besd-dense \n");
        }
        else if (0 == strcmp(option_str[i], "--besd-shrink")){
            make_besd_flag=true;
            besd_shrink_flag=true;
            LOGPRINTF("--besd-shrink \n");
        }
        else if (0 == strcmp(option_str[i], "--add-eff")){
            enveff=true;
            LOGPRINTF("--add-eff \n");
            LOGPRINTF("NOTE: Private use to add an environmental effect ( N(0,sd) ) to each probe. \n");
        }
        else if (strcmp(option_str[i], "--cor-mat")==0){
            corMatFName = option_str[++i];
            FLAG_VALID_CK("--cor-mat", corMatFName);
            LOGPRINTF("--cor-mat %s\n", corMatFName);
            FileExist(corMatFName);
        }
        else if (strcmp(option_str[i], "--eff-file")==0){
            efffname = option_str[++i];
            FLAG_VALID_CK("--eff-file", efffname);
            LOGPRINTF("--eff-file %s\n", efffname);
            FileExist(efffname);
        }
        else if(strcmp(option_str[i],"--eff-probe")==0){
            effproblstName=option_str[++i];
            FLAG_VALID_CK("--eff-probe", effproblstName);
            FileExist(effproblstName);
            LOGPRINTF("--eff-probe %s\n",effproblstName);
        }
        else if(0==strcmp(option_str[i],"--pcc-z")){
            zflag=true;
            LOGPRINTF("--pcc-z\n");
        }
        else if (0 == strcmp(option_str[i], "--prt-mid")){
            prt_mid_rlt=true;
            LOGPRINTF("--prt-mid \n");
        }
        else if (strcmp(option_str[i], "--std-probe") == 0) {
            stdprb = true;
            LOGPRINTF( "--std-probe \n");
        }
        else if(strcmp(option_str[i],"--freq-file")==0){
            freqFName=option_str[++i];
            FLAG_VALID_CK("--freq-file", freqFName);
            FileExist(freqFName);
            LOGPRINTF("--freq-file %s\n",freqFName);
        }
        else if(strcmp(option_str[i],"--var-file")==0){
            varFName=option_str[++i];
            FLAG_VALID_CK("--var-file", varFName);
            FileExist(varFName);
            LOGPRINTF("--var-file %s\n",varFName);
        }
        else if(0==strcmp(option_str[i],"--lambda-wind")){
            lambda_wind=atof(option_str[++i]);
            if (lambda_wind < 0 ) {
                LOGPRINTF("\nError: --lambda-wind should be over 0 .\n");
                TERMINATE();
            }
            LOGPRINTF("--lambda-wind %0.2e\n", lambda_wind);
        }
        else if(0==strcmp(option_str[i],"--fast-linear")){
            fastlinear=true;
            LOGPRINTF("--fast-linear \n");
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

    char tmpch[5]="osca";
    if(outfileName == NULL) outfileName=tmpch;
    
    if(merge_beed_flag) merge_beed( outfileName, befileFlstName, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, beta2m, m2beta);
    else if(make_beed_flag) make_beed( outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType, beta2m, m2beta, std_thresh,upperBeta,lowerBeta, dpvalfName,  thresh_det_pval,  thresh_prpt_prb,  thresh_prpt_spl,  filter_det_pval_mth,missing_ratio_prob, adjprb,  covfileName, qcovfileName,enveff, effproblstName, efffname);
    else if(make_efile_flag) make_efile(outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType, beta2m, m2beta,std_thresh,upperBeta,lowerBeta,transposedout,impute_mean_flag);
    else if(make_erm_flag) make_erm(outfileName, efileName,befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno,erm_bin_flag, erm_alg, beta2m, m2beta,std_thresh,upperBeta,lowerBeta, transposedin, efileType, no_fid_flag, valueType);
    else if(mlma_flag) mlma(outfileName, befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno, erm_bin_flag,  erm_alg,  covfileName, qcovfileName,  erm_file,  subtract_erm_file,  m_erm_flag,within_family,priors,priors_var,no_constrain,reml_mtd,MaxIter,reml_fixed_var_flag,reml_force_inv_fac_flag,reml_force_converge_flag, reml_no_converge_flag, mlma_no_adj_covar,percentage_out, lambda_wind,fastlinear);
    else if(mlma_loco_flag) mlma_loco(outfileName, befileName,problstName,problst2exclde,genelistName, chr,prbname, fromprbname, toprbname,prbWind,fromprbkb, toprbkb, prbwindFlag, genename, probe2rm, indilstName, indilst2remove, phenofileName, mpheno, covfileName, qcovfileName, MaxIter, priors, priors_var, no_constrain, mlma_no_adj_covar, reml_mtd, reml_fixed_var_flag, reml_force_inv_fac_flag, reml_force_converge_flag,  reml_no_converge_flag, autosome_num,percentage_out,lambda_wind,fastlinear);
    else if(pca_flag) pca(outfileName, erm_file, indilstName, indilst2remove,  erm_cutoff,erm_cutoff_2sides, m_erm_flag,  out_pc_num);
    else if(diffflag) diff(befileName,befileName2);
    else if(update_epi_file_flag) update_epifile(befileName, epifname);
    else if(refacotr_flag) getRefactor(outfileName, befileName,problstName, problst2exclde, genelistName,  chr, prbname, fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm,indilstName,indilst2remove, covfileName, qcovfileName, celltype_num,dmr_num,out_pc_num);
    else if(assoc_flag) assoc(outfileName, befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno,  covfileName, qcovfileName,std_thresh ,upperBeta,lowerBeta,estn_flag);
    else if(linear_flag) linear(outfileName, befileName, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove,phenofileName,mpheno,  covfileName, qcovfileName,std_thresh ,upperBeta,lowerBeta,tsk_ttl,tsk_id,fastlinear);
    else if(simu_qt_flag || simu_cc) EWAS_simu(outfileName, befileName, simu_rep, simu_causal,  simu_case_num,  simu_control_num,  simu_h2,  simu_K,  simu_seed, simu_eff_mod,simu_residual_only);
    else if(reml_flag) fit_reml(outfileName, phenofileName, mpheno, erm_bin_flag, grm_bin_flag, erm_alg, covfileName, qcovfileName, erm_file, grm_file, m_erm_flag, within_family, priors, priors_var, no_constrain, reml_mtd, MaxIter, reml_fixed_var_flag, reml_force_inv_fac_flag, reml_force_converge_flag, reml_no_converge_flag, mlma_no_adj_covar, pred_rand_eff, est_fix_eff, no_lrt, prevalence, mlma_flag,reml_drop,indilstName,indilst2remove,  NULL, erm_cutoff,erm_cutoff_2sides, -2.0,-2,prt_residiual);
    else if (ext_inde_flag)  extract_inden_probes(outfileName, befileName, inde_num, ldrsq,   simu_seed, problstName, problst2exclde, genelistName,  chr,prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName, indilst2remove);
    else if(getvariance_flag || getmean_flag) getPrbVarianceMean(outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType,getvariance_flag,getmean_flag);
    else if(blup_probe_flag) blup_probe(outfileName, efileName,  befileName, transposedin,  efileType, problstName, problst2exclde,genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, probe2rm, indilstName,indilst2remove, no_fid_flag,valueType, std_thresh,upperBeta,lowerBeta,blup_indi_file);
    else if(score_flag) scoreIndividuals(outfileName,befileName, score_file, col_prb, col_score, hasHeader, phenofileName,  mpheno, problstName, problst2exclde,genelistName, chr, prbname,  fromprbname, toprbname, prbWind,fromprbkb, toprbkb,prbwindFlag, genename, probe2rm, indilstName, indilst2remove,  std_thresh, impute_mean_flag);
    else if (vqtl) V_QTL(outfileName,  efileName, befileName,  bFileName,  transposedin, efileType, problstName, problst2exclde, genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,problst2exclde, indilstName, indilst2remove,  no_fid_flag, valueType, beta2m, m2beta,  std_thresh, upperBeta, lowerBeta, dpvalfName,  thresh_det_pval,  thresh_prpt_prb,  thresh_prpt_spl,  filter_det_pval_mth,  missing_ratio_prob, autosome_num,  maf, snplstName, snplst2exclde, tsk_ttl, tsk_id,vqtl_mtd,covfileName, qcovfileName,besd_snp_major, cis_flag,  cis_itvl);
    else if(eqtl) eQTL(outfileName,  efileName, befileName,  bFileName,  transposedin, efileType, problstName, problst2exclde, genelistName,  chr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,problst2exclde, indilstName, indilst2remove,  no_fid_flag, valueType, beta2m, m2beta,  std_thresh, upperBeta, lowerBeta, dpvalfName,  thresh_det_pval,  thresh_prpt_prb,  thresh_prpt_spl,  filter_det_pval_mth,  missing_ratio_prob, autosome_num,  maf, snplstName, snplst2exclde, tsk_ttl, tsk_id,covfileName, qcovfileName,besd_snp_major);
    else if(queryBesd) query_besd(outfileName, beqtlFileName,  snplstName,  snplst2exclde,  problstName,problst2exclde,  genelistName,  pQueryBesd,  chr,   prbchr, snpchr,  snprs,  fromsnprs, tosnprs,  prbname,  fromprbname,  toprbname,snpWind,  prbWind, genename, fromsnpkb,  tosnpkb,  fromprbkb,  toprbkb,  snpwindFlag,  prbwindFlag, cis_flag,  cis_itvl,  probe2rm, snp2rm);
    else if (metaflg) {
        if(besdflstName != NULL) meta( besdflstName ,outfileName,meta_mtd,pmecs, cis_flag,  cis_itvl);
        else if(gwasflstName != NULL) meta_gwas(gwasflstName, outfileName, meta_mtd, pmecs,mecs_mth,corMatFName,snplstName, zflag);
        else {
            LOGPRINTF("ERROR: please input a file list using --besd-flist or --gwas-flist.\n");
            TERMINATE();
        }
    }
    else if(make_besd_flag) make_besd(outfileName, beqtlFileName,  snplstName,  snplst2exclde,  problstName,problst2exclde,  genelistName,  pQueryBesd,  chr,   prbchr, snpchr,  snprs,  fromsnprs, tosnprs,  prbname,  fromprbname,  toprbname,snpWind,  prbWind, genename, fromsnpkb,  tosnpkb,  fromprbkb,  toprbkb,  snpwindFlag,  prbwindFlag, cis_flag,  cis_itvl,  probe2rm, snp2rm,save_dense_flag,to_smr_flag,besd_shrink_flag, stdprb, freqFName, varFName);

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

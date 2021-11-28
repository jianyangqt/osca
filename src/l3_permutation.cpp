#include "l3_permutation.hpp"
#include <time.h>
#include <fstream>
#include <sstream>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "cis_learn_beta.hpp"

namespace PERMU
{
    static int
    trans_eQTL_num(eInfo *einfo, bInfo *bdata, int trans_itvl, vector<vector<uint32_t>> &snpids, vector<int> &trans_num)
    {
        //bp as start and gd as end
        int maxr = 0;
        trans_num.resize(einfo->_epi_include.size());
        snpids.resize(einfo->_epi_include.size());
        for (int i = 0; i < einfo->_epi_include.size(); i++)
        {
            int prbchr = einfo->_epi_chr[einfo->_epi_include[i]];
            int prbstart = einfo->_epi_bp[einfo->_epi_include[i]];
            int prbend = einfo->_epi_gd[einfo->_epi_include[i]];
            int trans_left = (prbstart - trans_itvl * 1000 > 0) ? (prbstart - trans_itvl * 1000) : 0;
            int trans_right = prbend + trans_itvl * 1000;
            int r = 0;
            for (int kk = 0; kk < bdata->_include.size(); kk++)
            {
                int snpchr = bdata->_chr[bdata->_include[kk]];
                int snpbp = bdata->_bp[bdata->_include[kk]];
                if (snpchr != prbchr || snpbp <= trans_left || snpbp >= trans_right)
                {
                    snpids[i].push_back(kk);
                    r++;
                }
            }
            trans_num[i] = r;
            if (r > maxr)
            {
                maxr = r;
            }
        }
        return maxr;
    }
    
    
    static int
    cis_eQTL_num_2(eInfo *einfo, bInfo *bdata, int cis_itvl, vector<vector<uint32_t>> &snpids, vector<int> &cis_num)
    {
        //bp as start and gd as end
        int maxr = 0;
        cis_num.resize(einfo->_epi_include.size());
        snpids.resize(einfo->_epi_include.size());
        for (int i = 0; i < einfo->_epi_include.size(); i++)
        {
            int prbchr = einfo->_epi_chr[einfo->_epi_include[i]];
            int prbstart = einfo->_epi_bp[einfo->_epi_include[i]];
            int prbend = einfo->_epi_gd[einfo->_epi_include[i]];
            int cisstart = (prbstart - cis_itvl * 1000 > 0) ? (prbstart - cis_itvl * 1000) : 0;
            int cisend = prbend + cis_itvl * 1000;
            int r = 0;
            for (int kk = 0; kk < bdata->_include.size(); kk++)
            {
                int snpchr = bdata->_chr[bdata->_include[kk]];
                int snpbp = bdata->_bp[bdata->_include[kk]];
                if (snpchr == prbchr && snpbp >= cisstart && snpbp <= cisend)
                {
                    snpids[i].push_back(kk);
                    r++;
                }
            }
            cis_num[i] = r;
            if (r > maxr)
                maxr = r;
        }
        return maxr;
    }


    static void
    read_annofile(eqtlInfo *eqtlinfo, char *bedFileName)
    {
        FILE *epifile = NULL;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum = 6;
        if (fopen_checked(&epifile, bedFileName, "r"))
            TERMINATE();
        LOGPRINTF("Reading annotation information from %s ...\n", bedFileName);
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        eqtlinfo->_snpNum = 0;
        eqtlinfo->_probNum = 0;
        eqtlinfo->_valNum = 0;
        eqtlinfo->_sampleNum = 0;
        bool chrwarning = false;
        bool genewarning = false;
        bool orienwarning = false;
        while (fgets(Tbuf, MAX_LINE_SIZE, epifile))
        {
            split_str(Tbuf, strlist, 0);
            if (Tbuf[0] == '\0')
            {
                LOGPRINTF("ERROR: Line %u is blank.\n", line_idx);
                TERMINATE();
            }
            if (strlist.size() > colnum)
            {
                //LOGPRINTF("WARNING: Line %u has more than %d items. The first %d columns would be used. \n", line_idx,colnum,colnum);
            }
            if (strlist.size() < colnum)
            {
                LOGPRINTF("ERROR: Line %u has less than %d items.\n", line_idx, colnum);
                TERMINATE();
            }
            eqtlinfo->_probe_name_map.insert(pair<string, int>(strlist[3], line_idx));
            if (eqtlinfo->_probe_name_map.size() == line_idx)
            {
                LOGPRINTF("ERROR: Duplicate probe : %s.\n", strlist[1].c_str());
                TERMINATE();
            }
            if (strlist[0] == "X" || strlist[0] == "x")
                eqtlinfo->_epi_chr.push_back(23);
            else if (strlist[0] == "Y" || strlist[0] == "y")
                eqtlinfo->_epi_chr.push_back(24);
            else if (strlist[0] == "NA" || strlist[0] == "na")
            {
                eqtlinfo->_epi_chr.push_back(-9);
                if (!chrwarning)
                {
                    LOGPRINTF("WARNING: At least one probe chr is missing.\n");
                    chrwarning = true;
                }
            }
            else if (atoi(strlist[0].c_str()) == 0)
            {
                //LOGPRINTF("ERROR: unrecongized chromosome found:\n");
                //LOGPRINTF("WARNING: unrecongized chromosome found. This chromosome is set to 0:\n");
                //LOGPRINTF("%s\n",Tbuf);
                //TERMINATE();
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            }
            else if (atoi(strlist[0].c_str()) > 24 || atoi(strlist[0].c_str()) < 0)
            {
                //LOGPRINTF("WARNING: abmormal chromosome found:\n");
                //LOGPRINTF("%s\n",Tbuf);
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));
            }
            else
                eqtlinfo->_epi_chr.push_back(atoi(strlist[0].c_str()));

            if (strlist[1] == "NA" || strlist[1] == "na" || strlist[2] == "NA" || strlist[2] == "na")
            {
                LOGPRINTF("ERROR: NA start / end position found:\n");
                LOGPRINTF("%s\n", Tbuf);
                TERMINATE();
            }
            int startp = atoi(strlist[1].c_str());
            int endp = atoi(strlist[2].c_str());
            int pos = startp + (endp - startp) / 2;
            eqtlinfo->_epi_bp.push_back(startp);
            eqtlinfo->_epi_prbID.push_back(strlist[3]);
            eqtlinfo->_epi_gd.push_back(endp);

            if (strlist[5] == "NA" || strlist[5] == "na")
            {
                if (!genewarning)
                {
                    LOGPRINTF("WARNING: at least one gene id is missing.\n");
                    genewarning = true;
                }
            }
            eqtlinfo->_epi_gene.push_back(strlist[5].c_str());
            if (strlist[4] == "NA")
            {
                eqtlinfo->_epi_orien.push_back('*');
                if (!orienwarning)
                {
                    LOGPRINTF("WARNING: At least one gene strand is missing.\n");
                    orienwarning = true;
                }
            }
            else
                eqtlinfo->_epi_orien.push_back(strlist[4][0]);
            eqtlinfo->_include.push_back(line_idx);
            line_idx++;
        }
        eqtlinfo->_probNum = line_idx;
        fclose(epifile);
        LOGPRINTF("%llu probes to be included from  %s .\n", eqtlinfo->_probNum, bedFileName);
    }

    

    static void
    gene_check(eInfo *sqtlinfo, vector<vector<int>> &tranids, char *annofileName, bInfo *bdata, eInfo *einfo, bool cis_flag, bool trans_flag)
    {
        eqtlInfo tmpinfo;
        read_annofile(&tmpinfo, annofileName);

        int ids = 0;
        map<int, int> bchr_map;
        map<int, int>::iterator iter1;
        for (int i = 0; i < bdata->_include.size(); i++)
            bchr_map.insert(pair<int, int>(bdata->_chr[bdata->_include[i]], ids++));

        map<string, int> gene_count_map;
        map<string, int>::iterator iter;
        vector<string> gene;
        vector<vector<int>> idx;
        ids = 0;
        for (int i = 0; i < einfo->_epi_include.size(); i++)
        {
            string gn = einfo->_epi_gene[einfo->_epi_include[i]];
            to_upper(gn);
            iter = gene_count_map.find(gn);
            if (iter == gene_count_map.end())
            {
                gene_count_map.insert(pair<string, int>(gn, ids++));
                gene.push_back(gn);
                vector<int> tmpidx;
                tmpidx.push_back(einfo->_epi_include[i]);
                idx.push_back(tmpidx);
            }
            else
            {
                int curid = iter->second;
                idx[curid].push_back(einfo->_epi_include[i]);
            }
        }
        vector<int> inids;
        int tmpc = 0;
        for (int i = 0; i < idx.size(); i++)
            if (idx[i].size() > 1)
            {
                iter = tmpinfo._probe_name_map.find(gene[i]);
                if (iter != tmpinfo._probe_name_map.end())
                {
                    int tid = iter->second;
                    int tchr = tmpinfo._epi_chr[tid];
                    iter1 = bchr_map.find(tchr);
                    if (cis_flag && !trans_flag)
                    {
                        if (iter1 != bchr_map.end())
                        {
                            //printf("inter cis mode\n");
                            sqtlinfo->_epi_prb.push_back(tmpinfo._epi_prbID[tid]);
                            sqtlinfo->_epi_chr.push_back(tchr);
                            sqtlinfo->_epi_gd.push_back(tmpinfo._epi_gd[tid]);
                            sqtlinfo->_epi_gene.push_back(tmpinfo._epi_gene[tid]);
                            sqtlinfo->_epi_orien.push_back(tmpinfo._epi_orien[tid]);
                            sqtlinfo->_epi_bp.push_back(tmpinfo._epi_bp[tid]);
                            sqtlinfo->_epi_include.push_back(tmpc++);
                            inids.push_back(i);
                        }
                    }
                    else if (trans_flag && !cis_flag)
                    {
                        //printf("inter trans mode\n");
                        sqtlinfo->_epi_prb.push_back(tmpinfo._epi_prbID[tid]);
                        sqtlinfo->_epi_chr.push_back(tchr);
                        sqtlinfo->_epi_gd.push_back(tmpinfo._epi_gd[tid]);
                        sqtlinfo->_epi_gene.push_back(tmpinfo._epi_gene[tid]);
                        sqtlinfo->_epi_orien.push_back(tmpinfo._epi_orien[tid]);
                        sqtlinfo->_epi_bp.push_back(tmpinfo._epi_bp[tid]);
                        sqtlinfo->_epi_include.push_back(tmpc++);
                        inids.push_back(i);
                    }
                    else
                    {
                        fprintf(stderr, "cis or trans must choosed and choose one\n");
                    }
                }
            }

        for (int i = 0; i < inids.size(); i++)
        {
            vector<int> tmp;
            tmp.swap(idx[inids[i]]);
            tranids.push_back(tmp);
        }
        LOGPRINTF("%ld genes are retained after gene check with annaotation information, genotype information, and gene expression information.\n", tranids.size());
    }


    static void
    permute_trpv(vector < vector <double> > & trpv, uint32_t nindi, uint32_t permu_count)
    {
        if (permu_count == 0) {
            return;
        }
        int i = 0, j = 0;
        vector <int> random_seq;
        for (i = 0; i < nindi; i++) {
            random_seq.push_back(i);
        }
        srand(time(NULL));
        random_shuffle(random_seq.begin(), random_seq.end());
        vector < vector<double> > tmp;
        vector< double > tmp_row;
        tmp_row.resize(nindi);
        for (i = 0; i < trpv.size(); i++) {
            for (j = 0; j < nindi; j++) {
                //printf("%le\n", trpv[i][random_seq[j]]);
                tmp_row[j] = trpv[i][random_seq[j]]; 
            }
            tmp.push_back(tmp_row);
        }

        trpv = tmp;

        return;
    }


    static double *
    ajust_pval_permu(vector <double> & p_permuted)
    {
        uint32_t i = 0, j = 0;
        uint32_t permu_times = p_permuted.size() - 1;
        double best_p = p_permuted[0];
        vector <double> best_permuted_pv;
        double sum_perm_p = 0.0, mean_perm_p = 0.0, varient_perm_p = 0.0;
        double pbml = 0, pemp = 0;
        double beta_ml1 = 0.0, beta_ml2 = 0.0;
        double beta_mm1 = 0.0, beta_mm2 = 0.0;
        for (i = 1; i < permu_times + 1; i++) {
            sum_perm_p += p_permuted[i];
            best_permuted_pv.push_back(p_permuted[i]);
        }
        mean_perm_p = sum_perm_p / permu_times;
        for (i = 1; i < permu_times + 1; i++) {
            varient_perm_p += (p_permuted[i] - mean_perm_p) * (p_permuted[i] - mean_perm_p);
        }
        varient_perm_p /= (permu_times - 1);
        if (mean_perm_p != 1.0 && varient_perm_p != 0.0) {
            beta_mm1 = mean_perm_p * (mean_perm_p * (1 - mean_perm_p) / varient_perm_p - 1);
            beta_mm2 = beta_mm1 * (1 / mean_perm_p - 1);
            beta_ml1 = beta_mm1;
            beta_ml2 = beta_mm2;
            try {
                learnBetaParameters(best_permuted_pv, beta_ml1, beta_ml2);
            } catch (const std::exception & e) {
                fprintf(stderr, "Maximum Likelihood estimation failed, use Moment Matching instead!\n");
                beta_ml1 = beta_mm1;
                beta_ml2 = beta_mm2;
            }
        } else {
            fprintf(stderr, "mean of permuted p value equal 1 or its varient equal 0\n");
        }

        uint32_t perm_p_less_p = 0;
        for (i = 1; i < permu_times + 1; i++) {
            if (abs(p_permuted[i]) <= abs(best_p)) {
                perm_p_less_p++;
            }
        }

        pemp = ((double) perm_p_less_p + 1.0) / ((double) permu_times + 1.0);
        pbml = pbeta(best_p, beta_ml1, beta_ml2, 1, 0);
        double * dt_out = (double *)malloc(sizeof(double) * 4);
        //printf("permu_pval_less_than_best_nom%u, permu_time: %u\n", perm_p_less_p, permu_times);
        //printf("beat_mm1: %le, beta_mm2: %le\n", beta_ml1, beta_mm2);
        dt_out[0] = pemp;
        dt_out[1] = pbml;
        dt_out[2] = beta_ml1;
        dt_out[3] = beta_ml2;
        return dt_out;
    }


    static void
    print_res(output_data * dataout, FILE * fout)
    {
        fprintf(stdout, "Writing result file...\n");
        fprintf(fout,
            "#probe_id\tchrom\tprobe_pos\tgene_name\toritation\tsnp_contained"
            "\tbest_snp_id\tbest_snp_chrom\tbest_spn_pos\tp_nominal"
            "\tbeta_ml1\tbeta_ml2\tp_emp\tp_bml\n");
        output_data * tmp;
        tmp = dataout;
        while (tmp) {
            fprintf(fout, 
                "%s\t%d\t%u\t%s\t%c\t%u"
                "\t%s\t%d\t%u\t%le"
                "\t%le\t%le\t%le\t%le\n", 
                tmp -> probe_id, tmp -> chrom, tmp -> probe_pos,
                tmp -> gene_name, tmp -> orientation, tmp -> snp_contained,
                tmp -> top_snp_id, tmp -> top_snp_chrom, tmp -> top_snp_pos, tmp -> p_nominal,
                tmp -> beta_ml1, tmp -> beta_ml2, tmp -> pemp, tmp -> pbml
                );
                tmp = tmp -> next;
        }
        fclose(fout);
        tmp = dataout;
        while (dataout) {
            tmp = dataout -> next;
            free(dataout);
            dataout = tmp;
        }

        return;
    }


    void
    permu_sqtl (char *outFileName, char *efileName, char *befileName,
              char *bFileName, bool transposed, int efileType, char *problstName,
              char *problst2exclde, char *genelistName, int chr, char *prbname,
              char *fromprbname, char *toprbname, int prbWind, int fromprbkb,
              int toprbkb, bool prbwindFlag, char *genename, char *probe2exclde,
              char *indilstName, char *indilst2remove, bool no_fid_flag, int valueType,
              bool beta2m, bool m2beta, double std_thresh, double upperBeta,
              double lowerBeta, char *dpvalfName, double dp_thresh, double prb_thresh,
              double spl_thresh, int filter_mth, double mssratio_prob, int autosome_num,
              double maf, char *snplstName, char *snplst2exclde, int tsk_ttl,
              int tsk_id, char *covfileName, char *qcovfileName,
              bool nofastlinear, bool cis_flag, int cis_itvl, double zeroratio, double call,
              char *annofileName, char *covbodfileName, char *covefileName,
              bool transopse_ecov, bool use_top_p, bool not_user_top,
              bool trans_flag, int trans_itvl, uint32_t permu_times)
    {
        if (!use_top_p && !not_user_top) {
            use_top_p = true;
        } else if (use_top_p && not_user_top) {
            fprintf(stderr, "use_top and not_use_top can not both be true.\n");
            fprintf(stderr, "Set use_top_p to true.\n");
            not_user_top = false;
        }
        if (cis_flag && trans_flag)
        {
            fprintf(stderr, "The cis_flag and trans_flag can not be true at same time.");
            TERMINATE();
        }
        //by default, will set cis_flag = true
        if (!cis_flag && !trans_flag)
        {
            cis_flag = true;
        }

        setNbThreads(thread_num);
        LOGPRINTF("\nUsing %d thread(s) to conduct analysis ...\n", thread_num);

        eInfo einfo;
        bInfo bdata;
        eInfo eCov;
        char * phenofileName = NULL;
        int xqtlNO = 3;
        init_einfo(&einfo);
        init_einfo(&eCov);
        load_vqtl_workspace(&einfo, &bdata, efileName, befileName, phenofileName,
                            bFileName, transposed, efileType, problstName, problst2exclde, genelistName,
                            chr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb,
                            prbwindFlag, genename, probe2exclde, indilstName, indilst2remove, no_fid_flag,
                            valueType, beta2m, m2beta, std_thresh, upperBeta, lowerBeta, dpvalfName,
                            dp_thresh, prb_thresh, spl_thresh, filter_mth, mssratio_prob, autosome_num,
                            snplstName, snplst2exclde, tsk_ttl, tsk_id, covfileName, qcovfileName,
                            NULL, xqtlNO, zeroratio, &eCov, covbodfileName, covefileName,
                            transopse_ecov); // using _keep and _eii_include, the individuals are aligned.
        if (maf > 0)
            filter_snp_maf(&bdata, maf);
        if (call > 0)
            filter_snp_call(&bdata, call);
        char outputname[FNAMESIZE];
        outputname[0] = '\0';
        if (tsk_ttl > 1)
        {
            if (outFileName != NULL)
            {
                string tmp = string(outFileName) + "_" + atos(tsk_ttl) + "_" + atos(tsk_id) + "permutation.txt"; 
                strcpy(outputname, tmp.c_str());
                outFileName = outputname;
            }
        }
        FILE * fout = fopen(outFileName, "w");
        if (!fout) {
            fputs("Open outFile failed\n", stderr);
            exit(1);
        }
        eInfo sqtlinfo;
        vector<vector<int>> tranids; //smaller einfo._epi_include
        gene_check(&sqtlinfo, tranids, annofileName, &bdata, &einfo, cis_flag, trans_flag);

        //LOGPRINTF("\nPerforming sQTL analysis ...\n");

        if (covfileName != NULL || qcovfileName != NULL)
            adjprobe(&einfo);

        vector<vector<uint32_t>> snpids; //to save _include id not _include value
        if (cis_flag)
        {
            vector<int> cis_num;
            cis_eQTL_num_2(&sqtlinfo, &bdata, cis_itvl, snpids, cis_num);
        }
        else if (trans_flag)
        {
            vector<int> trans_num;
            trans_eQTL_num(&sqtlinfo, &bdata, trans_itvl, snpids, trans_num);
        }

        bool warned = false;
        int nindi = (int)einfo._eii_include.size();
        output_data * head = NULL, * tail = NULL;
        uint32_t probe_num_ok = sqtlinfo._epi_include.size();
        int probe_couter = 0;

        #pragma omp parallel for
        for (int jj = 0; jj < sqtlinfo._epi_include.size(); jj++)
        {   
            uint32_t prb_idx = sqtlinfo._epi_include[jj];
            string gene_name = sqtlinfo._epi_gene[prb_idx];
            string prbid = sqtlinfo._epi_prb[prb_idx];
            
            //printf("%s\n", gene_name.c_str());
            /*
            if (gene_name != "UFD1" && gene_name != "MOV10L1")
            {
                continue;
            }
            */
            int prb_chr = sqtlinfo._epi_chr[prb_idx];
            char prb_ori = sqtlinfo._epi_orien[prb_idx];
            output_data * grow = (output_data *)malloc(sizeof(output_data));
            grow -> next = NULL;
            strncpy(grow -> probe_id, prbid.c_str(), 1023);
            grow -> probe_pos = sqtlinfo._epi_bp[prb_idx];
            strncpy(grow -> gene_name, gene_name.c_str(), 1023);
            grow -> chrom = prb_chr;
            grow -> orientation = prb_ori;
            // snpids should indexed by jj not prb_idx.
            int numTrans = (int)tranids[jj].size();
            uint32_t snp_num = snpids[jj].size();
            grow->snp_contained = snp_num;
            #pragma omp critical
            {
                LOGPRINTF("\n\033[0;32m>\033[0mProcessing gene:%s, probeid:%s (%d/%d)\n", 
                    gene_name.c_str(), prbid.c_str(), ++probe_couter, probe_num_ok);
                LOGPRINTF("    This gene contain %d transcritps/isoform, and have %d SNPs\n",
                        numTrans, snpids[jj].size());
               
            }
            if (snpids[jj].size() == 0)
            {
                fprintf(stderr, "gene %s was passed, because not snp contained.\n", gene_name);
                free(grow);
                continue;
            }
            MatrixXd _X;
            make_XMat(&bdata, snpids[jj], _X);
            vector<vector<double>> trpv;

            trpv.resize(numTrans);
            for (int kk = 0; kk < numTrans; kk++)
                trpv[kk].resize(nindi);
            for (int kk = 0; kk < numTrans; kk++)
            {
                for (int ll = 0; ll < nindi; ll++) {
                    trpv[kk][ll] = einfo._val[tranids[jj][kk] * einfo._eii_num \
                        + einfo._eii_include[ll]];
                }
            }

            grow -> p_nominal = 1.0;
            (grow -> top_snp_id)[0] = '\0';
            grow -> top_snp_pos = 0;
            // do permutation loop
            //uint32_t permu_times = 100;
            vector <double> p_permuted;
            p_permuted.resize(permu_times + 1);
            for (int perm_count = 0; perm_count < permu_times + 1; perm_count++) {
                p_permuted[perm_count] = 1.0;
            }
            if (thread_num == 1) {
                printf("    Start Permutation...\n");
                printf("    Permutation Count:          ");
            }
         
            int print_cnt = 0;
            for (int perm_count = 0; perm_count < permu_times + 1; perm_count++ )
            {   
                if (thread_num == 1) {
                    for (print_cnt = 0; print_cnt < 10; print_cnt++)
                    {
                        printf("\b");
                    }
                    printf("%-10d", perm_count);
                    fflush(stdout);
                }

                permute_trpv(trpv, nindi, perm_count);
                vector<vector<double>> cor_null;
                cor_null.resize(numTrans);
                for (int kk = 0; kk < numTrans; kk++)
                    cor_null[kk].resize(numTrans);

                for (int kk = 0; kk < numTrans; kk++)
                {
                    for (int ll = kk + 1; ll < numTrans; ll++)
                    {
                        vector<double> y, x;
                        for (int mm = 0; mm < einfo._eii_include.size(); mm++)
                        {
                            if (trpv[kk][mm] < 1e9 && trpv[ll][mm] < 1e9)
                            {
                                y.push_back(trpv[kk][mm]);
                                x.push_back(trpv[ll][mm]);
                            }
                        }
                        cor_null[kk][ll] = cor_null[ll][kk] = cor(y, x);
                    }
                }

                //modified by fanghl
                //remove row and columns which contain value is 1
                vector<int> need_remove;
                int i = 0, j = 0, k = 0, l = 0;
                vector<int>::iterator it;

                for (i = 0; i < cor_null.size(); i++)
                {
                    for (j = 1 + i; j < cor_null[i].size(); j++)
                    {
                        if (abs(cor_null[i][j] - 1) < 1e-15)
                        {
                            need_remove.push_back(i);
                            break;
                        }
                    }
                }

                numTrans -= need_remove.size();
                //LOGPRINTF("    %d transcripts left after filter.", numTrans);
                vector<vector<double>> cor_null_clean;
                vector<vector<double>> trpv_clean;
                vector<double> tmp;
                if (need_remove.size() > 0 && numTrans > 1)
                {
                    cor_null_clean.resize(numTrans);
                    for (i = 0; i < numTrans; i++)
                    {
                        cor_null_clean[i].resize(numTrans);
                    }

                    k = 0;
                    for (i = 0; i < cor_null.size(); i++)
                    {
                        l = 0;
                        it = find(need_remove.begin(), need_remove.end(), i);
                        if (it != need_remove.end())
                        {
                            k++;
                        }
                        else
                        {
                            for (j = 0; j < cor_null[i].size(); j++)
                            {
                                it = find(need_remove.begin(), need_remove.end(), j);
                                if (it != need_remove.end())
                                {
                                    l++;
                                }
                                else
                                {
                                    cor_null_clean[i - k][j - l] = cor_null[i][j];
                                }
                            }
                        }
                    }

                    for (int kk = 0; kk < numTrans; kk++)
                        cor_null_clean[kk][kk] = 1.0;

                    for (i = 0; i < trpv.size(); i++)
                    {
                        tmp.clear();
                        it = find(need_remove.begin(), need_remove.end(), i);
                        if (it == need_remove.end())
                        {
                            for (j = 0; j < trpv[i].size(); j++)
                            {
                                tmp.push_back(trpv[i][j]);
                            }
                            trpv_clean.push_back(tmp);
                        }
                    }
                }
                else
                {
                    cor_null_clean = cor_null;
                    for (int kk = 0; kk < numTrans; kk++)
                        cor_null_clean[kk][kk] = 1.0;
                    trpv_clean = trpv;
                }

                //_X.cols() ==snpids[jj].size()
                for (int kk = 0; kk < _X.cols(); kk++)
                {
                    uint32_t snpid = snpids[jj][kk];
                    string snprs = bdata._snp_name[bdata._include[snpid]];

                    double snpfreq = bdata._mu[bdata._include[snpid]] / 2;
                    if (snpfreq == 0 || snpfreq == 1)
                    {
                        if (!warned)
                        {
                            LOGPRINTF("WARNING: MAF found 0 or 1 with SNP(s).\n");
                            warned = 1;
                        }
                        continue;
                    }
                    vector<double> beta(numTrans), se(numTrans);
                    for (int ll = 0; ll < numTrans; ll++)
                    {
                        vector<double> y, x, rst;
                        for (int mm = 0; mm < einfo._eii_include.size(); mm++)
                        {
                            double bval = _X(mm, kk), tval = trpv_clean[ll][mm];
                            if (bval < 1e5 && tval < 1e9)
                            {
                                y.push_back(tval);
                                x.push_back(bval);
                            }
                        }
                        reg(y, x, rst);

                        beta[ll] = rst[0];
                        se[ll] = rst[1];
                    }

                    int varnum = numTrans * (numTrans - 1) / 2;
                    int i = 0, j = 0, k = 0;
                    VectorXd d(varnum), vardev(varnum), chisq_dev(varnum);
                    MatrixXd vdev(varnum, varnum);
                    MatrixXd corr_dev(varnum, varnum);
                    double beta_mean = 0;
                    double var_mean = 0;
                    double tmp1 = 0;
                    double tmp2 = 0;
                    if (use_top_p)
                    {
                        d.resize(numTrans);
                        vardev.resize(numTrans);
                        chisq_dev.resize(numTrans);
                        vdev.resize(numTrans, numTrans);
                        corr_dev.resize(numTrans, numTrans);

                        for (i = 0; i < numTrans; i++)
                        {
                            beta_mean += beta[i];
                        }
                        beta_mean = beta_mean / numTrans;

                        tmp1 = 0;
                        tmp2 = 0;
                        for (i = 0; i < numTrans; i++)
                        {
                            tmp1 += se[i] * se[i];
                        }
                        for (i = 0; i < numTrans - 1; i++)
                        {
                            for (j = i + 1; j < numTrans; j++)
                            {
                                tmp2 += 2 * se[i] * se[j] * cor_null_clean[i][j];
                            }
                        }
                        var_mean = (tmp1 + tmp2) / (numTrans * numTrans);

                        for (i = 0; i < numTrans; i++)
                        {
                            tmp1 = 0;
                            d[i] = beta[i] - beta_mean;
                            for (j = 0; j < numTrans; j++)
                            {
                                tmp1 += se[i] * se[j] * cor_null_clean[i][j];
                            }
                            vardev[i] = se[i] * se[i] + var_mean - 2 * (tmp1 / numTrans);
                        }

                        for (i = 0; i < numTrans; i++)
                        {
                            chisq_dev[i] = d[i] * d[i] / vardev[i];
                        }

                        for (i = 0; i < numTrans; i++)
                        {
                            for (j = 0; j < numTrans; j++)
                            {
                                tmp1 = 0;
                                tmp2 = 0;
                                for (k = 0; k < numTrans; k++)
                                {
                                    tmp1 += se[i] * se[k] * cor_null_clean[i][k];
                                }
                                tmp1 = tmp1 / numTrans;

                                for (k = 0; k < numTrans; k++)
                                {
                                    tmp2 += se[j] * se[k] * cor_null_clean[i][k];
                                }
                                tmp2 = tmp2 / numTrans;

                                vdev(i, j) = se[i] * se[j] * cor_null_clean[i][j] -
                                            tmp1 - tmp2 + var_mean;
                            }
                        }

                        for (i = 0; i < numTrans; i++)
                        {
                            for (j = i; j < numTrans; j++)
                            {
                                corr_dev(i, j) = corr_dev(j, i) =
                                    vdev(i, j) / sqrt(vdev(i, i) * vdev(j, j));
                            }
                        }
                    }
                    else
                    {
                        k = 0;
                        for (int m1 = 0; m1 < numTrans - 1; m1++)
                        {
                            for (int m2 = m1 + 1; m2 < numTrans; m2++)
                            {
                                d[k] = beta[m1] - beta[m2];
                                vardev[k] = se[m1] * se[m1] + se[m2] * se[m2] -
                                            2 * cor_null_clean[m1][m2] * se[m1] * se[m2];
                                k++;
                            }
                        }

                        for (int m1 = 0; m1 < varnum; m1++)
                            chisq_dev[m1] = d[m1] * d[m1] / vardev[m1];

                        int mi = 0, mj = 0;
                        for (int m1 = 0; m1 < numTrans - 1; m1++)
                        {
                            for (int m2 = m1 + 1; m2 < numTrans; m2++)
                            {
                                mj = 0;
                                for (int m3 = 0; m3 < numTrans - 1; m3++)
                                {
                                    for (int m4 = m3 + 1; m4 < numTrans; m4++)
                                    {
                                        vdev(mi, mj) = se[m1] * se[m3] * cor_null_clean[m1][m3] -
                                                    se[m1] * se[m4] * cor_null_clean[m1][m4] -
                                                    se[m2] * se[m3] * cor_null_clean[m2][m3] +
                                                    se[m2] * se[m4] * cor_null_clean[m2][m4];
                                        mj++;
                                    }
                                }
                                mi++;
                            }
                        }

                        for (int m1 = 0; m1 < varnum; m1++)
                        {
                            for (int m2 = m1; m2 < varnum; m2++)
                            {
                                corr_dev(m1, m2) = corr_dev(m2, m1) = vdev(m1, m2) / sqrt(vdev(m1, m1) *
                                                                                        vdev(m2, m2));
                            }
                        }
                    }

                    VectorXd lambda;
                    #pragma omp critical
                    {
                        SelfAdjointEigenSolver<MatrixXd> es;
                        es.compute(corr_dev, EigenvaluesOnly);
                        lambda = es.eigenvalues();
                    }

                    double z = 0.0;
                    #pragma omp critical
                    {
                        double sumChisq_dev = chisq_dev.sum();
                        double pdev = 0.0;
                        pdev = pchisqsum(sumChisq_dev, lambda);
                        p_permuted[perm_count] = (p_permuted[perm_count] > pdev)? pdev: p_permuted[perm_count];
                        if (perm_count == 0) {
                            if (grow -> p_nominal > pdev) {
                                grow -> p_nominal = pdev;
                                strcpy(grow -> top_snp_id, snprs.c_str());
                                grow -> top_snp_pos = bdata._bp[snpid];
                                grow -> top_snp_chrom = bdata._chr[snpid];
                            }
                        }
                    }

                }
            }
            double * pemp_pbml_betaml1_betaml2 = ajust_pval_permu(p_permuted);
            grow -> pemp = pemp_pbml_betaml1_betaml2[0];
            grow -> pbml = pemp_pbml_betaml1_betaml2[1];
            grow -> beta_ml1 = pemp_pbml_betaml1_betaml2[2];
            grow -> beta_ml2 = pemp_pbml_betaml1_betaml2[3];
            
            #pragma omp critical
            {
                if (head) {
                    tail -> next = grow;
                    tail = grow;
                } else {
                    head = tail = grow;
                }
            }
        }
        print_res(head, fout);
        fprintf(stdout, "Done, result was saved to %s\n", outFileName);
    }
}


//#define DEBUG_MAIN
#ifdef DEBUG_MAIN
using namespace PERMU;
int
main(int argc, char * argv[])
{
    logfile = fopen("log.txt", "w");
    unsigned char *wkspace_ua;
    uint64_t mb = 2048;
    wkspace_ua = (unsigned char *)malloc(mb * 1048576 * sizeof(char));
    memset(wkspace_ua, 0, mb * 1048576 * sizeof(char));
    uint64_t llxx = getMemSize_Plink();
    mem_left = getAllocMB_Plink(llxx);
    permu_sqtl("out1", NULL, "rint", "chr22", false, 0, NULL, NULL, NULL, -9,
        NULL, NULL, NULL, 1000, -9, -9, false, NULL, NULL, NULL, NULL, false,
        2, false, false, 0, 1, 0, NULL, 0.05, 0.01, 0.01, 0, 1, 22, 0.01, NULL,
        NULL, 1, 1, NULL, "eigenvec.txt", false, true, 2000, 0.8, 0.85,
        "anno.txt", NULL, NULL, false, false, true, false, 5000, 100);

}
#endif

/*
    Utility: OSCA;
    Module: VQTL;
    Methods: drm, svlm;

    See Module_vqtl_drm_svlm.h for more information.

    Note:
    1. library was splited out from this src file.

    2. The function Module_vqtl_drm and Module_vqtl_svlm is introspective. This
    mean you just offer argc and argv to the function. The function would decide
    whether to run or not himself. This is good to seperate the code from other
    modules and easy to maintain.


    --Benjamin Fang, 20230215
*/


#define OSCA_VQTL_DRM_SVLM_SRC
#include <ctype.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_fit.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#include "../lib/plinklite.h"
#include "../lib/bodfile.h"
#include "../lib/sysinfo.h"
#include "../lib/besdfile.h"


#define MODULE_NAME "vqtl"
#define VQTL_DRM_METHOD "drm"
#define VQTL_SVLM_METHOD "svlm"


/*->>vqtl args*/
typedef struct {
    bool flag_vqtl;
    bool flag_help;
    const char *method;

    const char *arg_geno_file;    // arg plink bed file.
    const char *opt_pheno_bod;    // option phenotype bod file.
    const char *opt_pheno_txt;    // option phenotype text file.

    int opt_thread;         // thread number.
    int opt_start_variant;  // split the job to n trunck.
    int opt_end_variant;    // job trunck to run.
    int opt_start_probe;
    int opt_end_probe;
    int opt_tast_num;
    int opt_tast_id;

    bool flag_trans;
    int opt_trans_distance_bp;
    bool flag_cis;
    int opt_cis_window_bp;

    float pthresh;
    float opt_mem;

    const char *opt_outname;  // output file name.
    const char *opt_outformat;

} VQTL_ARGS, *VQTL_ARGS_ptr;
/*-<<vqtl args*/


/*->> DRM Thread argument structure*/
typedef struct {

    int thread_num;
    int thread_index;

    uint32_t probe_slice_start_index;
    uint32_t probe_slice_len;
    uint32_t probe_offset;

    uint32_t variant_slice_start_index;
    uint32_t variant_slice_len;

    uint32_t fam_num;
    uint32_t oii_num;
    uint32_t align_len;

    char not_need_align;

    uint32_t *fam_index_array;
    uint32_t *oii_index_array;

    char *variant_data;
    uint64_t variant_data_len_char;

    //pointer of mem need to allocate.
    double *probe_data;
    double *g0_array;
    double *g1_array;
    double *g2_array;
    double *geno_array;
    double *pheno_array;
    float *result;

} DRM_THREAD_ARGS, *DRM_THREAD_ARGS_ptr;
/*<<- DRM Thread argument structure*/


/*->> SVLM Thread argument structure*/
typedef struct {
    int thread_index;
    int thread_num;

    uint32_t probe_slice_start_index;
    uint32_t probe_slice_len;
    uint32_t probe_offset;

    uint32_t variant_slice_start_index;
    uint32_t variant_slice_len;

    uint32_t fam_num;
    uint32_t oii_num;
    uint32_t align_len;

    char not_need_align;

    uint32_t *fam_index_array;
    uint32_t *oii_index_array;

    char *variant_data;
    uint64_t variant_data_len_char;

    // pointer need to allocate.
    double *probe_data;
    double *geno_array;
    double *pheno_array;
    float *result;

} SVLM_THREAD_ARGS, *SVLM_THREAD_ARGS_ptr;
/*<<- SVLM Thread argument structure*/


static void help_legacy(void);
static VQTL_ARGS_ptr vqtl_parse_args_legacy(int argc, char *argv[],
    const char *method, VQTL_ARGS_ptr args);
static void get_logfilename(int argc, char *argv[], char *fname);

static int compare_uint32(const void *a, const void *b);
static double qmedian(double *array, int p, int r, int pos);
static unsigned int BKDRHash(char *str);
static int linner_regression(const double *x, const double *y, uint32_t array_len,
    double *beta1, double *beta1_se, double *p_beta1);

static int align_fam_oii_ids(FAM_LINE_ptr fam_lines, uint32_t fam_line_num,
                             OII_LINE_ptr oii_lines, uint32_t oii_line_num,
                             uint32_t *fam_index_array,
                             uint32_t *oii_index_array, uint32_t *aligned_len);

int Module_vqtl_drm(int argc, char *argv[]);
static DRM_THREAD_ARGS_ptr make_drm_threads_args(
    int thread_num, uint32_t indi_num_fam, uint32_t indi_num_oii,
    uint32_t align_len, char not_need_align, uint32_t *fam_index_array,
    uint32_t *oii_index_array, uint32_t probe_slice_start,
    uint32_t probe_slice_len, char *variant_data, uint64_t variant_data_len,
    uint32_t variant_load_len, DRM_THREAD_ARGS_ptr thread_args);
static void free_drm_threads_args_malloc(DRM_THREAD_ARGS_ptr args, int thread_num);
static void *drm_thread_worker(void *args);

int Module_vqtl_svlm(int argc, char *argv[]);
static SVLM_THREAD_ARGS_ptr make_svlm_threads_args(
    int thread_num, uint32_t indi_num_fam, uint32_t indi_num_oii,
    uint32_t align_len, char not_need_align, uint32_t *fam_index_array,
    uint32_t *oii_index_arrary, uint32_t probe_slice_start,
    uint32_t probe_slice_len, char *variant_data, uint64_t variant_data_len,
    uint32_t variant_load_len, SVLM_THREAD_ARGS_ptr thread_args);
static void free_svlm_threads_args_malloc(SVLM_THREAD_ARGS_ptr thread_args,
    int thread_num);
static void * svlm_thread_worker(void *args);

static void write_tmp_data(void *thread_args_ori, char *args_type,
                                int thread_num, FILE *fout, float pthresh,
                                uint32_t *varint_index_pass_thresh,
                                float *beta_value, float *se_value);


int
Module_vqtl_drm(int argc, char *argv[])
{
    
    VQTL_ARGS_ptr args = (VQTL_ARGS_ptr)malloc(sizeof(VQTL_ARGS));
    args = vqtl_parse_args_legacy(argc, argv, VQTL_DRM_METHOD, args);
    if (!(args->flag_vqtl) || !(args->method)) {
        return 0;
    } else if (args->flag_help) {
        help_legacy();
        return 1;
    }
    char logfname[1024];
    get_logfilename(argc, argv, logfname);
    FILE *flog = fopen(logfname, "a");
    printf("\n\033[0;32m>>>VQTL Module DRM method\033[0m Begin\n");
    
    const char *plinkf_filename = args->arg_geno_file;
    PLINKFILE plink_data = plinkopen(plinkf_filename);
    const char *bod_filename = args->opt_pheno_bod;
    BODFILE bod_data = bodfileopen(bod_filename);

    uint32_t indi_num_fam = plink_data.individual_num;
    uint32_t indi_num_oii = bod_data.individual_num;
    uint32_t vari_num = plink_data.variant_num;
    uint32_t probe_num = bod_data.probe_num;
    
    printf("fam len: %u;\noii len: %u;\nvariant num: %u;\nprobe num: %u;\n",
        indi_num_fam, indi_num_oii, vari_num, probe_num);
    fprintf(flog, "fam len: %u;\noii len: %u;\nvariant num: %u;\nprobe num: %u;\n",
           indi_num_fam, indi_num_oii, vari_num, probe_num);
    //aligne fam and oii ids
    uint32_t *fam_index_array =
        (uint32_t *)malloc(sizeof(uint32_t) * indi_num_fam);
    uint32_t *oii_index_array =
        (uint32_t *)malloc(sizeof(uint32_t) * indi_num_oii);
    FAM_LINE_ptr fam_lines =
        (FAM_LINE_ptr)malloc(sizeof(FAM_LINE) * indi_num_fam);
    OII_LINE_ptr oii_lines =
        (OII_LINE_ptr)malloc(sizeof(OII_LINE) * indi_num_oii);
    famreadlines(&plink_data, fam_lines, indi_num_fam);
    oiireadlines(&bod_data, oii_lines, indi_num_oii);

    uint32_t align_len = 0;
    align_fam_oii_ids(fam_lines, indi_num_fam, oii_lines, indi_num_oii,
        fam_index_array, oii_index_array, &align_len);
    // if oii fam already aligned, then do not align in future.
    char not_need_align = 0;
    if ((indi_num_fam == indi_num_oii) && (align_len == indi_num_fam)) {
        not_need_align = 1;
        for (int i = 0; i < indi_num_oii; i++) {
            if ((i != oii_index_array[i]) || (i != fam_index_array[i])) {
                not_need_align = 0;
                break;
            }
        }
    }

    printf("fam and oii ids was aligned, alinged length: %u.\n", align_len);
    fprintf(flog, "fam and oii ids was aligned, alinged length: %u.\n", align_len);
#if defined VQTL_DEBUG_INFO || DEBUG_INFO
    printf("fam and oii aligned information:\n");
    for (uint32_t i = 0; i < align_len; i++) {
        uint32_t fam_index = fam_index_array[i];
        uint32_t oii_index = oii_index_array[i];
        printf("    (fam)%s-%s (oii)%s-%s\n", (fam_lines[fam_index]).family_id,
               (fam_lines[fam_index]).within_famid,
               (oii_lines[oii_index]).family_id,
               (oii_lines[oii_index]).indiv_id);
    }
#endif

    uint32_t probe_start_offset = 0;
    uint32_t probe_end = probe_num;
    uint32_t variant_start_offset = 0;
    uint32_t variant_end = vari_num;
    
    int task_num = args->opt_tast_num;
    int task_id = args->opt_tast_id;
    
    char res_fname[1024];
    if (args->opt_outname) {
        strcpy(res_fname, args->opt_outname);
    } else {
        strcpy(res_fname, "out");
    }
    if (task_id > task_num) {
        fprintf(stderr, "task id should less than task number.\n");
        return 1;
    }
    int task_len = 0;
    if (task_num > 1) {
        char sufix[1000];
        task_len = ceil((double)probe_num / task_num);
        if (task_id > 1 && task_id <= task_num) {
            sprintf(sufix, "_%d_%d", task_num, task_id);
            strcat(res_fname, sufix);
            probe_start_offset = task_len * (task_id - 1);

        } else {
            sprintf(sufix, "_%d_%d", task_num, 1);
            strcat(res_fname, sufix);
            probe_start_offset = 0;
        }
        // probe_end has been assigned as probe_num previousely.
        if (probe_start_offset + task_len < probe_num) {
            probe_end = probe_start_offset + task_len;
        }
    }
    int res_fname_len = strlen(res_fname);

    printf("start probe offset: %u;\nend probe: %u;\n"
        "start variant offset: %u;\nend variant: %u;\n", probe_start_offset,
            probe_end, variant_start_offset, variant_end);

    fprintf(flog,
        "start probe offset: %u;\nend probe: %u;\n"
        "start variant offset: %u;\nend variant: %u;\n",
        probe_start_offset, probe_end, variant_start_offset, variant_end);

    //write epi
    res_fname[res_fname_len] = '\0';
    strcat(res_fname, ".epi");
    FILE *epi_fout = fopen(res_fname, "w");
    //pass lines
    for (int i = 0; i < probe_start_offset; i++) {
        OPI_LINE opi_line;
        opireadline(&bod_data, &opi_line);
    }
    for (int i = probe_start_offset; i < probe_end; i++) {
        OPI_LINE opi_line;
        char chrom[4];
        char ori[4];
        opireadline(&bod_data, &opi_line);
        if (opi_line.chrom == 201) {
            strcpy(chrom, "X");
        } else if (opi_line.chrom == 202) {
            strcpy(chrom, "Y");
        } else if (opi_line.chrom == 203) {
            strcpy(chrom, "MT");
        } else {
            sprintf(chrom, "%d", opi_line.chrom);
        }
        if (opi_line.ori == 0) {
            strcpy(ori, "NA");
        } else if (opi_line.ori == 1) {
            strcpy(ori, "+");
        } else {
            strcpy(ori, "-");
        }
        fprintf(epi_fout, "%s\t%s\t%u\t%s\t%s\n", chrom, opi_line.probe_id,
            opi_line.position, opi_line.gene_id, ori);
    }
    fclose(epi_fout);
    printf("epi file was writen.\n");
    fprintf(flog, "epi file was writen.\n");

    //grant memory
    uint64_t mem_size = 0;
    float mem = args->opt_mem;
    if (mem == 0.0) {
        SYSINFO sysinfo_dt;
        get_sysinfo(&sysinfo_dt);
        mem_size = sysinfo_dt.mem_size_byte * (3 / 4);
    } else {
        mem_size = (uint64_t)floor((double)mem * 1024 * 1024 * 1024);
    }
    printf("granted mem size %llu\n", mem_size);
    fprintf(flog, "granted mem size %llu\n", mem_size);

    uint32_t variant_load_len = mem_size / (indi_num_fam * sizeof(char));
    uint32_t variant_total_len = variant_end - variant_start_offset;

    char *variant_data = NULL;
    uint64_t variant_data_len = 0;
    if (variant_total_len > variant_load_len) {
        variant_data_len = variant_load_len * indi_num_fam;
        variant_data =
            (char *)malloc(sizeof(char) * variant_data_len);
    } else {
        variant_load_len = variant_total_len;
        variant_data_len = variant_total_len * indi_num_fam;
        variant_data =
            (char *)malloc(sizeof(char) * variant_data_len);
    }
    printf("variant total len: %u\n", variant_total_len);
    printf("variant load len: %u\n", variant_load_len);
    printf("variant_data_len: %llu\n", variant_data_len);
    fprintf(flog, "variant total len: %u\n", variant_total_len);
    fprintf(flog, "variant load len: %u\n", variant_load_len);
    fprintf(flog, "variant_data_len: %llu\n", variant_data_len);

    int thread_num = args->opt_thread;
    pthread_t *restrict thread_ids =
        (pthread_t *)malloc(sizeof(pthread_t) * thread_num);
    DRM_THREAD_ARGS_ptr threads_args = (DRM_THREAD_ARGS_ptr)malloc(
        sizeof(DRM_THREAD_ARGS) * thread_num);

    uint32_t probe_slice_start = probe_start_offset;
    uint32_t probe_slice_len = probe_end - probe_start_offset;
    make_drm_threads_args(thread_num, indi_num_fam, indi_num_oii, align_len,
        not_need_align,fam_index_array, oii_index_array,
        probe_slice_start, probe_slice_len,
        variant_data, variant_data_len, variant_load_len,
        threads_args);

    // creat tmp directory
    char tmp_dir_name[] = "oscatmp";
    char tmp_fname[1024];
    if (access(tmp_dir_name, F_OK)) {
        mkdir(tmp_dir_name, S_IRWXU);
    } else {
        DIR *dirp = opendir(tmp_dir_name);
        struct dirent *dp = NULL;
        while ((dp = readdir(dirp)) != NULL) {
            if (strcmp(".", dp->d_name) != 0 && strcmp("..", dp->d_name) != 0) {
                strcpy(tmp_fname, tmp_dir_name);
                strcat(tmp_fname, "/");
                strcat(tmp_fname, dp->d_name);
                if (unlink(tmp_fname)) {
                    fprintf(stderr, "%s remove failed.\n", tmp_fname);
                }
            }
        }
    }

    // malloc buffer for drm_write_tmp_data funtion.
    uint32_t *variant_index_pass_thresh = (uint32_t *)malloc(sizeof(uint32_t) *
        variant_load_len);
    float *beta_value_pass_thresh = (float *)malloc(sizeof(float) * variant_load_len);
    float *se_value_pass_thresh = (float *)malloc(sizeof(float) * variant_load_len);

    uint64_t variant_data_len_real = 0;
    uint32_t variant_slice_start = 0;
    uint32_t variant_slice_len = 0;
    char tmp_fname_new[256];
    uint32_t tmp_meta_arrary[4];
    float pthresh = args->pthresh;

    res_fname[res_fname_len] = '\0';
    strcat(res_fname, ".esi");
    FILE *esi_fout = fopen(res_fname, "w");
    //pass first variant_start_offset lines.
    for (int i = 0; i < variant_start_offset; i++) {
        BIM_LINE bim_line;
        bimreadline(&plink_data, &bim_line);
    }

    for (int i = variant_start_offset; i < variant_end; i += variant_load_len) {
        variant_slice_start = i;
        if (i + variant_load_len > variant_end) {
            variant_slice_len = variant_end - i;
        } else {
            variant_slice_len = variant_load_len;
        }
        printf(">>variant slice<< %u to %u\n", variant_slice_start,
            variant_slice_start + variant_slice_len);

        strcpy(tmp_fname, tmp_dir_name);
        strcat(tmp_fname, "/");
        sprintf(tmp_fname_new, "tmp_%u_%u_%u_%u", probe_slice_start,
            probe_slice_len, variant_slice_start, variant_slice_len);
        strcat(tmp_fname, tmp_fname_new);
        tmp_meta_arrary[0] = probe_slice_start;
        tmp_meta_arrary[1] = probe_slice_len;
        tmp_meta_arrary[2] = variant_slice_start;
        tmp_meta_arrary[3] = variant_slice_len;
        FILE *fout = fopen(tmp_fname, "w");
        if (!fout) {
            fprintf(stderr, "open OSCA tmp file failed.\n");
        }
        fwrite(tmp_meta_arrary, sizeof(uint32_t), 4, fout);

        variant_data_len_real = variant_slice_len * indi_num_fam;
        bedloaddata_n(&plink_data, variant_data, variant_data_len_real,
            variant_slice_start, variant_slice_len);
        char *variant_data_ptr = variant_data;
        for (int i = 0; i < variant_slice_len; i++) {
            float first_allel_freq = 0;
            uint32_t value_num = 0;
            uint32_t first_allel_count = 0;
            for (int j = 0; j < align_len; j++) {
                uint32_t index = fam_index_array[j];
                if (variant_data_ptr[index] != 4) {
                    first_allel_count += variant_data_ptr[index];
                    value_num++;
                }
            }

            first_allel_freq = (float)first_allel_count / (2.0 * value_num);
            BIM_LINE bim_line;
            char chrom[4];
            bimreadline(&plink_data, &bim_line);
            if (bim_line.chrom == 201) {
                strcpy(chrom, "X");
            } else if (bim_line.chrom == 202) {
                strcpy(chrom, "Y");
            } else if (bim_line.chrom == 203) {
                strcpy(chrom, "MT");
            } else {
                sprintf(chrom, "%d", bim_line.chrom);
            }
            char *allel1 = NULL;
            char *allel2 = NULL;
            uint32_t line_index = plink_data.current_bim_line_index - 1;
            struct Allel_trimed_stu *trimed_allel_tmp = NULL;
            uint32_t trimed_list_len = plink_data.trimed_allel_list_len;
            if (strlen(bim_line.allel1) == 0) {
                trimed_allel_tmp = plink_data.trime_allel_list;
                for (int i = 0; i < trimed_list_len; i++) {
                    if (trimed_allel_tmp->allel_choose == 1 &&
                        trimed_allel_tmp->line_index == line_index) {
                        allel1 = trimed_allel_tmp->allel_ptr;
                        break;
                    }
                    trimed_allel_tmp = trimed_allel_tmp->next;
                }
                if (!allel1) {
                    fprintf(stderr, "can not find trimed allel.\n");
                    return 1;
                }
            } else {
                allel1 = bim_line.allel1;
            }
            if (strlen(bim_line.allel2) == 0) {
                trimed_allel_tmp = plink_data.trime_allel_list;
                for (int i = 0; i < trimed_list_len; i++) {
                    if (trimed_allel_tmp->allel_choose == 2 &&
                        trimed_allel_tmp->line_index == line_index) {
                        allel2 = trimed_allel_tmp->allel_ptr;
                        break;
                    }
                    trimed_allel_tmp = trimed_allel_tmp->next;
                }
                if (!allel2) {
                    fprintf(stderr, "can not find trimed allel.\n");
                    return 1;
                }
            } else {
                allel2 = bim_line.allel2;
            }
            fprintf(esi_fout, "%s\t%s\t%f\t%u\t%s\t%s\t%f\n", chrom, bim_line.rsid,
                bim_line.phy_pos, bim_line.pos, allel1,
                allel2, first_allel_freq);            
            variant_data_ptr += indi_num_fam;
        }


        for (int k = 0; k < thread_num; k++) {
            threads_args[k].variant_data_len_char = variant_data_len_real;
            threads_args[k].variant_slice_start_index = variant_slice_start;
            threads_args[k].variant_slice_len = variant_slice_len;
        }
        
        uint32_t j_limit = probe_end - thread_num + 1;
        uint32_t j = 0;
        
        bodfileseek(&bod_data, probe_start_offset);
        printf("progress by probe of this variant slice: %10u/%-10u", 0, probe_slice_len);
        for (j = probe_start_offset; j < j_limit; j += thread_num) {
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("%10u/%-10u", j, probe_slice_len);
            fflush(stdout);

            for (int n = 0; n < thread_num; n++) {
                double *readdata;
                uint32_t readdata_len;
                readdata = threads_args[n].probe_data;
                readdata_len = threads_args[n].oii_num;
                bodreaddata(&bod_data, readdata, readdata_len);
                threads_args[n].probe_offset = j + n;
            }

            for (int m = 0; m < thread_num; m++) {
                int pthread_status = pthread_create(&(thread_ids[m]), NULL, drm_thread_worker,
                    &(threads_args[m]));
                if (pthread_status) {
                    fprintf(stderr, "creat thread failed.\n");
                    return 1;
                }
            }
            for (int p = 0; p < thread_num; p++) {
                int join_status = pthread_join(thread_ids[p], NULL);
                if (join_status) {
                    fprintf(stderr, "thread join failed.\n");
                    return 1;
                }
            }
            
            write_tmp_data(threads_args, VQTL_DRM_METHOD, thread_num, fout,
                               pthresh, variant_index_pass_thresh,
                               beta_value_pass_thresh, se_value_pass_thresh);
        }

        int left_probe_n = 0;
        uint32_t j_keep = j;
        for (; j < probe_end; j++) {


            double *readdata;
            uint32_t readdata_len;
            readdata = threads_args[left_probe_n].probe_data;
            readdata_len = threads_args[left_probe_n].oii_num;
            bodreaddata(&bod_data, readdata, readdata_len);
            threads_args[left_probe_n].probe_offset = j;
            left_probe_n++;
        }
        for (int m = 0; m < left_probe_n; m++) {
            int pthread_status = pthread_create(&(thread_ids[m]), NULL, drm_thread_worker,
                           &(threads_args[m]));
            if (pthread_status) {
                fprintf(stderr, "creat thread failed.\n");
                return 1;
            }
        }
        for (int p = 0; p < left_probe_n; p++) {
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("%10u/%-10u", j_keep + p, probe_slice_len);
            fflush(stdout);

            int join_status = pthread_join(thread_ids[p], NULL);
            if (join_status) {
                fprintf(stderr, "thread join failed.\n");
                return 1;
            }
        }

        write_tmp_data(threads_args, VQTL_DRM_METHOD, left_probe_n, fout, pthresh,
                           variant_index_pass_thresh, beta_value_pass_thresh,
                           se_value_pass_thresh);

        fclose(fout);
        printf("\n");
    }
    fclose(esi_fout);
    printf("\nesi file was writen.\n");
    fprintf(flog, "\nesi file was writen.\n");

    //merge result into besd or plain text file.
    uint32_t tmp_file_num = (uint32_t)ceil((double)(variant_end - variant_start_offset) /
        variant_load_len);
    FILE **tmp_files_fin = (FILE **)malloc(sizeof(FILE *) * tmp_file_num);
    uint32_t tmp_file_index = 0; 
    for (int i = variant_start_offset; i < variant_end; i+=variant_load_len) {
        variant_slice_start = i;
        if (i + variant_load_len > variant_end) {
            variant_slice_len = variant_end - i;
        } else {
            variant_slice_len = variant_load_len;
        }

        strcpy(tmp_fname, tmp_dir_name);
        strcat(tmp_fname, "/");
        sprintf(tmp_fname_new, "tmp_%u_%u_%u_%u", probe_slice_start,
                probe_slice_len, variant_slice_start, variant_slice_len);
        strcat(tmp_fname, tmp_fname_new);
        //printf("tmp file name : %s\n", tmp_fname);
        FILE *fin = fopen(tmp_fname, "r");
        if (!fin) {
            fprintf(stderr, "open OSCA tmp file failed.\n");
        }
        tmp_files_fin[tmp_file_index] = fin;
        tmp_file_index++;
    }

    //use besd sparse format
    int besd_file_format = 3;
    uint32_t besd_sample_num = align_len;
    uint32_t besd_esi_num = variant_end - variant_start_offset;
    uint32_t besd_epi_num = probe_end - probe_start_offset;
    uint64_t besd_value_num = 0;
    uint64_t *besd_offset = (uint64_t *)malloc(sizeof(uint64_t) * (besd_epi_num * 2 + 1));
    besd_offset[0] = 0;
    uint32_t *besd_index = (uint32_t *)malloc(sizeof(uint32_t) * besd_esi_num);
    float *besd_beta = (float *)malloc(sizeof(float) * besd_esi_num);
    float *besd_se = (float *)malloc(sizeof(float) * besd_esi_num);
    
    for (int i = 0; i < tmp_file_num; i++) {
        uint32_t first4[4];
        fread(first4, sizeof(uint32_t), 4, tmp_files_fin[i]);
    }
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_meta");
    FILE *besd_meta_fout = fopen(tmp_fname, "w");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_index");
    FILE *besd_index_fout = fopen(tmp_fname, "w");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_beta_se");
    FILE *besd_beta_se_fout = fopen(tmp_fname, "w");
    for (int i = 0; i < besd_epi_num; i++) {
        uint32_t data_num_probe = 0;
        uint32_t *besd_index_ptr = besd_index;
        float *besd_beta_ptr = besd_beta;
        float *besd_se_ptr = besd_se;
        for (int j = 0; j < tmp_file_num; j++) {
            uint32_t data_num_file = 0;
            fread(&data_num_file, sizeof(uint32_t), 1, tmp_files_fin[j]);
            data_num_probe += data_num_file;
            fread(besd_index_ptr, sizeof(uint32_t), data_num_file, tmp_files_fin[j]);
            fread(besd_beta_ptr, sizeof(float), data_num_file, tmp_files_fin[j]);
            fread(besd_se_ptr, sizeof(float), data_num_file, tmp_files_fin[j]);
            besd_index_ptr += data_num_file;
            besd_beta_ptr += data_num_file;
            besd_se_ptr += data_num_file;
        }

        besd_value_num += 2 * data_num_probe;
        besd_offset[i * 2 + 1] = besd_offset[i * 2] + data_num_probe;
        besd_offset[i * 2 + 2] = besd_offset[i * 2 + 1] + data_num_probe;
        besd_sparse_write_variant_index(besd_index, data_num_probe, besd_index_fout);
        besd_sparse_write_beta_se_data(besd_beta, besd_se, data_num_probe, besd_beta_se_fout);
        /*
            printf("probe: %d, besd_offset: %lu, value_num: %lu\n", i,
            besd_value_num, besd_offset[i * 2 + 2]);
        */
    }
    besd_sparse_write_meta(besd_file_format, besd_sample_num, besd_esi_num,
        besd_epi_num, besd_value_num, besd_offset, besd_meta_fout);
    for (int k = 0; k < tmp_file_num; k++) {
        fclose(tmp_files_fin[k]);
    }
    fclose(besd_index_fout);
    fclose(besd_beta_se_fout);
    fclose(besd_meta_fout);

    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_meta");
    FILE *besd_meta_fin = fopen(tmp_fname, "r");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_index");
    FILE *besd_index_fin = fopen(tmp_fname, "r");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_beta_se");
    FILE *besd_beta_se_fin = fopen(tmp_fname, "r");
    res_fname[res_fname_len] = '\0';
    strcat(res_fname, ".besd");
    FILE *besd_fout = fopen(res_fname, "w");

    //here borrow variant data as a file read buffer.
    uint64_t read_len = 0;
    while (read_len = fread(variant_data, sizeof(char), variant_data_len, besd_meta_fin)) {
        fwrite(variant_data, sizeof(char), read_len, besd_fout);
    }

    while (read_len = fread(variant_data, sizeof(char), variant_data_len,
                            besd_index_fin)) {
        fwrite(variant_data, sizeof(char), read_len, besd_fout);
    }

    while (read_len = fread(variant_data, sizeof(char), variant_data_len,
                            besd_beta_se_fin)) {
        fwrite(variant_data, sizeof(char), read_len, besd_fout);
    }
    fclose(besd_meta_fin);
    fclose(besd_index_fin);
    fclose(besd_beta_se_fin);
    fclose(besd_fout);
    printf("besd file was writen.\n");
    fprintf(flog, "besd file was writen.\n");
    fclose(flog);
    //remove tmp directory.
    if (!access(tmp_dir_name, F_OK)) {
        DIR *dirp = opendir(tmp_dir_name);
        struct dirent *dp = NULL;
        while ((dp = readdir(dirp)) != NULL) {
            if (strcmp(".", dp->d_name) != 0 && strcmp("..", dp->d_name) != 0) {
                strcpy(tmp_fname, tmp_dir_name);
                strcat(tmp_fname, "/");
                strcat(tmp_fname, dp->d_name);
                if (unlink(tmp_fname)) {
                    fprintf(stderr, "%s remove failed.\n", tmp_fname);
                }
            }
        }
        rmdir(tmp_dir_name);
    }

    free(args);
    plinkclose(&plink_data);
    bodfileclose(&bod_data);
    free(oii_index_array);
    free(fam_index_array);
    free(fam_lines);
    free(oii_lines);
    free(variant_data);
    free(thread_ids);
    free_drm_threads_args_malloc(threads_args, thread_num);
    free(threads_args);
    free(variant_index_pass_thresh);
    free(beta_value_pass_thresh);
    free(se_value_pass_thresh);
    free(tmp_files_fin);
    free(besd_index);
    free(besd_beta);
    free(besd_se);
    printf("\n\033[0;32m<<<VQTL Module DRM method\033[0m End\n");
    return 1;
}


int
Module_vqtl_svlm(int argc, char *argv[])
{
    VQTL_ARGS_ptr args = (VQTL_ARGS_ptr)malloc(sizeof(VQTL_ARGS));
    args = vqtl_parse_args_legacy(argc, argv, VQTL_SVLM_METHOD, args);
    if (!(args->flag_vqtl) || !(args->method)) {
        return 0;
    } else if (args->flag_help) {
        help_legacy();
        return 1;
    }
    char logfname[1024];
    get_logfilename(argc, argv, logfname);
    FILE *flog = fopen(logfname, "a");
    printf("\n\033[0;32m>>>VQTL Module SVLM method\033[0m Begin\n");

    const char *plinkf_filename = args->arg_geno_file;
    PLINKFILE plink_data = plinkopen(plinkf_filename);
    const char *bod_filename = args->opt_pheno_bod;
    BODFILE bod_data = bodfileopen(bod_filename);

    uint32_t indi_num_fam = plink_data.individual_num;
    uint32_t indi_num_oii = bod_data.individual_num;
    uint32_t vari_num = plink_data.variant_num;
    uint32_t probe_num = bod_data.probe_num;

    printf("fam len: %u;\noii len: %u;\nvariant num: %u;\nprobe num: %u;\n",
           indi_num_fam, indi_num_oii, vari_num, probe_num);
    fprintf(flog, "fam len: %u;\noii len: %u;\nvariant num: %u;\nprobe num: %u;\n",
           indi_num_fam, indi_num_oii, vari_num, probe_num);

    // aligne fam and oii ids
    uint32_t *fam_index_array =
        (uint32_t *)malloc(sizeof(uint32_t) * indi_num_fam);
    uint32_t *oii_index_array =
        (uint32_t *)malloc(sizeof(uint32_t) * indi_num_oii);
    FAM_LINE_ptr fam_lines =
        (FAM_LINE_ptr)malloc(sizeof(FAM_LINE) * indi_num_fam);
    OII_LINE_ptr oii_lines =
        (OII_LINE_ptr)malloc(sizeof(OII_LINE) * indi_num_oii);
    famreadlines(&plink_data, fam_lines, indi_num_fam);
    oiireadlines(&bod_data, oii_lines, indi_num_oii);
    uint32_t align_len = 0;
    align_fam_oii_ids(fam_lines, indi_num_fam, oii_lines, indi_num_oii,
                      fam_index_array, oii_index_array, &align_len);
    // if oii fam already aligned, then do not align in future.
    char not_need_align = 0;
    if ((indi_num_fam == indi_num_oii) && (align_len == indi_num_fam)) {
        not_need_align = 1;
        for (int i = 0; i < indi_num_oii; i++) {
            if ((i != oii_index_array[i]) || (i != fam_index_array[i])) {
                not_need_align = 0;
                break;
            }
        }
    }
    printf("fam and oii ids was aligned, alinged length: %u\n", align_len);
    fprintf(flog, "fam and oii ids was aligned, alinged length: %u\n", align_len);

#if defined VQTL_DEBUG_INFO || DEBUG_INFO
    printf("fam and oii aligned information:\n");
    for (uint32_t i = 0; i < align_len; i++) {
        uint32_t fam_index = fam_index_array[i];
        uint32_t oii_index = oii_index_array[i];
        printf("    (fam)%s-%s (oii)%s-%s\n", (fam_lines[fam_index]).family_id,
               (fam_lines[fam_index]).within_famid,
               (oii_lines[oii_index]).family_id,
               (oii_lines[oii_index]).indiv_id);
    }
#endif

    uint32_t probe_start_offset = 0;
    uint32_t probe_end = probe_num;
    uint32_t variant_start_offset = 0;
    uint32_t variant_end = vari_num;

    int task_num = args->opt_tast_num;
    int task_id = args->opt_tast_id;
    char res_fname[1024];
    if (args->opt_outname) {
        strcpy(res_fname, args->opt_outname);
    } else {
        strcpy(res_fname, "out");
    }
    if (task_id > task_num) {
        fprintf(stderr, "task id should less than task number.\n");
        return 1;
    }
    int task_len = 0;
    if (task_num > 1) {
        char sufix[1000];
        task_len = ceil((double)probe_num / task_num);
        if (task_id > 1 && task_id <= task_num) {
            sprintf(sufix, "_%d_%d", task_num, task_id);
            strcat(res_fname, sufix);
            probe_start_offset = task_len * (task_id - 1);

        } else {
            sprintf(sufix, "_%d_%d", task_num, 1);
            strcat(res_fname, sufix);
            probe_start_offset = 0;
        }
        //probe_end has been assigned as probe_num previousely.
        if (probe_start_offset + task_len < probe_num) {
            probe_end = probe_start_offset + task_len;
        }
    }
    int res_fname_len = strlen(res_fname);

    printf("start probe offset: %u\n", probe_start_offset);
    printf("end probe: %u\n", probe_end);
    printf("start variant offset: %u\n", variant_start_offset);
    printf("end variant: %u\n", variant_end);
    fprintf(flog, "start probe offset: %u\n", probe_start_offset);
    fprintf(flog, "end probe: %u\n", probe_end);
    fprintf(flog, "start variant offset: %u\n", variant_start_offset);
    fprintf(flog, "end variant: %u\n", variant_end);
    // write epi
    res_fname[res_fname_len] = '\0';
    strcat(res_fname, ".epi");
    FILE *epi_fout = fopen(res_fname, "w");
    // pass lines
    for (int i = 0; i < probe_start_offset; i++) {
        OPI_LINE opi_line;
        opireadline(&bod_data, &opi_line);
    }
    for (int i = probe_start_offset; i < probe_end; i++) {
        OPI_LINE opi_line;
        char chrom[4];
        char ori[4];
        opireadline(&bod_data, &opi_line);
        if (opi_line.chrom == 201) {
            strcpy(chrom, "X");
        } else if (opi_line.chrom == 202) {
            strcpy(chrom, "Y");
        } else if (opi_line.chrom == 203) {
            strcpy(chrom, "MT");
        } else {
            sprintf(chrom, "%d", opi_line.chrom);
        }
        if (opi_line.ori == 0) {
            strcpy(ori, "NA");
        } else if (opi_line.ori == 1) {
            strcpy(ori, "+");
        } else {
            strcpy(ori, "-");
        }
        fprintf(epi_fout, "%s\t%s\t%u\t%s\t%s\n", chrom, opi_line.probe_id,
                opi_line.position, opi_line.gene_id, ori);
    }
    fclose(epi_fout);
    printf("epi file was writen.\n");
    fprintf(flog, "epi file was writen.\n");
    //grant memory
    uint64_t mem_size = 0;
    float mem = args->opt_mem;
    if (mem == 0.0) {
        SYSINFO sysinfo_dt;
        get_sysinfo(&sysinfo_dt);
        mem_size = sysinfo_dt.mem_size_byte * (3 / 4);
    } else {
        mem_size = (uint64_t)floor((double)mem * 1024 * 1024 * 1024);
    }
    printf("granted mem size %llu\n", mem_size);

    uint32_t variant_load_len = mem_size / (indi_num_fam * sizeof(char));
    uint32_t variant_total_len = variant_end - variant_start_offset;

    char *variant_data = NULL;
    uint64_t variant_data_len = 0;
    if (variant_total_len > variant_load_len) {
        variant_data_len = variant_load_len * indi_num_fam;
        variant_data = (char *)malloc(sizeof(char) * variant_data_len);
    } else {
        variant_load_len = variant_total_len;
        variant_data_len = variant_total_len * indi_num_fam;
        variant_data = (char *)malloc(sizeof(char) * variant_data_len);
    }
    printf("variant load len: %u\n", variant_load_len);
    printf("variant total len: %u\n", variant_total_len);
    printf("variant_data_len: %llu\n", variant_data_len);
    fprintf(flog, "variant load len: %u\n", variant_load_len);
    fprintf(flog, "variant total len: %u\n", variant_total_len);
    fprintf(flog, "variant_data_len: %llu\n", variant_data_len);

    int thread_num = args->opt_thread;
    pthread_t *restrict thread_ids =
        (pthread_t *)malloc(sizeof(pthread_t) * thread_num);
    SVLM_THREAD_ARGS_ptr threads_args =
        (SVLM_THREAD_ARGS_ptr)malloc(sizeof(DRM_THREAD_ARGS) * thread_num);

    uint32_t probe_slice_start = probe_start_offset;
    uint32_t probe_slice_len = probe_end - probe_start_offset;
    make_svlm_threads_args(thread_num, indi_num_fam, indi_num_oii, align_len,
                          not_need_align, fam_index_array, oii_index_array,
                          probe_slice_start, probe_slice_len, variant_data,
                          variant_data_len, variant_load_len, threads_args);

    // creat tmp directory
    char tmp_dir_name[] = "oscatmp";
    char tmp_fname[1024];
    if (access(tmp_dir_name, F_OK)) {
        mkdir(tmp_dir_name, S_IRWXU);
    } else {
        DIR *dirp = opendir(tmp_dir_name);
        struct dirent *dp = NULL;
        while ((dp = readdir(dirp)) != NULL) {
            if (strcmp(".", dp->d_name) != 0 && strcmp("..", dp->d_name) != 0) {
                strcpy(tmp_fname, tmp_dir_name);
                strcat(tmp_fname, "/");
                strcat(tmp_fname, dp->d_name);
                if (unlink(tmp_fname)) {
                    fprintf(stderr, "%s remove failed.\n", tmp_fname);
                }
            }
        }
    }

    // malloc buffer for drm_write_tmp_data funtion.
    uint32_t *variant_index_pass_thresh =
        (uint32_t *)malloc(sizeof(uint32_t) * variant_load_len);
    float *beta_value_pass_thresh =
        (float *)malloc(sizeof(float) * variant_load_len);
    float *se_value_pass_thresh =
        (float *)malloc(sizeof(float) * variant_load_len);

    uint64_t variant_data_len_real = 0;
    uint32_t variant_slice_start = 0;
    uint32_t variant_slice_len = 0;
    char tmp_fname_new[256];
    uint32_t tmp_meta_arrary[4];
    float pthresh = args->pthresh;

    res_fname[res_fname_len] = '\0';
    strcat(res_fname, ".esi");
    FILE *esi_fout = fopen(res_fname, "w");
    // pass first variant_start_offset lines.
    for (int i = 0; i < variant_start_offset; i++) {
        BIM_LINE bim_line;
        bimreadline(&plink_data, &bim_line);
    }

    for (int i = variant_start_offset; i < variant_end; i += variant_load_len) {
        variant_slice_start = i;
        if (i + variant_load_len > variant_end) {
            variant_slice_len = variant_end - i;
        } else {
            variant_slice_len = variant_load_len;
        }
        printf(">>variant slice<< %u to %u\n", variant_slice_start,
               variant_slice_start + variant_slice_len);

        strcpy(tmp_fname, tmp_dir_name);
        strcat(tmp_fname, "/");
        sprintf(tmp_fname_new, "tmp_%u_%u_%u_%u", probe_slice_start,
                probe_slice_len, variant_slice_start, variant_slice_len);
        strcat(tmp_fname, tmp_fname_new);
        tmp_meta_arrary[0] = probe_slice_start;
        tmp_meta_arrary[1] = probe_slice_len;
        tmp_meta_arrary[2] = variant_slice_start;
        tmp_meta_arrary[3] = variant_slice_len;
        FILE *fout = fopen(tmp_fname, "w");
        if (!fout) {
            fprintf(stderr, "open OSCA tmp file failed.\n");
        }
        fwrite(tmp_meta_arrary, sizeof(uint32_t), 4, fout);

        variant_data_len_real = variant_slice_len * indi_num_fam;
        bedloaddata_n(&plink_data, variant_data, variant_data_len_real,
                      variant_slice_start, variant_slice_len);
        char *variant_data_ptr = variant_data;
        for (int i = 0; i < variant_slice_len; i++) {
            float first_allel_freq = 0;
            uint32_t value_num = 0;
            uint32_t first_allel_count = 0;
            for (int j = 0; j < align_len; j++) {
                uint32_t index = fam_index_array[j];
                if (variant_data_ptr[index] != 4) {
                    first_allel_count += variant_data_ptr[index];
                    value_num++;
                }
            }

            first_allel_freq = (float)first_allel_count / (2.0 * value_num);
            BIM_LINE bim_line;
            char chrom[4];
            bimreadline(&plink_data, &bim_line);
            if (bim_line.chrom == 201) {
                strcpy(chrom, "X");
            } else if (bim_line.chrom == 202) {
                strcpy(chrom, "Y");
            } else if (bim_line.chrom == 203) {
                strcpy(chrom, "MT");
            } else {
                sprintf(chrom, "%d", bim_line.chrom);
            }
            fprintf(esi_fout, "%s\t%s\t%f\t%u\t%s\t%s\t%f\n", chrom,
                    bim_line.rsid, bim_line.phy_pos, bim_line.pos,
                    bim_line.allel1, bim_line.allel2, first_allel_freq);
            variant_data_ptr += indi_num_fam;
        }

        for (int k = 0; k < thread_num; k++) {
            threads_args[k].variant_data_len_char = variant_data_len_real;
            threads_args[k].variant_slice_start_index = variant_slice_start;
            threads_args[k].variant_slice_len = variant_slice_len;
        }

        uint32_t j_limit = probe_end - thread_num + 1;
        uint32_t j = 0;
        
        bodfileseek(&bod_data, probe_start_offset);
        printf("progress by probe of this variant slice: %10u/%-10u", 0,
               probe_slice_len);
        for (j = probe_start_offset; j < j_limit; j += thread_num) {
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("%10u/%-10u", j, probe_slice_len);
            fflush(stdout);
#if defined VQTL_DEBUG_INFO || DEBUG_INFO
            printf("\n");
#endif    
            for (int n = 0; n < thread_num; n++) {
                double *readdata;
                uint32_t readdata_len;
                readdata = threads_args[n].probe_data;
                readdata_len = threads_args[n].oii_num;
                bodreaddata(&bod_data, readdata, readdata_len);
                threads_args[n].probe_offset = j + n;
            }

            for (int m = 0; m < thread_num; m++) {
                int pthread_status = pthread_create(&(thread_ids[m]), NULL, svlm_thread_worker,
                               &(threads_args[m]));
                if (pthread_status) {
                    fprintf(stderr, "creat thread failed.\n");
                    return 1;
                }
            }
            for (int p = 0; p < thread_num; p++) {
                int join_status = pthread_join(thread_ids[p], NULL);
                if (join_status) {
                    fprintf(stderr, "thread join failed.\n");
                    return 1;
                }
            }
            
            write_tmp_data(threads_args, VQTL_SVLM_METHOD, thread_num, fout, pthresh,
                               variant_index_pass_thresh,
                               beta_value_pass_thresh, se_value_pass_thresh);
        }

        int left_probe_n = 0;
        uint32_t j_keep = j;
        for (; j < probe_end; j++) {
            double *readdata;
            uint32_t readdata_len;
            readdata = threads_args[left_probe_n].probe_data;
            readdata_len = threads_args[left_probe_n].oii_num;
            bodreaddata(&bod_data, readdata, readdata_len);
            threads_args[left_probe_n].probe_offset = j;
            left_probe_n++;
        }
        for (int m = 0; m < left_probe_n; m++) {
            int pthread_status = pthread_create(&(thread_ids[m]), NULL, svlm_thread_worker,
                           &(threads_args[m]));
            if (pthread_status) {
                fprintf(stderr, "cread thread failed.\n");
                return 1;
            }
        }
        for (int p = 0; p < left_probe_n; p++) {
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("%10u/%-10u", j_keep + p, probe_slice_len);
            fflush(stdout);

            int join_status = pthread_join(thread_ids[p], NULL);
            if (join_status) {
                fprintf(stderr, "thread join failed.\n");
                return 1;
            }
            
        }

        write_tmp_data(threads_args, VQTL_SVLM_METHOD, left_probe_n, fout, pthresh,
                           variant_index_pass_thresh, beta_value_pass_thresh,
                           se_value_pass_thresh);
        printf("\n");
        fclose(fout);
    }
    printf("\nesi file was writen.\n");
    fprintf(flog, "\nesi file was writen.\n");
    // merge result into besd or plain text file.
    uint32_t tmp_file_num = (uint32_t)ceil(
        (double)(variant_end - variant_start_offset) / variant_load_len);
    FILE **tmp_files_fin = (FILE **)malloc(sizeof(FILE *) * tmp_file_num);
    uint32_t tmp_file_index = 0;
    for (int i = variant_start_offset; i < variant_end; i += variant_load_len) {
        variant_slice_start = i;
        if (i + variant_load_len > variant_end) {
            variant_slice_len = variant_end - i;
        } else {
            variant_slice_len = variant_load_len;
        }

        strcpy(tmp_fname, tmp_dir_name);
        strcat(tmp_fname, "/");
        sprintf(tmp_fname_new, "tmp_%u_%u_%u_%u", probe_slice_start,
                probe_slice_len, variant_slice_start, variant_slice_len);
        strcat(tmp_fname, tmp_fname_new);
        //printf("tmp file name %s\n", tmp_fname);
        FILE *fin = fopen(tmp_fname, "r");
        if (!fin) {
            fprintf(stderr, "open OSCA tmp file failed.\n");
        }
        tmp_files_fin[tmp_file_index] = fin;
        tmp_file_index++;
    }

    int besd_file_format = 3;
    uint32_t besd_sample_num = align_len;
    uint32_t besd_esi_num = variant_end - variant_start_offset;
    uint32_t besd_epi_num = probe_end - probe_start_offset;
    uint64_t besd_value_num = 0;
    uint64_t *besd_offset =
        (uint64_t *)malloc(sizeof(uint64_t) * (besd_epi_num * 2 + 1));
    besd_offset[0] = 0;
    uint32_t *besd_index = (uint32_t *)malloc(sizeof(uint32_t) * besd_esi_num);
    float *besd_beta = (float *)malloc(sizeof(float) * besd_esi_num);
    float *besd_se = (float *)malloc(sizeof(float) * besd_esi_num);

    for (int i = 0; i < tmp_file_num; i++) {
        uint32_t first4[4];
        fread(first4, sizeof(uint32_t), 4, tmp_files_fin[i]);
    }

    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_meta");
    FILE *besd_meta_fout = fopen(tmp_fname, "w");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_index");
    FILE *besd_index_fout = fopen(tmp_fname, "w");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_beta_se");
    FILE *besd_beta_se_fout = fopen(tmp_fname, "w");
    for (int i = 0; i < besd_epi_num; i++) {
        uint32_t data_num_probe = 0;
        uint32_t *besd_index_ptr = besd_index;
        float *besd_beta_ptr = besd_beta;
        float *besd_se_ptr = besd_se;
        for (int j = 0; j < tmp_file_num; j++) {
            uint32_t data_num_file = 0;
            fread(&data_num_file, sizeof(uint32_t), 1, tmp_files_fin[j]);
            data_num_probe += data_num_file;
            fread(besd_index_ptr, sizeof(uint32_t), data_num_file,
                  tmp_files_fin[j]);
            fread(besd_beta_ptr, sizeof(float), data_num_file,
                  tmp_files_fin[j]);
            fread(besd_se_ptr, sizeof(float), data_num_file, tmp_files_fin[j]);
            besd_index_ptr += data_num_file;
            besd_beta_ptr += data_num_file;
            besd_se_ptr += data_num_file;
        }

        besd_value_num += 2 * data_num_probe;
        besd_offset[i * 2 + 1] = besd_offset[i * 2] + data_num_probe;
        besd_offset[i * 2 + 2] = besd_offset[i * 2 + 1] + data_num_probe;
        /*
            printf("probe: %d, besd_offset: %lu, value_num: %lu\n", i,
               besd_value_num, besd_offset[i * 2 + 2]);
        */

        besd_sparse_write_variant_index(besd_index, data_num_probe,
                                        besd_index_fout);
        besd_sparse_write_beta_se_data(besd_beta, besd_se, data_num_probe,
                                       besd_beta_se_fout);
    }
    besd_sparse_write_meta(besd_file_format, besd_sample_num, besd_esi_num,
                           besd_epi_num, besd_value_num, besd_offset,
                           besd_meta_fout);
    for (int k = 0; k < tmp_file_num; k++) {
        fclose(tmp_files_fin[k]);
    }
    fclose(besd_index_fout);
    fclose(besd_beta_se_fout);
    fclose(besd_meta_fout);

    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_meta");
    FILE *besd_meta_fin = fopen(tmp_fname, "r");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_index");
    FILE *besd_index_fin = fopen(tmp_fname, "r");
    strcpy(tmp_fname, tmp_dir_name);
    strcat(tmp_fname, "/");
    strcat(tmp_fname, "besd_beta_se");
    FILE *besd_beta_se_fin = fopen(tmp_fname, "r");

    res_fname[res_fname_len] = '\0';
    strcat(res_fname, ".besd");
    FILE *besd_fout = fopen(res_fname, "w");

    // here borrow variant data as a file read buffer.
    uint64_t read_len = 0;
    while (read_len = fread(variant_data, sizeof(char), variant_data_len,
                            besd_meta_fin)) {
        fwrite(variant_data, sizeof(char), read_len, besd_fout);
    }

    while (read_len = fread(variant_data, sizeof(char), variant_data_len,
                            besd_index_fin)) {
        fwrite(variant_data, sizeof(char), read_len, besd_fout);
    }

    while (read_len = fread(variant_data, sizeof(char), variant_data_len,
                            besd_beta_se_fin)) {
        fwrite(variant_data, sizeof(char), read_len, besd_fout);
    }
    fclose(besd_meta_fin);
    fclose(besd_index_fin);
    fclose(besd_beta_se_fin);
    fclose(besd_fout);
    printf("besd file was writen.\n");
    fprintf(flog, "besd file was writen.\n");

    // remove tmp directory.
    if (!access(tmp_dir_name, F_OK)) {
        DIR *dirp = opendir(tmp_dir_name);
        struct dirent *dp = NULL;
        while ((dp = readdir(dirp)) != NULL) {
            if (strcmp(".", dp->d_name) != 0 && strcmp("..", dp->d_name) != 0) {
                strcpy(tmp_fname, tmp_dir_name);
                strcat(tmp_fname, "/");
                strcat(tmp_fname, dp->d_name);
                if (unlink(tmp_fname)) {
                    fprintf(stderr, "%s remove failed.\n", tmp_fname);
                }
            }
        }
        rmdir(tmp_dir_name);
    }

    free(args);
    plinkclose(&plink_data);
    bodfileclose(&bod_data);
    free(fam_index_array);
    free(oii_index_array);
    free(fam_lines);
    free(oii_lines);
    free(variant_data);
    free(thread_ids);
    free(variant_index_pass_thresh);
    free(beta_value_pass_thresh);
    free(se_value_pass_thresh);
    free_svlm_threads_args_malloc(threads_args, thread_num);
    free(threads_args);
    free(besd_index);
    free(besd_beta);
    free(besd_se);
    printf("\n\033[0;32m<<<VQTL Module SVLM method\033[0m End\n");
    return 1;
}


static void
help_legacy(void)
{
    printf(
        "\nHelp:\n"
        "--help,        flag, print this message and exit.\n"
        "--vqtl,        flag, use vqtl module.\n"
        "--method        STR, vqtl methods. 'drm' for method DRM, 'svlm' for SVLM.\n"
        "--geno          STR, plink genotype files.\n"
        "--pheno         STR, phenotype file in plain text.(not supported yet)\n"
        "--pheno-bod     STR, phenotype files in bod file format.\n"
        "--thread-num    INT, number of threads to use, default is 1.\n"
        "--start-var     INT, index of first variant to calculate. default is 1.(not supported yet)\n"
        "--end-var       INT, last variant to calculate. default is last one (not supported yet)\n"
        "--start-probe   INT, index of first probe.(not supported yet)\n"
        "--end-probe     INT, last probe, defaulte is last probe.(not supported yet)\n"
        "--tast-num      INT, task number would like to divide.\n"
        "--tast-id       INT, task id.\n"
        "--trans,       flag, only calculate trans region.(not supported yet)\n"
        "--trans-distance-bp\n"
        "                INT, distance from probe in basepair to define as trans.(not supported yet)\n"
        "--cis           INT, only calculate cis region.(not supported yet)\n"
        "--cis-window-pb\n"
        "                INT, widow width in basepare defined as cis.(not supported yet)\n"
        "--pthresh     FLOAT, p value to filter beta1 t test result.\n"
        "--mem         FLOAT, GB, memoray used by program, default is 3/4 all memoray.\n"
        "--out           STR, output file name.\n"
        "--outformat     STR, output format, default is besd.(not supported yet)\n\n"
    );
    return;
}


static  VQTL_ARGS_ptr
vqtl_parse_args_legacy(int argc, char *argv[], const char *method, VQTL_ARGS_ptr args) 
{
    args->flag_vqtl = false;
    args->flag_help = false;
    args->method = NULL;

    args->arg_geno_file = NULL;
    args->opt_pheno_bod = NULL;
    args->opt_pheno_txt = NULL;

    args->opt_thread = 1;
    args->opt_start_variant = -1;
    args->opt_end_variant = -1;
    args->opt_start_probe = -1;
    args->opt_end_probe = -1;
    args->opt_tast_num = 1;
    args->opt_tast_id = 1;

    args->flag_trans = false;
    args->opt_trans_distance_bp = 5000000;
    args->flag_cis = false;
    args->opt_cis_window_bp = 2000000;

    args->pthresh = 0.5;

    args->opt_outname = NULL;
    args->opt_outformat = NULL;
    args->opt_mem = 0.0;


    char out_tast_suffix[64];
    int  run_this_method = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--vqtl") == 0) {
            run_this_method++;
        }
        if (strcmp(argv[i], "--method") == 0) {
            if (i + 1 < argc && strcmp(argv[i + 1], method) == 0) {
                run_this_method++;
            }
        } 
    }
    char logfname[1024];
    get_logfilename(argc, argv, logfname);
    FILE *flog = fopen(logfname, "a");
    if (run_this_method == 2) {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--vqtl") == 0) {
                args->flag_vqtl = true;
                printf("--vqtl\n");
                fprintf(flog, "--vqtl\n");
                continue;
            }

            if (strcmp(argv[i], "--help") == 0) {
                args->flag_help = true;
                break;
            }

            if (strcmp(argv[i], "--method") == 0) {
                if ( i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->method = argv[++i];
                    printf("--method %s\n", argv[i]);
                    fprintf(flog, "--method %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--method need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }


            if (strcmp(argv[i], "--geno") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->arg_geno_file = argv[++i];
                    printf("--geno %s\n", argv[i]);
                    fprintf(flog, "--geno %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--geno need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--pheno") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_pheno_txt = argv[++i];
                    printf("--pheno %s\n", argv[i]);
                    fprintf(flog, "--pheno %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--pheno need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--pheno-bod") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_pheno_bod = argv[++i];
                    printf("--pheno-bod %s\n", argv[i]);
                    fprintf(flog, "--pheno-bod %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--pheno-bod need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--thread-num") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_thread = atoi(argv[++i]);
                    printf("--thread-num %s\n", argv[i]);
                    fprintf(flog, "--thread-num %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--thread-num need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--start-var") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_start_variant = atoi(argv[++i]);
                    printf("--start-var %s\n", argv[i]);
                    fprintf(flog, "--start-var %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--start-var need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--end-var") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_end_variant = atoi(argv[++i]);
                    printf("--end-var %s\n", argv[i]);
                    fprintf(flog, "--end-var %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--end-var need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--start-probe") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_start_probe = atoi(argv[++i]);
                    printf("--start-probe %s\n", argv[i]);
                    fprintf(flog, "--start-probe %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--start-probe need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--end-probe") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_end_probe = atoi(argv[++i]);
                    printf("--end-probe %s\n", argv[i]);
                    fprintf(flog, "--end-probe %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--end-probe need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--task-num") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_tast_num = atoi(argv[++i]);
                    printf("--task-num %s\n", argv[i]);
                    fprintf(flog, "--task-num %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--task-num need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--task-id") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_tast_id = atoi(argv[++i]);
                    printf("--task-id %s\n", argv[i]);
                    fprintf(flog, "--task-id %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--task-id need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--trans") == 0) {
                args->flag_trans = true;
                printf("--trans\n");
                fprintf(flog, "--trans");
                continue;
            }

            if (strcmp(argv[i], "--trans-distance-bp") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_trans_distance_bp = atoi(argv[++i]);
                    printf("--trans-distance-bp %s\n", argv[i]);
                    fprintf(flog, "--trans-distance-bp %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--trans-distance-bp need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--cis") == 0) {
                args->flag_cis = true;
                printf("--cis\n");
                fprintf(flog, "--cis\n");
                continue;
            }

            if (strcmp(argv[i], "--cis-window-bp") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_cis_window_bp = atoi(argv[++i]);
                    printf("--cis-window-bp %s\n", argv[i]);
                    fprintf(flog, "--cis-window-bp %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--cis-window-bp need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--pthresh") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->pthresh = atof(argv[++i]);
                    printf("--pthresh %s\n", argv[i]);
                    fprintf(flog, "--pthresh %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--pthresh need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--mem") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_mem = atof(argv[++i]);
                    printf("--mem %s\n", argv[i]);
                    fprintf(flog, "--mem %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--mem need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--out") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_outname = argv[++i];
                    printf("--out %s\n", argv[i]);
                    fprintf(flog, "--out %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--out need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--outformat") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_outformat = argv[++i];
                    printf("--outformat %s\n", argv[i]);
                    fprintf(flog, "--outformat %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--outformat need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            fprintf(stderr, "option %s not recgnized.\n", argv[i]);
            fprintf(flog, "option %s not recgnized.\n", argv[i]);
            args->flag_help = true;
            break;
        }
    }
    fclose(flog);
    return args;
}


static int
compare_uint32(const void *a, const void *b)
{
    return (*(uint32_t *)a - *(uint32_t *)b);
}


/*
     BKDR Hash Function
 */
static unsigned int
BKDRHash(char *str)
{
    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int hash = 0;
    while (*str)
    {
        hash = hash * seed + (*str++);
    }
    return (hash & 0x7FFFFFFF);
}


static int
linner_regression(const double *x, const double *y, uint32_t array_len, double *beta1,
    double *se_beta1, double *p_beta1)
{
    //printf("%lf %lf %u\n", x[0], y[0], array_len);
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(x, 1, y, 1, array_len, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

/*
    double stdev0 = sqrt(cov00);
    double t0 = c0 / stdev0;
    double pv0 = (t0 < 0)? 2 * (1 - gsl_cdf_tdist_P(-t0, array_len - 2)): 2 *
        (1 - gsl_cdf_tdist_P(t0, array_len - 2));
*/
    double stdev1 = sqrt(cov11);
    double t1 = c1 / stdev1;
    double pv1 = t1 < 0 ? 2 * (1 - gsl_cdf_tdist_P(-t1, array_len - 2)): 2 *
        (1 - gsl_cdf_tdist_P(t1, array_len - 2));
/*
    int i = 0;
    double dl = array_len - 2;
    double y_mean = 0;
    for (i = 0; i < array_len; i++) {
        y_mean += y[i];
    }
    y_mean = y_mean / array_len;

    double y_var = 0;
    for (i = 0; i < array_len; i++) {
        y_var += pow(y[i] - y_mean, 2);
    }

    double R2 = 1 - sumsq / y_var;
    double F = R2 * dl / (1 - R2);
    double p_value = 1 - gsl_cdf_fdist_P(F, 1, dl);
*/
    //printf("%le %le %le %le\n", c1, stdev1, t1, p_value);
    *beta1 = c1;
    *se_beta1 = stdev1;
    *p_beta1 = pv1;

    return 0;
}


static int
align_fam_oii_ids(FAM_LINE_ptr fam_lines, uint32_t fam_line_num,
    OII_LINE_ptr oii_lines, uint32_t oii_line_num, uint32_t *fam_index_array,
    uint32_t *oii_index_array, uint32_t *aligned_len)
{
    struct HASH_NODE {
        int32_t data_index;
        struct HASH_NODE *next;
    };

    struct HASH_NODE *hash_array = (struct HASH_NODE *)malloc(
        sizeof(struct HASH_NODE) * fam_line_num);
    int32_t i = 0, j = 0;
    for (i = 0; i < fam_line_num; i++) {
        hash_array[i].data_index = -1;
        hash_array[i].next = NULL;
    }

    uint32_t hash_value = 0;
    char str1[128];
    char str2[128];
    for (i = 0; i < fam_line_num; i++) {
        strcpy(str1, fam_lines[i].family_id);
        strcat(str1, fam_lines[i].within_famid);
        hash_value = BKDRHash(str1);
        hash_value %= fam_line_num;
        if (hash_array[hash_value].data_index == -1) {
            hash_array[hash_value].data_index = i;
        } else {
            struct HASH_NODE *hash_bucket_tail = NULL;
            struct HASH_NODE *new_hash_node = (struct HASH_NODE *)
                malloc(sizeof(struct HASH_NODE));
            new_hash_node->data_index = i;
            new_hash_node->next = NULL;
            hash_bucket_tail = &(hash_array[hash_value]);
            while(hash_bucket_tail->next) {
                hash_bucket_tail = hash_bucket_tail->next;
            }
            hash_bucket_tail->next = new_hash_node;
        }
    }

    // remove repeat
    struct HASH_NODE *this_ptr = NULL, *next_ptr = NULL;
    uint32_t duplicat_num = 0;

    for (i = 0; i < fam_line_num; i++) {
        this_ptr = &(hash_array[i]);
        while (this_ptr) {
            if (this_ptr->data_index != -1) {
                strcpy(str1, fam_lines[this_ptr->data_index].family_id);
                strcat(str1, fam_lines[this_ptr->data_index].within_famid);
                next_ptr = this_ptr->next;
                while (next_ptr && next_ptr->data_index != -1) {
                    strcpy(str2, fam_lines[next_ptr->data_index].family_id);
                    strcat(str2, fam_lines[next_ptr->data_index].within_famid);
                    if (strcmp(str1, str2) == 0) {
                        duplicat_num++;
                        next_ptr->data_index = -1;
                    }
                    next_ptr = next_ptr->next;
                }
            }
            this_ptr = this_ptr->next;
        }
    }
    if (duplicat_num > 0) {
        fprintf(stderr, "duplicated fam-indi id found\n");
    }

    // algin oii index
    this_ptr = NULL;
    int32_t oii_mapped_index_to_fam = -1;
    int32_t *oii_mapped_index_array_tmp = (int32_t *)malloc(sizeof(int32_t) *
        oii_line_num);
    for (i = 0; i < oii_line_num; i++) {
        oii_mapped_index_to_fam = -1;
        strcpy(str1, oii_lines[i].family_id);
        strcat(str1, oii_lines[i].indiv_id);
        hash_value = BKDRHash(str1);
        hash_value %= fam_line_num;
        this_ptr = &(hash_array[hash_value]);
        while (this_ptr) {
            if (this_ptr->data_index != -1) {
                strcpy(str2, fam_lines[this_ptr->data_index].family_id);
                strcat(str2, fam_lines[this_ptr->data_index].within_famid);
                if (strcmp(str1, str2) == 0) {
                    oii_mapped_index_to_fam = this_ptr->data_index;
                    break;
                }
            }
            this_ptr = this_ptr->next;
        }
        oii_mapped_index_array_tmp[i] = oii_mapped_index_to_fam;
    }

    //remove duplication of oii individuals.
    for (i = 0; i < oii_line_num; i++) {
        if (oii_mapped_index_array_tmp[i] != -1) {
            for (j = i + 1; j < oii_line_num; j++) {
                if (oii_mapped_index_array_tmp[i] == oii_mapped_index_array_tmp[j]) {
                    oii_mapped_index_array_tmp[j] = -1;
                }
            }
        }
    }

    //copy oii index to new array and count oii aligned length.
    uint32_t oii_index_mapped_num = 0;
    for (i = 0; i < oii_line_num; i++) {
        
        if (oii_mapped_index_array_tmp[i] != -1) {
            oii_index_array[oii_index_mapped_num] = i;
            fam_index_array[oii_index_mapped_num] = oii_mapped_index_array_tmp[i];
            oii_index_mapped_num++;
        }
    }
    free(oii_mapped_index_array_tmp);

    *aligned_len = oii_index_mapped_num;
    return 0;
}


static DRM_THREAD_ARGS_ptr
make_drm_threads_args(
    int thread_num, uint32_t indi_num_fam, uint32_t indi_num_oii,
    uint32_t align_len, char not_need_align, uint32_t *fam_index_array,
    uint32_t *oii_index_array,
    uint32_t probe_slice_start, uint32_t probe_slice_len,
    char *variant_data, uint64_t variant_data_len, uint32_t variant_load_len,
    DRM_THREAD_ARGS_ptr thread_args) {
        
    for (int i = 0; i < thread_num; i++) {
        thread_args[i].thread_index = i;
        thread_args[i].thread_num = thread_num;

        thread_args[i].fam_num = indi_num_fam;
        thread_args[i].oii_num = indi_num_oii;
        thread_args[i].align_len = align_len;
        thread_args[i].not_need_align = not_need_align;
        thread_args[i].fam_index_array = fam_index_array;
        thread_args[i].oii_index_array = oii_index_array;

        thread_args[i].probe_slice_start_index = probe_slice_start;
        thread_args[i].probe_slice_len = probe_slice_len;

        thread_args[i].variant_data = variant_data;
        thread_args[i].variant_data_len_char = variant_data_len;

        //printf("%u %u %u %u\n", indi_num_oii, align_len, vari_num, variant_load_num);
        //all need test malloc exit status.
        thread_args[i].probe_data = (double *)malloc(sizeof(double) * indi_num_oii);
        thread_args[i].g0_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].g1_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].g2_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].geno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].pheno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].result = (float *)malloc(sizeof(float) * variant_load_len * 3);

        if (!thread_args[i].probe_data ||
            !thread_args[i].g0_array || !thread_args[i].g1_array ||
            !thread_args[i].g2_array || !thread_args[i].geno_array ||
            !thread_args[i].pheno_array || !thread_args[i].result) {
                fprintf(stderr, "allocat for threads args failed.\n");
                return NULL;
        }
    }
    
    return thread_args;
}


static void
get_logfilename(int argc, char *argv[], char *fname)
{
    int task_num = -1;
    int task_id = -1;
    char *outname = NULL;
    for (int i = 0; i < argc; i++) {
        if (strcmp("--task-num", argv[i]) == 0) {
            task_num = atoi(argv[++i]);
            continue;
        }
        if (strcmp("--task-id", argv[i]) == 0) {
            task_id = atoi(argv[++i]);
            continue;
        }
        if (strcmp("--out", argv[i]) == 0) {
            outname = argv[++i];
            continue;
        }
    }
    if (task_num > 0 && task_id > 0) {
        sprintf(fname, "%s_%d_%d.log", outname, task_num, task_id);
    } else if (task_num > 0) {
        sprintf(fname, "%s_%d_1.log", outname, task_num);
    } else if (outname){
        sprintf(fname, "%s.log", outname);
    } else {
        sprintf(fname, "%s", "osca.log");
    }
    return;
}


static void
free_drm_threads_args_malloc(DRM_THREAD_ARGS_ptr thread_args, int thread_num)
{
    for (int i = 0; i < thread_num; i++) {
        if (thread_args[i].fam_index_array) {
            thread_args[i].fam_index_array = NULL;
        }

        if (thread_args[i].oii_index_array) {
            thread_args[i].oii_index_array = NULL;
        }
        
        if (thread_args[i].variant_data) {
            thread_args[i].variant_data = NULL;
        }

        if (thread_args[i].probe_data) {
            free(thread_args[i].probe_data);
            thread_args[i].probe_data = NULL;
        }

        if (thread_args[i].g0_array) {
            free(thread_args[i].g0_array);
            thread_args[i].g0_array = NULL;
        }

        if (thread_args[i].g1_array) {
            free(thread_args[i].g1_array);
            thread_args[i].g1_array = NULL;
        }

        if (thread_args[i].g2_array) {
            free(thread_args[i].g2_array);
            thread_args[i].g2_array = NULL;
        }

        if (thread_args[i].geno_array) {
            free(thread_args[i].geno_array);
            thread_args[i].geno_array = NULL;
        }

        if (thread_args[i].pheno_array) {
            free(thread_args[i].pheno_array);
            thread_args[i].pheno_array = NULL;
        }

        if (thread_args[i].result) {
            free(thread_args[i].result);
            thread_args[i].result = NULL;
        }

    }
    return;
}


static double
qmedian(double *array, int p, int r, int pos)
{
    double x = array[r];
    int i = p - 1;
    double tmp;
    for (int j = p; j < r; j++) {
        if (array[j] <= x) {
            i++;
            tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }
    }
    array[r] = array[i + 1];
    array[i + 1] = x;

    if (pos > i + 1) {
        return qmedian(array, i + 2, r, pos);
    } else if (pos < i + 1) {
        return qmedian(array, p, i, pos);
    } else {
        return array[i + 1];
    }
}


static void *
drm_thread_worker(void *args) {
    DRM_THREAD_ARGS_ptr args_in = (DRM_THREAD_ARGS_ptr)args;

#if defined VQTL_DEBUG_INFO || DEBUG_INFO
    printf(">%d\n", args_in->thread_index);
    clock_t t1 = clock();
#endif

    uint32_t variant_slice_len = args_in->variant_slice_len;

    uint32_t fam_num = args_in->fam_num;
    uint32_t align_len = args_in->align_len;
    char not_need_align = args_in->not_need_align;

    uint32_t *fam_index_array = args_in->fam_index_array;
    uint32_t *oii_index_array = args_in->oii_index_array;

    char *variant_data_loaded = args_in->variant_data;

    double *pheno_data = args_in->probe_data;
    double *g0_array = args_in->g0_array;
    double *g1_array = args_in->g1_array;
    double *g2_array = args_in->g2_array;

    double *geno_data_aligned = args_in->geno_array;
    double *pheno_data_aligned = args_in->pheno_array;
    float *result = args_in->result;

    for (int i = 0; i < variant_slice_len; i++) {
        uint32_t align_len_rm_missing = 0;
        char *current_geno_one = variant_data_loaded + i * fam_num;
        uint32_t g0_num = 0, g1_num = 0, g2_num = 0;

        for (uint32_t j = 0; j < align_len; j++) {
            uint32_t geno_index = 0;
            uint32_t pheno_index = 0;
            char geno_value;
            double pheno_value;
            geno_index = fam_index_array[j];
            geno_value = current_geno_one[geno_index];
            pheno_index = oii_index_array[j];
            pheno_value = pheno_data[pheno_index];
           
            // do I need use other way to compare float number?
            if (geno_value != 4 && pheno_value != -9.0) {
                if (geno_value == 0) {
                    g0_array[g0_num] = pheno_value;
                    g0_num++;
                } else if (geno_value == 1) {
                    g1_array[g1_num] = pheno_value;
                    g1_num++;
                } else if (geno_value == 2) {
                    g2_array[g2_num] = pheno_value;
                    g2_num++;
                } else {
                    fprintf(stderr, "geno value can not be others.\n");
                    return NULL;
                }

                align_len_rm_missing++;
            }
        }

        /*new method to get median.*/
        if (g0_num > 0) {
            double geno_0_median = 0.0;
            if (g0_num % 2) {
                geno_0_median = qmedian(g0_array, 0, g0_num - 1, g0_num / 2);
            } else {
                geno_0_median =
                    (qmedian(g0_array, 0, g0_num - 1, g0_num / 2 - 1) +
                     qmedian(g0_array, 0, g0_num - 1, g0_num / 2)) /
                    2;
            }
            uint32_t copy_len = 0;
            for (uint32_t i = 0; i < g0_num; i++) {
                double tmp = 0;
                geno_data_aligned[copy_len] = (double)0.0;
                tmp = g0_array[i] - geno_0_median;
                pheno_data_aligned[copy_len] = (tmp < 0)? -tmp: tmp;
                copy_len++;
            }
        }

        if (g1_num > 0) {
            double geno_1_median = 0.0;
            if (g1_num % 2) {
                geno_1_median = qmedian(g1_array, 0, g1_num - 1, g1_num / 2);
            } else {
                geno_1_median =
                    (qmedian(g1_array, 0, g1_num - 1, g1_num / 2 - 1) +
                     qmedian(g1_array, 0, g1_num - 1, g1_num / 2)) /
                    2;
            }
            uint32_t copy_len = g0_num;
            for (uint32_t i = 0; i < g1_num; i++) {
                double tmp = 0.0;
                geno_data_aligned[copy_len] = (double)1.0;
                tmp = g1_array[i] - geno_1_median;
                pheno_data_aligned[copy_len] = (tmp < 0)? -tmp: tmp;
                copy_len++;
            }
        }
        if (g2_num > 0) {
            double geno_2_median = 0.0;
            if (g2_num % 2) {
                geno_2_median = qmedian(g2_array, 0, g2_num - 1, g2_num / 2);
            } else {
                geno_2_median =
                    (qmedian(g2_array, 0, g2_num - 1, g2_num / 2 - 1) +
                     qmedian(g2_array, 0, g2_num - 1, g2_num / 2)) /
                    2;
            }
            uint32_t copy_len = g0_num + g1_num;
            for (uint32_t i = 0; i < g2_num; i++) {
                double tmp = 0.0;
                geno_data_aligned[copy_len] = (double)2.0;
                tmp = g2_array[i] - geno_2_median;
                pheno_data_aligned[copy_len] = (tmp < 0)? -tmp: tmp;
                copy_len++;
            }
        }

        double beta1 = 0.0, se_beta1 = 0.0, p_beta1 = 0.0;
        double beta0 = 0.0, se_beta0 = 0.0;
        double cov01, sumsq;
        gsl_fit_linear(geno_data_aligned, 1, pheno_data_aligned, 1,
            align_len_rm_missing, &beta0, &beta1, &se_beta0, &cov01, &se_beta1,
            &sumsq);
        se_beta1 = sqrt(se_beta1);
        double t1 = beta1 / se_beta1;
        p_beta1 = t1 < 0 ? 2 * (1 - gsl_cdf_tdist_P(-t1, align_len_rm_missing - 2))
                         : 2 * (1 - gsl_cdf_tdist_P(t1, align_len_rm_missing - 2));
        
        result[i * 3] = (float)beta1;
        result[i * 3 + 1] = (float)se_beta1;
        result[i * 3 + 2] = (float)p_beta1;
    }

#if defined VQTL_DEBUG_INFO || DEBUG_INFO
    clock_t t2 = clock();
    printf("<%d   %ld cpu ticks\n", args_in->thread_index, t2 - t1);
#endif

    return NULL;
}


static SVLM_THREAD_ARGS_ptr
make_svlm_threads_args(
    int thread_num, uint32_t indi_num_fam, uint32_t indi_num_oii,
    uint32_t align_len, char not_need_align, uint32_t *fam_index_array,
    uint32_t *oii_index_arrary, uint32_t probe_slice_start,
    uint32_t probe_slice_len, char *variant_data, uint64_t variant_data_len,
    uint32_t variant_load_len, SVLM_THREAD_ARGS_ptr thread_args) {

    for (int i = 0; i < thread_num; i++) {
        thread_args[i].thread_index = i;
        thread_args[i].thread_num = thread_num;
        thread_args[i].probe_slice_start_index = probe_slice_start;
        thread_args[i].probe_slice_len = probe_slice_len;
        
        // will assigned later
        thread_args[i].probe_offset = 0;
        thread_args[i].variant_slice_start_index = 0;
        thread_args[i].variant_slice_len = 0;

        thread_args[i].fam_num = indi_num_fam;
        thread_args[i].oii_num = indi_num_oii;
        thread_args[i].align_len = align_len;

        thread_args[i].not_need_align = not_need_align;

        thread_args[i].fam_index_array = fam_index_array;
        thread_args[i].oii_index_array = oii_index_arrary;

        thread_args[i].variant_data = variant_data;
        thread_args[i].variant_data_len_char = variant_data_len;

        thread_args[i].probe_data = (double *)malloc(sizeof(double) * indi_num_oii);
        thread_args[i].geno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].pheno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].result = (float *)malloc(sizeof(float) * variant_load_len * 3);
    }
    return thread_args;
}


static void
free_svlm_threads_args_malloc(SVLM_THREAD_ARGS_ptr thread_args, int thread_num)
{
    for (int i = 0; i < thread_num; i++) {
        if (thread_args[i].fam_index_array) {
            thread_args[i].fam_index_array = NULL;
        }

        if (thread_args[i].oii_index_array) {
            thread_args[i].oii_index_array = NULL;
        }

        if (thread_args[i].variant_data) {
            thread_args[i].variant_data = NULL;
        }

        if (thread_args[i].probe_data) {
            free(thread_args[i].probe_data);
            thread_args[i].probe_data = NULL;
        }

        if (thread_args[i].geno_array) {
            free(thread_args[i].geno_array);
            thread_args[i].geno_array = NULL;
        }

        if (thread_args[i].pheno_array) {
            free(thread_args[i].pheno_array);
            thread_args[i].pheno_array = NULL;
        }

        if (thread_args[i].result) {
            free(thread_args[i].result);
            thread_args[i].result = NULL;
        }
    }

    return;
}


static void *
svlm_thread_worker(void *args)
{
    SVLM_THREAD_ARGS_ptr args_in = (SVLM_THREAD_ARGS_ptr)args;

#if defined VQTL_DEBUG_INFO || DEBUG_INFO
    printf(">%d\n", args_in->thread_index);
    clock_t t1 = clock();
#endif

    uint32_t fam_num = args_in->fam_num;
    uint32_t align_len = args_in->align_len;
    uint32_t variant_slice_len = args_in->variant_slice_len;
    char *variant_data = args_in->variant_data;
    char not_need_align = args_in->not_need_align;

    uint32_t *fam_index_array = args_in->fam_index_array;
    uint32_t *oii_index_array = args_in->oii_index_array;

    double *probe_data = args_in->probe_data;
    double *geno_array = args_in->geno_array;
    double *pheno_array = args_in->pheno_array;
    float *result = args_in->result;

    for (int i = 0; i < variant_slice_len; i++) {
        char *varinat_data_one = variant_data + (i * fam_num);
        uint32_t align_len_rm_missing = 0;
        uint32_t geno_index = 0, pheno_index = 0;
        char geno_value = 0;
        double pheno_value = 0.0;

        for (int j = 0; j < align_len; j++) {

            geno_index = fam_index_array[j];
            pheno_index = oii_index_array[j];
            geno_value = varinat_data_one[geno_index];
            pheno_value = probe_data[pheno_index];
            
            if (geno_value != 4 && pheno_value != -9.0) {
                geno_array[align_len_rm_missing] = (double)geno_value;
                pheno_array[align_len_rm_missing] = pheno_value;
                align_len_rm_missing++;
            }            
        }

        double beta1 = 0.0, se_beta1 = 0.0, p_beta1 = 0.0;
        double beta0 = 0.0, se_beta0 = 0.0;
        double cov01, sumsq;

        gsl_fit_linear(geno_array, 1, pheno_array, 1, align_len_rm_missing,
            &beta0, &beta1, &se_beta0, &cov01, &se_beta1, &sumsq);

        // pheno_array turn into array of residual square.
        for (int k = 0; k < align_len_rm_missing; k++) {
            double residule;
            residule = pheno_array[k] - (beta0 + beta1 * geno_array[k]);
            pheno_array[k] = residule * residule;
        }
        gsl_fit_linear(geno_array, 1, pheno_array, 1, align_len_rm_missing,
                       &beta0, &beta1, &se_beta0, &cov01, &se_beta1, &sumsq);

        se_beta1 = sqrt(se_beta1);
        double t1 = beta1 / se_beta1;
        p_beta1 = t1 < 0
                      ? 2 * (1 - gsl_cdf_tdist_P(-t1, align_len_rm_missing - 2))
                      : 2 * (1 - gsl_cdf_tdist_P(t1, align_len_rm_missing - 2));

        result[i * 3] = (float)beta1;
        result[i * 3 + 1] = (float)se_beta1;
        result[i * 3 + 2] = (float)p_beta1;
        
    }

#if defined VQTL_DEBUG_INFO || DEBUG_INFO
    clock_t t2 = clock();
    printf("<%d    %lu cpu ticks\n", args_in->thread_index, t2 - t1);
#endif 

    return NULL;
}


static void
write_tmp_data(void *thread_args_ori, char *args_type,
    int thread_num, FILE *fout, float pthresh,
    uint32_t *varint_index_pass_thresh,
    float *beta_value, float *se_value)
{

    if (strcmp(args_type, VQTL_DRM_METHOD) == 0) {
        DRM_THREAD_ARGS_ptr thread_args = (DRM_THREAD_ARGS_ptr)thread_args_ori;
        for (int i = 0; i < thread_num; i++) {
            uint32_t variant_num_pass_thresh = 0;
            float *result = thread_args[i].result;
            uint32_t variant_slice_len = thread_args[i].variant_slice_len;
            uint32_t variant_slice_start =
                thread_args[i].variant_slice_start_index;
            for (int j = 0; j < variant_slice_len; j++) {
                uint32_t offset = j * 3;
                if (result[offset + 2] <= pthresh) {
                    varint_index_pass_thresh[variant_num_pass_thresh] =
                        j + variant_slice_start;
                    beta_value[variant_num_pass_thresh] = result[offset];
                    se_value[variant_num_pass_thresh] = result[offset + 1];
                    variant_num_pass_thresh++;
                }
            }
            /*
            printf("probe %u: variant %u pass thresh.\n", 
                thread_args[i].probe_offset, variant_num_pass_thresh);
            */
            fwrite(&variant_num_pass_thresh, sizeof(uint32_t), 1, fout);
            fwrite(varint_index_pass_thresh, sizeof(uint32_t),
                   variant_num_pass_thresh, fout);
            fwrite(beta_value, sizeof(float), variant_num_pass_thresh, fout);
            fwrite(se_value, sizeof(float), variant_num_pass_thresh, fout);
        }
    } else if (strcmp(args_type, VQTL_SVLM_METHOD) == 0) {
        SVLM_THREAD_ARGS_ptr thread_args = (SVLM_THREAD_ARGS_ptr)thread_args_ori;
        for (int i = 0; i < thread_num; i++) {
            uint32_t variant_num_pass_thresh = 0;
            float *result = thread_args[i].result;
            uint32_t variant_slice_len = thread_args[i].variant_slice_len;
            uint32_t variant_slice_start =
                thread_args[i].variant_slice_start_index;
            for (int j = 0; j < variant_slice_len; j++) {
                uint32_t offset = j * 3;
                if (result[offset + 2] <= pthresh) {
                    varint_index_pass_thresh[variant_num_pass_thresh] =
                        j + variant_slice_start;
                    beta_value[variant_num_pass_thresh] = result[offset];
                    se_value[variant_num_pass_thresh] = result[offset + 1];
                    variant_num_pass_thresh++;
                }
            }
            /*
            printf("probe %u: variant %u pass thresh.\n",
                   thread_args[i].probe_offset, variant_num_pass_thresh);
            */
            fwrite(&variant_num_pass_thresh, sizeof(uint32_t), 1, fout);
            fwrite(varint_index_pass_thresh, sizeof(uint32_t),
                   variant_num_pass_thresh, fout);
            fwrite(beta_value, sizeof(float), variant_num_pass_thresh, fout);
            fwrite(se_value, sizeof(float), variant_num_pass_thresh, fout);
        }
    } else {
        fprintf(stderr, "method type not recognized.\n");
        return;
    }

    return;
    
}
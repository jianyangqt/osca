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

#include "../lib/plinklite.h"
#include "../lib/bodfile.h"
#include "../lib/sysinfo.h"


/*
#define MODULE_LEVER_DB_INFO
#define METHOD_LEVER_DB_INFO
#define LEVER_3_DB_INFO
#define LEVER_4_DB_INFO
*/

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
    int opt_mem;

    const char *opt_outname;  // output file name.
    const char *opt_outformat;

} VQTL_ARGS, *VQTL_ARGS_ptr;
/*-<<vqtl args*/


/*->> DRM Thread argument structure*/
typedef struct {
    int thread_index;
    int thread_num;

    uint32_t fam_num;
    uint32_t oii_num;
    uint32_t vari_num;
    uint32_t probe_num;
    uint32_t align_len;
    char not_need_align;

    uint32_t *fam_index_array;
    uint32_t *oii_index_array;

    uint32_t variant_load_num;
    char *variant_data;
    uint64_t variant_data_len_char;

    //pointer of mem need to allocate.
    double *probe_data;
    double *g0_array;
    double *g1_array;
    double *g2_array;
    double *geno_array;
    double *pheno_array;
    double *result;

} DRM_THREAD_ARGS, *DRM_THREAD_ARGS_ptr;
/*<<- DRM Thread argument structure*/

/*->> SVLM Thread argument structure*/
typedef struct {
    int thread_index;

    uint32_t fam_num;
    uint32_t oii_num;
    uint32_t probe_num;
    uint32_t align_len;

    // pointer to memory have alloced.
    double *bod_data_all_probe;
    uint32_t *fam_index_array;
    uint32_t *oii_index_arrya;

    // pointer need to allocate.
    double *geno_data_this_variant;
    double *geno_array;
    double *pheno_array;
    double *result;

} SVLM_THREAD_ARGS, *SVLM_THREAD_ARGS_ptr;
/*<<- SVLM Thread argument structure*/


static void help_legacy(void);
static VQTL_ARGS_ptr vqtl_parse_args_legacy(int argc, char *argv[], const char *method,
                                     VQTL_ARGS_ptr args);


static int compare_uint32(const void *a, const void *b);
static int compare_double(const void *a, const void *b);
static unsigned int BKDRHash(char *str);
static int linner_regression2(const double *x, const double *y,
    const uint32_t array_len,
    double *beta0, double *beta0_se, double *beta0_p,
    double *beta1, double *beta1_se, double *beta1_p,
    double *ftest_p_val);
static int linner_regression(const double *x, const double *y, uint32_t array_len, double *c1_res,
    double *stdev1_res, double *t1_res, double *p_value_res);

static int align_fam_oii_ids(FAM_LINE_ptr fam_lines, uint32_t fam_line_num,
                             OII_LINE_ptr oii_lines, uint32_t oii_line_num,
                             uint32_t *fam_index_array,
                             uint32_t *oii_index_array, uint32_t *aligned_len);


int Module_vqtl_drm(int argc, char *argv[]);
static DRM_THREAD_ARGS_ptr make_drm_threads_args(
    int thread_num, uint32_t indi_num_fam, uint32_t indi_num_oii,
    uint32_t vari_num, uint32_t probe_num, uint32_t align_len, char not_need_align, uint32_t variant_load_num,
    char *variant_data, uint64_t variant_data_len, uint32_t *fam_index_array,
    uint32_t *oii_index_array, DRM_THREAD_ARGS_ptr thread_args);
static void free_drm_threads_args_malloc(DRM_THREAD_ARGS_ptr args, int thread_num);
static void *drm_thread_worker(void *args);
static void drm_print_res(DRM_THREAD_ARGS_ptr thread_args, int thread_num,
                      int probe_num, FILE *outfile);


int Module_vqtl_svlm(int argc, char *argv[]);
static SVLM_THREAD_ARGS_ptr make_svlm_threads_args(int thread_num,
    uint32_t fam_node_num, uint32_t oii_node_num,
    uint32_t probe_num, uint32_t align_len, double *bod_data_all,
    uint32_t *fam_index_array, uint32_t *oii_index_array,
    SVLM_THREAD_ARGS_ptr thread_args);
static void free_svlm_threads_args_malloc(SVLM_THREAD_ARGS_ptr thread_args,
    int thread_num);
static void * svlm_thread_worker(void *args);
static void svlm_print_res(SVLM_THREAD_ARGS_ptr thread_args, int thread_num,
    int probe_num, FILE *outfile);


enum PROGRESS {
    BEGIN,
    ARG_PARSED,
    GENO_FILE_MALLOC,
    BIM_READED,
    FAM_READED,
    PHENO_FILE_MALLOC,
    OPI_READED,
    OII_READED,
    FAM_OII_INDEX_ARRAY_MALLOC,
    BOD_DATA_INITED,
    ALL_BOD_DATA_LOADED,
    BED_DATA_INITED,
    GET_VAR_REANGE,
    THREAD_ID_MALLOC,
    THREAD_ARGS_MALLOC,
    THREAD_ARGS_MEMBER_MALLOC,
    OUT_FIEL_OPENED,
    END
};


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

    printf("\n\033[0;32m>>>VQTL Module DRM method\033[0m Begin\n");
    
    const char *plinkf_filename = args->arg_geno_file;
    PLINKFILE plink_data = plinkopen(plinkf_filename);
    const char *bod_filename = args->opt_pheno_bod;
    BODFILE bod_data = bodfileopen(bod_filename);

    uint32_t indi_num_fam = plink_data.individual_num;
    uint32_t indi_num_oii = bod_data.individual_num;
    uint32_t vari_num = plink_data.variant_num;
    uint32_t probe_num = bod_data.probe_num;
    printf("fam len: %u\noii len: %u\nvariant num: %u\nprobe num: %u\n",
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

#ifdef DEBUG_INFO
    for (int i = 0; i < indi_num_fam; i++) {
        printf("fam_line: %s %s\n", fam_lines[i].family_id, fam_lines[i].within_famid);
    }
    for (int i = 0; i < indi_num_oii; i++) {
        printf("oii_line: %s %s\n", oii_lines[i].family_id, oii_lines[i].indiv_id);
    }
#endif

    uint32_t align_len = 0;
    align_fam_oii_ids(fam_lines, indi_num_fam, oii_lines, indi_num_oii,
        fam_index_array, oii_index_array, &align_len);
    printf("aligne len: %u\n", align_len);
    
    char not_need_align = 0; //if oii fam already aligned, then do not align in future.
    if ((indi_num_fam == indi_num_oii) && (align_len == indi_num_fam)) {
        not_need_align = 1;
        for (int i = 0; i < indi_num_oii; i++) {
            if ((i != oii_index_array[i]) || (i != fam_index_array[i])) {
                not_need_align = 0;
                break;
            }
        }
    }
    printf("not need align: %d\n", not_need_align);

#ifdef DEBUG_INFO
    for (int i = 0; i < align_len; i++) {
        printf("%u\t%u\n", oii_index_array[i], fam_index_array[i]);
    }
#endif

    uint32_t probe_start_offset = 0;
    uint32_t probe_end = probe_num;
    uint32_t variant_start_offset = 0;
    uint32_t variant_end = vari_num;
    
    int task_num = args->opt_tast_num;
    int task_id = args->opt_tast_id;
    if (task_id > task_num) {
        fprintf(stderr, "task id should less than task number.\n");
        return 1;
    }
    int task_len = 0;
    if (task_num > 1) {
        task_len = ceil((double)probe_num / task_num);
        if (task_id > 1 && task_id <= task_num) {
            probe_start_offset = task_len * (task_id - 1);

        } else {
            probe_start_offset = 0;
        }

        if (probe_start_offset + task_len < probe_num) {
            probe_end = probe_start_offset + task_len;
        }

    }
    printf("start probe offset: %u\n", probe_start_offset);
    printf("end probe: %u\n", probe_end);
    printf("start variant offset: %u\n", variant_start_offset);
    printf("end variant: %u\n", variant_end);


    uint64_t mem_size = 0;
    int mem = args->opt_mem;
    printf("--mem %d\n", mem);
    if (mem == 0) {
        SYSINFO sysinfo_dt;
        get_sysinfo(&sysinfo_dt);
        mem_size = sysinfo_dt.mem_size_byte * (3 / 4);
    } else {
        mem_size = (uint64_t)mem * 1024 * 1024 * 1024;
    }
    printf("mem size %lu\n", mem_size);

 

    uint32_t variant_load_len = 0;
    variant_load_len = mem_size / (indi_num_fam * sizeof(char));
    uint32_t variant_process_len = variant_end - variant_start_offset;
    printf("variant_load_len %u; variant_process_len %u\n", variant_load_len,
        variant_process_len);

    char *variant_data = NULL;
    uint64_t variant_data_len = 0;
    if (variant_process_len > variant_load_len) {
        variant_data_len = variant_load_len * indi_num_fam;
        variant_data =
            (char *)malloc(sizeof(char) * variant_data_len);
    } else {
        variant_data_len = variant_process_len * indi_num_fam;
        variant_data =
            (char *)malloc(sizeof(char) * variant_data_len);
    }

    printf("%lu\n", variant_data_len);

    uint64_t variant_data_len_real = 0;
    uint32_t variant_load_len_real = 0;
    int thread = args->opt_thread;
    pthread_t *restrict thread_ids =
        (pthread_t *)malloc(sizeof(pthread_t) * thread);
    DRM_THREAD_ARGS_ptr threads_args = (DRM_THREAD_ARGS_ptr)malloc(
        sizeof(DRM_THREAD_ARGS) * thread);
    make_drm_threads_args(thread, indi_num_fam, indi_num_oii, vari_num,
                          probe_num, align_len, not_need_align, variant_load_len, variant_data,
                          variant_data_len, fam_index_array, oii_index_array,
                          threads_args);


    for (int i = variant_start_offset; i < variant_end; i += variant_load_len) {
        if (i + variant_load_len > variant_end) {
            variant_load_len_real = variant_end - i;
        } else {
            variant_load_len_real = variant_load_len;
        }
        variant_data_len_real = variant_load_len_real * indi_num_fam;
        bedloaddata_n(&plink_data, variant_data, variant_data_len_real, i, variant_load_len_real);
        for (int k = 0; k < thread; k++) {
            threads_args->variant_load_num = variant_load_len_real;
            threads_args->variant_data_len_char = variant_data_len_real;
        }
        
        uint32_t j_limit = probe_end - thread + 1;
        int j = 0;
        bodfileseek(&bod_data, probe_start_offset);
        for (j = probe_start_offset; j < j_limit; j += thread) {
            double *readdata;
            uint32_t readdata_len;
            for (int n = 0; n < thread; n++) {
                readdata = threads_args[n].probe_data;
                readdata_len = threads_args[n].oii_num;
                bodreaddata(&bod_data, readdata, readdata_len);
            }

            for (int m = 0; m < thread; m++) {
                pthread_create(&(thread_ids[m]), NULL, drm_thread_worker, &(threads_args[m]));
            }
            for (int p = 0; p < thread; p++) {
                pthread_join(thread_ids[p], NULL);
            }

        }

        int left_probe_n = 0;
        for (; j < probe_end; j++) {
            left_probe_n++;
            double *readdata;
            uint32_t readdata_len;
            readdata = threads_args[j].probe_data;
            readdata_len = threads_args[j].oii_num;
            bodreaddata(&bod_data, readdata, readdata_len);
        }
        for (int m = 0; m < left_probe_n; m++) {
            pthread_create(&(thread_ids[m]), NULL, drm_thread_worker,
                           &(threads_args[m]));
        }
        for (int p = 0; p < left_probe_n; p++) {
            pthread_join(thread_ids[p], NULL);
        }
    }

    plinkclose(&plink_data);
    bodfileclose(&bod_data);
    free(oii_index_array);
    free(fam_index_array);
    free(fam_lines);
    free(oii_lines);

    return 1;
}


int
Module_vqtl_svlm(int argc, char *argv[])
{
/*
    enum PROGRESS svlm_progress = BEGIN;
    VQTL_ARGS_ptr args = (VQTL_ARGS_ptr)malloc(sizeof(VQTL_ARGS));
    args = vqtl_parse_args_legacy(argc, argv, VQTL_SVLM_METHOD, args);
    if (!(args->flag_vqtl) || !(args->method)) {
        return 0;
    } else if (args->flag_help) {
        help_legacy();
        return 1;
    }
    svlm_progress = ARG_PARSED;
   

    printf("\n\033[0;32m>>>VQTL Module SVLM method\033[0m Begin\n");

    int geno_name_len = strlen(args->arg_geno_file);
    char *geno_file = (char *)malloc(sizeof(char) * (geno_name_len + 5));
    if (!geno_file) {
        fprintf(stderr, "malloc for geno file name failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = GENO_FILE_MALLOC;

    strcpy(geno_file, args->arg_geno_file);
    strcat(geno_file, ".bim");
    BIM_DATA_ptr bim_data_res = read_bim_file(geno_file);
    if (!bim_data_res) {
        fprintf(stderr, "read bim failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = BIM_READED;

    geno_file[geno_name_len] = '\0';
    strcat(geno_file, ".fam");
    FAM_DATA_ptr fam_data_res = read_fam_file(geno_file);
    if (!fam_data_res) {
        fprintf(stderr, "read fam file failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = FAM_READED;
    
    int pheno_name_len = strlen(args->opt_pheno_bod);
    char *pheno_file = (char *)malloc(sizeof(char) * (pheno_name_len + 5));
    if (!pheno_file) {
        fprintf(stderr, "malloc for pheno file failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = PHENO_FILE_MALLOC;

    strcpy(pheno_file, args->opt_pheno_bod);
    strcat(pheno_file, ".opi");
    OPI_DATA_ptr opi_data_res = read_opi_file(pheno_file);
    if (!opi_data_res) {
        fprintf(stderr, "read pheno file failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = OPI_READED;

    pheno_file[pheno_name_len] = '\0';
    strcat(pheno_file, ".oii");
    OII_DATA_ptr oii_data_res = read_oii_file(pheno_file);
    if (!oii_data_res) {
        fprintf(stderr, "read oii file failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = OII_READED;

    FAM_NODE_ptr *fam_node_array = (FAM_NODE_ptr *)make_node_ptr_array(
        (void *)fam_data_res, FAM_DATA_TYPE);
    if (!fam_node_array) {
        fprintf(stderr, "malloc for fam node array failed.\n");
        goto clean_before_exit;
    }

    OII_NODE_ptr *oii_node_array = (OII_NODE_ptr *)make_node_ptr_array(
        (void *)oii_data_res, OII_DATA_TYPE);
    if (!oii_node_array) {
        fprintf(stderr, "malloc for oii node arrary failed.\n");
        goto clean_before_exit;
    }


    uint32_t fam_node_num = fam_data_res->line_num;
    uint32_t oii_node_num = oii_data_res->line_num;
    uint32_t align_len = 0;
    uint32_t *fam_index_array = NULL, *oii_index_array = NULL;
    align_fam_oii_ids(fam_node_array, fam_node_num, oii_node_array,
                      oii_node_num, &fam_index_array, &oii_index_array,
                      &align_len);
    if (!fam_index_array || !oii_index_array) {
        fprintf(stderr, "make fam index array or oii index array failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = FAM_OII_INDEX_ARRAY_MALLOC;


    // load bod data to memory
    pheno_file[pheno_name_len] = '\0';
    strcat(pheno_file, ".bod");
    BOD_DATA_ptr bod_data_res = init_bod_data_struct(pheno_file);
    if (!bod_data_res) {
        fprintf(stderr, "init bod data failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = BOD_DATA_INITED;

    double *bod_data_all = load_bod_to_mem(bod_data_res);
    if (!bod_data_all) {
        fprintf(stderr, "load all bod data into mem failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = ALL_BOD_DATA_LOADED;

    geno_file[geno_name_len] = '\0';
    strcat(geno_file, ".bed");
    BED_DATA_ptr bed_data_res = init_bed_data_struct(
        geno_file, fam_data_res->line_num, bim_data_res->line_num);
    if (!bed_data_res) {
        fprintf(stderr, "init bed data failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = BED_DATA_INITED;

    uint32_t variant_num_bim = bim_data_res->line_num;
    uint32_t probe_num_opi = opi_data_res->line_num;
    int thread_num = args->opt_thread;
    int start_variant = 0;
    int end_variant = variant_num_bim;

    int task_num = args->opt_tast_num;
    int task_id = args->opt_tast_id;
    if (task_num > 1) {
        int task_width = ceil(variant_num_bim / task_num);
        if (task_id != task_num) {
            start_variant = (task_id - 1) * task_width;
            end_variant = task_id * task_width;
        } else {
            start_variant = (task_id - 1) * task_width;
            end_variant = variant_num_bim;
        }
    }
    if (args->opt_start_variant > 0) {
        start_variant = args->opt_start_variant - 1;
    }
    if (args->opt_end_variant > 0) {
        end_variant = args->opt_end_variant;
    }
    if (start_variant >= end_variant) {
        fprintf(
            stderr,
            "start variant index should less equal than end variant index.");
        goto clean_before_exit;
    }
    svlm_progress = GET_VAR_REANGE;


    pthread_t *restrict thread_ids =
        (pthread_t *)malloc(sizeof(pthread_t) * thread_num);
    if (!thread_ids) {
        fprintf(stderr, "malloc for threads id failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = THREAD_ID_MALLOC;

    uint32_t limit = end_variant - thread_num + 1;
    SVLM_THREAD_ARGS_ptr thread_args = 
        (SVLM_THREAD_ARGS_ptr)malloc(sizeof(SVLM_THREAD_ARGS) * thread_num);
    if (!thread_args) {
        fprintf(stderr, "malloc for theads args failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = THREAD_ARGS_MALLOC;

    SVLM_THREAD_ARGS_ptr status = make_svlm_threads_args(thread_num, fam_node_num, oii_node_num, probe_num_opi,
                      align_len, bod_data_all, fam_index_array, oii_index_array,
                      thread_args);
    if (!status) {
        fprintf(stderr, "malloc for thread args inner memory failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = THREAD_ARGS_MEMBER_MALLOC;
    
    char outname[1024];
    if (!args->opt_outname) {
        strcpy(outname, "out");
    } else {
        strcpy(outname, args->opt_outname);
    }
    FILE *outfile = fopen(outname, "w");
    if (!outfile) {
        fprintf(stderr, "open out file failed.\n");
        goto clean_before_exit;
    }
    svlm_progress = OUT_FIEL_OPENED;

    uint32_t i = 0;
    char *geno_buf_by_varian = NULL;
    for (i = start_variant; i < limit; i += thread_num ) {
        for (int j = 0; j < thread_num; j++) {
            geno_buf_by_varian = read_bed_by_variant(bed_data_res, NULL);
            for (uint32_t k = 0; k < fam_node_num; k++) {
                ((thread_args[j]).geno_data_this_variant)[k] =
                    (double)geno_buf_by_varian[k];
            }
        }

        for (uint32_t l = 0; l < thread_num; l++) {
            pthread_create(&(thread_ids[l]), NULL, svlm_thread_worker,
                           &(thread_args[l]));
        }
        for (uint32_t m = 0; m < thread_num; m++) {
            pthread_join(thread_ids[m], NULL);
        }
        svlm_print_res(thread_args, thread_num, probe_num_opi, outfile);
    }

    for (; i < end_variant; i++) {
        geno_buf_by_varian = read_bed_by_variant(bed_data_res, NULL);
        for (int j = 0; j < fam_node_num; j++) {
            (thread_args[0]).geno_data_this_variant[j] =
                (double)geno_buf_by_varian[j];
        }
        svlm_thread_worker(thread_args);
        svlm_print_res(thread_args, 1, probe_num_opi, outfile);
    }

    svlm_progress = END;

    
clean_before_exit:
    switch (svlm_progress) {
        case END:;
        case OUT_FIEL_OPENED:
            fclose(outfile);
        case THREAD_ARGS_MEMBER_MALLOC:
            free_svlm_threads_args_malloc(thread_args, thread_num);
        case THREAD_ARGS_MALLOC:
            free(thread_args);
        case THREAD_ID_MALLOC:
            free(thread_ids);
        case GET_VAR_REANGE:;
        case BED_DATA_INITED:
            finalize_bed_data_struct(bed_data_res);
        case ALL_BOD_DATA_LOADED:
            clean_load_bod(bod_data_all);
        case BOD_DATA_INITED:
            finalize_bod_data_struct(bod_data_res);
        case FAM_OII_INDEX_ARRAY_MALLOC:
            clean_fam_oii_align_array(oii_index_array, fam_index_array);
        case OII_READED:
            clean_oii_data(oii_data_res);
        case OPI_READED:
            clean_opi_data(opi_data_res);
        case PHENO_FILE_MALLOC:
            free(pheno_file);
        case FAM_READED:
            clean_fam_data(fam_data_res);
        case BIM_READED:
            clean_bim_data(bim_data_res);
        case GENO_FILE_MALLOC:
            free(geno_file);
        case ARG_PARSED:
            free(args);
        default:
            if (svlm_progress == END) {
                printf("Run successfully.\n");
            }
            printf("\033[0;32m<<<VQTL Module SVLM method\033[0m Returned\n\n");
            
    }
*/
    return 1;
}


static void
help_legacy(void)
{
    printf(
        "\nHelp:\n"
        "--help,        flag, print this help message.\n"
        "--vqtl,        flag, use vqtl module.\n"
        "--method        STR, vqtl method. 'drm' for method DRM, 'svlm' for SVLM.\n"
        "--geno          STR, plink genotype files.\n"
        "--pheno         STR, phenotyp file in plain text.(no support yet)\n"
        "--pheno-bod     STR, phenotyp file in bod file format.\n"
        "--thread-num    INT, number of threads to use, default is 1.\n"
        "--start-var     INT, index of first variant to calculate. default is 1.(not support yet)\n"
        "--end-var       INT, last variant to calculate. default is last one of bim "
                             "file.(not support yet)\n"
        "--start-probe   INT, index of first probe.(not support yet)\n"
        "--end-probe     INT, last probe, defaulte is last one.(not support yet)\n"
        "--tast-num      INT, tast number.\n"
        "--tast-id       INT, task id.\n"
        "--trans,       flag, only calculate trans region.(not support yet)\n"
        "--trans-distance-bp\n"
        "                INT, distance from probe in basepair to define as trans.(not support yet)\n"
        "--cis           INT, only calculate cis region.(not support yet)\n"
        "--cis-window-pb\n"
        "                INT, widow width in basepare defined as cis.(not support yet)\n"
        "--pthresh     FLOAT, p value to filter beta1 t test result.\n"
        "--mem           INT, GB, memoray used by program, default is 3/4 all memoray.\n"
        "--out           STR, output file name.\n"
        "--outformat     STR, output format, default is besd.(not support yet)\n\n"
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
    args->opt_mem = 0;


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

    if (run_this_method == 2) {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--vqtl") == 0) {
                args->flag_vqtl = true;
                printf("--vqtl\n");
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
                    printf("--pheno-bod %s", argv[i]);
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
                    continue;
                } else {
                    fprintf(stderr, "--task-id need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--trans") == 0) {
                args->flag_trans = true;
                printf("--trans");
                continue;
            }

            if (strcmp(argv[i], "--trans-distance-bp") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_trans_distance_bp = atoi(argv[++i]);
                    printf("--trans-distance-bp %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--trans-distance-bp need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--cis") == 0) {
                args->flag_cis = true;
                printf("--cis");
                continue;
            }

            if (strcmp(argv[i], "--cis-window-bp") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_cis_window_bp = atoi(argv[++i]);
                    printf("--cis-window-bp %s\n", argv[i]);
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
                    continue;
                } else {
                    fprintf(stderr, "--pthresh need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            if (strcmp(argv[i], "--mem") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_mem = atoi(argv[++i]);
                    printf("--mem %s\n", argv[i]);
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
                    continue;
                } else {
                    fprintf(stderr, "--outformat need a argument.\n");
                    args->flag_help = true;
                    break;
                }
            }

            fprintf(stderr, "option %s not recgnized.\n", argv[i]);
            args->flag_help = true;
            break;
        }
    }
    return args;
}


static int
compare_uint32(const void *a, const void *b)
{
    return (*(uint32_t *)a - *(uint32_t *)b);
}


static int
compare_double(const void *a, const void *b)
{
    return ((*((double *)a) - *((double *)b)) > 0) ? 1: -1;
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
linner_regression(const double *x, const double *y, uint32_t array_len, double *c1_res,
    double *stdev1_res, double *t1_res, double *p_value_res)
{
    //printf("%lf %lf %u\n", x[0], y[0], array_len);
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(x, 1, y, 1, array_len, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    //printf("%lf %lf %lf %lf %lf %lf\n", c0, c1, cov00, cov01, cov11, sumsq);

    // double stdev0 = sqrt(cov00);
    // double t0 = c0 / stdev0;
    // double pv0 = (t0 < 0)? 2 * (1 - gsl_cdf_tdist_P(-t0, n - 2)): 2 * (1 -
    // gsl_cdf_tdist_P(t0, n - 2));

    double stdev1 = sqrt(cov11);
    double t1 = c1 / stdev1;
    // double pv1 = t1 < 0 ? 2 * (1 - gsl_cdf_tdist_P(-t1, n - 2)): 2 * (1 -
    // gsl_cdf_tdist_P(t1, n - 2));

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

    // double ym = 0.2 * (y[0] + y[1] + y[2] + y[3] + y[4]);
    // double sct = pow(y[0] - ym, 2) + pow(y[1]-ym, 2) + pow(y[2] - ym, 2) +
    // pow(y[3] - ym, 2) + pow(y[4] - ym, 2);
    double R2 = 1 - sumsq / y_var;
    double F = R2 * dl / (1 - R2);
    double p_value = 1 - gsl_cdf_fdist_P(F, 1, dl);

    //printf("%le %le %le %le\n", c1, stdev1, t1, p_value);
    *c1_res = c1;
    *stdev1_res = stdev1;
    *t1_res = t1;
    *p_value_res = p_value;

    return 0;
}


static int
linner_regression2(const double *x, const double *y, const uint32_t array_len,
    double *beta0, double *beta0_se, double *beta0_p,
    double *beta1, double *beta1_se, double *beta1_p,
    double *ftest_p_val)
{
    /*
        y = beta0 + beta1 * x
        c0: beta0
        c1: beta1
        cov00: square of standard error of beta0
        cov11: square of standard error of beta1
        sumsq: sum of residule square(RSS)
    */

    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(x, 1, y, 1, array_len, &c0, &c1, &cov00, &cov01, &cov11,
                   &sumsq);

    /*
        t = beta / se
        then carry out t test.
    */
    double stdev0 = sqrt(cov00);
    double t0 = c0 / stdev0;
    double pv0 = (t0 < 0)? 2 * (1 - gsl_cdf_tdist_P(-t0, array_len - 2)):
        2 * (1 - gsl_cdf_tdist_P(t0, array_len - 2));

    double stdev1 = sqrt(cov11);
    double t1 = c1 / stdev1;
    double pv1 = t1 < 0 ? 2 * (1 - gsl_cdf_tdist_P(-t1, array_len - 2)):
        2 * (1 - gsl_cdf_tdist_P(t1, array_len - 2));

    /*
        calculate mean of y and Variance of y.
    */
    double dl = array_len - 2;
    double y_mean = 0;
    for (int i = 0; i < array_len; i++) {
        y_mean += y[i];
    }
    y_mean = y_mean / array_len;

    double y_var = 0;
    for (int i = 0; i < array_len; i++) {
        y_var += pow(y[i] - y_mean, 2);
    }

    /*
        Get f value and carry out f test.
    */
    double R2 = 1 - sumsq / y_var;
    double F = R2 * dl / (1 - R2);
    double ftest_p_value = 1 - gsl_cdf_fdist_P(F, 1, dl);

    *beta0 = c0;
    *beta0_se = stdev0;
    *beta0_p = pv0;
    *beta1 = c1;
    *beta1_se = stdev1;
    *beta1_p = pv1;
    *ftest_p_val = ftest_p_value;

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
/*
    uint32_t *fam_index_arrary_tmp = fam_index_array;
    uint32_t not_duplicated_fam_id_num = 0;
    for (i = 0; i < fam_line_num; i++) {
        this_ptr = &(hash_array[i]);
        while (this_ptr) {
            if (this_ptr->data_index != -1) {
                fam_index_arrary_tmp[not_duplicated_fam_id_num] =
                    (uint32_t)(this_ptr->data_index);
                not_duplicated_fam_id_num++;
            }
            this_ptr = this_ptr->next;
        }
    }
    qsort(fam_index_arrary_tmp, not_duplicated_fam_id_num, sizeof(uint32_t),
        compare_uint32);
*/
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
make_drm_threads_args(int thread_num, uint32_t indi_num_fam, uint32_t indi_num_oii,
    uint32_t vari_num, uint32_t probe_num, uint32_t align_len, char not_need_align, uint32_t variant_load_num,
    char *variant_data, uint64_t variant_data_len,
    uint32_t *fam_index_array, uint32_t *oii_index_array,
    DRM_THREAD_ARGS_ptr thread_args)
{
    for (int i = 0; i < thread_num; i++) {
        thread_args[i].thread_index = i;
        thread_args[i].thread_num = thread_num;

        thread_args[i].fam_num = indi_num_fam;
        thread_args[i].oii_num = indi_num_oii;
        thread_args[i].vari_num = vari_num;
        thread_args[i].probe_num = probe_num;
        thread_args[i].align_len = align_len;
        thread_args[i].not_need_align = not_need_align;
        thread_args[i].fam_index_array = fam_index_array;
        thread_args[i].oii_index_array = oii_index_array;
        thread_args[i].variant_load_num = variant_load_num;
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
        thread_args[i].result = (double *)malloc(sizeof(double) * variant_load_num * 3);

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
free_drm_threads_args_malloc(DRM_THREAD_ARGS_ptr thread_args, int thread_num)
{
    for (int i = 0; i < thread_num; i++) {
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

static void *drm_thread_worker(void *args) {
    DRM_THREAD_ARGS_ptr args_in = (DRM_THREAD_ARGS_ptr)args;
    uint32_t fam_num = args_in->fam_num;
    uint32_t oii_num = args_in->oii_num;
    uint32_t vari_num = args_in->vari_num;
    uint32_t probe_num = args_in->probe_num;
    uint32_t align_len = args_in->align_len;
    char not_need_align = args_in->not_need_align;

    uint32_t *fam_index_array = args_in->fam_index_array;
    uint32_t *oii_index_array = args_in->oii_index_array;

    uint32_t varin_num_loaded = args_in->variant_load_num;
    char *variant_data_loaded = args_in->variant_data;
    uint64_t variant_data_loaded_len = args_in->variant_data_len_char;

    double *pheno_data = args_in->probe_data;
    double *g0_array = args_in->g0_array;
    double *g1_array = args_in->g1_array;
    double *g2_array = args_in->g2_array;

    double *geno_data_aligned = args_in->geno_array;
    double *pheno_data_aligned = args_in->pheno_array;
    double *result = args_in->result;

    double c1_res = 0, stdev_res = 0, t1_res = 0, p_value_res = 0;
    double geno_0_median = 0.0, geno_1_median = 0.0, geno_2_median = 0.0;
    int g0_num = 0, g1_num = 0, g2_num = 0;

    double *g0_array_tmp = NULL;
    double *g1_array_tmp = NULL;
    double *g2_array_tmp = NULL;
    int align_len_rm_missing = 0;
    uint32_t geno_index = 0;
    uint32_t pheno_index = 0;
    double geno_value;
    double pheno_value;

    char *current_geno_one = NULL;
    for (int i = 0; i < varin_num_loaded; i++) {
        current_geno_one = variant_data_loaded + i * fam_num;

        align_len_rm_missing = 0;
        g0_array_tmp = g0_array;
        g1_array_tmp = g1_array;
        g2_array_tmp = g2_array;
        g0_num = 0;
        g1_num = 0;
        g2_num = 0;

        for (int j = 0; j < align_len; j++) {
            if (not_need_align) {
                geno_value = (double)current_geno_one[j];
                pheno_value = pheno_data[j];
            } else {
                geno_index = fam_index_array[j];
                geno_value = (double)current_geno_one[geno_index];
                pheno_index = oii_index_array[j];
                pheno_value = pheno_data[pheno_index];
            }
            // do I need use other way to compare float number?
            if (geno_value != 4.0 && pheno_value != -9.0) {
                geno_data_aligned[align_len_rm_missing] = geno_value;
                pheno_data_aligned[align_len_rm_missing] = pheno_value;
                align_len_rm_missing++;
            }
        }

        for (int k = 0; k < align_len_rm_missing; k++) {
                if (geno_data_aligned[k] == 0.0) {
                *g0_array_tmp = pheno_data_aligned[k];
                g0_num++;
                g0_array_tmp++;
                } else if (geno_data_aligned[k] == 1.0) {
                *g1_array_tmp = pheno_data_aligned[k];
                g1_num++;
                g1_array_tmp++;
                } else {
                *g2_array_tmp = pheno_data_aligned[k];
                g2_num++;
                g2_array_tmp++;
                }
        }

        qsort(g0_array, g0_num, sizeof(double), compare_double);
        qsort(g1_array, g1_num, sizeof(double), compare_double);
        qsort(g2_array, g2_num, sizeof(double), compare_double);

        geno_0_median =
            (g0_num % 2)
                ? g0_array[g0_num / 2]
                : (g0_array[g0_num / 2 - 1] + g0_array[g0_num / 2]) / 2;

        geno_1_median =
            (g1_num % 2)
                ? g1_array[g1_num / 2]
                : (g1_array[g1_num / 2 - 1] + g1_array[g1_num / 2]) / 2;

        geno_2_median =
            (g2_num % 2)
                ? g2_array[g2_num / 2]
                : (g2_array[g2_num / 2 - 1] + g2_array[g2_num / 2]) / 2;

        // substract conressponding median from pheno value. and use absulute
        // value.
        for (int i = 0; i < align_len_rm_missing; i++) {
                pheno_data_aligned[i] =
                    (geno_data_aligned[i] == 0.0)
                        ? fabs(pheno_data_aligned[i] - geno_0_median)
                        : ((geno_data_aligned[i] == 1.0)
                               ? fabs(pheno_data_aligned[i] - geno_1_median)
                               : fabs(pheno_data_aligned[i] - geno_2_median));
        }

        linner_regression(geno_data_aligned, pheno_data_aligned,
                          align_len_rm_missing, &c1_res, &stdev_res, &t1_res,
                          &p_value_res);
        // c1_res beta1 beta1_se t_value_of_beta1, f_test_p
        // printf("%lf %lf %lf %lf\n", c1_res, stdev_res, t1_res, p_value_res);
        result[i * 3] = c1_res;
        result[i * 3 + 1] = stdev_res;
        result[i * 3 + 2] = p_value_res;
    }

    return NULL;
}


static void
drm_print_res(DRM_THREAD_ARGS_ptr thread_args, int thread_num, int probe_num,
    FILE *outfile)
{
    double *res = NULL;
    for (int i = 0; i < thread_num; i++) {
        res = (thread_args[i]).result;
        for (int j = 0; j < probe_num; j++) {
            fprintf(outfile, "%le\t%le\t%le\t%le\t", res[j*4], res[j*4 + 1],
                res[j*4 + 2], res[j*4 + 3]);
        }
        fprintf(outfile, "\n");
    }
    return;
}



static SVLM_THREAD_ARGS_ptr make_svlm_threads_args(
    int thread_num, uint32_t fam_node_num, uint32_t oii_node_num,
    uint32_t probe_num, uint32_t align_len, double *bod_data_all,
    uint32_t *fam_index_array, uint32_t *oii_index_array,
    SVLM_THREAD_ARGS_ptr thread_args) {
    for (int i = 0; i < thread_num; i++) {
        thread_args[i].thread_index = i;
        thread_args[i].fam_num = fam_node_num;
        thread_args[i].oii_num = oii_node_num;
        thread_args[i].probe_num = probe_num;
        thread_args[i].align_len = align_len;

        thread_args[i].bod_data_all_probe = bod_data_all;
        thread_args[i].fam_index_array = fam_index_array;
        thread_args[i].oii_index_arrya = oii_index_array;

        // need test mallc exit status
        thread_args[i].geno_data_this_variant =
            (double *)malloc(sizeof(double) * fam_node_num);
        thread_args[i].pheno_array =
            (double *)malloc(sizeof(double) * align_len);
        thread_args[i].geno_array =
            (double *)malloc(sizeof(double) * align_len);
        thread_args[i].result =
            (double *)malloc(sizeof(double) * probe_num * 8);
        if (!thread_args[i].geno_data_this_variant ||
            !thread_args[i].pheno_array || !thread_args[i].geno_array ||
            !thread_args[i].result) {
                fprintf(stderr, "malloc for thread args number mem failed.\n");
                return NULL;
        }
    }
    return thread_args;
}


static void
free_svlm_threads_args_malloc(SVLM_THREAD_ARGS_ptr thread_args, int thread_num)
{
    for (int i = 0; i < thread_num; i++) {
        free(thread_args[i].geno_data_this_variant);
        free(thread_args[i].geno_array);
        free(thread_args[i].pheno_array);
        free(thread_args[i].result);
    }

    return;
}



static void *
svlm_thread_worker(void *args)
{
    SVLM_THREAD_ARGS_ptr args_in = (SVLM_THREAD_ARGS_ptr)args;
//    int thread_index = args_in->thread_index;
//    uint32_t fam_num = args_in->fam_num;
    uint32_t oii_num = args_in->oii_num;
    uint32_t probe_num = args_in->probe_num;
    uint32_t align_len = args_in->align_len;

    double * bod_data_all = args_in->bod_data_all_probe;
    uint32_t *fam_index_array = args_in->fam_index_array;
    uint32_t *oii_index_array = args_in->oii_index_arrya;

    double *geno_data_one_variant = args_in->geno_data_this_variant;
    double *pheno_data_one_probe = NULL;

    double *geno_array = args_in->geno_array;
    double *pheno_array = args_in->pheno_array;
    double *result = args_in->result;

    uint32_t geno_index = 0, pheno_index = 0;
    uint32_t align_len_rm_missing = 0;
    double geno_value = 0.0, pheno_value = 0.0;
    double beta0, beta0_se, beta0_p, beta1, beta1_se, beta1_p, ftest_p_val;
    uint32_t result_index = 0;
    for (int i = 0; i < probe_num; i++) {
        pheno_data_one_probe = bod_data_all + (i * oii_num);
        align_len_rm_missing = 0;
        for (int j = 0; j < align_len; j++) {
            geno_index = fam_index_array[j];
            pheno_index = oii_index_array[j];
            geno_value = geno_data_one_variant[geno_index];
            pheno_value = pheno_data_one_probe[pheno_index];

            if (geno_value != 4.0 && pheno_value != -9.0) {
                geno_array[align_len_rm_missing] = geno_value;
                pheno_array[align_len_rm_missing] = pheno_value;
                align_len_rm_missing++;
            }            
        }
        linner_regression2(geno_array, pheno_array, align_len_rm_missing, &beta0,
            &beta0_se, &beta0_p, &beta1, &beta1_se, &beta1_p, &ftest_p_val);
        result[result_index++] = beta0;
        result[result_index++] = beta0_se;
        result[result_index++] = beta1;
        result[result_index++] = beta1_se;
        // pheno_array turn into array of residual square.
        for (int k = 0; k < align_len_rm_missing; k++) {
            double residule;
            residule = pheno_array[k] - (beta0 + beta1 * geno_array[k]);
            pheno_array[k] = residule * residule;
        }
        linner_regression2(geno_array, pheno_array, align_len_rm_missing,
                           &beta0, &beta0_se, &beta0_p, &beta1, &beta1_se,
                           &beta1_p, &ftest_p_val);
        result[result_index++] = beta1;
        result[result_index++] = beta1_se;
        result[result_index++] = beta1_p;
        result[result_index++] = ftest_p_val;
        printf("%le\n", ftest_p_val);
    }

    return NULL;
}


static void
svlm_print_res(SVLM_THREAD_ARGS_ptr thread_args, int thread_num, int probe_num,
    FILE *outfile)
{
    /*
        res of every point: beta0, beta0_se, beta1, beta1_se, beta1_secode,
        beta1_se_seconde, beta1_p_seconde, ftest_p_val_seconde
    */

    double *res = NULL;
    for (int i = 0; i < thread_num; i++) {
        res = (thread_args[i]).result;
        for (int j = 0; j < probe_num; j++) {
            double chi2 = res[j * 4 + 4] * res[j * 4 + 5];
            chi2 = chi2 * chi2;
            fprintf(outfile, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t",
                res[j * 8], res[j * 8 + 1], res[j * 8 + 2], res[j * 8 + 3],
                res[j * 8 + 4], res[j * 8 + 5], res[j * 8 + 6], res[j * 8 + 7],
                chi2
            );
        }
        fprintf(outfile, "\n");
    }
    return;
}

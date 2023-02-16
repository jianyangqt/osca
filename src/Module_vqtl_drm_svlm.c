/*
    Utility: OSCA;
    Module: VQTL;
    Methods: drm, svlm;

    See Module_vqtl_drm_svlm.h for more information.

    Note:
    1. All related code was writed in this file for convience. May split into
    files in future for good coding prectice.

    2. The function Module_vqtl_drm and Module_vqtl_svlm is introspective. This
    mean you just offer argc and argv to the function. The function would decide
    whether to run or not himself. This is good to seperate the code from other
    modules and easy to maintain.


    --Benjamin Fang, 20230204
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

/*
#define MODULE_LEVER_DB_INFO
#define METHOD_LEVER_DB_INFO
#define LEVER_3_DB_INFO
#define LEVER_4_DB_INFO
*/

#define MODULE_NAME "vqtl"
#define VQTL_DRM_METHOD "drm"
#define VQTL_SVLM_METHOD "svlm"


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


/*->>plink file*/
#define CharM1 0x03
#define CharM2 0x0c
#define CharM3 0x30

#define CONVERT_GENO(x) ((x) == 0) ? 2 : (((x) == 2) ? 1 : (((x) == 3) ? 0 : 4))

#define BIM_DATA_TYPE 1
typedef struct bim_node {
    char chrom;
    char rsid[64];
    float phy_pos;
    uint32_t pos;
    char *allel1;
    char *allel2;
    struct bim_node *next;

} BIM_NODE, *BIM_NODE_ptr;

typedef struct {
    uint32_t line_num;
    uint32_t current_node_index;
    BIM_NODE_ptr next_iter;
    BIM_NODE_ptr bim_data;
    BIM_NODE_ptr *node_ptr_array;
} BIM_DATA, *BIM_DATA_ptr;

#define FAM_DATA_TYPE 2
typedef struct fam_node {
    char family_id[64];
    char within_famid[64];
    char father_id[64];
    char mother_id[64];
    char sex;
    char phenotype_value;
    struct fam_node *next;
} FAM_NODE, *FAM_NODE_ptr;

typedef struct {
    uint32_t line_num;
    FAM_NODE_ptr next_iter;
    uint32_t current_node_index;
    FAM_NODE_ptr fam_data;
    FAM_NODE_ptr *node_ptr_array;
} FAM_DATA, *FAM_DATA_ptr;

typedef struct {
    char magic_char[3];
    uint32_t indiv_num;
    uint32_t variant_num;
    FILE *filehandle;
    uint32_t current_variant;
    uint64_t current_file_pos; 
    uint32_t data_len;
    char *data;
    uint32_t read_buf_len;
    uint8_t *read_buf;
} BED_DATA, *BED_DATA_ptr;
/*<<-plink file*/


/*->>bod file*/
#define OII_DATA_TYPE 3
typedef struct oii_node {
    char family_id[64];
    char indiv_id[64];
    char parental_id[64];
    char maternal_id[64];
    char sex;
    struct oii_node *next;
} OII_NODE, *OII_NODE_ptr;

typedef struct {
    uint32_t line_num;
    OII_NODE_ptr next_iter;
    uint32_t current_node_index;
    OII_NODE_ptr oii_data;
    OII_NODE_ptr *node_ptr_array;
} OII_DATA, *OII_DATA_ptr;

#define OPI_DATA_TYPE 4
typedef struct opi_node {
    unsigned char chrom;  // 0 for NA, 200 for X, 201 for Y, 202 for MT
    char probe_id[64];
    uint64_t position;  // 0 for NA
    char gene_id[64];
    char ori;  // 0 for NA, 1 for +, 2 for -
    struct opi_node *next;
} OPI_NODE, *OPI_NODE_ptr;

typedef struct {
    uint32_t line_num;
    OPI_NODE_ptr next_iter;
    uint32_t current_node_index;
    OPI_NODE_ptr opi_data;
    OPI_NODE_ptr *node_ptr_array;
} OPI_DATA, *OPI_DATA_ptr;

typedef struct {
    char value_type;
    char data_type;
    uint32_t indiv_num;
    uint32_t probe_num;
    FILE *filehandle;
    uint32_t current_probe;
    uint64_t current_file_pos;
    uint64_t data_len;
    double *data;
} BOD_DATA, *BOD_DATA_ptr;

/*<<-bod file*/


/*->> DRM Thread argument structure*/
typedef struct {
    int thread_index;

    uint32_t fam_num;
    uint32_t oii_num;
    uint32_t probe_num;
    uint32_t align_len;
    
    //pointer to memory have alloced.
    double *bod_data_all_probe;
    uint32_t *fam_index_array;
    uint32_t *oii_index_arrya;
    
    //pointer of mem need to allocate.
    double *geno_data_this_variant;
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


static BIM_DATA_ptr read_bim_file(const char *bim_filename);
static int clean_bim_data(BIM_DATA_ptr data_in);
static FAM_DATA_ptr read_fam_file(const char *fam_filename);
static int clean_fam_data(FAM_DATA_ptr data_in);
static BED_DATA_ptr init_bed_data_struct(const char * bed_filename, uint32_t indiv_num,
    uint32_t variant_num);
static char *read_bed_by_variant(BED_DATA_ptr data_in, char *ext_data_array);
static int finalize_bed_data_struct(BED_DATA_ptr data_in);


static OII_DATA_ptr read_oii_file(const char * oii_filename);
static int clean_oii_data(OII_DATA_ptr data_in);
static OPI_DATA_ptr read_opi_file(const char * opi_filename);
static int clean_opi_data(OPI_DATA_ptr data_in);
static BOD_DATA_ptr init_bod_data_struct(const char * bod_filename);
static double *read_bod_by_probe(BOD_DATA_ptr data_in, double * ext_data_array);
static int finalize_bod_data_struct(BOD_DATA_ptr data_in);
static double *load_bod_to_mem(BOD_DATA_ptr data_in);
static int clean_load_bod(double *bod_data);


static void *next(void *data_in, char data_type);
static void *reset_iter(void *data_in, char data_type);
static void *make_node_ptr_array(void *data_in, char data_type);


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


static int align_fam_oii_ids(FAM_NODE_ptr *fam_node_array, uint32_t fam_node_num,
    OII_NODE_ptr *oii_node_array, uint32_t oii_node_num, uint32_t **fam_index_arrary,
    uint32_t **oii_index_array, uint32_t *aligned_len);
static int clean_fam_oii_align_array(uint32_t *oii_index_array,
    uint32_t *fam_index_array);

int Module_vqtl_drm(int argc, char *argv[]);
static DRM_THREAD_ARGS_ptr make_drm_threads_args(int thread_num, uint32_t fam_node_num,
    uint32_t oii_node_num, uint32_t probe_num, uint32_t align_len,
    double *bod_data_all, uint32_t *fam_index_array, uint32_t *oii_index_array,
    DRM_THREAD_ARGS_ptr thread_args);
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
    enum PROGRESS drm_progress = BEGIN;
    VQTL_ARGS_ptr args = (VQTL_ARGS_ptr)malloc(sizeof(VQTL_ARGS));
    args = vqtl_parse_args_legacy(argc, argv, VQTL_DRM_METHOD, args);
    if (!(args->flag_vqtl) || !(args->method)) {
        return 0;
    } else if (args->flag_help) {
        help_legacy();
        return 1;
    }
    drm_progress = ARG_PARSED;
    return 1;

    printf("\n\033[0;32m>>>VQTL Module DRM method\033[0m Begin\n");
    
    int geno_name_len = strlen(args->arg_geno_file);
    char *geno_file = (char *)malloc(sizeof(char) * (geno_name_len + 5));
    if (!geno_file) {
        fprintf(stderr, "malloc for geno file name failed.\n");
        goto clean_before_exit;
    }
    drm_progress = GENO_FILE_MALLOC;
    
    strcpy(geno_file, args->arg_geno_file);
    strcat(geno_file, ".bim");
    BIM_DATA_ptr bim_data_res = read_bim_file(geno_file);
    if (!bim_data_res) {
        fprintf(stderr, "read bim failed.\n");
        goto clean_before_exit;
    }
    drm_progress = BIM_READED;

    geno_file[geno_name_len] = '\0';
    strcat(geno_file, ".fam");
    FAM_DATA_ptr fam_data_res = read_fam_file(geno_file);
    if (!fam_data_res) {
        fprintf(stderr, "read fam file failed.\n");
        goto clean_before_exit;
    }
    drm_progress = FAM_READED;

    int pheno_name_len = strlen(args->opt_pheno_bod);
    char *pheno_file = (char *)malloc(sizeof(char) * (pheno_name_len + 5));
    if (!pheno_file) {
        fprintf(stderr, "malloc for pheno file failed.\n");
        goto clean_before_exit;
    }
    drm_progress = PHENO_FILE_MALLOC;

    strcpy(pheno_file, args->opt_pheno_bod);
    strcat(pheno_file, ".opi");
    OPI_DATA_ptr opi_data_res = read_opi_file(pheno_file);
    if (!opi_data_res) {
        fprintf(stderr, "read pheno file failed.\n");
        goto clean_before_exit;
    }
    drm_progress = OPI_READED;

    pheno_file[pheno_name_len] = '\0';
    strcat(pheno_file, ".oii");
    OII_DATA_ptr oii_data_res = read_oii_file(pheno_file);
    if (!oii_data_res) {
        fprintf(stderr, "read oii file failed.\n");
        goto clean_before_exit;
    }
    drm_progress = OII_READED;

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
    align_fam_oii_ids(fam_node_array, fam_node_num, oii_node_array, oii_node_num,
        &fam_index_array, &oii_index_array, &align_len);
    if (!fam_index_array || !oii_index_array) {
        fprintf(stderr, "make fam index array or oii index array failed.\n");
        goto clean_before_exit;
    }
    drm_progress = FAM_OII_INDEX_ARRAY_MALLOC;

#if defined METHOD_LEVER_DB_INFO
    for (int i = 0; i < align_len; i++) {
        printf("%d %u %u %s %s\n", i, fam_index_array[i], oii_index_array[i],
            fam_node_array[fam_index_array[i]]->within_famid,
            oii_node_array[oii_index_array[i]]->indiv_id);
    }
#endif

    //load bod data to memory
    pheno_file[pheno_name_len] = '\0';
    strcat(pheno_file, ".bod");
    BOD_DATA_ptr bod_data_res = init_bod_data_struct(pheno_file);
    if (!bod_data_res) {
        fprintf(stderr, "init bod data failed.\n");
        goto clean_before_exit;
    }
    drm_progress = BOD_DATA_INITED;

    double *bod_data_all = load_bod_to_mem(bod_data_res);
    if (!bod_data_all) {
        fprintf(stderr, "load all bod data into mem failed.\n");
        goto clean_before_exit;
    }
    drm_progress = ALL_BOD_DATA_LOADED;

#if defined METHOD_LEVER_DB_INFO
    double *first_probe_dt = read_bod_by_probe(bod_data_res, NULL);
    for (int i = 0; i < oii_node_num; i++) {
        printf("%lf %lf\n", bod_data_all[0 * oii_node_num + i], first_probe_dt[i]);
    }
#endif

    geno_file[geno_name_len] = '\0';
    strcat(geno_file, ".bed");
    BED_DATA_ptr bed_data_res = init_bed_data_struct(
        geno_file, fam_data_res->line_num, bim_data_res->line_num);
    if (!bed_data_res) {
        fprintf(stderr, "init bed data failed.\n");
        goto clean_before_exit;
    }
    drm_progress = BED_DATA_INITED;

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
        fprintf(stderr, "start variant index should less equal than end variant index.");
        goto clean_before_exit;
    }
    drm_progress = GET_VAR_REANGE;

    pthread_t * restrict thread_ids = (pthread_t *)malloc(sizeof(pthread_t) * thread_num);
    if (!thread_ids) {
        fprintf(stderr, "malloc for threads id failed.\n");
        goto clean_before_exit;
    }
    drm_progress = THREAD_ID_MALLOC;

    uint32_t limit = end_variant - thread_num + 1;
    DRM_THREAD_ARGS_ptr thread_args =
        (DRM_THREAD_ARGS_ptr)malloc(sizeof(DRM_THREAD_ARGS) * thread_num);
    if (!thread_args) {
        fprintf(stderr, "malloc for theads args failed.\n");
        goto clean_before_exit;
    }
    drm_progress = THREAD_ARGS_MALLOC;

    DRM_THREAD_ARGS_ptr status = make_drm_threads_args(thread_num, fam_node_num, oii_node_num, probe_num_opi, align_len,
        bod_data_all, fam_index_array, oii_index_array, thread_args);
    if (!status) {
        fprintf(stderr, "malloc for thread args inner memory failed.\n");
        goto clean_before_exit;
    }
    drm_progress = THREAD_ARGS_MEMBER_MALLOC;
    
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
    drm_progress = OUT_FIEL_OPENED;

    uint32_t i = 0;
    char *geno_buf_by_varian = NULL;
    for (i = start_variant; i < limit; i += thread_num) {
        for (uint32_t j = 0; j < thread_num; j++) {
            geno_buf_by_varian = read_bed_by_variant(bed_data_res, NULL);
            for (uint32_t k = 0; k < fam_node_num; k++) {
                ((thread_args[j]).geno_data_this_variant)[k] = 
                    (double)geno_buf_by_varian[k];
            }
        }

        for (uint32_t l = 0; l < thread_num; l++) {
            pthread_create(&(thread_ids[l]), NULL, drm_thread_worker, &(thread_args[l]));
        }

        for (uint32_t m = 0; m < thread_num; m++) {
            pthread_join(thread_ids[m], NULL);
        }
        drm_print_res(thread_args, thread_num, probe_num_opi, outfile);
    }

    for (; i < end_variant; i++) {
        geno_buf_by_varian = read_bed_by_variant(bed_data_res, NULL);
        for (int j = 0; j < fam_node_num; j++) {
            (thread_args[0]).geno_data_this_variant[j] =
                (double)geno_buf_by_varian[j];
        }
        drm_thread_worker(thread_args);
        drm_print_res(thread_args, 1, probe_num_opi, outfile);
    }

    drm_progress = END;

clean_before_exit:
    switch (drm_progress) {
        case END:;
        case OUT_FIEL_OPENED:
            fclose(outfile);
        case THREAD_ARGS_MEMBER_MALLOC:
            free_drm_threads_args_malloc(thread_args, thread_num);
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
            if (drm_progress == END) {
                printf("Run successfully.\n");
            }
            printf("\033[0;32m<<<VQTL Module DRM method\033[0m Returned\n\n");
    }
    return 1;
}


int
Module_vqtl_svlm(int argc, char *argv[])
{
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
    return 1;
}


static void
help_legacy(void)
{
    printf(
        "Help:\n"
        "--help,    flag, print this help message.\n"
        "--vqtl,    flag, use vqtl module.\n"
        "--method    STR, vqtl method. 'drm' for method DRM, 'svlm' for SVLM.\n"
        "--geno      STR, plink genotype files.\n"
        "--pheno     STR, phenotyp file in plain text.(no support yet)\n"
        "--pheno-bod STR, phenotyp file in bod file format.\n"
        "--threads   INT, number of threads to use, default is 1.\n"
        "--start-var INT, index of first variant to calculate. default is 1.(not support yet)\n"
        "--end-var   INT, last variant to calculate. default is last one of bim "
        "file.(not support yet)\n"
        "--start-probe INT, index of first probe.(not support yet)\n"
        "--end-probe INT, last probe, defaulte is last one.(not support yet)\n"
        "--tast-num  INT, tast number.\n"
        "--tast-id   INT, task id.\n"
        "--trans, flag, only calculate trans region.(not support yet)\n"
        "--trans-distance-bp\n"
        "            INT, distance from probe in basepair to define as trans.(not support yet)\n"
        "--cis       INT, only calculate cis region.(not support yet)\n"
        "--cis-window-pb\n"
        "            INT, widow width in basepare defined as cis.(not support yet)\n"
        "--pthresh FLOAT, p value to filter beta1 t test result.\n"
        "--mem       INT, memoray used by program, default is 3/4 all memoray.\n"
        "--out       STR, output file name.\n"
        "--outformat STR, output format, default is besd.(not support yet)\n\n"
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

            if (strcmp(argv[i], "--threads") == 0) {
                if (i + 1 < argc && strncmp(argv[i + 1], "--", 2) != 0) {
                    args->opt_thread = atoi(argv[++i]);
                    printf("--threads %s\n", argv[i]);
                    continue;
                } else {
                    fprintf(stderr, "--threads need a argument.\n");
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


static BIM_DATA_ptr
read_bim_file(const char *bim_filename)
{
#if defined DEBUG_INFO_abc || DEBUG_INFO
    printf("bim file: %s\n", bim_filename);
#endif
    FILE * fin = fopen(bim_filename, "r");
    if (!fin) {
        fprintf(stderr, "open bim file failed.\n");
        return NULL;
    }

    BIM_DATA_ptr data_out = (BIM_DATA_ptr)malloc(sizeof(BIM_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc bim data failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->current_node_index = 0;
    data_out->next_iter = NULL;
    data_out->bim_data = NULL;
    data_out->node_ptr_array = NULL;

    char line_buf[2048];
    line_buf[2047] = '\0';
    BIM_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;
    char chrom[16];
    chrom[15] = '\0';
    char rsid[64];
    rsid[63] = '\0';
    char phy_pos[16];
    phy_pos[15] = '\0';
    char pos[16];
    pos[15] = '\0';
    char allel1[1024];
    allel1[1023] = '\0';
    char allel2[1024];
    allel2[1023] = '\0';
    int str_len = 0;
    char *char_ptr = NULL;
    char *allel1_ptr;
    char *allel2_ptr;
    while (fgets(line_buf, 2048, fin)) {
        if (line_buf[2047] != '\0') {
            fprintf(stderr, "bim line overflow.\n");
            clean_bim_data(data_out);
            return NULL;
        }

        if (sscanf(line_buf, "%s %s %s %s %s %s", chrom, rsid, phy_pos, pos, allel1, allel2) != 6) {
            fprintf(stderr, "split bim field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        new_node = (BIM_NODE_ptr) malloc(sizeof(BIM_NODE));
        if (!new_node) {
            fprintf(stderr, "malloc for bim node failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        new_node->allel1 = NULL;
        new_node->allel2 = NULL;
        if (chrom[15] == '\0') {
            str_len = strlen(chrom);
            char_ptr = chrom;
            while (isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->chrom = (unsigned char)atoi(chrom);
            } else if (strcmp(chrom, "X") == 0) {
                new_node->chrom = 200;
            } else if (strcmp(chrom, "Y") == 0) {
                new_node->chrom = 201;
            } else if (strcmp(chrom, "MT") == 0) {
                new_node->chrom = 202;
            } else if (strcmp(chrom, "NA") == 0) {
                new_node->chrom = 0;
            } else {
                fprintf(stderr, "chrom not reganized.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        }

        if (rsid[63] == '\0') {
            strcpy(new_node->rsid, rsid);
        } else {
            fprintf(stderr, "rsid field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        if (phy_pos[15] == '\0') {
            new_node->phy_pos = (float)atof(phy_pos);
        } else {
            fprintf(stderr, "physical position field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        if (pos[15] == '\0') {
            str_len = strlen(pos);
            char_ptr = pos;
            while(isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->pos = (uint32_t)atol(pos);
            } else if (strcmp(pos, "NA") == 0) {
                new_node->pos = 0;
            } else {
                fprintf(stderr, "bim pos not recognized.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "bim position field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        if (allel1[1023] == '\0') {
            str_len = strlen(allel1);
            allel1_ptr = (char *)malloc(++str_len);
            if (allel1_ptr) {
                strcpy(allel1_ptr, allel1);
                new_node->allel1 = allel1_ptr;
            } else {
                fprintf(stderr, "malloc for allel1 failed.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        }
        if (allel2[1023] == '\0') {
            str_len = strlen(allel2);
            allel2_ptr = (char *)malloc(++str_len);
            if (allel2_ptr) {
                strcpy(allel2_ptr, allel2);
                new_node->allel2 = allel2_ptr;
            } else {
                fprintf(stderr, "malloc for allel2 failed.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        }

        if (!first_node) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->next_iter = new_node;
            data_out->bim_data = new_node;
            data_out->line_num++;
        } else {
            node_ptr->next = new_node;
            node_ptr = new_node;
            data_out->line_num++;
        }
    }
    fclose(fin);
    return data_out;
}


static int
clean_bim_data(BIM_DATA_ptr data_in)
{
    if (!data_in) {
        return 0;
    }
    if (data_in->node_ptr_array) {
        free(data_in->node_ptr_array);
        data_in->node_ptr_array = NULL;
    }
    if (data_in->next_iter) {
        data_in->next_iter = NULL;
    }
    BIM_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->bim_data;
    data_in->bim_data = NULL;

    while(last_node) {
        next_node = last_node->next;
        free(last_node);
        last_node = next_node;
    }
    
    free(data_in);
    return 0;
}


static FAM_DATA_ptr
read_fam_file(const char *fam_filename)
{
#if defined DEBUG_INFO_abc || DEBUG_INFO
    printf("fam file: %s\n", fam_filename);
#endif

    FILE *fin = fopen(fam_filename, "r");
    if (!fin) {
        fprintf(stderr, "open fam file failed.\n");
        return NULL;
    }

    FAM_DATA_ptr data_out = (FAM_DATA_ptr)malloc(sizeof(FAM_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc for fam data failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->next_iter = NULL;
    data_out->current_node_index = 0;
    data_out->fam_data = NULL;
    data_out->node_ptr_array = NULL;

    char line_buf[512];
    line_buf[511] = '\0';
    char fam_id[64];
    fam_id[63] = '\0';
    char within_famid[64];
    within_famid[63] = '\0';
    char father_id[64];
    father_id[63] = '\0';
    char mother_id[64];
    mother_id[63] = '\0';
    char sex[16];
    sex[15] = '\0';
    char pheno_val[16];
    pheno_val[15] = '\0';
    FAM_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;

    while (fgets(line_buf, 512, fin)) {
        if (line_buf[511] != '\0') {
            fprintf(stderr, "fam line buffer overflow.\n");
            clean_fam_data(data_out);
            return NULL;
        }

        if (sscanf(line_buf, "%s %s %s %s %s %s", fam_id, within_famid, 
            father_id, mother_id, sex, pheno_val) != 6) {
            fprintf(stderr, "split fam line field failed.\n");
            clean_fam_data(data_out);
            return NULL;   
        }
        new_node = (FAM_NODE_ptr)malloc(sizeof(FAM_NODE));
        if (fam_id[63] == '\0') {
            strcpy(new_node->family_id, fam_id);
        } else {
            fprintf(stderr, "family id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (within_famid[63] == '\0') {
            strcpy(new_node->within_famid, within_famid);
        } else {
            fprintf(stderr, "within family id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (father_id[63] == '\0') {
            strcpy(new_node->father_id, father_id);
        } else {
            fprintf(stderr, "father id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (mother_id[63] == '\0') {
            strcpy(new_node->mother_id, mother_id);
        } else {
            fprintf(stderr, "mother id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (sex[15] == '\0') {
            if ((strlen(sex) == 1) && (isdigit(sex[0]))) {
                new_node->sex = (char)atoi(sex);
            } else if (strcmp(sex, "NA") == 0) {
                new_node->sex = 0;
            } else {
                fprintf(stderr, "sex field is not recognized.\n");
                clean_fam_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "sex field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }

        if (pheno_val[15] == '\0') {
            if ((strlen(pheno_val) == 1) && (isdigit(pheno_val[0]))) {
                new_node->phenotype_value = (char)atoi(pheno_val);
            } else if (strcmp(pheno_val, "-9") == 0) {
                new_node->phenotype_value = 0;
            } else if (strcmp(pheno_val, "NA") == 0) {
                new_node->phenotype_value = 0;
            } else {
                fprintf(stderr, "phenotype value field failed.\n");
                clean_fam_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "phenotype value field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }

        if (!first_node) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->next_iter = new_node;
            data_out->fam_data = new_node;
            data_out->line_num++;
        } else {
            node_ptr->next = new_node;
            node_ptr = new_node;
            data_out->line_num++;
        }
    }
    fclose(fin);
    return data_out;
}


static int
clean_fam_data(FAM_DATA_ptr data_in)
{
    if (!data_in) {
        return 0;
    }
    if (data_in->node_ptr_array) {
        free(data_in->node_ptr_array);
        data_in->node_ptr_array = NULL;
    }
    if (data_in->next_iter) {
        data_in->next_iter = NULL;
    }
    FAM_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->fam_data;
    data_in->fam_data = NULL;
    while (last_node) {
        next_node = last_node->next;
        free(last_node);
        last_node = next_node;
    }
    free(data_in);
    return 0;
}


static BED_DATA_ptr
init_bed_data_struct(const char * bed_filename, uint32_t indiv_num, uint32_t variant_num)
{
    FILE *fin = fopen(bed_filename, "r");
    if (!fin) {
        fprintf(stderr, "open bed file failed.\n");
        return NULL;
    }

    BED_DATA_ptr data_out = (BED_DATA_ptr)malloc(sizeof(BED_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc for BED_DATA failed.\n");
        return NULL;
    }
    data_out->filehandle = NULL;
    data_out->data_len = 0;
    data_out->data = NULL;
    data_out->current_file_pos = 0;
    data_out->current_variant = 0;
    data_out->indiv_num = 0;
    data_out->variant_num = 0;
    data_out->read_buf_len = 0;
    data_out->read_buf = NULL;

    char buf[3];
    if (fread(buf, sizeof(char), 3, fin) != 3) {
        fprintf(stderr, "read bed first 3 char failed.\n");
        fclose(fin);
        return NULL;
    }
    memcpy(data_out->magic_char, buf, 3);
    data_out->indiv_num = indiv_num;
    data_out->variant_num = variant_num;
    data_out->filehandle = fin;
    data_out->current_file_pos = 3;
    data_out->current_variant = 0;
    data_out->data = (char *)malloc(ceil((double)indiv_num / 4) * 4);
    data_out->read_buf_len = ceil((double)indiv_num / 4);
    data_out->read_buf = (uint8_t *)malloc(ceil((double)indiv_num / 4));
    return data_out;
}


static char *
read_bed_by_variant(BED_DATA_ptr data_in, char *ext_data_array)
{
    if (data_in->current_variant >= data_in->variant_num) {
        return NULL;
    }

    data_in->data_len = 0;
    uint32_t read_len = data_in->read_buf_len;
    if (fread(data_in->read_buf, sizeof(uint8_t), data_in->read_buf_len, data_in->filehandle) != data_in->read_buf_len) {
        fprintf(stderr, "read bed file by variant failed.\n");
        return NULL;
    }
    data_in->current_file_pos += data_in->read_buf_len;
    data_in->current_variant += 1;

    char *store_array = NULL;
    if (ext_data_array) {
        store_array = ext_data_array;
    } else {
        store_array = data_in->data;
        data_in->data_len = data_in->indiv_num;
    }
    int i = 0, j = 0;
    uint8_t char_tmp;
    char Low_1_2b = 0,  Low_2_2b = 0, Low_3_2b = 0, Low_4_2b = 0;
    uint8_t *read_buf = data_in->read_buf;
    for (i = 0; i < read_len; i++) {
        char_tmp = read_buf[i];
        Low_1_2b = char_tmp & CharM1;
        Low_2_2b = (char_tmp & CharM2) >> 2;
        Low_3_2b = (char_tmp & CharM3) >> 4;
        Low_4_2b = char_tmp >> 6;
        store_array[j] = CONVERT_GENO(Low_1_2b);
        store_array[++j] = CONVERT_GENO(Low_2_2b);
        store_array[++j] = CONVERT_GENO(Low_3_2b);
        store_array[++j] = CONVERT_GENO(Low_4_2b);

        j++;
    }
    if (j != data_in->indiv_num) {
        fprintf(stderr, "decode genotype failed.\n");
        return NULL;
    }

    return store_array;
}


static int
finalize_bed_data_struct(BED_DATA_ptr data_in)
{
    if (data_in) {
        if (data_in->filehandle) {
            fclose(data_in->filehandle);
        }
        if (data_in->data) {
            free(data_in->data);
        }
        if (data_in->read_buf) {
            free(data_in->read_buf);
        }
        free(data_in);
    }
    return 0;
}


static OII_DATA_ptr
read_oii_file(const char * oii_filename)
{
    FILE *fin = fopen(oii_filename, "r");
    if (!fin) {
        fprintf(stderr, "open oii file failed.\n");
        return NULL;
    }

    OII_DATA_ptr data_out = (OII_DATA_ptr) malloc(sizeof(OII_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->next_iter = NULL;
    data_out->current_node_index = 0;
    data_out->oii_data = NULL;
    data_out->node_ptr_array = NULL;

    
    char family_id[64];
    family_id[63] = '\0';
    char indiv_id[64];
    indiv_id[63] = '\0';
    char parental_id[64];
    parental_id[63] = '\0';
    char maternal_id[64];
    maternal_id[63] = '\0';
    char sex[64];
    sex[63] = '\0';
    OII_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;
    char line_buf[512];
    line_buf[511] = '\0';
    while(fgets(line_buf, 512, fin)) {
        if (line_buf[511] != '\0') {
            fprintf(stderr, "read oii line buffer overflow.\n");
            return NULL;
        }
        if (strlen(line_buf) > 262) {
            fprintf(stderr, "oii line field overflow.\n");
            return NULL;
        }

        if (sscanf(line_buf, "%s %s %s %s %s", family_id, indiv_id, parental_id,
            maternal_id, sex) != 5) {
            fprintf(stderr, "sscanf failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }

        new_node = (OII_NODE_ptr) malloc(sizeof(OII_NODE));
        if (!new_node) {
            fprintf(stderr, "malloc for oii node failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        new_node->next = NULL;

        if (family_id[63] == '\0') {
            strcpy(new_node->family_id, family_id);
        } else {
            fprintf(stderr, "family id field overflow.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (indiv_id[63] == '\0') {
            strcpy(new_node->indiv_id, indiv_id);
        } else {
            fprintf(stderr, "oii individual field overflow.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (parental_id[63] == '\0') {
            strcpy(new_node->parental_id, parental_id);
        } else {
            fprintf(stderr, "parental id field failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (maternal_id[63] == '\0') {
            strcpy(new_node->maternal_id, maternal_id);
        } else {
            fprintf(stderr, "maternal id field failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (sex[63] == '\0') {
            if (strcmp(sex, "NA") == 0) {
                new_node->sex = 0;
            } else if (strlen(sex) == 1 && isdigit(*sex)) {
                new_node->sex = (char) atoi(sex);
            } else {
                fprintf(stderr, "sex field not recgnized.\n");
                clean_oii_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "sex field failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        
        if (first_node == NULL) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->line_num ++;
            data_out->next_iter = first_node;
            data_out->oii_data = first_node;
        } else {
            data_out->line_num ++;
            node_ptr->next = new_node;
            node_ptr = new_node;
        }
    }
    
    return data_out;
}


static int
clean_oii_data(OII_DATA_ptr data_in) 
{
    if (data_in->node_ptr_array) {
        free(data_in->node_ptr_array);
        data_in->node_ptr_array = NULL;
    }
    if (data_in->next_iter) {
        data_in->next_iter = NULL;
    }
    OII_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->oii_data;
    data_in->oii_data = NULL;
    while(last_node) {
        next_node = last_node -> next;
        free(last_node);
        last_node = next_node;
    }

    free(data_in);
    return 0;
}


static OPI_DATA_ptr
read_opi_file(const char * opi_filename)
{
    OPI_DATA_ptr data_out = (OPI_DATA_ptr) malloc(sizeof(OPI_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc for opi data failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->next_iter = NULL;
    data_out->current_node_index = 0;
    data_out->opi_data = NULL;
    data_out->node_ptr_array = NULL;

    FILE * fin = fopen(opi_filename, "r");
    if (!fin) {
        fprintf(stderr, "open opi file failed.\n");
        return NULL;
    }

    char chrom[16];
    chrom[15] = '\0';
    char probe_id[64];
    probe_id[63] = '\0';
    char position[16];
    position[15] = '\0';
    char gene_id[64];
    gene_id[63] = '\0';
    char oritation[16];
    oritation[15] = '\0';
    OPI_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;
    char line_buf[512];
    line_buf[511] = '\0';
    char * char_ptr = NULL;
    int str_len = 0;
    while(fgets(line_buf, 512, fin)) {
        if (line_buf[511] != '\0') {
            fprintf(stderr, "read oii line buffer overflow.\n");
            return NULL;
        }
        if (sscanf(line_buf, "%s %s %s %s %s", chrom, probe_id, position, 
            gene_id, oritation) != 5) {
            fprintf(stderr, "opi field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }

        new_node = (OPI_NODE_ptr) malloc(sizeof(OPI_NODE));
        if (!new_node) {
            fprintf(stderr, "malloc for new opi node failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        new_node->next = NULL;
        if (chrom[15] == '\0') {
            str_len = strlen(chrom);
            char_ptr = chrom;
            while (isdigit(*char_ptr)) { 
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->chrom = (unsigned char) atoi(chrom);
            } else if (strcmp(chrom, "X") == 0) {
                new_node->chrom = 200;
            } else if (strcmp(chrom, "Y") == 0) {
                new_node->chrom = 201;
            } else if (strcmp(chrom, "MT") == 0) {
                new_node->chrom = 202;
            } else if (strcmp(chrom, "NA") == 0){
                new_node->chrom = 0;
            } else {
                fprintf(stderr, "chrom not reganized.\n");
                clean_opi_data(data_out);
                return NULL;
            }
        }
        if (probe_id[63] == '\0') {
            strcpy(new_node->probe_id, probe_id);
        } else {
            fprintf(stderr, "probe id field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        if (position[15] == '\0') {
            char_ptr = position;
            str_len = strlen(position);
            while(isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->position = (uint64_t) atol(position);
            } else if (strcmp(position, "NA") == 0) {
                new_node->position = 0;
            } else {
                fprintf(stderr, "unknow position field.\n");
                clean_opi_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "opi position failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        if (gene_id[63] == '\0') {
            strcpy(new_node->gene_id, gene_id);
        } else {
            fprintf(stderr, "gene id field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        if (oritation[15] == '\0') {
            if (strcmp(oritation, "+") == 0) {
                new_node->ori = 1;
            } else if (strcmp(oritation, "-") == 0) {
                new_node->ori = 2;
            } else if (strcmp(oritation, "NA") == 0) {
                new_node->ori = 0;
            } else {
                fprintf(stderr, "oritation not recognized.\n");
                clean_opi_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "oritation field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }

        if (first_node == NULL) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->line_num ++;
            data_out->opi_data = first_node;
            data_out->next_iter = first_node;

        } else {
            data_out->line_num ++;
            node_ptr->next = new_node;
            node_ptr = new_node;
        }
    }

    return data_out;

}


static int
clean_opi_data(OPI_DATA_ptr data_in)
{
    if (data_in->node_ptr_array) {
        free(data_in->node_ptr_array);
        data_in->node_ptr_array = NULL;
    }
    data_in->next_iter = NULL;
    OPI_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->opi_data;
    data_in->opi_data = NULL;
    if (last_node) {
        next_node = last_node->next;
        free(last_node);
        last_node = next_node;
    }
    free(data_in);
    return 0;
}


static BOD_DATA_ptr
init_bod_data_struct(const char * bod_filename)
{
    FILE * fin = fopen(bod_filename, "rb");
    if (!fin){
        fprintf(stderr, "open file failed.\n");
        return NULL;
    }

    BOD_DATA_ptr data_out = (BOD_DATA_ptr)malloc(sizeof(BOD_DATA));
    if (!data_out) {
        fprintf(stderr, "alloc mem failed.\n");
        return NULL;
    }
    data_out->filehandle = NULL;
    data_out->data = NULL;
    data_out->data_len = 0;

    char buf[12];
    if (fread(buf, sizeof(char), 12, fin) != 12) {
        fprintf(stderr, "read bod first 12 bytes faild.\n");
        fclose(fin);
        return NULL;
    }

    data_out->value_type = buf[0];
    data_out->data_type = buf[1];
    data_out->indiv_num = ((uint32_t *)buf)[1];
    data_out->probe_num = ((uint32_t *)buf)[2];
    data_out->filehandle = fin;
    data_out->current_probe = 0;
    data_out->current_file_pos = 12;
    data_out->data = (double *)malloc(sizeof(double) * data_out->indiv_num);
    
    return data_out;
}


static double *
read_bod_by_probe(BOD_DATA_ptr data_in, double *ext_data_array)
{
    if (data_in->current_probe >= data_in->probe_num){
        return NULL;
    }

    double *data_out = NULL;
    data_in->data_len = 0;
    if (ext_data_array) {
        data_out = ext_data_array;
        if (fread(ext_data_array, sizeof(double), data_in->indiv_num,
                  data_in->filehandle) != data_in->indiv_num) {
            fprintf(stderr, "read probe data failed.\n");
            return NULL;
        }
    } else {
        data_out = data_in->data;
        if (fread(data_in->data, sizeof(double), data_in->indiv_num,
                  data_in->filehandle) != data_in->indiv_num) {
            fprintf(stderr, "read probe data failed.\n");
            return NULL;
        }
        data_in->data_len = data_in->indiv_num;
    }

    data_in->current_probe++;
    data_in->current_file_pos += sizeof(double) * data_in->indiv_num;
    return data_out;
}


static int
finalize_bod_data_struct(BOD_DATA_ptr data_in)
{
    if (data_in->data_len != 0){
        free(data_in->data);
    }
    if (data_in->filehandle) {
        fclose(data_in->filehandle);
    }
    free(data_in);
    return 0;
}


static double *
load_bod_to_mem(BOD_DATA_ptr data_in)
{
    fseek(data_in->filehandle, 12, SEEK_SET);
    uint32_t read_len = data_in->probe_num * data_in->indiv_num;
    double *bod_data_all = (double *)malloc(sizeof(double) * read_len);
    if (fread(bod_data_all, sizeof(double), read_len, data_in->filehandle) != read_len) {
        fprintf(stderr, "load bod to memory failed.\n");
        return NULL;
    }
    fseek(data_in->filehandle, data_in->current_file_pos, SEEK_SET);
    return bod_data_all;
}


static int
clean_load_bod(double * bod_data)
{
    if (bod_data) {
        free(bod_data);
    }
    return 0;
}


static void *
next(void *data_in, char data_type)
{
    if (!data_in) {
        return NULL;
    }
    void *data_out;
    //bim data_type: 1
    if (data_type == BIM_DATA_TYPE) {
        BIM_DATA_ptr data_bim = (BIM_DATA_ptr) data_in;
        if (data_bim->next_iter) {
            data_bim->current_node_index ++;
            data_out = (void *)data_bim->next_iter;
            data_bim->next_iter = data_bim->next_iter->next;
        } else {
            data_out = NULL;
        }
    } else if (data_type == FAM_DATA_TYPE) {
        FAM_DATA_ptr data_fam = (FAM_DATA_ptr)data_in;
        if (data_fam->next_iter) {
            data_fam->current_node_index ++;
            data_out = (void *)data_fam->next_iter;
            data_fam->next_iter = data_fam->next_iter->next;
        } else {
            data_out = NULL;
        }
    } else if (data_type == OII_DATA_TYPE) {
        OII_DATA_ptr data_oii = (OII_DATA_ptr)data_in;
        if (data_oii->next_iter) {
            data_oii->current_node_index ++;
            data_out = (void *)data_oii->next_iter;
            data_oii->next_iter = data_oii->next_iter->next;
        } else {
            data_out = NULL;
        }
    } else if (data_type == OPI_DATA_TYPE) {
        OPI_DATA_ptr data_opi = (OPI_DATA_ptr)data_in;
        if (data_opi->next_iter) {
            data_opi->current_node_index ++;
            data_out = (void *)data_opi->next_iter;
            data_opi->next_iter = data_opi->next_iter->next;
        } else {
            data_out = NULL;
        }
    } else {
        fprintf(stderr, "FILE TYPE not recognized.\n");
        return NULL;
    }

    return data_out;
}


static void *
reset_iter(void * data_in, char data_type)
{
    if (!data_in) {
        return NULL;
    }
    void *data_out;
    // bim data_type: 1
    if (data_type == BIM_DATA_TYPE) {
        BIM_DATA_ptr data_bim = (BIM_DATA_ptr)data_in;
        data_bim->next_iter = data_bim->bim_data;
        data_bim->current_node_index = 0;
        data_out = (void *)data_bim;
    } else if (data_type == FAM_DATA_TYPE) {
        FAM_DATA_ptr data_fam = (FAM_DATA_ptr)data_in;
        data_fam->next_iter = data_fam->fam_data;
        data_fam->current_node_index = 0;
        data_out = (void *)data_fam;
    } else if (data_type == OII_DATA_TYPE) {
        OII_DATA_ptr data_oii = (OII_DATA_ptr)data_in;
        data_oii->next_iter = data_oii->oii_data;
        data_oii->current_node_index = 0;
        data_out = (void *)data_oii;
    } else if (data_type == OPI_DATA_TYPE) {
        OPI_DATA_ptr data_opi = (OPI_DATA_ptr)data_in;
        data_opi->next_iter = data_opi->opi_data;
        data_opi->current_node_index = 0;
        data_out = (void *)data_opi;
    } else {
        fprintf(stderr, "FILE TYPE not recognized.\n");
        return NULL;
    }
    return data_out;
}


static void *
make_node_ptr_array(void * data_in, char data_type)
{
    if (!data_in) {
        return NULL;
    }
    reset_iter(data_in, data_type);
    void * data_out = NULL;
    uint32_t i = 0;
    if (data_type == BIM_DATA_TYPE) {
        BIM_DATA_ptr data_bim = (BIM_DATA_ptr)data_in;
        data_bim->node_ptr_array = (BIM_NODE_ptr *)malloc(sizeof(BIM_NODE_ptr) * 
            data_bim->line_num);
        BIM_NODE_ptr *node_ptr_array_bim = data_bim->node_ptr_array;
        BIM_NODE_ptr node_iter_bim = NULL;
        while ((node_iter_bim = (BIM_NODE_ptr)next((void *)data_bim, data_type))) {
            node_ptr_array_bim[i] = node_iter_bim;
            i++;
        }
        data_out = (void *)node_ptr_array_bim;
    } else if (data_type == FAM_DATA_TYPE) {
        FAM_DATA_ptr data_fam = (FAM_DATA_ptr)data_in;
        data_fam->node_ptr_array = (FAM_NODE_ptr *)malloc(sizeof(FAM_NODE_ptr) * 
            data_fam->line_num);
        FAM_NODE_ptr *node_ptr_array_fam = data_fam->node_ptr_array;
        FAM_NODE_ptr node_iter_fam = NULL;
        while ((node_iter_fam = (FAM_NODE_ptr)next((void *)data_fam, data_type))) {
            node_ptr_array_fam[i] = node_iter_fam;
            i++;
        }
        data_out = (void *)node_ptr_array_fam;
    } else if (data_type == OII_DATA_TYPE) {
        OII_DATA_ptr data_oii = (OII_DATA_ptr)data_in;
        data_oii->node_ptr_array = (OII_NODE_ptr *)malloc(sizeof(OII_NODE_ptr) *
            data_oii->line_num);
        OII_NODE_ptr *node_ptr_array_oii = data_oii->node_ptr_array;
        OII_NODE_ptr node_iter_oii = NULL;
        while ((node_iter_oii = (OII_NODE_ptr)next((void *)data_oii, data_type))) {
            node_ptr_array_oii[i] = node_iter_oii;
            i++;
        }
        data_out = node_ptr_array_oii;
    } else if (data_type == OPI_DATA_TYPE) {
        OPI_DATA_ptr data_opi = (OPI_DATA_ptr)data_in;
        data_opi->node_ptr_array = (OPI_NODE_ptr *)malloc(sizeof(OPI_NODE_ptr) *
            data_opi->line_num);
        OPI_NODE_ptr *node_ptr_array_opi = data_opi->node_ptr_array;
        OPI_NODE_ptr node_iter_opi = NULL;
        while ((node_iter_opi = (OPI_NODE_ptr)next((void *)data_opi, data_type))) {
            node_ptr_array_opi[i] = node_iter_opi;
            i++;
        }
        data_out = node_ptr_array_opi;
    } else {
        fprintf(stderr, "file type not recognized.\n");
        data_out = NULL;
    }
    reset_iter(data_in, data_type);
    return data_out;
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
align_fam_oii_ids(FAM_NODE_ptr * fam_node_array, uint32_t fam_node_num,
    OII_NODE_ptr * oii_node_array, uint32_t oii_node_num, uint32_t **fam_index_arrary,
    uint32_t **oii_index_array_mapped_to_fam, uint32_t *aligned_len)
{
    struct HASH_NODE {
        int32_t data_index;
        struct HASH_NODE *next;
    };

    struct HASH_NODE *hash_array = (struct HASH_NODE *)malloc(sizeof(struct HASH_NODE) *
        fam_node_num);
    int32_t i = 0, j = 0;
    for (i = 0; i < fam_node_num; i++) {
        hash_array[i].data_index = -1;
        hash_array[i].next = NULL;
    }

    uint32_t hash_value = 0;
    for (i = 0; i < fam_node_num; i++) {
        char *within_famid = NULL;
        within_famid = fam_node_array[i]->within_famid;
        hash_value = BKDRHash(within_famid);
        hash_value %= fam_node_num;
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
#if defined DEBUG_INFO_abc || DEBUG_INFO
printf("fam has table\n");
    for (i = 0; i < fam_node_num; i++) {
        this_ptr = &(hash_array[i]);
        while (this_ptr && (this_ptr->data_index != -1)) {
            printf("|%u -> %d ", i, this_ptr->data_index);
            this_ptr = this_ptr->next;
        }
        printf("\n");
    }
#endif
    char *str_a = NULL, *str_b = NULL;
    uint32_t duplicat_num = 0;
    for (i = 0; i < fam_node_num; i++) {  
        this_ptr = &(hash_array[i]);
        while (this_ptr && (this_ptr->data_index != -1)) {
            str_a = fam_node_array[this_ptr->data_index]->within_famid;
            next_ptr = this_ptr->next;
            while (next_ptr && next_ptr->data_index != -1) {
                str_b = fam_node_array[next_ptr->data_index]->within_famid;
                if (strcmp(str_a, str_b) == 0) {
                    duplicat_num++;
                    next_ptr->data_index = -1;    
                }
                next_ptr = next_ptr->next;
            }
            this_ptr = this_ptr->next;
        }     
    }
    if (duplicat_num > 0) {
        fprintf(stderr, "duplicated fam file within family id found\n");
    }

    *fam_index_arrary = (uint32_t *)malloc(sizeof(uint32_t) * fam_node_num);
    uint32_t *fam_index_arrary_tmp = *fam_index_arrary;
    uint32_t not_duplicated_fam_id_num = 0;
    for (i = 0; i < fam_node_num; i++) {
        this_ptr = &hash_array[i];
        while (this_ptr) {
            if (this_ptr->data_index != -1) {
                fam_index_arrary_tmp[not_duplicated_fam_id_num] = (uint32_t)(this_ptr->data_index);
                not_duplicated_fam_id_num++;
            }
            this_ptr = this_ptr->next;
        }
    }
    //printf("%u %u\n", duplicat_num, not_duplicated_fam_id_num);
    qsort(fam_index_arrary_tmp, not_duplicated_fam_id_num, sizeof(uint32_t), compare_uint32);
#if defined DEBUG_INFO_abc || DEBUG_INFO
    for (i = 0; i < not_duplicated_fam_id_num; i++) {
        printf("%u\n", fam_index_arrary_tmp[i]);
    }
#endif

    // algin oii index
    this_ptr = NULL;
    int32_t oii_mapped_index_to_fam = -1;
    int32_t *oii_mapped_index_array_tmp = (int32_t *)malloc(sizeof(int32_t) * oii_node_num);
    for (i = 0; i < oii_node_num; i++) {
        oii_mapped_index_to_fam = -1;
        char *oii_individ = oii_node_array[i]->indiv_id;
        hash_value = BKDRHash(oii_individ);
        hash_value %= fam_node_num;
        this_ptr = &(hash_array[hash_value]);
        while (this_ptr && this_ptr->data_index != -1) {
            if (strcmp(fam_node_array[this_ptr->data_index]->within_famid,
                oii_node_array[i]->indiv_id) == 0)
            oii_mapped_index_to_fam = this_ptr->data_index;
            this_ptr = this_ptr->next;
        }
        oii_mapped_index_array_tmp[i] = oii_mapped_index_to_fam;
    }

    //remove duplication of oii individuals.
    for (i = 0; i < oii_node_num; i++) {
        if (oii_mapped_index_array_tmp[i] != -1) {
            for (j = i + 1; j < oii_node_num; j++) {
                if (oii_mapped_index_array_tmp[i] == oii_mapped_index_array_tmp[j]) {
                    oii_mapped_index_array_tmp[j] = -1;
                }
            }
        }
    }
    //copy oii index to new array and count oii aligned length.
    uint32_t oii_index_mapped_num = 0;
    *oii_index_array_mapped_to_fam = (uint32_t *)malloc(sizeof(uint32_t) * oii_node_num);
    for (i = 0; i < oii_node_num; i++) {
        if (oii_mapped_index_array_tmp[i] != -1) {
            (*oii_index_array_mapped_to_fam)[oii_index_mapped_num] = oii_mapped_index_array_tmp[i];
            oii_index_mapped_num++;
        }
    }
    free(oii_mapped_index_array_tmp);
    
    if (not_duplicated_fam_id_num != oii_index_mapped_num) {
        printf("%u %u\n", not_duplicated_fam_id_num, oii_index_mapped_num);
    }
    *aligned_len = oii_index_mapped_num;
    return 0;
}


static int
clean_fam_oii_align_array(uint32_t *oii_index_array, uint32_t *fam_index_array)
{
    if (oii_index_array) {
        free(oii_index_array);
    }
    if (fam_index_array) {
        free(fam_index_array);
    }
    return 0;
}


static DRM_THREAD_ARGS_ptr
make_drm_threads_args(int thread_num, uint32_t fam_node_num, uint32_t oii_node_num,
    uint32_t probe_num, uint32_t align_len, double *bod_data_all,
    uint32_t *fam_index_array, uint32_t *oii_index_array,
    DRM_THREAD_ARGS_ptr thread_args)
{
#if defined DEBUG    
    printf("fam_node_num: %u, oii_node_num: %u, probe_num: %u, align_len: %u, "
        "bod_data_all: %p, fam_index_array: %p, oii_index_array: %p\n",
        fam_node_num, oii_node_num, probe_num, align_len, bod_data_all,
        fam_index_array, oii_index_array);
#endif
    for (int i = 0; i < thread_num; i++) {
        thread_args[i].thread_index = i;
        thread_args[i].fam_num = fam_node_num;
        thread_args[i].oii_num = oii_node_num;
        thread_args[i].probe_num = probe_num;
        thread_args[i].align_len = align_len;

        thread_args[i].bod_data_all_probe = bod_data_all;
        thread_args[i].fam_index_array = fam_index_array;
        thread_args[i].oii_index_arrya = oii_index_array;

        //all need test malloc exit status.
        thread_args[i].geno_data_this_variant = (double *)malloc(
            sizeof(double) * fam_node_num);
        thread_args[i].g0_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].g1_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].g2_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].geno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].pheno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].result = (double *)malloc(sizeof(double) * probe_num * 4);
        if (!thread_args[i].geno_data_this_variant || !thread_args[i].g0_array ||
            !thread_args[i].g1_array || !thread_args[i].g2_array ||
            !thread_args[i].geno_array || !thread_args[i].pheno_array ||
            !thread_args[i].result) {
            fprintf(stderr, "malloc for thread args number mem failed.\n");
            return NULL;
        }
    }
    
    return thread_args;
}


static SVLM_THREAD_ARGS_ptr
make_svlm_threads_args(int thread_num, uint32_t fam_node_num, uint32_t oii_node_num,
    uint32_t probe_num, uint32_t align_len, double *bod_data_all,
    uint32_t *fam_index_array, uint32_t *oii_index_array,
    SVLM_THREAD_ARGS_ptr thread_args)
{
    for (int i = 0; i < thread_num; i++) {
        thread_args[i].thread_index = i;
        thread_args[i].fam_num = fam_node_num;
        thread_args[i].oii_num = oii_node_num;
        thread_args[i].probe_num = probe_num;
        thread_args[i].align_len = align_len;

        thread_args[i].bod_data_all_probe = bod_data_all;
        thread_args[i].fam_index_array = fam_index_array;
        thread_args[i].oii_index_arrya = oii_index_array;

        //need test mallc exit status
        thread_args[i].geno_data_this_variant = (double *)malloc(
            sizeof(double) * fam_node_num);
        thread_args[i].pheno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].geno_array = (double *)malloc(sizeof(double) * align_len);
        thread_args[i].result = (double *)malloc(sizeof(double) * probe_num * 8);
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
free_drm_threads_args_malloc(DRM_THREAD_ARGS_ptr thread_args, int thread_num)
{
    for (int i = 0; i < thread_num; i++) {
        free(thread_args[i].geno_data_this_variant);
        free(thread_args[i].g0_array);
        free(thread_args[i].g1_array);
        free(thread_args[i].g2_array);
        free(thread_args[i].geno_array);
        free(thread_args[i].pheno_array);
        free(thread_args[i].result);
    }
    return;
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
drm_thread_worker(void *args)
{
    
    DRM_THREAD_ARGS_ptr args_in = (DRM_THREAD_ARGS_ptr)args;
//    int thread_index = args_in->thread_index;

//    uint32_t fam_num = args_in->fam_num;
    uint32_t oii_num = args_in->oii_num;
    uint32_t probe_num = args_in->probe_num;
    uint32_t align_len = args_in->align_len;

    double *pheno_data_all = args_in->bod_data_all_probe;
    uint32_t *fam_index_array = args_in->fam_index_array;
    uint32_t *oii_index_array = args_in->oii_index_arrya;

    double *geno_data_one_varaint = args_in->geno_data_this_variant;
    double *pheno_data_one_probe = NULL;

    double *geno_data_aligned = args_in->geno_array;
    double *pheno_data_aligned = args_in->pheno_array;
    double *result = args_in->result;

    double c1_res = 0, stdev_res = 0, t1_res = 0, p_value_res = 0;
    double geno_0_median = 0.0, geno_1_median = 0.0, geno_2_median = 0.0;
    int g0_num = 0, g1_num = 0, g2_num = 0;
    double *g0_array = args_in->g0_array;
    double *g1_array = args_in->g1_array;
    double *g2_array = args_in->g2_array;
    double *g0_array_tmp = NULL;
    double *g1_array_tmp = NULL;
    double *g2_array_tmp = NULL;
    int align_len_rm_missing = 0;
    int geno_index = 0;
    int pheno_index = 0;
    double geno_value;
    double pheno_value;
#if defined DEBUG
    printf("thread_index: %d, align_len: %u, pheno_data_all: %p, fam_index_array: %p, "
        "oii_index_array: %p, geno_data_one_variant: %p, geno_data_aligned: %p, "
        "pheno_data_aligned: %p, result: %p, g0_array: %p, g1_array: %p,"
        " g2_array: %p\n",
        thread_index, align_len, pheno_data_all, fam_index_array, oii_index_array,
        geno_data_one_varaint, geno_data_aligned, pheno_data_aligned, result,
        g0_array, g1_array, g2_array);
#endif
    for (int i = 0; i < probe_num; i++) {
        //printf("probe_index: %d\n", i);
        pheno_data_one_probe = pheno_data_all + i * oii_num;
        align_len_rm_missing = 0;
        g0_array_tmp = g0_array;
        g1_array_tmp = g1_array;
        g2_array_tmp = g2_array;
        g0_num = 0;
        g1_num = 0;
        g2_num = 0;

        for (int j = 0; j < align_len; j++) {
            geno_index = fam_index_array[j];
            geno_value = geno_data_one_varaint[geno_index];
            pheno_index = oii_index_array[j];
            pheno_value = pheno_data_one_probe[pheno_index];

            // do I need use other way to compare float number?
            if (geno_value != 4.0 && pheno_value != -9.0) {
                geno_data_aligned[align_len_rm_missing] = geno_value;
                pheno_data_aligned[align_len_rm_missing] = pheno_value;
                align_len_rm_missing ++;
            }
        }

#if defined DEBUG       
            for (int i = 0; i < align_len_rm_missing; i++) {
                printf("%d %d %lf %lf\n", thread_index, i, geno_data_aligned[i],
                       pheno_data_aligned[i]);
            }
#endif        

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
        
        //printf("%u %u %u\n", g0_num, g1_num, g2_num);

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
        //printf("%lf %lf %lf \n", geno_0_median, geno_1_median, geno_2_median);

        //substract conressponding median from pheno value. and use absulute value.
        for (int i = 0; i < align_len_rm_missing; i++) {
            pheno_data_aligned[i] = (geno_data_aligned[i] == 0.0)
                ? fabs(pheno_data_aligned[i] - geno_0_median)
                : ((geno_data_aligned[i] == 1.0)
                    ? fabs(pheno_data_aligned[i] - geno_1_median)
                    : fabs(pheno_data_aligned[i] - geno_2_median));
        }

        linner_regression(geno_data_aligned, pheno_data_aligned, align_len_rm_missing,
            &c1_res, &stdev_res, &t1_res, &p_value_res);
        //c1_res beta1 beta1_se t_value_of_beta1, f_test_p
        //printf("%lf %lf %lf %lf\n", c1_res, stdev_res, t1_res, p_value_res);
        result[i * 4] = c1_res;
        result[i * 4 + 1] = stdev_res;
        result[i * 4 + 2] = t1_res;
        result[i * 4 + 3] = p_value_res;
        //printf("%le\n", p_value_res);
    }
#if defined DEBUG
    for (int i = 0; i < probe_num; i++) {
        printf("%lf %lf %lf %lf\n", result[4*i], result[4*i+1], result[4*i+2], result[4*i+3]);
    }
#endif

    return NULL;
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

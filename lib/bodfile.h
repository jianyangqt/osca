#ifndef BODFILE_HEAD
#define BODFILE_HEAD

#ifndef BODFILE_SRC
#define BODFILE_EXTERN extern
#else
#define BODFILE_EXTERN
#endif

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define BOD_OPEN_SUCCESS 0
#define BOD_OPEN_FAIL 1
#define BOD_FILE_OPEN_FAIL 2
#define BOD_FILE_EMPTY 3
#define BOD_READ_FAIL 4
#define BOD_VALUE_TYPE_UNKNOW 5
#define BOD_DATA_TYPE_UNKNOW 6
#define BOD_NUMBER_NOT_CONSISTENT 7
#define BOD_MALLOC_FAIL 8
#define BOD_LINE_BUFFER_OVERFLOW 9
#define BOD_FIELD_BUFFER_OVERFLOW 10
#define BOD_FIELD_SPLIT_FAIL 11
#define BOD_FIELD_UNKNOW 12
#define BOD_LINE_NOT_EAQUAL -1
#define BOD_EOF -1


#define BOD_VALUE_TYPE_MET 0
#define BOD_VALUE_TYPE_MET_M 1
#define BOD_VALUE_TYPE_OTHER 2

#define BOD_DATA_TYPE_GENE_EXP 0
#define BOD_DATA_TYPE_MET 1
#define BOD_DATA_TYPE_OTHER 2

#define BOD_LINE_BUFFER_LEN 2048

#define BOD_OII_NA 0
#define BOD_OPI_NA 0

#define BOD_MISSING_CHROM 0
#define BOD_X_CHROM 201
#define BOD_Y_CHROM 202
#define BOD_MT_CHROM 203
#define BOD_CHROM_UNKNOW 13

#define BOD_ORIEN_POS 0
#define BOD_ORIEN_NEG 1
#define BOD_ORIEN_NA 2

typedef struct bodfile {
    uint32_t individual_num;
    uint32_t probe_num;

    FILE *oii_file;
    FILE *opi_file;
    FILE *bod_file;

    uint32_t current_oii_line_index;
    uint32_t current_opi_line_index;
    uint32_t current_bod_data_index;

    uint64_t oii_offset_byte;
    uint64_t opi_offset_byte;
    uint64_t bod_offset_byte;

    char bod_value_type;
    char bod_data_type;

    uint32_t bod_line_buf_len;
    char *bod_line_buf;

    char status;

} BODFILE, *BODFILE_ptr;


typedef struct oii_line_stu {
    char family_id[64];
    char indiv_id[64];
    char parental_id[64];
    char maternal_id[64];
    char sex;
} OII_LINE, *OII_LINE_ptr;


typedef struct opi_line_stu {
    unsigned char chrom;  // 0 for NA, 201 for X, 202 for Y, 203 for MT
    char probe_id[64];
    uint32_t position;  // 0 for NA
    char gene_id[64];
    char ori;  // 0 for NA, 1 for +, 2 for -
} OPI_LINE, *OPI_LINE_ptr;

BODFILE_EXTERN BODFILE bodfileopen(const char *filename);
BODFILE_EXTERN int bodfileclose(BODFILE_ptr bod_data);
BODFILE_EXTERN void bodfilerewind(BODFILE_ptr bod_data);
BODFILE_EXTERN int bodfileseek(BODFILE_ptr bod_data, uint32_t seek_len);

BODFILE_EXTERN int oiireadline(BODFILE_ptr bod_data, OII_LINE_ptr oii_line);
BODFILE_EXTERN int oiireadlines(BODFILE_ptr bod_data, OII_LINE_ptr oii_lines,
                 uint32_t oii_line_num);

BODFILE_EXTERN int opireadline(BODFILE_ptr bod_data, OPI_LINE_ptr opi_line);
BODFILE_EXTERN int opireadlines(BODFILE_ptr bod_data, OPI_LINE_ptr opi_lines,
                 uint32_t line_num);

BODFILE_EXTERN int bodreaddata(BODFILE_ptr bod_data, double *bod_readout, uint64_t readout_len);
BODFILE_EXTERN int bodloaddata_n(BODFILE_ptr bod_data, double *bod_readoutn,
                  uint64_t readout_len, uint32_t start_offset, uint32_t load_len);
BODFILE_EXTERN int bodloaddata_all(BODFILE_ptr bod_data, double *bod_readoutall,
                    uint64_t readout_len);

#endif


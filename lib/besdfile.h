#ifndef BESDFILE_HEAD
#define BESDFILE_HEAD
#include <stdio.h>
#include <stdint.h>


#ifndef BESDFILE_SRC
#define BESDFILE_EXTERN extern
#else
#define BESDFILE_EXTERN
#endif

BESDFILE_EXTERN int besd_sparse_write_meta(int file_format, int32_t sample_size,
                        uint32_t variant_num,
                        uint32_t probe_num,
                        uint64_t vaule_num, uint64_t *offset,
                        FILE *fout);
BESDFILE_EXTERN int besd_sparse_write_variant_index(uint32_t *index, uint32_t index_len,
                                    FILE *fout);

BESDFILE_EXTERN int besd_sparse_write_beta_se_data(float *beta, float *se,
                                                   uint32_t data_len,
                                                   FILE *fout);

#endif
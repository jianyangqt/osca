#define BESDFILE_SRC

#include "besdfile.h"

int
besdfileopen(const char *fname, void *besd_data)
{

    return 0;
}

void
besdfilerewind(void *besd_data)
{

    return;
}


int
besdfileseek(void *besd_data)
{

    return 0;
}


int
epireadline(void)
{
    return 0;
}


int
epireadlines(void)
{
    return 0;
}


int
esireadline(void)
{
    return 0;
}


int
esireadlines(void)
{
    return 0;
}


int
besd_sparse_write_meta(int file_format, int32_t sample_size, uint32_t variant_num,
    uint32_t probe_num, uint64_t vaule_num, uint64_t *offset, FILE *fout)
{
    int meta16[16];
    meta16[0] = file_format;
    meta16[1] = sample_size;
    meta16[2] = variant_num;
    meta16[3] = probe_num;

    for (int i = 4; i < 16; i++) {
        meta16[i] = -9;
    }

    fwrite(meta16, sizeof(int32_t), 16, fout);
    fwrite(&vaule_num, sizeof(uint64_t), 1, fout);
    fwrite(offset, sizeof(uint64_t), probe_num * 2 + 1, fout);
    return 0;
}


int
besd_sparse_write_variant_index(uint32_t *index, uint32_t index_len, FILE *fout)
{
    fwrite(index, sizeof(uint32_t), index_len, fout);
    fwrite(index, sizeof(uint32_t), index_len, fout);
    return 0;
}


int
besd_sparse_write_beta_se_data(float *beta, float *se, uint32_t data_len, FILE *fout)
{
    fwrite(beta, sizeof(float), data_len, fout);
    fwrite(se, sizeof(float), data_len, fout);
    return 0;
}


int
besd_sparse_write(void)
{

    /*
        method one,
        creat tmp files, and then concatnate them.
    */

    return 0;

}


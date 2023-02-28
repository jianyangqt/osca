#define PLINKLITE_SRC

/*
    plinklite is lightwight library to read and write plink bed file.

    -Benjamin Fang 20230214
*/

#include "plinklite.h"


PLINKFILE
plinkopen(const char *filename)
{
    PLINKFILE data_out;
    data_out.fam_file = NULL;
    data_out.bim_file = NULL;
    data_out.bed_file = NULL;
    data_out.individual_num = 0;
    data_out.variant_num = 0;
    data_out.current_fam_line_index = 0;
    data_out.current_bim_line_index = 0;
    data_out.current_bed_data_index = 0;

    data_out.bed_byte_offset = 0;

    data_out.status = PLINK_FAIL;
    data_out.raw_buf_len = 0;
    data_out.raw_buf = NULL;
    data_out.decode_buf_len = 0;
    data_out.decoded_buf = NULL;
    data_out.line_buf_len = PLINK_LINE_BUF_LEN;
    data_out.line_buf = NULL;

    int filename_len = strlen(filename);
    char *filename_full = (char *)malloc(sizeof(char) * (filename_len + 5));
    strcpy(filename_full, filename);

    //open fam file
    FILE *fin = NULL;
    strcat(filename_full, ".fam");
    fin = fopen(filename_full, "r");
    if (!fin) {
        fprintf(stderr, "open fam file failed.\n");
        data_out.status = PLINK_FILE_OPEN_FAIL;
        return data_out;
    }
    data_out.fam_file = fin;
    
    //open bim file
    fin = NULL;
    filename_full[filename_len] = '\0';
    strcat(filename_full, ".bim");
    fin = fopen(filename_full, "r");
    if (!fin) {
        fprintf(stderr, "open bim file failed.\n");
        data_out.status = PLINK_FILE_OPEN_FAIL;
        return data_out;
    }
    data_out.bim_file = fin;
    
    //open bed file
    fin = NULL;
    filename_full[filename_len] = '\0';
    strcat(filename_full, ".bed");
    fin = fopen(filename_full, "r");
    if (!fin) {
        fprintf(stderr, "open bed file failed.\n");
        data_out.status = PLINK_FILE_OPEN_FAIL;
        return data_out;
    }
    data_out.bed_file = fin;
    free(filename_full);

    //get  individual number.
    uint32_t line_counter = 0;
    int last_char = -1;
    int cc = 0;
    fin = data_out.fam_file;
    while ((cc = fgetc(fin)) != EOF) {
        if (cc == '\n') {
            line_counter++;
        }
        last_char = cc;
    }
    if (last_char != '\n') {
        fprintf(stderr, "Warning, file not end by new line character.\n");
        line_counter++;
    }
    if (line_counter == 0) {
        fprintf(stderr, "fam file is empty.\n");
        data_out.status = PLINK_FILE_EMPTY;
        return data_out;
    }
    data_out.individual_num = line_counter;
    rewind(fin);

    //get variant number.
    line_counter = 0;
    last_char = -1;
    fin = data_out.bim_file;
    while ((cc = fgetc(fin)) != EOF) {
        if (cc == '\n') {
            line_counter++;
        }
        last_char = cc;
    }

    if (last_char != '\n') {
        fprintf(stderr, "Warning, file not end by new line character.\n");
        line_counter++;
    }
    if (line_counter == 0) {
        fprintf(stderr, "empty bim file.\n");
        data_out.status = PLINK_FILE_EMPTY;
        return data_out;
    }
    data_out.variant_num = line_counter;
    rewind(fin);

    //handle bed file.
    char megic_num[3];
    fin = data_out.bed_file;
    if (fread(megic_num, sizeof(char), PLINK_MEGIC_NUM_LEN, fin) != 3) {
        fprintf(stderr, "read megic num failed.\n");
        data_out.status = PLINK_FILE_READ_FAIL;
        return data_out;
    }
    memcpy(data_out.bed_megic_num, megic_num, PLINK_MEGIC_NUM_LEN);
    data_out.bed_byte_offset = PLINK_MEGIC_NUM_LEN;

    int bed_raw_buf_len = 0;
    bed_raw_buf_len = ceil((float)data_out.individual_num / 4.0);
    int bed_decode_buf_len = bed_raw_buf_len * 4;
    
    char *buf = NULL;
    buf = (char *)malloc(bed_raw_buf_len * sizeof(char));
    if (!buf) {
        fprintf(stderr, "malloc plink internal buffer failed.\n");
        data_out.status = PLINK_MALLOC_BUF_FAIL;
        return data_out;
    }
    data_out.raw_buf = buf;
    data_out.raw_buf_len = bed_raw_buf_len;

    buf = NULL;
    buf = (char *)malloc(bed_decode_buf_len * sizeof(char));
    if (!buf) {
        fprintf(stderr, "malloc plink internal buffer failed.\n");
        data_out.status = PLINK_MALLOC_BUF_FAIL;
        return data_out;
    }
    data_out.decoded_buf = buf;
    data_out.decode_buf_len = bed_decode_buf_len;

    buf = NULL;
    buf = (char *)malloc(sizeof(char) * PLINK_LINE_BUF_LEN);
    if (!buf) {
        fprintf(stderr, "malloc plink intenal buffer failed.\n");
        data_out.status = PLINK_MALLOC_BUF_FAIL;
        return data_out;
    }
    data_out.line_buf = buf;
    
    data_out.status = PLINK_SUCCESS;
    return data_out;
}


/*
    Finalize the plink data.
*/
int
plinkclose(PLINKFILE_ptr plink_data)
{
    int status = 0;
    if (plink_data->fam_file) {
        status += fclose(plink_data->fam_file);
        plink_data->fam_file = NULL;
    }

    if (plink_data->bim_file) {
        status += fclose(plink_data->bim_file);
        plink_data->bim_file = NULL;
    }

    if (plink_data->bed_file) {
        status += fclose(plink_data->bed_file);
        plink_data->bed_file = NULL;
    }

    if (plink_data->raw_buf) {
        free(plink_data->raw_buf);
        plink_data->raw_buf = NULL;
    }

    if (plink_data->decoded_buf) {
        free(plink_data->decoded_buf);
        plink_data->decoded_buf = NULL;
    }

    if (plink_data->line_buf) {
        free(plink_data->line_buf);
        plink_data->line_buf = NULL;
    }
    return status;
}


/*
    Rewind the plink data to its begin.
*/
void
plinkrewind(PLINKFILE_ptr plink_data)
{

    rewind(plink_data->fam_file);
    plink_data->current_fam_line_index = 0;

    rewind(plink_data->bim_file);
    plink_data->current_bim_line_index = 0;

    fseek(plink_data->bed_file, PLINK_MEGIC_NUM_LEN, SEEK_SET);
    plink_data->current_bed_data_index = 0;
    plink_data->bed_byte_offset = PLINK_MEGIC_NUM_LEN;

    return;
}


/*
    Seek the plink data to a specific position from begin.
    0 returned for success,
    null 0 for fail.
*/
int
plinkseek(PLINKFILE_ptr plink_data, uint32_t seek_len)
{
    plinkrewind(plink_data);
    FAM_LINE fam_line;
    BIM_LINE bim_line;
    uint32_t indi_num = plink_data->individual_num;
    char *bed_data = (char *)malloc(sizeof(char) * indi_num);
    for (uint32_t i = 0; i < seek_len; i++) {
        famreadline(plink_data, &fam_line);
        bimreadline(plink_data, &bim_line);
        bedreaddata(plink_data, bed_data, indi_num);
    }
    free(bed_data);
    return 0;
}


int
famreadline(PLINKFILE_ptr plink_data, FAM_LINE_ptr fam_line)
{
   
    FILE *fin = plink_data->fam_file;
    char *line_buf = plink_data->line_buf;
    line_buf[PLINK_LINE_BUF_LEN - 1] = 1;

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


    if (fgets(line_buf, PLINK_LINE_BUF_LEN, fin)) {

        
        if (line_buf[PLINK_LINE_BUF_LEN - 1] != 1) {
            fprintf(stderr, "fam line buffer overflow.\n");
            return 1;
        }
        if (sscanf(line_buf, "%s %s %s %s %s %s", fam_id, within_famid,
                   father_id, mother_id, sex, pheno_val) != 6) {
            fprintf(stderr, "split fam line field failed.\n");
            return 1;
        }

        if (fam_id[63] == '\0') {
            strcpy(fam_line->family_id, fam_id);
        } else {
            fprintf(stderr, "family id field failed.\n");
            return 1;
        }

        if (within_famid[63] == '\0') {
            strcpy(fam_line->within_famid, within_famid);
        } else {
            fprintf(stderr, "within family id field failed.\n");
            return 1;
        }

        if (father_id[63] == '\0') {
            strcpy(fam_line->father_id, father_id);
        } else {
            fprintf(stderr, "father id field failed.\n");
            return 1;
        }

        if (mother_id[63] == '\0') {
            strcpy(fam_line->mother_id, mother_id);
        } else {
            fprintf(stderr, "mother id field failed.\n");
            return 1;
        }

        if (sex[15] == '\0') {
            if ((strlen(sex) == 1) && (isdigit(sex[0]))) {
                fam_line->sex = (char)atoi(sex);
            } else if (strcmp(sex, "NA") == 0) {
                fam_line->sex = PLINK_FAM_NA;
            } else {
                fprintf(stderr, "sex field is not recognized.\n");
                return 1;
            }
        } else {
            fprintf(stderr, "sex field failed.\n");
            return 1;
        }

        if (pheno_val[15] == '\0') {
            if ((strlen(pheno_val) == 1) && (isdigit(pheno_val[0]))) {
                fam_line->phenotype_value = (char)atoi(pheno_val);
            } else if (strcmp(pheno_val, "-9") == 0) {
                fam_line->phenotype_value = PLINK_FAM_NA;
            } else if (strcmp(pheno_val, "NA") == 0) {
                fam_line->phenotype_value = PLINK_FAM_NA;
            } else {
                fprintf(stderr, "phenotype value field failed.\n");
                return 1;
            }
        } else {
            fprintf(stderr, "phenotype value field failed.\n");
            return 1;
        }
    } else {
        return -1;
    }
    plink_data->current_fam_line_index += 1;
    return 0;
}


int
famreadlines(PLINKFILE_ptr plink_data, FAM_LINE_ptr fam_lines, uint32_t fam_line_num)
{

    FILE *fin = plink_data->fam_file;

    rewind(fin);
    plink_data->current_fam_line_index = 0;

    if (fam_line_num != plink_data->individual_num) {
        fprintf(stderr, "fam line num do not match.\n");
        return 1;
    }

    int line_count = 0;
    while (famreadline(plink_data, &(fam_lines[line_count])) == 0) {
        line_count++;
    }

    return line_count;
}


int
bimreadline(PLINKFILE_ptr plink_dara, BIM_LINE_ptr bim_line)
{

    FILE *fin = plink_dara->bim_file;
    char *line_buf = plink_dara->line_buf;
    line_buf[PLINK_LINE_BUF_LEN - 1] = 1;

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

    if (fgets(line_buf, PLINK_LINE_BUF_LEN, fin)) {
        
        if (line_buf[PLINK_LINE_BUF_LEN - 1] != 1) {
            fprintf(stderr, "bim file line buf overflow.\n");
            return 1;
        }
        if (sscanf(line_buf, "%s %s %s %s %s %s", chrom, rsid, phy_pos, pos,
            allel1, allel2) != 6) {
            fprintf(stderr, "split bim line fields failed.\n");
            return 1;
        }

        if (chrom[15] == '\0') {
            str_len = strlen(chrom);
            char_ptr = chrom;
            while (isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                bim_line->chrom = (unsigned char)atoi(chrom);
            } else if (strcmp(chrom, "X") == 0) {
                bim_line->chrom = PLINK_X_CHROM;
            } else if (strcmp(chrom, "Y") == 0) {
                bim_line->chrom = PLINK_Y_CHROM;
            } else if (strcmp(chrom, "MT") == 0) {
                bim_line->chrom = PLINK_MT_CHROM;
            } else if (strcmp(chrom, "NA") == 0) {
                bim_line->chrom = PLINK_MISSING_CHROM;
            } else {
                fprintf(stderr, "chrom not reganized.\n");
                return 1;
            }
        }

        if (rsid[63] == '\0') {
            strcpy(bim_line->rsid, rsid);
        } else {
            fprintf(stderr, "rsid field failed.\n");
            return 1;
        }
        if (phy_pos[15] == '\0') {
            bim_line->phy_pos = (float)atof(phy_pos);
        } else {
            fprintf(stderr, "physical position field failed.\n");
            return 1;
        }
        if (pos[15] == '\0') {
            str_len = strlen(pos);
            char_ptr = pos;
            while (isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                bim_line->pos = (uint32_t)atol(pos);
            } else if (strcmp(pos, "NA") == 0) {
                bim_line->pos = PLINE_BIM_POS_NA;
            } else {
                fprintf(stderr, "bim pos not recognized.\n");
                return 1;
            }
        } else {
            fprintf(stderr, "bim position field failed.\n");
            return 1;
        }
        if (allel1[1023] == '\0') {
            str_len = strlen(allel1);
            if (str_len >= PLINE_MAX_ALLEL_LEN) {
                //fprintf(stderr, "allel %s trimed\n", allel1);
                ;
            }
            (bim_line->allel1)[PLINE_MAX_ALLEL_LEN - 1] = '\0';
            strncpy(bim_line->allel1, allel1, PLINE_MAX_ALLEL_LEN - 1);
            
        } else {
            fprintf(stderr, "bim allel failed.\n");
            return 1;
        }

        if (allel2[1023] == '\0') {
            str_len = strlen(allel2);
            if (str_len >= PLINE_MAX_ALLEL_LEN) {
                //fprintf(stderr, "allel %s trimed\n", allel2);
                ;
            }
            bim_line->allel2[PLINE_MAX_ALLEL_LEN - 1] = '\0';
            strncpy(bim_line->allel2, allel2, PLINE_MAX_ALLEL_LEN - 1);
        } else {
            fprintf(stderr, "bim allel failed.\n");
            return 1;
        }

    } else {
        return -1;
    }

    plink_dara->current_bim_line_index++;
    return 0;
}


int
bimreadlines(PLINKFILE_ptr plink_data, BIM_LINE_ptr bim_lines, uint32_t bim_line_num)
{

    FILE *fin = plink_data->bim_file;

    rewind(fin);
    plink_data->current_bim_line_index = 0;

    if (bim_line_num != plink_data->variant_num) {
        fprintf(stderr, "bim line number do not match.\n");
        return 1;
    }

    int line_count = 0;
    while (bimreadline(plink_data, &(bim_lines[line_count])) == 0) {
        line_count++;
    }

    return line_count;
}


int
bedreaddata(PLINKFILE_ptr plink_data, char *bed_data, uint64_t bed_data_len)
{
    
    FILE *fin = plink_data->bed_file;
    if (bed_data_len != plink_data->individual_num) {
        fprintf(stderr, "bed data length error.\n");
        return 1;
    }

    char *raw_buf = plink_data->raw_buf;
    char *decode_buf = plink_data->decoded_buf;
    int raw_len = plink_data->raw_buf_len;
    int decode_len = plink_data->decode_buf_len;
    int readlen = 0;

    if ((readlen = fread(raw_buf, sizeof(char), raw_len, fin)) != raw_len) {
        if (readlen > 0 && readlen < raw_len) {
            fprintf(stderr, "readlen not equal to expected. %d\n", readlen);
            return 1;
        } else if (readlen == 0) {
            fprintf(stderr, "Warning, EOF reached. %d\n", readlen);
            return -1;
        }

    }

    int i = 0, j = 0;
    uint8_t char_tmp;
    char Low_1_2b = 0, Low_2_2b = 0, Low_3_2b = 0, Low_4_2b = 0;
    for (i = 0; i < raw_len; i++) {
        char_tmp = raw_buf[i];
        Low_1_2b = char_tmp & CharM1;
        Low_2_2b = (char_tmp & CharM2) >> 2;
        Low_3_2b = (char_tmp & CharM3) >> 4;
        Low_4_2b = char_tmp >> 6;
        decode_buf[j] = CONVERT_GENO(Low_1_2b);
        decode_buf[++j] = CONVERT_GENO(Low_2_2b);
        decode_buf[++j] = CONVERT_GENO(Low_3_2b);
        decode_buf[++j] = CONVERT_GENO(Low_4_2b);
        j++;
    }

    if (j != decode_len) {
        fprintf(stderr, "bed decode failed.\n");
        return 1;
    }

    for (int i = 0; i < bed_data_len; i++) {
        bed_data[i] = decode_buf[i];
    }

    plink_data->current_bed_data_index++;
    plink_data->bed_byte_offset += raw_len;

    return 0;
}


int
bedloaddata_n(PLINKFILE_ptr plink_data, char *bed_data, uint64_t bed_data_len,
    int start_offset, int load_length)
{
    FILE *fin = plink_data->bed_file;
    uint32_t indiv_num = plink_data->individual_num;
    uint64_t current_byte_offset = plink_data->bed_byte_offset;
    uint32_t current_data_index = plink_data->current_bed_data_index;
    uint32_t row_len = plink_data->raw_buf_len;

    if ((load_length * indiv_num) != bed_data_len) {
        fprintf(stderr, "bed data length error.\n");
        return 1;
    }
    uint64_t seek_len = start_offset * row_len + 3;
    fseek(fin, seek_len, SEEK_SET);
    
    plink_data->bed_byte_offset = seek_len;
    plink_data->current_bed_data_index = start_offset;
    for (int i = 0; i < load_length; i++) {
        int status = 0;
        if ((status = bedreaddata(plink_data, bed_data, indiv_num)) == 0) {
            bed_data += indiv_num;

        } else if (status == 1) {
            fprintf(stderr, "read bed failed.\n");
            return 1;
        } else if (status == -1) {
            fprintf(stderr, "Warning, read out of bunderay.\n");
            return -1;
        }
    }
    fseek(fin, current_byte_offset, SEEK_SET);
    plink_data->current_bed_data_index = current_data_index;
    plink_data->bed_byte_offset = current_byte_offset;

    return 0;
}


int
bedloaddata_all(PLINKFILE_ptr plink_data, char *bed_data, uint64_t bed_data_len)
{
    uint32_t variant_num = plink_data->variant_num;
    uint32_t indi_num = plink_data->individual_num;
    if (bed_data_len != (variant_num * indi_num)) {
        fprintf(stderr, "bed mem allocation may not correct.\n");
        return 1;
    }
    if (bedloaddata_n(plink_data, bed_data, bed_data_len, 0, variant_num) != 0) {
        return 1;
    }
    return 0;
}




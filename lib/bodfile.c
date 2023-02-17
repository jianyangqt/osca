#define BODFILE_SRC

#include "bodfile.h"

BODFILE
bodfileopen(const char *filename)
{

    BODFILE data_out;
    data_out.individual_num = 0;
    data_out.probe_num = 0;

    data_out.oii_file = NULL;
    data_out.opi_file = NULL;
    data_out.bod_file = NULL;

    data_out.current_bod_data_index = 0;
    data_out.current_oii_line_index = 0;
    data_out.current_opi_line_index = 0;

    data_out.oii_offset_byte = 0;
    data_out.opi_offset_byte = 0;
    data_out.bod_offset_byte = 0;

    data_out.bod_line_buf_len = 0;
    data_out.bod_line_buf = NULL;

    int filename_len = strlen(filename);
    char *filename_full = (char *)malloc(sizeof(char) * (filename_len + 5));
    strcpy(filename_full, filename);

    //open oii file
    FILE *fin = NULL;
    strcat(filename_full, ".oii");
    fin = fopen(filename_full, "r");
    if (!fin) {
        fprintf(stderr, "open oii file failed.\n");
        data_out.status = BOD_FILE_OPEN_FAIL;
        return data_out;
    }
    data_out.oii_file = fin;


    //open opi file
    fin = NULL;
    filename_full[filename_len] = '\0';
    strcat(filename_full, ".opi");
    fin = fopen(filename_full, "r");
    if (!fin) {
        fprintf(stderr, "open opi file failed.\n");
        data_out.status = BOD_FILE_OPEN_FAIL;
        return data_out;
    }
    data_out.opi_file = fin;

    //open bod file
    fin = NULL;
    filename_full[filename_len] = '\0';
    strcat(filename_full, ".bod");
    fin = fopen(filename_full, "r");
    if (!fin) {
        fprintf(stderr, "open bod file failed.\n");
        data_out.status = BOD_FILE_OPEN_FAIL;
        return data_out;
    }
    data_out.bod_file = fin;

    //get individual number
    uint32_t line_counter = 0;
    int last_char = -1;
    int cc = 0;
    fin = data_out.oii_file;
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
        data_out.status = BOD_FILE_EMPTY;
        return data_out;
    }
    data_out.individual_num = line_counter;
    rewind(fin);

    //get probe number
    fin = data_out.opi_file;
    line_counter = 0;
    last_char = -1;

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
        data_out.status = BOD_FILE_EMPTY;
        return data_out;
    }
    data_out.probe_num = line_counter;
    rewind(fin);

    //handle bod file
    fin = data_out.bod_file;
    char read_buf[12];
    if (fread(read_buf, sizeof(char), 12, fin) != 12) {
        fprintf(stderr, "read bod first 12 bytes failed.\n");
        data_out.status = BOD_READ_FAIL;
        return data_out;
    }

    if (read_buf[0] == BOD_VALUE_TYPE_MET) {
        data_out.bod_value_type = BOD_VALUE_TYPE_MET;
    } else if (read_buf[0] == BOD_VALUE_TYPE_MET_M) {
        data_out.bod_data_type = BOD_VALUE_TYPE_MET_M;
    } else if (read_buf[0] == BOD_VALUE_TYPE_OTHER) {
        data_out.bod_value_type = BOD_DATA_TYPE_OTHER;
    } else {
        fprintf(stderr, "value type unknow.\n");
        data_out.status = BOD_VALUE_TYPE_UNKNOW;
        return data_out;
    }

    if (read_buf[1] == BOD_DATA_TYPE_GENE_EXP) {
        data_out.bod_data_type = BOD_DATA_TYPE_GENE_EXP;
    } else if (read_buf[1] == BOD_DATA_TYPE_MET) {
        data_out.bod_data_type = BOD_DATA_TYPE_MET;
    } else if (read_buf[1] == BOD_DATA_TYPE_OTHER) {
        data_out.bod_data_type = BOD_DATA_TYPE_OTHER;
    } else {
        fprintf(stderr, "data type unknow.\n");
        data_out.status = BOD_DATA_TYPE_UNKNOW;
        return data_out;
    }

    uint32_t *int_ptr = NULL;
    int_ptr = (uint32_t *)(read_buf + 4);
    if (*int_ptr != data_out.individual_num) {
        fprintf(stderr, "individual number not consistent.\n");
        data_out.status = BOD_NUMBER_NOT_CONSISTENT;
        return data_out;
    }

    int_ptr = (uint32_t *)(read_buf + 8);
    if (*int_ptr != data_out.probe_num) {
        fprintf(stderr, "probe number not consistent.\n");
        data_out.status = BOD_NUMBER_NOT_CONSISTENT;
        return data_out;
    }

    data_out.bod_offset_byte = 12;


    //allocate buffer
    char *line_buf = (char *)malloc(sizeof(char) * BOD_LINE_BUFFER_LEN);
    if (!line_buf) {
        fprintf(stderr, "malloc for bod line buffer failed.\n");
        data_out.status = BOD_MALLOC_FAIL;
        return data_out;
    }
    data_out.bod_line_buf = line_buf;
    data_out.bod_line_buf_len = BOD_LINE_BUFFER_LEN;

    return data_out;
}


void
bodfileclose(BODFILE_ptr bod_data)
{
    if (bod_data->oii_file) {
        fclose(bod_data->oii_file);
        bod_data->oii_file = NULL;
    }

    if (bod_data->opi_file) {
        fclose(bod_data->opi_file);
        bod_data->opi_file = NULL;
    }

    if (bod_data->bod_file) {
        fclose(bod_data->bod_file);
        bod_data->bod_file = NULL;
    }


    if (bod_data->bod_line_buf) {
        free(bod_data->bod_line_buf);
        bod_data->bod_line_buf = NULL;
    }
    return;
}


int
oiireadline(BODFILE_ptr bod_data, OII_LINE_ptr oii_line)
{
    FILE *fin = bod_data->oii_file;
    char *line_buf = bod_data->bod_line_buf;
    line_buf[BOD_LINE_BUFFER_LEN - 1] = '\0';

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

    if (fgets(line_buf, BOD_LINE_BUFFER_LEN, fin)) {

        if (line_buf[BOD_LINE_BUFFER_LEN - 1] != '\0') {
            fprintf(stderr, "line buffer overflow.\n");
            return BOD_LINE_BUFFER_OVERFLOW;
        }

        if (sscanf(line_buf, "%s %s %s %s %s", family_id, indiv_id,
            parental_id, maternal_id, sex) != 6) {
            fprintf(stderr, "line buffer split error.\n");
            return BOD_FIELD_SPLIT_FAIL;
        }

        if (family_id[63] == '\0') {
            strcpy(oii_line->family_id, family_id);
        } else {
            fprintf(stderr, "family id field overflow.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }
        if (indiv_id[63] == '\0') {
            strcpy(oii_line->indiv_id, indiv_id);
        } else {
            fprintf(stderr, "oii individual field overflow.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }
        if (parental_id[63] == '\0') {
            strcpy(oii_line->parental_id, parental_id);
        } else {
            fprintf(stderr, "parental id field failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }
        if (maternal_id[63] == '\0') {
            strcpy(oii_line->maternal_id, maternal_id);
        } else {
            fprintf(stderr, "maternal id field failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }
        if (sex[63] == '\0') {
            if (strcmp(sex, "NA") == 0) {
                oii_line->sex = BOD_OII_NA;
            } else if (strlen(sex) == 1 && isdigit(*sex)) {
                oii_line->sex = (char)atoi(sex);
            } else {
                fprintf(stderr, "sex field not recgnized.\n");
                return BOD_FIELD_UNKNOW;
            }
        } else {
            fprintf(stderr, "sex field failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }

    } else {
        return BOD_EOF;
    }
    bod_data->current_oii_line_index++;
    return 0;
}


int
oiireadlines(BODFILE_ptr bod_data, OII_LINE_ptr oii_lines, uint32_t oii_line_num)
{
    FILE *fin = bod_data->oii_file;
    bod_data->current_oii_line_index = 0;
    rewind(fin);
    if (oii_line_num != bod_data->individual_num) {
        fprintf(stderr, "oii line number not match.\n");
        return BOD_LINE_NOT_EAQUAL;
    }

    int line_counter = 0;
    while (oiireadline(bod_data, &(oii_lines[line_counter])) == 0) {
        line_counter++;
    }
    return line_counter;
}


int
opireadline(BODFILE_ptr bod_data, OPI_LINE_ptr opi_line)
{
    FILE *fin = bod_data->opi_file;
    char *line_buf = bod_data->bod_line_buf;
    line_buf[BOD_LINE_BUFFER_LEN - 1] = '\0';

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

    char *char_ptr = NULL;
    int str_len = 0;

    if (fgets(line_buf, BOD_LINE_BUFFER_LEN, fin)) {
        if (line_buf[BOD_LINE_BUFFER_LEN - 1] != '\0') {
            fprintf(stderr, "read oii line buffer overflow.\n");
            return BOD_LINE_BUFFER_OVERFLOW;
        }
        if (sscanf(line_buf, "%s %s %s %s %s", chrom, probe_id, position,
                   gene_id, oritation) != 5) {
            fprintf(stderr, "opi field failed.\n");
            return BOD_FIELD_SPLIT_FAIL;
        }

        if (chrom[15] == '\0') {
            str_len = strlen(chrom);
            char_ptr = chrom;
            while (isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                opi_line->chrom = (unsigned char)atoi(chrom);
            } else if (strcmp(chrom, "X") == 0) {
                opi_line->chrom = BOD_X_CHROM;
            } else if (strcmp(chrom, "Y") == 0) {
                opi_line->chrom = BOD_Y_CHROM;
            } else if (strcmp(chrom, "MT") == 0) {
                opi_line->chrom = BOD_MT_CHROM;
            } else if (strcmp(chrom, "NA") == 0) {
                opi_line->chrom = BOD_MISSING_CHROM;
            } else {
                fprintf(stderr, "chrom not reganized.\n");
                return BOD_CHROM_UNKNOW;
            }
        } else {
            fprintf(stderr, "field overflow.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }

        if (probe_id[63] == '\0') {
            strcpy(opi_line->probe_id, probe_id);
        } else {
            fprintf(stderr, "probe id field failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }

        if (position[15] == '\0') {
            char_ptr = position;
            str_len = strlen(position);
            while (isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                opi_line->position = (uint64_t)atol(position);
            } else if (strcmp(position, "NA") == 0) {
                opi_line->position = BOD_OPI_NA;
            } else {
                fprintf(stderr, "unknow position field.\n");
                return BOD_FIELD_UNKNOW;
            }
        } else {
            fprintf(stderr, "opi position failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }

        if (gene_id[63] == '\0') {
            strcpy(opi_line->gene_id, gene_id);
        } else {
            fprintf(stderr, "gene id field failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }

        if (oritation[15] == '\0') {
            if (strcmp(oritation, "+") == 0) {
                opi_line->ori = BOD_ORIEN_POS;
            } else if (strcmp(oritation, "-") == 0) {
                opi_line->ori = BOD_ORIEN_NEG;
            } else if (strcmp(oritation, "NA") == 0) {
                opi_line->ori = BOD_ORIEN_NA;
            } else {
                fprintf(stderr, "oritation not recognized.\n");
                return BOD_FIELD_UNKNOW;
            }
        } else {
            fprintf(stderr, "oritation field failed.\n");
            return BOD_FIELD_BUFFER_OVERFLOW;
        }
    } else {
        return BOD_EOF;
    }
    bod_data->current_opi_line_index++;
    return 0;
}

int
opireadlines(BODFILE_ptr bod_data, OPI_LINE_ptr opi_lines, uint32_t line_num)
{
    FILE *fin = bod_data->opi_file;
    rewind(fin);
    bod_data->current_opi_line_index = 0;

    uint32_t line_counter = 0;
    while (opireadline(bod_data, &(opi_lines[line_counter])) == 0) {
        line_counter++;
    }

    return line_counter;
}


int
bodreaddata(BODFILE_ptr bod_data, double *bod_readout, uint32_t readout_len)
{
    FILE *fin = bod_data->bod_file;
    uint32_t indi_num = bod_data->individual_num;
    if (indi_num != readout_len) {
        fprintf(stderr, "alloced bod readline is not equal individual length.\n");
        return BOD_LINE_NOT_EAQUAL;
    }

    if (fread(bod_readout, sizeof(double), readout_len, fin) != readout_len) {
        fprintf(stderr, "read bod file failed.\n");
        return BOD_READ_FAIL;
    }

    bod_data->bod_offset_byte += readout_len;
    bod_data->current_bod_data_index++;

    return 0;
}


int
bodloaddata_n(BODFILE_ptr bod_data, double *bod_readoutn, uint32_t readout_len,
    int start, int end)
{
    FILE *fin = bod_data->bod_file;
    uint32_t indi_num = bod_data->individual_num;

    if (readout_len != indi_num * (end - start + 1)) {
        fprintf(stderr, "read len error.\n");
        return 1;
    }

    uint32_t current_data_index = bod_data->current_bod_data_index;
    uint64_t current_offset = bod_data->bod_offset_byte;

    uint64_t seek_len = sizeof(double) * indi_num * (start - 1);

    fseek(fin, seek_len, SEEK_SET);

    for (int i = start; i <= end; i++) {
        int status = 0;
        if ((status = bodreaddata(bod_data, bod_readoutn, indi_num)) == 0) {
            bod_readoutn += indi_num;
        } else if (status == -1) {
            return -1;
        } else {
            return 1;
        }
    }
    fseek(fin, current_offset, SEEK_SET);
    bod_data->current_bod_data_index = current_data_index;
    bod_data->bod_offset_byte = current_offset;

    return 0;
}


int
bodloaddata_all(BODFILE_ptr bod_data, double *bod_readoutall, uint32_t readout_len)
{
    FILE *fin = bod_data->bod_file;

    uint32_t end = bod_data->probe_num;
    uint32_t start = 1;

    uint32_t indi_num = bod_data->individual_num;
    
    if (readout_len != (end - start + 1) * indi_num) {
        fprintf(stderr, "bed malloc len error.\n");
        return 1;
    }

    if (bodloaddata_n(bod_data, bod_readoutall, readout_len, start, end) != 0) {
        fprintf(stderr, "read bod failed.\n");
        return 1;
    }

    return 0;
}

/*
 * Copyright 2023 Nick Conway; Copyright 2018, Nick Conway; Wyss Institute
 * Harvard University
 *
 * Utilities for DNA C-string to integer conversion, CIGAR format IO,
 * and ctypes-friendly SSW alignment wrapper.
 */
#include "str_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int8_t DNA_BASE_LUT[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    /*   A     C            G                          */
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    /*             T                                   */
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    /*   a     c            g                          */
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    /*             t                                   */
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void dna_to_int8(const char *c_str, int8_t *arr, uint32_t len) {
    const char *str_lim = c_str + len;
    while (c_str < str_lim) {
        unsigned char ch = (unsigned char)*c_str++;
        *arr++ = (ch < 128) ? DNA_BASE_LUT[ch] : 4;
    }
}

void ssw_write_cigar(const s_align *a) {
    int32_t c;
    if (a && a->cigar) {
        printf("CIGAR start index %d: ", a->ref_begin1);
        for (c = 0; c < a->cigarLen; ++c) {
            char letter = cigar_int_to_op(a->cigar[c]);
            uint32_t length = cigar_int_to_len(a->cigar[c]);
            printf("%u%c", length, letter);
        }
        printf("\n");
    } else {
        printf("no CIGAR\n");
    }
}

int ssw_align_dna(
    const char *read_seq,
    int32_t read_len,
    const char *ref_seq,
    int32_t ref_len,
    const int8_t *mat,
    int32_t mat_n,
    uint8_t gap_open,
    uint8_t gap_extend,
    uint8_t flag,
    ssw_ctypes_result *out
) {
    int8_t *read_num = NULL;
    int8_t *ref_num = NULL;
    s_profile *prof = NULL;
    s_align *aln = NULL;
    int32_t mask_len;

    if (!read_seq || !ref_seq || !mat || !out || read_len <= 0 || ref_len <= 0 || mat_n <= 0) {
        return -1;
    }

    memset(out, 0, sizeof(*out));
    out->ref_begin1 = -1;
    out->read_begin1 = -1;

    read_num = (int8_t *)malloc((size_t)read_len * sizeof(int8_t));
    ref_num = (int8_t *)malloc((size_t)ref_len * sizeof(int8_t));
    if (!read_num || !ref_num) {
        free(read_num);
        free(ref_num);
        return -2;
    }

    dna_to_int8(read_seq, read_num, (uint32_t)read_len);
    dna_to_int8(ref_seq, ref_num, (uint32_t)ref_len);

    prof = ssw_init(read_num, read_len, mat, mat_n, 2);
    if (!prof) {
        free(read_num);
        free(ref_num);
        return -3;
    }

    mask_len = (read_len / 2) < 15 ? 15 : (read_len / 2);
    aln = ssw_align(prof, ref_num, ref_len, gap_open, gap_extend, flag, 0, 0, mask_len);
    if (!aln) {
        init_destroy(prof);
        free(read_num);
        free(ref_num);
        return -4;
    }

    out->score1 = aln->score1;
    out->score2 = aln->score2;
    out->ref_begin1 = aln->ref_begin1;
    out->ref_end1 = aln->ref_end1;
    out->read_begin1 = aln->read_begin1;
    out->read_end1 = aln->read_end1;
    out->cigar_len = aln->cigarLen;

    if (aln->cigar && aln->cigarLen > 0) {
        size_t nbytes = (size_t)aln->cigarLen * sizeof(uint32_t);
        out->cigar = (uint32_t *)malloc(nbytes);
        if (!out->cigar) {
            align_destroy(aln);
            init_destroy(prof);
            free(read_num);
            free(ref_num);
            return -5;
        }
        memcpy(out->cigar, aln->cigar, nbytes);
    }

    align_destroy(aln);
    init_destroy(prof);
    free(read_num);
    free(ref_num);
    return 0;
}

void ssw_free_ctypes_result(ssw_ctypes_result *out) {
    if (out) {
        free(out->cigar);
        out->cigar = NULL;
        out->cigar_len = 0;
    }
}

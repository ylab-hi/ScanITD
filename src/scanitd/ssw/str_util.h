/*
 * Copyright 2023 Nick Conway; Copyright 2018, Nick Conway; Wyss Institute
 * Harvard University
 *
 * Utilities for DNA C-string to integer conversion, CIGAR format IO,
 * and ctypes-friendly SSW alignment wrapper.
 */
#ifndef SSWPY_STR_UTIL_H
#define SSWPY_STR_UTIL_H

#include <stdint.h>
#include "ssw.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint16_t score1;
    uint16_t score2;
    int32_t ref_begin1;
    int32_t ref_end1;
    int32_t read_begin1;
    int32_t read_end1;
    int32_t cigar_len;
    uint32_t *cigar;
} ssw_ctypes_result;

void dna_to_int8(const char *c_str, int8_t *arr, uint32_t len);
void ssw_write_cigar(const s_align *a);

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
);

void ssw_free_ctypes_result(ssw_ctypes_result *out);

#ifdef __cplusplus
}
#endif

#endif /* SSWPY_STR_UTIL_H */

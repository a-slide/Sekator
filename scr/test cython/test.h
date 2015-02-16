#ifndef TEST_H
#define TEST_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

typedef struct {
    uint16_t score1;
    uint16_t score2;
    int32_t ref_begin1;
    int32_t ref_end1;
    int32_t read_begin1;
    int32_t read_end1;
    int32_t ref_end2;
    uint32_t* cigar;
    int32_t cigarLen;
} align;


int8_t* score_matrix (const uint8_t match, const uint8_t mismatch, const uint8_t ambiguous);
int8_t* DNA_seq_to_int(char * DNA_seq, int seq_len);
int8_t* malloc_int8_t (int dim);
int8_t* calloc_int8_t (int dim);
void print_int8_t (int8_t* mat, int dim);
void free_int8_t (int8_t* mat);

#endif  // TEST_H

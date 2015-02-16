#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

/// structure of the alignment result profile
struct _alignment_end;
typedef struct _alignment_end s_alignment_end;

struct _align;
typedef struct _align s_align;

// Public functions
int8_t* score_matrix (
    const int8_t match,
    const int8_t mismatch,
    const int8_t ambiguous);

int8_t* DNA_seq_to_int(
    char* DNA_seq,
    int seq_len);

s_align ssw_align (
    const int8_t* query,
    int32_t queryLen,
    const int8_t* ref,
    int32_t refLen,
    const int8_t* mat,
    const int8_t gapO,
    const int8_t gapE);

// Private functions
static __m128i* qP_word (
    const int8_t* query_num,
    const int8_t* mat,
    const int32_t queryLen);

static s_alignment_end sw_sse2_word (
    const int8_t* ref,
    int8_t ref_dir,
    int32_t refLen,
    int32_t queryLen,
    const uint8_t gapO,
    const uint8_t gapE,
    const __m128i* vProfile,
    uint16_t terminate);

static int8_t* seq_reverse(
    const int8_t* seq,
    int32_t end);

// Memory management functions
void print_int8_t (int8_t* mat, int dim);
int8_t* malloc_int8_t (int dim);

#endif  // SSW_H

/*
@package    FastqFT
@brief      Simplified Fast C implementation of Smith Waterman based on SSW Library
* An SIMD Smith-Waterman C/C++ Library for Use in Genomic Applications
* from Mengyao Zhao & Wan-Ping Lee available
* Available at https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library

@copyright  The MIT licence
* The MIT License
* Copyright (c) 2012-1015 Boston College.
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:

* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
* BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
* ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
@author Adrien Leger
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
*/

#include <emmintrin.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ssw.h"

//#################################################################################################
// MACRO
//#################################################################################################

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* Convert the coordinate in the scoring matrix into the coordinate in one line of the band. */
#define set_u(u, w, i, j) { int x=(i)-(w); x=x>0?x:0; (u)=(j)-x+1; }

/* Convert the coordinate in the direction matrix into the coordinate in one line of the band. */
#define set_d(u, w, i, j, p) { int x=(i)-(w); x=x>0?x:0; x=(j)-x; (u)=x*3+p; }

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

//#################################################################################################
// STRUCTURES AND GLOBAL VARIABLES
//#################################################################################################

struct _alignment_end {
    uint16_t score;
    int32_t ref;     //0-based position
    int32_t query;    //alignment ending position on query, 0-based
};
typedef struct _alignment_end s_alignment_end;

/* structure of the alignment result profile  */
struct _align{
    uint16_t score;
    int32_t ref_begin;
    int32_t ref_end;
    int32_t query_begin;
    int32_t query_end;
};
typedef struct _align s_align;

/* This table is used to transform nucleotide letters into numbers.
 * A=a=0, C=c=1, G=g=2, T=t=3, all other char = 4*/
static const int8_t NT_TABLE[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//#################################################################################################
// TEST MAIN
//#################################################################################################

int main(void)
{
    int8_t* query;
    int8_t* ref;
    int8_t* score;
    s_align res;

    query = DNA_seq_to_int("ATCGN", 5);
    print_int8_t (query,5);
    ref = DNA_seq_to_int("AATCGTCAGATCNNAGTCTGC", 20);
    print_int8_t (ref,20);
    score = score_matrix (2, 2);
    print_int8_t (score,25);

    res = ssw_align (query, 5, ref, 20, score, 3, 1);
    printf("Score %d\n", res.score);
    printf("Ref begin %d\n", res.ref_begin);
    printf("Ref end %d\n", res.ref_end);
    printf("Query begin %d\n", res.query_begin);
    printf("Query end %d\n", res.query_end);

    return (1);
}

//#################################################################################################
// FUNCTIONS
//#################################################################################################

int8_t* score_matrix (
    const int8_t match, // Bonus for a match (POSITIVE)
    const int8_t mismatch) // Malus for a mismatch (POSITIVE)
{
    int32_t i, j, k;
    int8_t* mat;
    mat = malloc_int8_t(25);

    // a score of 0 is automatically applied for ambiguous bases
    for (i = k = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) mat[k++] = i == j ? match : - mismatch;
        mat[k++] = 0;
    }
    for (i = 0; i < 5; ++i) mat[k++] = 0;

    return mat;
}

int8_t* DNA_seq_to_int(char * DNA_seq, const int32_t seq_len)
{
    int32_t i;
    int8_t* mat;
    mat = malloc_int8_t(seq_len);

    for (i = 0; i < seq_len; ++i) mat[i] = NT_TABLE[(int)DNA_seq[i]];

    return mat;
}

s_align ssw_align (const int8_t* query, // Query sequence encoded by integers
                    int32_t queryLen, // Length of the query sequence
                    const int8_t* ref, // Reference sequence encoded by integers
                    int32_t refLen, // Length of the reference sequence
                    const int8_t* mat, // Score matrix produced by score_matrix()
                    const int8_t gapO, // Weight of gap open (POSITIVE)
                    const int8_t gapE) // Weight of gap extend (POSITIVE)
{
    s_alignment_end best;
    s_alignment_end best_reverse;
    __m128i* vProfile_for;
    __m128i* vProfile_rev;
    int8_t* query_reverse;
    s_align res;

    // Find the alignment scores and ending positions
    vProfile_for = qP_word(query, mat, queryLen);
    best = sw_sse2_word(ref, 0, refLen, queryLen, gapO, gapE, vProfile_for, -1);
    free(vProfile_for);
    res.score = best.score;
    res.ref_end = best.ref;
    res.query_end = best.query;

    // Find the beginning position of the best alignment.
    query_reverse = seq_reverse(query, res.query_end);
    vProfile_rev = qP_word(query_reverse, mat, res.query_end+1);
    best_reverse = sw_sse2_word(ref, 1, res.ref_end+1, res.query_end+1, gapO, gapE, vProfile_rev, res.score);
    free(vProfile_rev);
    free(query_reverse);
    res.ref_begin = best_reverse.ref;
    res.query_begin = res.query_end - best_reverse.query;

    return res;
}

static __m128i* qP_word (const int8_t* query_num,
                  const int8_t* mat,
                  const int32_t queryLen)
{
    int32_t segLen = (queryLen + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(5 * segLen * sizeof(__m128i));
    int16_t* t = (int16_t*)vProfile;
    int32_t nt, i, j;
    int32_t segNum;

    /* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
    for (nt = 0; LIKELY(nt < 5); nt ++) {
        for (i = 0; i < segLen; i ++) {
            j = i;
            for (segNum = 0; LIKELY(segNum < 8) ; segNum ++) {
                *t++ = j>= queryLen ? 0 : mat[nt * 5 + query_num[j]];
                j += segLen;
            }
        }
    }
    return vProfile;
}

static s_alignment_end sw_sse2_word (const int8_t* ref,
                             int8_t ref_dir,    // 0: forward ref; 1: reverse ref
                             int32_t refLen,
                             int32_t queryLen,
                             const uint8_t gapO, /* positive but will be used as - */
                             const uint8_t gapE, /* positive but will be used as - */
                             const __m128i* vProfile,
                             uint16_t terminate)
{

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
                    (m) = _mm_extract_epi16((vm), 0)

    uint16_t max = 0;                            /* the max alignment score */
    int32_t end_query = queryLen - 1;
    int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
    int32_t segLen = (queryLen + 7) / 8; /* number of segment */

    /* array to record the largest score of each reference position */
    uint16_t* maxColumn = (uint16_t*) calloc(refLen, 2);

    /* array to record the alignment query ending position of the largest score of each reference position */
    int32_t* end_query_column = (int32_t*) calloc(refLen, sizeof(int32_t));

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_set1_epi32(0);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

    int32_t i, j, k;
    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(gapO);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gapE);

    __m128i vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxMark = vZero; /* Trace the highest score till the previous column. */
    __m128i vTemp;
    int32_t begin = 0, end = refLen, step = 1;

    /* outer loop to process the reference sequence */
    if (ref_dir == 1) {
        begin = refLen - 1;
        end = -1;
        step = -1;
    }
    for (i = begin; LIKELY(i != end); i += step) {
        int32_t cmp;
        __m128i e, vF = vZero; /* Initialize F value to 0.
                               Any errors to vH values will be corrected in the Lazy_F loop.
                             */
        __m128i vH = pvHStore[segLen - 1];
        vH = _mm_slli_si128 (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;

        __m128i vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */

        const __m128i* vP = vProfile + ref[i] * segLen; /* Right part of the vProfile */
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (j = 0; LIKELY(j < segLen); j ++) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + j));

            /* Get max from vH, vE and vF. */
            e = _mm_load_si128(pvE + j);
            vH = _mm_max_epi16(vH, e);
            vH = _mm_max_epi16(vH, vF);
            vMaxColumn = _mm_max_epi16(vMaxColumn, vH);

            /* Save vH values. */
            _mm_store_si128(pvHStore + j, vH);

            /* Update vE value. */
            vH = _mm_subs_epu16(vH, vGapO); /* saturation arithmetic, result >= 0 */
            e = _mm_subs_epu16(e, vGapE);
            e = _mm_max_epi16(e, vH);
            _mm_store_si128(pvE + j, e);

            /* Update vF value. */
            vF = _mm_subs_epu16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + j);
        }

        /* Lazy_F loop: has been revised to disallow adjacent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        for (k = 0; LIKELY(k < 8); ++k) {
            vF = _mm_slli_si128 (vF, 2);
            for (j = 0; LIKELY(j < segLen); ++j) {
                vH = _mm_load_si128(pvHStore + j);
                vH = _mm_max_epi16(vH, vF);
                _mm_store_si128(pvHStore + j, vH);
                vH = _mm_subs_epu16(vH, vGapO);
                vF = _mm_subs_epu16(vF, vGapE);
                if (UNLIKELY(! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH)))) goto end;
            }
        }

end:
        vMaxScore = _mm_max_epi16(vMaxScore, vMaxColumn);
        vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
        cmp = _mm_movemask_epi8(vTemp);
        if (cmp != 0xffff) {
            uint16_t temp;
            vMaxMark = vMaxScore;
            max8(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (LIKELY(temp > max)) {
                max = temp;
                end_ref = i;
                for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
            }
        }

        /* Record the max score of current column. */
        max8(maxColumn[i], vMaxColumn);
        if (maxColumn[i] == terminate) break;
    }

    /* Trace the alignment ending position on query. */
    uint16_t *t = (uint16_t*)pvHmax;
    int32_t column_len = segLen * 8;
    for (i = 0; LIKELY(i < column_len); ++i, ++t) {
        int32_t temp;
        if (*t == max) {
            temp = i / 8 + i % 8 * segLen;
            if (temp < end_query) end_query = temp;
        }
    }

    free(pvHmax);
    free(pvE);
    free(pvHLoad);
    free(pvHStore);

    /* Find the best alignment. */
    s_alignment_end best;
    best.score = max;
    best.ref = end_ref;
    best.query = end_query;

    free(maxColumn);
    free(end_query_column);

    return best;
}

static int8_t* seq_reverse(const int8_t* seq, int32_t end)  /* end is 0-based alignment ending position */
{
    int8_t* reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));
    int32_t start = 0;
    while (LIKELY(start <= end)) {
        reverse[start] = seq[end];
        reverse[end] = seq[start];
        ++ start;
        -- end;
    }
    return reverse;
}


int8_t* malloc_int8_t (int dim)
{
    int8_t* mat;

    mat = malloc (sizeof(int8_t) * dim);
    if (mat == NULL)
    {
        fprintf (stderr, "Memory allocation error\n\n");
        exit (EXIT_FAILURE);
    }

    return mat;
}

void print_int8_t (int8_t* mat, int dim)
{
    int i;

    for (i = 0 ; i < dim ; i++)
        printf("%d ", mat [i]);

    printf("\n");
    return;
}

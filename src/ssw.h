/*
@package    FastqFT
@brief      Simplified Fast C implementation of Striped Smith Waterman based on SSW Library
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

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

//// Structures ////

struct _alignment_end;
typedef struct _alignment_end s_alignment_end;
/*! @typedef    structure of the alignment result profile*/

struct _align;
typedef struct _align s_align;
/*! @typedef    structure of the alignment result
    @field  score       the best alignment score
    @field  ref_begin   0-based best alignment beginning position on reference; ref_begin = -1 when the best alignment beginning
                        position is not available
    @field  ref_end     0-based best alignment ending position on reference
    @field  read_begin  0-based best alignment beginning position on read; read_begin = -1 when the best alignment beginning
                        position is not available
    @field  read_end    0-based best alignment ending position on read
*/

//// Public functions ////

int8_t* score_matrix (const int8_t match, const int8_t mismatch);

/*! @function   Generate a 1D integer array representing a 5*5 score matrix from the provided match, mismatch and ambiguous penalties
    @param match        Bonus in case of match (POSITIVE)
    @param mismatch     Penalty in case of mismatch (POSITIVE)
    @return 8 bits integer array encoding a score matrix (see bellow)
    @note   If the penalty for match and mismatch are 2, 2, respectively, the score matrix will be:
             A  C  G  T  N  (or other ambiguous code)
             2 -2 -2 -2  0   A
            -2  2 -2 -2  0   C
            -2 -2  2 -2  0   G
            -2 -2 -2  2  0   T
             0  0  0  0  0   N (or other ambiguous code)
            The corresponding array is {2,-2,-2,-2,0,-2,2,-2,-2,0,-2,-2,2,-2,0,-2,-2,-2,2,0,0,0,0,0,0}
*/

int8_t* DNA_seq_to_int(char* DNA_seq, const int32_t seq_len);

/*! @function       Transform a DNA sequence in a int8_t* integer array with A,a=0; C,c=1; G,g=2; T,t=3; other char = 4
    @param  seq_len Length of the DNA sequence
    @param  DNA_seq DNA sequence in char
    @return 8 bits integer array encoding the sequence
    @note   Example if the query sequence is: ACGTN, the array will be : {0, 1, 2, 3, 4}
*/

s_align ssw_align (const int8_t* query, int32_t queryLen, const int8_t* ref, int32_t refLen, const int8_t* mat, const int8_t gapO, const int8_t gapE);

/*! @function   Do Striped Smith-Waterman alignment.
    @param query    Pointer to the query sequence; the query sequence needs to be numbers
    @param queryLen Length of the query sequence
    @param ref      Pointer to the reference sequence; the reference sequence needs to be numbers
    @param refLen   Length of the reference sequence
    @param mat      Score matrix produced by score_matrix
    @param gapO     Penalty in case of gap opening (POSITIVE)
    @param gapE     Penalty in case of gap extension (POSITIVE)
    @return Alignment result structure
*/

//// Private functions ////

static __m128i* qP_word (const int8_t* query_num, const int8_t* mat, const int32_t queryLen);
// Generate query profile rearrange query sequence & calculate the weight of match/mismatch

static s_alignment_end sw_sse2_word (const int8_t* ref, int8_t ref_dir, int32_t refLen,
    int32_t queryLen, const uint8_t gapO, const uint8_t gapE, const __m128i* vProfile,
    uint16_t terminate);
/* Striped Smith-Waterman
   Record the highest score of each reference position. Return the alignment score and ending position of the best alignment.
   Gap begin and gap extension are different. with_match > 0, all other weights < 0. The returned positions are 0-based.
 */

static int8_t* seq_reverse(const int8_t* seq, int32_t end);
// Reverse the order of "letters" in a int8_t* sequence

void print_int8_t (int8_t* mat, int dim);
// Print an int8_t* mat on one line

int8_t* malloc_int8_t (int dim);
// Allocate a 1D int8_t mat and verify its allocation

#endif  // SSW_H

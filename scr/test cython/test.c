#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "test.h"

int8_t* score_matrix (const uint8_t match, const uint8_t mismatch, const uint8_t ambiguous)
{
    int32_t i, j, k;
    int8_t* mat;
    mat = malloc_int8_t(25);

    for (i = k = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) mat[k++] = i == j ? match : - mismatch;
        mat[k++] = ambiguous;
    }
    for (i = 0; i < 5; ++i) mat[k++] = ambiguous;

    return mat;
}

int8_t* DNA_seq_to_int (char * DNA_seq, int seq_len)
{
    int32_t i;
    int8_t* mat;
    mat = malloc_int8_t(seq_len);

    for (i = 0; i < seq_len; ++i) {

        switch (DNA_seq[i])
        {
          case 'A':
            mat[i]=0;
            break;
          case 'C':
            mat[i]=1;
            break;
          case 'G':
            mat[i]=2;
            break;
          case 'T':
            mat[i]=3;
            break;
          case 'a':
            mat[i]=0;
            break;
          case 'c':
            mat[i]=1;
            break;
          case 'g':
            mat[i]=2;
            break;
          case 't':
            mat[i]=3;
            break;
          default:
            mat[i]=4;
            break;
        }
    }

    return mat;
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

    printf("\n");
    for (i = 0 ; i < dim ; i++)
        printf("%d ", mat [i]);

    printf("\n");
    return;
}

void free_int8_t (int8_t* mat)
{
    free (mat);
    return;
}

int main (int argc, char * const argv[]) {

    int8_t* mat;
    int8_t* seq;

    mat = score_matrix (1, 2, 3);
    print_int8_t (mat, 25);

    seq = DNA_seq_to_int("ACGTNX", 6);
    print_int8_t (seq, 6);

    free_int8_t (mat);
    free_int8_t (seq);

    return(0);
}

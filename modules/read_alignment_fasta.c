/*=============================================================================
 *
 * read_alignment_fasta.c reads a fasta alignment file
 *
 * Input: file name, max_gap_fraction
 *
 * Output: length, sequences, alignment in numeric format
 *
 * This is a MEX-file for MATLAB.
 *
 * This code accompanies the paper "Fast and accurate multivariate Gaussian
 * modeling of protein families: Predicting residue contacts and
 * protein-interaction partners" by Carlo Baldassi, Marco Zamparo, Christoph
 * Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weight and Andrea
 * Pagnani, submitted to PLOS ONE (2013).
 *
 * The code is released under the terms of the GNU General Public License
 * version 3 or, at your option, any later version.
 *
 *============================================================================*/

#include <mex.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* full AA alphabet */
const double l2n_alphabet[] = { 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20};
                             /* A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y */

double letter2number(char l)
{
    int i = l - 0x41; /* 0x41 == 'A' */
    if (i >= 0 && i < 25) {
        return l2n_alphabet[i];
    } else {
        return 21;
    }
}

void convert_seq(char * s, int * inds, double ** Z, int i)
{
    int j, k, l;
    l = strlen(s);
    k = 0;
    for (j = 0; j < l; ++j) if (inds[j]) {
        Z[k][i] = letter2number(s[j]);
        ++k;
    }
}

/* Parser */
void parse_seq_pass1(kseq_t * seq, int ** inds_ptr, int ** zinds_ptr, int * N_ptr, int * M_ptr, double max_gap_fraction)
{
    int t;
    int l0 = 0;
    int N = 0;
    int M = 0;
    int * inds = NULL;
    int * zinds = NULL;
    int zl = 1000;
    int c;
    int ngaps;
    int seqn = 0;

    zinds = malloc(zl * sizeof(int));

    while ((t = kseq_read(seq)) >= 0) {
        char * s = seq->seq.s;
        int l = strlen(s);
        ngaps = 0;
        if (M == 0) {
            inds = malloc(l * sizeof(int));
            for (c = 0; c < l; ++c) {
                inds[c] = (s[c] != '.' && s[c] == toupper(s[c]));
                ngaps += (s[c] == '-');
            }
            l0 = l;
            for (c = 0; c < l0; ++c) {
                N += inds[c];
            }
        } else {
            if (l != l0) {
                mexErrMsgIdAndTxt("read_alignment_fasta:input", "input data is unaligned");
            }
            for (c = 0; c < l; ++c) {
                if (inds[c] != (s[c] != '.' && s[c] == toupper(s[c]))) {
                    mexErrMsgIdAndTxt("read_alignment_fasta:input", "input data is unaligned?");
                }
                ngaps += (s[c] == '-');
            }
        }
        if (seqn > zl) {
            zl *= 2;
            zinds = realloc(zinds, zl * sizeof(int));
        }
        if ((double) ngaps / N <= max_gap_fraction) {
            zinds[seqn] = 1;
            ++M;
        } else {
            zinds[seqn] = 0;
        }
        ++seqn;
    }
    *inds_ptr = inds;
    *zinds_ptr = zinds;
    *N_ptr = N;
    *M_ptr = M;
    return;
}

/* Parser */
void parse_seq_pass2(kseq_t * seq, double ** Z, int * inds, int * zinds)
{
    int l;
    int i = 0;
    int j = 0;
    while ((l = kseq_read(seq)) >= 0) {
        if (!zinds[j++]) {
            continue;
        }
        convert_seq(seq->seq.s, inds, Z, i++);
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char * filename;
    double max_gap_fraction;
    gzFile fp;
    kseq_t *seq;
    int N, M;
    double * N_ptr;
    double * M_ptr;
    double * Z_ptr;
    double ** Z;
    int * inds;
    int * zinds;

    /* check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("read_alignemnt_fasta:nrhs", "Two inputs required: filename, max_gap_fraction.");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("read_alignemnt_fasta:nlhs", "Three outputs required: N, M, Z.");
    }

    /* get the value of the frequence matrices  */
    filename = mxArrayToString(prhs[0]);

    fp = gzopen(filename, "r");

    if (fp == Z_NULL) {
        mexErrMsgIdAndTxt("read_alignemnt_fasta:open_file", "Error opening file");
    }

    seq = kseq_init(fp);

    /* get the max_gap_fraction value */
    max_gap_fraction = mxGetScalar(prhs[1]);

    /* create the outputs N, M */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

    N_ptr = mxGetPr(plhs[0]);
    M_ptr = mxGetPr(plhs[1]);

    parse_seq_pass1(seq, &inds, &zinds, &N, &M, max_gap_fraction);

    *N_ptr = (double) N;
    *M_ptr = (double) M;

    /* create the output matrix Z */
    plhs[2] = mxCreateDoubleMatrix(M, N, mxREAL);

    Z_ptr = mxGetPr(plhs[2]);
    Z = malloc(N * sizeof(double));
    {
        int i;
        for (i = 0; i < N; ++i) {
            Z[i] = Z_ptr;
            Z_ptr += M;
        }
    }

    gzrewind(fp);
    kseq_rewind(seq);

    parse_seq_pass2(seq, Z, inds, zinds);

    /* release memory */
    kseq_destroy(seq);
    gzclose(fp);

    mxFree(filename);

    free(Z);

    free(inds);
}

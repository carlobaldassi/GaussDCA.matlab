/*=============================================================================
 *
 * compute_true_frequencies.c calculates sequence weights for MSA
 *
 * Input: MSA in numerical form, matrix dimensions, q, reweighting flag
 *
 * Output: frequencies matrices, weight vector
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

#if !(PACKBITS == 128 || PACKBITS == 64 || PACKBITS == 32)
#  warning "Invalid PACKBITS value; it must be 32, 64 or 128: using default value (64)"
#  undef PACKBITS
#  define PACKBITS 64
#endif

#if PACKBITS == 128
#  define packed __uint128_t
const packed _u = ((packed) 0x2108421084210842 << 64) | ((packed) 0x1084210842108421);
const packed _alt = ((packed) 0x1f07c1f07c1f07c1 << 64) | ((packed) 0xf07c1f07c1f07c1f);
#  define collapse10(t) \
    t += t >> 10; \
    t += t >> 20; \
    t += t >> 40; \
    t += t >> 80; \
    t &= 0x3ff;
#elif PACKBITS == 64
#  define packed unsigned long long
const packed _u = 0x0084210842108421;
const packed _alt = 0x007c1f07c1f07c1f;
#  define collapse10(t) \
    t += t >> 10; \
    t += t >> 20; \
    t += t >> 40; \
    t &= 0x3ff
#else /* PACKBITS == 32 */
#  define packed unsigned int
const packed _u = 0x02108421;
const packed _alt = 0x01f07c1f;
#  define collapse10(t) \
    t += t >> 10; \
    t += t >> 20; \
    t &= 0x3ff
#endif

#define packfactor (PACKBITS / 5)
#define packrest (PACKBITS % 5)

const packed _msk = ~((packed) 0);

#define clength(l) (((l)-1)/packfactor + 1)
#define crest(l)   (((l)-1)%packfactor + 1)

#define nnz_aux(x) \
    (((x) | ((x) >> 1) | ((x) >> 2) | ((x) >> 3) | ((x) >> 4)) & _u)

#define nz_aux(nx) \
    ((nx) & ((nx) >> 1) & ((nx) >> 2) & ((nx) >> 3) & ((nx) >> 4) & _u)

#define nz_aux2(nx, s) nz_aux(nx) & (_msk >> (s))

#define collapse(z, t) \
    t = (z & _alt) + ((z >> 5) & _alt); \
    collapse10(t);

packed ** compressedalign(double **align, int m, int length)
{
    packed * calignmem;
    packed ** calign;
    int i;
    int cl;

    if (sizeof(packed) * CHAR_BIT < PACKBITS) {
        fprintf(stderr, "ERROR unsupported system\n");
        exit(1);
    }

    cl = clength(length);

    calignmem = calloc(m * cl, sizeof(packed));
    calign = malloc(m * sizeof(packed*));
    for (i = 0; i < m; ++i) {
        calign[i] = calignmem;
        calignmem += cl;
    }

#ifdef OMP
#  pragma omp parallel for
#endif
    for (i = 0; i < m; ++i) {
        int k;
        for (k = 0; k < length; ++k) {
            int k0 = k / packfactor;
            int k1 = k % packfactor;
            calign[i][k0] |= ((packed)align[i][k] << 5*k1);
        }
    }

    return calign;
}

double compute_theta(packed **calign, int m, int length)
{
    int i;
    double meanfracid;
    double theta = 0.;

    const int cl = clength(length);
    const int cr = 5 * (packfactor - crest(length)) + packrest;

    const int kmax = (cl - 1) / 31;
    const int rmax = (cl - 1) % 31;

    unsigned long long nids = 0;

    /* count seqs below theta dist */
#ifdef OMP
#  pragma omp parallel for schedule(dynamic) reduction (+ : nids)
#endif
    for (i = 0; i < m - 1; ++i) {
        int j;
        for (j = i + 1; j < m; ++j) {
            int k, r;
            packed * caligni = calign[i];
            packed * calignj = calign[j];
            packed z, ny, t;
            for (k = 0; k < kmax; ++k) {
                z = 0;
                for (r = 0; r < 31; ++r, ++caligni, ++calignj) {
                    ny = ~(*caligni) ^ *calignj;
                    z += nz_aux(ny);
                }
                collapse(z, t);
                nids += (int) t;
            }

            z = 0;
            for (r = 0; r < rmax; ++r, ++caligni, ++calignj) {
                ny = ~(*caligni) ^ *calignj;
                z += nz_aux(ny);
            }
            ny = ~(*caligni) ^ *calignj;
            z += nz_aux2(ny, cr);
            collapse(z, t);
            nids += (int) t;
        }
    }
    meanfracid = (2.0 * nids) / (m * (m - 1.0) * length);
    theta = 0.38 * 0.32 / meanfracid;

    mexPrintf("theta = %g\n", theta);
    return theta;
}

/* Weight vector computation */
double compute_weights(packed ** calign, double theta, int m, int length, double *W)
{
    int i;
    const double thresh = (int)(theta * length);
    const int cl = clength(length);
    const int kmax = (cl - 1) / 31;
    const int rmax = (cl - 1) % 31 + 1;
    double Meff = 0.;

    mexPrintf("length = %i\n", length);
    mexPrintf("threshold = %g\n", thresh);

    /* initialize W to one */

#ifdef OMP
#  pragma omp parallel for
#endif
    for (i = 0; i < m; ++i) {
        W[i] = 1.;
    }

    if (theta == 0.) {
        mexPrintf("M = %i N = %i Meff = %g\n", m, length, (double) m);
        return (double) m;
    }

    /* count seqs below theta dist */
#ifdef OMP
#  pragma omp parallel for schedule(dynamic)
#endif
    for (i = 0; i < m - 1; ++i) {
        int j;
        for (j = i + 1; j < m; ++j) {
            int k, r;
            int dist = 0;
            packed * caligni = calign[i];
            packed * calignj = calign[j];
            packed z, y, t;
            for (k = 0; k < kmax && dist < thresh; ++k) {
                z = 0;
                for (r = 0; r < 31; ++r, ++caligni, ++calignj) {
                    y = *caligni ^ *calignj;
                    z += nnz_aux(y);
                }
                collapse(z, t);
                dist += (int) t;
            }
            if (dist < thresh) {
                z = 0;
                for (r = 0; r < rmax; ++r, ++caligni, ++calignj) {
                    y = *caligni ^ *calignj;
                    z += nnz_aux(y);
                }
                collapse(z, t);
                dist += (int) t;
            }

            if (dist < thresh)
#ifdef OMP
#  pragma omp critical
#endif
            {
                ++W[i];
                ++W[j];
            }
        }
    }

    /* weight = 1 / number of sequences */

#ifdef OMP
#  pragma omp parallel for reduction (+ : Meff)
#endif
    for (i = 0; i < m; ++i) {
        W[i] = 1.0 / W[i];
        Meff += W[i];
    }

    mexPrintf("M = %i N = %i Meff = %g\n", m, length, Meff);
    return Meff;
}

/* Frequencies computation */
void compute_frequencies(double **** Pij, double ** Pi, double ** align,
        int M, int N, int s, double * W, double Meff)
{
    int i, j;

#ifdef OMP
#  pragma omp parallel for
#endif
    for (i = 0; i < N; ++i) {
        int j, a, b;
        for (a = 0; a < s; ++a) {
            for (j = 0; j < N; ++j) {
                for (b = 0; b < s; ++b) {
                    Pij[i][a][j][b] = 0.;
                }
            }
        }
    }
#ifdef OMP
#  pragma omp parallel for
#endif
    for (j = 0; j < N; ++j) {
        int a;
        for (a = 0; a < s; ++a) {
            Pi[j][a] = 0.;
        }
    }

#ifdef OMP
#  pragma omp parallel for
#endif
    for (j = 0; j < N; ++j) {
        int i, a;
        for (i = 0; i < M; ++i) {
            int aj = align[i][j] - 1;
            double Wi = W[i];
            if (aj == s) continue;
            Pi[j][aj] += Wi;
        }
        for (a = 0; a < s; ++a) {
            Pi[j][a] /= Meff;
        }
    }

#ifdef OMP
#  pragma omp parallel for schedule(dynamic)
#endif
    for (j = 0; j < N; ++j) {
        int i, k, a, b;
        for (i = 0; i < M; ++i) {
            double Wi = W[i];
            a = align[i][j] - 1;
            if (a == s) continue;
            for (k = j + 1; k < N; ++k) {
                b = align[i][k] - 1;
                if (b == s) continue;
                Pij[j][a][k][b] += Wi;
            }
        }
        for (a = 0; a < s; ++a) {
            for (k = j + 1; k < N; ++k) {
                for (b = 0; b < s; ++b) {
                    Pij[j][a][k][b] /= Meff;
                    Pij[k][b][j][a] = Pij[j][a][k][b];
                }
            }
        }
    }

#ifdef OMP
#  pragma omp parallel for
#endif
    for (j = 0; j < N; ++j) {
        int a;
        for (a = 0; a < s; ++a) {
            Pij[j][a][j][a] = Pi[j][a];
        }
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *align_ptr;   /* alignment pointer */
    double **align;      /* alignment matrix */
    int m;               /* number of sequences */
    int length;          /* number of AAs = length of seq */
    int q;               /* number of AAs in the alphabet */
    double ****Pij;      /* output matrix Pij */
    double **Pi;         /* output matrix Pi */
    double *W;           /* output matrix W */
    double theta;        /* scalar threshold */
    packed **calign = NULL;
    double *Pi_ptr;
    double *Pij_ptr;
    int s;
    int sN;
    double Meff;

    /* check for proper number of arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("Compute_True_Frequencies:nrhs", "Five inputs required: align, M, N, q, theta.");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("Compute_True_Frequencies:nlhs", "Three outputs required: Pij, Pi, W.");
    }

    /* get the value of the alignment  */
    align_ptr = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    m = (int) (mxGetScalar(prhs[1]) + 1e-12);
    length = (int) (mxGetScalar(prhs[2]) + 1e-12);
    q = (int) (mxGetScalar(prhs[3]) + 1e-12);
    s = q - 1;

    if (q > 32) {
        mexErrMsgIdAndTxt("Compute_True_Frequencies:q", "Parameter q is too big");
    }

    /* get theta (a value of -1 means 'auto') */
    theta = mxGetScalar(prhs[4]);

    /* set up the alignment matrix */
    align = malloc(m * sizeof(double*));
    {
        int i;
        for (i = 0; i < m; ++i) {
            align[i] = align_ptr;
            align_ptr += length;
        }
    }

    calign = compressedalign(align, m, length);

    if (theta == -1) {
        theta = compute_theta(calign, m, length);
    }

    /* create the output matrix W */
    plhs[2] = mxCreateDoubleMatrix(1, m, mxREAL);

    /* get a pointer to the real data in W */
    W = mxGetPr(plhs[2]);

    /* compute the weights */
    Meff = compute_weights(calign, theta, m, length, W);

    /* Pij matrix size */
    sN = length * s;

    /* create the output matrices Pij, Pi */
    plhs[0] = mxCreateDoubleMatrix(sN, sN, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(sN, 1, mxREAL);

    /* get pointers to the real data in Pij, Pi */
    Pij_ptr = mxGetPr(plhs[0]);
    Pi_ptr = mxGetPr(plhs[1]);

    /* set up the Pij, Pi matrices */
    Pij = malloc(length * sizeof(double***));
    Pi = malloc(length * sizeof(double*));
    {
        int i, j, a;
        for (i = 0; i < length; ++i) {
            Pij[i] = malloc(s * sizeof(double**));
            for (a = 0; a < s; ++a) {
                Pij[i][a] = malloc(length * sizeof(double*));
                for (j = 0; j < length; ++j) {
                    Pij[i][a][j] = Pij_ptr;
                    Pij_ptr += s;
                }
            }
        }
        for (i = 0; i < length; ++i) {
            Pi[i] = Pi_ptr;
            Pi_ptr += s;
        }
    }

    compute_frequencies(Pij, Pi, align, m, length, s, W, Meff);

    free(align);

    free(calign);

    {
        int i, a;
        for (i = 0; i < length; ++i) {
            for (a = 0; a < s; ++a) {
                free(Pij[i][a]);
            }
            free(Pij[i]);
        }
        free(Pij);
    }
    free(Pi);
}

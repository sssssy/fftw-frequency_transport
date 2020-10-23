#include <stdlib.h>

#include <fftw3.h>

inline double amp(fftw_complex c)
{
    return sqrt(c[0] * c[0] + c[1] * c[1]);
}

inline double* ijo(fftw_complex* c, int N, int i, int j, int offset)
{
    i = (i + offset) % N;
    j = (j + offset) % N;
    return c[j + N*i];
}

inline double* ijklo(fftw_complex* c, int N, int i, int j, int k, int l, int offset)
{
    i = (i + offset) % N;
    j = (j + offset) % N;
    k = (k + offset) % N;
    l = (l + offset) % N;
    return c[l + N*(k + N*(j + N*i))];
}

inline double* ijklo(double* c, int N, int i, int j, int k, int l, int offset)
{
    i = (i + offset) % N;
    j = (j + offset) % N;
    k = (k + offset) % N;
    l = (l + offset) % N;
    return &c[l + N*(k + N*(j + N*i))];
}

void real2complex(double* dpMatrix, fftw_complex *cpIn, int N)
{
    double *dp = dpMatrix;
    fftw_complex *cp;
    fftw_complex *cpEnd = cpIn + (int)pow(1.0*N, 4);

    for (cp = cpIn; cp < cpEnd; cp++, dp++)
    {
        *cp[0] = *dp;
        *cp[1] = 0.0;
    }
}

void fftscaling(fftw_complex* c, int N, int D)
{
    int fftSize = pow(1.0*N, D);
    double scale = fftSize / 2.0;
    //printf("scale = %lf\n", scale);
    fftw_complex* cp, *cpEnd = c + fftSize;
    /*for (cp = c; cp < cpEnd; cp++)
    {
        *cp[0] /= scale;
        *cp[1] /= scale;
    }*/
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++)
                {
                    ijklo(c, N, i, j, k, l, 0)[0] /= scale;
                    ijklo(c, N, i, j, k, l, 0)[1] /= scale;
                    //c[l + N*(k + N*(j + N*i))][0] /= scale;
                    //c[l + N*(k + N*(j + N*i))][1] /= scale;
                }
}

void fftshift2d(fftw_complex* in, fftw_complex* out, int N)
{
    int offset = floor(N / 2.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            ijo(out, N, i, j, offset)[0] = ijo(in, N, i, j, 0)[0];
            ijo(out, N, i, j, offset)[1] = ijo(in, N, i, j, 0)[1];
            //out[l + N*(k + N*(j + N*i))][0] = in[l + N*(k + N*(j + N*i))][0];
            //out[l + N*(k + N*(j + N*i))][1] = in[l + N*(k + N*(j + N*i))][1];
        }
}

void fftshift4d(fftw_complex* in, fftw_complex* out, int N)
{
    int offset = floor(N / 2.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++)
                {
                    ijklo(out, N, i, j, k, l, 0)[0] = ijklo(in, N, i+offset, j, k+offset, l, 0)[0];
                    ijklo(out, N, i, j, k, l, 0)[1] = ijklo(in, N, i+offset, j, k+offset, l, 0)[1];
                    //out[l + N*(k + N*(j + N*i))][0] = in[l + N*(k + N*(j + N*(i)))][0];
                    //out[l + N*(k + N*(j + N*i))][1] = in[l + N*(k + N*(j + N*(i)))][1];
                }
}

void shift4d(double* in, double *out, int N)
{
    int offset = floor(N / 2.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++)
                {
                    *ijklo(out, N, i, j, k, l, 0) = *ijklo(in, N, i, j, k, l, 5);
                }
}
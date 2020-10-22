#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <windows.h>

#include <fftw3.h>

#include "bmpFile.h"
#include "utils.h"

using namespace std;
typedef unsigned char BYTE;

const double PI = 3.1415926;
const int UNDEFINED = 0;

class Spectrum
{
public:
    Spectrum(int iSamplingRate, int iDims = 4, string sName = "unnamed")
        : iTime(0), iDims(iDims), iSamplingRate(iSamplingRate), iRadius(iSamplingRate/2), 
        sName(sName), sLastOp("NULL"), bFakeBilinear(true)
    {
        dpMatrix = (double*)malloc(
            sizeof(double)*(int)pow(1.0*iSamplingRate, iDims));
        cpFmatrix = (fftw_complex*)fftw_malloc(
            sizeof(fftw_complex)*(int)pow(1.0*iSamplingRate, iDims));
    }
    ~Spectrum(){}

    int iTime;
    int iDims;
    int iSamplingRate;
    int iRadius;
    string sName;
    string sLastOp;
    bool bFakeBilinear;

    /* matrice */
    // x * y * theta * phi
    double *dpMatrix; 
    fftw_complex *cpFmatrix;

    void setZero()
    { 
        double *dp, *dpEnd = dpMatrix + (int)pow(float(iSamplingRate), 4);
        for (dp = dpMatrix; dp < dpEnd; dp++)
        {
            *dp = 1.0;
        }
        sLastOp = "setZero";
        printf("setZero.\n");
    }

    void setRect()
    {
        for (int i = 0; i < iSamplingRate; i++)
            for (int j = 0; j < iSamplingRate; j++)
                for (int k = 0; k < iSamplingRate; k++)
                    for (int l = 0; l < iSamplingRate; l++)
                    {
                        if (i >= 0.4*(iSamplingRate-1) && i <= 0.6*(iSamplingRate-1)
                            && k >= 0.3*(iSamplingRate-1) && k <= 0.7*(iSamplingRate-1))
                        {
                            *ijklo(dpMatrix, iSamplingRate, i, j, k, l, 0) = 1.0;
                        }
                        else
                            *ijklo(dpMatrix, iSamplingRate, i, j, k, l, 0) = 0.0;
                    }
        sLastOp = "setRect";
        printf("setRect.\n");
    }

    void fourier()
    {
        int N = 5;
        fftw_complex *in, *out, *shift;
        fftw_plan p1d, p2d;

        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
        shift = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
        p1d = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        p2d = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        /* init in array */
        // 必须在创建plan之后再初始化数据，因为FFTW_MEASURE会覆盖你的数据。
        srand((unsigned)time(NULL));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                in[i*N + j][0] = i * 5 + j + rand() / (double)(RAND_MAX + 1);
                in[i*N + j][1] = 0;
            }

        /* show array */
        printf("Primary: \n\n");
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                printf("%.3f+%.3fi,\t", in[i*N + j][0], in[i*N + j][1]);
            }
            printf("\n");
        }

        fftw_execute(p2d);
        fftshift2d(out, shift, N);
        printf("\nFourier: \n\n");

        /* show spectrum array */
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                printf("|%.3f|,\t", amp(shift[i*N + j]));
            }
            printf("\n");
        }

        fftw_destroy_plan(p1d);
        fftw_destroy_plan(p2d);
        fftw_free(in);
        fftw_free(out);
    }

    void print(double* m)
    {
        double* dp;
        for (int i = 0; i < iSamplingRate; i++)
        {
            for (int k = 0; k < iSamplingRate; k++)
            {
                printf("%.3f, ", *ijklo(m, iSamplingRate, i, 0, k, 0, 0));
            }
            printf("\n");
        }
    }

    void save()
    {
        BYTE *pSaveImg = (unsigned char *)malloc(
            sizeof(unsigned char)*(int)(pow(float(iSamplingRate), 2)));
        BYTE *pSave = pSaveImg;

        // Relative paths are relative to the current working directory, 
        // not the path of the executable. 
        // The current working directory is the directory from which you started the program.
        string sDir = "E:/repositories/cpp-prefiltering/frequency_transport_test/visualstudio/freqTransTest/Debug/images/";
        string sExt = ".bmp";
        string sFileName = sDir + sLastOp + sExt;
        for(int i=0; i<iSamplingRate; i++)
            for (int k = 0; k < iSamplingRate; k++)
            {
                *(pSave++) = (BYTE)*ijklo(dpMatrix, iSamplingRate, i, 0, k, 0, 0);
            }
        if (Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFileName.c_str()))
        {
            printf("primary image saved.\n");
        }
        else
        {
            printf("primary image saved failed.\n");
            return;
        }

        pSave = pSaveImg;
        string sFExt = "_F.bmp";
        string sFFileName = sDir + sLastOp + sFExt;
        double real, imag, N = pow(1.0*iSamplingRate, iDims);
        for (int i = 0; i<iSamplingRate; i++)
            for (int k = 0; k < iSamplingRate; k++)
            {
                *(pSave++) = (BYTE)amp(ijklo(cpFmatrix, iSamplingRate, i, 0, k, 0, 0));
            }
        bool bSuc1 = Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFFileName.c_str());
        
        if (bSuc1)
        {
            printf("saved.\n");
        }
        else
        {
            printf("fourier save failed.\n");
        }

        delete pSaveImg;
    }
};
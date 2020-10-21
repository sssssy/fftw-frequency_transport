#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <windows.h>

#include <fftw3.h>

#include "bmpFile.h"

using namespace std;
typedef unsigned char BYTE;

double PI = 3.1415926;
int UNDEFINED = 0;

class Spectrum
{
    public:
        Spectrum(int iSamplingRate) : iDims(4)
    {
        this->iSamplingRate = iSamplingRate;
        pMatrix = (double*)calloc(sizeof(double), (int)(pow(float(iSamplingRate), iDims)));
        pFmatrix = fftw_alloc_complex(sizeof(fftw_complex) * (int)(pow(float(iSamplingRate), iDims)));
        fpMatrix = fftw_alloc_complex(sizeof(fftw_complex) * (int)(pow(float(iSamplingRate), iDims)));
        fpFmatrix = fftw_alloc_complex(sizeof(fftw_complex) * (int)(pow(float(iSamplingRate), iDims)));
        for (int i = 0; i < iDims;)
        {
            pDimsIndex[i++] = iDims;
        }
    }
    ~Spectrum(){}

    int iTime = 0;
    string sName = "unamed";
    int iDims;
    int iSamplingRate;
    int iRadius = iSamplingRate / 2;
    double *pMatrix; // x * y * theta * phi
    fftw_complex *pFmatrix;
    string sLastOp = "NULL";
    bool bFakeBilinear = true;

    /* fftw stuff */
    fftw_complex *fpMatrix, *fpFmatrix;
    int* pDimsIndex = new int[iDims];
    fftw_plan fft;

    void setZero()
    { 
        double *pData;
        for (pData = pMatrix; pData < pMatrix + (int)pow(float(iSamplingRate), 4); pData++){
            *(pData) = 255;
        }
        sLastOp = "setZero";
        printf("setZero.\n");
        printM(pMatrix);
    }

    void setRect()
    {
        for (int i = 0; i < iSamplingRate; i++)
            for (int j = 0; j < iSamplingRate; j++)
                for (int k = 0; k < iSamplingRate; k++)
                    for (int l = 0; l < iSamplingRate; l++){
                        if(i>0.4*iSamplingRate && i<0.6*iSamplingRate 
                            && k>0.3*iSamplingRate && k<0.7*iSamplingRate)
                            pMatrix[((i*iSamplingRate+j)*iSamplingRate+k)*iSamplingRate+l] = 255;
                        else
                            pMatrix[((i * iSamplingRate + j) * iSamplingRate + k) * iSamplingRate + l] = 0;
                    }

        sLastOp = "setRect";
        printf("setRect.\n");
        printM(pMatrix);
    }

    void printM(double* m)
    {
        for (int i = 0; i < pDimsIndex[0]; i++)
            for (int j = 0; j < pDimsIndex[1]; j++)
            {
                for (int k = 0; k < pDimsIndex[2]; k++)
                {
                    for (int l = 0; l < pDimsIndex[3]; l++)
                    {
                        printf("%f ", m[((i*iSamplingRate + j)*iSamplingRate + k)*iSamplingRate + l]);
                        if (l == pDimsIndex[3] - 1) printf("\n");
                    }
                    if (k == pDimsIndex[2] - 1) printf("\n\n");
                }
                if (j == pDimsIndex[1] - 1) printf("========\n\n");
            }
    }

    void fourier()
    {
        auto fp = fpMatrix;
        auto fpEnd = (fftw_complex*)(fp + (int)pow(1.0*iSamplingRate, iDims));
        auto p = pMatrix;
        for (; fp < fpEnd; fp++)
        {
            (*fp)[0] = *(p++);
            (*fp)[1] = 0.0;
        }
        fft = fftw_plan_dft_r2c(iDims, pDimsIndex, (double*)pMatrix, fpFmatrix, FFTW_ESTIMATE);
        fftw_execute(fft);
        
        fftw_destroy_plan(fft);
    }

    void fourier2()
    {

    }

    void save()
    {
        BYTE *pSaveImg = (unsigned char *)malloc(
            sizeof(unsigned char)*(int)(pow(float(iSamplingRate), 2)));
        BYTE *pSave = pSaveImg;

        string sDir = "images/";
        string sExt = ".bmp";
        string sFileName = sDir + sLastOp + sExt;
        for(int i=0; i<iSamplingRate; i++)
            for (int j=0; j < iSamplingRate; j++)
                *(pSave++) = (BYTE)pMatrix[((i * iSamplingRate) * iSamplingRate + j) * iSamplingRate];
        Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFileName.c_str());

        pSave = pSaveImg;
        string sFExt = "_F.bmp";
        string sFFileName = sDir + sLastOp + sFExt;
        double real, imag, N = pow(1.0*iSamplingRate, iDims);
        for (int i = 0; i<iSamplingRate; i++)
            for (int j = 0; j < iSamplingRate; j++)
            {
                real = pFmatrix[((i * iSamplingRate) * iSamplingRate + j) * iSamplingRate][0];
                imag = pFmatrix[((i * iSamplingRate) * iSamplingRate + j) * iSamplingRate][1];
                *(pSave++) = 2 * sqrt(pow(real, 2) + pow(imag, 2)) / N;
            }
        Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFFileName.c_str());

        printf("saved.\n");

        delete pSaveImg;
    }
};
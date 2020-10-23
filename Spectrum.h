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
        : iTime(0), iDims(iDims), iSamplingRate(iSamplingRate), iRadius(floor(iSamplingRate/2.0)), 
        sName(sName), sLastOp("NULL"), bFakeBilinear(true)
    {
        for (int i = 0; i < iDims; i++)
        {
            ipDimsLength[i] = iSamplingRate;
        }

        dpMatrix = (double*)malloc(
            sizeof(double)*(int)pow(1.0*iSamplingRate, iDims));
        cpFmatrix = (fftw_complex*)fftw_malloc(
            sizeof(fftw_complex)*(int)pow(1.0*iSamplingRate, iDims));
        cpIn = (fftw_complex*)fftw_malloc(
            sizeof(fftw_complex)*(int)pow(1.0*iSamplingRate, iDims));
        cpShift = (fftw_complex*)fftw_malloc(
            sizeof(fftw_complex)*(int)pow(1.0*iSamplingRate, iDims));
        fftn = fftw_plan_dft(iDims, ipDimsLength, cpIn, cpFmatrix, FFTW_FORWARD, FFTW_MEASURE);
    }
    ~Spectrum()
    {
        delete dpMatrix;
        fftw_free(cpFmatrix);
        fftw_free(cpIn);
        fftw_free(cpShift);
    }

    int iTime;
    int iDims;
    int iSamplingRate;
    int iRadius;
    string sName;
    string sLastOp;
    bool bFakeBilinear;
    fftw_plan fftn;
    int ipDimsLength[4];

    /* matrice */
    // x * y * theta * phi
    double *dpMatrix;
    fftw_complex *cpFmatrix, *cpIn, *cpShift;


    void setConst(double n)
    { 
        /*
        * output double value in [0., n.]
        */
        double *dp, *dpEnd = dpMatrix + (int)pow(float(iSamplingRate), 4);
        for (dp = dpMatrix; dp < dpEnd; dp++)
        {
            *dp = n;
        }
        sLastOp = "_setConst";
        printf("setConst.\n");
    }

    void setCos()
    {
        /*
        * output double value in [0., 1.]
        */   
        double *dist = (double*)malloc(sizeof(double)*iSamplingRate);
        double *dp, *dpEnd = dist + iSamplingRate;
        double i = 0.0;

        for (dp = dist; dp < dpEnd; dp++, i++)
        {
            *dp = -cos(2*PI*i/iSamplingRate)*0.5+0.5;
        }

        for (int i = 0; i < iSamplingRate; i++)
            for (int j = 0; j < iSamplingRate; j++)
                for (int k = 0; k < iSamplingRate; k++)
                    for (int l = 0; l < iSamplingRate; l++)
                        *ijklo(dpMatrix, iSamplingRate, i, j, k, l, 0) = dist[k];
        
        sLastOp = "_setCos";
        printf("setCos. \n");
    }

    void setRect()
    {
        /*
        * output double value in [0., 1.]
        */
        for (int i = 0; i < iSamplingRate; i++)
            for (int j = 0; j < iSamplingRate; j++)
                for (int k = 0; k < iSamplingRate; k++)
                    for (int l = 0; l < iSamplingRate; l++)
                    {
                        if (i >= 0.3*(iSamplingRate-1) && i <= 0.7*(iSamplingRate-1)
                            && k >= 0.4*(iSamplingRate-1) && k <= 0.6*(iSamplingRate-1))
                        {
                            *ijklo(dpMatrix, iSamplingRate, i, j, k, l, 0) = 1.0;
                        }
                        else
                            *ijklo(dpMatrix, iSamplingRate, i, j, k, l, 0) = 0.0;
                    }
        sLastOp = "_setRect";
        printf("setRect.\n");
    }

    void fourier()
    {
        //printf("dpMatrix: \n");
        real2complex(dpMatrix, cpIn, iSamplingRate);  
        //print(dpMatrix);
        //printf("cpIn: \n");
        //print(cpIn);

        fftw_execute(fftn);
        //printf("cpFmatrix: \n");
        //print(cpFmatrix);

        fftscaling(cpFmatrix, iSamplingRate, iDims);
        //printf("cpFmatrix after scaling: \n");
        //print(cpFmatrix);
        
        fftshift4d(cpFmatrix, cpShift, iSamplingRate);
        //printf("cpShift: \n");
        //print(cpShift);
    }

    void print(double* m)
    {
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                printf("%.3f, ", *ijklo(m, iSamplingRate, i, 0, k, 0, 0));
            }
            printf("\n");
        }
        printf("\n");
    }

    void print(fftw_complex* c)
    {
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                printf("|%.3f|, ", amp(ijklo(c, iSamplingRate, i, 0, k, 0, 0)));
            }
            printf("\n");
        }
        printf("\n");
    }

    void save()
    {
        // Relative paths are relative to the current working directory, 
        // not the path of the executable. 
        // The current working directory is the directory from which you started the program.
        string sDir = "E:/repositories/cpp-prefiltering/frequency_transport_test/visualstudio/freqTransTest/Debug/images/";
        
        BYTE *pSaveImg = (unsigned char *)malloc(
            sizeof(unsigned char)*(int)(pow(float(iSamplingRate), 2)));
        BYTE *pEnd = pSaveImg + iSamplingRate * iSamplingRate;

        BYTE *pSave = pSaveImg;
        double *dp;
        string sExt = ".bmp";
        string sFileName = sDir + sName + sLastOp + sExt;
        double max, min, pixel;
        bool bSuc = true;


        /* save dpMatrix */

        max = 0.0, min = 0.0;
        // theta first (for y axis)
        for(int k=0; k<iSamplingRate; k++)
            for (int i = 0; i < iSamplingRate; i++)
            {
                pixel = *ijklo(dpMatrix, iSamplingRate, i, 0, k, 0, 0);
                min = pixel < min ? pixel : min;
                max = pixel > max ? pixel : max;
            }

        max -= min;
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                *(pSave++) = (BYTE)(*ijklo(dpMatrix, iSamplingRate, i, 0, k, 0, 0) *255.0 / max);
            }
        }

        /*pSave = pSaveImg;
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                printf("%d ", *pSave++);
            }
            printf("\n");
        }*/

        bSuc &= Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFileName.c_str());


        /* save cpFmatrix */

        pSave = pSaveImg;
        sExt = "_F.bmp";
        sFileName = sDir + sName + sLastOp + sExt;
        max = 0.0, min = 0.0;
        // theta first (for y axis)
        for (int k = 0; k<iSamplingRate; k++)
            for (int i = 0; i < iSamplingRate; i++)
            {
                pixel = amp(ijklo(cpFmatrix, iSamplingRate, i, 0, k, 0, 0));
                min = pixel < min ? pixel : min;
                max = pixel > max ? pixel : max;
            }

        max -= min;
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                *(pSave++) = (BYTE)(amp(ijklo(cpFmatrix, iSamplingRate, i, 0, k, 0, 0)) *255.0 / max);
            }
        }

       /* pSave = pSaveImg;
        printf("cpFmatrix: \n");
        for (int k = 0; k < iSamplingRate; k++)
            {
                for (int i = 0; i < iSamplingRate; i++)
                {
                   printf("%d ", *pSave++);
                }
            printf("\n");
        }*/

        bSuc &= Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFileName.c_str());


        /* save cpShift */

        pSave = pSaveImg;
        sExt = "_F_shift.bmp";
        sFileName = sDir + sName + sLastOp + sExt;
        max = 0.0, min = 0.0;
        // theta first (for y axis)
        for (int k = 0; k<iSamplingRate; k++)
            for (int i = 0; i < iSamplingRate; i++)
            {
                pixel = amp(ijklo(cpShift, iSamplingRate, i, 0, k, 0, 0));
                min = pixel < min ? pixel : min;
                max = pixel > max ? pixel : max;
            }

        max -= min;
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                *(pSave++) = (BYTE)(amp(ijklo(cpShift, iSamplingRate, i, 0, k, 0, 0)) *255.0 / max);
            }
        }

        /*pSave = pSaveImg;
        printf("cpShift: \n");
        for (int k = 0; k < iSamplingRate; k++)
        {
            for (int i = 0; i < iSamplingRate; i++)
            {
                printf("%d ", *pSave++);
            }
            printf("\n");
        }*/

        bSuc &= Rmw_Write8BitImg2BmpFile(pSaveImg, iSamplingRate, iSamplingRate, sFileName.c_str());
        
        delete pSaveImg;
        if (bSuc)
            printf("saved.\n");
        else
            printf("error when saving.\n");
    }
};
#include "Spectrum.h"

int main()
{
    Spectrum s(100);
    s.setZero();
    //s.setRect();
    //s.print(s.dpMatrix);
    s.save();

#if 0
    fftw_cleanup();

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
            in[i*N + j][0] = i * 5 + j +rand() / (double)(RAND_MAX + 1);
            in[i*N + j][1] = 0;
        }

    /* show array */
    printf("Primary: \n\n");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%.3f+%.3fi,\t", in[i*N+j][0], in[i*N+j][1]);
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

#endif
    //getchar();
    return 0;
}


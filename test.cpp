#include "Spectrum.h"

int main()
{
    fftw_cleanup();

    int N = 5;
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    /* init in array */
    // 必须在创建plan之后再初始化数据，因为FFTW_MEASURE会覆盖你的数据。
    for (int i = 0; i < N; i++)
    {
        in[i][0] = i;
        in[i][1] = 0;
    }

    /* show array */
    printf("Primary: \n\n");
    for (int i = 0; i < N; i++)
    {
        printf("%f+%fi \n", in[i][0], in[i][1]);
    }

    fftw_execute(p);
    printf("\nFourier: \n\n");

    /* show spectrum array */
    for (int i = 0; i < N; i++)
    {
        printf("%f+%fi \n", out[i][0], out[i][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    //getchar();
    return 0;
}

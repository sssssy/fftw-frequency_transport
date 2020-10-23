#include "Spectrum.h"

int main()
{
    fftw_cleanup();
    
    Spectrum s(71, 4, "test");
    //s.setConst(1.0);
    s.setCos();
    //s.setRect();
    s.fourier();
    s.save();

    //getchar();
    return 0;
}

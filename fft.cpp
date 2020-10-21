#include <complex.h>
#include <math.h>

using namespace std;
typedef unsigned char BYTE;

//   FFT implementation
void fft(complex<BYTE> * Data, unsigned int N, bool Inverse = false)
{
   const double pi = Inverse ? 3.14159265358979323846 : -3.14159265358979323846;
   //   Iteration through dyads, quadruples, octads and so on...
   for (unsigned int Step = 1; Step < N; Step <<= 1)
   {
      //   Jump to the next entry of the same transform factor
      const unsigned int Jump = Step << 1;
      //   Angle increment
      const double delta = pi / double(Step);
      //   Auxiliary sin(delta / 2)
      const double Sine = sin(delta * .5);
      //   Multiplier for trigonometric recurrence
      const complex<BYTE> Multiplier(-2. * Sine * Sine, sin(delta));
      //   Start value for transform factor, fi = 0
      complex<BYTE> Factor(1.);
      //   Iteration through groups of different transform factor
      for (unsigned int Group = 0; Group < Step; ++Group)
      {
         //   Iteration within group 
         for (unsigned int Pair = Group; Pair < N; Pair += Jump)
         {
            //   Match position
            const unsigned int Match = Pair + Step;
            //   Second term of two-point transform
            const complex<BYTE> Product(Factor * Data[Match]);
            //   Transform for fi + pi
            Data[Match] = Data[Pair] - Product;
            //   Transform for fi
            Data[Pair] += Product;
         }
         //   Successive transform factor via trigonometric recurrence
         Factor = Multiplier * Factor + Factor;
      }
   }
}


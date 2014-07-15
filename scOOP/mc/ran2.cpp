#include "ran2.h"


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double Ran2::operator ()() {
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (seed <= 0) {
        if (-(seed) < 1) seed=1;
        else seed = -(seed);
        idum2=(seed);
        for (j=NTAB+7;j>=0;j--) {
            k=(seed)/IQ1;
            seed=IA1*(seed-k*IQ1)-k*IR1;
            if (seed < 0) seed += IM1;
            if (j < NTAB) iv[j] = seed;
        }
        iy=iv[0];
    }
    k=(seed)/IQ1;
    seed=IA1*(seed-k*IQ1)-k*IR1;
    if (seed < 0) seed += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = seed;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

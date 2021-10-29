#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <time.h> 

using namespace NTL;

/*
 * Ring Regev SPHF
 */
class rregsphf{
    Vec<ZZ_pE> hk;
    
    ZZ_pE sp;
    Vec<ZZ_pE> s;    
  public:
    Vec<ZZ_pE> ph;

    ZZ q;
    long d, N, mu, lambda;
    Vec<ZZ_pE> Ap, b;
    Mat<ZZ_pE> A;
    rregsphf(long N, long d, ZZ q, Vec<ZZ_pE> Ap, Vec<ZZ_pE> b);
    ZZ_pE Round(ZZ_pE x);
    void HashKG();
    void ProjKG();
    ZZ_pE Hash(Vec<ZZ_pE> c, ZZ_pE m);
    ZZ_pE ProjHash(Vec<ZZ_pE> c, ZZ_pE m, Vec<ZZ_pE> w);
};

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>
#include <time.h> 
#include <mp.h> 

using namespace NTL;

/*
 * Plain Regev SPHF
 */
class mpsphf{
    Vec<ZZ_p> hk;
    
    Vec<ZZ_p> sp;
    Vec<ZZ_p> s;    
  public:
    mp *scheme;
    Vec<ZZ_p> ph;

    ZZ q;
    long m, N, mu, lambda;
    Vec<ZZ_p> b;
    Mat<ZZ_p> Ap;
    Mat<ZZ_p> A;
    mpsphf(long N, long d, ZZ q, Mat<ZZ_p> A, mp *scheme);
    ZZ_p Round(ZZ_p x);
    void HashKG();
    void ProjKG();
    ZZ_p Hash(Vec<ZZ_p> c, Vec<ZZ_p> m);
    ZZ_p ProjHash(Vec<ZZ_p> c, Vec<ZZ_p> m, Vec<ZZ_p> w);
};

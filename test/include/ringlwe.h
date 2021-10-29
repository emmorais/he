#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <time.h> 

using namespace NTL;

class ringlwe{
    double delta;
    ZZ_pE sp;
    Vec<ZZ_pE> s;    
  public:
    ZZ q;
    long d, N, mu, lambda;
    Vec<ZZ_pE> Ap, b;
    Mat<ZZ_pE> A;
    ringlwe(long N, long d, ZZ q, double delta);
    Vec<ZZ_pE> Mult(Mat<ZZ_pE> a, Vec<ZZ_pE> b);
    Mat<ZZ_pE> Transpose(Mat<ZZ_pE> a);
    ZZ_pE InnerProduct(Vec<ZZ_pE> a, Vec<ZZ_pE> b);
    ZZ_pE SampleMessage();	    
    Vec<ZZ_pE> SampleRandom();
    Vec<ZZ_pE> SampleR();
    Vec<ZZ_pE> SampleNoiseArray();
    ZZ_pE SampleNoise();
    ZZ_pE SampleKey();
    void KeyGen();
    Vec<ZZ_pE> EncryptR(ZZ_pE m, Vec<ZZ_pE> r);
    Vec<ZZ_pE> Encrypt(ZZ_pE m);
    ZZ_pE Decrypt(Vec<ZZ_pE> c);
    ZZ_pE Mod2(ZZ_pE a);
};


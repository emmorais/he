#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <time.h> 

using namespace NTL;

class plainlwe{
    double delta;
    Vec<ZZ_p> sp;
    Vec<ZZ_p> s;    
  public:
    ZZ q;
    long m, N, mu, lambda;
    Vec<ZZ_p> b;
    Mat<ZZ_p> Ap;
    Mat<ZZ_p> A;
    plainlwe(long N, long m, ZZ q);
    Vec<ZZ_p> Add(Vec<ZZ_p> a, Vec<ZZ_p> b);
    Vec<ZZ_p> Mult(Mat<ZZ_p> a, Vec<ZZ_p> b);
    Mat<ZZ_p> Transpose(Mat<ZZ_p> a);
    ZZ_p InnerProduct(Vec<ZZ_p> a, Vec<ZZ_p> b);
    ZZ_p SampleMessage();	    
    Mat<ZZ_p> SampleRandom();
    Vec<ZZ_p> SampleR();
    Vec<ZZ_p> SampleNoiseArray(long n);
    ZZ_p SampleNoise();
    ZZ_p SampleKey();
    Vec<ZZ_p> Double(Vec<ZZ_p> x);
    void KeyGen();
    Vec<ZZ_p> EncryptR(ZZ_p m, Vec<ZZ_p> r);
    Vec<ZZ_p> Encrypt(ZZ_p m);
    ZZ_p Decrypt(Vec<ZZ_p> c);
    ZZ_p Mod2(ZZ a);
};


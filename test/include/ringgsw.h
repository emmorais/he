#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <time.h> 

using namespace NTL;

class ringgsw{
    double delta;
    Vec<ZZ_p> sp;
    Vec<ZZ_p> s;    
  public:
    Mat<ZZ_p> G;
    ZZ q;
    long Q, m, N, n, mu, lambda;
    Vec<ZZ_p> b;
    Mat<ZZ_p> Ap;
    Mat<ZZ_p> A;
    ringgsw(long N, long n, long m, ZZ q, long Q);
    Vec<ZZ_p> Add(Vec<ZZ_p> a, Vec<ZZ_p> b);
    Mat<ZZ_p> Add(Mat<ZZ_p> a, Mat<ZZ_p> b);
    Mat<ZZ_p> Mult(Mat<ZZ_p> a, Mat<ZZ_p> b);
    Vec<ZZ_p> Mult(Mat<ZZ_p> a, Vec<ZZ_p> b);
    Mat<ZZ_p> Transpose(Mat<ZZ_p> a);
    ZZ_p SampleMessage();	    
    Mat<ZZ_p> SampleRandom();
    Vec<ZZ_p> SampleRandomArray(int len);
    Mat<ZZ_p> SampleRMatrix();
    Vec<ZZ_p> SampleR();
    Vec<ZZ_p> SampleNoiseArray(long n);
    Vec<ZZ_p> SampleShortNoise(int len);
    ZZ_p SampleNoise();
    ZZ_p SampleKey();
    void KeyGen();
    //Vec<ZZ_p> EncryptR(ZZ_p m, Vec<ZZ_p> r);
    Mat<ZZ_p> Encrypt(ZZ_p m);
    ZZ_p Decrypt(Mat<ZZ_p> c);
    ZZ_p Mod2(ZZ a);
    void ComputeG();
};


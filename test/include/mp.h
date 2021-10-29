#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <time.h> 
#include <NTL/mat_ZZ.h>

using namespace NTL;

class mp{
    double delta;
    Vec<ZZ_p> sp;
    Vec<ZZ_p> s;    
    Mat<ZZ> invS; 
    Mat<ZZ> HNF_SS; 
    Mat<ZZ> invHNF_SS; 
    Mat<ZZ> SS; 
    Mat<ZZ> invSS; 
    Mat<ZZ> tSS; 
    Mat<ZZ> tinvSS; 
    Mat<ZZ_p> RvID; 
  public:
    Mat<ZZ_p> R; 
    Mat<ZZ> S; 
    ZZ q;
    long m, Q, N, mu, lambda;
    Vec<ZZ_p> b;
    Mat<ZZ_p> Ap, A, G, H;
    mp(long N, long m, ZZ q, long Q);
    Vec<ZZ_p> Add(Vec<ZZ_p> a, Vec<ZZ_p> b);
    Vec<ZZ_p> Sub(Vec<ZZ_p> a, Vec<ZZ_p> b);
    Vec<ZZ_p> MultT(Mat<ZZ_p> a, Vec<ZZ_p> b);
    Vec<ZZ_p> Mult(Mat<ZZ_p> a, Vec<ZZ_p> b);
    Mat<ZZ_p> Transpose(Mat<ZZ_p> a);
    Mat<ZZ> Transpose(Mat<ZZ> a);
    ZZ_p InnerProduct(Vec<ZZ_p> a, Vec<ZZ_p> b);
    Vec<ZZ_p> SampleMessage();	    
    ZZ_p SampleM();
    Mat<ZZ_p> SampleRandom();
    Vec<ZZ_p> SampleS();
    Vec<ZZ_p> SampleNoiseArray(long n);
    ZZ_p SampleNoise();
    ZZ_p SampleKey();
    Vec<ZZ_p> Double(Vec<ZZ_p> x);
    void KeyGen();
    Vec<ZZ_p> g(Vec<ZZ_p> s, Vec<ZZ_p> e);
    Vec<ZZ_p> EncryptR(Vec<ZZ_p> m, Vec<ZZ_p> s, Vec<ZZ_p> e);
    Vec<ZZ_p> Encrypt(Vec<ZZ_p> m);
    Vec<ZZ_p> DecryptR(Vec<ZZ_p> b, Vec<ZZ_p> e);
    Vec<ZZ_p> Decrypt(Vec<ZZ_p> c);
    ZZ_p Mod2(ZZ a);

    Mat<ZZ> SetIDZZ(int n);
    Mat<ZZ_p> SetID(int n);
    Mat<ZZ_p> SampleR(int n, int m);
    Vec<ZZ_p> concat(Vec<ZZ_p> a, Vec<ZZ_p> b);
    Mat<ZZ_p> concat(Mat<ZZ_p> a, Mat<ZZ_p> b);
    Mat<ZZ_p> concatv(Mat<ZZ_p> a, Mat<ZZ_p> b);
    Mat<ZZ_p> mul(Mat<ZZ_p> a, Mat<ZZ_p> b);
    Vec<ZZ> mulv(Mat<ZZ> a, Vec<ZZ> b);
    Mat<ZZ_p> add(Mat<ZZ_p> a, Mat<ZZ_p> b);
    Mat<ZZ_p> sub(Mat<ZZ_p> a, Mat<ZZ_p> b);
    void ComputeG();
    void ComputeRvID();
    int D();
    void GenTrap(Mat<ZZ_p> Ab, Mat<ZZ_p> H);
    void InvGA(Vec<ZZ_p> &s, Vec<ZZ_p> &e, Vec<ZZ_p> b); 
    void InvGOracle(ZZ_p &s, Vec<ZZ_p> &e, Vec<ZZ_p> b);
    void ComputeHermiteS();
    void ComputeS();
    void ComputeSS();
    Vec<ZZ> Encode(Vec<ZZ> m);
    Vec<ZZ> Decode(Vec<ZZ> m);
    Vec<ZZ_p> Zeros(int n);
    Vec<ZZ_p> ScalarMul(Vec<ZZ_p> x, ZZ_p c);
    Vec<ZZ> ConvVecToZZ(Vec<ZZ_p> x);
    Vec<ZZ_p> ConvVecToZZp(Vec<ZZ> x);
};

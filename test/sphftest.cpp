#include <rregsphf.h>
#include <pregsphf.h>
#include <ringlwe.h>
#include <plainlwe.h>
#include <colors.h>
#include <mpsphf.h>
#include <iostream>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

#define NUM 100
#define N 20
#define M 12
#define Q 19
#define D 128

bool rregsphfTest(){
  int count=0, ok=0, nok=0;
  rregsphf *scheme;
  ringlwe *scheme_ringlwe;
  long n = N;
  long d = D;
  ZZ_pE b, bb;
  ZZ q;
  Vec<ZZ_pE> w;
  Vec<ZZ_pE> Ap;

  GenPrime(q, Q);
  ZZ_p::init(conv<ZZ>(q)); 
  ZZ_pX P;
  SetCoeff(P, 0, 1);
  SetCoeff(P, d, 1);
  ZZ_pE::init(P); 
  scheme_ringlwe = new ringlwe(n, d, q, 1);
  scheme_ringlwe->KeyGen();
  scheme = new rregsphf(n, d, q, scheme_ringlwe->Ap, scheme_ringlwe->b);
  
  while(count < NUM){
    scheme->HashKG();
    scheme->ProjKG();
    
    w = scheme_ringlwe->SampleR();
    ZZ_pE m = scheme_ringlwe->SampleMessage();
    Vec<ZZ_pE> c = scheme_ringlwe->EncryptR(m, w);
    ZZ_pE mm = scheme_ringlwe->Decrypt(c);
    
    b = scheme->Hash(c, m); 
    bb = scheme->ProjHash(c, m, w);
    
    if(b == bb)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}

bool pregsphfTest(){
  int count=0, ok=0, nok=0;
  pregsphf *scheme;
  plainlwe *scheme_plainlwe;
  long n = N;
  long m = M;
  ZZ_p b, bb;
  ZZ q;
  Vec<ZZ_p> w;
  Vec<ZZ_p> Ap;

  GenPrime(q, Q);
  ZZ_p::init(conv<ZZ>(q)); 
  scheme_plainlwe = new plainlwe(n, m, q);
  scheme_plainlwe->KeyGen();
  scheme = new pregsphf(n, m, q, scheme_plainlwe->Ap, scheme_plainlwe->b);
  
  while(count < NUM){
    scheme->HashKG();
    scheme->ProjKG();
    
    w = scheme_plainlwe->SampleR();
    ZZ_p m = scheme_plainlwe->SampleMessage();
    Vec<ZZ_p> c = scheme_plainlwe->EncryptR(m, w);
    ZZ_p mm = scheme_plainlwe->Decrypt(c);
    
    b = scheme->Hash(c, m); 
    bb = scheme->ProjHash(c, m, w);
    
    if(b == bb)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}

bool mpsphfTest(){
  int count=0, ok=0, nok=0;
  mpsphf *scheme;
  mp *scheme_mp;
  long n = 1;
  long m = 1*9;
  ZZ_p b, bb;
  ZZ q;
  Vec<ZZ_p> w;
  Vec<ZZ_p> Ap;
  Mat<ZZ_p> s;
  Vec<ZZ_p> e;

  GenPrime(q, 9);
  q = conv<ZZ>(509);
  ZZ_p::init(conv<ZZ>(q)); 
  scheme_mp = new mp(n, m, q, 9);
  scheme_mp->KeyGen();
  scheme = new mpsphf(n, m, q, scheme_mp->A, scheme_mp);
  
  while(count < NUM){
    scheme->HashKG();
    scheme->ProjKG();
    
    s = scheme_mp->SampleR(n, m);
    s[0][0] = 1;
    cout << "s:     " << s << ";\n";
    e = scheme_mp->SampleNoiseArray(scheme_mp->A.NumCols());
    Vec<ZZ_p> m = scheme_mp->SampleMessage();
    m[8] = 0;
    m[7] = 0;
    m[6] = 0;

    Vec<ZZ_p> c = scheme_mp->EncryptR(m, s[0], e);
    Vec<ZZ_p> mm = scheme_mp->Decrypt(c);
    cout << "c:  " << c << ";\n";
    cout << "m:  " << m << ";\n";
    cout << "mm: " << mm << ";\n";
    
    b = scheme->Hash(c, m); 
    bb = scheme->ProjHash(c, m, s[0]);
    
    cout << "Hash:     " << b << ";\n";
    cout << "ProjHash: " << bb << ";\n";
    if(b == bb)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}


int main(int argc, char **argv){
  srand (time(NULL));
  bool t1 = rregsphfTest();
  cout << FGRN("Ring Regev SPHF test result: ") << (t1 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  bool t2 = pregsphfTest();
  cout << FGRN("Plain Regev SPHF test result: ") << (t2 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  //bool t3 = mpsphfTest();
  //cout << FGRN("Micciancio-Peikert SPHF test result: ") << (t3 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  return 0;
}

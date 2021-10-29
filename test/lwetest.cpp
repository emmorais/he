#include <iostream>
#include <NTL/ZZ.h>
#include <colors.h>
#include <ringlwe.h>
#include <plainlwe.h>

using namespace std;
using namespace NTL;

#define NUM 100
#define N 240
#define M 12
#define Q 79
#define D 128

/*
 * Encryption/Decryption of plain LWE scheme.
 */
bool plainLWETest(){
  int count=0, ok=0, nok=0;
  plainlwe *scheme;
  double delta = 1;
  long n = N;
  long m = M;
  ZZ q;
  GenPrime(q, Q);
  ZZ_p::init(conv<ZZ>(q)); 
  ZZ_pX P;
  scheme = new plainlwe(n, m, q);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_p m = scheme->SampleMessage();
    Vec<ZZ_p> c1 = scheme->Encrypt(m);
    ZZ_p mm = scheme->Decrypt(c1);
    if(m == mm)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}

/*
 * Encryption/Decryption of ring LWE scheme.
 */
bool ringLWETest(){
  int count=0, ok=0, nok=0;
  ringlwe *scheme;
  double delta = 1;
  long n = N;
  long d = D;
  ZZ q;
  GenPrime(q, Q);
  ZZ_p::init(conv<ZZ>(q)); 
  ZZ_pX P;
  SetCoeff(P, 0, 1);
  SetCoeff(P, d, 1);
  ZZ_pE::init(P); 
  scheme = new ringlwe(n, d, q, delta);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_pE m = scheme->SampleMessage();
    Vec<ZZ_pE> c1 = scheme->Encrypt(m);
    ZZ_pE mm = scheme->Decrypt(c1);
    if(m == mm)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}

/*
 * Addition homomorphism of ring LWE scheme.
 */
bool ringLWEAddTest(){
  int count=0, ok=0, nok=0;
  ringlwe *scheme;
  double delta = 1;
  long n = N;
  long d = D;
  ZZ q;
  GenPrime(q, Q);
  ZZ_p::init(conv<ZZ>(q)); 
  ZZ_pX P;
  SetCoeff(P, 0, 1);
  SetCoeff(P, d, 1);
  ZZ_pE::init(P); 
  scheme = new ringlwe(n, d, q, delta);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_pE m1 = scheme->SampleMessage();
    ZZ_pE m2 = scheme->SampleMessage();
    Vec<ZZ_pE> c1 = scheme->Encrypt(m1);
    Vec<ZZ_pE> c2 = scheme->Encrypt(m2);
    ZZ_pE mm = scheme->Decrypt(c1+c2);
    if(scheme->Mod2(m1+m2) == mm)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}


int main(){
  srand (time(NULL));
  bool t1 = ringLWETest();
  cout << FGRN("Ring BGV Encryption/Decryption test result: ") << (t1 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  bool t2 = plainLWETest();
  cout << FGRN("Plain BGV Encryption/Decryption test result: ") << (t2 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  bool t3 = ringLWEAddTest();
  cout << FGRN("Addition Homomorphism test result: ") << (t3 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  return 0;
}

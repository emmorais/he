#include <iostream>
#include <NTL/ZZ.h>
#include <util.h>
#include <sample.h>
#include <gsw.h>
#include <ringgsw.h>
#include <colors.h>

using namespace std;
using namespace NTL;

#define NUM 100
#define N 64
#define M 50
#define Q 29

/*
 * Tests Encryption and Decryption of ring GSW.
 */
bool ringgswTest(){
  int count=0, ok=0, nok=0;
  ringgsw *scheme;
  double delta = 1;
  long n = 1;
  long m = M;
  ZZ q;
  q = conv<ZZ>("536870912");
  ZZ_p::init(conv<ZZ>(q)); 
  scheme = new ringgsw(N, n, m, q, Q);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_p m = sample::SampleMessage();
    Mat<ZZ_p> c1 = scheme->Encrypt(m);
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
 * Tests Encryption and Decryption of GSW.
 */
bool gswTest(){
  int count=0, ok=0, nok=0;
  gsw *scheme;
  double delta = 1;
  long n = 1;
  long m = M;
  ZZ q;
  q = conv<ZZ>("536870912");
  ZZ_p::init(conv<ZZ>(q)); 
  scheme = new gsw(N, n, m, q, Q);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_p m = sample::SampleMessage();
    Mat<ZZ_p> c1 = scheme->Encrypt(m);
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
 * Tests addition homomorphism.
 */
bool gswAddTest(){
  int count=0, ok=0, nok=0;
  gsw *scheme;
  double delta = 1;
  long n = 1;
  long m = M;
  ZZ q;
  q = conv<ZZ>("536870912");
  ZZ_p::init(conv<ZZ>(q)); 
  scheme = new gsw(N, n, m, q, Q);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_p m1 = sample::SampleMessage();
    ZZ_p m2 = sample::SampleMessage();
    Mat<ZZ_p> c1 = scheme->Encrypt(m1);
    Mat<ZZ_p> c2 = scheme->Encrypt(m2);
    ZZ_p mm1 = scheme->Decrypt(c1);
    ZZ_p mm2 = scheme->Decrypt(c2);
    ZZ_p mm1mm2 = scheme->Decrypt(util::Add(c1, c2));
    if(conv<ZZ>(m1+m2)%2 == mm1mm2)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}

/*
 * Tests multiplication homomorphism.
 */
bool gswMulTest(){
  int count=0, ok=0, nok=0;
  gsw *scheme;
  double delta = 1;
  long n = 1;
  long m = M;
  ZZ q;
  q = conv<ZZ>("536870912");
  ZZ_p::init(conv<ZZ>(q)); 
  scheme = new gsw(N, n, m, q, Q);
  
  while(count < NUM){
    scheme->KeyGen();
    ZZ_p m1 = sample::SampleMessage();
    ZZ_p m2 = sample::SampleMessage();
    Mat<ZZ_p> c1 = scheme->Encrypt(m1);
    Mat<ZZ_p> c2 = scheme->Encrypt(m2);
    ZZ_p mm1 = scheme->Decrypt(c1);
    ZZ_p mm2 = scheme->Decrypt(c2);
    ZZ_p mm1mm2 = scheme->Decrypt(util::Mult(c1,c2));
    if(conv<ZZ>(m1*m2)%2 == mm1*mm2)
      ok++;
    else
      nok++;
    count++;
  }
  int result = (ok == count); 
  return result; 
}

/*
 * Runs all tests.
 */
int main(){
  srand (time(NULL));
  bool t1 = gswTest();
  cout << FGRN("GSW Encryption/Decryption test result: ") << (t1 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  bool t2 = gswAddTest();
  cout << FGRN("GSW Addition Homomorphism test result: ") << (t2 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  bool t3 = gswMulTest();
  cout << FGRN("GSW Multiplication Homomorphism test result: ") << (t3 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";
  bool t4 = ringgswTest();
  cout << FGRN("Ring GSW Encryption/Decryption test result: ") << (t4 ? BOLD(FGRN("True")) : BOLD(FRED("False"))) << "\n";

  return 0;
}

#include <iostream>
#include <NTL/ZZ.h>
#include <mp.h>


using namespace std;
using namespace NTL;

#define NUM 100
#define N 10
#define Q 16

void mpEncTest(){
  int i, count=0, ok=0, nok=0;
  mp *scheme_mp;
  long n = N;
  long m = N*Q;
  ZZ_p b, bb;
  ZZ q;
  Vec<ZZ_p> w;
  Vec<ZZ_p> Ap;
  Vec<ZZ_p> s;
  Vec<ZZ_p> e;

  q = conv<ZZ>(power(conv<ZZ>(2), Q));

  ZZ_p::init(conv<ZZ>(q)); 
  cout << "Micciancio-Peikert Encryption scheme test:\n";
  scheme_mp = new mp(n, m, q, Q);
  
  int achou = 0;
  while(!achou) {
    count = 0;
    ok = 0; nok = 0;
    scheme_mp->KeyGen();
  
    while(count < NUM){
      s = scheme_mp->SampleS();
      cout << "s: " << s << ";\n";
      
      // change mod to 2q
      ZZ_p::init(conv<ZZ>(q)*conv<ZZ>(2)); 
      e = scheme_mp->SampleNoiseArray(scheme_mp->A.NumCols());
      cout << "e: " << e << "\n";
      // change mod back to q
      ZZ_p::init(conv<ZZ>(q)); 
      
      Vec<ZZ_p> m = scheme_mp->SampleMessage();
      Vec<ZZ_p> c = scheme_mp->EncryptR(m, s, e);
      Vec<ZZ_p> mm = scheme_mp->Decrypt(c);
      cout << "c:  " << c << ";\n";
      cout << "m:  " << m << ";\n";
      cout << "mm: " << mm << ";\n";
      
      if(m == mm)
        ok++;
      else { nok++; break;}
      count++;
    }
    if (nok>0) { achou = 1; }
    cout << "ok,nok: " << ok << "," << nok << "\n";
  }
}

int main(){
  srand (time(NULL));
  mpEncTest();
  return 0;
}

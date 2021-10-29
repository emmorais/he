#include <mpsphf.h>
#include <math.h>

using namespace NTL;
using namespace std;

mpsphf::mpsphf(long N, long m, ZZ q, Mat<ZZ_p> A, mp *scheme){
  this->N = N;
  this->m = m;
  this->q = q;
  this->A = A;
  this->scheme = scheme;
}

void mpsphf::HashKG(){
  int i;
  Vec<ZZ_p> res;
  res.SetLength(this->m);
  for(i=0;i<this->m;i++){
    res[i] = random_ZZ_p();
  }
  this->hk = res;
  return;
}

void mpsphf::ProjKG(){
  int i,j;
  Vec<ZZ_p> res;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    for(j=0;j<this->m;j++){
      res[i] += this->A[i][j]*this->hk[j];
    }
  }
  this->ph = res;
  return; 
}

ZZ_p mpsphf::Round(ZZ_p x){
  int i;
  ZZ_p res = x;
  if(conv<ZZ>(res)<this->q/conv<ZZ>(2)){
    res = conv<ZZ_p>(conv<ZZ>(res)%conv<ZZ>(2));
  }
  else{
    res = conv<ZZ_p>(conv<ZZ>(1) - conv<ZZ>(res)%conv<ZZ>(2));
  }
  return res; 
}
 
ZZ_p mpsphf::Hash(Vec<ZZ_p> c, Vec<ZZ_p> m){
  int i;
  ZZ_p res;
  Vec<ZZ> mm;
  Vec<ZZ_p> cc;
  cc.SetLength(c.length());
  mm.SetLength(m.length());
  res = conv<ZZ_p>(conv<ZZ>(0));
  mm = this->scheme->Encode(this->scheme->ConvVecToZZ(m));
  cout << "mm: " << mm << "\n";
  for(i=0;i<cc.length();i++){
    cc[i] = c[i];
  }
  cc = this->scheme->Sub(cc, this->scheme->concat(this->scheme->Zeros(this->N*this->scheme->Q), this->scheme->ConvVecToZZp(mm)));
  for(i=0;i<m.length();i++){
    res += (cc[i])*this->hk[i];
  }
  for(i=m.length();i<this->m;i++){
    res += c[i]*this->hk[i];
  }  
  res = Round(res);
  return res; 
}

ZZ_p mpsphf::ProjHash(Vec<ZZ_p> c, Vec<ZZ_p> m, Vec<ZZ_p> w){
  int i;
  ZZ_p res;
  res = conv<ZZ_p>(conv<ZZ>(0));
  for (i=0;i<this->N;i++){
    res += this->ph[i]*w[i];
  }
  res = Round(res);
  return res;
}

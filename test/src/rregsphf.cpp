#include <rregsphf.h>
#include <math.h>

using namespace NTL;
using namespace std;

rregsphf::rregsphf(long N, long d, ZZ q, Vec<ZZ_pE> Ap, Vec<ZZ_pE> b){
  this->N = N;
  this->d = d;
  this->q = q;
  this->Ap = Ap;
  this->b = b;
}

void rregsphf::HashKG(){
  int i;
  Vec<ZZ_pE> res;
  res.SetLength(2);
  for(i=0;i<2;i++){
    random(res[i]);
  }
  this->hk = res;
  return;
}

void rregsphf::ProjKG(){
  int i;
  Vec<ZZ_pE> res;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = this->b[i]*this->hk[0] - this->Ap[i]*this->hk[1];
  }
  this->ph = res;
  return; 
}

ZZ_pE rregsphf::Round(ZZ_pE x){
  int i;
  ZZ_pE res;
  ZZ_pX pol;
  ZZ_p a;
  pol = conv<ZZ_pX>(x);
  for(i=0;i<this->d;i++){
    a = coeff(pol, i);
    if(conv<ZZ>(a)<this->q/conv<ZZ>(2)){
      SetCoeff(pol, i, conv<ZZ_p>(conv<ZZ>(a)%conv<ZZ>(2)));
    }
    else{
      SetCoeff(pol, i, conv<ZZ_p>(conv<ZZ>(1) - conv<ZZ>(a)%conv<ZZ>(2)));
    }
  }
  res = conv<ZZ_pE>(pol);
  return res; 
}
 
ZZ_pE rregsphf::Hash(Vec<ZZ_pE> c, ZZ_pE m){
  int i;
  ZZ_pE res;
  res = conv<ZZ_pE>(conv<ZZ>(0));
  res = (c[0]-m)*this->hk[0]+c[1]*this->hk[1];  
  res = Round(res);
  return res; 
}

ZZ_pE rregsphf::ProjHash(Vec<ZZ_pE> c, ZZ_pE m, Vec<ZZ_pE> w){
  int i;
  ZZ_pE res;
  res = conv<ZZ_pE>(conv<ZZ>(0));
  for (i=0;i<this->N;i++){
    res += this->ph[i]*w[i];
  }
  res = Round(res);
  return res;
}

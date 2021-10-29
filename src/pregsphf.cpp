#include <pregsphf.h>
#include <math.h>

using namespace NTL;
using namespace std;

pregsphf::pregsphf(long N, long m, ZZ q, Mat<ZZ_p> Ap, Vec<ZZ_p> b){
  this->N = N;
  this->m = m;
  this->q = q;
  this->Ap = Ap;
  this->b = b;
}

void pregsphf::HashKG(){
  int i;
  Vec<ZZ_p> res;
  res.SetLength(this->m+1);
  for(i=0;i<this->m+1;i++){
    res[i] = random_ZZ_p();
  }
  this->hk = res;
  return;
}

void pregsphf::ProjKG(){
  int i,j;
  Vec<ZZ_p> res;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = this->b[i]*this->hk[0];
    for(j=0;j<this->m;j++){
      res[i] -= this->Ap[i][j]*this->hk[1+j];
    }
  }
  this->ph = res;
  return; 
}

ZZ_p pregsphf::Round(ZZ_p x){
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
 
ZZ_p pregsphf::Hash(Vec<ZZ_p> c, ZZ_p m){
  int i;
  ZZ_p res;
  res = conv<ZZ_p>(conv<ZZ>(0));
  res = (c[0]-m)*this->hk[0];
  for(i=0;i<this->m;i++){
    res += c[i+1]*this->hk[i+1];
  }  
  res = Round(res);
  return res; 
}

ZZ_p pregsphf::ProjHash(Vec<ZZ_p> c, ZZ_p m, Vec<ZZ_p> w){
  int i;
  ZZ_p res;
  res = conv<ZZ_p>(conv<ZZ>(0));
  for (i=0;i<this->N;i++){
    res += this->ph[i]*w[i];
  }
  res = Round(res);
  return res;
}

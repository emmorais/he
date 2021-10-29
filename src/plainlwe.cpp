#include <plainlwe.h>
#include <math.h>

using namespace NTL;
using namespace std;

/*
 * Constructor
 */
plainlwe::plainlwe(long N, long m, ZZ q){
  this->N = N;
  this->m = m;
  this->q = q;
}

ZZ_p plainlwe::SampleMessage(){
  int i, mm;
  ZZ_p res;
  res = rand()%2;
  return res;
}
    
Mat<ZZ_p> plainlwe::SampleRandom(){
  int i,j;
  Mat<ZZ_p> res;
  res.SetDims(this->N, this->m);
  for(i=0;i<this->N;i++){
    for(j=0;j<this->m;j++){
      res[i][j] = random_ZZ_p();
    }
  }
  return res;
}

ZZ_p plainlwe::SampleNoise(){
  ZZ_p res;
  double randNormal, randStdNormal;
  double u1 = ((double) rand()) / (RAND_MAX); 
  double u2 = ((double) rand()) / (RAND_MAX); 
  u2 *= 6.2831;
  u1 = sqrt (-2 * log(u1));
  randNormal = u1 * sin(u2); 
  randNormal = 8 * randNormal; 
  randNormal = u1 * cos(u2); 
  randNormal = randNormal; 
  res = randNormal;
  return res;
}

Vec<ZZ_p> plainlwe::SampleNoiseArray(long n){
  Vec<ZZ_p> res;
  long i;
  res.SetLength(n);
  for(i=0;i<n;i++){
    res[i] = SampleNoise();
  }
  return res;
}

Vec<ZZ_p> plainlwe::SampleR(){
  Vec<ZZ_p> res;
  long i;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = conv<ZZ_p>(SampleMessage());
  }
  return res;
}

Vec<ZZ_p> plainlwe::Double(Vec<ZZ_p> x){
  long i;
  Vec<ZZ_p> res;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = 2*x[i];
  }
  return res;
}

void plainlwe::KeyGen(){	
  int i,j;
  Vec<ZZ_p> e;
  this->sp = SampleNoiseArray(this->m);
  this->s.SetLength(this->m+1);
  this->s[0] = conv<ZZ_p>(1);
  for(i=0;i<this->m;i++){
    this->s[i+1] = this->sp[i];
  }
  this->Ap = SampleRandom();  
  e.SetLength(this->N);
  e = SampleNoiseArray(this->N);
  e = Double(e);
  this->b = Add(Mult(this->Ap,this->sp), e);
  this->A.SetDims(this->N, this->m+1);
  for(i=0;i<this->N;i++){
    this->A[i][0] = this->b[i];
    for(j=0;j<this->m;j++){
      this->A[i][j+1] = -this->Ap[i][j];
    }
  }
}
    
Mat<ZZ_p> plainlwe::Transpose(Mat<ZZ_p> a){
  Mat<ZZ_p> res;
  int i,j;
  res.SetDims(a.NumCols(), a.NumRows());
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[j][i] = a[i][j];
    }
  }
  return res;
}

Vec<ZZ_p> plainlwe::Add(Vec<ZZ_p> a, Vec<ZZ_p> b){
  int i;
  Vec<ZZ_p> res;
  res.SetLength(a.length());
  for(i=0;i<a.length();i++){
    res[i] = a[i] + b[i];
  }
  return res;
}

Vec<ZZ_p> plainlwe::Mult(Mat<ZZ_p> a, Vec<ZZ_p> b){
  int i,j,k;
  Vec<ZZ_p> res;
  res.SetLength(a.NumRows());
  for(i=0;i<a.NumRows();i++){
    res[i] = conv<ZZ_p>(conv<ZZ>(0));
    for(k=0;k<b.length();k++){
      res[i] += a[i][k]*b[k];
    }
  }
  return res;
}

Vec<ZZ_p> plainlwe::EncryptR(ZZ_p m, Vec<ZZ_p> r){
  Vec<ZZ_p> c;
  c.SetLength(this->m+1);
  c = Mult(Transpose(this->A),r);
  c[0]+=m;
  return c;
}

Vec<ZZ_p> plainlwe::Encrypt(ZZ_p m){
  Vec<ZZ_p> c;
  Vec<ZZ_p> r;
  c.SetLength(this->m+1);
  r = SampleR();
  c = Mult(Transpose(this->A),r);
  c[0]+=m;
  return c;
}

ZZ_p plainlwe::InnerProduct(Vec<ZZ_p> a, Vec<ZZ_p> b){
  int i;
  ZZ_p res = conv<ZZ_p>(conv<ZZ>(0));
  for(i=0;i<a.length();i++){
    res += a[i]*b[i];
  }
  return res;
}

ZZ_p plainlwe::Decrypt(Vec<ZZ_p> c){
  int i;
  ZZ_p res;
  ZZ_p a;
  a = InnerProduct(c,this->s);
  if(conv<ZZ>(a)<this->q/conv<ZZ>(2)){
      res = conv<ZZ_p>(conv<ZZ>(a)%conv<ZZ>(2));
  }
  else{
      res = conv<ZZ_p>(conv<ZZ>(1) - conv<ZZ>(a)%conv<ZZ>(2));
  }
  return res; 
}


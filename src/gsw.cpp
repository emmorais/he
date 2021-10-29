#include <gsw.h>
#include <sample.h>
#include <util.h>
#include <math.h>

using namespace NTL;
using namespace std;

/*
 * Constructor
 */
gsw::gsw(long N, long n, long m, ZZ q, long Q){
  this->N = N;
  this->n = n;
  this->m = m;
  this->q = q;
  this->Q = Q;
}


/* 
 * Key Generation algorithm.
 */ 
void gsw::KeyGen(){	
  int i,j;
  Vec<ZZ_p> e;
  this->sp = sample::SampleRandomArray(this->n);
  this->s.SetLength(this->n+1);
  this->s[0] = conv<ZZ_p>(1);
  for(i=0;i<this->n;i++){
    this->s[i+1] = this->sp[i];
  }
  this->Ap = sample::SampleRandom(this->m, this->n);  
  e.SetLength(this->N);
  e = sample::SampleNoiseArray(this->N);
  this->b = util::Add(util::Mult(this->Ap,this->sp), e);
  this->A.SetDims(this->m, this->n+1);
  for(i=0;i<this->m;i++){
    this->A[i][0] = this->b[i];
    for(j=0;j<this->n;j++){
      this->A[i][j+1] = -this->Ap[i][j];
    }
  }
  ComputeG();
}
    
ZZ_p gsw::Mod2(ZZ a){
  ZZ_p res;
  res = conv<ZZ_p>(a%conv<ZZ>(2));
  return res;
}

/* 
 * Encryption algorithm.
 */ 
Mat<ZZ_p> gsw::Encrypt(ZZ_p m){
  Mat<ZZ_p> c;
  Mat<ZZ_p> r;
  c.SetDims(this->n+1, this->N);
  r = sample::SampleRMatrix(this->m, this->N);
  c = util::Mult(util::Transpose(this->A),r);

  // if m not zero, then add G
  if (m != 0) {
    c = util::Add(c, this->G);
  }
  return c;
}

void gsw::ComputeG(){
  int i,j;
  int power2 = 1;
  int *two;
  two = (int *) calloc(this->Q, sizeof(int));
  for(j=0;j<this->Q;j++){
    two[j] = power2;
    power2 *= 2;
  }
  this->G.SetDims(this->n+1, this->N);
  for(i=0;i<this->n+1;i++){
    for(j=0;j<this->Q;j++){
      this->G[i][i*this->Q+j] = two[j]; 
    }
  }
}

/* 
 * Decryption algorithm.
 */ 
ZZ_p gsw::Decrypt(Mat<ZZ_p> c){
  int i;
  ZZ_p res;
  ZZ_p a;
  Vec<ZZ_p> columnI;
  columnI.SetLength(this->n+1);
  // take I-th column of C
  for (i=0; i<this->n+1; i++) {
    columnI[i] = c[i][this->Q-1];
  }
  a = util::InnerProduct(columnI, this->s);
  ZZ aa = conv<ZZ>(a);
  if(conv<ZZ>(a)>this->q/conv<ZZ>(2)) {
    aa = aa - this->q;
  }
  if(abs(aa)<=this->q/conv<ZZ>(4)) {
      res = conv<ZZ_p>(0);
  }
  else {
      res = conv<ZZ_p>(1);
  }
  return res; 
}


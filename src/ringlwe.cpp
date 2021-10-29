#include <ringlwe.h>
#include <math.h>

using namespace NTL;
using namespace std;

/*
 * Constructor
 */
ringlwe::ringlwe(long N, long d, ZZ q, double delta){
  this->N = N;
  this->d = d;
  this->q = q;
  this->delta = delta;
}

/*
 * SampleMessage is responsible for generating a binary polynomial in ZZ_p[x].
 */
ZZ_pE ringlwe::SampleMessage(){
  int i, mm;
  ZZ_pE res;
  ZZ_pX pol;
  for(i=0;i<this->d;i++){
    mm = rand()%2;
    SetCoeff(pol, i, conv<ZZ_p>(mm));
  }
  res = conv<ZZ_pE>(pol);
  return res;
}
    
/*
 * SampleRandom is used to generate the public key A elements.
 */
Vec<ZZ_pE> ringlwe::SampleRandom(){
  int i;
  Vec<ZZ_pE> res;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    random(res[i]);
  }
  return res;
}

/*
 * SampleNoise is responsible for the Gaussian sampling algorithm.
 * This is a very simple implementation used for our experimentation purposes.
 * In order to a secure implementation of LWE to be used in production, this
 * noise must be generated in a secure way. 
 */
ZZ_pE ringlwe::SampleNoise(){
  ZZ_pE res;
  ZZ_pX pol;
  long i;
  double randNormal, randStdNormal;
  for(i = 0; i < this->d; i++){
    double u1 = ((double) rand()) / (RAND_MAX); 
    double u2 = ((double) rand()) / (RAND_MAX); 
    u2 *= 6.2831;
    u1 = sqrt (-2 * log(u1));
    randNormal = u1 * sin(u2); 
    randNormal = this->delta*8 * randNormal; 
    SetCoeff(pol, i, conv<ZZ_p>((long)randNormal));
    i++;
    if (i <= this->d){
      randNormal = u1 * cos(u2); 
      randNormal = this->delta * randNormal; 
      SetCoeff(pol, i, conv<ZZ_p>((long)randNormal));
    }
  }
  res = conv<ZZ_pE>(pol);
  return res;
}

/*
 * SampleNoiseArray is used to compute the error vector is the LWE
 * encryption scheme.
 */
Vec<ZZ_pE> ringlwe::SampleNoiseArray(){
  Vec<ZZ_pE> res;
  long i;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = SampleNoise();
  }
  return res;
}

/*
 * SampleR is used to generate the randomization factor from the 
 * encryption algorithm.
 */
Vec<ZZ_pE> ringlwe::SampleR(){
  Vec<ZZ_pE> res;
  long i;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = SampleMessage();
  }
  return res;
}

/*
 *
 */
void ringlwe::KeyGen(){	
  int i;
  Vec<ZZ_pE> e;
  this->sp = SampleNoise();
  this->s.SetLength(2);
  this->s[0] = conv<ZZ_pE>(1);
  this->s[1] = this->sp;
  this->Ap = SampleRandom();  
  e.SetLength(this->N);
  e = SampleNoiseArray();
  this->b = this->Ap*this->sp + 2*e;
  this->A.SetDims(this->N, 2);
  for(i=0;i<this->N;i++){
    this->A[i][0] = this->b[i];
    this->A[i][1] = -this->Ap[i];
  }
}
    
/*
 *
 */
ZZ_pE ringlwe::Mod2(ZZ_pE a){
  ZZ_pE res;
  ZZ_pX pol;
  int i;
  pol = conv<ZZ_pX>(a);
  for(i=0;i<this->d;i++){
    ZZ aux = conv<ZZ>(coeff(pol, i));
    aux %= conv<ZZ>(2);
    SetCoeff(pol, i, conv<ZZ_p>(aux));
  } 
  res = conv<ZZ_pE>(pol);
  return res;
}

/*
 *
 */
Mat<ZZ_pE> ringlwe::Transpose(Mat<ZZ_pE> a){
  Mat<ZZ_pE> res;
  int i,j;
  res.SetDims(a.NumCols(), a.NumRows());
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[j][i] = a[i][j];
    }
  }
  return res;
}

/*
 *
 */
Vec<ZZ_pE> ringlwe::Mult(Mat<ZZ_pE> a, Vec<ZZ_pE> b){
  int i,j,k;
  Vec<ZZ_pE> res;
  res.SetLength(2);
  for(i=0;i<a.NumRows();i++){
    res[i] = conv<ZZ_pE>(conv<ZZ>(0));
    for(k=0;k<b.length();k++){
      res[i] += a[i][k]*b[k];
    }
  }
  return res;
}

/*
 *
 */
Vec<ZZ_pE> ringlwe::EncryptR(ZZ_pE m, Vec<ZZ_pE> r){
  Vec<ZZ_pE> c;
  c.SetLength(2);
  c = Mult(Transpose(this->A),r);
  c[0]+=m;
  return c;
}

/*
 *
 */
Vec<ZZ_pE> ringlwe::Encrypt(ZZ_pE m){
  Vec<ZZ_pE> c;
  Vec<ZZ_pE> r;
  c.SetLength(2);
  r = SampleR();
  c = Mult(Transpose(this->A),r);
  c[0]+=m;
  return c;
}

/*
 *
 */
ZZ_pE ringlwe::InnerProduct(Vec<ZZ_pE> a, Vec<ZZ_pE> b){
  int i;
  ZZ_pE res = conv<ZZ_pE>(conv<ZZ>(0));
  for(i=0;i<a.length();i++){
    res += a[i]*b[i];
  }
  return res;
}

/*
 *
 */
ZZ_pE ringlwe::Decrypt(Vec<ZZ_pE> c){
  int i;
  ZZ_pE res;
  ZZ_pX pol;
  ZZ_p a;
  res = InnerProduct(c,this->s);
  pol = conv<ZZ_pX>(res);
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


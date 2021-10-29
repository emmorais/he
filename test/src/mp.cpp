#include <math.h>
#include <stdlib.h>

#include <mp.h>
using namespace NTL;
using namespace std;

/*
 * Constructor
 */
mp::mp(long N, long m, ZZ q, long Q){
  this->N = N;
  this->m = m;
  this->q = q;
  this->Q = Q;
}

////////////////////////////////////////////////////////////////////////////////
//  SAMPLING
////////////////////////////////////////////////////////////////////////////////

Vec<ZZ_p> mp::SampleMessage(){
  int i, mm;
  Vec<ZZ_p> res;
  res.SetLength(this->N*this->Q);
  for(i=0;i<res.length();i++){
    res[i] = rand()%2;
  }
  return res;
}

ZZ_p mp::SampleM(){
  int i, mm;
  ZZ_p res;
  res = rand()%2;
  return res;
}

Mat<ZZ_p> mp::SampleRandom(){
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

ZZ_p mp::SampleNoise(){
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

Vec<ZZ_p> mp::SampleNoiseArray(long n){
  Vec<ZZ_p> res;
  long i;
  res.SetLength(n);
  for(i=0;i<n;i++){
    res[i] = SampleNoise();
  }
  return res;
}

Vec<ZZ_p> mp::SampleS(){
  Vec<ZZ_p> res;
  long i;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = conv<ZZ_p>(SampleM());
  }
  return res;
}

Mat<ZZ_p> mp::SampleR(int n, int m){
  Mat<ZZ_p> res;
  res.SetDims(n, m);
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      res[i][j] = conv<ZZ_p>(D());
    }
  }
  return res;
}

////////////////////////////////////////////////////////////////////////////////
//  UTILS
////////////////////////////////////////////////////////////////////////////////

Vec<ZZ_p> mp::Double(Vec<ZZ_p> x){
  long i;
  Vec<ZZ_p> res;
  res.SetLength(this->N);
  for(i=0;i<this->N;i++){
    res[i] = 2*x[i];
  }
  return res;
}

ZZ_p mp::Mod2(ZZ a){
  ZZ_p res;
  res = conv<ZZ_p>(a%conv<ZZ>(2));
  return res;
}

Mat<ZZ_p> mp::Transpose(Mat<ZZ_p> a){
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

Mat<ZZ> mp::Transpose(Mat<ZZ> a){
  Mat<ZZ> res;
  int i,j;
  res.SetDims(a.NumCols(), a.NumRows());
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[j][i] = a[i][j];
    }
  }
  return res;
}

Vec<ZZ_p> mp::Add(Vec<ZZ_p> a, Vec<ZZ_p> b){
  int i;
  Vec<ZZ_p> res;
  res.SetLength(a.length());
  for(i=0;i<a.length();i++){
    res[i] = a[i] + b[i];
  }
  return res;
}

Vec<ZZ_p> mp::Sub(Vec<ZZ_p> a, Vec<ZZ_p> b){
  int i;
  Vec<ZZ_p> res;
  res.SetLength(a.length());
  for(i=0;i<a.length();i++){
    res[i] = a[i] - b[i];
  }
  return res;
}

Vec<ZZ_p> mp::MultT(Mat<ZZ_p> a, Vec<ZZ_p> b){
  int i,j,k;
  Vec<ZZ_p> res;
  res.SetLength(a.NumCols());
  for(k=0;k<a.NumCols();k++){
    res[k] = conv<ZZ_p>(conv<ZZ>(0));
    for(i=0;i<a.NumRows();i++){
      res[k] += a[i][k]*b[i];
    }
  }
  return res;
}


Vec<ZZ_p> mp::Mult(Mat<ZZ_p> a, Vec<ZZ_p> b){
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

Vec<ZZ_p> mp::Zeros(int n){
  Vec<ZZ_p> res;
  res.SetLength(n);
  int i;
  for(i=0;i<n;i++){
    res[i] = conv<ZZ_p>(0);
  }
  return res;
}

Vec<ZZ_p> mp::ScalarMul(Vec<ZZ_p> x, ZZ_p c){
  Vec<ZZ_p> res;
  res.SetLength(x.length());
  int i;
  for(i=0;i<x.length();i++){
    res[i] = x[i]*c;
  }
  return res;
}

Vec<ZZ_p> mp::ConvVecToZZp(Vec<ZZ> x){
  Vec<ZZ_p> res;
  res.SetLength(x.length());
  int i;
  for(i=0;i<x.length();i++){
    res[i] = conv<ZZ_p>(x[i]);
  }
  return res;
}

Vec<ZZ> mp::ConvVecToZZ(Vec<ZZ_p> x){
  Vec<ZZ> res;
  res.SetLength(x.length());
  int i;
  for(i=0;i<x.length();i++){
    res[i] = conv<ZZ>(x[i]);
  }
  return res;
}

int mp::D() {
  int r = rand() % 4;
  if (r % 2  == 0)
    return 0;
  else {
    if (r/2 % 2 == 0)
      return -1;
    else 
      return 1;
  }
}

Mat<ZZ> mp::SetIDZZ(int n){
  int i,j;
  Mat<ZZ> res;
  res.SetDims(n, n);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (i==j)
        res[i][j] = 1;
      else
        res[i][j] = 0;
    }
  }
  return res;
}

Mat<ZZ_p> mp::SetID(int n){
  int i,j;
  Mat<ZZ_p> res;
  res.SetDims(n, n);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (i==j)
        res[i][j] = 1;
      else
        res[i][j] = 0;
    }
  }
  return res;
}

Vec<ZZ_p> mp::concat(Vec<ZZ_p> a, Vec<ZZ_p> b){
  Vec<ZZ_p> res;
  res.SetLength(a.length()+b.length());
  int i;
  for(i=0;i<a.length();i++){
    res[i] = a[i];
  }
  for(i=a.length();i<a.length()+b.length();i++){
    res[i] = b[i-a.length()];
  }
  return res;
}

Mat<ZZ_p> mp::concat(Mat<ZZ_p> a, Mat<ZZ_p> b){
  Mat<ZZ_p> res;
  res.SetDims(a.NumRows(), a.NumCols()+b.NumCols());
  int i,j;
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[i][j] = a[i][j];
    }
    for(j=a.NumCols();j<a.NumCols()+b.NumCols();j++){
      res[i][j] = b[i][j-a.NumCols()];
    }
  }
  return res;
}

Mat<ZZ_p> mp::concatv(Mat<ZZ_p> a, Mat<ZZ_p> b){
  Mat<ZZ_p> res;
  res.SetDims(a.NumRows()+b.NumRows(), a.NumCols());
  int i,j;
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[i][j] = a[i][j];
    }
  }
  for(i=a.NumRows();i<a.NumRows()+b.NumRows();i++){
    for(j=0;j<b.NumCols();j++){
      res[i][j] = b[i-a.NumRows()][j];
    }
  }
  return res;
}

Mat<ZZ_p> mp::mul(Mat<ZZ_p> a, Mat<ZZ_p> b){
  int i,j,k;
  Mat<ZZ_p> res;
  res.SetDims(a.NumRows(), b.NumCols());
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<b.NumCols();j++){
      res[i][j] = 0;
      for(k=0;k<a.NumCols();k++){
        res[i][j] += a[i][k]*b[k][j]; 
      }
    }
  }
  return res;
}

Vec<ZZ> mp::mulv(Mat<ZZ> a, Vec<ZZ> b){
  int i,j;
  Vec<ZZ> res;
  res.SetLength(a.NumRows());
  for(i=0;i<a.NumRows();i++){
    res[i] = 0;
    for(j=0;j<a.NumCols();j++){
      res[i] += a[i][j]*b[j]; 
    }
  }
  return res;
}


Mat<ZZ_p> mp::add(Mat<ZZ_p> a, Mat<ZZ_p> b){
  int i,j;
  Mat<ZZ_p> res;
  res.SetDims(a.NumRows(), a.NumCols());
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[i][j] = a[i][j] + b[i][j]; 
    }
  }
  return res;
}

Mat<ZZ_p> mp::sub(Mat<ZZ_p> a, Mat<ZZ_p> b){
  int i,j;
  Mat<ZZ_p> res;
  res.SetDims(a.NumRows(), a.NumCols());
  for(i=0;i<a.NumRows();i++){
    for(j=0;j<a.NumCols();j++){
      res[i][j] = a[i][j] - b[i][j]; 
    }
  }
  return res;
}

////////////////////////////////////////////////////////////////////////////////
//  SCHEME
////////////////////////////////////////////////////////////////////////////////

void mp::ComputeRvID(){
  Mat<ZZ_p> ID = this->SetID(this->N*this->Q);
  this->RvID.SetDims(this->R.NumRows()+this->N*this->Q, this->R.NumCols());
  this->RvID = concatv(this->R, ID);
}
    
void mp::ComputeG(){
  int i,j;
  int power2 = 1;
  int *two;
  two = (int *) calloc(this->Q, sizeof(int));
  for(j=0;j<this->Q;j++){
    two[j] = power2;
    power2 *= 2;
  }
  this->G.SetDims(this->N, this->N*this->Q);
  for(i=0;i<this->N;i++){
    for(j=0;j<this->Q;j++){
      this->G[i][i*this->Q+j] = two[j]; 
    }
  }
}

Vec<ZZ_p> mp::g(Vec<ZZ_p> s, Vec<ZZ_p> e){
  Vec<ZZ_p> b;
  b.SetLength(e.length());
  b = MultT(this->A,s);
  b = Add(b, e);
  return b;
}

Vec<ZZ_p> mp::EncryptR(Vec<ZZ_p> m, Vec<ZZ_p> s, Vec<ZZ_p> e){
  Vec<ZZ_p> b;
  b.SetLength(e.length());
  b = Mult(Transpose(this->A),s);
  // change mod to 2q
  ZZ_p::init(conv<ZZ>(this->q)*conv<ZZ>(2)); 
  b = ScalarMul(b, conv<ZZ_p>(2));
  b = Add(b, e);
  b = Add(b, concat(Zeros((2*this->N*this->Q)-this->N*this->Q), ConvVecToZZp(Encode(ConvVecToZZ(m)))));
  // change mod back to q
  ZZ_p::init(conv<ZZ>(this->q)); 
  return b;
}

Vec<ZZ_p> mp::Encrypt(Vec<ZZ_p> m){
  Vec<ZZ_p> s;
  s.SetLength(this->N);
  s = SampleS();
  Vec<ZZ_p> e;
  e.SetLength(this->N*this->Q*2);
  e = SampleNoiseArray(e.length());
  Vec<ZZ_p> b;
  b.SetLength(e.length());
  b = Mult(Transpose(this->A),s);
  // change mod to 2q
  ZZ_p::init(conv<ZZ>(this->q)*conv<ZZ>(2)); 
  b = ScalarMul(b, conv<ZZ_p>(2));
  b = Add(b, e);
  b = Add(b, concat(Zeros(this->N*this->Q), ConvVecToZZp(Encode(ConvVecToZZ(m)))));
  // change mod back to q
  ZZ_p::init(conv<ZZ>(this->q)); 
  return b;
}

Vec<ZZ_p> mp::DecryptR(Vec<ZZ_p> b, Vec<ZZ_p> e){
  int i;
  Vec<ZZ_p> res, v, vt;
  v.SetLength(2*this->N*this->Q);
  vt.SetLength(this->N);
  // change mod to 2q
  ZZ_p::init(conv<ZZ>(this->q)*conv<ZZ>(2)); 
  v = b - e;
  vt = MultT(this->RvID, v);
  res= ConvVecToZZp(Decode(ConvVecToZZ(vt)));
  // change mod back to q
  ZZ_p::init(conv<ZZ>(this->q)); 
  return res; 
}

Vec<ZZ_p> mp::Decrypt(Vec<ZZ_p> b){
  Vec<ZZ_p> bb, bbb;
  bb.SetLength(b.length());
  bbb.SetLength(b.length());
  int i;
  for(i=0;i<b.length();i++){
    bb[i] = conv<ZZ_p>(conv<ZZ>(b[i]));
  }
  Vec<ZZ_p> res, v, vt, et, eet;
  res.SetLength(this->N*this->Q);
  v.SetLength(2*this->N*this->Q);
  vt.SetLength(this->N);

  Vec<ZZ_p> e, s, ee;
  e.SetLength(b.length());
  ee.SetLength(e.length());
  s.SetLength(this->N);
  InvGA(s, e, bb);
  et = MultT(this->RvID, e);
  // change mod to 2q
  ZZ_p::init(conv<ZZ>(this->q)*conv<ZZ>(2)); 
  for(i=0;i<b.length();i++){
    if (conv<ZZ>(e[i]) > q/2) {
      ee[i] = conv<ZZ_p>(conv<ZZ>(e[i]) - q);
    } else {
      ee[i] = e[i];
    }
  }
  eet = MultT(this->RvID, ee);
  v = b - ee;
  vt = MultT(this->RvID, v);
  // Fix some positions
  for(i=0;i<eet.length();i++){
    if ((conv<ZZ>(eet[i]) > this->q/2) && (conv<ZZ>(eet[i]) < (2*this->q - q/2))) {
      vt[i] = conv<ZZ_p>(conv<ZZ>(vt[i]) - this->q);
    } 
  }
  res = ConvVecToZZp(Decode(ConvVecToZZ(vt)));
  // change mod back to q
  ZZ_p::init(conv<ZZ>(this->q)); 
  return res; 
}


void mp::KeyGen(){	
  int i,j;
  Vec<ZZ_p> e;
 
  this->Ap = SampleRandom();  
  ComputeS();
  ComputeSS();
  ComputeG();
  this->H = SetID(this->N);
  GenTrap(this->Ap, this->H);
  ComputeRvID();
}
    

void mp::GenTrap(Mat<ZZ_p> Ab, Mat<ZZ_p> H){
  long nq = this->N*this->Q;
  this->A.SetDims(this->N, nq+this->m);
  Mat<ZZ_p> suf;
  suf.SetDims(this->N, nq);
  Mat<ZZ_p> abr;
  abr.SetDims(this->N, nq);
  
  // Choose R from distribution D: -1 -> 1/4; 0 -> 1/2; 1 -> 1/2
  this->R = SampleR(nq, nq);
  this->R = SetID(nq);
  cout << "R: " << this->R << "\n";
  
  // A = [Ab | HG - AbR]
  suf = mul(H, this->G);
  abr = mul(Ab, this->R);
  suf = sub(suf, abr);
  this->A = concat(Ab, suf);
  cout << "A: " << this->A << "\n";
  return; 
}
    

/*
 * Recover the components of s and e individually, given the corresponding part of b.
 */ 
void mp::InvGOracle(ZZ_p &s, Vec<ZZ_p> &e, Vec<ZZ_p> b){
  int i,k,condition, p2, pi2;
  ZZ_p res;

  k = b.length();
  e.SetLength(k);
  p2 = 1;
  pi2 = pow(2, k-1);
  s = 0;
  for(i=k-1;i>=0;i--){
    condition = (abs(conv<ZZ>(b[i] - conv<ZZ_p>(pi2)*s)) < q/4);
    s = s + p2*condition;
    e[i] = b[i] - conv<ZZ_p>(pi2*s);

    pi2 /= 2;
    p2 *= 2;
  }

  return;
}

/*
 * This method is used to recover the whole vectors s and e, and not only its components. 
 * The components can be recovered in parallel using InvOracle method. 
 */
void mp::InvGA(Vec<ZZ_p> &s, Vec<ZZ_p> &e, Vec<ZZ_p> b){
  Vec<ZZ_p> ee;
  ZZ val;
  Vec<ZZ_p> bt = MultT(this->RvID, b);
   
  int i,j,condition, p2, pi2;

  ee.SetLength(b.length());
  s.SetLength(this->N);
  for(j=0;j<this->N;j++){
    s[j] = 0;
    p2 = 1;
    pi2 = pow(2, this->Q-1);
    for(i=this->Q-1;i>=0;i--){
      val = conv<ZZ>(bt[j*this->Q + i] - conv<ZZ_p>(pi2)*s[j]);
      if (val > q/2){
        val -= q; 
      }
      condition = (abs(val) > q/4);
      s[j] = s[j] + p2*condition;
      ee[j*this->Q + i] = bt[j*this->Q + i] - conv<ZZ_p>(pi2*s[j]);
  
      pi2 /= 2;
      p2 *= 2;
    }
  }

  // adjust since modulus is prime (not a power of 2)
  /*for(i=0;i<this->N;i++){
    if(conv<ZZ>(s[i]) > q/2){
      s[i] = s[i] - p2;
    }
  }*/

  // e = b - A^t*s 
  Vec<ZZ_p> Ats = MultT(this->A, s);
  e = Sub(b, Ats);
   
}


void mp::ComputeSS(){
  int i,j,k;
  ZZ det;
  this->SS.SetDims(this->N*this->Q, this->N*this->Q);
  this->tSS.SetDims(this->N*this->Q, this->N*this->Q);
  for(i=0;i<this->Q;i++){
    for(j=0;j<this->Q;j++){
      for(k=0;k<this->N;k++){
        this->SS[i+k*this->Q][j+k*this->Q] = PowerMod(conv<ZZ>(2), conv<ZZ>(i+j), this->q);
        this->tSS[i+k*this->Q][j+k*this->Q] = 2*this->SS[i][j];
      }
    }
  }
  inv(det, this->invSS, this->SS);
  inv(det, this->tinvSS, this->tSS);
  // In some cases it is necessary to adjust the matrix, like Q=2
  /*for(i=0;i<this->Q;i++){
    for(j=0;j<this->Q;j++){
      this->tinvSS[i][j] *= -1;
    }
  }*/
}

void mp::ComputeS(){
  int i,j,k;
  ZZ_p d;
  long nq = this->N*this->Q;

  long qbits;
  ZZ qnext;
  ZZ det;
  // compute S
  this->S.SetDims(nq,nq); 
  for(i=0;i<nq;i++){
    if(i%this->Q != this->Q-1){
      this->S[i][i] = 2;
      this->S[i+1][i] = -1;
    }
    else {
      this->S[i][i] = 2;
    }
  }
  // compute invS
  inv(det, this->invS, this->S);
}

Vec<ZZ> mp::Encode(Vec<ZZ> m){
  int i;
  Vec<ZZ> res;
  res = mulv(this->SS, m);
  return res;
}

Vec<ZZ> mp::Decode(Vec<ZZ> m){
  Vec<ZZ> res, red, prod;
  int i;
  // Reduce to find coset
  prod = mulv(this->tinvSS, m);
  for(i=0;i<prod.length();i++){
    div(prod[i], conv<ZZ>(prod[i]), power(this->q, this->N*this->Q));
  }
  prod = mulv(this->tSS, prod);
  red = m - prod;
  
  // Decode
  res = mulv(this->invSS, red);
  for(i=0;i<m.length();i++){
    div(res[i], conv<ZZ>(res[i]), power(this->q, this->N*(this->Q-1)));
  }
  return res;
}	

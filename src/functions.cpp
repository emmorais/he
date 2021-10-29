#include "functions.h"

#define NUM 100
#define N 1024
#define Q 157
#define T 1024
#define W 2
#define E 0
#define LWQ 158

ZZ_pE Floor(ZZ_pE x, ltv *scheme){
  ZZ_pE res;
  ZZ_pX p10;
  ZZ_pE r10;  
  ZZ_pE half = Encoding::EncodeRZZ(conv<RR>(power(conv<RR>(0.5), 16)), scheme);
  ZZ_pE p2e = Encoding::EncodeRZZ(conv<RR>(power(conv<RR>(2), 16)), scheme);
  ZZ_pE chalf = scheme->Encrypt(half);
  ZZ_pE cp2e = scheme->Encrypt(p2e);
  res = scheme->Mul(cp2e, x, scheme->b->q);
  res = scheme->Mul(res, chalf, scheme->b->q);
  return res;
}


void lheMaxTest(){
  int count=0, ok=0, nok=0;
  long lwq = Q+1;
  double delta = 1;
  long n = N;
  long t = T;
  long w = W;
  long e = E;
  ltv *fnord;
  ZZ q;
  GenPrime(q, Q);
  cout << "q: " << q << "\n";
  ZZ_p::init(conv<ZZ>(q));
  ZZ_pX P, pol;
  SetCoeff(P, 0, 1);
  SetCoeff(P, n, 1);
  ZZ_pE::init(P); 
  cout << "Leveled homomorphic encryption test:\n";
  fnord = new ltv();
  fnord->ParamsGen(t, w, e, n, q, delta, lwq);
  while(count < NUM){
    fnord->KeyGen();
    ZZ_pE m1 = Encoding::EncodeRZZ(conv<RR>(47), fnord);
    ZZ_pE mm1 = Encoding::EncodeRZZ(conv<RR>(-47), fnord);
    ZZ_pE m2 = Encoding::EncodeRZZ(conv<RR>(4), fnord);
    ZZ_pE mm2 = Encoding::EncodeRZZ(conv<RR>(-4), fnord);
    ZZ_pE m3 = Encoding::EncodeRZZ(conv<RR>(1), fnord);
    ZZ_pE c1 = fnord->Encrypt(m1);
    ZZ_pE mc1 = fnord->Encrypt(mm1);
    ZZ_pE c2 = fnord->Encrypt(m2);
    ZZ_pE mc2 = fnord->Encrypt(mm2);
    //ZZ_pE c3 = fnord->Encrypt(m3);
    //ZZ_pE cmul = fnord->Mul(c1,c3,q);
    ZZ_pE cmax = Functions::Max(c1, c2, fnord);
    ZZ_pE mm = fnord->Decrypt(cmax);
    cout << "mm: " << mm << "\n";
    cout << "Decode(mm): " << Encoding::DecodeRZZ(mm, fnord) << "\n";
    if (fnord->ModN(m1*m1, t) == mm){
      ok++; 
      cout << "ok!\n";
    }
    else {
      nok++; 
      cout << "nok!\n";
    }
    count++;
  }
  cout << "ok,nok: " << ok << "," << nok << "\n";
}



ZZ_pE Functions::RotLeft(ZZ_pE a, int k, ltv *scheme){
    ZZ_pE res;
    ZZ_pX pol;
    int i;
    for(i=0;i<scheme->b->n;i++){
        if (i == k)
            SetCoeff(pol, i, 1);
        else 
            SetCoeff(pol, i, 0);
    }
    res = conv<ZZ_pE>(pol);
    res = inv(res);
    res = scheme->Encrypt(res);
    res = scheme->Mul(res,a,scheme->b->q);
    return res;
}

ZZ_pE Functions::RotRight(ZZ_pE a, int k, ltv *scheme){
    ZZ_pE res;
    ZZ_pX pol;
    int i;
    for(i=0;i<scheme->b->n;i++){
        if (i == k)
            SetCoeff(pol, i, 1);
        else 
            SetCoeff(pol, i, 0);
    }
    res = conv<ZZ_pE>(pol);
    res = scheme->Encrypt(res);
    res = scheme->Mul(res,a,scheme->b->q);
    return res;
}

ZZ_pE Functions::Equal(ZZ_pE a, ZZ_pE b, ltv *scheme){
    ZZ_pE res;
    return res;
}

ZZ_pE Functions::GeqZero(ZZ_pE a, ltv *scheme){
    ZZ_pE res;
    ZZ_pE um = Encoding::EncodeRZZ(conv<RR>(1), scheme);
    ZZ_pE m2 = Encoding::EncodeRZZ(conv<RR>(pow(2, 31)), scheme);
    m2 = inv(m2);
    res = scheme->Encrypt(m2);
    um = scheme->Encrypt(um);
    res = scheme->Mul(res,a,scheme->b->q);
    res = um-res;
    return res;
}

ZZ_pE Functions::LeqZero(ZZ_pE a, ltv *scheme){
    ZZ_pE res;
    ZZ_pE mum = Encoding::EncodeRZZ(conv<RR>(1), scheme);
    mum = scheme->Encrypt(mum);
    res = a - mum;
    ZZ_pE m = scheme->Decrypt(res);
    //cout << "a + mum: " << m << "\n";
    //cout << "dec(a + mum): " << Encoding::DecodeRZZ(m, scheme) << "\n";
    res = LZero(res,scheme);
    res = RotLeft(res, scheme->b->n-31, scheme);
    res = GeqZero(res, scheme);
    //cout << "LeqZero(a): " << scheme->Decrypt(res) << "\n";
    //cout << "Decode(LeqZero(a)): " << Encoding::DecodeRZZ(scheme->Decrypt(res), scheme) << "\n";
    res = res - mum;
    //cout << "LeqZero(a)-1: " << scheme->Decrypt(res) << "\n";
    //cout << "Decode(LeqZero(a)-1): " << Encoding::DecodeRZZ(scheme->Decrypt(res), scheme) << "\n";
    return res;
}

ZZ_pE Functions::GZero(ZZ_pE a, ltv *scheme){
    ZZ_pE res;
    ZZ_pE um = Encoding::EncodeRZZ(conv<RR>(1), scheme);
    um = scheme->Encrypt(um);
    res = LeqZero(a,scheme);
    return um - res;
}

ZZ_pE Functions::LZero(ZZ_pE a, ltv *scheme){
    ZZ_pE res;
    ZZ_pE m2 = Encoding::EncodeRZZ(conv<RR>(pow(2, 31)), scheme);
    m2 = inv(m2);
    res = scheme->Encrypt(m2);
    res = scheme->Mul(res,a,scheme->b->q);
    return res;
}
        
ZZ_pE Functions::Max(ZZ_pE a, ZZ_pE b, ltv *scheme){
    ZZ_pE res;
    ZZ_pE aux = GeqZero(a-b, scheme);
    ZZ_pE m = scheme->Decrypt(aux);
    cout << "GeqZero(a-b): " << Encoding::DecodeRZZ(m, scheme) << "\n";
    res = scheme->Mul(a, GeqZero(a-b, scheme), scheme->b->q); 
    m = scheme->Decrypt(b-a);
    cout << "Dec(b-a): " << m << "\n";
    cout << "b-a: " << Encoding::DecodeRZZ(m, scheme) << "\n";
    aux = GeqZero(b-a, scheme);
    m = scheme->Decrypt(aux);
    cout << "GeqZero(b-a): " << Encoding::DecodeRZZ(m, scheme) << "\n";
    ZZ_pE aux2 = Encoding::EncodeRZZ(conv<RR>(pow(2, 33)-1), scheme);
    aux2 += conv<ZZ_pE>(1);
    cout << "Encode(pow(2, 33)-1): " << aux2 << "\n";
    aux2 = scheme->Encrypt(aux2);

    cout << "GeqZero(b-a+aux2): " << Encoding::DecodeRZZ(scheme->Decrypt(GeqZero(b-a+aux2, scheme)), scheme) << "\n";

    res += scheme->Mul(b, GeqZero(b-a+aux2, scheme), scheme->b->q); 
    return res;
}

ZZ_pE Functions::Max(ZZ_pE a, ZZ_pE ma, ZZ_pE b, ZZ_pE mb, ltv *scheme){
    ZZ_pE res;
    //ZZ_pE aux = GeqZero(a-b, scheme);
    //ZZ_pE m = scheme->Decrypt(aux);
    //cout << "GeqZero(a-b): " << Encoding::DecodeRZZ(m, scheme) << "\n";
    cout << "a: " << Encoding::DecodeRZZ(scheme->Decrypt(a), scheme) << "\n";
    cout << "Dec(a-b): " << scheme->Decrypt(a-b) << "\n";
    cout << "a-b: " << Encoding::DecodeRZZ(scheme->Decrypt(a-b), scheme) << "\n";
    cout << "geq(a-b): " << Encoding::DecodeRZZ(scheme->Decrypt(GeqZero(a-b, scheme)), scheme) << "\n";
    res = scheme->Mul(a, GeqZero(a-b, scheme), scheme->b->q); 
    cout << "res: " << Encoding::DecodeRZZ(scheme->Decrypt(res), scheme) << "\n";

    ZZ_pE x = Encoding::EncodeRZZ(conv<RR>(2), scheme);
    x = scheme->Encrypt(x);
    ZZ_pE aa = scheme->Decrypt(scheme->Mul(a, x, scheme->b->q));
    cout << "(2a): " << scheme->Decrypt(aa) << "\n";
    cout << "(2a-a): " << scheme->Decrypt(aa+ma) << "\n";
    aa = aa+ma;
    cout << "GeqZero((2a)-a): " << Encoding::DecodeRZZ(scheme->Decrypt(GeqZero(aa, scheme)), scheme) << "\n";
    ZZ_pE bb = scheme->Decrypt(scheme->Mul(b, x, scheme->b->q));
    cout << "(2b): " << scheme->Decrypt(bb) << "\n";
    cout << "(2b-b): " << scheme->Decrypt(bb+mb) << "\n";
    bb = bb+mb;
    cout << "GeqZero((2b)-b): " << Encoding::DecodeRZZ(scheme->Decrypt(GeqZero(bb, scheme)), scheme) << "\n";
    cout << "((2a)-a)-((2b)-b): " << scheme->Decrypt((aa)-(bb)) << "\n";

    //cout << "GeqZero(a+mb): " << Encoding::DecodeRZZ(scheme->Decrypt(GeqZero(a+mb, scheme)), scheme) << "\n";
    //cout << "b+ma: " << scheme->Decrypt(b+ma) << "\n";
    //cout << "GeqZero(b+ma): " << Encoding::DecodeRZZ(scheme->Decrypt(GeqZero(b+ma, scheme)), scheme) << "\n";
    //ZZ_pE m = scheme->Decrypt(b-mb+ma-a);
    //cout << "Dec(2b-2a): " << m << "\n";
    //cout << "2(b-a): " << Encoding::DecodeRZZ(m, scheme) << "\n";

    cout << "(a-b)+(a+mb): " << scheme->Decrypt((a-b)+(a+mb)) << "\n";
    cout << "Decode((a-b)+(a+mb)): " << Encoding::DecodeRZZ(scheme->Decrypt((a-b)+(a+mb)), scheme) << "\n";
    cout << "(b-a)+(b+ma): " << scheme->Decrypt((b-a)+(b+ma)) << "\n";
    cout << "Decode((b-a)+(b+ma)): " << Encoding::DecodeRZZ(scheme->Decrypt((b-a)+(b+ma)), scheme) << "\n";
    
    ZZ_pE aux = GeqZero(b-a, scheme) + GeqZero(b+ma, scheme);
    ZZ_pE m = scheme->Decrypt(aux);
    cout << "aux: " << m << "\n";
    cout << "GeqZero(b-a)+GeqZero(b+ma): " << Encoding::DecodeRZZ(m, scheme) << "\n";
    //ZZ_pE aux2 = Encoding::EncodeRZZ(conv<RR>(pow(2, 33)-1), scheme);
    //aux2 += conv<ZZ_pE>(1);
    //cout << "Encode(pow(2, 33)-1): " << aux2 << "\n";
    //aux2 = scheme->Encrypt(aux2);

    cout << "b: " << Encoding::DecodeRZZ(scheme->Decrypt(b), scheme) << "\n";
    ZZ_pE r = scheme->Mul(aux, b, scheme->b->q); 
    cout << "r: " << scheme->Decrypt(r) << "\n";
    cout << "Decode(r): " << Encoding::DecodeRZZ(scheme->Decrypt(r), scheme) << "\n";
    res += r;

    return res;
}

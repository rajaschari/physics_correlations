#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <bits/stdc++.h>


using namespace std;

const double pi = M_PI;

//level 0--vectors
double dot(vector<double> a,vector<double> b){
    return a[0]*b[0]+a[1]*b[1];
}

vector<double> sumVec(vector<double> a,vector<double> b){


    vector<double> c=a; c[0]+=b[0];c[1]+=b[1];
    return c;
};

vector<double> negVec(vector<double> a){

    vector<double> b(2); b[0]=-a[0];b[1]=-a[1];

    return b;

}


vector<double> zeroVec (2,0.0);

double n1p[] = {1.0/2.0,sqrt(3.0)/2.0};
double n2p[] = {-1.0/2.0,sqrt(3.0)/2.0};

double b1p[] = {2.0*pi,2.0*pi/sqrt(3.0)};
double b2p[] = {-2.0*pi,2.0*pi/sqrt(3.0)};

//double rxp[] = {1.0,1.0/sqrt(3.0)};
//double ryp[] = {-1.0,1.0/sqrt(3.0)};
//double rzp[] = {0.0,-2.0/sqrt(3.0)};
//
//double rxyp[] = {-3.0/2.0,-1.0/(2.0*sqrt(3.0))};
//double rxzp[] = {-1.0,-2.0/sqrt(3.0)};
//double ryxp[] = {-3.0/2.0,-1.0/(2.0*sqrt(3.0))};
//double ryzp[] = {1.0,-2.0/sqrt(3.0)};
//double rzxp[] = {1.0/2.0,5.0/2*sqrt(3.0)};
//double rzyp[] = {-1.0/2.0,5.0/2*sqrt(3.0)};

vector<double> n1 (n1p, n1p + sizeof(n1p) / sizeof(double) );
vector<double> n2 (n2p, n2p + sizeof(n2p) / sizeof(double) );

vector<double> b1 (b1p, b1p + sizeof(b1p) / sizeof(double) );
vector<double> b2 (b2p, b2p + sizeof(b2p) / sizeof(double) );

//vector<double> rx (rxp, rxp + sizeof(rxp) / sizeof(double) );
//vector<double> ry (ryp, ryp + sizeof(ryp) / sizeof(double) );
//vector<double> rz (rzp, rzp + sizeof(rzp) / sizeof(double) );
//
//vector<double> rxy (rxyp, rxyp + sizeof(rxyp) / sizeof(double) );
//vector<double> rxz (rxzp, rxzp + sizeof(rxzp) / sizeof(double) );
//vector<double> ryx (ryxp, ryxp + sizeof(ryxp) / sizeof(double) );
//vector<double> ryz (ryzp, ryzp + sizeof(ryzp) / sizeof(double) );
//vector<double> rzx (rzxp, rzxp + sizeof(rzxp) / sizeof(double) );
//vector<double> rzy (rzyp, rzyp + sizeof(rzyp) / sizeof(double) );

vector<double> ex = n2;
vector<double> ey = negVec(n1);
vector<double> ez = sumVec(n1,negVec(n2));
//
//ex = sumVec(n2,zeroVec,ex);
//ey = sumVec(negVec(n1),zeroVec,ey);
//ez = sumVec(n1,negVec(n2),ez);

double Ep[] = {1.0/sqrt(6.0),-2.0/sqrt(6.0),1.0/sqrt(6.0)};
double Bp[] = {1.0/sqrt(6.0),1.0/sqrt(6.0),-2.0/sqrt(6.0)};
//level 1-- field declarations

//level 3-- real scalar energy scales 1

double Rxx(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ex + (Ey + Ez)/sqrt(2.0) ;
}

double Rxy(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ey+ (Ex)/sqrt(2.0) ;
}

double Rxz(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ez+ (Ex)/sqrt(2.0) ;
}

double Ryx(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ex+ (Ey)/sqrt(2.0) ;
}

double Ryy(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ey+ (Ex + Ez)/sqrt(2.0) ;
}

double Ryz(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ez+ (Ey)/sqrt(2.0) ;
}

double Rzx(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ex+ (Ez)/sqrt(2.0) ;
}

double Rzy(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ey+ (Ez)/sqrt(2.0) ;
}

double Rzz(double t){
    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  Ez+ (Ex+Ey)/sqrt(2.0) ;
}


double tx(double t){

    return  (2.0/0.27)*Ryx(t)*Rzx(t) ;
}

double ty(double t){

    return  (2.0/0.27)*Rxy(t)*Rzy(t) ;
}

double tz(double t){

    return  (2.0/0.27)*Rxz(t)*Ryz(t) ;
}


double Jx(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    return  1.0+pow(Bx,2)/0.27+pow(Rxx(t),2)/0.47 ;
}

double Jy(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    return  1.0+pow(By,2)/0.27+pow(Ryy(t),2)/0.47 ;
}

double Jz(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    return  1.0+pow(Bz,2)/0.27+pow(Rzz(t),2)/0.47 ;
}


double K(double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    return (1.0/(0.27*0.27))*Bx*By*Bz;
}

//2*Bx*By*Bz/((double)pow(0.27,2.0))
double Qx(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  (2.0/0.27)*(Bz*(Ryz(t))-By*(Rzy(t))) ;
}

double Qy(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  (2.0/0.27)*(Bx*(Rzx(t))-Bz*(Rxz(t)));
}

double Qz(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  (2.0/0.27)*(By*(Rxy(t)) - Bx*(Ryx(t)));
}


//level 3-- complex scalar hamiltonian elements
//complex<double> z2 = exp((complex<double>)1i*pi); // imaginary unit squared

complex<double> f(vector<double> k,double t, double r){
//    double kx=k[0];
//    double ky=k[1];

    return  2.0*(double(Jx(t,r))*exp((complex<double>) 1i*dot(negVec(k),n1)) + double(Jy(t,r))*exp((complex<double>) 1i*dot(negVec(k),n2)) + double(Jz(t,r)));
}

double Delta(vector<double> k, double r){
//    double kx=k[0];
//    double ky=k[1];
    vector<double> n3 = sumVec(n2,negVec(n1));
    return  4.0*(double(K(r)))*(sin(dot(k,n1))+sin(dot(negVec(k),n2))+sin(dot(k,n3)));
}

double G(vector<double> k, double t, double r){
//    double kx=k[0];
//    double ky=k[1];

    return  4.0*((double(Qx(t,r)))*sin(dot(k,ex))+(double(Qy(t,r)))*sin(dot(k,ey))+(double(Qz(t,r)))*sin(dot(k,ez)));
}

complex<double> T(vector<double> k, double t){
//    double kx=k[0];
//    double ky=k[1];

    vector<double> rx= sumVec(n1,negVec(n2));
    vector<double> ry= sumVec(n2,negVec(n1));
    vector<double> rz= negVec(sumVec(n1,n2));

    return  (double(tx(t)))*exp((complex<double>) 1i*dot(k,rx))+(double(ty(t)))*exp((complex<double>) 1i*dot(k,ry))+(double(tz(t)))*exp((complex<double>) 1i*dot(k,rz));
}

complex<double> P(vector<double> k, double t){
//    double kx=k[0];
//    double ky=k[1];

    vector<double> rxy= sumVec(sumVec(n2,negVec(n1)),negVec(n1));
    vector<double> rxz= sumVec(negVec(n1),negVec(n1));
    vector<double> ryx= sumVec(sumVec(n1,negVec(n2)),negVec(n2));
    vector<double> ryz= sumVec(negVec(n2),negVec(n2));
    vector<double> rzx=n1;
    vector<double> rzy=n2;


    return  ((double(pow(Rxy(t),2)))*exp((complex<double>) 1i*dot(k,rxy))+(double(pow(Rxz(t),2)))*exp((complex<double>) 1i*dot(k,rxz))+
    (double(pow(Ryx(t),2)))*exp((complex<double>) 1i*dot(k,ryx))+(double(pow(Ryz(t),2)))*exp((complex<double>) 1i*dot(k,ryz))+
    (double(pow(Rzx(t),2)))*exp((complex<double>) 1i*dot(k,rzx))+(double(pow(Rzy(t),2)))*exp((complex<double>) 1i*dot(k,rzy)))/0.27;
}


complex<double> f1(vector<double> k,double t,double r){

    return f(k,t,r)+P(k,t)+T(k,t);
}

complex<double> f1e(vector<double> k,double t,double r){

    return (f1(k,t,r)+f1(negVec(k),t,r))/2.0;
}

complex<double> f1o(vector<double> k,double t,double r){

    return (f(k,t,r)-f(negVec(k),t,r))/2.0;
}

//level 4--energy

double e1(vector<double> k, double t, double r){
    const complex<double> X = f(k,t,r)+P(k,t)+T(k,t);
    const complex<double> Y = Delta(k,r);
  return real(G(k,t,r)-sqrt(X*conj(X)+Y*conj(Y)));
}

double e2(vector<double> k, double t, double r){
    const complex<double> X = f(k,t,r)+P(k,t)+T(k,t);
    const complex<double> Y = Delta(k,r);
  return real(G(k,t,r)+sqrt(X*conj(X)+Y*conj(Y)));
}

//level 5 -- kernel and Cp

double kernel(vector<double> k,double t,double r,double TT,double e11,double e21){

    return  (
    pow(e11/TT,2.0)*exp(e11/TT)/pow((1+exp(e11/TT)),2)+
    pow(e21/TT,2.0)*exp(e21/TT)/pow((1+exp(e21/TT)),2));
}

int main(int argc,char** argv) {

  const int L = atoi(argv[1]);
  const double t = atof(argv[2]);
  const double r = atof(argv[3]);
  const double Tmin=atof(argv[4]);
  const double Tmax=atof(argv[5]);

  double en1[L][L];
  double en2[L][L];
  vector<double> ks[L][L];

  double De[L][L];
  double Ge[L][L];

  for (int i1 = 0; i1 < L; i1++) {
    for (int i2 = 0; i2 < L; i2++) {
      vector<double> k(2);
      k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
      k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;
      en1[i1][i2] = e1(k,t,r);
      en2[i1][i2] = e2(k,t,r);
//      De[i1][i2] = Delta(k,r);
//      Ge[i1][i2] = G(k,t,r);
//      ks[i1][i2] = k;
//	  cout<<k[0]<<" "<<k[1]<<" "<<en1[i1][i2]<<" "<<en2[i1][i2]<<" "<< 1.0*De[i1][i2] <<endl;
    }
  }
//return 0;
  const int nT = 128;

  for (int iT = 0; iT < nT; iT++) {
    const double T = ((double)iT/(double)nT)*(Tmax-Tmin) + Tmin;

    double C = 0.0;
    for (int i1 = 0; i1 < L; i1++) {
      for (int i2 = 0; i2 < L; i2++) {
        C += kernel(ks[i1][i2],t,r,T,en1[i1][i2],en2[i1][i2]);
      }
    }
    C /= (double)(L*L*2.0);
    cout<<T<<" "<<C<<endl;

  }
  return 0;
}


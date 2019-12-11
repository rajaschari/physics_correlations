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

complex<double> comDiv(complex<double> aa,complex<double> bb){

    double a=real(aa);
    double b=imag(aa);
    double c=real(bb);
    double d=imag(bb);

    double c1=(a*c-b*d)/(c*c+d*d);
    double c2=(b*c+a*d)/(c*c+d*d);

    return c1+1i*c2;

}


vector<double> zeroVec (2,0.0);

double n1p[] = {1.0/2.0,sqrt(3.0)/2.0};
double n2p[] = {-1.0/2.0,sqrt(3.0)/2.0};

double b1p[] = {2.0*pi,2.0*pi/sqrt(3.0)};
double b2p[] = {-2.0*pi,2.0*pi/sqrt(3.0)};

double rxp[] = {1.0,1.0/sqrt(3.0)};
double ryp[] = {-1.0,1.0/sqrt(3.0)};
double rzp[] = {0.0,-2.0/sqrt(3.0)};

double rxyp[] = {-3.0/2.0,-1.0/(2.0*sqrt(3.0))};
double rxzp[] = {-1.0,-2.0/sqrt(3.0)};
double ryxp[] = {-3.0/2.0,-1.0/(2.0*sqrt(3.0))};
double ryzp[] = {1.0,-2.0/sqrt(3.0)};
double rzxp[] = {1.0/2.0,5.0/2*sqrt(3.0)};
double rzyp[] = {-1.0/2.0,5.0/2*sqrt(3.0)};

vector<double> n1 (n1p, n1p + sizeof(n1p) / sizeof(double) );
vector<double> n2 (n2p, n2p + sizeof(n2p) / sizeof(double) );

vector<double> b1 (b1p, b1p + sizeof(b1p) / sizeof(double) );
vector<double> b2 (b2p, b2p + sizeof(b2p) / sizeof(double) );

vector<double> rx (rxp, rxp + sizeof(rxp) / sizeof(double) );
vector<double> ry (ryp, ryp + sizeof(ryp) / sizeof(double) );
vector<double> rz (rzp, rzp + sizeof(rzp) / sizeof(double) );

vector<double> rxy (rxyp, rxyp + sizeof(rxyp) / sizeof(double) );
vector<double> rxz (rxzp, rxzp + sizeof(rxzp) / sizeof(double) );
vector<double> ryx (ryxp, ryxp + sizeof(ryxp) / sizeof(double) );
vector<double> ryz (ryzp, ryzp + sizeof(ryzp) / sizeof(double) );
vector<double> rzx (rzxp, rzxp + sizeof(rzxp) / sizeof(double) );
vector<double> rzy (rzyp, rzyp + sizeof(rzyp) / sizeof(double) );

vector<double> ex = n2;
vector<double> ey = negVec(n1);
vector<double> ez = sumVec(n1,negVec(n2));
//
//ex = sumVec(n2,zeroVec,ex);
//ey = sumVec(negVec(n1),zeroVec,ey);
//ez = sumVec(n1,negVec(n2),ez);

double Ep[] = {1.0/sqrt(2.0),-1.0/sqrt(2.0),0};
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

    return  2*Ryx(t)*Rzx(t) ;
}

double ty(double t){

    return  2*Rxy(t)*Rzy(t) ;
}

double tz(double t){

    return  2*Rxz(t)*Ryx(t) ;
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

    return 2*Bx*By*Bz/pow(0.27,2.0);
}


double Qx(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  7.40741*(By*(-Ex - 2*Ey) + Bz*(Ex + 2*Ez)) ;
}

double Qy(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  7.40741*(Bz*(-Ey - 2*Ez) + Bx*(2*Ex + Ez)) ;
}

double Qz(double t,double r){

    double Bx=r*Bp[0];
    double By=r*Bp[1];
    double Bz=r*Bp[2];

    double Ex=t*Ep[0];
    double Ey=t*Ep[1];
    double Ez=t*Ep[2];

    return  7.40741*(Bz*(-2*Ex - Ez) + By*(2*Ey + Ez));
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

complex<double> G(vector<double> k, double t, double r){
//    double kx=k[0];
//    double ky=k[1];

    return  4.0*((double(Qx(t,r)))*sin(dot(k,ex))+(double(Qy(t,r)))*sin(dot(k,ey))+(double(Qz(t,r)))*sin(dot(k,ez)));
}

complex<double> T(vector<double> k, double t){
//    double kx=k[0];
//    double ky=k[1];

    return  (double(tx(t)))*exp((complex<double>) 1i*dot(k,rx))+(double(ty(t)))*exp((complex<double>) 1i*dot(k,ry))+(double(tz(t)))*exp((complex<double>) 1i*dot(k,rz));
}

complex<double> P(vector<double> k, double t){
//    double kx=k[0];
//    double ky=k[1];

    return  (double(pow(Rxy(t),2)))*exp((complex<double>) 1i*dot(k,rxy))+(double(pow(Rxz(t),2)))*exp((complex<double>) 1i*dot(k,rxz))+
    (double(pow(Ryx(t),2)))*exp((complex<double>) 1i*dot(k,ryx))+(double(pow(Ryz(t),2)))*exp((complex<double>) 1i*dot(k,ryz))+
    (double(pow(Rzx(t),2)))*exp((complex<double>) 1i*dot(k,rzx))+(double(pow(Rzy(t),2)))*exp((complex<double>) 1i*dot(k,rzy));
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

complex<double> ras1(vector<double> k, double t, double r){

    return comDiv((real(f1e(k,t,r))+sqrt(pow(real(f1e(k,t,r)),2.0)+pow(imag(f1o(k,t,r)),2.0)+pow(Delta(k,r),2.0))),(f1o(k,t,r)-Delta(k,r)));
}

complex<double> ras2(vector<double> k, double t, double r){

    return comDiv((real(f1e(k,t,r))-sqrt(pow(real(f1e(k,t,r)),2.0)+pow(imag(f1o(k,t,r)),2.0)+pow(Delta(k,r),2.0))),(f1o(k,t,r)-Delta(k,r)));
}

vector<double> psi1(double c1,double c2,double t,double r, int L){

    vector<double> k(2);
    k[0] = b1[0]*(double)c1/(double)L+b2[0]*(double)c2/(double)L;
    k[1] = b1[1]*(double)c1/(double)L+b2[1]*(double)c2/(double)L;

    vector<double> ps1(4);

    ps1[0]=real(ras1(k,t,r));
    ps1[1]=imag(ras1(k,t,r));
    ps1[2]=1.0;
    ps1[3]=0.0;

    return ps1;

}

vector<double> psi2(double c1,double c2,double t,double r, int L){

    vector<double> k(2);
    k[0] = b1[0]*(double)c1/(double)L+b2[0]*(double)c2/(double)L;
    k[1] = b1[1]*(double)c1/(double)L+b2[1]*(double)c2/(double)L;

    vector<double> ps2(4);

    ps2[0]=real(ras2(k,t,r));
    ps2[1]=imag(ras2(k,t,r));
    ps2[2]=1.0;
    ps2[3]=0.0;

    return ps2;
}

//level 5 -- phi, Pp , F , kernel(berry phase)

complex<double> phi1(double c1,double c2,double c1p,double c2p,double t, double r,int L){

    complex<double> h1=psi1(c1,c2,t,r,L)[0]+1i*psi1(c1,c2,t,r,L)[1];
    complex<double> h2=psi1(c1,c2,t,r,L)[2]+1i*psi1(c1,c2,t,r,L)[3];

    complex<double> s11=conj(h1);
    complex<double> s12=conj(h2);

    complex<double> s21=psi1(c1p,c2p,t,r,L)[0]+1i*psi1(c1p,c2p,t,r,L)[1];
    complex<double> s22=psi1(c1p,c2p,t,r,L)[2]+1i*psi1(c1p,c2p,t,r,L)[3];


    return (s11*s21)+(s12*s22);
}

complex<double> phi2(double c1,double c2,double c1p,double c2p,double t, double r,int L){

    complex<double> h1=psi2(c1,c2,t,r,L)[0]+1i*psi1(c1,c2,t,r,L)[1];
    complex<double> h2=psi2(c1,c2,t,r,L)[2]+1i*psi1(c1,c2,t,r,L)[3];

    complex<double> s11=conj(h1);
    complex<double> s12=conj(h2);

    complex<double> s21=psi2(c1p,c2p,t,r,L)[0]+1i*psi1(c1p,c2p,t,r,L)[1];
    complex<double> s22=psi2(c1p,c2p,t,r,L)[2]+1i*psi1(c1p,c2p,t,r,L)[3];


    return (s11*s21)+(s12*s22);
}

complex<double> Pp1(double wf00_0,double wf10_0,double wf11_0,double wf01_0,double wf00_1,double wf10_1,double wf11_1,double wf01_1, double t,double r,int L){

    complex<double> pr0010= (1.0 + ((wf00_0*wf10_0)+(wf00_1*wf10_1)))+ 1i*(wf00_0*wf10_1 - wf00_1*wf10_0);
    complex<double> pr1011= (1.0 + ((wf10_0*wf11_0)+(wf10_1*wf11_1)))+ 1i*(wf10_0*wf11_1 - wf10_1*wf11_0);
    complex<double> pr1101= (1.0 + ((wf11_0*wf01_0)+(wf11_1*wf01_1)))+ 1i*(wf11_0*wf01_1 - wf11_1*wf11_0);
    complex<double> pr0100= (1.0 + ((wf01_0*wf00_0)+(wf01_1*wf00_1)))+ 1i*(wf01_0*wf00_1 - wf01_1*wf00_0);

//    cout << pr0010 << " " << pr1011 << " " << pr1101 << " " << pr0100 << endl;
    cout << pr0010*pr1011*pr1101*pr0100 << endl;
    cout << arg(pr0010*pr1011*pr1101*pr0100) << endl;

    return pr0010*pr1011*pr1101*pr0100;
}

complex<double> Pp2(double c1,double c2, double t,double r,int L){

    return phi2(c1,c2,c1+1,c2,t,r,L)*phi2(c1+1,c2,c1+1,c2+1,t,r,L)*phi2(c1+1,c2+1,c1,c2+1,t,r,L)*phi2(c1,c2+1,c1,c2,t,r,L);
}

double kernel1(double e11,double wf00_0,double wf10_0,double wf11_0,double wf01_0,double wf00_1,double wf10_1,double wf11_1,double wf01_1,double t,double r,int L){


    complex<double> pr0010= (1.0 + ((wf00_0*wf10_0)+(wf00_1*wf10_1)))+ 1i*(wf00_0*wf10_1 - wf00_1*wf10_0);
    complex<double> pr1011= (1.0 + ((wf10_0*wf11_0)+(wf10_1*wf11_1)))+ 1i*(wf10_0*wf11_1 - wf10_1*wf11_0);
    complex<double> pr1101= (1.0 + ((wf11_0*wf01_0)+(wf11_1*wf01_1)))+ 1i*(wf11_0*wf01_1 - wf11_1*wf11_0);
    complex<double> pr0100= (1.0 + ((wf01_0*wf00_0)+(wf01_1*wf00_1)))+ 1i*(wf01_0*wf00_1 - wf01_1*wf00_0);

    cout << pr0010 << " " << pr1011 << " " << pr1101 << " " << pr0100 << endl;
    cout << wf11_0 << " " << wf11_1 << endl;



    if(e11<0 && pr0010*pr1011*pr1101*pr0100!=(0.0,0.0)){
        return arg(pr0010*pr1011*pr1101*pr0100);
    }
    else{
        return 0;
    }
}

double kernel2(double c1,double c2,double t,double r,double e21,int L){

    if(e21<0){
        return arg(Pp2(c1,c2,t,r,L));
    }
    else{
        return 0;
    }
}

int main(int argc,char** argv) {

//  const int L = atoi(argv[1]);
//  const double r = atof(argv[2]);
//  const double tmin=atof(argv[3]);
//  const double tmax=atof(argv[4]);

    const int L=16;
    const double r= 0.1;
    const double tmin=0;
    const double tmax=1;


//real code starts here

  double en1[L][L];
  double en2[L][L];

  double psit10[L][L];
  double psit11[L][L];
//  double psit12[L][L];
//  double psit13[L][L];

  double psit20[L][L];
  double psit21[L][L];
//  double psit22[L][L];
//  double psit23[L][L];

  vector<double> ks[L][L];

//
////  for (int i1 = 0; i1 < L; i1++) {
////    for (int i2 = 0; i2 < L; i2++) {
////      vector<double> k(2);
////      k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
////      k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;
////      en1[i1][i2] = e1(k,t,r);
////      en2[i1][i2] = e2(k,t,r);
////      ks[i1][i2] = k;
////	//cout<<k[0]<<" "<<k[1]<<" "<<en1[i1][i2]<<" "<<en2[i1][i2]<<" "<<real(G(k,t,r))<<endl;
////    }
////  }
//
//  const int nt = 16;
//
//  for (int it = 0; it < nt; it++) {
//    const double tt = ((double)it/(double)nt)*(tmax-tmin) + tmin;
//
//    for (int i1 = 0; i1 < L; i1++) {
//     for (int i2 = 0; i2 < L; i2++) {
//       vector<double> k(2);
//       k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
//       k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;
//       en1[i1][i2] = e1(k,tt,r);
//       en2[i1][i2] = e2(k,tt,r);
//
//       psit10[i1][i2]=psi1(i1,i2,tt,r,L)[0];
//       psit11[i1][i2]=psi1(i1,i2,tt,r,L)[1];
////       psit12[i1][i2]=psi1(i1,i2,tt,r,L)[2];
////       psit13[i1][i2]=psi1(i1,i2,tt,r,L)[3];
//
//       psit20[i1][i2]=psi2(i1,i2,tt,r,L)[0];
//       psit21[i1][i2]=psi2(i1,i2,tt,r,L)[1];
////       psit22[i1][i2]=psi2(i1,i2,tt,r,L)[2];
////       psit23[i1][i2]=psi2(i1,i2,tt,r,L)[3];
//
//       ks[i1][i2] = k;
//	//cout<<k[0]<<" "<<k[1]<<" "<<en1[i1][i2]<<" "<<en2[i1][i2]<<" "<<real(G(k,t,r))<<endl;
//      }
//    }
//
//
//
//    double C = 0.0;
//    for (int i1 = 1; i1 < L; i1++) {
//      for (int i2 = 1; i2 < L; i2++) {
//
//        vector<double> k(2);
//        k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
//        k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;
////
////        C += kernel1(i1,i2,tt,r,e1(k,tt,r),L) + kernel2(i1,i2,tt,r,e2(k,tt,r),L);
//        if((i1%L+i2%L!=0) || ((i1+1)%L+(i2+1)%L)!=0){
//            C += kernel1(en1[i1][i2],psit10[i1][i2],psit10[i1+1][i2],psit10[i1+1][i2+1],psit10[i1][i2+1],
//                  psit11[i1][i2],psit11[i1+1][i2],psit11[i1+1][i2+1],psit11[i1][i2+1],tt,r,L) ;
//        }
////          C += kernel2(i1,i2,tt,r,e2(k,tt,r),L);
//      }
//    }
////    C /= (double)(L*L*2.0);
//
//    cout<<tt<<" "<<C<<endl;
//}

//real code ends here

  int i1=0;
//  int i2=0;
  double tt=0.1;
//  vector<double> k(2);
//  k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
//  k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;

  for (int i1 = 0; i1 < L; i1++) {
     for (int i2 = 0; i2 < L; i2++) {
       vector<double> k(2);
       k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
       k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;
       en1[i1][i2] = e1(k,tt,r);
       en2[i1][i2] = e2(k,tt,r);

       psit10[i1][i2]=psi1(i1,i2,tt,r,L)[0];
       psit11[i1][i2]=psi1(i1,i2,tt,r,L)[1];
//       psit12[i1][i2]=psi1(i1,i2,tt,r,L)[2];
//       psit13[i1][i2]=psi1(i1,i2,tt,r,L)[3];

       psit20[i1][i2]=psi2(i1,i2,tt,r,L)[0];
       psit21[i1][i2]=psi2(i1,i2,tt,r,L)[1];
//       psit22[i1][i2]=psi2(i1,i2,tt,r,L)[2];
//       psit23[i1][i2]=psi2(i1,i2,tt,r,L)[3];

       ks[i1][i2] = k;
	//cout<<k[0]<<" "<<k[1]<<" "<<en1[i1][i2]<<" "<<en2[i1][i2]<<" "<<real(G(k,t,r))<<endl;
      }
    }


//
//  cout << e1(k,0.1,0.1) <<  " " << e2(k,0.1,0.1) <<  " " << kernel1(1,1,0.1,0.1,e1(k,0.1,0.1),128) << endl;
  double C = 0.0;
    for (int i1 = 0; i1 < L; i1++) {
      for (int i2 = 0; i2 < L; i2++) {

        vector<double> k(2);
        k[0] = b1[0]*(double)i1/(double)L+b2[0]*(double)i2/(double)L;
        k[1] = b1[1]*(double)i1/(double)L+b2[1]*(double)i2/(double)L;
//
//        C += kernel1(i1,i2,tt,r,e1(k,tt,r),L) + kernel2(i1,i2,tt,r,e2(k,tt,r),L);


        if((i1%L+i2%L)!=0 && ((i1+1)%L+(i2)%L)!=0  && ((i1+1)%L+(i2+1)%L)!=0 && ((i1)%L+(i2+1)%L)!=0 ){

            C += kernel1(en1[i1][i2],psit10[i1][i2],psit10[i1+1][i2],psit10[i1+1][i2+1],psit10[i1][i2+1],
                        psit11[i1][i2],psit11[i1+1][i2],psit11[i1+1][i2+1],psit11[i1][i2+1],tt,r,L) ;


             cout << i1 <<" "<< i2 <<  " " << C << endl;
        }
      }
    }

//
//
//  double tt=0.1;
//  int i11=2;
//  int i12=3;
//  vector<double> k1(2);
//  k1[0] = b1[0]*(double)i11/(double)L+b2[0]*(double)i12/(double)L;
//  k1[1] = b1[1]*(double)i11/(double)L+b2[1]*(double)i12/(double)L;
//
//  cout << i11 << " "<< i12 << " " <<kernel1(en1[i11][i12],psit10[i11][i12],psit10[i11+1][i12],psit10[i11+1][i12+1],psit10[i11][i12+1],
//        psit11[i11][i12],psit11[i11+1][i12],psit11[i11+1][i12+1],psit11[i11][i12+1],tt,r,L) <<endl;
//  cout << e1(k1,0.1,0.1) <<  " " << kernel1(en1[i11][i12],psit10[i11][i12],psit10[i11+1][i12],psit10[i11+1][i12+1],psit10[i11][i12+1],
//        psit11[i11][i12],psit11[i11+1][i12],psit11[i11+1][i12+1],psit11[i11][i12+1],tt,r,L) << endl;
//
//
//  cout << Pp1(psit10[i11][i12],psit10[i11+1][i12],psit10[i11+1][i12+1],psit10[i11][i12+1],
//        psit11[i11][i12],psit11[i11+1][i12],psit11[i11+1][i12+1],psit11[i11][i12+1],tt,r,L)  << endl;
//
//
//  cout << psit10[i11][i12] << " " << psit10[i11+1][i12] << " " << psit10[i11+1][i12+1] << " " << psit10[i11][i12+1]
//        <<" "<<psit11[i11][i12]<<" "<<psit11[i11+1][i12]<<" "<<psit11[i11+1][i12+1]<<" "<<psit11[i11][i12+1] << endl;
//
//
//  cout << psi1(0,0,0.1,0.1,4)[0] <<  " " << psi1(0,0,0.1,0.1,4)[1]<<  " " << psi1(0,0,0.1,0.1,4)[2]<< " " << psi1(0,0,0.1,0.1,4)[3]<< endl;
//  cout << ras1(k1,0.1,0.1) <<  " " << 1 <<  " " << 1 << " " << 1 << endl;
//  cout << f1o(k1,0.1,0.1) <<" "<< Delta(k1,0.1) << " "<< K(0.1) <<endl;

  return 0;
}

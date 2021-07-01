#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"


using namespace std;


double error(double  av, double av2, int k){
  if(k==0)
    return 0;
  else
    return sqrt((av2-av*av)/double (k));
}

void RW3Ddisc(double *a, double random){
  int q=random+1;
  switch(q){
    case 1:
    a[0]++;
    break;
    case 2:
    a[0]--;
    break;
    case 3:
    a[1]++;
    break;
    case 4:
    a[1]--;
    break;
    case 5:
    a[2]++;
    break;
    case 6:
    a[2]--;
    break;
  }
}

double thetaHM(double random1, double random2){
  if (random2<0){
    return 2*M_PI-acos(random1/sqrt(random1*random1+random2*random2));
  }
  else
    return acos(random1/sqrt(random1*random1+random2*random2));
}

double phiHM(double random1, double random2){
  if(random1>0)
    return asin(random2/sqrt(random1*random1+random2*random2));
  else
    return M_PI-asin(random2/sqrt(random1*random1+random2*random2));
}

void makevet(double theta, double phi, double *c){
  c[0]+=sin(phi)*cos(theta);
  c[1]+=sin(phi)*sin(theta);
  c[2]+=cos(phi);
}

void makevet2(double theta, double phi, double *c){
  c[0]+=sin(phi)*cos(theta)*sin(phi)*cos(theta);
  c[1]+=sin(phi)*sin(theta)*sin(phi)*sin(theta);
  c[2]+=cos(phi)*cos(phi);
}

int main(){

  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;


int M=10000;
int N=100;
int A=M/N,i,j,k;
//double ave[N], ave2[N], sum[N], sum2[N], err[N];
//double aveis[N], ave2is[N], sumis[N], sum2is[N], erris[N];
int v[3]= {0};
double c[3]={0};
double mm[N][3], mm2[N][3];
double cl2[N][3];

/*
for(i=0; i<3; i++){
  cout << v[i]<<endl;
}

for(i=0; i<N; i++){
  cout << "Blocco numero: "<< i+1<<endl;
  for( j=0; j<A; j++){
    RW3Ddisc(v,rnd.Rannyu(0,6));
  }

  for(k=0; k<3; k++){
    cout << v[k]<<endl;
    v[k]=0;
  }
  cout << endl;
}

*/

/*
for(i=0; i<3; i++){
  cout << c[i]<<endl;
}

for(i=0; i<N; i++){
  cout << "Blocco numero: "<< i+1<<endl;
  for( j=0; j<A; j++){
    double r1=rnd.Rannyu(-1,1);
    double r2=rnd.Rannyu(-1,1);
    double r3=rnd.Rannyu(-1,1);
    double r4=rnd.Rannyu();

    while (r1*r1+r2*r2>1){
      r1=rnd.Rannyu(-1,1);
      r2=rnd.Rannyu(-1,1);
    }

    while (r3*r3+r4*r4>1){
      r3=rnd.Rannyu(-1,1);
      r4=rnd.Rannyu();
    }

    double theta=thetaHM(r1,r2);
    double phi = phiHM(r3,r4);
    //cout << "theta: "<< theta<<"  phi: "<<phi<<endl;
    makevet(theta, phi, c);
  }

  for(k=0; k<3; k++){
    cout << c[k]<<endl;
    c[k]=0;
  }
  cout << endl;
}

*/

cout << "t: 0"<<endl;


for(i=0; i<N; i++){
  cout << "t: "<< i+1<<endl;
  mm[i][0]=0;
  mm[i][1]=0;
  mm[i][2]=0;
  for(j=0; j<M; j++){
    double r1=rnd.Rannyu(-1,1);
    double r2=rnd.Rannyu(-1,1);
    double r3=rnd.Rannyu(-1,1);
    double r4=rnd.Rannyu();

    while (r1*r1+r2*r2>1){
      r1=rnd.Rannyu(-1,1);
      r2=rnd.Rannyu(-1,1);
    }

    while (r3*r3+r4*r4>1){
      r3=rnd.Rannyu(-1,1);
      r4=rnd.Rannyu();
    }

    double theta=thetaHM(r1,r2);
    double phi = phiHM(r3,r4);
    //cout << "theta: "<< theta<<"  phi: "<<phi<<endl;
    makevet(theta, phi, mm[i]);
    //cout << mm[i][0]<<"   "<< mm[i][0]<< "   " << mm[i][0]<<endl;

  }
  if(i!=0){
    mm[i][0]+=mm[i-1][0];
    mm[i][1]+=mm[i-1][1];
    mm[i][2]+=mm[i-1][2];
  }

  cout << mm[i][0]<<"   "<< mm[i][1]<< "   " << mm[i][2]<<endl;

  cout << endl;
}

for(i=0; i<N; i++){
  cout << "t: "<<i+1<<endl;
  cout << mm[i][0]/double (M)<< "  "<< mm[i][1]/double (M)<<"   "<<mm[i][2]/double (M)<<endl;
}





cout << "t: 0"<<endl;


for(i=0; i<N; i++){
  cout << "t: "<< i+1<<endl;
  mm2[i][0]=0;
  mm2[i][1]=0;
  mm2[i][2]=0;
  for(j=0; j<M; j++){
    double r1=rnd.Rannyu(-1,1);
    double r2=rnd.Rannyu(-1,1);
    double r3=rnd.Rannyu(-1,1);
    double r4=rnd.Rannyu();

    while (r1*r1+r2*r2>1){
      r1=rnd.Rannyu(-1,1);
      r2=rnd.Rannyu(-1,1);
    }

    while (r3*r3+r4*r4>1){
      r3=rnd.Rannyu(-1,1);
      r4=rnd.Rannyu();
    }

    double theta=thetaHM(r1,r2);
    double phi = phiHM(r3,r4);
    //cout << "theta: "<< theta<<"  phi: "<<phi<<endl;
    makevet2(theta, phi, mm2[i]);
    //cout << mm[i][0]<<"   "<< mm[i][0]<< "   " << mm[i][0]<<endl;

  }
  if(i!=0){
    mm2[i][0]+=mm2[i-1][0];
    mm2[i][1]+=mm2[i-1][1];
    mm2[i][2]+=mm2[i-1][2];
  }

  cout << mm2[i][0]<<"   "<< mm2[i][1]<< "   " << mm2[i][2]<<endl;

  cout << endl;
}

for(i=0; i<N; i++){
  cout << "t: "<<i+1<<endl;
  cout << sqrt(mm2[i][0]/double (M))<< "  "<< sqrt(mm2[i][1]/double (M))<<"   "<<sqrt(mm2[i][2]/double (M))<<endl;
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

cout << "t: 0"<<endl;


for(i=0; i<N; i++){
  cout << "t: "<< i+1<<endl;
  cl2[i][0]=0;
  cl2[i][1]=0;
  cl2[i][2]=0;
  for(j=0; j<M; j++){
    //cout << "theta: "<< theta<<"  phi: "<<phi<<endl;
    RW3Ddisc(cl2[i],rnd.Rannyu(0,6));
    //cout << mm[i][0]<<"   "<< mm[i][0]<< "   " << mm[i][0]<<endl;

  }

  if(i!=0){
    cl2[i][0]+=cl2[i-1][0];
    cl2[i][1]+=cl2[i-1][1];
    cl2[i][2]+=cl2[i-1][2];
  }

  //cout << cl2[i][0]<<"   "<< cl2[i][1]<< "   " << cl2[i][2]<<endl;

  cout << endl;
}

for(i=0; i<N; i++){
  cout << "t: "<<i+1<<endl;
  //cout << sqrt(cl2[i][0]/double (M))<< "  "<< sqrt(cl2[i][1]/double (M))<<"   "<<sqrt(cl2[i][2]/double (M))<<endl;
  cout << cl2[i][0]<<"   "<< cl2[i][1]<< "   " << cl2[i][2]<<endl;
  double q0=cl2[i][0]/double(M);
  double q1=cl2[i][1]/double(M);
  double q2=cl2[i][2]/double(M);
  cout << q0<<"   "<< q1<< "   " << q2<<endl;
  //cout << sqrt(q0/double (M))<< "   "<< sqrt(q1/double (M))<<"   "<<sqrt(q2/double (M))<<endl;
}


  rnd.SaveSeed();

  return 0;
}

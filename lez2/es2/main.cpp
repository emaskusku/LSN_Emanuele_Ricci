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

double modulo (double *v){
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

int main(){

  Random rnd;
  int seed[5];
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


  int M=100000;
  int N=100;
  int L=M/N, i,j,k;
  double mc[N][3], mod[N], modcl[N];
  ofstream fileout;
  //double mc2[N][3];

  for(i=0; i<N; i++){
    mod[i]=0;
  }

  for(i=0; i<M; i++){
    for(j=0; j<N; j++){
      mc[j][0]=0;
      mc[j][1]=0;
      mc[j][2]=0;

      if(j!=0){
        mc[j][0]+=mc[j-1][0];
        mc[j][1]+=mc[j-1][1];
        mc[j][2]+=mc[j-1][2];
      }

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
      makevet(theta,phi,mc[j]);
    }

    for(k=0; k<N; k++){
      mod[k]+=modulo(mc[k]);
    }
  }

  fileout.open("RWCont.txt");

  for(i=0; i<N; i++){
    mod[i]=sqrt(mod[i]/double (M));
    cout <<mod[i]<<endl;
    fileout<<mod[i]<<endl;
  }

  fileout.close();

  /******************************************************/

  for(i=0; i<N; i++){
    modcl[i]=0;
  }

  for(i=0; i<M; i++){
    for(j=0; j<N; j++){
      mc[j][0]=0;
      mc[j][1]=0;
      mc[j][2]=0;

      if(j!=0){
        mc[j][0]+=mc[j-1][0];
        mc[j][1]+=mc[j-1][1];
        mc[j][2]+=mc[j-1][2];
      }

      RW3Ddisc(mc[j], rnd.Rannyu(0,6));
    }

    for(k=0; k<N; k++){
      modcl[k]+=modulo(mc[k]);
    }
  }

  fileout.open("RWCL.txt");

  for(i=0; i<N; i++){
    modcl[i]=sqrt(modcl[i]/double (M));
    cout <<modcl[i]<<endl;
    fileout<<modcl[i]<<endl;
  }

  fileout.close();












  rnd.SaveSeed();

  return 0;
}

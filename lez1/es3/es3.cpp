#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;

double module (double a, double b, double a1, double b1){
  return sqrt((a-a1)*(a-a1)+(b-b1)*(b-b1));
}

double error(double *av, double *av2, int k){
  if (k==0)
    return 0;
  else{
    return sqrt((av2[k]-av[k]*av[k])/double(k));
  }
}

double thetaHM(double random1, double random2){
  if (random2<0){
    return 2*M_PI-acos(random1/sqrt(random1*random1+random2*random2));
  }
  else
    return acos(random1/sqrt(random1*random1+random2*random2));
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


  double x,y,x1,y1;
  int M=100000;
  int N=100;
  int L=M/N;
  double min=0.;
  double max=10000.;
  double d=2.;
  double l=1.5;
  int v[N],v2[N];
  double ave[N],ave2[N], err[N];
  ofstream fileout;

  for(int i=0; i<N; i++){
    v[i]=0;
    for(int j=0; j<L; j++){
      x=rnd.Rannyu(min,max);
      y=rnd.Rannyu(min,max);

      x1=rnd.Rannyu(-1,1);
      y1=rnd.Rannyu(-1,1);

      while (x1*x1+y1*y1>1){
        x1=rnd.Rannyu(-1,1);
        y1=rnd.Rannyu(-1,1);
      }

      double theta=thetaHM(x1,y1);
      double xf=x+l*cos(theta);
      double yf=y+l*sin(theta);

      for(int k=0; k<max/d-1; k++){
        if((yf-d*(k+1))*(y-d*(k+1))<0)
          v[i]++;

        if(y1==0. || y==0.)
          v[i]++;
      }
    }
    v2[i]=v[i]*v[i];
  }

  for(int i=0; i<N; i++){
    ave[i]=0;
    ave2[i]=0;
    for(int j=0; j<i+1; j++){
      double accu = double (v[j]*v[j]);
      double appo =double (i+1);
      //cout <<accu<<endl;
      ave[i]+=double (v[j])/double(i+1);
      ave2[i]+=accu/appo;
    }

    err[i]=error(ave,ave2,i);
  }

  fileout.open("misurePi.txt");

  for(int i=0; i<N; i++) {
      int x=L*i+L;
      double pi=2*l*L/(ave[i]*d);
      double erpi=err[i]*2*l*L/(ave[i]*ave[i]*d);
      cout << pi<< "  "<< erpi <<endl;
      fileout << x << "  "<<pi << "  "<< erpi<<endl;
  }


  fileout.close();
  rnd.SaveSeed();

  return 0;
}

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

double dprob (double random){
  double x=asin(random)*2/M_PI;
  return x;
}

double dprob2(double random){
  return 1.+sqrt(1-random);
}

double g (double x){
  double coseno=M_PI*cos(M_PI*x/2.)/2.;
  double d=(1.-x)*2.;
  return coseno/d;
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
int A=M/N,i,j;
double ave[N], ave2[N], sum[N], sum2[N], err[N];
double aveis[N], ave2is[N], sumis[N], sum2is[N], erris[N];

for(i=0; i<N; i++){
  ave[i]=0;
  ave2[i]=0;
  for(j=0; j<A; j++){
    ave[i]+=M_PI/2*cos(M_PI*rnd.Rannyu()/2);
  }
  ave[i]=ave[i]/double (A);
  ave2[i]=ave[i]*ave[i];
  //cout << ave[i]<<"    "<< ave2[i]<<endl;

}

for(i=0; i<N; i++){
  sum[i]=0;
  sum2[i]=0;
  for(j=0; j<i+1; j++){
    sum[i]+=ave[j];
    sum2[i]+=ave2[j];
  }
  sum[i]=sum[i]/double (i+1);
  sum2[i]=sum2[i]/double (i+1);
  cout <<sum[i]<< "  " << sum2[i]<<endl;
}



cout << endl<<"Importance sampling"<<endl<<endl;

for(i=0; i<N; i++){
  aveis[i]=0;
  ave2is[i]=0;
  for(j=0; j<A; j++){
    aveis[i]+=g(dprob2(rnd.Rannyu()));
  }
  aveis[i]=aveis[i]/double (A);
  ave2is[i]=aveis[i]*aveis[i];
  cout << aveis[i]<<"    "<< ave2is[i]<<endl;

}

for(i=0; i<N; i++){
  sumis[i]=0;
  sum2is[i]=0;
  for(j=0; j<i+1; j++){
    sumis[i]+=aveis[j];
    sum2is[i]+=ave2is[j];
    //cout <<sumis[i]<< "  " << sum2is[i]<<endl;
  }
  sumis[i]=sumis[i]/double (i+1);
  sum2is[i]=sum2is[i]/double (i+1);
  cout <<sumis[i]<< "  " << sum2is[i]<<endl;
}

cout << "varianza"<<endl;

for(i=0; i<N; i++){
  err[i]=error(sum[i], sum2[i], i);
  erris[i]=error(sumis[i], sum2is[i], i);
  cout << err[i]<< "  "<<erris[i]<<"  "<< A+A*i<<endl;
  //cout << sum2[i]-sum[i]*sum[i]<< "  "<<sum2is[i]-sumis[i]*sumis[i]<<endl ;
}



ofstream fileout;

fileout.open("IntUnif.txt");
for(i=0; i<N; i++){
  fileout<<A+A*i<<"  "<<sum[i]<< "  "<<err[i]<<endl;
}
fileout.close();



fileout.open("IntIS.txt");
for(i=0; i<N; i++){
  fileout<<A+A*i<<"  "<<sumis[i]<< "  "<<erris[i]<<endl;
}
fileout.close();





  rnd.SaveSeed();

  return 0;
}

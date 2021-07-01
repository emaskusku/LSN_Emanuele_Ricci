/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#define I 100.

using namespace std;


double N(double x){
    return 0.5 * (1. + erf(x / sqrt(2.)));
  }

double black_scholes_C( double S0,double K,double T,double r, double sigma){
  double d1 = 1./(sigma * sqrt(T)) * (log(S0 / K) + (r + (sigma*sigma) / 2.) * T);
  double d2 = d1 - sigma * sqrt(T);
  double C = S0 * N(d1) - K * exp(-r * T) * N(d2);
  return C;
}

double black_scholes_P(double S0,double K,double T,double r, double sigma){
  double d1 = 1./(sigma * sqrt(T)) * (log(S0 / K) + (r + (sigma*sigma) / 2.) * T);
  double d2 = d1 - sigma * sqrt(T);
  double P = S0 *(N(d1) - 1.) - K * exp(-r * T) * (N(d2)-1.);
  return P;
}


double AP(double S0,double T,double r, double sigma, double gauss){
  return S0*exp((r-0.5*sigma*sigma)*T+sigma*gauss);
}

double AP_discr(double r, double sigma, double gauss, double passo){
  return exp((r-0.5*sigma*sigma)*passo+sigma*gauss*sqrt(passo));
}

double error(double *av, double *av2, int k){
  if (k==0)
    return 0;
  else{
    return sqrt((av2[k]-av[k]*av[k])/double(k));
  }
}


int main (int argc, char *argv[]){

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


  double S0 = 100.;
  double K = 100.;
  double T = 1.;
  double r = 0.1;
  double sigma = 0.25;
  double p;
  int M=100000;
  int N=100;
  int L=int(M/N);
  double ave[N];
  double av2[N];
  double sum_prog[N];
  double su2_prog[N];
  double err_prog[N];
  double sum;
  int i,j;
  double passo=T/I;

  ofstream fileout;


  for(i=0; i<N; i++){
    ave[i]=0;
    sum=0;
    for(j=0; j<L; j++){
      p=(AP(S0, T, r, sigma, rnd.Gauss(0,1))-K)*exp(-r*T);
      if(p<0)
        p=0;
        sum+=p;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
  }

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("1formulaC.txt");

  for(int k=0; k<N; k++){
    fileout<<k*L<<"  "<< sum_prog[k]<<"  "<<err_prog[k]<<endl;
  }

  fileout.close();



  for(i=0; i<N; i++){
    ave[i]=0;
    sum=0;
    for(j=0; j<L; j++){
      p=(K-AP(S0, T, r, sigma, rnd.Gauss(0,1)))*exp(-r*T);
      if(p<0)
        p=0;
        sum+=p;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
  }

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("1formulaP.txt");

  for(int k=0; k<N; k++){
    fileout<<k*L <<"  "<< sum_prog[k]<<"  "<<err_prog[k]<<endl;
  }

  fileout.close();

  /*****************************************************/
  /*****************************************************/
  /*****************************************************/

  for(i=0; i<N; i++){
    ave[i]=0;
    sum=0;
    for(j=0; j<L; j++){
      p=S0;
      for(int u=0; u<I; u++){
        p=p*AP_discr(r, sigma, rnd.Gauss(0,1), passo);
      }
        p=(p-K)*exp(-r*T);
      if(p<0)
        p=0;

      //cout << p<<endl;
      sum+=p;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
  }

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  for(int k=0; k<N; k++){
    //cout << sum_prog[k]<<"  "<<err_prog[k]<<endl;
  }


  fileout.open("1formulaCDiscr.txt");

  for(int k=0; k<N; k++){
    fileout<<k*L<<"  " << sum_prog[k]<<"  "<<err_prog[k]<<endl;
  }

  fileout.close();




  for(i=0; i<N; i++){
    ave[i]=0;
    sum=0;
    for(j=0; j<L; j++){
      p=S0;
      for(int u=0; u<I; u++){
        p=p*AP_discr(r, sigma, rnd.Gauss(0,1), passo);
      }
        p=(K-p)*exp(-r*T);
      if(p<0)
        p=0;

      //cout << p<<endl;
      sum+=p;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
  }

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  for(int k=0; k<N; k++){
    //cout << sum_prog[k]<<"  "<<err_prog[k]<<endl;
  }


  fileout.open("1formulaPDiscr.txt");

  for(int k=0; k<N; k++){
    fileout<<k*L<<"  "  << sum_prog[k]<<"  "<<err_prog[k]<<endl;
  }

  fileout.close();

   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <cmath>
#include "random.h"
#include <fstream>
#include <string>

#define NN 100


using namespace std;

int main(){

  Random rnd;
  int seed[12];
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
           input >> seed[0] >> seed[5] >> seed[7] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;


  int M=100,i,j,k,pos;
  int n=10000;
  float nAsp=float (n/M);
  int v[M];
  float chiq;

  ofstream fileout;
  fileout.open("chiq.txt");

  //cout << nAsp<<endl;

  for(i=0; i<M; i++){
    v[i]=0;
  }

  for(j=0; j<NN; j++){
    for(i=0; i<n; i++){
      float r=rnd.Rannyu();
      //cout << rnd.Rannyu()<<endl;
      pos=r*100.;
      //cout << pos << endl;
      v[pos]++;
      //cout << v[pos]<<endl;
    }
    chiq=0;

    for(k=0; k<M; k++){
      float accu= float (v[k])-nAsp;
      //cout << accu<<endl;
      chiq+=accu*accu/nAsp;
    }

    fileout << chiq<<endl;

    for(i=0; i<M; i++){
      //cout << v[i]<<endl;
      v[i]=0;
    }
  }











  return 0;
}

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

#define NN 10000

using namespace std;

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


  rnd.SaveSeed();

  double par=1;
  double mean=0;
  double par1=1;

  ofstream fileoutU1;
  ofstream fileoutE1;
  ofstream fileoutCL1;
  ofstream fileoutU2;
  ofstream fileoutE2;
  ofstream fileoutCL2;
  ofstream fileoutU10;
  ofstream fileoutE10;
  ofstream fileoutCL10;
  ofstream fileoutU100;
  ofstream fileoutE100;
  ofstream fileoutCL100;

  fileoutU1.open("U1.txt");
  fileoutE1.open("E1.txt");
  fileoutCL1.open("CL1.txt");
  fileoutU2.open("U2.txt");
  fileoutE2.open("E2.txt");
  fileoutCL2.open("CL2.txt");
  fileoutU10.open("U10.txt");
  fileoutE10.open("E10.txt");
  fileoutCL10.open("CL10.txt");
  fileoutU100.open("U100.txt");
  fileoutE100.open("E100.txt");
  fileoutCL100.open("CL100.txt");

  double sunif=0, sexp=0, scl=0;

  for(int i=0; i<NN; i++){
    for(int j=0; j<2; j++){
      sunif+=rnd.Rannyu()/2.;
      sexp+=rnd.Exp(par)/2.;
      scl+=rnd.CL(mean, par1)/2.;
    }
    fileoutU2<<sunif<<endl;
    fileoutE2<<sexp<<endl;
    fileoutCL2<<scl<<endl;
    sunif=0;
    sexp=0;
    scl=0;

    for(int j=0; j<10; j++){
      sunif+=rnd.Rannyu()/10.;
      sexp+=rnd.Exp(par)/10.;
      scl+=rnd.CL(mean, par1)/10.;
    }
    fileoutU10<<sunif<<endl;
    fileoutE10<<sexp<<endl;
    fileoutCL10<<scl<<endl;
    sunif=0;
    sexp=0;
    scl=0;

    for(int j=0; j<100; j++){
      sunif+=rnd.Rannyu()/100.;
      sexp+=rnd.Exp(par)/100.;
      scl+=rnd.CL(mean, par1)/100.;
    }
    fileoutU100<<sunif<<endl;
    fileoutE100<<sexp<<endl;
    fileoutCL100<<scl<<endl;
    sunif=0;
    sexp=0;
    scl=0;

    fileoutU1<<rnd.Rannyu()<<endl;
    fileoutE1<<rnd.Exp(par)<<endl;
    fileoutCL1<<rnd.CL(mean, par1)<<endl;
  }

  fileoutU1.close();
  fileoutE1.close();
  fileoutCL1.close();
  fileoutU2.close();
  fileoutE2.close();
  fileoutCL2.close();
  fileoutU10.close();
  fileoutE10.close();
  fileoutCL10.close();
  fileoutU100.close();
  fileoutE100.close();
  fileoutCL100.close();

  return 0;
}

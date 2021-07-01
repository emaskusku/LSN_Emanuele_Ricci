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

using namespace std;


double psi1 (double x, double y, double z){
  return exp(-2*sqrt(x*x+y*y+z*z))/M_PI;
}

double psi2 (double x, double y, double z){
  return sqrt(2)*sqrt(x*x+y*y+z*z)*exp(-sqrt(x*x+y*y+z*z)/2)*x/(8*sqrt(M_PI)*sqrt(x*x+y*y))*sqrt(2)*sqrt(x*x+y*y+z*z)*exp(-sqrt(x*x+y*y+z*z)/2)*x/(8*sqrt(M_PI)*sqrt(x*x+y*y));
}

double acc(double p){
  if(p>1)
    return 1;
  else
    return p;
}

void trans_r(double *v, double d, Random rnd){
  v[0]=v[0]+rnd.Rannyu(-1,1)*d;
  v[1]=v[1]+rnd.Rannyu(-1,1)*d;
  v[2]=v[2]+rnd.Rannyu(-1,1)*d;
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


int M=60000;
int N=100;
int L=6000;
double D=1.2;
double x[3]={1,0,0};  //con x[0]=374 cambia totalmente; -> probabilmente per la precisione double
double x_try[3]={0,0,0};
double ave_r=0,r;
int accepted=0,tries=0;
double p_old,p_new;
double A;
double sum;
double ave[N], av2[N];
double sum_prog[N], su2_prog[N];
double err_prog[N];
int i,j;
ofstream fileout;


for(int i=0; i<N; i++){
  sum=0;
  for(int j=0; j<L; j++){
  x_try[0]=x[0]+D*rnd.Rannyu(-1,1);
  x_try[1]=x[1]+D*rnd.Rannyu(-1,1);
  x_try[2]=x[2]+D*rnd.Rannyu(-1,1);

  p_old=psi1(x[0],x[1],x[2]);
  p_new=psi1(x_try[0],x_try[1],x_try[2]);

  A=p_new/p_old;
  A=acc(A);

  if(rnd.Rannyu()<A){
    x[0] = x_try[0];
    x[1] = x_try[1];
    x[2] = x_try[2];
    accepted++;
  }

  tries++;
  r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  ave_r+=r;
  sum+=r;
  //cout << x[0]<< " "<<x[1]<< " "<<x[2]<<endl;
}
ave[i]=sum/double(L);
av2[i]=ave[i]*ave[i];
}

//cout << ave_r/M<<endl;
//cout << accepted<< "  "<<tries<<endl;

for(i=0; i<N; i++){
  sum_prog[i]=0;
  su2_prog[i]=0;
  for(j=0; j<i+1; j++){
    sum_prog[i]+=ave[j];
    su2_prog[i]+=av2[j];
  }
  sum_prog[i]/=double(i+1);
  //cout << sum_prog[i]<<endl;
  su2_prog[i]/=double(i+1);
  err_prog[i]= error(sum_prog, su2_prog, i);
}


fileout.open("psi1.txt");
for(i=0; i<N; i++){
fileout << i*L+L<<" "<< sum_prog[i]<<" "<<err_prog[i]<<endl;
}
fileout.close();



double x2[3]={1,0,0};
double x2_try[3]={0,0,0};
ave_r=0;
accepted=0;
tries=0;
D=3.;
M=35000;
L=35000;


for(int i=0; i<N; i++){
  sum=0;
  for(int j=0; j<L; j++){
  x2_try[0]=x2[0]+D*rnd.Rannyu(-1,1);
  x2_try[1]=x2[1]+D*rnd.Rannyu(-1,1);
  x2_try[2]=x2[2]+D*rnd.Rannyu(-1,1);

  p_old=psi2(x2[0],x2[1],x2[2]);
  p_new=psi2(x2_try[0],x2_try[1],x2_try[2]);

  A=p_new/p_old;
  A=acc(A);

  if(rnd.Rannyu()<A){
    x2[0] = x2_try[0];
    x2[1] = x2_try[1];
    x2[2] = x2_try[2];
    accepted++;
  }

  tries++;
  r=sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
  ave_r+=r;
  sum+=r;
  //cout << x[0]<< " "<<x[1]<< " "<<x[2]<<endl;
}
ave[i]=sum/double(L);
av2[i]=ave[i]*ave[i];
}

//cout << ave_r/M<<endl;
//cout << accepted<< "  "<<tries<<endl;

for(i=0; i<N; i++){
  sum_prog[i]=0;
  su2_prog[i]=0;
  for(j=0; j<i+1; j++){
    sum_prog[i]+=ave[j];
    su2_prog[i]+=av2[j];
  }
  sum_prog[i]/=double(i+1);
  //cout << sum_prog[i]<<endl;
  su2_prog[i]/=double(i+1);
  err_prog[i]= error(sum_prog, su2_prog, i);
}


fileout.open("psi2.txt");
for(i=0; i<N; i++){
  fileout << i*L+L<<" "<< sum_prog[i]<<" "<<err_prog[i]<<endl;
}
fileout.close();



/************Gauss gauss gauss********************/



double x3[3]={1,0,0};
double x3_try[3]={0,0,0};
ave_r=0;
accepted=0;
tries=0;
D=0.7;
M=35000;
L=35000;


for(int i=0; i<N; i++){
  sum=0;
  for(int j=0; j<L; j++){
  x3_try[0]=x3[0]+D*rnd.Gauss(0,1);
  x3_try[1]=x3[1]+D*rnd.Gauss(0,1);
  x3_try[2]=x3[2]+D*rnd.Gauss(0,1);

  p_old=psi1(x3[0],x3[1],x3[2]);
  p_new=psi1(x3_try[0],x3_try[1],x3_try[2]);

  A=p_new/p_old;
  A=acc(A);

  if(rnd.Rannyu()<A){
    x3[0] = x3_try[0];
    x3[1] = x3_try[1];
    x3[2] = x3_try[2];
    accepted++;
  }

  tries++;
  r=sqrt(x3[0]*x3[0] + x3[1]*x3[1] + x3[2]*x3[2]);
  ave_r+=r;
  sum+=r;
  //cout << x[0]<< " "<<x[1]<< " "<<x[2]<<endl;
}
ave[i]=sum/double(L);
av2[i]=ave[i]*ave[i];
}

cout << ave_r/M/100<<endl;
cout << accepted<< "  "<<tries<<endl;

for(i=0; i<N; i++){
  sum_prog[i]=0;
  su2_prog[i]=0;
  for(j=0; j<i+1; j++){
    sum_prog[i]+=ave[j];
    su2_prog[i]+=av2[j];
  }
  sum_prog[i]/=double(i+1);
  //cout << sum_prog[i]<<endl;
  su2_prog[i]/=double(i+1);
  err_prog[i]= error(sum_prog, su2_prog, i);
}

fileout.open("psi1gauss.txt");
for(i=0; i<N; i++){
fileout << i*L+L<<" "<< sum_prog[i]<<" "<<err_prog[i]<<endl;
}
fileout.close();





double x4[3]={1,0,0};
double x4_try[3]={0,0,0};
ave_r=0;
accepted=0;
tries=0;
D=1.2;
M=35000;
L=35000;


for(int i=0; i<N; i++){
  sum=0;
  for(int j=0; j<L; j++){
  x4_try[0]=x4[0]+D*rnd.Gauss(0,1);
  x4_try[1]=x4[1]+D*rnd.Gauss(0,1);
  x4_try[2]=x4[2]+D*rnd.Gauss(0,1);

  p_old=psi2(x4[0],x4[1],x4[2]);
  p_new=psi2(x4_try[0],x4_try[1],x4_try[2]);

  A=p_new/p_old;
  A=acc(A);

  if(rnd.Rannyu()<A){
    x4[0] = x4_try[0];
    x4[1] = x4_try[1];
    x4[2] = x4_try[2];
    accepted++;
  }

  tries++;
  r=sqrt(x4[0]*x4[0] + x4[1]*x4[1] + x4[2]*x4[2]);
  ave_r+=r;
  sum+=r;
  //cout << x[0]<< " "<<x[1]<< " "<<x[2]<<endl;
}
ave[i]=sum/double(L);
av2[i]=ave[i]*ave[i];
}

cout << ave_r/M/100<<endl;
cout << accepted<< "  "<<tries<<endl;

for(i=0; i<N; i++){
  sum_prog[i]=0;
  su2_prog[i]=0;
  for(j=0; j<i+1; j++){
    sum_prog[i]+=ave[j];
    su2_prog[i]+=av2[j];
  }
  sum_prog[i]/=double(i+1);
  //cout << sum_prog[i]<<endl;
  su2_prog[i]/=double(i+1);
  err_prog[i]= error(sum_prog, su2_prog, i);
}

fileout.open("psi2gauss.txt");
for(i=0; i<N; i++){
fileout << i*L+L<<" "<< sum_prog[i]<<" "<<err_prog[i]<<endl;
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

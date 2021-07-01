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
#include "popolazione.h"

#define N 32

using namespace std;

double prob(double T, SimulAnn p){
  return exp(-p.GetLen()/T);
}

bool Metropolis(double p1, double p2){
  double metro_p=p1/p2;

  if(metro_p>1){
    metro_p=1.;
  }

  if(rnd.Rannyu()<metro_p){
    return true;
  }
  else {
    return false;
  }
}

int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p2 >> p1 ;
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


   InitRandom();
   coord v[N];

   ofstream fileout;
   double T:
   double accu=2*M_PI/N;
   for(int i=0; i<N; i++){
     v[i].x=cos(accu*i);
     v[i].y=sin(accu*i);
   }

  fileout.open("cities.txt");
  for(int i=0; i<N; i++){
    fileout<<v[i].x<<" "<<v[i].y<<endl;
  }
  fileout.close();

   double l=0;
   for(int i=0; i<N-1; i++){
     l+=(v[i].x-v[i+1].x)*(v[i].x-v[i+1].x);
     l+=(v[i].y-v[i+1].y)*(v[i].y-v[i+1].y);
   }
   l+=(v[N-1].x-v[0].x)*(v[N-1].x-v[0].x);
   l+=(v[N-1].y-v[0].y)*(v[N-1].y-v[0].y);
   cout<< l<<endl;
   //cout << sumx << " "<<sumy<<" "<<prod<<endl;


   SimulAnn myann,appo;

   myann.SetN(N);
   myann.InitPath();
   myann.CreaPath(v);
   myann.CalcolaLen();

   appo.SetN(N);
   appo.InitPath();
   appo.CopiaPath1(myann);
   appo.CalcolaLen();

   fileout.open("path.txt");
   fileout << myann.GetLen()<<endl;

   double dt=0.01;
   cout<<T/dt<<endl;
   double
   for(int i=0; i<10000; i++){
      if(rnd.Rannyu(0,10)>=9){
         appo.PairPermutation();
       }

       if(rnd.Rannyu(0,10)>=9){
         appo.InvertOrder();
       }

       if(rnd.Rannyu(0,10)>=9){
         appo.Shiftn();
       }

       if(rnd.Rannyu(0,10)>=9){
         appo.Mirror();
       }

       appo.CalcolaLen();
       if(i%100==0 && T!=0){
        //mygen[j+1].m_gen[i].SetT(mygen[j+1].m_gen[i].GetT()-dt);
        //myann[j+1].m_gen[i].SetT(mygen[j+1].m_gen[i].GetT()*0.5);
        T*=0.9;
      }
       if(Metropolis(prob(T,appo),prob(T,myann))==true){
         myapp.CopiaPath1(appo);
         myapp.CalcolaLen();
       }
       fileout<<myann.GetLen()<<endl;
     }

     //cout<<mygen[j+1].m_gen[0].GetT()<<endl;
   fileout.close();
  fileout.open("bestpath.txt");
  for(int i=0; i<N; i++){
    fileout<<myann.GetCity(i).x<<" "<<myann.GetCity(i).y<<endl;
  }
  fileout<< myann.GetCity(0).x<<" "<<myann.GetCity(0).y<<endl;
  fileout.close();













//CITTàààààààààààààààààààà RANDOMICHEEEEEEEEEEEE




  for(int i=0; i<N; i++){
    //double a=rnd.Rannyu(0,2*M_PI);
    v[i].x=rnd.Rannyu();
    sumx+=v[i].x;
    v[i].y=rnd.Rannyu();
    sumy+=v[i].y;
    prod+=v[i].x*v[i].y;
  }

/*
  for(int i=0; i<N ;i++){
    cout <<v[i].x <<"  "<<v[i].y<<endl;
  }
*/


 fileout.open("citiesrandom.txt");
 for(int i=0; i<N; i++){
   fileout<<v[i].x<<" "<<v[i].y<<endl;
 }
 fileout.close();

  l=0;
  for(int i=0; i<N-1; i++){
    l+=(v[i].x-v[i+1].x)*(v[i].x-v[i+1].x);
    l+=(v[i].y-v[i+1].y)*(v[i].y-v[i+1].y);
  }
  l+=(v[N-1].x-v[0].x)*(v[N-1].x-v[0].x);
  l+=(v[N-1].y-v[0].y)*(v[N-1].y-v[0].y);
  cout<< l<<endl;
  //cout << sumx << " "<<sumy<<" "<<prod<<endl;


  Generation mygen1[GG];

  for(int j=0; j<GG; j++){
    mygen1[j].SetDim(DIM);
    mygen1[j].CreaGen();
    for(int i=0; i<DIM; i++){
      mygen1[j].m_gen[i].SetN(N);
      mygen1[j].m_gen[i].InitPath();
      //mygen[j].m_gen[i].CreaPath(v);
      //mygen[j].m_gen[i].CalcolaLen();
    }
  }

  for(int i=0; i<DIM; i++){
    mygen1[0].m_gen[i].CreaPath(v);
    mygen1[0].m_gen[i].CalcolaLen();
    //cout <<mygen[0].m_gen[i].GetLen()<<endl;
  }

  mygen1[0].merge_sort(0,DIM-1);
  fileout.open("bestpathgenrandom.txt");
  fileout << mygen1[0].m_gen[0].GetLen()<<endl;
  avelen.open("avebestpathrandom.txt");
  ave=0;
  for(int i=0; i<DIM/2; i++){
    ave+=mygen1[0].m_gen[i].GetLen();
  }
  ave/=DIM/2;
  avelen<<ave<<endl;

  for(int j=0; j<GG-1; j++){
    for(int i=0; i<DIM/2; i++){
      mygen1[j].Crossover(mygen1[j+1],i);
      //mygen[j+1].m_gen[i].CalcolaLen();
    }

    if(DIM%2==1){
      mygen1[j].m_gen[DIM-1].CopiaPath(mygen1[j+1].m_gen[DIM-1]);
    }


    for(int i=0; i<DIM; i++){
      if(rnd.Rannyu(0,10)>=9){
        mygen1[j+1].m_gen[i].PairPermutation();
       // mygen[j+1].m_gen[i].PrintPath();
      }
      else {
        //mygen[j+1].m_gen[i].CopiaPath(mygen[j+1].m_gen[i]);
        //mygen[j+1].m_gen[i].PrintPath();
      }

      if(rnd.Rannyu(0,10)>=9){
        mygen1[j+1].m_gen[i].InvertOrder();
      }
      else {
        //mygen[j+1].m_gen[i].CopiaPath(mygen[j+1].m_gen[i]);
      }

      if(rnd.Rannyu(0,10)>=9){
        mygen1[j+1].m_gen[i].Shiftn();
      }

     if(rnd.Rannyu(0,10)>=9){
       mygen1[j+1].m_gen[i].Mirror();
     }

    }

    for(int i=0; i<DIM; i++){
      mygen1[j+1].m_gen[i].CalcolaLen();
    }

    mygen1[j+1].merge_sort(0,DIM-1);
    double ave=0;
    for(int i=0; i<DIM/2; i++){
      ave+=mygen1[j+1].m_gen[i].GetLen();
    }
    ave/=DIM/2;
    avelen<<ave<<endl;
    //cout << mygen[j+1].m_gen[0].GetLen()<<endl;
    //cout<< mygen.m_gen[1].GetLen()<<endl;
    fileout<<mygen1[j+1].m_gen[0].GetLen()<<endl;
  }
  fileout.close();
  avelen.close();


/*
double avel=0;
  for(int i=0; i<GG; i++){
    for(int j=0; j<DIM; j++){
    avel+=mygen[i].m_gen[j].GetLen();
   }
   cout << avel/DIM<<endl;
   avel=0;
 }
*/
  //mygen.merge_sort(0,DIM-1);

  cout << "culo"<<endl;
  cout << mygen1[GG-1].m_gen[0].GetLen()<<endl;
  cout << mygen1[GG-1].m_gen[1].GetLen()<<endl;
 cout << mygen1[GG-1].m_gen[2].GetLen()<<endl;
   cout << mygen1[GG-1].m_gen[3].GetLen()<<endl;

 //mygen.Crossover();
 mygen1[GG-1].m_gen[0].PrintPath();
 mygen1[GG-1].m_gen[1].PrintPath();

 fileout.open("bestpathrandom.txt");
 for(int i=0; i<N; i++){
   fileout<<mygen1[GG-1].m_gen[0].GetCity(i).x<<" "<<mygen1[GG-1].m_gen[0].GetCity(i).y<<endl;
 }
 fileout<< mygen1[GG-1].m_gen[0].GetCity(0).x<<" "<<mygen1[GG-1].m_gen[0].GetCity(0).y<<endl;
 fileout.close();

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

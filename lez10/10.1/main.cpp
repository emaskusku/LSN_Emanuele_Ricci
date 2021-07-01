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
#define NN 100000

using namespace std;

double prob(double T, SimulAnn p){
  return exp(-p.GetLen()/T);
}

bool Metropolis(double p1, double p2, double rannyu){
  double metro_p=p1/p2;

  if(metro_p>1){
    metro_p=1.;
  }

  if(rannyu<metro_p){
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

   ofstream fileout, tout, t_acc_tot;
   double accu=2*M_PI/N;
   double r=5.;
   for(int i=0; i<N; i++){
     v[i].x=r*cos(accu*i);
     v[i].y=r*sin(accu*i);
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
   //myann.CreaPerfect(v);
   myann.CalcolaLen();

   appo.SetN(N);
   appo.InitPath();
   appo.CopiaPath1(myann);
   appo.CalcolaLen();

   fileout.open("path.txt");
   fileout << myann.GetLen()<<endl;

   t_acc_tot.open("metropolisAccTot.txt");
   tout.open("temp.txt");
   double T=200;
   tout<<T<<endl;

   int tot=0, acc=0;
   double dt=T/NN;
   cout<<T/dt<<endl;
   for(int i=1; i<NN+1; i++){
     //for(int j=0; j<5; j++){
       if(rnd.Rannyu(0,10)>=5){
         appo.PairPermutation();
         appo.CalcolaLen();
         if(Metropolis(prob(T,appo),prob(T,myann),rnd.Rannyu())==true){
           myann.CopiaPath1(appo);
           myann.CalcolaLen();
           acc++;
         } else {
           appo.CopiaPath1(myann);
         }

         tot++;
       }

       if(rnd.Rannyu(0,10)>=5){
         appo.InvertOrder();
         appo.CalcolaLen();
         if(Metropolis(prob(T,appo),prob(T,myann),rnd.Rannyu())==true){
           myann.CopiaPath1(appo);
           myann.CalcolaLen();
           acc++;
         } else {
           appo.CopiaPath1(myann);
         }
         tot++;
       }

       if(rnd.Rannyu(0,10)>=5){
         appo.Shiftn();
         appo.CalcolaLen();
         if(Metropolis(prob(T,appo),prob(T,myann),rnd.Rannyu())==true){
           myann.CopiaPath1(appo);
           myann.CalcolaLen();
           acc++;
         } else {
           appo.CopiaPath1(myann);
         }
         tot++;
       }

       if(rnd.Rannyu(0,10)>=5){
         appo.Mirror();
         appo.CalcolaLen();
         if(Metropolis(prob(T,appo),prob(T,myann),rnd.Rannyu())==true){
           myann.CopiaPath1(appo);
           myann.CalcolaLen();
           acc++;
         } else {
           appo.CopiaPath1(myann);
         }
         tot++;
       }

       //appo.CalcolaLen();
       /*
       if(Metropolis(prob(T,appo),prob(T,myann),rnd.Rannyu())==true){
         myann.CopiaPath1(appo);
         myann.CalcolaLen();
       }
       */

     //}
       if(i%1000==0 && T!=0){
        //mygen[j+1].m_gen[i].SetT(mygen[j+1].m_gen[i].GetT()-dt);
        //myann[j+1].m_gen[i].SetT(mygen[j+1].m_gen[i].GetT()*0.5);
        t_acc_tot<<acc<<" "<<tot<<" "<<T<<endl;
        acc=0;
        tot=0;
        T*=0.93;
      }

       fileout<<myann.GetLen()<<endl;
       tout<<T<<endl;
     }

     //cout<<mygen[j+1].m_gen[0].GetT()<<endl;
   fileout.close();
  fileout.open("bestpath.txt");
  for(int i=0; i<N; i++){
    fileout<<myann.GetCity(i).x<<" "<<myann.GetCity(i).y<<endl;
  }
  fileout<< myann.GetCity(0).x<<" "<<myann.GetCity(0).y<<endl;
  fileout.close();

  //t_acc_tot<<acc<<" "<<tot<<" "<<T<<endl;
  t_acc_tot.close();
  tout.close();







//CITTàààààààààààààààààààà RANDOMICHEEEEEEEEEEEE






for(int i=0; i<N; i++){
  //double a=rnd.Rannyu(0,2*M_PI);
  v[i].x=rnd.Rannyu(0,5);
  v[i].y=rnd.Rannyu(0,5);
}

fileout.open("citiesrandom.txt");
for(int i=0; i<N; i++){
  fileout<<v[i].x<<" "<<v[i].y<<endl;
}
fileout.close();

SimulAnn myann1,appo1;

myann1.SetN(N);
myann1.InitPath();
myann1.CreaPath(v);
//myann.CreaPerfect(v);
myann1.CalcolaLen();

appo1.SetN(N);
appo1.InitPath();
appo1.CopiaPath1(myann1);
appo1.CalcolaLen();

fileout.open("pathrandom.txt");
fileout << myann1.GetLen()<<endl;

t_acc_tot.open("metropolisAccTotrandom.txt");
tout.open("temprandom.txt");
T=75;
tout<<T<<endl;

tot=0;
acc=0;
dt=T/NN;
cout<<T/dt<<endl;
for(int i=1; i<NN+1; i++){
  //for(int j=0; j<5; j++){
    if(rnd.Rannyu(0,10)>=5){
      appo1.PairPermutation();
      appo1.CalcolaLen();
      if(Metropolis(prob(T,appo1),prob(T,myann1),rnd.Rannyu())==true){
        myann1.CopiaPath1(appo1);
        myann1.CalcolaLen();
        acc++;
      } else {
        appo1.CopiaPath1(myann1);
      }
      tot++;
    }

    if(rnd.Rannyu(0,10)>=5){
      appo1.InvertOrder();
      appo1.CalcolaLen();
      if(Metropolis(prob(T,appo1),prob(T,myann1),rnd.Rannyu())==true){
        myann1.CopiaPath1(appo1);
        myann1.CalcolaLen();
        acc++;
      } else {
        appo1.CopiaPath1(myann1);
      }
      tot++;
    }

    if(rnd.Rannyu(0,10)>=5){
      appo1.Shiftn();
      appo1.CalcolaLen();
      if(Metropolis(prob(T,appo1),prob(T,myann1),rnd.Rannyu())==true){
        myann1.CopiaPath1(appo1);
        myann1.CalcolaLen();
        acc++;
      } else {
        appo1.CopiaPath1(myann1);
      }
      tot++;
    }

    if(rnd.Rannyu(0,10)>=5){
      appo1.Mirror();
      appo1.CalcolaLen();
      if(Metropolis(prob(T,appo1),prob(T,myann1),rnd.Rannyu())==true){
        myann1.CopiaPath1(appo1);
        myann1.CalcolaLen();
        acc++;
      } else {
        appo1.CopiaPath1(myann1);
      }
      tot++;
    }

    //appo.CalcolaLen();
    /*
    if(Metropolis(prob(T,appo),prob(T,myann),rnd.Rannyu())==true){
      myann.CopiaPath1(appo);
      myann.CalcolaLen();
    }
    */

  //}
    if(i%1000==0 && T!=0){
     //mygen[j+1].m_gen[i].SetT(mygen[j+1].m_gen[i].GetT()-dt);
     //myann[j+1].m_gen[i].SetT(mygen[j+1].m_gen[i].GetT()*0.5);
     t_acc_tot<<acc<<" "<<tot<<" "<<T<<endl;
     acc=0;
     tot=0;
     T*=0.93;
   }

    fileout<<myann1.GetLen()<<endl;
    tout<<T<<endl;
  }

  //cout<<mygen[j+1].m_gen[0].GetT()<<endl;
fileout.close();
fileout.open("bestpathrandom.txt");
for(int i=0; i<N; i++){
 fileout<<myann1.GetCity(i).x<<" "<<myann1.GetCity(i).y<<endl;
}
fileout<< myann1.GetCity(0).x<<" "<<myann1.GetCity(0).y<<endl;
fileout.close();

//t_acc_tot<<acc<<" "<<tot<<" "<<T<<endl;
t_acc_tot.close();

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

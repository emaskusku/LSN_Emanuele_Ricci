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
#include "mpi.h"
#include "popolazione.h"

#define N 32
#define DIM 1500
#define GG 1000



//importante da fare: togliere la paralelizzazione del random del main e metterla
// nella "popolazione", quindi cambiare anche il makefile

using namespace std;

int main (int argc, char *argv[]){

  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat1, stat2;
  MPI_Request req;
  string c= to_string(rank);
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

   InitRandom(rank);
   coord v[N];

   ofstream fileout;
   ofstream avelen;
   double sumx=0;
   double sumy=0;
   double prod=0;
   int itag=1; int itag2=2;

   double accu=2*M_PI/N;
   for(int i=0; i<N; i++){
     //double a=rnd.Rannyu(0,2*M_PI);
     v[i].x=cos(accu*i);
     sumx+=v[i].x;
     v[i].y=sin(accu*i);
     sumy+=v[i].y;
     prod+=v[i].x*v[i].y;
   }

/*
   for(int i=0; i<N ;i++){
     cout <<v[i].x <<"  "<<v[i].y<<endl;
   }
*/


  fileout.open("cities_"+c+".txt");
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


   Generation mygen[GG];

   for(int j=0; j<GG; j++){
     mygen[j].SetDim(DIM);
     mygen[j].CreaGen();
     for(int i=0; i<DIM; i++){
       mygen[j].m_gen[i].SetN(N);
       mygen[j].m_gen[i].InitPath();
       //mygen[j].m_gen[i].CreaPath(v);
       //mygen[j].m_gen[i].CalcolaLen();
     }
   }

   for(int i=0; i<DIM; i++){
     mygen[0].m_gen[i].CreaPath(v);
     mygen[0].m_gen[i].CalcolaLen();
     //cout <<mygen[0].m_gen[i].GetLen()<<endl;
   }

   mygen[0].merge_sort(0,DIM-1);
   fileout.open("bestpathgen_"+c+".txt");
   fileout << mygen[0].m_gen[0].GetLen()<<endl;
   avelen.open("avebestpath_"+c+".txt");
   double ave=0;
   for(int i=0; i<DIM/2; i++){
     ave+=mygen[0].m_gen[i].GetLen();
   }
   ave/=DIM/2;
   avelen<<ave<<endl;

   for(int j=0; j<GG-1; j++){
     for(int i=0; i<DIM/2; i++){
       mygen[j].Crossover(mygen[j+1],i);
       //mygen[j+1].m_gen[i].CalcolaLen();
     }

     if(DIM%2==1){
       mygen[j].m_gen[DIM-1].CopiaPath(mygen[j+1].m_gen[DIM-1]);
     }


     for(int i=0; i<DIM; i++){
       if(rnd.Rannyu(0,10)>=9){
         mygen[j+1].m_gen[i].PairPermutation();
        // mygen[j+1].m_gen[i].PrintPath();
       }
       else {
         //mygen[j+1].m_gen[i].CopiaPath(mygen[j+1].m_gen[i]);
         //mygen[j+1].m_gen[i].PrintPath();
       }

       if(rnd.Rannyu(0,10)>=9){
         mygen[j+1].m_gen[i].InvertOrder();
       }
       else {
         //mygen[j+1].m_gen[i].CopiaPath(mygen[j+1].m_gen[i]);
       }

       if(rnd.Rannyu(0,10)>=9){
         mygen[j+1].m_gen[i].Shiftn();
       }

      if(rnd.Rannyu(0,10)>=9){
        mygen[j+1].m_gen[i].Mirror();
      }

     }

     for(int i=0; i<DIM; i++){
       mygen[j+1].m_gen[i].CalcolaLen();
     }

     mygen[j+1].merge_sort(0,DIM-1);

     int ind;
      mygen[j+1].m_gen[0].CreaVettori();
      ind=(int)rnd.Rannyu(0,DIM);
      /*
      for(int kk=0; kk<N; kk++){
        ll[kk]=mygen[j+1].m_gen[0].m_x[0];
        mm[kk]=mygen[j+1].m_gen[ind].m_x[0];
      }
      */
     while(ind==0 || ind==DIM){
       ind=(int)rnd.Rannyu(0,DIM);
     }
     //METTI LA paralelizzazione QUI
     if(rank==1){
       //MPI_Isend(ll,n,MPI_DOUBLE_PRECISION,0,itag, MPI_COMM_WORLD,&req);
       MPI_Isend(&mygen[j+1].m_gen[0].m_x[0],N,MPI_DOUBLE_PRECISION,0,itag, MPI_COMM_WORLD,&req);
       MPI_Isend(&mygen[j+1].m_gen[0].m_y[0],N,MPI_DOUBLE_PRECISION,0,itag, MPI_COMM_WORLD,&req);
       MPI_Recv(&mygen[j+1].m_gen[ind].m_x[0],N,MPI_DOUBLE_PRECISION,0,itag2, MPI_COMM_WORLD,&stat2);
       MPI_Recv(&mygen[j+1].m_gen[ind].m_y[0],N,MPI_DOUBLE_PRECISION,0,itag2, MPI_COMM_WORLD,&stat2);
       mygen[j+1].m_gen[ind].VettoriToPath();
     } else if(rank==0){
        MPI_Send(&mygen[j+1].m_gen[0].m_x[0],N, MPI_DOUBLE_PRECISION,1,itag2,MPI_COMM_WORLD);
        MPI_Send(&mygen[j+1].m_gen[0].m_y[0],N, MPI_DOUBLE_PRECISION,1,itag2,MPI_COMM_WORLD);
        //MPI_Recv(mm,n,MPI_DOUBLE_PRECISION,1,itag,MPI_COMM_WORLD, &stat1);
        MPI_Recv(&mygen[j+1].m_gen[ind].m_x[0],N,MPI_DOUBLE_PRECISION,1,itag,MPI_COMM_WORLD, &stat1);
        MPI_Recv(&mygen[j+1].m_gen[ind].m_y[0],N,MPI_DOUBLE_PRECISION,1,itag,MPI_COMM_WORLD, &stat1);
        mygen[j+1].m_gen[ind].VettoriToPath();
      }

      mygen[j+1].m_gen[ind].CalcolaLen();
      mygen[j+1].merge_sort(0,DIM-1);


     double ave=0;
     for(int i=0; i<DIM/2; i++){
       ave+=mygen[j+1].m_gen[i].GetLen();
     }
     ave/=DIM/2;
     avelen<<ave<<endl;
     //cout << mygen[j+1].m_gen[0].GetLen()<<endl;
     //cout<< mygen.m_gen[1].GetLen()<<endl;
     fileout<<mygen[j+1].m_gen[0].GetLen()<<endl;
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


   //cout << mygen[GG-1].m_gen[0].GetLen()<<endl;
  // cout << mygen[GG-1].m_gen[1].GetLen()<<endl;
  //cout << mygen[GG-1].m_gen[2].GetLen()<<endl;
    //cout << mygen[GG-1].m_gen[3].GetLen()<<endl;

  //mygen.Crossover();
  //mygen[GG-1].m_gen[0].PrintPath();
  //mygen[GG-1].m_gen[1].PrintPath();

  fileout.open("bestpath_"+c+".txt");
  for(int i=0; i<N; i++){
    fileout<<mygen[GG-1].m_gen[0].GetCity(i).x<<" "<<mygen[GG-1].m_gen[0].GetCity(i).y<<endl;
  }
  fileout<< mygen[GG-1].m_gen[0].GetCity(0).x<<" "<<mygen[GG-1].m_gen[0].GetCity(0).y<<endl;
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


 fileout.open("citiesrandom_"+c+".txt");
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
  fileout.open("bestpathgenrandom_"+c+".txt");
  fileout << mygen1[0].m_gen[0].GetLen()<<endl;
  avelen.open("avebestpathrandom_"+c+".txt");
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


    int ind;
     mygen1[j+1].m_gen[0].CreaVettori();
     ind=(int)rnd.Rannyu(0,DIM);
     /*
     for(int kk=0; kk<N; kk++){
       ll[kk]=mygen[j+1].m_gen[0].m_x[0];
       mm[kk]=mygen[j+1].m_gen[ind].m_x[0];
     }
     */
    while(ind==0 || ind==DIM){
      ind=(int)rnd.Rannyu(0,DIM);
    }
    //METTI LA paralelizzazione QUI
    if(rank==1){
      //MPI_Isend(ll,n,MPI_DOUBLE_PRECISION,0,itag, MPI_COMM_WORLD,&req);
      MPI_Isend(&mygen1[j+1].m_gen[0].m_x[0],N,MPI_DOUBLE_PRECISION,0,itag, MPI_COMM_WORLD,&req);
      MPI_Isend(&mygen1[j+1].m_gen[0].m_y[0],N,MPI_DOUBLE_PRECISION,0,itag, MPI_COMM_WORLD,&req);
      MPI_Recv(&mygen1[j+1].m_gen[ind].m_x[0],N,MPI_DOUBLE_PRECISION,0,itag2, MPI_COMM_WORLD,&stat2);
      MPI_Recv(&mygen1[j+1].m_gen[ind].m_y[0],N,MPI_DOUBLE_PRECISION,0,itag2, MPI_COMM_WORLD,&stat2);
      mygen1[j+1].m_gen[ind].VettoriToPath();
    } else if(rank==0){
       MPI_Send(&mygen1[j+1].m_gen[0].m_x[0],N, MPI_DOUBLE_PRECISION,1,itag2,MPI_COMM_WORLD);
       MPI_Send(&mygen1[j+1].m_gen[0].m_y[0],N, MPI_DOUBLE_PRECISION,1,itag2,MPI_COMM_WORLD);
       //MPI_Recv(mm,n,MPI_DOUBLE_PRECISION,1,itag,MPI_COMM_WORLD, &stat1);
       MPI_Recv(&mygen1[j+1].m_gen[ind].m_x[0],N,MPI_DOUBLE_PRECISION,1,itag,MPI_COMM_WORLD, &stat1);
       MPI_Recv(&mygen1[j+1].m_gen[ind].m_y[0],N,MPI_DOUBLE_PRECISION,1,itag,MPI_COMM_WORLD, &stat1);
       mygen1[j+1].m_gen[ind].VettoriToPath();
     }

     mygen1[j+1].m_gen[ind].CalcolaLen();
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

  //cout << "culo"<<endl;
  //cout << mygen1[GG-1].m_gen[0].GetLen()<<endl;
  //cout << mygen1[GG-1].m_gen[1].GetLen()<<endl;
 //cout << mygen1[GG-1].m_gen[2].GetLen()<<endl;
   //cout << mygen1[GG-1].m_gen[3].GetLen()<<endl;

 //mygen.Crossover();
 //mygen1[GG-1].m_gen[0].PrintPath();
 //mygen1[GG-1].m_gen[1].PrintPath();

 fileout.open("bestpathrandom_"+c+".txt");
 for(int i=0; i<N; i++){
   fileout<<mygen1[GG-1].m_gen[0].GetCity(i).x<<" "<<mygen1[GG-1].m_gen[0].GetCity(i).y<<endl;
 }
 fileout<< mygen1[GG-1].m_gen[0].GetCity(0).x<<" "<<mygen1[GG-1].m_gen[0].GetCity(0).y<<endl;
 fileout.close();


   MPI_Finalize();
   //rnd.SaveSeed();
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

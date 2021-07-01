/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "popolazione.h"
//#include "random.h"
#include <cmath>
#include <fstream>

using namespace std;

Random rnd;

void InitRandom(int rank){

  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if(rank==1){
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
 }else{
  if (Primes.is_open()){
     Primes >> p2 >> p1 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
 }
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

  /*
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
  */
}


Path::~Path(){}; //distruttore

void Path::InitRandom(){
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
}

void Path::CreaPath(coord v[]) {//va bene ma da testare con numeri più grandi
  m_path[0]=v[0];
  coord w[m_N];
  int j=1, *index;
  index = new int;
  *index=0;
  for(int i=0; i<m_N; i++){
    w[i].x=v[i].x;
    w[i].y=v[i].y;
  }
  w[0].x=0;
  w[0].y=0;

  while(j<(m_N)){
    while(w[*index].x==0 && w[*index].y==0){
       *index= (int)rnd.Rannyu(1,m_N);
      //cout<<rnd.Rannyu(1,m_N)<<endl;
    }
    m_path[j]=w[*index];
    m_x[j]=w[*index].x;
    m_y[j]=w[*index].y;
    w[*index].x=0;
    w[*index].y=0;
    j++;
  }

  delete index;
}

void Path::Controllo(){
  double sumx=0;
  double sumy=0;
  double prod=0;

  for(int i=0; i<m_N; i++){
    sumx+=m_path[i].x;
    sumy+=m_path[i].y;
    prod+=m_path[i].x*m_path[i].y;

  }
  cout<< sumx<< " "<<sumy<<" "<<prod<<endl;
}


void Path::CalcolaLen(){
  double l=0;
  for(int i=0; i<m_N-1; i++){
    l+=(m_path[i].x-m_path[i+1].x)*(m_path[i].x-m_path[i+1].x);
    l+=(m_path[i].y-m_path[i+1].y)*(m_path[i].y-m_path[i+1].y);
  }
  l+=(m_path[m_N-1].x-m_path[0].x)*(m_path[m_N-1].x-m_path[0].x);
  l+=(m_path[m_N-1].y-m_path[0].y)*(m_path[m_N-1].y-m_path[0].y);
  m_len=l;
}


void Path::Mirror(){
  int n=rnd.Rannyu(1,(m_N)/2+m_N%2);
  coord appo=m_path[n];
  m_path[n]=m_path[m_N-n];
  m_path[m_N-n]=appo;
}

void Path::Shiftn(){
  int m=rnd.Rannyu(0,m_N);
  coord appo[m_N-1];
  int accu;

  for(int i=0; i<m_N-1; i++){
    appo[i].x=m_path[i+1].x;
    appo[i].y=m_path[i+1].y;
  }

  for(int i=1; i<m_N; i++){
    accu=i+m;
    if(accu>m_N-1){
      accu-=(m_N-1);
    }
    m_path[i]=appo[accu-1];
  }
}

void Path::PairPermutation(){
  int a=(int)rnd.Rannyu(1,m_N);
  int b=(int)rnd.Rannyu(1,m_N);

  while (a==b){
    b=(int)rnd.Rannyu(1,m_N);
  }

/*
  for(int i=0; i<m_N; i++){
    p.SetCity(GetCity(i),i);
  }

  p.SetCity(m_path[a],b);
  p.SetCity(m_path[b],a);
*/
  coord appo=m_path[a];
  m_path[a]=m_path[b];
  m_path[b]=appo;
}


void Path::InvertOrder(){
  int a=(int)rnd.Rannyu(1,m_N);
  int b=(int)rnd.Rannyu(1,m_N);

  while (a==b){
    b=(int)rnd.Rannyu(1,m_N);
  }

  if(b<a){
    int c=a;
    a=b;
    b=c;
  }

  /*
  for(int i=0; i<m_N; i++){
    pp.SetCity(m_path[i],i);
  }

  for(int i=0; i<(sqrt((a-b)*(a-b))+1)/2; i++){
    pp.SetCity(m_path[a+i],b-i);
    pp.SetCity(m_path[b-i],a+i);
  }
  */

  coord appo[m_N];
  for(int i=0; i<m_N; i++){
    appo[i]=m_path[i];
  }

  for(int i=a; i<=(b-a)/2; i++){
    m_path[a+i]=appo[b-i];
    m_path[b-i]=appo[a+i];
  }

}


void Path::CopiaPath(Path p){
  for(int i=0; i<m_N; i++){
    p.SetCity(m_path[i],i);
  }
}

void Path::PrintPath (){
  for(int i=0; i<m_N; i++){
    cout << "Città "<< i+1<<":"<<endl;
    cout<< "x: "<< m_path[i].x<<"   y: "<<m_path[i].y<<endl;
  }
}


void Path::CreaVettori(){
  for(int i=0; i<m_N; i++){
    m_x[i]=m_path[i].x;
    m_y[i]=m_path[i].y;
  }
}

void Path::VettoriToPath(){
  for(int i=0; i<m_N; i++){
    m_path[i].x=m_x[i];
    m_path[i].y=m_y[i];
  }
}






void Generation::InitRandom(){
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
}





int Generation::Selector(){//p=1.7 mi pare buono
  double a=pow(rnd.Rannyu(),m_p);
  return (int)m_dim*a;
}

  //Questa procedura implementa il mergeSort (per interi)

void Generation::merge_sort(int low, int high){//Da ricontrollare
      int mid;
      if(low<high) {
          mid = low + (high-low)/2; //This avoids overflow when low, high are too large
          merge_sort(low,mid);
          merge_sort(mid+1,high);
          merge(low,mid,high);
      }
}


void Generation::merge(int low,int mid,int high){
  //low: indice piu` piccolo dell'array da fondere
  //mid: mezzo dell'array da fondere
  //high: indice piu` grande dell'array da fondere

      int h,i,j,k;
      Path b[m_dim]; //Abbiamo bisogno di un vettore di appoggio dove copiare i dati
      h=low; //Indice libero per scandire vettore di sinistra
      i=low; //Indica la prima posizione libera dell'array in cui stiamo copiando i dati
      j=mid+1; //Indice libero per scandire vettore di destra

      while((h<=mid)&&(j<=high)){ //Mentre non hai esaurito uno dei due sottovettori
          if(m_gen[h].GetLen()<=m_gen[j].GetLen()){
              b[i]=m_gen[h];
              h++; //Avanza di uno con l'indice libero dell'array di sinistra
          }
          else{
              b[i]=m_gen[j];
              j++; //Avanza di uno con l'indice libero dell'array di destra
          }
          i++; //
      }

      if(h>mid){ //Se hai esaurito il sottovettore di sinistra, riversa in b quanto rimane del sottovettore di destra
          for(k=j;k<=high;k++){
              b[i]=m_gen[k];
              i++;
          }
      }
      else{
          for(k=h;k<=mid;k++){ //Se hai esaurito il sottovettore di destra, riversa in b quanto rimane del sottovettore di sinistra
              b[i]=m_gen[k];
              i++;
          }
      }
      for(k=low;k<=high;k++) m_gen[k]=b[k]; //Ricopia in a il vettore ordinato b. Questo passo si puo` ottimizzare usando vettori allocati dinamicamente.
  }


void Generation::Crossover(Generation g, int pos){

  int a=Selector();
  int b=Selector();

  while (a==b){
    b=Selector();
  }

  if (rnd.Rannyu(0,9)>4){

  bool flag;
  int k=m_gen[0].GetN()/2+1;


  //coord appo[m_gen[0].GetN()];
  //coord appo1[m_gen[0].GetN()];

  for(int i=0; i<m_gen[0].GetN(); i++){
    g.m_gen[pos*2].SetCity(m_gen[a].GetCity(i),i);
    //perc1[i].y=m_gen[a].GetCity(i).y;
    //perc2[i].x=m_gen[b].GetCity(i).x;
    //perc2[i].y=m_gen[b].GetCity(i).y;
    g.m_gen[pos*2+1].SetCity(m_gen[b].GetCity(i),i);
  }

  //m_gen[a].PrintPath();
  //m_gen[b].PrintPath();

  for(int i=1; i<m_gen[0].GetN(); i++){
    flag=true;
    for(int j=1; j<=m_gen[0].GetN()/2; j++){
      if(g.m_gen[pos*2].GetCity(j).x==m_gen[b].GetCity(i).x && g.m_gen[pos*2].GetCity(j).y==m_gen[b].GetCity(i).y){
        flag=false;
        break;
      }
    }
    if(flag==true){
      g.m_gen[pos*2].SetCity(m_gen[b].GetCity(i),k);
      k++;
    }
  }

  k=m_gen[0].GetN()/2+1;

  for(int i=1; i<m_gen[0].GetN(); i++){
    flag=true;
    for(int j=1; j<=m_gen[0].GetN()/2; j++){
      if(g.m_gen[pos*2+1].GetCity(j).x==m_gen[a].GetCity(i).x && g.m_gen[pos*2+1].GetCity(j).y==m_gen[a].GetCity(i).y){
        flag=false;
        break;
      }
    }
    if(flag==true){
      g.m_gen[pos*2+1].SetCity(m_gen[a].GetCity(i),k);
      k++;
    }
  }

  //m_gen[a].PrintPath();
  //m_gen[b].PrintPath();
  //return
}
else{
  m_gen[a].CopiaPath(g.m_gen[pos*2]);
  m_gen[b].CopiaPath(g.m_gen[pos*2+1]);
}
}

bool Generation::provacrossover(int a, int b){
  for(int i=0; i<m_gen[0].GetN()/2+m_gen[0].GetN()%2; i++){
    for(int j=0; j<m_gen[0].GetN()/2+m_gen[0].GetN()%2; j++){
      if(m_gen[a].GetCity(i).x==m_gen[b].GetCity(j).x && m_gen[a].GetCity(i).y==m_gen[b].GetCity(j).y){
        return false;
      }
    }
  }
  return true;
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

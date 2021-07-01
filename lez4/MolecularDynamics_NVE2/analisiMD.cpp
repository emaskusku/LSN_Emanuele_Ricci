#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <fstream>


using namespace std;


double error(double *av, double *av2, int k){
  if (k==0)
    return 0;
  else{
    return sqrt((av2[k]-av[k]*av[k])/double(k));
  }

}


int main(){

  int M=10000;
  int N=100;
  int L=int (M/N);
  double ave[N];
  double av2[N];
  double sum_prog[N];
  double su2_prog[N];
  double err_prog[N];
  int i,j;
  double sum, accu;
  string end="solid";
  //string end="solid_equilib";
  //double sigma=0.00000000034;


  ofstream fileout;
  ifstream filein;


  filein.open("output_ekin.dat");

  for(i=0; i<N; i++){
    sum=0;
    for(j=0; j<L; j++){
      filein>>accu;
      sum+=accu;
      //cout << sum<<endl;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
    //cout << ave[i]<<"  "<< av2[i] <<endl;
  }

  filein.close();

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

  fileout.open("ave_ekin" + end + ".out");

  for(i=0; i<N; i++)
    fileout<<i*L+L<<" "<< sum_prog[i]<< " " <<err_prog[i]<<endl;

  fileout.close();

/********************************************************/
/********************************************************/

  filein.open("output_epot.dat");


  for(i=0; i<N; i++){
    sum=0;
    for(j=0; j<L; j++){
      filein>>accu;
      sum+=accu;
      //cout << sum<<endl;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
    //cout << ave[i]<<"  "<< av2[i] <<endl;
  }

  filein.close();

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    cout << sum_prog[i]<<endl;
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_epot" + end + ".out");

  for(i=0; i<N; i++)
    fileout<<i*L+L<<" "<< sum_prog[i]<< " " <<err_prog[i]<<endl;

  fileout.close();

  /********************************************************/
  /********************************************************/


  filein.open("output_etot.dat");


  for(i=0; i<N; i++){
    sum=0;
    for(j=0; j<L; j++){
      filein>>accu;
      sum+=accu;
      //cout << sum<<endl;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
    //cout << ave[i]<<"  "<< av2[i] <<endl;
  }

  filein.close();

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    cout << sum_prog[i]<<endl;
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_etot" + end + ".out");

  for(i=0; i<N; i++)
    fileout<<i*L+L<<" "<< sum_prog[i]<< " " <<err_prog[i]<<endl;

  fileout.close();


  /********************************************************/
  /********************************************************/


  filein.open("output_temp.dat");


  for(i=0; i<N; i++){
    sum=0;
    for(j=0; j<L; j++){
      filein>>accu;
      sum+=accu;
      //cout << sum<<endl;
    }
    ave[i]=sum/double(L);
    av2[i]=ave[i]*ave[i];
    //cout << ave[i]<<"  "<< av2[i] <<endl;
  }

  filein.close();

  for(i=0; i<N; i++){
    sum_prog[i]=0;
    su2_prog[i]=0;
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=double(i+1);
    cout << sum_prog[i]<<endl;
    su2_prog[i]/=double(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_temp" + end + ".out");

  for(i=0; i<N; i++)
    fileout<<i*L+L<<" "<< sum_prog[i]<< " " <<err_prog[i]<<endl;

  fileout.close();


  return 0;
}

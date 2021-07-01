#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <fstream>


using namespace std;


float error(float *av, float *av2, int k){
  if (k==0)
    return 0;
  else{
    return sqrt((av2[k]-av[k]*av[k])/float(k));
  }

}


int main(){

  int M=1000;
  int N=100;
  int L=int (M/N);
  float ave[N];
  float av2[N];
  float sum_prog[N];
  float su2_prog[N];
  float err_prog[N];
  int i,j;
  float sum, accu;

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
    ave[i]=sum/float(L);
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
    sum_prog[i]/=float(i+1);
    su2_prog[i]/=float(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_ekin.out");

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
    ave[i]=sum/float(L);
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
    sum_prog[i]/=float(i+1);
    cout << sum_prog[i]<<endl;
    su2_prog[i]/=float(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_epot.out");

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
    ave[i]=sum/float(L);
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
    sum_prog[i]/=float(i+1);
    cout << sum_prog[i]<<endl;
    su2_prog[i]/=float(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_etot.out");

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
    ave[i]=sum/float(L);
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
    sum_prog[i]/=float(i+1);
    cout << sum_prog[i]<<endl;
    su2_prog[i]/=float(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }

  fileout.open("ave_temp.out");

  for(i=0; i<N; i++)
    fileout<<i*L+L<<" "<< sum_prog[i]<< " " <<err_prog[i]<<endl;

  fileout.close();



/*
  for(i=0; i<N; i++){
    fileout<< sum_prog[i]<<endl;
  }

  fileout.close();
  fileout.open("err_prog1.txt");

  for(i=0; i<N; i++){
    fileout<< err_prog[i]<<endl;
  }



  for(i=0; i<N; i++){
    ave[i]=0;
    av2[i]=0;
    sum_prog[i]=0;
    su2_prog[i]=0;
    err_prog[i]=0;
  }

  for(i=0; i<N; i++){
    sum=0;
    for(j=0; j<L; j++){
      sum+=(rnd.Rannyu()-0.5)*(rnd.Rannyu()-0.5);
      //cout << rnd.Rannyu() <<endl;
    }

    ave[i]=sum/float(L);
    av2[i]=ave[i]*ave[i];
    //cout << ave[i]<<"  "<< av2[i] <<endl;
  }


  for(i=0; i<N; i++){
    for(j=0; j<i+1; j++){
      sum_prog[i]+=ave[j];
      su2_prog[i]+=av2[j];
    }
    sum_prog[i]/=float(i+1);
    //cout << sum_prog[i]<<endl;
    su2_prog[i]/=float(i+1);
    err_prog[i]= error(sum_prog, su2_prog, i);
  }


  for(i=0; i<N; i++) {
      int x=L*i+L;
      //cout<<sum_prog[i]<<endl;
      sum_prog[i]=sum_prog[i]-1/12;
      myGraph1.SetPoint(i,x,sum_prog[i]);
      myGraph1.SetPointError(i,0,err_prog[i]);
  }

  */

/*
  fileout.close();
  fileout.open("sum_prog2.txt");

  for(i=0; i<N; i++){
    fileout<< sum_prog[i]<<endl;
  }

  fileout.close();
  fileout.open("err_prog2.txt");

  for(i=0; i<N; i++){
    fileout<< err_prog[i]<<endl;
  }

  fileout.close();

  */



  return 0;
}

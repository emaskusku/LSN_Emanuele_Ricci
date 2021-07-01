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
#include <cmath>
#include "random.h"

using namespace std;


double psiT (double mu, double sigma, double x) {
  double psi=exp( - (x-mu)*(x-mu) / (2.*sigma*sigma) ) + exp( - (x+mu)*(x+mu) / (2.*sigma*sigma) ) ;
  return psi*psi;
}

double H_psi (double mu, double sigma, double x) {

    double psi_T = exp( - (x-mu)*(x-mu) / (2.*sigma*sigma) ) + exp( - (x+mu)*(x+mu) / (2.*sigma*sigma) ) ;
    //double D2_psi_T = -exp( - (x+m_mu)*(x+m_mu) / (2.*m_sigma2) ) * ( exp(2.*m_mu*x/m_sigma2) * ( (m_mu-x)*(m_mu-x) - m_sigma2 ) + (m_mu+x)*(m_mu+x) - m_sigma2 ) / m_sigma2*m_sigma2 ;
    //double D2_psi_T = exp( - (x+m_mu)*(x+m_mu) / (2.*m_sigma2) ) * (x+m_mu)*(x+m_mu)/m_sigma2*m_sigma2 + exp( - (x-m_mu)*(x-m_mu) / (2.*m_sigma2) ) * (x-m_mu)*(x-m_mu)/m_sigma2*m_sigma2 - exp( - (x-m_mu)*(x-m_mu) / (2.*m_sigma2) )/m_sigma2 - exp( - (x+m_mu)*(x+m_mu) / (2.*m_sigma2) )/m_sigma2;
    //double D2_psi_T = exp( - (x-m_mu)*(x-m_mu) / (2.*m_sigma2) ) * (x*x - 2.*m_mu*x + m_mu*m_mu + m_sigma2) /m_sigma2*m_sigma2 + exp( - (x+m_mu)*(x+m_mu) / (2.*m_sigma2) ) * (x*x + 2.*m_mu*x + m_mu*m_mu + m_sigma2) / m_sigma2*m_sigma2 ;
    double add_minus = exp( - (x-mu)*(x-mu) / (2.*sigma*sigma) ) * ( (mu-x)*(mu-x) - sigma*sigma );
    double add_plus = exp( - (x+mu)*(x+mu) / (2.*sigma*sigma) ) * ( (mu+x)*(mu+x) - sigma*sigma );
    double D2_psi_T = ( add_minus + add_plus ) / (sigma*sigma*sigma*sigma);
    //remember we are implying that m=1 and h=2pi
    //using potential: v=x^4-5/2x^2
    double H_psi_T = -0.5*D2_psi_T + ( x*x*x*x - 5./2.*x*x )*psi_T ;
    return H_psi_T/psi_T;
}


double acc(double p){
  if(p>1)
    return 1;
  else
    return p;
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

   rnd.SaveSeed();



   double x_try, x=1;
   double mu=0.2, sigma=1,H;
   double mu_try, sigma_try;
   double D=2.4;
   double F=0.5;
   int L=10000, M=10000, N=100, K=100000, tries=0, accepted=0;
   double A, p_old, p_new, sum=0;
   double E=0.46,E_new;
   double ave[N], av2[N];
   double sum_prog[N], su2_prog[N];
   double err_prog[N];
   ofstream fileout;
   ofstream GS;
   ofstream param;


   for(int i=0; i<M; i++){
     sum=0;
     x=1;
     mu_try=mu+F*rnd.Rannyu(-1,1);
     sigma_try=sigma+F*rnd.Rannyu(-1,1);

     for(int j=0; j<L; j++){
       x_try=x+D*rnd.Rannyu(-1,1);

       p_old=psiT(mu_try, sigma_try, x);
       p_new=psiT(mu_try, sigma_try, x_try);

       A=p_new/p_old;
       A=acc(A);

       if(rnd.Rannyu()<A){
         x = x_try;
       }

       H=H_psi(mu_try, sigma_try, x);
       //ave_H+=H;
       sum+=H;
     }

     E_new=sum/(double) L;

     if(E_new<=E){
       E=E_new;
       mu=mu_try;
       sigma=sigma_try;
     }
     //cout<<"mu: "<<mu << "  sigma: "<<sigma<<"  E: "<<E<<endl;
   }

   cout<<sigma<<endl<<mu<<endl<<E<<endl;
   param.open("parameters.txt");
   param<<sigma<<endl;
   param<<mu<<endl;
   param.close();

   for(int i=0; i<N; i++){
     sum=0;
     x=1;
     for(int j=0; j<L; j++){
       x_try=x+D*rnd.Rannyu(-1,1);

       p_old=psiT(mu, sigma, x);
       p_new=psiT(mu, sigma, x_try);

       A=p_new/p_old;
       A=acc(A);

       if(rnd.Rannyu()<A){
         x = x_try;
         accepted++;
       }

       tries++;
       H=H_psi(mu, sigma, x);
       sum+=H;
       //cout << x[0]<< " "<<x[1]<< " "<<x[2]<<endl;
     }
     ave[i]=sum/double(L);
     av2[i]=ave[i]*ave[i];
   }

   //cout << (double)accepted/tries<<endl;
   //cout << sum/double(L)<<endl;


   for(int i=0; i<N; i++){
     sum_prog[i]=0;
     su2_prog[i]=0;
     for(int j=0; j<i+1; j++){
       sum_prog[i]+=ave[j];
       su2_prog[i]+=av2[j];
     }
     sum_prog[i]/=double(i+1);
     //cout << sum_prog[i]<<endl;
     su2_prog[i]/=double(i+1);
     err_prog[i]= error(sum_prog, su2_prog, i);
   }


   fileout.open("E0.txt");
   for(int i=0; i<N; i++){
   fileout << i*L+L<<" "<< sum_prog[i]<<" "<<err_prog[i]<<endl;
   }
   fileout.close();


GS.open("conf.txt");
for(int i=0; i<K; i++){
  sum=0;
  x=1;
  for(int j=0; j<L; j++){
    x_try=x+D*rnd.Rannyu(-1,1);

    p_old=psiT(mu, sigma, x);
    p_new=psiT(mu, sigma, x_try);

    A=p_new/p_old;
    A=acc(A);

    if(rnd.Rannyu()<A){
      x = x_try;
    }
    //cout << x[0]<< " "<<x[1]<< " "<<x[2]<<endl;
  }
  GS<<x<<endl;
}
GS.close();


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

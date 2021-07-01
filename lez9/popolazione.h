/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __path_h__
#define __path_h__

#include <cmath>
#include <iostream>
#include "random.h"

void InitRandom();

struct coord {
  double x;
  double y;
};

class Path{
public:

  void SetN(int N) { m_N = N; }
  int GetN() {return m_N;}
  void InitPath () { m_path = new coord [m_N]; }
  void InitRandom();
  void CreaPath(coord v[]);
  void CalcolaLen();
  void PairPermutation();
  void InvertOrder();
  double GetLen(){return m_len;}
  void SetCity(coord c, int i){m_path[i]=c;}
  coord GetCity(int pos) {return m_path[pos];}
  void PrintPath();
  void Controllo();
  void CopiaPath (Path p);
  void Shiftn();
  void Mirror();

  ~Path();

private:
  coord *m_path;
  int m_N;
  //Random rnd; // <- se lo fai agire i numeri random sono correlati
  double m_len;

};


class Generation {
public:
  void InitRandom();//va
  void SetDim (int dim) { m_dim=dim; }//va
  int GetDim () {return m_dim;}//va
  int Selector ();//va
  void CreaGen () {m_gen=new Path [m_dim];}
  void merge(int low,int mid,int high);
  void merge_sort(int low, int high); //usa questa
  void Crossover(Generation g, int pos);
  bool provacrossover(int a, int b);
  Path *m_gen;
private:
  int m_dim;
  double m_p=2.;
  //Random rnd; <- se lo fai agire i numeri random sono correlati

};

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

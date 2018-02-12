#ifndef OPTIMISEUR_H
#define OPTIMISEUR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#include "vecteur.hpp"
#include "matrice.hpp"

using namespace std;


class Probleme_abstract
{
public:
  virtual double value(const Vecteur& x) = 0;
  virtual Vecteur gradient(const Vecteur& x) = 0;
};


class Optimiseur
{
private:
  Probleme_abstract* prob_;

public:
  Optimiseur(Probleme_abstract &Pb) {
    prob_ = &Pb;
  };

  double f(const Vecteur& x) {
    return prob_->value(x);
  };

  Vecteur g(const Vecteur& x) {
    return prob_->gradient(x);
  };

  Vecteur gradientFixe(const Vecteur& x0, double eps, int itMax, bool verbose=false);

  Vecteur gradientPasVar(const Vecteur& x0, double eps, int itMax, bool verbose=false);

  Vecteur gradientConj(const Vecteur& x0, double eps, int itMax, bool verbose=false);

  Vecteur quasiNewtonBFGS(const Vecteur& x0, double eps, int itMax, bool verbose=false);

  // Recherche linéaire de Wolfe, au point, gradient et pas initial spécifié
  double rechercheLineaireWolfe(const Vecteur&, const Vecteur&, double, bool);
};

#endif

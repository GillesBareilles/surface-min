#ifndef FCTOPT_H
#define FCTOPT_H

#include <cmath>

#include "vecteur.hpp"
#include "optimiseur.hpp"

class Rosenbrock : public Probleme_abstract
{
public:
  double value(const Vecteur& x) {
    double t1 = (1 - x(1));
    double t2 = (x(2) - x(1) * x(1));
    return t1 * t1 + 100 * t2 * t2;
  };
  Vecteur gradient(const Vecteur& x) {
    Vecteur grad(2, 0.);
    grad(1)  = -2 * (1 - x(1)) + 200 * (x(2) - x(1) * x(1)) * (-2 * x(1));
    grad(2)  =                   200 * (x(2) - x(1) * x(1));
    return grad;
  };
};

class StyblinskiTang : public Probleme_abstract
{
public:
  double value(const Vecteur& x) {
    double r = 0;
    for (int i=1;i<=x.size();i++) r+=pow(x(i),4)-16*pow(x(i),2)+5*x(i);
    return r/2;
  };
  Vecteur gradient(const Vecteur& x) {
    Vecteur grad(x.size());
    for (int i=1;i<=x.size();i++) grad(i) = 2*pow(x(i),3)-16*x(i)+5/2;
    return grad;
  };
};

class Booth : public Probleme_abstract
{
public:
  double value(const Vecteur& x) {
    return pow(x(1)+2*x(2)-7,2) + pow(2*x(1)+x(2)-5,2);
  };
  Vecteur gradient(const Vecteur& x) {
    Vecteur grad(x.size());
    grad(1) = 2*(x(1)+2*x(2)-7) + 4*(2*x(1)+x(2)-5);
    grad(2) = 4*(x(1)+2*x(2)-7) + 2*(2*x(1)+x(2)-5);
    return grad;
  };
};


double hyperplan(const Vecteur &x) {
  return x(1);
}

double unite(const Vecteur &x) {
  return 1;
}

double norme_quad(const Vecteur &x) {
  return pow(x(1),2) + pow(x(2),2);
}

double norme_quad_opp(const Vecteur &x) {
  return pow(x(1),2) - pow(x(2),2);
}

double norme_tri(const Vecteur &x) {
  return pow(2*x(1),3) + pow(2*x(2),3);
}

double linSin(const Vecteur &x) {
  return -x(1) + sin(5*x(2));
}

#endif

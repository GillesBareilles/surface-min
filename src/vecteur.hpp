#ifndef VECTEUR_H
#define VECTEUR_H

#include <iostream>
#include <cmath>
#include <vector>
#include "matrice.hpp"

using namespace std;


class Vecteur : public Matrice {
public:
    Vecteur();
    Vecteur(int, double v=0.);
    Vecteur(const Vecteur&);
    Vecteur(const std::vector<double>&);

    int size() const;
    double operator() (int) const;
    double& operator() (int);
    double operator[] (int) const;
    double& operator[] (int);

    Vecteur& operator+=(double);
    Vecteur& operator-=(double);
    Vecteur& operator*=(double);
    Vecteur& operator/=(double);
    Vecteur& operator+=(const Vecteur&);
    Vecteur& operator-=(const Vecteur&);

    Matrice t();
};

Vecteur operator+(const Vecteur&);                  // +U
Vecteur operator-(const Vecteur&);                  // -U
Vecteur operator+(const Vecteur&, double);          // U+x
Vecteur operator-(const Vecteur&, double);          // U-x
Vecteur operator*(const Vecteur&, double);          // U*x
Vecteur operator/(const Vecteur&, double);          // U/x
Vecteur operator+(double, const Vecteur&);          // x+U
Vecteur operator-(double, const Vecteur&);          // x-U
Vecteur operator*(double, const Vecteur&);          // x*U
Vecteur operator+(const Vecteur&, const Vecteur&);  // U+V
Vecteur operator-(const Vecteur&, const Vecteur&);  // U-V

ostream& operator<<(ostream& os, const Vecteur&);

double ps(const Vecteur &, const Vecteur &);
double norme(const Vecteur &);

Vecteur operator*(const Matrice&, const Vecteur&);
Matrice operator*(const Vecteur&, const Matrice&);
#endif

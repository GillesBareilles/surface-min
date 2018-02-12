#ifndef MATRICE_H
#define MATRICE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class Matrice {
protected:
  int dim_l_, dim_c_;
  vector<double> val_;

public:
  Matrice();
  Matrice (int, int, double v=0.);
  Matrice(const Matrice&);

  Matrice& operator=(const Matrice &);

  double& operator() (int, int);
  double operator() (int, int) const;

  int dim_l() const;
  int dim_c() const;
  const double* data() const;

  Matrice& operator+=(double);
  Matrice& operator-=(double);
  Matrice& operator*=(double);
  Matrice& operator/=(double);
  Matrice& operator+=(const Matrice&);
  Matrice& operator-=(const Matrice&);

  bool operator==(const Matrice &) const;
  bool operator!=(const Matrice &) const;

  Matrice& identite();
  double det() const;
  double norme() const;

  void saveToFile(string fichier = "Vecteur.m");
};

Matrice operator+(const Matrice&);                  // +A
Matrice operator-(const Matrice&);                  // -A
Matrice operator+(const Matrice&, double);          // A+x
Matrice operator-(const Matrice&, double);          // A-x
Matrice operator*(const Matrice&, double);          // A*x
Matrice operator/(const Matrice&, double);          // A/x
Matrice operator+(double, const Matrice&);          // x+A
Matrice operator-(double, const Matrice&);          // x-A
Matrice operator*(double, const Matrice&);          // x*A
Matrice operator+(const Matrice&, const Matrice&);  // A+B
Matrice operator-(const Matrice&, const Matrice&);  // A-B
Matrice operator*(const Matrice&, const Matrice&);  // A*B

ostream& operator<<(ostream& os, const Matrice&);
istream& operator>>(istream& is, Matrice&);

Matrice t(const Matrice&);
Matrice inv(const Matrice&);

#endif

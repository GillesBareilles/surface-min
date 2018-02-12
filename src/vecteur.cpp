#include "vecteur.hpp"
#include <stdlib.h>

Vecteur::Vecteur() : Matrice() {}

Vecteur::Vecteur(int d, double v) : Matrice(d,1,v) {}

Vecteur::Vecteur(const Vecteur& V) : Matrice(V.size(), 1) {
  for (int i=0;i<dim_l_; i++) val_[i] = V.val_[i];
}

Vecteur::Vecteur(const std::vector<double>& V) : Matrice(V.size(), 1) {
  for (int i=0;i<dim_l_; i++) val_[i] = V[i];
}

// Accesseurs
int Vecteur::size() const {
  return dim_l_;
}

double Vecteur::operator() (int i) const {
  return val_.at(i-1);
}

double& Vecteur::operator() (int i) {
  return val_.at(i-1);
}

double Vecteur::operator[] (int i) const {
  return val_.at(i);
}

double& Vecteur::operator[] (int i) {
  return val_.at(i);
}


// Opérateurs algébriques internes
Vecteur& Vecteur::operator+=(double x) {
  double* pU=val_.data();
  for (int i=0;i<dim_l_; i++, pU++) (*pU)+=x;
  return *this;
}

Vecteur& Vecteur::operator-=(double x) {
  double* pU=val_.data();
  for (int i=0;i<dim_l_; i++, pU++) (*pU)-=x;
  return *this;
}

Vecteur& Vecteur::operator*=(double x) {
  double* pU=val_.data();
  for (int i=0;i<dim_l_; i++, pU++) (*pU)*=x;
  return *this;
}

Vecteur& Vecteur::operator/=(double x) {
  if (x==0) {
    cout << "\nDivision par zéro dans V/=x\n";
    exit(-1);
  }
  double* pU=val_.data();
  for (int i=0;i<dim_l_; i++, pU++) (*pU)/=x;
  return *this;
}

Vecteur& Vecteur::operator+=(const Vecteur& V) {
  if (dim_l_!=V.dim_l_) {
    cout << "\nAjout entre vecteurs incompatibles\n";
    exit(-1);
  }
  double* pU=val_.data();
  const double* pV=V.val_.data();
  for (int i=0;i<dim_l_; i++, pU++, pV++) (*pU)+=(*pV);
  return *this;
}

Vecteur& Vecteur::operator-=(const Vecteur& V) {
  if (dim_l_!=V.dim_l_) {
    cout << "Différence entre vecteurs incompatibles\n";
    exit(-1);
  }
  double* pU=val_.data();
  const double* pV=V.val_.data();
  for (int i=0;i<dim_l_; i++, pU++, pV++) (*pU)-=(*pV);
  return *this;
}

// Transposition
Matrice Vecteur::t() {
  dim_c_ = dim_l_;
  dim_l_ = 1;
  Matrice R(*this);
  dim_l_ = dim_c_;
  dim_c_ = 1;
  return R;
}

ostream& operator<<(ostream &os, const Vecteur &V)
{
  os << "[ ";
  for (int i=0;i<V.size();i++) os << V(i+1) << " ";
  os << "]";
  return os;
}

// Fonctions externes
Vecteur operator+(const Vecteur& U) {                    // +U
  return U;}
Vecteur operator-(const Vecteur& U) {                    // -U
  Vecteur R(U.size(), 0.); return R-=U;}

Vecteur operator+(const Vecteur &U, double x) {          // U+x
  Vecteur R(U); return R+=x;}
Vecteur operator-(const Vecteur &U, double x) {          // U-x
  Vecteur R(U); return R-=x;}
Vecteur operator*(const Vecteur &U, double x) {          // U*x
  Vecteur R(U); return R*=x;}
Vecteur operator/(const Vecteur &U, double x) {          // U/x
  Vecteur R(U); return R/=x;}

Vecteur operator+(double x, const Vecteur &U) {          // x+U
  return U+x;}
Vecteur operator-(double x, const Vecteur &U) {          // x-U
  Vecteur R(U.size(), x); return R-=U;}
Vecteur operator*(double x, const Vecteur &U) {          // x*U
  return U*x;}

Vecteur operator+(const Vecteur &U, const Vecteur &V) {  // U+V
  Vecteur R(U); return R+=V;}
Vecteur operator-(const Vecteur &U, const Vecteur &V) {  // U-V
  Vecteur R(U); return R-=V;}


double ps(const Vecteur &V1, const Vecteur &V2) {
  if (V1.size()!=V2.size()) {
    cout << "\nProduit scalaire entre vecteurs incompatibles\n";
    exit(-1);
  }
  double r=0;
  const double* p1=V1.data();
  const double* p2=V2.data();
  for (int i=1;i<=V1.size();i++,p1++,p2++) r += (*p1) * (*p2);
  return r;
}

double norme(const Vecteur &U) {
  return sqrt(ps(U,U));
}

// Produit Matrice*Vecteur -> Vecteur
Vecteur operator*(const Matrice &A, const Vecteur &V) {
    if (A.dim_c()!=V.size()) {
      cout << "\nErreur : operator*(Matrice, Vecteur), dimensions non compatibles\n";
      exit(-1);
    }
    Vecteur R(A.dim_l());
    const double *pA = A.data();
    const double *pV = V.data();
    double *pR = &R(1);
    for (int j=1;j<=A.dim_c();j++,pV++) {
      pR = &R(1);
      for (int i=1;i<=A.dim_l();i++,pA++,pR++) (*pR)+=((*pA)*(*pV));
    }
    return R;
}

// Produit Vecteur*Matrice -> Matrice (matrice 1 ligne, ~ forme linéaire)
Matrice operator*(const Vecteur &V, const Matrice &A) {
    if ((A.dim_l()!=1) and (A.dim_l()!=V.size())) {
      cout << "\nErreur : operator*(Vecteur, Matrice), dimensions non compatibles\n";
      exit(-1);
    }
    Matrice R(V.size(),A.dim_c());
    for (int i=1;i<=V.size();i++) {
      for (int j=1;j<=A.dim_c();j++) {
        R(i,j) = V(i) * A(1,j);
      }
    }
    return R;
}

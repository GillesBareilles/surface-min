#include "matrice.hpp"

Matrice::Matrice () : dim_l_(0), dim_c_(0), val_(vector<double>()) {}

Matrice::Matrice (int nbL, int nbC, double v) : dim_l_(nbL), dim_c_(nbC), val_(vector<double>(dim_l_*dim_c_, v)) {}

Matrice::Matrice(const Matrice &M) : dim_l_(M.dim_l_), dim_c_(M.dim_c_), val_(vector<double>(dim_l_*dim_c_, 0)) {
  int d = dim_l_*dim_c_;
  for (int i=0;i<d;i++) val_[i] = M.val_[i];
}

// Opérateur d'assignation
Matrice& Matrice::operator=(const Matrice& M) {
  if (this==(&M)) return *this;
  dim_l_=M.dim_l_;
  dim_c_=M.dim_c_;
  val_.resize(dim_l_*dim_c_);
  for (int i=0;i<dim_l_*dim_c_; i++) val_[i] = M.val_[i];
  return *this;
}

// Acces a un element
double& Matrice::operator() (int i, int j) {
  return val_.at((j-1)*dim_l_+i-1);
}

double Matrice::operator() (int i, int j) const {
  return val_.at((j-1)*dim_l_+i-1);
}

int Matrice::dim_l() const {
  return dim_l_;
}
int Matrice::dim_c() const {
  return dim_c_;
}
const double* Matrice::data() const {
  return val_.data();
}

// Opérateurs algébriques internes
Matrice& Matrice::operator+=(double x) {
  double *pA=val_.data();
  int d=dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++) (*pA)+=x;
  return *this;
}

Matrice& Matrice::operator-=(double x) {
  double *pA=val_.data();
  int d=dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++) (*pA)-=x;
  return *this;
}

Matrice& Matrice::operator*=(double x) {
  double *pA=val_.data();
  int d=dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++) (*pA)*=x;
  return *this;
}

Matrice& Matrice::operator/=(double x) {
  if (x==0) {
    cout << "\nErreur : operator/=(Matrice, double), division par zéro\n";
    exit(-1);
  }
  double *pA=val_.data();
  int d=dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++) (*pA)/=x;
  return *this;
}

Matrice& Matrice::operator+=(const Matrice& B) {
  if ((dim_l_!=B.dim_l_) or (dim_c_!=B.dim_c_)) {
    cout << "\nErreur : operator+=(Matrice, Matrice), Matrices incompatibles\n";
    exit(-1);
  }
  double *pA=val_.data();
  const double *pB=B.val_.data();
  int d=dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++, pB++) (*pA)+=(*pB);
  return *this;
}

Matrice& Matrice::operator-=(const Matrice& B) {
  if ((dim_l_!=B.dim_l_) or (dim_c_!=B.dim_c_)) {
    cout << "\nErreur : operator-=(Matrice, Matrice), Matrices incompatibles\n";
    exit(-1);
  }
  double *pA=val_.data();
  const double *pB=B.val_.data();
  int d=dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++, pB++) (*pA)-=(*pB);
  return *this;
}

// Opérateur d'égalité
bool Matrice::operator==(const Matrice &B) const {
  if ((dim_l_!=B.dim_l_) or (dim_c_!=B.dim_c_)) return false;
  const double *pA=val_.data();
  const double *pB=B.val_.data();
  int d = dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++, pB++) if (*(pA)!=*(pB)) return false;
  return true;
}

bool Matrice::operator!=(const Matrice &B) const {
  return !((*this)==B);
}

// Définit la matrice courante à l'identité (si bonnes dimensions)
Matrice& Matrice::identite() {
  if (dim_l_!=dim_c_) {
    cout << "\nErreur : identite(), Matrice non carrée\n";
    exit(-1);
  }
  double* pA=val_.data();
  int d = dim_l_*dim_c_;
  for (int i=0;i<d; i++, pA++) *(pA) = 0;
  pA = val_.data();
  for (int i=0;i<dim_l_;i++) (*(pA+(dim_l_+1)*i))=1;
  return *this;
}

// Déterminant des matrices 2x2
double Matrice::det() const {
  if (dim_l_!=dim_c_) {
    cout << "\nErreur : inv(), Matrice non carrée\n";
    exit(-1);
  }
  if (dim_l_!=2) {
    cout << "\nErreur : inv(), Matrice de dim > 2 ...\n";
    exit(-1);
  }
  return (val_[0]*val_[3]-val_[1]*val_[2]);
}

// Transpostion
Matrice t(const Matrice &M) {
  Matrice T(M.dim_c(), M.dim_l());
  const double *pM = M.data();
  double *pT = &T(1,1);
  for (int j=1;j<=T.dim_c();j++) {
    for (int i=1;i<=T.dim_l();i++, pT++) *(pT) = *(pM +(i-1)*M.dim_l()+j-1);
  }
  return T;
}


// Inversion des matrices 2x2
Matrice inv(const Matrice& M) {
  if (M.dim_l()!=M.dim_c()) {
    cout << "\nErreur : inv(), Matrice non carrée\n";
    exit(-1);
  }
  if (M.dim_l()!=2) {
    cout << "\nErreur : inv(), Matrice de dim > 2 ...\n";
    exit(-1);
  }
  if ((M(1,1)*M(2,2)-M(1,2)*M(2,1))==0) {
    cout << "\nErreur : inv(), matrice non inversible ...\n";
    exit(-1);
  }
  Matrice inverse(2,2);
  inverse(1,1) =  M(2,2);
  inverse(2,2) =  M(1,1);
  inverse(1,2) = -M(1,2);
  inverse(2,1) = -M(2,1);

  return inverse/(M(1,1)*M(2,2)-M(1,2)*M(2,1));
}

double Matrice::norme() const {
  double res = 0;
  const double* pU = val_.data();
  for (int i=0;i<dim_l_*dim_c_;i++,pU++) res+=*(pU);
  return sqrt(res);
}

void Matrice::saveToFile(string fichier) {
  fstream fs;
  fs.open (fichier.c_str(), fstream::in | fstream::out | fstream::trunc);

  const double *pU=val_.data();
  for (int i=0;i<dim_l_;i++) {
    for (int j=0;j<dim_c_;j++,pU++) fs << *(pU) << " ";
    fs << endl;
  }

  fs.close();
}


// Fonctions externes à la classe
Matrice operator+(const Matrice& A) {                    // +A
  return A;}
Matrice operator-(const Matrice& A) {                    // -A
  Matrice R(A.dim_l(), A.dim_c(), 0.); return R-=A;}     // Non optimal

Matrice operator+(const Matrice &A, double x) {          // A+x
  Matrice R(A); return R+=x;}
Matrice operator-(const Matrice &A, double x) {          // A-x
  Matrice R(A); return R-=x;}
Matrice operator*(const Matrice &A, double x) {          // A*x
  Matrice R(A); return R*=x;}
Matrice operator/(const Matrice &A, double x) {          // A/x
  Matrice R(A); return R/=x;}

Matrice operator+(double x, const Matrice &A) {          // x+A
  return A+x;}
Matrice operator-(double x, const Matrice &A) {          // x-A
  Matrice R(A.dim_l(), A.dim_c(), x); return R-=A;}
Matrice operator*(double x, const Matrice &A) {          // x*A
  return A*x;}

Matrice operator+(const Matrice &A, const Matrice &B) {  // A+B
  Matrice R(A); return R+=B;}
Matrice operator-(const Matrice &A, const Matrice &B) {  // A-B
  Matrice R(A); return R-=B;}
Matrice operator*(const Matrice &A, const Matrice &B) {  // A*B
  Matrice R(A.dim_l(), B.dim_c());
  for (int i=1;i<=R.dim_l();i++) {
    for (int j=1;j<=R.dim_c();j++) {
      for (int k=1;k<=A.dim_c();k++) R(i,j) += A(i,k)*B(k,j);
    }
  }
  return R;
}

//Operateur de flux
ostream& operator<<(ostream& os, const Matrice &A) {
  os << "Matrice de dimensions (" << A.dim_l() << "," << A.dim_c() << ")" << endl;
  for (int i=1;i<=A.dim_l();i++) {
    for (int j=1;j<=A.dim_c();j++) os << A(i,j) << " ";
    os << endl;
  }
  return os;
}

istream& operator>>(istream& is, Matrice &A) {
  for (int i=1;i<=A.dim_l();i++) {
    for (int j=1;j<=A.dim_c();j++) is >> A(i,j);
  }
  return is;
}

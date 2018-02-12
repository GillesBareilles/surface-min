#include "maillage.hpp"

//---------------------- Maillage
void Maillage::summary()
{
  if (points_.empty()) cout << "Pas de points_ !" << endl;
  else
  {
    cout << "nombre de points_ : " << points_.size() << endl;
    for (Point p:points_) p.print();

    cout << endl << "nombre de segments_ : " << segments_.size() << endl;
    for (Segment s:segments_) s.print();

    cout << endl << "nombre de triangles_ : " << triangles_.size() << endl;
    for (Triangle t:triangles_) t.print();

    cout << endl << "nombre de nodes_ : " << nodes_.size() << endl;
    // for (vector<double> p:nodes_) cout << "x, y, z -- " << p[0] << " " << p[1] << " " << p[2] << endl;

    cout << endl << "maillage d'ordre : " << ordre_elts_geo_ << endl;
  }
}

Maillage::Maillage(string fichier, int ordre, double a, double b)
{
    alpha_ = a;
    beta_ = b;

    //Test de l'existence du fichier .geo
    int res;
    string maCommande = "ls " + fichier;
    const char *cmd = maCommande.c_str();
    res = system(cmd);
    if (res != 0) {
      cout << "chemin specifie pour le fichier .geo invalide" << endl;
      exit(-1);
    }

    // Création du fichier .msh
    maCommande = "gmsh -2 -order " + to_string(ordre) + " " + fichier;
    cout << maCommande << endl << endl;
    const char *cmd1 = maCommande.c_str();
    res = system(cmd1);
    if (res != 0) {
      cout << "l'execution de gmsh sur le fichier fourni a échouée" << endl;
      exit(-1);
    }

    cout << "... done creating " << fichier.substr(0,fichier.length()-4) << ".msh" << endl;
    ifstream geomFile;
    geomFile.open(fichier.substr(0,fichier.length()-4) + ".msh", ifstream::in);

    string buffer;
    for(int k=0;k<4;k++) getline(geomFile,buffer); //skip head of file

    //Récup nb_nodes_
    int nb_nodes_;
    geomFile >> nb_nodes_;

    vector<int> id;

    int i;
    while (geomFile>>i) {
      id.push_back(i);
      double w;
      vector<double> P;
      geomFile>>w;
      P.push_back(w);
      geomFile>>w;
      P.push_back(w);
      geomFile>>w;
      P.push_back(w);
      nodes_.push_back(P);
    }

    // cout << "nodes_ size : " << nodes_.size() << endl;
    cout << "Node parsing done" << endl;


    //String de nouveau valide ; saut de '$Endnodes_\n$Elements'
    geomFile.clear();
    for(int k=0;k<2;k++) getline(geomFile,buffer);


    int nb_elts;
    geomFile >> nb_elts;
    // cout << "Parsing "<<nb_elts<<" elements"<<endl;

    id.clear();
    int type_elem, nb_flags, flag;
    int idNode1, idNode2, idNode3, idNode4, idNode5, idNode6;

    // Tant que le prochain elt du flux est un int...
    // /!/ Les indices gmsh commencent à 1, sont donc décrémentés
    while (geomFile >> i) {
      id.push_back(i);
      geomFile >> type_elem;
      geomFile >> nb_flags;
      for (int i=0; i<nb_flags; i++) geomFile >> flag; //Flags inutiles
      if (type_elem==15) { //points_
        geomFile >> idNode1;
        idNode1--;
        points_.push_back(Point(idNode1));
      }
      if (type_elem==1) { //segments_ P1
        geomFile >> idNode1 >> idNode2;
        idNode1--; idNode2--;
        segments_.push_back(Segment(idNode1, idNode2));

        // insertion des points dans le bord...
        insereDansTrie(points_bord_, idNode1);
        insereDansTrie(points_bord_, idNode2);
      }
      if (type_elem==2) { //triangles_ P1
        geomFile >> idNode1 >> idNode2 >> idNode3;
        idNode1--; idNode2--; idNode3--;
        triangles_.push_back(Triangle(idNode1, idNode2, idNode3));
      }
      if (type_elem==8) { //segments_ P2
        geomFile >> idNode1 >> idNode2 >> idNode3;
        idNode1--; idNode2--; idNode3--;
        segments_.push_back(Segment(idNode1, idNode2, idNode3));

        // insertion des points dans le bord...
        insereDansTrie(points_bord_, idNode1);
        insereDansTrie(points_bord_, idNode2);
        insereDansTrie(points_bord_, idNode3);
      }
      if (type_elem==9) { //Triangle P2
        geomFile >> idNode1 >> idNode2 >> idNode3 >> idNode4 >> idNode5 >> idNode6;
        idNode1--; idNode2--; idNode3--;idNode4--; idNode5--; idNode6--;
        triangles_.push_back(Triangle(idNode1, idNode2, idNode3, idNode4, idNode5, idNode6));
      }
    }

    ordre_elts_geo_ = triangles_[0].ordre();

    cout << "Element parsing done" << endl;
    cout << "... " << points_.size() << " points_" << endl;
    cout << "... " << segments_.size() << " segments_" << endl;
    cout << "... " << triangles_.size() << " triangles_" << endl;
    cout << "... " << nodes_.size() << " nodes_" << endl;
    cout << "... " << points_bord_.size() << " points_bord_" << endl;
    cout << "... " << "Elements géométriques d'ordre " << ordre_elts_geo_ << endl << endl;

    geomFile.close();
}



Point Maillage::points(int i) const {
  return points_[i];
}

Segment Maillage::segments(int i) const {
  return segments_[i];
}

Triangle Maillage::triangles(int i) const {
  return triangles_[i];
}

vector<double> Maillage::nodes(int i) const {
  return nodes_[i];
}

int Maillage::points_bord(int i) const {
  return points_bord_[i];
}

int Maillage::points_size() const {
  return points_.size();
}

int Maillage::segments_size() const {
  return segments_.size();
}

int Maillage::triangles_size() const {
  return triangles_.size();
}

int Maillage::nodes_size() const {
  return nodes_.size();
}

int Maillage::points_bord_size(int) const {
  return points_bord_.size();
}

bool Maillage::estAuBord(int ind) const {
  if (points_bord_.empty()) return false;
  if ((points_bord_.front() == ind) or (points_bord_.back() == ind)) return true;

  int n1 = 0;
  int n2 = points_bord_.size();
  while ((n2-n1)>1) {
    if (points_bord_[(n1+n2)/2] == ind) return true;
    else if (points_bord_[(n1+n2)/2] < ind) n1 = (n1+n2)/2;
    else n2 = (n1+n2)/2;
  }
  return false;
}

double Maillage::mesure(const Forme &T) const {
  vector<double> P0 = nodes_[T[0]];
  vector<double> P1 = nodes_[T[1]];
  vector<double> P2 = nodes_[T[2]];
  return abs( (P1[0]-P0[0])*(P2[1]-P0[1]) - (P1[1]-P0[1])*(P2[0]-P0[0]))/2;
}

vector<Vecteur> Maillage::gradients(const Forme &T) const {
  vector<Vecteur> G;
  Matrice J(2,2);
  J(1,1) =  nodes_[T[1]][0] - nodes_[T[0]][0];
  J(2,1) =  nodes_[T[1]][1] - nodes_[T[0]][1];
  J(1,2) =  nodes_[T[2]][0] - nodes_[T[0]][0];
  J(2,2) =  nodes_[T[2]][1] - nodes_[T[0]][1];
  Matrice invtJ = inv(t(J));
  Vecteur g(2);
  g(1) = -1; g(2) = -1;
  G.push_back(invtJ*g);
  g(1) = 1; g(2) = 0;
  G.push_back(invtJ*g);
  g(1) = 0; g(2) = 1;
  G.push_back(invtJ*g);
  return G;
}

// Calcul de surface explicite avec des polynômes d'ordre 1.
double Maillage::valueExpl(const Vecteur &U) {
  double J = 0;
  for (Triangle T:triangles_) {
    vector<Vecteur> G = this->gradients(T);
    double mesure = this->mesure(T);
    double S = 0;
    Vecteur sum_grad(2);
    for (int k=0;k<3;k++) {
      sum_grad += U[T[k]]*G[k];
    }
    S = ps(sum_grad, sum_grad);
    J += mesure*(alpha_*sqrt(1+S)+beta_*S);
  }
  return J;
}

// Calcul de gradient explicite avec des polynômes d'ordre 1.
Vecteur Maillage::gradientExpl(const Vecteur &U) {
  Vecteur gradient(U.size());

  for (Triangle T:triangles_) {
    vector<Vecteur> G = this->gradients(T);
    double mesure = this->mesure(T);

    Vecteur sum_grad(2);
    for (int k=0;k<3;k++) {
      sum_grad += U[T[k]]*G[k];
    }
    for (int i=0;i<3;i++) {
      double g2 = ps(sum_grad, sum_grad);

      if (!this->estAuBord(T[i])) {
        double prodScal = ps(G[i], sum_grad);
        gradient[T[i]] += mesure*(alpha_*prodScal/sqrt(1+g2)+2*beta_*prodScal);
      }
    }
  }
  return gradient;
}


// Methode de quadrature de gaussLegendre
void Maillage::gaussLegendre(vector<Vecteur> &ptsQuadT, Vecteur &pdsQuadT, vector< vector< Vecteur> > &gradsTRef) const {
  // Calcul des 7 points de quadrature
  double s0 = 1./6;
  double s1 = 2./3;
  int nbq = 3;

  ptsQuadT = vector<Vecteur>();
  Vecteur temp(2);
  temp(1)=s0;temp(2)=s0;
  ptsQuadT.push_back(temp);
  temp(1)=s1;temp(2)=s0;
  ptsQuadT.push_back(temp);
  temp(1)=s0;temp(2)=s1;
  ptsQuadT.push_back(temp);

  // Calcul des poids de quadrature.
  double pds = 1./6;
  pdsQuadT = Vecteur(nbq, pds);

  // Calcul des gradients des 3 fonctions de base (triangulation P1) sur les 3 points de quadrature
  gradsTRef = vector< vector< Vecteur> >(nbq);
  for (int k=0; k<nbq; k++) {
    temp(1) = -1; temp(2) = -1;
    gradsTRef[k].push_back(temp);
    temp(1) = 1; temp(2) = 0;
    gradsTRef[k].push_back(temp);
    temp(1) = 0; temp(2) = 1;
    gradsTRef[k].push_back(temp);
  }
}

// Méthode de quadrature de Gauss Lobatto
void Maillage::gaussLobatto(vector<Vecteur> &ptsQuadT, Vecteur &pdsQuadT, vector< vector< Vecteur> > &gradsTRef) const {
  // Calcul des 7 points de quadrature
  double os = sqrt(15);
  int nbq=7;
  double s3=1/3;
  double pp1=(6-os)/21;   double pp2=(6+os)/21;
  double pp3=(9+2*os)/21; double pp4=(9-2*os)/21;

  ptsQuadT = vector<Vecteur>();
  Vecteur temp(2);
  temp(1)=s3;temp(2)=s3;
  ptsQuadT.push_back(temp);
  temp(1)=pp1;temp(2)=pp1;
  ptsQuadT.push_back(temp);
  temp(1)=pp1;temp(2)=pp3;
  ptsQuadT.push_back(temp);
  temp(1)=pp3;temp(2)=pp1;
  ptsQuadT.push_back(temp);
  temp(1)=pp2;temp(2)=pp2;
  ptsQuadT.push_back(temp);
  temp(1)=pp2;temp(2)=pp4;
  ptsQuadT.push_back(temp);
  temp(1)=pp4;temp(2)=pp2;
  ptsQuadT.push_back(temp);

  // Calcul des poids de quadrature.
  pp1 = (155.-os)/2400; pp2 = (155.+os)/2400;
  pdsQuadT = Vecteur(7, pp2);
  pdsQuadT(1) = 9./80;
  pdsQuadT(2) = pp1;
  pdsQuadT(3) = pp1;
  pdsQuadT(4) = pp1;

  // Calcul des gradients des 6 fonctions de base (triangulation P2) sur les 7 points de quadrature
  gradsTRef = vector< vector< Vecteur> >(nbq);
  for (int k=0; k<nbq; k++) {
    double x = ptsQuadT[k](1);
    double y = ptsQuadT[k](2);
    temp(1) = 4*(x+y)-3;
    temp(2) = 4*(x+y)-3;
    gradsTRef[k].push_back(temp);
    temp(1) = 4*x-1;
    temp(2) = 0;
    gradsTRef[k].push_back(temp);
    temp(1) = 0;
    temp(2) = 4*y-1;
    gradsTRef[k].push_back(temp);
    temp(1) = -4*(-1+y+2*x);
    temp(2) = -4*x;
    gradsTRef[k].push_back(temp);
    temp(1) = 4*y;
    temp(2) = 4*x;
    gradsTRef[k].push_back(temp);
    temp(1) = -4*y;
    temp(2) = -4*(-1+x+2*y);
    gradsTRef[k].push_back(temp);
  }
}



// définition des valeurs de la fonction au bord
void Maillage::placeBord(Vecteur &U, double f(const Vecteur&)) {
  for (int i:points_bord_) {
    Vecteur X(nodes_[i]);
    U[i] = f(X);
  }
}

void Maillage::initFonction(Vecteur &U, double f(const Vecteur&)) {
  for (unsigned int i=0; i<nodes_.size(); i++) {
    Vecteur X(nodes_[i]);
    U[i] = f(X);
  }
}

// Calcul de surface avec des polynômes d'ordre 2, approx quadr 7pts
double Maillage::value(const Vecteur &U) {
  vector<Vecteur> ptsQuadT;
  Vecteur pdsQuadT;
  vector< vector< Vecteur> > gradsTRef;

  if (ordre_elts_geo_ == 1) {
    this->gaussLegendre(ptsQuadT, pdsQuadT, gradsTRef);
  } else {
    this->gaussLobatto(ptsQuadT, pdsQuadT, gradsTRef);
  }

  int nbq = ptsQuadT.size();
  int nbPts = gradsTRef[0].size();

  double Surf = 0;

  for (Triangle T:triangles_) {
    Matrice J(2,2);
    J(1,1) =  nodes_[T[1]][0] - nodes_[T[0]][0];
    J(2,1) =  nodes_[T[1]][1] - nodes_[T[0]][1];
    J(1,2) =  nodes_[T[2]][0] - nodes_[T[0]][0];
    J(2,2) =  nodes_[T[2]][1] - nodes_[T[0]][1];
    Matrice invtJ = inv(t(J));
    double Jl = 0;
    for (int k=0;k<nbq;k++) {

      Vecteur sum_grad(2);
      for (int i=0; i<nbPts; i++) {
        sum_grad += U[T[i]] * gradsTRef[k][i];
      }
      sum_grad = invtJ * sum_grad;

      double g2 = ps(sum_grad,sum_grad);
      Jl+=pdsQuadT[k]*(alpha_*sqrt(1+g2)+beta_*g2);
    }
    Surf += Jl*J.det();
  }
  return Surf;
}


// Calucl de surface avec des polynômes d'ordre 2, approx quadr 7pts
Vecteur Maillage::gradient(const Vecteur &U) {
  vector<Vecteur> ptsQuadT;
  Vecteur pdsQuadT;
  vector< vector< Vecteur> > gradsTRef;

  if (ordre_elts_geo_ == 1) {
    this->gaussLegendre(ptsQuadT, pdsQuadT, gradsTRef);
  } else {
    this->gaussLobatto(ptsQuadT, pdsQuadT, gradsTRef);
  }

  int nbq = ptsQuadT.size(); // nb de points de quadrature
  int nbPts = gradsTRef[0].size(); // nb de pts physiques pour 1 triangle

  Vecteur gradient(U.size());

  for (Triangle T:triangles_) {
      Matrice J(2,2);
      J(1,1) =  nodes_[T[1]][0] - nodes_[T[0]][0];
      J(2,1) =  nodes_[T[1]][1] - nodes_[T[0]][1];
      J(1,2) =  nodes_[T[2]][0] - nodes_[T[0]][0];
      J(2,2) =  nodes_[T[2]][1] - nodes_[T[0]][1];
      Matrice invtJ = inv(t(J));
      double mesure = this->mesure(T);

      for (int s=0;s<3*ordre_elts_geo_;s++) {
        if (!(this->estAuBord(T[s]))) {
          for (int k=0;k<nbq;k++) {

            Vecteur sum_grad(2);
            for (int i=0; i<nbPts; i++) {
              sum_grad += U[T[i]] * gradsTRef[k][i];
            }
            sum_grad = invtJ * sum_grad;

            double g2 = ps(sum_grad,sum_grad);
            double prodScal = ps(invtJ * gradsTRef[k][s], sum_grad);

            gradient[T[s]] += 2*mesure*pdsQuadT[k]*(alpha_*prodScal/sqrt(1+g2) + 2*beta_*prodScal);
          }
        }
      }
  }
  return gradient;
}


void Maillage::saveToFile(string fichier) {
  fstream fs;
  fs.open ((fichier + "_ptsX.m").c_str(), fstream::in | fstream::out | fstream::trunc);
  for (vector<double> X:nodes_) fs << X[0] << endl;
  fs.close();

  fs.open ((fichier + "_ptsY.m").c_str(), fstream::in | fstream::out | fstream::trunc);
  for (vector<double> X:nodes_) fs << X[1] << endl;
  fs.close();

  fs.open ((fichier + "_triangles.m").c_str(), fstream::in | fstream::out | fstream::trunc);
  if (ordre_elts_geo_ == 1) {
    for (Triangle T:triangles_) fs << T[0]+1 << " " << T[1]+1 << " " << T[2]+1 << endl;
    fs.close();
  } else {
    for (Triangle T:triangles_) {
      fs << T[0]+1 << " " << T[3]+1 << " " << T[5]+1 << endl;
      fs << T[3]+1 << " " << T[1]+1 << " " << T[4]+1 << endl;
      fs << T[5]+1 << " " << T[3]+1 << " " << T[4]+1 << endl;
      fs << T[5]+1 << " " << T[4]+1 << " " << T[2]+1 << endl;
    }
    fs.close();
  }
}


// instertion d'un elt dans une liste triee
void insereDansTrie(vector<int>& vect, const int x) {
  bool done = false;
  for (vector<int>::iterator it = vect.begin() ; it!=vect.end()&&!done; ++it) {
    if (*it == x) done = true;
    if (*it > x) {
      vect.emplace(it, x);
      done = true;
    }
  }
  if (!done) vect.push_back(x);
}

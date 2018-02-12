#ifndef MAILLAGE_H_INCLUDED
#define MAILLAGE_H_INCLUDED

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "formes.hpp"
#include "vecteur.hpp"
#include "matrice.hpp"
#include "optimiseur.hpp"

using namespace std;

class Maillage : public Probleme_abstract
{
private:
  vector<Point> points_;
  vector<Segment> segments_;
  vector<Triangle> triangles_;
  vector< vector<double> > nodes_;
  vector<int> points_bord_;
  int ordre_elts_geo_;
  double alpha_;
  double beta_;

public:
  Maillage(string, int ordre = 1, double a = 1, double b = 0.01);
  void summary();

  // Accesseurs
  Point points(int) const;
  Segment segments(int) const;
  Triangle triangles(int) const;
  vector<double> nodes(int) const;
  int points_bord(int) const;
  int points_size() const;
  int segments_size() const;
  int triangles_size() const;
  int nodes_size() const;
  int points_bord_size(int) const;

  // Test d'appartenance d'un noeud au bord
  bool estAuBord(int) const;

  // définition des valeurs de la fonction au bord
  void placeBord(Vecteur&, double(const Vecteur&));
  void initFonction(Vecteur&, double(const Vecteur&));

  // Methodes de quadrature
  void gaussLegendre(vector<Vecteur> &ptsQuadT, Vecteur &pdsQuadT, vector< vector< Vecteur> > &gradsTRef) const;
  void gaussLobatto(vector<Vecteur> &ptsQuadT, Vecteur &pdsQuadT, vector< vector< Vecteur> > &gradsTRef) const;


  // Fonctions sur triangle
  double mesure(const Forme&) const;
  vector<Vecteur> gradients(const Forme&) const;

  // Fonctionnelle et gradient
  double valueExpl(const Vecteur &U);
  Vecteur gradientExpl(const Vecteur &U);

  double value(const Vecteur &U);
  Vecteur gradient(const Vecteur &U);

  // Enregistrement du maillage (pts_x, pts_y, triangles)
  void saveToFile(string fichier = "Vecteur.m");
};

// Utilitaire d'insertion sur un vector<int> *trié*
void insereDansTrie(vector<int>&, const int);

#endif // MAILLAGE_H_INCLUDED

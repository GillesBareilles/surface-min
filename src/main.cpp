#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "fctOpt.hpp"
#include "maillage.hpp"
#include "vecteur.hpp"
#include "matrice.hpp"
#include "optimiseur.hpp"

using namespace std;



int main(int argc, char const *argv[]) {
  string monFichier;
  if (argc==1) {
    monFichier += "./geom_descr/geomTrou.geo";
  } else if (argc==2) {
    monFichier += argv[1];
  }
  Maillage M(monFichier);

  Vecteur U(M.nodes_size(), 0);
  Vecteur G(M.nodes_size(), 0);
  Vecteur GExpl(M.nodes_size(), 0);

  cout << "Fonctionnelle " << M.value(U) << endl;
  cout << "FonctionnelleExpl " << M.valueExpl(U) << endl;
  cout << "gradient " << M.gradient(U).norme() << endl;
  cout << "gradientExpl " << M.gradientExpl(U).norme() << endl;

  M.initFonction(U, norme_quad);
  G = M.gradient(U);
  GExpl = M.gradientExpl(U);


  U.saveToFile("out/vecteur_init.m");
  G.saveToFile("out/vecteur_init_g.m");
  GExpl.saveToFile("out/vecteur_init_gExpl.m");
  M.saveToFile("out/maillage");

  cout << "Fonctionnelle " << M.value(U) << endl;
  cout << "FonctionnelleExpl " << M.valueExpl(U) << endl;
  cout << "gradient " << norme(M.gradient(U)) << endl;
  cout << "gradientExpl " << norme(M.gradientExpl(U)) << endl;

  Optimiseur opt(M);

  Vecteur SolBFGS = opt.quasiNewtonBFGS(U, 0.01,500,true);
  SolBFGS.saveToFile("out/vecteur_fin.m");


  Vecteur Sol = opt.gradientConj(U, 0.01,500,true);
  Sol.saveToFile("out/vecteur_fin.m");

  return 0;
}

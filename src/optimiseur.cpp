#include "optimiseur.hpp"


Vecteur Optimiseur::gradientFixe(const Vecteur& x0, double eps, int itMax, bool verbose) {

  cout << "------------- Optimisation : gradient à pas fixe -------------" << endl;
  cout << "--------------------------------------------------------------"<<endl;

  int iteration = 0;
  int iterationMax = itMax;

  Vecteur x = x0;
  Vecteur g = this->prob_->gradient(x0);
  Vecteur dx(x0.size());

  cout << "Norme du gradient initial "<< norme(g) << "  ; eps : "<< eps << endl ;

  while ((norme(g) > eps) and (iteration < iterationMax))  {
    iteration +=1;
    g = this->prob_->gradient(x);

    //Choix de la direction de descente
    dx = -g/norme(g);

    //Descente
    double alpha = 1;
    x += alpha*dx;
    if ((iteration % 1 == 0) and verbose) {
      cout << "it " << iteration << " (|dx| = " << norme(alpha*dx) << ")" << "  alpha " << alpha << "  grad " << norme(g) << endl;
    }
  }

  cout << "\n----------------- Optimisation terminée ----------------------" << endl;
  cout << "Optimisation en "<<iteration<<" iterations (itMax="<<iterationMax<< ")." << endl;
  cout << " (|g|="<<norme(g)<<")." << endl << endl;

  return x;
}


Vecteur Optimiseur::gradientPasVar(const Vecteur& x0, double eps, int itMax, bool verbose) {
  cout << "----------- Optimisation : gradient à pas variable -----------" << endl;
  cout << "---------------- (recherche linéaire de Wolfe) ---------------"<<endl;

  int iteration=0;
  int iterationMax = itMax;
  double alpha=0;

  Vecteur x = x0;
  Vecteur g = this->prob_->gradient(x0);
  Vecteur d(x0.size());

  cout << "Norme du gradient initial "<< norme(g) << "  ; eps : "<< eps << endl;

  while ((norme(g) > eps) and (iteration < iterationMax))  {
    iteration +=1;
    g = this->prob_->gradient(x);

    //Calcul de la direction de descente
    d = -g;

    //Calcul de la longueur du pas de gradient
    alpha = rechercheLineaireWolfe(x, d, 1, false); //

    //Mise à jour des variables
    x += alpha*d;

    if ((iteration % 1 == 0) and verbose) {
      cout << "it " << iteration << " (|dx| = " << norme(alpha*d) << ")" << "  alpha " << alpha << "  grad " << norme(g) << "  f " << this->prob_->value(x) << endl;
    }
  }

  cout << "\n----------------- Optimisation terminée ----------------------" << endl;
  cout << "Optimisation en "<<iteration<<" iterations (itMax="<<iterationMax<< ")." << endl;
  cout << " f_final: " << this->prob_->value(x) << " (|dx| = " << norme(alpha*d) << ", |g|="<<norme(g)<<")." << endl;
  return x;
}

Vecteur Optimiseur::gradientConj(const Vecteur& x0, double eps, int itMax, bool verbose) {
  cout << "----------- Optimisation : gradient à pas conjugué -----------" << endl;
  cout << "---------------- (Ploak-Ribière - RL de Wolfe) ---------------"<<endl;

  int iteration=0;
  int iterationMax=itMax;
  double beta_k; double alpha=0;
  clock_t start_opt;
  start_opt = clock();

  Vecteur x = x0;
  Vecteur g = this->prob_->gradient(x0);
  Vecteur gp = g;
  Vecteur d(x0.size());

  cout << "Norme du gradient initial: "<< norme(g) << " eps : "<< eps << endl;

  while ((norme(g) > eps) and (iteration < iterationMax))  {
    iteration +=1;
    gp = g;
    g = this->prob_->gradient(x);

    //Calcul de la direction de descente
    if (iteration == 1) d = -g;
    else {
      beta_k = ps(g, g-gp)/ps(gp,gp);
      d = -g + beta_k*d;
    }

    //Calcul de la longueur du pas de gradient
    alpha = rechercheLineaireWolfe(x, d, 1, false);

    //Mise à jour des variables
    x += alpha*d;

    if ((iteration % 1 == 0) and verbose) {
      cout << "it " << iteration << " (|dx| = " << norme(alpha*d) << ")" << "  alpha " << alpha << "  grad " << norme(g) << "  f " << this->prob_->value(x) << endl;
    }
  }

  cout << "\n----------------- Optimisation terminée ----------------------" << endl;
  cout << "Optimisation en "<<iteration<<" iterations (itMax="<<iterationMax<< ") et " << ( std::clock() - start_opt ) / (double) CLOCKS_PER_SEC << "s  à  " << ( std::clock() - start_opt ) / (double) (CLOCKS_PER_SEC * iteration)   << " s/it." << endl;
  cout << " f_final: " << this->prob_->value(x) << " (|dx| = " << norme(alpha*d) << ", |g|="<<norme(g)<<")." << endl;
  return x;
}


Vecteur Optimiseur::quasiNewtonBFGS(const Vecteur& x0, double eps, int itMax, bool verbose) {
  cout << "----------- Optimisation : Quasi Newton (BFGS) -----------" << endl;
  cout << "------------ (et recherche linéaire de Wolfe) ------------"<<endl;

  int iteration=1;
  int iterationMax=itMax;
  double alpha = 0;
  clock_t start_opt;
  start_opt = clock();

  Vecteur x = x0;
  Vecteur x_p = x0;
  Vecteur grad = this->prob_->gradient(x0);
  Vecteur gradPrev = this->prob_->gradient(x0);

  Vecteur d(x0.size());
  Vecteur s(x0.size());  // x_k+1 - x_k
  Vecteur s_bar(x0.size());
  Vecteur y(x0.size());  // g_k+1 - g_k
  Vecteur y_bar(x0.size());

  Matrice W_k(x0.size(), x0.size());
  Matrice Id(x0.size(), x0.size());
  Id.identite();
  W_k.identite();

  cout << "Norme du gradient initial "<< norme(grad) << "  ; eps : "<< eps << endl;

  while ((norme(grad) > eps) and (iteration < iterationMax))  {

    //Calcul de la direction de descente
    d = -W_k * grad;

    if (ps(d,grad) >= 0) {
      cout << "La direction de descente ne convient pas. Reset H" << endl;
      W_k.identite();
      d = -grad;
    }

    //Calcul de la longueur du pas de gradient
    alpha = rechercheLineaireWolfe(x, d, 1, false);

    //Mise à jour des variables
    s = alpha*d;
    x += s;
    gradPrev = grad;
    grad = this->prob_->gradient(x);
    y = grad-gradPrev;

    //Mise a jour de W_k :
    double rho = 1.0 / ps(y,s);


    clock_t start; double duration;
    start = clock();

    // W_k = W_k - rho * (matL1(s, y)*W_k + matL1((W_k * y), s)) + rho * rho * (ps(y, W_k * y) + 1.0 / rho) * matL1(s,s); // \!/ produit matrice matrice, pas opti
    W_k = W_k - rho * (s * y.t() *W_k + W_k * y * s.t()) + rho * rho * (ps(y, W_k * y) + 1.0 / rho) * (s * s.t());
    if (verbose) duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;


    if (iteration > 1) {
    } else {
      W_k *= ps(y,s)/ps(y,y); // Mise à l'échelle à l'étape 1
    }

    iteration +=1;

    if ((iteration % 1 == 0) and verbose) {
      cout << "it " << iteration << " (|dx| = " << norme(alpha*d) << ")" << "  alpha " << alpha << "  grad " << norme(grad) << "  f " << this->prob_->value(x) << " temps de calcul : "<< duration <<  endl;
    }
  }

  cout << "\n----------------- Optimisation terminée ----------------------" << endl;
  cout << "Optimisation en "<<iteration<<" iterations (itMax="<<iterationMax<< ") et " << ( std::clock() - start_opt ) / (double) CLOCKS_PER_SEC << "s  à  " << ( std::clock() - start_opt ) / (double) (CLOCKS_PER_SEC * iteration)   << " s/it." << endl;
  cout << " f_final: " << this->prob_->value(x) << " (|dx| = " << norme(alpha*d) << ", |g|="<<norme(grad)<<")." << endl;
  return x;
}

/**
  Recherche linéaire de Wolfe : détermination du pas optimal pour le point et la direction de descente d spécifiée. Le pas initial alpha est donné.
  Algo de Fletcher-Lemarechal
*/
double Optimiseur::rechercheLineaireWolfe(const Vecteur &x, const Vecteur &d, double alpha0, bool verbose) {
  double omega1 = 0.1;
  double omega2 = 0.9;
  double dltx = 0.000001; //écart minimal entre les points donnes par deux pas successifs de l'algorithme.

  double infPos = 1.0e+15;
  double alphaMin = 0.; double alphaMax = infPos;

  // Initialisation de l'algorithme
  double alpha = alpha0;
  double alphaPrev = alpha0 + 2*dltx/norme(d);

  double fn; double f=this->prob_->value(x);
  Vecteur xn, gn; Vecteur g=this->prob_->gradient(x);
  bool wolfe1, wolfe2;

  if (verbose) cout << "\n---- debut de la recherche linéaire de Wolfe ----" << endl;

  do {
    xn = x + (alpha*d);
    fn = this->prob_->value(xn);
    gn = this->prob_->gradient(xn);
    wolfe1 = (fn <= f + omega1*alpha*(ps(g,d))); // décroissance ?
    wolfe2 = (ps(gn,d) >= omega2*(ps(g,d))); // pas pas trop petit ?

    if (verbose) cout << "   Wolfe1,wolfe2 " << wolfe1 << "," << wolfe2 << " alphaMin " << alphaMin << "   alphaMax " << alphaMax << endl;

    alphaPrev=alpha;

    if (!wolfe1) {
      alphaMax=alpha;
      alpha=0.5*(alphaMin+alphaMax);
    } else {
      if (!wolfe2) {
        alphaMin=alpha;
        if  (alphaMax == infPos) alpha=2*alphaMin;
        else alpha=0.5*(alphaMin+alphaMax);
      }
    }
  } while ( !(wolfe1 && wolfe2) && (abs(alpha-alphaPrev)*norme(d) > dltx));

  if (verbose) {
    cout << "------------- Fin de recherche de Wolfe ---------" << endl;
    cout << "alphaMin " << alphaMin << "   alphaMax " << alphaMax << endl;
  }
  return alpha;
}

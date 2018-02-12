%% Affichage Vecteur Solution

clear all
nom_maillage = '../geom_descr/geomCarre.msh'; % Attention au fichier géométrie

% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes] = lecture_msh(nom_maillage);


V = importdata('vecteur_grad_init.m');
% Affichage
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),V);
% shading interp
shading faceted %interp
% shading flat
colorbar;

% ajouter eventuellement un titre
title('toto');

%% Affichage
    
X = importdata('maillage_ptsX.m');
Y = importdata('maillage_ptsY.m');
Triangles = importdata('maillage_triangles.m');
W = importdata('vecteur_init.m');


% Affichage
trisurf(Triangles,X,Y,W);

% shading interp
shading faceted
% shading flat
colorbar;

% ajouter eventuellement un titre
title('toto');
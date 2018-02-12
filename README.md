# surface-min

When pulling a wire from soap water, a thin soapy surface creates. Its shape results from the equilibrium between superficial tension and gravity. This project aims at computing such a surface given a description of the wire(s).


From a mathematical standpoint, given an open set ![equation](https://latex.codecogs.com/gif.latex?%5COmega%20%5Csubset%20%5Cmathbb%7BR%7D%5E2) of boundary ![equation](https://latex.codecogs.com/gif.latex?%5CGamma), and a function ![equation](https://latex.codecogs.com/gif.latex?g%20%3A%20%5CGamma%20%5Crightarrow%20%5Cmathbb%7BR%7D), the exact formulation consists in finding a solution ![equation](https://latex.codecogs.com/gif.latex?u) s.t. :

<center>
![equation](https://latex.codecogs.com/gif.latex?%5Cmin_%7Bu%5Cin%20H%5E1%28%5COmega%29%2C%20u%7C_%5CGamma%20%3D%20g%7D%20%5Cint_%5COmega%20%281&plus;%7C%5Cnabla%20u%7C%5E2%29%5E%7B%5Cfrac%7B1%7D%7B2%7D%7D%20d%5COmega)
</center>

However, this problem is ill-formulated, so we consider, for

<center>
![equation](https://latex.codecogs.com/gif.latex?%5Cmin_%7Bu%5Cin%20H%5E1%28%5COmega%29%2C%20u%7C_%5CGamma%20%3D%20g%7D%20%5Calpha%20%5Cint_%5COmega%20%281&plus;%7C%5Cnabla%20u%7C%5E2%29%5E%7B%5Cfrac%7B1%7D%7B2%7D%7D%20d%5COmega%20&plus;%20%5Cbeta%20%5Cint_%5COmega%20%7C%5Cnabla%20u%20%7C%5E2%20d%5COmega)
</center>

Classes implémentées et fonctionnalités :
---

Calcul:
- classe Vecteur :
- classe Matrice :
- classe Optimiseur : gradient à pas variable (RL de Wolfe), gradient conjugué, _BFGS (en cours)_


Structure de données:
- classe Point, Segment, Triangle: héritent de Shape
- classe Maillage: contient tous les nodes, points segments et triangles de la surface considérée. Le constructeur se base sur un fichier msh.

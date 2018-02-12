# surface-min

When pulling a wire from soap water, a thin soapy surface creates. Its shape results from the equilibrium between superficial tension and gravity. This project aims at computing such a surface given a description of the wire(s).

From a mathematical standpoint, given an open set ![equation](https://latex.codecogs.com/gif.latex?%5COmega%20%5Csubset%20%5Cmathbb%7BR%7D%5E2) of boundary $\Gamma$, and a function $g : \Gamma \rightarrow \mathbb{R}$, the problem consists in finding a solution $u$ s.t. :


Classes implémentées et fonctionnalités :
---

Calcul:
- classe Vecteur :
- classe Matrice :
- classe Optimiseur : gradient à pas variable (RL de Wolfe), gradient conjugué, _BFGS (en cours)_


Structure de données:
- classe Point, Segment, Triangle: héritent de Shape
- classe Maillage: contient tous les nodes, points segments et triangles de la surface considérée. Le constructeur se base sur un fichier msh.

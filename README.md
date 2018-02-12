# surface-min

Classes implémentées et fonctionnalités :
---

Calcul:
- classe Vecteur :
- classe Matrice :
- classe Optimiseur : gradient à pas variable (RL de Wolfe), gradient conjugué, _BFGS (en cours)_


Structure de données:
- classe Point, Segment, Triangle: héritent de Shape
- classe Maillage: contient tous les nodes, points segments et triangles de la surface considérée. Le constructeur se base sur un fichier msh.

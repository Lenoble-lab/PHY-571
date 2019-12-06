# Implémentation de l'algorithme de Barn_Hut

L'algorithme se trouve dans le fichier array_node et array_node_2nd_order, selon l'ordre du schéma.

Dans notre implémentation, la seule classe est la classe Node qui correspond à un noeud de l'arbre. 
Si ce noeud est virtuel, elle contient la position du centre de masse et la masse notamment. 
Si est le noeud est une particule réelle, elle contient la masse, la force, et l'energie potentielle de la particule. 
Nous construisons l'arbre en initialisant la racine directement (nous n'effectuons que des descentes de l'arbre et pas d'opérations sur son architecture)

Enfin, les particules ont leurs informations de stockées dans différents tableaux numpy : 
    positions pour la position
    velocities pour la vitess
    masse pour la masse
    force/energie,...
A chaque particule est associé un identifiant unique qui correspond à la postition de la particule dans cette liste. 

# PHY-571

# Architecture
class Particule 
- position, masse, accel, velocity


class Node 
- particule, nb_children, box_size, id, [NW, NE, SW, SE]
- init
- create_tree (simul)
- calculate_centers, calculate_masses


class simulation 
- delta_t, nb_point, size, G, particules[]
- init
- init_gaussian (et plus selon affinites)
- calculate_accel("cree l'arbre et calcule les accel")
- update_velocity
- update_position

main : conditions initiales et fait avancer le pas de temps

class exploitation


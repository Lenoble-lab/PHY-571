# PHY-571

# Architecture
class Particule 
- position, masse, force, velocity


class Node 
- particule, nb_children, box_size, id, [NW, NE, SW, SE]
- init
- create_tree (simul)
- calculate_centers, calculate_masses


class simulation 
- delta_t, nb_point, size, G, particules[], root
- init
- init_gaussian (et plus selon affinites)
- set_tree ("met root a jour")
- calculate_force(calcule force, potentiel)
- calculate_accel
- update_velocity
- update_position

main : conditions initiales et fait avancer le pas de temps

class exploitation
???

import numpy as np

class Particule :

    def __init__ (self, pos, mass, force, velocity, id):
        self.position = pos     ##position vector of the particule (ex : r)
        self.mass = mass        ## mass of the particule
        self.force = force      ## force applied to the particule
        self.velocity = velocity    ## velocity of the particule
        self.id = id            ## number of the particule (-1 if the paticule doesnt really exist)
    

class Node :

    def __init__ (self, virtual_particule, box_size, box_center):
        self.virtual_particule = virtual_particule          ## particule in the node == real particule if leaf
        self.nb_children = 0          ## nb of children
        self.box_size = box_size                ## size of the box if leaf
        self.children = np.zeros(4)               ## array of the children (nodes)
        self.box_center = box_center            ## coordinates of the center of the box

    def calculate_mass_COM() : #update mass and center of mass, implicit recursive function with create_tree
        mass = 0
        position = np.zeros_like(self.virtual_particule.position)

        for i in range (self.children):
            mass += self.children[i].virtual_particule.mass

            position += self.children[i].virtual_particule.mass * self.children[i].virtual_particule.position

        self.virtual_particule.mass += mass
        self.virtual_particule.position = position / self.virtual_particule.mass


class Tree :

    def __init__(self, node) :
        self.current_node = node  ## node on which the pointeur is

    
    def create_tree (self, particules):
    


        if len(particules)==1 :     ##ie there is only one particule
            self.current_node = Node (particules[0], self.current_node.box_size/2 , np.array([0,0]))  #we don't care box_center
        
        else :
            self.generate_children(particules)          ##function which creates the tree
            current_node.calculate_mass_COM()               ## calculates the total mass and the center of mass

    
    def generate_children(self, particules):
        x, y = self.current_node.box_center             
        size = self.current_node.box_size
        
        particules_nw = []
        particules_ne = []
        particules_sw = []
        particules_se = []

        for part in particules :
            part_x, part_y = part.position

            if part_x > x :
                if part_y > y :
                    particules_ne.append(part)
                else :
                    partiucles_se.append(part)
            else :
                if part_y > y :
                    particules_nw.append(part)
                else :
                    particules_sw.append(part)

        
        ## create the nw child
        if len(particules_nw)>0 :
            nwnode = Node(Particule((x-size/2, y+size/2),0,0,0,-1), size/2, (x-size/2, y+size/2))
            nwnode.create_tree(particules_nw)
            self.current_node.children[0] = nwnode
            self.current_node.nb_children += 1

        ## creates the ne child
        if len(particules_ne)>0 :
            nenode = Node(Particule((x+size/2, y+size/2),0,0,0,-1), size/2, (x+size/2, y+size/2))
            nenode.create_tree(particules_ne)
            self.current_node.children[1] = nenode
            self.current_node.nb_children += 1


         ## creates the se child
        if len(particules_se)>0 :
            senode = Node(Particule((x+size/2, y-size/2),0,0,0,-1), size/2, (x+size/2, y-size/2))
            senode.create_tree(particules_se)
            self.current_node.children[1] = senode
            self.current_node.nb_children += 1

         ## creates the sw child
        if len(particules_sw)>0 :
            swnode = Node(Particule((x-size/2, y-size/2),0,0,0,-1), size/2, (x-size/2, y-size/2))
            swnode.create_tree(particules_sw)
            self.current_node.children[1] = swnode
            self.current_node.nb_children += 1






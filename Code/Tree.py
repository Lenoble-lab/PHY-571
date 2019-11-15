import numpy as np

class Particule :

    def __init__ (self, pos, mass, force, velocity, id):
        self.position = pos     ##position vector of the particule (ex : r)
        self.mass = mass        ## mass of the particule
        self.force = force      ## force applied to the particule
        self.velocity = velocity    ## velocity of the particule
        self.id = id            ## number of the particule (-1 if the paticule doesnt really exist)
        self.potential = 0        #For calculation of total energy
    
    def print_particule(self):
        print("particule")
        print(self.position, ' position')
        print(self.velocity, ' velocity')
        print(self.force, ' force')
        print(self.mass, " mass")
        print(" ")

class Node :

    def __init__ (self, box_size, box_center, particules): #créé l'arbre donné pour la liste "particules" de particules

        self.nb_children = 0          ## nb of children
        self.box_size = box_size                ## size of the box if leaf
        self.children = np.empty(4, dtype = Node)               ## array of the children (nodes)
        self.box_center = box_center            ## coordinates of the center of the box

        if len(particules)==1 :     ##if there is only one particule
            self.virtual_particule = particules[0]
            
        else :
            self.virtual_particule = Particule(self.box_center, 0, 0, np.array([0,0]), -1)
            self.generate_children(particules)          ##function which creates the tree
            self.calculate_mass_COM()               
    
    
    def calculate_mass_COM(self) : #update mass and center of mass, implicit recursive function with create_tree
        mass = 0
        position = np.zeros_like(self.virtual_particule.position)
        

        for i in range (len(self.children)):
            if self.children[i] != None:
                self.children[i].virtual_particule.print_particule()
                mass += self.children[i].virtual_particule.mass
                position += self.children[i].virtual_particule.mass * self.children[i].virtual_particule.position

        self.virtual_particule.mass += mass
        self.virtual_particule.position = position / self.virtual_particule.mass

    
    def generate_children(self, particules):
        
        x, y = self.box_center             
        size = self.box_size
        
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
                    particules_se.append(part)
            else :
                if part_y > y :
                    particules_nw.append(part)
                else :
                    particules_sw.append(part)

        if len(particules_nw)>0 :
            self.children[0] = Node(self.box_size/2, np.array([x-size/2, y+size/2]), particules_nw)
            self.nb_children += 1

        ## creates the ne child
        if len(particules_ne)>0 :
            self.children[1] = Node(self.box_size/2, np.array([x+size/2, y+size/2]), particules_ne)
            self.nb_children += 1


         ## creates the se child
        if len(particules_se)>0 :
            self.children[2] = Node(self.box_size/2, np.array([x+size/2, y-size/2]), particules_se)
            self.nb_children += 1

         ## creates the sw child
        if len(particules_sw)>0 :
            self.children[3] = Node(self.box_size/2, np.array([x-size/2, y-size/2]), particules_sw)
            self.nb_children += 1






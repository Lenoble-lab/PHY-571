from Tree import *

class Simulation :
    def __init__(self, delta_t, size, center, G=1, theta = 1):
        self.delta_t = delta_t      #pas temporel
        self.nb_point = 0       #number of points of the simulation
        self.size = size            #size of the simulation
        self.center = center        #center of the grid
        self.G = G                  #gravitation constante
        self.particules = []        #particules, needs to be update
        self.root = None        #root of the tree
        self.theta = theta          #parameter of the simulation
    
    def step(self):
        self.nb_point = len(self.particules)
        self.set_tree()                 #create the tree with the current state of the particules
        print("set_tree finish")
        self.calculate_force()          #calcul the interaction on each particule
        self.update_velocity()          #update the velocity of each particule
        self.update_position()          #update the position of each particule


    def set_tree(self):
        self.root = Node(self.size, self.center, self.particules)

    def calculate_force(self):
        for i in range (self.nb_point):
            self.particules[i].force = np.array([0.,0.]) # re-set the force for the next step
            print("self.root")
            self.root.virtual_particule.print_particule()
            self.calculate_force_target(self.particules[i], self.root)
            
    
    def calculate_force_target(self, target_particule, temp_node):
        
        particule2 = temp_node.virtual_particule  #particule (or virtual particule) in question

        vect_r = particule2.position - target_particule.position
        r = np.sqrt(np.sum(vect_r**2))

        d = temp_node.box_size
        
        if target_particule.id == 2 : print(r, ' r')

        if r>0:  #if the node has only one particule
            
            if target_particule.id ==1 : 
                print(d/r, " d/r")
                print(target_particule.force, 'target part force')

            if (temp_node.nb_children == 0) or (d/r <=  self.theta): 
                
                #if the node is far enough then we approximate
                #or termination condition

                target_particule.force += self.G * particule2.mass * target_particule.mass / r**3 * vect_r
                
            
            else :  
                #not favorable case, split the node and repeat
                print("not favorable case")
                for i in range(4):
                    if temp_node.children[i] != None :
                        """
                        print("target particule : ")
                        target_particule.print_particule()
                        print("noeud enfant")
                        temp_node.children[i].virtual_particule.print_particule()
                        print(self.calculate_force_target(target_particule, temp_node.children[i]), " force to add")
                        """
                        self.calculate_force_target(target_particule, temp_node.children[i])
                        
    

    def update_velocity (self) :
        for i in range(self.nb_point):
            self.particules[i].velocity += (self.particules[i].force/self.particules[i].mass )*self.delta_t   #1st order of integration


    def update_position (self):
        for i in range(self.nb_point):
            self.particules[i].position += self.particules[i].velocity*self.delta_t   # 1st order of integration





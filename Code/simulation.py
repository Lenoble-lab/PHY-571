from Tree import *

class Simulation :
    def __init__(self, delta_t, nb_point, size, center, G=1, theta = 1):
        self.delta_t = delta_t      #pas temporel
        self.nb_point = nb_point    #number of points of the simulation
        self.size = size            #size of the simulation
        self.center = center        #center of the grid
        self.G = G                  #gravitation constante
        self.particules = []        #particules, needs to be update
        self.root                   #root of the tree
        self.theta = theta          #parameter of the simulation
    
    def step(self):
        self.set_tree()                 #create the tree with the current state of the particules
        self.calculate_force()          #calcul the interaction on each particule
        self.calculate_acceleration()   #calculate acceleration
        self.update_velocity()          #update the velocity of each particule
        self.update_position()          #update the position of each particule


    def set_tree(self):
        root_particule = Particule(self.center, 0, 0, 0, -1)
        root_node = Node(root_particule, self.size, self.center)

        tree = Tree(root_node)
        tree.create_tree(self.particules)   #create the tree
        self.root = tree.current_node       #extract the top node for calculation

    def calculate_force(self):
        for i in range (self.nb_point):
            self.particules[i].force = 0 # re-set the force for the next step
            self.calculate_force_target(self.particules[i])
    
    def calculate_force_target(self, target_particule, temp_node = self.root):
        
        particule2 = temp_node.particule  #particule (or virtual particule) in question

        vect_r = particule2.position - target_particule.position
        r = np.sqrt(np.sum(vect_r**2))

        d = temp_node.box_size

        if (d/r < self.theta) or (temp_node.nb_children == 0):  
            
            #if the node is far enough then we approximate
            #or termination condition

            target_particule.force += self.G * particule2.mass * target_particule.mass / r**2
        
        else :  
            #not favorable case

            for i in range(4):
                if temp_node.children[i] != 0 :
                    target_particule.force += calculate_force_target(target_particule, temp_node.children[i])
    


    def calculate_velocity (self) :
        for i in range(self.nb_point):
            self.particules[i].velocity += (-self.particules[i].force/self.particules[i].mass )*self.delta_t   #1st order of integration


    def update_position (self):
        for i in range(self.nb_point):
            self.particules[i].position += self.particules[i].velocity*self.delta_t   # 1st order of integration





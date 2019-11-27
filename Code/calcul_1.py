import numpy as np
import numpy.random as rnd
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numba import jit


def init_syst_soleil(N_part):
    M_soleil = 10**3
    R_max = 50.  
    P_soleil = np.array([0.,0.])
    V_soleil = np.array([0.,0.])



    positions = np.zeros((N_part, 2))
    masses = np.ones(N_part) 
    velocities = np.zeros((N_part, 2))

    positions[0] = P_soleil
    masses[0] = M_soleil
    velocities[0] = V_soleil
    
    
    i = 1
    while i < N_part : 

        theta = rnd.random() * 2 * np.pi
        r = rnd.random() * R_max
        #r = np.abs(rnd.normal(0,10))
        if r!=0 : 
            positions[i] = np.array([r*np.cos(theta), r*np.sin(theta)])
            M_tot = M_soleil + (r/R_max)**2 * N_part
            velocities[i] = np.sqrt(M_tot/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])
            i += 1
  
    return positions, masses, velocities




def calculate_force(temp_node, target_node, eps = 5, thetamax=0.7, G=1.0):
    """
    Calculate the interaction/energy between target_node and temp_node.
    Then : update the values in target_node.
    To compute everty interaction for a current node, initialise with temp_node = summit of the tree
    """

    vect_r = temp_node.COM - target_node.COM    # vector between nodes' centres of mass, for force calculation
    r = np.sqrt(np.sum(vect_r**2))   # distance between them
    
    if r>0:
        #checking that the calculate force is valid (no division by 0)
        # if the node only has one particle or theta is small enough,
        #  add the field contribution to value stored in node.force
        if (len(temp_node.children)==0) or (temp_node.size/r < thetamax):
            target_node.force += G * temp_node.mass * vect_r/r /(r**2 + eps**2)
            target_node.energy += -G * temp_node.mass * target_node.mass /(r + eps)
        else:
            # otherwise split up the node and repeat
            #not favorable case, split the node and repeat

            for children in temp_node.children: 
                calculate_force(children, target_node, eps, thetamax, G)



class Node:
    def __init__(self, center, size, masses, positions, ids, leaves=[]):
    
        """
        take as parameter the list (np.array) of masses, positions and ids of the particules
        center, size for the box
        """

        self.center = center                    # center of the node's box
        self.size = size                        # size of the box
        self.children = []                      # childrens, empty for now
 
        N_points = len(positions)
 
        if N_points == 1:
            # if there is one point, then the node is a real particule
            leaves.append(self)
            self.COM = positions[0]
            self.mass = masses[0]
            self.id = ids[0]
            self.force = np.zeros(2)        # at each point, we initailise the gravitational field
            self.energy = 0                  #and the enrgy potential
            self.der_force = 0          #for second order integration shema
        else:
            self.GenerateChildren(positions, masses, ids, leaves)     # if we have at least 2 points in the node, we generate the children
                                                             
            # updating the COM and mass of each node, we will calculate recursively on the node's children
            com_total = np.zeros(2) # running total for mass moments to get COM
            m_total = 0.            # running total for masses
            for c in self.children:
                m, com = c.mass, c.COM
                m_total += m
                com_total += com * m   # add the moments of each child
            self.mass = m_total
            self.COM = com_total / self.mass  
 
    def GenerateChildren(self, positions, masses, ids, leaves):
        """Generates the node's children"""
        tree_pos = (positions > self.center)  #does all comparisons needed to determine points' octants
        
        for i in range(2): #looping over the 8 octants
            for j in range(2):

                in_tree = np.all(tree_pos == np.bool_([i,j]), axis=1)
                if not np.any(in_tree): continue           # if no particles, don't make a node
                dx = 0.5*self.size*(np.array([i,j])-0.5)   # offset between parent and child box centers
                self.children.append(Node(self.center+dx,
                                                self.size/2,
                                                masses[in_tree],
                                                positions[in_tree],
                                                ids[in_tree],
                                                leaves))

@jit
def GravAccel(positions, masses, sec_order = False, thetamax=0.7, G=1.):
    center = (np.max(positions,axis=0)+np.min(positions,axis=0))/2       #center of bounding box is the mean of the max and min particule
    topsize = np.max(np.max(positions,axis=0)-np.min(positions,axis=0))  #size of bounding box
    leaves = []  # want to keep track of leaf nodes
    topnode = Node(center, topsize, masses, positions, np.arange(len(masses)), leaves) #build the tree
    
    """
    if sec_order : der_force = np.empty_like(positions)
    """
    force = np.empty_like(positions)
    energy = np.empty_like(positions)
    

    for i,leaf in enumerate(leaves):
        calculate_force(topnode, leaf, thetamax, G)  # update energy and force of every particules
        force[leaf.id] = leaf.force  # store force and accéleration in order to update the velocity and the position later
        energy[leaf.id] = leaf.energy
        if sec_order : der_force[leaf.id] = leaf.der_force
    
    if sec_order : return force, energy, der_force
    return force, energy
    


def step_1st_order(positions, masses, velocities, delta_t):
    """
    put it together : update the tree, calculate the force/energy and finaly update position/velocity with a first order shema
    """
    force, energy = GravAccel(positions, masses)
    N_points = len(positions)

    velocities += delta_t * force
    positions += delta_t * velocities
    return positions, velocities, np.sum(energy)

def step_leap_frog(positions, masses, velocities, force_i, delta_t):
    """
    put it together : update the tree, calculate the force/energy and finaly update position/velocity with a leapfrog shema
    """
    
    force, energy = GravAccel(positions, masses)
    velocities += 0.5*(force_i + force) * delta_t 
    positions += velocities * delta_t + 0.5*force_i * delta_t**2
    return positions, velocities, force, np.sum(energy)






sys.setrecursionlimit(10**5) 

t = time.clock()



positions, masses, velocities = init_syst_soleil(20)
N_cycle = 1000
N_part = len(positions)

pos = np.zeros((N_cycle+1, N_part, 2))
energy_pot = np.zeros(N_cycle)
energy_cin = np.zeros(N_cycle)
cintetic_momentum = np.zeros(N_cycle)

pos[0] = positions
force, energy = GravAccel(positions, masses)

for i in range (N_cycle):
    if i % 100 == 0:
        print (i, ' t ')
    
    positions, velocities, force, energy_pot_t = step_leap_frog(positions, masses, velocities, force, 0.005)
    # postitions, velocities, energy_pot_t = step_1st_order(positions, masses, velocities, 0.005)
    # positions, velocities, force, energy_pot_t = step_2nd_order(positions, masses, velocities, force, 0.005)

    pos [i+1] = positions
    energy_pot[i] = energy_pot_t/2
    energy_cin[i] = 0.5 * np.sum(masses * [np.sum(velocities[i,:]**2) for i in range (len(velocities))])
    cintetic_momentum[i] =  np.sum(masses * [positions[i][0] * velocities[i][1] - positions[i][1] * velocities[i][0] for i in range (len(positions))])
    


np.save("../analyse/data", np.array([pos, energy_pot, energy_cin, cintetic_momentum]))

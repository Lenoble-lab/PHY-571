import numpy as np
import numpy.random as rnd
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numba import jit
from Init import *

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
    N_points = len(positions)
    velocities += 0.5*(force_i + force) * delta_t 
    positions += velocities * delta_t + 0.5*force_i * delta_t**2
    return positions, velocities, force, np.sum(energy)


sys.setrecursionlimit(10**5) 

t = time.clock()



fig = plt.figure()

positions, masses, velocities = init_terr_soleil()
N_cycle = 1500
N_part = len(positions)

pos = np.zeros((N_cycle+1, N_part, 2))
energy_pot = np.zeros(N_cycle)
energy_cin = np.zeros(N_cycle)
cintetic_momentum = np.zeros(N_cycle)

pos[0] = positions

fig = plt.figure()
frames_size = 400
ax1 = fig.add_subplot(1,2,1)
ax1.set_xlim(-frames_size, frames_size)
ax1.set_ylim(-frames_size, frames_size)
line, = ax1.plot([], [], 'o', markersize=2)

ax2 = fig.add_subplot(1,2,2)

line_ene_pot, = ax2.plot([0], [0], 'r', label = 'energy potential', linewidth = 2)
line_ene_cin, = ax2.plot([0], [0], 'b', label = 'energy cinetic', lw = 2)
line_ene_tot, = ax2.plot([0], [0], 'g', label = 'energy total', lw = 2)
line_mom, = ax2.plot([0], [0], 'black', label = 'cinetic momentum', lw = 2)

ax2.set_xlim(0, 200)

energy_cin_0 = 0.5 * np.sum(masses * [np.sum(velocities[i,:]**2) for i in range (len(velocities))])

ax2.set_ylim(5 * -energy_cin_0, energy_cin_0 *1.5)
ax2.legend()

nt = 200

"""
t = time.clock()
for i in range (40):
    positions, velocities = step(positions, masses, velocities, 0.005)

print(time.clock()-t)
"""
force, energy = GravAccel(positions, masses)

def make_frame(i):
    if i % 10 == 0:
        print (i, ' t ')
    global positions
    global velocities
    global force
    
    positions, velocities, force, energy_pot_t = step(positions, masses, velocities, force, 0.005)
    # postitions, velocities, energy_pot_t = step_1st_order(positions, masses, velocities, 0.005)

    pos [i+1] = positions
    energy_pot[i] = energy_pot_t/2
    energy_cin[i] = 0.5 * np.sum(masses * [np.sum(velocities[i,:]**2) for i in range (len(velocities))])
    cintetic_momentum[i] =  0.1 * np.sum(masses * [positions[i][0] * velocities[i][1] - positions[i][1] * velocities[i][0] for i in range (len(positions))])
    
    line.set_data(pos[i,:,0], pos[i,:,1])
    if i<nt and i>0: 
        line_ene_cin.set_data(range(0,i), energy_cin[0:i])
        line_ene_pot.set_data(range(0,i), energy_pot[0:i])
        line_ene_tot.set_data(range(0,i), energy_pot[0:i] + energy_cin[0:i])
        line_mom.set_data(range(0,i), cintetic_momentum[0:i])

    elif i>=nt : 
        line_ene_cin.set_data(range(0,nt), energy_cin[-nt +i:i])
        line_ene_pot.set_data(range(0,nt), energy_pot[-nt + i:i])
        line_ene_tot.set_data(range(0,nt), energy_pot[-nt+i:i] + energy_cin[-nt+i:i])
        line_mom.set_data(range(0,nt), cintetic_momentum[-nt + i:i])


    return (line,  line_ene_tot, line_ene_cin, line_ene_pot, line_mom)

#ani = animation.FuncAnimation(fig, make_frame, frames=N_cycle, interval=30, repeat = False)
#ani.save("1000_part_0.005_deltat_0.005_gaussian.mp4")
print(time.clock() - t, 'temps')
#plt.show()

def calcul(N_cycle):
    global positions
    global velocities
    global force
    for i in range (N_cycle):
        if i % 100 == 0:
            print (i, ' t ')
        
        positions, velocities, force, energy_pot_t = step_leap_frog(positions, masses, velocities, force, 0.007)
        # positions, velocities, energy_pot_t = step_1st_order(positions, masses, velocities, 0.007)

        pos [i+1] = positions
        energy_pot[i] = energy_pot_t/2
        energy_cin[i] = 0.5 * np.sum(masses * [np.sum(velocities[i,:]**2) for i in range (len(velocities))])
        cintetic_momentum[i] =  - np.sum(masses * [positions[i][0] * velocities[i][1] - positions[i][1] * velocities[i][0] for i in range (len(positions))])
t = time.clock()
calcul(N_cycle)
print(time.clock()-t)

name = "trace_terre_soleil_0.007_leapfrog"

plt.figure()
plt.title("Trajectoire du système")
for i in range (N_part):
    plt.plot(pos[:,i,0], pos[:,i,1], 'o', markersize = 1)
plt.axis('equal')
plt.savefig("../rapport/" + name + "_trajectoire.jpg")

plt.figure()
plt.title("energie et moment cinétique")
plt.plot(range(0,N_cycle), energy_pot,  'r', label = 'energy potential', linewidth = 2)
plt.plot(range(0,N_cycle), energy_cin, 'b', label = 'energy cinetic', lw = 2)
plt.plot(range(0, N_cycle), energy_cin + energy_pot , 'g', label = 'energy total', lw = 2)
plt.plot(range(0,N_cycle), cintetic_momentum, 'black', label = 'cinetic momentum', lw = 2)
plt.legend()
#plt.axis('equal')
plt.savefig("../rapport/" + name + "_energy.jpg")
plt.show()


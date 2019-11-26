import numpy as np
import numpy.random as rnd
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numba import jit


def calculate_force(temp_node, target_node, thetamax=0.7, G=1.0):
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
            target_node.force += G * temp_node.mass * vect_r/r /(r**2 + 5)
        else:
            # otherwise split up the node and repeat
            #not favorable case, split the node and repeat

            for c in temp_node.children: calculate_force(c, target_node, thetamax, G)



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
def GravAccel(positions, masses, thetamax=0.7, G=1.):
    center = (np.max(positions,axis=0)+np.min(positions,axis=0))/2       #center of bounding box is the mean of the max and min particule
    topsize = np.max(np.max(positions,axis=0)-np.min(positions,axis=0))  #size of bounding box
    leaves = []  # want to keep track of leaf nodes
    topnode = Node(center, topsize, masses, positions, np.arange(len(masses)), leaves) #build the tree
 
    accel = np.empty_like(positions)
    for i,leaf in enumerate(leaves):
        calculate_force(topnode, leaf, thetamax, G)  # do field summation
        accel[leaf.id] = leaf.force  # get the stored acceleration
 
    return accel

def update_position(positions, accel, velocities, delta_t):
    N_points = len(positions)
    for i in range (N_points):
        velocities[i] += delta_t * accel[i]
        positions[i] += delta_t * velocities[i]
    return positions, velocities

def step(positions, masses, velocities, delta_t):
    accel = GravAccel(positions, masses)
    update_position(positions, accel, velocities, delta_t)
    return positions, velocities

def init_terr_soleil():
    M_soleil = 10**4
    M_terre = 1
    R = 10.
    
    P_soleil = np.array([0.,0.])
    P_terre = np.array([0., R])

    V_soleil = np.array([0.,0.])
    V_terre =np.array([np.sqrt(M_soleil/R), 0.])
    
    return np.array([P_soleil, P_terre]), np.array([M_soleil, M_terre]), np.array([V_soleil, V_terre])

def init_syst_soleil(N_part):
    M_soleil = 10**3
    R_max = 50.  
    P_soleil = np.array([0.,0.])
    V_soleil = np.array([0.,0.])



    positions = np.zeros((N_part, 2))
    masses = np.zeros(N_part)
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
            masses[i] = 1
            velocities[i] = np.sqrt(M_soleil/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])
            i += 1
  
    return positions, masses, velocities



sys.setrecursionlimit(10**5) 

t = time.clock()



fig = plt.figure()

positions, masses, velocities = init_syst_soleil(1000)
N_cycle = 1000
N_part = len(positions)

pos = np.zeros((N_cycle+1, N_part, 2))

pos[0] = positions

# pos_sol = [positions[0].tolist()]
# pos_terre = [positions[1].tolist()]
ax1 = plt.axes()
frames_size = 100
# ax1 = fig.subplot(1,2,1)
ax1.set_xlim(-frames_size, frames_size)
ax1.set_ylim(-frames_size, frames_size)
line, = ax1.plot([], [], 'o', markersize=3)

"""
t = time.clock()
for i in range (40):
    positions, velocities = step(positions, masses, velocities, 0.005)

print(time.clock()-t)
"""
def make_frame(i):
    if i % 10 == 0:
        print (i, ' t ')
    global positions
    global velocities
    positions, velocities = step(positions, masses, velocities, 0.005)
    
    pos [i+1] = positions
    # pos_sol += [positions[0].tolist()]
    # pos_terre += [positions[1].tolist()]
    line.set_data(pos[i,:,0], pos[i,:,1])
    return line,

ani = animation.FuncAnimation(fig, make_frame, frames=N_cycle, interval=30, repeat = False)
ani.save("1000_part_0.005_deltat_.mp4")
print(time.clock() - t, 'temps')
plt.show()

plt.figure()
for i in range (N_part):
    plt.plot(pos[:,i,0], pos[:,i,1], 'o')
plt.axis('equal')
plt.show()

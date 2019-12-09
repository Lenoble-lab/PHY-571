
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy.random as rnd
import sys
import time
import scipy.special
import scipy.integrate
import matplotlib.pyplot as plt


def init_carre_random(N_part, R):   #random initialisation of N_part particules in a square
    positions = np.zeros((N_part, 2))           #array for the positions
    masses = np.ones(N_part)                    #array for the masses
    velocities = np.zeros((N_part, 2))          #array for the velocities

    for i in range (N_part):
        positions[i] = np.array([R/2 * rnd.random() - R/2, R/2 * rnd.random() - R/2])

    return positions, masses, velocities

def init_syst_2_corps():            #initialisation of a two-particules system
    M_1 = 100                 #mass of the 1st particule
    M_2 = 100                  #mass of the 2nd particule
    R = 5.                    #initial distance between the particules
    


    #T = np.sqrt(a**3 * 4*np.pi**2 / (M_soleil + M_terre))

    P_1 = np.array([M_2/M_1 * R/(1 + M_2/M_1),0.])  #position of the 1st particule
    P_2 = np.array([M_1/M_2 * R/(1 + M_1/M_2), R])  #position of the 2nd particle

    # T_1 = 

    V_soleil = np.array([2*np.pi*R/T,0.])           #velocity of the 1st particule (Sun)
    V_terre =np.array([-2*np.pi*R/T, 0.])           #velocity of the 2nd particule (Earth)
 

    return np.array([P_soleil, P_terre]), np.array([M_soleil, M_terre]), np.array([V_soleil, V_terre])
def init_terr_soleil():                 #initialisation of an Earth-Sun system
    M_soleil = 10**4  #mass of the Sun
    M_terre = 1         #mass of the Earth
    R = 5.             #radius between them
    
    P_soleil = np.array([0.,0.])       #initial positions
    P_terre = np.array([0., R])

    V_soleil = np.array([0.,0.])             #initial velocities
    V_terre =np.array([np.sqrt(M_soleil/R), 0.])
    
    return np.array([P_soleil, P_terre]), np.array([M_soleil, M_terre]), np.array([V_soleil, V_terre])

def init_syst_soleil(N_part, R_max = 50.):         #initialisation of a solar system with N_part particules
    M_soleil = 10**6                    #mass of the Sun
    P_soleil = np.array([0.,0.])        #position of the Sun
    V_soleil = np.array([0.,0.])        #velocity of the Sun



    positions = np.zeros((N_part, 2))       #positions
    masses = np.ones(N_part)                #and masses
    velocities = np.zeros((N_part, 2))      #and velocities of the particules

    positions[0] = P_soleil         #initialisation 
    masses[0] = M_soleil
    velocities[0] = V_soleil
    
    
    i = 1
    while i < N_part :      #creation of N-part-1 random particules

        theta = rnd.random() * 2 * np.pi
        r = rnd.random() * R_max
        #r = np.abs(rnd.normal(0,10))
        if r > R_max/30 + 7 : 
            positions[i] = np.array([r*np.cos(theta), r*np.sin(theta)])
            M_tot = M_soleil + (r/R_max)**2 * N_part
            velocities[i] = np.sqrt(M_tot/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])
            i += 1
  
    return positions, masses, velocities

def init_galaxy(N_part, R_max = 100):
    """
    initialisation of a galaxy with the same proportion as the milky way
    we use a simplified density, with an exponenitial decrease
    """

    M_BH = 10**3              #mass of the black hole
    M_tot = 10**6             #total mass
    M_disk = M_tot - M_BH     #mass of the disk
    R_max = 100.              #maximal radius
    h = 1
    gamma = 10 #core radius
    r_c = 10**3 #cut of radius
    P_BH = np.array([0.,0.])    #position of the black hole
    V_BH = np.array([0.,0.])    #velocity of the black hole



    positions = np.zeros((N_part, 2))
    masses = np.ones(N_part) * M_disk/N_part
    velocities = np.zeros((N_part, 2))

    positions[0] = P_BH
    masses[0] = M_BH
    velocities[0] = V_BH

    plt.figure()
    x = np.linspace(0,10,1000)
    plt.plot(x, gamma * np.exp(-(x/r_c)**2)/(x**2 + gamma))

    plt.show()
    
    
    i = 1
    while i < N_part :      #creation of N_part-1 random particules

        theta = rnd.random() * 2 * np.pi
        r = -h * np.log(R_max * rnd.random()/M_tot)             #distance to the center of the galaxy

        if r > R_max/10 : 
            positions[i] = np.array([r*np.cos(theta), r*np.sin(theta)])

            M_int = M_BH + 2*np.pi * M_disk * (-1 + (1-r)*np.exp(-r))  #density of mass at the center
            velocities[i] = np.sqrt(M_int/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])      #velocity to have a circular trajectoire
            i += 1
  
    return positions, masses, velocities


def init_collision_galaxies(N_part, R_max = 300):       #initialisation of two galaxies of N_part particules together

    pos_1, masses_1, vel_1 = init_syst_soleil( 8* N_part//9, R_max)         #initialisation of the 1st galaxy

    pos_2, masses_2, vel_2 = init_syst_soleil(N_part//9, R_max/3)           #initialisation of the 2nd galaxy

    pos_2 = pos_2 + np.ones_like(pos_2) * R_max 

    masses_1[0] = 10**6
    masses_2[0] = masses_1[0] 


    for i in range (len(vel_2)):
        vel_2[i] = vel_2[i] + np.sqrt(np.sum(masses_1)/np.sqrt(R_max**2 * 2)) * np.array([-1., 1.])/2
    
    return np.concatenate((pos_1, pos_2), axis = 0), np.concatenate((masses_1, masses_2), axis = 0), np.concatenate((vel_1, vel_2), axis = 0)


# positions, masses, velocities = init_carre_random(4000, 400)

"""
plt.figure()
plt.plot(positions[:,0], positions[:,1], 'o', markersize = 1)
plt.show()
"""




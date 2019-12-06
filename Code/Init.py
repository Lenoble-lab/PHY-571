import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy.random as rnd
import sys
import time
import scipy.special
import scipy.integrate
import matplotlib.pyplot as plt

"""
contains every different initialisation that we have use
"""


def init_carre_random(N_part, R):
    """
    dispose de façon aléatoire des points dans un carré
    """
    positions = np.zeros((N_part, 2))
    masses = np.ones(N_part) 
    velocities = np.zeros((N_part, 2))

    for i in range (N_part):
        positions[i] = np.array([R/2 * rnd.random() - R/2, R/2 * rnd.random() - R/2])

    return positions, masses, velocities

def init_syst_2_corps():
    """
    initiate the 2 corps system with trajectories
    """

    M_1 = 100
    M_2 = 100
    R = 5.
    


    #T = np.sqrt(a**3 * 4*np.pi**2 / (M_soleil + M_terre))

    P_1 = np.array([M_2/M_1 * R/(1 + M_2/M_1),0.])
    P_2 = np.array([M_1/M_2 * R/(1 + M_1/M_2), R])

    # T_1 = 

    V_soleil = np.array([2*np.pi*R/T,0.])
    V_terre =np.array([-2*np.pi*R/T, 0.])
 

    return np.array([P_soleil, P_terre]), np.array([M_soleil, M_terre]), np.array([V_soleil, V_terre])

def init_terr_soleil():
    """
    initiate the 2 corps system with circular trajectories and one corps heavier than the second
    """
    M_soleil = 10**4
    M_terre = 1
    R = 5.
    
    P_soleil = np.array([0.,0.])
    P_terre = np.array([0., R])

    V_soleil = np.array([0.,0.])
    V_terre =np.array([np.sqrt(M_soleil/R), 0.])
    
    return np.array([P_soleil, P_terre]), np.array([M_soleil, M_terre]), np.array([V_soleil, V_terre])

def init_syst_soleil(N_part, R_max = 50.):
    """
    distribution in 1/R around a massive corps with a speed to ensure stability
    we use this for the simulation of galaxies
    """

    M_soleil = 10**6
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

    M_BH = 10**3
    M_tot = 10**6
    M_disk = M_tot - M_BH
    R_max = 100.
    h = 1
    gamma = 10 #core radius
    r_c = 10**3 #cut of radius
    P_BH = np.array([0.,0.])
    V_BH = np.array([0.,0.])



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
    while i < N_part : 

        theta = rnd.random() * 2 * np.pi
        r = -h * np.log(R_max * rnd.random()/M_tot)             #distance to the center of the galaxy

        if r > R_max/10 : 
            positions[i] = np.array([r*np.cos(theta), r*np.sin(theta)])

            M_int = M_BH + 2*np.pi * M_disk * (-1 + (1-r)*np.exp(-r))  #density of mass at the center
            velocities[i] = np.sqrt(M_int/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])      #velocity to have a circular trajectoire
            i += 1
  
    return positions, masses, velocities


def init_collision_galaxies(N_part, R_max = 300):
    """
    initialisation for a collision between galaxies
    the smaler galaxy contains 1/3 of the total points and is in rotation around the other galaxie
    """
    pos_1, masses_1, vel_1 = init_syst_soleil( 8* N_part//9, R_max)

    pos_2, masses_2, vel_2 = init_syst_soleil(N_part//9, R_max/3)

    pos_2 = pos_2 + np.ones_like(pos_2) * R_max 

    masses_1[0] = 10**6
    masses_2[0] = masses_1[0] 


    for i in range (len(vel_2)):
        vel_2[i] = vel_2[i] + np.sqrt(np.sum(masses_1)/np.sqrt(R_max**2 * 2)) * np.array([-1., 1.])/2
    
    return np.concatenate((pos_1, pos_2), axis = 0), np.concatenate((masses_1, masses_2), axis = 0), np.concatenate((vel_1, vel_2), axis = 0)


def init_milkyWay(N_part) :

    G = 1
    h = 1
    Md = 1
    zo = 0.2
    Mh = 5.8
    gamma = 1
    rc = 0.1
    Mb = 1.0/3
    a = 0.2
    c = 0.1
    Q = 1.5
    Ro = 8.5/3.5
    alpha = (1-np.sqrt(np.pi)*gamma*np.exp((gamma/rc)**2)/rc*(1-scipy.special.erf(gamma/rc)))**(-1)
    
    v_disk = 10  
    v_halos = 0.1
    v_bord = 50
    
    positions = np.zeros((N_part,2))
    masses = np.zeros(N_part)
    velocities = np.zeros((N_part, 2))

    for i in range(N_part//3) :
        #initialisation for the disk
        theta = rnd.random() * 2 * np.pi        
        r = np.abs(rnd.normal(0, h))
        positions [i] = np.array([r*np.cos(theta), r*np.sin(theta)])
        
        masses [i] = r*Md*np.exp(-r/h)/(4*np.pi*h**2)
        
        

        v = np.abs(rnd.normal(0,v_disk))
        velocities[i] = v *np.array([-np.sin(theta), np.cos(theta)])


        # initialisation for the dark holes
        theta = rnd.random() * 2 * np.pi        
        r = np.abs(rnd.normal(0, h))
        positions [N_part//3+i] = np.array([r*np.cos(theta), r*np.sin(theta)])
        masses[i+N_part//3] = r*Mh*alpha*np.exp(-r**2/rc**2)/(2*np.pi**(3/2)*rc*(r**2+gamma**2))

        """
        vr_h = G*(Mh*alpha/np.pi)**2/rc*scipy.integrate.quad(lambda x : np.exp(-(x/rc)**2)/(x**2+gamma**2)*scipy.integrate.quad(lambda y : y**2*np.exp(-y**2)/(y**2+(gamma/rc)**2), 0, x/rc), r,np.inf)
        v = 2*vr_h*np.scipy.special.erfinv(vr_h**4*0.5*np.random())
        """

        v = np.abs(rnd.normal(0, v_halos))
        velocities[i+N_part//3] = v*np.array([-np.sin(theta), np.cos(theta)])
        
        
        #initialisation for the bulges
        theta = rnd.random() * 2 * np.pi         
        x,y = np.abs(rnd.normal(0, a)), np.abs(rnd.normal(0, a))
        m = (x**2+y**2)/a**2
        positions[2*N_part//3+i] = np.array([x, y])
        masses[i+2*N_part//3] = np.sqrt(x**2+y**2)*Mb/(2*np.pi*a*m*(m+1))

        """
        v_rhd = (1/(Mh*alpha*np.exp(-r**2/rc**2)/(2*np.pi**(3/2)*rc*(r**2+gamma**2))+Md*np.exp(-r/h)/(4*np.pi*h**2)))*scipy.integrate.quad(lambda x: (Mh*alpha*np.exp(-x**2/rc**2)/(2*np.pi**(3/2)*rc*(x**2+gamma**2))+Md*np.exp(-x/h)/(4*np.pi*h**2))*G*(Md+scipy.integrate.quad(lambda y : y**2*np.exp(-y**2)/(y**2+(gamma/rc)**2), 0, x/rc))/x**2,r, np.inf) 
        v_rb = Mb/(2*np.pi*a*m*(m+1))*np.scipy.integrate.quad(lambda x : Mb/(2*np.pi*a*(x/a)**2*((x/a)**2+1))*G*Mb*(x/a)**4/(1+(x/a)**2)**2, r, np.inf)
        vr_b= np.sqrt(v_rhd+v_rb)
        v = 2*np.sqrt(2/np.pi)*np.scipy.special.erfinv(vr_b**4*0.5*np.random())/(vr_b**2)
        """
        v = np.abs(rnd.normal(0, v_bord))
        velocities[i+2*N_part//3] = v*np.array([-np.sin(theta), np.cos(theta)])

    return positions, masses, velocities

"""
positions, masses, velocities = init_milkyWay(5000)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(positions[:,0], positions[:,1], 'o', markersize = 1)
plt.show()
"""
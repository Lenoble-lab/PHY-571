from Tree import *
from simulation import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy.random as rnd
import sys
import time
import scipy.special
import scipy.integrate


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
    M_soleil = 10**1
    R_max = 2.  
    P_soleil = np.array([0.,0.])
    V_soleil = np.array([0.,0.])

    masses = np.zeros(N_part)
    velocities = np.zeros((N_part, 2))

    positions[0] = P_soleil
    masses[0] = M_soleil
    velocities[0] = V_soleil
    
    
    i = 1
    while i < N_part : 

        theta = rnd.random() * 2 * np.pi
        # r = rnd.random() * R_max
        r = np.abs(rnd.normal(0,10))
        if r!=0 : 
            positions[i] = np.array([r*np.cos(theta), r*np.sin(theta)])
            masses[i] = 1
            velocities[i] = np.sqrt(M_soleil/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])
            i += 1
  
    return positions, masses, velocities


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
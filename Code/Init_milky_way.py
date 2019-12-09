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
import random as rd


#fonction à intégrer dans la masse du halos
def f(x,q) :
    return x*np.exp(-x**2)/(x**2+q**2)

#intégration de f 
def primitive_f (r,rc,q) :
    return scipy.integrate.quad(f,0, r/rc,args=(q,))[0]

#masse totale
def M(r, Mh,Md,Mb,alpha,rc,q,a,h): 
    return Md/2+0.5*np.exp(-r/h)*Md*(r/h-1)+Mh*alpha*primitive_f(r,rc,q)/(np.sqrt(np.pi)*rc)+Mb*r*(r/a+2)/(2*a*(1+r/a)**2)

#densité du halos
def rho_h(r,Mh,alpha,rc,gamma):
    return Mh*alpha*np.exp(-r**2/rc**2)/(2*np.pi**(1.5)*rc*(r**2+gamma**2))

#densité du bulbe
def rho_b(r,Mb,a):
    return Mb/(2*np.pi*r*(1+r/a)**3)

#densité du disque
def rho_d(r,Md,h):
    return Md*np.exp(-r/h)/(4*np.pi*h**2)

#fonction à intégrer pour trouver vrh
def int_h(r,Mh,Md,Mb,alpha,rc,q,a,h,gamma):
    return rho_h(r,Mh,alpha,rc,gamma)*M(r, Mh,Md,Mb,alpha,rc,q,a,h)/r**2

#carré de la vitesse moyenne du halos
def vr_h(r,Mh,Md,Mb,alpha,rc,q,a,h,gamma,G):
    return G*scipy.integrate.quad(int_h,r,np.inf,args=(Mh,Md,Mb,alpha,rc,q,a,h,gamma,))[0]/rho_h(r,Mh,alpha,rc,gamma)

#fonction à intégrer pour trouver vr_b
def int_b(r,Mh,Md,Mb,alpha,rc,q,a,h,gamma):
    return rho_b(r,Mb,a)*M(r,Mh,Md,Mb,alpha,rc,q,a,h)/r**2

#carré de la vitesse moyenne du bulbe
def vr_b(r,Mh,Md,Mb,alpha,rc,q,a,h,gamma,G):
    return G*scipy.integrate.quad(int_b,r,np.inf,args=(Mh,Md,Mb,alpha,rc,q,a,h,gamma,))[0]/rho_b(r,Mb,a)

#carré de la vitesse du disque 
def vr_d(r,h,Q):
    return Q*np.exp(-r/h)

#densité de la vitesse
def vitesse(v,vr) :
    return np.sqrt(2/np.pi)*v**2*np.exp(-v**2/(2*vr**2))/vr**3


def init_milkyWay(N_part) :
    #constantes du problème telles que dans l'article
    G = 1
    h = 1
    Md = 1
    Mh = 5.8
    gamma = 1
    rc = 0.1
    Mb = 1.0/3
    a = 0.2
    Q = 1.5
    alpha = (1-np.sqrt(np.pi)*gamma*np.exp((gamma/rc)**2)/rc*(1-scipy.special.erf(gamma/rc)))**(-1)
    q = gamma/rc
    
    #initialisation des tableaux
    positions = np.zeros((N_part,2))
    masses = np.zeros(N_part)
    velocities = np.zeros((N_part, 2))

    for i in range(N_part//3) :
        #initialisation for the disk
        borne_sup = Md/(4*np.pi*h**2)
        X,Y = (rd.random())*100000, (rd.random())*borne_sup
        while Y >= rho_d(X,Md,h) : 
            X,Y = (rd.random())*100000, (rd.random())*borne_sup
        r = X
        theta = rnd.random() * 2 * np.pi 
        positions [i] = np.array([r*np.cos(theta), r*np.sin(theta)])
        
        masses [i] = 3*Md/N_part        
        
        v_d = np.sqrt(vr_d(r,h,Q))
        borne_sup = np.sqrt(2/np.pi)*2*np.exp(-1)/v_d
        X,Y = (rd.random())*100000, (rd.random())*borne_sup
        while Y >= vitesse(X,v_d): 
            X,Y = (rd.random())*100000, (rd.random())*borne_sup
        v = X
        velocities[i] = v*np.array([-np.sin(theta), np.cos(theta)])

        # initialisation for the halos
        #position aléatoire selon la densité rho_h grâce à la méthode du rejet
        borne_sup = Mh*alpha/(2*np.pi**1.5*rc)
        X,Y = (rd.random())*100000, (rd.random())*borne_sup
        while Y >= rho_h(X,Mh,alpha,rc,gamma): 
            X,Y = (rd.random())*100000, (rd.random())*borne_sup
        
        theta = rnd.random() * 2 * np.pi 
        r = X 
        positions [i + N_part//3] = np.array([r*np.cos(theta), r*np.sin(theta)])

        masses[i+N_part//3] = 3*Mh/N_part

        v_h = np.sqrt(vr_h(r,Mh,Md,Mb,alpha,rc,q,a,h,gamma,G))
        borne_sup = np.sqrt(2/np.pi)*2*np.exp(-1)/v_h
        X,Y = (rd.random())*100000, (rd.random())*borne_sup
        while Y >= vitesse(X,v_h) : 
            X,Y = (rd.random())*100000, (rd.random())*borne_sup
        v = X
        velocities[i] = v*np.array([-np.sin(theta), np.cos(theta)])      
        
        #initialisation for the bulges
        borne_sup = Mb/(2*np.pi)
        X,Y = (rd.random())*100000, (rd.random())*borne_sup
        while Y >= rho_b(X,Mb,a) : 
            X,Y = (rd.random())*100000, (rd.random())*borne_sup
        
        theta = rnd.random() * 2 * np.pi 
        r = X 
        positions [i + N_part//3] = np.array([r*np.cos(theta), r*np.sin(theta)])

        masses[i+2*N_part//3] = 3*Mb/N_part
        
        v_b = np.sqrt(vr_b(r,Mh,Md,Mb,alpha,rc,q,a,h,gamma,G))
        borne_sup = np.sqrt(2/np.pi)*2*np.exp(-1)/v_b
        X,Y = (rd.random())*100000, (rd.random())*borne_sup
        while Y >= vitesse(X,v_b) : 
            X,Y = (rd.random())*100000, (rd.random())*borne_sup
        v = X
        velocities[i] = v*np.array([-np.sin(theta), np.cos(theta)]) 

    return positions, masses, velocities

N_point = 3
positions, masses, velocities = init_milkyWay(N_point)


filename = "milky_way_" + str(N_point)

np.savez(filename , positions, masses, velocities)

data = []
data = np.load(filename + ".npy", allow_pickle = True)

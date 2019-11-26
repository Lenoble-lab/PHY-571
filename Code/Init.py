from Tree import *
from simulation import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy.random as rnd
import sys
import time


def init_terr_soleil(s):
    M_soleil = 10**4
    M_terre = 1
    R = 10.
    
    P_soleil = np.array([0.,0.])
    P_terre = np.array([0., R])

    V_soleil = np.array([0.,0.])
    V_terre =np.array([np.sqrt(s.G * M_soleil/R), 0.])
    s.particules = [Particule(P_soleil, M_soleil, 0, V_soleil, 0), 
                            Particule(P_terre, M_terre, 0,V_terre, 0)]

    R = 5.
    
    P_terre = np.array([R,0. ])
    V_terre =np.array([0., np.sqrt(s.G * M_soleil/R)])
    s.particules.append(Particule(P_terre, M_terre, 0,-V_terre, 0))

def init_syst_soleil(s):
    M_soleil = 10**1
    R_max = 2.  
    P_soleil = np.array([0.,0.])
    V_soleil = np.array([0.,0.])
    s.particules = [Particule(P_soleil, M_soleil, 0, V_soleil, 0)]

    N_part = 20

    for i in range (N_part):
        theta = rnd.random() * 2 * np.pi
        #r = rnd.random() * R_max
        r = np.abs(rnd.normal(0,10))
        if r!=0 : 
            pos = np.array([r*np.cos(theta), r*np.sin(theta)])

            v = np.sqrt(s.G * M_soleil/np.abs(r)) *np.array([-np.sin(theta), np.cos(theta)])
            s.particules = s.particules + [Particule(pos, 1, 0, v, i)]


def init_milkyWay(s) :

    G = 1
    h = 1
    Md = 1
    zo = 0.2
    Mh = 5.8
    gamma = 1
    rc = 10
    Mb = 1.0/3
    a = 0.2
    c = 0.1
    Q = 1.5
    Ro = = 8.5/3.5

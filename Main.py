from Tree import *
from simulation import *
import numpy as np
import matplotlib.pyplot as plt


    
simul = Simulation (0.0001, 15, np.array([0.,0.]))

simul.particules = [Particule(np.array([3,3.]), 1, 0,np.array([0.,0.]),0), 
                            Particule(np.array([-3,3.]), 1, 0,np.array([0.,0.]),1),
                            Particule(np.array([-3., -3.]), 1, 0,np.array([0.,0.]),2),
                            Particule(np.array([3., -3.]), 1, 0,np.array([0.,0.]),3)]


def init_terr_soleil(s):
    M_soleil = 10**3
    M_terre = 1
    R = 10.
    
    P_soleil = np.array([0.,0.])
    P_terre = np.array([0., R])

    V_soleil = np.array([0.,0.])
    V_terre = np.array([np.sqrt(s.G * M_soleil/R), 0.])
    s.particules = [Particule(P_soleil, M_soleil, 0, V_soleil, 0), 
                            Particule(P_terre, M_terre, 0,V_terre, 0)]


init_terr_soleil(simul)
N_part = len(simul.particules)
N_cycle = 4000

pos_x = [[]]*N_part
pos_y = [[]]*N_part


enery_cin = np.zeros(N_cycle)
enery_pot = np.zeros(N_cycle)

for j in range(N_cycle) :
    parti = simul.particules


    for i in range (N_part):
        pos_x[i].append(parti[i].position[0])
        pos_y[i].append(parti[i].position[1])

    simul.step()
    parti = simul.particules

    for i in range (N_part):
        enery_cin[j] += 1/2 * np.sum(parti[i].velocity**2) * parti[i].mass
        enery_pot[j] += parti[i].potential /2


plt.figure()
plt.subplot(1,2,1)
plt.axis('equal')
for i in range (N_part):
    plt.plot(pos_x[i], pos_y[i], 'o')

plt.subplot(1,2,2)
plt.plot(enery_cin, label = 'energy cinetic')
plt.plot(enery_pot, label = 'energy potential')
plt.plot(enery_cin + enery_pot, label = 'energy tot')

plt.legend()
plt.show()
    


from Tree import *
from simulation import *
import numpy as np
import matplotlib.pyplot as plt


    
simul = Simulation (1, 2, 15, np.array([0.,0.]))

simul.particules = [Particule(np.array([-5.,1.]), 1, 0,np.array([0.,0.]),0), 
                                Particule(np.array([5.,1.]), 1, 0,np.array([0.,0.]),1)]




pos_x = [[], []]
pos_y = [[], []]
for j in range(10) :
    print(j, 'j')
    simul.step()
    parti = simul.particules
    
    for i in range (2):
    #parti.print_particule()
        print(parti[i].force, 'force')
        pos_x[i].append(parti[i].position[0])
        pos_y[i].append(parti[i].position[1])


plt.figure()
for i in range (2):
    plt.plot(pos_x[i], pos_y[i], '-o')
plt.show()
    


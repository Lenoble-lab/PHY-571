from Tree import *
from simulation import *
import numpy as np
import matplotlib.pyplot as plt


    
simul = Simulation (1, 15, np.array([0.,0.]))

simul.particules = [Particule(np.array([0.1,3.]), 1, 0,np.array([0.,0.]),0), 
                            Particule(np.array([-0.1,3.]), 1, 0,np.array([0.,0.]),1),
                            Particule(np.array([0., -4]), 1, 0,np.array([0.,0.]),2)]


N_part = len(simul.particules)
pos_x = [[]]*N_part
pos_y = [[]]*N_part
for j in range(1) :
    parti = simul.particules

    for i in range (N_part):
        pos_x[i].append(parti[i].position[0])
        pos_y[i].append(parti[i].position[1])

    simul.step()
    parti = simul.particules
    
    for i in range (N_part):
    #parti.print_particule()
        print(parti[i].force, 'force ' + str(i))
        print(parti[i].position[0], 'pos x')


plt.figure()
for i in range (N_part):
    plt.plot(pos_x[i], pos_y[i], 'o')
plt.show()
    


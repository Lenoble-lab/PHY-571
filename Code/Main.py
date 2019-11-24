from Tree import *
from simulation import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy.random as rnd
import sys

sys.setrecursionlimit(10**5) 


simul = Simulation (0.005, 100, np.array([0.,0.]))
N_part = len(simul.particules)
N_cycle = 1000

simul.particules = [Particule(np.array([0.1,0.1]), 1, 0,np.array([0.,0.]),0), 
                            Particule(np.array([-0.1,0.1]), 1, 0,np.array([0.,0.]),1),
                            Particule(np.array([-0.1, -0.1]), 1, 0,np.array([0.,0.]),2),
                            Particule(np.array([0.1, -0.1]), 1, 0,np.array([0.,0.]),3)]


simul.particules = [Particule(np.array([0., 0.]), 1, 0,np.array([0.,0.]),0), 
                            Particule(np.array([0,1/2**2001]), 1, 0,np.array([0.,0.]),1)]

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





#init_terr_soleil(simul)
init_syst_soleil(simul)

N_part = len(simul.particules)
simul.step_second_order()
enery_pot_0, enery_cin_0 = 0, 0

for i in range (N_part):
        enery_cin_0 += 1/2 * np.sum(simul.particules[i].velocity**2) * simul.particules[i].mass
        enery_pot_0 += simul.particules[i].potential /2

pos_x = [[]]*N_part
pos_y = [[]]*N_part


enery_cin = np.zeros(N_cycle)
enery_pot = np.zeros(N_cycle)


fig = plt.figure()
frames_size = 100
ax1 = fig.add_subplot(1,2,1)
ax1.set_xlim(-frames_size, frames_size)
ax1.set_ylim(-frames_size, frames_size)
line, = ax1.plot([], [], 'o', markersize=5)

ax2 = fig.add_subplot(1,2,2)

line_ene_pot, = ax2.plot([0], [0], 'r', label = 'energy potential', linewidth = 2)
line_ene_cin, = ax2.plot([0], [0], 'b', label = 'energy cinetic', lw = 2)
line_ene_tot, = ax2.plot([0], [0], 'g', label = 'energy total', lw = 2)
ax2.set_xlim(0, 200)
ax2.set_ylim(1.2 * enery_pot_0, enery_cin_0 *1.2)
ax2.legend()

steps = []
ene_pot = []
ene_cin = []
ene_tot =[]
nt = 200
def make_frame(t):
    N_part = len(simul.particules)

    for i in range (N_part):
        parti = simul.particules

        pos_x[i] = pos_x[i] + [parti[i].position[0]]
        pos_y[i] = pos_y[i] + [parti[i].position[1]]

    if __name__ == '__main__':
        simul.step_second_order()
    parti = simul.particules
    N_part = len(simul.particules)

    for i in range (N_part):
        enery_cin[t] += 1/2 * np.sum(parti[i].velocity**2) * parti[i].mass
        enery_pot[t] += parti[i].potential /2
    
    ene_cin.append(enery_cin[t])
    ene_pot.append(enery_pot[t])
    ene_tot.append(enery_cin[t] + enery_pot[t])

    if t<=nt : steps.append(t)

    if t>nt : 
        ene_cin.pop(0)
        ene_pot.pop(0)
        ene_tot.pop(0)
    line.set_data(np.array(pos_x)[:,t], np.array(pos_y)[:,t])
    line_ene_cin.set_data(steps, ene_cin)
    line_ene_pot.set_data(steps, ene_pot)
    line_ene_tot.set_data(steps, ene_tot)
    


    print(' ')
    print(t, ' t')
    print(enery_pot[t])
    '''
    print(max(np.array(pos_x)[:,t]), "max x")
    print(min(np.array(pos_x)[:,t]), "min x")
    print(max(np.array(pos_y)[:,t]), "max y")
    print(min(np.array(pos_y)[:,t]), "min y")

    sort_pos = np.sort(np.array(pos_x)[:,t])
    print(max( [sort_pos[i] - sort_pos[i+1] for i in range (len(pos_x)-1)]), "min distane") 
    '''
        


    return (line,  line_ene_tot, line_ene_cin, line_ene_pot)

ani = animation.FuncAnimation(fig, make_frame, frames=N_cycle, interval=30, repeat = False)
ani.save('video.mp4')
plt.show()


data = [pos_x, pos_y]

#np.save("../results/data", data)

plt.figure()
plt.subplot(1,2,1)
plt.axis('equal')
for i in range (N_part):
    plt.plot(pos_x[i], pos_y[i], 'o', markersize=2)

plt.subplot(1,2,2)
plt.plot(enery_cin, label = 'energy cinetic')
plt.plot(enery_pot, label = 'energy potential')
plt.plot(enery_cin + enery_pot, label = 'energy tot')

plt.legend()
plt.savefig('test.jpg')
#plt.savefig("../results/graph_energy")

plt.show()  
    


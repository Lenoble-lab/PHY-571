import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#creation of the video


#data = np.load("../results/data.npy")
data = []

filename = "data_3000_part_0.005_deltat_collision_1129"
data = np.load(filename + ".npy", allow_pickle = True)

[pos, energy_pot, energy_cin, cintetic_momentum] = data



fig = plt.figure()

N_cycle = len(pos)
N_part = len(pos[0])


fig = plt.figure()
frames_size = 600
center = []
ax1 = fig.add_subplot(1,2,1)
ax1.set_xlim(-frames_size, frames_size)
ax1.set_ylim(-frames_size, frames_size)
line, = ax1.plot([], [], 'o', markersize=1)

ax2 = fig.add_subplot(1,2,2)

line_ene_pot, = ax2.plot([0], [0], 'r', label = 'energy potential', linewidth = 2)
line_ene_cin, = ax2.plot([0], [0], 'b', label = 'energy cinetic', lw = 2)
line_ene_tot, = ax2.plot([0], [0], 'g', label = 'energy total', lw = 2)
line_mom, = ax2.plot([0], [0], 'black', label = 'cinetic momentum', lw = 2)

ax2.set_xlim(0, 200)
# ax2.set_ylim(1,2 * (energy_cin[1] + energy_pot[1]), (energy_cin[1] + energy_pot[1]) *-12)
ax2.set_ylim(1.5 * energy_pot[0],1.5 * energy_cin[0])
ax2.legend()

nt = 200

def make_frame(i):

    if i % 100 == 0:
        print(i) 

    line.set_data(pos[i,:,0], pos[i,:,1])
    if i<nt and i>0: 
        line_ene_cin.set_data(range(0,i), energy_cin[0:i])
        line_ene_pot.set_data(range(0,i), energy_pot[0:i])
        line_ene_tot.set_data(range(0,i), energy_pot[0:i] + energy_cin[0:i])
        line_mom.set_data(range(0,i), cintetic_momentum[0:i])

    elif i>=nt : 
        line_ene_cin.set_data(range(0,nt), energy_cin[-nt +i:i])
        line_ene_pot.set_data(range(0,nt), energy_pot[-nt + i:i])
        line_ene_tot.set_data(range(0,nt), energy_pot[-nt+i:i] + energy_cin[-nt+i:i])
        line_mom.set_data(range(0,nt), cintetic_momentum[-nt + i:i])


    return (line,  line_ene_tot, line_ene_cin, line_ene_pot, line_mom)

ani = animation.FuncAnimation(fig, make_frame, frames=N_cycle, interval=30, repeat = False)
# plt.show()
ani.save(filename + "_energy.mp4")


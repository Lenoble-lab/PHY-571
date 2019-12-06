import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"""
prgm pour convertir un fichier .npy sous la forme enregristrée en vidéo avec uniquement la position des particules
"""


data = []

filename = "0312_collision_1st_order"
data = np.load(filename + ".npy", allow_pickle = True)

[pos, energy_pot, energy_cin, cintetic_momentum] = data



fig = plt.figure()

N_cycle = len(pos)
N_part = len(pos[0])


fig = plt.figure()
frames_size = 700
center = [0, 0]
plt.axis('equal')
plt.axes( xlim = (-frames_size + center[0] , +center[0] + frames_size), ylim = (-frames_size + center[1], center[1] + frames_size))
line, = plt.plot([], [], 'o', markersize=1)



def make_frame(i):
    if i%100 == 0:
        print(i)
    line.set_data(pos[i,:,0], pos[i,:,1])
    

    return (line)

ani = animation.FuncAnimation(fig, make_frame, frames=30, interval=30, repeat = False)
#plt.show()
ani.save(filename + "_star.mp4")


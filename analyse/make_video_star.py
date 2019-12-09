import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation




data = []

filename = "0612_collision_1st_order"

data = np.load(filename + ".npy", allow_pickle = True)

[pos, energy_pot, energy_cin, cintetic_momentum] = data




N_cycle = len(pos)
N_part = len(pos[0])

print(N_cycle, " N_cyles")
print(N_part, " N_part")





fig = plt.figure()

frames_size = 1000
center = [0, 200]
plt.axis('equal')
plt.axes( xlim = (-frames_size + center[0] , +center[0] + frames_size), ylim = (-frames_size + center[1], center[1] + frames_size))
line, = plt.plot([], [], 'o', markersize=1)


def make_frame(i):
    if i%100 == 0:
        print(i)
    line.set_data(pos[i,:,0], pos[i,:,1])
    

    return (line)

ani = animation.FuncAnimation(fig, make_frame, frames=len(pos), interval=30, repeat = False)
#plt.show()
ani.save(filename + "_star.mp4")


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation




data = np.load(r'C:\Users\romain Lenoble\Documents\PHY575\Projet\PHY-571\Data\data1.npy')

pos_x, pos_y = np.array(data)

N_part = len(pos_x)
N_cycle = len(pos_x[0])

pos_x = np.array(pos_x)
pos_y  = np.array(pos_y) 



fig = plt.figure()

ax = plt.axes(xlim=(-13, 13), ylim=(-13, 13))
line, = ax.plot([], [], 'o')
#line2, = ax.plot([], [], lw=3)

def make_frame(t):
    line.set_data(pos_x[:,t], pos_y[:,t])
    print(pos_x[:,t])
    return line,    

ani = animation.FuncAnimation(fig, make_frame, frames=40, interval=200)

ani.save('dynamic_images.mp4')

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"""
prgm pour convertir un fichier .npy et obetnir le graphe de l'énergie/moment cinétique en fonction du temps
"""


data = []

filename = "results/collision_5/collision_0112_1_2000"

data = np.load(filename + ".npy", allow_pickle = True)

[pos, energy_pot, energy_cin, cintetic_momentum] = data




N_cycle = len(energy_pot)
N_part = len(pos[0])

print(N_cycle, 'n_cyles')
print(N_part, "N_part")

plt.figure()
plt.title("energie et moment cinétique")
plt.plot(range(0,N_cycle), energy_pot,  'r', label = 'energy potential', linewidth = 2)
plt.plot(range(0,N_cycle), energy_cin, 'b', label = 'energy cinetic', lw = 2)
plt.plot(range(0, N_cycle), energy_cin + energy_pot , 'g', label = 'energy total', lw = 2)
plt.plot(range(0,N_cycle), 0.1 * cintetic_momentum, 'black', label = 'cinetic momentum x 10%', lw = 2)
plt.legend()
plt.savefig(filename + "_graph_energy.jpg")
plt.show()

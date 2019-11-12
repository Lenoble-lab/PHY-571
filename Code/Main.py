from Tree import *
from simultation import *
import numpy as np



simul = Simulation (0.1, 1, 10, np.array([0,0]))

simul.particules = [Particule(np.array([0,0]), 1, 0,np.array([1,1]),0)]


for i in range(10) :
    simul.step()


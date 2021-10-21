from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

x = [11, 22, 33]
y = [22, 44, 66]

file_x = open("x_axis.txt","w")
file_y = open("y_axis.txt","w")

for element in x:
    file_x.write(str(element) + "\n")
file_x.close()

for element in y:
    file_y.write(str(element) + "\n")
file_y.close()
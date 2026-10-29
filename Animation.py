import random
from itertools import count
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import MDAnalysis as mda
import argparse
plt.style.use("fivethirtyeight")

def = Main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-xy", help="Generates animated plot of xy-movement for the Lipids." )
	arg = parser.parse_args()
	

#CHL = genfromtxt("/Users/tomunterleiter/Desktop/PythonKurs/CHL.csv", delimiter=",")
u = mda.Universe("/Users/tomunterleiter/Desktop/mem.pdb")
CHL = u.select_atoms("resname CHL")

x = []
y = []
def animate(i):

	u.trajectory[i]
		
	x.append(CHL.center_of_mass(compound="residues")[:, 0])
	y.append(CHL.center_of_mass(compound="residues")[:, 1])
	
	plt.cla()

	plt.plot(x, y, label="xymove", )
	plt.title("xy Movement of CHL")
	if i == 200:
		exit()
    
ani = FuncAnimation(plt.gcf(), animate, interval=500)
    
plt.tight_layout()
plt.show()

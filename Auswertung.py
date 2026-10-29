import numpy as np
from numpy import genfromtxt
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def mk_array():
	
	CHL = genfromtxt("CHL.csv", delimiter=",")
	return CHL

def plotMSD(CHL):
	
	plt.style.use('ggplot')
	
	
#	for i in range(0, len(CHL[i, :])):
#		plt.plot(CHL[:, i], linewidth = 1.0)
#		plt.xlabel("MSD with different tau in 100 ps")
#		plt.ylabel("MSD in nm**2")
#		plt.show()

	CHL = np.delete(CHL, -1, 0) 
	x = np.arange(200)*100
	D_Fit, c = np.polyfit(x, np.mean(CHL, 1), 1)
	MSD, = plt.plot(x, np.mean(CHL, 1), label='MSD')
	Slope, = plt.plot(x, D_Fit*x+c, label='Slope')
	
	plt.legend([MSD, Slope])
	Gromacs = np.loadtxt("msd_CHOL.xvg", skiprows=21)
	plt.plot(x, Gromacs[:, 1])
	D_Fit_G, c_G = np.polyfit(x, Gromacs[:,1], 1)
	plt.plot(x, D_Fit_G*x+c_G)

	plt.title('MSD for CHL for different time steps with fitted slope')
	plt.xlabel('Tau in ps')
	plt.ylabel('MSD in cm^2/s')
	plt.legend()
#	plt.savefig('MSD.png')	 
	
	print("Diff Coeff: ", D_Fit*10**(-2)*10**5, "E-5 cm^2/s")
	plt.show()
	
CHL = mk_array()
plotMSD(CHL)


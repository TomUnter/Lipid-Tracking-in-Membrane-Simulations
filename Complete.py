import argparse
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt
from matplotlib.animation import FuncAnimation

cliparser = argparse.ArgumentParser()

cliparser.add_argument("routine", help="calculates msd for lips and saves it to a .csv file", choices={'data', 'plots', 'xy', 'xz'})
cliargs = cliparser.parse_args()

if __name__ == '__main__':
	
	if cliargs.routine == 'data':
		
#MDAnalysis topology, selecting atoms in respect to Lipid	
		
		mem = input("Please enter absolute path to the trajectory")
		u = mda.Universe(mem.pdb)
		Mem = u.select_atoms('all')
		CHL = u.select_atoms('resname CHL')
		DAPC = u.select_atoms('resname DAPC')
		DSPC = u.select_atoms('resname DSPC')

#Grouping atoms in respect to residue

		Lipid_CHL_MSD_Mean = np.zeros((201, len(CHL.residues)))
		Lipid_DAPC_MSD_Mean = np.zeros((201, len(DAPC.residues)))
		Lipid_DSPC_MSD_Mean = np.zeros((201, len(DSPC.residues)))	
	
#For MSD a sliding window is needed tau sets different start time step, 
#while win will iterate over all possible timesteps

		for tau in range(1, 201,  1):

			MSD_CHL = np.zeros((200-tau+1, len(CHL.residues)))
			MSD_DAPC = np.zeros((200-tau+1, len(DAPC.residues)))
			MSD_DSPC = np.zeros((200-tau+1, len(DSPC.residues)))

			for win in range(0, 200-tau+1):

#Calculates center of mass for every residue and subtracts com of membrane to remove com motion
#Divides by ten to get coords in nm
			
				u.trajectory[win]
				COM_CHL0 = (CHL.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
				COM_DAPC0 = (DAPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
				COM_DSPC0 = (DSPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			
				u.trajectory[win+tau]
				COM_CHL = (CHL.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
				COM_DAPC = (DAPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
				COM_DSPC = (DSPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			
# MSD_* writes generates an array with MSD-values with coloumns corresponding to every lipid molecule
# and the respective row to the timestep
			
				MSD_CHL[win,:] = (COM_CHL[:, 0] - COM_CHL0[:, 0])**2 + (COM_CHL[:, 1] - COM_CHL0[:, 1])**2
				MSD_DAPC[win,:] = (COM_DAPC[:, 0] - COM_DAPC0[:, 0])**2 + (COM_DAPC[:, 1] - COM_DAPC0[:, 1])**2
				MSD_DSPC[win,:] = (COM_DSPC[:, 0] - COM_DSPC0[:, 0])**2 + (COM_DSPC[:, 1] - COM_DSPC0[:, 1])**2
			
#Lipid_*_MSD_Mean writes to mean MSD for every lipid for every timestep
			
			print(tau)
			Lipid_CHL_MSD_Mean[tau,:] = MSD_CHL.mean(axis=0)
			Lipid_DAPC_MSD_Mean[tau,:] = MSD_DAPC.mean(axis=0)
			Lipid_DSPC_MSD_Mean[tau,:] = MSD_DSPC.mean(axis=0)

#save data in *.csv files for further needed calculation in working directory		
	
		np.savetxt("CHL.csv", Lipid_CHL_MSD_Mean, delimiter=",")
		np.savetxt("DAPC.csv", Lipid_DAPC_MSD_Mean, delimiter=",")
		np.savetxt("DSPC.csv", Lipid_DSPC_MSD_Mean, delimiter=",")
	
	if cliargs.routine == 'plots':
		
		def mk_array():
	
			DAPC = genfromtxt("CHL.csv", delimiter=",")
			DSPC = genfromtxt("DAPC.csv", delimiter=",")
			CHL = genfromtxt("DSPC.csv", delimiter=",")
			return CHL, DAPC, DSPC
		
		def plotMSD(CHL, DAPC, DSPC):
	
			plt.style.use('ggplot')
			numrows1, numcols1 = CHL.shape
			numrows2, numcols2 = DAPC.shape
			numrows3, numcols3 = DSPC.shape
			
			fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1)
			
			x = np.arange(201)*100
			
			D_Fit1, c1 = np.polyfit(x, np.mean(CHL, 1), 1)
			MSD1, = ax1.plot(x, np.mean(CHL, 1), linewidth=3, label='MSD - CHL')
			Slope1, = ax1.plot(x, D_Fit1*x+c1, label='Slope - CHL')
			
			D_Fit2, c2 = np.polyfit(x, np.mean(DAPC, 1), 1)
			MSD2, = ax2.plot(x, np.mean(DAPC, 1), linewidth=3, label='MSD - DAPC')
			Slope2, = ax2.plot(x, D_Fit2*x+c2, label='Slope - DAPC')
			
			D_Fit3, c3 = np.polyfit(x, np.mean(DSPC, 1), 1)
			MSD3, = ax3.plot(x, np.mean(DSPC, 1), linewidth=3, label='MSD - DSPC')
			Slope3, = ax3.plot(x, D_Fit3*x+c3, label='Slope - DSPC')
			
			print(numcols1)
			
	
#			for i in range(0, numcols1):
#				ax1.plot(x, CHL[:, i], "k", linewidth = 0.1)
#				
#			for j in range(0, numcols2):
#				ax2.plot(x, DAPC[:, i], "k", linewidth = 0.1)
#			
#			for k in range (0, numcols3):
#				ax3.plot(x, DSPC[:, i], "k", linewidth = 0.1)


			ax1.set_title('MSD for different time steps with fitted slope')
			ax1.set_xlabel('Tau in ps')
			ax1.set_ylabel('MSD in cm^2/s')
			ax1.legend()
			ax1.set_ylim(0, 0.2)
#			plt.savefig('MSD.png')	 
			
			print("Diff Coeff: ", D_Fit1*10**(-2)*10**5, "E-5 cm^2/s")
			plt.show()
		
		CHL, DAPC, DSPC = mk_array()
		plotMSD(CHL, DAPC, DSPC)
		
	if cliargs.routine == 'xy':
		
		Residue = input('Please enter the the residue you want to track (CHL, DAPC, or DSPC)')
		
		if Residue == 'CHL':
			leaflet = input('Enter lower, upper, or both to select the requested leaflet')
			mem = input("Please enter absolute path to the trajectory")
			u = mda.Universe(mem)
			CHL = u.select_atoms('resname CHL')
			
		

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
			
		if Residue == 'DAPC':
			mem = input("Please enter absolute path to the trajectory")
			u = mda.Universe(mem)
			DAPC = u.select_atoms('resname DAPC')
			
		

			x = []
			y = []
			def animate(i):

				u.trajectory[i]
				
				x.append(DAPC.center_of_mass(compound="residues")[:, 0])
				y.append(DAPC.center_of_mass(compound="residues")[:, 1])
			
				plt.cla()

				plt.plot(x, y, label="xymove", )
				plt.title("xy-movement of DAPC")
				if i == 200:
					exit()
			
			ani = FuncAnimation(plt.gcf(), animate, interval=500)
			
			plt.tight_layout()
			plt.show()
			
		if Residue == 'DSPC':
			mem = input("Please enter absolute path to the trajectory")
			u = mda.Universe(mem)
			DSPC = u.select_atoms('resname DAPC')
			
		

			x = []
			y = []
			def animate(i):

				u.trajectory[i]
				
				x.append(DSPC.center_of_mass(compound="residues")[:, 0])
				y.append(DSPC.center_of_mass(compound="residues")[:, 1])
			
				plt.cla()

				plt.plot(x, y, label="xymove", )
				plt.title("xy-movement of DAPC")
				if i == 200:
					exit()
			
			ani = FuncAnimation(plt.gcf(), animate, interval=500)
			
			plt.tight_layout()
			plt.show()
			
	if cliargs.routine == 'xz':
		mem = input("Please enter absolute path to the trajectory")
		u = mda.Universe(mem)
		CHL = u.select_atoms('resname CHL')
			
		

		x = []
		z = []
		def animate(i):

			u.trajectory[i]
			x.append(CHL.center_of_mass(compound="residues")[:, 0])
			z.append(CHL.center_of_mass(compound="residues")[:, 2])
			
			plt.cla()

			plt.plot(x, z, label="xymove", )
			plt.title("xz Movement of CHL")
			if i == 200:
				exit()
			
		ani = FuncAnimation(plt.gcf(), animate, interval=500)
			
		plt.tight_layout()
		plt.show() 
		
		

			
			
		
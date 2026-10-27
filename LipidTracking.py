import argparse
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt
from matplotlib.animation import FuncAnimation

cliparser = argparse.ArgumentParser()

cliparser.add_argument("routine", help="data: MSD for every lipid\n plot: plots MSD and calculates Diff. Coeff. \n xy: Tracking animation for Lipids in xy direction \n xz: Tracking Animation for Lipids in xz direction", choices={'data', 'plots', 'xy', 'xz'})
cliargs = cliparser.parse_args()

if __name__ == '__main__':
	
	if cliargs.routine == 'data':
		
#MDAnalysis topology, selecting atoms in respect to Lipid	
		
		mem = input("Please enter absolute path to the trajectory  ")
		u = mda.Universe(mem)
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

		for tau in range(1, 2,  1):

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

#Generates a plot for every lipid type containing the corresponding MSD for every time step and linear regression
		
		def mk_array():
		
#load data from .csv into an array
			CHL = genfromtxt("CHL.csv", delimiter=",")
			DAPC = genfromtxt("DAPC.csv", delimiter=",")
			DSPC = genfromtxt("DSPC.csv", delimiter=",")
			return CHL, DAPC, DSPC

#plots the curves in a figure with the linear fit through np.polyfit
		
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
			
			D_singlefit_CHL = np.zeros(numcols1)
			D_singlefit_DAPC = np.zeros(numcols2)
			D_singlefit_DSPC = np.zeros(numcols3)
			c_singlefit_CHL = np.zeros(numcols1)
			c_singlefit_DAPC = np.zeros(numcols2)
			c_singlefit_DSPC = np.zeros(numcols3)
			
			for i in range(0, numcols1):
				D_singlefit_CHL[i], c_singlefit_CHL[i] = np.polyfit(x, CHL[:, i], 1)
			
			for j in range (0, numcols2):
				D_singlefit_DAPC[j], c_singlefit_DAPC[j] = np.polyfit(x, DAPC[:, j], 1)
				
			for k in range (0, numcols3):
				D_singlefit_DSPC[k], c_singlefit_DSPC[k] = np.polyfit(x, DSPC[:, k], 1)
			
			print("Mean Diff. Coeff. of all CHLs: ", np.mean(D_singlefit_CHL)*10**(-2)*10**5, "E-5 cm^2/s")
			print("Mean Diff. Coeff. of all CHLs: ", np.mean(D_singlefit_DAPC)*10**(-2)*10**5, "E-5 cm^2/s")
			print("Mean Diff. Coeff. of all CHLs: ", np.mean(D_singlefit_DSPC)*10**(-2)*10**5, "E-5 cm^2/s")
			
			
#			for i in range(0, numcols1):
#				ax1.plot(x, CHL[:, i], "k", linewidth = 0.1)
#				
#			for j in range(0, numcols2):
#				ax2.plot(x, DAPC[:, j], "k", linewidth = 0.1)
#			
#			for k in range (0, numcols3):
#				ax3.plot(x, DSPC[:, k], "k", linewidth = 0.1)

# Those outcommented for loops would plot the MSDs for every lipid, but the plots look confusing, so they were not plotted
# additional ylimits would be needed since there are lipids that move around alot more than the average



			ax1.set_title('MSD for different time steps with fitted slope')
			ax3.set_xlabel('tau in [ps]')
			ax2.set_ylabel('MSD in [nm^2]')
			ax1.legend()
			ax2.legend()
			ax3.legend()
			
#			ax1.set_ylim(0, 0.2)

			
			 
			
			print("Diff. Coeff. for CHL with mean MSDs: ", D_Fit1*10**(-2)*10**5, "E-5 cm^2/s")
			print("Diff. Coeff. for DAPC with mean MSDs: ", D_Fit2*10**(-2)*10**5, "E-5 cm^2/s")
			print("Diff. Coeff. for DSPC with mean MSDs: ", D_Fit3*10**(-2)*10**5, "E-5 cm^2/s")
			plt.savefig('MSD.png')
			plt.show()
			
			
		CHL, DAPC, DSPC = mk_array()
		plotMSD(CHL, DAPC, DSPC)
		
		
	if cliargs.routine == 'xy':
	
#returns a tracking animation for the chosen lipid type in xy direction in both, upper, or lower leaflet
		
		Residue = input('Please enter the the residue you want to track (CHL, DAPC, or DSPC)  ')
		
		if Residue == 'CHL':
			leaflet = input('Enter lower, upper, or both to select the requested molecules  ')
			mem = input("Please enter absolute path to the trajectory  ")
			u = mda.Universe(mem)
			CHL = u.select_atoms('resname CHL')
			
		
			if leaflet == 'both':
				x = []
				y = []
				def animate(i):

					u.trajectory[i]
					
					x.append(CHL.center_of_mass(compound="residues")[:, 0])
					y.append(CHL.center_of_mass(compound="residues")[:, 1])
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of CHL-molecules in both leaflets")
					plt.xlabel('y in [angstrom]')
					plt.ylabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
				
			if leaflet == 'upper':
				
				COM_CHL = CHL.center_of_mass(compound='residues')
				mean_z = np.mean(COM_CHL[:, 2])
				upper_condition = COM_CHL[:, 2] >= mean_z
				upper = COM_CHL[upper_condition]
				
				x = []
				y = []
				
				def animate(i):
					
					
					u.trajectory[i]
					COM_CHL = CHL.center_of_mass(compound='residues')
					upperx = COM_CHL[:, 0][upper_condition]
					x.append(upperx)
					uppery = COM_CHL[:, 1][upper_condition]
					y.append(uppery)
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of CHL-molecules in upper leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
			
			if leaflet == 'lower':
				
				COM_CHL = CHL.center_of_mass(compound='residues')
				mean_z = np.mean(COM_CHL[:, 2])
				lower_condition = COM_CHL[:, 2] <= mean_z
				lower = COM_CHL[lower_condition]
				
				x = []
				y = []
				
				def animate(i):
					
					
					u.trajectory[i]
					COM_CHL = CHL.center_of_mass(compound='residues')
					lowerx = COM_CHL[:, 0][lower_condition]
					x.append(lowerx)
					lowery = COM_CHL[:, 1][lower_condition]
					y.append(lowery)
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of CHL-molecules in lower leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
				
				
			
		if Residue == 'DAPC':
			leaflet = input('Enter lower, upper, or both to select the requested molecules')
			mem = input("Please enter absolute path to the trajectory")
			u = mda.Universe(mem)
			DAPC = u.select_atoms('resname DAPC')
			
		
			if leaflet == 'both':
				x = []
				y = []
				def animate(i):

					u.trajectory[i]
					
					x.append(DAPC.center_of_mass(compound="residues")[:, 0])
					y.append(DAPC.center_of_mass(compound="residues")[:, 1])
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of DAPC-molecules")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
			
			if leaflet == 'upper':
			
				COM_DAPC = DAPC.center_of_mass(compound='residues')
				mean_z = np.mean(COM_DAPC[:, 2])
				upper_condition = COM_DAPC[:, 2] >= mean_z
				upper = COM_DAPC[upper_condition]
				
				x = []
				y = []
				
				def animate(i):
					
					
					u.trajectory[i]
					COM_DAPC = DAPC.center_of_mass(compound='residues')
					upperx = COM_DAPC[:, 0][upper_condition]
					x.append(upperx)

					uppery = COM_DAPC[:, 1][upper_condition]
					y.append(uppery)
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of DAPC-molecules in upper leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
			
			if leaflet == 'lower':
				
				COM_DAPC = DAPC.center_of_mass(compound='residues')
				mean_z = np.mean(COM_DAPC[:, 2])
				lower_condition = COM_DAPC[:, 2] <= mean_z
				lower = COM_DAPC[lower_condition]
				
				x = []
				y = []
				
				def animate(i):
					
					
					u.trajectory[i]
					COM_DAPC = DAPC.center_of_mass(compound='residues')
					lowerx = COM_DAPC[:, 0][lower_condition]
					x.append(lowerx)
					lowery = COM_DAPC[:, 1][lower_condition]
					y.append(lowery)
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of DAPC-molecules in lower leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
			
		if Residue == 'DSPC':
			leaflet = input('Enter lower, upper, or both to select the requested molecules')
			mem = input("Please enter absolute path to the trajectory")
			u = mda.Universe(mem)
			DSPC = u.select_atoms('resname DSPC')
			
		
			if leaflet == 'both':
				x = []
				y = []
				def animate(i):

					u.trajectory[i]
					
					x.append(DSPC.center_of_mass(compound="residues")[:, 0])
					y.append(DSPC.center_of_mass(compound="residues")[:, 1])
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of DSPC-molecules in upper leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
				
			if leaflet == 'upper':
							
				COM_DSPC = DSPC.center_of_mass(compound='residues')
				mean_z = np.mean(COM_DSPC[:, 2])
				upper_condition = COM_DSPC[:, 2] >= mean_z
				upper = COM_DSPC[upper_condition]
				
				x = []
				y = []
				
				def animate(i):
					
					
					u.trajectory[i]
					COM_DSPC = DSPC.center_of_mass(compound='residues')
					upperx = COM_DSPC[:, 0][upper_condition]
					x.append(upperx)

					uppery = COM_DSPC[:, 1][upper_condition]
					y.append(uppery)
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of DSPC-molecules in upper leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
			
			if leaflet == 'lower':
				
				COM_DSPC = DSPC.center_of_mass(compound='residues')
				mean_z = np.mean(COM_DSPC[:, 2])
				lower_condition = COM_DSPC[:, 2] <= mean_z
				lower = COM_DSPC[lower_condition]
				
				x = []
				y = []
				
				def animate(i):
					
					
					u.trajectory[i]
					COM_DSPC = DSPC.center_of_mass(compound='residues')
					lowerx = COM_DSPC[:, 0][lower_condition]
					x.append(lowerx)
					lowery = COM_DSPC[:, 1][lower_condition]
					y.append(lowery)
				
					plt.cla()

					plt.plot(x, y, label="xymove", )
					plt.title("xy-movement of DSPC-molecules in lower leaflet")
					plt.ylabel('y in [angstrom]')
					plt.xlabel('x in [angstrom]')
					if i == 200:
						exit()
				
				ani = FuncAnimation(plt.gcf(), animate, interval=500)
				
				plt.tight_layout()
				plt.show()
			
			
	if cliargs.routine == 'xz':
	
#returns a tracking animation for CHL in xz direction

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
			plt.title("xz Movement of CHL-molecules")
			plt.ylabel('z in [nm]')
			plt.xlabel('x in [nm]')
			if i == 200:
				exit()
			
		ani = FuncAnimation(plt.gcf(), animate, interval=500)

#Answer to question if needed. Yes on the left side one can see that a few CHLs move between leaflets
			
		plt.tight_layout()
		plt.show() 
		
		

			
			
		

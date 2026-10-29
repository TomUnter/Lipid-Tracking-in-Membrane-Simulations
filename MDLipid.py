import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import csv

if __name__ == '__main__':
	
	u = mda.Universe('mem.pdb')
	Mem = u.select_atoms('all')
	CHL = u.select_atoms('resname CHL')
	DAPC = u.select_atoms('resname DAPC')
	DSPC = u.select_atoms('resname DSPC')

	Lipid_CHL_MSD_Mean = np.zeros((201, len(CHL.residues)))
	Lipid_DAPC_MSD_Mean = np.zeros((201, len(DAPC.residues)))
	Lipid_DSPC_MSD_Mean = np.zeros((201, len(DSPC.residues)))	

	for tau in range(1, 201,  1):

		MSD_CHL = np.zeros((200-tau+1, len(CHL.residues)))
		MSD_DAPC = np.zeros((200-tau+1, len(DAPC.residues)))
		MSD_DSPC = np.zeros((200-tau+1, len(DSPC.residues)))

		for win in range(0, 200-tau+1):
			
			u.trajectory[win]
			COM_CHL0 = (CHL.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			COM_DAPC0 = (DAPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			COM_DSPC0 = (DSPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			
			u.trajectory[win+tau]
			COM_CHL = (CHL.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			COM_DAPC = (DAPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			COM_DSPC = (DSPC.center_of_mass(compound="residues") - Mem.center_of_mass()) / 10
			
			MSD_CHL[win,:] = (COM_CHL[:, 0] - COM_CHL0[:, 0])**2 + (COM_CHL[:, 1] - COM_CHL0[:, 1])**2
			MSD_DAPC[win,:] = (COM_DAPC[:, 0] - COM_DAPC0[:, 0])**2 + (COM_DAPC[:, 1] - COM_DAPC0[:, 1])**2
			MSD_DSPC[win,:] = (COM_DSPC[:, 0] - COM_DSPC0[:, 0])**2 + (COM_DSPC[:, 1] - COM_DSPC0[:, 1])**2
			
			
		
		print(tau)
		Lipid_CHL_MSD_Mean[tau,:] = MSD_CHL.mean(axis=0)
		Lipid_DAPC_MSD_Mean[tau,:] = MSD_DAPC.mean(axis=0)
		Lipid_DSPC_MSD_Mean[tau,:] = MSD_DSPC.mean(axis=0)							
	
	np.savetxt("CHL.csv", Lipid_CHL_MSD_Mean, delimiter=",")
	np.savetxt("DAPC.csv", Lipid_DAPC_MSD_Mean, delimiter=",")
	np.savetxt("DSPC.csv", Lipid_DSPC_MSD_Mean, delimiter=",")



			


		
		
		
			
			
			
		
			
	
		
#	print("diffusion coeff:", D_Fit*(10**-2)*10**5, "E-5 cm^2/s")

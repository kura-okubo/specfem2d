# make 'STATIONS' file for forward modeling with point source
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

fname="../OUTPUT_FILES_grid/externalsource.txt"
fi = open(fname,'r')
temp_df = pd.read_csv(fi,  header=None, engine='python', sep =',', comment='#')
fi.close()
temp = temp_df.as_matrix()
iele = temp[:,0]
DT = temp[:,1]
cx = temp[:,2]
cz = temp[:,3]
NumofCE = len(iele) #number of coupling elements

# add stations at locations where seismograms are validated
R_validate = 40e3 # radius of validating locations
theta_span = 10 #deg
theta = np.arange(-90,90,theta_span-1e-6) * np.pi / 180.
val_cx = R_validate*np.cos(theta - np.pi/2)
val_cz = R_validate*np.sin(theta - np.pi/2)

number_of_station = NumofCE + len(theta)
cx_all=np.append(cx, val_cx)
cz_all=np.append(cz, val_cz)

# write STATION file
with open("../DATA/STATIONS", 'w') as fo:
	for j in range(number_of_station):
		fo.write("S%04d  AA, %20.8f, %20.8f 0.0 0.0\n" % (j, cx_all[j], cz_all[j]));

# plot locations
plt.plot(cx,cz,'bv')
plt.plot(val_cx,val_cz,'rv')
plt.axes().set_aspect('equal', 'datalim')
plt.show()

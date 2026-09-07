# make 'STATIONS' file for forward modeling with point source
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

fname="../OUTPUT_FILES_grid/externalsource.txt"
fi = open(fname,'r')
temp_df = pd.read_csv(fi,  header=None, engine='python', sep =',', comment='#')
fi.close()
temp = temp_df.values
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
# Azimuth of each validation receiver, measured CLOCKWISE FROM NORTH, with
# north = +z and east = +x. theta above is only the parameter that generates the
# ring; the receivers themselves sit at (theta - 90 deg) from +x, so their true
# azimuth is 90 deg minus that. Positions are unchanged -- this is the label only.
azimuth = (90.0 - np.degrees(theta - np.pi/2)) % 360.0


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
ax=plt.gca()
ax.set_aspect('equal', 'datalim')
plt.show()

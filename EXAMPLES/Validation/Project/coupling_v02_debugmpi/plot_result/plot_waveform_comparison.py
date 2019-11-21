# plot waveform comparison between point source and acceleration injection
import numpy as np
import pandas as pd
from scipy.signal import butter, lfilter, filtfilt
import matplotlib.pyplot as plt
import os

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    #y = lfilter(b, a, data)
    y = filtfilt(b, a, data)
    return y


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
R_validate = 20e3 # radius of validating locations
theta_span = 30 #deg
theta = np.arange(0,360,theta_span) * np.pi / 180.
val_cx = R_validate*np.cos(theta)
val_cz = R_validate*np.sin(theta)

number_of_station = NumofCE + len(theta)
cx_all=np.append(cx, val_cx)
cz_all=np.append(cz, val_cz)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8), dpi=80, sharey=True)

# read waveform
itheta = 0
amp = 5e5
lw=1.5
cutoff = 50.0
t0 = 3.0

# Applying lowpass filter to eliminate high-frequency numerical oscilation
order = 6
dt = DT[0]
fs = 1/dt# sample rate, Hz
cutoff = 10.0  # desired cutoff frequency of the filter, Hz

for i in range(NumofCE, number_of_station):

	plottheta = theta[itheta] * 180 / np.pi
	itheta = itheta + 1

	with open("../OUTPUT_FILES_P/AA.S%04d.BXX.sema"%i, 'r') as fi:
		temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')

	temp	= temp_df.as_matrix()
	tx 		= temp[:,0]
	ax_P 		= temp[:,1]

	with open("../OUTPUT_FILES_P/AA.S%04d.BXZ.sema"%i, 'r') as fi:
		temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')

	temp	= temp_df.as_matrix()
	tz 		= temp[:,0]
	az_P 		= temp[:,1]

	with open("../OUTPUT_FILES_C/AA.S%04d.BXX.sema"%i, 'r') as fi:
		temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')

	temp	= temp_df.as_matrix()
	ax_C 		= temp[:,1]

	with open("../OUTPUT_FILES_C/AA.S%04d.BXZ.sema"%i, 'r') as fi:
		temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')

	temp	= temp_df.as_matrix()
	az_C 		= temp[:,1]


	# apply lowpass filter to remove numerical oscilation
	ax_P = butter_lowpass_filter(ax_P, cutoff, fs, order)
	ax_C = butter_lowpass_filter(ax_C, cutoff, fs, order)
	az_P = butter_lowpass_filter(az_P, cutoff, fs, order)
	az_C = butter_lowpass_filter(az_C, cutoff, fs, order)


	ax1.plot(tx+t0, plottheta+amp*ax_P, 'k-', lw=lw)
	ax1.plot(tx+t0, plottheta+amp*ax_C, 'r:', lw=lw)

	ax2.plot(tx+t0, plottheta+amp*az_P, 'k-', lw=lw)
	ax2.plot(tx+t0, plottheta+amp*az_C, 'r:', lw=lw)


#ax1.set_xlim(0.0, right)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Azimuth ($^{\circ}$)')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Azimuth ($^{\circ}$)')
ax2.legend(['true model (point source)', 'external source coupling'], loc=1)

ax1.set_title('x-acceleration (m/s/s)')
ax2.set_title('z-acceleration (m/s/s)')

fig.suptitle('Validation of acceleration injection: Low-Pass with %4.2f (Hz)'%(cutoff), fontsize=12)

plt.savefig("validation_full.png", dpi=300)
plt.show()

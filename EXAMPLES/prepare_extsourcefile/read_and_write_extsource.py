import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
def ricker(f0=1.0, st=0.0, et=10.0, dt=0.1, init_t=0.0, amp=1.0):
	"""
	ricker(f0, st, et, init_t,amp)
	returns time series of Ricker wavelet
	"""
	# f0=1.0
	# st=0.0
	# et=10.0
	# dt=0.1
	# init_t=5.0
	# amp=1.0
	T = et - st
	N = T/dt
	t = np.linspace(st, et, endpoint=True, num=N+1)
	rt0 = round(3 * np.sqrt(6)/(np.pi*f0)) # duration of ricker wavelet
	Nr = np.round(rt0/dt)
	if np.mod(Nr, 2) != 0:
		Nr = Nr+1
	rt = np.linspace(-rt0, rt0, num=2*Nr+1)
	a = (np.pi*f0)**2
	ricker1 = amp*(1-2*a*rt**2)*np.exp(-a*rt**2)
	stid = round(init_t/dt)
	stshift = round(rt0/dt)
	ricker = np.zeros(int(N+1))
	ricker[int(stid-stshift):int(stid-stshift+2*Nr+1)] = ricker1
	return ricker, t


fname="externalsource.txt"
fi = open(fname,'r')
temp_df = pd.read_csv(fi,  header=None, engine='python', sep =',', comment='#')
fi.close()
temp = temp_df.as_matrix()
iele = temp[:,0]
DT = temp[:,1]
cx = temp[:,2]
cz = temp[:,3]
NumofCE = len(iele) #number of coupling elements

# Parameters of sources
dt=DT[0]
st= 0.0
et= 10.0

f0=1.0
init_t=3.0 # test source init
amp=1.0

output = "./extsource"
if not os.path.isdir(output):
	os.makedirs(output)

x, t = ricker(f0, st, et, dt, init_t, amp)
# compute acceleration of x
ax = np.gradient(x, t)
az = np.gradient(x, t)

# make files
for i in range(NumofCE):
	fname = output+"/EXT%05d.dat"%iele[i]
	print(fname)
	with open(fname, 'w') as fo:
		for j in range(len(t)):
			fo.write("%20.8f, %20.8f, %20.8f\n" % (t[j], ax[j], az[j]));



plt.plot(t,x)
plt.plot(t,ax)
plt.show()
#
# reftime=UTCDateTime("2019-11-01T00")
# st = obspy.Trace()
# st.data = x
#
#
# net = "AA"
# sta = "HOSS"
# loc = "--"
# cha = "BHN"
#
# xcoord = 1000.0
# ycoord = -2000.0
# foname = "test1.sac"
#
# # st.stats['sampling_rate'] =  1.0/dt
# # st.stats['delta'] =  dt
# # st.stats['starttime'] =  reftime
# # st.stats['npts'] =  len(x)
# # st.stats['calib'] =  1.0
# # st.stats['network'] =  "AA"
# # st.stats['station'] =  "HOSS"
# # st.stats['location'] =  "--"
# # st.stats['channel'] =  "BHN"
#
# # stats.sampling_rate = 1/f0
# # st.write(foname, format="SAC")
#
# header = {'knetwk': net, 'kstnm': sta, 'kcmpnm': cha, 'stla': xcoord, 'stlo': ycoord,
#          'evla': 0.0, 'evlo': 0.0, 'evdp': 0, 'nzyear': reftime.year,
#            'nzjday': reftime.julday, 'nzhour': reftime.hour, 'nzmin': reftime.minute, 'nzsec': reftime.second,
#            'nzmsec': reftime.microsecond, 'delta': dt}
#
# sac = SACTrace(data=x, **header)
#
# sac.write(foname)
#
# sac1 = SACTrace.read(foname, debug_strings=True, headonly=False)
#
# tr1 = sac1.to_obspy_trace()

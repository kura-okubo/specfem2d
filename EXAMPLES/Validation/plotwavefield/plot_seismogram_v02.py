# Plot seismogram from large model with HOSS_SPECFEM2D coupling
# 11/2019 Kurama Okubo
#
import os
import sys
import numpy as np
import obspy
from obspy.core.trace import Trace
from obspy.core.trace import Stats
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
import time
import pandas as pd
import matplotlib.pyplot as plt

def divergence(field):
    "return the divergence of a n-D field"
    return np.sum(np.gradient(field),axis=0)

### Parameters ###
#--- File Paths ---#

inputdir = '../Project/_coupling_v05_discontinuum_largedomain/'
outputdir = inputdir + 'seismograms'

#--- HOSS non-dimensionalization info ---#
R0 = 551.92583897 #[m]
cs = 3333.3333333 #[m/s]

#--- Finite fault model geometry info ---#
L = 40000 #[m]

#--- station id info
s_sta = 0 # start id of stations
R_validate_list = [20e3, 40e3, 60e3, 80e3] # radius of validating locations
theta_span = 10.0 #deg

#--- channel parameters ---#
starttime = UTCDateTime(2019, 11, 22, 0, 0)
net = 'HOSS'
sta = 'H'
cha = 'HH'

#--- filtering parameters ---#
trim_end_time = 40.0 #6.5 discontinuum 8.0 # [s]: trimming the end of trace because it's reflection
freqmin = 0.001 #[Hz]
freqmax =5 #[Hz]

amp = 24.0 # plot amplitude
lw=1.0 #plot line width
azimuthticks = np.arange(0, 360.1, 30)
PlotUnit = 'ACC' # you can save with 'VEL', but it has to be 'ACC' for coupling

PlotType = 'raw' #'raw': acc or vel x and z components, 'divandcurl' divergence and curl of vector

fontfamily='arial'
fontsize = 20

IfSaveFig = True
PlotFigure = True # This must be false if using MPI

#######################################

#--- make output file
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

# add stations at locations where seismograms are recorded

cx_all = []
cz_all = []
number_of_station = []

theta = np.arange(0,360,theta_span) * np.pi / 180.
for R_validate in R_validate_list:
	val_cx = R_validate*np.cos(theta - np.pi/2)
	val_cz = R_validate*np.sin(theta - np.pi/2)
	number_of_station.append(len(theta))
	cx_all.append(val_cx)
	cz_all.append(val_cz)


#--- load stations

for rid in range(len(number_of_station)):
#for rid in range(1):
	# each R_validation has a stream
	stemp = Stream()
	sid = s_sta + rid*number_of_station[rid]
	eid = (rid+1)*number_of_station[rid]
	thetacount = 0
	for i in range(sid, eid):
		with open(inputdir+"/OUTPUT_FILES_C/AA.S%04d.BXX.sema"%i, 'r') as fi:
			temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')
		temp	= temp_df.as_matrix()
		tx 		= temp[:,0]
		ax 		= temp[:,1]
		sampling_rate_x=tx[1] - tx[0]

		with open(inputdir+"/OUTPUT_FILES_C/AA.S%04d.BXZ.sema"%i, 'r') as fi:
			temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')
		temp	= temp_df.as_matrix()
		tz 		= temp[:,0]
		sampling_rate_z=tz[1] - tz[0]
		az		= temp[:,1]
		#--- store traces in obspy.trace
		#------------------------------------------#
		#store data to obspy channel
		#------------------------------------------#
		# store data
		traceN = Trace(ax)
		traceE = Trace(az)
		if PlotUnit == 'VEL':
			#integrate accerelation if plotting velocity
			traceN.integrate(method='cumtrapz')
			traceE.integrate(method='cumtrapz')
		plottheta = theta[thetacount] * 180 / np.pi
		thetacount = thetacount + 1
		# store stats
		stationname = sta+'%04d'%i
		channelnameN = cha+'%s'%'N'
		channelnameE = cha+'%s'%'E'

		# for NS components
		statsN = Stats()
		statsN.sampling_rate = 1.0/sampling_rate_x
		statsN.delta =sampling_rate_x
		statsN.starttime = starttime
		statsN.npts = len(traceN.data)
		statsN.network = net
		statsN.station = stationname
		statsN.location = ''
		statsN.channel = channelnameN
		traceN.stats=statsN
		traceN.stats.sac = obspy.core.AttribDict()
		traceN.stats.sac.back_azimuth = plottheta # use this as azimuth of station

		#---applying filters---#
		traceN.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
		tN = traceN.stats.starttime
		traceN.trim(starttime=tN, endtime=tN+trim_end_time)
		traceN.taper(0.05, side='right')

		#----------------------#
		stemp += Stream(traceN)
		# for EW components
		statsE = Stats()
		statsE.sampling_rate = 1.0/sampling_rate_z
		statsE.delta = sampling_rate_z
		statsE.starttime = starttime
		statsE.npts = len(traceN.data)
		statsE.network = net
		statsE.station = stationname
		statsE.location = ''
		statsE.channel = channelnameE
		traceE.stats=statsE
		traceE.stats.sac = obspy.core.AttribDict()
		traceE.stats.sac.back_azimuth = plottheta # use this as azimuth of station
		#---applying filters---#
		traceE.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
		tE = traceE.stats.starttime
		traceE.trim(starttime=tE, endtime=tE+trim_end_time)
		traceE.taper(0.05, side='right')
		#----------------------#
		stemp += Stream(traceE)
		#print(len(stemp_N))
		print(traceN)

	# 1. plot all azimuth
	widths = [1 ,1]
	heights = [1]
	fig = plt.figure(figsize=(14, 10), dpi=80, constrained_layout=True)
	spec = fig.add_gridspec(ncols=2, nrows=1, width_ratios=widths,
	                          height_ratios=heights)
	ax = []
	for row in range(1):
		for col in range(2):
			ax.append(fig.add_subplot(spec[row, col]))

	for i in range(number_of_station[rid]):
		Nid = 2*i
		Eid = 2*i+1
		ba = stemp[Nid].stats.sac.back_azimuth

		if PlotType == 'raw':
			ax[0].plot(stemp[Eid].times(), ba+amp*stemp[Eid].data, 'k-', lw=lw)
			ax[1].plot(stemp[Nid].times(), ba+amp*stemp[Nid].data, 'k-', lw=lw)
			ax[0].set_title('ax',fontsize=fontsize, family=fontfamily)
			ax[1].set_title('az',fontsize=fontsize, family=fontfamily)


		elif PlotType == 'divandcurl':
			print("divergence and curl not implemented.")
			sys.exit(1)
			# ax[0].plot(stemp[Eid].times(), ba+amp*stemp[Eid].data, 'k-', lw=lw)
			# ax[1].plot(stemp[Nid].times(), ba+amp*stemp[Nid].data, 'k-', lw=lw)
			# ax[0].set_title(r"div $\vec a$",fontsize=fontsize, family=fontfamily)
			# ax[1].set_title(r"|curl $\vec a$|",fontsize=fontsize, family=fontfamily)


	ax[0].yaxis.set_ticks(azimuthticks)
	ax[1].yaxis.set_ticks(azimuthticks)

	for j in range(2):
		ax[j].set_xlabel('Time (s)',fontsize=fontsize, family=fontfamily)
		ax[j].set_ylabel('Azimuth ($^\circ$)',fontsize=fontsize, family=fontfamily)
		plt.setp(ax[j].get_yticklabels(), fontsize=fontsize, color="black", family=fontfamily)
		plt.setp(ax[j].get_xticklabels(), fontsize=fontsize, color="black", family=fontfamily)


	fig.suptitle('Far-field acceleration (m/s/s)',fontsize=fontsize, family=fontfamily)

	if IfSaveFig:
		plt.savefig(outputdir+"/azimuthplot_R%dkm.png"%(R_validate_list[rid]/1e3), dpi=300)

	if PlotFigure:
		plt.show()


	#output stream
	foname = outputdir+'/%s_stream.pkl'%(net)
	stemp.write(foname, format="PICKLE")




print("p convert done.")



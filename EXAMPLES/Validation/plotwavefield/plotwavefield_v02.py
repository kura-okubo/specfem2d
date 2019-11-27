"""
Plot wavefield
11/2019 Kurama Okubo
"""
# make external source file for injecting acceleration at coupling elements
import numpy as np
import math
import pandas as pd
import os
import sys
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline
import gc
import time
#from mpi4py import MPI
import matplotlib.pyplot as plt
from importParaviewColormap import importParaviewColormap

from mpi4py import MPI
from numba import jit

#----------------------------------------------------------------------------#
@jit
def plot_wavefield(fname):
	print(fname)
	tt1 = time.time()

	framenumber = int(fname[9:16])

	if ('minframe' in vars()) and ('maxframe' in vars()):
		# plot only frames within minframe and maxframe
		if (framenumber < minframe) or (framenumber > maxframe):
			return;

	# avoid overwrite
	if not Overwrite:
		foverwrite = Inputdir+"/../snapshot_wavefield/"+"wavefield%0.7d.png"%framenumber
		if os.path.isfile(foverwrite):
			return;

	# read vectors
	fi = open(Inputdir+"/"+fname, 'rb')
	vec = np.fromfile(file=fi, dtype=np.float32)
	fi.close()
	num_of_nodes = round(coord.size/2)

	# correct artifacts in velocity due to discontinuity in the source region
	vec_ref = [] # to avoid error with JIT
	if Is_correct_discountinuity:
		fname_ref="wavefield%0.7d_01.bin"%correct_discountinuity_frame
		if not os.path.isfile(Inputdir+"/"+fname_ref):
			print("error: correct artifacts reference frame does not exist.")
			sys.exit(1)
		fi = open(Inputdir+"/"+fname_ref, 'rb')
		vec_ref = np.fromfile(file=fi, dtype=np.float32)
		fi.close()


	nid = 0
	mag = np.zeros(num_of_nodes)

	for i in range(num_of_nodes):
		xid = 2*i
		zid = 2*i+1

		if plotmode == 'mag':
			mag[nid]   = math.sqrt(vec[xid]**2 + vec[zid]**2)
		elif plotmode == 'divergence':
			print("not implemented.")
		elif plotmode == 'curl':
			print("not implemented.")
		else:
			print("plotmode should be mag, divergence or curl.")
			sys.exit(1)

		if Is_correct_discountinuity and (DT*framenumber >= correct_starttime):
			if plotmode == 'mag':
				mag_ref   = math.sqrt(vec_ref[xid]**2 + vec_ref[zid]**2)
			elif plotmode == 'divergence':
				print("not implemented.")
			elif plotmode == 'curl':
				print("not implemented.")
			else:
				print("plotmode should be mag, divergence or curl.")
				sys.exit(1)
			# subtract residual velocity field caused by discontinuity in source regioin
			if (mag_ref > 0.1*cmin):
				mag[nid] = abs(mag[nid] - mag_ref)

				#mag_ref = 0.1*cmin # avoid invalid subtraction

			#mag[nid] = abs(mag[nid] - mag_ref)
			#mag[nid] = abs(mag_ref)

		if zeropad_sourceregion:

			cox = coord_x[i]
			coz = coord_z[i]

			if COUPLING_SHAPE == "circle":
				if (math.sqrt((cox-extori_x)**2 + (coz-extori_z)**2)) <= R_ext:
					mag[nid] = 0.9*cmin

			elif COUPLING_SHAPE == "rectangle":
				if (rec_xmin<=cox and rec_xmax>=cox and
					rec_zmin<=coz and rec_zmax>=coz):
					mag[nid] = 0.9*cmin

		#debug
		# eps = 1e-3;
		# cox = coord_x[i]
		# coz = coord_z[i]
		# if (math.sqrt((cox - 12e3)**2) < eps) and (math.sqrt((coz - 0.0)**2) < eps):
		# 	print("debug: %12.8e, %12.8e, %12.8e, %12.8e, %12.8e"%(cox, coz,vec_x[nid], vec_z[nid],mag[nid]))

		nid = nid + 1

	# avoid log10(0)
	mag[mag < 1e-8] = 1e-8

	#make interpolate grid
	PBox = [min(coord_x), max(coord_x), min(coord_z), max(coord_z)]
	L = PBox[1] - PBox[0]
	W = PBox[3] - PBox[2]
	xqgrid = np.linspace(PBox[0], PBox[1], plotresolution)
	zqgrid = np.linspace(PBox[2], PBox[3], plotresolution)
	xq, zq= np.meshgrid(xqgrid, zqgrid)
	tt2 = time.time()
	#zq1 = griddata((coord_x, coord_z), np.log10(mag), (xq,yq), method=method); # This is very slow (~7-8 times)
	#zq = np.nan_to_num(zq1)

	#---assembling matrix from x,y,z data
	arr = np.vstack((coord_x, coord_z, np.log10(mag)))
	rowIndex = 0
	# Sort 2D numpy array by 2nd Column
	arr1 = arr[:, arr[1].argsort()] # sort on z
	# search number of points on z axis
	for i in range(len(coord_x)):
		if i == 0:
			z1 = arr1[1,i]

		elif arr1[1,i] != z1:
			num_x = i
			break;

	xqcoord = []
	zqcoord = []
	vq	= []
	zcount = 0
	while True:
		v_vectemp = []
		# slice arr1 to sort in z
		sliceObj = slice(zcount*num_x,(zcount+1)*num_x)
		subarr = arr1[:, sliceObj]
		arr2 = subarr[:, subarr[0].argsort()] # sort on x
		if not arr2.any():
			# this is the end of matrix
			num_z = zcount
			break

		# create x and z coord for interpolation
		zqcoord.append(arr2[1,0])
		if zcount == 0:
			xqcoord.append(arr2[0,:])
			xqcoord = np.squeeze(xqcoord)

		# create value matrix
		v_vectemp = arr2[2, :]
		if zcount == 0:
			vq = v_vectemp
		else:
			vq = np.vstack((v_vectemp, vq))

		zcount = zcount + 1

	# make interpolation
	vq = np.array(vq)
	interp_spline = RectBivariateSpline(xqcoord, zqcoord, vq,  kx=3, ky=3)
	vq_interp = interp_spline(xqgrid, zqgrid)
	tt3 = time.time()
	cm = importParaviewColormap(colormapname, IscolormapInverse)

	fig = plt.figure(num=None, figsize=(12.0, 12.0), dpi=80)
	#---linear scale plot---#
	#plt.contourf(xq/1e3,yq/1e3,zq, contourlevels, cmap=cm, locator=ticker.LogLocator())
	#---logscale plot---#
	clogmin = np.log10(cmin);
	clogmax = np.log10(cmax);
	#p1 = plt.contourf(xq, yq, zq, contourlevels, cmap=cm, extend='both', norm=Normalize(vmin=clogmin, vmax=clogmax))

	levels=np.linspace(cmin, cmax, contourlevels)
	p1 = plt.contourf(xq/1e3, zq/1e3, vq_interp, np.log10(levels), cmap=cm, extend='both', vmin=clogmin, vmax=clogmax)
	#---------------------#
	#set axis
	ax = plt.gca()
	ax.set_aspect('equal', adjustable='box')
	plt.xlabel('x (km)', fontsize=fontsize, family=fontfamily)
	plt.ylabel('y (km)', fontsize=fontsize, family=fontfamily)
	t1 = DT*framenumber
	ax.set_title('Far-field radiation: T = %6.2f'%t1, fontsize=fontsize, color="black", family=fontfamily)
	plt.setp(plt.gca().get_yticklabels(), fontsize=fontsize, color="black", family=fontfamily)
	plt.setp(plt.gca().get_xticklabels(), fontsize=fontsize, color="black", family=fontfamily)

	#colorbar
	#L = np.linspace(cmin, cmax, 5) #(m/s)
	L = [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 1.0, 1.5, 2.0] #(m/s)
	Lstr = [str(round(s, 2))for s in L]
	Ltick = np.log10(L);

	#p1.set_clim(clogmin, clogmax)
	m1 = plt.cm.ScalarMappable(cmap=cm)
	m1.set_array(zq)
	m1.set_clim(clogmin, clogmax)
	cbar = fig.colorbar(m1, orientation='horizontal', shrink=0.7, format="%.2f", pad=0.1)
	cbar.set_label(r"Particle velocity magnitude (m/s)", fontsize=fontsize, family=fontfamily)
	cbar.ax.tick_params(labelsize=0.8*fontsize)
	cbar.set_ticks(Ltick)  # vertically oriented colorbar
	cbar.set_ticklabels(Lstr)  # vertically oriented colorbar

	tt4 = time.time()

	if SaveFigure:
		fodir = Inputdir+"/../snapshot_wavefield"
		if not os.path.isdir(fodir):
			os.mkdir(fodir);

		foname = fodir+"/wavefield%0.7d.png"%framenumber
		fig.savefig(foname, dpi=300, format='png', transparent=False, frameon=False)


	if PlotFigure:
		plt.show()

	#print("%s done by process %d of %d on %s." % (fname, rank, NumOfProcessors, name))

	gc.collect()
	print("Computation time: %8.4f, %8.4f, %8.4f"%(tt2-tt1, tt3-tt2, tt4-tt3))

	return 0

#----------------------------------------------------------------------------#
# MAIN PROGRAM
#----------------------------------------------------------------------------#

#---parameters---#
Inputdir = "../Project/_coupling_v05_discontinuum_largedomain/OUTPUT_FILES_C" # path to where wavefield_grid_for_dumps is located

plotresolution = 1000 # number of grid par space: large number takes long time for output
method = 'linear' #'cubic'

colormapname = 'paracooltowarm.json' #paracoolwarm pararainbow
IscolormapInverse = 0
contourlevels=400
fontfamily='arial'
fontsize = 20
cmax = 2.0 #0.4
cmin = 0.4 #0.025 * cmax

DT = 2.5e-3 # same with Par_file

# zeropadding source region
zeropad_sourceregion = True # zero padding source region
COUPLING_SHAPE       = "rectangle"      # rectangle or circle
extori_x             = 0.0            # x origin of coupling source region
extori_z             = 0.0       # z origin of coupling source region
R_ext                = 5000.0       # Radius of coupling source elements: Need to compute source within this radius by the other numerical software
dR_ext               = 250.0        # Threshold of diference in distance (R-norm(cx, cz)) to detect source element: normally dx/2 works
rec_xmin             = -11.0e3 # -10.48e3 # bottom-left corner of box
rec_zmin             = -4.0e3 #-3.5e3  #-3.1e3 # bottom-left corner of box
rec_xmax             = 11.0e3 #10.48e3  # top-right corner of box
rec_zmax             = 4.0e3 #3.5e3 #3.1e3  # top-right corner of box
rec_dx               = 120.0   # Threshold of diference in distance from box line.

# correct artifacts due to discontinuity in the source region
Is_correct_discountinuity = True
correct_discountinuity_frame = 19000 # reference frame where velocity field is equillibrium but non_zero due to fractures
correct_starttime = 6.0 # start correcting the discontinuity effect: should be after the faulting process
plotmode = 'mag' #'mag', 'divergence' and 'curl'
minframe = 0 #4380 #2340 #0
maxframe = 20 #4500 #2440 #16000

SaveFigure = True
Overwrite = True
FileFormat = 'png'
PlotFigure = False # This must be false if using MPI
#----------------#

#1. read grids
fd = open(Inputdir+"/wavefield_grid_for_dumps.bin", 'rb')
coord = np.fromfile(file=fd, dtype=np.float32)
fd.close()

num_of_nodes = round(coord.size/2)

nid = 0
coord_x = np.zeros(num_of_nodes)
coord_z = np.zeros(num_of_nodes)
for i in range(num_of_nodes):
    xid = 2*i
    zid = 2*i+1
    coord_x[nid] = coord[xid]
    coord_z[nid] = coord[zid]
    nid = nid + 1

#2. plot velocity field
#search dumped file name
flist = os.listdir(Inputdir)
fnamelist = [f for f in flist if (f[0:9]=='wavefield') and (f[9]!='_')]

start = time.time()

comm = MPI.COMM_WORLD
NumOfProcessors = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
mpiitrcount = 0

for sind, fname in enumerate(fnamelist):
	processID = sind - (NumOfProcessors * mpiitrcount)
	if rank == processID:
		plot_wavefield(fname)
		mpiitrcount = mpiitrcount + 1

elapsed_time = time.time() - start
comm.Barrier()
#Script for decompression
if rank == 0:
	print("Total_time:{0}".format(elapsed_time) + "[sec] by process %d of %d on %s.\n" % (rank, NumOfProcessors, name))
	print('All process has been successfully done.')


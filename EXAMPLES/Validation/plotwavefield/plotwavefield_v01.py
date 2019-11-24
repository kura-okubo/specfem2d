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

	# read vectors
	fi = open(Inputdir+"/"+fname, 'rb')
	vec = np.fromfile(file=fi, dtype=np.float32)
	fi.close()
	num_of_nodes = round(coord.size/2)

	nid = 0
	vec_x = np.zeros(num_of_nodes)
	vec_z = np.zeros(num_of_nodes)
	mag = np.zeros(num_of_nodes)
	for i in range(num_of_nodes):
		xid = 2*i
		zid = 2*i+1
		vec_x[nid] = vec[xid]
		vec_z[nid] = vec[zid]
		#mag[nid] = np.linalg.norm([vec[xid], vec[zid]])
		mag[nid]   = math.sqrt(vec[xid]**2 + vec[zid]**2)
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
	ax.set_title('Farfield radiation: T = %6.2f'%t1, fontsize=fontsize, color="black", family=fontfamily)
	plt.setp(plt.gca().get_yticklabels(), fontsize=fontsize, color="black", family=fontfamily)
	plt.setp(plt.gca().get_xticklabels(), fontsize=fontsize, color="black", family=fontfamily)

	#colorbar
	L = np.linspace(cmin, cmax, 5) #(m/s)
	Lstr = [str(round(s, 2))for s in L]
	Ltick = np.log10(L);

	#p1.set_clim(clogmin, clogmax)
	m1 = plt.cm.ScalarMappable(cmap=cm)
	m1.set_array(zq)
	m1.set_clim(clogmin, clogmax)
	cbar = fig.colorbar(m1, orientation='horizontal', shrink=0.7, format="%.2f")
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

	print("Computation time: %8.4f, %8.4f, %8.4f"%(tt2-tt1, tt3-tt2, tt4-tt3))

	return 0

#----------------------------------------------------------------------------#
# MAIN PROGRAM
#----------------------------------------------------------------------------#

#---parameters---#
Inputdir = "../Project/_coupling_v04_discontinuum_plotwaveform/OUTPUT_FILES_C" # path to where wavefield_grid_for_dumps is located

plotresolution = 400 # number of grid par space
method = 'linear' #'cubic'

colormapname = 'pararainbow.json'
IscolormapInverse = 0
contourlevels=400
fontfamily='arial'
fontsize = 20
cmax = 0.6
cmin = 0.2 * cmax

DT = 2.5e-3 # same with Par_file

SaveFigure = True
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


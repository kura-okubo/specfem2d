# make external source file for injecting acceleration at coupling elements
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

t0 = 6.0
# read receiver files at coupling elements and write to external source format

if not os.path.exists("../extsource"):
	os.makedirs("../extsource")

for i in range(NumofCE):
	with open("../OUTPUT_FILES_P/AA.S%04d.BXX.sema"%i, 'r') as fi:
		temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')

	temp	= temp_df.as_matrix()
	tx 		= temp[:,0]
	ax 		= temp[:,1]

	with open("../OUTPUT_FILES_P/AA.S%04d.BXZ.sema"%i, 'r') as fi:
		temp_df = pd.read_csv(fi,  header=None, engine='python', delim_whitespace=True, comment='#')

	temp	= temp_df.as_matrix()
	tz 		= temp[:,0]
	az 		= temp[:,1]


	with open("../extsource/EXT%08d.dat"%iele[i], 'w') as fo:
		for j in range(len(tx)):
			fo.write("%20.8f, %20.8e, %20.8e\n" % (tx[j], ax[j], az[j]));

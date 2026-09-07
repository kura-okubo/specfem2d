def importParaviewColormap(colormapname, IscolormapInverse):
	"""importParaviewColormap(colormapname, IscolormapInverse)
	Import colormap from paraview xml colormap file
	1. Export json file of colormap with Paraview
	2. choose colormapname
	"""

	with open('./paracolor/%s'%colormapname, 'r') as fi:

		from matplotlib.colors import LinearSegmentedColormap
		import numpy as np

		tline = fi.readline()
		samplepoints = []
		samplepointsref = 0;
		tempmap = []

		while tline:

			if 'RGBPoints' in tline:
				tline = fi.readline()

				while not ']' in tline:
					if float(tline.replace(",", "")) > samplepointsref or not samplepoints:
						samplepoints.append(float(tline.replace(",", "")))
						samplepointsref = float(tline.replace(",", ""))
						tline = fi.readline()
						r = float(tline.replace(",", ""))
						tline = fi.readline()
						g = float(tline.replace(",", ""));
						tline = fi.readline()
						b = float(tline.replace(",", ""));
						tempmap.append([r, g, b]);
						tline = fi.readline()
					else:
						tline = fi.readline()
						tline = fi.readline()
						tline = fi.readline()
						tline = fi.readline()

			tline = fi.readline()

		if IscolormapInverse:
			tempmap = np.flipud(tempmap)


		values = range(len(samplepoints))
		vmax = np.ceil(np.max(values))
		color_list = []

		for v, c in zip(values, tempmap):
			color_list.append( ( v/ vmax, c) )

		return LinearSegmentedColormap.from_list('custom_cmap', color_list)


if __name__ == "__main__":

	import matplotlib.pyplot as plt
	import numpy as np

	plt.style.use(u'seaborn-talk')
	# make these smaller to increase the resolution
	dx, dy = 0.15, 0.05

	# generate 2 2d grids for the x & y bounds
	y, x = np.mgrid[slice(-3, 3 + dy, dy),
	                slice(-3, 3 + dx, dx)]
	z = (1 - x / 2. + x ** 5 + y ** 3) * np.exp(-x ** 2 - y ** 2)
	# x and y are bounds, so z should be the value *inside* those bounds.
	# Therefore, remove the last value from the z array.
	z = z[:-1, :-1]
	z_min, z_max = -np.abs(z).max(), np.abs(z).max()

	colormapname = "pararainbow.json"
	IscolormapInverse = 0

	cm = importParaviewColormap(colormapname, IscolormapInverse)

	plt.figure
	plt.pcolor(x, y, z, cmap=cm, vmin=z_min, vmax=z_max)
	plt.title('pcolor')
	# set the limits of the plot to the limits of the data
	plt.axis([x.min(), x.max(), y.min(), y.max()])
	plt.colorbar()

	plt.show()




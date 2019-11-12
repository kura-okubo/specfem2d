module Vonkarman
# cross-correlation module
export von_karman_velocity

using Plotly, ORCA, Random, LinearAlgebra, Printf, FFTW, Distributions
#gr()

"""
    von_karman_velocity(xmin, xmax, zmin, zmax, gridsize, randomseed)
Calculate heterogeneous velocity model following Von Karman autocorrelation
# Arguments
v:: average velocity [m/s]
"""

function von_karman_velocity(xmin::Real, xmax::Real, zmin::Real, zmax::Real, gridsize::Real, v_avrg::Real;
a::Real=40.0, sig::Real=10, randomseed::Int=1234, IsPlotFig::Bool=false)

if xmin >= xmax || zmin >= zmax
	error("x,z max must be larger than x, z min")
end

if mod((xmax-xmin)/gridsize, 1) != 0.0 || mod((zmax-zmin)/gridsize, 1) != 0.0
	error("(x, zmax-x, zmin)/gridsize is not integer.")
end

NX = round(Int, (xmax - xmin)/gridsize) + 1;
NZ = round(Int, (zmax - zmin)/gridsize) + 1;

xgrid = collect(range(xmin, stop=xmax, length=NX))
zgrid = collect(range(zmin, stop=zmax, length=NZ))

rng = MersenneTwister(randomseed)
randv = rand(rng, Float64, NX, NZ)

v0 = zeros(Float64, NX, NZ)
v1 = zeros(Float64, NX, NZ)

for i = 1:NX
	for j = 1:NZ
		v0[i,j] = v_avrg
	end
end

#calculate Nyquist wave number
dkx = 1/(NX*gridsize)
dkz = 1/(NZ*gridsize)
kx_nyq = pi/(gridsize);
kz_nyq = pi/(gridsize);
kx = collect(-kx_nyq:2*pi*dkx:kx_nyq-dkx)
kz = collect(-kz_nyq:2*pi*dkz:kz_nyq-dkz)

#Von Karman Fourier spectrum

VK1=zeros(Float64, NX, NZ);
VK2=zeros(Float64, NX, NZ);
VK3=zeros(Float64, NX, NZ);

for i = 1:NX
	for j=1:NZ
		kx0 = kx[i]
		kz0 = kz[j]
		kr = sqrt(kx0^2 + kz0^2)
		#From Obermann et al. (2013) eq. 4
		#VK1[i, j] = (a^2/(1+kr^2*a^2)) * (sig/100)^2
		VK1[i, j] = (a^2/(1+kr^2*a^2))
	end
end

#Normalize VK
VK2 =  VK1 ./ maximum(VK1)
VK3 = sqrt(VK2)

data=(randv.-0.5);

fdata=fft(data);
# MULTIPLYING fdata by sqrt(C_hh)
#println("Applying filter")
newfdata=fdata.*VK3;

#FROM FOURIER TO SPACE DOMAIN
#println("Going to spacedomain (ifft)")
randdata=real(ifft(newfdata));
rdata=randdata;

#PDF method
#assume gauss velocity variation with given sigma

N_vel = 50
vmin = (1 - 3 * sig/100) * v_avrg
vmax = (1 + 3 * sig/100) * v_avrg

vel_sample = collect(range(vmin, stop=vmax, length=N_vel))
frac = zeros(Float64, N_vel)
s1 = (sig/100) * v_avrg

for i = 1:N_vel
	frac[i] = (1/(s1*sqrt(2*pi))) * exp(-0.5*((vel_sample[i]-v_avrg)/s1)^2)
	#frac[i] = rand(rng, Float64)
end

#------#

frac=cumsum(frac./sum(frac)); # normalize and cumsum
#println("Sorting data")

sdata=sort(hcat(randdata)[:], rev=false)


NN=NX*NZ;

# Calculate critical values in randdata to resemble fraction
fraclimit=zeros(length(frac)+1,1);
fraclimit[1]=sdata[1];

for n=1:length(frac)

	sid = round(Int, NN*frac[n])
	
	if sid == 0
		sid = 1
	end

	fraclimit[n+1]=sdata[sid];
end

fraclimit[1] = fraclimit[1]- 0.1*abs(fraclimit[1]);
fraclimit[end] = fraclimit[end] + 0.1*abs(fraclimit[end]);

rdata=zeros(size(randdata));

for n=1:length(frac)
	#println("Assigning velocities to fraction ", string(n))

	# OLD SLOW SOLUTION
	#[x1,z1]=find(randdata>fraclimit(n) & randdata<fraclimit(n+1) ;
	#ixz =    find(randdata>fraclimit(n) & randdata<fraclimit(n+1) ); 

	mask = zeros(Float64, NX, NZ)

	for i = 1:NX
		for j = 1:NZ
			#println(i, ",", j)
			if randdata[i,j] > fraclimit[n] && randdata[i,j] <= fraclimit[n+1]
				#println("if test")
				mask[i,j] = 1.0
			end
    	end
    end

	#global rdata = rdata + mask .* vel_sample[n];
	rdata = rdata + mask .* vel_sample[n];

end

v1 = rdata

#-------------------#



#plot result
if IsPlotFig
	data_VK = contour(x=kx, y = kz, z=real(VK3);
		colorscale="Picnic",
		contours=attr(coloring="fill",	
		showlines=false),
		ncontours=60,
		colorbar=attr(;title="Power spectral", titleside="right"),
		)

	layout = Layout(title="Von Karman",
		hight=600,
		width=600)

	plot(data_VK, layout)
	
	plot_v1 = contour(x=xgrid/1e3, y=zgrid/1e3, z=v1;
		colorscale="Portland",
		xlabel="x[km]",
		ylabel="y[km]",
		contours=attr(coloring="fill",	
		showlines=false),
		ncontours=60,
		colorbar=attr(;title="Vp [m/s]", titleside="right"),
		)
	layout = Layout(title="Vp init",
		hight=600,
		width=600)

	p = plot(plot_v1, layout)

	savefig(p, "./velocitymodel.png")

end

return v1

end

end
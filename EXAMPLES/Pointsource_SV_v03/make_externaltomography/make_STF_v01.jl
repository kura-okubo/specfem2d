"""
Make external STF file for specfem2d
05/06/2019
Kurama Okubo
"""

using Plotly, Printf, DelimitedFiles, Interpolations, ORCA

#---input parameters---
dt = 0.8e-3 #This should be same with SPECFEM2D Par_file
NSTEP = 40000 #This should be same with SPECFEM2D Par_file
Mfactor = 1e10

input_filename = "stf_withdamage.dat"
output_filename = "source.dat"

#------------------------------------------------#

T = dt * NSTEP

tout = collect(0:dt:T-dt) 


#read STF
A = readdlm(input_filename, ',', Float64, '\n')
t0 = A[:, 1]
stf0 = A[:, 2] 

#interpolate with tout
stf1 = 

# linear interpolation
interp_cubic =  LinearInterpolation(t0, stf0, extrapolation_bc = 0)
stf1 = interp_cubic(tout) # exactly log(3)

#Normalize STF
stfout = stf1 ./ maximum(stf1)

trace1 = scatter(x=tout, y=stfout, mode="lines",
	name="STF",
	linewidth=1.5,
	line_color="rgb(205, 12, 24)"
	)

layout = Layout(;title="Normalized STF",
				#xaxis_range=[0.75, 5.25],
                yaxis_range=[0, 1.2],
               )
p = plot(trace1, layout)
savefig(p, "./STF.png")

#Output external tomography file

open(output_filename, "w") do f

	for i = 1:NSTEP

		write(f,  @sprintf("%.8fd0 %.8fd0\n", tout[i], Mfactor * stfout[i]))

	end
end





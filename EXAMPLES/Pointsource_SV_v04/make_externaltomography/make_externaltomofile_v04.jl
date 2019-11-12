"""
Make external tomography file for specfem2d
05/06/2019
Kurama Okubo

---following from src of specfem2d---
!
!-------------------------------------------------------------------------------------------
!

  subroutine read_tomo_file()

! ----------------------------------------------------------------------------------------
! This subroutine reads the external ASCII tomo file TOMOGRAPHY_FILE (path to which is
! given in the Par_file).
! This file format is not very clever however it is the one used in specfem3D hence
! we chose to implement it here as well
! The external tomographic model is represented by a grid of points with assigned material
! properties and homogeneous resolution along each spatial direction x and z. The xyz file
! TOMOGRAPHY_FILE that describe the tomography should be located in the TOMOGRAPHY_PATH
! directory, set in the Par_file. The format of the file, as read from
! define_external_model_from_xyz_file.f90 looks like :
!
! ORIGIN_X ORIGIN_Z END_X END_Z
! SPACING_X SPACING_Z
! NX NZ
! VP_MIN VP_MAX VS_MIN VS_MAX RHO_MIN RHO_MAX
! x(1) z(1) vp vs rho
! x(2) z(1) vp vs rho
! ...
! x(NX) z(1) vp vs rho
! x(1) z(2) vp vs rho
! x(2) z(2) vp vs rho
! ...
! x(NX) z(2) vp vs rho
! x(1) z(3) vp vs rho
! ...
! ...
! x(NX) z(NZ) vp vs rho
!
! Where :
! _x and z must be increasing
! _ORIGIN_X, END_X are, respectively, the coordinates of the initial and final tomographic
!  grid points along the x direction (in meters)
! _ORIGIN_Z, END_Z are, respectively, the coordinates of the initial and final tomographic
!  grid points along the z direction (in meters)
! _SPACING_X, SPACING_Z are the spacing between the tomographic grid points along the x
!  and z directions, respectively (in meters)
! _NX, NZ are the number of grid points along the spatial directions x and z,
!  respectively; NX is given by [(END_X - ORIGIN_X)/SPACING_X]+1; NZ is the same as NX, but
!  for z direction.
! _VP_MIN, VP_MAX, VS_MIN, VS_MAX, RHO_MIN, RHO_MAX are the minimum and maximum values of
!  the wave speed vp and vs (in m.s-1) and of the density rho (in kg.m-3); these values
!  could be the actual limits of the tomographic parameters in the grid or the minimum
!  and maximum values to which we force the cut of velocity and density in the model.
! _After these first four lines, in the file file_name the tomographic grid points are
!  listed with the corresponding values of vp, vs and rho, scanning the grid along the x
!  coordinate (from ORIGIN_X to END_X with step of SPACING_X) for each given z (from ORIGIN_Z
!  to END_Z, with step of SPACING_Z).
! ----------------------------------------------------------------------------------------
"""

using Plotly, LinearAlgebra, Printf

include("./vonkarman.jl")

using .Vonkarman

#---input parameters---

xmin = 0 #[m]
xmax = 60000 #[m]
zmin = -60000 #[m]
zmax = 0 #[m]


vp_average = 3000; #[m/s]
rho = 2700 #[kg/m3]

#gridsize = 600;

gridsize = 600;


output_filename = "tomo_file.xyz"
Figdirroot = "./"
Filetype = "pdf"
#------------------------------------------------#

if xmin >= xmax || zmin >= zmax
	error("x,z max must be larger than x, z min")
end

spacing_x = gridsize; #[m]
spacing_z = gridsize; #[m]

NX = round(Int, (xmax - xmin)/spacing_x) + 1;
NZ = round(Int, (zmax - zmin)/spacing_z) + 1;

#calculate von karman velocity model
vp = von_karman_velocity(xmin, xmax, zmin, zmax, gridsize, vp_average; a=10.0, sig=5, randomseed=507, IsPlotFig=true)
vs = vp ./ sqrt(3)

vpmin = minimum(vp)
vpmax = maximum(vp)
vsmin = minimum(vs)
vsmax = maximum(vs)
#assume rho is constant
rhomin = rho
rhomax = rho

xgrid = collect(range(xmin, stop=xmax, length=NX))
zgrid = collect(range(zmin, stop=zmax, length=NZ))


#Output external tomography file

open(output_filename, "w") do f
	write(f, "$xmin $zmin $xmax $zmax\n")
	write(f, "$spacing_x $spacing_z\n")
	write(f, "$NX $NZ\n")
	write(f, @sprintf("%.2fd0 %.2fd0 %.2fd0 %.2fd0 %.2fd0 %.2fd0\n ", vpmin, vpmax, vsmin, vsmax, rhomin, rhomax))

	for j = 1:NZ
		for i = 1:NX

			#if (i == 3 || i == 4) && (j == 3 || j==4)
			#	write(f,  @sprintf("%.2fd0 %.2fd0 %.2fd0 %.2fd0 %.2fd0\n", xgrid[i], zgrid[j], vpmin, vsmin, rhomin))
			#else
			#	write(f,  @sprintf("%.2fd0 %.2fd0 %.2fd0 %.2fd0 %.2fd0\n", xgrid[i], zgrid[j], vpmax, vsmax, rhomax))
			#end
			write(f,  @sprintf("%.2fd0 %.2fd0 %.2fd0 %.2fd0 %.2fd0\n", xgrid[i], zgrid[j], vp[i,j], vs[i,j], rho))

		end
	end
end








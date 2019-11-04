!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================
! This subroutine is developped for SPECFEM2D-HOSS coupling
! 11/2019 Kurama Okubo (Harvard University)
!========================================================================

  subroutine determine_external_source_elements()
    !----------------------------------------------------------------------------------------------
    ! Determine which element is element where external source in imposed and which is not when the mesh has been made by the internal
    ! mesher.
    ! The elements are numbered like this :
    !
    !                              nxread
    ! <----------------------------------------------------->
    !      ______________________________________________________
    !     |          |          |    |           | nxread*nzread | ^
    !     |__________|__________|____|___________|_______________| |
    !     |   ...    |    ...   |....|    ...    |      ...      | |
    !     |__________|__________|____|___________|_______________| |
    !     | nxread+1 | nxread+2 |....|2*nxread-1 |   2*nxread    | | nzread
    !     |__________|__________|____|___________|_______________| |
    !     |     1    |    2     |....| nxread-1  |    nxread     | |
    !     |__________|__________|____|___________|_______________| v
    !
    !     !WARNING! is_pml array starts from 0
    !
    !----------------------------------------------------------------------------------------------

    use constants, only: IMAIN, OUTPUT_FILES

    use part_unstruct_par, only: elmnts, nodes_coords,nx,nz,nxread,nzread

    use shared_parameters, only: ngnod, extori_x, extori_z, R_ext, dR_ext, DT

    implicit none

    ! external functions
    integer, external :: num_4, num_9

    ! local parameters
    integer :: i,j
    integer :: num_node
    integer :: ier
    integer :: iele
    integer :: num_elmnt
    integer :: num_sourceele = 0
    integer :: nxele, nzele
    integer, dimension (ngnod) :: npele
    double precision, dimension (2, 4) :: elecoords
    double precision :: cx, cz

    !temporal parameters in input file
    ! double precision, dimension (2) :: origin_ext = (/ 0.d0, -50000.d0 /)
    ! double precision                :: R_ext = 15000.d0
    ! double precision                :: dR_ext = 5000.d0

    write(IMAIN,*)
    write(IMAIN,*) "Test determine_external_source_elements"
    write(IMAIN,*)
    write(IMAIN,*) "Size of elmnts: ", size(elmnts)
    write(IMAIN,*)

    if (ngnod == 4) then
      num_elmnt = size(elmnts)/ngnod
    else if (ngnod == 9) then
      num_elmnt = size(elmnts)/ngnod * 4
    else
      call stop_the_code('ngnod should be either 4 or 9')
    endif

    write(IMAIN,*) "Number of elements: ", num_elmnt
    write(IMAIN,*)

    write(IMAIN,*) "nx, nz ", nx, nz
    write(IMAIN,*)

    if (ngnod == 4) then
      do j = 0, nz
        do i = 0, nx
          num_node = num_4(i,j,nxread)
          write(IMAIN,*) num_node, nodes_coords(1, num_node), nodes_coords(2, num_node)
        enddo
      enddo
    else if (ngnod == 9) then
      do j = 0, nz
        do i = 0, nx
          num_node = num_9(i,j,nxread,nzread)
          write(IMAIN,*) num_node, nodes_coords(1, num_node), nodes_coords(2, num_node)
        enddo
      enddo
    else
      call stop_the_code('ngnod should be either 4 or 9')
    endif

    !output for gnuplot and more
    open(unit=20,file=trim(OUTPUT_FILES)//'gridfile_externalsource.gnu',status='unknown',iostat=ier)
    if (ier /= 0 ) then
      print *,'Error opening gnuplot file for writing: ',trim(OUTPUT_FILES)//'gridfile_externalsource.gnu'
      print *,'Please make sure directory ',trim(OUTPUT_FILES),' exists...'
      call stop_the_code('Error saving gnuplot file')
    endif

    open(unit=30,file=trim(OUTPUT_FILES)//'externalsource.txt',status='unknown',iostat=ier)
    if (ier /= 0 ) then
      print *,'Error opening gnuplot file for writing: ',trim(OUTPUT_FILES)//'gridfile_externalsource.gnu'
      print *,'Please make sure directory ',trim(OUTPUT_FILES),' exists...'
      call stop_the_code('Error saving gnuplot file')
    endif

    write(30, *) "#-----------------------------------------------#"
    write(30, *) "# List of location for External source"
    write(30, *) "# index of element, dt, cx, cy"
    write(30, *) "# How to coupling:"
    write(30, *) "# You need to obtain accerelation time series at cx, cy with dt, from other numerical solution."
    write(30, *) "# Then, please make input files named as 'EXTXXXXX.dat' (XXXXX is index of element with zero padding)."
    write(30, *) "# (Ex. EXT00012.dat)"
    write(30, *) "# The contents should be t, ax, az"
    write(30, *) "#-----------------------------------------------#"

    ! compute 4 nodes on element
    do iele = 1,num_elmnt

      ! find four nodes id
      if (ngnod == 4) then
        nzele = (iele-1)/nxread + 1
        nxele = mod((iele-1), nxread) + 1
        npele(1) = num_4(nxele-1,nzele-1,nxread)
        npele(2) = num_4(nxele,nzele-1,nxread)
        npele(3) = num_4(nxele,nzele,nxread)
        npele(4) = num_4(nxele-1,nzele,nxread)
      else if (ngnod == 9) then
        nzele = (iele-1)/(2*nxread) + 1
        nxele = mod((iele-1), 2*nxread) + 1
        npele(1) = num_9(nxele-1,nzele-1,nxread,nzread)
        npele(2) = num_9(nxele,nzele-1,nxread,nzread)
        npele(3) = num_9(nxele,nzele,nxread,nzread)
        npele(4) = num_9(nxele-1,nzele,nxread,nzread)
      endif

      write(IMAIN,*) iele, nxele, nzele, ": ", npele(1), npele(2), npele(3), npele(4)

      ! obtain coordinates
      do i = 1,4
        elecoords(1,i) = nodes_coords(1, npele(i))
        elecoords(2,i) = nodes_coords(2, npele(i))
      enddo

      cx = sum(elecoords(1,:))/4.d0
      cz = sum(elecoords(2,:))/4.d0

      write(IMAIN,*) "cx, cz: ", cx, cz

      if (abs(R_ext-sqrt((cx-extori_x)**2+(cz-extori_z)**2)) < dR_ext) then
        write(IMAIN,*) "This element goes to source element.", cx, cz, elecoords(1,1), elecoords(2,1)
        num_sourceele = num_sourceele + 1

        write(20, "(A,I6,A,1F20.8,A,1F20.8,A,1F20.8,A,1F20.8,A)") "set object ", num_sourceele, &
        " rect from ", elecoords(1,1), ",", elecoords(2,1), " to ", elecoords(1,3), ",", elecoords(2,3), &
        " fc rgb 'red'"

        write(30, "(I6, 1F20.8, 1F20.8,1F20.8)") iele, DT, cx, cz

      endif
    enddo

    close(20)
    close(30)

  end subroutine determine_external_source_elements

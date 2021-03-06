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

    use shared_parameters, only: ngnod, extori_x, extori_z, R_ext, dR_ext, &
                              rec_xmin, rec_zmin, rec_xmax, rec_zmax, rec_dx, DT, COUPLING_SHAPE

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
    double precision :: a, b, c, d


    ! write(IMAIN,*)
    ! write(IMAIN,*) "Test determine_external_source_elements"
    ! write(IMAIN,*)
    ! write(IMAIN,*) "Size of elmnts: ", size(elmnts)
    ! write(IMAIN,*)

    if (ngnod == 4) then
      num_elmnt = size(elmnts)/ngnod
    else if (ngnod == 9) then
      num_elmnt = size(elmnts)/ngnod * 4
    else
      call stop_the_code('ngnod should be either 4 or 9')
    endif

    ! write(IMAIN,*) "Number of elements: ", num_elmnt
    ! write(IMAIN,*)
    !
    ! write(IMAIN,*) "nx, nz ", nx, nz
    ! write(IMAIN,*)

    if (ngnod == 4) then
      do j = 0, nz
        do i = 0, nx
          num_node = num_4(i,j,nxread)
          !write(IMAIN,*) num_node, nodes_coords(1, num_node), nodes_coords(2, num_node)
        enddo
      enddo
    else if (ngnod == 9) then
      do j = 0, nz
        do i = 0, nx
          num_node = num_9(i,j,nxread,nzread)
          !write(IMAIN,*) num_node, nodes_coords(1, num_node), nodes_coords(2, num_node)
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
    write(30, *) "#"
    write(30, *) "# index of element, dt, cx, cz"
    write(30, *) "#"
    write(30, *) "# How to coupling:"
    write(30, *) "# You need to obtain accerelation time series at"
    write(30, *) "# cx, cz with dt, from other numerical solution."
    write(30, *) "# Then, please make input files named as"
    write(30, *) "# 'EXTXXXXXXXX.dat' (XXXXXXXX is index of element"
    write(30, *) "# with zero padding)."
    write(30, *) "# (Ex. EXT00000012.dat)"
    write(30, *) "# The contents should be t, ax, az"
    write(30, *) "# Then, please make 'extsource' directory at the same"
    write(30, *) "# directory with 'run_***.sh' (where you run the specfem)"
    write(30, *) "# and locate all EXT********.dat files in the directory."
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

      !write(IMAIN,*) iele, nxele, nzele, ": ", npele(1), npele(2), npele(3), npele(4)

      ! obtain coordinates
      do i = 1,4
        elecoords(1,i) = nodes_coords(1, npele(i))
        elecoords(2,i) = nodes_coords(2, npele(i))
      enddo

      cx = sum(elecoords(1,:))/4.d0
      cz = sum(elecoords(2,:))/4.d0

      !write(IMAIN,*) "cx, cz: ", cx, cz

      select case (trim(COUPLING_SHAPE))
      case ('circle')
          ! select coupling elements on circle
          if (abs(R_ext-sqrt((cx-extori_x)**2+(cz-extori_z)**2)) < dR_ext) then
            !write(IMAIN,*) "This element goes to source element.", cx, cz, elecoords(1,1), elecoords(2,1)
            num_sourceele = num_sourceele + 1

            write(20, "(A,I8,A,1F20.8,A,1F20.8,A,1F20.8,A,1F20.8,A)") "set object ", num_sourceele, &
            " rect from ", elecoords(1,1), ",", elecoords(2,1), " to ", elecoords(1,3), ",", elecoords(2,3), &
            " fc rgb 'red'"

            write(30, "(I8,',',1F20.8,',',1F20.8,',',1F20.8)") iele, DT, cx, cz

          endif

      case ('rectangle')
        ! select coupling elements on rectangle
        ! loop four edges
        do i = 1,4
          if (i == 1) then
            a = 0
            b= 1
            c = -rec_zmin
          elseif (i == 2) then
            a = 1
            b= 0
            c = -rec_xmax
          elseif (i == 3) then
            a = 0
            b= 1
            c = -rec_zmax
          elseif (i == 4) then
            a = 1
            b= 0
            c = -rec_xmin
          endif

          d = abs(a*cx + b*cz + c) / (sqrt(a**2+b**2))
          ! select coupling elements on circle
          if (d < rec_dx .and.&
              cx > rec_xmin-rec_dx .and.&
              cx < rec_xmax+rec_dx .and.&
              cz > rec_zmin-rec_dx .and.&
              cz < rec_zmax+rec_dx) then
            !write(IMAIN,*) "This element goes to source element.", cx, cz, elecoords(1,1), elecoords(2,1)
            num_sourceele = num_sourceele + 1

            write(20, "(A,I8,A,1F20.8,A,1F20.8,A,1F20.8,A,1F20.8,A)") "set object ", num_sourceele, &
            " rect from ", elecoords(1,1), ",", elecoords(2,1), " to ", elecoords(1,3), ",", elecoords(2,3), &
            " fc rgb 'red'"
            write(30, "(I8,',',1F20.8,',',1F20.8,',',1F20.8)") iele, DT, cx, cz
            exit
          endif
        enddo


      case default
          print *,"Error: unrecognized COUPLING_SHAPE = ",trim(COUPLING_SHAPE)
          print *,'COUPLING_SHAPE should be "circle" or "rectangle"'
          call stop_the_code('Invalid MODEL parameter')
      end select

    enddo

    close(20)
    close(30)

  end subroutine determine_external_source_elements

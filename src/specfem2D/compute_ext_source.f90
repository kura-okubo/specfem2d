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

  subroutine read_ext_source()

  ! reads in time series based on external source files
  use constants, only: EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX
  use shared_parameters, only: iele, extsource
  use shared_input_parameters, only: DT

  ! local parameters
  integer :: i, j, tid
  character*200:: line
  character*200:: extsource_filename, source_timeseries_name

  ! double precision :: DT
  ! double precision :: t_temp, ax_temp, az_temp
  ! character*200:: line
  ! character*200:: extsource_filename, source_timeseries_name
  ! integer :: ier, NumSource
  ! integer, dimension(100) :: iele
  ! integer :: i, j, tid
  ! double precision, dimension(100, 3, 100000) :: extsource

  implicit none

  ! read external source element

  extsource_filename = "./externalsource.txt"
  open(unit=31,file=trim(extsource_filename),status='old',action='read',iostat=ier)
  !if (ier /= 0) call stop_the_code('Error opening externalsource.txt, please make sure file exists...')

  NumSource = 0
  write(*,*) ier

  do while(ier == 0)

    call read_extsource_id(31, line, ier)

    write(*,*) line
    if (ier == 0) then
      NumSource = NumSource + 1
      read (line, *) iele(NumSource)
    endif

    !if (NumSource == 10000)then
    !endif
  enddo

  write(*,*) iele(1:NumSource)
  write(*,*) NumSource

  close (31)

  ! read source time series

  do i = 1, NumSource
    write(source_timeseries_name, '("extsource/EXT",I0.5,".dat")') iele(i)
    write(*,*) source_timeseries_name
    open(unit=31,file=trim(source_timeseries_name),status='old',action='read',iostat=ier)
    write(*,*) ier

    !if (ier /= 0) call stop_the_code('Error opening source file, please make sure file exists...')

    ier = 0
    tid = 0
    do while(ier == 0)
      read (31, *, iostat=ier) t_temp, ax_temp, az_temp
      if (ier == 0) then
        ! nan check
        if ((ax_temp /= ax_temp) .or. (az_temp /= az_temp)) then
          !if (ier /= 0) call stop_the_code('External source has NaN value. Please check the external source files.')
        endif
        tid = tid + 1
        extsource(i, 1, tid) = t_temp
        extsource(i, 2, tid) = ax_temp
        extsource(i, 3, tid) = az_temp
      endif
    enddo

    write(*,*) extsource(i, 3, 3000:3010)
    close (31)

  enddo


  end subroutine read_ext_source

  !========================================================================

  subroutine add_ext_source()

  ! inject the "source" from external source files

  ! use constants, only: EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX
  ! use shared_parameters, only: iele, extsource
  ! use shared_input_parameters, only: DT
  use constants, only: IMAIN

  implicit none

  ! local parameters
  ! integer :: ier, NumSource
  ! double precision :: t_temp, ax_temp, az_temp
  ! integer :: i,j, tid
  ! character*200:: line
  ! character*200:: extsource_filename, source_timeseries_name

  !allocate(iele(EXT_SOURCE_NUM_MAX),stat=ier)
  !if (ier /= 0) call stop_the_code('Error allocating array iele')
  ! double precision :: DT
  ! double precision :: t_temp, ax_temp, az_temp
  ! character*200:: line
  ! character*200:: extsource_filename, source_timeseries_name
  ! integer :: ier, NumSource
  ! integer, dimension(100) :: iele
  ! integer :: i, j, tid
  ! double precision, dimension(100, 3, 100000) :: extsource

  write(IMAIN,*)
  write(IMAIN,*) "add_ext_source Test"
  write(IMAIN,*)


  end subroutine add_ext_source


  subroutine read_extsource_id(chan, line, ier)
     integer, intent(in):: chan
     integer, intent(inout):: ier
     character*200, intent(inout):: line

  10 continue
     read (chan, *, iostat=ier) line
     if (line(1:1) .eq. '#') goto 10
     return
  end subroutine read_extsource_id

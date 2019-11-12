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
  use constants, only: IMAIN, OUTPUT_FILES, CUSTOM_REAL, EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX
  use specfem_par, only: NSTEP
  use shared_parameters, only: iele, extsource, number_of_extsource
  use shared_input_parameters, only: DT

  implicit none

  ! local parameters
  integer :: i, tid
  integer :: ier
  double precision :: t_temp, ax_temp, az_temp
  character(len=200):: line
  character(len=200):: source_timeseries_name
  ! double precision :: DT
  ! double precision :: t_temp, ax_temp, az_temp
  ! character*200:: line
  ! character*200:: extsource_filename, source_timeseries_name
  ! integer :: ier, NumSource
  ! integer, dimension(100) :: iele
  ! integer :: i, j, tid
  ! double precision, dimension(100, 3, 100000) :: extsource


  ! read external source element

  if (EXT_SOURCE_TRACE_MAX < NSTEP) then
    call stop_the_code('EXT_SOURCE_TRACE_MAX is smaller than NSTEP: increase EXT_SOURCE_TRACE_MAX.')
  endif

  allocate(iele(EXT_SOURCE_NUM_MAX))
  allocate(extsource(3,EXT_SOURCE_TRACE_MAX,EXT_SOURCE_NUM_MAX))

  ! initialize external source traces
  extsource(:,:,:) = 0._CUSTOM_REAL

  ! externalsource filename is fixed.

  open(unit=32,file=trim(OUTPUT_FILES)//'externalsource.txt',status='old',action='read',iostat=ier)
  if (ier /= 0) call stop_the_code('Error opening externalsource.txt, please make sure file exists...')

  number_of_extsource = 0

  do while(ier == 0)

    call read_extsource_id(32, line, ier)

    write(*,*) line
    if (ier == 0) then
      number_of_extsource = number_of_extsource + 1
      read (line, *) iele(number_of_extsource)
    endif

    if (number_of_extsource == EXT_SOURCE_NUM_MAX) then
      call stop_the_code('Number of external source is too many. Abort.')
    endif
  enddo

  write(*,*) iele(1:number_of_extsource)
  write(*,*) number_of_extsource

  close (32)

  ! read source time series

  do i = 1, number_of_extsource

    write(source_timeseries_name, '(A,"/EXT",I0.5,".dat")') './extsource', iele(i)
    write(*,*) source_timeseries_name
    open(unit=32,file=trim(source_timeseries_name),status='old',action='read',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening external source file, please make sure file exists...')
    write(*,*) ier

    ier = 0
    tid = 0
    do while(ier == 0)
      read (32, *, iostat=ier) t_temp, ax_temp, az_temp
      if (ier == 0) then
        ! nan check
        if ((ax_temp /= ax_temp) .or. (az_temp /= az_temp)) then
          call stop_the_code('External source has NaN value. Please check the external source files.')
        endif
        tid = tid + 1
        extsource(1, tid, i) = t_temp
        extsource(2, tid, i) = ax_temp
        extsource(3, tid, i) = az_temp
      endif
    enddo

    ! timestep check
    if ((extsource(1, 2, i) - extsource(1, 1, i)) /= DT) then
      call stop_the_code('Time step of external source is diferent with input file. Please check the time step.')
    endif

    write(*,*) extsource(2, 3000:3010, i)
    close (32)

  enddo


  write(IMAIN,*)
  write(IMAIN,*) '*****************************************************'
  write(IMAIN,*) '*** External source for injection of acceleration ***'
  write(IMAIN,*) '*****************************************************'
  write(IMAIN,*)
  write(IMAIN,*) '*** Number of external source = ',number_of_extsource
  write(IMAIN,*) '*** Time series length of external source = ',tid
  write(IMAIN,*)

  end subroutine read_ext_source

  !========================================================================

  subroutine add_ext_source(accel_elastic)

  ! inject the "source" from external source files

  ! use constants, only: EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX
  ! use shared_parameters, only: iele, extsource
  ! use shared_input_parameters, only: DT


  use constants, only: CUSTOM_REAL, IMAIN, NDIM, EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX, &
                      NGLLX,NGLLZ
  use specfem_par, only: nglob, P_SV, it, ibool, jacobian, NSTEP, deltat, &
                        wxgll,wzgll
  use shared_parameters, only: iele, extsource, number_of_extsource

  implicit none

  ! local parameters
  integer :: i, j, k, iglob, eleid, length_of_timeseries
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic
  integer, dimension(2,number_of_extsource*NGLLZ*NGLLX) :: iglob_add_accel_smooth_idlist
  integer :: iadd, ioverlap
  logical :: is_iglob_overlap
  !integer, dimension(2, number_of_extsource*NGLLZ*NGLLX) :: iglob_add_accel_smooth
  double precision :: t0


  t0 = ((NSTEP-1)/2.0_CUSTOM_REAL) * deltat
  length_of_timeseries = size(extsource(1,:,1))

  if (it == 1) then
    write(IMAIN,*) "DEBUG: t0 IS ", t0
    write(IMAIN,*) "DEBUG: length_of_timeseries IS ", length_of_timeseries
  endif

  !timeval = (it-1) * deltat
  iadd = 0
  ioverlap = 0

  do eleid = 1, number_of_extsource
    if (P_SV) then
      ! P-SV calculation
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,iele(eleid))

          is_iglob_overlap = .false.
          ! find if the acceleration is already added to the iglob
          ! this process is essentially only needed for once in the beginning, but
          ! due to the readablility and not so efficient even if implementing subroutine
          ! for the initialization of iglob_add_accel_smooth_idlist, we find it every time step.
          ! iglob_add_accel_smooth_idlist (1,:) = iglob id
          ! iglob_add_accel_smooth_idlist (2,:) = number of adding acceleration at the iglob

          do k = 1, iadd
              if (iglob_add_accel_smooth_idlist(1, k) .eq. iglob) then
                  is_iglob_overlap = .true.
                  iglob_add_accel_smooth_idlist(2, k) = iglob_add_accel_smooth_idlist(2, k) + 1
                  exit
              endif
          end do

          if (.not. is_iglob_overlap) then
            ! this iglob is first time to add
            accel_elastic(1,iglob) = extsource(1,it,eleid) * wxgll(i)*wzgll(j)*jacobian(i,j,iele(eleid))
            accel_elastic(2,iglob) = extsource(2,it,eleid) * wxgll(i)*wzgll(j)*jacobian(i,j,iele(eleid))
            iadd = iadd + 1
            iglob_add_accel_smooth_idlist(1, iadd) = iglob
            iglob_add_accel_smooth_idlist(2, iadd) = 1
          else
            ! this iglob already has value: it will be averaged
            accel_elastic(1,iglob) = accel_elastic(1,iglob) + extsource(1,it,eleid) * wxgll(i)*wzgll(j)*jacobian(i,j,iele(eleid))
            accel_elastic(2,iglob) = accel_elastic(2,iglob) + extsource(2,it,eleid) * wxgll(i)*wzgll(j)*jacobian(i,j,iele(eleid))
          endif

        enddo
      enddo
    else
      ! SH (membrane) calculation
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,iele(eleid))
          is_iglob_overlap = .false.

          do k = 1, iadd
              if (iglob_add_accel_smooth_idlist(1, k) .eq. iglob) then
                  is_iglob_overlap = .true.
                  iglob_add_accel_smooth_idlist(2, k) = iglob_add_accel_smooth_idlist(2, k) + 1
                  exit
              endif
          end do

          if (.not. is_iglob_overlap) then
            ! this iglob is first time to add
            accel_elastic(1,iglob) = extsource(1,it,eleid) * wxgll(i)*wzgll(j)*jacobian(i,j,iele(eleid))
            iadd = iadd + 1
            iglob_add_accel_smooth_idlist(1, iadd) = iglob
            iglob_add_accel_smooth_idlist(2, iadd) = 1
          else
            ! this iglob already has value: it will be averaged
            accel_elastic(1,iglob) = accel_elastic(1,iglob) + extsource(1,it,eleid) * wxgll(i)*wzgll(j)*jacobian(i,j,iele(eleid))
          endif

        enddo
      enddo
    endif
  enddo

  ! take an average based on the number of adding acceleration
  do  k = 1,iadd
    if (P_SV) then
      !write(IMAIN,*) iglob_add_accel_smooth_idlist(1, k), iglob_add_accel_smooth_idlist(2, k)
      accel_elastic(1,iglob_add_accel_smooth_idlist(1, k)) = &
      accel_elastic(1,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
      accel_elastic(2,iglob_add_accel_smooth_idlist(1, k)) = &
      accel_elastic(2,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
    else
      accel_elastic(1,iglob_add_accel_smooth_idlist(1, k)) = &
      accel_elastic(1,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
    endif
  enddo

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

  ! write(IMAIN,*)
  ! write(IMAIN,*) "add_ext_source Test"
  ! write(IMAIN,*)

  ! if (P_SV) then
  !   ! P-SV simulation
  !   do j = 1,NGLLZ
  !     do i = 1,NGLLX
  !       ! iglob = ibool(i,j,ispec_noise)
  !       source_array_noise(1,i,j,:) = time_function_noise(:) * hxi(i) * hgamma(j)
  !       source_array_noise(2,i,j,:) = time_function_noise(:) * hxi(i) * hgamma(j)
  !     enddo
  !   enddo
  ! else
  !   ! SH (membrane) simulation
  !   do j = 1,NGLLZ
  !     do i = 1,NGLLX
  !       ! iglob = ibool(i,j,ispec_noise)
  !       source_array_noise(1,i,j,:) = time_function_noise(:) * hxi(i) * hgamma(j)
  !     enddo
  !   enddo
  ! endif
  !
  !


  ! do ispec = 1, nspec
  !   do j = 1, NGLLZ
  !     do i = 1, NGLLX
  !       iglob = ibool(i,j,ispec)
  !       if (P_SV) then
  !         ! P-SV calculation
  !         !accel_elastic(1,iglob) = accel_elastic(1,iglob) + surface_movie_x_noise(iglob) * &
  !         !                         mask_noise(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
  !         ! only vertical for now...
  !         accel_elastic(2,iglob) = accel_elastic(2,iglob) + surface_movie_y_or_z_noise(iglob) * &
  !                                  mask_noise(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
  !
  !       else
  !         ! SH (membrane) calculation
  !         accel_elastic(1,iglob) = accel_elastic(1,iglob) + surface_movie_y_or_z_noise(iglob) * &
  !                                  mask_noise(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
  !       endif
  !     enddo
  !   enddo
  ! enddo


  !local
  ! integer :: i,j,iglob
  !
  ! if (P_SV) then
  !   ! P-SV calculation
  !   do j = 1,NGLLZ
  !     do i = 1,NGLLX
  !       iglob = ibool(i,j,ispec_noise)
  !       accel_elastic(1,iglob) = accel_elastic(1,iglob) + sin(angle_noise)*source_array_noise(1,i,j,it)
  !       accel_elastic(2,iglob) = accel_elastic(2,iglob) - cos(angle_noise)*source_array_noise(2,i,j,it)
  !     enddo
  !   enddo
  ! else
  !   ! SH (membrane) calculation
  !   do j = 1,NGLLZ
  !     do i = 1,NGLLX
  !       iglob = ibool(i,j,ispec_noise)
  !       accel_elastic(1,iglob) = accel_elastic(1,iglob) - source_array_noise(1,i,j,it)
  !     enddo
  !   enddo
  ! endif

  end subroutine add_ext_source


  subroutine read_extsource_id(chan, line, ier)
     integer, intent(in):: chan
     integer, intent(inout):: ier
     character(len=200), intent(inout):: line

  10 continue
     read (chan, *, iostat=ier) line
     if (line(1:1) .eq. '#') goto 10
     return
  end subroutine read_extsource_id

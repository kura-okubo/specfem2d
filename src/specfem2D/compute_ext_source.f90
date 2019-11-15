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
! added this subroutine by Kurama Okubo 11/2019
!========================================================================

  subroutine read_ext_source_num()

    ! reads in time series based on external source files
    use constants, only: IMAIN, OUTPUT_FILES, EXT_SOURCE_NUM_MAX, MAX_STRING_LEN
    use specfem_par, only: myrank, NSTEP
    use shared_parameters, only: iele, number_of_extsource, NPROC

    implicit none

    ! local parameters
    integer :: i, j, tid
    integer :: ier
    double precision :: t_temp, ax_temp, az_temp
    character(len=MAX_STRING_LEN):: line
    integer, dimension(1) :: ibuf_num

    number_of_extsource = 0

    allocate(iele(EXT_SOURCE_NUM_MAX))

    if (myrank==0) then
      ! read external source element
      ! externalsource filename is fixed.
      open(unit=32,file=trim(OUTPUT_FILES)//'externalsource.txt',status='old',action='read',iostat=ier)
      if (ier /= 0) call stop_the_code('Error opening externalsource.txt, please make sure file exists...')

      do while(ier == 0)

        call read_extsource_id(32, line, ier)

        if (ier == 0) then
          number_of_extsource = number_of_extsource + 1
          read (line, *) iele(number_of_extsource)
        endif

        if (number_of_extsource == EXT_SOURCE_NUM_MAX) then
          call stop_the_code('Number of external source is too many. Abort.')
        endif

      enddo

      close (32)

      write(IMAIN,*)
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*** External source for injection of acceleration ***'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*)
      write(IMAIN,*) '*** Number of external source = ',number_of_extsource
      write(IMAIN,*) '*** Time series length of external source = ',tid
      write(IMAIN,*)
    endif

    if (myrank == 0) then
      ibuf_num(1) = number_of_extsource
      do i = 1, NPROC-1
        call send_singlei(ibuf_num, i, 1)
        do j = 1, number_of_extsource
          call send_i(iele, j, i, 2)
        enddo
      enddo
    else
      call recv_singlei(ibuf_num, 0, 1)
      number_of_extsource = ibuf_num(1)
      do j = 1, number_of_extsource
        call recv_i(iele,j, 0, 2)
      enddo
    endif

    write(IMAIN, *) "myrank, number_of_extsource: ", myrank, number_of_extsource

  end subroutine read_ext_source_num

  !========================================================================

  subroutine add_ext_source(accel_elastic)

  ! inject the "source" from external source files
  use constants, only: CUSTOM_REAL, IMAIN, NDIM, EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX, &
                      NGLLX,NGLLZ
  use specfem_par, only: nglob, P_SV, it, ibool, NSTEP, deltat, myrank, nspec,&
                    nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic, &
                    glob2loc_table

  use shared_parameters, only: iele, extsource, number_of_extsource

  implicit none

  ! local parameters
  integer :: i, j, k, iglob, eleid, length_of_timeseries
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic
  integer, dimension(2,number_of_extsource*NGLLZ*NGLLX) :: iglob_add_accel_smooth_idlist
  integer :: iadd, ioverlap
  logical :: is_iglob_overlap
  double precision :: t0

  !double precision :: dummy
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  !integer :: iphase
  integer :: iphase, num_elements
  integer, allocatable, dimension(:) :: iphase_elelist
  integer :: ispec_p, num_modified_ele
  logical :: is_ispec_in
  integer :: loc_rank, loc_eleid


  !write(IMAIN,*) "test 1"

  ! Initialized full extsource array
  if (it == 1) then
    allocate(extsource(3,EXT_SOURCE_TRACE_MAX,EXT_SOURCE_NUM_MAX))
    extsource(:,:,:) = 0._CUSTOM_REAL
  endif


  t0 = ((NSTEP-1)/2.0_CUSTOM_REAL) * deltat
  length_of_timeseries = size(extsource(1,:,1))

  !timeval = (it-1) * deltat
  iadd = 0
  ioverlap = 0
  !
  ! do eleid = 1, number_of_extsource
  !   write(IMAIN, *) iele(eleid)
  ! enddo
  num_modified_ele = 0

  !write(IMAIN, *) " myrank, iphase, number_of_extsource: ", myrank, iphase, number_of_extsource
  !
  ! write(IMAIN, *) "#my rank = 0 -----------------------------#"
  ! if (myrank == 0) then
  !   do i = 1, 40000
  !     write(IMAIN, *) glob2loc_table(i, 1), glob2loc_table(i, 2), glob2loc_table(i, 3)
  !   enddo
  ! endif
  !

  if (it == 1) call setup_glob2loctable()

  ! if (myrank == 1) then
  ! write(IMAIN, *) "#my rank = 1 -----------------------------#"
  !
  !   do i = 1, nspec
  !     write(IMAIN, *) glob2loc_table(i, 1), glob2loc_table(i, 2), glob2loc_table(i, 3)
  !   enddo
  ! endif

  !do iphase = 1, 2
  do iphase = 1, 2
    write(IMAIN, *) "test000"

    !write(IMAIN,*) "test -1", myrank

    ! choses inner/outer elements
    if (iphase == 1) then
      num_elements = nspec_outer_elastic
    else
      num_elements = nspec_inner_elastic
    endif

    allocate(iphase_elelist(num_elements))

    ! make list of outer/inner elements
    do ispec_p = 1,num_elements
      iphase_elelist(ispec_p) = phase_ispec_inner_elastic(ispec_p,iphase)
    enddo

    write(IMAIN, *) "myrank, number_of_extsource: ", myrank, number_of_extsource

    do eleid = 1, number_of_extsource

      ! check if iele(eleid) is in the ispec_plist
      ! is_ispec_in = .false.
      ! do ispec_p = 1,num_elements
      !   if (iele(eleid) == phase_ispec_inner_elastic(ispec_p,iphase)) then
      !     is_ispec_in = .true.
      !     write(IMAIN, *) "myrank, iphase, modifiediele: ", myrank, iphase, iele(eleid)
      !     exit
      !   endif
      ! enddo

      ! if (is_ispec_in) then
      !write(IMAIN, *) iele(eleid)

      ! find local element id in each MPI domain from glob2loc_table
      is_ispec_in = .false.

      do i = 1, nspec
        if (glob2loc_table(i, 1) == iele(eleid)) then
          loc_rank = glob2loc_table(i, 2)
          loc_eleid = glob2loc_table(i, 3)
          ! check if this local element is inner or outer
          do ispec_p = 1,num_elements
            if (iphase_elelist(ispec_p) == loc_eleid) then
              is_ispec_in = .true.
              num_modified_ele = num_modified_ele + 1
              if (it == 1) then
                call read_ext_source(extsource, eleid, iele(eleid))
                write(IMAIN, *) "EXT", iele(eleid), ".dat is read."
              endif
              !write(IMAIN, *) extsource(2, 500:600, eleid)
              !write(*,*) loc_rank, loc_eleid
              exit
            endif
          enddo
        endif
      enddo

       !write(IMAIN,*) "test 0", myrank

      if (is_ispec_in) then
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX

              ! write(IMAIN,*) "test 1", myrank

              iglob = ibool(i,j,loc_eleid)

               !write(IMAIN,*) "test 2", myrank

              !write(IMAIN,*) "test 2",myrank,eleid,i,j

              is_iglob_overlap = .false.
              ! find if the acceleration is already added to the iglob
              ! this process is essentially only needed for once in the beginning, but
              ! due to the readablility and not so efficient even if implementing subroutine
              ! for the initialization of iglob_add_accel_smooth_idlist, we find it every time step.
              ! iglob_add_accel_smooth_idlist (1,:) = iglob id
              ! iglob_add_accel_smooth_idlist (2,:) = number of adding acceleration at the iglob
              !write(IMAIN,*) "test 3", myrank

              do k = 1, iadd
                  !write(IMAIN,*) "test 2.1", myrank, iele(eleid)
                  if (iglob_add_accel_smooth_idlist(1, k) .eq. iglob) then
                      is_iglob_overlap = .true.
                      iglob_add_accel_smooth_idlist(2, k) = iglob_add_accel_smooth_idlist(2, k) + 1
                      !write(IMAIN,*) "test 2.3", myrank
                      exit
                  endif
                  !write(IMAIN,*) "test 2.2", myrank
              end do
               !write(IMAIN,*) "test 4", myrank

              ! write(IMAIN,*) "test 3", myrank

              if (.not. is_iglob_overlap) then
                ! this iglob is first time to add


                !write(IMAIN, *) accel_elastic(1,iglob) !, extsource(2,it,eleid)



                !write(IMAIN, *) accel_elastic(2,iglob) !, extsource(3,it,eleid)
                 !write(IMAIN,*) "test 5", myrank
                accel_elastic(1,iglob) = extsource(2,it,eleid)
                 !write(IMAIN,*) "test 6", myrank

                accel_elastic(2,iglob) = extsource(3,it,eleid)
                iadd = iadd + 1
                ! write(IMAIN,*) "test 5", myrank
                iglob_add_accel_smooth_idlist(1, iadd) = iglob
                iglob_add_accel_smooth_idlist(2, iadd) = 1
                ! write(IMAIN,*) "test 6", myrank

              else
                ! this iglob already has value: it will be averaged
                ! write(IMAIN,*) "test 7", myrank
                 !write(IMAIN,*) "test 7", myrank

                accel_elastic(1,iglob) = accel_elastic(1,iglob) + extsource(2,it,eleid)
                 !write(IMAIN,*) "test 8", myrank

                accel_elastic(2,iglob) = accel_elastic(2,iglob) + extsource(3,it,eleid)
                 !write(IMAIN,*) "test 9", myrank

              endif

              !write(IMAIN,*) "test 4"

            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,loc_eleid)
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
                accel_elastic(1,iglob) = extsource(2,it,eleid)
                iadd = iadd + 1
                iglob_add_accel_smooth_idlist(1, iadd) = iglob
                iglob_add_accel_smooth_idlist(2, iadd) = 1

              else
                ! this iglob already has value: it will be averaged
                accel_elastic(1,iglob) = accel_elastic(1,iglob) + extsource(2,it,eleid)
              endif

            enddo
          enddo
        endif
      endif
    enddo

    ! write(IMAIN,*) "test sync", myrank
     !write(IMAIN,*) "test 10", myrank

    ! take an average based on the number of adding acceleration
    ! do  k = 1,iadd
    !   if (P_SV) then
    !     ! write(IMAIN,*) "test 9", myrank
    !      !write(IMAIN,*) "test 11", myrank
    !
    !     !write(IMAIN,*) iglob_add_accel_smooth_idlist(1, k), iglob_add_accel_smooth_idlist(2, k)
    !     accel_elastic(1,iglob_add_accel_smooth_idlist(1, k)) = &
    !     accel_elastic(1,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
    !     accel_elastic(2,iglob_add_accel_smooth_idlist(1, k)) = &
    !     accel_elastic(2,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
    !      !write(IMAIN,*) "test 12", myrank
    !
    !   else
    !     accel_elastic(1,iglob_add_accel_smooth_idlist(1, k)) = &
    !     accel_elastic(1,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
    !   endif
    ! enddo

    !write(IMAIN,*) "test 13", myrank

    deallocate(iphase_elelist)

  enddo ! end do iphase
  !write(IMAIN, *) "myrank :", myrank," total modified ele:", num_modified_ele
  !  write(IMAIN,*) "test 10", myrank

  end subroutine add_ext_source


  subroutine read_ext_source(extsource, extsource_eleid, extsource_globalid)

    ! reads in time series based on external source files
    use constants, only: CUSTOM_REAL, EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX
    use shared_input_parameters, only: DT

    implicit none

    ! local parameters
    integer, intent(in) :: extsource_eleid, extsource_globalid
    double precision, dimension(3,EXT_SOURCE_TRACE_MAX,EXT_SOURCE_NUM_MAX), intent(inout) :: extsource
    integer :: tid, ier
    double precision :: t_temp, ax_temp, az_temp
    character(len=200):: line
    character(len=200):: source_timeseries_name

    ! read external source element
    write(source_timeseries_name, '(A,"/EXT",I0.8,".dat")') './extsource', extsource_globalid
    !write(*,*) source_timeseries_name
    open(unit=32,file=trim(source_timeseries_name),status='old',action='read',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening external source file, please make sure file exists...')

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
        extsource(1, tid, extsource_eleid) = t_temp
        extsource(2, tid, extsource_eleid) = ax_temp
        extsource(3, tid, extsource_eleid) = az_temp
      endif
    enddo

    ! timestep check
    if (abs((extsource(1, 2, extsource_eleid) - extsource(1, 1, extsource_eleid)) - DT) > 1.d-4*DT) then
      call stop_the_code('Time step of external source is diferent with input file. Please check the time step.')
    endif

    !write(*,*) extsource(2, 3000:3010, i)
    close (32)

  end subroutine read_ext_source
  !
  !-----------------------------------------------------------------------------------
  !
  subroutine read_extsource_id(chan, line, ier)
     integer, intent(in):: chan
     integer, intent(inout):: ier
     character(len=200), intent(inout):: line

  10 continue
     read (chan, *, iostat=ier) line
     if (line(1:1) .eq. '#') goto 10
     return
  end subroutine read_extsource_id

  !
  !-----------------------------------------------------------------------------------
  !
  subroutine setup_glob2loctable()
    ! read glob2loctable for acceleration injection at coupling elements
    use constants, only: IMAIN, IIN, OUTPUT_FILES, MAX_STRING_LEN
    use specfem_par, only: myrank, glob2loc_table, nspec
    use shared_input_parameters, only: NPROC

    implicit none

    integer :: i, ier
    integer, dimension(:, :), allocatable :: glob2loc_table_root
    character(len=MAX_STRING_LEN) :: prname

    !if (myrank==0) then
    allocate(glob2loc_table(nspec, 3))

    write(prname,"(a,i5.5,a)") trim(OUTPUT_FILES)//'glob2loc_table',myrank,'.bin'
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Failed to open glob2loc_table.bin with COUPLING_IN = .true. Please check the file.')


    do i = 1, nspec
      write(IMAIN, *) myrank

      ! read(IIN) glob2loc_table_root(i, 1), glob2loc_table_root(i, 2), glob2loc_table_root(i, 3)
      ! call bcast_all_i(glob2loc_table_root(i, 1))
      ! call bcast_all_i(glob2loc_table_root(i, 2))
      ! call bcast_all_i(glob2loc_table_root(i, 3))

      read(IIN) glob2loc_table(i, 1), glob2loc_table(i, 2), glob2loc_table(i, 3)
      !write(IMAIN, *) glob2loc_table(i, 1), glob2loc_table(i, 2), glob2loc_table(i, 3)
    enddo

      !   ! call bcast_all_singlei(gl1)
      !   ! call bcast_all_singlei(gl2)
      !   ! call bcast_all_singlei(gl3)
      !   call MPI_BCAST(gl1,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
      !   call MPI_BCAST(gl1,2,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
      !   call MPI_BCAST(gl1,3,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
      !
      !   !write(IMAIN,*) glob2loc_table(i, 1), glob2loc_table(i, 2), glob2loc_table(i, 3)
      !   write(IMAIN, *) "test 1"
      !   write(IMAIN, *) gl1, gl2, gl3
      !   call send_singlei(dummy, 1, i)
      ! else
      !   write(IMAIN, *) "test 2"
      !   write(IMAIN, *) "myrank is:", myrank
      !   write(IMAIN, *) "NPROC is:", NPROC
      !
      !   call recv_singlei(dummy, myrank-1, i)
      !   if (myrank < NPROC-1) call send_singlei(1, myrank+1, 1)
      ! endif
      !
      ! write(IMAIN, *) "test 3"
      !
      ! glob2loc_table(i, 1) = gl1
      ! glob2loc_table(i, 2) = gl2
      ! glob2loc_table(i, 3) = gl3
      !
      !  write(IMAIN, *) gl1, gl2, gl3

    ! enddo
    ! !call synchronize_all()
    !
    ! do i = 1, nspec_total
    !
    !
    !    write(IMAIN, *) gl1, gl2, gl3
    !   ! glob2loc_table(i, 1) = glob2loc_table_root(i, 1)
    !   ! glob2loc_table(i, 2) = glob2loc_table_root(i, 2)
    !   ! glob2loc_table(i, 3) = glob2loc_table_root(i, 3)
    ! enddo
    close(IIN)

    ! if (myrank == 0) close(IIN)

  end subroutine setup_glob2loctable

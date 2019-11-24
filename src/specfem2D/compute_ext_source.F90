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
    use specfem_par, only: myrank
    use shared_parameters, only: iele, number_of_extsource, NPROC

    implicit none

    ! local parameters
    integer :: i, j
    integer :: ier
    character(len=MAX_STRING_LEN):: line
    integer, dimension(2) :: ibuf_num

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
      write(IMAIN,*)
    endif

    if (myrank == 0) then
      ibuf_num(1) = number_of_extsource
      ibuf_num(2) = 0 ! # avoid error with ifort
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
  use specfem_par, only: nglob, P_SV, it, ibool, NSTEP, deltat, nspec,&
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
  !double precision :: cx, cz
  ! integer, dimension(num_iglob_interface) :: buffer_iglob
  ! double precision, dimension(num_iglob_interface, 2) :: buffer_accel

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

  num_modified_ele = 0

  if (it == 1) then
    !1. setup table for global id to MPI local id
    call setup_glob2loctable()
    !2. setup table for MPI interface igrob
    call setup_iglob_interface()
  endif

  do iphase = 1, 2

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

    do eleid = 1, number_of_extsource

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
              exit
            endif
          enddo
        endif
      enddo

      if (is_ispec_in) then
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,loc_eleid)
              ! find if the acceleration is already added to the iglob
              ! this process is essentially only needed for once in the beginning, but
              ! due to the readablility and not so efficient even if implementing subroutine
              ! for the initialization of iglob_add_accel_smooth_idlist, we find it every time step.
              ! iglob_add_accel_smooth_idlist (1,:) = iglob id
              ! iglob_add_accel_smooth_idlist (2,:) = number of adding acceleration at the iglob
              !write(IMAIN,*) "test 3", myrank
              is_iglob_overlap = .false.

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

              if (.not. is_iglob_overlap) then
                ! this iglob is first time to add

                accel_elastic(1,iglob) = extsource(2,it,eleid)
                accel_elastic(2,iglob) = extsource(3,it,eleid)
                iadd = iadd + 1
                iglob_add_accel_smooth_idlist(1, iadd) = iglob
                iglob_add_accel_smooth_idlist(2, iadd) = 1

              else
                ! this iglob already has value: it will be averaged
                accel_elastic(1,iglob) = accel_elastic(1,iglob) + extsource(2,it,eleid)
                accel_elastic(2,iglob) = accel_elastic(2,iglob) + extsource(3,it,eleid)
              endif ! if the iglob overlap in a single MPI domain
            enddo ! do loop all gll point on an element
          enddo ! do loop all gll point on an element
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

    ! take an average based on the number of adding acceleration
    do  k = 1,iadd
      if (P_SV) then
        accel_elastic(1,iglob_add_accel_smooth_idlist(1, k)) = &
        accel_elastic(1,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
        accel_elastic(2,iglob_add_accel_smooth_idlist(1, k)) = &
        accel_elastic(2,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
      else
        accel_elastic(1,iglob_add_accel_smooth_idlist(1, k)) = &
        accel_elastic(1,iglob_add_accel_smooth_idlist(1, k))/dble(iglob_add_accel_smooth_idlist(2, k))
      endif
    enddo

    deallocate(iphase_elelist)

  enddo ! end do iphase

  ! smooth between MPI domain
#ifdef USE_MPI
  call smooth_MPI_interface(accel_elastic)
#endif

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

    implicit none

    integer :: i, ier
    character(len=MAX_STRING_LEN) :: prname

    allocate(glob2loc_table(nspec, 3))

    write(prname,"(a,i5.5,a)") trim(OUTPUT_FILES)//'glob2loc_table',myrank,'.bin'
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Failed to open glob2loc_table.bin with COUPLING_IN = .true. Please check the file.')


    do i = 1, nspec
      read(IIN) glob2loc_table(i, 1), glob2loc_table(i, 2), glob2loc_table(i, 3)
    enddo

    close(IIN)

  end subroutine setup_glob2loctable
  !
  !-----------------------------------------------------------------------------------
  !
  subroutine smooth_MPI_interface(accel_elastic)
    ! inject the "source" from external source files

    use my_mpi_communicator
    use mpi
    use constants, only: CUSTOM_REAL, IMAIN, NDIM, EXT_SOURCE_NUM_MAX, EXT_SOURCE_TRACE_MAX, &
                        NGLLX,NGLLZ
    use specfem_par, only: myrank, nglob, nspec, P_SV, ibool, coord, nspec_outer_elastic,phase_ispec_inner_elastic, &
                      glob2loc_table, num_iglob_interface, coord_interface, iglob_interface_table

    use shared_parameters, only: iele, number_of_extsource

    implicit none

    ! local parameters
    integer :: i, j, k, iglob, eleid, ier
    integer :: iphase, num_elements
    real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic
    integer, allocatable, dimension(:) :: iphase_elelist
    integer :: ispec_p
    logical :: is_ispec_in
    integer :: loc_rank, loc_eleid
    double precision :: cx, cz
    integer, dimension(num_iglob_interface) :: buffer_iglob
    double precision, dimension(num_iglob_interface, 2) :: buffer_accel
    double precision, dimension(num_iglob_interface, 2) :: sum_accel

    ! find gll points on interface in outer elements

    buffer_iglob(:) = 0
    buffer_accel(:,:) = 0._CUSTOM_REAL
    sum_accel(:,:) = 0._CUSTOM_REAL

    iphase = 1
    num_elements = nspec_outer_elastic
    allocate(iphase_elelist(num_elements))

    ! make list of outer/inner elements
    do ispec_p = 1,num_elements
      iphase_elelist(ispec_p) = phase_ispec_inner_elastic(ispec_p,iphase)
    enddo

    do eleid = 1, number_of_extsource

      is_ispec_in = .false.

      do i = 1, nspec
        if (glob2loc_table(i, 1) == iele(eleid)) then
          loc_rank = glob2loc_table(i, 2)
          loc_eleid = glob2loc_table(i, 3)
          ! check if this local element is inner or outer
          do ispec_p = 1,num_elements
            if (iphase_elelist(ispec_p) == loc_eleid) then
              is_ispec_in = .true.
              exit
            endif
          enddo
        endif
      enddo

      if (is_ispec_in) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,loc_eleid)
            cx    = coord(1, iglob) ! this iglob is local id in each MPI domain
            cz    = coord(2, iglob) ! this iglob is local id in each MPI domain

            do k = 1, num_iglob_interface
              if (coord_interface(k, 1) .eq. cx .and. &
                coord_interface(k, 2) .eq. cz) then
                ! This GLL point is on the MPI interface. Process separately using the table.
                !if (myrank==0) write(*,*) coord_interface(i1, 1) , coord_interface(i1, 2) , cx, cz
                buffer_iglob(k) = iglob
                buffer_accel(k, 1) = accel_elastic(1, iglob)
                buffer_accel(k, 2) = accel_elastic(2, iglob)
                exit
              endif
            enddo
          enddo
        enddo
      endif
    enddo

    ! do i = 1, num_iglob_interface
    !   if (myrank==0) write(IMAIN, *) myrank, ":",i, buffer_iglob(i), buffer_accel(i, 1), buffer_accel(i, 2)
    ! enddo
    ! reduce buffer_accel

    call MPI_ALLREDUCE(buffer_accel,sum_accel,num_iglob_interface*2,MPI_DOUBLE_PRECISION, &
                        MPI_SUM,my_local_mpi_comm_world,ier)


    do i = 1, num_iglob_interface
      !write(IMAIN, *) myrank, ":",i, buffer_iglob(i), sum_accel(i, 1), sum_accel(i, 2)
      if (iglob_interface_table(i, myrank+2) == 1) then
        ! This MPI domain has the shared point on MPI interface.
        if (P_SV) then
          accel_elastic(1, buffer_iglob(i)) = sum_accel(i, 1) / dble(iglob_interface_table(i, 1))
          accel_elastic(2, buffer_iglob(i)) = sum_accel(i, 2) / dble(iglob_interface_table(i, 1))
        else
          accel_elastic(1, buffer_iglob(i)) = sum_accel(i, 1) / dble(iglob_interface_table(i, 1))
        endif
      endif
    enddo

    deallocate(iphase_elelist)

  end subroutine smooth_MPI_interface

  !-----------------------------------------------------------------------------------
  ! The subroutines below are for the smoothing on iglob located at MPI interface.
  ! We assemble a table which indicate the iglob id, number of overlap and processor rank
  !-----------------------------------------------------------------------------------

  subroutine setup_iglob_interface()
    ! read glob2loctable for acceleration injection at coupling elements
    use mpi
    use constants, only: IMAIN, EXT_SOURCE_NUM_MAX, CUSTOM_REAL
    use specfem_par, only: myrank, NPROC
    implicit none

    ! iglob_interface_table_temp contains:
    ! (iglob, num of overlap on MPI interface, contain flag on the MPI Domain)
    integer, dimension(EXT_SOURCE_NUM_MAX*3,NPROC+1) :: iglob_interface_table_temp
    double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2) :: coord_interface_temp
    integer :: num_iglob_temp

    ! evaluate only outer elements in MPI domain

    !initialize array
    iglob_interface_table_temp(:,:) = 0._CUSTOM_REAL
    num_iglob_temp = 0

    if (NPROC > 1) then

      write(IMAIN, *) "start searching iglob on MPI interface"
      if (myrank == 0) then

        call recv_iglob_on_interface(1, iglob_interface_table_temp, coord_interface_temp,num_iglob_temp)
        write(IMAIN, *) "receive at root"

        call find_iglob_on_interface(iglob_interface_table_temp,coord_interface_temp, num_iglob_temp)
        write(IMAIN, *) "num_iglob_temp: ", num_iglob_temp
        write(IMAIN,*) "MPI circulation of find_iglob_on_interface done."

      else
        if (myrank == NPROC-1) then
          ! start finding iglob froom this MPI domain
          write(IMAIN, *) "start ", myrank
          call find_iglob_on_interface(iglob_interface_table_temp,coord_interface_temp, num_iglob_temp)
          ! send iglob_interface_table_temp and num_iglob_temp to next MPI domain
          call send_iglob_on_interface(myrank - 1, iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)
          write(IMAIN, *) "num_iglob_temp: ", num_iglob_temp, "send to", myrank - 1

        else
          call recv_iglob_on_interface(myrank + 1, iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)
          write(IMAIN, *) "receive from", myrank + 1
          call find_iglob_on_interface(iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)
          call send_iglob_on_interface(myrank - 1, iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)
          write(IMAIN, *) "num_iglob_temp: ", num_iglob_temp, "send to", myrank - 1

        endif
      endif

      ! synchronizes processes
      call synchronize_all()

      write(IMAIN, *) "start broad cast ", myrank

      call bcast_iglob_on_interface(iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)

      !-------!
      !Debug
      !call debug_dump_table(iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)
      !-------!

      ! allocate and fill out iglob_interface_table
      call assemble_iglob_interface_table(iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)

      !-------!
      !Debug
      !call debug_assemble_dump_table()
      !-------!

    else
      ! serial run case: in this case num_iglob_temp remains to be 0.
    endif

  end subroutine setup_iglob_interface
  !
  !-----------------------------------------------------------------------------------
  !
  subroutine find_iglob_on_interface(iglob_interface_table_temp, coord_interface_temp, num_iglob_temp)

    ! read glob2loctable for acceleration injection at coupling elements
    use constants, only: IMAIN, EXT_SOURCE_NUM_MAX, NGLLX,NGLLZ
    use specfem_par, only: myrank, nspec, NPROC, glob2loc_table, nspec_outer_elastic, &
                          ibool, phase_ispec_inner_elastic, coord
    use shared_parameters, only: iele, number_of_extsource

    implicit none

    integer :: i, i1, j1, k, k1, eleid
    integer, dimension(EXT_SOURCE_NUM_MAX * 3, NPROC + 1), intent(inout) :: iglob_interface_table_temp
    double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2), intent(inout)  :: coord_interface_temp
    integer, intent(inout) :: num_iglob_temp
    integer, dimension(:), allocatable :: iphase_elelist
    double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2) :: MPIdomain_duplication_list
    integer :: ispec_p, id_iglob_overlap, iglob, num_elements, num_MPIdomain_duplication
    integer :: loc_rank, loc_eleid
    double precision :: cx, cz
    logical :: is_MPIdomain_duplication


    ! choses outer elements
    num_elements = nspec_outer_elastic

    allocate(iphase_elelist(num_elements))

    ! make list of outer elements
    do ispec_p = 1,num_elements
      iphase_elelist(ispec_p) = phase_ispec_inner_elastic(ispec_p,1)
    enddo

    ! Debug
    !write(IMAIN, *) myrank, ": coord size:" ,size(coord,2)
    do ispec_p = 1,num_elements
      do j1 = 1,NGLLZ
        do i1 = 1,NGLLX
          iglob = ibool(i1,j1,iphase_elelist(ispec_p))
          !write(IMAIN, *) coord(1,iglob), coord(2,iglob)
        enddo
      enddo
    enddo

    num_MPIdomain_duplication = 0
    do eleid = 1, number_of_extsource
      do i = 1, nspec
        if (glob2loc_table(i, 1) == iele(eleid)) then
          loc_rank = glob2loc_table(i, 2)
          loc_eleid = glob2loc_table(i, 3)
          ! do search on outer elements
          do ispec_p = 1,num_elements
            ! check if this  element is outer
            if (iphase_elelist(ispec_p) == loc_eleid) then
              ! loop GLL points in coupling element
              do j1 = 1,NGLLZ
                do i1 = 1,NGLLX
                  iglob = ibool(i1,j1,loc_eleid)
                  cx = coord(1,iglob)
                  cz = coord(2,iglob)
                  id_iglob_overlap = 0
                  ! search iglob in the present domain

                  is_MPIdomain_duplication = .false.
                  do k1 = 1, num_MPIdomain_duplication
                    if (MPIdomain_duplication_list(k1, 1) .eq. cx .and. &
                      MPIdomain_duplication_list(k1, 2) .eq. cz) then
                      id_iglob_overlap = -1
                      is_MPIdomain_duplication = .true.
                    endif
                  enddo

                  if (.not. is_MPIdomain_duplication) then
                    do k = 1, num_iglob_temp
                      !write(IMAIN,*) "test 2.1", myrank, iele(eleid)
                      if (coord_interface_temp(k, 1) .eq. cx .and. &
                        coord_interface_temp(k, 2) .eq. cz) then
                        ! this gll point is overlapped on MPI interfacec
                        id_iglob_overlap = k
                        exit
                      endif
                    enddo
                  endif

                  ! modify coord_interface_temp and iglob_interface_table_temp
                  if (id_iglob_overlap == -1) then
                    ! write(IMAIN,*) "test avoid duplication"
                    ! write(IMAIN, *) cx, cz

                    continue
                  elseif (id_iglob_overlap == 0) then
                    ! this iglob is first time to add
                    !write(IMAIN, *) "find id_iglob_overlap"
                    num_iglob_temp = num_iglob_temp + 1
                    ! coordinates
                    coord_interface_temp(num_iglob_temp, 1) = cx
                    coord_interface_temp(num_iglob_temp, 2) = cz
                    ! num of overlap
                    iglob_interface_table_temp(num_iglob_temp, 1) = &
                        iglob_interface_table_temp(num_iglob_temp, 1) + 1
                    ! MPI domain flag
                    iglob_interface_table_temp(num_iglob_temp, myrank+2) = 1
                    ! save iglob to avoid counting overlap within this MPIDomain
                    num_MPIdomain_duplication = num_MPIdomain_duplication + 1
                    MPIdomain_duplication_list(num_MPIdomain_duplication, 1) = cx
                    MPIdomain_duplication_list(num_MPIdomain_duplication, 2) = cz

                  else
                    ! this iglob already has value: it will be averaged
                    iglob_interface_table_temp(id_iglob_overlap, 1) = &
                        iglob_interface_table_temp(id_iglob_overlap, 1) + 1
                    iglob_interface_table_temp(id_iglob_overlap, myrank+2) = 1

                    ! avoid duplication in this MPIDomain
                    num_MPIdomain_duplication = num_MPIdomain_duplication + 1
                    MPIdomain_duplication_list(num_MPIdomain_duplication, 1) = cx
                    MPIdomain_duplication_list(num_MPIdomain_duplication, 2) = cz

                  endif !id_iglob_overlap
                enddo ! i1 NGLLX
              enddo ! j1 NGLLZ
            endif ! if outer
          enddo ! do search within iphase_elelist

          exit ! exit nspec loop because this eleid is already found and processed

        endif ! if this iele is in MPIDomain
      enddo ! do search all elements within MPIDomain
    enddo ! loop all coupling element

    deallocate(iphase_elelist)

    end subroutine find_iglob_on_interface

    !
    !-----------------------------------------------------------------------------------
    !
    subroutine send_iglob_on_interface(dest, iglob_interface_table_temp, coord_interface_temp, num_iglob_temp)

      use constants, only: EXT_SOURCE_NUM_MAX
      use specfem_par, only: NPROC

      ! send iglob_interface_table_temp and num_iglob_temp to dest
      use my_mpi_communicator
      use mpi

      implicit none

      integer :: dest
      integer, dimension(2) :: sendbuf
      integer :: ier

      integer, dimension(EXT_SOURCE_NUM_MAX*3,NPROC+1), intent(in) :: iglob_interface_table_temp
      double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2), intent(in)  :: coord_interface_temp

      integer, intent(in) :: num_iglob_temp

      !send iglob_interface_table_temp
      call MPI_SEND(iglob_interface_table_temp,(EXT_SOURCE_NUM_MAX*3)*(NPROC+1), &
                  MPI_INTEGER,dest,1,my_local_mpi_comm_world,ier)

      !send coord_interface_temp
      call MPI_SEND(coord_interface_temp,(EXT_SOURCE_NUM_MAX*3)*2, &
                  MPI_DOUBLE_PRECISION,dest,1,my_local_mpi_comm_world,ier)

      !send num_iglob_temp
      sendbuf(1) = num_iglob_temp
      sendbuf(2) = 0
      call MPI_SEND(sendbuf,1, MPI_INTEGER,dest,2,my_local_mpi_comm_world,ier)

    end subroutine send_iglob_on_interface
    !
    !-----------------------------------------------------------------------------------
    !
    subroutine recv_iglob_on_interface(dest, iglob_interface_table_temp, coord_interface_temp, num_iglob_temp)

      ! synchronuous/blocking receive
      use constants, only: EXT_SOURCE_NUM_MAX
      use specfem_par, only: NPROC

      use my_mpi_communicator
      use mpi

      implicit none

      integer :: dest
      integer, dimension(2) :: recvbuf
      integer :: ier

      integer, dimension(EXT_SOURCE_NUM_MAX*3,NPROC+1), intent(inout) :: iglob_interface_table_temp
      double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2), intent(inout)  :: coord_interface_temp

      integer, intent(inout) :: num_iglob_temp

      !receive iglob_interface_table_temp
      call MPI_RECV(iglob_interface_table_temp,(EXT_SOURCE_NUM_MAX*3)*(NPROC+1), &
                  MPI_INTEGER,dest,1,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

      !receive coord_interface_temp
      call MPI_RECV(coord_interface_temp,(EXT_SOURCE_NUM_MAX*3)*2, &
                  MPI_DOUBLE_PRECISION,dest,1,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

      !receive num_iglob_temp
      call MPI_RECV(recvbuf, 1, MPI_INTEGER,dest,2,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

      num_iglob_temp = recvbuf(1)

    end subroutine recv_iglob_on_interface
    !
    !-----------------------------------------------------------------------------------
    !
    subroutine bcast_iglob_on_interface(iglob_interface_table_temp, coord_interface_temp, num_iglob_temp)

      ! synchronuous/blocking receive
      use constants, only: EXT_SOURCE_NUM_MAX
      use specfem_par, only: NPROC

      use my_mpi_communicator
      use mpi

      implicit none

      integer, dimension(2) :: bcastbuf
      integer :: ier

      integer, dimension(EXT_SOURCE_NUM_MAX*3,NPROC+1), intent(inout) :: iglob_interface_table_temp
      double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2), intent(inout)  :: coord_interface_temp
      integer, intent(inout) :: num_iglob_temp

      call MPI_BCAST(iglob_interface_table_temp,(EXT_SOURCE_NUM_MAX*3)*(NPROC+1), &
                    MPI_INTEGER,0,my_local_mpi_comm_world,ier)

      call MPI_BCAST(coord_interface_temp,(EXT_SOURCE_NUM_MAX*3)*2, &
                    MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

      bcastbuf(1) = num_iglob_temp
      bcastbuf(2) = 0

      call MPI_BCAST(bcastbuf,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

      num_iglob_temp = bcastbuf(1)

    end subroutine bcast_iglob_on_interface
    !
    !-----------------------------------------------------------------------------------
    !

    subroutine assemble_iglob_interface_table(iglob_interface_table_temp, coord_interface_temp, num_iglob_temp)

      ! read glob2loctable for acceleration injection at coupling elements
      use constants, only: IMAIN, EXT_SOURCE_NUM_MAX, CUSTOM_REAL
      use specfem_par, only: NPROC, iglob_interface_table, coord_interface, num_iglob_interface

      implicit none

      integer :: i, j, numoverlap
      integer, dimension(EXT_SOURCE_NUM_MAX*3,NPROC+1), intent(in) :: iglob_interface_table_temp
      double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2), intent(in)  :: coord_interface_temp

      integer, intent(in) :: num_iglob_temp

      ! count the number of overlapped iglob at MPI interface
      numoverlap = 0
      do i = 1, num_iglob_temp
        if (iglob_interface_table_temp(i, 1) > 1) then
          ! this iglob has overlap between MPI domains
          numoverlap = numoverlap + 1
        endif
      enddo

      ! initializing iglob_interface_table
      allocate(iglob_interface_table(numoverlap, NPROC+1))
      allocate(coord_interface(numoverlap, 2))

      num_iglob_interface = numoverlap

      ! reassemble iglob_interface_table from iglob_interface_table_temp on every processor
      numoverlap = 0
      do i = 1, num_iglob_temp
        if (iglob_interface_table_temp(i, 1) > 1) then
          !write(*,*) "test find overlap:", iglob_interface_table_temp(i, 1)
          numoverlap = numoverlap + 1
          ! this iglob has overlap between MPI domains
          do j = 1,NPROC+1
            !write(*,*) "test find numoverlap:", numoverlap
            iglob_interface_table(numoverlap, j) = iglob_interface_table_temp(i, j)
            !write(IMAIN,*)  "test find numoverlap in table:", iglob_interface_table(numoverlap, j)
          enddo
          !write(IMAIN,*) iglob_interface_table(numoverlap,1:NPROC+1)
          coord_interface(numoverlap, 1) = coord_interface_temp(i, 1)
          coord_interface(numoverlap, 2) = coord_interface_temp(i, 2)

          ! write(IMAIN,*) myrank, "test find coord1:", coord_interface_temp(i, 1), coord_interface_temp(i, 2)
          ! write(IMAIN,*) myrank, "test find coord2:", coord_interface(numoverlap, 1), coord_interface(numoverlap, 2)

        endif
      enddo

      !write(IMAIN,*) "myrank: ", myrank
      !write(IMAIN,*) coord_interface
    end subroutine assemble_iglob_interface_table
    !
    !-----------------------------------------------------------------------------------
    !
    subroutine debug_dump_table(iglob_interface_table_temp,coord_interface_temp,num_iglob_temp)

      use constants, only: IMAIN, OUTPUT_FILES, EXT_SOURCE_NUM_MAX
      use specfem_par, only: myrank, NPROC

      implicit none

      integer, dimension(EXT_SOURCE_NUM_MAX*3,NPROC+1), intent(in) :: iglob_interface_table_temp
      double precision, dimension(EXT_SOURCE_NUM_MAX*3, 2), intent(in)  :: coord_interface_temp

      integer, intent(in) :: num_iglob_temp
      integer :: i, ier
      double precision :: cx, cz
      character(len=200) :: foname

      write(foname, "(A23,I0.3,A4)") "debug_iglob_interface_table_temp_", myrank,".csv"
      !write(IMAIN, *) foname
      open(unit=1002,file=trim(OUTPUT_FILES)//foname,status='unknown',iostat=ier)
      if (ier /= 0 ) then
        call stop_the_code('Error saving debug iglob_interface_table_temp.csv')
      endif

      do i = 1, num_iglob_temp

        cx = coord_interface_temp(i, 1)
        cz = coord_interface_temp(i, 2)

        write(1002, '(1(I5,","), 2(F20.8,","), 4(I5,","))') iglob_interface_table_temp(i, 1), &
                    cx, cz, iglob_interface_table_temp(i, 2:NPROC+1)
      enddo

      close(1002)

    end subroutine debug_dump_table
    !
    !-----------------------------------------------------------------------------------
    !
    subroutine debug_assemble_dump_table()
      ! read glob2loctable for acceleration injection at coupling elements
      use constants, only: OUTPUT_FILES
      use specfem_par, only: myrank, NPROC, iglob_interface_table, coord_interface, num_iglob_interface

      implicit none

      integer :: i, ier
      double precision :: cx, cz
      character(len=200) :: foname

      write(*,*) myrank, "test"
      write(foname, "(A22,I0.3,A4)") "debug_assembled_table_", myrank,".csv"
      open(unit=1003,file=trim(OUTPUT_FILES)//foname,status='unknown',iostat=ier)
      if (ier /= 0 ) then
        call stop_the_code('Error saving debug iglob_interface_table_temp.csv')
      endif

      do i = 1, num_iglob_interface

        cx = coord_interface(i, 1)
        cz = coord_interface(i, 2)
        write(*,*) myrank, "test", cx, cz

        write(1003, '(1(I5,","), 2(F20.8,","), 4(I5,","))') iglob_interface_table(i, 1), &
                    cx, cz, iglob_interface_table(i, 2:NPROC+1)
      enddo

      close(1003)

    end subroutine debug_assemble_dump_table
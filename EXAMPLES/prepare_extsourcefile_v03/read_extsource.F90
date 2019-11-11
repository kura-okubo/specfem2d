!scratch to test reading external source
  program read_extsource

  implicit none

  ! local parameters

  double precision :: DT
  double precision :: t_temp, ax_temp, az_temp
  character*200:: line
  character*200:: extsource_filename, source_timeseries_name
  integer :: ier, NumSource
  integer, dimension(100) :: iele
  integer :: i, j, tid
  double precision, dimension(100, 3, 100000) :: extsource


  DT = 0.1
  write(*,*) "test", DT

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

  ! external source sanity check


  !
  !
  ! ! counts number of lines
  ! icounter = 0
  ! do while(ier == 0)
  !   read(read,"(a)",iostat=ier) string_read
  !
  !   if (ier == 0) then
  !     ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
  !     if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)
  !
  !     ! suppress leading and trailing white spaces, if any
  !     string_read = adjustl(string_read)
  !     string_read = string_read(1:len_trim(string_read))
  !
  !     ! if the line is not empty and is not a comment, count it
  !     if (len_trim(string_read) > 0 .and. (index(string_read,'#') == 0 .or. index(string_read,'#') > 1)) icounter = icounter + 1
  !   endif
  ! enddo
  ! close(31)


  end program read_extsource

  subroutine read_extsource_id(chan, line, ier)
     integer, intent(in):: chan
     integer, intent(inout):: ier
     character*200, intent(inout):: line

  10 continue
     read (chan, *, iostat=ier) line
     if (line(1:1) .eq. '#') goto 10
     return
  end subroutine read_extsource_id

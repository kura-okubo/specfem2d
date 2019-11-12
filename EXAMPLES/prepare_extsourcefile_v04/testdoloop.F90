!scratch to test reading external source
  program read_extsource

  implicit none
  integer :: i

  do i = 1,3
    write(*,*) dble(i)
  enddo

  end program read_extsource

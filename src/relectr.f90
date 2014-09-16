subroutine relectr(iounit,filename)
! read electrode node numbers from file

! Andreas Kemna                                            11-Oct-1993
! last change                                              24-Oct-1996

! ..............................................................................

  USE electrmod
  USE elemmod
  USE errmod

  IMPLICIT none

! ..............................................................................
! IO unit
  INTEGER           iounit
! filename
  CHARACTER (80)    filename
! ..............................................................................
! internals
! index vars
  INTEGER           i,ifp
! ..............................................................................
! open file
  fetxt = filename
  errnr = 1
  open(iounit,file=TRIM(fetxt),status='old',err=999)
  errnr = 3
! read number of electrodes
  read(iounit,*,end=1001,err=1000) eanz
  ALLOCATE (enr(eanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation enr'
     errnr = 94
     goto 1000
  END IF
! read electrode node numbers
  do i=1,eanz
     read(iounit,*,end=1001,err=1000) enr(i)
     WRITE (ifp,*)i,sx(snr(enr(i))),sy(snr(enr(i)))
     if (enr(i).gt.sanz) then
        fetxt = ' '
        errnr = 29
        goto 1000
     end if
  end do

! write electrode positions to file
  CALL get_unit(ifp)
  OPEN (ifp,FILE='inv.elecpositions',STATUS='replace')
  WRITE (ifp,*)eanz
  do i=1,eanz
     WRITE (ifp,*)i,sx(snr(enr(i))),sy(snr(enr(i)))
  end do

  write(*,'(a,I7,a)') ' read',eanz,' electrode positions'
  CLOSE (ifp)
  close(iounit)
  errnr = 0
  return
! ..............................................................................
! error messages
999 return
1000 close(iounit)
  return
1001 close(iounit)
  errnr = 2
  return
end subroutine relectr

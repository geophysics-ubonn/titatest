subroutine rdatm(iounit,filename)
! read electrode configurations and current strength from file
!
! Andreas Kemna                                            11-Oct-1993
! last change   22-Feb-2006
! ..............................................................................

  USE datmod
  USE electrmod
  USE errmod

  IMPLICIT none

! ..............................................................................
! IO unit
  INTEGER           iounit
! filename
  CHARACTER (80)    filename
! ..............................................................................
! internals
! index var
  INTEGER           i
! electrode numbers
  INTEGER           elec1, elec2, elec3, elec4
! crtomo format compliance
  LOGICAL           ::    crtf
! ..............................................................................
! open file
  fetxt = filename
  errnr = 1
  open( iounit, file = TRIM( fetxt ), status = 'old', err = 999 )
  errnr = 3

! read number of configurations
  read( iounit, *, end = 1001, err = 1000 ) nanz
! check for crtomo format compliance
  read( iounit, *, end = 1001, err = 1000 ) elec1
  BACKSPACE( iounit )
  elec3 = elec1 - 10000
  crtf = ( elec3 .ge. 0 ) ! crtomo compliance

  ALLOCATE (strnr(nanz),strom(nanz),volt(nanz),sigmaa(nanz),&
       kfak(nanz),vnr(nanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'error memory allocation volt '
     errnr = 94
     goto 1000
  END IF

! read current and voltage electrodes
  do i=1,nanz
     IF ( crtf ) THEN
        read( iounit, *, end = 1001, err = 1000 ) strnr(i), vnr(i)
     ELSE
        read( iounit, *, end = 1001, err = 1000 ) elec1, elec2, elec3, elec4
        strnr(i) = elec1*10000 + elec2
        vnr(i)   = elec3*10000 + elec4
     END IF
! override current strength with 1 A
     strom(i) = 1d0
! compute current electrodes
     elec1 = mod(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000
! compute voltage electrodes
     elec3 = mod(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000
! errors
     if (elec1.lt.0.or.elec1.gt.eanz.or. &
          elec2.lt.0.or.elec2.gt.eanz.or. &
          elec3.lt.0.or.elec3.gt.eanz.or. &
          elec4.lt.0.or.elec4.gt.eanz) then
        WRITE (fetxt,'(a,I5,a)')'electrode pair ',i,'not correct '
        errnr = 46
        goto 1000
     end if
! >> RM
! plausibility check of possible electrode intersection
! devide the strnr and vnr into elec{1,2,3,4}
     IF ((elec1.eq.elec2).OR.(elec3.eq.elec4).OR.&
          &((((elec1.eq.elec3).or.(elec1.eq.elec4)).and.(elec1.ne.0)).or.&
          (((elec2.eq.elec3).or.(elec2.eq.elec4)).and.(elec2.ne.0)))) THEN
        WRITE (fetxt,'(a,I7)')' duplicate electrodes for reading ',i
        errnr = 73
        GOTO 1000
     END IF
! << RM
  end do
write(*,'(a,I7,a)') ' read',nanz,' configurations'
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

end subroutine rdatm

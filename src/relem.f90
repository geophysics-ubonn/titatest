SUBROUTINE relem(iounit,filename)
! read grid nodes and element edges/nodes from file 
!
! Andreas Kemna                                            11-Oct-1993
! last change                                              24-Oct-1996
! ..............................................................................

  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT NONE

! ..............................................................................
! IO unit
  INTEGER           iounit
! filename
  CHARACTER (80)    filename
! ..............................................................................
! internal
! index vars
  INTEGER           i, j, k
! dummy vars
  INTEGER           idum, ifln, iflnr
  LOGICAL           my_check, failed(2)
! border to element connection (rnr) dummy vars
  INTEGER           ik1, ik2, jk1, jk2, ic, l
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: my_nrel
! ..............................................................................
! open file
  fetxt = filename
  errnr = 1
  OPEN(iounit,file=TRIM(fetxt),status='old',err=999)
  errnr = 3

! read number of nodes, number of element types and matrix bandwidth
  READ(iounit,*,END=1001,err=1000) sanz,typanz,mb

  ALLOCATE (sx(sanz),sy(sanz),snr(sanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation sx failed'
     errnr = 94
     GOTO 999
  END IF

  ALLOCATE (typ(typanz),nelanz(typanz),selanz(typanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation selanz failed'
     errnr = 94
     GOTO 999
  END IF

! read element types, number of elements of each type and number of nodes of 
! each element type
  READ(iounit,*,END=1001,err=1000)(typ(i),nelanz(i),selanz(i),i=1,typanz)

! set maximum number of nodes for all elements
  smaxs = MAXVAL(selanz)

!     Anzahl der Elemente (ohne Randelemente) und Anzahl der Randelemente
!     bestimmen
! get number of area and boundary elements 
  elanz  = 0
  relanz = 0
  my_check = .FALSE.

  DO i=1,typanz
     IF (typ(i).GT.10) THEN
        relanz = relanz + nelanz(i)
     ELSE
        elanz  = elanz  + nelanz(i)
     END IF
  END DO

  ALLOCATE (nrel(elanz+relanz,smaxs),my_nrel(elanz+relanz,smaxs),rnr(relanz),&
       stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'error memory allocation nrel,rnr failed'
     errnr = 94
     GOTO 999
  END IF

  ALLOCATE (espx(elanz),espy(elanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'error memory allocation espx failed'
     errnr = 94
     GOTO 999
  END IF
! area element midpoints
  espx = 0.;espy = 0.

! read nodes
  READ(iounit,*,END=1001,err=1000) (snr(i),sx(i),sy(i),i=1,sanz)
  write(*,'(a,I7,a)') ' read',sanz,' nodes'

! read area element node numbers
  idum = 0 
  ifln = 0
  iflnr = 0
  DO i=1,typanz
     DO j=1,nelanz(i)
        READ(iounit,*,END=1001,err=1000)(nrel(idum+j,k),k=1,selanz(i))
        IF (typ(i) < 10) THEN ! set area element midpoints
           ifln = ifln + 1
           DO k = 1,selanz(i)
              espx(ifln) = espx(ifln) + sx(snr(nrel(idum+j,k)))
              espy(ifln) = espy(ifln) + sy(snr(nrel(idum+j,k)))
           END DO
           espx(ifln) = espx(ifln) / selanz(i)
           espy(ifln) = espy(ifln) / selanz(i)
        END IF
     END DO
  write(*,'(a,I7,a,I2,a)') ' read',nelanz(i),' type ',typ(i),' elements'
  idum = idum + nelanz(i)

  END DO
! read boundary element node numbers
  READ(iounit,*,END=1001,err=1000) (rnr(i),i=1,relanz)

  CLOSE(iounit)

! ..............................................................................
! RM: write new grid if one of the checks failed
  IF (ANY(failed)) THEN

101  FORMAT (I9)
102  FORMAT (2(I9,1X))
103  FORMAT (3(I9,1X))
104  FORMAT (4(I9,1X))
105  FORMAT (I9,2(1X,G12.4))

     fetxt = TRIM(filename)//'_new'
     PRINT*,'+++ WRITING OUT IMPROVED GRID --> Writing',TRIM(fetxt)
     errnr = 1
     OPEN(iounit,file=TRIM(fetxt),status='replace')
     
     errnr = 3
! write elem.dat header
     WRITE(iounit,103) sanz,typanz,mb
     WRITE(iounit,103) (typ(i),nelanz(i),selanz(i),i=1,typanz)
! write nodes
     DO i=1,sanz
        WRITE(iounit,105) snr(i),sx(i),sy(i)
     END DO
! write element node numbers
     idum = 0;ifln = 0;iflnr = 0
     DO i=1,typanz
        DO j=1,nelanz(i)
           IF (typ(i) == 8)  THEN
              WRITE(iounit,104)(my_nrel(idum+j,k),k=1,selanz(i))
           ELSE IF (typ(i) == 4)  THEN
              WRITE(iounit,103)(my_nrel(idum+j,k),k=1,selanz(i))
           ELSE
              WRITE(iounit,102)(my_nrel(idum+j,k),k=1,selanz(i))
           END IF
        END DO
        idum = idum + nelanz(i)
     END DO
     WRITE(iounit,101) (rnr(i),i=1,relanz)
     CLOSE (iounit)
  END IF
  DEALLOCATE (my_nrel)
! << RM
  errnr = 0
  RETURN

! ..............................................................................
! error messages

999 RETURN

1000 CLOSE(iounit)
  RETURN

1001 CLOSE(iounit)
  errnr = 2
  RETURN

END SUBROUTINE relem

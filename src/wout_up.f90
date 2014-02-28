SUBROUTINE wout_up(kanal,it,itr,switch)

!!!$     Unterprogramm zum Schreiben der Widerstandsverteilung und der
!!!$     modellierten Daten inkl. Elektrodenkennungen.

!!!$     Andreas Kemna                                            28-Sep-1994
!!!$     Letzte Aenderung   10-Mar-2007

!!!$.....................................................................
use alloci, only: prec
  USE datmod, ONLY : nanz, strnr, vnr, sigmaa
  USE invmod, ONLY : m0, wdfak
  USE sigmamod, ONLY : sigma
  USE modelmod, ONLY : mnr
  USE elemmod, ONLY : elanz, espx, espy 
  USE errmod, ONLY : errnr, fetxt
  USE konvmod, ONLY: pharms, betrms, lverb, nrmsd, lam
  USE pathmod, ONLY : ramd, lnramd, slash, mkdir, rmdir
  USE femmod, ONLY : ldc
  IMPLICIT NONE

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND=4),INTENT (IN) ::  kanal
!!!$     Iterationsnummer
  INTEGER (KIND=4),INTENT (IN) ::  it,itr
!!!$     Switching between parameter or data model output
  LOGICAL ,INTENT (IN) ::  switch

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:
!!!$     Indexvariablen
  INTEGER (KIND=4) ::  i,i_f
  INTEGER          :: ifp
!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::  idum,idum2
  CHARACTER (256)   ::  file,itdir
  COMPLEX(prec)  ::  dum

!!!$     diff+<
  REAL(prec)  ::   dum2,dum3
!!!$     diff+>
  CHARACTER (12)   ::   c_i ! iteration number as string for directory like
!!!$ IT_<number> with number in two digits (compatible with old CRTomo versions..)
  CHARACTER(2)    ::   ci ! iteration number as string
  CHARACTER(3)    ::   cu ! model update number as string
!!$  CHARACTER (6)   ::   c_i ! iteration number as string for directory like
!!$!!!$ IT_<number> with number in three digits..
!!$  CHARACTER(3)    ::   ci ! iteration number as string

!!!$     Dateinamen
  CHARACTER (80)   ::  foutfn

  LOGICAL :: crtf
!!!$.....................................................................

!!$ formatstrings


  foutfn = 'update'
  WRITE (ci,'(I2.2)')it
  WRITE (cu,'(I3.3)')itr

  c_i = 'IT_'//ci//'_UP_'//cu

  itdir = TRIM(ramd(1:lnramd))//slash//c_i

  IF (switch) THEN

     PRINT*,'overwrite existing '//TRIM(itdir)
     CALL get_unit(ifp)
     !  print*,TRIM(itdir)
!!$! workaround for compability issues with ifort..
     file = TRIM(itdir)//slash//'tmp.check'
     !  PRINT*,TRIM(file)
     crtf = .FALSE.
!!$ test if you can open a file in the directory..
     OPEN (ifp,FILE=TRIM(file),STATUS='replace',ERR=97)
!!$ if you can, you can, the directory exits and you can remove it safely
     CLOSE(ifp,STATUS='delete')
!!$ set this switch to circumvent mkdir
!!$  PRINT*,'Inversion directory exists, removing content '
!!$ TODO: remove content of invdir?
     CALL SYSTEM (rmdir//TRIM(itdir))
     crtf = .TRUE.
97   CONTINUE
!!$  PRINT*,'Creating inversion directory '//TRIM(itdir)
     CALL SYSTEM (mkdir//TRIM(itdir))

!!!$ write data to iteration directory
!!$ extracting the 'rho' from foutfn path 
     file = TRIM(itdir)//slash//TRIM(foutfn)//'.mag'
     !  print*,'Trying to write '//TRIM(file)
     OPEN (kanal,FILE=TRIM(file),STATUS='replace',ERR=1000)
     WRITE (kanal,*,err=1000) elanz,betrms,lam
     WRITE (kanal,'(3(G12.4,2x))',err=1000) (REAL(espx(i)),REAL(espy(i)),&
          REAL(LOG10(ABS(1d0/sigma(i)))),i=1,elanz)
     CLOSE (kanal)

     !  print*,TRIM(foutfn),TRIM(dvolt)
     file = TRIM(itdir)//slash//TRIM(foutfn)//'.pha'
     !  print*,'Trying to write '//TRIM(file)
     errnr = 1
     OPEN (kanal,FILE=TRIM(file),STATUS='replace',ERR=1000)
     errnr = 4
     WRITE (kanal,*,err=1000) elanz,pharms,lam
     WRITE (kanal,'(3(G12.4,2x))',err=1000)(REAL(espx(i)),REAL(espy(i)),&
          REAL(1d3*ATAN2(AIMAG(1./sigma(i)),REAL(1./sigma(i)))),i=1,elanz)

     CLOSE (kanal)

     file = TRIM(itdir)//slash//TRIM(foutfn)//'.modl'
     !  print*,'Trying to write '//TRIM(file)
     errnr = 1
     OPEN (kanal,FILE=TRIM(file),STATUS='replace',ERR=1000)
     errnr = 4
     WRITE(kanal,*,err=1000) elanz,nrmsd,lam
     WRITE (kanal,'(2(G12.4,2x))',err=1000)(1./REAL(sigma(i)),&
          REAL(1d3*ATAN2(AIMAG(1./sigma(i)),REAL(1./sigma(i)))),&
          i=1,elanz)

     CLOSE (kanal)

  ELSE
     file = TRIM(itdir)//slash//TRIM(foutfn)//'.dat'
     PRINT*,'Trying to write '//TRIM(file)
     errnr = 1
     OPEN (kanal,FILE=TRIM(file),STATUS='replace',ERR=1000)
     errnr = 4
     WRITE(kanal,*,err=1000) nanz

!!!     Stromelektrodennummern, Spannungselektrodennummern und scheinbare
     !!     Widerstandswerte (Betrag und Phase (in mrad)) schreiben
     IF (ldc) THEN
        DO i=1,nanz
           WRITE(kanal,*,err=1000)strnr(i),vnr(i),&
                !         real(1d0/EXP(REAL(sigmaa(i))))
!!$c     diff+<
                REAL(1d0/EXP(REAL(sigmaa(i)))),wdfak(i)
!!$c     diff+>
        END DO
     ELSE
        DO i=1,nanz
           WRITE(kanal,*,err=1000)strnr(i),vnr(i),&
                REAL(1d0/EXP(REAL(sigmaa(i)))),&
!!$c     diff-     1                      real(-1d3*aimag(sigmaa(i)))
!!$c     diff+<
                REAL(-1d3*aimag(sigmaa(i))),wdfak(i)
!!$c     diff+>
        END DO
     END IF

     CLOSE(kanal)
  END IF
  errnr = 0
  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 RETURN

1000 CLOSE(kanal)
  RETURN

END SUBROUTINE wout_up

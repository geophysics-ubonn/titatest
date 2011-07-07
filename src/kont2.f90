subroutine kont2(lsetup)

!!!$     Unterprogramm zur Ausgabe der Kontrollvariablen.

!!!$     Andreas Kemna                                         16-Apr-1996
!!!$     Letzte Aenderung   11-Nov-1997

!!!$.....................................................................

  USE datmod
!!!$     mw
  USE invmod
  USE cjgmod,ONLY:cgres,ncg
  USE errmod
  USE konvmod
  USE pathmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Hilfsschalter
  LOGICAL  ::  lsetup

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND=4)  ::  i,k
!!!$     Platzhalter
  INTEGER,PARAMETER :: ncdump=110
  CHARACTER(ncdump) :: cdump
!!!$.....................................................................
100 FORMAT (t1,a3, t5,i3,t11,g10.4,t69,g10.4,t81,g10.4,t93,i4,t105,g9.3)
101 FORMAT (t1,a3,t5,i3,t11,g10.4,t69,g10.4,t81,g10.4,t93,i4)

110 FORMAT (t1,a3,t5,i3,t11,g10.4,t23,g10.4,t34,g10.4,t46,g10.4,t58,&
       i4,t69,g10.4,t81,g10.4,t93,i4,t105,g9.3,t117,f5.3)
111 FORMAT (t1,a3,t5,i3,t11,g10.4,t23,g10.4,t34,g10.4,t46,g10.4,t58,&
       i4,t69,g10.4,t81,g10.4,t93,i4,t105,f5.3)

105 FORMAT (t1,a3,t5,i3,t11,g10.4,t23,g9.3,t34,g10.4,t46,g10.4,t58,&
       i4,t105,f5.3)

!!!$     'inv.ctr' oeffnen
  DO i=1,ncdump-1
     cdump(i:i+1)='*'
  END DO
  errnr = 1
  fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
  open(fpinv,file=TRIM(fetxt),status='old',err=1000,position='append')
  fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
  open(fpcjg,file=TRIM(fetxt),status='old',err=1000,position='append')
  errnr = 4
  ncg = INT(cgres(1))
!!!$     Erste Iteration (?)
  if (lsetup) then
     IF (lip) THEN          ! FPI?
        write(fpinv,'(a)',err=1000)cdump
!!!$     Robuste Inversion
        if (lrobust) then
           write(fpinv,100,err=1000)'PIT',it,nrmsd,betrms,pharms,npol,l1rat
        else
           write(fpinv,101,err=1000)'PIT',it,nrmsd,betrms,pharms,npol
        end if
        write(fpinv,'(a)',err=1000)cdump
     ELSE
        write(fpinv,'(a)',err=1000)cdump
!!!$     Robuste Inversion
        if (lrobust) then
           write(fpinv,100,err=1000)'IT',it,nrmsd,betrms,pharms,npol,l1rat
        else
           write(fpinv,101,err=1000)'IT',it,nrmsd,betrms,pharms,npol
        end if
        write(fpinv,'(a)',err=1000)cdump
     END IF
  else
!!!$     
!!!$     Hauptiterationen
     if (llam.and..not.lstep) then
        write(fpinv,'(a)',err=1000)cdump
        if (lip) then
           if (lrobust) then
              write(fpinv,110,err=1000)'PIT',it,nrmsd,bdpar,lam,rough,&
                   ncg,betrms,pharms,npol,l1rat,step
           else
              write(fpinv,111,err=1000)'PIT',it,nrmsd,bdpar,lam,rough,&
                   ncg,betrms,pharms,npol,step
           end if
        else
!!!$     kein FPI
           if (lrobust) then
              write(fpinv,110,err=1000)'IT',it,nrmsd,bdpar,lam,rough,&
                   ncg,betrms,pharms,npol,l1rat,step
           else
              write(fpinv,111,err=1000)'IT',it,nrmsd,bdpar,lam,rough,&
                   ncg,betrms,pharms,npol,step
           end if

        end if

        write(fpinv,'(a)',err=1000)cdump

        write(fpcjg,*,err=1000)
        write(fpcjg,*,err=1000) it
        do k=1,ncg
           write(fpcjg,*,err=1000) cgres(k+1)
        end do
!!!$     lambda search and steplength search
     else
        IF (lip) THEN
           write(fpinv,105,err=1000)'PUP',itr,nrmsd,bdpar,lam,rough,ncg,step
        ELSE 
           write(fpinv,105,err=1000)'UP',itr,nrmsd,bdpar,lam,rough,ncg,step
        END IF
     end if
  end if

!!!$     'inv.ctr' schliessen
  close(fpinv)
  close(fpcjg)

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine kont2

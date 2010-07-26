      subroutine kont2(lsetup)
      
c     Unterprogramm zur Ausgabe der Kontrollvariablen.

c     Andreas Kemna                                            16-Apr-1996
c     Letzte Aenderung   11-Nov-1997

c.....................................................................
      USE datmod

      IMPLICIT none

      INCLUDE 'err.fin'
      INCLUDE 'path.fin'
      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
c     mw
      INCLUDE 'inv.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Hilfsschalter
      logical         * 4     lsetup

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i,k
c     Platzhalter
      INTEGER,PARAMETER :: ncdump=110
      CHARACTER(ncdump) :: cdump
c.....................................................................
 100  FORMAT (t1,a3, t5,i3,t11,g10.4,t69,g10.4,t81,g10.4,
     1     t93,i4,t105,g9.3)
 101  FORMAT (t1,a3,t5,i3,t11,g10.4,t69,g10.4,t81,g10.4,
     1     t93,i4)

 110  FORMAT (t1,a3,t5,i3,t11,g10.4,t23,g9.3,t34,g10.4,
     1     t46,g10.4,t58,i4,t69,g10.4,t81,g10.4,
     1     t93,i4,t105,g9.3,t117,f5.3)
 111  FORMAT (t1,a3,t5,i3,t11,g10.4,t23,g9.3,t34,g10.4,
     1     t46,g10.4,t58,i4,t69,g10.4,t81,g10.4,
     1     t93,i4,t105,f5.3)

 105  FORMAT (t1,a3,t5,i3,t11,g10.4,t23,g9.3,t34,g10.4,
     1     t46,g10.4,t58,i4,t105,f5.3)

c     'inv.ctr' oeffnen
      DO i=1,ncdump-1
         cdump(i:i+1)='*'
      END DO
      errnr = 1
      fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
      open(fpinv,file=fetxt,status='old',err=1000,position='append')
      fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
      open(fpcjg,file=fetxt,status='old',err=1000,position='append')
      errnr = 4

c     Erste Iteration (?)
      if (lsetup) then
         IF (lip) THEN          ! FPI?
            write(fpinv,'(a)',err=1000)cdump
c     Robuste Inversion
            if (lrobust) then
               write(fpinv,100,err=1000)
     1              'PIT',it,nrmsd,betrms,pharms,npol,l1rat
            else
               write(fpinv,101,err=1000)
     1              'PIT',it,nrmsd,betrms,pharms,npol
            end if
            write(fpinv,'(a)',err=1000)cdump
         ELSE
            write(fpinv,'(a)',err=1000)cdump
c     Robuste Inversion
            if (lrobust) then
               write(fpinv,100,err=1000)
     1              'IT',it,nrmsd,betrms,pharms,npol,l1rat
            else
               write(fpinv,101,err=1000)
     1              'IT',it,nrmsd,betrms,pharms,npol
            end if
            write(fpinv,'(a)',err=1000)cdump
         END IF
      else
c     
         ncg = int(cgres(1))

c     Hauptiterationen
         if (llam.and..not.lstep) then
            write(fpinv,'(a)',err=1000)cdump
            if(lip) then
               if (lrobust) then
                  write(fpinv,110,err=1000)
     1                 'PIT',it,nrmsd,bdpar,lam,rough,ncg,betrms,
     1                 pharms,npol,l1rat,step
               else
                  write(fpinv,111,err=1000)
     1                 'PIT',it,nrmsd,bdpar,lam,rough,ncg,betrms,
     1                 pharms,npol,step
               end if
            else
c     kein FPI
               if (lrobust) then
                  write(fpinv,110,err=1000)
     1                 'IT',it,nrmsd,bdpar,lam,rough,ncg,betrms,
     1                 pharms,npol,l1rat,step
               else
                  write(fpinv,111,err=1000)
     1                 'IT',it,nrmsd,bdpar,lam,rough,ncg,betrms,
     1                 pharms,npol,step
               end if

            end if

            write(fpinv,'(a)',err=1000)cdump

            write(fpcjg,*,err=1000)
            write(fpcjg,*,err=1000) it
            do k=1,ncg
               write(fpcjg,*,err=1000) cgres(k+1)
            end do
c     lambda search and steplength search
         else
            IF (lip) THEN
               write(fpinv,105,err=1000)
     1              'PUP',itr,nrmsd,bdpar,lam,rough,ncg,step
            ELSE 
               write(fpinv,105,err=1000)
     1              'UP',itr,nrmsd,bdpar,lam,rough,ncg,step
            END IF
         end if
      end if

c     'inv.ctr' schliessen
      close(fpinv)
      close(fpcjg)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

      subroutine kont2(lsetup)
      
c     Unterprogramm zur Ausgabe der Kontrollvariablen.

c     Andreas Kemna                                            16-Apr-1996
c     Letzte Aenderung   11-Nov-1997

c.....................................................................

      IMPLICIT none
      INCLUDE 'err.fin'
      INCLUDE 'path.fin'
      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'dat.fin'

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

      if (lsetup) then

         write(fpinv,'(a)',err=1000)cdump
         if (lrobust) then
            write(fpinv,'(t1,i3,t7,g10.4,t65,g10.4,t77,g10.4,
     1t89,i4,t101,g9.3)',err=1000)
     1           it,nrmsd,betrms,pharms,npol,l1rat
         else
            write(fpinv,'(t1,i3,t7,g10.4,t65,g10.4,t77,g10.4,
     1t89,i4)',err=1000)
     1           it,nrmsd,betrms,pharms,npol
         end if
         write(fpinv,'(a)',err=1000)cdump
      else
         ncg = int(cgres(1))

         if (llam.and..not.lstep) then
            write(fpinv,'(a)',err=1000)cdump
            if (lrobust) then
               write(fpinv,'(t1,i3,t7,g10.4,t19,g9.3,t30,g10.4,
     1t42,g10.4,t54,i4,t65,g10.4,t77,g10.4,
     1t89,i4,t101,g9.3,t113,f5.3)',err=1000)
     1              it,nrmsd,bdpar,lam,rough,ncg,betrms,
     1              pharms,npol,l1rat,step
            else
               write(fpinv,'(t1,i3,t7,g10.4,t19,g9.3,t30,g10.4,
     1t42,g10.4,t54,i4,t65,g10.4,t77,g10.4,
     1t89,i4,t101,f5.3)',err=1000)
     1              it,nrmsd,bdpar,lam,rough,ncg,betrms,
     1              pharms,npol,step
            end if
            write(fpinv,'(a)',err=1000)cdump

            write(fpcjg,*,err=1000)
            write(fpcjg,*,err=1000) it
            do k=1,ncg
               write(fpcjg,*,err=1000) cgres(k+1)
            end do
         else

            write(fpinv,'(t1,i3,t7,g10.4,t19,g9.3,t30,g10.4,
     1t42,g10.4,t54,i4,t101,f5.3)',err=1000)
     1           itr,nrmsd,bdpar,lam,rough,ncg,step
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

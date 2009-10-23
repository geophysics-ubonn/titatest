        subroutine kont2(lsetup)
     
c Unterprogramm zur Ausgabe der Kontrollvariablen.

c Andreas Kemna                                            16-Apr-1996
c                                       Letzte Aenderung   11-Nov-1997

c.....................................................................

        INCLUDE 'err.fin'
        INCLUDE 'path.fin'
        INCLUDE 'parmax.fin'
        INCLUDE 'konv.fin'
        INCLUDE 'dat.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Hilfsschalter
        logical         * 4     lsetup

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Indexvariable
        integer         * 4     k
c Platzhalter
        INTEGER,PARAMETER :: ncdump=97
        CHARACTER(ncdump) :: cdump
c.....................................................................

c 'inv.ctr' oeffnen
        DO i=1,ncdump-1
           cdump(i:i+1)='*'
        END DO
        errnr = 1
        fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
        open(13,file=fetxt,status='old',err=1000,position='append')
c$$$1       read(13,*,end=2)
c$$$            goto 1
c$$$2       backspace(13)
        errnr = 4

        if (lsetup) then

           write(13,'(a)',err=1000)cdump
            if (lrobust) then
                write(13,'(t1,i3,t7,g10.4,t65,g10.4,t77,g10.4,
     1                     t89,i4,t101,g9.3)',err=1000)
     1                   it,nrmsd,betrms,pharms,npol,l1rat
            else
                write(13,'(t1,i3,t7,g10.4,t65,g10.4,t77,g10.4,
     1                     t89,i4)',err=1000)
     1                   it,nrmsd,betrms,pharms,npol
            end if
            write(13,'(a)',err=1000)cdump
        else
            ncg = int(cgres(1))

            if (llam.and..not.lstep) then
                write(13,'(a)',err=1000)cdump
                if (lrobust) then
                   write(13,'(t1,i3,t7,g10.4,t19,g9.3,t30,g10.4,
     1                  t42,g10.4,t54,i4,t65,g10.4,t77,g10.4,
     1                  t89,i4,t101,g9.3)',err=1000)
     1                  it,nrmsd,step,lam,rough,ncg,betrms,
     1                  pharms,npol,l1rat,bdpar
                else
                   write(13,'(t1,i3,t7,g10.4,t19,g9.3,t30,g10.4,
     1                  t42,g10.4,t54,i4,t65,g10.4,t77,g10.4,
     1                  t89,i4,t95,g9.3)',err=1000)
     1                  it,nrmsd,step,lam,rough,ncg,betrms,
     1                  pharms,npol,bdpar
                end if
                write(13,'(a)',err=1000)cdump

                fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
                write(14,*,err=1000)
                write(14,*,err=1000) it
                do k=1,ncg
                    write(14,*,err=1000) cgres(k+1)
                end do
            else
                write(13,'(t1,i3,t7,g10.4,t19,g9.3,t30,g10.4,
     1                     t42,g10.4,t54,i4)',err=1000)
     1                   it,nrmsd,step,lam,rough,ncg
            end if
        end if

c 'inv.ctr' schliessen
        close(13)

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

1000    return

        end

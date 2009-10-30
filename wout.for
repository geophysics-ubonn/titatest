      subroutine wout(kanal,dsigma,dvolt)

c     Unterprogramm zum Schreiben der Widerstandsverteilung und der
c     modellierten Daten inkl. Elektrodenkennungen.

c     Andreas Kemna                                            28-Sep-1994
c     Letzte Aenderung   10-Mar-2007
      
c.....................................................................
      
      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'path.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'fem.fin'
      INCLUDE 'konv.fin'
c     diff+<
      INCLUDE 'model.fin'
      INCLUDE 'inv.fin'
c     diff+>

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Dateinamen
      character       * 80    dsigma,
     1     dvolt

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c     Indexvariablen
      integer         * 4     i,j,k

c     Hilfsfelder
      real            * 8     xkoord(elmax),
     1     ykoord(elmax)

c     Hilfsvariablen
      real            * 8     xdum,ydum
      integer         * 4     idum,idum2,lit
      character       * 80    htxt
      character       * 12    htxt2
      complex         * 16    dum

c     Hilfsfunctions
      character       * 12    intcha
      character       * 80    filpat

c     diff+<
      real            * 8     dum2,dum3
c     diff+>
      logical         * 4     exi
c.....................................................................

c     Schwerpunktkoordinaten der Elemente bestimmen
      iel = 0

      do i=1,typanz
         if (typ(i).gt.10) goto 10

         nkel = selanz(i)

         do j=1,nelanz(i)
            iel  = iel + 1

            xdum = 0d0
            ydum = 0d0

            do k=1,nkel
               xdum = xdum + sx(snr(nrel(iel,k)))
               ydum = ydum + sy(snr(nrel(iel,k)))
            end do

            xkoord(iel) = xdum/dble(nkel)
            ykoord(iel) = ydum/dble(nkel)
         end do
      end do

c     'dsigma' modifizieren
 10   lit  = int(log10(real(itmax)))+1
      htxt = filpat(dsigma,idum2,1,slash(1:1))
      idum = idum2+index(dsigma(idum2+1:80),'.')-1

      if ((idum-idum2-1).gt.(8-lit)) then
         fetxt = dsigma
         errnr = 15
         goto 999
      end if

      htxt2 = intcha(it,lit)
      htxt  = dsigma(1:idum)//htxt2(1:lit)//dsigma(idum+1:idum+4)

c     Betraege ausgeben
      idum  = index(htxt,' ')
      fetxt = htxt(1:idum-4)//'mag'
      OPEN (kanal,FILE='tmp.lastmod',STATUS='replace')
      WRITE (kanal,*)TRIM(fetxt)
      CLOSE (kanal)
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4
      write(kanal,*,err=1000) elanz

      do i=1,elanz
c     diff+<
         if (.not.ldiff) then
c     diff+>
            dum = dcmplx(1d0)/sigma(i)
c     ak            write(kanal,*,err=1000) real(xkoord(i)),real(ykoord(i)),
c     ak     1                              real(cdabs(dum))
            write(kanal,*,err=1000) real(xkoord(i)),real(ykoord(i)),
     1           real(dlog10(cdabs(dum)))
c     ro            write(kanal,*,err=1000) real(ykoord(i)),real(xkoord(i)),
c     ro     1                              real(cdabs(dum))
c     diff+<
         else
            dum3 = cdabs(dcmplx(1d0)/sigma(i))
            dum2 = cdabs(dcmplx(1d0)/cdexp(m0(mnr(i))))
            write(kanal,'(7(f10.4,2x))',err=1000)
     1           real(xkoord(i)),real(ykoord(i)),
     1           real(dlog10(dum3)),
     1           real(1d2*(dum3/dum2-1d0)),
     1           real(dum3-dum2),
     1           real(1d2*(dum2/dum3-1d0)),
     1           real(1d3/dum3-1d3/dum2)
         end if
c     diff+>
      end do
      close(kanal)
      fetxt = htxt(1:idum-4)//'modl'

c$$$  
c$$$  INQUIRE (FILE='tmp.matinfo',EXIST=exi)
c$$$  IF (exi) THEN
c$$$  OPEN (kanal,FILE='tmp.matinfo',STATUS='old',
c$$$  1        POSITION='append')
c$$$  WRITE (kanal,*)TRIM(fetxt)
c$$$  CLOSE (kanal)
c$$$  ELSE
c$$$  OPEN (kanal,FILE='tmp.matinfo',STATUS='new')
c$$$  WRITE (kanal,*)TRIM(fetxt)
c$$$  CLOSE (kanal)
c$$$  END IF

      OPEN (kanal,FILE=fetxt,status='replace')
      WRITE (kanal,'(I7)',err=1000) elanz
      DO i=1,elanz
         WRITE (kanal,'(2(1x,G12.4))')1./REAL(sigma(i)),0.0
      END DO
      CLOSE (kanal)
c     Ggf. Phasen ausgeben
      if (.not.ldc) then
         fetxt = htxt(1:idum-4)//'pha'
         errnr = 1
         open(kanal,file=fetxt,status='replace',err=999)
         errnr = 4
         write(kanal,*,err=1000) elanz

         do i=1,elanz
            dum = dcmplx(1d0)/sigma(i)
            write(kanal,*,err=1000)
c     ak Default
     1           real(xkoord(i)),real(ykoord(i)),
     1           real(1d3*datan2(dimag(dum),dble(dum)))
c     ak MMAJ
c     ak     1                      real(xkoord(i)),real(ykoord(i)),
c     ak     1                      -real(1d3*datan2(dimag(dum),dble(dum)))
c     ro     1                      real(ykoord(i)),real(xkoord(i)),
c     ro     1                      -real(1d3*datan2(dimag(dum),dble(dum)))
         end do
         close(kanal)
      end if

c     'dvolt' modifizieren
      htxt = filpat(dvolt,idum2,1,slash(1:1))
      idum = idum2+index(dvolt(idum2+1:80),'.')-1

      if ((idum-idum2-1).gt.(8-lit)) then
         fetxt = dvolt
         errnr = 15
         goto 999
      end if

      htxt2 = intcha(it,lit)
      htxt  = dvolt(1:idum)//htxt2(1:lit)//dvolt(idum+1:idum+4)

      fetxt = htxt
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)

      errnr = 4
      write(kanal,*,err=1000) nanz

c     Stromelektrodennummern, Spannungselektrodennummern und scheinbare
c     Widerstandswerte (Betrag und Phase (in mrad)) schreiben
      if (ldc) then
         do i=1,nanz
            write(kanal,*,err=1000)
     1           strnr(i),vnr(i),
c     diff-     1                      real(1d0/dexp(dble(sigmaa(i))))
c     diff+<
     1           real(1d0/dexp(dble(sigmaa(i)))),wdfak(i)
c     diff+>
         end do
      else
         do i=1,nanz
            write(kanal,*,err=1000)
     1           strnr(i),vnr(i),
     1           real(1d0/dexp(dble(sigmaa(i)))),
c     diff-     1                      real(-1d3*dimag(sigmaa(i)))
c     diff+<
     1           real(-1d3*dimag(sigmaa(i))),wdfak(i)
c     diff+>
         end do
      end if

      close(kanal)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end

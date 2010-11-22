      subroutine wout(kanal,dsigma,dvolt)

!!!$     Unterprogramm zum Schreiben der Widerstandsverteilung und der
!!!$     modellierten Daten inkl. Elektrodenkennungen.

!!!$     Andreas Kemna                                            28-Sep-1994
!!!$     Letzte Aenderung   10-Mar-2007
      
!!!$.....................................................................
      
      USE datmod
      USE femmod
      USE invmod
      USE sigmamod
      USE modelmod
      USE elemmod
      USE errmod
      USE konvmod
      USE pathmod

      IMPLICIT none

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
      integer         * 4     kanal

!!!$     Dateinamen
      character       * 80    dsigma,
     1     dvolt

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Aktuelle Elementnummer
      integer         * 4     iel

!!!$     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

!!!$     Indexvariablen
      integer         * 4     i,j,k

!!!$     Hilfsvariablen
      real            * 8     xdum,ydum
      integer         * 4     idum,idum2,lit
      character       * 80    htxt
      character       * 12    htxt2
      complex         * 16    dum

!!!$     Hilfsfunctions
      character       * 12    intcha
      character       * 80    filpat

!!!$     diff+<
      real            * 8     dum2,dum3
!!!$     diff+>
      logical         * 4     exi
      CHARACTER       * 2     ci
!!!$.....................................................................

!!!$     'dsigma' modifizieren

 10   IF (it>9) THEN
         WRITE (ci,'(I2)')it
      ELSE 
         WRITE (ci,'(a,I1)')'0',it
      END IF
      htxt = filpat(dsigma,idum2,1,slash(1:1))
      idum = idum2+index(dsigma(idum2+1:80),'.')-1
      htxt  = dsigma(1:idum)//ci//dsigma(idum+1:idum+4)

!!!$     Betraege ausgeben
      idum  = index(htxt,' ')
      fetxt = htxt(1:idum-4)//'mag'
      OPEN (kanal,FILE='inv.lastmod',STATUS='replace')
      WRITE (kanal,*)TRIM(fetxt)
      CLOSE (kanal)
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4
      write(kanal,*,err=1000) elanz,ACHAR(9),betrms

      do i=1,elanz
!!!$     diff+<
         if (.not.ldiff.AND..NOT.lprior) then
c$$$  if (.not.ldiff) then
!!!$     diff+>
            dum = dcmplx(1d0)/sigma(i)
!!!$     ak            write(kanal,*,err=1000) real(espx(i)),real(espy(i)),
!!!$     ak     1                              real(cdabs(dum))
            write(kanal,'(3F10.4,2X)',err=1000) 
     1           real(espx(i)),real(espy(i)),
     1           real(dlog10(cdabs(dum)))
!!!$     ro            write(kanal,*,err=1000) real(espy(i)),real(espx(i)),
!!!$     ro     1                              real(cdabs(dum))
!!!$     diff+<
         ELSE IF (lprior) THEN
            dum3 = CDABS(dcmplx(1d0)/sigma(i))
            dum2 = CDABS(dcmplx(1d0)/cdexp(m0(mnr(i))))
!     dum3 = REAL(CDLOG(sigma(i))/m0(mnr(i)))
            write(kanal,'(7(f10.4,2x))',err=1000)
     1           REAL(espx(i)),
     1           REAL(espy(i)),
     1           real(dlog10(dum3)),
     1           (DLOG10(dum3) - DLOG10(dum2)),
     1           real(dlog10(dum2)),
     1           real(1d2*(1d0-dum2/dum3)),
     1           real(1d2*(1d0-dum3/dum2))
         else
            dum3 = cdabs(dcmplx(1d0)/sigma(i))
            dum2 = cdabs(dcmplx(1d0)/cdexp(m0(mnr(i))))
            write(kanal,'(7(f10.4,2x))',err=1000)
     1           real(espx(i)),real(espy(i)),
     1           real(dlog10(dum3)),
     1           real(1d2*(dum3/dum2-1d0)),
     1           real(dum3-dum2),
     1           real(1d2*(dum2/dum3-1d0)),
     1           real(1d3/dum3-1d3/dum2)
         end if
!!!$     diff+>
      end do
      close(kanal)
      fetxt = htxt(1:idum-4)//'modl'
      OPEN (kanal,FILE=fetxt,status='replace')

      WRITE (kanal,'(I7)',err=1000) elanz
      IF (.NOT.lprior) then
         DO i=1,elanz
            dum = dcmplx(1d0)/sigma(i)
            dum2 = real(1d3*datan2(dimag(dum),dble(dum)))
            WRITE (kanal,'(2(1x,G12.4))')1./REAL(sigma(i)),dum2
         END DO
      ELSE
         DO i=1,elanz
            dum = dcmplx(1d0)/sigma(i) ! - DCMPLX(1d0)/CDEXP(m0(mnr(i)))
            dum2 = real(1d3*datan2(dimag(dum),dble(dum)))
            WRITE (kanal,'(2(1x,G12.4))')1./REAL(sigma(i)),dum2
         END DO
      END IF
      CLOSE (kanal)
!!!$     Ggf. Phasen ausgeben
      if (.not.ldc) then
         fetxt = htxt(1:idum-4)//'pha'
         errnr = 1
         open(kanal,file=fetxt,status='replace',err=999)
         errnr = 4
         write(kanal,*,err=1000) elanz,ACHAR(9),pharms

         do i=1,elanz
            dum = dcmplx(1d0)/sigma(i)
            write(kanal,*,err=1000)
!!!$     ak Default
     1           real(espx(i)),real(espy(i)),
     1           real(1d3*datan2(dimag(dum),dble(dum)))
!!!$     ak MMAJ
!!!$     ak     1                      real(espx(i)),real(espy(i)),
!!!$     ak     1                      -real(1d3*datan2(dimag(dum),dble(dum)))
!!!$     ro     1                      real(espy(i)),real(espx(i)),
!!!$     ro     1                      -real(1d3*datan2(dimag(dum),dble(dum)))
         end do
         close(kanal)
      end if

!!!$     'dvolt' modifizieren
      htxt = filpat(dvolt,idum2,1,slash(1:1))
      idum = idum2+index(dvolt(idum2+1:80),'.')-1

      htxt  = dvolt(1:idum)//ci//dvolt(idum+1:idum+4)

      fetxt = htxt
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)

      errnr = 4
      write(kanal,*,err=1000) nanz

!!!$     Stromelektrodennummern, Spannungselektrodennummern und scheinbare
!!!$     Widerstandswerte (Betrag und Phase (in mrad)) schreiben
      if (ldc) then
         do i=1,nanz
            write(kanal,*,err=1000)
     1           strnr(i),vnr(i),
!!!$     diff-     1                      real(1d0/dexp(dble(sigmaa(i))))
!!!$     diff+<
     1           real(1d0/dexp(dble(sigmaa(i)))),wdfak(i)
!!!$     diff+>
         end do
      else
         do i=1,nanz
            write(kanal,*,err=1000)
     1           strnr(i),vnr(i),
     1           real(1d0/dexp(dble(sigmaa(i)))),
!!!$     diff-     1                      real(-1d3*dimag(sigmaa(i)))
!!!$     diff+<
     1           real(-1d3*dimag(sigmaa(i))),wdfak(i)
!!!$     diff+>
         end do
      end if

      close(kanal)

      errnr = 0
      return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end

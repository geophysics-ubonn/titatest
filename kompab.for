      subroutine kompab(nelec,ki)

c     Unterprogramm zur Kompilation der FE-Matrix 'a' in Bandform
c     (vorgegebene Bandbreite 'mb') und des Konstantenvektors 'b'
c     ( A * x + b = 0 ).

c     Andreas Kemna                                            17-Dec-1993
c     Letzte Aenderung   15-Jul-2007
      
c.....................................................................

      USE alloci
      USE femmod
      USE sigmamod
      USE electrmod
      USE elemmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'waven.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Aktuelle Elektrodennummer
      integer         * 4     nelec

c     Aktueller Wellenzahlindex
      integer         * 4     ki

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Aktuelle Randelementnummer
      integer         * 4     rel

c     Aktueller Elementtyp
      integer         * 4     ntyp

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c     Hilfsvariablen
      real            * 8     dum
      complex         * 16    dum2
      integer         * 4     im,imax,imin
      integer         * 4     nzp,nnp,idif,ikl,idum

c     Indexvariablen
      integer         * 4     i,j,k,l

c.....................................................................

c     Gesamtsteifigkeitsmatrix und Konstantenvektor auf Null setzen
      im = (mb+1)*sanz

      do i=1,im
         a(i) = dcmplx(0d0)
      end do

      do i=1,sanz
         b(i) = dcmplx(0d0)
      end do

      iel = 0

      do i=1,typanz
         ntyp = typ(i)
         nkel = selanz(i)

         do 10 j=1,nelanz(i)
            iel = iel + 1
            ikl = 0

            if (ntyp.gt.11) goto 10

            do k=1,nkel
               nzp = nrel(iel,k)

               do l=1,k
                  nnp  = nrel(iel,l)
                  idif = iabs(nzp-nnp)

c     Ggf. Fehlermeldung
                  if (idif.gt.mb) then
                     WRITE (fetxt,*)
     1                    'kompab idif',idif,' iel ',iel
                     fetxt = ' '
                     errnr = 19
                     goto 1000

                  else

c     Aufbau der Gesamtsteifigkeitsmatrix und ggf. des Konstantenvektors
                     ikl = ikl + 1

                     imax = max0(nzp,nnp)
                     imin = min0(nzp,nnp)
                     im   = imax*mb + imin

                     if (ntyp.eq.11) then
                        rel  = iel - elanz
                        dum  = relbg(rel,ikl) * kg(rel,nelec,ki)
                        idum = rnr(rel)
                     else
                        dum  = elbg(iel,ikl,ki)
                        idum = iel
                     end if

                     a(im) = a(im) + dcmplx(dum) * sigma(idum)

                     if (lsr) then
                        dum2   = dcmplx(dum) * (sigma(idum)-sigma0)
                        b(nzp) = b(nzp) + dum2 * pota(nnp)
                        if (nnp.ne.nzp)
     1                       b(nnp) = b(nnp) + dum2 * pota(nzp)
                     end if
                  end if
               end do
            end do
 10      continue
      end do

c     Ggf. Konstantenvektor belegen
      if (.not.lsr) b(enr(nelec)) = dcmplx(-1d0)

c     akc BAW-Tank
c     ak        b(211) = dcmplx(1d0)
c     akc Model EGS2003
c     ak        b(1683) = dcmplx(1d0)
c     akc Lysimeter hor_elem\normal
c     ak        b(129) = dcmplx(1d0)
c     akc Lysimeter hor_elem\fine
c     ak        b(497) = dcmplx(1d0)
c     akc Simple Tucson Model
c     ak        b(431) = dcmplx(1d0)
c     akc TU Berlin Mesokosmos
c     ak        b(201) = dcmplx(1d0)
c     akc Andy
c     ak        b(2508) = dcmplx(1d0)
c     akc Sandra (ele?_anom)
c     ak        b(497) = dcmplx(1d0)
c     akc Adrian (Tank)
c     ak        b(1660) = dcmplx(1d0)

      if (lsink) b(nsink) = dcmplx(1d0)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

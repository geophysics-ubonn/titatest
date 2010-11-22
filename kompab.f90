      subroutine kompab(nelec,ki)

!!!$     Unterprogramm zur Kompilation der FE-Matrix 'a' in Bandform
!!!$     (vorgegebene Bandbreite 'mb') und des Konstantenvektors 'b'
!!!$     ( A * x + b = 0 ).

!!!$     Andreas Kemna                                            17-Dec-1993
!!!$     Letzte Aenderung   15-Jul-2007
      
!!!$.....................................................................

      USE alloci
      USE femmod
      USE sigmamod
      USE electrmod
      USE elemmod
      USE wavenmod
      USE errmod

      IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Aktuelle Elektrodennummer
      integer         * 4     nelec

!!!$     Aktueller Wellenzahlindex
      integer         * 4     ki

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Aktuelle Elementnummer
      integer         * 4     iel

!!!$     Aktuelle Randelementnummer
      integer         * 4     rel

!!!$     Aktueller Elementtyp
      integer         * 4     ntyp

!!!$     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

!!!$     Hilfsvariablen
      real            * 8     dum
      complex         * 16    dum2
      integer         * 4     im,imax,imin
      integer         * 4     nzp,nnp,idif,ikl,idum

!!!$     Indexvariablen
      integer         * 4     i,j,k,l

!!!$.....................................................................

!!!$     Gesamtsteifigkeitsmatrix und Konstantenvektor auf Null setzen
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

!!!$     Ggf. Fehlermeldung
                  if (idif.gt.mb) then
                     WRITE (fetxt,*)
     1                    'kompab idif',idif,' iel ',iel
                     fetxt = ' '
                     errnr = 19
                     goto 1000

                  else

!!!$     Aufbau der Gesamtsteifigkeitsmatrix und ggf. des Konstantenvektors
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

!!!$     Ggf. Konstantenvektor belegen
      if (.not.lsr) b(enr(nelec)) = dcmplx(-1d0)

!!!$     akc BAW-Tank
!!!$     ak        b(211) = dcmplx(1d0)
!!!$     akc Model EGS2003
!!!$     ak        b(1683) = dcmplx(1d0)
!!!$     akc Lysimeter hor_elem\normal
!!!$     ak        b(129) = dcmplx(1d0)
!!!$     akc Lysimeter hor_elem\fine
!!!$     ak        b(497) = dcmplx(1d0)
!!!$     akc Simple Tucson Model
!!!$     ak        b(431) = dcmplx(1d0)
!!!$     akc TU Berlin Mesokosmos
!!!$     ak        b(201) = dcmplx(1d0)
!!!$     akc Andy
!!!$     ak        b(2508) = dcmplx(1d0)
!!!$     akc Sandra (ele?_anom)
!!!$     ak        b(497) = dcmplx(1d0)
!!!$     akc Adrian (Tank)
!!!$     ak        b(1660) = dcmplx(1d0)

      if (lsink) b(nsink) = dcmplx(1d0)

      errnr = 0
      return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 1000 return

      end

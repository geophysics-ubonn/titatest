subroutine kompab(nelec,ki,a_komp,b_komp)

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

  COMPLEX (KIND (0D0)),DIMENSION(*) ::     a_komp
  COMPLEX (KIND (0D0)),DIMENSION(*) ::     b_komp

!!!$     Aktuelle Elektrodennummer
  INTEGER (KIND = 4) ::     nelec

!!!$     Aktueller Wellenzahlindex
  INTEGER (KIND = 4) ::     ki

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Aktuelle Elementnummer
  INTEGER (KIND = 4) ::     iel

!!!$     Aktuelle Randelementnummer
  INTEGER (KIND = 4) ::     rel

!!!$     Aktueller Elementtyp
  INTEGER (KIND = 4) ::     ntyp

!!!$     Anzahl der Knoten im aktuellen Elementtyp
  INTEGER (KIND = 4) ::     nkel

!!!$     Hilfsvariablen
  REAL (KIND(0D0))   ::     dum
  COMPLEX(KIND(0D0)) ::    dum2
  INTEGER (KIND = 4) ::     im,imax,imin
  INTEGER (KIND = 4) ::     nzp,nnp,idif,ikl,idum

!!!$     Indexvariablen
  INTEGER (KIND = 4) ::     i,j,k,l

!!!$.....................................................................

!!!$     Gesamtsteifigkeitsmatrix und Konstantenvektor auf Null setzen
  im = (mb+1)*sanz
  
! a_komp = DCMPLX(0D0)

  do i=1,im
     a_komp(i) = dcmplx(0d0)
  end do

! b_komp = DCMPLX(0D0)

  do i=1,sanz
     b_komp(i) = dcmplx(0d0)
  end do

  iel = 0

  do i=1,typanz
     ntyp = typ(i)
     nkel = selanz(i)

     do j=1,nelanz(i)
        iel = iel + 1
        ikl = 0

        if (ntyp.gt.11) CYCLE

        do k=1,nkel
           nzp = nrel(iel,k)

           do l=1,k
              nnp  = nrel(iel,l)
              idif = iabs(nzp-nnp)

!!!$     Ggf. Fehlermeldung
              if (idif.gt.mb) then
                 WRITE (fetxt,*)'kompab idif',idif,' iel ',iel
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

                 a_komp(im) = a_komp(im) + dcmplx(dum) * sigma(idum)

                 if (lsr) then
                    dum2   = dcmplx(dum) * (sigma(idum)-sigma0)
                    b_komp(nzp) = b_komp(nzp) + dum2 * pota(nnp)
                    if (nnp.ne.nzp) b_komp(nnp) = b_komp(nnp) + dum2 * pota(nzp)
                 end if
              end if
           end do
        end do
     END do
  end do

!!!$     Ggf. Konstantenvektor belegen
  if (.not.lsr) b_komp(enr(nelec)) = dcmplx(-1d0)

!!!$     akc BAW-Tank
!!!$     ak        b_komp(211) = dcmplx(1d0)
!!!$     akc Model EGS2003
!!!$     ak        b_komp(1683) = dcmplx(1d0)
!!!$     akc Lysimeter hor_elem\normal
!!!$     ak        b_komp(129) = dcmplx(1d0)
!!!$     akc Lysimeter hor_elem\fine
!!!$     ak        b_komp(497) = dcmplx(1d0)
!!!$     akc Simple Tucson Model
!!!$     ak        b_komp(431) = dcmplx(1d0)
!!!$     akc TU Berlin Mesokosmos
!!!$     ak        b_komp(201) = dcmplx(1d0)
!!!$     akc Andy
!!!$     ak        b_komp(2508) = dcmplx(1d0)
!!!$     akc Sandra (ele?_anom)
!!!$     ak        b_komp(497) = dcmplx(1d0)
!!!$     akc Adrian (Tank)
!!!$     ak        b_komp(1660) = dcmplx(1d0)

  if (lsink) b_komp(nsink) = dcmplx(1d0)

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine kompab

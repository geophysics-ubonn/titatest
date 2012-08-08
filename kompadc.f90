subroutine kompadc(nelec,ki)

!!!$     Unterprogramm zur Kompilation der FE-Matrix 'adc' in Bandform
!!!$     (vorgegebene Bandbreite 'mb') und des Konstantenvektors 'bdc'
!!!$     ( A * x + b = 0 ).

!!!$     Andreas Kemna                                            17-Dec-1993
!!!$     Letzte Aenderung   16-Jul-2007

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
  REAL (KIND (0D0)) ::     dum
  REAL (KIND (0D0)) ::     dum2
  INTEGER (KIND = 4) ::     im,imax,imin
  INTEGER (KIND = 4) ::     nzp,nnp,idif,ikl,idum

!!!$     Indexvariablen
  INTEGER (KIND = 4) ::     i,j,k,l

!!!$.....................................................................

!!!$     Gesamtsteifigkeitsmatrix und Konstantenvektor auf Null setzen
  im = (mb+1)*sanz

  adc = 0D0
  bdc = 0D0

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
                    dum2 = DBLE(sigma(rnr(rel)))
                    dum2 = DBLE(sigma0)
                    dum2 = 0d0 ! which removes the influence
                 else
                    dum  = elbg(iel,ikl,ki)
                    dum2 = DBLE(sigma(iel))
                 end if
!!!$ BUG
!!$ THE PROBLEM IS HERE, IF THERE ARE TOO FEW OUTER GRID CELLS
!!!$ SIGMA IS NO MORE DENFINED!!!

                 a(im) = a(im) + dum * dum2

                 if (lsr) then
                    dum2   = dum * DBLE(dum2 - sigma0)
                    bdc(nzp) = bdc(nzp) + dum2 * dble(pota(nnp))
                    if (nnp.ne.nzp) bdc(nnp) = bdc(nnp) + dum2 * &
                         dble(pota(nzp))
                 end if

              end if

           end do ! l=1,k
        end do ! k=1,nkel
     END do ! j=1,nelanz(i)
  end do ! i=1,typanz
  
!!!$     Ggf. Konstantenvektor belegen
  if (.not.lsr) bdc(enr(nelec)) = -1d0

!!!$     akc BAW-Tank
!!!$     ak        bdc(211) = 1d0
!!!$     akc Model EGS2003
!!!$     ak        bdc(1683) = 1d0
!!!$     akc Lysimeter hor_elem\normal
!!!$     ak        bdc(129) = 1d0
!!!$     akc Lysimeter hor_elem\fine
!!!$     ak        bdc(497) = 1d0
!!!$     akc Simple Tucson Model
!!!$     ak        bdc(431) = 1d0
!!!$     akc TU Berlin Mesokosmos
!!!$     ak        bdc(201) = 1d0
!!!$     akc Andy
!!!$     ak        bdc(2508) = 1d0
!!!$     akc Sandra (ele?_anom)
!!!$     ak        bdc(497) = 1d0

  if (lsink) bdc(nsink) = 1d0

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine kompadc

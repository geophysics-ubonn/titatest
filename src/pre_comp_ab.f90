subroutine pre_comp_ab(ki,my_a_mat_band)

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

  COMPLEX (prec),Dimension(2*mb+1,sanz) :: my_a_mat_band

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
  REAL (prec)   ::     dum
  COMPLEX(prec) ::     dum2,value
  INTEGER (KIND = 4) ::     im,imax,imin
  INTEGER (KIND = 4) ::     nzp,nnp,idif,ikl,idum
  INTEGER :: la_Kd,la_N,la_i,la_j

!!!$     Indexvariablen
  INTEGER (KIND = 4) ::     i,j,k,l,index_i

!!!$.....................................................................

!!!$     Gesamtsteifigkeitsmatrix und Konstantenvektor auf Null setzen
  my_a_mat_band = dCMPLX(0D0)

  iel = 0
  do i=1,typanz
     ntyp = typ(i)
     nkel = selanz(i)
     do j=1,nelanz(i)
        iel = iel + 1
        ikl = 0
        if (ntyp.gt.8) CYCLE
        do k=1,nkel
           nzp = nrel(iel,k)
           do l=1,nkel
              nnp  = nrel(iel,l)
              idif = iabs(nzp-nnp)  

!!!$     Aufbau der Gesamtsteifigkeitsmatrix und ggf. des Konstantenvektors
                 ikl = ikl + 1
                 dum  = elbg(iel,ikl,ki)
                 dum2 = sigma(iel)
                 value = dcmplx(dum)*dum2
                 ! band matrix index (see above)
                call assign_zgbsvx(my_a_mat_band,nzp,nnp,mb,sanz,value)
!                call assign_zpbsvx(my_a_mat_band,nzp,nnp,mb,sanz,dum2)
           end do
        end do
     END do
  end do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine pre_comp_ab

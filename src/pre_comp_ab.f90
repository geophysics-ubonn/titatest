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

  COMPLEX (prec),Dimension(3*mb+1,sanz) :: my_a_mat_band
!  COMPLEX (prec),DIMENSION(sanz) ::     my_b

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
  COMPLEX(prec) ::     dum2
  INTEGER (KIND = 4) ::     im,imax,imin
  INTEGER (KIND = 4) ::     nzp,nnp,idif,ikl,idum

!!!$     Indexvariablen
  INTEGER (KIND = 4) ::     i,j,k,l,index_i

!!!$.....................................................................

!!!$     Gesamtsteifigkeitsmatrix und Konstantenvektor auf Null setzen
  my_a_mat_band = CMPLX(0D0)

!BAND STORAGE
!An m-by-n band matrix with kl subdiagonals and ku superdiagonals may be stored compactly in a two-dimensional array with kl+ku+1 rows and n columns. Columns of the matrix are stored in corresponding columns of the array, and diagonals of the matrix are stored in rows of the array. This storage scheme should be used in practice only if $kl, ku \ll \min(m,n)$, although LAPACK routines work correctly for all values of kl and ku. In LAPACK, arrays that hold matrices in band storage have names ending in `B'.

!To be precise, aij is stored in AB(ku+1+i-j,j) for $\max(1,j-ku) \leq i \leq \min(m,j+kl)$. For example, when m = n = 5, kl = 2 and ku = 1:
!source: http://www.netlib.org/lapack/lug/node124.html

  iel = 0
  do i=1,typanz
     ntyp = typ(i)
     nkel = selanz(i)
     do j=1,nelanz(i)
        iel = iel + 1
!         print*,'precomp',sigma(iel),iel,ntyp
        ikl = 0
        if (ntyp.gt.8) CYCLE
        do k=1,nkel
           nzp = nrel(iel,k)
           do l=1,nkel
              nnp  = nrel(iel,l)
              idif = iabs(nzp-nnp)

!!!$     Aufbau der Gesamtsteifigkeitsmatrix und ggf. des Konstantenvektors
                 ikl = ikl + 1
!if ((nzp.ge.max(1,nnp-mb)).and.(nzp.le.nnp)) then
!                 imax = max0(nzp,nnp)
!                 imin = min0(nzp,nnp)
!                 im   = imax*mb + imin
                 dum  = elbg(iel,ikl,ki)
                 dum2 = sigma(iel)
                 ! band matrix index (see above)
!               index_i = mb+1+nzp-nnp
!               my_a_mat_band(index_i,nnp) = my_a_mat_band(index_i,nnp) + dum * dum2
!end if
         index_i = mb+mb+1+nzp-nnp
               my_a_mat_band(index_i,nnp) = my_a_mat_band(index_i,nnp) + CMPLX(dum) * dum2 
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

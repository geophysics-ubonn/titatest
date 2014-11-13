!> \file scalab.f90
!> \brief scale the FE equations for better numerical accuracy
!> \details Kemna (2000) Appendix C: The procedure used to evaluate the inverse Fourier integral \f[ \phi_i = \frac{1}{\pi} \int_0^\infty \tilde \phi_i(k)dk. \f] is based on the approach of LaBreque et al. (1996a). Since the asymptotic behavior of the transformed potential \f$\tilde \phi (k)\f$ for small and large arguments is different, the integral is split into two parts to which, respectively, adequate numerical integration techniques are then applied, i.e.,
!> \f[\int_0^\infty \tilde \phi_i(k)dk = \int_0^{k_0} \tilde \phi_i(k)dk + \int_{k_0}^\infty \tilde \phi_i(k)dk, \f]
!> with some characteristic wavenumber \f$ k_0 \f$.
!> According to the analytic solution of the 2D Helmholtz equation,
!> \f[ \tilde \phi_p = \frac{I}{2 \pi \sigma_0} \left( K-= (kr_-) + K_0 (kr_+) \right) ,\f] 
!> the (primary) transformed potential is proportional to the modified Bessel function \f$ K_0(u) \f$, with \f$ u= kr\f$. Here, \f$ r\f$ denotes the representative radial distance from the (image) source. The corresponding asymptotic behavior is given by (Abramowitz and Stegun, 1984)
!> \f[\begin{array}{cc} u \rightarrow 0 : & K_0(u) \propto -\ln u, \\ u \rightarrow \infty : & K_0 \propto e^{-u} / \sqrt{u}  \end{array}\f]
!> To overcome the singularity in the integrand at zero, the change of variable \f$ k' = (k/k_0)^{1/2}\f$ is made in the first integral of the right-hand side of the first equation. The resulting nonsingular integral, with integrand \f$ g(k') = 2 k_0 k' \tilde \phi (k) \f$, is evaluated by conventional \f$N_g-\f$point Gaussian quadrature (Press et al., 1992), yielding appropriate abscissa \f$k'_n\f$ with corresponding weights \f$ w'_n \f$. Summing up both steps, it is
!> \f[ \int_0^{k_0} \tilde \phi_i(k)dk = \int_0^1 g(k')dk' = \sum_{n=1}^{N_G} w'_n g(k'_n) = \pi \sum_{n=1}^{N_G} w_n \tilde \phi (k_n)\f]
!> where \f$ k_n = k_0 k^{'2}_0 \f$ and \f$ w_n = 2 k_0 k'_n w'_n / \pi \f$.
!> From the relations in the 3rd equation, it is seen that the behavior of \f$ \tilde\phi(k) \f$ for large arguments is characterized by an exponential decrease. Therefore, the upper integration in the 2nd equation is performed using a \f$N_L-\f$point Laguerre-type formula (Press et al., 1992). Analogous to the 4th equation, one finds
!> \f[ \int_{k_0}^\infty \tilde \phi(k) dk = \int_0^\infty e^{-k'} g(k') dk' = \sum_{n=1}^{N_L} w'_n g(k'_n) = \pi \sum_{n=1}^{N_L} w_n \tilde \phi (k_n) \f]
!> with the rescaled abscissa and weights, respectively, \f$ k_n = k_0 (k'_n +1) \f$ and \f$ w_n = k_0 e^{k'_n} w'_n /\pi \f$. Note that in the last equation, after substitution of \f$ k' = k/k_0 -1 \f$, the new function \f$ g(k') = k_0 e^{k'} \tilde \phi (k) \f$ is temporarily defined.
!> The integration bound \f$ k_0 \f$ has to be specified in relation to the spatial scale of the considered problem. A characteristic quantity in this regard is given by the minimum distance between corresponding transmitting and receiving electrodes within the underlying set of measurement configurations, denoted by \f$ r_{min} \f$. The number of employed abscissa in the integration formulas determines the valid range of the integration in terms of the distance from the source. From a numerical analysis for a homogeneous half-space, it was found that the variable choice
!> \f[ N_G = \int \[ 6 \log (r_{max}/r_{min}) \] \f]
!> together with \f$ N_L = 4 \f$ and \f$ k_0 = 1/(2r_{min}) \f$, guarantees an error of less than 0.5 % in the inverse Fourier integration for distances ranging from \f$ r_{min}\f$ to \f$r_{max}\f$ (see Figure C.1). Herein,  denotes the maximum distance between source and receiver electrodes as being of interest within the survey.
!> @author Andreas Kemna 
!> @date 10/11/1993

subroutine scalab(a_scal,b_scal,fak_scal)

!     Unterprogramm skaliert 'a' und 'b' und liefert die Skalierungs-
!     faktoren im Vektor 'fak'.

!     ( Vgl. Subroutine 'SCALBNDN' in Schwarz (1991) )

!     Andreas Kemna                                            11-Oct-1993
!     Letzte Aenderung   07-Mar-2003

!.....................................................................

  USE alloci
  USE femmod
  USE elemmod
  USE errmod

  IMPLICIT none


!.....................................................................
!!$Gesamtsteifigkeitsmatrix
  COMPLEX (KIND(0D0)),DIMENSION(*):: a_scal
!!$ Berechnete Potentialwerte (bzw. Loesungsverktor)
  COMPLEX (KIND(0D0)),DIMENSION(*):: b_scal
!!$ Skalirerungsfaktor
  REAL (KIND(0D0)),DIMENSION(*)  :: fak_scal

!     Hilfsvariablen
  INTEGER (KIND=4) ::    idi,i0
  INTEGER (KIND=4) ::     ja
  REAL(KIND(0D0))  ::     dum

!     Indexvariablen
  INTEGER (KIND=4) ::     i,j

!.....................................................................

  do  i=1,sanz

     idi = i*(mb+1)
     dum = cdabs(a_scal(idi))

     if (dum.le.0d0) then
        WRITE (fetxt,*)'scalab idi',idi,'i',i
        errnr = 27
        goto 1000
     end if

     a_scal(idi) = a_scal(idi) / dcmplx(dum)
     fak_scal(i) = 1d0 / dsqrt(dum)
     b_scal(i)   = b_scal(i) * dcmplx(fak_scal(i))

     if (i.eq.1) CYCLE

     i0 = i*mb
     ja = max0(1,i-mb)

     do j=ja,i-1
        a_scal(i0+j) = a_scal(i0+j) * dcmplx(fak_scal(i)*fak_scal(j))
     END do

  END do

  errnr = 0
  return

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     Fehlermeldungen

1000 return

end subroutine scalab

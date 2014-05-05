Minimum Gradient Support
========================

Available Versions
------------------

* (=5)	- MGS regularization (pure MGS after Zhdanov & Portniaguine):
  :math:`C_m^-1 \approx \int
  \frac{(\nabla m_{ij})^2}{(\nabla m_{ij})^2+\beta^2}\;dA`

* (=6)	- MGS with sensitivity weighted beta: = beta / f_{ij} (from Blaschek
  2008) where f_{ij} = (1 + g(i) + g(k))^2 and g(i,k) = log_{10} coverage
  (m_{i,k})

* (=7)	- MGS with sensitivity weighted beta (as with 6) but normalized: f_{ij}
  = ((1 + g(i) + g(k))/mean{coverage})^2

* (=8)	- MGS as in 6, but here a modified version of Roland Blaschek (which I
  apparently didn't understood, but was given to me from Fred...).  For more
  details please have a look into bsmatm_mod.f90 -> SUBROUTINE bsmamtmmgs

* (=9)- Same as 8, but with mean coverage normalization


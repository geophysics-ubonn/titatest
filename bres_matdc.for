      SUBROUTINE bres_matdc
      
      USE alloci
      USE femmod
      
      IMPLICIT none
      
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'


!     help variable  
      INTEGER :: i,j,k,l

      
      IF (.NOT.ALLOCATED(atadc)) ALLOCATE (atadc(manz,manz))
      
      atadc=MATMUL(TRANSPOSE(sensdc),sensdc)
      
      
      IF (.NOT.ALLOCATED(inv_atadc)) 
     1     ALLOCATE (inv_atadc(manz,manz))
      
      END SUBROUTINE bres_matdc

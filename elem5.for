      subroutine elem5()

c     Unterprogramm liefert die Elementmatrizen 'elmas(4,4)' und 'elmam(4,4)'
c     sowie den Elementvektor 'elve(4)' fuer ein Parallelogrammelement mit
c     bilinearem Ansatz ( Elementtyp Nr.5 ).

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   10-Nov-1997
      
c.....................................................................

      USE elemmod,only:xk,yk,elmam,elmas,elve
      USE errmod

      IMPLICIT none


c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Grundelementmatrizen
      integer         * 4     s1(4,4),s2(4,4),s3(4,4),
     1     s4(4,4)

c     Grundelementvektor
      integer         * 4     sb(4)

c     Hilfsvariablen
      real            * 8     x21,x31,y21,y31,
     1     det,a,b,c

c     Indexvariablen
      integer         * 4     i,j

c.....................................................................

      data s1/2,-2,-1,1,-2,2,1,-1,-1,1,2,-2,1,-1,-2,2/
      data s2/1,0,-1,0,0,-1,0,1,-1,0,1,0,0,1,0,-1/
      data s3/2,1,-1,-2,1,2,-2,-1,-1,-2,2,1,-2,-1,1,2/
      data s4/4,2,1,2,2,4,2,1,1,2,4,2,2,1,2,4/
      data sb/1,1,1,1/

      x21 = xk(2) - xk(1)
      x31 = xk(4) - xk(1)
      y21 = yk(2) - yk(1)
      y31 = yk(4) - yk(1)
      det = x21*y31 - x31*y21

c     Ggf. Fehlermeldung
      if (det.le.0d0) then
         fetxt = ' '
         errnr = 26
         return
      end if

      a =   (x31*x31 + y31*y31) / det
      b = - (x31*x21 + y31*y21) / det
      c =   (x21*x21 + y21*y21) / det

      do i=1,4
         elve(i) = det * dble(sb(i)) / 4d0

         do j=1,4
            elmas(i,j) =     (a*dble(s1(i,j)) +
     1           3d0*b*dble(s2(i,j)) +
     1           c*dble(s3(i,j))) / 6d0
            elmam(i,j) = det * dble(s4(i,j)) / 3.6d1
         end do
      end do

      errnr = 0

      return
      end

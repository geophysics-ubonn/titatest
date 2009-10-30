      subroutine elem3(xk,yk,elmas,elmam,elve)

c     Unterprogramm liefert die Elementmatrizen 'elmas(3,3)' und 'elmam(3,3)'
c     sowie den Elementvektor 'elve(3)' fuer ein Dreieckelement mit linearem
c     Ansatz ( Elementtyp Nr.3 ).

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   05-Nov-1997

c.....................................................................

      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      
c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     x-Koordinaten der Eckknotenpunkte
      real            * 8     xk(selmax)

c     y-Koordinaten der Eckknotenpunkte
      real            * 8     yk(selmax)

c     Elementmatrizen
      real            * 8     elmas(selmax,selmax),
     1     elmam(selmax,selmax)

c     Elementvektor
      real            * 8     elve(selmax)

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Grundelementmatrizen
      integer         * 4     s1(3,3),s2(3,3),s3(3,3),
     1     s4(3,3)

c     Grundelementvektor
      integer         * 4     sb(3)

c     Hilfsvariablen
      real            * 8     x21,x31,y21,y31,
     1     det,a,b,c

c     Indexvariablen
      integer         * 4     i,j

c.....................................................................

      data s1/1,-1,0,-1,1,0,0,0,0/
      data s2/2,-1,-1,-1,0,1,-1,1,0/
      data s3/1,0,-1,0,0,0,-1,0,1/
      data s4/2,1,1,1,2,1,1,1,2/
      data sb/1,1,1/

      x21 = xk(2) - xk(1)
      x31 = xk(3) - xk(1)
      y21 = yk(2) - yk(1)
      y31 = yk(3) - yk(1)
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

      do i=1,3
         elve(i) = det * dble(sb(i)) / 6d0

         do j=1,3
            elmas(i,j) = (a*dble(s1(i,j)) +
     1           b*dble(s2(i,j)) +
     1           c*dble(s3(i,j))) / 2d0
            elmam(i,j) = det * dble(s4(i,j)) / 2.4d1
         end do
      end do

      errnr = 0

      return
      end

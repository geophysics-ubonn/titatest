      subroutine elem4()

c     Unterprogramm liefert die Elementmatrizen 'elmas(6,6)' und 'elmam(6,6)'
c     sowie den Elementvektor 'elve(6)' fuer ein Dreieckelement mit quadra-
c     tischem Ansatz ( Elementtyp Nr.4 ).

c     ( Vgl. Subroutine 'DRQELL' in Schwarz (1991) )

c     Andreas Kemna                                            11-Feb-1993
c     Letzte Aenderung   25-Jul-2003

c.....................................................................


      USE elemmod,only:xk,yk,elmam,elmas,elve
      USE errmod

      IMPLICIT none
      
      
c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Grundelementmatrizen
      integer         * 4     s1(6,6),s2(6,6),s3(6,6),
     1     s4(6,6)

c     Grundelementvektor
      integer         * 2     sb(6)

c     Hilfsvariablen
      real            * 8     x21,x31,y21,y31,
     1     det,a,b,c

c     Indexvariablen
      integer         * 4     i,j

c.....................................................................

      data s1/3,1,0,-4,0,0,1,3,0,-4,0,0,6*0,
     1     -4,-4,0,8,0,0,4*0,8,-8,4*0,-8,8/

      data s2/6,1,1,-4,0,-4,1,0,-1,-4,4,0,1,-1,0,0,4,-4,
     1     -4,-4,0,8,-8,8,0,4,4,-8,8,-8,-4,0,-4,8,-8,8/

      data s3/3,0,1,0,0,-4,6*0,1,0,3,0,0,-4,3*0,8,-8,0,
     1     3*0,-8,8,0,-4,0,-4,0,0,8/

      data s4/6,-1,-1,0,-4,0,-1,6,-1,0,0,-4,-1,-1,6,-4,0,0,
     1     0,0,-4,32,16,16,-4,0,0,16,32,16,0,-4,0,16,16,32/

      data sb/0,0,0,1,1,1/

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

      do i=1,6
         elve(i) = det * dble(sb(i)) / 6d0

         do j=1,6
            elmas(i,j) = (a*dble(s1(i,j)) +
     1           b*dble(s2(i,j)) +
     1           c*dble(s3(i,j))) / 6d0
            elmam(i,j) = det * dble(s4(i,j)) / 3.6d2
         end do
      end do

      errnr = 0

      return
      end

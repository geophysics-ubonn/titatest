      subroutine precal()
      
c     Unterprogramm zur Berechnung der Element- und Randelementbeitraege
c     sowie der Konfigurationsfaktoren zur Berechnung der gemischten
c     Randbedingung.

c     Andreas Kemna                                            21-Dec-1995
c     Letzte Aenderung   22-Sep-1998

c.....................................................................

      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'waven.fin'
      INCLUDE 'fem.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsfunction
      real            * 8     beta

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Aktuelle Randelementnummer
      integer         * 4     rel

c     Aktueller Elementtyp
      integer         * 4     ntyp

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c     x-Koordinaten der Eckknotenpunkte
      real            * 8     xk(selmax)

c     y-Koordinaten der Eckknotenpunkte
      real            * 8     yk(selmax)

c     Elementmatrizen
      real            * 8     elmas(selmax,selmax),
     1     elmam(selmax,selmax)

c     Elementvektor
      real            * 8     elve(selmax)

c     Indexvariablen
      integer         * 4     i,j,
     1     imn,m,n,
     1     l,k

c.....................................................................

      lbeta  = .false.
      lrandb = .false.
      iel    = 0

      do i=1,typanz
         ntyp = typ(i)
         nkel = selanz(i)

         do 10 j=1,nelanz(i)
            iel = iel + 1
            do m=1,nkel
               xk(m) = sx(snr(nrel(iel,m)))
               yk(m) = sy(snr(nrel(iel,m)))
            end do
            WRITE (fetxt,'(a,I7,2F10.2)')'Elementnr',iel
c     Randelement, linearer Ansatz
            if (ntyp.eq.11) then

c     Ggf. Fehlermeldung
               if (lrandb) then
                  fetxt = ' '
                  errnr = 101
                  goto 1000
               end if

               lbeta = .true.
               call elem1(xk,yk,elmam,elve)

c     Randelementbeitraege berechnen
               imn = 0
               rel = iel - elanz

c     Ggf. Fehlermeldung
               if (rel.le.0) then
                  fetxt = ' '
                  errnr = 36
                  goto 1000
               end if

               do m=1,nkel
                  do n=1,m
                     imn = imn + 1
                     relbg(rel,imn) = elmam(m,n)
                  end do
               end do

c     Konfigurationsfaktoren zur Berechnung der gemischten Randbedingung
c     berechnen
               do l=1,eanz
                  do k=1,kwnanz
                     kg(rel,l,k) = beta(l,k,xk,yk)
                     if (errnr.ne.0) goto 1000
                  end do
               end do

            else if (ntyp.eq.12) then
               goto 10

            else if (ntyp.eq.13) then

c     Ggf. Fehlermeldung
               if (lbeta) then
                  fetxt = ' '
                  errnr = 101
                  goto 1000
               end if

               lrandb = .true.
               goto 10
               
c     Zusammengesetztes Viereckelement
c     (vier Teildreiecke mit linearem Ansatz)
            else if (ntyp.eq.8) then

               do k=1,kwnanz
                  call elem8(xk,yk,elmas,elve,kwn(k))
                  if (errnr.ne.0) goto 1000

c     Elementbeitraege berechnen
                  imn = 0

                  do m=1,nkel
                     do n=1,m
                        imn = imn + 1
                        elbg(iel,imn,k) = elmas(m,n)
                     end do
                  end do
               end do

            else

c     Dreieckelement, linearer Ansatz
               if (ntyp.eq.3) then

                  call elem3(xk,yk,elmas,elmam,elve)
                  if (errnr.ne.0) goto 1000

c     Parallelogrammelement, bilinearer Ansatz
               else if (ntyp.eq.5) then

                  call elem5(xk,yk,elmas,elmam,elve)
                  if (errnr.ne.0) goto 1000

c     Fehlermeldung
               else
                  fetxt = ' '
                  errnr = 18
                  goto 1000
               end if

c     Elementbeitraege berechnen
               do k=1,kwnanz
                  imn = 0

                  do m=1,nkel
                     do n=1,m
                        imn = imn + 1
                        elbg(iel,imn,k) = elmas(m,n) +
     1                       elmam(m,n)*kwn(k)*kwn(k)
                     end do
                  end do
               end do

            end if
 10      continue
      end do

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

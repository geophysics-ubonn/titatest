      subroutine precal()
      
!!!$     Unterprogramm zur Berechnung der Element- und Randelementbeitraege
!!!$     sowie der Konfigurationsfaktoren zur Berechnung der gemischten
!!!$     Randbedingung.

!!!$     Andreas Kemna                                            21-Dec-1995
!!!$     Letzte Aenderung   22-Sep-1998

!!!$.....................................................................

      USE femmod
      USE electrmod
      USE elemmod
      USE wavenmod
      USE errmod
      USE konvmod, ONLY : lsytop
      IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfunction
      real            * 8     beta

!!!$     Aktuelle Elementnummer
      integer         * 4     iel

!!!$     Aktuelle Randelementnummer
      integer         * 4     rel

!!!$     Aktueller Elementtyp
      integer         * 4     ntyp

!!!$     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel


!!!$     Indexvariablen
      integer         * 4     i,j,imn,m,n,l,k

!!!$.....................................................................

      lbeta  = .false.
      lrandb = .false.
      iel    = 0
      
      ALLOCATE (xk(smaxs),yk(smaxs),elmas(smaxs,smaxs),
     1     elmam(smaxs,smaxs),elve(smaxs),stat=errnr)
      IF (errnr /= 0) then
         fetxt = 'allocation problem elbg elmam'
         errnr = 97 
         GOTO 1000
      END IF
      ALLOCATE (elbg(elanz,(smaxs*(smaxs+1))/2,kwnanz),stat=errnr)
      IF (errnr /= 0) then
         fetxt = 'allocation problem elbg elbg'
         errnr = 97 
         GOTO 1000
      END IF
      ALLOCATE (relbg(relanz,(smaxs*(smaxs+1))/2),stat=errnr)
      ALLOCATE (kg(relanz,eanz,kwnanz),stat=errnr)
      IF (errnr /= 0) then
         fetxt = 'allocation problem elbg relbg or kg'
         errnr = 97 
         GOTO 1000
      END IF

      IF (lsytop) THEN
         CALL bsytop
      ELSE
         sytop = 0.
      END IF

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
!!!$     Randelement, linearer Ansatz
            if (ntyp.eq.11) then

!!!$     Ggf. Fehlermeldung
               if (lrandb) then
                  fetxt = ' '
                  errnr = 101
                  goto 1000
               end if

               lbeta = .true.
               call elem1

!!!$     Randelementbeitraege berechnen
               imn = 0
               rel = iel - elanz

!!!$     Ggf. Fehlermeldung
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

!!!$     Konfigurationsfaktoren zur Berechnung der gemischten Randbedingung
!!!$     berechnen
               do l=1,eanz
                  do k=1,kwnanz
                     kg(rel,l,k) = beta(l,k)
                     if (errnr.ne.0) goto 1000
                  end do
               end do

            else if (ntyp.eq.12) then
               goto 10

            else if (ntyp.eq.13) then

!!!$     Ggf. Fehlermeldung
               if (lbeta) then
                  fetxt = ' '
                  errnr = 101
                  goto 1000
               end if

               lrandb = .true.
               goto 10
               
!!!$     Zusammengesetztes Viereckelement
!!!$     (vier Teildreiecke mit linearem Ansatz)
            else if (ntyp.eq.8) then

               do k=1,kwnanz
                  call elem8(elmas,elve,kwn(k),smaxs)
                  if (errnr.ne.0) goto 1000

!!!$     Elementbeitraege berechnen
                  imn = 0

                  do m=1,nkel
                     do n=1,m
                        imn = imn + 1
                        elbg(iel,imn,k) = elmas(m,n)
                     end do
                  end do
               end do

            else

!!!$     Dreieckelement, linearer Ansatz
               if (ntyp.eq.3) then

                  call elem3
                  if (errnr.ne.0) goto 1000

!!!$     Parallelogrammelement, bilinearer Ansatz
               else if (ntyp.eq.5) then

                  call elem5
                  if (errnr.ne.0) goto 1000

!!!$     Fehlermeldung
               else
                  fetxt = ' '
                  errnr = 18
                  goto 1000
               end if

!!!$     Elementbeitraege berechnen
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

      IF (ALLOCATED (xk)) DEALLOCATE (xk,yk,elmas,elmam,elve)

      errnr = 0

      return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 1000 return

      end

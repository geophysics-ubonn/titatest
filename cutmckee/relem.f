      subroutine relem(kanal,io,datei)

c     Unterprogramm zum Einlesen der FEM-Parameter aus 'datei'.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   24-Oct-1996

c.....................................................................

      USE elemmod

      IMPLICIT none

      INCLUDE 'err.fin'
      
c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal,io

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariablen
      integer         * 4     i,j,k

c     Hilfsvariable
      integer         * 4     idum

c.....................................................................

c     'datei' oeffnen
      fetxt = datei

      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)
      OPEN(io,file='ctmck.info',STATUS='replace')
      errnr = 3

c     Anzahl der Knoten (bzw. Knotenvariablen), Anzahl der Elementtypen
c     sowie Bandbreite der Gesamtsteifigkeitsmatrix einlesen
      read(kanal,*,end=1001,err=1000) sanz,typanz,mb
!     !$ now get some memory for the fields..
!     !$ first the sanz fields
      ALLOCATE (sx(sanz),sy(sanz),snr(sanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation sx failed'
         errnr = 94
         GOTO 999
      END IF

      ALLOCATE (typ(typanz),nelanz(typanz),selanz(typanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation selanz failed'
         errnr = 94
         GOTO 999
      END IF

c     Elementtypen, Anzahl der Elemente eines bestimmten Typs sowie
c     Anzahl der Knoten in einem Elementtyp einlesen
      read(kanal,*,end=1001,err=1000)
     1     (typ(i),nelanz(i),selanz(i),i=1,typanz)

!     !$ set number of node points for regular elements
      smaxs = MAXVAL(selanz)

c     Anzahl der Elemente (ohne Randelemente) und Anzahl der Randelemente
c     bestimmen
      relanz = 0
      elanz  = 0
      do i=1,typanz
         if (typ(i).gt.10) then
            relanz = relanz + nelanz(i)
         else
            elanz  = elanz  + nelanz(i)
         end if
      end do
!     !$ get memory for the element integer field      
      ALLOCATE (nrel(elanz+relanz,smaxs),rnr(relanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation nrel failed'
         errnr = 94
         GOTO 999
      END IF

      PRINT*,'check dis out ',4+sanz+elanz+relanz*2

c     Zeiger auf Koordinaten, x-Koordinaten sowie y-Koordinaten der Knoten
c     einlesen
      print*,'knoten::',sanz
      read(kanal,*,end=1001,err=1000) (snr(i),sx(i),sy(i),i=1,sanz)
      i = 0;j = 0
      DO k=1,sanz
         IF (sx(snr(k))==sx(snr(1))) j = j+1 !counts vertical grid nodes
         IF (sy(snr(k))==sy(snr(1))) i = i+1 !counts horizontal grid nodes
      END DO
      IF ((i-1) == 0 .OR. (j-1) == 0) THEN
         PRINT*,ACHAR(9)//'## GRID not regular ##'
         WRITE (io,'(a)')'## Grid not regular'
      ELSE
         WRITE (*,*)'Regular grid info NX=',i-1,
     1        '  NY=',j-1,'  NXY=',(i-1)*(j-1)
         WRITE (*,*)'x0=',MINVAL(sx(1:sanz)),'y0=',
     1        MINVAL(ABS(sy(1:sanz)))
         WRITE (*,'(2(a,2x,F10.4))')'Linear dx',sx((2))-sx((1)),
     1        'dy',sy((i+1))-sy((1))
         WRITE (io,'(a)')'## Regular grid'
         WRITE (io,'(A,2X,I4)')'NX',i-1
         WRITE (io,'(a,2x,F10.3)')'X0',MINVAL(sx(1:sanz))         
         WRITE (io,'(a,2x,F10.3)')'DX',sx((2))-sx((1))
         WRITE (io,'(A,2X,I4)')'NY',j-1
         WRITE (io,'(a,2x,F10.3)')'Y0',MINVAL(ABS(sy(1:sanz)))
         WRITE (io,'(a,2x,F10.3)')'DY',sy((i+1))-sy((1))
         WRITE (io,'(A,2X,I4)')'NXY',(i-1)*(j-1)
      END IF
c     Knotennummern der Elemente einlesen
      print*,'post knoten position',4+sanz
      idum = 0
      print*,'elemente::',elanz,relanz
      do i=1,typanz
         do j=1,nelanz(i)
            read(kanal,*,end=1001,err=1000)
     1           (nrel(idum+j,k),k=1,selanz(i))
         end do
         idum = idum + nelanz(i)
      end do
      print*,'post elemente position',4+sanz+elanz+relanz
      print*,'randpunkte::',relanz
c     Zeiger auf Werte der Randelemente einlesen
      read(kanal,*,end=1001,err=1000) (rnr(i),i=1,relanz)
      print*,'post randpunkte position',4+sanz+elanz+relanz*2
c     'datei' schliessen
      close(kanal)
      CLOSE(io)
      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

 1001 close(kanal)
      errnr = 2
      return

      end

      subroutine welem(kanal,datei)

c     Unterprogramm zum Schreiben der FEM-Parameter in 'datei'.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   21-Jan-2003

c.....................................................................

      USE elemmod

      IMPLICIT none

      INCLUDE 'err.fin'
c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

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
      open(kanal,file=fetxt,status='replace',err=999)

      errnr = 4

c     Anzahl der Knoten (bzw. Knotenvariablen), Anzahl der Elementtypen
c     sowie Bandbreite der Gesamtsteifigkeitsmatrix schreiben
      write(kanal,*,err=1000) sanz,typanz,mb

c     Elementtypen, Anzahl der Elemente eines bestimmten Typs sowie
c     Anzahl der Knoten in einem Elementtyp schreiben
      do i=1,typanz
         write(kanal,*,err=1000) typ(i),nelanz(i),selanz(i)
      end do

c     Zeiger auf Koordinaten, x-Koordinaten sowie y-Koordinaten der Knoten
c     schreiben
      do i=1,sanz
         write(kanal,'(I8,2F12.3)',err=1000) 
     1        snr(i),real(sx(i)),real(sy(i))
      end do

c     Knotennummern der Elemente schreiben
      idum = 0
      do i=1,typanz
         do j=1,nelanz(i)
            write(kanal,*,err=1000) (nrel(idum+j,k),k=1,selanz(i))
         end do
         idum = idum + nelanz(i)
      end do

c     Zeiger auf Werte der Randelemente schreiben
      do i=1,relanz
         write(kanal,*,err=1000) rnr(i)
      end do

c     'datei' schliessen
      close(kanal)
      
      fetxt='ctm_triangles.txt'
      OPEN (kanal,FILE=fetxt,STATUS='replace',ACCESS='sequential')
      do j=1,nelanz(1)
         write(kanal,*,err=1000) (nrel(j,k),k=1,selanz(1))
      end do
      CLOSE (kanal)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end

      program fem

c     Hauptprogramm zum Complex-Resistivity-2.5D-Modelling.

c     Belegte Kanaele:  9 - error.dat
c     11 - in-/output
c     12 - crmod.cfg

c     Andreas Kemna                                            15-Dec-1994
c     Letzte Aenderung   02-May-2008

c.....................................................................

      USE alloci
      USE tic_toc
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'waven.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'fem.fin'
      INCLUDE 'randb.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     Kanalnummer
      integer         * 4     kanal


c     Dateinamen
      character       * 80    delem,
     1     delectr,
     1     dsigma,
     1     dstrom,
     1     dkpot,
     1     dpot,
     1     dvolt,
     1     dsens,
     1     drandb

c     Schalter ob transformierte Potentialwerte ausgegeben werden sollen
      logical         * 4     lkpot

c     Schalter ob Potentialwerte ausgegeben werden sollen
      logical         * 4     lpot

c     Schalter ob Spannungswerte ausgegeben werden sollen
      logical         * 4     lvolt

c     Schalter ob weiterer Datensatz modelliert werden soll
      logical         * 4     lagain

c     Schalter ob mit K-Faktor modelliert werden soll (default ohne)
      logical         * 4     wkfak

c     Schalter ob nur analytisch modelliert werden soll (default ohne)
      logical         * 4     lana

c     Indexvariablen
      integer         * 4     j,k,l,c1

      character       *256    ftext
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      CALL tic(c1)
c     'crmod.cfg' oeffnen
      fetxt = 'crmod.cfg'
      errnr = 1
      open(12,file=fetxt,status='old',err=999)

c     Allgemeine Parameter setzen
      kanal = 11

c     "singularity removal" ?
c     ak        lsr = .true.
      lsr = .false.

c     Art der Ruecktransformation
c     ak        swrtr = 1
c     ak        swrtr = 0

c     Sonstiges
      lkpot = .false.
      wkfak = .false.
      lana = .FALSE.
      lprior = .false. ! kann prior modell verrauschen
      lnsepri = .false.
c     ak        lkpot = .true.
c     ak        dkpot = '..\tmp\kpot.ctr'

c     Kontrollausgabe
 5    write(*,*)
      write(*,'(a20)') ' Reading Input-Files'

c     'crmod.cfg' einlesen
      fetxt = 'crmod.cfg'
      errnr = 3
      mswitch = 0
      read(12,*,end=1001,err=999)
      read(12,'(a80)',end=1001,err=999) delem
      read(12,'(a80)',end=1001,err=999) delectr
      read(12,'(a80)',end=1001,err=999) dsigma
      read(12,'(a80)',end=1001,err=999) dstrom
      read(12,*,end=1001,err=999) lpot
      read(12,'(a80)',end=1001,err=999) dpot
      read(12,*,end=1001,err=999) lvolt
      read(12,'(a80)',end=1001,err=999) dvolt
      read(12,*,end=1001,err=999) lsens
      read(12,'(a80)',end=1001,err=999) dsens
      read(12,*,end=1001,err=999) lagain
      read(12,*,end=1001,err=999) swrtr
      read(12,*,end=1001,err=999) lsink
      read(12,*,end=1001,err=999) nsink
      read(12,*,end=1001,err=999) lrandb2
      read(12,'(a80)',end=1001,err=999) drandb
      fetxt = 'reading mswitch'
      read(12,'(I2)',end=101,err=100) mswitch

      GOTO 101

 100  BACKSPACE (12)

 101  lana = BTEST(mswitch,0)
      wkfak = BTEST(mswitch,1)
      lsr = BTEST(mswitch,2)
      


c     Alles einlesen
      call relem(kanal,delem)
      if (errnr.ne.0) goto 999

      call relectr(kanal,delectr)
      if (errnr.ne.0) goto 999

      call rdatm(kanal,dstrom)
      if (errnr.ne.0) goto 999

      if (swrtr.eq.0) then
         lana = .FALSE.
         lsr = .FALSE.
         kwnanz = 1
         kwn(1) = 0d0
      else
         call rwaven()
         if (errnr.ne.0) goto 999
      end if
      print*,''
      IF (wkfak) PRINT*,'Modelling with K-Faktor!'
      IF (lana) PRINT*,'Analytical solution only'
      IF (lsr) THEN
         PRINT*,'Singularity removal'
         lana = .false.
      END IF
      lsr = lana
      print*,''

      call rsigma(kanal,dsigma)
      if (errnr.ne.0) goto 999

      if (lrandb2) then
         call rrandb(kanal,drandb)
         if (errnr.ne.0) goto 999
      end if

c     Ggf. Referenzleitfaehigkeit bestimmen
      if (lsr) call refsig()
      
      IF (lsink) WRITE(6,'(/A,I5,2F12.3/)')
     1     'Fictious sink @ node ',nsink,sx(snr(nsink)),sy(snr(nsink))
c     Element- und Randelementbeitraege sowie Konfigurationsfaktoren zur
c     Berechnung der gemischten Randbedingung bestimmen
      call precal()
      if (errnr.ne.0) goto 999

      if (.not.lbeta) lsr=.false.

c     'a', 'hpot' und 'kpot' zuweisen
      ALLOCATE(a((mb+1)*sanz),hpot(sanz,eanz),
     1     kpot(sanz,eanz,kwnanz),stat=errnr)
      if (errnr.ne.0) then
         errnr = 97 
         goto 999
      end if

c     POTENTIALWERTE BERECHNEN
      do k=1,kwnanz
c     Kontrollausgabe
         if (swrtr.eq.0) then
            write(*,'(a)') ' Calculating Potentials'
         else
            write(*,'(a,i5)',ADVANCE='no')ACHAR(13)//
     1           ' Calculating Potentials : Wavenumber ',k
         end if

         do l=1,eanz
            if (lsr.or.lbeta.or.l.eq.1) then
               
c     Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
               if (lsr) call potana(l,k)
               
c     Kompilation des Gleichungssystems (fuer Einheitsstrom !)
               call kompab(l,k)
               if (errnr.ne.0) goto 999
                  
c     Ggf. Randbedingung beruecksichtigen
               if (lrandb) call randb()
               if (lrandb2) call randb2()
               
c     Gleichungssystem skalieren
               call scalab()
               if (errnr.ne.0) goto 999
               
c     Cholesky-Zerlegung der Matrix
               call chol()
               if (errnr.ne.0) goto 999
            else
               
c     Stromvektor modifizieren
               call kompb(l)
            end if
            
c     Gleichungssystem loesen
            IF (.NOT.lana) call vre()
c     Potentialwerte zurueckskalieren und umspeichern sowie ggf.
c     analytische Loesung addieren
            do j=1,sanz
               IF (lana) THEN
                  kpot(j,l,k) = pota(j)
               ELSE
                  kpot(j,l,k) = pot(j) * dcmplx(fak(j))
                  if (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
               END IF
c     ak (fuer Testzwecke)
c     ak                    kpot(j,l,k) = pota(j)
               if (swrtr.eq.0) hpot(j,l) = kpot(j,l,k)
            end do
         end do
      end do
      print*,''
c     Ggf. Ruecktransformation der Potentialwerte
      if (swrtr.eq.1) call rtrafo()

c     Ggf. transformierte Potentialwerte ausgeben
      if (lkpot) then
         call wkpot(kanal,dkpot)
         if (errnr.ne.0) goto 999
      end if

c     Ggf. Potentialwerte berechnen und ausgeben
      if (lpot) then
         call bpot(kanal,dpot)
         if (errnr.ne.0) goto 999
      end if

c     Ggf. Spannungswerte berechnen und ausgeben
c     (bzw. scheinbaren Widerstandswerte)
      if (lvolt) then
         IF (wkfak.AND.lbeta) then
            call bkfak()
            if (errnr.ne.0) goto 999
         else
c     nur echte Spannungen ausgeben...
            kfak=1.0 
         end if

         call bvolt()
         if (errnr.ne.0) goto 999

         call wdatm(kanal,dvolt)
         if (errnr.ne.0) goto 999
      end if

c     'a' und 'hpot' freigeben
      DEALLOCATE(a,hpot)

c     Ggf. Sensitivitaeten aller Messungen berechnen und ausgeben
      if (lsens) then

c     Modelleinteilung gemaess Elementeinteilung belegen
         manz = elanz

         if (manz.gt.mmax) then
            fetxt = ' '
            errnr = 63
            goto 999
         end if

         do j=1,elanz
            mnr(j) = j
         end do

c     'sens' zuweisen
         ALLOCATE(sens(nanz,manz),stat=errnr)
         if (errnr.ne.0) then
            errnr = 97 
            goto 999
         end if

         call bsens()

         call wsens(kanal,dsens)
         if (errnr.ne.0) goto 999

c     'sens' freigeben
         DEALLOCATE(sens)
      end if

c     'kpot' freigeben
      DEALLOCATE(kpot)

c     Ggf. weiteren Datensatz modellieren
      if (lagain) goto 5

c     'crmod.cfg' schliessen
      close(12)

c     Kontrollausgabe

      fetxt = 'solution time'
      CALL TOC(c1,fetxt)

      write(*,*)
      write(*,'(a)',ADVANCE='no')' Modelling completed'

      STOP '0'
      
c.....................................................................

c     (Fehler-) Meldung schreiben
 999  open(9,file='error.dat',status='replace')
      errflag = 2
      CALL get_error(ftext,errnr,errflag,fetxt)
      write(9,'(a80,i3,i1)') fetxt,errnr,errflag
      write(9,*)ftext
      close(9)
      stop '-1'

 1001 open(9,file='error.dat',status='replace')
      errnr   = 2
      errflag = 2
      CALL get_error(ftext,errnr,errflag,fetxt)
      write(9,'(a80,i3,i1)') fetxt,errnr,errflag
      write(9,*)ftext
      close(9)
      stop '-2'

      end

      program cutmck

c     Programm zur optimalen Numerierung der Knotenpukte nach dem
c     Algorithmus von Cuthill-McKee. Ist die Anzahl der Startpunkte 'spanz'
c     kleiner Null, werden die Startpunkte mit minimalem Grad automatisch
c     bestimmt. Andernfalls werden die vorgegebenen Nummern der Startpunkte
c     aus 'dstart' gelesen.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung, RM   11-Sep-2009

c     ( Vgl. entspr. Hauptprogramm in Schwarz (1991) )

c.....................................................................

      USE elemmod
      USE electrmod

      IMPLICIT none

      INCLUDE 'err.fin'

c.....................................................................

c     Schalter ob Kontrolldateien ('*.ctr') ausgegeben werden sollen
      logical         * 1     kont

c     Dateinamen
      character       * 80    delem,delectr
      CHARACTER(256)      ::  ftext
c     Maximaler Grad der Knotenpunkte
      INTEGER,PARAMETER   :: grmax=500
c     Maximale Anzahl vorgegebener Startpunkte im Cuthill-McKee-Algorithmus
      INTEGER,PARAMETER   :: spmax=100
c     Maximale Anzahl der Stufen im Cuthill-McKee-Algorithmus
      INTEGER,PARAMETER   :: stmax=100

c     Permutationsvektor der Umnumerierung
      INTEGER,DIMENSION(:),ALLOCATABLE   :: perm

c     Knotennummern der Startpunkte
      INTEGER,DIMENSION(:),ALLOCATABLE   :: start

c     Hilfsvariablen
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: graph
      INTEGER,DIMENSION(:),ALLOCATABLE   :: grad,neu,neuin,level

      integer         * 4   gradzp,fcm,kbdm,nstart,spanz,
     1     i,j,k,l,m,idum,is,nzp,nnp,maxgd,mingd,
     1     minbd,mmin,mingr,levs,leve,nlev

      LOGICAL,DIMENSION(:),ALLOCATABLE   :: num
      logical         * 1     exi1,exi2,exi3
      integer         * 4     c1,c2,se,mi,st,ta,fp,fp2
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 4    FORMAT(/'  minimum grade =',I6/'  maximum grade =',I6,
     1     '  minimum bandwidth =',I6)
 6    FORMAT(/'Ergebnisse der Nuenummerierungen '/,'Startpunkt',3x,
     1     'Bandbreite',3x,'Profil (normal)',3x,'Profil (reverse)')
 8    FORMAT(3(I8,4X),I8)
 9    FORMAT(//'  minimum bandwidth =',i5,'  for init.node',
     1     i5//'  vector of permutation ='/)
 11   FORMAT((3x,10i5))
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CALL SYSTEM_CLOCK (c1,i)
      kont   = .FALSE.
      
      fp = 11
      fp2 = 12
c     Fehlerdatei oeffnen
      open(10,file='error.dat',status='unknown')

c     Fehlerflag setzen
      errflag = 2

c     'CutMck.cfg' einlesen
      fetxt = 'cutmck.cfg'
      INQUIRE (FILE=fetxt,EXIST=exi1)
      IF (exi1) THEN
         errnr = 1
         open(fp,file=fetxt,status='old',err=1000)
         
         errnr = 3
         read(fp,'(a,/,a)',end=1001,err=1000)
     1        delem,delectr
         close(fp)
      ELSE                      !check for defaults
         delem='elem.dat'
         delectr='elec.dat'
         INQUIRE (FILE=delem,EXIST=exi2)
         INQUIRE (FILE=delectr,EXIST=exi3)
         IF (exi2.AND.exi3) THEN
            PRINT*,'trying default::',TRIM(delem),'  ',TRIM(delectr)
         ELSE
            IF (.NOT.exi2) THEN
               fetxt=TRIM(delem)
               PRINT*,'no default grid file::',TRIM(delem)
            END IF
            IF (.NOT.exi3)THEN
               fetxt=TRIM(fetxt)//' '//TRIM(delectr)
               PRINT*,'no default electrode file::',TRIM(delectr)
            END IF
            errnr = 3
            GOTO 1000
         END IF
      END IF
      print*,'read in elements...'
c     'delem' einlesen
      call relem(fp,fp2,delem)
      if (errnr.ne.0) goto 1000
      print*,'read in electrodes'
c     'delectr' einlesen
      call relectr(fp,delectr)
      if (errnr.ne.0) goto 1000
      print*,'ok'

      ALLOCATE (graph(grmax,sanz),perm(sanz),grad(sanz),neu(sanz),
     1     neuin(sanz),level(sanz),start(sanz),num(sanz),stat=errnr)
      

c     Aufbau des Graphen aufgrund der Knotennummern der Elemente
      grad=0;graph=0

      idum = 0
      grad=0.;graph=0.
      DO l = 1 , typanz

         DO m = 1 , nelanz(l)

            DO i = 1 , selanz(l)-1

               nzp = nrel(idum+m,i)

               DO j=i+1,selanz(l)

                  nnp = nrel(idum+m,j)

                  DO k=1,grmax

                     IF (graph(k,nzp) == nnp) GOTO 10
                     IF (graph(k,nzp) > 0) CYCLE

                     graph(k,nzp)=nnp
                     grad(nzp)=grad(nzp)+1

                     EXIT
                  END DO
                  
                  grad(nnp) = grad(nnp)+1
                  IF (grad(nnp) <= grmax) graph(grad(nnp),nnp)=nzp
                  
 10               CONTINUE

               END DO

            END DO

         END DO

         idum = idum + nelanz(l)

      END DO

      PRINT*,'graph erstellt'
      
      mingd = MINVAL(grad(1:sanz))
      maxgd = MAXVAL(grad(1:sanz))

      IF (mingd==0) THEN
         PRINT*,'existance of zero nodes.. aborting'
         STOP
      ENDIF

      minbd = (maxgd+1)/2
c     Ggf. Kontrolldatei oeffnen

      IF (kont) THEN
         fetxt = 'cutmck.ctr'
         OPEN(fp,file=TRIM(fetxt),status='unknown')
         WRITE(fp,4) mingd,maxgd,minbd
      END IF

c     rm
      spanz = INT(elanz / 10)
      spanz = MAX(spanz,1000)
      k = 0
      spanz = MIN(spanz,spmax)
      WRITE(*,'(a,I8)',ADVANCE='no')'mingd::',mingd
 110  do 120 i=1,sanz
         if (grad(i).eq.mingd) then
            k = k+1
            start(k) = i
            if (k.ge.spanz) goto 130
         end if
 120  continue
      IF (mingd>elanz/10) GOTO 130
      mingd = mingd+1
      WRITE(*,'(a,I8)',ADVANCE='no')ACHAR(13)//'mingd::',mingd
      goto 110

 130  continue

      spanz=k-1
      print*,'errechnete startpunkte ',spanz
      if (kont) then

         errnr = 4
         write(11,7,err=1000)
 7       format(//'  init.node   bandwidth'/)
      end if

c     Neunumerierung der Knotenpunkte fuer alle Startpunkte
      mmin = sanz
      kbdm = 0

      do is=1,spanz
         nstart        = start(is)
         
         WRITE (*,'(2(a,1X,I6,2X))',ADVANCE='no')
     1        ACHAR(13)//'startknoten',nstart,'bandbreite',mmin
         
         neu(1)        = nstart
         neuin(nstart) = 1
         
         do i=1,sanz
            num(i) = .false.
         end do
         
         num(nstart) = .true.
         level(1)    = 1
         levs        = 1
         leve        = 1
         nlev        = 1
         l           = 1

 150     do 180 j=levs,leve
            nzp    = neu(j)
            gradzp = grad(nzp)
 160        mingr  = grmax
            k      = 0

            do 170 i=1,gradzp
               nnp = graph(i,nzp)

               if (num(nnp).or.grad(nnp).gt.mingr) goto 170

               mingr = grad(nnp)
               k     = nnp
 170        continue

            if (k.eq.0) goto 180

            l        = l+1
            neu(l)   = k
            neuin(k) = l
            num(k)   = .true.

            goto 160
 180     continue
         levs        = levs+level(nlev)
         nlev        = nlev+1
         level(nlev) = l-levs+1
         leve        = leve+level(nlev)


         if (leve.lt.sanz) goto 150

c     Bandbreite 'mb'
         mb = 0

         do i=1,sanz
            nzp    = neuin(i)
            fcm    = nzp
            gradzp = grad(i)

            do j=1,gradzp
               k   = neuin(graph(j,i))
               mb  = max0(mb,iabs(k-nzp))
               fcm = min0(fcm,k)
            end do
         end do

         if (kont) then
            write(fp,8,err=1000) nstart,mb
         end if

         if (mb.lt.mmin) then
            
            mmin = mb
            kbdm = is
            do i=1,sanz
               perm(i) = neuin(i)
            end do
         end if
         
      end do                    ! is

      if (kont) then
         write(fp,9,err=1000) mmin,start(kbdm)
         write(fp,11,err=1000) (perm(i),i=1,sanz)
         close(fp)
      end if

c     Minimale Bandbreite setzen
      mb = mmin

c     UMNUMERIERUNG

c     Knotennummern der Elemente umspeichern
      idum = 0

      do i=1,typanz
         do j=1,nelanz(i)
            do k=1,selanz(i)
               nrel(idum+j,k) = perm(nrel(idum+j,k))
            end do
         end do
         idum = idum+nelanz(i)
      end do

c     Zeiger auf Koordinaten der Knoten umspeichern ('grad' als Hilfsfeld)
      do i=1,sanz
         grad(i) = snr(i)
      end do

      do i=1,sanz
         snr(perm(i)) = grad(i)
      end do

c     Knotennummern der Elektroden umspeichern
      do i=1,eanz
         enr(i) = perm(enr(i))
      end do

c     Startwerte umspeichern (zur Kontrolle)
      do i=1,spanz
         start(i) = perm(start(i))
      end do
      PRINT*,' writing out new values '
c     Elementeinteilung und Elektrodenverteilung schreiben
c     delem=TRIM(ADJUSTL(delem))//'_ctm'
      fetxt='cp -f '//TRIM(delem)//' '//TRIM(delem)//'.orig'
      CALL SYSTEM (fetxt)
      call welem(11,delem)
      if (errnr.ne.0) goto 1000

      fetxt='cp -f '//TRIM(delectr)//' '//TRIM(delectr)//'.orig'
      CALL SYSTEM (fetxt)
c     delectr=TRIM(ADJUSTL(delectr))//'_ctm'
      call welectr(11,delectr)
      if (errnr.ne.0) goto 1000

c     Fehlerdatei loeschen
      close(10,status='delete')

      CALL SYSTEM_CLOCK (c2,i)
      k=(c2-c1)/j               ! Sekunden
      mi=INT(k/60)              ! Minuten
      st=INT(k/60/60)           ! Stunden
      ta=INT(k/60/60/24)        ! Tage
      se=k-mi*60-st*60*60-ta*60*60*24 ! Sekunden
 3    FORMAT(I2,'d/',1X,I2,'h/',1X,I2,'m/',1X,I2,'s')
      WRITE (*,3)ta,st,mi,se

      RETURN

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 CALL get_error(ftext,errnr,errflag,fetxt)
      write(10,'(a80,i2,i1)') fetxt,errnr,errflag
      write(10,*)ftext
      close(10)
      stop 'abbruch '

 1001 errnr = 2
      CALL get_error(ftext,errnr,errflag,fetxt)
      write(10,'(a80,i2,i1)') fetxt,errnr,errflag
      write(10,*)ftext
      close(10)
      stop 'abbruch '

      end

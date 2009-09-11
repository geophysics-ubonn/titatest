      program cutmck

c     Programm zur optimalen Numerierung der Knotenpukte nach dem
c     Algorithmus von Cuthill-McKee. Ist die Anzahl der Startpunkte 'spanz'
c     kleiner Null, werden die Startpunkte mit minimalem Grad automatisch
c     bestimmt. Andernfalls werden die vorgegebenen Nummern der Startpunkte
c     aus 'dstart' gelesen.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   21-Jan-2003

c     ( Vgl. entspr. Hauptprogramm in Schwarz (1991) )

c.....................................................................

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'cutmck.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'

c.....................................................................

c     Schalter ob Kontrolldateien ('*.ctr') ausgegeben werden sollen
      logical         * 1     kont

c     Dateinamen
      character       * 80    delem,
     1     delectr
      CHARACTER(256)      ::  ftext
c     Permutationsvektor der Umnumerierung
      integer         * 4     perm(smax)

c     Hilfsvariablen
      integer         * 4     graph(grmax,smax),
     1     grad(smax),
     1     neu(smax),
     1     neuin(smax),
     1     level(stmax),
     1     gradzp,
     1     fcm,kbdm,nstart,
     1     i,j,k,l,m,
     1     idum,is,nzp,nnp,
     1     maxgd,mingd,
     1     minbd,mmin,mingr,
     1     levs,leve,nlev

      logical         * 1     num(smax)
      logical         * 1     exi1,exi2,exi3
      integer         * 4     c1,c2,se,mi,st,ta
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CALL SYSTEM_CLOCK (c1,i)
      kont   = .FALSE.

c     Fehlerdatei oeffnen
      open(10,file='error.dat',status='unknown')

c     Fehlerflag setzen
      errflag = 2

c     'CutMck.cfg' einlesen
      fetxt = 'cutmck.cfg'
      INQUIRE (FILE=fetxt,EXIST=exi1)
      IF (exi1) THEN
         errnr = 1
         open(11,file=fetxt,status='old',err=1000)
         
         errnr = 3
         read(11,'(a,/,a)',end=1001,err=1000)
     1        delem,delectr
         close(11)
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
      call relem(11,delem)
      if (errnr.ne.0) goto 1000
      print*,'read in electrodes'
c     'delectr' einlesen
      call relectr(11,delectr)
      if (errnr.ne.0) goto 1000
      print*,'ok'
c     Aufbau des Graphen aufgrund der Knotennummern der Elemente
      grad=0;graph=0

      idum = 0

      do l=1,typanz
         do m=1,nelanz(l)

            do 90 i=1,selanz(l)-1
               nzp = nrel(idum+m,i)

               do 80 j=i+1,selanz(l)
                  nnp = nrel(idum+m,j)
                  
                  do 50 k=1,grmax
                     if (graph(k,nzp).eq.nnp) goto 80
                     if (graph(k,nzp).gt.0) goto 50

                     graph(k,nzp) = nnp
                     grad(nzp)    = grad(nzp)+1

                     goto 60
 50               continue

c     Fehlermeldung
                  fetxt = ' '
                  errnr = 22
                  goto 1000

 60               grad(nnp) = grad(nnp)+1
                  if (grad(nnp).le.grmax) goto 70

c     Fehlermeldung
                  fetxt = ' '
                  errnr = 22
                  goto 1000

 70               graph(grad(nnp),nnp) = nzp

 80            continue
 90         continue

         end do
         idum = idum+nelanz(l)
      end do
      PRINT*,'graph erstellt'
      maxgd = grad(1)
      mingd = grad(1)
      do i=2,sanz
         maxgd = max0(maxgd,grad(i))
         mingd = min0(mingd,grad(i))
      end do
c$$$  DO i=1,sanz
c$$$  PRINT*,'Knoten ',i,'grad ',grad(i)
c$$$  ENDDO

      IF (mingd==0) THEN
         PRINT*,'existance of zero nodes.. aborting'
         STOP
      ENDIF

      minbd = (maxgd+1)/2
c     Ggf. Kontrolldatei oeffnen
      
      if (kont) then

         fetxt = 'cutmck.ctr'

         errnr = 1
         open(11,file=fetxt,status='unknown',err=1000)

         errnr = 4
         write(11,4,err=1000) mingd,maxgd,minbd
 4       format(/'  minimum grade =',i4/'  maximum grade =',i4/
     1        '  minimum bandwidth =',i4)
      end if

c     ak
      spanz = sanz
      k = 0
      spanz = MIN(sanz,spmax)
      PRINT*,'mingd_1::',mingd
 110  do 120 i=1,sanz
         if (grad(i).eq.mingd) then
            k = k+1
            start(k) = i
            if (k.ge.spanz) goto 130
         end if
 120  continue

      mingd = mingd+1
      PRINT*,'mingd::',mingd
      goto 110

 130  continue

      spanz=k-1
      print*,'errechnete startpkte ',spanz
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
            write(11,8,err=1000) nstart,mb
 8          format(i8,8x,i4)
         end if

         if (mb.lt.mmin) then
            
            mmin = mb
            kbdm = is
            do i=1,sanz
               perm(i) = neuin(i)
            end do
         end if
         
      end do ! is

      if (kont) then
         write(11,9,err=1000) mmin,start(kbdm)
 9       format(//'  minimum bandwidth =',i5,'  for init.node',
     1        i5//'  vector of permutation ='/)

         write(11,11,err=1000) (perm(i),i=1,sanz)
 11      format((3x,10i5))

         close(11)
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
      delem=TRIM(ADJUSTL(delem))//'_ctm'
      call welem(11,delem)
      if (errnr.ne.0) goto 1000

      delectr=TRIM(ADJUSTL(delectr))//'_ctm'
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

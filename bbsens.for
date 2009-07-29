        subroutine bbsens(kanal,datei)

c Unterprogramm zur Berechnung der Summe der Sensitivitaeten aller
c Messungen (normiert)
c ('kpot' als Hilfsfeld benutzt).

c Andreas Kemna                                            02-Mar-1995
c                                       Letzte Aenderung   10-Mar-2007

c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'model.fin'
        INCLUDE 'inv.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Kanalnummer
        integer         * 4     kanal

c Datei
        character       * 80    datei

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvariablen
        real            * 8     xdum,ydum,dum

c Indexvariablen
        integer         * 4     i,j,k,l,i2

c Aktuelle Elementnummer
        integer         * 4     iel

c Anzahl der Knoten im aktuellen Elementtyp
        integer         * 4     nkel

c.....................................................................

c Werte berechnen
        dum = 0d0

        do j=1,manz
            l  = j/smax + 1
            i2 = mod(j,smax)

            kpot(i2,l,1) = dcmplx(0d0)
            kpot(i2,l,2) = dcmplx(0d0)
            kpot(i2,l,3) = dcmplx(0d0)
            kpot(i2,l,4) = dcmplx(0d0)

            do i=1,nanz
                kpot(i2,l,1) = kpot(i2,l,1) +
     1                         dconjg(sens(i,j))*sens(i,j)*
     1                         dcmplx(wmatd(i)*dble(wdfak(i)))
            end do

            dum = dmax1(dum,dble(kpot(i2,l,1)))
        end do

c Summe der Sensitivitaeten normieren
        if (dum.gt.1d-12) then
            do j=1,manz
                l = j/smax + 1
                i = mod(j,smax)
                kpot(i,l,1) = kpot(i,l,1) / dcmplx(dum)
            end do
        end if

c Schwerpunktkoordinaten der (ggf. zusammengesetzten) Elemente bestimmen
        iel = 0

        do i=1,typanz
            if (typ(i).gt.10) goto 10

            nkel = selanz(i)

            do j=1,nelanz(i)
                iel  = iel + 1

                xdum = 0d0
                ydum = 0d0

                do k=1,nkel
                    xdum = xdum + sx(snr(nrel(iel,k)))
                    ydum = ydum + sy(snr(nrel(iel,k)))
                end do

                l  = mnr(iel)/smax + 1
                i2 = mod(mnr(iel),smax)

                kpot(i2,l,2) = kpot(i2,l,2) + dcmplx(xdum/dble(nkel))
                kpot(i2,l,3) = kpot(i2,l,3) + dcmplx(ydum/dble(nkel))
                kpot(i2,l,4) = kpot(i2,l,4) + dcmplx(1d0)
            end do
        end do

10      continue

        do j=1,manz
            l = j/smax + 1
            i = mod(j,smax)
            kpot(i,l,2) = kpot(i,l,2) / kpot(i,l,4)
            kpot(i,l,3) = kpot(i,l,3) / kpot(i,l,4)
        end do

c 'datei' oeffnen
        fetxt = datei
        errnr = 1
        open(kanal,file=fetxt,status='replace',err=999)
        errnr = 4

cak
c Maximale Sensitivitaet schreiben
        write(kanal,*,err=1000) dum

c Koordinaten und Sensitivitaetsbetraege schreiben
c (logarithmierter (Basis 10) normierter Betrag)
        do j=1,manz
            l   = j/smax + 1
            i   = mod(j,smax)
            dum = dble(kpot(i,l,1))

            if (dum.gt.0d0) then
                write(kanal,*,err=1000) real(dble(kpot(i,l,2))),
     1                                  real(dble(kpot(i,l,3))),
     1                                  real(dlog10(dum))
            else
                write(kanal,*,err=1000) real(dble(kpot(i,l,2))),
     1                                  real(dble(kpot(i,l,3))),
     1                                  -1.e2
            end if
        end do

c 'datei' schliessen
        close(kanal)

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

999     return

1000    close(kanal)
        return

        end

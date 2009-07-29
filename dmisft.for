        subroutine dmisft(lsetup)

c Unterprogramm zum Bestimmen des Misfits der Daten.

c Andreas Kemna                                            01-Mar-1995
c                                       Letzte Aenderung   15-Jan-2001
         
c.....................................................................

        INCLUDE 'err.fin'
        INCLUDE 'path.fin'
        INCLUDE 'parmax.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'inv.fin'
        INCLUDE 'konv.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Hilfsschalter
        logical         * 4     lsetup

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsfelder
        real            * 8     eps2(nmax),psi(nmax)
        integer         * 4     wdlok(nmax)

c Hilfsvariablen
        integer         * 4     i,idum
        complex         * 16    cdum,cdat,csig
        real            * 8     dum,norm,norm2

c.....................................................................

c Ggf. Ausgabe
        if ((llam.and..not.lstep).or.lsetup) then
            fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
            errnr = 4
            write(15,*,err=1000)
            write(15,*,err=1000) it
        end if

c RMS-WERTE BERECHNEN
        nrmsd  = 0d0
        betrms = 0d0
        pharms = 0d0
        idum   = 0

        do i=1,nanz
            wdlok(i) = 1

c Phasen lokal korrigieren
            call chkpo2(dat(i),sigmaa(i),cdat,csig,wdlok(i),lpol)

cak Ggf. Daten mit Phase betraglich groesser 200 mrad nicht beruecksichtigen
cak (Standardabweichung auf 1d4 hochsetzen)
cak           if (dabs(1d3*dimag(csig)).gt.200d0) wmatd(i)=1d-8

cdiff-            cdum = cdat - csig
cdiff+<
            if (.not.ldiff) then
                cdum = cdat - csig
            else
                cdum = cdat - csig - (d0(i) - fm0(i))
            end if
cdiff+>

            if (lip) then
                psi(i) = dsqrt(wmatd(i))*dabs(dimag(cdum))
            else
                psi(i) = dsqrt(wmatd(i))*cdabs(cdum)
            end if

c Ggf. 'eps_i', 'psi_i' und Hilfsfeld ausgeben
            if ((llam.and..not.lstep).or.lsetup)
     1          write(15,*,err=1000) real(1d0/dsqrt(wmatd(i))),
     1                               real(psi(i)),wdlok(i)
     
            idum   = idum   + wdlok(i)
            nrmsd  = nrmsd  + psi(i)*psi(i)*dble(wdlok(i))
            betrms = betrms + wmatd(i)*dble(wdlok(i))*
     1                        dble(cdum)*dble(cdum)
            pharms = pharms + wmatd(i)*dble(wdlok(i))*
     1                        dimag(cdum)*dimag(cdum)
        end do
        
c Ggf. Fehlermeldung
        if (idum.eq.0) then
            fetxt = ' '
            errnr = 99
            goto 1000
        end if

        npol   = nanz-idum
        rmssum = nrmsd
        nrmsd  = dsqrt(nrmsd /dble(idum))
        betrms = dsqrt(betrms/dble(idum))
        pharms = dsqrt(pharms/dble(idum))

c Ggf. ROBUST INVERSION (nach Doug' LaBrecque)
        if (lrobust) then

c 'estimated weights' und 1-Normen berechnen
            norm  = 0d0
            norm2 = 0d0

            do i=1,nanz
                dum     = 1d0/dsqrt(wmatd(i))
                eps2(i) = dum*dsqrt(psi(i))
                norm    = norm  + psi(i)*dble(wdlok(i))
                norm2   = norm2 + psi(i)*dble(wdlok(i))*dum/eps2(i)
            end do

c 'estimated weights' normieren
            do i=1,nanz
                eps2(i) = eps2(i) * norm2/norm
            end do

c Kleinere Standardabweichung ausschliessen und neue 1-Norm berechnen
            norm2 = 0d0

            do i=1,nanz
                dum     = 1d0/dsqrt(wmatd(i))
                eps2(i) = dmax1(dum,eps2(i))
                norm2   = norm2 + psi(i)*dble(wdlok(i))*dum/eps2(i)
            end do

            l1rat = norm/norm2

c Ggf. neue Wichtungsfaktoren belegen
            if (l1rat.gt.l1min) then
                do i=1,nanz
                    dum = 1d0/eps2(i)/eps2(i)

c Ausgabe, falls 'eps_neu' > 1.1 * 'eps_alt'
                    if (dum.lt.0.83d0*wmatd(i).and.
     1                  ((llam.and..not.lstep).or.lsetup)) then
                        open(10,file=ramd(1:lnramd)//slash(1:1)
     1                                 //'run.ctr',status='old')
1                       read(10,*,end=2)
                            goto 1
2                       backspace(10)
                        write(10,'(i4,a32,f5.1)')
     1                           i,' : increase standard deviation *',
     1                           real(eps2(i)*dsqrt(wmatd(i)))
                        close(10)
                    end if

                    wmatd(i) = dum
                end do
            end if
        end if

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen
1000    close(14)
        close(15)
        return

        end

c*********************************************************************

        subroutine chkpo2(dati,sigi,cdat,csig,wdlok,ldum)

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Eingabe
        complex         * 16    dati,sigi

c Ausgabe
        complex         * 16    cdat,csig
        
c Schalter
        integer         * 4     wdlok
        logical         * 4     ldum

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvariablen
        integer         * 4     idat,isig

c Real-, Imaginaerteile
        real            * 8     redat,imdat,
     1                          resig,imsig

c Pi
        real            * 8     pi

c.....................................................................

        pi = dacos(-1d0)

c Logarithmierte Betraege in den Realteilen,
c Phasen (in rad) in den Imaginaerteilen
        redat = dble(dati)
        imdat = dimag(dati)
        resig = dble(sigi)
        imsig = dimag(sigi)

c Phasenbereich checken
        if (imdat.gt.pi/2d0) then
            idat = -1
        else if (imdat.le.-pi/2d0) then
            idat = 1
        else
            idat = 0
        end if

        if (imsig.gt.pi/2d0) then
            isig = -1
        else if (imsig.le.-pi/2d0) then
            isig = 1
        else
            isig = 0
        end if

        if (idat.eq.0.and.isig.ne.0) then

c Falls ldum=.true., angenommene Polaritaet des Messdatums falsch,
c ggf. Korrektur; auf jeden Fall Polaritaetswechsel
            imsig = imsig + dble(isig)*pi
            if (.not.ldum) imdat=imdat-dsign(pi,imdat)
            
            wdlok = 0

        else if (idat.ne.0.and.isig.eq.0) then

c Falls ldum=.true., angenommene Polaritaet des Messdatums falsch,
c ggf. Korrektur
            if (ldum) imdat=imdat+dble(idat)*pi
            
            wdlok = 0

        else if (idat.ne.0.and.isig.ne.0) then

c Polaritaetswechsel
            imsig = imsig + dble(isig)*pi
            imdat = imdat + dble(idat)*pi

        end if

c 'cdat' und 'csig' speichern
        cdat = dcmplx(redat,imdat)
        csig = dcmplx(resig,imsig)

        return
        end

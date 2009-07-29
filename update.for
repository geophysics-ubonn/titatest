        subroutine update(dpar2,cgres2)

c Unterprogramm zum Bestimmen und Anbringen der Modellverbesserung
c mittels 'Smoothness Least Squares Method' und konjugierten
c Gradienten.

c Andreas Kemna                                            01-Mar-1996
c                                       Letzte Aenderung   15-Jan-2001
       
c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'sigma.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'model.fin'
        INCLUDE 'fem.fin'
        INCLUDE 'inv.fin'
        INCLUDE 'konv.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Felder zum Zwischenspeichern
        complex         * 16    dpar2(mmax)
        real            * 4     cgres2(mmax+1)

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvariablen
        complex         * 16    cdum
        real            * 8     dum

c Hilfsfelder
        logical         * 4     lfeld(mmax)
        complex         * 16    bvec(mmax)

c Indexvariablen
        integer         * 4     i,j,k

c.....................................................................

c Parametervektor belegen
        do j=1,manz
            lfeld(j) = .false.
        end do

        do k=1,elanz
            j = mnr(k)

            if (.not.lfeld(j)) then
                lfeld(j) = .true.
                par(j)   = cdlog(sigma(k))
            end if
        end do

        if (.not.llam) then

c Felder speichern
            do j=1,manz
                dpar2(j) = dpar(j)
            end do

            i = int(cgres(1))+1
            do k=1,i
                cgres2(k) = cgres(k)
            end do

c Smoothnessvektor berechnen
            do i=1,manz
                cdum = dcmplx(0d0)

cdiff+<
              if (.not.ldiff) then
cdiff+>
                if (i.gt.1)
     1              cdum = dcmplx(smatm(i-1,2))*par(i-1)
                if (i.lt.manz)
     1              cdum = cdum + dcmplx(smatm(i,2))*par(i+1)
                if (i.gt.nx)
     1              cdum = cdum + dcmplx(smatm(i-nx,3))*par(i-nx)
                if (i.lt.manz-nx+1)
     1              cdum = cdum + dcmplx(smatm(i,3))*par(i+nx)

                bvec(i) = cdum + dcmplx(smatm(i,1))*par(i)
cdiff+<
              else
                if (i.gt.1)
     1            cdum = dcmplx(smatm(i-1,2))*(par(i-1)-m0(i-1))
                if (i.lt.manz)
     1            cdum = cdum + dcmplx(smatm(i,2))*(par(i+1)-m0(i+1))
                if (i.gt.nx)
     1            cdum = cdum + dcmplx(smatm(i-nx,3))*
     1                                             (par(i-nx)-m0(i-nx))
                if (i.lt.manz-nx+1)
     1            cdum = cdum + dcmplx(smatm(i,3))*(par(i+nx)-m0(i+nx))

                bvec(i) = cdum + dcmplx(smatm(i,1))*(par(i)-m0(i))
              end if
cdiff+>
            end do

c Skalierungsfaktoren bestimmen
            do j=1,manz
                dum = 0d0

                if (ldc) then
                    do i=1,nanz
                        dum = dum + sensdc(i,j)*sensdc(i,j)*
     1                              wmatd(i)*dble(wdfak(i))
                    end do
                else if (lip) then
                    do i=1,nanz
                        dum = dum + dble(sens(i,j))*dble(sens(i,j))*
     1                              wmatd(i)*dble(wdfak(i))
                    end do
                else
                    do i=1,nanz
                        dum = dum + dble(dconjg(sens(i,j))*sens(i,j))*
     1                              wmatd(i)*dble(wdfak(i))
                    end do
                end if

                dum    = dum + lam*smatm(j,1)
                fak(j) = 1d0/dsqrt(dum)
            end do

c Konstantenvektor berechen und skalieren
            do j=1,manz
                cdum = dcmplx(0d0)

cdiff+<
              if (.not.ldiff) then
cdiff+>
                if (ldc) then
                    do i=1,nanz
                        cdum = cdum + dcmplx(sensdc(i,j)*wmatd(i)*
     1                                dble(wdfak(i)))*(dat(i)-sigmaa(i))
                    end do
                else if (lip) then
                    do i=1,nanz
                        cdum = cdum + dcmplx(dble(sens(i,j))*wmatd(i)*
     1                                dble(wdfak(i)))*(dat(i)-sigmaa(i))
                    end do
                else
                    do i=1,nanz
                        cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)*
     1                                dble(wdfak(i)))*(dat(i)-sigmaa(i))
                    end do
                end if
cdiff+<
              else
                if (ldc) then
                    do i=1,nanz
                        cdum = cdum + dcmplx(sensdc(i,j)*wmatd(i)*
     1                                dble(wdfak(i)))*(dat(i)-sigmaa(i)-
     1                                                 (d0(i)-fm0(i)))
                    end do
                else if (lip) then
                    do i=1,nanz
                        cdum = cdum + dcmplx(dble(sens(i,j))*wmatd(i)*
     1                                dble(wdfak(i)))*(dat(i)-sigmaa(i)-
     1                                                 (d0(i)-fm0(i)))
                    end do
                else
                    do i=1,nanz
                        cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)*
     1                                dble(wdfak(i)))*(dat(i)-sigmaa(i)-
     1                                                 (d0(i)-fm0(i)))
                    end do
                end if
              end if
cdiff+>

                bvec(j) = cdum - dcmplx(lam)*bvec(j)
                bvec(j) = bvec(j)*dcmplx(fak(j))
            end do

c Modellverbesserung mittels konjugierter Gradienten bestimmen
            if (ldc.or.lip) then
                call cjggdc(bvec)
            else
                call cjggra(bvec)
            end if

c Ggf. Verbesserung umspeichern
            if (lip) then
                do j=1,manz
                    dpar(j) = dcmplx(0d0,dble(dpar(j)))
                end do
            end if

c Verbesserung skalieren
            do j=1,manz
                dpar(j) = dpar(j)*dcmplx(fak(j))
            end do
        else

c Felder zuruecksetzen
            do j=1,manz
                dpar(j) = dpar2(j)
            end do

            i = int(cgres2(1))+1
            do k=1,i
                cgres(k) = cgres2(k)
            end do
        end if

        do j=1,manz
        
c Verbesserung anbringen        
            par(j) = par(j) + dpar(j)*dcmplx(step)

c Ggf. (Leitfaehigkeits-)Phasen < 0 mrad korrigieren
            if (lphi0.and.dimag(par(j)).lt.0d0)
     1          par(j) = dcmplx(dble(par(j)))
cakc Ggf. (Leitfaehigkeits-)Phasen < 1 mrad korrigieren
cak            if (lphi0.and.dimag(par(j)).lt.1d-3)
cak     1          par(j) = dcmplx(dble(par(j)),1d-3)
        end do

c Betrag des Verbesserungsvektors bestimmen
        bdpar = 0d0

        do j=1,manz
            bdpar = bdpar + dble(dpar(j)*dconjg(dpar(j)))
        end do

        bdpar = dsqrt(bdpar/dble(manz))

        return
        end

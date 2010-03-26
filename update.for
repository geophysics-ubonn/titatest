      subroutine update(dpar2,cgres2)
c     
c     Unterprogramm zum Bestimmen und Anbringen der Modellverbesserung
c     mittels 'Smoothness Least Squares Method' und konjugierten
c     Gradienten.
c     Fuer beliebige Triangulierung und Stochastische Regularisierung
c     
c     Andreas Kemna                                            01-Mar-1996
c     
c     Letzte Aenderung                                         03-Aug-2009
c     
c.....................................................................

      USE alloci
      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'fem.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Felder zum Zwischenspeichern
      complex         * 16    dpar2(mmax)
      real            * 4     cgres2(mmax+1)

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      complex         * 16    cdum
      real            * 8     dum

c     Hilfsfelder
      complex         * 16    bvec(mmax)

c     Indexvariablen
      integer         * 4     i,j,k,ij,in
      integer         * 4     smaxs !
c.....................................................................

      smaxs=MAXVAL(selanz)

      if (.not.llam) then

c     Felder speichern
         dpar2 = dpar

         i = int(cgres(1))+1
         cgres2(1:i) = cgres(1:i)

c     Smoothnessvektor berechnen
c     triang>
         if (ltri==0) then
c     triang<
            do i=1,manz
               cdum = dcmplx(0d0)
c     diff+<
               if (.not.lprior) then
c     diff+>
                  if (i.gt.1) cdum =
     1                 dcmplx(smatm(i-1,2))*par(i-1)
                  if (i.lt.manz) cdum = cdum +
     1                 dcmplx(smatm(i,2))*par(i+1)
                  if (i.gt.nx) cdum = cdum +
     1                 dcmplx(smatm(i-nx,3))*par(i-nx)
                  if (i.lt.manz-nx+1) cdum = cdum +
     1                 dcmplx(smatm(i,3))*par(i+nx)
                  
                  bvec(i) = cdum + dcmplx(smatm(i,1))*par(i)
c     diff+<
               else
                  if (i.gt.1) cdum =
     1                 dcmplx(smatm(i-1,2))*(par(i-1)-m0(i-1))
                  if (i.lt.manz) cdum = cdum +
     1                 dcmplx(smatm(i,2))*(par(i+1)-m0(i+1))
                  if (i.gt.nx) cdum = cdum +
     1                 dcmplx(smatm(i-nx,3))*(par(i-nx)-m0(i-nx))
                  if (i.lt.manz-nx+1) cdum = cdum +
     1                 dcmplx(smatm(i,3))*(par(i+nx)-m0(i+nx))
                  
                  bvec(i)=cdum+dcmplx(smatm(i,1))*(par(i)-m0(i))
               end if
c     diff+>
            end do
c Damping--
         else if (ltri == 3.OR.ltri == 4) THEN

            WRITE(*,'(a)',ADVANCE='no')
     1           'update:: damping has no part in the gradient'

c     triang>
         else if (ltri == 1.OR.ltri == 2.OR.
     1           (ltri > 4 .AND. ltri < 15)) then
            do i=1,manz
               cdum = dcmplx(0d0)
               DO ij=1,smaxs
                  in=nachbar(i,ij)
                  IF (in/=0) then
                     if (.not.lprior) then
                        cdum = cdum + DCMPLX(smatm(i,ij))*
     1                       par(in)
                     else
                        cdum = cdum + DCMPLX(smatm(i,ij))* 
     1                       (par(in) - m0(in))
                     end if
                  end if
               END DO 
               if (.not.lprior) then
                  bvec(i) = cdum + DCMPLX(smatm(i,smaxs+1))*
     1                 par(i)
               else
                  bvec(i) = cdum + DCMPLX(smatm(i,smaxs+1))*
     1                 (par(i) - m0(i))
               end if
            end do

         else if (ltri == 15) THEN

            if (.not.lprior) then
               bvec(1:manz) = 
     1              MATMUL(DCMPLX(smatm),par(1:manz))
            else
               bvec(1:manz) = 
     1              MATMUL(dcmplx(smatm),(par(1:manz) - m0(1:manz)))
            end if
         END IF
c     triang<
         
c     Skalierungsfaktoren bestimmen
         do j=1,manz
            dum = 0d0

            if (ldc) then
               do i=1,nanz
                  dum = dum + sensdc(i,j)*
     1                 sensdc(i,j)*wmatd(i)*dble(wdfak(i))
               end do
            else if (lip) then
               do i=1,nanz
                  dum = dum + dble(sens(i,j))*dble(sens(i,j))*
     1                 wmatd(i)*dble(wdfak(i))
               end do
            else
               do i=1,nanz
                  dum = dum + dble(dconjg(sens(i,j))*
     1                 sens(i,j))*wmatd(i)*dble(wdfak(i))
               end do
            end if

c     triang< 
            IF (ltri==0) THEN
               dum    = dum + lam*smatm(j,1)

            ELSE IF (ltri == 1.OR.ltri == 2.OR.
     1              (ltri > 4 .AND. ltri < 15)) THEN
               dum    = dum + lam*smatm(j,smaxs+1)

            ELSE IF (ltri == 3.OR.ltri == 4) THEN
               dum = dum + lam * smatm(j,1)

            ELSE IF (ltri == 15) THEN
               dum    = dum + lam*smatm(j,j)

            END IF
c     triang> 
            
            fak(j) = 1d0/dsqrt(dum)

         end do

c     Konstantenvektor berechen und skalieren (RHS)
         do j=1,manz

            cdum = dcmplx(0d0)

c     diff+<
            if (.not.ldiff) then
c     diff+>
               if (ldc) then
                  do i=1,nanz
                     cdum = cdum + dcmplx(sensdc(i,j)*wmatd(i)*
     1                    dble(wdfak(i)))*(dat(i)-sigmaa(i))
                  end do
               else if (lip) then
                  do i=1,nanz
                     cdum = cdum + dcmplx(dble(sens(i,j))*wmatd(i)*
     1                    dble(wdfak(i)))*(dat(i)-sigmaa(i))
                  end do
               else
                  do i=1,nanz
                     cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)*
     1                    dble(wdfak(i)))*(dat(i)-sigmaa(i))
                  end do
               end if
c     diff+<
            else
               if (ldc) then
                  do i=1,nanz
                     cdum = cdum + dcmplx(sensdc(i,j)*wmatd(i)*
     1                    dble(wdfak(i)))*(dat(i)-sigmaa(i)-
     1                    (d0(i)-fm0(i)))
                  end do
               else if (lip) then
                  do i=1,nanz
                     cdum = cdum + dcmplx(dble(sens(i,j))*wmatd(i)*
     1                    dble(wdfak(i)))*(dat(i)-sigmaa(i)-
     1                    (d0(i)-fm0(i)))
                  end do
               else
                  do i=1,nanz
                     cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)*
     1                    dble(wdfak(i)))*(dat(i)-sigmaa(i)-
     1                    (d0(i)-fm0(i)))
                  end do
               end if
            end if
c     diff+>
            IF (ltri == 3.OR.ltri == 4) THEN

               bvec(j) = cdum

            ELSE

               bvec(j) = cdum - dcmplx(lam)*bvec(j)
            END IF

            bvec(j) = bvec(j)*dcmplx(fak(j))

         end do

c     Modellverbesserung mittels konjugierter Gradienten bestimmen
         if (ldc.or.lip) then
            call cjggdc(bvec)
         else
            call cjggra(bvec)
         end if

c     Ggf. Verbesserung umspeichern
         if (lip) then
            do j=1,manz
               dpar(j) = dcmplx(0d0,dble(dpar(j)))
            end do
         end if

c     Verbesserung skalieren
         do j=1,manz
            dpar(j) = dpar(j)*dcmplx(fak(j))
         end do

      else                      !(llam)


c     Felder zuruecksetzen
         do j=1,manz
            dpar(j) = dpar2(j)
         end do

         i = int(cgres2(1))+1

         do k=1,i
            cgres(k) = cgres2(k)
         end do
         
      end if
      do j=1,manz
         
c     Verbesserung anbringen
         par(j) = par(j) + dpar(j)*dcmplx(step)
c     Ggf. (Leitfaehigkeits-)Phasen < 0 mrad korrigieren
c     if (lphi0.and.dimag(par(j)).lt.0d0)
c     1        par(j) = dcmplx(dble(par(j)))
c     akc Ggf. (Leitfaehigkeits-)Phasen < 1 mrad korrigieren
c     ak            if (lphi0.and.dimag(par(j)).lt.1d-3)
c     ak     1          par(j) = dcmplx(dble(par(j)),1d-3)
      end do

c     i.e Stepsize = ||\delta m||
      bdpar = 0d0

      do j=1,manz
         bdpar = bdpar + dble(dpar(j)*dconjg(dpar(j)))
      end do
      bdpar = bdpar * step
c$$$      bdpar = dsqrt(bdpar/dble(manz))

      end

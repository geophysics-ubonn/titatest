      subroutine update()
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
      USE femmod
      USE datmod
      USE invmod
      USE cjgmod
      USE sigmamod
      USE modelmod
      USE elemmod
      USE cg_mod
      USE errmod
      USE konvmod

      IMPLICIT none


c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      complex         * 16    cdum
      real            * 8     dum

c     Indexvariablen
      integer         * 4     i,j,k,ij,in
c.....................................................................

      if (.not.llam) then

c     Felder speichern
         dpar2 = dpar
         cgres2 = cgres

C     felder allozieren
         IF (.NOT.ALLOCATED (bvec)) 
     1        ALLOCATE (bvec(manz),stat=errnr)
         IF (errnr /= 0) THEN
            fetxt = 'Error memory allocation bvec in update'
            errnr = 94
            RETURN
         END IF
 
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
         ELSE IF (ltri == 3.OR.ltri == 4) THEN

            WRITE(*,'(a)',ADVANCE='no')
     1           'update:: damping has no part in the gradient'

c     triang>
         ELSE IF (ltri == 1.OR.ltri == 2.OR.
     1           (ltri > 4 .AND. ltri < 15)) then
            DO i=1,manz
               cdum = dcmplx(0d0)
               DO ij=1,smaxs
                  in = nachbar(i,ij)
                  IF (in /= 0) THEN
                     IF (.NOT. lprior) THEN
                        cdum = cdum + DCMPLX(smatm(i,ij)) *
     1                       par(in)
                     ELSE
                        cdum = cdum + DCMPLX(smatm(i,ij)) * 
     1                       (par(in) - m0(in))
                     END IF
                  END IF
               END DO 
               IF (.NOT. lprior) THEN
                  bvec(i) = cdum + DCMPLX(smatm(i,smaxs+1)) *
     1                 par(i)
               ELSE
                  bvec(i) = cdum + DCMPLX(smatm(i,smaxs+1)) *
     1                 (par(i) - m0(i))
               END IF
            END DO

         ELSE IF (ltri == 15) THEN
            IF (.NOT. lprior) THEN
               bvec = MATMUL( DCMPLX(smatm),par)
            ELSE
               bvec = MATMUL( DCMPLX(smatm),( par - m0 ) )
            END IF
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
         CALL cjg

c     Ggf. Verbesserung umspeichern
         if (lip) then
            dpar = DCMPLX(0D0,DBLE(dpar))
         end if

c     Verbesserung skalieren
         do j=1,manz
            dpar(j) = dpar(j)*dcmplx(fak(j))
         end do

      else                      !(llam)


c     Felder zuruecksetzen
         dpar = dpar2

         i = int(cgres2(1))+1

         cgres(1:i) = cgres2(1:i)
         
      end if

      IF (ALLOCATED (bvec)) DEALLOCATE (bvec)

c     Ggf. (Leitfaehigkeits-)Phasen < 0 mrad korrigieren
c     if (lphi0.and.dimag(par(j)).lt.0d0)
c     1        par(j) = dcmplx(dble(par(j)))
c     akc Ggf. (Leitfaehigkeits-)Phasen < 1 mrad korrigieren
c     ak            if (lphi0.and.dimag(par(j)).lt.1d-3)
c     ak     1          par(j) = dcmplx(dble(par(j)),1d-3)

c     i.e Stepsize = ||\delta m||
      bdpar = 0d0
      in = 0
      do j=1,manz
         
         par(j) = par(j) + DCMPLX(step) * dpar(j) ! model update
         
c$$$  ! eventually correct for phase < 0 mrad
         IF (lphi0 .AND. dimag(par(j)) < 0D0) THEN 
            par(j) = DCMPLX(DBLE(par(j))) 
            in = in + 1 
         END IF 
         
         bdpar = bdpar + dble(dpar(j)*dconjg(dpar(j)))
         
      end do
      
      IF (in > 0) WRITE (*,'(a,I9,a)',ADVANCE = 'no')
     1     ' forcing zero ',in,' times'
      bdpar = bdpar * step

      end

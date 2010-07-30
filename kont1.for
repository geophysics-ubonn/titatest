cdiff-        subroutine kont1(delem,delectr,dstrom,drandb)
cdiff+<
      subroutine kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0)
c     diff+>

c     Unterprogramm zur Ausgabe der Kontrollvariablen.

c     Andreas Kemna                                            16-Apr-1996
c     Letzte Aenderung   16-Jul-2007
      
c.....................................................................

      USE variomodel
      USE femmod
      USE datmod
      USE cjgmod
      USE sigmamod
      USE modelmod
      USE elemmod
      USE wavenmod
      USE randbmod
      USE errmod
      USE konvmod

      IMPLICIT none

      INCLUDE 'path.fin'
c.....................................................................
c     EIN-/AUSGABEPARAMETER:

c     Dateinamen
      character       * 80    delem,
     1     delectr,
     1     dstrom,
c     diff+<
     1     dd0,
     1     dm0,
     1     dfm0,
c     diff+>
     1     drandb
      REAL(KIND(0D0)) :: Ix,Iy
 
      fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
      OPEN(fpinv,file=fetxt,status='old',POSITION='append',err=999)
      
c     HEADER AUSGEBEN
      write(fpinv,'(a11)',err=999) '***FILES***'
      write(fpinv,'(a80)',err=999) delem
      write(fpinv,'(a80)',err=999) delectr
      write(fpinv,'(a80)',err=999) dstrom
      write(fpinv,'(a60)',err=999) ramd
c     diff+<
      write(fpinv,'(l1,t18,a24)',err=999) ldiff,
     1     '! difference inversion ?'
      write(fpinv,'(a80)',err=999) dd0
      write(fpinv,'(l1,t18,a24)',err=999) lprior,
     1     '! smooth (m - m_{prior}) ?'
      write(fpinv,'(a80)',err=999) dm0
      write(fpinv,'(a80)',err=999) dfm0
c     diff+>
      write(fpinv,'(a16)',err=999) '***PARAMETERS***'
      IF (ltri == 0) THEN
         write(fpinv,'(i4,t18,a24)',err=999) nx,
     1        '! # cells in x-direction'
         write(fpinv,'(i4,t18,a24)',err=999) nz,
     1        '! # cells in z-direction'
      END IF
      write(fpinv,'(g11.5,t18,a36)',err=999) alfx,
     1     '! smoothing parameter in x-direction'
      write(fpinv,'(g11.5,t18,a36)',err=999) alfz,
     1     '! smoothing parameter in z-direction'
      write(fpinv,'(i2,t18,a29)',err=999) itmax,
     1     '! max. # inversion iterations'
c     ak        write(fpinv,'(g11.5,t18,a15)',err=999) nrmsdm,'! min. data RMS'
      write(fpinv,'(l1,t18,a16)',err=999) ldc,'! DC inversion ?'
c     ak        write(fpinv,'(l1,t18,a23)',err=999) lsr,'! singularity removal ?'
      write(fpinv,'(l1,t18,a20)',err=999) lrobust,'! robust inversion ?'
c     ak        write(fpinv,'(l1,t18,a33)',err=999) lpol,
c     ak     1           '! automatic polarity adjustment ?'
      write(fpinv,'(l1,t18,a27)',err=999) lfphai,
     1     '! final phase improvement ?'
c     ak        write(fpinv,'(l1,t18,a20)',err=999) lindiv,'! individual error ?'
      write(fpinv,'(g11.5,t18,a76)',err=999) stabw0,
     1     '! rel. resistance error level (%)'//
     1     '  (parameter A1 in err(R) = A1*abs(R) + A2)'
      write(fpinv,'(g11.5,t18,a76)',err=999) stabm0,
     1     '! min. abs. resistance error (ohm)'//
     1     ' (parameter A2 in err(R) = A1*abs(R) + A2)'
      write(fpinv,'(g11.5,t18,a93)',err=999) stabpA1,
     1     '! phase error model parameter A1 (mrad/ohm^B) '//
     1     '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
      write(fpinv,'(g11.5,t18,a93)',err=999) stabpB,
     1     '! phase error model parameter B  (-)          '//
     1     '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
      write(fpinv,'(g11.5,t18,a93)',err=999) stabpA2,
     1     '! phase error model parameter A2 (%)          '//
     1     '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
      write(fpinv,'(g11.5,t18,a93)',err=999) stabp0,
     1     '! phase error model parameter A3 (mrad)       '//
     1     '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
      write(fpinv,'(a,1X,l1)',err=999)
     1     '! (NEW) restart final phase with homogenous phase model?',
     1     lffhom
      write(fpinv,'(l1,t18,a38)',err=999) lrho0,
     1     '! homogeneous background resistivity ?'
      write(fpinv,'(g11.5,t18,a30)',err=999) bet0,
     1     '! background magnitude (ohm*m)'
      write(fpinv,'(g11.5,t18,a25)',err=999) pha0,
     1     '! background phase (mrad)'
      write(fpinv,'(i1,t18,a22)',err=999) swrtr,
     1     '! 2D (=0) or 2.5D (=1)'
      write(fpinv,'(l1,t18,a19)',err=999) lsink,
     1     '! fictitious sink ?'
      write(fpinv,'(i6,t18,a29)',err=999) nsink,
     1     '! fictitious sink node number'
      write(fpinv,'(l1,t18,a19)',err=999) lrandb2,
     1     '! boundary values ?'
      write(fpinv,'(a80)',err=999) drandb
      write(fpinv,'(/a)',err=999)      '***Model stats***'
      write(fpinv,*,err=999)'# Model parameters  : ',manz
      write(fpinv,*,err=999)'# Data points       : ',nanz
      write(fpinv,*,err=999)'Add data noise ?    : ',lnse
      write(fpinv,*,err=999)'Couple to Err. Modl?: ',.NOT.lnse2
      write(fpinv,*,err=999)'    seed            : ',iseed
      write(fpinv,*,err=999)'    Variance        : ',nstabw0
      write(fpinv,*,err=999)'Add model noise ?   : ',lnsepri
      write(fpinv,*,err=999)'    seed            : ',iseedpri
      write(fpinv,*,err=999)'    Variance        : ',modl_stdn
      write(fpinv,'(/a)',err=999)
     1     '******** Regularization Part *********'
      write(fpinv,*,err=999)'Regularization      : ',ltri
      write(fpinv,*,err=999)'Regular grid smooth : ',(ltri==0)
      write(fpinv,*,err=999)'Triangular regu     : ',(ltri==1)
      write(fpinv,*,err=999)'Triangular regu2    : ',(ltri==2)
      write(fpinv,*,err=999)'Levenberg damping   : ',(ltri==3)
      write(fpinv,*,err=999)'Marquardt damping   : ',(ltri==4)
      write(fpinv,*,err=999)'Minimum grad supp   : ',(ltri==5)
      write(fpinv,*,err=999)'MGS beta/sns1 (RM)  : ',(ltri==6)
      write(fpinv,*,err=999)'MGS beta/sns2 (RM)  : ',(ltri==7)
      write(fpinv,*,err=999)'MGS beta/sns1 (RB)  : ',(ltri==8)
      write(fpinv,*,err=999)'MGS beta/sns2 (RB)  : ',(ltri==9)
      write(fpinv,*,err=999)'TV (Huber)          : ',(ltri==10)

      IF (ltri>4.AND.ltri<15)
     1     write(fpinv,*,err=999)'  Stabilizer beta  : ',betamgs
      write(fpinv,*,err=999)'Stochastic regu     : ',(ltri==15)

      IF (lvario) THEN
         WRITE (fpinv,'(a)',err=999)'Experimental Variogram::'
         WRITE (fpinv,'(a,I4)',err=999)ACHAR(9)//
     1        'nx-switch  : ',nx
         CALL get_vario (Ix,Iy,fetxt,0) ! get korrelation lengths
         WRITE (fpinv,'(2(a,F5.2))',ERR=999)ACHAR(9)//
     1        'Integral lengths Ix/Iy',Ix,'/',Iy
         WRITE (*,'(/2(a,F5.2))')ACHAR(9)//
     1        'Integral lengths Ix/Iy',Ix,'/',Iy
         WRITE (fpinv,'(a)',err=999)ACHAR(9)//
     1        'Variogram('//TRIM(fetxt)//')'
         WRITE (*,'(a)')ACHAR(9)//
     1        'Variogram('//TRIM(fetxt)//')'
         IF (ltri == 15) THEN
            CALL get_vario (Ix,Iy,fetxt,1) ! get covariance..
            WRITE (fpinv,'(a)',err=999)ACHAR(9)//
     1           'Covariance('//TRIM(fetxt)//')'
            WRITE (*,'(a)')ACHAR(9)//
     1           'Covariance('//TRIM(fetxt)//')'
         END IF
      END IF

      write(fpinv,*,err=999)'Fixed lambda       : ',llamf,lamfix
      write(fpinv,'(/a)',err=999)
     1     '******** Additional output *********'
      write(fpinv,*,err=999)'Read start model   : ',lstart
      write(fpinv,*,err=999)'Write coverage     : ',lsens
      write(fpinv,*,err=999)'Write MCM 1        : ',lcov1
      write(fpinv,*,err=999)'Write resolution   : ',lres
      write(fpinv,*,err=999)'Write MCM 2        : ',lcov2
      write(fpinv,*,err=999)'Using Gauss ols    : ',lgauss
      IF (nz<0) THEN
         write(fpinv,'(1x,a)',err=999,ADVANCE='no')
     1        'taking easy lam_0 : '
         IF (nz<-1) write(fpinv,*,err=999) -REAL(nz)
         IF (nz==-1) write(fpinv,*,err=999)MAX(REAL(manz),REAL(nanz))
      END IF

      write(fpinv,'(/a)',err=999) '***FIXED***'
      if (swrtr.eq.1) then
         write(fpinv,'(a,t50,i2)',err=999) ' # wavenumbers :',kwnanz
         write(fpinv,'(a,t50,g11.5,t62,a1)',err=999)
     1        ' Inverse Fourier transform range :',amin,'m'
         write(fpinv,'(t50,g11.5,t62,a1)',err=999) amax,'m'
      end if
      if (.not.lrho0.and..not.lstart) then
         bet0 = cdabs(dcmplx(1d0)/sigma0)
         pha0 = 1d3*datan2(dimag(dcmplx(1d0)/sigma0),
     1        dble(dcmplx(1d0)/sigma0))
         write(fpinv,'(a,t50,g11.5,t62,a5)',err=999)
     1        ' Background resistivity :',bet0,'ohm*m'
         write(fpinv,'(t50,g11.5,t62,a4)',err=999)
     1        pha0,'mrad'
      end if
      write(fpinv,'(a,t50,l1)',err=999)
     1     ' Force negative phase ?',lphi0
      write(fpinv,'(a16,t50,l1)',err=999) ' Ratio dataset ?',lratio
      if (lrobust.or.lfphai) write(fpinv,'(a13,t50,g11.5)',err=999)
     1     ' Min. L1 norm',l1min
      write(fpinv,'(a,t50,g11.5)',err=999)
     1     ' Min. rel. decrease of data RMS :',mqrms
      write(fpinv,'(a,t50,g11.5)',err=999)
     1     ' Min. steplength              :',stpmin
      write(fpinv,'(a,t50,g11.5)',err=999)
     1     ' Min. stepsize (||\delta m||) :',bdmin
      write(fpinv,'(a,t50,g11.5)',err=999)
     1     ' Min. error in relaxation :',eps
      write(fpinv,'(a,t50,i5)',err=999)
     1     ' Max. # relaxation iterations :',ncgmax
      write(fpinv,'(a,t50,i3)',err=999)
     1     ' Max. # regularization steps :',nlam
      write(fpinv,'(a,t50,g11.5)',err=999)
     1     ' Initial step factor :',fstart
      write(fpinv,'(a,t50,g11.5)',err=999)
     1     ' Final   step factor :',fstop
      write(fpinv,*,err=999)
      write(fpinv,'(a48,a48,a13)',err=999)
     1     '------------------------------------------------',
     1     '------------------------------------------------',
     1     '-------------'
      write(fpinv,*,err=999)
c     Robuste Inversion
      if (lrobust) then
         write(fpinv,'(t1, a3, t5,a3,t11,a8,t23,a8,t34,a8,t46,a8,t58,a8,
     1t69,a8,t81,a8,t93,a8,t105,a8,t117,a10)',err=999)
     1        'ID','it.','data RMS','stepsize',' lambda ',' roughn.',
     1        'CG-steps',' mag RMS',' pha RMS','- # data',
     1        'L1-ratio','steplength'
      else
         write(fpinv,'(t1, a3, t5,a3,t11,a8,t23,a8,t34,a8,t46,a8,t58,a8,
     1t69,a8,t81,a8,t93,a8,t105,a10)',err=999)
     1        'ID','it.','data RMS','stepsize',' lambda ',' roughn.',
     1        'CG-steps',' mag RMS',' pha RMS','- # data',
     1        'steplength'
      end if
      write(fpinv,*,err=999)

c     'inv.ctr' schliessen
      close(fpinv)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

      end

cdiff-        subroutine kont1(delem,delectr,dstrom,drandb)
cdiff+<
        subroutine kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0)
cdiff+>

c Unterprogramm zur Ausgabe der Kontrollvariablen.

c Andreas Kemna                                            16-Apr-1996
c                                       Letzte Aenderung   16-Jul-2007
 
c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'path.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'sigma.fin'
        INCLUDE 'waven.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'fem.fin'
        INCLUDE 'konv.fin'
        INCLUDE 'randb.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Dateinamen
        character       * 80    delem,
     1                          delectr,
     1                          dstrom,
cdiff+<
     1                          dd0,
     1                          dm0,
     1                          dfm0,
cdiff+>
     1                          drandb

c.....................................................................

c Kontrolldateien oeffnen
        errnr = 1
        fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
        open(13,file=fetxt,status='replace',err=999)
      
        fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
        open(14,file=fetxt,status='replace',err=999)

        fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
        open(15,file=fetxt,status='replace',err=999)

        errnr = 4

c HEADER AUSGEBEN
        write(13,'(a11)',err=999) '***FILES***'
        write(13,'(a80)',err=999) delem
        write(13,'(a80)',err=999) delectr
        write(13,'(a80)',err=999) dstrom
        write(13,'(a60)',err=999) ramd
cdiff+<
        write(13,'(l1,t18,a24)',err=999) ldiff,
     1           '! difference inversion ?'
        write(13,'(a80)',err=999) dd0
        write(13,'(a80)',err=999) dm0
        write(13,'(a80)',err=999) dfm0
cdiff+>
        write(13,'(a16)',err=999) '***PARAMETERS***'
        write(13,'(i4,t18,a24)',err=999) nx,'! # cells in x-direction'
        write(13,'(i4,t18,a24)',err=999) nz,'! # cells in z-direction'
        write(13,'(g11.5,t18,a36)',err=999) alfx,
     1           '! smoothing parameter in x-direction'
        write(13,'(g11.5,t18,a36)',err=999) alfz,
     1           '! smoothing parameter in z-direction'
        write(13,'(i2,t18,a29)',err=999) itmax,
     1           '! max. # inversion iterations'
cak        write(13,'(g11.5,t18,a15)',err=999) nrmsdm,'! min. data RMS'
        write(13,'(l1,t18,a16)',err=999) ldc,'! DC inversion ?'
cak        write(13,'(l1,t18,a23)',err=999) lsr,'! singularity removal ?'
        write(13,'(l1,t18,a20)',err=999) lrobust,'! robust inversion ?'
cak        write(13,'(l1,t18,a33)',err=999) lpol,
cak     1           '! automatic polarity adjustment ?'
        write(13,'(l1,t18,a27)',err=999) lfphai,
     1           '! final phase improvement ?'
cak        write(13,'(l1,t18,a20)',err=999) lindiv,'! individual error ?'
        write(13,'(g11.5,t18,a76)',err=999) stabw0,
     1           '! rel. resistance error level (%)'//
     1           '  (parameter A1 in err(R) = A1*abs(R) + A2)'
        write(13,'(g11.5,t18,a76)',err=999) stabm0,
     1           '! min. abs. resistance error (ohm)'//
     1           ' (parameter A2 in err(R) = A1*abs(R) + A2)'
        write(13,'(g11.5,t18,a93)',err=999) stabpA1,
     1           '! phase error model parameter A1 (mrad/ohm^B) '//
     1           '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
        write(13,'(g11.5,t18,a93)',err=999) stabpB,
     1           '! phase error model parameter B  (-)          '//
     1           '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
        write(13,'(g11.5,t18,a93)',err=999) stabpA2,
     1           '! phase error model parameter A2 (%)          '//
     1           '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
        write(13,'(g11.5,t18,a93)',err=999) stabp0,
     1           '! phase error model parameter A3 (mrad)       '//
     1           '(in err(pha) = A1*abs(R)**B + A2*abs(pha) + A3)'
        write(13,'(l1,t18,a38)',err=999) lrho0,
     1           '! homogeneous background resistivity ?'
        write(13,'(g11.5,t18,a30)',err=999) bet0,
     1           '! background magnitude (ohm*m)'
        write(13,'(g11.5,t18,a25)',err=999) pha0,
     1           '! background phase (mrad)'
        write(13,'(i1,t18,a22)',err=999) swrtr,
     1           '! 2D (=0) or 2.5D (=1)'
        write(13,'(l1,t18,a19)',err=999) lsink,
     1           '! fictitious sink ?'
        write(13,'(i6,t18,a29)',err=999) nsink,
     1           '! fictitious sink node number'
        write(13,'(l1,t18,a19)',err=999) lrandb2,
     1           '! boundary values ?'
        write(13,'(a80)',err=999) drandb
        write(13,'(a11)',err=999) '***FIXED***'
        if (swrtr.eq.1) then
            write(13,'(a16,t50,i2)',err=999) ' # wavenumbers :',kwnanz
            write(13,'(a34,t50,g11.5,t62,a1)',err=999)
     1               ' Inverse Fourier transform range :',amin,'m'
            write(13,'(t50,g11.5,t62,a1)',err=999) amax,'m'
        end if
        if (.not.lrho0.and..not.lstart) then
            bet0 = cdabs(dcmplx(1d0)/sigma0)
            pha0 = 1d3*datan2(dimag(dcmplx(1d0)/sigma0),
     1                         dble(dcmplx(1d0)/sigma0))
            write(13,'(a25,t50,g11.5,t62,a5)',err=999)
     1               ' Background resistivity :',bet0,'ohm*m'
            write(13,'(t50,g11.5,t62,a4)',err=999)
     1               pha0,'mrad'
        end if
        write(13,'(a23,t50,l1)',err=999)
     1           ' Force negative phase ?',lphi0
        write(13,'(a16,t50,l1)',err=999) ' Ratio dataset ?',lratio
        if (lrobust.or.lfphai) write(13,'(a13,t50,g11.5)',err=999)
     1                                  ' Min. L1 norm',l1min
        write(13,'(a33,t50,g11.5)',err=999)
     1           ' Min. rel. decrease of data RMS :',mqrms
        write(13,'(a19,t50,g11.5)',err=999)
     1           ' Min. step-length :',stpmin
        write(13,'(a27,t50,g11.5)',err=999)
     1           ' Min. error in relaxation :',eps
        write(13,'(a31,t50,i5)',err=999)
     1           ' Max. # relaxation iterations :',ncgmax
        write(13,'(a30,t50,i3)',err=999)
     1           ' Max. # regularization steps :',nlam
        write(13,'(a22,t50,g11.5)',err=999)
     1           ' Initial step factor :',fstart
        write(13,'(a22,t50,g11.5)',err=999)
     1           ' Final   step factor :',fstop
        write(13,*,err=999)
        write(13,'(a48,a48,a12)',err=999)
     1           '------------------------------------------------',
     1           '------------------------------------------------',
     1           '------------'
        write(13,*,err=999)
        if (lrobust.or.lfphai) then
            write(13,'(t1,a3,t7,a8,t19,a8,t30,a8,t42,a8,t54,a8,t65,a8,
     1                 t77,a8,t89,a8,t101,a8)',err=999)
     1               'it.','data RMS','stepsize',' lambda ',' roughn.',
     1               'CG-steps',' mag RMS',' pha RMS','- # data',
     1               'L1-ratio'
        else
            write(13,'(t1,a3,t7,a8,t19,a8,t30,a8,t42,a8,t54,a8,t65,a8,
     1                 t77,a8,t89,a8)',err=999)
     1               'it.','data RMS','stepsize',' lambda ',' roughn.',
     1               'CG-steps',' mag RMS',' pha RMS','- # data'
        end if
        write(13,*,err=999)

c 'inv.ctr' schliessen
        close(13)

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

999     return

        end

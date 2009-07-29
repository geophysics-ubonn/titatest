        subroutine parfit(fa,fb,fc,fmin,smin)

c Unterprogramm fittet Parabel durch die drei Punkte (0,fa), (0.5,fb)
c und (1,fc) und liefert Abszisse des Minimums in 'step', falls Minimum
c zwischen 0 und 1 existiert und Minimum < fmin. Sonst geeignete lineare
c Interpolation auf fmin.
c (Formel aus 'Numerical Recipes' (S. 395, eq. (10.2.1)))

c Andreas Kemna                                            28-May-1996
c                                       Letzte Aenderung   22-Sep-1998
        
c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'konv.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Funktionswerte, Grenzwerte
        real            * 8     fa,fb,fc,fmin,smin

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c x-Werte
        real            * 8     a,b,c

c Hilfsvariablen
        real            * 8     bma,bmc,
     1                          fbmfc,fbmfa,
     1                          zaehler,nenner

c.....................................................................

        a = 0d0
        b = 0.5d0
        c = 1d0

        if (fa.gt.fb) then

            if (fb.le.fc) then

                if (fb.le.fmin) then
                    if (fa.lt.fmin.and.fc.gt.fa) then

c                     Zwischen 'fb' und 'fc' linear auf fmin interpolieren
                        step = b + (c-b)*(fmin-fb)/(fc-fb)
                    else 

c                     Zwischen 'fa' und 'fb' linear auf fmin interpolieren
                        step = a + (b-a)*(fmin-fa)/(fb-fa)
                    end if
                else                    	 

c                 Parabolische Interpolation auf Minimum
                    bma     = b-a
                    fbmfc   = fb-fc
                    bmc     = b-c
                    fbmfa   = fb-fa
                    zaehler = bma*bma*fbmfc - bmc*bmc*fbmfa
                    nenner  = bma*fbmfc - bmc*fbmfa
                    step    = b - zaehler/(2d0*nenner)
                end if

            else if (fb.gt.fc) then

                if (fb.le.fmin) then

c                 Zwischen 'fa' und 'fb' linear auf fmin interpolieren
                    step = a + (b-a)*(fmin-fa)/(fb-fa)
                else if (fc.le.fmin) then

c                 Zwischen 'fb' und 'fc' linear auf fmin interpolieren
                    step = b + (c-b)*(fmin-fb)/(fc-fb)
                else

c                 Full step-length
                    step = c
                end if

            end if

        else if (fa.lt.fb) then

            if (fb.ge.fc) then

                if (fb.ge.fmin) then
                    if (fa.gt.fmin.and.fc.lt.fa) then

c                     Zwischen 'fb' und 'fc' linear auf fmin interpolieren
                        step = b + (c-b)*(fmin-fb)/(fc-fb)
                    else 

c                     Zwischen 'fa' und 'fb' linear auf fmin interpolieren
                        step = a + (b-a)*(fmin-fa)/(fb-fa)
                    end if
                else                    	 

c                 Parabolische Interpolation auf Maximum
                    bma     = b-a
                    fbmfc   = fb-fc
                    bmc     = b-c
                    fbmfa   = fb-fa
                    zaehler = bma*bma*fbmfc - bmc*bmc*fbmfa
                    nenner  = bma*fbmfc - bmc*fbmfa
                    step    = b - zaehler/(2d0*nenner)
                end if

            else if (fb.lt.fc) then

                if (fb.ge.fmin) then

c                 Zwischen 'fa' und 'fb' linear auf fmin interpolieren
                    step = a + (b-a)*(fmin-fa)/(fb-fa)
                else if (fc.ge.fmin) then

c                 Zwischen 'fb' und 'fc' linear auf fmin interpolieren
                    step = b + (c-b)*(fmin-fb)/(fc-fb)
                else

c                 Full step-length
                    step = c
                end if

            end if

        else if (fa.eq.fb) then

            if (fb.gt.fc) then

                if (fb.lt.fmin) then

c                 Parabolische Interpolation auf Maximum
                    bma     = b-a
                    fbmfc   = fb-fc
                    bmc     = b-c
                    fbmfa   = fb-fa
                    zaehler = bma*bma*fbmfc - bmc*bmc*fbmfa
                    nenner  = bma*fbmfc - bmc*fbmfa
                    step    = b - zaehler/(2d0*nenner)
                else if (fc.le.fmin) then

c                 Zwischen 'fb' und 'fc' linear auf fmin interpolieren
                    step = b + (c-b)*(fmin-fb)/(fc-fb)
                else

c                 Full step-length
                    step = c
                end if

            else if (fb.lt.fc) then

                if (fb.gt.fmin) then

c                 Parabolische Interpolation auf Minimum
                    bma     = b-a
                    fbmfc   = fb-fc
                    bmc     = b-c
                    fbmfa   = fb-fa
                    zaehler = bma*bma*fbmfc - bmc*bmc*fbmfa
                    nenner  = bma*fbmfc - bmc*fbmfa
                    step    = b - zaehler/(2d0*nenner)
                else if (fc.ge.fmin) then

c                 Zwischen 'fb' und 'fc' linear auf fmin interpolieren
                    step = b + (c-b)*(fmin-fb)/(fc-fb)
                else

c                 Full step-length
                    step = c
                end if

            else if (fb.eq.fc) then

c             Mindest-step-length
                step = smin
            end if

        end if

c Genzen einhalten (wegen Möglichkeit der linearen "Extrapolation")
        step = dmax1(smin,step)
        step = dmin1(c,step)

        return
        end

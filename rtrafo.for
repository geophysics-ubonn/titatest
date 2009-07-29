        subroutine rtrafo()
  
c Unterprogramm zur Ruecktransformation.

c Andreas Kemna                                            20-Dec-1993
c                                       Letzte Aenderung   13-Nov-1997

c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'electr.fin'
        INCLUDE 'waven.fin'
        INCLUDE 'fem.fin'

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Pi
        real            * 8     pi

c Hilfsvariablen
        complex         * 16    summe
        real            * 8     summdc

c Indexvariablen
        integer         * 4     j,k,l

c.....................................................................

        pi = dacos(-1d0)

        if (ldc) then

            do l=1,eanz
                do j=1,sanz
                    summdc = 0d0

                    do k=1,kwnanz
                        summdc = summdc + kpotdc(j,l,k)*kwnwi(k)
                    end do

                    hpotdc(j,l) = summdc / pi
                end do
            end do

	  else

            do l=1,eanz
                do j=1,sanz
                    summe = dcmplx(0d0)

                    do k=1,kwnanz
                        summe = summe + kpot(j,l,k)*dcmplx(kwnwi(k))
                    end do

                    hpot(j,l) = summe / dcmplx(pi)
                end do
            end do

	  end if

        return
        end

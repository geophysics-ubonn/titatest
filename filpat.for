      character * (*) function filpat(disfile,ln,sw,slash)

c     Die Funktion separiert aus einem kompletten Pfad 'disfile' den
c     File-Namen (sw = 0) bzw. den Ordner (sw = 1) der Laenge 'ln'.

c     Andreas Kemna                                            19-Jan-1993
c     Letzte Aenderung   24-Oct-1996

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

      IMPLICIT none
      character       * (*)   disfile
      integer         * 4     ln
      integer         * 4     sw
      character       * 1     slash

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

      integer         * 4     id1,id2

c.....................................................................

      id2 = 0
 10   id1 = index(disfile(id2+1:len(disfile)),slash(1:1))
      if (id1.ne.0) then
         id2 = id1+id2
         goto 10
      end if

      if (sw.eq.0) then
         filpat = disfile(id2+1:len(disfile))
      else
         if (id2.ge.2) then
            filpat = disfile(1:id2-1)
         else
            filpat = ' '
         end if
      end if

      ln = index(filpat,' ')-1
      ln = max(ln,1)

      return
      end

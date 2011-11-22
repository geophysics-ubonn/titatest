MODULE pathmod
!!$c 'path.fin'
!!$
!!$c Andreas Kemna                                            25-Jun-1996
!!$c                                       Letzte Aenderung   24-Oct-1996
!!$
!!$c.....................................................................

!!$c (Back-) Slash
  CHARACTER(1),PUBLIC :: slash

!!$c RAM-Disk-Pfad
  CHARACTER (60),PUBLIC ::    ramd

!!$c Laenge des RAM-Disk-Pfades
  INTEGER(KIND = 4),PUBLIC ::     lnramd

!!$c UNIX COMMANDS
  CHARACTER(6),PUBLIC,PARAMETER :: mkdir='mkdir '
  CHARACTER(6),PUBLIC,PARAMETER :: rmdir='rm -R '

END MODULE pathmod

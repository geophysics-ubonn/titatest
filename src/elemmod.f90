MODULE elemmod
!!$c 'elem.fin'
!!$
!!$c Andreas Kemna                                            24-Nov-1993
!!$c                                       Letzte Aenderung   20-Nov-2009
!!$         
!!$c.....................................................................

!!!$Anzahl der Knoten (bzw. Knotenvariablen)
  INTEGER(KIND = 4),PUBLIC                            :: sanz
!!!$Anzahl der Elementtypen
  INTEGER(KIND = 4),PUBLIC                            :: typanz
!!!$Bandbreite der Gesamtsteifigkeitsmatrix 'a'
  INTEGER(KIND = 4),PUBLIC                            ::  mb

!!!$Elementtypen
!!!$(Randelemente (ntyp > 10) am Schluss !)
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE   :: typ

!!!$Anzahl der Elemente eines bestimmten Typs
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE   :: nelanz

!!!$Anzahl der Knoten (bzw. Knotenvariablen) in einem Elementtyp
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE   :: selanz

!!!$Zeiger auf Koordinaten der Knoten
!!!$(Inverser Permutationsvektor der Umnumerierung)
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE   :: snr

!!!$x-Koordinaten der Knoten
  REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: sx

!!!$y-Koordinaten der Knoten
  REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: sy

!!!$ Elementschwerpunktkoordinaten (ESP) der Flaechenelemente
  REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: espx,espy

!!!$Knotennummern der Elemente (Reihenfolge !)
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:,:),ALLOCATABLE :: nrel

!!!$Anzahl der Elemente (ohne Randelemente)
  INTEGER(KIND = 4),PUBLIC                            :: elanz

!!!$Anzahl der Randelemente
  INTEGER(KIND = 4),PUBLIC                            :: relanz

!!!$Zeiger auf Werte der Randelemente
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE   :: rnr

!!!$ Groeste Anzahl der Knoten der Flaechenelemente
  INTEGER(KIND = 4),PUBLIC                            :: smaxs

!!!$Gitter statistiken..
!!!$Min/Maximaler Abstand zwischen (Flaechen) Elementschwerpunkten
  REAL(KIND(0D0)),PUBLIC                              :: esp_min,&
       esp_max
!!!$Mittelwert/Median und Standardabweichung der ESP
  REAL(KIND(0D0)),PUBLIC                              :: esp_mit,&
       esp_med,esp_std

!!!$Min/Maximale  Gitterabstaende
  REAL(KIND(0D0)),PUBLIC                              :: grid_min,&
       grid_max
!!!$Min/Maximaler Gitterabstand in x-,y-richtung
  REAL(KIND(0D0)),PUBLIC                              :: grid_minx,&
       grid_maxx
  REAL(KIND(0D0)),PUBLIC                              :: grid_miny,&
       grid_maxy
!!!$switch/number fictitious sink node (only for 2D)
  LOGICAL,PUBLIC                                      :: lsink
!!!$ number of grid node for sink
  INTEGER(KIND = 4),PUBLIC                            :: nsink
!!$ switch boundary values
  LOGICAL,PUBLIC                                      :: lrandb2

!!!$mittlere y-Koordinate aller Randelemente vom Typ 12 ("no flow")
  REAL(KIND(0D0)),PUBLIC                              :: sytop 

!!!$    x-Koordinaten der Eckknotenpunkte
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: xk

!!!$    y-Koordinaten der Eckknotenpunkte
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: yk

!!!$    Elementmatrizen
  REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: elmas,elmam

!!!$    Elementvektor
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: elve

END MODULE elemmod

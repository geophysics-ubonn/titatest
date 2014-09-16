MODULE get_ver
! get git version information from file 'my_git_version.h' which is created
! at each ./autogen.sh run

  IMPLICIT none

  CHARACTER (256)  ::   version(5)

CONTAINS

  SUBROUTINE get_git_ver(version)
  CHARACTER (256)  ::   version(5)

  INCLUDE 'my_git_version.h'
  
  version(1)=TRIM(ADJUSTL(my_git_version(1)))
  version(2)=TRIM(ADJUSTL(my_git_version(2)))
  version(3)=TRIM(ADJUSTL(my_git_version(3)))
  version(4)=TRIM(ADJUSTL(my_git_version(4)))
  version(5)=TRIM(ADJUSTL(my_git_version(5)))
  
  PRINT*,'git branch        ',TRIM(version(1))
  PRINT*,'commit ID         ',TRIM(version(2))
  PRINT*,'created on        ',TRIM(version(3))
  PRINT*,'compiler          ',TRIM(version(4))
  PRINT*,'opertating system ',TRIM(version(5))
  
END SUBROUTINE get_git_ver
END MODULE get_ver

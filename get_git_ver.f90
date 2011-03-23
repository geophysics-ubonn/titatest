MODULE get_ver

  IMPLICIT none

  CHARACTER (256),PUBLIC  ::   version(3)

  PUBLIC :: get_git_ver
  
CONTAINS

  SUBROUTINE get_git_ver

  INCLUDE 'my_git_version.h'
  
  version(1)=TRIM(ADJUSTL(my_git_version(1)))
  version(2)=TRIM(ADJUSTL(my_git_version(2)))
  version(3)=TRIM(ADJUSTL(my_git_version(3)))
  
  PRINT*,'Git-Branch  ',TRIM(version(1))
  PRINT*,'Commit-ID   ',TRIM(version(2))
  PRINT*,'Created     ',TRIM(version(3))
  
END SUBROUTINE get_git_ver
END MODULE get_ver
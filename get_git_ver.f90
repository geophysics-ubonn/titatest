SUBROUTINE get_git_ver(buff)

  CHARACTER(256),INTENT (OUT) :: buff(3)
  
  INCLUDE 'my_git_version.h'
  
  buff(1)=TRIM(ADJUSTL(my_git_version(1)))
  buff(2)=TRIM(ADJUSTL(my_git_version(2)))
  buff(3)=TRIM(ADJUSTL(my_git_version(3)))
  
  
END SUBROUTINE get_git_ver

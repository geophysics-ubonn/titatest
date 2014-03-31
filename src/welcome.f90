subroutine welcome
use errmod
use get_ver
use alloci
implicit none
integer pid
pid = getpid()
  fetxt = 'crtomo.pid'
  print*,'THIS IS CRTOMO'
  print*,'Process_ID ::',pid
  open (fprun,FILE=trim(fetxt),STATUS='replace',err=999)
  write (fprun,*)pid
  close (fprun)

  print*,'GIT version information:'
  ! Get git version and store it in 'version'
  call get_git_ver(version)
    ! Display stiffness matrix precision
  write(*,'(a,i3,a)')' Precision  ',precision(a),' digits' 
  write(*,*) '------------------------------------------'
return
999 print*,'Error: could not create crtomo.pid'
end subroutine welcome


program main
  use Photochem, only: Atmosphere, err_len, real_kind
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  logical :: success
  character(len=:), allocatable :: reaction
  integer :: i
  real(real_kind) :: tn
  
  err = ""
  call pc%init("../data", &
               "../data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  ! pc%var%verbose = 0
  
  call pc%initialize_stepper(pc%var%usol_init,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  do i = 1,10
    tn = pc%step(err)
    if (len(trim(err)) > 0) then
      print*,trim(err)
      stop 1
    endif
  enddo
  
  call pc%destroy_stepper(err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  call pc%photochemical_equilibrium(success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  ! do i = 1, pc%dat%nrT
  !   print*,i,trim(pc%dat%reaction_equations(i))
  ! enddo

end program

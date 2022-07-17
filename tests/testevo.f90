program testevo
  use photochem, only: EvoAtmosphere, version, dp
  implicit none
  character(:), allocatable :: err
  type(EvoAtmosphere) :: pc
  logical :: success
  integer :: i, j
  real(dp), allocatable :: t_eval(:)
  real(dp) :: from, to
  
  print*,'photochem version == ',trim(version)

  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarthEvo.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarthEvo.txt", &
               err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

  allocate(t_eval(100))
  from = 0
  to = 1
  do i = 1,size(t_eval)
    t_eval(i) = from + (to - from)*(i-1)/(size(t_eval)-1)
  enddo
  t_eval = 10.0_dp**t_eval

  
  success = pc%evolve('../test1.dat',0.0_dp,pc%var%usol_init, t_eval, .true., err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif



end program
program wavetest
  use theorywaves
  
  implicit none

  real*8   :: u10      ! 10 meter wind (m/s)
  real*8   :: ustar    ! water-side surface friction velocity (m/s)
  real*8   :: hbl      ! boundary layer depth (m) 

  real*8   :: EF       ! Enhancement Factor
  real*8   :: ustokes  ! Surface layer averaged Stokes drift

  write(*,*) 'Enter wind speed at 10m (m/s): '
  read(*,*) u10

  write(*,*) 'Enter water-side surface friction velocity (m/s): '
  read(*,*) ustar

  write(*,*) 'Enter boundary layer depth (m): '
  read(*,*) hbl

  write(*,*) 'Wind Speed (m/s): ', u10
  write(*,*) 'Friction Velocity (m/s): ', ustar
  write(*,*) 'Boundary Layer Depth (m): ', hbl

  EF = EFactor_model(u10, ustar, hbl)
  
  write(*,*) 'Enhancement Factor: ', EF

  ustokes = ustokes_SL_model(u10, hbl)

  write(*,*) 'Surface Layer Averaged Stokes Drift: ', ustokes 
 
end program wavetest

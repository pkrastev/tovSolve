!=====================================================================
! MAIN PROGRAM
!=====================================================================
program tov2019
  use eos
  use constants
  use terminate
  use profile
  implicit none

  real(8)           :: rhoc
  real(8)           :: rho_start
  real(8)           :: rho_end
  real(8)           :: rho_tmp
  real(8)           :: rho_h
  integer(4)        :: nsteps
  integer(4)        :: i
  character(len=30) :: eos_file

  external tov_derivs, y_deriv, rkqc

  ! +++ Parameters +++
  eos_file  = 'eos_MDI_x0.0.in' ! EOS  
  rho_start = 0.09d0            ! Starting density
  rho_end   = 1.5d0             ! Ending density
  nsteps    = 100               ! Number of density steps
  ! +++ End Parameters +++

  rho_h = (rho_end - rho_start) / dfloat( nsteps - 1 )

  rho_tmp = rho_start

100 format(1x,t10,a,t22,a,t35,a,t46,a,t62a,t73,a,t86,a)
  write(6,100) 'M', 'R', 'k2', 'lambda', 'I', 'beta', 'rhoc'

  ! +++ Loop over central density +++
  do i = 1, nsteps
    rhoc = rho_tmp * fm3cm3
    call solve_tov(eos_file, rhoc)
    rho_tmp = rho_tmp + rho_h
  enddo

end program tov2019

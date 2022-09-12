!=====================================================================
! Modules
!.....................................................................
module constants
  implicit none
  ! Mathematical constants............................................
  real(8), parameter :: pi       = 3.1415926535897932384d0
  real(8), parameter :: eulercon = 0.577215664901532861d0   
  real(8), parameter :: a2rad    = pi/180.0d0
  real(8), parameter :: rad2a    = 180.0d0/pi
  ! Physical constants (cgs)..........................................
  real(8), parameter :: c       = 2.99792458d10
  real(8), parameter :: g       = 6.67259d-8
  real(8), parameter :: h       = 6.6260755d-27 
  real(8), parameter :: hbar    = 0.5d0 * h / pi
  real(8), parameter :: qe      = 4.8032068d-10   
  real(8), parameter :: avo     = 6.0221367d23 
  real(8), parameter :: kerg    = 1.380658d-16
  real(8), parameter :: kev     = 8.617385d-5 
  real(8), parameter :: amu     = 1.6605402d-24
  real(8), parameter :: mn      = 1.6749286d-24   
  real(8), parameter :: mp      = 1.6726231d-24   
  real(8), parameter :: me      = 9.1093897d-28   
  real(8), parameter :: rbohr   = hbar * hbar / ( me * qe * qe )
  real(8), parameter :: fine    = qe * qe / (hbar * c )  
  real(8), parameter :: hion    = 13.605698140d0  
  real(8), parameter :: ev2erg  = 1.602d-12
  real(8), parameter :: ssol    = 5.67051d-5
  real(8), parameter :: asol    = 4.0d0 * ssol / c 
  real(8), parameter :: weinlam = h * c / ( kerg * 4.965114232d0 ) 
  real(8), parameter :: weinfre = 2.821439372d0 * kerg / h   
  real(8), parameter :: rhonuc  = 2.342d14 
  ! Astrophysical constants...........................................
  real(8), parameter :: r_sun   = 6.95997d10
  real(8), parameter :: m_sun   = 1.9892d33
  real(8), parameter :: lsol    = 3.8268d33 
  real(8), parameter :: mearth  = 5.9764d27 
  real(8), parameter :: rearth  = 6.37d8 
  real(8), parameter :: ly      = 9.460528d17   
  real(8), parameter :: pc      = 3.261633d0 * ly 
  real(8), parameter :: au      = 1.495978921d13
  real(8), parameter :: secyer  = 3.1558149984d7
  ! Conversion factors................................................
  real(8), parameter :: fm3cm3      = 1.0d+39
  real(8), parameter :: mevfm3gcm3  = ( 1.0d0/1.7827d0 )*1.0d-12
  real(8), parameter :: dycm2mevfm3 = ( 197.3269788d0/3.1616d0 )*1.0d-35
  real(8), parameter :: hbar_c      = 197.3269788d0
  real(8), parameter :: pG          = (g / (c*c*c*c))*1.0d4  ! Eden, P in [1/m2]
  real(8), parameter :: mG          = (g/(c*c))*1.0d-2       ! Mass in [m]
end module constants

! MODULE EOS..........................................................
module eos
  implicit none
  integer(4)           :: np       ! Number of rows in EOS file
  real(8), allocatable :: eray(:)  ! Energy density
  real(8), allocatable :: pray(:)  ! Pressure
  real(8), allocatable :: hray(:)  ! Enthalpy
  real(8), allocatable :: xnray(:) ! Number density
end module eos

! Termination energy density and preassure............................
module terminate
  implicit none
  real(8) :: eden_term
  real(8) :: p_term
end module terminate

! NS profile..........................................................
module profile
  implicit none
  integer(4)           :: nr       ! Dimension of profile arrays
  real(8), allocatable :: ns_ed(:) ! Energy density
  real(8), allocatable :: ns_pr(:) ! Pressure
  real(8), allocatable :: ns_gm(:) ! Gravitational mass
  real(8), allocatable :: ns_r(:)  ! R
end module profile

!=====================================================================
! SUBROUTINE SOLVE_TOV: Compute spherical NS model
!                       {M, R, k2, lambda, I, beta, rhoc}
!=====================================================================
subroutine solve_tov(file_in, rhoc)
  use eos
  use constants
  use terminate
  use profile
  implicit none

  integer(4)            :: i
  integer(4)            :: kount
  integer(4)            :: ngood
  integer(4)            :: nbad
  integer(4), parameter :: ydim = 2     ! Number of equations to solve
  integer(4), parameter :: xdim = 1000  ! Number of rows for storage array
  real(8)               :: p_c          ! Central pressure
  real(8)               :: rhoc         ! Central number density
  real(8)               :: edenc        ! Central energy density
  real(8)               :: rho
  real(8)               :: eden
  real(8)               :: xp(xdim)     ! Storage for R
  real(8)               :: yp(ydim,xdim)! Storage for M, P
  real(8)               :: bc(ydim)
  real(8)               :: h_try
  real(8)               :: h_min
  real(8)               :: h_start
  real(8)               :: h_stop
  real(8)               :: eps
  real(8)               :: dxsav
  character(len=30)     :: file_in      ! Input file
  character(len=30)     :: op_eos
  character(len=30)     :: op_den
  character(len=30)     :: str_den

  real(8) :: yr
  real(8) :: cs2
  real(8) :: beta
  real(8) :: k2
  real(8) :: M
  real(8) :: R
  real(8) :: lambda

  real(8) :: fr   ! phi at r = R (Eq. (12) from arxiv: 1810.10992)
  real(8) :: i_ns ! Moment of inertia I

  external tov_derivs, y_deriv, rkqc
  
  ! Read EOS file.....................................................
  call load_eos(file_in)

  ! Get central pressure and energy density...........................
  call eosinv(rhoc, xnray, pray, np, p_c)
  call eosinv(p_c, pray, eray, np, edenc)

  ! Integrate in radius (in cm) from h_start to h_stop................
  h_start = 1.0d2
  h_try   = h_start
  h_min   = 1.0d0
  h_stop  = r_sun

  ! Initial conditions................................................ 
  bc(1) = (4.0d0/3.0d0)*pi*(h_start**3)*edenc                ! Gravitational mass
  bc(2) = p_c-0.5d0*(4.0d0/3.0d0*pi)*g*(h_start**2)*edenc**2 ! Central pressure

  ! Get termination density and pressure..............................
  eden_term = 10.d0**eray(2)
  call eosinv(eden_term, eray, pray, np, p_term)

  ! Set parameters for integration of ODEs............................
  eps   = 1.0e-6
  dxsav = 0.0d0

  ! Integrate the TOV equations in tov_derivs.........................
  call tovint(bc, ydim, xdim, h_start, h_stop, xdim, kount, dxsav, &
       xp, yp, eps, h_try, h_min, ngood, nbad, tov_derivs, rkqc)

  ! Allocate NS profile arrays.......................................
  if ( .not. allocated (ns_ed) )  allocate(ns_ed(kount) )
  if ( .not. allocated (ns_pr) )  allocate(ns_pr(kount) )
  if ( .not. allocated (ns_gm) )  allocate(ns_gm(kount) )
  if ( .not. allocated (ns_r) )   allocate(ns_r(kount) )

  ! Fill in NS profile arrays ........................................

  ! convert to geometric units
  ! This is required for computing k2 and the moment of inertia

  do i = 1, kount

     call eosinv( yp(2,i), pray, eray, np, eden ) ! get eden
     call eosinv( yp(2,i), pray, xnray, np, rho ) ! get rho
     
     ns_ed(i)  = log10(eden*(c*c*pG)) ! Energy density (1/m2)
     ns_pr(i)  = log10(yp(2,i)*pG)    ! Pressure (1/m2)
     ns_gm(i)  = log10(yp(1,i)*mG)    ! Gravitational mass (m)
     ns_r(i)   = log10(xp(i)/1.0d2)   ! Radius (m)

  end do

  nr = kount

  ! Calculate yR ....................................................
  call solve_y(yr)

  ! Calculate fR ....................................................
  call solve_f(fr)
  
  ! Calculate k2, lambda and I ......................................
  M    = yp(1,kount)                                  ! Gravitational mass (g)
  R    = xp(kount)                                    ! Radius (cm)
  beta = ( g*M ) / ( R*c*c )                          ! Compactness parameter (dimentionless)
  call solve_k2( yr, beta, k2 )                       ! k2 (dimentionless)
  lambda = (2.0d0 / (3.0d0*g))*k2*(R**5.0d0)          ! lambda (cm2 g s2)
  i_ns = (((R**3.0d0)*fr)/(6.0d0+2.0d0*fr)) * (c*c/g) ! Moment of inertia

  ! Write out results ................................................
  write(6,'(7(2x,f11.6))') &
        M/m_sun,           & ! Mass, M 
        R/1.0d5,           & ! Radius, R
        k2,                & ! Love number, k2
        lambda/1.0d36,     & ! Tidal deformability, lambda
        i_ns/1.0d45,       & ! Moment of inertia, I
        beta,              & ! Compactness, beta
        rhoc/fm3cm3          ! Central number density, rho_c

  ! Allocate NS profile arrays.......................................
  if ( allocated (ns_ed) )  deallocate(ns_ed)
  if ( allocated (ns_pr) )  deallocate(ns_pr)
  if ( allocated (ns_gm) )  deallocate(ns_gm)
  if ( allocated (ns_r) )   deallocate(ns_r)


  return
end subroutine solve_tov

!=====================================================================
! Load EOS file
! Format: energy_density pressure in CGS units (after RNS)
!=====================================================================
subroutine load_eos(eos_file)
  use constants
  use eos
  implicit none
  integer(4)        :: i
  real(8)           :: p        ! Pressure
  real(8)           :: eden     ! Energy density
  real(8)           :: h0       ! Enthalpy
  real(8)           :: n0       ! Number density
  character(len=30) :: eos_file ! EOS file  

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Read EOS file
  ! {P, Energy_Density, H, Number_Density}
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  open(unit=1,file=eos_file)
  read(1,*) np
  if ( .not. allocated(pray) ) allocate( pray(np) )
  if ( .not. allocated(eray) ) allocate( eray(np) )
  if ( .not. allocated(hray) ) allocate( hray(np) )
  if ( .not. allocated(xnray) ) allocate( xnray(np) )
  do i = 1, np
     read(1,*) eden, p, h0, n0
     eray(i)  = log10(eden)
     pray(i)  = log10(p)
     hray(i)  = log10(h0)
     xnray(i) = log10(n0)
  end do
  close(unit=1)

  return
end subroutine load_eos

!=====================================================================
! Calculate speed of sound
! cs2
!=====================================================================
subroutine speed_of_sound(r, cs2)
  use constants
  use eos
  use profile
  implicit none
  integer(4) :: i
  real(8)    :: r
  real(8)    :: rp
  real(8)    :: rm
  real(8)    :: pp
  real(8)    :: pm
  real(8)    :: edp
  real(8)    :: edm
  real(8)    :: cs2
  real(8), parameter :: rh = 1.0d0 ! step in meters

  rp = r + rh
  rm = r - rh

  call eosinv(rp, ns_r, ns_ed, nr, edp)
  call eosinv(rm, ns_r, ns_ed, nr, edm)
  call eosinv(rp, ns_r, ns_pr, nr, pp)
  call eosinv(rm, ns_r, ns_pr, nr, pm)

  ! Simple 2-point derivative approxiamtion...........................
  cs2 = ( pp - pm ) / (edp - edm )

  return
end subroutine speed_of_sound

!=====================================================================
! SUBROUTINE EOSINV: Interface to LANGINT
!                    Returns pressure, energy density, or 
!                    number density
!=====================================================================
subroutine eosinv(xx, xi, yi, ni, yy)
  implicit none
  integer(4), parameter :: order = 2 ! Order of interpolation
  integer(4)            :: ni        ! Size of arrays xi() and yi()
  real(8)               :: xx        ! Abscissa at which the interpolation is to be evaluated
  real(8)               :: yy        ! Interpolated value
  real(8)               :: xi(ni)    ! Arrays of data abscissas
  real(8)               :: yi(ni)    ! Arrays of data ordinates
  real(8), external     :: lagint
  real(8)               :: xl
  real(8)               :: yl
  xl = log10(xx)
  yl = lagint(xl, xi, yi, ni, order+1)
  yy = 10.0d0**yl
  return
end subroutine eosinv
 
!=====================================================================
! FUNCTION LAGINT: Lagrange interpolation
!
! Input:
! xx    -- abscissa at which the interpolation is to be evaluated
! xi()  -- arrays of data abscissas
! yi()  -- arrays of data ordinates
! ni    -- size of the arrays xi() and yi()
! n     -- number of points for interpolation (order of interp. = n-1)
!
! Output:
! lagint- interpolated value
!
! Comments:
! if ( n > ni ) n = ni
! Program works for both equally and unequally spaced xi()
!=====================================================================
function lagint(xx, xi, yi, ni, n)
  implicit none
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: js
  integer(4) :: jl
  integer(4) :: ni
  integer(4) :: n
  real(8)    :: xi(ni)
  real(8)    :: yi(ni)
  real(8)    :: lambda(ni)
  real(8)    :: y
  real(8)    :: lagint
  real(8)    :: xx

  ! check order of interpolation
  if ( n > ni ) n = ni

  ! if x is ouside the xi(1)-xi(ni) interval take a boundary value
  if ( xx <= xi(1) ) then
     lagint = yi(1)
     return
  end if
  if ( xx >= xi(ni) ) then
     lagint = yi(ni)
     return
  end if

  ! a binary (bisectional) search to find i so that xi(i) < x < xi(i+1)
  i = 1
  j = ni
  do while ( j > i+1 )
     k = ( i + j ) / 2
     if ( xx < xi(k) ) then
        j = k
     else
        i = k
     end if
  end do

  ! shift i that will correspond to n-th order of interpolation
  ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
  i = i + 1 - n/2

  ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
  if ( i < 1 ) i = 1
  if ( i + n > ni ) i = ni - n + 1

  ! Lagrange interpolation
  y = 0.0
  do js =i, i + n - 1
     lambda(js) = 1.0
     do jl = i, i + n - 1
        if( jl /= js ) lambda(js)=lambda(js)*(xx-xi(jl))/(xi(js)-xi(jl))
     end do
     y = y + yi(js)*lambda(js)
  end do
  lagint = y
  return
end function lagint

!=====================================================================
! SUBROUTINE TOV_DERIVS: Setting up the Tolman-Oppenheimer-Volkoff 
!                        equations for spherically symmetric object 
!                        in hydrostatic equilibrium
!                        (cgs)
!=====================================================================
subroutine tov_derivs(r, y, dydr)
  use constants
  use eos
  use terminate
  implicit none
  real(8) :: r
  real(8) :: y(2)
  real(8) :: dydr(2)
  real(8) :: grmass
  real(8) :: eden
  real(8) :: p
  real(8) :: relcor
  real(8), parameter :: c2 = c*c

  !...................................................................
  ! y(1): gravitational mass
  ! y(2): pressure
  !...................................................................

  grmass = y(1)
  p      = y(2)

  ! Get energy density and number density.............................
  call eosinv(p, pray, eray, np, eden)

  ! Gravitational mass................................................
  dydr(1) = (4.0d0*pi)*r**2*eden

  ! Pressure.......................................................... 
  relcor = (1.0d0+p/(eden*c2))*                       &
       ( 1.0d0+((4.0d0*pi)*p*r**3 )/( grmass*c2) )/   &
       ( 1.0d0-(2.0d0*g*grmass )/( r*c2) )
  dydr(2) = -( g*grmass/(r**2) ) * eden * relcor

  return
end subroutine tov_derivs

!=====================================================================
! SUBROUTINE Y_DERIV: Set up equation for y to compute k2
!                     ( c = G = 1 )
!=====================================================================
subroutine y_deriv(r, yr, dydr)
  use constants
  use profile
  implicit none
  real(8) :: r
  real(8) :: gm
  real(8) :: p
  real(8) :: ed
  real(8) :: cs2
  real(8) :: y
  real(8) :: yr(1)
  real(8) :: dydr(1)
  real(8) :: F
  real(8) :: Q
  real(8) :: L
  real(8) :: r2, r3, r4, gm2
  y = yr(1)
  call eosinv(r, ns_r, ns_ed, nr, ed)
  call eosinv(r, ns_r, ns_pr, nr, p)
  call eosinv(r, ns_r, ns_gm, nr, gm)
  call speed_of_sound(r, cs2)
  r2 = r*r     ! r2
  r3 = r*r*r   ! r3
  r4 = r*r*r*r ! r4
  gm2 = gm*gm  ! gm2
  L = 1.0d0 / ( 1.0d0 - 2.0d0*gm/r )
  F = ( 1.0d0 - (4.0d0*pi*r2)*( ed - p ) ) * L
  Q = (4.0d0*pi)*( 5.0d0*ed + 9.0d0*p + ((ed+p)/cs2) ) * L &
       -(6.0d0/r2)* L - (4.0d0*gm2/r4)                     &
       *( ( 1.0d0 + 4.0d0*pi*r3*p/gm )**2.0d0 ) * (L**2.0d0)
  dydr(1) = -(y*y/r) -(y*F/r) -(r*Q)
  return
end subroutine y_deriv

!=====================================================================
! SUBROUTINE SOLVE_Y: Integrate the ODE in y_deriv and compute yR
!                     Use ode.f90 library
!=====================================================================
subroutine solve_y(yr) 
  use constants 
  use profile
  implicit none
  integer(4), parameter :: neqn = 1 
  integer(4) :: i 
  integer(4) :: iflag
  integer(4) :: iwork(5) 
  real(8)    :: yr 
  real(8)    :: relerr 
  real(8)    :: abserr 
  real(8)    :: ri 
  real(8)    :: rf 
  real(8)    :: work(100+21*neqn) 
  real(8)    :: y(neqn) 
  external y_deriv

  abserr = 0.00001d0
  relerr = 0.00001d0

  iflag = 1

  ri = 10.0d0**ns_r(1)  ! Initial radius
  rf = 10.0d0**ns_r(nr) ! Final radius

  y(1) = 2.0d0          !BC

  ! Integrate from ri to rf ( center to edge of NS )..................
  call ode( y_deriv, neqn, y, ri, rf, relerr, abserr, iflag, work, iwork )

  yr = y (1) ! yR 

  return
end subroutine solve_y

!=====================================================================
! SUBROUTINE SOLVE_K2: Evaluate Love number k2
!=====================================================================
subroutine solve_k2(yr, beta, k2)
  implicit none
  real(8) :: beta
  real(8) :: yR
  real(8) :: k2
  real(8) :: num
  real(8) :: den
  num = (8.0d0/5.0d0)*(beta**5)*((1.0d0-2.0d0*beta)**2)                      &
       *(2.0d0-yR+2.0d0*beta*(yR-1.0d0))
  den =  2.0d0*beta*(6.0d0-3.0d0*yR+3.0d0*beta*(5.0d0*yR-8.0d0))+            &
         4.0d0*beta*beta*beta*(13.0d0-11.0d0*yR+beta*(3.0d0*yR-2.0d0)+2.0d0* &
         beta*beta*(1.0d0+yR))+                                              &
         3.0d0*((1.0d0-2.0d0*beta)**2)*(2.0-yR+2.0*beta*(yR-1.0))*           &
         log(1.0-2.0*beta)
  k2 = num / den
  return
end subroutine solve_k2

!=====================================================================
! SUBROUTINE F_DERIV
!=====================================================================
subroutine f_deriv(r, fr, dfdr)
  use constants
  use profile
  implicit none
  real(8) :: r
  real(8) :: gm
  real(8) :: p
  real(8) :: ed
  real(8) :: f
  real(8) :: fr(1)
  real(8) :: dfdr(1)
  real(8) :: r2
  f = fr(1)
  call eosinv(r, ns_r, ns_ed, nr, ed)
  call eosinv(r, ns_r, ns_pr, nr, p)
  call eosinv(r, ns_r, ns_gm, nr, gm)
  r2 = r*r ! r2
  dfdr(1) = -(f/r)*(f+3.0d0) + (4.0d0+f)*(4.0d0*pi*r2)*(ed+p)/(r-2.0d0*gm) 
  return
end subroutine f_deriv

!=====================================================================
! SUBROUTINE: SOLVE_F
!=====================================================================
subroutine solve_f(fr)
  use constants
  use profile
  implicit none
  integer(4), parameter :: neqn = 1
  integer(4) :: i
  integer(4) :: iflag
  integer(4) :: iwork(5)
  real(8)    :: fr
  real(8)    :: relerr
  real(8)    :: abserr
  real(8)    :: ri
  real(8)    :: rf
  real(8)    :: work(100+21*neqn)
  real(8)    :: f(neqn)
  external f_deriv

  abserr = 0.00001d0
  relerr = 0.00001d0

  iflag = 1

  ri = 10.0d0**ns_r(1)  ! Initial radius
  rf = 10.0d0**ns_r(nr) ! Final radius

  f(1) = 0.000001d0     ! BC

  ! Integrate from ri to rf ( center to edge of NS )..................
  call ode( f_deriv, neqn, f, ri, rf, relerr, abserr, iflag, work, iwork )

  fr = f(1) ! fR 
  return
end subroutine solve_f

!=====================================================================
! SUBROUTINE RKQC: Adapted from Numerical Recipes
!                  (C) Copr. 1986-92 Numerical Recipes Software
!=====================================================================
subroutine rkqc(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)   
  implicit none    
  integer(4)            :: n
  integer(4)            :: i
  integer(4), parameter :: nmax = 2000     
  real(8)               :: x
  real(8)               :: htry
  real(8)               :: eps
  real(8)               :: hdid
  real(8)               :: hnext
  real(8)               :: h
  real(8)               :: hh
  real(8)               :: xsav
  real(8)               :: errmax
  real(8)               :: y(n)
  real(8)               :: dydx(n)
  real(8)               :: yscal(n)  
  real(8)               :: ytemp(nmax)
  real(8)               :: ysav(nmax)
  real(8)               :: dysav(nmax)
  real(8), parameter    :: fcor   = 1.0d0/15.0d0
  real(8), parameter    :: safety = 0.9d0
  real(8), parameter    :: errcon = 6.0d-4
  real(8), parameter    :: pgrow  = -0.2d0
  real(8), parameter    :: pshrnk = -0.25d0
  external derivs
 
  h = htry 
  xsav =  x 
  do i = 1, n    
     ysav(i)  = y(i) 
     dysav(i) = dydx(i) 
  end do
   
1 hh = 0.5d0 * h   
  call rk4(ysav, dysav, n, xsav, hh, ytemp, derivs)    
  x = xsav + hh  
  call derivs(x, ytemp, dydx)  
  call rk4(ytemp, dydx, n, x, hh, y, derivs)   
  x = xsav + h   
  if ( x == xsav ) stop 'stepsize not significant in rkqc'  
   
  call rk4(ysav, dysav, n, xsav, h, ytemp, derivs) 
  
  errmax = 0.0d0 
  do i = 1, n    
     ytemp(i) = y(i) - ytemp(i)   
     errmax = max(errmax,abs(ytemp(i)/yscal(i)))    
  end do
  errmax = errmax/eps  
 
  if (errmax > 1.0) then 
     h = safety * h * (errmax**pshrnk)  
     go to  1   
  else   
     hdid = h   
     if ( errmax > errcon ) then 
        hnext = safety * h * (errmax**pgrow) 
     else 
        hnext = 4.0d0 * h 
     end if
  end if
  do i = 1, n    
     y(i) = y(i) + ytemp(i)*fcor  
  end do
  return 
end subroutine rkqc

!=====================================================================
! SUBROUTINE RK4: Adapted from Numerical Recipes
!                 (C) Copr. 1986-92 Numerical Recipes Software
!=====================================================================
subroutine rk4(y, dydx, n, x, h, yout, derivs)  
  implicit none
  integer(4)            :: i
  integer(4)            :: n
  integer(4), parameter :: nmax = 200
  real(8)               :: x
  real(8)               :: h
  real(8)               :: y(n)
  real(8)               :: dydx(n)
  real(8)               :: yout(n)
  real(8)               :: yt(nmax)
  real(8)               :: dyt(nmax)
  real(8)               :: dym(nmax)
  real(8)               :: hh
  real(8)               :: h6
  real(8)               :: xh    
  external derivs

  hh = h * 0.5d0 
  h6 = h / 6.0d0  
  xh = x + hh   

  ! Step 1............................................................    
  do i = 1, n   
     yt(i) = y(i) + hh*dydx(i)   
  end do

  ! Step 2............................................................
  call derivs(xh, yt, dyt)    
  do i = 1, n   
     yt(i) = y(i) + hh*dyt(i)    
  end do

  ! Step 3............................................................
  call derivs(xh, yt, dym)    
  do i = 1, n   
     yt(i)  = y(i) + h*dym(i) 
     dym(i) = dyt(i) + dym(i)    
  end do

  ! Step 4............................................................
  call derivs(x+h, yt, dyt)   
  do i = 1, n   
     yout(i) = y(i) + h6*(dydx(i) +dyt(i) + 2.0d0*dym(i))  
  end do
  return    
end subroutine rk4

!=====================================================================
! SUBROUTINE TOVINT: Adapted from Numerical Recipes
!                    (C) Copr. 1986-92 Numerical Recipes Software
!=====================================================================
subroutine tovint(ystart, nvar, xdim, x1, x2, kmax, kount, dxsav, &
     xp, yp, eps, h1, hmin, nok, nbad, derivs, steper) 

  use terminate
  implicit none

  integer(4)            :: nok
  integer(4)            :: nbad
  integer(4)            :: kmax
  integer(4)            :: kount
  integer(4)            :: i
  integer(4)            :: nstp
  integer(4)            :: nvar
  integer(4)            :: xdim
  integer(4), parameter :: nmax = 20
  integer(4), parameter :: maxstp = 20000
  real(8)               :: x1
  real(8)               :: x2
  real(8)               :: h1
  real(8)               :: hmin
  real(8)               :: eps
  real(8)               :: dxsav
  real(8)               :: x
  real(8)               :: xsav
  real(8)               :: h
  real(8)               :: hdid
  real(8)               :: hnext
  real(8)               :: xp(xdim)
  real(8)               :: yp(nvar,xdim)
  real(8)               :: yscal(nmax)
  real(8)               :: y(nmax)
  real(8)               :: ystart(nvar)
  real(8)               :: dydx(nmax)
  real(8), parameter    :: zero = 0.0d0
  real(8), parameter    :: one = 1.0d0
  real(8), parameter    :: tiny = 1.0d-15

  external derivs, steper 

  x     = x1 
  h     = sign(h1,x2-x1) 
  nok   = 0 
  nbad  = 0 
  kount = 0

  ! Store the first step..............................................
  do i = 1, nvar 
     y(i) = ystart(i) 
  end do
  xsav = x - 2.0d0 * dxsav

  ! Do at most maxstp.................................................
  do nstp = 1, maxstp 
     call derivs(x, y, dydx)
     
     ! Scaling vector used to monitor accuracy........................
     do i = 1, nvar 
        yscal(i) = abs(y(i)) + tiny 
     end do

     ! Store intermediate results.....................................
     if ( kmax > 0 ) then 
        if ( (abs(dxsav) - abs(x-xsav)) <= tiny ) then
           if ( kount < (kmax-1) )  then 
              kount = kount + 1 
              xp(kount) = x 
              do i = 1, nvar 
                 yp(i,kount) = y(i) 
              end do
              xsav = x 
           end if
        end if
     end if

     ! If the step can overshoot the stop point or the dxsav increment then cut it
     if ( (x+h-x2)*(x+h-x1) > zero ) h = x2 - x
     if ( dxsav .ne. zero .and. h > (xsav-x+dxsav) ) h = xsav + dxsav-x
 
     ! Take an integration step.......................................
     call steper(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs) 
     if ( hdid == h ) then 
        nok = nok + 1 
     else 
        nbad = nbad + 1 
     end if

     ! Exit point: save the final step................................
     if ( (x-x2)*(x2-x1) >= zero .or. y(2) <= p_term ) then 
        do i = 1, nvar 
           ystart(i) = y(i) 
        end do
        if ( kmax /= 0 ) then 
           kount = kount + 1 
           xp(kount) = x 
           do i = 1, nvar 
              yp(i,kount) = y(i) 
           end do
        end if
        return 
     end if

     !Set the step size for the next iteration........................
     h = hnext
     if ( abs(hnext) < hmin ) return 
     
  end do

  return 
end subroutine tovint

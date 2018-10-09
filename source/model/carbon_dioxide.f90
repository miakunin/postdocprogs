MODULE CARBON_DIOXIDE_MOD

contains
!> Subroutine CARBON_DIOXIDE calculates :
!! 1. Diffusion of dissolved inorganic carbon (DIC), DOC, POCL, POCD, in the water column
!! 2. Settlement of POCD
!! 3. Gas exchange of CO_2 between solution and bubbles.
!! Note, thate \f$DIC\f$ array stands for DIC concentration (mol/m**3 of C atomes in DIC).
SUBROUTINE CARBON_DIOXIDE &
& (gs, Twater, fbbleflx_co2_sum, fbbleflx_co2, lamcarbdi, gsp, bathymwater, lsc, fracb0, &
& pressure, wind10, co2_pres_atm, febul0, ls, dt, eps_surf, sodbot, gas, soilswitch, Flux_atm)

use LAKE_DATATYPES, only : ireals, iintegers

use PHYS_CONSTANTS, only : Kelvin0

use ARRAYS, only : gas_type, workc, &
& water_carbdi_indic, water_carbon_indic
use ARRAYS_GRID, only : gridspacing_type, gridsize_type
use ARRAYS_BATHYM, only : bathym, layers_type
use ARRAYS_SOIL, only : lsh_type

use METH_OXYG_CONSTANTS, only : CO2O2_sod, ngasb, &
& Henry_const0_co2, Henry_temp_dep_co2, Henry_temp_ref, &
& pH, dpart

use PHYS_FUNC, only : &
& HENRY_CONST, &
& GAS_WATATM_FLUX, HC_CORR_CARBEQUIL, &
& W_SEDIM

use T_SOLVER_MOD, only : &
& DIFF_COEF

use NUMERICS, only : PROGONKA

use NUMERIC_PARAMS, only : vector_length

use DRIVING_PARAMS, only : carbon_model

implicit none

!> Grid spacing group
type(gridspacing_type), intent(in) :: gsp
!> Gridsizes group
type(gridsize_type),    intent(in) :: gs
!> Depths of physical layers 
type(layers_type),      intent(in) :: ls
!> Bathymetry characteristics of water body
type(bathym),           intent(in) :: bathymwater(1:gs%M+1)
!> Source term due to marginal vertical \f$CO_2\f$ flux into sediments
type(lsh_type),         intent(in) :: lsc
!> Group of dissolved gases
type(gas_type),         intent(inout) :: gas
!> Water temperature profile, deg.Celsius
real(kind=ireals), intent(in) :: Twater(1:gs%M+1) 
!> The bubble carbon-dioxide flux, horizontally integrated, 
!! normalized by bottom value, n/d
real(kind=ireals), intent(inout) :: fbbleflx_co2_sum(0:gs%M+1)
!> The bubble carbon-dioxide flux above the deepest soil column, 
!! normalized by bottom value, n/d
real(kind=ireals), intent(in) :: fbbleflx_co2(0:gs%M+1) 
!> Vertical diffusivity for \f$CO_2\f$
real(kind=ireals), intent(in) :: lamcarbdi(1:gs%M)
!> Molar fractions of gases in a bubble
real(kind=ireals), intent(in) :: fracb0(1:ngasb) 
!> Atmospheric pressure, Pa
real(kind=ireals), intent(in) :: pressure 
!> Wind speed at 10 m above the surface
real(kind=ireals), intent(in) :: wind10 
!> Carbon dioxide partial pressure in the atmosphere, Pa
real(kind=ireals), intent(in) :: co2_pres_atm 
!> Methane bottom bubble flux
real(kind=ireals), intent(in) :: febul0 
!> Time step, s
real(kind=ireals), intent(in) :: dt
!> TKE disipation rate at the surface, m**2/s**3
real(kind=ireals), intent(in) :: eps_surf
!> Sedimentary oxygen demand at the bottom
real(kind=ireals), intent(in) :: sodbot
!> Soil switch
integer(kind=iintegers), intent(in) :: soilswitch
!> Diffusive flux of \f$CO_2\f$ to the atmosphere
real(kind=ireals), intent(out) :: Flux_atm

! Local variables
real(kind=ireals), parameter :: small_value = 1.d-20

real(kind=ireals), allocatable :: a(:), b(:), c(:), f(:), y(:)
real(kind=ireals), allocatable, target :: conc(:)

real(kind=ireals) :: Fcarbdi1 ! the flux of carbon dioxide at the bottom

real(kind=ireals) :: x, xx, xxx, xxxx
integer(kind=iintegers) :: i, nspec

allocate (a(1:vector_length),b(1:vector_length),c(1:vector_length), &
&         f(1:vector_length),y(1:vector_length))
allocate (conc(1:gs%M+1))

Fcarbdi1 = - sodbot*CO2O2_sod

if (ls%h1 == 0) return ! Nothing to do if there is no water layer

if     (carbon_model%par == 1) then
  nspec = 1 ! No need to calculate DOC, POCL, POCD
elseif (carbon_model%par == 2) then
  nspec = 5
endif

!Vertical diffusion of DIC, DOC, POCL, POCD.
!Diffusivities are assumed to be the same.
do i = 1, nspec

  select case(i)
  case (1)
    conc = gas%DIC(:,1)
  case (2)
    conc = gas%DOC(:,1) !autochtonous DOC
  case (3)
    conc = gas%POCL
  case (4)
    conc = gas%POCD
  case (5)
    conc = gas%DOC(:,2) !allochtonous DOC
  end select

  !i = 1: DIC
  !i = 2: DOC (autochtonous)
  !i = 3: POCL
  !i = 4: POCD
  !i = 5: DOC (allochtonous)

  ! 1-st step of splitting-up scheme - diffusion
  if (ls%l1 == 0 .and. i == 1 ) then !No diffusion to the atmosphere for DOC, POCL, POCD
    !c(1) = 1.
    !b(1) = 0.
    !f(1) = co2_pres_atm * &
    !& HENRY_CONST(Henry_const0_co2, Henry_temp_dep_co2, Henry_temp_ref, Twater(1)+273.15) !co2_atm0 
    x = HENRY_CONST(Henry_const0_co2, Henry_temp_dep_co2, Henry_temp_ref, Twater(1)+Kelvin0)
    xx = gas%DIC(1,1)/HC_CORR_CARBEQUIL(Twater(1)+Kelvin0,pH)! Converting DIC to CO_2(aq)
    Flux_atm = GAS_WATATM_FLUX &
    & (Twater(1),wind10,xx,co2_pres_atm,x,water_carbdi_indic,eps_surf)
  else
    Flux_atm = 0.
  endif
  ! Tridiagonal system coefficients for the top layer
  if (i == 1) then
    ! Morphometry effects for DIC
    xxx  = bathymwater(1)%area_half/bathymwater(1)%area_int
    xxxx = dt*lsc%water(1)
  else
    ! No morphometry effects for DOC, POCL, POCD
    xxx  = 1.
    xxxx = 0.
  endif
  xx     = 0.5*( - xxx * lamcarbdi(1)/(ls%h1*gsp%ddz(1)) + 0.5*ls%dhw0/dt )
  x      = 0.5*gsp%ddz(1)*ls%h1/dt
  c(1)   = xx - x
  b(1)   = xx
  f(1)   = - x*conc(1) + Flux_atm + b(1)*conc(2) - xx*conc(1) - x*xxxx
  ! Specifying coefficients for internal layers
  if (i == 1) then
    call DIFF_COEF(a,b,c,f,2,gs%M,2,gs%M,water_carbdi_indic,dt)
  else
    workc => conc
    call DIFF_COEF(a,b,c,f,2,gs%M,2,gs%M,water_carbon_indic,dt)
    nullify(workc)
  endif
  ! Tridiagonal system coefficients for the bottom layer
  if (i == 1) then
    ! Morphometry effects only for DIC, otherwise xxx = 1, specified above
    xxx  = bathymwater(gs%M)%area_half/bathymwater(gs%M+1)%area_int
  endif
  xx     = 0.5*( - xxx * &
  & lamcarbdi(gs%M)/(gsp%ddz(gs%M)*ls%h1) + 0.5*(ls%dhw - ls%dhw0)/dt )
  x      = 0.5*gsp%ddz(gs%M)*ls%h1/dt
  c(gs%M+1) = xx - x
  a(gs%M+1) = xx
  f(gs%M+1) = - x*conc(gs%M+1) + a(gs%M+1)*conc(gs%M) - xx*conc(gs%M+1)
  ! Bottom diffusive flux only for DIC
  if (ls%ls1 == 0. .and. i == 1) f(gs%M+1) = f(gs%M+1) + Fcarbdi1
  call PROGONKA (a,b,c,f,y,1,gs%M+1)
  y(1:gs%M+1) = max(y(1:gs%M+1),0.d0)

  select case(i)
  case (1)
    gas%DIC   (1:gs%M+1,2) = y(1:gs%M+1)
  case (2)
    gas%DOC   (1:gs%M+1,1) = y(1:gs%M+1)
  case (3)
    gas%POCL  (1:gs%M+1)   = y(1:gs%M+1)
  case (4)
    gas%POCD  (1:gs%M+1)   = y(1:gs%M+1)
  case (5)
    gas%DOC   (1:gs%M+1,2) = y(1:gs%M+1)
  end select

enddo

if (carbon_model%par == 2) then
  !Settlement of POCD, explicit (Euler), central-differences scheme
  !> @todo: do different sedimentation speed depending on turbulence regime
  !> @todo: check scheme for sedimentation term in MyLake model
  x = W_SEDIM(dpart,1_iintegers) ! the same laminar-limit sedimentation speed at all depths
  conc(1) = gas%POCD(1) - dt*x*(gas%POCD(2) - gas%POCD(1))/(ls%h1*gsp%ddz(1))
  do i = 2, gs%M
    conc(i) = gas%POCD(i) - dt*x*(gas%POCD(i+1) - gas%POCD(i-1))/(ls%h1*gsp%ddz2(i-1))
  enddo
  conc(gs%M+1) = gas%POCD(gs%M+1) - dt*x*(gas%POCD(gs%M+1) - gas%POCD(gs%M))/(ls%h1*gsp%ddz(gs%M))
  gas%POCD(1:gs%M+1) = conc(1:gs%M+1)
endif

! Horizontally integrated CO_2 bubble flux in case with single soil column
if (gs%nsoilcols == 1) then
  fbbleflx_co2_sum(0)      = febul0*fracb0(2)/fracb0(1)*fbbleflx_co2(0)     *bathymwater(1)     %area_int
  fbbleflx_co2_sum(1:gs%M) = febul0*fracb0(2)/fracb0(1)*fbbleflx_co2(1:gs%M)*bathymwater(1:gs%M)%area_half
  fbbleflx_co2_sum(gs%M+1) = febul0*fracb0(2)/fracb0(1)*fbbleflx_co2(gs%M+1)*bathymwater(gs%M+1)%area_int
endif

if (soilswitch == 1) then
  ! 2-d step of splitting-up scheme - DIC source due to bubble dissolution
  !x = dt*febul0*fracb0(2)/fracb0(1)
  gas%DIC(1,2) = gas%DIC(1,2) + &
  & 2.d0*dt/(ls%h1*gsp%ddz(1))*(fbbleflx_co2_sum(1) - fbbleflx_co2_sum(0)) / &
  & bathymwater(1)%area_int
  gas%DIC(1,2) = max(gas%DIC(1,2),small_value)
  do i = 2, gs%M
    gas%DIC(i,2) = gas%DIC(i,2) + &
    & dt/(ls%h1*gsp%ddz05(i-1))*(fbbleflx_co2_sum(i) - fbbleflx_co2_sum(i-1)) / &
    & bathymwater(i)%area_int
    gas%DIC(i,2) = max(gas%DIC(i,2),small_value)
  enddo
  gas%DIC(gs%M+1,2) = gas%DIC(gs%M+1,2) + 2.d0*dt/(ls%h1*gsp%ddz(gs%M)) * &
  & (fbbleflx_co2_sum(gs%M+1) - fbbleflx_co2_sum(gs%M)) / &
  & bathymwater(gs%M+1)%area_int
  gas%DIC(gs%M+1,2) = max(gas%DIC(gs%M+1,2),small_value)
endif


deallocate (a, b, c, f, y)
deallocate(conc)

END SUBROUTINE CARBON_DIOXIDE

END MODULE CARBON_DIOXIDE_MOD

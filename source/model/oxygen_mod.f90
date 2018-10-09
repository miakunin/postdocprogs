MODULE OXYGEN_MOD

use NUMERICS, only : PROGONKA, STEP
use NUMERIC_PARAMS, only : vector_length
use LAKE_DATATYPES, only : ireals, iintegers
use ARRAYS_GRID, only : gridsize_type
use ARRAYS_BATHYM, only : layers_type, bathym 
use ARRAYS_GRID, only : gridspacing_type
use ARRAYS_SOIL, only : lsh_type
use METH_OXYG_CONSTANTS
use DRIVING_PARAMS, only : carbon_model
use UNITS, only : kg_to_g

real(4), parameter :: k_PA = 1.7 ! n/d, coefficient in Poole&Atkins formula
real(4), parameter :: extwat1 = k_PA/2., extwat2 = k_PA/3.5 ! Extinction coefficients
                                                          ! identifying eutrophic,
                                                          ! mesotrophic and oligotrophic lakes,
                                                          ! derived from Poole and Atkins (1929) 
                                                          ! formula and data from Stefan and Fang, 1994
real(4), parameter :: chla_aver_oligo = 2.E-3 !2.E-3 ! mg/l, Mean chlorophyll density in oligotrophic lake, Stefan and Fang, 1994
real(4), parameter :: chla_aver_meso  = 6.E-3 !6.E-3 ! mg/l, Mean chlorophyll density in mesotrophic lake, Stefan and Fang, 1994
real(4), parameter :: chla_aver_eu    = 15.E-3 !15.E-3 ! mg/l, Mean chlorophyll density in eutrophic lake, Stefan and Fang, 1994

contains
!> Subroutine OXYGEN calculates vertical diffusion of dissolved oxygen and interaction with bubbles
SUBROUTINE OXYGEN &
& (gs, Twater, fbbleflx_o2_sum, fbbleflx_o2, lamoxyg, gsp, fracb0, pressure, wind10, &
& o2_pres_atm, ls, bathymwater, lso, dt, febul0, eps_surf, sodbot, oxyg, soilswitch, Flux_atm)

use PHYS_FUNC, only : &
& HENRY_CONST, &
& GAS_WATATM_FLUX

use PHYS_CONSTANTS, only : &
& Kelvin0

use T_SOLVER_MOD, only : &
& DIFF_COEF

implicit none

!> Grid size group
type(gridsize_type),    intent(in) :: gs
!> Grid spacing group
type(gridspacing_type), intent(in) :: gsp
!> Group containing depths of physical layers
type(layers_type),      intent(in) :: ls
!> Bathymetry characteristics of a water body
type(bathym),           intent(in) :: bathymwater(1:gs%M+1)
!> Oxygen flux at the water-sediments interface
type(lsh_type),         intent(in) :: lso
!> Water temperature profile, Celsius
real(kind=ireals), intent(in) :: Twater(1:gs%M+1) 
!> The bubble oxygen flux, horizontally integrated, normalized by bottom value, n/d
real(kind=ireals), intent(inout) :: fbbleflx_o2_sum(0:gs%M+1)
!> The bubble oxygen flux above the deepest soil column, normalized by bottom value, n/d
real(kind=ireals), intent(in) :: fbbleflx_o2(0:gs%M+1) 
!> Vertical diffusivity for oxygen, m**2/s
real(kind=ireals), intent(in) :: lamoxyg(1:gs%M)
!> Molar fractions of gases in a bubble
real(kind=ireals), intent(in) :: fracb0(1:ngasb) 
!> Bottom methane flux, mol/m**2/s
real(kind=ireals), intent(in) :: febul0
!> Atmospheric pressure, Pa
real(kind=ireals), intent(in) :: pressure 
!> Wind speed at 10 m above the surface
real(kind=ireals), intent(in) :: wind10 
!> Oxygen partial pressure in the atmosphere, Pa
real(kind=ireals), intent(in) :: o2_pres_atm 
!> Time step, s
real(kind=ireals), intent(in) :: dt
!> TKE dissipation rate at the water surface, m**2/s**3
real(kind=ireals), intent(in) :: eps_surf 
!> Sedimentary oxygen demand (oxygen flux at the bottom), mol/(m**2*s)
real(kind=ireals), intent(in) :: sodbot 
!> Oxygen concentration in water, mol/m**3
real(kind=ireals), intent(inout) :: oxyg(1:gs%M+1,1:2)
!> Switch for soil (sediments)
integer(kind=iintegers), intent(in) :: soilswitch
!> Flux to the atmosphere, mol/m**2/s
real(kind=ireals), intent(out) :: Flux_atm

! Local variables
integer(kind=iintegers), parameter :: water_oxygen_indic = 9

real(kind=ireals), parameter :: small_value = 1.d-20

real(kind=ireals), allocatable :: a(:), b(:), c(:), f(:), y(:)

real(kind=ireals), save :: Foxyg1 = 0.d0 ! Bottom oxygen flux, mol/(m**2*s)

real(kind=ireals) :: x, xx
real(kind=ireals), save :: botO2 = 0.    ! Bottom oxygen concentration, mol/m**3

integer(kind=iintegers) :: i
integer(kind=iintegers), parameter :: switchbotbc = 1 !1 - Neumann, 2 - Dirichlet

allocate (a(1:vector_length),b(1:vector_length),c(1:vector_length), &
&         f(1:vector_length),y(1:vector_length))
a(:) = 0.; b(:) = 0.; c(:) = 0.; f(:) = 0.; y(:) = 0.

Foxyg1 = sodbot

if (ls%h1 > 0) then

! 1-st step of splitting-up scheme - diffusion
  if (ls%l1 == 0) then
!    c(1) = 1.
!    b(1) = 0.
!    f(1) = o2_pres_atm * &
!    & HENRY_CONST(Henry_const0_o2, Henry_temp_dep_o2, Henry_temp_ref, Twater(1)+273.15) !o2_atm0 
    x = HENRY_CONST(Henry_const0_o2, Henry_temp_dep_o2, Henry_temp_ref, &
    & Twater(1) + Kelvin0) 
    Flux_atm = GAS_WATATM_FLUX &
    & (Twater(1),wind10,oxyg(1,1),o2_pres_atm,x,water_oxygen_indic,eps_surf)
  else
    Flux_atm = 0.
  endif
  xx = 0.5* ( - bathymwater(1)%area_half/bathymwater(1)%area_int * &
  & lamoxyg(1)/(ls%h1*gsp%ddz(1)) + 0.5*ls%dhw0/dt)
  x  = 0.5*gsp%ddz(1)*ls%h1/dt
  c(1)   = xx - x
  b(1)   = xx
  f(1)   = - x*oxyg(1,1) + Flux_atm + b(1)*oxyg(2,1) - xx*oxyg(1,1) - &
  & x*dt*lso%water(1)
  call DIFF_COEF(a,b,c,f,2,gs%M,2,gs%M,water_oxygen_indic,dt)
  if (ls%ls1 > 0) then  
    ! Switching bottom b.c.
    select case (switchbotbc)
    case(1)
      xx = 0.5*( - bathymwater(gs%M)%area_half/bathymwater(gs%M+1)%area_int * &
      & lamoxyg(gs%M)/(gsp%ddz(gs%M)*ls%h1) + 0.5*(ls%dhw - ls%dhw0)/dt)
      x = 0.5*gsp%ddz(gs%M)*ls%h1/dt
      c(gs%M+1) = xx - x
      a(gs%M+1) = xx
      f(gs%M+1) = - x*oxyg(gs%M+1,1) + a(gs%M+1)*oxyg(gs%M,1) - xx*oxyg(gs%M+1,1)
    case(2)
      c(gs%M+1) = 1.
      a(gs%M+1) = 0.
      f(gs%M+1) = botO2
    end select
    call PROGONKA (a,b,c,f,y,1,gs%M+1)
    y(1:gs%M+1) = max(y(1:gs%M+1),0.d0)
    oxyg(1:gs%M+1,2) = y(1:gs%M+1)
  else  
    ! Switching bottom b.c.
    select case (switchbotbc)
    case(1)
      xx = 0.5*( - bathymwater(gs%M)%area_half/bathymwater(gs%M+1)%area_int * &
      & lamoxyg(gs%M)/(gsp%ddz(gs%M)*ls%h1) + 0.5*(ls%dhw-ls%dhw0)/dt)
      x  = 0.5*gsp%ddz(gs%M)*ls%h1/dt
      a(gs%M+1) = xx
      c(gs%M+1) = xx - x
      f(gs%M+1) = - oxyg(gs%M+1,1)*x + Foxyg1 + &
      & a(gs%M+1)*oxyg(gs%M,1) - xx*oxyg(gs%M+1,1)
    case(2)
      c(gs%M+1) = 1.
      a(gs%M+1) = 0.
      f(gs%M+1) = botO2
    end select
    !print*, lamoxyg
    !print*, b
    !print*, c
    !print*, f
    call PROGONKA (a,b,c,f,y,1_iintegers,gs%M+1)
    y(1:gs%M+1) = max(y(1:gs%M+1),0.d0)
    oxyg(1:gs%M+1,2) = y(1:gs%M+1)
  endif

  if (gs%nsoilcols == 1) then
    fbbleflx_o2_sum(0)      = febul0*fracb0(3)/fracb0(1)*fbbleflx_o2(0)     *bathymwater(1)     %area_int
    fbbleflx_o2_sum(1:gs%M) = febul0*fracb0(3)/fracb0(1)*fbbleflx_o2(1:gs%M)*bathymwater(1:gs%M)%area_half
    fbbleflx_o2_sum(gs%M+1) = febul0*fracb0(3)/fracb0(1)*fbbleflx_o2(gs%M+1)*bathymwater(gs%M+1)%area_int
  endif

  if (soilswitch == 1) then
    ! 2-d step of splitting-up scheme - source due to bubble dissolution
    !x = dt*febul0*fracb0(3)/fracb0(1)
    oxyg(1,2) = oxyg(1,2) + 2.d0*dt/(ls%h1*gsp%ddz(1))*(fbbleflx_o2_sum(1) - fbbleflx_o2_sum(0)) / &
    & bathymwater(1)%area_int
    oxyg(1,2) = max(oxyg(1,2),small_value)
    do i = 2, gs%M
      oxyg(i,2) = oxyg(i,2) + dt/(ls%h1*gsp%ddz05(i-1))*(fbbleflx_o2_sum(i) - fbbleflx_o2_sum(i-1)) / &
      & bathymwater(i)%area_int
      oxyg(i,2) = max(oxyg(i,2),small_value)
    enddo
    oxyg(gs%M+1,2) = oxyg(gs%M+1,2) + 2.d0*dt/(ls%h1*gsp%ddz(gs%M)) * &
    & (fbbleflx_o2_sum(gs%M+1) - fbbleflx_o2_sum(gs%M)) / &
    & bathymwater(gs%M+1)%area_int
    oxyg(gs%M+1,2) = max(oxyg(gs%M+1,2),small_value)
  endif

!  x = oxyg(1,2)*ddz(1)/2.
!  x = x + oxyg(M+1,2)*ddz(M)/2.
!  do i = 2, M
!    x = x + oxyg(i,2)*(ddz(i-1)+ddz(i))/2.
!  enddo
  

endif 


deallocate (a, b, c, f, y)

END SUBROUTINE OXYGEN 


!> Calculates the rate of oxygen sources and sinks in the water column,
!! generally following Stefan and Fang, Ecological Modelling, 71 (1994) 37-68
!! V.S., 06.2014
SUBROUTINE OXYGEN_PRODCONS(gs,gsp,wst,bathymsoil,gas,area_int,area_half,ddz05, &
&                          h,dt,Tsoil,Chl_a,dzs,por,oxygsoil, &
&                          itroph,prodox,resp,bod,sod,sodbot)

use DRIVING_PARAMS, only : soilcolconjtype

use PHYS_CONSTANTS, only : &
& Wm22Em2, short2PAR, z0_bot

use TIMEVAR, only : &
& hour_sec, day_sec

use METH_OXYG_CONSTANTS, only : &
& molmass_o2, &
& k_O2_SOD, mubeta0, kc0, &
& thetaC_SOD, thetamu_SOD, &
& tortuosity_coef, &

& T0, &   ! Reference temperature, Celsius
& T00, &  ! Another reference temperature, Celsius
& c1_pmax, c2_pmax, &! Constants for temperature dependence,
                     ! (Stefan and Fang, 1994; Megard et al., 1984)
& k1_c1, k1_c2, &    ! (Megard et al., 1984)
& k2_c1, k2_c2, &
& YCHO2, &   ! mass ratio of chlorophyll-a to oxygen utilized in respiration
& k_r, &     ! day**(-1), (Brown and Barnwell, QUAL2E, 1987; Riley, 1988)
& theta_r, & ! temperature dependence of respiration (Ambrose et al., 1988)
& k_b, &     ! day**(-1), 1-st order decay coefficient, (Riley, 1988; Brown and Barnwell, 1987)
& theta_b1, theta_b2, &   ! temperature-dependence coefficient for biochemical oxygen demand
& mO2C_dec, mCChla_dec, & ! Stoichiometry constants for organics decay in the water column, 
                          ! for BOD, basing on carbonaceous BOD and neglecting 
                          ! nitropgenous BOD, see details in Stefan and Fang, 1994
& Sb20, & ! estimates: 0.5E+3 for oligotrophic, 1.5E+3 for eutrophic lakes, Stefan and Fang, 1994,
          ! mg/(m**2*day), reference sedimentary oxygen demand (SOD)
& theta_s1, theta_s2, &! temperature dependence coefficients for SOD, Stefan and Fang, 1994
& T_md, &              ! the temperature of maximal density

& alpha_POCL, tau_DOC, tau_POCD !constants of Hanson et al. 2004 model

use ARRAYS_GRID, only : &
& gridsize_type, gridspacing_type

use ARRAYS_WATERSTATE, only : waterstate_type
use RADIATION, only : RadWater
use ARRAYS_BATHYM, only : bathym
use ARRAYS, only : gas_type

use PHYS_FUNC, only : DIFF_WATER_OXYGEN, LOGFLUX

use RADIATION, only : nbands

implicit none

! Input variables
!> Grid sizes
type(gridsize_type), intent(in) :: gs
!> Grid spacings
type(gridspacing_type), intent(in) :: gsp
!> Water state variables
type(waterstate_type), intent(in) :: wst
!> Bathymetry of soil columns
type(bathym), intent(in) :: bathymsoil(1:gs%nsoilcols+1,1:gs%nx,1:gs%ny)
!> Carbon species in water and oxygen
type(gas_type), intent(in) :: gas

!> Lake cross-section area at cell interfaces
real(kind=ireals),    intent(in) :: area_int(1:gs%M+1) 
!> Lake cross-section area at cell centers
real(kind=ireals),    intent(in) :: area_half(1:gs%M)  
!> Spacing between cell centers, n/d
real(kind=ireals),    intent(in) :: ddz05(0:gs%M) 
!> Soil temperature, deg. Celsius
real(kind=ireals),    intent(in) :: Tsoil(1:gs%ns,1:gs%nsoilcols)  
!> Chlorophyll-a concentration, mg * l**(-1)
real(kind=ireals),    intent(in) :: Chl_a(1:gs%M+1) 
!> Soil grid spacing, m
real(kind=ireals),    intent(in) :: dzs  (1:gs%ns)  
!> Soil porosity, n/d
real(kind=ireals),    intent(in) :: por  (1:gs%ns)  
!> Mean oxygen concentration in top soil layer
real(kind=ireals),    intent(inout) :: oxygsoil (1:gs%nsoilcols)  
!> Lake depth, m
real(kind=ireals),    intent(in) :: h            
!> Timestep, s
real(kind=ireals),    intent(in) :: dt           
!> Trophic status (1 - oligotrophic, 2 - mesotrophic, 3 - eutrophic)
integer(kind=iintegers), intent(in) :: itroph 

! Output variables
!> Production of O2 by photosynthesis
real(kind=ireals), intent(out) :: prodox(1:gs%M+1) 
!> Sink due to respiration
real(kind=ireals), intent(out) :: resp(1:gs%M+1)   
!> Biochemical oxygen demand
real(kind=ireals), intent(out) :: bod(1:gs%M+1)    
!> Sediment oxygen demand (sink)
real(kind=ireals), intent(out) :: sod(1:gs%M+1)    
!> Sediment oxygen demand at the bottom, positive downwards
real(kind=ireals), intent(out) :: sodbot        

! Local variables

integer(kind=iintegers), parameter :: nSOD = 2 ! The switch for SOD models

real(4) :: Pmax  ! Maximum oxygen production by photosynthesis, h**(-1)
real(4) :: minL  ! Light limiting factor
real(4) :: PAR   ! Photosynthetically active radiation, Einstein /(m*m*h)
real(4) :: BOD_   ! detritus as oxygen equivalent, mg/l


real(4) :: k1_minL, k2_minL, kc, mubeta, oxysurf
real(kind=ireals) :: x, step10, step20, chla_aver, flux, diff
real(kind=ireals), allocatable :: PARz(:)

integer(kind=iintegers) :: i, j, k
integer(kind=iintegers), allocatable :: kcompl(:)

if (nSOD == 3) then
  allocate (kcompl(1:gs%nsoilcols))
  kcompl(:) = 0
endif

! Calculating mean chlorophyll-a density
if     (itroph == 3) then 
! Eutrophic lake
  chla_aver = chla_aver_eu
elseif (itroph == 2) then
! Mesotrophic lake
  chla_aver = chla_aver_meso
elseif (itroph == 1) then
! Oligotrophic lake
  chla_aver = chla_aver_oligo
endif

allocate(PARz(0:gs%M+1))
if (nbands > 1) then
  PARz(0:gs%M+1) = RadWater%flux(0:gs%M+1,2) !2-d band assumed to be PAR
else
  PARz(0:gs%M+1) = RadWater%flux(0:gs%M+1,1)
endif

do i = 1, gs%M
  ! Oxygen production by photosynthesis
  Pmax = c1_pmax*c2_pmax**(wst%Tw2(i) - T0)
  k1_minL = k1_c1*k1_c2**(wst%Tw2(i) - T0)
  step10 = STEP(wst%Tw2(i) - T00) ! step function around 10 Celsius
  k2_minL = step10*k2_c2 + (1 - step10)*k2_c1
  if (i > 1) then
    PAR = 0.5*(PARz(i-1) + PARz(i))*Wm22Em2*hour_sec*short2PAR !converting to Einstein /(m**2*hour)
  else
    PAR = PARz(0)*Wm22Em2*hour_sec*short2PAR !converting to Einstein /(m**2*hour)
  endif
  minL = PAR*(1. + 2.*sqrt(k1_minL/k2_minL)) / (PAR + k1_minL + PAR*PAR/k2_minL)
  prodox(i) = Pmax*minL*Chl_a(i)
  prodox(i) = prodox(i)/(hour_sec*molmass_o2) ! Converting from mg/(l*h) to mol/(m**3*s)
  if     (carbon_model%par == 1) then ! Stefan and Fang model for respiration and BOD
    ! Oxygen sink by respiration
    resp(i) = 1./YCHO2*k_r*theta_r**(wst%Tw2(i) - T0)*Chl_a(i)
    resp(i) = resp(i)/(day_sec *molmass_o2) ! Converting from mg/(l*day) to mol/(m**3*s)
    ! Biochemical oxygen demand (BOD)
    BOD_ = Chla_aver*mO2C_dec*mCChla_dec ! stoichiometry details see in Stefan and Fang, 1994
    step20 = STEP(wst%Tw2(i) - T0) ! step function around 20 Celsius
    x = (theta_b1*step20 + theta_b2*(1-step20))*STEP(wst%Tw2(i) - T_md) ! using step function around 4 Celsius
    if (x == 0.) then
      bod(i) = 0.
    else
      bod(i) = k_b*x**(wst%Tw2(i) - T0)*BOD_
    endif
    bod(i) = bod(i)/(day_sec *molmass_o2) ! --||-- (Converting from mg/(l*day) to mol/(m**3*s))
  elseif (carbon_model%par == 2) then !Hanson et al. model for respiration and BOD
    resp(i) = alpha_POCL*prodox(i)
    !> @todo: introduce different timescale for degradation of allochtonous DOC
    bod(i) = kg_to_g/molmass_ch2o*(gas%POCD(i)/tau_POCD + &
    & gas%DOC(i,1)/tau_DOC + gas%DOC(i,2)/tau_DOC) !Check units!
  endif
  ! Sedimentary oxygen demand (SOD)
  select case (nSOD)
  case(1) ! The model from Stefan and Fang (1994)
    x = step10*theta_s1 + (1 - step10)*theta_s2 !using step function
    sod(i) = Sb20*x**(wst%Tw2(i) - T0)/(day_sec *molmass_o2 * 1.E+3)
  case(2) ! The model from Walker and Sondgrass (1986)
    kc     = kc0*thetaC_SOD**(wst%Tw2(i) - 20.)
    mubeta = mubeta0*thetamu_SOD**(wst%Tw2(i) - 25.)
    sod(i) = kc*gas%oxyg(i,1) + mubeta*gas%oxyg(i,1)/(k_O2_SOD + gas%oxyg(i,1))
  case(3) ! New SOD model 
    k = gsp%ksoilcol(i) !The number of soil column intersecting i-th layer
    if (kcompl(k) == 0) then !checking if SOD calculated for k-th column
      if     (soilcolconjtype == 1) then
        j = bathymsoil(k,gs%ix,gs%iy)%icent
        x = bathymsoil(k,gs%ix,gs%iy)%dzSLc
      elseif (soilcolconjtype == 2) then
        j = bathymsoil(k,gs%ix,gs%iy)%itop
        x = bathymsoil(k,gs%ix,gs%iy)%dzSL
      endif
      ! Oxygen concentration at the soil surface from soil side
      diff = tortuosity_coef*DIFF_WATER_OXYGEN(Tsoil(1,k))
      oxysurf = oxygsoil(k) !OXSURF(Tsoil(1,k),0.5*dzs(1),diff,oxygsoil(k))
      ! Oxygen flux is calculated from logarithmic law
      sod(i) = LOGFLUX(sqrt(wst%u2(j)*wst%u2(j) + wst%v2(j)*wst%v2(j)), &
      & oxysurf/por(1) - gas%oxyg(j,1), x, z0_bot, z0_bot, 1._ireals, 1)
      ! Oxygen concentration is updated to the next time step
      call OXYGEN_TOPSOIL(Tsoil(1,k),sod(i),0.5*dzs(1),dt,oxygsoil(k))
      kcompl(k) = 1
    else
      sod(i) = sod(i-1) ! SOD flux is the same for all water layer intesecting the same soil column
    endif
  end select
  if (i > 1) then
    sod(i) = - sod(i) / area_int(i)*(area_half(i) - area_half(i-1)) / &
    & (h*ddz05(i-1))
  else
    sod(i) = - sod(i) / area_int(i)*(area_half(i) - area_int(i)) / &
    & (h*ddz05(i-1))
  endif
enddo
!sod   (:) = sod   (:)/(day_sec *molmass_o2 * 1.E+3) ! Converting from mg/(m*m*day) to mol/(m*m*s)

! Sedimentary oxygen demand at the bottom (flux)
select case (nSOD)
case(1) ! The model from Stefan and Fang (1994)
  step10 = STEP(wst%Tw2(gs%M+1) - 10.) ! step function around 10 Celsius
  x = step10*theta_s1 + (1 - step10)*theta_s2 ! using step function
  sodbot = Sb20*x**(wst%Tw2(gs%M+1) - T0)/(day_sec *molmass_o2 * 1.E+3) ! Converting from mg/(m*m*day) to mol/(m*m*s)
case(2) ! The model from Walker and Sondgrass (1986)
  kc = kc0*thetaC_SOD**(wst%Tw2(gs%M+1) - 20.)
  mubeta = mubeta0*thetamu_SOD**(wst%Tw2(gs%M+1) - 25.)
  sodbot = kc*gas%oxyg(gs%M+1,1) + mubeta*gas%oxyg(gs%M+1,1)/(k_O2_SOD + gas%oxyg(gs%M+1,1))
case(3) ! The new SOD model
  diff = tortuosity_coef*DIFF_WATER_OXYGEN(Tsoil(1,gs%nsoilcols))
  oxysurf = oxygsoil(gs%nsoilcols) !OXSURF(Tsoil(1,gs%nsoilcols),0.5*dzs(1),diff,oxygsoil(gs%nsoilcols))
  sodbot = LOGFLUX(sqrt(wst%u2(gs%M+1)*wst%u2(gs%M+1) + wst%v2(gs%M+1)*wst%v2(gs%M+1)), &
  & oxysurf/por(1) - gas%oxyg(gs%M+1,1), 0.25*gsp%ddz(gs%M)*h, z0_bot, z0_bot, 1._ireals, 1)
  call OXYGEN_TOPSOIL(Tsoil(1,gs%nsoilcols),sodbot,0.5*dzs(1),dt,oxygsoil(gs%nsoilcols))
end select


deallocate(PARz)

END SUBROUTINE OXYGEN_PRODCONS


SUBROUTINE CHLOROPHYLLA(z_full, H_mixed_layer, H_photic_zone, extwat, M, Chl_a, itroph)
! This subroutine specifies the chlorophyll profile, following
! Stefan and Fang (1994) with modifications:
! - no annual cycle


implicit none

! Input variables
real(kind=ireals), intent(in) :: z_full(1:M+1) ! Depth at cell interfaces, m
real(kind=ireals), intent(in) :: H_mixed_layer ! Mixed layer depth, m
real(kind=ireals), intent(in) :: H_photic_zone ! Photic zone depth, m
real(kind=ireals), intent(in) :: extwat        ! Extinction coefficient, m**(-1)

integer(kind=iintegers), intent(in) :: M ! Number of grid layers (cells)

! Output variables
real(kind=ireals), intent(out) :: Chl_a(1:M+1) ! Chlorophyll-a density, mg/l
integer(kind=iintegers), intent(out) :: itroph ! Trophic status

! Local variables
real(4) :: depthact
real(4) :: chla_aver

integer(kind=iintegers) :: i

if (extwat > extwat1) then 
! Eutrophic lake
  chla_aver = chla_aver_eu
  itroph = 3
elseif (extwat <= extwat1 .and. extwat > extwat2) then
! Mesotrophic lake
  chla_aver = chla_aver_meso
  itroph = 2
else
! Oligotrophic lake
  chla_aver = chla_aver_oligo
  itroph = 1
endif

depthact = max(H_mixed_layer,H_photic_zone)
do i = 1, M+1
  Chl_a(i) = chla_aver*STEP(depthact - z_full(i)) ! step function
enddo

END SUBROUTINE CHLOROPHYLLA


!> Subroutine ADDOXPRODCONS adds to O_2 and CO_2 tendencies associated
!! with organics dynamics in a lake (photosynthesis, respiration,
!! BOD and SOD), and the same for DOC, POCL, POCD
SUBROUTINE ADDOXPRODCONS(M,i_maxN,prodox,resp,bod,sod,dt,ls,gas)

use ARRAYS, only : gas_type
use ARRAYS_BATHYM, only : layers_type

use METH_OXYG_CONSTANTS, only : &
& alpha_POCL, tau_POCD, tau_DOC, &
& beta_POCL, tau_Dh_epi, tau_Dh_hypo, &
& molmass_ch2o

use DRIVING_PARAMS, only : &
& carbon_model

implicit none

! Input variables
!> Production of O2 by photosynthesis, mol/(m**3*s)
real(kind=ireals), intent(in) :: prodox(1:M+1) 
!> Sink due to respiration, mol/(m**3*s)
real(kind=ireals), intent(in) :: resp(1:M+1)   
!> Biochemical oxygen demand, mol/(m**3*s)
real(kind=ireals), intent(in) :: bod(1:M+1)    
!> Sediment oxygen demand (sink), mol/(m**3*s)
real(kind=ireals), intent(in) :: sod(1:M+1)    
!> Time step, s
real(kind=ireals), intent(in) :: dt 
!> Number of model layers (cells)
integer(kind=iintegers), intent(in) :: M 
!> The computational level of maximal Brunt-Vaisala frequency
integer(kind=iintegers), intent(in) :: i_maxN
!> Physical layers' thicknesses
type(layers_type), intent(in) :: ls

!Input/ouput variables
!> Gas concentrations in water, carbon species content
type(gas_type), intent(inout) :: gas

! Local variables

real(kind=ireals) :: x
real(kind=ireals) :: E_POCL, D_DOC1, D_DOC2, P_POCL, R_POCL, Dh_POCL, D_POCD, tau_Dh

integer(kind=iintegers) :: i

do i = 1, M
  x = min(gas%DIC(i,1),dt*CO2O2_prod*prodox(i))
  gas%DIC(i,1)    = gas%DIC(i,1) - x
  gas%oxyg(i,1)   = gas%oxyg(i,1)   + x/CO2O2_prod

  x = min(gas%oxyg(i,1),dt*resp(i))
  gas%oxyg(i,1)   = gas%oxyg(i,1)   - x
  gas%DIC(i,1)    = gas%DIC(i,1) + x*CO2O2_resp

  x = min(gas%oxyg(i,1),dt*bod(i))
  gas%oxyg(i,1)   = gas%oxyg(i,1)   - x
  gas%DIC(i,1)    = gas%DIC(i,1) + x*CO2O2_bod

  x = min(gas%oxyg(i,1),dt*sod(i))
  gas%oxyg(i,1)   = gas%oxyg(i,1)   - x
  gas%DIC(i,1)    = gas%DIC(i,1) + x*CO2O2_sod
enddo

! Hanson et al. model tendencies for carbon species
if (carbon_model%par == 2) then
  do i = 1, M
    ! The death timescale assumed to be the same in hypolimion and during ice period
    if (ls%l1 > 0 .or. i >= i_maxN) then
      tau_Dh = tau_Dh_hypo
    else
      tau_Dh = tau_Dh_epi
    endif
    P_POCL = prodox(i)*molmass_ch2o / kg_to_g
    R_POCL = alpha_POCL*P_POCL
    E_POCL = beta_POCL *P_POCL
    Dh_POCL = gas%POCL(i)/tau_Dh
    D_POCD = gas%POCD(i)/tau_POCD
    D_DOC1 = gas%DOC(i,1)/tau_DOC
    !> @todo: introduce different timescale for degradation of allochtonous DOC
    D_DOC2 = gas%DOC(i,2)/tau_DOC

    gas%DOC(i,1) = gas%DOC(i,1) + dt*(E_POCL - D_DOC1)
    gas%DOC(i,2) = gas%DOC(i,2) + dt*(       - D_DOC2)
    gas%POCL(i) = gas%POCL(i) + dt*(P_POCL - R_POCL - E_POCL - Dh_POCL)
    gas%POCD(i) = gas%POCD(i) + dt*(Dh_POCL - D_POCD)
  enddo
endif


END SUBROUTINE ADDOXPRODCONS


SUBROUTINE OXYGEN_TOPSOIL(Tsoil,influx,h,dt,oxygsoil)

!A bulk model for oxygen in top of soil column,
!calculates oxygen concentration in an aerobic top
!soil layer at the next timestep

use METH_OXYG_CONSTANTS, only : &
k_O2_SOD, mubeta0, thetamu_SOD 

implicit none

!Input variables
real(kind=ireals), intent(in) :: Tsoil  !Temperature of aerobic layer, C
real(kind=ireals), intent(in) :: influx !Oxygen flux from water column, positive downwards
real(kind=ireals), intent(in) :: h      !Thickness of top aerobic soil layer, m
real(kind=ireals), intent(in) :: dt     !Timestep, s

!Input/output variables
real(kind=ireals), intent(inout) :: oxygsoil !Mean oxygen concentration in the top soil layer

!Local variables
real(kind=ireals) :: mubeta

! Michaelis-Menthen (Monod) oxygen depletion according to Walker and Snodgrass, 1986, 
! neglecting "chemical sediment oxygen demand"
mubeta = mubeta0*thetamu_SOD**(Tsoil - 25.)/h
oxygsoil = (oxygsoil + influx/h*dt)/(1. + dt*mubeta/(k_O2_SOD + oxygsoil))

END SUBROUTINE OXYGEN_TOPSOIL


FUNCTION OXSURF(Tsoil,h,diff,oxmean)

!Calculates surface oxygen concentration from
!mean concentration in anaerobic layer. Mean is calculated
!from analytical solution of problem for oxygen concentration
!d^2 C/d^2 z - alpha*C=0 with Dirichlet b.c.s

use METH_OXYG_CONSTANTS, only : &
k_O2_SOD, mubeta0, thetamu_SOD 

implicit none

!Input variables
real(kind=ireals), intent(in) :: Tsoil  !Temperature of aerobic layer, C
real(kind=ireals), intent(in) :: h      !Thickness of anaerobic layer, m
real(kind=ireals), intent(in) :: diff   !Oxygen diffusivity in anaerobic layer
real(kind=ireals), intent(in) :: oxmean !Mean oxygen concentration in anaerobic layer

!Output variables
real(kind=ireals) :: OXSURF

!Local variables
real(kind=ireals) :: alpha, mubeta

mubeta = mubeta0*thetamu_SOD**(Tsoil - 25.)/h
alpha = sqrt(mubeta/(k_O2_SOD*diff))

OXSURF = oxmean*alpha*h*(exp(alpha*h) + 1)/(exp(alpha*h) - 1.)

END FUNCTION OXSURF


END MODULE OXYGEN_MOD

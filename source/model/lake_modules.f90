
MODULE TIMEVAR

use LAKE_DATATYPES, only : ireals, iintegers

real(kind=ireals), parameter :: hour_sec = 60.*60.
real(kind=ireals), parameter :: day_sec = 24.*60.*60.
real(kind=ireals), parameter :: year_sec  = 60.*60.*24.*365.
real(kind=ireals), parameter :: month_sec  = 60.*60.*24.*365.25/12.
real(kind=ireals), parameter :: hour_min = 60.

END MODULE TIMEVAR



MODULE UNITS

use LAKE_DATATYPES, only : ireals

real(kind=ireals), parameter :: kg_to_g = 1.E+3

END MODULE UNITS



MODULE PHYS_PARAMETERS

use LAKE_DATATYPES, only : ireals, iintegers

real(kind=ireals) :: whc, hcond, sncr

SAVE

contains
SUBROUTINE DEFINE_PHYS_PARAMETERS()
implicit none

!cccccccccccccccccccccc SNOW cccccccccccccccccccccccc
whc = 0.04   ! water holding capacity of snow
hcond = 0.01 ! hydraulic conductivity in snow (m/s)
SNCR = 0.01  ! M OF WATER, FROM WHICH SNOW COVERS ALL THE CELL

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
END SUBROUTINE DEFINE_PHYS_PARAMETERS
END MODULE PHYS_PARAMETERS

 
MODULE PHYS_CONSTANTS

use LAKE_DATATYPES, only : ireals, iintegers

real(kind=ireals) :: g
real(kind=ireals) :: omega
real(kind=ireals) :: roa0, row0, rofrsn
real(kind=ireals), target :: roi
real(kind=ireals) :: cw, csnow
real(kind=ireals), target :: ci
real(kind=ireals), target :: Lwi
real(kind=ireals) :: Lwv, Liv
real(kind=ireals) :: kappa
real(kind=ireals) :: sigma
real(kind=ireals) :: cp
real(kind=ireals) :: aMag, bMag
real(kind=ireals) :: Rd, Rwv, R_univ
real(kind=ireals) :: Avogadro_number
real(kind=ireals) :: Planck_const
real(kind=ireals) :: clight
real(kind=ireals) :: lambda_PAR0
real(kind=ireals) :: short2PAR
real(kind=ireals) :: lami, lamw0, lamair
real(kind=ireals) :: lamsal0
real(kind=ireals) :: niu_atm, niu_wat
real(kind=ireals) :: surf_tension_wat
real(kind=ireals) :: roughness0
real(kind=ireals) :: charnock_const
real(kind=ireals) :: albedoofice, albedoofsnow, albedoofwater, albedoofsoil
real(kind=ireals) :: albedoofice_lw, albedoofsnow_lw, albedoofwater_lw, albedoofsoil_lw
real(kind=ireals), target :: emissivityofice, emissivityofwater 
real(kind=ireals), target :: emissivityofsnow, emissivityofsoil
real(kind=ireals) :: aMagw, bMagw
real(kind=ireals) :: aMagi, bMagi
real(kind=ireals) :: alsal, almeth, aloxyg, alcarbdi
real(kind=ireals) :: sabs
real(kind=ireals) :: sabs0
real(kind=ireals) :: Kelvin0
real(kind=ireals) :: pref
real(kind=ireals) :: esatsurf0
real(kind=ireals) :: salice0, poricebot, asalice
real(kind=ireals) :: snmeltwat, icetopfr

real(kind=ireals) :: z0_bot, c_bot

!Derived constants
real(kind=ireals) :: ci_m_roi, cw_m_row0
real(kind=ireals) :: row0_m_Lwi, roi_m_Lwi
real(kind=ireals) :: row0_d_roi, roi_d_row0
real(kind=ireals) :: roa0_d_row0, roa0_d_roi
real(kind=ireals) :: row0_d_roa0
real(kind=ireals) :: Rd_d_cp, Rd_d_Rwv
real(kind=ireals) :: Wm22Em2
real(kind=ireals) :: g_d_Kelvin0


SAVE

contains
SUBROUTINE DEFINE_PHYS_CONSTANTS 
implicit none

g     = 9.814    ! acceleration due to gravity,   m/s**2
omega = 7.29d-5  ! angular velocity of the Earth's rotation,  1/s
sigma = 5.67d-8  ! Stefan-Boltzman constant,W/(m**2*K**4)

roa0  = 1.273    ! reference density of air,kg/m**3
row0  = 1000.    ! reference density of water,    kg/m**3
rofrsn  = 100.   ! reference density of fresh snow,    kg/m**3
roi   = 917.     ! density of ice,    kg/m**3

cw    = 3990.    ! heat capacity of water,  J/(kg*K)
ci    = 2150.    ! heat capacity of ice,    J/(kg*K)
cp    = 1005.    ! heat capacity of air at constant pressure, J/(kg*K)
csnow = 2100.    ! the snow specific heat content,J/(kg K)

lamair = 0.024   ! molecular conductivity of air, J/(m*s*K)
lami  = 2.2      ! molecular conductivity of ice, J/(m*s*K)
lamw0 = 0.561    ! molecular conductivity of water,     J/(m*s*K)
lamsal0 = lamw0/(cw*row0)*1.E-2 !molecular diffusivity coefficient for salinity, m**2/s
niu_atm = 1.8d-5 ! molecular viscosity of air,    m**2/s
niu_wat = 1.307d-6 ! molecular viscosity of water at T = 10 C,m**2/s
surf_tension_wat = 72.d-3 ! surface tension of water at 25 C, N/m

Lwi   = 333500. !267900.  ! latent heat of freezing, J/kg
Lwv   = 2501000. ! latent heat of evaporation,    J/kg
Liv   = 2834600. ! latent heat of sublimation,    J/kg

kappa = 0.38     ! Karman constant,   n/d 
Rd    = 287.05   ! gas constant for dry air,J/(kg*K)
Rwv   = 461.91   ! gas constant for water vapor,  J/(kg*K)
R_univ = 8.31441 ! universal gas constant,  J/(mol*K)
Avogadro_number = 6.022141d+23 !mol**(-1)
Planck_const = 6.62607d-34     !J*s
clight = 3.d+8   ! light speed in vacuum,   m / s
lambda_PAR0 = 550.d-9 ! ~the wavelength of the PAR spectrum center, m
short2PAR = 0.48 ! the global to PAR radiation ratio,   n/d

aMagw = 7.6326   ! coefficient for Magnus formula for water surface
bMagw = 241.9    ! coefficient for Magnus formula for water surface
aMagi = 9.5! coefficient for Magnus formula for ice   surface
bMagi = 265.5    ! coefficient for Magnus formula for ice   surface

albedoofwater = 0.06 ! 0.07 ! albedo of water surface for visible radiation, n/d ! 0.07
albedoofice   = 0.40 ! 0.5! albedo of ice   surface,    n/d ! 0.35
albedoofsnow  = 0.85 ! 0.8! albedo of snow  surface,    n/d ! 0.62
albedoofsoil  = 0.3  ! albedo of soil  surface,   n/d

!     albedo_lw ~ 1 - emissivity
albedoofwater_lw = 0.0 ! albedo of water surface for atmospheric radiation, n/d
albedoofice_lw   = 0.03  ! albedo of ice surface for atmospheric radiation,   n/d
albedoofsnow_lw  = 0.01  ! albedo of snow surface for atmospheric radiation,  n/d
albedoofsoil_lw  = 0.0  ! albedo of soil surface for atmospheric radiation,  n/d

emissivityofwater = 0.98 !emissivity of water surface, n/d ! 0.95
emissivityofice   = 0.97 ! emissivity of ice   surface, n/d ! 0.95
emissivityofsnow  = 0.99 ! emissivity of snow  surface, n/d
emissivityofsoil  = 0.9  ! emissivity of soil  surface, n/d

roughness0 = 1.d-4 ! the reference roughness of water, m
charnock_const = 0.0144 ! Charnock constant, n/d, 0.0144 - classical value
      
Kelvin0 = 273.15
pref    = 1.d+5 ! Reference atmospheric pressure, Pa
esatsurf0 = 610.7 ! Saturation partial water vapor pressure at 0 Celsius, Pa

sabs   = 0. !0.35 !0.4 ! the part of solar radiation, absorbed at the water surface = fraction of NIR 
sabs0  = sabs
alsal  = 1.  ! the ratio of eddy diffusivity for salinity to that of heat
almeth = 1.  ! the ratio of eddy diffusivity for dissolved methane to that of heat
aloxyg = 1.  ! the ratio of eddy diffusivity for dissolved oxygen to that of heat
alcarbdi = 1.  ! the ratio of eddy diffusivity for dissolved carbon dioxode to that of heat

salice0 = 0.! 4.d-3 ! reference salinity of ice
poricebot = 0.05 ! porosity of ice at the bottom of ice layer
asalice = 1./(86400.*365.) !Characteristic time of salt removal from ice

snmeltwat = 1. ! The fraction of snowmelt transferred to water through cracks in ice, 
               ! the rest is partly frozen on a top of the ice layer

icetopfr = 0. ! The switch for freezing of water at the ice top

z0_bot = 1.d-3 ! Roughness parameter of bottom surface, m
c_bot = 1.d-2 !2.5d-3 ! Momentum exchange coefficient at the sloping bottom

ci_m_roi  = ci * roi
cw_m_row0 = cw * row0
row0_m_Lwi = row0 * Lwi
roi_m_Lwi = roi * Lwi

row0_d_roi = row0/roi
roi_d_row0 = roi/row0
roa0_d_row0 = roa0/row0
roa0_d_roi = roa0/roi
row0_d_roa0 = row0/roa0

Rd_d_cp = Rd/cp
Rd_d_Rwv = Rd/Rwv

Wm22Em2 = lambda_PAR0/(Avogadro_number*Planck_const*clight) ! converting multiplyer from W/m2 to Einstein/(m*m*s)
                                                            ! for PAR region

g_d_Kelvin0 = g/Kelvin0

END SUBROUTINE DEFINE_PHYS_CONSTANTS
END MODULE PHYS_CONSTANTS


MODULE NUMERIC_PARAMS

use LAKE_DATATYPES, only : ireals, iintegers

integer(kind=iintegers), parameter :: KL = 1, NT = 52592, &
& Num_Soil = 11, Num_Veget = 13, ML = 131, MS = 100 
integer(kind=iintegers), target, save :: ML_ = ML, MS_ = MS
real(kind=ireals), parameter :: dznorm = 0.04, dzmin = 0.015, &
& UpperLayer = 0.02 ! UpperLayer = 0.1, UpperLayer = 0.008 

integer(kind=iintegers), parameter :: vector_length = 150
real(kind=ireals),    parameter :: pi = 4.*datan(1.d0)

real(kind=ireals), parameter :: small_value = 1.d-20

real(kind=ireals), parameter :: min_ice_thick = 0.02 
real(kind=ireals), parameter :: min_water_thick = 0.02
real(kind=ireals), parameter :: T_phase_threshold = 1.d-5

 !   ML --- total number of layers
 !   MS --- maximal number of layers in snow
 !   dznorm --- standard layer depth in snow (m)
 !   dzmin --- min layer depth in snow (m)
 !   depth --- depth of the lowest layer in soil (m)
 !   UpperLayer --- upper soil layer depth (m)

END MODULE NUMERIC_PARAMS


MODULE ATMOS

use LAKE_DATATYPES, only : ireals, iintegers

real(kind=ireals) :: tempair
real(kind=ireals) :: relhum
real(kind=ireals) :: humair
real(kind=ireals) :: pressure
real(kind=ireals) :: hw
real(kind=ireals) :: xlew
real(kind=ireals) :: cdmw
real(kind=ireals) :: hw_input, xlew_input, cdmw_input, surfrad_input
real(kind=ireals) :: uwind, vwind
real(kind=ireals) :: wind, wind10, windwork
real(kind=ireals) :: longwave
real(kind=ireals) :: Radbal
real(kind=ireals) :: hflux
real(kind=ireals) :: surfrad
real(kind=ireals) :: hbal
real(kind=ireals) :: zref
real(kind=ireals) :: botflux
real(kind=ireals) :: cdmw2
real(kind=ireals) :: shortwave, sabspen
real(kind=ireals) :: precip
real(kind=ireals) :: Sflux0
real(kind=ireals) :: Elatent
real(kind=ireals) :: Radbal_surf
real(kind=ireals) :: eflux, eflux0, eflux0_kinem
real(kind=ireals) :: hskin
real(kind=ireals) :: tau, tau_wav
real(kind=ireals) :: velfrict, velfrict_prev
real(kind=ireals) :: cdm1
real(kind=ireals) :: turb_density_flux
real(kind=ireals) :: cloud

! Former common-block wind
real(kind=ireals) :: urel, vrel, u, v
     
SAVE 

END MODULE ATMOS


!> Radiation data structures and operations on them
MODULE RADIATION

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

type, public :: rad_type
  sequence
  integer(kind=iintegers) :: nbands, nlevels
  real(kind=ireals), allocatable       :: extinct(:,:), flux(:,:), integr(:)
  procedure(RAD_UPDATE), pointer, pass :: RAD_UPDATE => RAD_UPDATE
  procedure(RAD_INIT)  , pointer, pass :: RAD_INIT   => RAD_INIT
endtype rad_type
type(rad_type) :: RadWater, RadIce, RadDeepIce


integer(kind=iintegers), parameter :: nbands = 4 !wavebands: UV, PAR, NIR, IR
real(kind=ireals), parameter :: extwatbands(1:nbands) = (/3.13,0.81,1.73,1.087E+3/) !Extinction coefficients for bands
real(kind=ireals), parameter :: fracbands(1:nbands) = (/0.046,0.429,0.208,0.317/) !The energy fraction of bands in incident shortwave radiation

contains

!>Subroutine integrates Beer-Lambert equation in all spectral bands
SUBROUTINE RAD_UPDATE(this,thickness,flux0)
implicit none
type(rad_type), intent(inout) :: this
!Input variables
real(kind=ireals), intent(in) :: thickness(0:this%nlevels-1)
real(kind=ireals), intent(in) :: flux0(this%nbands)
!Local variables
integer(kind=iintegers) :: i, j

this%flux(0,1:this%nbands) = flux0(1:this%nbands)
this%integr(0) = sum(this%flux(0,1:this%nbands))
do j = 1, this%nbands
  do i = 1, this%nlevels
    this%flux(i,j) = this%flux(i-1,j) * exp(-this%extinct(i,j)*thickness(i-1))
  enddo
enddo
forall(i=1:this%nlevels) this%integr(i) = sum(this%flux(i,1:this%nbands))
END SUBROUTINE RAD_UPDATE


!> Initializes the radiation variable structure
SUBROUTINE RAD_INIT(this,nbands,nlevels)
type(rad_type), intent(inout) :: this
integer(kind=iintegers), intent(in) :: nbands, nlevels
this%nbands = nbands
this%nlevels = nlevels
allocate(this%flux   (0:nlevels,this%nbands))
allocate(this%extinct(1:nlevels,this%nbands))
allocate(this%integr (0:nlevels            ))
END SUBROUTINE RAD_INIT


END MODULE RADIATION


MODULE ARRAYS_WATERSTATE

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

 ! Water state variables

 real(kind=ireals), allocatable, target :: Tw1(:),Tw2(:),Ti1(:),Ti2(:),Tis1(:), &
 & Tis2(:),RS(:),lamw(:),lamw_back(:),lamsal(:),Sal1(:),Sal2(:), &
 & preswat(:), salice(:), porice(:), ci_m_roi_v(:), lami_v(:)
 !real(kind=ireals), pointer :: SR(:),SRi(:),SRdi(:)
 type, public :: waterstate_type
   sequence
   real(kind=ireals), pointer :: Tw1(:),Tw2(:),Ti1(:),Ti2(:),Tis1(:), &
   & Tis2(:),RS(:),SR(:),lamw(:),lamsal(:),Sal1(:),Sal2(:),row(:), &
   & SRi(:),SRdi(:),preswat(:),u1(:),v1(:),u2(:),v2(:),extwatarr(:), &
   & salice(:), porice(:), ci_m_roi_v(:), lami_v(:), lamw_back(:)
 endtype waterstate_type
 type(waterstate_type) :: wst

 type, public :: intbal_type
   sequence
   real(kind=ireals), allocatable :: terms(:), termsdt(:)
   real(kind=ireals) :: storage, storage0, dstorage, resid
   integer(kind=iintegers) :: nterms
   !contains
   procedure(intbal_update), pointer, pass :: intbal_update => INTBAL_UPDATE
 endtype intbal_type
 
 contains
 

 SUBROUTINE INTBAL_UPDATE(this,terms,storage1,storage2,dt)
 implicit none
 !Input/output variables
 type(intbal_type) :: this
 real(kind=ireals), intent(in), dimension(:) :: terms
 real(kind=ireals), intent(in) :: storage1, storage2, dt
 !Local variables
 logical, save :: firstcall = .true.
 integer(kind=iintegers) :: i

 if (firstcall) then
   this%storage0 = storage1
 endif

 !do i = 1, this%nterms
 this%terms(:) = terms(:)
 this%termsdt(:) = this%termsdt(:) + terms(:)*dt
 this%storage = storage2
 this%dstorage = storage2 - this%storage0
 this%resid = this%dstorage - sum(this%termsdt(:))
 !enddo

 if (firstcall) firstcall = .false.
 END SUBROUTINE INTBAL_UPDATE


END MODULE ARRAYS_WATERSTATE



MODULE ARRAYS_GRID

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

! Grid size group 
integer(kind=iintegers), save, target :: nsoilcols ! Number of soil columns per one lake
type, public :: gridsize_type
  sequence
  integer(kind=iintegers), pointer :: M, Mice, ns, ms, ml, nsoilcols, nx, ny
  integer(kind=iintegers) :: ix, iy, isoilcol
end type
type(gridsize_type) :: gs

! Grid spacing group
real(kind=ireals), allocatable, target :: dz_full(:), z_full(:), z_half(:), zsoilcols(:,:,:)
real(kind=ireals), allocatable, target :: ddz(:), ddz2(:), ddz05(:), ddz052(:), ddz054(:)
real(kind=ireals), allocatable, target :: ddzi(:), ddzi05(:)
real(kind=ireals), allocatable, target :: dzeta_int(:), dzeta_05int(:)
real(kind=ireals), allocatable, target :: dzetai_int(:), dzetai_05int(:)      
integer(kind=iintegers), allocatable, target :: isoilcolc(:), ksoilcol(:)
type, public :: gridspacing_type
  sequence
  real(kind=ireals), pointer :: dz_full(:), z_full(:), z_half(:), zsoilcols(:,:,:)
  real(kind=ireals), pointer :: ddz(:), ddz2(:), ddz05(:), ddz052(:), ddz054(:)
  real(kind=ireals), pointer :: ddzi(:), ddzi05(:)
  real(kind=ireals), pointer :: dzeta_int(:), dzeta_05int(:)
  real(kind=ireals), pointer :: dzetai_int(:), dzetai_05int(:) 
  integer(kind=iintegers), pointer :: isoilcolc(:), ksoilcol(:)
end type
type(gridspacing_type) :: gsp

END MODULE ARRAYS_GRID



MODULE ARRAYS_METHANE

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

!Methane group
 real(kind=ireals) :: veg
 real(kind=ireals) :: fplant, fdiff_lake_surf 
 real(kind=ireals), allocatable :: fdiffbot(:)
 real(kind=ireals), allocatable :: febul0(:)
! real(kind=ireals) :: febul2(1:2) ! two-meth
! real(kind=ireals) :: fdiff2(1:2), fdiff_lake_surf2(1:2) ! two-meth
 real(kind=ireals), target :: febultot(1:2), febulbottot(1:2)
 real(kind=ireals) :: fdifftot(1:2), fdiff_lake_surftot(1:2)
 real(kind=ireals) :: methoxidwat(1:2), ice_meth_oxid
 real(kind=ireals) :: qwateroxidtot, qwateroxidML
 real(kind=ireals) :: ice_meth_oxid_total
 real(kind=ireals) :: meth_cont_wat_total0
! real(kind=ireals) :: febultot2(1:2,1:2), fdifftot2(1:2,1:2), fdiff_lake_surftot2(1:2,1:2) ! two-meth
 real(kind=ireals) :: plant_sum, bull_sum, oxid_sum, rprod_sum
 real(kind=ireals), allocatable :: rprod_total_oldC(:), rprod_total_newC(:)
 real(kind=ireals) :: rprod_total_newC_integr(1:2), rprod_total_oldC_integr(1:2)
! real(kind=ireals) :: rprod_total_newC_integr2(1:2,1:2), rprod_total_oldC_integr2(1:2,1:2) ! two-meth
 real(kind=ireals) :: h_talik
 real(kind=ireals) :: anox
 real(kind=ireals), allocatable :: rootss(:), TgrAnn(:)
 real(kind=ireals), allocatable, target :: qsoil(:,:), qwater(:,:), qmethane(:)
 real(kind=ireals), allocatable :: fbbleflx_ch4(:,:), fbbleflx_ch4_sum(:)
 real(kind=ireals), allocatable :: fracb0(:)
 real(kind=ireals), allocatable :: methgenmh(:) 
! real(kind=ireals), allocatable :: qsoil2(:,:), qwater2(:,:,:)  ! two-meth
 real(kind=ireals), allocatable :: lammeth(:)
 real(kind=ireals) :: tot_ice_meth_bubbles
! real(kind=ireals) :: tot_ice_meth_bubbles2(1:2)  ! two-meth

END MODULE ARRAYS_METHANE



MODULE ARRAYS_OXYGEN

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

!Oxygen group      
real(kind=ireals), allocatable, target :: oxyg(:,:)
real(kind=ireals), allocatable :: lamoxyg(:)
real(kind=ireals), allocatable :: fbbleflx_o2(:,:), fbbleflx_o2_sum(:)
real(kind=ireals), allocatable :: Chl_a(:) ! Chlorophyll-a concentration, mg/l
real(kind=ireals), allocatable :: prodox(:) ! O_2 production due to photosynthesis
real(kind=ireals), allocatable :: resp(:) ! O_2 sink by respiration
real(kind=ireals), allocatable :: bod(:) ! Biochemical oxygen demand 
real(kind=ireals), allocatable :: sod(:) ! Sedimentary oxygen demand (at the lake margins)
real(kind=ireals) :: sodbot ! Sedimentary oxygen demand at the lake bottom
real(kind=ireals), target :: H_photic_zone ! The photic zone depth, m
real(kind=ireals) :: fo2 ! Diffusive flux of oxygen to the atmosphere, mol/(m**2*s)
integer(kind=iintegers) :: itroph ! Trophic status of a lake

real(kind=ireals), allocatable :: oxygsoil(:) ! Mean oxygen concentration in the uppermost (aerobic)
                                           ! soil layer

!Carbon dioxide group 
real(kind=ireals), allocatable, target :: DIC(:,:)
real(kind=ireals), allocatable, target :: DOC(:,:) !DOC(:,1) stands for autochtonous DOC, DOC(:,2) -- for allochtonous one
real(kind=ireals), allocatable, target :: POCL(:), POCD(:)
real(kind=ireals), allocatable :: lamcarbdi(:)
real(kind=ireals) :: fco2 ! Diffusive flux of oxygen to the atmosphere, mol/(m**2*s)
real(kind=ireals), allocatable :: fbbleflx_co2(:,:), fbbleflx_co2_sum(:)

END MODULE ARRAYS_OXYGEN



MODULE ARRAYS_TURB

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

real(kind=ireals), allocatable :: S(:),Gen(:),F(:),TKE_turb_trans(:),Gen_seiches(:)
real(kind=ireals), allocatable :: KT(:),k2(:),E1(:),eps1(:)
real(kind=ireals), allocatable :: KC(:),KLengT(:),RSR(:) 
real(kind=ireals), allocatable :: E2(:),eps2(:), &
& k1(:),k3(:),k4(:),k5(:),u3(:), &
& v3(:),C1aup(:),Re(:),Ri(:),E_it1(:), &
& E_it2(:),C_num(:),E_it3(:),Feps(:),E_it21(:), &
& eps_it21(:),eps_it1(:),eps_it2(:), &
& l(:),k2_mid(:),k3_mid(:),k4_mid(:),k5_mid(:), &
& E12(:),eps12(:),knum(:),k2t(:), &
& Eeps(:),Ecent(:),Epscent(:),res_E(:), &
& res_eps(:),dresd(:,:,:,:),dres2dE(:),dres2deps(:), &
& veg_sw(:) 
real(kind=ireals), allocatable, target :: row(:), row2(:), rowc(:)
real(kind=ireals), allocatable :: WU_(:),WV_(:),GAMT(:),GAMU(:),GAMV(:), TF(:),KLengu(:)
real(kind=ireals), allocatable :: PEMF    (:) , PDENS_DOWN (:), &
&                       PT_DOWN (:) , PSAL_DOWN  (:), &
&                       pt_down_f(:)
real(kind=ireals), allocatable :: k_turb_T_flux(:), T_massflux(:)

integer(kind=iintegers) :: i_maxN
integer(kind=iintegers), allocatable :: itherm(:)

real(kind=ireals) :: Buoyancy0,H_mixed_layer,maxN,w_conv_scale,u_star,T_conv_scale,H_entrainment, &
& signwaveheight = 0.d0, Ri_bulk, Wedderburn, Rossby_rad, LakeNumber, ThermThick, ReTherm, RiTherm
real(kind=ireals) :: S_integr_positive, S_integr_negative, Gen_integr, &
& eps_integr, E_integr, TKE_turb_trans_integr 
real(kind=ireals) :: Seps_integr_positive, Seps_integr_negative, &
& Geneps_integr, epseps_integr
real(kind=ireals) :: TKE_balance, eps_balance
real(kind=ireals) :: Eseiches

! Data structure for turbulent characteristics
type, public :: turb_type
  sequence
  real(kind=ireals), allocatable :: Rp(:), Rpdens(:)
  real(kind=ireals), pointer :: Ri(:)
  real(kind=ireals), pointer :: Buoyancy0, H_mixed_layer, maxN, w_conv_scale, &
  & u_star, T_conv_scale, H_entrainment, signwaveheight, Ri_bulk, Wedderburn, &
  & Rossby_rad, LakeNumber, ThermThick, ReTherm, RiTherm
end type
type(turb_type) :: trb


END MODULE ARRAYS_TURB



MODULE ARRAYS_BATHYM

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

!Bathymetry group
type, public :: bathym
  sequence
  real(kind=ireals) :: area_int, area_half, Lx, Ly, Lx_half, Ly_half, dzSL, dzSLc, rad_int
  real(kind=ireals), allocatable :: sarea_int(:), sarea_half(:)
  integer(kind=iintegers) :: itop, ibot, icent
end type
type(bathym), allocatable :: bathymwater(:), bathymice(:) 
type(bathym), allocatable :: bathymdice(:), bathymsoil(:,:,:)
real(kind=ireals), allocatable :: area_int(:), area_half(:), Lx(:), Ly(:)
real(kind=ireals), allocatable :: aki(:,:), bki(:,:)

! Layers' thicknesses group
real(kind=ireals), target :: h1, l1, hs1, ls1, vol, botar
real(kind=ireals), target :: hx1, hx2, hy1, hy2, hx1t, hx2t, hy1t, hy2t
real(kind=ireals), allocatable, target :: hx1ml(:), hx2ml(:), hy1ml(:), hy2ml(:)
type, public :: layers_type
  sequence
  real(kind=ireals), pointer :: h1, hx1, hx2, hy1, hy2, l1, hs1, ls1, vol, &
  & dhw, dhw0, dhi, dhi0, dls, dls0, dhiimp, snmeltice, botar, H_photic_zone, &
  & hx1t, hx2t, hy1t, hy2t
  real(kind=ireals), pointer :: dhwhigh, dhwlow, dhwp, dhwe, dhwtrib, dhwls, dhwsnmelt
  real(kind=ireals), pointer :: dhihigh, dhilow, dhip, dhis, dhif
  real(kind=ireals), pointer :: hx1ml(:), hx2ml(:), hy1ml(:), hy2ml(:)
end type
type(layers_type) :: ls

real(kind=ireals) :: voldef = 0.d0
real(kind=ireals), target :: dhw,dhw0,dhi,dhi0,dls,dls0,dhiimp,snmeltice
real(kind=ireals), target :: dhwhigh, dhwlow, dhwp, dhwe, dhwtrib, dhwls, dhwsnmelt
real(kind=ireals), target :: dhihigh, dhilow, dhip, dhis, dhif

real(kind=ireals) :: dhwfsoil = 0.

integer(kind=iintegers), parameter :: form_ellipse = 1, form_rectangle = 2

END MODULE ARRAYS_BATHYM



MODULE ARRAYS_SOIL

use LAKE_DATATYPES, only : ireals, iintegers
!use DRIVING_PARAMS

implicit none

 ! Soil group
 real(kind=ireals), allocatable:: rosoil(:), csoil(:), &
 & lamsoil(:), dzs(:), dzss(:), zsoil(:), wsoil(:)
 real(kind=ireals), allocatable :: lammoist(:), rosdry(:), filtr(:), wa(:)
 real(kind=ireals), allocatable :: Tsoil1(:,:), Tsoil2(:), Tsoil3(:,:)
 real(kind=ireals), allocatable :: wl1(:,:), wl2(:), wl3(:), wl4(:,:)
 real(kind=ireals), allocatable :: wi1(:,:), wi2(:,:)
 real(kind=ireals), allocatable :: Sals1(:,:), Sals2(:,:)
 real(kind=ireals), allocatable :: soilflux(:,:,:)
 real(kind=ireals), allocatable :: methsoilflux(:), co2soilflux(:), o2soilflux(:)
 type, public :: lsh_type
   sequence
   real(kind=ireals), allocatable :: water(:), ice(:), dice(:)
 end type lsh_type
 type(lsh_type) :: lsh, lsm, lsc, lso, lsu, lsv

 real(kind=ireals), allocatable :: WLM0(:),WLM7(:),bH(:),PSIMAX(:),POR(:),FLWMAX(:),DLMAX(:)

END MODULE ARRAYS_SOIL



 MODULE ARRAYS
 
!MODULE ARRAYS contains allocatable arrays of the model  

 use LAKE_DATATYPES, only : ireals, iintegers
! use DRIVING_PARAMS
 
 integer(kind=iintegers), parameter :: soil_indic = 1
 integer(kind=iintegers), parameter :: ice_indic = 2
 integer(kind=iintegers), parameter :: ice_sal_indic = 12
 integer(kind=iintegers), parameter :: snow_indic = 3
 integer(kind=iintegers), parameter :: water_indic = 4
 integer(kind=iintegers), parameter :: deepice_indic = 5
 integer(kind=iintegers), parameter :: water_salinity_indic = 6
 integer(kind=iintegers), parameter :: soil_salinity_indic = 7
 integer(kind=iintegers), parameter :: water_methane_indic = 8
 integer(kind=iintegers), parameter :: water_oxygen_indic = 9
 integer(kind=iintegers), parameter :: water_carbdi_indic = 10
 integer(kind=iintegers), parameter :: water_carbon_indic = 11


 integer(kind=iintegers) :: snow,ice,water,deepice,nstep = 0
 real(kind=ireals) :: time, Erad, dep_av
 real(kind=ireals) :: zinv
 integer(kind=iintegers) :: time_int,time_int_old,man,par,flag_snow_init

!Computational time group
 integer(kind=iintegers), parameter :: ncomp = 10
 real(kind=ireals) :: comptime(1:ncomp) = 0.d0
 
!Parameter calibration group      
 real(kind=ireals), pointer, save :: calibr_par1, calibr_par2, cost_function(:)
 logical, save :: calibrate_parameter = .false.
  
 type, public :: gas_type
   sequence
   real(kind=ireals), pointer :: qmethane(:),qsoil(:,:),qwater(:,:), &
   & oxyg(:,:),DIC(:,:),DOC(:,:),POCL(:),POCD(:)
 end type gas_type
 type(gas_type), target :: gas

!Surface characteristics group
 real(kind=ireals) :: roughness,emissivity,albedo,albedo_lw,aM,bM,relhums
 
 real(kind=ireals) totalevaps,totalmelts,totalprecips
 real(kind=ireals), allocatable :: Tskin(:)
 real(kind=ireals) :: T_0dim

 real(kind=ireals) :: SR_botsnow
 real(kind=ireals) :: dt_keps
 real(kind=ireals) :: dt05, dt_inv, dt_inv05
 real(kind=ireals) :: deltaskin


 real(kind=ireals), allocatable :: z_watersoil(:)
 real(kind=ireals), allocatable, target :: u1(:),v1(:),u2(:),v2(:),um(:),vm(:),uv(:)

 real(kind=ireals), allocatable :: dep_2d(:,:)
 integer(kind=iintegers), allocatable :: init(:,:), num(:)

 real(kind=ireals), allocatable :: utend(:), vtend(:), Twtend(:), Saltend(:)
 real(kind=ireals), allocatable :: ueffl(:), veffl(:), Tweffl(:), Saleffl(:)

 real(kind=ireals), pointer :: workc(:)

 !SAVE

 END MODULE ARRAYS

      

 MODULE TURB_CONST

 use LAKE_DATATYPES, only : ireals, iintegers
 

!PHYSICAL CONSTANTS
  real(kind=ireals), parameter:: kar       = 0.4d0
  !real(kind=ireals), parameter:: niu       = 1.007d-6

!DIMENSIONLESS CONSTANTS IN E-EPS PARAMETERIZATION, according to Goudsmit et al., 2002
  real(kind=ireals), parameter:: CE0       = 0.09d0 
  real(kind=ireals), parameter:: CEt0      = 0.072d0
  real(kind=ireals), parameter:: ceps1     = 1.44d0
  real(kind=ireals), parameter:: ceps2     = 1.92d0
  real(kind=ireals), parameter:: ceps3_unstable = 1.14d0
  real(kind=ireals), parameter:: ceps3_stable   = -0.4d0
!        ceps3 = 1.14d0
!         ceps3 < -0.4 should be used for stable stratification (Burchard, 2002)
!         ceps3 = 1. ! following (Kochergin and Sklyar, 1992)
!         ceps3 = 1.14d0 ! Baum and Capony (1992)
!         Other values for ceps3:
!         0.<ceps3<0.29 for stable and ceps3 = 1.44 for unstable conditions (Rodi, 1987);
!         ceps3 >= 1 (Kochergin and Sklyar, 1992)
  real(kind=ireals), parameter:: sigmaE    = 1.d0
!  real(kind=ireals), parameter:: sigmaeps  = 1.3d0
  real(kind=ireals), save :: sigmaeps
  
  real(kind=ireals), save :: CL 
  real(kind=ireals), save :: alpham
                                          
  real(kind=ireals), parameter :: C_wstar = 0.4d0   ! Deardorff, 1970

  real(kind=ireals), parameter:: lam_T     = 1.d0 
  real(kind=ireals), parameter:: lam_gen   = 1.d0 
  real(kind=ireals), parameter:: lamTM     = 0.14d0

! real(kind=ireals), parameter:: cmu_0     = 0.094d0
! Other values for cmu_0, proposed in literature, are
! 0.121, 0.077 - see Burchard, 2002

  real(kind=ireals), parameter:: lam_E0     = CE0/sigmaE   
  real(kind=ireals), save :: lam_eps0

! The coefficients in formulas for stability functions (Canuto et al., 2001)
  real(kind=ireals), parameter:: CEcoef01   = 0.127d0
  real(kind=ireals), parameter:: CEcoef02   = 0.1526d-1
  real(kind=ireals), parameter:: CEcoef03   = 0.16d-3

  real(kind=ireals), parameter:: CEtcoef01  = 0.1190d0
  real(kind=ireals), parameter:: CEtcoef02  = 0.4294d-2
  real(kind=ireals), parameter:: CEtcoef03  = 0.66d-3

  real(kind=ireals), parameter:: CEcoef11   = 1.d0
  real(kind=ireals), parameter:: CEcoef12   = 0.2d0
  real(kind=ireals), parameter:: CEcoef13   = 0.315d-1
  real(kind=ireals), parameter:: CEcoef14   = 0.58d-2
  real(kind=ireals), parameter:: CEcoef15   = 0.4d-2
  real(kind=ireals), parameter:: CEcoef16   = 0.4d-4
 
! The parameters of Mellor and Yamada model (1982)      

  real(kind=ireals), parameter:: A1  =  0.92d0
  real(kind=ireals), parameter:: B1  =  1.66d+1
  real(kind=ireals), parameter:: C1  =  0.8d-1
  real(kind=ireals), parameter:: A2  =  0.74d0
  real(kind=ireals), parameter:: B2  =  1.01d+1
  real(kind=ireals), save :: CL_K_KL_MODEL

!Co is according to Satyanarayana et al., 1999
  real(kind=ireals), parameter:: Co        = 1.9d0

!CONSTANTS IN CONVECTION PARAMETERIZATION
  real(kind=ireals), parameter:: KC0       = 0.d0 !1.d+5 !500.
  real(kind=ireals), parameter:: dtdz0     = 2.d0

!CONSTANTS OF SIMOES'S PARAMETERIZATION
  real(kind=ireals), parameter:: Cs        = 1.5d0  !0.15
  real(kind=ireals), parameter:: Cs1       = 100.d0    
 
!lo=0.4*ddz*h1*3.                                                         VS,06.2007
  real(kind=ireals), parameter:: L0        = 0.1d0 !0.1d0

  real(kind=ireals), parameter:: CON0      = 0.05d0
  real(kind=ireals), parameter:: CON1      = 0.11d0
  real(kind=ireals), parameter:: CON2      = 1.56d0
  real(kind=ireals), parameter:: CONUV     = 1.d0

  real(kind=ireals), parameter :: kar_tilde = 2.d-1 ! Burchard (2002), p.80

  real(kind=ireals), parameter :: grav_wave_Gill_const = 0.7 !Burchard(2002), Gill(1982)

  !real(kind=ireals), parameter :: min_turb = 1.d-10 !Minimal turbulent viscosity/diffusivity coefficient, m**2/s
  real(kind=ireals), parameter :: min_diff = 1.d-10 !2.d-7 !Minimal turbulent diffusivity coefficient, m**2/s
  real(kind=ireals), parameter :: Pr_stable = 1. !1.E+1 !Turbulent Prandtl number at very stable stratification
  real(kind=ireals), parameter :: min_visc = min_diff*Pr_stable

  contains

  SUBROUTINE TURB_CONST_DERIVED_SET

  use DRIVING_PARAMS, only : kepsbc

  implicit none

  CL = CE0**(0.75d0) ! This value results from assumption
                     ! that shear production is balanced by dissipation
                     ! (Burchard, 2002)

  CL_K_KL_MODEL = 2.d0**(1.5d0)/B1

  alpham = (1.5*CE0**0.5*sigmaE*kar_tilde**(-2))**0.5

  select case(kepsbc%par)
    case(1)
      sigmaeps  = kar*kar/(CE0**0.5 * (ceps2 - ceps1) ) ! according to Burchard, 2002,
                                                        ! fulfills the law of the wall
                                                        ! in k-eps model
    case(2)
      sigmaeps = (4./3.*alpham + 1.)*(alpham + 1.) * &
      & kar_tilde*kar_tilde/(ceps2*CE0**0.5)  ! according to Burchard, 2002,
                                              ! fulfills the analytical solution for 
                                              ! unstratified shear-free flow with 
                                              ! wave breaking TKE input
    case default
      sigmaeps  = 1.3d0
  end select
!  sigmaeps = 1.3d0

  lam_eps0   = CE0/sigmaeps 


 END SUBROUTINE TURB_CONST_DERIVED_SET

 END MODULE TURB_CONST


!> Module EVOLUTION_VARIABLES contains variables
!! that have to be saved for the next time step (~prognostic variables)
MODULE EVOLUTION_VARIABLES

use LAKE_DATATYPES, only : ireals, iintegers

real(kind=ireals), allocatable, public :: l1_2d(:,:)
real(kind=ireals), allocatable, public :: h1_2d(:,:)
real(kind=ireals), allocatable, public :: hx1_2d(:,:), hx2_2d(:,:)
real(kind=ireals), allocatable, public :: hy1_2d(:,:), hy2_2d(:,:)
real(kind=ireals), allocatable, public :: hx1t_2d(:,:), hx2t_2d(:,:)
real(kind=ireals), allocatable, public :: hy1t_2d(:,:), hy2t_2d(:,:)
real(kind=ireals), allocatable, public :: hx1ml_2d(:,:,:), hx2ml_2d(:,:,:)
real(kind=ireals), allocatable, public :: hy1ml_2d(:,:,:), hy2ml_2d(:,:,:)
real(kind=ireals), allocatable, public :: ls1_2d(:,:)
real(kind=ireals), allocatable, public :: hs1_2d(:,:)
real(kind=ireals), allocatable, public :: time_2d(:,:) 
real(kind=ireals), allocatable, public :: cdm_2d(:,:)
real(kind=ireals), allocatable, public :: u_2d(:,:,:)
real(kind=ireals), allocatable, public :: v_2d(:,:,:)
real(kind=ireals), allocatable, public :: Tsoil1_2d(:,:,:,:)
real(kind=ireals), allocatable, public :: wi1_2d(:,:,:,:)
real(kind=ireals), allocatable, public :: wl1_2d(:,:,:,:)
real(kind=ireals), allocatable, public :: Tw1_2d(:,:,:)
real(kind=ireals), allocatable, public :: Ti1_2d(:,:,:)
real(kind=ireals), allocatable, public :: Tis1_2d(:,:,:)
real(kind=ireals), allocatable, public :: dz_2d(:,:,:)
real(kind=ireals), allocatable, public :: T_2d(:,:,:)
real(kind=ireals), allocatable, public :: wl_2d(:,:,:)
real(kind=ireals), allocatable, public :: dens_2d(:,:,:)
real(kind=ireals), allocatable, public :: E_2d(:,:,:)
real(kind=ireals), allocatable, public :: eps_2d(:,:,:)
real(kind=ireals), allocatable, public :: dhwfsoil_2d(:,:) 
real(kind=ireals), allocatable, public :: Sal1_2d(:,:,:)
real(kind=ireals), allocatable, public :: Sals1_2d(:,:,:,:)
real(kind=ireals), allocatable, public :: ueffl_2d(:,:,:)
real(kind=ireals), allocatable, public :: veffl_2d(:,:,:)
real(kind=ireals), allocatable, public :: Tweffl_2d(:,:,:)
real(kind=ireals), allocatable, public :: Saleffl_2d(:,:,:)
real(kind=ireals), allocatable, public :: Elatent_2d(:,:)
real(kind=ireals), allocatable, public :: dhw_2d(:,:)
real(kind=ireals), allocatable, public :: dhw0_2d(:,:)
real(kind=ireals), allocatable, public :: dhi_2d(:,:)
real(kind=ireals), allocatable, public :: dhi0_2d(:,:)
real(kind=ireals), allocatable, public :: dls0_2d(:,:)
real(kind=ireals), allocatable, public :: velfrict_2d(:,:)
real(kind=ireals), allocatable, public :: roughness_2d(:,:)
real(kind=ireals), allocatable, public :: eflux0_kinem_2d(:,:)
real(kind=ireals), allocatable, public :: lamw_2d(:,:,:)
real(kind=ireals), allocatable, public :: snmelt_2d(:,:)
real(kind=ireals), allocatable, public :: snowmass_2d(:,:)
real(kind=ireals), allocatable, public :: Tskin_2d(:,:)
real(kind=ireals), allocatable, public :: k_turb_T_flux_2d(:,:,:)
real(kind=ireals), allocatable, public :: qwater_2d(:,:,:)
real(kind=ireals), allocatable, public :: qsoil_2d(:,:,:,:)
real(kind=ireals), allocatable, public :: oxyg_2d(:,:,:)
real(kind=ireals), allocatable, public :: oxygsoil_2d(:,:,:)
real(kind=ireals), allocatable, public :: carbdi_2d(:,:,:)
real(kind=ireals), allocatable, public :: DOC_2d(:,:,:,:)
real(kind=ireals), allocatable, public :: POCL_2d(:,:,:)
real(kind=ireals), allocatable, public :: POCD_2d(:,:,:)
real(kind=ireals), allocatable, public :: tot_ice_meth_bubbles_2d(:,:)
real(kind=ireals), allocatable, public :: febul0_2d(:,:,:)
real(kind=ireals), allocatable, public :: Eseiches_2d(:,:)
real(kind=ireals), allocatable, public :: salice_2d(:,:,:)
real(kind=ireals), allocatable, public :: porice_2d(:,:,:)

integer(kind=iintegers), allocatable, public :: fl_sn_2d(:,:)
integer(kind=iintegers), allocatable, public :: fl_sn_init_2d(:,:) 
integer(kind=iintegers), allocatable, public :: itop_2d(:,:)
integer(kind=iintegers), allocatable, public :: nstep_2d(:,:)
integer(kind=iintegers), allocatable, public :: i_maxN_2d(:,:)
integer(kind=iintegers), allocatable, public :: itherm_2d(:,:,:)

logical, public :: arrays_allocated = .false.

SAVE

contains
!> Subroutine ALLOCATE_ARRAYS allocates arrays of EVOLUTION_VARIABLES module
SUBROUTINE ALLOCATE_ARRAYS(nx, ny, M, Mice, ns, ms, ml, nsoilcols)

use DRIVING_PARAMS, only : &
& dyn_pgrad, &
& carbon_model

implicit none

! Input variables
 integer(kind=iintegers), intent(in) :: nx, ny
 integer(kind=iintegers), intent(in) :: M, Mice, ns, ms, ml, nsoilcols

 allocate ( l1_2d(nx,ny) ) ; l1_2d = 0
 allocate ( h1_2d(nx,ny) ) ; h1_2d = 0
 allocate ( hx1_2d(nx,ny), hx2_2d(nx,ny) ) ; hx1_2d = 0 ; hx2_2d = 0
 allocate ( hy1_2d(nx,ny), hy2_2d(nx,ny) ) ; hy1_2d = 0 ; hy2_2d = 0
 allocate ( hx1t_2d(nx,ny), hx2t_2d(nx,ny) ) ; hx1t_2d = 0 ; hx2t_2d = 0
 allocate ( hy1t_2d(nx,ny), hy2t_2d(nx,ny) ) ; hy1t_2d = 0 ; hy2t_2d = 0
 if (dyn_pgrad%par == 3 .or. dyn_pgrad%par == 4) then ! Multilayer dynamic pressure gradient model
   allocate ( hx1ml_2d(1:M+1,nx,ny), hx2ml_2d(1:M+1,nx,ny) ) ; hx1ml_2d = 0. ; hx2ml_2d = 0.
   allocate ( hy1ml_2d(1:M+1,nx,ny), hy2ml_2d(1:M+1,nx,ny) ) ; hy1ml_2d = 0. ; hy2ml_2d = 0.
 endif
 allocate ( ls1_2d(nx,ny) ) ; ls1_2d = 0
 allocate ( cdm_2d(nx,ny) ) ; cdm_2d = 0
 allocate ( u_2d(M+1,nx,ny) ) ; u_2d = 0
 allocate ( v_2d(M+1,nx,ny) ) ; v_2d = 0
 allocate ( E_2d(M+1,nx,ny) ) ; E_2d = 0
 allocate ( eps_2d(M+1,nx,ny) ) ; eps_2d = 0
 allocate ( Tsoil1_2d(ns,nsoilcols,nx,ny) ) ; Tsoil1_2d = 0
 allocate ( wi1_2d(ns,nsoilcols,nx,ny) ) ; wi1_2d = 0
 allocate ( wl1_2d(ns,nsoilcols,nx,ny) ) ; wl1_2d = 0
 allocate ( Tw1_2d(M+1,nx,ny) ) ; Tw1_2d = 0
 allocate ( Sal1_2d(M+1,nx,ny) ) ; Sal1_2d = 0
 allocate ( Ti1_2d(Mice+1,nx,ny) ) ; Ti1_2d = 0
 allocate ( Tis1_2d(Mice+1,nx,ny) ) ; Tis1_2d = 0
 allocate ( k_turb_T_flux_2d(M+1,nx,ny) ) ; k_turb_T_flux_2d = 0
 allocate ( Sals1_2d(ns,nsoilcols,nx,ny) ) ; Sals1_2d = 0
 allocate ( ueffl_2d(M+1,nx,ny) ) ; ueffl_2d = 0
 allocate ( veffl_2d(M+1,nx,ny) ) ; veffl_2d = 0
 allocate ( Tweffl_2d(M+1,nx,ny) ) ; Tweffl_2d = 0
 allocate ( Saleffl_2d(M+1,nx,ny) ) ; Saleffl_2d = 0
 allocate ( hs1_2d(nx,ny) ) ; hs1_2d = 0
 allocate ( dz_2d(ms,nx,ny) ) ; dz_2d = 0
 allocate ( T_2d(ml,nx,ny) ) ; T_2d = 0
 allocate ( wl_2d(ml,nx,ny) ) ; wl_2d = 0
 allocate ( dens_2d(ms,nx,ny) ) ; dens_2d = 0
 allocate ( time_2d(nx,ny) ) ; time_2d = 0
 allocate ( dhwfsoil_2d(nx,ny) ) ; dhwfsoil_2d = 0
 allocate ( Elatent_2d(nx,ny) ) ; Elatent_2d = 0
 allocate ( dhw_2d(nx,ny) ) ; dhw_2d = 0
 allocate ( dhw0_2d(nx,ny) ) ; dhw0_2d = 0
 allocate ( dhi_2d(nx,ny) ) ; dhi_2d = 0
 allocate ( dhi0_2d(nx,ny) ) ; dhi0_2d = 0
 allocate ( dls0_2d(nx,ny) ) ; dls0_2d = 0
 allocate ( velfrict_2d(nx,ny) ) ; velfrict_2d = 0
 allocate ( roughness_2d(nx,ny) ) ; roughness_2d = 0
 allocate ( eflux0_kinem_2d(nx,ny) ) ; eflux0_kinem_2d = 0
 allocate ( lamw_2d(1:M,nx,ny) ) ; lamw_2d = 0
 allocate ( snmelt_2d(nx,ny) ) ; snmelt_2d = 0
 allocate ( snowmass_2d(nx,ny) ) ; snowmass_2d = 0
 allocate ( fl_sn_2d(nx,ny) ) ; fl_sn_2d = 0
 allocate ( fl_sn_init_2d(nx,ny) ) ; fl_sn_init_2d = 0
 allocate ( itop_2d(nx,ny) ) ; itop_2d = 0
 allocate ( nstep_2d(nx,ny) ) ; nstep_2d = 0
 allocate ( Tskin_2d(nx,ny) ) ; Tskin_2d = 0
 allocate ( qwater_2d(M+1,nx,ny) ) ; qwater_2d = 0.
 allocate ( qsoil_2d(ns,nsoilcols,nx,ny) ) ; qsoil_2d = 0.
 allocate ( oxyg_2d(M+1,nx,ny) ) ; oxyg_2d = 0.
 allocate ( oxygsoil_2d(nsoilcols,nx,ny) ) ; oxygsoil_2d = 0.
 allocate ( carbdi_2d(M+1,nx,ny) ) ; carbdi_2d = 0. 
 if (carbon_model%par == 2) then ! Model for DIC/DOC/POCL/POCD from Hanson et al., 2004
   allocate ( DOC_2d(M+1,nx,ny,2) ) ; DOC_2d = 0. 
   allocate ( POCL_2d(M+1,nx,ny) ) ; POCL_2d = 0. 
   allocate ( POCD_2d(M+1,nx,ny) ) ; POCD_2d = 0. 
 endif
 allocate ( tot_ice_meth_bubbles_2d(nx,ny) ) ; tot_ice_meth_bubbles_2d = 0.
 allocate ( febul0_2d(nsoilcols,nx,ny) ) ; febul0_2d = 0.
 allocate ( Eseiches_2d(nx,ny) ) ; Eseiches_2d = 0.
 allocate ( salice_2d(Mice+1,nx,ny) ) ; salice_2d = 0.
 allocate ( porice_2d(Mice+1,nx,ny) ) ; porice_2d = 0.
 allocate ( i_maxN_2d(nx,ny) ) ; i_maxN_2d = 0
 if (dyn_pgrad%par == 4) then
   allocate ( itherm_2d(M+2,nx,ny) ) 
   itherm_2d = 0
 endif

 arrays_allocated = .true.

END SUBROUTINE ALLOCATE_ARRAYS


!> Subroutine UPDATE_CURRENT_TIMESTEP gets the values of prognostic variables saved at the 
!! previous timestep
SUBROUTINE UPDATE_CURRENT_TIMESTEP ( &
& ix, iy, nx, ny, M, Mice, ns, ms, ml, nsoilcols, &
& l1, h1, hx1, hx2, hy1, hy2, ls1, hs1, &
& hx1t, hx2t, hy1t, hy2t, &
& hx1ml, hx2ml, hy1ml, hy2ml, &
& u1, v1, &
& E1, eps1, k_turb_T_flux, &
& Tsoil1, Sals1, wi1, wl1, &
& Tw1, Sal1, lamw, &
& Tskin, &
& Ti1, Tis1, &
& ueffl, veffl, Tweffl, Saleffl, &
& dz, T, wl, dens, &
& qwater, qsoil, &
& oxyg, oxygsoil, &
& DIC,DOC,POCL,POCD, &
& snmelt, snowmass, &
& cdmw, &
& time, &
& dhwfsoil, &
& Elatent, &
& dhw, dhw0, &
& dhi, dhi0, dls0, &
& velfrict, &
& roughness, &
& eflux0_kinem, &
& tot_ice_meth_bubbles, &
& febul0, Eseiches, salice, porice, &

& flag_snow, flag_snow_init, &
& itop, &
& nstep, i_maxN, itherm )

implicit none

! Input variables
integer(kind=iintegers), intent(in) :: ix, iy, nx, ny
integer(kind=iintegers), intent(in) :: M, Mice, ns, ms, ml, nsoilcols

! Output variables
real(kind=ireals), intent(out) :: l1, h1, hx1, hx2, hy1, hy2, ls1, hs1
real(kind=ireals), intent(out) :: hx1t, hx2t, hy1t, hy2t
real(kind=ireals), intent(inout), allocatable :: hx1ml(:), hx2ml(:), hy1ml(:), hy2ml(:)
real(kind=ireals), intent(out) :: u1(1:M+1), v1(1:M+1)
real(kind=ireals), intent(out) :: E1(1:M+1), eps1(1:M+1), k_turb_T_flux(1:M)
real(kind=ireals), intent(out) :: Tsoil1(1:ns,1:nsoilcols), Sals1(1:ns,1:nsoilcols)
real(kind=ireals), intent(out) :: wi1(1:ns,1:nsoilcols), wl1(1:ns,1:nsoilcols)
real(kind=ireals), intent(out) :: Tw1(1:M+1), Sal1(1:M+1), lamw(1:M)
real(kind=ireals), intent(out) :: Tskin
real(kind=ireals), intent(out) :: Ti1(1:Mice+1), Tis1(1:Mice+1)
real(kind=ireals), intent(out) :: ueffl(1:M+1), veffl(1:M+1), Tweffl(1:M+1), Saleffl(1:M+1)
real(kind=ireals), intent(out) :: dz(1:ms), T(1:ml), wl(1:ml), dens(1:ms)
real(kind=ireals), intent(out) :: qwater(1:M+1), qsoil(1:ns,1:nsoilcols)
real(kind=ireals), intent(out) :: oxyg(1:M+1), DIC(1:M+1)
real(kind=ireals), intent(inout), allocatable :: DOC(:,:), POCL(:), POCD(:)
real(kind=ireals), intent(out) :: oxygsoil(1:nsoilcols)
real(kind=ireals), intent(out) :: snmelt, snowmass
real(kind=ireals), intent(out) :: cdmw
real(kind=ireals), intent(out) :: time
real(kind=ireals), intent(out) :: dhwfsoil
real(kind=ireals), intent(out) :: Elatent
real(kind=ireals), intent(out) :: dhw, dhw0
real(kind=ireals), intent(out) :: dhi, dhi0, dls0
real(kind=ireals), intent(out) :: velfrict
real(kind=ireals), intent(out) :: roughness
real(kind=ireals), intent(out) :: eflux0_kinem
real(kind=ireals), intent(out) :: tot_ice_meth_bubbles
real(kind=ireals), intent(out) :: febul0(1:nsoilcols)
real(kind=ireals), intent(out) :: Eseiches
real(kind=ireals), intent(out) :: salice(1:Mice+1)
real(kind=ireals), intent(out) :: porice(1:Mice+1)

integer(kind=iintegers), intent(out) :: flag_snow, flag_snow_init
integer(kind=iintegers), intent(out) :: itop
integer(kind=iintegers), intent(out) :: nstep, i_maxN 
integer(kind=iintegers), intent(inout), allocatable :: itherm(:)

! Local variables
integer(kind=iintegers) :: i

if (.not.arrays_allocated) &
& call ALLOCATE_ARRAYS(nx, ny, M, Mice, ns, ms, ml, nsoilcols)

l1 = l1_2d(ix,iy) 
h1 = h1_2d(ix,iy)
hx1 = hx1_2d(ix,iy)
hx2 = hx2_2d(ix,iy)
hy1 = hy1_2d(ix,iy)
hy2 = hy2_2d(ix,iy)
ls1 = ls1_2d(ix,iy) 
hs1 = hs1_2d(ix,iy)

! Thermocline displacement
hx1t = hx1t_2d(ix,iy)
hx2t = hx2t_2d(ix,iy)
hy1t = hy1t_2d(ix,iy)
hy2t = hy2t_2d(ix,iy)

! Displacements of each layer
if (allocated(hx1ml)) then 
  hx1ml(1:M+1) = hx1ml_2d(1:M+1,ix,iy)
  hx2ml(1:M+1) = hx2ml_2d(1:M+1,ix,iy)
  hy1ml(1:M+1) = hy1ml_2d(1:M+1,ix,iy)
  hy2ml(1:M+1) = hy2ml_2d(1:M+1,ix,iy)
endif

do i = 1, M+1
  u1(i) = u_2d(i,ix,iy)
  v1(i) = v_2d(i,ix,iy)
enddo

do i = 1, M
  E1(i) = E_2d(i,ix,iy)
  eps1(i) = eps_2d(i,ix,iy)
  k_turb_T_flux(i) = k_turb_T_flux_2d(i,ix,iy)
enddo
    
do i = 1, ns 
  Tsoil1(i,1:nsoilcols) = Tsoil1_2d(i,1:nsoilcols,ix,iy)
  Sals1 (i,1:nsoilcols) = Sals1_2d (i,1:nsoilcols,ix,iy)
  wi1   (i,1:nsoilcols) = wi1_2d   (i,1:nsoilcols,ix,iy)
  wl1   (i,1:nsoilcols) = wl1_2d   (i,1:nsoilcols,ix,iy)
  qsoil (i,1:nsoilcols) = qsoil_2d (i,1:nsoilcols,ix,iy)
enddo
    
do i = 1, M+1 
  Tw1(i) = Tw1_2d(i,ix,iy)
  Sal1(i) = Sal1_2d(i,ix,iy)
  qwater(i) = qwater_2d(i,ix,iy)
  oxyg(i) = oxyg_2d(i,ix,iy)
  DIC(i) = carbdi_2d(i,ix,iy)
enddo
if (allocated(DOC)) then
  DOC (1:M+1,1:2) = DOC_2d (1:M+1,ix,iy,1:2)
  POCL(1:M+1)     = POCL_2d(1:M+1,ix,iy)
  POCD(1:M+1)     = POCD_2d(1:M+1,ix,iy)
endif

oxygsoil(1:nsoilcols) = oxygsoil_2d(1:nsoilcols,ix,iy)
  
Tskin = Tskin_2d(ix,iy)
 
do i = 1, Mice + 1 
  Ti1(i) = Ti1_2d(i,ix,iy)
  Tis1(i) = Tis1_2d(i,ix,iy)
enddo

ueffl(1:M+1) = ueffl_2d(1:M+1,ix,iy)
veffl(1:M+1) = veffl_2d(1:M+1,ix,iy)
Tweffl(1:M+1) = Tweffl_2d(1:M+1,ix,iy)
Saleffl(1:M+1) = Saleffl_2d(1:M+1,ix,iy)

flag_snow = fl_sn_2d(ix,iy)
flag_snow_init = fl_sn_init_2d(ix,iy)

itop = itop_2d(ix,iy)

do i = max(1,itop), ms
  dz(i) = dz_2d(i,ix,iy)
  T(i) = T_2d(i,ix,iy)
  wl(i) = wl_2d(i,ix,iy)
  dens(i) = dens_2d(i,ix,iy)
enddo

snmelt = snmelt_2d(ix,iy)
snowmass = snowmass_2d(ix,iy)
 
cdmw = cdm_2d(ix,iy)
  
time = time_2d(ix,iy)
nstep = nstep_2d(ix,iy)
i_maxN = i_maxN_2d(ix,iy)
if (allocated(itherm)) itherm(1:M+2) = itherm_2d(1:M+2,ix,iy)

dhwfsoil = dhwfsoil_2d(ix,iy)
Elatent = Elatent_2d(ix,iy)
dhw = dhw_2d(ix,iy)
dhw0 = dhw0_2d(ix,iy)
dhi = dhi_2d(ix,iy)
dhi0 = dhi0_2d(ix,iy)


dls0 = dls0_2d(ix,iy)
velfrict = velfrict_2d(ix,iy)
roughness = roughness_2d(ix,iy)
eflux0_kinem = eflux0_kinem_2d(ix,iy)

tot_ice_meth_bubbles = tot_ice_meth_bubbles_2d(ix,iy)
febul0(1:nsoilcols) = febul0_2d(1:nsoilcols,ix,iy)

Eseiches = Eseiches_2d(ix,iy)
salice(1:Mice+1) = salice_2d(1:Mice+1,ix,iy)
porice(1:Mice+1) = porice_2d(1:Mice+1,ix,iy)

lamw(1:M) = lamw_2d(1:M,ix,iy)
  
END SUBROUTINE UPDATE_CURRENT_TIMESTEP


!> Subroutine UPDATE_NEXT_TIMESTEP saves the values of prognostic variables for the 
!! next timestep
SUBROUTINE UPDATE_NEXT_TIMESTEP ( &
& ix, iy, nx, ny, M, Mice, ns, ms, ml, nsoilcols, &
& l1, h1, hx1, hx2, hy1, hy2, ls1, hs1, &
& hx1t, hx2t, hy1t, hy2t, &
& hx1ml, hx2ml, hy1ml, hy2ml, &
& u1, v1, &
& E1, eps1, k_turb_T_flux, &
& Tsoil1, Sals1, wi1, wl1, &
& Tw1, Sal1, lamw, &
& Tskin, &
& Ti1, Tis1, &
& ueffl, veffl, Tweffl, Saleffl, &
& dz, T, wl, dens, &
& qwater, qsoil, &
& oxyg, oxygsoil, &
& DIC, DOC, POCL, POCD, &
& snmelt, snowmass, &
& cdmw, &
& time, &
& dhwfsoil, &
& Elatent, &
& dhw,dhw0, &
& dhi, dhi0, dls0, &
& velfrict, &
& roughness, &
& eflux0_kinem, &
& tot_ice_meth_bubbles, &
& febul0, Eseiches, salice, porice, &

& flag_snow, flag_snow_init, &
& itop, &
& nstep, i_maxN, itherm)

implicit none

! Input variables
integer(kind=iintegers), intent(in) :: ix, iy, nx, ny
integer(kind=iintegers), intent(in) :: M, Mice, ns, ms, ml, nsoilcols

real(kind=ireals), intent(in) :: l1, h1, hx1, hx2, hy1, hy2, ls1, hs1
real(kind=ireals), intent(in) :: hx1t, hx2t, hy1t, hy2t
real(kind=ireals), intent(in), allocatable :: hx1ml(:), hx2ml(:), hy1ml(:), hy2ml(:)
real(kind=ireals), intent(in) :: u1(1:M+1), v1(1:M+1)
real(kind=ireals), intent(in) :: E1(1:M+1), eps1(1:M+1), k_turb_T_flux(1:M)
real(kind=ireals), intent(in) :: Tsoil1(1:ns,1:nsoilcols), Sals1(1:ns,1:nsoilcols) 
real(kind=ireals), intent(in) :: wi1(1:ns,1:nsoilcols), wl1(1:ns,1:nsoilcols)
real(kind=ireals), intent(in) :: Tw1(1:M+1), Sal1(1:M+1), lamw(1:M)
real(kind=ireals), intent(in) :: Tskin
real(kind=ireals), intent(in) :: Ti1(1:Mice+1), Tis1(1:Mice+1)
real(kind=ireals), intent(in) :: ueffl(1:M+1), veffl(1:M+1), Tweffl(1:M+1), Saleffl(1:M+1)
real(kind=ireals), intent(in) :: dz(1:ms), T(1:ml), wl(1:ml), dens(1:ms)
real(kind=ireals), intent(in) :: qwater(1:M+1), qsoil(1:ns,1:nsoilcols)
real(kind=ireals), intent(in) :: oxyg(1:M+1), DIC(1:M+1)
real(kind=ireals), intent(inout), allocatable :: DOC(:,:), POCL(:), POCD(:)
real(kind=ireals), intent(in) :: oxygsoil(1:nsoilcols)
real(kind=ireals), intent(in) :: snmelt, snowmass
real(kind=ireals), intent(in) :: cdmw
real(kind=ireals), intent(in) :: time
real(kind=ireals), intent(in) :: dhwfsoil
real(kind=ireals), intent(in) :: Elatent
real(kind=ireals), intent(in) :: dhw, dhw0
real(kind=ireals), intent(in) :: dhi, dhi0, dls0
real(kind=ireals), intent(in) :: velfrict
real(kind=ireals), intent(in) :: roughness
real(kind=ireals), intent(in) :: eflux0_kinem
real(kind=ireals), intent(in) :: tot_ice_meth_bubbles
real(kind=ireals), intent(in) :: febul0(1:nsoilcols)
real(kind=ireals), intent(in) :: Eseiches
real(kind=ireals), intent(in) :: salice(1:Mice+1)
real(kind=ireals), intent(in) :: porice(1:Mice+1)

integer(kind=iintegers), intent(in) :: flag_snow, flag_snow_init
integer(kind=iintegers), intent(in) :: itop
integer(kind=iintegers), intent(in) :: nstep, i_maxN 
integer(kind=iintegers), intent(inout), allocatable :: itherm(:)

! Local variables
integer(kind=iintegers) :: i

if (.not.arrays_allocated) &
& call ALLOCATE_ARRAYS(nx, ny, M, Mice, ns, ms, ml, nsoilcols)

l1_2d(ix,iy)  = l1 
h1_2d(ix,iy)  = h1
hx1_2d(ix,iy) = hx1
hx2_2d(ix,iy) = hx2
hy1_2d(ix,iy) = hy1
hy2_2d(ix,iy) = hy2
ls1_2d(ix,iy) = ls1
hs1_2d(ix,iy) = hs1

! Thermocline displacement
hx1t_2d(ix,iy) = hx1t
hx2t_2d(ix,iy) = hx2t
hy1t_2d(ix,iy) = hy1t
hy2t_2d(ix,iy) = hy2t

! Displacements of each layer
if (allocated(hx1ml)) then
  hx1ml_2d(1:M+1,ix,iy) = hx1ml(1:M+1)
  hx2ml_2d(1:M+1,ix,iy) = hx2ml(1:M+1)
  hy1ml_2d(1:M+1,ix,iy) = hy1ml(1:M+1)
  hy2ml_2d(1:M+1,ix,iy) = hy2ml(1:M+1)
endif

do i = 1, M+1
  u_2d(i,ix,iy) = u1(i)
  v_2d(i,ix,iy) = v1(i)
enddo

do i = 1, M
  E_2d(i,ix,iy) = E1(i)
  eps_2d(i,ix,iy) = eps1(i)
  k_turb_T_flux_2d(i,ix,iy) = k_turb_T_flux(i)
enddo
    
do i = 1, ns 
  Tsoil1_2d(i,1:nsoilcols,ix,iy) = Tsoil1(i,1:nsoilcols)
  Sals1_2d (i,1:nsoilcols,ix,iy) = Sals1 (i,1:nsoilcols) 
  wi1_2d   (i,1:nsoilcols,ix,iy) = wi1   (i,1:nsoilcols)
  wl1_2d   (i,1:nsoilcols,ix,iy) = wl1   (i,1:nsoilcols)
  qsoil_2d (i,1:nsoilcols,ix,iy) = qsoil (i,1:nsoilcols)
enddo
    
do i = 1, M+1 
  Tw1_2d(i,ix,iy) = Tw1(i)
  Sal1_2d(i,ix,iy) = Sal1(i) 
  qwater_2d(i,ix,iy) = qwater(i)
  oxyg_2d(i,ix,iy) = oxyg(i)
  carbdi_2d(i,ix,iy) = DIC(i)
enddo
if (allocated(DOC)) then
  DOC_2d (1:M+1,ix,iy,1:2) = DOC (1:M+1,1:2)
  POCL_2d(1:M+1,ix,iy)     = POCL(1:M+1)
  POCD_2d(1:M+1,ix,iy)     = POCD(1:M+1)
endif 
oxygsoil_2d(1:nsoilcols,ix,iy) = oxygsoil(1:nsoilcols)

Tskin_2d(ix,iy) = Tskin

do i = 1, Mice+1
  Ti1_2d(i,ix,iy) = Ti1(i)
  Tis1_2d(i,ix,iy) = Tis1(i)
enddo

ueffl_2d(1:M+1,ix,iy) = ueffl(1:M+1)
veffl_2d(1:M+1,ix,iy) = veffl(1:M+1)
Tweffl_2d(1:M+1,ix,iy) = Tweffl(1:M+1)
Saleffl_2d(1:M+1,ix,iy) = Saleffl(1:M+1)
 
fl_sn_2d(ix,iy) = flag_snow
fl_sn_init_2d(ix,iy) = flag_snow_init

itop_2d(ix,iy) = itop

do i = max(1,itop), ms
  dz_2d(i,ix,iy) = dz(i)
  T_2d(i,ix,iy) = T(i)
  wl_2d(i,ix,iy) = wl(i)
  dens_2d(i,ix,iy) = dens(i)
enddo

snmelt_2d(ix,iy) = snmelt
snowmass_2d(ix,iy) = snowmass

cdm_2d(ix,iy) = cdmw

time_2d(ix,iy) = time
nstep_2d(ix,iy) = nstep
i_maxN_2d(ix,iy) = i_maxN
if (allocated(itherm)) itherm_2d(1:M+2,ix,iy) = itherm(1:M+2)


dhwfsoil_2d(ix,iy) = dhwfsoil
Elatent_2d(ix,iy) = Elatent
dhw_2d(ix,iy) = dhw
dhw0_2d(ix,iy) = dhw0
dhi_2d(ix,iy) = dhi
dhi0_2d(ix,iy) = dhi0
dls0_2d(ix,iy) = dls0
velfrict_2d(ix,iy) = velfrict
roughness_2d(ix,iy) = roughness
eflux0_kinem_2d(ix,iy) = eflux0_kinem

tot_ice_meth_bubbles_2d(ix,iy) = tot_ice_meth_bubbles
febul0_2d(1:nsoilcols,ix,iy) = febul0(1:nsoilcols)

Eseiches_2d(ix,iy) = Eseiches
salice_2d(1:Mice+1,ix,iy) = salice(1:Mice+1)
porice_2d(1:Mice+1,ix,iy) = porice(1:Mice+1)

lamw_2d(1:M,ix,iy) = lamw(1:M)

END SUBROUTINE UPDATE_NEXT_TIMESTEP

END MODULE EVOLUTION_VARIABLES



MODULE METH_OXYG_CONSTANTS

use LAKE_DATATYPES, only : ireals, iintegers
use TIMEVAR, only : day_sec
!use NUMERIC_PARAMS, only : small_value

real(kind=ireals), parameter :: pH = 6.0 !The pH of the water in a lake, assumed to be constant

real(kind=ireals), parameter :: row0_ = 1000. !The same as row0 in PHYS_CONSTANTS

real(kind=ireals), parameter :: molmass_ch4 = 16. ! methane molar mass, g/mol
real(kind=ireals), parameter :: molmass_co2 = 44. ! co2 molar mass, g/mol
real(kind=ireals), parameter :: molmass_o2  = 32. ! o2 molar mass, g/mol        
real(kind=ireals), parameter :: molmass_n2  = 28. ! n2 molar mass, g/mol
real(kind=ireals), parameter :: molmass_ar  = 40. ! argon molar mass, g/mol
real(kind=ireals), parameter :: molmass_h2o = 18. ! water molar mass, g/mol
real(kind=ireals), parameter :: molmass_so4 = 96. ! sulfate molar mass, g/mol
real(kind=ireals), parameter :: molmass_ch2o = 30. ! sugar molar mass, g/mol

real(kind=ireals), parameter :: molm3tonM = 1.d+6 ! multiplier converting mol/m**3 to nanomol/l
real(kind=ireals), parameter :: molm3tomcM = 1.d+3 ! multiplier converting mol/m**3 to micromol/l
real(kind=ireals), parameter :: molm3toppm_co2 = molmass_co2*1.d+3/row0_
real(kind=ireals), parameter :: molm3toppm_ch4 = molmass_ch4*1.d+3/row0_
real(kind=ireals), parameter :: molm3toppm_o2  = molmass_o2 *1.d+3/row0_
real(kind=ireals), parameter :: molm3tomkgl_ch4 = molmass_ch4*1.d+3
real(kind=ireals), parameter :: molm3tomgl_co2 = molmass_co2
real(kind=ireals), parameter :: molm3tomgl_o2 = molmass_o2
real(kind=ireals), parameter :: molm3tomgl_ch4 = molmass_ch4

real(kind=ireals), parameter :: tortuosity_coef = 0.66 ! the tortuosity of soil pores

real(kind=ireals), parameter :: diff_unsat = 2.d-7 ! Diffusion constant in unsaturated soil, m**2/s
real(kind=ireals), parameter :: diff_water = 2.d-9 ! Diffusion constant in water, m**2/s

real(kind=ireals), parameter :: ch4_atm0 = 2.28d-6 ! concentration of methane in water, in equilibrium with
                                         ! atmospheric partial methane pressure (0.175 Pa), 
                                         ! according to Henry law with reference constant, mol/m**3
                                         ! 7.6d-5 - atmospheric concentration of methane, mol/m**3
real(kind=ireals), parameter :: ch4_pres_atm0 = 1.75d-1 ! atmospheric partial methane pressure, Pa
real(kind=ireals), parameter :: o2_atm0 = 2.73d-1  ! concentration of oxygen in water, in equilibrium with
                                         ! atmospheric partial oxygen pressure (21278 Pa), 
                                         ! according to Henry law with reference constant, mol/m**3
                                         ! 8. - atmospheric concentration of oxygen, mol/m**3
real(kind=ireals), parameter :: o2_pres_atm0 = 21278. ! atmospheric partial oxygen pressure, Pa
real(kind=ireals), parameter :: n2_atm0 = 4.96d-1  ! concentration of nitrogen in water, in equilibrium with
                                         ! atmospheric partial nitrogen pressure (80046 Pa), 
                                         ! according to Henry law with reference constant, mol/m**3
real(kind=ireals), parameter :: n2_pres_atm0 = 80046. ! atmospheric partial nitrogen pressure, Pa
real(kind=ireals), parameter :: co2_pres_atm0 = 39.01 ! atmospheric partial carbon dioxide pressure, Pa
real(kind=ireals), parameter :: co2_atm0 = 1.3d-2  ! concentration of carbon dioxide in water, in equilibrium with
                                         ! atmospheric partial carbon dioxide pressure (39.01 Pa), 
                                         ! according to Henry law with reference constant, mol/m**3
                                        
real(kind=ireals), parameter :: Henry_const0_ch4 = 1.3d-5 ! Henry constant for methane at the reference temperature, mol/(m**3*Pa)
real(kind=ireals), parameter :: Henry_const0_o2 = 1.3d-5 ! Henry constant for oxygen at the reference temperature, mol/(m**3*Pa)
real(kind=ireals), parameter :: Henry_const0_n2 = 6.2d-6 ! Henry constant for nitrogen at the reference temperature, mol/(m**3*Pa)
real(kind=ireals), parameter :: Henry_const0_co2 = 3.36d-4 ! Henry constant for carbon doixide at the reference temperature, mol/(m**3*Pa)
real(kind=ireals), parameter :: Henry_const0_ar = 1.4d-5 ! Henry constant for argon at the reference temperature, mol/(m**3*Pa)
real(kind=ireals), parameter :: Henry_temp_ref = 298.15 ! The reference temperature for Henry constants, K
real(kind=ireals), parameter :: Henry_temp_dep_ch4 = 1.7d+3 ! The temperature dependence of Henry constant for methane, K
real(kind=ireals), parameter :: Henry_temp_dep_o2 = 1.7d+3 ! The temperature dependence of Henry constant for oxygen, K
real(kind=ireals), parameter :: Henry_temp_dep_n2 = 1.3d+3 ! The temperature dependence of Henry constant for nitrogen, K
real(kind=ireals), parameter :: Henry_temp_dep_co2 = 2.4d+3 ! The temperature dependence of Henry constant for carbon dioxide, K
real(kind=ireals), parameter :: Henry_temp_dep_ar = 1.1d+3 ! The temperature dependence of Henry constant for argon, K


real(kind=ireals), parameter :: dens_ch4 = 4.2262d+2 ! liquid methane density at boiling point, kg/m**3
real(kind=ireals), parameter :: dens_co2 = 1.032d+3  ! liquid co2 density at boiling point, kg/m**3
real(kind=ireals), parameter :: dens_o2  = 1.141d+3  ! liquid o2 density at boiling point, kg/m**3
real(kind=ireals), parameter :: dens_n2  = 8.08607d+2 ! liquid n2 density at boiling point, kg/m**3
real(kind=ireals), parameter :: dens_ar  = 1.3928d+3 ! liquid argon density at boiling point, kg/m**3

real(kind=ireals), target, save :: molvol_ch4 = molmass_ch4*1.d-3/dens_ch4 ! methane molar volume, m**3/mol
real(kind=ireals), target, save :: molvol_co2 = molmass_co2*1.d-3/dens_co2 ! co2 molar volume, m**3/mol
real(kind=ireals), target, save :: molvol_o2  = molmass_o2 *1.d-3/dens_o2  ! o2 molar volume, m**3/mol
real(kind=ireals), target, save :: molvol_n2  = molmass_n2 *1.d-3/dens_n2  ! n2 molar volume, m**3/mol
real(kind=ireals), target, save :: molvol_ar  = molmass_ar *1.d-3/dens_ar  ! argon molar volume, m**3/mol

real(kind=ireals), parameter :: mrat_so4_sw = 0.077 !the mass contribution of sulfate to total
                                                    !seawater salinity, kg/kg
real(kind=ireals), parameter :: kgkg_sal_to_molm3_so4 = mrat_so4_sw*1.e+3*row0_/molmass_so4 !multiplier converting
                                                                                           !seawater salinity, kg/kg, to 
                                                                                           !sulfate concentration, mol/m**3

real(kind=ireals), parameter :: nch4_d_nh2o = 1./5.75 ! Molar ratio CH4/H2O in methane hydrate
real(kind=ireals), parameter :: molmass_methhydr = &
& (molmass_h2o + molmass_ch4*nch4_d_nh2o)/(1. + nch4_d_nh2o) ! methane hydrate molar mass, g/mol
real(kind=ireals), save, target :: cpmethhydr = 2160. ! Specific heat at constant pressure, J/(kg*K)
real(kind=ireals), parameter :: methhydrdiss_ = 5.53d+4 ! Methane hydrate dissociation heat, J/mol (Nakagawa et al. 2008)
real(kind=ireals), save, target :: methhydrdiss = methhydrdiss_*1.d+3/(molmass_methhydr) ! Methane hydrate dissociation heat, J/kg
real(kind=ireals), parameter :: alphamh = nch4_d_nh2o*molmass_ch4/molmass_h2o ! parameter used in seperation of water and methane
                                                                    ! after methane hydrate dissociation
real(kind=ireals), save, target :: densmh = 9.d+2 ! Methane hydrate density

real(kind=ireals), parameter :: n2_exp_decay = 50. ! the decay rate in exponential law for nitrogen conecntration in soil, m**(-1)
                                                   ! after (Bazhin, 2001)

real(kind=ireals), parameter :: n2_ocean_ref = 4.46d-1 ! reference nitrogen concentration in seawater, mol/m**3, 
                                             ! according to http://www.seafriends.org.nz/oceano/seawater.htm
real(kind=ireals), parameter :: ar_ocean_ref = 1.d-2   ! reference argon concentration in seawater, mol/m**3, 
                                             ! according to http://www.seafriends.org.nz/oceano/seawater.htm

! Parameters for sediment oxygen demand model by Walker and Snodgrass (1986, J.Env.Eng.)
real(kind=ireals), parameter :: k_O2_SOD = 1.4/molm3tomgl_o2 ! Half-saturation constant for oxygen 
real(kind=ireals), parameter :: mubeta0 = 0.6/(molmass_o2*day_sec) ! maximum aerobic oxidation rate in Monod eq. at 25 deg.C,
                                                                ! range 0.58-5.52
real(kind=ireals), parameter :: kc0 = 0.045/day_sec !reference diffusivity at 20 deg.C, range 0.031 - 0.060 
real(kind=ireals), parameter :: thetaC_SOD = 1.103
real(kind=ireals), parameter :: thetamu_SOD = 1.085

!*               parameters for production
!*    =============================================================
!      real rnpp,rnppmax,rnroot,q100,r0,
!     *     vmax,rkm,q10,
!     *     rke,cmin,punveg,cthresh,
!     *     rkp,tveg,pox,rlmin,rl,rlmax
!*    =============================================================

real(kind=ireals), parameter :: rnpp = 100.          ! 100. - previous value ! NPP (gC/(m**2*month))
real(kind=ireals), parameter :: rnppmax = 180.       ! the maximum value of the NPP
real(kind=ireals), parameter :: q100 = 2.3           ! 2.3 (Liikanen et al., Biogeochemistry, 2002) 
real(kind=ireals), parameter :: r0 = 1.67d-7         ! the constant rate factor, mol/(m**3*s)
!real(kind=ireals), parameter :: r0_oliglake = 1.7d-8  !1.d-8 ! the constant rate factor for oligotrophic lake, mol/(m**3*s)
real(kind=ireals), target, save :: r0_oliglake = 3d-8  ! the constant for CH_4 production below photic zone, mol/(m**3*s)
                                                ! 2.55d-8 - optimal for Shuchi lake (Stepanenko et al., 2011)
real(kind=ireals), target, save :: r0_oliglake1 = 3.d-8  ! the constant for CH_4 production in the photic zone, mol/(m**3*s)
!real(kind=ireals), target, save :: r0_oliglake = 16.31d-10  !1.d-8 ! the constant rate factor for oligotrophic lake, mol/(m**3*s)
real(kind=ireals), parameter :: forg0 = 0.857        ! the constant in depth dependence of methane production rate
real(kind=ireals), parameter :: lambda_new_org = 3. !5. in Walter & Heimann model, 
                                           ! the parameter describing the decrease rate 
                                           ! of methane production term 
                                           ! (due to "new" organics decomposition ~ NPP)
                                           ! with soil depth, m**(-1)

real(kind=ireals), parameter :: r0_oldorg = 7.3d-7 !1.67d-7  ! the constant rate factor for "old" organics decomposition,  
                                                   ! according to formulation based on first-order kinetics, mol/(m**3*s)
real(kind=ireals), parameter :: r0_oldorg_star = 4.5d-8 !1.67d-7  ! the constant rate factor for "old" organics decomposition,  
                                                   ! according to formulation based on first-order kinetics, mol/(kg*s)
real(kind=ireals), parameter :: r0_oldorg2 = 7.3d-8          ! the constant rate factor for "old" organics decomposition, 
                                                   ! according to Michaelis-Menten based formulation, mol/(m**3*s)
!real(kind=ireals), parameter :: r0_oldorg2_star = 1.d-10 !2.d-10 ! the constant rate factor for "old" organics decomposition, 
!                                                   ! according to Michaelis-Menten based formulation, mol/(kg*s)
real(kind=ireals), target, save :: r0_oldorg2_star = 6.9d-11 !2.d-10 ! the constant rate factor for "old" organics decomposition, 
                                                   ! according to Michaelis-Menten based formulation, mol/(kg*s)                                                   
!real(kind=ireals), target, save :: r0_oldorg2_star = 1.758d-10 !2.d-10 ! the constant rate factor for "old" organics decomposition, 
                                                   ! according to Michaelis-Menten based formulation, mol/(kg*s)         
real(kind=ireals), parameter :: alpha_old_org = 2.d-1 !2.d-1 ! the constant of decomposition, year**(-1), value needs to be verified (!)
real(kind=ireals), parameter :: k_oldorg = 3.d-1 ! the constant in the denominator of Michaelis-Menten equation when deriving 
                                       ! the equation for methane generation from old organics decomposition, kg/m**3
                                       ! varies from 0.1 to 0.3 (rough interval)
real(kind=ireals), parameter :: V_oldorg = 0.2d-2 !0.2d-2 ! the constant in the numinator of Michaelis-Menten equation when deriving 
                                       ! the equation for methane generation from old organics decomposition, kg/m**3/year

real(kind=ireals), parameter :: O2inhib_conc = 10./molm3toppm_o2 !10 ppm is a critical oxygen concentration at which
                                                                 !methanogenesis is completely suppressed (Borrel et al., 2011) 
real(kind=ireals), parameter :: O2inhib_const = (100. - 1.)/O2inhib_conc !Complete inhibition of methanogenesis in the model
                                                                         !by O2 is simulated as 100 times decrease of production
real(kind=ireals), parameter :: O2_exp_decay = 2.d+2*log(10.) ! the decay rate in exponential law for oxygen concentration 
                                                              ! in soil, m**(-1), ensuring 100 times decrease in 1 centimeter

real(kind=ireals), parameter :: SO4inhib_conc = 30./molm3tomcM ! 30 mcmol/l is a critical sulfate concentration at which
                                                               ! methanogenesis is completely outcompeted by 
                                                               ! sulfate reducers  (Lovley and Klug, 1986) 
real(kind=ireals), parameter :: SO4_exp_decay = 2.d+2*log(10.) ! the decay rate in exponential law for sulfate concentration 
                                                               ! in soil, m**(-1), ensuring 100 times decrease in 1 centimeter

real(kind=ireals), parameter :: CH4_exp_growth = 2.d+2*log(10.) ! the growth rate in exponential law for methane concentration 
                                                                ! in soil, m**(-1), ensuring 100 times increase in 1 centimeter

real(kind=ireals), parameter :: C0_oldorg = 18.  ! the density of organics sequestered under the talik, kg/m**3
real(kind=ireals), parameter :: ice_trap_bubbl_meth_decr = 0.625 ! the fraction of methane that remains in bubbles
                                                       ! when they have been trapped by the ice cover;
                                                       ! the decrease of methane happens due to oxidation and 
                                                       ! the exchange with gases dissolved in water
real(kind=ireals), parameter :: methox2sod = 0. ! a fraction of SOD, utilized for methane oxidation in the top of bottom sediments,
                                       ! 0.1 - 0.64 (Liikanen et al., 2002 and references therein)

integer(kind=iintegers), parameter :: ngasb = 5 ! Number of gases considered in a bubble
       
!*    =============================================================
!*               parameters for oxidation
!*             =============================

real(kind=ireals), parameter :: vmax = 45./3600./1000. ! Michaelis-Menten, mol/(m**3*s) !(muM/s)
real(kind=ireals), parameter :: rkm = 5./1000.         ! coefficients, mol/(m**3*s) !(muM)
real(kind=ireals), parameter :: q10 = 2.               ! the observed value for oxidation, n/d

real(kind=ireals), parameter :: vq_max = 6.d-4 ! maximum methane oxidation potential (10 deg C), 
                                     ! after Arah&Stephen (1998), Watson et al. (1997)
real(kind=ireals), parameter :: temp0 = 283. ! reference temperature in Arrhenius equation, K
real(kind=ireals), parameter :: k_ch4 = 0.6 / molm3tomgl_ch4 !(Liikanen et al. 2002; Lofton et al. 2013)
                                   !0.44 ! Michaelis constant for methane in methane oxidation,
                                   ! after Arah&Stephen (1998), Nedwell&Watson (1995) 
real(kind=ireals), parameter :: k_o2 = 0.672 / molm3tomgl_o2 ! (Lidstrom and Somers, 1984)
                                   !0.33  ! Michaelis constant for oxygen in methane oxidation, 
                                   ! after Arah&Stephen (1998), Nedwell&Watson (1995)
real(kind=ireals), parameter :: delta_Eq = 5.d+4 ! activation energy for methane oxidation, J/mol, 
                                       ! after Arah&Stephen (1998), Dunfield et al. (1993),
                                       ! Nedwell&Watson (1995)
real(kind=ireals), parameter :: Vmaxw = 1.d-1/86400. ! reaction potential in oxygen-saturated Michaelis-Menten kinetics, 
                                          ! after (Liikanen et al., 2002)
real(kind=ireals), parameter :: k_ch40 = 1.d-2      ! half-saturation constant in oxygen-saturated Michaelis-Menten kinetics, 
                                          ! after (Liikanen et al., 2002)

real(kind=ireals), parameter :: koxyg = 0.38/86400. ! Methane oxidation constant in 1-st order kinetics, (Striegl et al., 1998)

real(kind=ireals), parameter :: k_ch4_soil = 9.5/molm3tomcM !CH4 half-saturation constant in M-M kinetics for top 1 cm sediments
                                                            !in Lidstrom and Sommers, 1984

real(kind=ireals), parameter :: k_o2_soil  = 21./molm3tomcM !O2 half-saturation constant in M-M kinetics for top 1 cm sediments
                                                            !in Lidstrom and Sommers, 1984
real(kind=ireals), parameter :: Vmax_soil = 40./molm3tomcM/3600. ! reaction potential in Michaelis-Menten methane oxidation in soil, 
                                                                 ! after (Lidstrom and Sommers, 1984)      

real(kind=ireals), parameter :: CO2O2_prod = 1., CO2O2_resp = 1. ! Stoichiometry ratios
real(kind=ireals), parameter :: CO2O2_bod  = 1., CO2O2_sod  = 1. ! Stoichiometry ratios

!*     =============================================================
!*                parameters for ebullition
!*              =============================

real(kind=ireals), parameter :: rke = 1./3600.      ! a rate constant of the unit, 1/s
real(kind=ireals), parameter :: cmin = 500./1000.   ! the concentration at which bubble formation occurs, mol/(m**3*s) !(muM)
real(kind=ireals), parameter :: punveg = 0.         ! the percentage of unvegetated, bare soil
real(kind=ireals), parameter :: cthresh = cmin*(1.+punveg/100.) !the threshold concentration for bubble formation
real(kind=ireals), parameter :: rel_conc_ebul_crit = 0.4 !0.4, 0.45 the threshold relative concentration of methane (Wania, 2007)
                                          ! 1. is taken to prevent in-ground diffusion flux in East-Siberian sea experiments
                                          ! for ebullition; 
                                          ! relative concentration = concentration of bubble formation divided by 
                                          ! maximal concentration according to Henry law at a given
                                          ! pressure
real(kind=ireals), parameter :: meth_ebul_ice_intercept = 0.9 ! The portion of gas (methane) bubbles
                                                    ! intercepted by the ice cover during
                                                    ! a wintertime, n/d

!*     =============================================================
!*            parameters for plant-mediated transport
!*          ===========================================

real(kind=ireals), parameter :: rkp = 0.01/3600.    ! a rate constant of the unit 0.01/s
real(kind=ireals), parameter :: tveg = 15.          ! a factor describing the quality of plant-mediated transport at a site, n/d
real(kind=ireals), parameter :: pox = 0.5           ! a certain fraction of methane oxidized in
                                          ! the oxic zone around the roots of plants, n/d
real(kind=ireals), parameter :: rlmin = 0.          ! the parameter for fgrow(t), n/d
real(kind=ireals), parameter :: rl = 4.             ! the parameter for fgrow(t), n/d
real(kind=ireals), parameter :: rlmax = 4.          ! rlmax=rlmin+rl, n/d

real(kind=ireals), parameter :: mfs = molmass_ch4*8.64d+7 ! the transform multiplier for methane flux
                                                ! from mol/(m**2*s) to mg/(m**2*day)
real(kind=ireals), parameter :: mol2mcmol = 1.d+6 !mol to micromol 
real(kind=ireals), parameter :: mol2mmol  = 1.d+3 !mol to millimol 
real(kind=ireals), parameter :: mf = molmass_ch4*1.d+3    ! the transform multiplier for methane amount
                                                ! from mol/(m**2) to mg/(m**2)

real(kind=ireals), parameter :: photic_threshold = 0.2 ! The fraction of surface irradiance that defines the 
                                                       ! depth of photic zone


! ===================================================
! Parameters of Stefan and Fang model for DO and CO_2
! ===================================================

real(4), parameter :: T0 = 20. ! Reference temperature, Celsius
real(4), parameter :: T00 = 10. ! Another reference temperature, Celsius
real(4), parameter :: c1_pmax = 9.6, c2_pmax = 1.036 ! Constants for temperature dependence,
                                                     ! (Stefan and Fang, 1994; Megard et al., 1984)
real(4), parameter :: k1_c1 = 0.687, k1_c2 = 1.086   ! (Megard et al., 1984)
real(4), parameter :: k2_c1 = 5., k2_c2 = 15.
real(4), parameter :: YCHO2 = 8.E-3 ! mass ratio of chlorophyll-a to oxygen utilized in respiration
real(4), parameter :: k_r = 0.1 ! day**(-1), (Brown and Barnwell, QUAL2E, 1987; Riley, 1988)
real(4), parameter :: theta_r = 1.045 ! temperature dependence of respiration (Ambrose et al., 1988)
real(4), parameter :: k_b = 0.1 ! day**(-1), 1-st order decay coefficient, (Riley, 1988; Brown and Barnwell, 1987)
real(4), parameter :: theta_b1 = 1.047, theta_b2 = 1.13 ! temperature-dependence coefficient for biochemical oxygen demand
real(4), parameter :: mO2C_dec = 8./3., mCChla_dec = 30. ! Stoichiometry constants for organics decay in the water column, 
                                                         ! for BOD, basing on carbonaceous BOD and neglecting 
                                                         ! nitropgenous BOD, see details in Stefan and Fang, 1994
real(4), parameter :: Sb20 = 0.5E+3 ! estimates: 0.5E+3 for oligotrophic, 1.5E+3 for eutrophic lakes, Stefan and Fang, 1994,
                                   ! mg/(m**2*day), reference sedimentary oxygen demand (SOD)
real(4), parameter :: theta_s1 = 1.065, theta_s2 = 1.13 ! temperature dependence coefficients for SOD, Stefan and Fang, 1994
real(4), parameter :: T_md = 4. ! The temperature of maximal water density, Celsius


! =====================================
! Parameters of Hanson et al. 2004 model for DIC, DOC, POC
! =====================================

real(kind=ireals), parameter :: alpha_POCL = 0.8 !respiration scaling in respect to photosynthesis, n/d
real(kind=ireals), parameter :: tau_POCD = 20.*day_sec !timescale for dead POC decay, s
real(kind=ireals), parameter :: tau_DOC = 200.*day_sec !timescale for DOC decay, s
real(kind=ireals), parameter :: dpart = 5.E-6 !Diameter of particles, m (Wetzel, 2001)
                                              !have to be checked for "living particles"
real(kind=ireals), parameter :: beta_POCL = 0.03 !the scaling factor of exudation rate to photosynthesis
real(kind=ireals), parameter :: tau_Dh_epi = 33.*day_sec !time scale for death of POCL in epilimnion
real(kind=ireals), parameter :: tau_Dh_hypo = 1.1*day_sec !time scale for death of POCL in hypolimnion

END MODULE METH_OXYG_CONSTANTS


MODULE BL_MOD_LAKE

use LAKE_DATATYPES, only : ireals
use NUMERIC_PARAMS, only : ms

implicit none

real(kind=ireals) :: ST
real(kind=ireals) :: PGR
real(kind=ireals) :: TGROLD, QGROLD, RADIAT, WSOLD, SNOLD
real(kind=ireals) :: ZS
real(kind=ireals) :: thsoil, whsoil
real(kind=ireals) :: bOLD
real(kind=ireals) :: RF1, RF2
real(kind=ireals) :: SNMELT
real(kind=ireals) :: HSNOW
real(kind=ireals) :: HS, ES
real(kind=ireals) :: TGRNEW, QGRNEW, WSNEW, SNNEW
real(kind=ireals) :: RUNOFF
real(kind=ireals) :: ElatOld, HFold, PRSold
real(kind=ireals), allocatable :: extinct(:)

contains
SUBROUTINE BL_MOD_LAKE_ALLOC
implicit none
allocate(extinct(1:ms))
END SUBROUTINE BL_MOD_LAKE_ALLOC

SUBROUTINE BL_MOD_LAKE_DEALLOC
implicit none
deallocate(extinct)
END SUBROUTINE BL_MOD_LAKE_DEALLOC

END MODULE BL_MOD_LAKE



MODULE SNOWSOIL_MOD

use LAKE_DATATYPES, only : ireals, iintegers
use NUMERIC_PARAMS, only : ms, ML

implicit none

real(kind=ireals), allocatable :: lams(:), q(:) !former common-block watericesnowarr
real(kind=ireals), allocatable :: AL(:),DLT(:),DVT(:),ALLL(:),DL(:) !former common-block soilsol
real(kind=ireals), allocatable :: ALV(:),DV(:),Z(:),T(:),WL(:),WV(:),WI(:),dens(:) !former common-block soilsol
real(kind=ireals), allocatable :: dz(:)!former common-block soildat
real(kind=ireals), allocatable :: Tsn(:),cs(:)!former common-block snow_char

integer(kind=iintegers) :: itop!former common-block soildat


contains
SUBROUTINE SNOWSOIL_MOD_ALLOC
implicit none
allocate(lams(1:ms), q(1:110))
allocate(AL(1:ML),DLT(1:ML),DVT(1:ML),ALLL(1:ML),DL(1:ML) )
allocate(ALV(1:ML),DV(1:ML),Z(1:ML),T(1:ML),WL(1:ML),WV(1:ML),WI(1:ML),dens(1:ms))
allocate(dz(1:ms))
allocate(Tsn(1:ms),cs(1:ms))
END SUBROUTINE SNOWSOIL_MOD_ALLOC

SUBROUTINE SNOWSOIL_MOD_DEALLOC
implicit none
deallocate(lams, q)
deallocate(AL, DLT, DVT, ALLL, DL )
deallocate(ALV, DV, Z, T, WL, WV, WI, dens)
deallocate(dz)
deallocate(Tsn, cs)
END SUBROUTINE SNOWSOIL_MOD_DEALLOC

END MODULE SNOWSOIL_MOD

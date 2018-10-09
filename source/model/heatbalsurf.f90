FUNCTION HEATBALANCE(Tsurf,surftyp,RadWater,RadIce,fetch,dt)

use LAKE_DATATYPES, only : ireals, iintegers

 use ATMOS,  only : &
 & Radbal      , &
 & Radbal_surf , &
 & hflux       , &
 & Elatent     , &
 & eflux       , &
 & eflux0      , &
 & eflux0_kinem, &
 & hskin, sabspen

 use ARRAYS_WATERSTATE, only: Tw1, Tw2, Ti1, Ti2, Sal2, preswat
 use ARRAYS, only : nstep
 use ARRAYS_BATHYM, only : h1, l1

 use DRIVING_PARAMS, only : skin

 use PHYS_CONSTANTS, only : &
 & cw, &
 & row0, &
 & cw_m_row0
 
 use SKIN_MOD, only : TSKIN_SOLVER

 use RADIATION, only : rad_type
 
 implicit none

 real(kind=ireals) :: HEATBALANCE

 integer(kind=iintegers), parameter :: water_indic = 4

 real(kind=ireals), intent(in) :: Tsurf
 real(kind=ireals), intent(in) :: dt
 !real(kind=ireals), intent(in) :: extwat
 !real(kind=ireals), intent(in) :: extice
 type(rad_type), intent(in) :: RadWater
 type(rad_type), intent(in) :: RadIce

 real(kind=ireals), intent(in) :: fetch
 
 integer(kind=iintegers), intent(in) :: surftyp

 real(kind=ireals) :: dHdt
 real(kind=ireals) :: convect_flux

 real(kind=ireals), external :: SHORTBAL_TOPHALFLAYER

 if     (skin%par == 1 .and. surftyp == water_indic) then
   call TSKIN_SOLVER(Tsurf,RadWater,RadIce,fetch,dt,hskin,sabspen)
   call HEATFLUXSURF(surftyp,dt,dHdt,convect_flux)

   HEATBALANCE = dHdt - (SHORTBAL_TOPHALFLAYER(RadWater,sabspen) + hskin) + &
   & eflux + convect_flux
   eflux0 = hskin
   eflux0_kinem = eflux0 /(cw_m_row0)
 else ! if (skin%par == 0) then
   call HEATFLUXSURF     (surftyp,dt,dHdt,convect_flux)
   call SENSLATMOM_FLUXES(Tsurf,fetch)
   call RADBALANCE       (Tsurf,RadWater,RadIce,surftyp)

   HEATBALANCE = dHdt - (Radbal - hflux - Elatent) + &
   & eflux + convect_flux
   eflux0      = Radbal_surf - hflux - Elatent
   eflux0_kinem = eflux0/(cw_m_row0)
!   if (nstep >= 10) then
!     print*, 'Bal ', dHdt, Radbal, hflux, Elatent, &
!     & eflux, convect_flux
!     print*, 'Tw1 ', Tw1
!     print*, 'Tw2 ', Tw2
!     print*, 'Ti1 ', Ti1
!     print*, 'Ti2 ', Ti2
!     print*, 'h1 ', h1, 'l1 ', l1
!     print*, 'Sal2 ', Sal2(1), 'preswat ', preswat(1)
!     read*
!   endif
 endif

 RETURN
 END FUNCTION HEATBALANCE


 SUBROUTINE HEATFLUXSURF(surftyp,dt,dHdt,convect_flux)

 use LAKE_DATATYPES, only : ireals, iintegers

 use ATMOS, only: &
 & botflux, &
 & eflux

 use ARRAYS, only : dt_inv
 use ARRAYS_SOIL, only : Tsoil2,Tsoil1,csoil,rosoil,lamsoil,dzs
 use ARRAYS_WATERSTATE, only : Tw2,Tw1,Ti2,Ti1,lamw,lami_v,ci_m_roi_v
 use ARRAYS_TURB, only : PEMF, pt_down_f
 use ARRAYS_GRID, only : nsoilcols,ddz,ddzi
 use ARRAYS_BATHYM, only : h1,l1,hs1,ls1,dhw0,dhi0,bathymwater, bathymice

 use PHYS_CONSTANTS, only: &
 & lami,ci,roi,cw,row0, &
 & ci_m_roi, cw_m_row0

 use NUMERIC_PARAMS, only: &
 & ms,ML

 use DRIVING_PARAMS, only: M

 use SNOWSOIL_MOD, only : Tsn, T, cs, lams, dz, itop, dens

 implicit none

!Parameters
 integer(kind=iintegers), parameter :: soil_indic = 1
 integer(kind=iintegers), parameter :: ice_indic = 2
 integer(kind=iintegers), parameter :: snow_indic = 3
 integer(kind=iintegers), parameter :: water_indic = 4

!Input variables
 integer(kind=iintegers), intent(in) :: surftyp
 real(kind=ireals), intent(in) :: dt

!Output variables
 real(kind=ireals), intent(out) :: dHdt
 real(kind=ireals), intent(out) :: convect_flux

 !integer(kind=iintegers) itop 
 !real(kind=ireals) AL,DLT,DVT,ALLL,DL, &
 !& ALV,DV,Z,T,WL,WV,WI,dens,dz
 !common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML), &
 !& ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
 !common /SOILDAT/ dz(ms),itop

 SELECT CASE (surftyp)   
  case (soil_indic)
    eflux = - (lamsoil(1) + lamsoil(2))*0.25d0 * &
    & (Tsoil2(2) + Tsoil1(2,nsoilcols) - &
    &  Tsoil2(1) - Tsoil1(1,nsoilcols))/dzs(1)
  case (ice_indic)
    eflux = - 0.5d0*lami_v(1)/bathymice(1)%area_int * bathymice(1)%area_half * &
    & (Ti2(2) + Ti1(2) - Ti2(1) - Ti1(1))/(ddzi(1)*l1) 
    ! if (l1<0.1) eflux = -lami*(Ti2(2)-Ti2(1))/(ddz*0.1) 
  case (snow_indic)
    eflux = - (lams(itop) + lams(itop+1))*0.25d0 * &
    & (Tsn(itop+1) + T(itop+1) - Tsn(itop) - T(itop))/dz(itop)
  case (water_indic)  
    eflux = - 0.5d0*lamw(1)/bathymwater(1)%area_int * bathymwater(1)%area_half * &
    & (Tw2(2) + Tw1(2) - Tw2(1) - Tw1(1))/(ddz(1)*h1) 
    ! if (h1<0.1) eflux = -lamw(1)*(Tw2(2)-Tw2(1))/(ddz*0.1) 
 END SELECT

 SELECT CASE (surftyp)   
  case (soil_indic)
    dHdt = 0.5d0*dt_inv*csoil(1)*rosoil(1)*dzs(1) * &
    & (Tsoil2(1) - Tsoil1(1,nsoilcols))
  case (ice_indic)
    dHdt = 0.5d0*dt_inv*ci_m_roi_v(1)*(ddzi(1)*l1*(Ti2(1) - Ti1(1)) + &
    &    0.5d0*dhi0*(Ti2(2) + Ti1(2) - Ti2(1) - Ti1(1)) )
  case (snow_indic)
    dHdt = 0.5d0*dt_inv*cs(itop)*dens(itop)*dz(itop)*(Tsn(itop) - T(itop))
  case (water_indic)  
    dHdt = 0.5d0*dt_inv*cw_m_row0*(ddz(1)*h1*(Tw2(1) - Tw1(1)) + &
    &    0.5d0*dhw0*(Tw2(2) + Tw1(2) - Tw2(1) - Tw1(1)) )
 END SELECT

 SELECT CASE (surftyp)
   case(soil_indic)
     convect_flux = 0.d0
   case(ice_indic)
     convect_flux = 0.d0
   case(snow_indic)
     convect_flux = 0.d0
   case(water_indic)
     convect_flux = 0.d0 !cw_m_row0*PEMF(1)*(pt_down_f(1)-0.5d0*(Tw2(2)+Tw2(1)) )
 END SELECT

 botflux = - 0.5d0*lamw(M)/bathymwater(M+1)%area_int * bathymwater(M)%area_half * &
 & (Tw2(M+1) + Tw1(M+1) - Tw2(M) - Tw1(M))/(ddz(M)*h1)

 RETURN
 END SUBROUTINE HEATFLUXSURF


 SUBROUTINE SENSLATMOM_FLUXES(Tsurf,fetch)
 
!Subroutine SENSLATMOM_FLUXES calculates
!sensible, latent heat and momentum fluxes

use LAKE_DATATYPES, only : ireals, iintegers

 use ATMOS, only : &
 & hw_input, xlew_input, cdmw_input, &
 & hw, xlew, cdmw, velfrict, &
 & uwind, vwind, wind10, zref, &
 & tempair, pressure, humair, &
 & hflux, Elatent, tau, &
 & u, v, urel, vrel

 use DRIVING_PARAMS, only : &
 & missing_value , &
 & relwind , & 
 & PBLpar , &
 & waveenh , &
 & c_d, &
 & sensflux0, momflux0

 use PHYS_CONSTANTS, only : &
 & Rd,cp,Lwv,Liv,g

 use ARRAYS, only : roughness,aM,bM,relhums
 use ARRAYS_BATHYM, only : h1,l1

 use SURF_SCHEME1_MOD, only : DRAGVL
 use SURF_SCHEME_INM, only : DRAG3_LAKE

 use SFCFLX, only : SFCFLX_MOMSENLAT

 implicit none

 ! Local variables  
 real(kind=ireals), parameter :: minwind = 1.d-2

 real(kind=ireals), intent(in) :: Tsurf
 real(kind=ireals), intent(in) :: fetch
 
 real(kind=ireals) :: TET2,TET1,ro,esatsurf,humsurf,bx(7),bix(11), &
 & c_u,c_t,ksw,hwave,SHF,LHF
 real(kind=4) :: bx4(7),bix4(11)
 real(kind=ireals) :: wr
 real(kind=8) :: xx, yy, zz, xx1 ! For interface with FLake routines
 real(kind=ireals) :: deltaR, deltaT, wr_st, wr_n
 
 data ksw /2.d0/
 
 integer(kind=iintegers) itdrag
 
 SAVE
 
!(urel,vrel) is wind, relative to lake currents

 if (relwind%par == 1) then
   urel = uwind - u
   vrel = vwind - v
 elseif (relwind%par == 0) then
   urel = uwind
   vrel = vwind
 endif

 wr       = max(sqrt(urel**2+vrel**2), minwind)
 wind10   = wr*log(10./roughness)/log(zref/roughness) 
 TET2     = (tempair+273.15)*(100000./pressure)**0.286
 TET1     = (Tsurf  +273.15)*(100000./pressure)**0.286  
 ro       = pressure/(Rd*(tempair+273.15))  
 esatsurf = 610.7*10.**(aM*Tsurf/(bM+Tsurf))
 humsurf  = 0.622/pressure*esatsurf*relhums
 

 if (PBLpar%par == -1) then
!Sensible, latent heat and momentum fluxes are prescribed outside the lake model
   hflux   = hw_input
   Elatent = xlew_input
   tau = cdmw_input*ro*wr
 endif
 
 if (PBLpar%par == 0) then
!In this case the fluxes are set to constants, specified in input file
   hflux   = sensflux0%par
   Elatent = 0.d0
   tau     = momflux0%par
 endif

 if (PBLpar%par == 1) then
 
!Businger-Dayer, Beljaars parameterization for exchange coefficients

   itdrag= 10
   bx(1) = wr !wind
   bx(2) = TET2
   bx(3) = TET1
   bx(4) = humair 
   bx(5) = humsurf
   bx(6) = zref
   bx(7) = roughness 
   call DRAGVL (bx, bix, itdrag)
   c_u = bix(8) !*sqrt(0.5_ireals)! - sensitivity experiment!
   c_t = bix(9) !*sqrt(0.5_ireals)

   !deltaT=max(0.,-(TET2-TET1))
   !deltaR=max(0.,-(humair-humsurf))
   !wr_st=max(0.01,(zref*c_u*c_t*wr*(deltaT+0.61*(Tsurf+273.15)*deltaR)))
   !wr_n = sqrt(wr**2 + 9*(0.333*wr_st)**0.66)

   hflux   = -cp*ro*c_u*c_t*wr*(TET2-TET1)
   if (l1 == 0) then
     ! No ice and snow cover above lake: 
     ! the latent heat of evaporation is used
     Elatent = - Lwv*ro*c_u*c_t*wr*(humair-humsurf) 
   else
     ! There is ice and (probably) snow cover above lake: 
     ! the latent heat of sublimation is used
     Elatent = - Liv*ro*c_u*c_t*wr*(humair-humsurf)
   endif
   tau = ro*c_u**2*wr**2
 endif

 ! Parameterization of exchange coefficients acording to (Louis, 1979)

 if (PBLpar%par == 2) then
   call RichPBL(hflux,Elatent,tau,Tsurf,humsurf,roughness)
 endif

 if (PBLpar%par == 3) then
 ! The parameterization of fluxes from lake model Flake (Mironov et al., 2006)
   call SFCFLX_MOMSENLAT(dble(zref)    , dble(zref)    , dble(fetch)  , &
   &    dble(wr) , dble(tempair + 273.15) , dble(humair)  , dble(Tsurf + 273.15)  , &
   &    dble(pressure) , dble(l1)       , xx1     , yy           , &
   &    zz  , xx)
   hflux = yy
   Elatent = zz
   tau = - xx1
 endif 
 
 if (PBLpar%par == 4) then
   call SURF_SCHEME3(wr,TET2,TET1,humair,humsurf,roughness,zref,ro,l1,  &
   & hflux, Elatent,tau)
 endif

 if (PBLpar%par == 5) then

   !Surface flux scheme from INM RAS climate model
   itdrag= 10
   bx4(1) = wr !wind
   bx4(2) = TET2 - TET1
   bx4(3) = 0.5*(TET1 + TET2)
   bx4(4) = humair - humsurf 
   bx4(5) = zref
   bx4(6) = roughness
   call DRAG3_LAKE(bx4,itdrag,bix4)
   c_u = bix4(8) !*sqrt(0.5_ireals)! - sensitivity experiment!
   c_t = bix4(9) !*sqrt(0.5_ireals)

   hflux   = -cp*ro*c_u*c_t*wr*(TET2-TET1)
   if (l1 == 0) then
     ! No ice and snow cover above lake: 
     ! the latent heat of evaporation is used
     Elatent = - Lwv*ro*c_u*c_t*wr*(humair-humsurf) 
   else
     ! There is ice and (probably) snow cover above lake: 
     ! the latent heat of sublimation is used
     Elatent = - Liv*ro*c_u*c_t*wr*(humair-humsurf)
   endif
   tau = ro*c_u**2*wr**2

 endif

 ! Constant drag coefficient
 if (c_d%par /= missing_value) then
   tau = ro*c_d%par*wr**2
 endif
 
!SHALLOW WATER EFFECT ON HEAT, MOISTURE AND MOMENTUM FLUXES (PANIN ET. AL., 2006)
 
 if (l1 == 0.and.h1 /= 0.and.waveenh%par == 1.and.(.not.PBLpar%par == 0)) then
   hwave   = 0.07*wind10**2*(g*h1/max(wind10,1.d0)**2)**(3./5.)/g
   hflux   = hflux   + hflux   * ksw*hwave/h1
   Elatent = Elatent + Elatent * ksw*hwave/h1
   tau     = tau     + tau     * ksw*hwave/h1
 endif
 cdmw = tau/(ro*wr)
 xlew = Elatent
 hw   = hflux
     
 velfrict = sqrt(tau/ro)    
      
! CALCULATED sensible and latent heat fluxes are used in heat balance

!SHF=hflux
!LHF=Elatent
      
! MEASURED sensible and latent heat fluxes are used in heat balance
!SHF = Hm
!LHF = LEm

 RETURN  
 END SUBROUTINE SENSLATMOM_FLUXES



 SUBROUTINE RADBALANCE(Tsurf,RadWater,RadIce,surftyp)

use LAKE_DATATYPES, only : ireals, iintegers

 use ATMOS, only: &
 & longwave    , &
 & Radbal      , &
 & Radbal_surf , &
 & surfrad     , &
 & shortwave   , &
 & humair      , &
 & pressure    , &
 & tempair     , &
 & cloud

 use DRIVING_PARAMS, only: skin, ifrad, missing_value

 use PHYS_CONSTANTS, only: &
 & sigma       , &
 & sabs        , &
 & Rd_d_Rwv    , &
 & Kelvin0

 use ARRAYS, only: roughness,emissivity,albedo,albedo_lw,bM,relhums
 use ARRAYS_GRID, only : ddz,ddzi
 use ARRAYS_BATHYM, only : h1,l1,hs1,ls1

 use PHYS_FUNC, only : &
 & EXTINCT_SNOW, &
 & NETLWRAD

 use RADIATION, only : rad_type

 implicit none

!Parameters
 integer(kind=iintegers), parameter :: soil_indic = 1
 integer(kind=iintegers), parameter :: ice_indic = 2
 integer(kind=iintegers), parameter :: snow_indic = 3
 integer(kind=iintegers), parameter :: water_indic = 4

 real(kind=ireals),    intent(in) :: Tsurf
 !real(kind=ireals),    intent(in) :: extwat
 !real(kind=ireals),    intent(in) :: extice
 type(rad_type), intent(in) :: RadWater
 type(rad_type), intent(in) :: RadIce
 
 integer(kind=iintegers), intent(in) :: surftyp

 real(kind=ireals) :: longwave1, epsa

 
 if     (ifrad%par == 1) then
   surfrad = emissivity*sigma*(Tsurf+273.15)**4
 elseif (ifrad%par == 0) then
   surfrad = 0.d0
 endif

 if (longwave == missing_value) then
!  Longwave radiation is missing in atmospheric forcing
   epsa = humair*pressure/Rd_d_Rwv
   longwave1 = NETLWRAD(Tsurf+Kelvin0,tempair+Kelvin0,epsa,cloud,emissivity) + surfrad
 else
   longwave1 = longwave
 endif

 SELECT CASE(surftyp)
 CASE(soil_indic)
   Radbal = shortwave*(1.d0-albedo)+longwave1*(1.d0-albedo_lw)-surfrad
 CASE(ice_indic)
   Radbal = & !shortwave*(1.d0-albedo)*(1.d0-exp(-extice*0.5d0*ddzi(1)*l1)) + &
   & RadIce%integr(0) - RadIce%integr(1) + &
   & longwave1*(1.d0-albedo_lw)-surfrad
 CASE(snow_indic)
   Radbal = longwave1*(1.d0-albedo_lw) - surfrad
!  & + shortwave*(1-albedo)*
!  & (1-EXTINCT_SNOW(dens(itop))**(0.5d0*dz(itop)) ) +
 CASE(water_indic)
   if     (skin%par == 1) then
!  In this case Radbal is the net radiation balance of the
!  surface skin layer
!    Radbal      = shortwave*(1.d0-albedo)*sabs + & ! previous treatment of albedo
     Radbal      = shortwave*sabs + longwave1*(1.d0-albedo_lw) - surfrad
!    In this case Radbal_surf is the net radiation balance at the
!    top surface of skin layer
!    Radbal_surf = shortwave*(1.d0-albedo) + & ! not used if skin is included
!    & longwave1*(1.d0-albedo_lw) - surfrad
   elseif (skin%par == 0) then
!    In this case Radbal is the net radiation balance of the
!    top half layer of the water coloumn grid, assuming that
!    the sabs part of shortwave radiation is absorbed at the top of it
!    Radbal = shortwave*(1.d0-albedo)*(1.d0-(1.d0-sabs) * & ! previous treatment of albedo
!    & exp(-extwat*0.5d0*ddz(1)*h1))+longwave1-surfrad
    Radbal = shortwave*sabs + &
    !& shortwave*(1.d0-sabs)*(1.d0-albedo)*(1.d0 - exp(-extwat*0.5d0*ddz(1)*h1)) + &
    RadWater%integr(0) - RadWater%integr(1) + &
    & longwave1*(1.d0-albedo_lw)-surfrad
!    In this case Radbal_surf is the net radiation balance at the
!    outer boundary of top half layer
!    Radbal_surf = shortwave*(1.d0-albedo)*sabs + & ! previous treatment of albedo
    Radbal_surf = shortwave*sabs + longwave1*(1.d0-albedo_lw) - surfrad
   endif
 END SELECT

 RETURN
 END SUBROUTINE RADBALANCE

 
 FUNCTION SHORTBAL_TOPHALFLAYER(RadWater,sabspen)

!The function SHORTBAL_TOPHALFLAYER computes
!the shortwave radiative balance of the top half layer
!of water coloumn, omitting the sabs fraction,
!absorbed by the skin layer aloft

use LAKE_DATATYPES, only : ireals, iintegers

 use ATMOS,  only: &
 & shortwave

 use ARRAYS, only: roughness,emissivity,albedo,aM,bM,relhums
 use ARRAYS_GRID, only : ddz
 use ARRAYS_BATHYM, only: h1

 use PHYS_CONSTANTS, only: &
 & sabs,sabs0

 use RADIATION, only : rad_type

 implicit none
 
 real(kind=ireals) :: SHORTBAL_TOPHALFLAYER
 
 real(kind=ireals), intent(in) :: sabspen

 type(rad_type), intent(in) :: RadWater

 SHORTBAL_TOPHALFLAYER = &
 !& shortwave*(1._ireals-albedo)*(1._ireals-sabs0)*(1._ireals-sabspen)* &
 !& (1._ireals-exp(-extwat*0.5_ireals*ddz(1)*h1))
 & (RadWater%integr(0) - RadWater%integr(1))*(1._ireals-sabspen)

 END FUNCTION SHORTBAL_TOPHALFLAYER

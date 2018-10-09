MODULE SKIN_MOD

use LAKE_DATATYPES, only : ireals, iintegers
use INOUT, only : CHECK_UNIT

contains
SUBROUTINE TSKIN_SOLVER(Tsurf,RadWater,RadIce,fetch,dt,hskin,sabspen)

!The subroutine TSKIN_SOLVER performs iterations to
!find the skin temperature from heat balance equation
 use DRIVING_PARAMS, only : path
 use ARRAYS, only: &
 & Tskin, &
 & deltaskin
 use PHYS_CONSTANTS, only : &
 & sabs,sabs0,albedoofwater
 use INOUT_PARAMETERS, only : &
 & lake_subr_unit_min, &
 & lake_subr_unit_max
 use RADIATION, only : rad_type

 implicit none

 real(kind=ireals), intent(in) :: Tsurf
 !real(kind=ireals), intent(in) :: extwat
 !real(kind=ireals), intent(in) :: extice
 type(rad_type), intent(in) :: RadWater
 type(rad_type), intent(in) :: RadIce
 real(kind=ireals), intent(in) :: fetch
 real(kind=ireals), intent(in) :: dt
 
 real(kind=ireals), intent(out) :: hskin
 real(kind=ireals), intent(out) :: sabspen

 real(kind=ireals) :: Tskin1,Tskin2,Tskin3,Balskin1,Balskin2,Balskin3
 real(kind=ireals) :: interval2
 
 integer(kind=iintegers) :: iter 
 integer(kind=iintegers) :: maxiter = 100
 integer(kind=iintegers) :: nunit_ = lake_subr_unit_min

! real(kind=ireals) :: HEATBALANCESKIN
! real(kind=ireals) :: SKIN_THICKNESS
! real(kind=ireals) :: FRAC_SOLAR_ABSORB_SKIN

 logical :: firstcall = .true.

 if (firstcall) then
   call CHECK_UNIT(lake_subr_unit_min,lake_subr_unit_max,nunit_)
   open (nunit_,file=path(1:len_trim(path))//'results/debug/iter_Tskin.dat')
 endif

 deltaskin = SKIN_THICKNESS()
 !The fraction of incident shortwave radiation absorbed in the skin layer
 !consists of its near ifrared part (sabs0) and the penetrated 
 !shortwave radiation absorbed inside the skin
 sabspen = FRAC_SOLAR_ABSORB_SKIN(deltaskin)
 sabs = sabs0 + (1. - sabs0)*(1.-albedoofwater)*sabspen

 interval2 = max((abs(Tsurf-Tskin(1)) + 0.1d0),10.d0)

 Tskin1 = Tsurf - interval2 !10.
 Tskin2 = Tsurf + interval2 !10.

 Balskin1 = HEATBALANCESKIN(Tskin1,Tskin(1),Tsurf, &
 & RadWater,RadIce,fetch,deltaskin,dt,hskin)
 Balskin2 = HEATBALANCESKIN(Tskin2,Tskin(1),Tsurf, &
 & RadWater,RadIce,fetch,deltaskin,dt,hskin)
 if (Balskin1*Balskin2>0) then
   STOP 'Chorde method is not applicable for skin temperature equation: STOP'
 endif

!print*, Balskin1, Balskin2

 iter=0
 Balskin3=1.d0
 
!CHORDE METHOD TO FIND SURFACE TEMPERATURE

 cyskin: do while (abs(Balskin3)>0.1d0)
   iter=iter+1
   Balskin1 = HEATBALANCESKIN(Tskin1,Tskin(1),Tsurf, &
   & RadWater,RadIce,fetch,deltaskin,dt,hskin)
   Balskin2 = HEATBALANCESKIN(Tskin2,Tskin(1),Tsurf, &
   & RadWater,RadIce,fetch,deltaskin,dt,hskin)
   
   Tskin3=(Tskin1*Balskin2-Tskin2*Balskin1)/(Balskin2-Balskin1)
   
   Balskin3 = HEATBALANCESKIN(Tskin3,Tskin(1),Tsurf, &
   & RadWater,RadIce,fetch,deltaskin,dt,hskin)
   if (iter>maxiter) then
     if (abs(Tskin1-Tskin2)<0.01d0) then
       write (nunit_,*) Balskin3
       exit cyskin
     else
       STOP 'Skin temperature iterations do not converge'
     endif        
   endif
   if     (Balskin1*Balskin3<0.) then
     Tskin2=Tskin3
   elseif (Balskin2*Balskin3<0.) then
     Tskin1=Tskin3
   endif
 enddo cyskin

!print*, 'Skin disbalance', Balskin3, Tskin3, hskin
!read*

 Tskin(2) = Tskin3

 if (firstcall) firstcall=.false.
 RETURN
 END SUBROUTINE TSKIN_SOLVER



 FUNCTION HEATBALANCESKIN(Tskin_iter,Tskin_prev,Tsurf, &
 & RadWater,RadIce,fetch,deltaskin,dt,hskin)
 
!The function HEATBALANCESKIN computes the heat balance
!of the skin at the top of water coloumn

 use PHYS_CONSTANTS, only: &
 & cw, &
 & row0, &
 & lamw0, &
 & sabs
 use ATMOS, only: &
 & hflux, &
 & Elatent, &
 & Radbal   
 use RADIATION, only : rad_type

 implicit none
 
 real(kind=ireals) :: HEATBALANCESKIN

 real(kind=ireals), intent(in) :: Tskin_iter
 real(kind=ireals), intent(in) :: Tskin_prev
 real(kind=ireals), intent(in) :: Tsurf
 !real(kind=ireals), intent(in) :: extwat
 !real(kind=ireals), intent(in) :: extice
 type(rad_type), intent(in) :: RadWater
 type(rad_type), intent(in) :: RadIce
 real(kind=ireals), intent(in) :: fetch
 real(kind=ireals), intent(in) :: deltaskin      
 real(kind=ireals), intent(in) :: dt
 
 real(kind=ireals), intent(out):: hskin

! real(kind=ireals) :: FRAC_SOLAR_ABSORB_SKIN
! real(kind=ireals) :: SKIN_THICKNESS
! real(kind=ireals) :: HEATCONTENTSKIN

 integer(kind=iintegers), save :: surftyp
 data                surftyp /4/
 

 call SENSLATMOM_FLUXES(Tskin_iter,fetch)

!deltaskin = SKIN_THICKNESS()
!    sabs = FRAC_SOLAR_ABSORB_SKIN(deltaskin)

 call RADBALANCE       (Tskin_iter,RadWater,RadIce,surftyp)

 hskin = - lamw0*(Tsurf-Tskin_iter)/deltaskin

 HEATBALANCESKIN = HEATCONTENTSKIN(Tskin_iter,Tskin_prev,deltaskin,dt) - &
 & Radbal + hflux + Elatent + hskin
 
 END FUNCTION HEATBALANCESKIN
 
 
 FUNCTION HEATCONTENTSKIN(Tskin_iter,Tskin_prev,deltaskin,dt)
 
!Function HEATCONTENTSKIN calculates the rate of change of
!skin heat content

 use PHYS_CONSTANTS, only: &
 & cw_m_row0
 
 implicit none
 real(kind=ireals) :: HEATCONTENTSKIN
 
 real(kind=ireals), intent(in) :: Tskin_iter
 real(kind=ireals), intent(in) :: Tskin_prev
 real(kind=ireals), intent(in) :: deltaskin
 real(kind=ireals), intent(in) :: dt
 
 HEATCONTENTSKIN = deltaskin*cw_m_row0*(Tskin_iter-Tskin_prev)/dt
 
 END FUNCTION HEATCONTENTSKIN
 
 
 FUNCTION SKIN_THICKNESS()
 
!Function SKIN_THICKNESS calculates the thickness of cool skin
!at the water surface
 
 use PHYS_CONSTANTS, only : &
 & roa0, &
 & row0, niu_wat
 use ATMOS, only : &
 & velfrict, &
 & velfrict_prev

 implicit none

 real(kind=ireals), parameter :: const_Saunders = 6.d0 ! Dimensionless
 real(kind=ireals), parameter :: SKIN_THICKNESS_min = 0.001d0
 real(kind=ireals), parameter :: SKIN_THICKNESS_max = 0.01d0

 real(kind=ireals) :: SKIN_THICKNESS

! SKIN_THICKNESS = 0.001d0
 if (velfrict.ne.0) then
   SKIN_THICKNESS = const_Saunders*niu_wat/(sqrt(roa0/row0)*velfrict_prev)
 else
   SKIN_THICKNESS = 0.01d0
 endif


 SKIN_THICKNESS = max(SKIN_THICKNESS, SKIN_THICKNESS_min)
 SKIN_THICKNESS = min(SKIN_THICKNESS, SKIN_THICKNESS_max)
 
 END FUNCTION SKIN_THICKNESS



 FUNCTION FRAC_SOLAR_ABSORB_SKIN(skin_thickness)
 
!Function FRAC_SOLAR_ABSORB_SKIN calculates the fraction of solar radiation
!absorbed by cool skin at the water surface according to Fairall et al., 1996
 
 implicit none
 
 real(kind=ireals) :: FRAC_SOLAR_ABSORB_SKIN
 
 real(kind=ireals), parameter :: c1 = 6.5d-2
 real(kind=ireals), parameter :: c2 = 1.1d+1
 real(kind=ireals), parameter :: c3 = 6.6d-5
 real(kind=ireals), parameter :: c4 = 8.0d-4
 
!Input variables      
 real(kind=ireals), intent(in) :: skin_thickness
 
 FRAC_SOLAR_ABSORB_SKIN = c1 + c2*skin_thickness - c3/skin_thickness * &
 & (1 - exp(-skin_thickness/c4))
 
 END FUNCTION FRAC_SOLAR_ABSORB_SKIN

END MODULE SKIN_MOD

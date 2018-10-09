MODULE TRIBUTARIES

use LAKE_DATATYPES, only : ireals, iintegers

use ARRAYS_GRID, only : gridspacing_type
use ARRAYS_WATERSTATE, only : waterstate_type

integer(kind=iintegers), save :: itriblev

contains
!>Subroutine TRIBUTARIES computes the change of the
!!horizontally average lake temperature, salinity and velocities
!!due to tributaries (inflows and outflows) advection
!!
!!INPUT variables:
!!dt --- timestep, sec
!!
!!INTPUT/OUTPUT variables:
!!Tw --- the temperature profile in lake, C
!!Sal --- salinity, kg/kg
!!(u,v) --- velocity components, m/s
SUBROUTINE TRIBTEMP(time,dt,h1,z_full,area_int,gsp,gas,wst, &
 & ueffl,veffl,Tweffl,Saleffl,spinup_done)

 use DRIVING_PARAMS , only : &
 & N_tribin, & 
 & N_triblev, & 
 & N_tribout, &
 & U_tribin, &
 & U_tribout, &
 & T_tribin, &
 & Sal_tribin, &
 & DIC_tribin, &
 & DOC_tribin, &
 & POC_tribin, &
 & meth_tribin, &
 & Ux_tribin, &
 & Uy_tribin, &
 & width_tribin, &
 & width_tribout, &
 & disch_tribin, &
 & disch_tribout, &
 & M, &
 & dttribupdate, &
 & missing_value, &
 & outflpar, &
 & tribheat, &
 & carbon_model

 use ARRAYS, only : gas_type

 implicit none

 real(kind=ireals), intent(in)    :: time, dt, h1
 real(kind=ireals), intent(in)    :: z_full(1:M+1)
 real(kind=ireals), intent(in)    :: area_int(1:M+1)
 real(kind=ireals), intent(in)    :: ueffl(1:M+1), veffl(1:M+1)
 real(kind=ireals), intent(in)    :: Tweffl(1:M+1), Saleffl(1:M+1)

 type(gridspacing_type), intent(in)   :: gsp
 type(gas_type)        , intent(inout):: gas
 type(waterstate_type) , intent(inout):: wst

 logical, intent(in) :: spinup_done

!Local variables
 real(kind=ireals), parameter :: day_sec = 24.*60.*60.

 real(kind=ireals) :: heatinflow, salinflow, uinflow, vinflow, &
 & methinflow, DICinflow, POCinflow, DOCinflow
 real(kind=ireals) :: meanTin, meansalin, meanuin, meanvin, &
 & meanmethin, meanDICin, meanDOCin, meanPOCin
 real(kind=ireals) :: outflow, work
 real(kind=ireals) :: invdt

 real(kind=ireals), allocatable :: U_tribin_(:,:), T_tribin_(:,:), Sal_tribin_(:,:), &
 & Ux_tribin_(:,:), Uy_tribin_(:,:), U_tribout_(:,:), &
 & meth_tribin_(:,:), DIC_tribin_(:,:),DOC_tribin_(:,:),POC_tribin_(:,:), &
 & width_tribin_(:,:), width_tribout_(:,:), z_(:)

 integer(kind=iintegers) :: i, j, k

 logical :: flag
 logical, save :: firstcall = .true.

 if (firstcall) then
   if ((tribheat%par == 2 .or. tribheat%par == 3) &
     & .and. dttribupdate == missing_value) then
     write(*,*) 'Missing tributaries update interval: STOP'
     STOP
   endif
   ! Allocation of inflow and outflow profiles
   if (tribheat%par == 1 .or. tribheat%par == 2) then
     itriblev = N_triblev%par
   elseif (tribheat%par == 3) then
     itriblev = M+1
   endif
   allocate (width_tribin  (1:N_tribin%par ,1:itriblev) )
   allocate (U_tribin      (1:N_tribin%par ,1:itriblev) )
   allocate (T_tribin      (1:N_tribin%par ,1:itriblev) )
   allocate (Sal_tribin    (1:N_tribin%par ,1:itriblev) )
   allocate (meth_tribin   (1:N_tribin%par ,1:itriblev) )
   allocate (DIC_tribin    (1:N_tribin%par ,1:itriblev) )
   allocate (DOC_tribin    (1:N_tribin%par ,1:itriblev) )
   allocate (POC_tribin    (1:N_tribin%par ,1:itriblev) )
   allocate (Ux_tribin     (1:N_tribin%par ,1:itriblev) )
   allocate (Uy_tribin     (1:N_tribin%par ,1:itriblev) )
   allocate (disch_tribin  (1:N_tribin%par ,1:M+1) )
   allocate (disch_tribout (1:N_tribout    ,1:M+1) )
   allocate (width_tribout (1:N_tribout    ,1:itriblev) )
   allocate (U_tribout     (1:N_tribout    ,1:itriblev) ) 
 endif 

 ! If tribheat == 1, tributaries' and effluent's data are read once, and kept constant throughout model integration time
 if (tribheat%par == 1 .and. firstcall) call TRIBUPDATE(spinup_done,gsp,wst)

 !print*, 'time', time, dt
 if (tribheat%par == 2 .or. tribheat%par == 3) then 
   work = mod(time-dt,dttribupdate*day_sec) ! time is for next time step, thus lateral inflow
                                            ! values are updated to the current time step
   !print*, 'work',work,firstcall
   if ((work >= 0.d0 .and. work < dt) .or. firstcall .or. spinup_done) then
     ! Updating tributary data
     call TRIBUPDATE(spinup_done,gsp,wst)
   endif
 endif

 allocate( U_tribin_(1:N_tribin%par,1:M+1), T_tribin_(1:N_tribin%par,1:M+1), &
 & Sal_tribin_(1:N_tribin%par,1:M+1), Ux_tribin_(1:N_tribin%par,1:M+1), &
 & Uy_tribin_(1:N_tribin%par,1:M+1), U_tribout_(1:N_tribout,1:M+1), &
 & meth_tribin_(1:N_tribin%par,1:M+1), DOC_tribin_(1:N_tribin%par,1:M+1), &
 & DIC_tribin_(1:N_tribin%par,1:M+1), POC_tribin_(1:N_tribin%par,1:M+1), &
 & width_tribin_(1:N_tribin%par,1:M+1), width_tribout_(1:N_tribout,1:M+1))

 !Interpolation of tributaries' data on the model grid
 if (itriblev /= M+1) then 
   allocate(z_(1:itriblev))
   forall (i=1:itriblev) z_(i) = (i-1)*h1/real(itriblev-1)
   do k = 1, N_tribin%par
     !call LININTERPOL(z_,U_tribin(k,1:itriblev),itriblev,z_full,U_tribin_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,T_tribin(k,1:itriblev),itriblev,z_full,T_tribin_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,Sal_tribin(k,1:itriblev),itriblev,z_full,Sal_tribin_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,Ux_tribin(k,1:itriblev),itriblev,z_full,Ux_tribin_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,Uy_tribin(k,1:itriblev),itriblev,z_full,Uy_tribin_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,U_tribout(k,1:itriblev),itriblev,z_full,U_tribout_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,width_tribin(k,1:itriblev),itriblev,z_full,width_tribin_(k,1:M+1),M+1,flag)
     !call LININTERPOL(z_,width_tribout(k,1:itriblev),itriblev,z_full,width_tribout_(k,1:M+1),M+1,flag)
     call PIECEWISECONSTINT(z_,U_tribin    (k,1:itriblev),itriblev,z_full,U_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,T_tribin    (k,1:itriblev),itriblev,z_full,T_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,Sal_tribin  (k,1:itriblev),itriblev,z_full,Sal_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,meth_tribin (k,1:itriblev),itriblev,z_full,meth_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,DIC_tribin  (k,1:itriblev),itriblev,z_full,DIC_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,DOC_tribin  (k,1:itriblev),itriblev,z_full,DOC_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,POC_tribin  (k,1:itriblev),itriblev,z_full,POC_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,Ux_tribin   (k,1:itriblev),itriblev,z_full,Ux_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,Uy_tribin   (k,1:itriblev),itriblev,z_full,Uy_tribin_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,width_tribin(k,1:itriblev),itriblev,z_full,width_tribin_(k,1:M+1),M+1,2_iintegers)
   enddo
   do k = 1, N_tribout
     call PIECEWISECONSTINT(z_,U_tribout(k,1:itriblev),itriblev,z_full,U_tribout_(k,1:M+1),M+1,2_iintegers)
     call PIECEWISECONSTINT(z_,width_tribout(k,1:itriblev),itriblev,z_full,width_tribout_(k,1:M+1),M+1,2_iintegers)
     !print*, U_tribout, U_tribout_
   enddo 
   deallocate(z_)
 else
   U_tribin_ = U_tribin
   T_tribin_ = T_tribin
   Sal_tribin_ = Sal_tribin
   meth_tribin_ = meth_tribin
   DOC_tribin_ = DOC_tribin
   DIC_tribin_ = DIC_tribin
   POC_tribin_ = POC_tribin
   Ux_tribin_ = Ux_tribin
   Uy_tribin_ = Uy_tribin
   U_tribout_ = U_tribout
   width_tribin_ = width_tribin
   width_tribout_ = width_tribout
 endif

 ! Calculating discharges of all tributaries and an effluent
 forall (i=1:M+1, j=1:N_tribin%par) disch_tribin(j,i) = &
 & width_tribin_(j,i)*U_tribin_(j,i)*gsp%ddz05(i-1)*h1
 !disch_tribout = 0. ! Assuming N_tribout == 1
 forall (i=1:M+1, j=1:N_tribout)   disch_tribout(j,i) = &
 & width_tribout_(j,i)*U_tribout_(j,i)*gsp%ddz05(i-1)*h1

 !! Diagnostic calculation of total heat influx by tributaries
 !work = 0.
 !do k = 1, N_tribin%par
 !  do i = 1, M+1
 !    work = work + U_tribin_(k,i)*T_tribin_(k,i)*width_tribin_(k,i)*gsp%ddz05(i-1)*h1
 !  enddo
 !enddo
 !print*, 'Total heat influx ', work
 !! Diagnostic calculation of total ouflow discharge
 !work = 0.
 !do k = 1, N_tribout
 !  do i = 1, M+1
 !    work = work + U_tribout_(k,i)*width_tribout_(k,i)*gsp%ddz05(i-1)*h1
 !  enddo
 !enddo
 !print*, 'Total outflow discharge ', work
 !stop

 invdt = 1./dt

 select case(outflpar%par)

 case(2)
!  Values at outflow are calculated before following "averaged Lagrangian" approach
!  Explicit scheme
   do i = 2, M
     wst%Tw1(i) = wst%Tw1(i) + dt*(U_tribin_(1,i)*T_tribin_(1,i)*width_tribin_(1,i) - &
     & U_tribout_(1,i)*Tweffl(i)*width_tribout_(1,i))/area_int(i)
     wst%Sal1(i) = wst%Sal1(i) + dt*(U_tribin_(1,i)*Sal_tribin_(1,i)*width_tribin_(1,i) - &
     ! Missing DIC, POC, DOC and ch4 inflows/outflows
     & U_tribout_(1,i)*Saleffl(i)*width_tribout_(1,i))/area_int(i)
     wst%u1(i) = wst%u1(i) + dt*(U_tribin_(1,i)*Ux_tribin_(1,i)*width_tribin_(1,i) - &
     & U_tribout_(1,i)*ueffl(i)*width_tribout_(1,i))/area_int(i)
     wst%v1(i) = wst%v1(i) + dt*(U_tribin_(1,i)*Uy_tribin_(1,i)*width_tribin_(1,i) - &
     & U_tribout_(1,i)*veffl(i)*width_tribout_(1,i))/area_int(i)
   enddo
 
 case default ! outflpar%par == 0 or 1
 
 do i = 1, M+1 !1, M+1

   heatinflow = 0.
   salinflow = 0.
   DICinflow = 0.
   DOCinflow = 0.
   POCinflow = 0.
   methinflow = 0.
   uinflow = 0.
   vinflow = 0.
   do k = 1, N_tribin%par
     heatinflow = heatinflow + U_tribin_(k,i)*T_tribin_(k,i)*width_tribin_(k,i)
     salinflow = salinflow + U_tribin_(k,i)*Sal_tribin_(k,i)*width_tribin_(k,i)
     methinflow = methinflow + U_tribin_(k,i)*meth_tribin_(k,i)*width_tribin_(k,i)
     DICinflow = DICinflow   + U_tribin_(k,i)*DIC_tribin_(k,i) *width_tribin_(k,i)
     DOCinflow = DOCinflow   + U_tribin_(k,i)*DOC_tribin_(k,i) *width_tribin_(k,i)
     POCinflow = POCinflow   + U_tribin_(k,i)*POC_tribin_(k,i) *width_tribin_(k,i)
     uinflow = uinflow + U_tribin_(k,i)*Ux_tribin_(k,i)*width_tribin_(k,i)
     vinflow = vinflow + U_tribin_(k,i)*Uy_tribin_(k,i)*width_tribin_(k,i)
   enddo
   heatinflow = heatinflow/area_int(i)
   salinflow = salinflow/area_int(i)
   methinflow = methinflow/area_int(i)
   DICinflow = DICinflow/area_int(i)
   DOCinflow = DOCinflow/area_int(i)
   POCinflow = POCinflow/area_int(i)
   uinflow = uinflow/area_int(i)
   vinflow = vinflow/area_int(i)

   outflow = 0.
   do k = 1, N_tribout
     outflow = outflow + U_tribout_(k,i)*width_tribout_(k,i)
   enddo
   outflow = outflow/area_int(i)
   !print*, 'outflow', outflow, U_tribin_(1,i),T_tribin_(1,i),width_tribin_(1,i), area_int(i), invdt

   meanTin = 0.
   meansalin = 0.
   meanmethin = 0.
   meanDICin = 0.
   meanDOCin = 0.
   meanPOCin = 0.
   meanuin = 0.
   meanvin = 0.
   if (outflpar%par == 1) then
     do k = 1, N_tribin%par
       meanTin = meanTin + T_tribin_(k,i)*width_tribin_(k,i)
       meansalin = meansalin + Sal_tribin_(k,i)*width_tribin_(k,i)
       meanmethin = meanmethin + meth_tribin_(k,i)*width_tribin_(k,i)
       meanDICin = meanDICin + DIC_tribin_(k,i)*width_tribin_(k,i)
       meanDOCin = meanDOCin + DOC_tribin_(k,i)*width_tribin_(k,i)
       meanPOCin = meanPOCin + POC_tribin_(k,i)*width_tribin_(k,i)
       meanuin = meanuin + Ux_tribin_(k,i)*width_tribin_(k,i)
       meanvin = meanvin + Uy_tribin_(k,i)*width_tribin_(k,i)
     enddo
     work = sum(width_tribin(1:N_tribin%par,i))
     meanTin = meanTin/work
     meansalin = meansalin/work
     meanmethin = meanmethin/work
     meanDICin = meanDICin/work
     meanDOCin = meanDOCin/work
     meanPOCin = meanPOCin/work
     meanuin = meanuin/work
     meanvin = meanvin/work
   endif

   !Implicit scheme for the outflow term
   !Check: if for any of tributaries the data on a variable is missing,
   !this variable is not updated with inlet/outlet sources
   work = invdt + outflow*(1. + outflpar%par)

   !print*, 'Tw1', Tw
   !print*, heatinflow, outflow, meanTin
   if (.not. any(T_tribin_(:,i) == missing_value)) &
   & wst%Tw1(i) = (invdt*wst%Tw1(i) + heatinflow + outflpar%par*outflow*meanTin)/work
   !print*, 'Tw2', Tw
   !print*, heatinflow, outflow, meanTin
   !read*

   if (.not. any(Sal_tribin_(:,i) == missing_value)) &
   & wst%Sal1(i) = (invdt*wst%Sal1(i) + salinflow + outflpar%par*outflow*meansalin)/work

   if (.not. any(meth_tribin_(:,i) == missing_value)) &
   & gas%qwater(i,1) = (invdt*gas%qwater(i,1) + methinflow + outflpar%par*outflow*meanmethin)/work

   if (carbon_model%par == 2) then
     ! Inlet DOC is added to the allochtonous DOC pool
     if (.not. any(DOC_tribin_(:,i) == missing_value)) &
     & gas%DOC(i,2) = (invdt*gas%DOC(i,2) + DOCinflow + outflpar%par*outflow*meanDOCin)/work
       gas%DOC(i,1) = (invdt*gas%DOC(i,1)                                             )/work
     if (.not. any(POC_tribin_(:,i) == missing_value)) &
     & gas%POCL(i) = (invdt*gas%POCL(i) + POCinflow + outflpar%par*outflow*meanPOCin)/work
   endif

   if (.not. any(DIC_tribin_(:,i) == missing_value)) &
   & gas%DIC(i,1) = (invdt*gas%DIC(i,1) + DICinflow + outflpar%par*outflow*meanDICin)/work

   !Sal(i) = (invdt*Sal(i) + salinflow + outflpar%par*outflow*meansalin)/work
   if (.not. any(Ux_tribin_(:,i) == missing_value)) &
   & wst%u1(i) = (invdt*wst%u1(i) + uinflow + outflpar%par*outflow*meanuin)/work

   if (.not. any(Uy_tribin_(:,i) == missing_value)) &
   & wst%v1(i) = (invdt*wst%v1(i) + vinflow + outflpar%par*outflow*meanvin)/work
   ! qwater(i) = (invdt*qwater(i) + 0.d0)/(invdt + outflow) ! Implying zero concentration of methane in the inflow

 enddo

 end select 

!print*, 'I am in TRIBUTARIES!' ! Debug string

 deallocate(U_tribin_, T_tribin_, Sal_tribin_, Ux_tribin_, Uy_tribin_, U_tribout_, &
 & width_tribin_, width_tribout_, DIC_tribin_, DOC_tribin_, &
 & POC_tribin_, meth_tribin_)

 if (firstcall) firstcall = .false.
 RETURN
 END SUBROUTINE TRIBTEMP


!> Subroutine TRIBUPDATE updates tributary data
SUBROUTINE TRIBUPDATE(spinup_done,gsp,wst)


use INOUT_PARAMETERS, only : &
& lake_misc_unit_min, &
& lake_misc_unit_max

use DRIVING_PARAMS, only : &
& N_tribin, & 
& N_triblev, & 
& N_tribout, &
& U_tribin, &
& U_tribout, &
& T_tribin, &
& Sal_tribin, &
& DIC_tribin, &
& DOC_tribin, &
& POC_tribin, &
& meth_tribin, &
& Ux_tribin, &
& Uy_tribin, &
& width_tribin, &
& width_tribout, &
& disch_tribin, &
& disch_tribout, &
& M, &
& eos,lindens, &
& missing_value, &
& tribheat, &
& fileinflow, &
& fileoutflow, &
& READPROFILE

use INOUT, only : CHECK_UNIT
use PHYS_CONSTANTS, only : pref
use WATER_DENSITY, only : DENSITY_W

implicit none

! Input variables

logical, intent(in) :: spinup_done

type(gridspacing_type), intent(in) :: gsp
type(waterstate_type) , intent(in) :: wst

! Local variables
integer(kind=iintegers), parameter :: mixtype = 0 !> Indicates the type of mixing of tributaries water at the
                                                  !! entrance to the lake, relevant it tribheat$par == 3

real(kind=ireals), allocatable :: work(:,:)
real(kind=ireals) :: disch, depth, width, temp, sal, &
& dirprojX, dirprojY, DOC, POC, DIC, CH_4, xx, rowtrib
integer(kind=iintegers), save :: nunit1 = lake_misc_unit_min, &
& nunit2 = lake_misc_unit_min + 1
integer(kind=iintegers) :: ncol, i, j, idep

character(len=20) :: workchar

logical, save :: firstcall = .true.
!print*, 'entered', tribheat%par
if (firstcall) then
  call CHECK_UNIT(lake_misc_unit_min,lake_misc_unit_max,nunit1)
  call CHECK_UNIT(lake_misc_unit_min,lake_misc_unit_max,nunit2)
  open (unit=nunit1,file='./data/'//fileinflow(1:len_trim(fileinflow)),status='old')
  open (unit=nunit2,file='./data/'//fileoutflow(1:len_trim(fileoutflow)),status='old')
endif

if (spinup_done) then
  rewind nunit1
  rewind nunit2
endif

if_tribform : if (tribheat%par == 2) then
  !Reading inflow profiles from file
  ncol = 11
  allocate (work(itriblev,ncol))
  do i = 1, N_tribin%par
    call READPROFILE(nunit1,itriblev,ncol,work) 
    width_tribin(i,:) = work (:,2)
    U_tribin    (i,:) = work (:,3)
    T_tribin    (i,:) = work (:,4)
    Sal_tribin  (i,:) = work (:,5)
    Ux_tribin   (i,:) = work (:,6)
    Uy_tribin   (i,:) = work (:,7)
    DOC_tribin  (i,:) = work (:,8)
    POC_tribin  (i,:) = work (:,9)
    DIC_tribin  (i,:) = work (:,10)
    meth_tribin (i,:) = work (:,11)
  enddo
  deallocate(work)

  ncol = 3
  allocate (work(itriblev,ncol))
  call READPROFILE(nunit2,itriblev,ncol,work) 
  ! In current version there is only ONE outflow 
  width_tribout(1,:) = work (:,2)
  U_tribout    (1,:) = work (:,3)
  deallocate(work)
else if (tribheat%par == 3) then
  ! Reading bulk characteristics of inlets and an outlet.
  ! The data for each inlet in an input file should be given in rows containing:
  ! timestring, discharge, depth, width, temperature, salinity, riverflow direction
  ! projection on X axis (-1,+1), riverflow direction projection on Y axis (-1,+1),
  ! DOC, POC, DIC, CH_4
  !print*, 'entered2'
  do i = 1, N_tribin%par

    read(nunit1,*) xx, disch, depth, width, temp, sal, &
    & dirprojX, dirprojY, DOC, POC, DIC, CH_4
    !print*, xx, disch, depth, width, temp, sal, &
    !& dirprojX, dirprojY, DOC, POC, DIC, CH_4
    !disch = 0.
    if (mixtype == 1 .and. temp /= missing_value) then
      ! Assuming that inlet properties are advected at levels from surface to the lake depth 
      ! where the water density equals to that of an inlet
      allocate(work(1:M+1,1))
      xx = sal
      if (sal == missing_value) xx = 0.
      call DENSITY_W(0_iintegers,eos%par,lindens%par,(/temp/),(/xx/),(/pref/),work(1,1))
      rowtrib = work(1,1)
      work(1:M+1,1) = abs(wst%row(1:M+1) - rowtrib)
      j = minloc(work(1:M+1,1),1)
      depth = gsp%z_half(j)
      !print*, work
      !print*, 'riv depth', j
      !read*
      deallocate(work)
    endif

    idep = 1
    do while (gsp%z_half(idep)/depth < 1.) 
      idep = idep + 1
    enddo
    if (idep /= 1_iintegers) idep = idep - 1
    depth = gsp%z_half(idep)

    width_tribin(i,1:idep) = width                         ; width_tribin(i,idep+1:M+1) = 0.
    U_tribin    (i,1:idep) = disch/(depth*width)           ; U_tribin    (i,idep+1:M+1) = 0.
    T_tribin    (i,1:idep) = temp                          ; T_tribin    (i,idep+1:M+1) = 0.
    Sal_tribin  (i,1:idep) = sal                           ; Sal_tribin  (i,idep+1:M+1) = 0.
    Ux_tribin   (i,1:idep) = disch/(depth*width)*dirprojX  ; Ux_tribin   (i,idep+1:M+1) = 0.
    Uy_tribin   (i,1:idep) = disch/(depth*width)*dirprojY  ; Uy_tribin   (i,idep+1:M+1) = 0.
    DOC_tribin  (i,1:idep) = DOC                           ; DOC_tribin  (i,idep+1:M+1) = 0.
    POC_tribin  (i,1:idep) = POC                           ; POC_tribin  (i,idep+1:M+1) = 0.
    DIC_tribin  (i,1:idep) = DIC                           ; DIC_tribin  (i,idep+1:M+1) = 0.
    meth_tribin (i,1:idep) = CH_4                          ; meth_tribin (i,idep+1:M+1) = 0.
    !if (idep > 10) then
    !print*, idep
    !read*
    !endif
  enddo
  !print*, 'in', width_tribin, U_tribin, T_tribin
  !read*

  ! In current version there is only ONE outflow 
  read(nunit2,*) workchar, disch, depth, width
  !disch = 0.
  idep = 1
  do while (gsp%z_half(idep)/depth < 1.) 
    idep = idep + 1
  enddo
  if (idep /= 1_iintegers) idep = idep - 1
  depth = gsp%z_half(idep)

  width_tribout(1,1:idep) = width               ; width_tribout(1,idep+1:M+1) = 0.
  U_tribout    (1,1:idep) = disch/(depth*width) ; U_tribout    (1,idep+1:M+1) = 0.
endif if_tribform

if (firstcall) firstcall = .false.
END SUBROUTINE TRIBUPDATE 

!> Calculates artificial abstraction rate (m**3/s) 
!! for the new year dependent
!! on reseroir volume (m**3) at the end of previous year.
FUNCTION ABSTR(vol)

implicit none

! Input variables
real(kind=ireals), intent(in) :: vol

! Output variables
real(kind=ireals) :: ABSTR

!Local variables
real(kind=ireals), parameter :: coef0 = 0.2770399597
real(kind=ireals), parameter :: coef1 = -0.0344416893
real(kind=ireals), parameter :: coef2 = 0.0153926477

ABSTR = coef2*vol*vol + coef1*vol + coef0

END FUNCTION ABSTR

!> Subroutine TENDCALC calculates tendencies of variables
SUBROUTINE TENDCALC(dt, utend, vtend, Twtend, Saltend)

use ARRAYS, only : u1, u2, v1, v2
use ARRAYS_GRID, only : ddz, dzeta_int
use ARRAYS_WATERSTATE, only : Tw1,Tw2,Sal1,Sal2
use DRIVING_PARAMS, only : M
use ARRAYS_BATHYM, only : dhw, h1


implicit none

! Input variables
real(kind=ireals), intent(in) :: dt ! Timestep, s

! Output variables
real(kind=ireals), intent(out) :: utend(1:M+1) ! Tendency of x-component of speed 
real(kind=ireals), intent(out) :: vtend(1:M+1) ! Tendency of y-component of speed
real(kind=ireals), intent(out) :: Twtend(1:M+1) ! Tendency of temperature
real(kind=ireals), intent(out) :: Saltend(1:M+1) ! Tendency of salinity

! Local variables
real(kind=ireals) :: xx
integer(kind=iintegers) :: i ! Loop index

! Tendencies are calculated in dzeta-coordinates
do i = 2, M
  xx = dhw*dzeta_int(i)/(dt*h1*(ddz(i-1) + ddz(i)))
  utend(i)    = (u2(i)   - u1(i))/dt - xx*(u2(i+1) - u2(i-1))
  vtend(i)    = (v2(i)   - v1(i))/dt - xx*(v2(i+1) - v2(i-1))
  Twtend(i)   = (Tw2(i)  - Tw1(i))/dt - xx*(Tw2(i+1) - Tw2(i-1))
  Saltend(i)  = (Sal2(i) - Sal1(i))/dt - xx*(Sal2(i+1) - Sal2(i-1))
enddo
utend(1)    = (u2(1)   - u1(1))/dt
vtend(1)    = (v2(1)   - v1(1))/dt
Twtend(1)   = (Tw2(1)  - Tw1(1))/dt
Saltend(1)  = (Sal2(1) - Sal1(1))/dt
xx = dhw/(dt*h1*ddz(M))
utend(M+1)    = (u2(M+1)   - u1(M+1))/dt - xx*(u2(M+1) - u2(M))
vtend(M+1)    = (v2(M+1)   - v1(M+1))/dt - xx*(v2(M+1) - v2(M))
Twtend(M+1)   = (Tw2(M+1)  - Tw1(M+1))/dt - xx*(Tw2(M+1) - Tw2(M))
Saltend(M+1)  = (Sal2(M+1) - Sal1(M+1))/dt - xx*(Sal2(M+1) - Sal2(M))

END SUBROUTINE TENDCALC

!> Subroutine LAGREFFL calculates model variables at the outflow based on
!! inflow values and the tendency of variables along Largangian trajectory.
!! Assumptions:
!! 1. Open boundaries ("water volume, surrounded by water")
!! 2. Currently assumes square shape of cross section!
SUBROUTINE LAGREFFL(utrib,vtrib,Twtrib,Saltrib, &
& utend,vtend,Twtend,Saltend, &
& u,v,area_int,dt,M, &
& ueffl,veffl,Tweffl,Saleffl)


implicit none

! Input variables
real(kind=ireals), intent(in) :: utrib(1:M+1), vtrib(1:M+1) ! Velocity components at inflow
real(kind=ireals), intent(in) :: Twtrib(1:M+1), Saltrib(1:M+1) ! Scalars at inflow
real(kind=ireals), intent(in) :: utend(1:M+1), vtend(1:M+1) ! Tendencies of velocity components, m/s**2
real(kind=ireals), intent(in) :: Twtend(1:M+1), Saltend(1:M+1) ! Tendencies of scalars
real(kind=ireals), intent(in) :: u(1:M+1), v(1:M+1) ! Velocity components, m/s
real(kind=ireals), intent(in) :: area_int(1:M+1) ! Area of horizontal cross-section of water body, m**2

real(kind=ireals), intent(in) :: dt ! Timestep, s

integer(kind=iintegers), intent(in) :: M ! Number od water layers

! Input/output variables
real(kind=ireals), intent(inout) :: ueffl(1:M+1), veffl(1:M+1) ! Velocity components at outflow
real(kind=ireals), intent(inout) :: Tweffl(1:M+1), Saleffl(1:M+1) ! Scalars at outflow

! Local variables
real(kind=ireals) :: L, tlagr
real(kind=ireals) :: xu, xv, xTw, xSal
integer(kind=iintegers) :: i ! Loop index

do i = 2, M
  L = sqrt(area_int(i)) ! An estimate of mean Lagrangian trajectory path length
  tlagr = L/sqrt(u(i)**2 + v(i)**2)
  xu   = utrib(i)   + utend(i)*tlagr
  xv   = vtrib(i)   + vtend(i)*tlagr
  xTw  = Twtrib(i)  + Twtend(i)*tlagr
  xSal = Saltrib(i) + Saltend(i)*tlagr
  if (tlagr > dt) then
!   Interpolation in time
    ueffl(i)   = ueffl(i)  *(1. - dt/tlagr) + xu*dt/tlagr
    veffl(i)   = veffl(i)  *(1. - dt/tlagr) + xv*dt/tlagr
    Tweffl(i)  = Tweffl(i) *(1. - dt/tlagr) + xTw*dt/tlagr
    Saleffl(i) = Saleffl(i)*(1. - dt/tlagr) + xSal*dt/tlagr
  else
!   Extrapolation in time
    ueffl(i)   = ueffl(i)   + dt/tlagr*(xu   - ueffl(i))
    veffl(i)   = veffl(i)   + dt/tlagr*(xv   - veffl(i))
    Tweffl(i)  = Tweffl(i)  + dt/tlagr*(xTw  - Tweffl(i))
    Saleffl(i) = Saleffl(i) + dt/tlagr*(xSal - Saleffl(i))
  endif
enddo

END SUBROUTINE LAGREFFL

END MODULE TRIBUTARIES

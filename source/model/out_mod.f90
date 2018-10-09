MODULE OUT_MOD

use LAKE_DATATYPES, only : ireals, iintegers
use INOUT, only : CHECK_UNIT
use DRIVING_PARAMS, only : missing_value

character, save :: outpath*60 

contains

SUBROUTINE MON_OUT(ix,iy,nx,ny,year,month,day,hour,time,z_out,nout)

use ARRAYS_WATERSTATE, only : Tw1, Sal1
use ARRAYS_TURB, only : E1, E2, eps1, S, Gen, Gen_seiches, row, TKE_turb_trans, TF, KT
use ARRAYS, only : u1, v1
use ARRAYS_GRID, only : dzeta_int, z_full, z_half
use ARRAYS_BATHYM, only : h1, l1, voldef, vol
use DRIVING_PARAMS, only : M
use PHYS_CONSTANTS, only : cw_m_row0

use INOUT_PARAMETERS, only : &
& lake_mon_out_unit_min, &
& lake_mon_out_unit_max

implicit none

! Input variables
! Reals
real(kind=ireals), intent(in) :: hour
real(kind=ireals), intent(in) :: time
real(kind=ireals), intent(in) :: z_out(1:nout) ! The output grid

! Integers
integer(kind=iintegers), intent(in) :: year, month, day
integer(kind=iintegers), intent(in) :: ix, iy
integer(kind=iintegers), intent(in) :: nx, ny
integer(kind=iintegers), intent(in) :: nout

! Local variables
! Reals
real(kind=ireals), allocatable :: accum_var (:,:,:,:)
real(kind=ireals), allocatable :: var       (:,:,:,:)
real(kind=ireals), allocatable :: accum_var_scalar (:,:,:)
real(kind=ireals), allocatable :: var_scalar       (:,:,:)
real(kind=ireals), allocatable :: tsteps    (:,:)
real(kind=ireals), allocatable :: Profile_out(:,:)
real(kind=ireals) :: voldef0, voldef0y
real(kind=ireals) :: work1, work2
real(kind=ireals), allocatable :: valmax(:)

real(kind=ireals), external:: DZETA

! Integers
integer(kind=iintegers), parameter :: n_var = 11, n_var_scalar = 5, n_var_max = 3
integer(kind=iintegers), parameter :: n_var1 = 4 ! Number of vector variables defined at cell interfaces
integer(kind=iintegers), allocatable :: month_old(:,:)
integer(kind=iintegers), allocatable :: out_unita(:,:)
integer(kind=iintegers) :: out_unit = lake_mon_out_unit_min
integer(kind=iintegers) :: i, j ! Loop indices

! Characters
character :: month1*2
character :: year1*4
character :: day1*2
character :: hour1*2
character :: timestring*6
character :: coords_point*6
character :: formatline*50

! Logicals
logical :: firstcall, flag
logical, allocatable :: firstcallixiy(:,:)

data firstcall /.true./

SAVE

if (firstcall) then
  allocate (firstcallixiy(1:nx,1:ny))
  allocate (out_unita(1:nx,1:ny))
  allocate (month_old(1:nx, 1:ny))
  allocate (tsteps    (1:nx, 1:ny))  
  allocate (var      (1:n_var, 1:nx, 1:ny, 1:M+1) )
  allocate (accum_var(1:n_var, 1:nx, 1:ny, 1:M+1) )
  allocate (var_scalar      (1:n_var_scalar, 1:nx, 1:ny) )
  allocate (accum_var_scalar(1:n_var_scalar, 1:nx, 1:ny) )
  allocate (valmax(1:n_var_max))
  firstcallixiy(:,:) = .true.
  out_unita(:,:) = lake_mon_out_unit_min
  month_old(:,:) = month
  tsteps   (:,:) = 0.d0  
  accum_var(:,:,:,:) = 0.d0
  accum_var_scalar(:,:,:) = 0.d0
  voldef0 = 0.
  voldef0y = 0.
  valmax(:) = 0.
endif

if (firstcallixiy(ix,iy)) then
  call CHECK_UNIT(lake_mon_out_unit_min,lake_mon_out_unit_max,out_unita(ix,iy))
  write (coords_point,'(2i3)') ix, iy
  open(out_unita(ix,iy),file=outpath(1:len_trim(outpath))//'monthly/'// &
  & 'depvol'//coords_point//'.dat', status='unknown')
  write (out_unita(ix,iy),*) '1 - year' 
  write (out_unita(ix,iy),*) '2 - month'
  write (out_unita(ix,iy),*) '3 - depth, m'
  write (out_unita(ix,iy),*) '4 - ice thickness, m' 
  write (out_unita(ix,iy),*) '5 - volume, m**3'
  write (out_unita(ix,iy),*) '6 - volume deficit, m**3'
  write (out_unita(ix,iy),*) '7 - accumulated volume deficit, m**3'
  write (out_unita(ix,iy),*) '8 - maximal monthly depth, m'
  write (out_unita(ix,iy),*) '9 - maximal monthly ice thickness, m'
  write (out_unita(ix,iy),*) '10 - maximal monthly volume, m**3'
endif

! Variables located at cell interfaces
var(1,ix,iy,1:M+1) = Tw1 (1:M+1)
var(2,ix,iy,1:M+1) = Sal1(1:M+1)
var(3,ix,iy,1:M+1) = row (1:M+1)
var(4,ix,iy,1:M+1) = sqrt(u1 (1:M+1)*u1 (1:M+1) + v1 (1:M+1)*v1 (1:M+1) )
! Variables located at cell centers
var(5,ix,iy,1:M)   = E1  (1:M)
var(6,ix,iy,1:M)   = TF(1:M)*E2(1:M) !eps1(1:M)
var(7,ix,iy,1:M)   = S   (1:M)
var(8,ix,iy,1:M)   = Gen (1:M)
var(9,ix,iy,1:M)   = Gen_seiches (1:M)
var(10,ix,iy,1:M)  = TKE_turb_trans (1:M)
var(11,ix,iy,1:M)  = KT (1:M)/cw_m_row0 !Kinematic heat conductance, m**2/s

var_scalar(1,ix,iy) = h1
var_scalar(2,ix,iy) = l1
var_scalar(3,ix,iy) = vol
var_scalar(4,ix,iy) = voldef - voldef0
var_scalar(5,ix,iy) = voldef - voldef0
voldef0 = voldef

valmax(1) = max(h1,valmax(1)) ! Calculating monthly maximum
valmax(2) = max(l1,valmax(2)) 
valmax(3) = max(vol,valmax(3))

accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) + var(:,ix,iy,:)
accum_var_scalar(:,ix,iy) = accum_var_scalar(:,ix,iy) + &
& var_scalar(:,ix,iy)
tsteps(ix,iy) = tsteps(ix,iy) + 1._ireals

if (month_old(ix,iy) /= month) then
  !print*, time/(365./12.*24.*3600.)
  !read*
  accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:)/tsteps(ix,iy)
  accum_var_scalar(1:3,ix,iy) = accum_var_scalar(1:3,ix,iy)/tsteps(ix,iy) ! excepting volume deficit

  call DATEMINUS(1,year,month,day,hour,year1,month1,day1,hour1)
  call TIMESTR(6,year1,month1,day1,hour1,timestring)

  call CHECK_UNIT(lake_mon_out_unit_min,lake_mon_out_unit_max,out_unit)
  write (coords_point,'(2i3)') ix, iy
  open(out_unit,file=outpath(1:len_trim(outpath))//'monthly/'// &
  & 'Profiles'//coords_point//timestring//'.dat', status='unknown')
  write (out_unit,*) '1 - depth, m' 
  write (out_unit,*) '2 - temperature, C'
  write (out_unit,*) '3 - salinity, kg/kg' 
  write (out_unit,*) '4 - density, kg/m**3'
  write (out_unit,*) '5 - absolute value of velocity, m/s'
  write (out_unit,*) '6 - turbulent kinetic energy, m**2/s**2'
  write (out_unit,*) '7 - disspation rate, m**2/s**3'
  write (out_unit,*) '8 - TKE production by buoyancy, m**2/s**3'
  write (out_unit,*) '9 - TKE production by shear   , m**2/s**3'
  write (out_unit,*) '10 - TKE production by seiches (Goudsmit paramaterization)  , m**2/s**3'
  write (out_unit,*) '11 - TKE turbulent transport  , m**2/s**3'
  write (out_unit,*) '12 - kinematic heat conductance, m**2/s'
  formatline = '(f12.6,f9.3,e12.4,f9.3,8e12.4)'
  if (nout > 0) then
    allocate (Profile_out(1:nout,1:n_var))
    ! Interpolating to given output levels
    do j = 1, n_var1
      call LININTERPOL (z_full,accum_var(j,ix,iy,:),M+1,z_out,Profile_out(:,j),nout,flag)
    enddo
    do j = n_var1+1, n_var
      call LININTERPOL (z_half,accum_var(j,ix,iy,:),M  ,z_out,Profile_out(:,j),nout,flag)
    enddo
    j = nout
    if (z_out(nout) > z_full(M+1)) then
      do while (z_out(j) > z_full(M+1))
        j = j - 1
      enddo
    endif
    ! Writing to file
    do i = 1, j 
      write (unit = out_unit, fmt = formatline) &
      & -z_out(i), Profile_out(i,1:n_var)
    enddo
  else
    allocate (Profile_out(1:M+1,1:n_var))
    forall (i = 1:M+1, j = 1:n_var1) Profile_out(i,j) = accum_var(j,ix,iy,i)
    ! Interpolating from cell centers to cell interfaces
    do j = n_var1+1, n_var
      call LININTERPOL (z_half,accum_var(j,ix,iy,:),M  ,z_full,Profile_out(:,j),M+1,flag)
    enddo
    ! Writing to file
    do i = 1, M+1 
      write (unit = out_unit, fmt = formatline) &
      & -dzeta_int(i)*h1, Profile_out(i,1:n_var)
    enddo
  endif
  deallocate(Profile_out)
  close(out_unit)


  read(year1,*) work1
  read(month1,*) work2
  write (out_unita(ix,iy),fmt=*) & !, fmt = formatline) &
  & int(work1), int(work2), accum_var_scalar(1:n_var_scalar,ix,iy), &
  & valmax(:)

  
  tsteps(ix,iy) = 0.d0
  accum_var(:,ix,iy,:) = 0.d0
  accum_var_scalar(1:4,ix,iy) = 0.d0
  if (month_old(ix,iy) == 12) accum_var_scalar(5,ix,iy) = 0.d0
  valmax(:) = 0.d0

  month_old(ix,iy) = month
endif

if (firstcall) firstcall = .false.
if (firstcallixiy(ix,iy)) firstcallixiy(ix,iy) = .false.
END SUBROUTINE MON_OUT


SUBROUTINE DAY_OUT(ix,iy,nx,ny,year,month,day,hour)

use ARRAYS_GRID, only : dzeta_int
use ARRAYS_WATERSTATE, only : Tw1, Sal1
use ARRAYS_TURB, only : E1, eps1, KT, Ri
use ARRAYS_BATHYM, only : h1
use ARRAYS, only : u1, v1
use DRIVING_PARAMS, only : M
use PHYS_CONSTANTS, only : cw_m_row0

use INOUT_PARAMETERS, only : &
& lake_day_out_unit_min, &
& lake_day_out_unit_max

implicit none

! Input variables
! Reals
real(kind=ireals), intent(in) :: hour

! Integers
integer(kind=iintegers), intent(in) :: year, month, day
integer(kind=iintegers), intent(in) :: ix, iy
integer(kind=iintegers), intent(in) :: nx, ny

! Local variables
! Reals
real(kind=ireals), allocatable :: tsteps   (:,:)
real(kind=ireals), allocatable :: accum_var(:,:,:,:)
real(kind=ireals), allocatable :: var      (:,:,:,:)

real(kind=ireals), external :: DZETA

! Integers
integer(kind=iintegers), allocatable :: day_old(:,:)
integer(kind=iintegers) :: i ! Loop index
integer(kind=iintegers) :: out_unit = lake_day_out_unit_min
integer(kind=iintegers) :: nvars = 7

character :: month1*2
character :: year1*4
character :: day1*2
character :: hour1*2
character :: timestring*8
character :: coords_point*6

logical :: firstcall

data firstcall /.true./
      
SAVE

if (firstcall) then
  allocate (day_old   (1:nx, 1:ny))
  allocate (tsteps    (1:nx, 1:ny))  
  allocate (var       (1:nvars, 1:nx, 1:ny, 1:M+1) )
  allocate (accum_var (1:nvars, 1:nx, 1:ny, 1:M+1) )
  day_old  (:,:) = day
  tsteps   (:,:) = 0.d0  
  accum_var(:,:,:,:) = 0.d0
  var      (:,:,:,:) = 0.d0
endif

var(1,ix,iy,1:M+1) = Tw1 (1:M+1)
var(2,ix,iy,1:M+1) = Sal1(1:M+1)
var(3,ix,iy,1:M)   = E1  (1:M)
var(4,ix,iy,1:M)   = eps1(1:M)
var(5,ix,iy,1:M)   = KT(1:M)/cw_m_row0
var(6,ix,iy,1:M+1) = sqrt(u1(1:M+1)**2 + v1(1:M+1)**2)
var(7,ix,iy,1:M)   = Ri(1:M)

accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) + var(:,ix,iy,:)
tsteps(ix,iy) = tsteps(ix,iy) + 1.d0

if (day_old(ix,iy)/=day) then
  accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:)/tsteps(ix,iy)
  call DATEMINUS (2,year,month,day,hour,year1,month1,day1,hour1)
  call TIMESTR (8,year1,month1,day1,hour1,timestring)
  call CHECK_UNIT(lake_day_out_unit_min,lake_day_out_unit_max,out_unit)
  write (coords_point, '(2i3)') ix, iy
  open(out_unit, file=outpath(1:len_trim(outpath))//'daily/'// &
  & 'Profiles'//coords_point//timestring//'.dat',status='unknown')
  write (out_unit,*) '1 - depth, m' 
  write (out_unit,*) '2 - temperature, C'
  write (out_unit,*) '3 - salinity, kg/kg' 
  write (out_unit,*) '4 - turbulent kinetic energy, m**2/s**2'
  write (out_unit,*) '5 - dissipation rate, m**2/s**3'
  write (out_unit,*) '6 - kinematic heat conductance, m**2/s'
  write (out_unit,*) '7 - modulus of horizontal velocity, m/s'
  write (out_unit,*) '8 - Richardson number, n/d'
  do i = 1, M+1 
    write (out_unit,'(f12.6,f11.5,6e12.4)') &
    &  -dzeta_int(i)*h1, accum_var(1:nvars,ix,iy,i)
  enddo
  close(out_unit)
  day_old(ix,iy) = day
  accum_var(:,ix,iy,:) = 0.d0
  tsteps(ix,iy) = 0.d0
endif

if (firstcall) firstcall=.false.
END SUBROUTINE DAY_OUT


SUBROUTINE HOUR_OUT(ix,iy,nx,ny,year,month,day,hour,time)

use ARRAYS_WATERSTATE, only : Tw1, Sal1
use ARRAYS_TURB!, only : E1, eps1, row, H_mixed_layer, H_entrainment, Wedderburn
use ARRAYS_BATHYM, only : h1
use ARRAYS_GRID, only : dzeta_int, dzeta_05int
use ARRAYS, only : u1, v1
use TIMEVAR, only : hour_sec

use DRIVING_PARAMS, only : M, scale_output

use INOUT_PARAMETERS, only : &
& lake_hour_out_unit_min, &
& lake_hour_out_unit_max

implicit none

! Input variables
! Reals
real(kind=ireals), intent(in) :: hour
real(kind=ireals), intent(in) :: time

! Integers
integer(kind=iintegers), intent(in) :: ix, iy
integer(kind=iintegers), intent(in) :: nx, ny
integer(kind=iintegers), intent(in) :: year, month, day

! Local variables
! Reals
real(kind=ireals), allocatable :: accum_var (:,:,:,:)
real(kind=ireals), allocatable :: var       (:,:,:,:)
real(kind=ireals), allocatable :: var_scalar(:,:,:)
real(kind=ireals), allocatable :: accum_var_scalar(:,:,:)

real(kind=ireals), external :: DZETA

! Integers
integer(kind=iintegers), parameter :: n_var = 23
integer(kind=iintegers), parameter :: n_var_scalar = 10
integer(kind=iintegers), allocatable :: hour_old(:,:)
integer(kind=iintegers), allocatable :: tsteps(:,:)
integer(kind=iintegers), allocatable :: n_unit(:,:)
integer(kind=iintegers) :: i
integer(kind=iintegers) :: out_unit = lake_hour_out_unit_min

character :: month1*2
character :: year1*4
character :: day1*2
character :: hour1*2
character :: timestring*10
character :: coords_point*6
character :: format_char*100

logical :: firstcall = .true.
logical, allocatable :: firstcall_ixiy(:,:)

      
SAVE

if (firstcall) then
  allocate (accum_var (n_var,1:nx,1:ny,1:M+1) )
  allocate (var       (n_var,1:nx,1:ny,1:M+1) )
  allocate (accum_var_scalar (n_var_scalar,1:nx,1:ny) )
  allocate (var_scalar       (n_var_scalar,1:nx,1:ny) )  
  allocate (tsteps    (1:nx,1:ny) )
  allocate (hour_old  (1:nx,1:ny) )
  allocate (n_unit(1:nx,1:ny))
  allocate (firstcall_ixiy(1:nx,1:ny))
  accum_var(:,:,:,:) = 0.d0
  accum_var_scalar(:,:,:) = 0.d0
  var(:,:,:,:) = 0.d0
  var_scalar(:,:,:) = 0.d0
  tsteps(:,:) = 0
  hour_old(:,:) = int(hour)
  firstcall_ixiy(:,:) = .true.
  n_unit(:,:) = out_unit + 1
endif

var(1,ix,iy,1:M)   = PEMF          (1:M)
var(2,ix,iy,1:M+1) = PDENS_DOWN    (1:M+1)
var(3,ix,iy,1:M+1) = PT_DOWN       (1:M+1)
var(4,ix,iy,1:M+1) = PSAL_DOWN     (1:M+1)
var(5,ix,iy,1:M)   = k_turb_T_flux (1:M)
var(6,ix,iy,1:M)   = T_massflux    (1:M)
var(7,ix,iy,1:M+1) = row           (1:M+1)
var(8,ix,iy,1:M+1) = Tw1           (1:M+1)
var(9,ix,iy,1:M+1) = Sal1          (1:M+1)
var(10,ix,iy,1:M)  = E1            (1:M)
var(11,ix,iy,1:M)  = eps1          (1:M)

if (scale_output%par == 1) then
  var(12,ix,iy,1)    = H_mixed_layer ! Scale
  var(13,ix,iy,1)    = Buoyancy0     ! Scale
  var(14,ix,iy,1)    = w_conv_scale  ! Scale
  var(15,ix,iy,1)    = T_conv_scale  ! Scale
elseif (scale_output%par == 0) then
  var(12,ix,iy,1)    = 1.d0 ! No scaling
  var(13,ix,iy,1)    = 1.d0 ! No scaling
  var(14,ix,iy,1)    = 1.d0 ! No scaling
  var(15,ix,iy,1)    = 1.d0 ! No scaling
else
  print*, 'Scale_output is ', scale_output%par, &
  & ', must be 0 or 1: STOP'  
  STOP
endif
      
var(16,ix,iy,1:M)  = S              (1:M)
var(17,ix,iy,1:M)  = Gen            (1:M)
var(18,ix,iy,1:M)  = Gen_seiches    (1:M)
var(19,ix,iy,1:M)  = TKE_turb_trans (1:M)
      
var(20,ix,iy,2:M)  = k3_mid         (2:M)

var(21,ix,iy,1:M+1)  = u1(1:M+1)
var(22,ix,iy,1:M+1)  = v1(1:M+1)
var(23,ix,iy,1:M+1)  = sqrt(u1(1:M+1)*u1(1:M+1) + v1(1:M+1)*v1(1:M+1))

var_scalar(1,ix,iy) = S_integr_positive
var_scalar(2,ix,iy) = S_integr_negative
var_scalar(3,ix,iy) = Gen_integr
var_scalar(4,ix,iy) = eps_integr
var_scalar(5,ix,iy) = TKE_balance
var_scalar(6,ix,iy) = H_entrainment
var_scalar(7,ix,iy) = TKE_turb_trans_integr
var_scalar(8,ix,iy) = Wedderburn
var_scalar(9,ix,iy) = LakeNumber
var_scalar(10,ix,iy) = Rossby_rad

accum_var_scalar(:,ix,iy) = accum_var_scalar(:,ix,iy) + var_scalar(:,ix,iy)
accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) + var(:,ix,iy,:)
tsteps(ix,iy) = tsteps(ix,iy) + 1

if (int(hour)/=hour_old(ix,iy)) then
  accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) / real(tsteps(ix,iy))
  accum_var_scalar(:,ix,iy) = accum_var_scalar(:,ix,iy) / real(tsteps(ix,iy))
  call DATEMINUS(3,year,month,day,hour,year1,month1,day1,hour1)
  call TIMESTR(10,year1,month1,day1,hour1,timestring)
  call CHECK_UNIT(lake_hour_out_unit_min,lake_hour_out_unit_max,out_unit)
  write (coords_point, '(2i3)') ix, iy
  
! Writing to the file Profiles<yyyymmddhh>.dat  
  open(out_unit, file=outpath(1:len_trim(outpath))// &
  & 'hourly/'//'Profiles'//coords_point//timestring//'.dat', &
  &  status='unknown')
  write (out_unit,*) '1 - depth, normalized' 
  write (out_unit,*) '2 - temperature, normalized'
  write (out_unit,*) '3 - salinity, kg/kg' 
  write (out_unit,*) '4 - water density, kg/m**3'  
  write (out_unit,*) '5 - turbulent kinetic energy, normalized'
  write (out_unit,*) '6 - disspation rate, normalized'
  write (out_unit,*) '7 - eddy viscosity for TKE, m**2/s'
  write (out_unit,*) '8  - mass flux, m/s'
  write (out_unit,*) '9  - downdraft temperature, C' 
  write (out_unit,*) '10 - downdraft salinity, kg/kg' 
  write (out_unit,*) '11 - downdraft density, kg/m**3'
  write (out_unit,*) '12 - x-component of speed, m/s'
  write (out_unit,*) '13 - y-component of speed, m/s'
  write (out_unit,*) '14 - modulus of speed, m/s'
  do i=1,M 
    write (out_unit,'(f14.6,f12.5,e10.3,f10.3,3e10.3,4e14.7,3f11.5)')      &
    & -dzeta_int(i)*h1                         /accum_var(12,ix,iy,1),  &
    & (accum_var(8,ix,iy,i) - scale_output%par*maxval(accum_var(8,ix,iy,:)) ) &
    &  /accum_var(15,ix,iy,1), &
    & accum_var(9,ix,iy,i),                                               &
    & accum_var(7,ix,iy,i),                                               &
    & accum_var(10,ix,iy,i)    / (accum_var(13,ix,iy,1)*accum_var(12,ix,iy,1) /       &
    &  accum_var(14,ix,iy,1) ),                                           &
    & accum_var(11,ix,iy,i)    / accum_var(13,ix,iy,1),                   &
    & accum_var(20,ix,iy,i),                                              & !accum_var(10,i)**2 / accum_var(11,i),
    & accum_var(1,ix,iy,i),                                               &
    & accum_var(3,ix,iy,i),                                               &
    & accum_var(4,ix,iy,i),                                               &
    & accum_var(2,ix,iy,i),                                               &
    & accum_var(21,ix,iy,i),                                              &
    & accum_var(22,ix,iy,i),                                              &
    & accum_var(23,ix,iy,i)
  enddo 
  close (out_unit)
  
! Writing to the file EDMF_profiles<yyyymmddhh>.dat
  open(out_unit, file=outpath(1:len_trim(outpath))//  &
  & 'hourly/'//'EDMF_profiles'//coords_point//timestring//'.dat',   &
  & status='unknown')
  write (out_unit,*) '1 - depth, normalized by mixed layer depth' 
  write (out_unit,*) '2 - "k - flux" of temperature, normalized'
  write (out_unit,*) '3 - mass flux of temperature,  normalized'
  write (out_unit,*) '4 - total flux of temperature, normalized'
  write (out_unit,'( 4(i10,6x) )') 1,2,3,4
  do i=1, M
    write (out_unit,'(f18.6,3e16.7)')                   &
    & -dzeta_05int(i)*h1/ accum_var(12,ix,iy,1),          &
    & accum_var(5,ix,iy,i)/(accum_var(14,ix,iy,1)*accum_var(15,ix,iy,1)), &
    & accum_var(6,ix,iy,i)/(accum_var(14,ix,iy,1)*accum_var(15,ix,iy,1)), &
    & (accum_var(5,ix,iy,i) + accum_var(6,ix,iy,i) ) /              &
    & (accum_var(14,ix,iy,1)*accum_var(15,ix,iy,1))
  enddo
  close (out_unit)
  
! Writing to the file TKE_budget<yyyymmddhh>.dat
  open(out_unit, file=outpath(1:len_trim(outpath))//    &
  & 'hourly/'//'TKE_budget'//coords_point//timestring//'.dat',        &
  &  status='unknown')
  write (out_unit,*) '1 - depth, normalised with mixed layer depth'
  write (out_unit,*) '2 - shear production, normalised with surface buoyancy flux'
  write (out_unit,*) '3 - seiches production, normalised with surface buoyancy flux'
  write (out_unit,*) '4 - buoyancy source,  normalised with surface buoyancy flux'
  write (out_unit,*) '5 - dissipation rate, normalised with surface buoyancy flux'
  write (out_unit,*) '6 - turbulent transport, normalised with surface buoyancy flux'
  write (out_unit,'( 6(i10,6x) )') 1,2,3,4,5,6
  do i = 1, M
    write (out_unit,'(f18.6,5e16.7)')          &
    & -dzeta_05int(i)*h1/ accum_var(12,ix,iy,1), &
    & accum_var(17,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(18,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(16,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(11,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(19,ix,iy,i) / accum_var(13,ix,iy,1)
  enddo
!  write (out_unit, '(a)')
  close (out_unit)
  
  if (firstcall_ixiy(ix,iy)) then
    call CHECK_UNIT(lake_hour_out_unit_min,lake_hour_out_unit_max,n_unit(ix,iy))
    open(n_unit(ix,iy), file=outpath(1:len_trim(outpath))//           &
    & 'hourly/'//'TKE_integr_conv'//coords_point//timestring//'.dat', &
    &  status='unknown')
    write (n_unit(ix,iy),*) 'Col. 1 - year'
    write (n_unit(ix,iy),*) 'Col. 2 - month'
    write (n_unit(ix,iy),*) 'Col. 3 - day'
    write (n_unit(ix,iy),*) 'Col. 4 - hour'
    write (n_unit(ix,iy),*) 'Col. 5 - the time from the start of integration, hours'
    write (n_unit(ix,iy),*) 'Col. 6 - H_entrainment, m    '
    write (n_unit(ix,iy),*) 'Col. 7 - B_integr+,   m**3/s**3'
    write (n_unit(ix,iy),*) 'Col. 8 - B_integr-,   m**3/s**3'
    write (n_unit(ix,iy),*) 'Col. 9 - S_integr,    m**3/s**3'
    write (n_unit(ix,iy),*) 'Col. 10 - eps_integr,  m**3/s**3'
    write (n_unit(ix,iy),*) 'Col. 11 - TKE turb trans integr,  m**3/s**3'
    write (n_unit(ix,iy),*) 'Col. 12 - TKE_balance, m**3/s**3'
    write (n_unit(ix,iy),*) 'Col. 13 - Wedderburn number, n/d'
    write (n_unit(ix,iy),*) 'Col. 14 - Lake number, n/d'
    write (n_unit(ix,iy),*) 'Col. 15 - Internal Rossby radius, m'
    firstcall_ixiy(ix,iy) = .false.
  endif
  
  format_char = '(3i7,f8.2,f13.2,f10.2,9e15.6)'
  write (n_unit(ix,iy), format_char) &
  &                    year,month,day,hour,time/hour_sec, &
  &                    accum_var_scalar(6,ix,iy), & 
  &                    accum_var_scalar(1,ix,iy), &
  &                    accum_var_scalar(2,ix,iy), &
  &                    accum_var_scalar(3,ix,iy), &
  &                    accum_var_scalar(4,ix,iy), &
  &                    accum_var_scalar(7,ix,iy), &
  &                    accum_var_scalar(5,ix,iy), &
  &                    accum_var_scalar(8,ix,iy), &
  &                    accum_var_scalar(9,ix,iy), &
  &                    accum_var_scalar(10,ix,iy)
  
  hour_old(ix,iy) = int(hour)
  accum_var(:,ix,iy,:) = 0.d0
  accum_var_scalar(:,ix,iy) = 0.d0
  tsteps(ix,iy)  = 0
   
endif

if (firstcall) firstcall = .false.
!if (firstcall_ixiy(ix,iy)) firstcall_ixiy(ix,iy) = .false.
END SUBROUTINE HOUR_OUT          
      


SUBROUTINE EVERYSTEP_OUT(ix,iy,nx,ny)

use ARRAYS!, only : time, nstep

use ARRAYS_WATERSTATE!, only : Tw1, Sal1 
use ARRAYS_BATHYM!, only : h1
use ARRAYS_GRID!, only : dzeta_int, M

use ARRAYS_TURB!, only : &
! & S_integr_positive, &
! & S_integr_negative, &
! & Gen_integr, &
! & eps_integr, &
! & TKE_balance, &
! & E_integr, &
! & H_entrainment, &
! & E1, eps1, &
! & k_turb_T_flux

use DRIVING_PARAMS, only : everystep, M

use INOUT_PARAMETERS, only : &
& lake_everystep_out_unit_min, &
& lake_everystep_out_unit_max

implicit none

! Input variables
! Integers
integer(kind=iintegers), intent(in) :: ix, iy
integer(kind=iintegers), intent(in) :: nx, ny

! Local variables
! Reals
real(kind=ireals), external :: DZETA
real(kind=ireals), parameter :: ACC = 1.d-20
real(kind=ireals), parameter :: hour_sec = 60.*60.

! Integers
integer(kind=iintegers), allocatable :: n_unit(:,:,:)
integer(kind=iintegers) :: i ! Loop index

! Characters
character :: coords_point*6
character :: format_char*100

! Logicals
logical, allocatable :: firstcall(:,:)

      
SAVE

if (.not.allocated(firstcall)) then
  allocate (firstcall(1:nx, 1:ny) )
  allocate (n_unit   (1:2, 1:nx, 1:ny) )
  firstcall = .true.
  n_unit = lake_everystep_out_unit_min
endif  
      
if (firstcall(ix,iy)) then
  write (coords_point, '(2i3)') ix, iy
  
  if (everystep%par /= 2) then
    call CHECK_UNIT(lake_everystep_out_unit_min,lake_everystep_out_unit_max, &
    & n_unit(1,ix,iy))
    open (n_unit(1,ix,iy), file=outpath(1:len_trim(outpath))//'everystep/'// &
    & 'Profiles'//coords_point//'.dat',  status='unknown')
    write (n_unit(1,ix,iy),*) '1 - depth, m' 
    write (n_unit(1,ix,iy),*) '2 - temperature, C'
    write (n_unit(1,ix,iy),*) '3 - salinity, kg/kg' 
    write (n_unit(1,ix,iy),*) '4 - turbulent kinetic energy, m**2/s**2'
    write (n_unit(1,ix,iy),*) '5 - disspation rate, m**2/s**3'
    write (n_unit(1,ix,iy),*) '6 - eddy diffusivity (TKE**2/dissipation), m**2/s'
    write (n_unit(1,ix,iy),*) '7 - k-flux of temperature, K*m/s'
  endif
    
  call CHECK_UNIT(lake_everystep_out_unit_min,lake_everystep_out_unit_max, &
  & n_unit(2,ix,iy))
  open (n_unit(2,ix,iy), file=outpath(1:len_trim(outpath))//'everystep/'// &
  & 'TKE_integr'//coords_point//'.dat',  status='unknown')
  write (n_unit(2,ix,iy),*) '1 - timestep              '
  write (n_unit(2,ix,iy),*) '2 - time,        hours    '
  write (n_unit(2,ix,iy),*) '3 - B_integr+,   m**3/s**3'
  write (n_unit(2,ix,iy),*) '4 - B_integr-,   m**3/s**3'
  write (n_unit(2,ix,iy),*) '5 - S_integr,    m**3/s**3'
  write (n_unit(2,ix,iy),*) '6 - eps_integr,  m**3/s**3'
  write (n_unit(2,ix,iy),*) '7 - TKE_balance, m**3/s**3'
  write (n_unit(2,ix,iy),*) '8 - E_integr,    m**3/s**2'
endif

if (everystep%par /= 2) then      
  write (n_unit(1,ix,iy),*) 'nstep = ', nstep
  format_char = '(f10.6,f9.3,5e12.4)'
  do i=1,M 
    write (n_unit(1,ix,iy), format_char) &
    & -dzeta_int(i)*h1,Tw1(i),Sal1(i),E1(i),eps1(i),  &
    & E1(i)**2/(eps1(i)+ACC), k_turb_T_flux(i)
  enddo
endif

format_char = '(i10, 2f10.2, 6e15.6)'
write(n_unit(2,ix,iy), format_char) &
& nstep, time/hour_sec, H_entrainment, S_integr_positive, S_integr_negative, &
& Gen_integr, eps_integr, TKE_balance, E_integr
      
if (firstcall(ix,iy)) firstcall(ix,iy)=.false.
END SUBROUTINE EVERYSTEP_OUT
      
      
SUBROUTINE SERIES_OUT(ix,iy,nx,ny,year,month,day,hour,tsw)

use DRIVING_PARAMS, only : M, ns, Mice, dt_out
use NUMERIC_PARAMS , only : &
& ms, small_value

use ARRAYS_TURB, only : row, H_mixed_layer, H_entrainment, &
& signwaveheight, w_conv_scale, u_star, Ri_bulk, maxN, i_maxN, &
& ThermThick, ReTherm, RiTherm

use ARRAYS_WATERSTATE, only : Tw1, Ti1, lamw, salice

use ARRAYS_GRID, only : nsoilcols

use ARRAYS_BATHYM, only : &
& h1,hx1,hx2,hy1,hy2, &
& hx1t,hx2t,hy1t,hy2t, &
& l1,hs1,ls1,  &
& voldef, vol,     &
& area_int, bathymwater

use ARRAYS_METHANE, only : &
& fbbleflx_ch4_sum, &
& fbbleflx_ch4,    &
& febul0, fplant,   &
& fdiff_lake_surf, &
& fdiffbot, &
& qwater, qsoil, &
& rprod_total_newC, &
& rprod_total_oldC, &
& h_talik, lammeth, &
& qwateroxidtot, qwateroxidML

use ARRAYS_OXYGEN, only : &
& fbbleflx_co2_sum, &
& fbbleflx_co2,    &
& fbbleflx_o2_sum, &
& fbbleflx_o2, &
& oxyg, &
& fco2, fo2

use ARRAYS_SOIL, only : Tsoil1, lsm

use ARRAYS, only:  &
& Tskin,       &
& zinv,            &
& time,            &
& T_0dim,           &
& roughness,emissivity,albedo,aM,bM,relhums

use ATMOS,  only:    &
& eflux0_kinem,      &
& eflux0,            &
& turb_density_flux, &
& hw, xlew, botflux, &
& tau, windwork, surfrad

use INOUT_PARAMETERS, only : &
& lake_series_out_unit_min, &
& lake_series_out_unit_max

use METH_OXYG_CONSTANTS, only : &
& molmass_ch4

use TIMEVAR, only : &
hour_sec, day_sec, hour_min

use PHYS_CONSTANTS, only : &
Kelvin0, cw_m_row0, roa0

use DZETA_MOD, only : VARMEAN

use ARRAYS_GRID, only : gsp

use SNOWSOIL_MOD, only : Tsn, itop

implicit none

! Input variables
! Reals
real(kind=ireals), intent(in) :: tsw
real(kind=ireals), intent(in) :: hour

! Integers
integer(kind=iintegers), intent(in) :: ix, iy
integer(kind=iintegers), intent(in) :: nx, ny
integer(kind=iintegers), intent(in) :: year, month, day

! Local variables
! Reals

real(kind=ireals), parameter :: mfs = molmass_ch4*8.64d+7 ! the transform multiplier for methane flux
                                                          ! from mol/(m**2*s) to mg/(m**2*day)
real(kind=ireals), parameter :: small_number = 1.d-5

real(kind=ireals) :: T_mean, voldef0 = 0.
real(kind=ireals) :: methfluxML, methfluxshalsed

real(kind=ireals), allocatable :: tsteps(:,:)
real(kind=ireals), allocatable :: var_scalar (:,:,:) 
real(kind=ireals), allocatable :: accum_var_scalar(:,:,:)

! Integers
integer(kind=iintegers), parameter :: nfiles = 5
integer(kind=iintegers), parameter :: numaccum = 5 ! The number of scalar variables being accumulated (averaged)
integer(kind=iintegers), allocatable :: n_unit(:,:,:)
integer(kind=iintegers), allocatable :: count_out(:,:)

! Characters
character :: coords_point*6
character :: format_char*100, workchar*10

! Logicals
logical, allocatable :: firstcall(:,:)


!real(kind=ireals) :: dz(ms)
integer(kind=iintegers) :: i
!common /SOILDAT/ dz,itop
!real(kind=ireals) :: Tsn(1:ms), cs(1:ms)
!common /SNOW_CHAR/ Tsn,cs

      
SAVE

if (.not.allocated(firstcall)) then
  allocate (firstcall(1:nx, 1:ny) )
  allocate (count_out(1:nx, 1:ny) )
  allocate (n_unit   (1:nfiles, 1:nx, 1:ny) )
  allocate (tsteps           (1:nx, 1:ny))  
  allocate (var_scalar       (1:numaccum, 1:nx, 1:ny) )
  allocate (accum_var_scalar (1:numaccum, 1:nx, 1:ny) )

  tsteps(:,:) = 0.d0
  var_scalar(:,:,:) = 0.d0
  accum_var_scalar(:,:,:) = 0.d0

  firstcall(:,:) = .true.
  n_unit = lake_series_out_unit_min
endif  
            
if (firstcall(ix,iy)) then
  write (coords_point, '(2i3)') ix, iy

  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(1,ix,iy))
  open (n_unit(1,ix,iy),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'layers'//coords_point//'.dat',  status='unknown')
  write (n_unit(1,ix,iy),*)'Col. 1 - year'
  write (n_unit(1,ix,iy),*)'Col. 2 - month'
  write (n_unit(1,ix,iy),*)'Col. 3 - day'
  write (n_unit(1,ix,iy),*)'Col. 4 - hour'
  write (n_unit(1,ix,iy),*)'Col. 5 - the time from the start of integration, hours'
  write (n_unit(1,ix,iy),*)'Col. 6 - water layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 7 - W mixed layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 8 - E mixed layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 9 - S mixed layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 10 - N mixed layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 11 - W lower layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 12 - E lower layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 13 - S lower layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 14 - N lower layer thickness, m'
  write (n_unit(1,ix,iy),*)'Col. 15 - ice layer thickness,   m'
  write (n_unit(1,ix,iy),*)'Col. 16 - snow layer thickness,  m'
  write (n_unit(1,ix,iy),*)'Col. 17 - bottom ice thickness,  m'
  write (n_unit(1,ix,iy),*)'Col. 18 - reservoir volume,  m**3'
  write (n_unit(1,ix,iy),*)'Col. 19 - volume deficit (accumulated),  m**3'
    
  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(2,ix,iy))
  open (n_unit(2,ix,iy),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'T_fluxes'//coords_point//'.dat',status='unknown')
  write (n_unit(2,ix,iy),*)'Col. 1 - year'
  write (n_unit(2,ix,iy),*)'Col. 2 - month'
  write (n_unit(2,ix,iy),*)'Col. 3 - day'
  write (n_unit(2,ix,iy),*)'Col. 4 - hour'
  write (n_unit(2,ix,iy),*)'Col. 5 - the time from the start of integration, hours'
  write (n_unit(2,ix,iy),*)'Col. 6 - surface temperature, C'
  write (n_unit(2,ix,iy),*)'Col. 7 - water skin temperature, C'
  write (n_unit(2,ix,iy),*)'Col. 8 - water surface temperature, C'
  write (n_unit(2,ix,iy),*)'Col. 9 - mean temperature of water coloumn, C'
  write (n_unit(2,ix,iy),*)'Col. 10 - maximal temperature in the water coloumn, C'
  write (n_unit(2,ix,iy),*)'Col. 11 - zero-dimensional model temperature, C'  
  write (n_unit(2,ix,iy),*)'Col. 12 - upper ice surface temperature, C'
  write (n_unit(2,ix,iy),*)'Col. 13 - upper snow surface temperature, C'
  write (n_unit(2,ix,iy),*)'Col. 14 - sensible heat flux,    W/m**2'
  write (n_unit(2,ix,iy),*)'Col. 15 - latent heat flux,      W/m**2'
  write (n_unit(2,ix,iy),*)'Col. 16 - downward heat flux at the upper lake surface, W/m**2'
  write (n_unit(2,ix,iy),*)'Col. 17 - downward heat flux at the lake bottom, W/m**2'
  write (n_unit(2,ix,iy),*)'Col. 18 - friction velocity in the surface air layer, m/s'
  write (n_unit(2,ix,iy),*)'Col. 19 - wind work at the water surface, W/m**2'
  write (n_unit(2,ix,iy),*)'Col. 20 - albedo of the lake-atmosphere interface, n/d'
  write (n_unit(2,ix,iy),*)'Col. 21 - significant wave height, m'
  write (n_unit(2,ix,iy),*)'Col. 22 - bottom ice salinity, kg/kg'
  write (n_unit(2,ix,iy),*)'Cols. 23... - top soil columns temperature, m'
  
  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(3,ix,iy))
  open (n_unit(3,ix,iy),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'conv_series'//coords_point//'.dat',status='unknown')
  write (n_unit(3,ix,iy),*)'Col. 1 - year'
  write (n_unit(3,ix,iy),*)'Col. 2 - month'
  write (n_unit(3,ix,iy),*)'Col. 3 - day'
  write (n_unit(3,ix,iy),*)'Col. 4 - hour'
  write (n_unit(3,ix,iy),*)'Col. 5 - the time from the start of integration, hours'
  write (n_unit(3,ix,iy),*)'Col. 6 - surface temperature, C'
  write (n_unit(3,ix,iy),*)'Col. 7 - surface density, kg/m**3'
  write (n_unit(3,ix,iy),*)'Col. 8 - heat flux downwards at the surface, K*m/s'
  write (n_unit(3,ix,iy),*)'Col. 9 - turbulent density flux at the surface, kg/m**2/s'
  write (n_unit(3,ix,iy),*)'Col. 10 - inversion depth (mass flux diagnostics),    m'
  write (n_unit(3,ix,iy),*)'Col. 11 - mixed layer depth,                          m'
  write (n_unit(3,ix,iy),*)'Col. 12 - entrainment depth,                          m'
  write (n_unit(3,ix,iy),*)'Col. 13 - convective velocity scale,                m/s'
  write (n_unit(3,ix,iy),*)'Col. 14 - friction velocity,                        m/s'
  write (n_unit(3,ix,iy),*)'Col. 15 - -w_star/u_star,                            m/s'
  write (n_unit(3,ix,iy),*)'Col. 16 - bulk Richardson number,                    n/d'
  write (n_unit(3,ix,iy),*)'Col. 17 - maximal Brunt-Vaisala frequency,       s**(-1)'
  write (n_unit(3,ix,iy),*)'Col. 18 - minimal thermal conductance,            m**2/s'
  write (n_unit(3,ix,iy),*)'Col. 19 - thermocline thickness,                       m'
  write (n_unit(3,ix,iy),*)'Col. 20 - Reynolds number in thermocline ,           n/d'
  write (n_unit(3,ix,iy),*)'Col. 21 - bulk Richardson number in thermocline,     n/d'
  
! The output of LakeMIP file 1 format   
  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(4,ix,iy))
  open (n_unit(4,ix,iy),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'LakeMIP_file1'//coords_point//'.dat',status='unknown')
!  write (n_unit(4,ix,iy),*)'Col. 1 - time, days'
!  write (n_unit(4,ix,iy),*)'Col. 1 - year'
!  write (n_unit(4,ix,iy),*)'Col. 2 - month'
!  write (n_unit(4,ix,iy),*)'Col. 3 - day in a month'
!  write (n_unit(4,ix,iy),*)'Col. 4 - hour'
!  write (n_unit(4,ix,iy),*)'Col. 5 - minute'
!  write (n_unit(4,ix,iy),*)'Col. 6 - mixed-layer temperature &
!  &(equals to the temperature of the upper water surface), degrees Celsius'
!  write (n_unit(4,ix,iy),*)'Col. 7 - mean temperature of the water column, degrees Celsius'
!  write (n_unit(4,ix,iy),*)'Col. 8 - bottom temperature, degrees Celsius'
!  write (n_unit(4,ix,iy),*)'Col. 9 - mixed-layer depth, meters'
!  write (n_unit(4,ix,iy),*)'Col. 10 - ice thickness, meters'
!  write (n_unit(4,ix,iy),*)'Col. 11 - snow thickness, meters'
!  write (n_unit(4,ix,iy),*)'Col. 12 - temperature at the ice upper surface, degrees Celsius'
!  write (n_unit(4,ix,iy),*)'Col. 13 - temperature at the snow upper surface, degrees Celsius'
!  write (n_unit(4,ix,iy),*)'Col. 14 - sensible heat flux at the lake-atmosphere &
!  &interface, averaged over output interval, upwards, W/m**2'
!  write (n_unit(4,ix,iy),*)'Col. 15 - latent heat flux at the lake-atmosphere &
!  &interface, averaged over output interval, upwards, W/m**2'
!  write (n_unit(4,ix,iy),*)'Col. 16 - momentum flux at the lake-atmosphere &
!  &interface, averaged over output interval, positive, N/m**2'
!  write (n_unit(4,ix,iy),*)'Col. 17 - upward long-wave radiation flux &
!  &at the lake-atmosphere interface, averaged over output interval, W/m**2'
!  write (n_unit(4,ix,iy),*)'Col. 18 - downward heat flux at the lake-atmosphere &
!  &interface, averaged over output interval, W/m**2'
!  write (n_unit(4,ix,iy),*)'Col. 19 - surface albedo, n/d'
      
  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(5,ix,iy))
  open (n_unit(5,ix,iy),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'methane_series'//coords_point//'.dat',status='unknown')
  write (n_unit(5,ix,iy),*)'Col. 1 - year'
  write (n_unit(5,ix,iy),*)'Col. 2 - month'
  write (n_unit(5,ix,iy),*)'Col. 3 - day'
  write (n_unit(5,ix,iy),*)'Col. 4 - hour'
  write (n_unit(5,ix,iy),*)'Col. 5 - the time from the start of integration, hours'
  write (n_unit(5,ix,iy),*)'Col. 6 - the talik depth, m'
  write (n_unit(5,ix,iy),*)'Col. 7 - lake surface methane concentration, mol/m**3'
  write (n_unit(5,ix,iy),*)'Col. 8 - lake bottom methane concentration, mol/m**3'
  write (n_unit(5,ix,iy),*)'Col. 9 - soil bottom methane concentration, mol/m**3'
  write (n_unit(5,ix,iy),*)'Col. 10 - lake surface oxygen concentration, mol/m**3'
  write (n_unit(5,ix,iy),*)'Col. 11 - lake bottom oxygen concentration, mol/m**3'
  write (n_unit(5,ix,iy),*)'Col. 12 - total methane production due to young C decomposition, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 13 - total methane production due to old C decomposition, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 14 - methane ebullition flux at the surface, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 15 - methane plant-mediated flux at the lake bottom, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 16 - methane diffusion flux at the lake bottom, upwards, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 17 - methane turbulent flux at the lake surface, upwards, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 18 - methane ebullition flux at the surface, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 19 - methane plant-mediated flux at the lake bottom, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 20 - methane diffusion flux at the lake bottom, upwards, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 21 - methane turbulent flux at the lake surface, upwards, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 22 - methane turbulent flux at the bottom of mixed layer normalized by &
  & surface area, upwards, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 23 - methane flux from sediments in the mixed layer normalized by &
  & surface area, upwards, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 24 - methane bubble flux at the bottom of the mixed layer normalized by &
  & surface area, upwards, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 25 - total methane oxidation in water normalized by &
  & surface area, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 26 - methane oxidation in mixed layer normalized by &
  & surface area, mg/(m**2*day)'
  write (n_unit(5,ix,iy),*)'Col. 27 - co2 turbulent flux at the lake surface, upwards, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 28 - co2 ebullition flux at the surface, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 29 - oxygen turbulent flux at the lake surface, upwards, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 30 - oxygen ebullition flux at the surface, mol/(m**2*s)'
  write (n_unit(5,ix,iy),*)'Col. 31..31+ns-1 - methane ebullition flux at the surface from different soil columns, mg/(m**2*day)'
  
  count_out(ix,iy) = int(time/(dt_out%par*hour_sec))
endif
  
! Accumulating variables
var_scalar(1,ix,iy) = hw
var_scalar(2,ix,iy) = xlew
var_scalar(3,ix,iy) = tau
var_scalar(4,ix,iy) = surfrad
var_scalar(5,ix,iy) = eflux0

accum_var_scalar(:,ix,iy) = accum_var_scalar(:,ix,iy) + var_scalar(:,ix,iy)
tsteps(ix,iy) = tsteps(ix,iy) + 1.d0

! Writing to files
if (int(time/(dt_out%par*hour_sec))>count_out(ix,iy) .or. &
  & abs(time - (count_out(ix,iy) + 1)*dt_out%par*hour_sec) < small_number) then

  T_mean = VARMEAN(Tw1,bathymwater,1_iintegers)
  format_char = '(3i7,f10.4,f15.4,f10.4,8e12.4,3f10.4,1x,e12.4,1x,f10.1)'
  write (n_unit(1,ix,iy),format_char) year,month,day ,hour,  &
  &                                   time/hour_sec ,        &
  &                                   h1, hx1, hx2, hy1, hy2,&
  &                                   hx1t, hx2t, hy1t, hy2t,&
  &                                   l1, hs1, ls1, vol,     &
  &                                   voldef - voldef0
  voldef0 = voldef

  format_char = '(3i7,f8.2,f13.2,14f11.4,f11.2,f11.4,f6.2'
  ! Formatting for possible multiple soil columns in a lake
  write(workchar,'(i2)') nsoilcols
  format_char = format_char(1:len_trim(format_char))//','// &
  &             workchar   (1:len_trim(workchar)   )//'f11.4)'
  write (n_unit(2,ix,iy),format_char) year,month,day,hour, &
  &                                   time/hour_sec ,      &
  &                                   tsw - Kelvin0,Tskin(1), &   
  &                                   Tw1(1),              &
  &                                   T_mean,maxval(Tw1(1:M+1)), &
  &                                   T_0dim,Ti1(1),Tsn(itop), &
  &                                   hw,xlew,eflux0,botflux,sqrt(tau/roa0), &
  &                                   windwork, &
  &                                   albedo,signwaveheight,      &
  &                                   salice(Mice+1),  &
  &                                   Tsoil1(1,1:nsoilcols)

  format_char = '(3i7,f8.2,f13.2,2f10.4,2e15.5,5f10.4,7e15.5)'
  write (n_unit(3,ix,iy),format_char) year,month,day,hour,        &
  &                                   time/hour_sec,              &
  &                                   Tw1(1),row(1),eflux0_kinem, &
  &                                   turb_density_flux, zinv,    &
  &                                   H_mixed_layer, H_entrainment, &
  &                                   w_conv_scale, u_star,        &
  &                                   min(max(- w_conv_scale/(u_star + small_value), &
  &                                   -1.e+3_ireals),1.e+3_ireals), &
  &                                   min(max(Ri_bulk, -1.e+3_ireals),1.e+3_ireals), &
  &                                   maxN, minval(lamw(max(i_maxN-3,1_iintegers):min(i_maxN+3,M)))/cw_m_row0, &
  &                                   ThermThick, ReTherm, RiTherm

  format_char = '(15(e15.8,1x))'! '(19(e15.8,1x))' !'(f10.5,10f8.2,f10.4,3f8.2)'
  accum_var_scalar(:,ix,iy) = accum_var_scalar(:,ix,iy)/tsteps(ix,iy)
  write (n_unit(4,ix,iy),format_char) & !float(year), float(month), float(day), hour, &
  &                                   time/day_sec, & !max((hour - floor(hour+small_number))*hour_min,0.d0), &
  &                                   Tw1(1), T_mean, &   
  &                                   Tw1(M+1), H_mixed_layer, &
  &                                   l1, hs1, Ti1(1), Tsn(itop), & !missing_value, missing_value, missing_value, missing_value,
  &                                   accum_var_scalar(1:numaccum,ix,iy), &
  &                                   albedo
  accum_var_scalar(:,ix,iy) = 0.d0
  tsteps(ix,iy) = 0.d0

  ! Methane flux at the bottom of mixed layer
  if (h1 > 0.) then
    methfluxML = lammeth(i_maxN)*(qwater(i_maxN+1,2) - qwater(i_maxN,2))/(gsp%ddz(i_maxN)*h1)*&
    & bathymwater(i_maxN)%area_half/bathymwater(1)%area_int
    ! Methane flux from shallow sediments (i.e. sediments in the mixed layer)
    methfluxshalsed = 0.
    do i = i_maxN, 1, -1
      methfluxshalsed = methfluxshalsed + bathymwater(i)%area_int * lsm%water(i) * gsp%ddz05(i-1)*h1
    enddo
    methfluxshalsed = methfluxshalsed/bathymwater(1)%area_int
  else
    methfluxML = 0.
    methfluxshalsed = 0.
  endif
  

  format_char = '(3i7,f8.2,f13.2,f8.2,30e15.5)'
  write (n_unit(5,ix,iy),format_char) year,month,day,hour,        &
  &                                   time/hour_sec, h_talik,     &
  &                                   qwater(1,1), qwater(M+1,1), &
  &                                   qsoil(ns,nsoilcols), oxyg(1,1),       &
  &                                   oxyg(M+1,1),                &
  &                                   rprod_total_newC(nsoilcols), &
  &                                   rprod_total_oldC(nsoilcols), &
  &                                   fbbleflx_ch4_sum(0)/area_int(1), fplant,   &
  &                                 - fdiffbot(nsoilcols), - fdiff_lake_surf, &
  &                                   fbbleflx_ch4_sum(0)/area_int(1)*mfs, fplant*mfs,      &
  &                                 - fdiffbot(nsoilcols)*mfs,    &
  &                                 - fdiff_lake_surf*mfs,        &
  &                                   methfluxML*mfs, methfluxshalsed*mfs, &
  &                                   fbbleflx_ch4_sum(i_maxN)/area_int(1)*mfs, &
  &                                   qwateroxidtot*mfs, qwateroxidML*mfs, &
  &                                   fco2, fbbleflx_co2_sum(0)/area_int(1), &
  &                                   fo2,  fbbleflx_o2_sum(0) /area_int(1), &      
  &                                   fbbleflx_ch4(0,1:nsoilcols)*febul0(1:nsoilcols)*mfs
!  format_char = '(3i7,f8.2,f13.2,f8.2,24e15.5)' ! two-meth
!  write (n_unit(5,ix,iy),format_char) year,month,day,hour,        &  ! two-meth
!  &                                   time/hour_sec, h_talik,     &  ! two-meth
!  &                                   qwater2(1,1,1:2),            &  ! two-meth
!  &                                   qwater2(M+1,1,1:2),          &  ! two-meth
!  &                                   qsoil2(ns,1:2),              &  ! two-meth
!  &                                   oxyg(1,1),                  &  ! two-meth
!  &                                   oxyg(M+1,1),                &  ! two-meth
!  &                                   rprod_total_newC,           &  ! two-meth
!  &                                   rprod_total_oldC,           &  ! two-meth
!  &                                   febul2(1:2), fplant,        &  ! two-meth
!  &                                   - fdiff2(1:2), - fdiff_lake_surf2(1:2), &  ! two-meth
!  &                                   febul2(1:2)*mfs, fplant*mfs,&  ! two-meth
!  &                                   - fdiff2(1:2)*mfs,          &  ! two-meth
!  &                                   - fdiff_lake_surf2(1:2)*mfs    ! two-meth
  
  count_out(ix,iy) = count_out(ix,iy) + 1 
endif 
      
if (firstcall(ix,iy)) firstcall(ix,iy)=.false.
END SUBROUTINE SERIES_OUT


SUBROUTINE PROFILE_OUTPUT &
& (ix, iy, nx, ny, &
&  year, month, day, hour, &
&  time, dt_out, &
&  Profile, z_prof, nprof, &
&  z_out, nout, &
&  outpath1, filename, &
&  unitindex, ndec, last, Profile2, Profile3)

! The subroutine WATER_TEMP_OUT writes the water temperature profiles
! with the timestep dt_out

use INOUT_PARAMETERS, only : &
& lake_series_out_unit_min, &
& lake_series_out_unit_max

use TIMEVAR, only : &
hour_sec, day_sec, hour_min

implicit none

! Input variables
! Reals
real(kind=ireals), intent(in) :: Profile(1:nprof) ! The profile
real(kind=ireals), intent(in) :: z_prof(1:nprof) ! The levels of the profile elements
real(kind=ireals), intent(in) :: hour ! The hour
real(kind=ireals), intent(in) :: z_out(1:nout) ! The output grid
real(kind=ireals), intent(in) :: time, dt_out ! The time from beginning of integration
                                    ! The timestep of output
! Integers
integer(kind=iintegers), intent(in) :: ix, iy ! The current coordinates of lake point
integer(kind=iintegers), intent(in) :: nx, ny ! The maximal coordinates of lake points
integer(kind=iintegers), intent(in) :: year, month, day ! Year, month and day
integer(kind=iintegers), intent(in) :: nprof ! The number of layers in water, (M+1) - number of levels
integer(kind=iintegers), intent(in) :: nout
integer(kind=iintegers), intent(in) :: unitindex
integer(kind=iintegers), intent(in) :: ndec ! Number of decimal digits

character(len=*), intent(in) :: outpath1
character(len=*), intent(in) :: filename

logical, intent(in) :: last

real(kind=ireals), intent(in), optional :: Profile2(1:nprof), Profile3(1:nprof) ! Optional profiles

! Local variables
! Reals

real(kind=ireals), parameter :: small_number = 1.d-5
real(kind=ireals), allocatable :: Profile_out(:,:)


! Integers
integer(kind=iintegers), parameter :: unitindex_max = 100
integer(kind=iintegers), allocatable, save :: n_unit(:,:,:)
integer(kind=iintegers), allocatable, save :: count_out(:,:)
integer(kind=iintegers) :: i, k !Loop index
integer(kind=iintegers) :: nvars
!integer(kind=iintegers) :: nform(1:2)

! Characters
character(len=6) :: coords_point
character(len=100) :: formats(1:2)

! Logicals
logical, allocatable, save :: firstcall(:,:)
logical :: flag


if (.not.allocated(firstcall)) then
  allocate (firstcall(1:nx, 1:ny) )
  allocate (count_out(1:nx, 1:ny) )
  allocate (n_unit   (1:unitindex_max, 1:nx, 1:ny) )
  firstcall(:,:) = .true.
  n_unit = lake_series_out_unit_min
endif  

nvars = 1
if (present(Profile2)) nvars = 2
if (present(Profile2) .and. present(Profile3)) nvars = 3
            
if (firstcall(ix,iy)) then
  write (coords_point, '(2i3)') ix, iy

  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(unitindex,ix,iy))
  open (n_unit(unitindex,ix,iy), &
  & file=outpath1(1:len_trim(outpath1))//'time_series/'// &
  & filename(1:len_trim(filename))//coords_point//'.dat',  status='unknown')
  write (n_unit(unitindex,ix,iy),*)'Col. 1 - year'
  write (n_unit(unitindex,ix,iy),*)'Col. 2 - month'
  write (n_unit(unitindex,ix,iy),*)'Col. 3 - day'
  write (n_unit(unitindex,ix,iy),*)'Col. 4 - hour'
  write (n_unit(unitindex,ix,iy),*) &
  & 'Col. 5 - the time from the start of integration, days' !& 'Col. 5 - the time from the start of integration, days'
  write (n_unit(unitindex,ix,iy),*) &
  & 'Col. [6...] - Depth1, Temp1, Depth2, Temp2, ..., DepthN, TempN'
 
  call CHECK_UNIT(lake_series_out_unit_min,lake_series_out_unit_max, &
  & n_unit(unitindex+1,ix,iy))
  open (n_unit(unitindex+1,ix,iy), &
  & file=outpath1(1:len_trim(outpath1))//'time_series/'// &
  & filename(1:len_trim(filename))//coords_point//'f2'//'.dat',  status='unknown')
  write (n_unit(unitindex+1,ix,iy),*)'Col. 1 - year'
  write (n_unit(unitindex+1,ix,iy),*)'Col. 2 - month'
  write (n_unit(unitindex+1,ix,iy),*)'Col. 3 - day'
  write (n_unit(unitindex+1,ix,iy),*)'Col. 4 - hour'
  write (n_unit(unitindex+1,ix,iy),*) &
  & 'Col. 5 - the time from the start of integration, days' !& 'Col. 5 - the time from the start of integration, days'
  write (n_unit(unitindex+1,ix,iy),*) &
  & 'Col. 6 - Depth, m'
  write (n_unit(unitindex+1,ix,iy),*) &
  & 'Col. 7 - Temp, m'
  
  count_out(ix,iy) = int(time/(dt_out*hour_sec))
endif

      
if (int(time/(dt_out*hour_sec))>count_out(ix,iy) .or. &
  & abs(time - (count_out(ix,iy) + 1)*dt_out*hour_sec) < small_number) then
  if (nout > 0) then
    allocate (Profile_out(1:nout,1:nvars))
    call LININTERPOL (z_prof,Profile,nprof,z_out,Profile_out(1,1),nout,flag)
    if (nvars >= 2) call LININTERPOL (z_prof,Profile2,nprof,z_out,Profile_out(1,2),nout,flag)
    if (nvars == 3) call LININTERPOL (z_prof,Profile3,nprof,z_out,Profile_out(1,3),nout,flag)
    k = nout
    call DEFFORMAT
    write (unit=n_unit(unitindex,ix,iy),fmt=formats(1)) year,month,day,hour, &
    &                                   time/day_sec,       & !max((hour - floor(hour+small_number))*hour_min,0.d0), & 
    &                                   (z_out(i),Profile_out(i,1:nvars), i = 1, nout)
    do i = 1, nout
      write (unit=n_unit(unitindex+1,ix,iy),fmt=formats(2)) year,month,day,hour, &
      &                                   time/day_sec,       & !max((hour - floor(hour+small_number))*hour_min,0.d0), & 
      &                                   -z_out(i),Profile_out(i,1:nvars)
    enddo
    deallocate (Profile_out)
  else
    allocate (Profile_out(1:nprof,1:nvars))
    Profile_out(1:nprof,1) = Profile(1:nprof)
    if (nvars >= 2) Profile_out(1:nprof,2) = Profile2(1:nprof)
    if (nvars == 3) Profile_out(1:nprof,3) = Profile3(1:nprof)
    k = nprof
    call DEFFORMAT
    write (unit=n_unit(unitindex,ix,iy),fmt=formats(1)) year,month,day,hour, &
    &                                   time/day_sec,       &  ! max((hour - floor(hour+small_number))*hour_min,0.d0), & 
    &                                   (z_prof(i),Profile_out(i,1:nvars), i = 1, nprof)
    do i = 1, nprof
      write (unit=n_unit(unitindex+1,ix,iy),fmt=formats(2)) year,month,day,hour, &
      &                                   time/day_sec,       & ! max((hour - floor(hour+small_number))*hour_min,0.d0), & 
      &                                   -z_prof(i),Profile_out(i,1:nvars)
    enddo
  endif
  if (last) count_out(ix,iy) = count_out(ix,iy) + 1 
endif 
      
if (firstcall(ix,iy) .and. last) firstcall(ix,iy) = .false.

contains
SUBROUTINE DEFFORMAT

implicit none

type charl
  character(len=20) :: ch
  integer(kind=iintegers) :: len
end type charl

type(charl) :: nc1, nc2, nc3, nv, ymdh

ymdh%ch = '3i7,f8.2,' !'3i7,f8.2,' ! year,month,day,hour
ymdh%len = len_trim(ymdh%ch)
write(nc1%ch,'(i3)') k     ; nc1%len = len_trim(nc1%ch)
write(nv%ch, '(i3)') nvars ; nv%len  = len_trim(nv%ch)
if (ndec >= 0) then
  write(nc2%ch,'(i3)') 6 + ndec ; nc2%len = len_trim(nc2%ch)
  write(nc3%ch,'(i3)') ndec     ; nc3%len = len_trim(nc3%ch)
  formats(1) = '('//ymdh%ch(1:ymdh%len)//'f13.4,'//nc1%ch(1:nc1%len)// &
  & '(f10.2, '//nv%ch(1:nv%len)//'f'//nc2%ch(1:nc2%len)//'.'//nc3%ch(1:nc3%len)//'))'
  formats(2) = '('//ymdh%ch(1:ymdh%len)//'f13.4,f10.2,'//nv%ch(1:nv%len)//'f'//nc2%ch(1:nc2%len)//'.'// &
  & nc3%ch(1:nc3%len)//')'
else
  formats(1) = '('//ymdh%ch(1:ymdh%len)//'f13.4,'//nc1%ch(1:nc1%len)//'(f10.2, '//nv%ch(1:nv%len)//'E15.7))'
  formats(2) = '('//ymdh%ch(1:ymdh%len)//'f13.4,f10.2,'//nv%ch(1:nv%len)//'E15.7)'
endif

END SUBROUTINE DEFFORMAT

!100 format (3i7,f8.2,f13.4,<k>(f10.2, f<6+ndec>.<ndec>))
!101 format (3i7,f8.2,f13.4,f10.2,f<6+ndec>.<ndec>)
!102 format (3i7,f8.2,f13.4,<k>(f10.2, E15.7))
!103 format (3i7,f8.2,f13.4,f10.2,E15.7)
END SUBROUTINE PROFILE_OUTPUT


SUBROUTINE ACCUM_VAR &
& (dt, l1, febul, febul0, fdiff, fdiff_lake_surf, qwateroxidtot, & 
& ice_meth_oxid_total, &
& rprod_total_newC, rprod_total_oldC, &
& febultot, febulbottot, fdifftot, fdiff_lake_surftot, &
& rprod_total_newC_integr, rprod_total_oldC_integr, &
& methoxidwat, ice_meth_oxid, add_to_winter)

use METH_OXYG_CONSTANTS, only : &
& molmass_ch4, mfs, mf

implicit none

! Input/output variables

real(kind=ireals), intent(in) :: dt
real(kind=ireals), intent(in) :: l1
real(kind=ireals), intent(in) :: febul, febul0
real(kind=ireals), intent(in) :: fdiff
real(kind=ireals), intent(in) :: fdiff_lake_surf
!real(kind=ireals), intent(in) :: febul(1:2) ! two-meth
!real(kind=ireals), intent(in) :: fdiff(1:2) ! two-meth
!real(kind=ireals), intent(in) :: fdiff_lake_surf(1:2) ! two-meth
real(kind=ireals), intent(in) :: rprod_total_newC
real(kind=ireals), intent(in) :: rprod_total_oldC
real(kind=ireals), intent(in) :: qwateroxidtot
real(kind=ireals), intent(in) :: ice_meth_oxid_total

logical, intent(inout) :: add_to_winter
real(kind=ireals), intent(inout) :: febultot(1:2), febulbottot(1:2)
real(kind=ireals), intent(inout) :: fdifftot(1:2)
real(kind=ireals), intent(inout) :: fdiff_lake_surftot(1:2)
real(kind=ireals), intent(inout) :: rprod_total_newC_integr(1:2)
real(kind=ireals), intent(inout) :: rprod_total_oldC_integr(1:2)
real(kind=ireals), intent(inout) :: methoxidwat(1:2)
real(kind=ireals), intent(inout) :: ice_meth_oxid

!real(kind=ireals), intent(inout) :: febultot(1:2,1:2) ! two-meth
!real(kind=ireals), intent(inout) :: fdifftot(1:2,1:2) ! two-meth
!real(kind=ireals), intent(inout) :: fdiff_lake_surftot(1:2,1:2) ! two-meth
!real(kind=ireals), intent(inout) :: rprod_total_newC_integr(1:2,1:2) ! two-meth
!real(kind=ireals), intent(inout) :: rprod_total_oldC_integr(1:2,1:2) ! two-meth

! Local variables

real(kind=ireals), parameter :: day_sec = 24.*60.*60.

integer(kind=iintegers) :: i ! loop index

if (l1 == 0. .and. (.not. add_to_winter)) then
  i = 1
else
  i = 2
  if (add_to_winter) add_to_winter = .false.
endif

febultot(i) = febultot(i) + febul*mfs*dt/day_sec
febulbottot(i) = febulbottot(i) + febul0*mfs*dt/day_sec
fdifftot(i) = fdifftot(i) - fdiff*mfs*dt/day_sec
fdiff_lake_surftot(i) = fdiff_lake_surftot(i) - fdiff_lake_surf*mfs*dt/day_sec
methoxidwat(i) = methoxidwat(i) - qwateroxidtot*mfs*dt/day_sec
rprod_total_newC_integr(i) = rprod_total_newC_integr(i) + rprod_total_newC*mfs*dt/day_sec 
rprod_total_oldC_integr(i) = rprod_total_oldC_integr(i) + rprod_total_oldC*mfs*dt/day_sec 

ice_meth_oxid = ice_meth_oxid + ice_meth_oxid_total*mf

!do i = 1, 2 ! two-meth
!  if (l1 == 0.) then ! two-meth
!    febultot(1,i) = febultot(1,i) + febul(i)*mfs*dt/day_sec ! two-meth
!    fdifftot(1,i) = fdifftot(1,i) - fdiff(i)*mfs*dt/day_sec ! two-meth
!    fdiff_lake_surftot(1,i) = fdiff_lake_surftot(1,i) - fdiff_lake_surf(i)*mfs*dt/day_sec ! two-meth
!    rprod_total_newC_integr(1,i) = rprod_total_newC_integr(1,i) + rprod_total_newC*mfs*dt/day_sec  ! two-meth
!    rprod_total_oldC_integr(1,i) = rprod_total_oldC_integr(1,i) + rprod_total_oldC*mfs*dt/day_sec  ! two-meth
!  else ! two-meth
!    febultot(2,i) = febultot(2,i) + febul(i)*mfs*dt/day_sec ! two-meth
!    fdifftot(2,i) = fdifftot(2,i) - fdiff(i)*mfs*dt/day_sec ! two-meth
!    fdiff_lake_surftot(2,i) = fdiff_lake_surftot(2,i) - fdiff_lake_surf(i)*mfs*dt/day_sec ! two-meth
!    rprod_total_newC_integr(2,i) = rprod_total_newC_integr(2,i) + rprod_total_newC*mfs*dt/day_sec  ! two-meth
!    rprod_total_oldC_integr(2,i) = rprod_total_oldC_integr(2,i) + rprod_total_oldC*mfs*dt/day_sec  ! two-meth
!  endif ! two-meth
!enddo ! two-meth

END SUBROUTINE ACCUM_VAR    


!> Subroutine TEMPLOC calculates the temperature and other scalars in a specified location
!! of a lake given the vertical profile of horizontally-averaged field 
!! and deviations of heights of constant-density layers
SUBROUTINE TEMPLOC(M,nvars,nvar,itherm,h,dt_out,time,hx1ml,hy1ml,Tw,gsp,gs,bathymwater,rtemp, &
&                  basename,firstcall_)

!Currently implemented for single-lake output only.

use NUMERIC_PARAMS, only : pi
use ARRAYS_GRID, only : gridspacing_type, gridsize_type
use ARRAYS_BATHYM, only : bathym
use INOUT_PARAMETERS, only : lake_misc_unit_min, lake_misc_unit_max
use TIMEVAR, only : hour_sec
use DRIVING_PARAMS, only : grarr1

implicit none

!Input variables
integer(kind=iintegers), intent(in) :: M !> Number of grid layers (number of levels - (M+1) )
integer(kind=iintegers), intent(in) :: nvars !> Number of scalar fields
integer(kind=iintegers), intent(in) :: nvar !> Number of scalar field
integer(kind=iintegers), intent(in) :: itherm(1:M+2) !> Grid levels between constant-density layers

real(kind=ireals), intent(in) :: h !> Lake depth, m
real(kind=ireals), intent(in) :: dt_out !> Output interval, s
real(kind=ireals), intent(in) :: time !> Time, s
real(kind=ireals), intent(in) :: hx1ml(1:M+1), hy1ml(1:M+1) !> Deviations of thickness of 
                                                            !! constant-density layers
real(kind=ireals), intent(in) :: Tw(1:M+1) !> Profile of horizontally-averaged scalar

type(gridspacing_type), intent(in) :: gsp !> Grid spacing group
type(gridsize_type), intent(in) :: gs !> Grid size group
type(bathym), intent(in) :: bathymwater(1:M+1) !> Bathymetry group
type(grarr1), intent(in) :: rtemp !> Coordinates of a location in a lake

character(len=*), intent(in) :: basename
logical, intent(in) :: firstcall_ != .true.

!Output variables
!real(kind=ireals), intent(out) :: Temp !> Temperature in the point set by coordinates rtemp

!Local variables
real(kind=ireals), allocatable :: Temp(:) ! Scalar in the point set by coordinates rtemp
real(kind=ireals) :: amplx, amply, sumh, sumh1, dsumh, larthick, dlarthick, alpha, Lx, Ly, h_ 
real(kind=ireals), save :: tcount = 0.

integer(kind=iintegers) :: i, j, k !Loop indices
integer(kind=iintegers), allocatable, save :: out_unit(:)


character :: coords_point*6

if (.not. allocated(out_unit)) then
  allocate(out_unit(1:nvars)) 
  out_unit = lake_misc_unit_min
endif
if (firstcall_) then
  call CHECK_UNIT(lake_misc_unit_min,lake_misc_unit_max,out_unit(nvar))
  write (coords_point,'(2i3)') gs%ix, gs%iy
  open(out_unit(nvar),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & basename(1:len_trim(basename))//coords_point//'.dat', status='unknown')
endif

if (int(time/(dt_out*hour_sec)) > tcount) then

  allocate (Temp(1:rtemp%numarr))
  Temp(:) = missing_value

  do k = 1, rtemp%numarr
    sumh = 0.
    cyc1 : do i = M+1, 1, -1
      if1 : if (itherm(i+1) > 0) then
        Lx = bathymwater(itherm(2))%Lx!bathymwater(itherm(i+1))%Lx ! Assuming cross-section 
                                                                   ! to be constant with depth
        Ly = bathymwater(itherm(2))%Ly!bathymwater(itherm(i+1))%Ly
        !if (abs(rtemp%arr1(1,k)) > 0.5*Lx .or. abs(rtemp%arr1(2,k)) > 0.5*Ly) exit if1
        amplx = - 0.5*hx1ml(i)*pi !/(Lx*Ly)
        amply = - 0.5*hy1ml(i)*pi !/(Lx*Ly)
        larthick = h*sum(gsp%ddz05(itherm(i):itherm(i+1)-1)) ! Equilibrium thickness of constant-density layer
        dlarthick = amplx*sin(pi*rtemp%arr1(1,k)/Lx) + amply*sin(pi*rtemp%arr1(2,k)/Ly)
        sumh = sumh + larthick + dlarthick
        !print*, 'ay', hy1ml
      endif if1
      if (sumh >= h - rtemp%arr1(3,k)) then
        !print*, itherm(i)+1, amplx
        h_ = h - rtemp%arr1(3,k)
        sumh1 = sumh
        cyclrs : do j = itherm(i)+1, itherm(i+1)
          dsumh = h*gsp%ddz05(j-1)*(1. + dlarthick/larthick) 
          if ( h_ <= sumh1 .and. h_ > sumh1 - 0.5*dsumh ) then
             alpha = (h_ - sumh1 + 0.5*dsumh)/dsumh !Assuming regular grid and the same thickness deviation
                                                    !in adjacent layers
             Temp(k) = Tw(max(j-1,1))*alpha + Tw(j)*(1. - alpha)
             !print*, '1', j, nvar
             exit cyclrs
          elseif ( h_ <= sumh1 - 0.5*dsumh .and. h_ > sumh1 - dsumh ) then
             alpha = (sumh1 - 0.5*dsumh - h_)/dsumh !Assuming regular grid and the same thickness deviation
                                                    !in adjacent layers
             Temp(k) = Tw(min(j+1,M+1))*alpha + Tw(j)*(1. - alpha) 
             !print*, '2', j, nvar
             exit cyclrs
          endif
          sumh1 = sumh1 - dsumh
        enddo cyclrs
        exit cyc1
      endif
    enddo cyc1

  enddo
  
  write(out_unit(nvar),*) time/hour_sec, Temp(1:rtemp%numarr)
  
  if (nvar == nvars) tcount = tcount + 1.

  deallocate(Temp)
endif
!write(*,*) Temp

!if (firstcall) firstcall = .false.
END SUBROUTINE TEMPLOC


END MODULE OUT_MOD

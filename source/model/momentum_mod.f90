MODULE MOMENTUM

use NUMERICS, only : PROGONKA, KRON, MATRIXPROGONKA
use TURB, only : CE_CANUTO, SMOMENT_GALPERIN
use NUMERIC_PARAMS, only : vector_length, pi, small_value

use LAKE_DATATYPES, only : ireals, iintegers

contains
SUBROUTINE MOMENTUM_SOLVER(ix, iy, nx, ny, ndatamax, year, month, day, hour, &
& kor, a_veg, c_veg, h_veg, &
& alphax, alphay, dt, b0, tau_air, tau_i, tau_gr, tau_wav, fetch, depth_area)

! Subroutine MOMENT_SOLVER solves the momentum equations
! for two horizontal components od speed

use ATMOS, only : &
& cdmw2, &
& zref, &
& uwind, vwind, wind, wind10, windwork, &
& u, v, urel, vrel

use DRIVING_PARAMS, only : &
& momflxpart, &
& dyn_pgrad, &
& Turbpar, &
& stabfunc, &
& relwind, &
& M, eos, lindens, &
& cuette, momflux0, &
& PBLpar, &
& tribheat, N_tribin, itribloc, N_tribout, iefflloc, &
& disch_tribin, disch_tribout

use ARRAYS_WATERSTATE, only : &
& Tw1, Sal1

use ARRAYS_BATHYM, only : &
& h1, hx1, hx2, hy1, hy2, &
& hx1t, hx2t, hy1t, hy2t, &
& hx1ml, hx2ml, hy1ml, hy2ml, &
& l1, ls1, &
& dhw, dhw0, &
& area_half, area_int, &
& Lx, Ly, &
& bathymwater, &
& aki, bki

use ARRAYS_GRID, only : &
& ddz, ddz05, dzeta_int, dzeta_05int

use ARRAYS_TURB, only : &
& E1, eps1, k2, &
& row, Ri, &
& knum, veg_sw, u_star, &
& i_maxN, itherm

use ARRAYS_SOIL, only : &
& lsu, lsv

use ARRAYS, only : &
& u1, u2, v1, v2, uv, &
& nstep

use PHYS_CONSTANTS, only: &
& row0, roa0, &
& lamw0, &
& cw, &
& g, &
& cw_m_row0, &
& roughness0, &
& z0_bot, c_bot, niu_wat

use PHYS_FUNC, only : &
& DOMWAVESPEED, &
& LOGFLUX

use TURB_CONST, only :  &
& min_visc, CE0, Cs, Cs1, kar, L0, &
& CL_K_KL_MODEL

use WATER_DENSITY, only: &
& water_dens_t_unesco, &
& water_dens_ts, &
& DDENS_DTEMP0, DDENS_DSAL0, &
& SET_DDENS_DTEMP, SET_DDENS_DSAL

implicit none


! Input variables

! Reals
real(kind=ireals), intent(in) :: dt
real(kind=ireals), intent(in) :: kor
real(kind=ireals), intent(in) :: h_veg, a_veg, c_veg
real(kind=ireals), intent(in) :: alphax, alphay
real(kind=ireals), intent(in) :: hour
real(kind=ireals), intent(in) :: fetch

real(kind=ireals), intent(in) :: depth_area(1:ndatamax,1:2)

! Integers
integer(kind=iintegers), intent(in) :: ix, iy
integer(kind=iintegers), intent(in) :: nx, ny
integer(kind=iintegers), intent(in) :: ndatamax
integer(kind=iintegers), intent(in) :: year, month, day

! Output variables

real(kind=ireals), intent(out) :: b0
real(kind=ireals), intent(out) :: tau_air, tau_gr, tau_i, tau_wav
        
! Local variables     
! Reals

real(kind=ireals) :: CE(M)

real(kind=ireals) :: a(vector_length)
real(kind=ireals) :: b(vector_length)
real(kind=ireals) :: c(vector_length)
real(kind=ireals) :: d(vector_length)

real(kind=ireals) :: am_(vector_length,2,2)
real(kind=ireals) :: bm_(vector_length,2,2)
real(kind=ireals) :: cm_(vector_length,2,2)
real(kind=ireals) :: ym_(vector_length,2)
real(kind=ireals) :: dm_(vector_length,2)

!real(kind=ireals) :: work(1:M+1,1:2) ! Work array

real(kind=ireals) :: row1, row2
real(kind=ireals) :: wr

real(kind=ireals) :: ufr
real(kind=ireals) :: ACC2, ACCk
real(kind=ireals) :: dudz2dvdz2
real(kind=ireals) :: taux, tauy
real(kind=ireals) :: kor2
real(kind=ireals) :: Cz_gr, Cz_i, C10
real(kind=ireals) :: coef1, coef2
real(kind=ireals) :: rp
real(kind=ireals) :: tau_sbl
real(kind=ireals) :: x, y, xx, yy, xx1, yy1, xx2, yy2, xx12, yy12
real(kind=ireals) :: lm
real(kind=ireals) :: ext_lamw
real(kind=ireals) :: month_sec, day_sec, hour_sec
real(kind=ireals) :: AL
real(kind=ireals) :: zup
real(kind=ireals) :: urels, vrels
real(kind=ireals) :: minwind
real(kind=ireals) :: lambda_M, lambda_N
real(kind=ireals) :: lam_force
real(kind=ireals), allocatable :: Force_x(:), Force_y(:)
real(kind=ireals), allocatable :: um(:), vm(:)
real(kind=ireals), allocatable :: rowlars(:)
real(kind=ireals), allocatable :: Amatrix1(:,:), Hlars(:), ulars(:), vlars(:)
real(kind=4), allocatable :: Amatrix(:,:), rhs(:) !to be passed to lapack routines
real(kind=ireals) :: tau_ix, tau_iy
real(kind=ireals) :: Lxi, Lxj, Lyi, Lyj
real(kind=ireals) :: tau_grx, tau_gry
real(kind=ireals), allocatable :: row0_(:)

! Integers
integer(kind=iintegers) :: i, j, k, nlevs, i0 = 0, ierr
integer(kind=iintegers), allocatable :: itherm_new(:)
integer(kind=iintegers), allocatable :: intwork(:)

! Logicals
logical, parameter :: impl_seiches = .true.
logical :: indstab
logical :: ind_bound
logical :: firstcall

! Characters
character :: tp_name*10
character :: numt*1


data firstcall /.true./

! Externals
real(kind=ireals), external :: DZETA
real(kind=4), external :: FindDet

SAVE

firstcall_if : if (firstcall) then
       
  month_sec   = 30.*24.*60.*60.
  day_sec     = 24*60.*60.
  hour_sec    = 60.*60.


  AL    = g/row0

! Parameters of numerical scheme
  ACC2    = 1.d-20 !1.d-20
  ACCk    = 1.d-20 !1.d-20
  knum    = 0.
  minwind = 1.d-2 !1.0d0
  
  lam_force = 4.d-1 !1.d0 !4.d-1
       
! Friction by vegetation

  if (h_veg>0) then
    print*, 'Vegetation friction parameterization &
    & is not operational currently: STOP'
    STOP
  endif

  veg_sw = 0.
  do i = 1, M+1
    if (h1-dzeta_int(i)*h1 < h_veg) veg_sw(i) = 1.
  enddo

  !allocate(row0_(1:M+1)); row0_(:) = row(:)

endif firstcall_if


allocate(Force_x(1:M+1), Force_y(1:M+1))
allocate(um(1:M+1),vm(1:M+1))

u = u1(1)
v = v1(1)


if (relwind%par == 1) then
  urel=uwind-u
  vrel=vwind-v
elseif (relwind%par == 0) then
  urel=uwind
  vrel=vwind
else
  print*, 'relwind=',relwind%par,' - incorrect meaning: STOP' 
  stop
endif

if (abs(urel)<minwind) then
  urels=sign(minwind,urel)
else
  urels=urel
endif

if (abs(vrel)<minwind) then
  vrels=sign(minwind,vrel)
else   
  vrels=vrel
endif

kor2  = 0.5d0*kor*dt 

wr   = sqrt(urels**2+vrels**2)
tau_air = roa0*cdmw2*wr
ufr  = sqrt(tau_air/row0)    

! Sensitivity test
! tau_air = 1.d-15

! PARAMETERIZATION OF WIND STRESS SPLIT-UP INTO TWO PARTS:
! Momentum exchange coefficient for wind speed at 10 m, Wuest and Lorke (2003)
wind10 = wr*log(10./roughness0)/log(zref/roughness0)
!if (wind10 > 5.) then
!  C10 = 1.26*10**(-3.)
!else
!  C10 = 0.0044*wind10**(-1.15)
!endif
C10 = cdmw2*wr/(wind10*wind10)
! a.PARAMETERIZATION USING SIGNIFICANT WAVE HEIGHT (NB1, pp. 8)
! co1=4.*kws/sqrt(C10*1000.*g) !1000 m is "average fetch" of lake Syrdakh
! kw=min(max(co1*wind10,kws),0.75)
! b.PARAMETERIZATION USING WAVE AGE (agw) (NB1, pp. 9)
! agw=0.14*(g*1000.)**(1./3.)*C10**(1./6.)*wind10**(-2./3.) 
! kw=min(max(kws*1.15/agw,kws),0.75)
! KW=0.75

!Partition of momentum flux between wave development and currents
!followng Lin et al. (2002, J. Phys. Ocean.)
tau_wav = momflxpart%par * &
& roa0*C10*(wind10 - DOMWAVESPEED(sqrt(tau_air/roa0),fetch)/1.2)**2 ! 1.2 - the ratio of dominant wave speed to mean wave speed
tau_sbl = max(0._ireals,tau_air - tau_wav)

u_star = sqrt(tau_sbl/row0)

! tau_wav = tau_air !*kw
! tau_sbl = tau_air !*(1.-kw)

! Calculation of Shezy coefficient
if (ls1 > 0) then
  rp = 0.01 ! height of roughness elements on ice bottom
else
  rp = 0.1 ! height of roughness elements on the ground bottom
endif
if (l1 > 0) then
  Cz_gr = 7.7*sqrt(g)*(0.5*h1/rp)**(1./6.) !50. !ground Shezy coefficient
else
  Cz_gr = 7.7*sqrt(g)*(h1/rp)**(1./6.)
endif

if (l1 > 0.) then
  rp = 0.01 
  Cz_i = 7.7*sqrt(g)*(0.5*h1/rp)**(1./6.)
endif

! Friction at the bottom
 tau_gr = row0*g*Cz_gr**(-2)*(u1(M+1)**2+v1(M+1)**2)
 tau_grx = 0. !row0*g*Cz_gr**(-2)*u1(M+1)*sqrt((u1(M+1)**2+v1(M+1)**2))
 tau_gry = 0. !row0*g*Cz_gr**(-2)*v1(M+1)*sqrt((u1(M+1)**2+v1(M+1)**2))
 
 coef2 = g*Cz_gr**(-2)*sqrt(u1(M+1)**2+v1(M+1)**2) !*0.005 !*row0

! Friction at the water - upper ice interface
 if (l1 > 0.) then
   tau_i = row0*g*Cz_i**(-2)*(u1(1)**2+v1(1)**2)
   tau_ix = 0. !- row0*g*Cz_i**(-2)*u1(1)*sqrt(u1(1)**2+v1(1)**2)
   tau_iy = 0. !- row0*g*Cz_i**(-2)*v1(1)*sqrt(u1(1)**2+v1(1)**2)
   
   coef1 = -g*Cz_i**(-2)*sqrt(u1(1)**2+v1(1)**2)
 endif
 
! Friction at the lateral boundaries
if (depth_area(1,2) >= 0.) then
  do i = 2, M

    !xx = c_bot !linear drag
    xx = c_bot*sqrt(u1(i)**2 + v1(i)**2) !quadratic law

    !xx = LOGFLUX(sqrt(u1(i)*u1(i) + v1(i)*v1(i)), - u1(i), &
    !& 0.5*ddz(i-1)*h1, z0_bot, z0_bot, 1._ireals, 2_iintegers)
    !lsu%water(i) = - 1./bathymwater(i)%area_int * &
    !& (bathymwater(i)%area_half - bathymwater(i-1)%area_half)/(ddz05(i-1)*h1)*xx

    !xx = LOGFLUX(sqrt(u1(i)*u1(i) + v1(i)*v1(i)), - v1(i), &
    !& 0.5*ddz(i-1)*h1, z0_bot, z0_bot, 1._ireals, 2_iintegers)
    lsv%water(i) = - 1./bathymwater(i)%area_int * &
    & (bathymwater(i)%area_half - bathymwater(i-1)%area_half)/(ddz05(i-1)*h1)*xx
  enddo

  !xx = c_bot !linear drag
  xx = c_bot*sqrt(u1(1)**2 + v1(1)**2) !quadratic law

  ! Top layer
  !xx = LOGFLUX(sqrt(u1(1)*u1(1) + v1(1)*v1(1)), - u1(1), &
  !& 0.5*ddz(1)*h1, z0_bot, z0_bot, 1._ireals, 2_iintegers)
  lsu%water(1) = - 2./bathymwater(1)%area_int * &
  & (bathymwater(1)%area_half - bathymwater(1)%area_int)/(ddz(1)*h1)*xx
  !xx = LOGFLUX(sqrt(u1(1)*u1(1) + v1(1)*v1(1)), - v1(1), &
  !& 0.5*ddz(1)*h1, z0_bot, z0_bot, 1._ireals, 2_iintegers)
  lsv%water(i) = - 2./bathymwater(1)%area_int * &
  & (bathymwater(1)%area_half - bathymwater(1)%area_int)/(ddz(1)*h1)*xx 
  lsu%water(M+1) = 0.
  lsv%water(M+1) = 0.
endif
!print*, xx
      
! Eddy viscosity parameterization
do i=1,M+1
  uv(i)=sqrt(u1(i)**2+v1(i)**2)
enddo
do i=1,M
! Ri(i)=min(max(g/row0/(((uv(i+1)-uv(i))/
! & (ddz*h1))**2)*(row(i+1)-row(i))/(ddz*h1),-10.d0),10.d0)
  dudz2dvdz2 = max( ( (u1(i+1)-u1(i))/(ddz(i)*h1) )**2 + &
  & ( (v1(i+1)-v1(i))/(ddz(i)*h1) )**2, small_value)
  Ri(i) = g/row0*(row(i+1)-row(i))/(ddz(i)*h1)/dudz2dvdz2
enddo

SELECT CASE (Turbpar%par)
       
! 1. "Empirical" parametrization: Stepanenko, Lykossov (2005)
  CASE (1)
    tp_name='SL_Empir' 
    k2(1) =(lamw0*10.+(wind10*2./zref/20.)*(lamw0*1000.-lamw0*10.))/ &
    & (cw_m_row0)
    ext_lamw = log(k2(1)*cw_m_row0/(lamw0*10.))/h1
    do i=2,M+1
      k2(i) = k2(1)*exp(-ext_lamw*dzeta_05int(i)*h1)+niu_wat
    enddo
    k2(1)=k2(1)+niu_wat

! 2. "E-epsilon" parameterization: k=CE*E**2/eps with 
!    prognostic equations for E and eps
  CASE (2)
    do i=1,M+1
      if (E1(i)<=0) E1(i)=10.**(-16)
      if (eps1(i)<=0) eps1(i)=10.**(-18)
    enddo
    do i = 1, M
      if (stabfunc%par == 1) then
        CE(i) = CE0
      else
        lambda_N = E1(i)*E1(i)/((eps1(i) + ACC2)*(eps1(i) + ACC2))* &
        & ( SET_DDENS_DTEMP(eos%par,lindens%par,Tw1(i),Sal1(i))*(Tw1(i+1) - Tw1(i) ) + &
        &   SET_DDENS_DSAL(eos%par,lindens%par,Tw1(i),Sal1(i))*(Sal1(i+1) - Sal1(i) ) ) / &
        & (h1*ddz(i)) *AL
        lambda_M = E1(i)*E1(i)/((eps1(i)+ACC2)*(eps1(i)+ACC2))* &
        & ( (u1(i+1)-u1(i))*(u1(i+1)-u1(i)) + &
        &   (v1(i+1)-v1(i))*(v1(i+1)-v1(i)) ) / &
        &    (h1*h1*ddz(i)*ddz(i))
        if (stabfunc%par == 2) then
          CE(i)  = CE_CANUTO (lambda_M, lambda_N)
        elseif (stabfunc%par == 3)  then
          CE(i)  = sqrt(2.d0)*CL_K_KL_MODEL*SMOMENT_GALPERIN (lambda_N)
        endif
      endif
    enddo
    do i = 1, M
      k2(i) = max(CE(i)*E1(i)**2/(eps1(i) + ACCk),min_visc) + niu_wat + knum(i)
    end do

! 3. Nickuradze (NICK) formulation: Rodi (1993)
  CASE (3)
    tp_name='NICK'
    do i=1,M
      zup=h1-dzeta_05int(i)*h1
      lm=h1*(0.14-0.08*(1-zup/h1)**2-0.06*(1-zup/h1)**4)
      k2(i)=lm**2*abs((uv(i+1)-uv(i))/(ddz(i)*h1))*exp(-Cs*Ri(i)) &
      & + niu_wat
    enddo
      
! 4. Parabolic (PARAB) formulation: Engelund (1976)
  CASE (4)
    tp_name='PARAB'
    do i=1,M
      zup=h1-dzeta_05int(i)*h1
      k2(i)=kar*ufr*zup*(1-zup/h1)*exp(-Cs*Ri(i)) &
      & +niu_wat
    enddo 

! 5. W2 (used in Version 2 of CE-QUAL-W2 model): Cole and Buchak (1995)   
  CASE (5)
    print*, 'Turbpar = 5 is not operational setting: STOP'
    STOP

! 6. W2N (W2 with mixing length of Nickuradze): Cole and Buchak (1995) and Rodi (1993)
  CASE (6)
    print*, 'Turbpar = 6 is not operational setting: STOP'
    STOP
   
! 7. RNG (re-normalization group) formulation: Simoes (1998)
  CASE (7)
    tp_name='RNG'
    do i=1,M
      zup=h1-dzeta_05int(i)*h1
      k2(i)=niu_wat*(1+max(3*kar*(zup/niu_wat*ufr)**3* &
      & (1-zup/h1)**3-Cs1,0.d0))**(1./3.)*exp(-Cs*Ri(i))+niu_wat  
    enddo
END SELECT
        
! Velocity components at the previous or intermediate timestep
! (the latter is in case of using splitting-up scheme for momentum equations) 
um(:) = u1(:)
vm(:) = v1(:)

bctop : if (cuette%par == 0 .or. cuette%par == 11) then

  ! General case of momentum b.c.s
  
  ! Equations for horizontal velocities (u and v)
  
  k2(1) = area_half(1)/area_int(1)*k2(1)
  
  xx = 0.5*(1. - dhw0*h1*ddz(1)/(k2(1)*dt))
  yy = (ddz(1)*h1)**2/(k2(1)*dt)
  
  ! The top boundary conditions for u and v
  ! 1-st case: water surface is free from ice:
  ! momentum flux = momentum flux from atmosphere

  if (l1 == 0.) then
  
    ! Constant momentum flux is imposed at the top boundary if Cuette == 11
    if (cuette%par == 11 .or. PBLpar%par == 0) then
      taux = momflux0%par/row0
      tauy = 0.
    else
      taux = urels/wr*tau_sbl/row0
      tauy = vrels/wr*tau_sbl/row0
    endif
    
    ! Seiche models
    seiche_model : if (dyn_pgrad%par == 1) then
      !Force_x = lam_force * &
      !& 1./h1*(tau_grx - taux - dhw/dt*u1(M+1) + dhw0/dt*(u1(M+1)-u1(1)))
      !Force_y = lam_force * &
      !& 1./h1*(tau_gry - tauy - dhw/dt*v1(M+1) + dhw0/dt*(v1(M+1)-v1(1)))
      call DHDXDHDY(M,u1,v1,area_int,Lx,Ly,ddz05,h1,dt,hx1,hx2,hy1,hy2)
      Force_x(1:M+1) = - g * pi * 0.5 * (hx2 - hx1)/Lx(1)
      Force_y(1:M+1) = - g * pi * 0.5 * (hy2 - hy1)/Ly(1)
    elseif (dyn_pgrad%par == 2 .and. (i_maxN > 1 .and. i_maxN < M+1) ) then
      !Two-layer approximation for pressure gradient
      !if (i0 == 0) i0 = i_maxN
      call DHDXDHDY2(i_maxN,M,u1,v1,area_int,area_half,Lx,Ly,ddz05,h1,dt, &
      & hx1,hx2,hy1,hy2,hx1t,hx2t,hy1t,hy2t)

      Force_x(1:i_maxN-1) = - g * pi * 0.5 * ( (hx2 - hx1)/Lx(1) + &
      & (hx2t - hx1t)/bathymwater(i_maxN-1)%Lx_half )
      Force_y(1:i_maxN-1) = - g * pi * 0.5 * ( (hy2 - hy1)/Ly(1) + &
      & (hy2t - hy1t)/bathymwater(i_maxN-1)%Ly_half )

      ! Variable frid spacing to be included
      row1 = sum(row(1:i_maxN-1))/real(i_maxN-1)
      row2 = sum(row(i_maxN:M+1))/real(M+2-i_maxN)

      Force_x(i_maxN:M+1) = - g * pi * 0.5 * ( row1/(row2*Lx(1))*(hx2 - hx1) + &
      & (hx2t - hx1t)/bathymwater(i_maxN-1)%Lx_half )
      Force_y(i_maxN:M+1) = - g * pi * 0.5 * ( row1/(row2*Ly(1))*(hy2 - hy1) + &
      & (hy2t - hy1t)/bathymwater(i_maxN-1)%Ly_half )

      !print*,hx1,hx2,hy1,hy2,hx1t,hx2t,hy1t,hy2t
      !read*
    elseif (dyn_pgrad%par == 3) then
      !Multi-layer approximation for pressure gradient
      call DHDXDHDYML(M,u1,v1,area_int,area_half,Lx,Ly,ddz05,h1,dt,hx1ml,hx2ml,hy1ml,hy2ml)
      !call DENSLAYERS(M,bathymwater,ddz05,row,h1,itherm)
      !x-pressure gradient
      xx1 = row(1)*(hx2ml(1) - hx1ml(1))/Lx(1)
      yy1 = 0.
      do i = 2, M+1
        yy1 = yy1 + (hx2ml(i) - hx1ml(i))/bathymwater(i-1)%Lx_half
      enddo
      do i = 1, M+1
        Force_x(i) = - 0.5/row(i)*pi*g*(xx1 + yy1*row(i))
        if (i /= M+1) then
          xx1 = xx1 + row(i+1)*(hx2ml(i+1) - hx1ml(i+1))/bathymwater(i)%Lx_half
          yy1 = yy1 - (hx2ml(i+1) - hx1ml(i+1))/bathymwater(i)%Lx_half
        endif
      enddo
      !y-pressure gradient
      xx1 = row(1)*(hy2ml(1) - hy1ml(1))/Ly(1)
      yy1 = 0.
      do i = 2, M+1
        yy1 = yy1 + (hy2ml(i) - hy1ml(i))/bathymwater(i-1)%Ly_half
      enddo
      do i = 1, M+1
        Force_y(i) = - 0.5/row(i)*pi*g*(xx1 + yy1*row(i))
        if (i /= M+1) then
          xx1 = xx1 + row(i+1)*(hy2ml(i+1) - hy1ml(i+1))/bathymwater(i)%Ly_half
          yy1 = yy1 - (hy2ml(i+1) - hy1ml(i+1))/bathymwater(i)%Ly_half
        endif
      enddo
      !print*, Force_x
      !print*, Force_y
      if (i_maxN > 1 .and. i_maxN < M+1) then
        hx1 = sum(hx1ml(1:i_maxN-1))
        hx2 = sum(hx2ml(1:i_maxN-1))
        hy1 = sum(hy1ml(1:i_maxN-1))
        hy2 = sum(hy2ml(1:i_maxN-1))
        hx1t = sum(hx1ml(i_maxN:M+1))
        hx2t = sum(hx2ml(i_maxN:M+1))
        hy1t = sum(hy1ml(i_maxN:M+1))
        hy2t = sum(hy2ml(i_maxN:M+1))
      else
        hx1 = sum(hx1ml(1:M+1))
        hx2 = sum(hx2ml(1:M+1))
        hy1 = sum(hy1ml(1:M+1))
        hy2 = sum(hy2ml(1:M+1))
        hx1t = 0.
        hx2t = 0.
        hy1t = 0.
        hy2t = 0.
      endif
    elseif (dyn_pgrad%par == 4) then
      !Multi-layer approximation for pressure gradient
      !if (.not. allocated(rowlars)) allocate (rowlars(1:M+1))
      !if (firstcall) call DENSLAYERS(M,i_maxN,bathymwater,ddz05,row,h1,nlevs,itherm,rowlars)
      !if (i0 == 0) i0 = i_maxN
      allocate (rowlars(1:M+1))
      allocate (itherm_new(1:M+2))
      itherm_new = itherm
      call DENSLAYERS(M,i_maxN,nstep,bathymwater,ddz05,ddz,row,h1,nlevs,itherm_new,rowlars)
      !print*, nlevs, itherm, rowlars

      if ( .not. (all(itherm_new == itherm)) .and. nstep > 1) then
        ! The layers of homogeneous density change at this timestep:
        call RELAYER(M,itherm,itherm_new,hx1ml,hx2ml,hy1ml,hy2ml)
        !print*, itherm, itherm_new
        !read*
      endif
      itherm = itherm_new

      seiche_scheme : if (nlevs == 1) then
        !Single layer (barotropic case), implicit scheme
        xx1 = 0.; xx2 = 0.
        yy1 = 0.; yy2 = 0.
        do j = 1, M+1
          xx1 = xx1 + ddz05(j-1)*Ly(j)*u1(j)
          xx2 = xx2 + ddz05(j-1)*Ly(j)
          yy1 = yy1 + ddz05(j-1)*Lx(j)*v1(j)
          yy2 = yy2 + ddz05(j-1)*Lx(j)
        enddo
        xx1 = xx1/xx2 !Mean x-velocity 
        yy1 = yy1/yy2 !Mean y-velocity
        x = 1. + 0.25*dt*dt*pi*pi*g*h1/(bathymwater(1)%Lx*bathymwater(1)%Lx)
        y = 1. + 0.25*dt*dt*pi*pi*g*h1/(bathymwater(1)%Ly*bathymwater(1)%Ly)
        xx12 = ( xx1 - 0.25*dt*pi*g/bathymwater(1)%Lx* &
        & (2.*(hx2ml(1)-hx1ml(1)) + dt*pi*h1*xx1/bathymwater(1)%Lx) )/x !Updated x-velocity
        yy12 = ( yy1 - 0.25*dt*pi*g/bathymwater(1)%Ly* &
        & (2.*(hy2ml(1)-hy1ml(1)) + dt*pi*h1*yy1/bathymwater(1)%Ly) )/y !Updated y-velocity
        hx2ml(1) = hx2ml(1) + 0.5*dt*pi*h1/bathymwater(1)%Lx*(xx1 + xx12)
        hx1ml(1) = hx1ml(1) - 0.5*dt*pi*h1/bathymwater(1)%Lx*(xx1 + xx12)
        hy2ml(1) = hy2ml(1) + 0.5*dt*pi*h1/bathymwater(1)%Ly*(yy1 + yy12)
        hy1ml(1) = hy1ml(1) - 0.5*dt*pi*h1/bathymwater(1)%Ly*(yy1 + yy12)
        um(:) = u1(:) + xx12 - xx1
        vm(:) = v1(:) + yy12 - yy1
        Force_x(:) = 0.
        Force_y(:) = 0. 

        !call DHDXDHDY(M,u1,v1,area_int,Lx,Ly,ddz05,h1,dt,hx1ml(1),hx2ml(1),hy1ml(1),hy2ml(1))
        !um(:) = u1(:) - dt * g * pi * 0.5 * (hx2ml(1) - hx1ml(1))/Lx(1)
        !vm(:) = v1(:) - dt * g * pi * 0.5 * (hy2ml(1) - hy1ml(1))/Ly(1)
        !print*, xx-xx1, yy-yy1, hx2ml(1), hy2ml(1)
        !read*
        !um(:) = u1(:)
        !vm(:) = v1(:)
      elseif (.not. impl_seiches) then
        ! Explicit scheme for constant-density layers' thicknesses and 
        ! tendency of speed due to horizontal pressure gradient      
        call DHDXDHDY3(itherm,M,u1,v1,area_int,area_half,Lx,Ly,ddz05,h1,dt, &
        & hx1ml,hx2ml,hy1ml,hy2ml)
        !x-pressure gradient
        xx1 = rowlars(1)*(hx2ml(1) - hx1ml(1))/Lx(1)
        yy1 = 0.
        do i = 2, nlevs
          yy1 = yy1 + (hx2ml(i) - hx1ml(i))/bathymwater(itherm(i))%Lx_half
        enddo
        do i = 1, nlevs
          xx2 = - 0.5/rowlars(i)*pi*g*(xx1 + yy1*rowlars(i))
          forall (j = itherm(i)+1:itherm(i+1)) Force_x(j) = xx2
          if (i /= nlevs) then
            xx1 = xx1 + rowlars(i+1)*(hx2ml(i+1) - hx1ml(i+1))/bathymwater(itherm(i+1))%Lx_half
            yy1 = yy1 - (hx2ml(i+1) - hx1ml(i+1))/bathymwater(itherm(i+1))%Lx_half
          endif
        enddo
        !y-pressure gradient
        xx1 = rowlars(1)*(hy2ml(1) - hy1ml(1))/Ly(1)
        yy1 = 0.
        do i = 2, nlevs
          yy1 = yy1 + (hy2ml(i) - hy1ml(i))/bathymwater(itherm(i))%Ly_half
        enddo
        do i = 1, nlevs
          yy2 = - 0.5/rowlars(i)*pi*g*(xx1 + yy1*rowlars(i))
          forall (j = itherm(i)+1:itherm(i+1)) Force_y(j) = yy2
          if (i /= nlevs) then
            xx1 = xx1 + rowlars(i+1)*(hy2ml(i+1) - hy1ml(i+1))/bathymwater(itherm(i+1))%Ly_half
            yy1 = yy1 - (hy2ml(i+1) - hy1ml(i+1))/bathymwater(itherm(i+1))%Ly_half
          endif
        enddo
        ! End of explicit scheme
      elseif (impl_seiches) then
        !Crank-Nicolson scheme for constant-density layers' thicknesses and 
        !tendency of speed due to horizontal pressure gradient
        allocate(Hlars(1:nlevs), ulars(1:nlevs), vlars(1:nlevs))
        allocate(Amatrix(1:nlevs,1:nlevs), rhs(1:nlevs))
        forall (i = 1:nlevs) Hlars(i) = sum(ddz05(itherm(i):itherm(i+1)-1))*h1
        ! Mean velocities over the layers of constant density
        do i = 1, nlevs
          xx1 = 0.; xx2 = 0.
          yy1 = 0.; yy2 = 0.
          do j = itherm(i)+1, itherm(i+1)
            !xx1 = xx1 + ddz05(j-1)*Ly(j)*u1(j)
            !xx2 = xx2 + ddz05(j-1)*Ly(j)
            !yy1 = yy1 + ddz05(j-1)*Lx(j)*v1(j)
            !yy2 = yy2 + ddz05(j-1)*Lx(j)
            xx1 = xx1 + ddz05(j-1)*bathymwater(j)%area_int*u1(j)
            xx2 = xx2 + ddz05(j-1)*bathymwater(j)%area_int
            yy1 = yy1 + ddz05(j-1)*bathymwater(j)%area_int*v1(j)
            yy2 = yy2 + ddz05(j-1)*bathymwater(j)%area_int
          enddo
          ulars(i) = xx1/xx2
          vlars(i) = yy1/yy2
        enddo
        ! Solve for u, hx2ml, hx1ml
        do i = 1, nlevs
          if (i == 1) then
            Lxi = bathymwater(1_iintegers)%Lx
          else
            Lxi = bathymwater(itherm(i))%Lx_half
          endif
          xx1 = 0.25*g*(pi*dt)**2/(rowlars(i)*Lxi)
          xx2 = 0.25*g*(pi*dt)   /(rowlars(i)*Lxi)
          yy1 = dt*pi
          j = 1
          Amatrix(i,j) = xx1*rowlars(min(i,j))*aki(j,i)*Hlars(j)/bathymwater(1_iintegers)%Lx
          do j = 2, nlevs
            Amatrix(i,j) = xx1*rowlars(min(i,j))*aki(j,i)*Hlars(j)/bathymwater(itherm(j))%Lx_half
          enddo
          Amatrix(i,i) = Amatrix(i,i) + 1.
          rhs(i) = ulars(i)
          j = 1
          rhs(i) = rhs(i) - xx2*rowlars(min(i,j))*aki(j,i)*(2.*(hx2ml(j) - hx1ml(j)) + yy1*Hlars(j)*ulars(j)/ &
          & bathymwater(1_iintegers)%Lx)
          do j = 2, nlevs
            rhs(i) = rhs(i) - xx2*rowlars(min(i,j))*aki(j,i)*(2.*(hx2ml(j) - hx1ml(j)) + yy1*Hlars(j)*ulars(j)/ &
            & bathymwater(itherm(j))%Lx_half) 
          enddo
        enddo
        ! Solving linear equations for u
        allocate(intwork(1:nlevs))
        !print*, FindDet(Amatrix, nlevs)
        call SGESV( nlevs, 1_4, Amatrix, nlevs, intwork, rhs, nlevs, ierr )
        ! rhs now stores a solution (updated ulars)
        ! Updating hx2ml, hx1ml and u
        i = 1
        xx1 = 0.5*dt*pi*Hlars(i)/bathymwater(1_iintegers)%Lx*(rhs(i) + ulars(i))
        hx2ml(i) = hx2ml(i) + xx1
        hx1ml(i) = hx1ml(i) - xx1
        forall ( j = itherm(i)+1:itherm(i+1) ) um(j) = u1(j) + rhs(i) - ulars(i)
        do i = 2, nlevs
          xx1 = 0.5*dt*pi*Hlars(i)/bathymwater(itherm(i))%Lx_half*(rhs(i) + ulars(i))
          hx2ml(i) = hx2ml(i) + xx1
          hx1ml(i) = hx1ml(i) - xx1
          forall ( j = itherm(i)+1:itherm(i+1) ) um(j) = u1(j) + rhs(i) - ulars(i)
        enddo

        ! Solve for v, hy2ml, hy1ml
        do i = 1, nlevs
          if (i == 1) then
            Lyi = bathymwater(1_iintegers)%Ly
          else
            Lyi = bathymwater(itherm(i))%Ly_half
          endif
          xx1 = 0.25*g*(pi*dt)**2/(rowlars(i)*Lyi)
          xx2 = 0.25*g*(pi*dt)   /(rowlars(i)*Lyi)
          yy1 = dt*pi
          j = 1
          Amatrix(i,j) = xx1*rowlars(min(i,j))*bki(j,i)*Hlars(j)/bathymwater(1_iintegers)%Ly
          do j = 2, nlevs
            Amatrix(i,j) = xx1*rowlars(min(i,j))*bki(j,i)*Hlars(j)/bathymwater(itherm(j))%Ly_half
          enddo
          Amatrix(i,i) = Amatrix(i,i) + 1.
          rhs(i) = vlars(i)
          j = 1
          rhs(i) = rhs(i) - xx2*rowlars(min(i,j))*bki(j,i)*(2.*(hy2ml(j) - hy1ml(j)) + yy1*Hlars(j)*vlars(j)/ &
          & bathymwater(1_iintegers)%Ly)
          do j = 2, nlevs
            rhs(i) = rhs(i) - xx2*rowlars(min(i,j))*bki(j,i)*(2.*(hy2ml(j) - hy1ml(j)) + yy1*Hlars(j)*vlars(j)/ &
            & bathymwater(itherm(j))%Ly_half)
          enddo
        enddo
        ! Solving linear equations for v
        call SGESV( nlevs, 1_4, Amatrix, nlevs, intwork, rhs, nlevs, ierr )
        ! rhs now stores a solution (updated vlars)
        ! Updating hy2ml, hy1ml and v
        i = 1
        yy1 = 0.5*dt*pi*Hlars(i)/bathymwater(1_iintegers)%Ly*(rhs(i) + vlars(i))
        hy2ml(i) = hy2ml(i) + yy1
        hy1ml(i) = hy1ml(i) - yy1
        forall ( j = itherm(i)+1:itherm(i+1) ) vm(j) = v1(j) + rhs(i) - vlars(i)
        do i = 2, nlevs
          yy1 = 0.5*dt*pi*Hlars(i)/bathymwater(itherm(i))%Ly_half*(rhs(i) + vlars(i))
          hy2ml(i) = hy2ml(i) + yy1
          hy1ml(i) = hy1ml(i) - yy1
          forall ( j = itherm(i)+1:itherm(i+1) ) vm(j) = v1(j) + rhs(i) - vlars(i)
        enddo

        deallocate(Hlars, ulars, vlars)
        deallocate(Amatrix, rhs)
        deallocate(intwork)

        Force_x(:) = 0.d0
        Force_y(:) = 0.d0

      endif seiche_scheme

      ! Adding tributaries inflow to the top constant-density layer
      if (N_tribin%par > 0 .and. tribheat%par > 0) then
        do j = 1, nlevs
          if (j == 1) then
            xx2 = bathymwater(1)%area_int
          else
            xx2 = bathymwater(itherm(j))%area_half
          endif
          do i = 1, N_tribin%par
            ! itribloc: location of tributary on a lake:
            ! Y
            ! -----------------------
            ! |    1     |     2    |
            ! |          |          |
            ! -----------------------
            ! |          |          |
            ! |    4     |     3    |
            ! ----------------------- X
            xx1 = 0.              
            do k = itherm(j)+1, itherm(j+1)
              xx1 = xx1 + 2.*dt*disch_tribin(i,k)/xx2
            enddo
            !print*, 'trib', disch_tribin(i)
            select case (itribloc(i))
              case(1)
                hx1ml(j) = hx1ml(j) + xx1
                hy2ml(j) = hy2ml(j) + xx1
              case(2)
                hx2ml(j) = hx2ml(j) + xx1
                hy2ml(j) = hy2ml(j) + xx1
              case(3)
                hx2ml(j) = hx2ml(j) + xx1
                hy1ml(j) = hy2ml(j) + xx1
              case(4)
                hx1ml(j) = hx1ml(j) + xx1
                hy1ml(j) = hy1ml(j) + xx1
            end select
          enddo
        enddo
      endif
  
      ! Subtracting effluent outflow from the top constant-density layer
      if (N_tribout > 0 .and. tribheat%par > 0) then
        do j = 1, nlevs
          if (j == 1) then
            xx2 = bathymwater(1)%area_int
          else
            xx2 = bathymwater(itherm(j))%area_half
          endif
          ! Single outflow is assumed
          ! iefflloc: location of tributary on a lake:
          ! Y
          ! -----------------------
          ! |    1     |     2    |
          ! |          |          |
          ! -----------------------
          ! |          |          |
          ! |    4     |     3    |
          ! ----------------------- X
          xx1 = 0.
          do k = itherm(j)+1, itherm(j+1)
            xx1 = xx1 + 2.*dt*disch_tribout(1,k)/xx2
          enddo
          !print*, 'effl', disch_tribout
          select case (iefflloc%par)
            case(1)
              hx1ml(j) = hx1ml(j) - xx1
              hy2ml(j) = hy2ml(j) - xx1
            case(2)
              hx2ml(j) = hx2ml(j) - xx1
              hy2ml(j) = hy2ml(j) - xx1
            case(3)
              hx2ml(j) = hx2ml(j) - xx1
              hy1ml(j) = hy2ml(j) - xx1
            case(4)
              hx1ml(j) = hx1ml(j) - xx1
              hy1ml(j) = hy1ml(j) - xx1
          end select
        enddo
      endif

      !print*, Force_x
      !print*, Force_y

      if (i_maxN > 1 .and. i_maxN < M+1) then
        hx1 = 0. ; hx2 = 0. ; hy1 = 0. ; hy2 = 0. 
        hx1t = 0.; hx2t = 0.; hy1t = 0.; hy2t = 0. 

        hx1 = hx1ml(1)
        hx2 = hx2ml(1)
        hy1 = hy1ml(1)
        hy2 = hy2ml(1)
        do i = 2, nlevs
          hx1t = hx1t + hx1ml(i)
          hx2t = hx2t + hx2ml(i)
          hy1t = hy1t + hy1ml(i)
          hy2t = hy2t + hy2ml(i)
        enddo

        !do i = 1, nlevs
        !  if (itherm(i+1) < i_maxN) then
        !    hx1 = hx1 + hx1ml(i)
        !    hx2 = hx2 + hx2ml(i)
        !    hy1 = hy1 + hy1ml(i)
        !    hy2 = hy2 + hy2ml(i)
        !  elseif (itherm(i) + 1 >= i_maxN) then
        !    hx1t = hx1t + hx1ml(i)
        !    hx2t = hx2t + hx2ml(i)
        !    hy1t = hy1t + hy1ml(i)
        !    hy2t = hy2t + hy2ml(i)
        !  endif
        !enddo
      else
        hx1 = sum(hx1ml(1:M+1))
        hx2 = sum(hx2ml(1:M+1))
        hy1 = sum(hy1ml(1:M+1))
        hy2 = sum(hy2ml(1:M+1))
        hx1t = 0.
        hx2t = 0.
        hy1t = 0.
        hy2t = 0.
      endif
      deallocate(rowlars)
      deallocate(itherm_new)
    else
      Force_x(:) = 0.d0
      Force_y(:) = 0.d0
    endif seiche_model
    
    do i = 1, 2
      cm_(1,i,i) = (xx + yy + 0.5*yy*dt*a_veg*c_veg * &
      & sqrt(um(1)**2 + vm(1)**2)*veg_sw(1))
      bm_(1,i,i) = xx
    enddo
    ! Adding friction at lateral boundaries
    cm_(1,1,1) = cm_(1,1,1) - yy*dt*lsu%water(1) !0.5*
    cm_(1,2,2) = cm_(1,2,2) - yy*dt*lsv%water(1) !0.5*
  
    bm_(1,1,2) = 0.
    bm_(1,2,1) = 0.
    cm_(1,1,2) = -kor*(h1*ddz(1))**2/(2*k2(1))
    cm_(1,2,1) = kor*(h1*ddz(1))**2/(2*k2(1))
  
    dm_(1,1)=(taux*ddz(1)*h1/k2(1)+yy*um(1)) + &
    & (um(2)-um(1))*xx + &
    & kor*(h1*ddz(1))**2*vm(1)/(2*k2(1)) + &
    & yy*dt*g*tan(pi*alphax/180.) + &
    & yy*dt*(Force_x(1)) ! + 0.5*lsu%water(1)*um(1)
    dm_(1,2)=(tauy*ddz(1)*h1/k2(1)+yy*vm(1)) + &
    & (vm(2)-vm(1))*xx - &
    & kor*(h1*ddz(1))**2*um(1)/(2*k2(1)) + &
    & yy*dt*g*tan(pi*alphay/180.) + &
    & yy*dt*(Force_y(1)) ! + 0.5*lsv%water(1)*v1(1)
  endif
  
  ! 2-nd case: thin ice:
  ! momentum flux = 
  ! weighted momentum flux from atmosphere plus weighted friction 
  ! at the water-ice interface
  
  if (l1 > 0 .and. l1 <= L0) then 
  
    b0 = l1/L0
    taux = urels/wr*tau_sbl/row0 !wdir+pi
    tauy = vrels/wr*tau_sbl/row0 !wdir+pi
    
    !if (dyn_pgrad%par > 0) then
    !  Force_x(1:M+1) = lam_force * 1./h1*(tau_grx - (1.-b0)*taux - b0*tau_ix - &
    !  & dhw/dt*u1(M+1) + dhw0/dt*(u1(M+1)-u1(1)))
    !  Force_y(1:M+1) = lam_force * 1./h1*(tau_gry - (1.-b0)*tauy - b0*tau_iy - &
    !  & dhw/dt*v1(M+1) + dhw0/dt*(v1(M+1)-v1(1)))
    !else
      Force_x(:) = 0.d0
      Force_y(:) = 0.d0
    !endif
    
    do i = 1, 2
      cm_(1,i,i) = (xx + yy - b0*coef1*ddz(1)*h1/k2(1) + &
      & 0.5*yy*dt*a_veg*c_veg*sqrt(um(1)**2 + vm(1)**2)*veg_sw(1))
      bm_(1,i,i) = xx
    enddo
    ! Adding friction at lateral boundaries
    cm_(1,1,1) = cm_(1,1,1) - yy*dt*lsu%water(1) !0.5*
    cm_(1,2,2) = cm_(1,2,2) - yy*dt*lsv%water(1) !0.5*
  
    bm_(1,1,2) = 0.
    bm_(1,2,1) = 0.
    cm_(1,1,2) = -kor*(h1*ddz(1))**2/(2*k2(1))
    cm_(1,2,1) = kor*(h1*ddz(1))**2/(2*k2(1))
  
    dm_(1,1) = yy*um(1) !-tau_ix*ddz*h1/k2(1)
    dm_(1,1) = dm_(1,1) + (um(2) - um(1))*xx + 2*taux*ddz(1)*h1/k2(1)*(1-b0) + &
    & kor*(h1*ddz(1))**2*vm(1)/(2*k2(1)) + yy*dt*g*tan(pi*alphax/180.) + &
    & yy*dt*(Force_x(1)) ! + 0.5*lsu%water(1)*um(1)
    dm_(1,2) = yy*vm(1) !-tau_iy*ddz*h1/k2(1)
    dm_(1,2) = dm_(1,2) + (vm(2) - vm(1))*xx + 2*tauy*ddz(1)*h1/k2(1)*(1-b0) - &
    & kor*(h1*ddz(1))**2*um(1)/(2*k2(1)) + yy*dt*g*tan(pi*alphay/180.) + &
    & yy*dt*(Force_y(1)) ! + 0.5*lsv%water(1)*vm(1)
  endif
  
  ! 3-rd case: thick ice:
  ! momentum flux = water-ice friction 
  
  if (l1 > L0) then
  
    !if (dyn_pgrad%par > 0) then
    !  Force_x(1:M+1) = lam_force * &
    !  & 1./h1*(tau_grx - tau_ix - dhw/dt*um(M+1) + dhw0/dt*(um(M+1)-um(1)))
    !  Force_y(1:M+1) = lam_force * &
    !  & 1./h1*(tau_gry - tau_iy - dhw/dt*vm(M+1) + dhw0/dt*(vm(M+1)-vm(1)))
    !else
      Force_x(:) = 0.d0
      Force_y(:) = 0.d0
    !endif
  
    do i = 1, 2
      cm_(1,i,i) = (xx + yy - coef1*ddz(1)*h1/k2(1) + &
      & 0.5*yy*dt*a_veg*c_veg*sqrt(um(1)**2 + vm(1)**2)*veg_sw(1))
      bm_(1,i,i) = xx
    enddo
    ! Adding friction at lateral boundaries
    cm_(1,1,1) = cm_(1,1,1) - yy*dt*lsu%water(1) !0.5*
    cm_(1,2,2) = cm_(1,2,2) - yy*dt*lsv%water(1) !0.5*
  
    bm_(1,1,2) = 0.
    bm_(1,2,1) = 0.
    cm_(1,1,2) = -kor*(h1*ddz(1))**2/(2*k2(1))
    cm_(1,2,1) = kor*(h1*ddz(1))**2/(2*k2(1))
  
    dm_(1,1) = yy*um(1) !-tau_ix*ddz*h1/k2(1)
    dm_(1,1) = dm_(1,1)+(um(2)-um(1))*xx + &
    & kor*(h1*ddz(1))**2*vm(1)/(2*k2(1)) + &
    & yy*dt*g*tan(pi*alphax/180.) + &
    & yy*dt*(Force_x(1)) ! + 0.5*lsu%water(1)*um(1)
    dm_(1,2) = yy*vm(1) !-tau_iy*ddz*h1/k2(1)
    dm_(1,2) = dm_(1,2)+(vm(2)-vm(1))*xx - &
    & kor*(h1*ddz(1))**2*um(1)/(2*k2(1)) + &
    & yy*dt*g*tan(pi*alphay/180.) + &
    & yy*dt*(Force_y(1)) ! + 0.5*lsv%water(1)*vm(1)
  endif
  
  k2(1) = area_int(1)/area_half(1)*k2(1)

elseif (cuette%par == 1) then

  ! B.c.s for Cuette flow
  ! Top boundary
  cm_(1,:,:) = 0.
  bm_(1,:,:) = 0.

  cm_(1,1,1) = 1.
  cm_(1,2,2) = 1.
  dm_(1,1) = u1(1)
  dm_(1,2) = v1(1)

endif bctop

bcbot : if (cuette%par == 0) then
  ! Boundary conditions for u and v at the lake bottom
  
  k2(M) = area_half(M)/area_int(M+1)*k2(M)
  
  xx = 0.5*(1 - ddz(M)*h1*(dhw - dhw0)/(k2(M)*dt))
  yy = (ddz(M)*h1)**2/(k2(M)*dt)
  
  do i = 1, 2
    cm_(M+1,i,i) = (xx + yy + coef2*ddz(M)*h1/k2(M) + &
    & 0.5*yy*dt*a_veg*c_veg*sqrt(um(M+1)**2 + vm(M+1)**2)*veg_sw(M+1))
    am_(M+1,i,i) = xx
  enddo
  
  am_(M+1,1,2) = 0.
  am_(M+1,2,1) = 0.
  cm_(M+1,1,2) = -kor*(h1*ddz(M))**2/(2*k2(M))
  cm_(M+1,2,1) = kor*(h1*ddz(M))**2/(2*k2(M))
  
  dm_(M+1,1) = yy*um(M+1) !-tau_grx*ddz*h1/k2(M)
  dm_(M+1,1) = dm_(M+1,1)-(um(M+1)-um(M))*xx + &
  & kor*(h1*ddz(M))**2*vm(M+1)/(2*k2(M)) + &
  & yy*dt*g*tan(pi*alphax/180.) + &
  & yy*dt*(Force_x(M+1)) ! + lsu%water(M+1)
  dm_(M+1,2) = yy*vm(M+1) !-tau_gry*ddz*h1/k2(M)
  dm_(M+1,2) = dm_(M+1,2)-(vm(M+1)-vm(M))*xx - &
  & kor*(h1*ddz(M))**2*um(M+1)/(2*k2(M)) + &
  & yy*dt*g*tan(pi*alphay/180.) + &
  & yy*dt*(Force_y(M+1)) ! + lsu%water(M+1)
  
  k2(M) = area_int(M+1)/area_half(M)*k2(M)

elseif (cuette%par == 1 .or. cuette%par == 11) then

  ! Bottom boundary
  cm_(M+1,:,:) = 0.
  am_(M+1,:,:) = 0.

  cm_(M+1,1,1) = 1.
  cm_(M+1,2,2) = 1.
  dm_(M+1,1) = u1(M+1)
  dm_(M+1,2) = v1(M+1)

endif bcbot



! The coefficients of equation for u and v in the internal points
! of the mesh

do k = 2, M

  do i = 1, 2
    do j = 1, 2
      am_(k,i,j) = KRON(i,j)*(dzeta_int(k)*dhw/(2.*(ddz(k-1)+ddz(k))*h1) - &
      & dhw0/(2.*(ddz(k-1)+ddz(k))*h1) - &
      & area_half(k-1)/area_int(k)*k2(k-1)*dt/(h1**2*(ddz(k)+ddz(k-1))*ddz(k-1) ))
      bm_(k,i,j) = -KRON(i,j)*(dzeta_int(k)*dhw/(2*(ddz(k-1)+ddz(k))*h1) + &
      & dhw0/(2*(ddz(k-1)+ddz(k))*h1) + &
      & area_half(k)/area_int(k)*k2(k)*dt/(h1**2*(ddz(k)+ddz(k-1))*ddz(k) ))
    enddo
  enddo

  cm_(k,1,1) = - &
  & (area_half(k)/area_int(k)*k2(k)/ddz(k) + &
  &  area_half(k-1)/area_int(k)*k2(k-1)/ddz(k-1))*dt/ &
  & (h1**2*(ddz(k)+ddz(k-1)) ) - 1. - &
  & 0.5*dt*c_veg*a_veg*sqrt(um(k)**2+vm(k)**2)*veg_sw(k)
  cm_(k,2,2) = cm_(k,1,1)
  ! Adding friction at lateral boundaries
  cm_(k,1,1) = cm_(k,1,1) - dt*lsu%water(k) !0.5*
  cm_(k,2,2) = cm_(k,2,2) - dt*lsv%water(k) !0.5*

  cm_(k,1,2) = kor2 !Krank-Nickolson
  cm_(k,2,1) = -cm_(k,1,2)

  dm_(k,1) = -um(k)-kor2*vm(k)-dt*g*tan(pi*alphax/180.) &
  & - (Force_x(k))*dt & ! + 0.5*lsu%water(k)*um(k)
  & - dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k) )**(-1)* &
  & (area_half(k)/area_int(k)*k2(k)*(um(k+1)-um(k))) &
  & + dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k-1) )**(-1)* &
  & (area_half(k-1)/area_int(k)*k2(k-1)*(um(k)-um(k-1))) &
  & - dzeta_int(k)*dhw*(2*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(um(k+1)-um(k-1)) &
  & + dhw0*(2.*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(um(k+1)-um(k-1))

  dm_(k,2) = -vm(k)+kor2*um(k)-dt*g*tan(pi*alphay/180.) &
  & - (Force_y(k))*dt & ! + 0.5*lsv%water(k)*vm(k)
  & - dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k) )**(-1)* &
  & (area_half(k)/area_int(k)*k2(k)*(vm(k+1)-vm(k))) &
  & + dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k-1) )**(-1)* &
  & (area_half(k-1)/area_int(k)*k2(k-1)*(vm(k)-vm(k-1))) &
  & - dzeta_int(k)*dhw*(2*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(vm(k+1)-vm(k-1)) &
  & + dhw0*(2.*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(vm(k+1)-vm(k-1))
enddo

! ind_bound=.false.
! call ind_stab_fact_db (a,b,c,1,M+1,indstab,ind_bound)
! if (indstab==.false.) then
! do i=2,M
!   a(i)=-k2(i)*dt/(h1**2*ddz**2)
!   b(i)=-k2(i+1)*dt/(h1**2*ddz**2)
! enddo
! endif  


call MATRIXPROGONKA(am_,bm_,cm_,dm_,ym_,M+1)
do k = 1, M+1
  u2(k) = ym_(k,1)
  v2(k) = ym_(k,2)
enddo

!print*, '(us, vs)', u2(1), v2(1)
     
windwork = row0*(taux*u2(1) + tauy*v2(1))

deallocate(Force_x, Force_y, um, vm)

if (firstcall) firstcall=.false.

END SUBROUTINE MOMENTUM_SOLVER


SUBROUTINE DHDXDHDY(M,u,v,area_int,Lx,Ly,ddz05,h1,dt,hx1,hx2,hy1,hy2)

! Calculates average levels of "quarters" of lake surface,
! assuming its elliptic or rectangular shape. Needed for average pressure
! gradient in momentum equations. The scheme of time integration is simple
! explicit. 

implicit none

!Input variables

integer(kind=iintegers), intent(in ) :: M ! Number of gridcells

real(kind=ireals), intent(in) :: u(1:M+1), v(1:M+1)   ! Horizontal velocity components, m/s
real(kind=ireals), intent(in) :: area_int(1:M+1)      ! Cross-section area, depth dependent, m**2
real(kind=ireals), intent(in) :: Lx(1:M+1), Ly(1:M+1) ! Length of X and Y central vertical cross-section of the lake body,
                                            ! depth dependent, m
real(kind=ireals), intent(in) :: ddz05(0:M)           ! Grid spacing, n/d
real(kind=ireals), intent(in) :: h1                   ! Lake depth, m
real(kind=ireals), intent(in) :: dt                   ! Timestep, s

!Input/output variables
real(kind=ireals), intent(inout) :: hx1, hx2, hy1, hy2 ! Average levels (heights above the lake bottom), of NW, NE, SE, SW
                                             ! quarters of the lake surface
!Ouput variables

!Local variables
real(kind=ireals) :: rint
integer(kind=iintegers) :: i

! Linearization: the change of cross-section area around 
! the cross-section of mean level in each half domain 
! (S,N,E and W) is neglected

rint = 0.
do i = 1, M+1
  rint = rint + h1*ddz05(i-1)*Lx(i)*v(i)
enddo
rint = pi*dt*rint/area_int(1)
hy2 = hy2 + rint
hy1 = hy1 - rint

rint = 0.
do i = 1, M+1
  rint = rint + h1*ddz05(i-1)*Ly(i)*u(i)
enddo
rint = pi*dt*rint/area_int(1)
hx2 = hx2 + rint
hx1 = hx1 - rint

END SUBROUTINE DHDXDHDY


SUBROUTINE DHDXDHDY2(itherm,M,u,v,area_int,area_half,Lx,Ly,ddz05,h1,dt, &
& hx1,hx2,hy1,hy2,hx1t,hx2t,hy1t,hy2t)

! Calculates average levels of "quarters" of lake surface and thermocline,
! assuming its elliptic or rectangular shape. Needed for average pressure
! gradient in momentum equations. The scheme of time integration is simple
! explicit. 

implicit none

!Input variables

integer(kind=iintegers), intent(in ) :: itherm ! The numerical level between two layers
integer(kind=iintegers), intent(in ) :: M ! Number of gridcells

real(kind=ireals), intent(in) :: u(1:M+1), v(1:M+1)   ! Horizontal velocity components, m/s
real(kind=ireals), intent(in) :: area_int   (1:M+1)   ! Cross-section area at cell interfaces, depth dependent, m**2
real(kind=ireals), intent(in) :: area_half  (1:M)     ! Cross-section area at cell centers, depth dependent, m**2
real(kind=ireals), intent(in) :: Lx(1:M+1), Ly(1:M+1) ! Length of X and Y central vertical cross-section of the lake body,
                                            ! depth dependent, m
real(kind=ireals), intent(in) :: ddz05(0:M)           ! Grid spacing, n/d
real(kind=ireals), intent(in) :: h1                   ! Lake depth, m
real(kind=ireals), intent(in) :: dt                   ! Timestep, s

!Input/output variables
real(kind=ireals), intent(inout) :: hx1, hx2, hy1, hy2 ! Average deviations of mixed-layer thickness, of NW, NE, SE, SW
                                                       ! quarters of the lake surface
real(kind=ireals), intent(inout) :: hx1t, hx2t, hy1t, hy2t ! Average deviations of layer below the mixed layer, of NW, NE, SE, SW
                                                           ! quarters of the lake
!Ouput variables

!Local variables
real(kind=ireals) :: rint
integer(kind=iintegers) :: i

! Linearization: the change of cross-section area around 
! the cross-section of mean level in each half domain 
! (S,N,E and W) is neglected


! Surface layer deviations
rint = 0.
do i = 1, itherm-1
  rint = rint + h1*ddz05(i-1)*Lx(i)*v(i)
enddo
rint = pi*dt*rint/area_int(1)
hy2 = hy2 + rint
hy1 = hy1 - rint

rint = 0.
do i = 1, itherm-1
  rint = rint + h1*ddz05(i-1)*Ly(i)*u(i)
enddo
rint = pi*dt*rint/area_int(1)
hx2 = hx2 + rint
hx1 = hx1 - rint

! Thermocline surface deviations
rint = 0.
do i = itherm, M+1
  rint = rint + h1*ddz05(i-1)*Lx(i)*v(i)
enddo
rint = pi*dt*rint/area_half(itherm-1)
hy2t = hy2t + rint
hy1t = hy1t - rint

rint = 0.
do i = itherm, M+1
  rint = rint + h1*ddz05(i-1)*Ly(i)*u(i)
enddo
rint = pi*dt*rint/area_half(itherm-1)
hx2t = hx2t + rint
hx1t = hx1t - rint

END SUBROUTINE DHDXDHDY2


SUBROUTINE DHDXDHDY3(itherm,M,u,v,area_int,area_half,Lx,Ly,ddz05,h1,dt, &
& hx1ml,hx2ml,hy1ml,hy2ml)

! Calculates average levels of "quarters" of lake surface and thermocline,
! assuming its elliptic or rectangular shape. Needed for average pressure
! gradient in momentum equations. The scheme of time integration is simple
! explicit. 

implicit none

!Input variables

integer(kind=iintegers), intent(in ) :: itherm(1:M+2) ! The numerical level between two layers
integer(kind=iintegers), intent(in ) :: M ! Number of gridcells

real(kind=ireals), intent(in) :: u(1:M+1), v(1:M+1)   ! Horizontal velocity components, m/s
real(kind=ireals), intent(in) :: area_int   (1:M+1)   ! Cross-section area at cell interfaces, depth dependent, m**2
real(kind=ireals), intent(in) :: area_half  (1:M)     ! Cross-section area at cell centers, depth dependent, m**2
real(kind=ireals), intent(in) :: Lx(1:M+1), Ly(1:M+1) ! Length of X and Y central vertical cross-section of the lake body,
                                            ! depth dependent, m
real(kind=ireals), intent(in) :: ddz05(0:M)           ! Grid spacing, n/d
real(kind=ireals), intent(in) :: h1                   ! Lake depth, m
real(kind=ireals), intent(in) :: dt                   ! Timestep, s

!Input/output variables
real(kind=ireals), intent(inout) :: hx1ml(1:M+1), hx2ml(1:M+1), hy1ml(1:M+1), hy2ml(1:M+1) 

!Local variables
real(kind=ireals) :: rint, areatop
integer(kind=iintegers) :: i, j

! Linearization: the change of cross-section area around 
! the cross-section of mean level in each half domain 
! (S,N,E and W) is neglected

main_do : do i = 1, M+1

  if (itherm(i+1) > 0) then

    if (i == 1) then
      areatop = area_int(1)
    else
      areatop = area_half(itherm(i))
    endif

    rint = 0.
    do j = itherm(i)+1, itherm(i+1)
      rint = rint + h1*ddz05(j-1)*Lx(j)*v(j)
    enddo
    rint = pi*dt*rint/areatop
    hy2ml(i) = hy2ml(i) + rint
    hy1ml(i) = hy1ml(i) - rint
    
    rint = 0.
    do j = itherm(i)+1, itherm(i+1)
      rint = rint + h1*ddz05(j-1)*Ly(j)*u(j)
    enddo
    rint = pi*dt*rint/areatop
    hx2ml(i) = hx2ml(i) + rint
    hx1ml(i) = hx1ml(i) - rint

  else
    exit main_do
  endif

enddo main_do

END SUBROUTINE DHDXDHDY3


SUBROUTINE DHDXDHDYML(M,u,v,area_int,area_half,Lx,Ly,ddz05,h1,dt,hx1ml,hx2ml,hy1ml,hy2ml)

! Calculates average displacements of "quarters" of each computational layer, 
! assuming their elliptic or rectangular shape. Needed for average pressure
! gradient in momentum equations. The scheme of time integration is simple
! explicit. 

implicit none

!Input variables

integer(kind=iintegers), intent(in ) :: M ! Number of gridcells

real(kind=ireals), intent(in) :: u(1:M+1), v(1:M+1)   ! Horizontal velocity components, m/s
real(kind=ireals), intent(in) :: area_int   (1:M+1)   ! Cross-section area at cell interfaces, depth dependent, m**2
real(kind=ireals), intent(in) :: area_half  (1:M)     ! Cross-section area at cell centers, depth dependent, m**2
real(kind=ireals), intent(in) :: Lx(1:M+1), Ly(1:M+1) ! Length of X and Y central vertical cross-section of the lake body,
                                                      ! depth dependent, m
real(kind=ireals), intent(in) :: ddz05(0:M)           ! Grid spacing, n/d
real(kind=ireals), intent(in) :: h1                   ! Lake depth, m
real(kind=ireals), intent(in) :: dt                   ! Timestep, s

!Input/output variables
real(kind=ireals), intent(inout) :: hx1ml(1:M+1), hx2ml(1:M+1), hy1ml(1:M+1), hy2ml(1:M+1) ! Average deviations of layer's thickness, of NW, NE, SE, SW
                                                                                           ! quarters their surface
!Ouput variables

!Local variables
real(kind=ireals) :: rint
integer(kind=iintegers) :: i

! Linearization: the change of cross-section area around 
! the cross-section of mean level in each half domain 
! (S,N,E and W) is neglected


! Layers' thickness deviations

rint = h1*ddz05(0)*Ly(1)*u(1)
rint = pi*dt*rint/area_int(1)
hx2ml(1) = hx2ml(1) + rint
hx1ml(1) = hx1ml(1) - rint
do i = 2, M+1
  rint = h1*ddz05(i-1)*Ly(i)*u(i)
  rint = pi*dt*rint/area_half(i-1)
  hx2ml(i) = hx2ml(i) + rint
  hx1ml(i) = hx1ml(i) - rint
enddo

rint = h1*ddz05(0)*Lx(1)*v(1)
rint = pi*dt*rint/area_int(1)
hy2ml(1) = hy2ml(1) + rint
hy1ml(1) = hy1ml(1) - rint
do i = 2, M+1
  rint = h1*ddz05(i-1)*Lx(i)*v(i)
  rint = pi*dt*rint/area_half(i-1)
  hy2ml(i) = hy2ml(i) + rint
  hy1ml(i) = hy1ml(i) - rint
enddo

END SUBROUTINE DHDXDHDYML


!> Subroutine seeks layers of quasiconstant density in continuous density profile
SUBROUTINE DENSLAYERS(M,i_maxN,nstep,bathymwater,ddz05,ddz,row,h1,nlevs,itherm,rowlars)

use ARRAYS_BATHYM, only : &
& bathym, aki, bki

use ARRAYS_TURB, only : rowc

use PHYS_CONSTANTS, only : g, row0

implicit none

!Input variables
integer(kind=iintegers), intent(in) :: M
integer(kind=iintegers), intent(in) :: i_maxN
integer(kind=iintegers), intent(in) :: nstep !> The number of timestep

type(bathym), intent(in) :: bathymwater(1:M+1)

real(kind=ireals), intent(in) :: ddz05(0:M)
real(kind=ireals), intent(in) :: ddz(1:M)
real(kind=ireals), intent(in) :: row(1:M+1)
real(kind=ireals), intent(in) :: h1

!Output variables
integer(kind=iintegers), intent(out) :: nlevs
integer(kind=iintegers), intent(inout) :: itherm(1:M+2)
real(kind=ireals), intent(out) :: rowlars(1:M+1)

!Local variables
real(kind=ireals), parameter :: a = 1.e-1, u = 0.001
real(kind=ireals) :: omega, dro, dro_, lnumx, lnumy,ldenx,ldeny !, maxN
real(kind=ireals), allocatable :: Hlars(:), work(:)
!integer(kind=iintegers) :: iwork(1:M+2)
real(kind=ireals), parameter :: Nmin = 1.E-1
integer(kind=iintegers), allocatable :: itherm_test(:)
integer(kind=iintegers), parameter :: layerapp = 3
integer(kind=iintegers) :: i, j, k

logical, save :: firstcall = .true.

!iwork(:) = - 1
if (layerapp == 1) then
  itherm(:) = - 1
  ! General algorithm
  j = 1
  itherm(j) = 0
  do i = 2, M+1
    omega = sqrt(g*ddz05(i-2)*ddz05(i-1)*h1*h1*max(row(i) - row(i-1),0._ireals)/ &
    & (row(i)*h1*(ddz05(i-2) + ddz05(i-1))))
    !print*, 0.5*omega*a/(u*ddz05(i-1)*h1)
    !if (0.5*omega*a/(u*ddz05(i-1)*h1) > 1._ireals) then
    !print*, omega*a/u
    if (TESTAMPL(ddz05(i-2)*h1,ddz05(i-1)*h1,row(i-1),row(i))) then
      itherm(j+1) = i-1
      j = j + 1
    endif
  enddo
  nlevs = j
elseif (layerapp == 2) then
  ! Three-layer approximation
  !maxN = sqrt( max( g/row0*(row(i_maxN+1) - row(i_maxN-1))/ &
  !&               ( h1*(ddz(i_maxN-1) + ddz(i_maxN)) ), 0._ireals) )
  !i = max(i_maxN,nint(0.1*M)); j = max(i_maxN,nint(0.1*M))
  !do
  !  if (sqrt( max(g/row0*(row(i+1) - row(i))/(h1*ddz(i)), 0._ireals))/maxN < 9.E-1) exit
  !  i = i - 1
  !enddo
  !do
  !  if (sqrt( max(g/row0*(row(j) - row(j-1))/(h1*ddz(j-1)), 0._ireals))/maxN < 1.E-1) exit
  !  j = j + 1
  !enddo

  !j = 1
  !iwork(j) = 0
  !do i = 2, M+1
  !  omega = sqrt(g*ddz05(i-2)*ddz05(i-1)*h1*h1*max(row(i) - row(i-1),0._ireals)/ &
  !  & (row(i)*h1*(ddz05(i-2) + ddz05(i-1))))
  !  !print*, 0.5*omega*a/(u*ddz05(i-1)*h1)
  !  !if (0.5*omega*a/(u*ddz05(i-1)*h1) > 1._ireals) then
  !  !print*, omega*a/u
  !  if (omega*a/u > 1._ireals) then
  !    iwork(j+1) = i-1
  !    j = j + 1
  !  endif
  !enddo

  !itherm(1) = 0
  !itherm(2) = iwork(2) !i+1
  !itherm(3) = iwork(j)

  !itherm(1) = 0
  !itherm(2) = 11
  !!itherm(3) = 21
  !itherm(3) = 31

  !if (mod(nstep-1,10000_iintegers) == 0) then
  if (nstep == 1 .or. nstep == 50000_iintegers) then  
    itherm(:) = - 1
    itherm(1) = 0
    dro = row(M+1) - row(1)
    do i = 2, M+1
      if (row(i) - sum(row(1:i-1))/real(i-1) > 0.2*dro) then
        itherm(2) = i-1
        exit
      endif
    enddo
    do i = M, 1, -1
      if (sum(row(i+1:M+1))/real(M-i+1) - row(i) > 0.1*dro) then
        itherm(3) = i
        exit
      endif
    enddo
  endif
  !print*, itherm(1:3)
  nlevs = 3
elseif (layerapp == 3) then
  nlevs = 1
  itherm(:) = - 1
  itherm(1) = 0
  itherm(2) = M+1
  dro = row(M+1) - row(1)
  if (dro > 0.) then
    allocate(itherm_test(1:M+2), Hlars(1:M+1), work(1:M+1))
    seekl : do
      itherm_test(:) = - 1
      itherm_test(1) = 0
      itherm_test(nlevs+2) = M+1
      dro_ = dro/real(nlevs+1)
      do j = 1, nlevs
        work(:) = (row(1) + j*dro_) - row(:)
        itherm_test(j+1) = minloc(work, dim = 1, mask = work > 0. )
        if (itherm_test(j+1) <= itherm_test(j)) exit seekl
        Hlars(j) = sum(ddz05(itherm_test(j):itherm_test(j+1)-1))*h1
      enddo
      Hlars(nlevs+1) = sum(ddz05(itherm_test(nlevs+1):itherm_test(nlevs+2)-1))*h1
      call AVERROW(nlevs+1,itherm_test)
      do j = 1, nlevs
        if (.not. TESTAMPL(Hlars(j),Hlars(j+1),rowlars(j),rowlars(j+1))) exit seekl
      enddo
      nlevs = nlevs + 1
      itherm = itherm_test
    enddo seekl 
    deallocate(itherm_test,Hlars,work)
  endif
endif
!itherm(2) = i_maxN-1
!nlevs=2

itherm(nlevs+1) = M+1

call AVERROW(nlevs,itherm)

!rowlars(1) = 997.90623796270870        
!rowlars(2) = 998.27504425868335        
!rowlars(3) = 998.63226709221465        
!rowlars(4) = 999.05361626201648        
!rowlars(5) = 999.43047189124843        
!rowlars(6) = 999.94501214548336
!print*, nlevs, rowlars
!stop
!itherm(1) = 0; itherm(2)=10; itherm(3)=15; itherm(4)=19; itherm(5)=121
!rowlars(1) = 998.04; rowlars(2)=998.54; rowlars(3)=999.13; rowlars(4)=999.94;

do i = 1, nlevs
  forall (j = itherm(i)+1:itherm(i+1)) rowc(j) = rowlars(i)
enddo

if (allocated(aki)) deallocate(aki)
if (allocated(bki)) deallocate(bki)

allocate(aki(1:nlevs,1:nlevs),bki(1:nlevs,1:nlevs))

!Pressure gradient transfer matrix
do j = 1, nlevs
  do i = 1, nlevs
    !if (i >= j) then
    !  aki(i,j) = 1.
    !  bki(i,j) = 1.
    !else
      if (i == 1) then
        ldenx = bathymwater(1_iintegers)%Lx
        lnumx = bathymwater(itherm(j))%Lx_half
        ldeny = bathymwater(1_iintegers)%Ly
        lnumy = bathymwater(itherm(j))%Ly_half
        aki(i,j) = sqrt(ldenx*ldeny/(lnumx*lnumy)) !sin(pi*0.5*lnum/lden)
        bki(i,j) = sqrt(ldenx*ldeny/(lnumx*lnumy)) !sin(pi*0.5*lnum/lden)
      else
        ldenx = bathymwater(itherm(i))%Lx_half
        ldeny = bathymwater(itherm(i))%Ly_half
        lnumx = bathymwater(itherm(j))%Ly_half
        lnumy = bathymwater(itherm(j))%Lx_half
        aki(i,j) = sqrt(ldenx*ldeny/(lnumx*lnumy)) !sin(pi*0.5*lnum/lden)
        bki(i,j) = sqrt(ldenx*ldeny/(lnumx*lnumy)) !sin(pi*0.5*lnum/lden)
      endif
    !endif
  enddo
enddo

!print*, 'aki', aki
!read*
!print*, 'bki', bki
!read*

!print*, nlevs, itherm
!print*, rowlars
!read*

if (firstcall) firstcall = .false. 

contains
!> Testing the smallness of layers thickness deviation in two-layer linear model
FUNCTION TESTAMPL(H1,H2,ro1,ro2)
implicit none
real(kind=ireals), intent(in) :: H1, H2, ro1 ,ro2
logical :: TESTAMPL
omega = sqrt(g * H1 * H2 * max(ro2 - ro1,0._ireals) / (row(i) * (H1 + H2)))
TESTAMPL = (omega*a/u > 1._ireals)
!print*, H1, H2, ro2-ro1, omega
END FUNCTION TESTAMPL

!> Averaging of density over layers of constant density
SUBROUTINE AVERROW(nlev,ithermm)
implicit none
integer(kind=iintegers), intent(in) :: nlev, ithermm(1:M+2)
! Have to be extended to include bathymetry
!print*, nlev, ithermm
rowlars(:) = -1.
do i = 1, nlev
  rowlars(i) = 0.
  do k = ithermm(i)+1,ithermm(i+1)
    rowlars(i) = rowlars(i) + ddz05(k-1)*row(k)
  enddo
  rowlars(i) = rowlars(i)/sum(ddz05(ithermm(i):ithermm(i+1)-1))
enddo
END SUBROUTINE AVERROW


END SUBROUTINE DENSLAYERS


!> Subroutine RELAYER distributes homogeneous-density-layers thickness deviations
!! induced by seiches from one set of layers to the other. NOTE: bathymetry effects not included
SUBROUTINE RELAYER(M,itherm,itherm_new,hx1ml,hx2ml,hy1ml,hy2ml)

implicit none

!Input variables
integer(kind=iintegers), intent(in ) :: M !> The number of computational layers
integer(kind=iintegers), intent(in ) :: itherm(1:M+2) !> The numerical level between two layers
integer(kind=iintegers), intent(in ) :: itherm_new(1:M+2) !> The numerical level between two layers, at the new timestep

!Input/output variables
real(kind=ireals), intent(inout) :: hx1ml(1:M+1), hx2ml(1:M+1), hy1ml(1:M+1), hy2ml(1:M+1) !> Average deviations of layer's thickness, of NW, NE, SE, SW
                                                                                           !> quarters their surface

!Local variables
integer(kind=iintegers) :: i
real(kind=ireals), allocatable :: hx1(:), hx2(:), hy1(:), hy2(:)


allocate (hx1(1:M+1), hx2(1:M+1), hy1(1:M+1), hy2(1:M+1))

! Distribution of const-dens layers thickness deviations between computational layers
do i = 1, M+1
  if (itherm(i)+1 /= 0) then ! itherm(i) /= -1
    hx1(itherm(i)+1:itherm(i+1)) = hx1ml(i)/real(itherm(i+1) - itherm(i))
    hx2(itherm(i)+1:itherm(i+1)) = hx2ml(i)/real(itherm(i+1) - itherm(i))
    hy1(itherm(i)+1:itherm(i+1)) = hy1ml(i)/real(itherm(i+1) - itherm(i))
    hy2(itherm(i)+1:itherm(i+1)) = hy2ml(i)/real(itherm(i+1) - itherm(i))
  endif
enddo

! Constructing the thickness deviations for new const-dens layers
do i = 1, M+1
  if (itherm_new(i)+1 /= 0) then ! itherm_new(i) /= -1
    hx1ml(i) = sum(hx1(itherm_new(i)+1:itherm_new(i+1)))
    hx2ml(i) = sum(hx2(itherm_new(i)+1:itherm_new(i+1)))
    hy1ml(i) = sum(hy1(itherm_new(i)+1:itherm_new(i+1)))
    hy2ml(i) = sum(hy2(itherm_new(i)+1:itherm_new(i+1)))
  endif
enddo

deallocate(hx1, hx2, hy1, hy2)

END SUBROUTINE RELAYER


END MODULE MOMENTUM





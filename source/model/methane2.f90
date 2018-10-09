    SUBROUTINE METHANE2 &
    & (pressure,wind10,zsoil,Tsoil,wl,wi,wa,rootss,rhogr,por,veg,qsoil,TgrAnn, &
    & ddz, Twater, lammeth, qwater, h1, l1, ls1, dhw, dhw0, &
    & fplant, febul, fdiff, ftot, fdiff_lake_surf, &
    & plant_sum,bull_sum,oxid_sum,rprod_sum, &
    & anox,M,ns,dt,eps_surf,rprod_total_oldC,rprod_total_newC, &
    & h_talik,tot_ice_meth_bubbles)
 
    use LAKE_DATATYPES, only : ireals, iintegers
    use DRIVING_PARAMS !, only : &
!    & tricemethhydr 
    use PHYS_CONSTANTS, only : &
    & row0, g, row0, roa0, &
    & row0_d_roi, row0_d_roa0, &
    & roa0_d_roi, roa0_d_row0
    use METH_OXYG_CONSTANTS
    use PHYS_PARAMETERS
    use NUMERIC_PARAMS
    use PHYS_FUNC, only : &
    & WL_MAX, &
    & HENRY_CONST, &
    & DIFF_WATER_METHANE, &
    & DIFF_AIR_METHANE, &
    & GAS_WATATM_FLUX
    use NUMERICS, only : &
    PROGONKA
    
    implicit none
    
!*  ===============================================================
!*  prod + oxid
!*	solution 1-D nonlineal equation
!*	describing production, oxidation
!*	and transport of methane from wetlands
!*  ===============================================================
!*	input parameters:
!*	=====================
!*    z(ml)    :	height of level in soil from soil surface (cm)
!*    Tsoil(ml)    :	the soil temperature (C degrees)
!*
!*	output parameters:
!*	=====================
!*    fplant   : the ch4 flux due to plant-mediated transport
!*    febul    : the ch4 flux due to ebullition
!*    fdiff    : the diffusive ch4 flux at soil/water=atm boundary
!*    ftot	   : the total ch4 emission
!*  =============================================================
!*	parameters for production
!*	=============================
!*    rnpp     ! NPP (gC/(m**2*mo))
!*	  rnppmax  ! the maximum value of the NPP
!*    rnroot   ! the rooting depth (cm)
!*    q100     ! value lying within the range of observed
!*             ! q100 values ranging from 1.7 to 16
!*    r0       ! the constant rate factor (muM/s)
!*	=============================================================
!*	parameters for oxidation
!*	=============================
!*	  vmax	        ! Michaelis-Menten (muM/s)
!*	  rkm	        ! coefficients (muM)
!*	  q10	        ! the observed value for oxidation
!*	  TgrAnn(ml)	! the annual mean soil temperature (deg C)
!*	=============================================================
!*	parameters for ebullition
!*	=============================
!*    rke	!	a rate constant of the unit 1/s
!*    cmin	!	the concentration at which bubble
!*          !   formation occurs (muM)
!*    punveg	!	the percentage of unvegetated, bare soil
!*    cthresh	!	the threshold concentration
!*              !   for bubble formation
!*	=============================================================
!*	parameters for plant-mediated transport
!*	===========================================
!*     rkp	!	a rate constant of the unit 0.01/s
!*     tveg	!	a factor describing the quality of
!*          !   plant-mediated transport at a site
!*     pox	!	a certain fraction of methane oxidized in
!*          !   the oxic zone around the roots of plants
!*     rlmin   ! the parameter for the function of fgrow(t)
!*     rl      ! the parameter for the function of fgrow(t)
!*     rlmax   ! rlmax=rlmin+rl
!*	===============================================================

    real(kind=ireals), intent(in) :: pressure ! Atmospheric pressure, Pa
    real(kind=ireals), intent(in) :: wind10 ! Wind speed at 10 m above the surface    
    real(kind=ireals), intent(in) :: zsoil(1:ns) ! numerical levels in soil, meters
    real(kind=ireals), intent(in) :: Tsoil(1:ns) ! Celsius
    real(kind=ireals), intent(in) :: wl(1:ns) ! liquid water content in a soil, kg/kg
    real(kind=ireals), intent(in) :: wi(1:ns) ! ice content in a soil, kg/kg
    real(kind=ireals), intent(in) :: wa(1:ns) ! air content in a soil, kg/kg
      
!      COMMON /VEINIT/ 
    real(kind=ireals), intent(in) :: rootss(1:ns) ! meters
      
!      COMMON /SOILDAT/ 
    real(kind=ireals), intent(in) :: rhogr(1:ns) ! density of dry soil, kg/m**3
    real(kind=ireals), intent(in) :: por(1:ns) ! soil porosity, n/d
    
    real(kind=ireals), intent(in) :: ddz(1:M) ! the thickness of dzeta-layers in water, n/d
    real(kind=ireals), intent(in) :: Twater(1:M+1) ! the water temperature, degrees Celsius
    real(kind=ireals), intent(in) :: lammeth(1:M) ! the heat turbulent diffusivity in water, m**2/s
    real(kind=ireals), intent(in) :: h1, l1, ls1 ! the thickness (depth) of the water column, m
                                       ! ice and bottom ice
    real(kind=ireals), intent(in) :: dhw, dhw0   ! the time increment of water column thickness
                                       ! and those at its top, m
    
    real(kind=ireals), intent(in) :: dt ! Timestep, sec
    real(kind=ireals), intent(in) :: eps_surf ! TKE dissipation rate, m**2/s**3

    integer(kind=iintegers), intent(in) :: M, ns ! Number of levels in water and soil
      
    real(kind=ireals), intent(inout) :: veg
            
!      common/ch4/ 
    real(kind=ireals), intent(inout) :: qsoil(1:ns,1:2) ! Methane concentration in soil
    real(kind=ireals), intent(inout) :: qwater(1:M+1,1:2,1:2) ! Methane concentration in water
    real(kind=ireals), intent(inout) :: TgrAnn(1:ns)
    real(kind=ireals), intent(inout) :: fplant, ftot, fdiff_lake_surf(1:2)
    real(kind=ireals), intent(inout) :: febul(1:2), fdiff(1:2)
    real(kind=ireals), intent(inout) :: plant_sum, bull_sum,oxid_sum,rprod_sum
    real(kind=ireals), intent(inout) :: anox
    real(kind=ireals), intent(inout) :: tot_ice_meth_bubbles(1:2)
    
    real(kind=ireals), intent(out) :: rprod_total_oldC
    real(kind=ireals), intent(out) :: rprod_total_newC
    real(kind=ireals), intent(out) :: h_talik ! talik depth, m
      
!     Local variables      

    integer(kind=iintegers), parameter :: water_methane_indic = 8
    
    real(kind=ireals), parameter :: C_depth_age = 0.25 ! from 0.25 to 0.30, according to West and Plug, 2007
    real(kind=ireals), parameter :: C_talik_age = 0.5  ! from 0.5 to 0.7, according to West and Plug, 2007
    real(kind=ireals), parameter :: Tmelt_pnt = 0. ! melting point temperature, degrees Celsius
    
    real(kind=ireals), save :: fcoarse = 0.5 !the relative volume of the coarse pores
    real(kind=ireals), save :: scf = 16.04/(1./864.)! the scale factor  mgCH4/m**2/day
!	 scf = 12./(24.*1./864)! the scale factor  mgC/m**2/h

    real(kind=ireals) :: diff(ns), forg(ns)
    real(kind=ireals) :: qplant(ns),qebul(1:ns,1:2),q_old(ns),bull(1:ns,1:2)
    real(kind=ireals) :: ch4_crit(1:ns,1:2)
    real(kind=ireals) :: oxid(ns),dz_soil(ns),z2(ns)
    real(kind=ireals) :: fcch4(ns,1:2)
    real(kind=ireals) :: bp(ns), ap(ns)
    real(kind=ireals) :: rprod(ns), rprod_new(ns), rprod_old(ns), plant(ns)
    real(kind=ireals) :: n2(ns)
    real(kind=ireals) :: water_porvolrat(1:ns), air_porvolrat(1:ns)
    
    real(kind=ireals), allocatable :: a(:), b(:), c(:), f(:), y(:)
!    real(kind=ireals), allocatable, save :: mean_rprod_oldc0(:)

    real(kind=ireals) :: rnroot
    real(kind=ireals) :: froot(ns), anox_crit(ns)
    
    real(kind=ireals) :: pow_oldorg
    real(kind=ireals) :: lake_age
      
    real(kind=ireals) :: tmp_veg, t_grow, t_mature, fnpp
    real(kind=ireals) :: fin, ft, ft1, t50, fgrow, q_sum_old, fdiff1, q_sum_new
    real(kind=ireals) :: balance
    real(kind=ireals) :: Flux_atm(1:2)
    real(kind=ireals) :: ratio_c13_c14
    real(kind=ireals) :: xx, yy ! working variables
    
    real(kind=ireals), allocatable :: z(:)
    real(kind=ireals), allocatable :: roots(:)
    
    real(kind=ireals), external :: DZETA
      
    integer(kind=iintegers) :: i, j ! Loop indices
    integer(kind=iintegers) :: i_talik

    real(kind=ireals), external :: MEAN_RPROD_OLDC, MEAN_RPROD_OLDC2
    
    allocate (z(1:ns))
    allocate (roots(1:ns))
    allocate (a(1:vector_length),b(1:vector_length),c(1:vector_length), &
    &         f(1:vector_length),y(1:vector_length))

   
    z(:) = zsoil(:)
    roots(:) = rootss(:)

!    TgrAnn(ns) = Tsoil(ns) !ms
    TgrAnn(ns) = 0. ! exludes the factor of annual mean temperature
    
!   Calculation of pore volume ratios for liquid water and air
!   They are needed for diffusivity calculations 
    do i = 1, ns
      water_porvolrat(i) = wl(i) / &
      & (por(i)*(row0/rhogr(i) + wl(i) + row0_d_roi*wi(i) + row0_d_roa0*wa(i)) )
      air_porvolrat(i) = wa(i) / &
      & (por(i)*(roa0/rhogr(i) + wa(i) + roa0_d_roi*wi(i) + roa0_d_row0*wl(i)) )
    enddo

!   Assuming exponential dependence of nitrogen concentration on depth
    n2(1) = 2.*n2_atm0*(1.-exp(-n2_exp_decay*0.5*z(2)) )/(z(2)*n2_exp_decay)
    n2(2:ns) = 0. ! nitrogen concentration in the layers below the surface is
                  ! set to zero, since it decays rapidly with depth (Bazhin, 2001)
                  
!   Calculation of the talik depth and of the power in production term due to 
!   old talik organics decomposition 
    i_talik = 0    
    do i = 1, ns-1
      if (Tsoil(i) > Tmelt_pnt .and. Tsoil(i+1) <= Tmelt_pnt) then
        h_talik = z(i) + (z(i+1) - z(i)) * &
        & (Tmelt_pnt - Tsoil(i))/(Tsoil(i+1) - Tsoil(i))
        if (h_talik < 0.5*(z(i) + z(i+1))) then
          i_talik = i
        else
          i_talik = i + 1
        endif
        exit
      endif
    enddo

!    if (.not.allocated(mean_rprod_oldc0)) then
!      allocate (mean_rprod_oldc0(1:ns))
!      do i = 2, ns-1
!        mean_rprod_oldc0(i) = &
!        & MEAN_RPROD_OLDC(0.5*(z(i-1) + z(i)), 0.5*(z(i+1) + z(i)), &
!        & h_talik, C_talik_age)
!      enddo
!      mean_rprod_oldc0(1) = &
!      & MEAN_RPROD_OLDC(z(1), 0.5*(z(2) + z(1)), &
!      & h_talik, C_talik_age)
!      mean_rprod_oldc0(ns) = &
!      & MEAN_RPROD_OLDC(0.5*(z(ns-1) + z(ns)), z(ns), &
!      & h_talik, C_talik_age)
!    endif
    
!   Calculation of approximate THERMOKARST lake age, according to results by West and Plug, 2007
!    lake_age = (h1/C_depth_age)**2 ! years
!    lake_age = (h_talik/C_talik_age)**2 ! years
!    pow_oldorg = 2.*alpha_old_org*lake_age
                  
    do i = 1, ns
      xx = HENRY_CONST(Henry_const0_ch4, Henry_temp_dep_ch4, Henry_temp_ref, Tsoil(i)+273.15)
      yy = HENRY_CONST(Henry_const0_n2, Henry_temp_dep_n2, Henry_temp_ref, Tsoil(i)+273.15)
      ch4_crit(i,1) = qsoil(i,1)*rel_conc_ebul_crit*por(i)*xx*(pressure + row0*g*h1 - n2(i)/yy) / &
      & (qsoil(i,1) + qsoil(i,2))
      ch4_crit(i,2) = ch4_crit(i,1)*qsoil(i,2)/qsoil(i,1)
      ! pressure -> 1.d+5
    enddo
    
    do i = 1, ns !ms+1, ml
!      anox_crit(i) = (por(i)*rhow/rhogr(i)-wi(i)*rhow/rhoi)*0.9
      anox_crit(i) = WL_MAX(por(i),rhogr(i),wi(i),tricemethhydr%par)*0.9
    end do

!   Calculation of rooting depth
    rnroot = 0.
    do i = 1, ns !ms, ml
      if (roots(i).gt.0.) rnroot = z(i)
    end do
    
!*	============================================================
!*	tgrow is the temperature at which plants start to grow
!*	============================================================

    if (Tsoil(2) .lt. 5.) then ! 5 deg. Celsius
      t_grow = 2. ! deg. Celsius
    else
      t_grow = 7. ! deg. Celsius
    end if
    t_mature = t_grow + 10. ! the temperature at which plants reach maturity, deg. Celsius
    
!*  =============================================================
!*	parameters for diffusion
!*	=============================

    anox = 0.
    do i = ns-1, 2, -1 !ml-1, ms+1, -1
      if (z(i) .le. 1.) then  ! 100 cm
        if (wl(i) .ge. anox_crit(i)) then
          anox = anox + (z(i+1)-z(i-1))*0.5
        end if
      end if
    end do

    if (rnroot > 0.) then
      do i = 1, ns !ms, ml
        if(z(i) .le. rnroot) then
          froot(i) = 2.*(rnroot-z(i))/rnroot
        else
          froot(i) = 0.
        end if
      end do
    else
      froot(:) = 0.
    endif

    do i = 1, ns
      if (z(i) .gt. rnroot) then                  ! for vegetated
        forg(i) = veg * exp((rnroot - z(i))*10.) ! soil
      else
!      if (z(i).le.rnroot) then
        forg(i) = veg * 1.
      endif
    end do
    
    do i = 1, ns                                                      ! for unvegetated
      forg(i) = forg(i) + (1.-veg) * forg0*exp(-z(i)*lambda_new_org) ! soils
!      forg(i) = forg(i) + (1.-veg) * forg0 ! the option for thermokarst lakes
                                            ! for which it is assumed the organic
                                            ! matter is homogeneously distributed 
                                            ! in soil
    end do

    do i = 2, ns 
!      if (wl(i) .lt. anox_crit(i)) then
!        diff(i) = diff_unsat*tortuosity_coef*fcoarse  ! in the unsaturated
!        
!      else          ! soil layers
!        diff(i) = diff_water*tortuosity_coef*fcoarse     ! in the water saturated
!        diff(i) = DIFF_WATER_METHANE(Tsoil(i)) * &
!        & tortuosity_coef*fcoarse
!      end if                                             ! soil layers
      diff(i) = 0.25*( por(i) + por(i-1) )*tortuosity_coef * &
      & ( (water_porvolrat(i) + water_porvolrat(i-1)) * &
      & DIFF_WATER_METHANE( 0.5*(Tsoil(i) + Tsoil(i-1)) ) + &
      & (air_porvolrat(i) + air_porvolrat(i-1)) * &
      & DIFF_AIR_METHANE( 0.5*(Tsoil(i) + Tsoil(i-1)) ) )
    end do

    fnpp = rnpp             !function of the variation of NPP with t
    fin = 1. + fnpp/rnppmax !variation of substrate availability with t

    rprod_total_newC = 0.
    rprod_total_oldC = 0.
    
    do i = 1, ns 
    
      if (wl(i) .ge. anox_crit(i)) then
        ft1 = 1. ! Note, methane generation through old organics decomposition
                 ! is zero below talik automatically, see the code below
        if (Tsoil(i) .gt. 0) then
          if (i /= i_talik) then
            ft = 1.
          else
            ft = (h_talik - 0.5*(z(i-1) + z(i)) ) / &
            & (0.5*(z(i+1) - z(i-1)))
          endif
        else
          if (i /= i_talik) then
            ft = 0.
          else
            ft = (h_talik - 0.5*(z(i-1) + z(i)) ) / &
            & (0.5*(z(i+1) - z(i-1)))
          endif
        endif
      else
        ft = 0.
        if (i /= ns .and. i /= 1) then
          if ((Tsoil(i) > 0. .and. Tsoil(i+1) < 0.) .or. &
          &   (Tsoil(i) < 0. .and. Tsoil(i-1) > 0.) ) then
!           Regularization since there is always minimum of soil moisture content
!           near the talik depth (due to diffusion from water saturated layers above
!           to frozen layers below)
            ft1 = 1. - (anox_crit(i) - wl(i))/anox_crit(i)
!            ft1 = 1.
!            ft = ft1
!            if (i == i_talik) &
!            & ft = ft * (h_talik - 0.5*(z(i-1) + z(i)) ) / &
!            & (0.5*(z(i+1) - z(i-1)))
          else
            ft1 = 0.
          endif
        else
          ft1 = 0.
        endif
      endif
      
!      xx = r0_oldorg*(min(zsoil(i)/h_talik,1.d0))**pow_oldorg * &
!      & q100**((Tsoil(i)-TgrAnn(ns))*0.1)*ft
!      xx = r0_oldorg*exp(-alpha_old_org*C_talik_age**(-2) * &
!      & (h_talik**2-(min(z(i),h_talik))**2))*ft * &
!      & q100**((Tsoil(i)-TgrAnn(ns))*0.1)

      ratio_c13_c14 = min(z(i)/h_talik, 1.d0)
!     Note that for the performance optimization MEAN_RPROD_OLDCs may
!     be computed only once - at the first timestep

!     Note that calculating MEAN_RPROD_OLDC the temperature step
!     function weighting is already implemented, so ft multiplication
!     is not needed
      if (i == 1) then
!        xx = MEAN_RPROD_OLDC(0.d0, 0.5*(z(2) + z(1)), h_talik, C_talik_age)
        xx = MEAN_RPROD_OLDC2(0.d0, 0.5*(z(2) + z(1)), h_talik, C_talik_age)        
!        xx = exp(-alpha_old_org*C_talik_age**(-2)*h_talik**2) * &
!        & mean_rprod_oldc0(1)
        rprod_old(i) = xx*q100**((Tsoil(1)-TgrAnn(ns))*0.1)*ft1*ratio_c13_c14
        rprod_new(i) = xx*q100**((Tsoil(1)-TgrAnn(ns))*0.1)*ft1*(1. - ratio_c13_c14)
      elseif (i == ns) then
!        xx = &
!        & ( r0_oldorg * &
!        & exp(-alpha_old_org*C_talik_age**(-2) * &
!        & (h_talik**2-min(0.5*(z(ns)+z(ns-1)),h_talik)**2 ) ) + &
!        & r0_oldorg * &
!        & exp(-alpha_old_org*C_talik_age**(-2) * &
!        & (h_talik**2-(min(z(ns),h_talik))**2 ) ) )* 0.5
!        xx = MEAN_RPROD_OLDC(0.5*(z(ns) + z(ns-1)), z(ns), h_talik, C_talik_age)
        xx = MEAN_RPROD_OLDC2(0.5*(z(ns) + z(ns-1)), z(ns), h_talik, C_talik_age)        
!        xx = exp(-alpha_old_org*C_talik_age**(-2)*h_talik**2) * &
!        & mean_rprod_oldc0(ns)
        rprod_old(i) = xx*q100**((Tsoil(ns)-TgrAnn(ns))*0.1)*ft1*ratio_c13_c14
        rprod_new(i) = xx*q100**((Tsoil(ns)-TgrAnn(ns))*0.1)*ft1*(1. - ratio_c13_c14)
      else
!        xx = &
!        & ( r0_oldorg * &
!        & exp(-alpha_old_org*C_talik_age**(-2) * &
!        & (h_talik**2-(min(0.5*(z(i)+z(i-1)),h_talik))**2)) + &
!        & r0_oldorg * &
!        & exp(-alpha_old_org*C_talik_age**(-2) * &
!        & (h_talik**2-(min(z(i),h_talik))**2)) ) * 0.5 * &
!        & (z(i)-z(i-1))/(z(i+1)-z(i-1)) + &
!        ( r0_oldorg * &
!        & exp(-alpha_old_org*C_talik_age**(-2) * &
!        & (h_talik**2-(min(0.5*(z(i+1)+z(i)),h_talik))**2)) + &
!        & r0_oldorg * &
!        & exp(-alpha_old_org*C_talik_age**(-2) * &
!        & (h_talik**2-(min(z(i),h_talik))**2)) ) * 0.5 * &
!        & (z(i+1)-z(i))/(z(i+1)-z(i-1))
!        xx = MEAN_RPROD_OLDC(0.5*(z(i) + z(i-1)), 0.5*(z(i+1) + z(i)), h_talik, C_talik_age)
        xx = MEAN_RPROD_OLDC2(0.5*(z(i) + z(i-1)), 0.5*(z(i+1) + z(i)), h_talik, C_talik_age)        
!        xx = exp(-alpha_old_org*C_talik_age**(-2)*h_talik**2) * &
!        & mean_rprod_oldc0(i)
        rprod_old(i) = xx*q100**((Tsoil(i)-TgrAnn(ns))*0.1)*ft1*ratio_c13_c14
        rprod_new(i) = xx*q100**((Tsoil(i)-TgrAnn(ns))*0.1)*ft1*(1. - ratio_c13_c14)
      endif

      rprod_new(i) = rprod_new(i) + &
      & r0_oliglake*forg(i)*fin * &
      & q100**((Tsoil(i)-TgrAnn(ns))*0.1)*ft
      
      if (i == 1) then
        rprod_total_oldC = rprod_total_oldC + rprod_old(1)*0.5*(z(2)-z(1))
        rprod_total_newC = rprod_total_newC + rprod_new(1)*0.5*(z(2)-z(1))
      elseif (i == ns) then
        rprod_total_oldC = rprod_total_oldC + rprod_old(ns)*0.5*(z(ns)-z(ns-1))
        rprod_total_newC = rprod_total_newC + rprod_new(ns)*0.5*(z(ns)-z(ns-1))
      else
        rprod_total_oldC = rprod_total_oldC + rprod_old(i)*0.5*(z(i+1)-z(i-1))
        rprod_total_newC = rprod_total_newC + rprod_new(i)*0.5*(z(i+1)-z(i-1))
      endif
      
!      rprod(i) = rprod_new(i) + rprod_old(i)
      
!      if(qsoil(i).ge.cthresh) then
      do j = 1, 2 
        if (qsoil(i,j) .ge. ch4_crit(i,j)) then
          fcch4(i,j) = 1.
        else
          fcch4(i,j) = 0.
        end if
      enddo
      
    end do
    
!*	===============================================================
!*	the function fgrow describes the growing state of the plants
!*	===============================================================
!*	t50 - the soil temperature at 50cm depth below soil surface

    do i = 1, ns !ms, ml
      if (z(i) .le. 0.5) then 
        t50 = Tsoil(i)
      else
        exit
      endif
    end do
    
    if (t50.lt.t_grow) fgrow = rlmin
    if (t50 .ge. t_grow .and. t50 .le. t_mature) then
      fgrow = rlmin + rl*(1.-((t_mature-t50)/(t_mature-t_grow))**2)
    endif
    if (t50.gt.t_mature) fgrow = rlmax
    
!	=============================================================

    do i = 1, ns !ms, ml
!*	oxid(i) - coefficient for oxidation
!*	methane oxidation occurs only in the unsaturated soil
!*	layers
      if (wl(i).lt.anox_crit(i)) then
!		if(z(i).lt.w) then
!        oxid(i) = vmax*q10**((Tsoil(i) - TgrAnn(ns))*0.1)/(rkm + qsoil(i)) !ml

!       Currently oxidation is turned off in soil layers.
!       It should be taken into account when oxygen concentration will
!       be somehow extimated (calculated).
        oxid(i) = 0.
      else
        oxid(i) = 0.
      end if

!*	  plant(i) - coefficient for plant-mediated transport
!	  plant(i) = rkp*tveg*froot(i)*fgrow*(1.-0.5)*veg

!      plant(i) = rkp*tveg*froot(i)*fgrow*veg
      
	  plant(i) = 0.
!*
!*	  bull(i) - coefficient for ebullition
!		if(z(i).ge.w) then

      do j = 1, 2
        if (wl(i) .ge. anox_crit(i)) then
          bull(i,j) = rke*fcch4(i,j)
        else
          bull(i,j) = 0.
        end if
      enddo
      
    end do
    
!    write (*,*) bull(:,2)
!    write (*,*) rprod_old(:)
!    write (*,*) ch4_crit(:,2)
    
!    oxid_sum = 0.
!    rprod_sum = 0.
!    bull_sum = 0.
!    plant_sum = 0.
!    q_sum_old = 0.
!    do i = 2, ns-1 !ms+1, ml-1
!      q_sum_old = q_sum_old + qsoil(i)*0.5*(z(i+1)-z(i-1))
!    end do
!
!    q_sum_old = q_sum_old + qsoil(ns)*0.5*(z(ns)-z(ns-1))
!    q_sum_old = q_sum_old + qsoil(1)*(0.5*(z(2)-z(1)) + 0.5*ddz(M)*h1)

!    do i = 1, ns !1, ml
!      q_old(i) = qsoil(i,j)
!    end do

!*	=============================================================
!*	set up the tridiagonal system of equations
!*	=============================================================

    meth2 : do j = 1, 2

    if (j == 1) rprod(:) = rprod_old(:)
    if (j == 2) rprod(:) = rprod_new(:)

    if (h1 > 0) then
      if (l1 == 0) then
!        c(1) = 1.
!        b(1) = 0.
!        f(1) = ch4_pres_atm0 * &
!        & HENRY_CONST(Henry_const0_ch4, Henry_temp_dep_ch4, Henry_temp_ref, Twater(1)+273.15) 
        ! ch4_atm0 ! atmospheric concentration
        xx = HENRY_CONST(Henry_const0_ch4, Henry_temp_dep_ch4, Henry_temp_ref, Twater(1)+273.15)
        Flux_atm(j) = GAS_WATATM_FLUX &
        & (Twater(1),wind10,qwater(1,1,j),ch4_pres_atm0,xx,water_methane_indic,eps_surf)        
      else
        Flux_atm(j) = 0.
      endif
      c(1)   = - lammeth(1)/(h1*ddz(1)) + dhw0/(2.d0*dt) - & 
      & ddz(1)*h1/(2.d0*dt)
      b(1)   = - lammeth(1)/(h1*ddz(1)) + dhw0/(2.d0*dt)
      f(1)   = - ddz(1)*h1*qwater(1,1,j)/(2.d0*dt) + Flux_atm(j)
      call DIFF_COEF(a,b,c,f,2,M,2,M,water_methane_indic,dt)
      if (ls1 > 0) then
!       Case water, deepice and soil; upper ice and snow are allowed       
!       Methane diffusion in water       
        c(M+1) = - lammeth(M)/(ddz(M)*h1)-ddz(M)*h1/(2.d0*dt) + & ! bug corrected: - (dhw-dhw0)/(2.d0*dt)
        & (dhw-dhw0)/(2.d0*dt) 
        a(M+1) = - lammeth(M)/(ddz(M)*h1)+(dhw-dhw0)/(2.d0*dt)
        f(M+1) = - qwater(M+1,1,j)*ddz(M)*h1/(2.d0*dt) ! + qflux1
        call PROGONKA (a,b,c,f,y,1,M+1)
        y(1:M+1) = max(y(1:M+1),0.d0)
        qwater(1:M+1,2,j) = y(1:M+1)
!       Methane evolution in soil under bottom ice
        b(1) = 2.*dt*diff(2)/((z(2)-z(1))*(z(2)-z(1)))
        c(1) = 1. + b(1) + dt*(bull(1,j) + plant(1) + oxid(1))
!       No diffusive flux is assumed at the ice-soil interface      
        f(1) = qsoil(1,j) + dt*(rprod(1) + bull(1,j)*ch4_crit(1,j)) ! cthresh -> ch4_crit(1)
        do i = 2, ns-1 !ms+1, ml-1
          ap(i) = 2.*dt/((z(i+1)-z(i-1))*(z(i+1)-z(i)))
          bp(i) = 2.*dt/((z(i+1)-z(i-1))*(z(i)-z(i-1)))
          a(i) = bp(i) * diff(i)
          c(i) = 1. + ap(i)*diff(i+1) + bp(i)*diff(i) + &
          & dt * (bull(i,j)+plant(i)+oxid(i))
          f(i) = qsoil(i,j) + dt*rprod(i) + dt*bull(i,j)*ch4_crit(i,j) ! cthresh -> ch4_crit(i)
          b(i) = ap(i)*diff(i+1)
        end do
        bp(ns) = 2.*dt/( (z(ns)-z(ns-1))*(z(ns)-z(ns-1)) )
        a(ns) =  bp(ns) * diff(ns) ! diff(ns-1)
        c(ns) = 1. + a(ns) + & ! bp(ns)*diff(ns)
        & dt * (bull(ns,j) + plant(ns) + oxid(ns))
        f(ns) = qsoil(ns,j) + dt*(rprod(ns) + bull(ns,j)*ch4_crit(ns,j) ) ! cthresh -> ch4_crit(ns)
        call PROGONKA(a, b, c, f, y, 1, ns)
        y(1:ns) = max(y(1:ns),0.d0)
        qsoil(1:ns,j) = y(1:ns)     
      else
!       Case water, soil; upper ice and snow are allowed       
!---------------------WATER-SOIL INTERFACE------------------
        a(M+1)= lammeth(M)/(ddz(M)*h1) - (dhw-dhw0)/(2.*dt)
        b(M+1)= diff(2)/( (z(2)-z(1)) )
        c(M+1) = a(M+1) + b(M+1) + (z(2)-z(1))/(2.*dt) + &
        & ddz(M)*h1/(2.*dt) + 0.5*(z(2)-z(1))*(bull(1,j) + plant(1) + oxid(1))
        f(M+1)= qwater(M+1,1,j)*((z(2)-z(1))/(2.*dt) + & 
        & ddz(M)*h1/(2.*dt)) + 0.5*(z(2)-z(1))*(rprod(1) + bull(1,j)*ch4_crit(1,j)) ! cthresh -> ch4_crit(1)
!-----------------------------------------------------------
        do i = 2, ns-1 !ms+1, ml-1
          ap(i) = 2.*dt/((z(i+1)-z(i-1))*(z(i+1)-z(i)))
          bp(i) = 2.*dt/((z(i+1)-z(i-1))*(z(i)-z(i-1)))
          a(i+M) = bp(i) * diff(i)
          c(i+M) = 1. + ap(i)*diff(i+1) + bp(i)*diff(i) + &
          & dt * (bull(i,j)+plant(i)+oxid(i))
          f(i+M) = qsoil(i,j) + dt*rprod(i) + dt*bull(i,j)*ch4_crit(i,j) ! cthresh -> ch4_crit(i)
          b(i+M) = ap(i)*diff(i+1)
        end do
        bp(ns) = 2.*dt/( (z(ns)-z(ns-1))*(z(ns)-z(ns-1)) )
        a(M+ns) =  bp(ns) * diff(ns) ! diff(ns-1)
        c(M+ns) = 1. + a(M+ns) + & ! bp(ns)*diff(ns)
        & dt * (bull(ns,j) + plant(ns) + oxid(ns))
        f(M+ns) = qsoil(ns,j) + dt*(rprod(ns) + bull(ns,j)*ch4_crit(ns,j)) ! cthresh -> ch4_crit(ns)
        call PROGONKA (a,b,c,f,y,1,M+ns)
        y(1:M+ns) = max(y(1:M+ns),0.d0)
        qwater(1:M+1,2,j) = y(1:M+1)
        qsoil(1:ns,j) = y(M+1:M+ns)
      endif
    else
!     Case of the soil and the ice above     
!     Methane evolution in soil under the ice
      b(1) = 2.*dt*diff(2)/((z(2)-z(1))*(z(2)-z(1)))
      c(1) = 1. + b(1) + dt*(bull(1,j) + plant(1) + oxid(1))
!     No diffusive flux is assumed at the ice-soil interface
      f(1) = qsoil(1,j) + dt*(rprod(1) + bull(1,j)*ch4_crit(1,j)) ! cthresh -> ch4_crit(1)
      do i = 2, ns-1 !ms+1, ml-1
        ap(i) = 2.*dt/((z(i+1)-z(i-1))*(z(i+1)-z(i)))
        bp(i) = 2.*dt/((z(i+1)-z(i-1))*(z(i)-z(i-1)))
        a(i) = bp(i) * diff(i)
        c(i) = 1. + ap(i)*diff(i+1) + bp(i)*diff(i) + &
        & dt * (bull(i,j)+plant(i)+oxid(i))
        f(i) = qsoil(i,j) + dt*rprod(i) + dt*bull(i,j)*ch4_crit(i,j) ! cthresh -> ch4_crit(i)
        b(i) = ap(i)*diff(i+1)
      end do
      bp(ns) = 2.*dt/( (z(ns)-z(ns-1))*(z(ns)-z(ns-1)) )
      a(ns) =  bp(ns) * diff(ns) ! diff(ns-1)
      c(ns) = 1. + a(ns) + &  ! bp(ns)*diff(ns)
      & dt * (bull(ns,j) + plant(ns) + oxid(ns))
      f(ns) = qsoil(ns,j) + dt*rprod(ns) + dt*bull(ns,j)*ch4_crit(ns,j) ! cthresh -> ch4_crit(ns)
      call PROGONKA(a, b, c, f, y, 1, ns)
      y(1:ns) = max(y(1:ns),0.d0)
      qsoil(1:ns,j) = y(1:ns)         
    endif 

    enddo meth2
    

!*	============================================================
!*	the methane flux due to plant-mediated transport
!*	=================fplant(t)==================================

!    fplant = 0.
!    do i = 2, ns-1 !ms+1, ml-1
!      fplant = fplant + 0.5*(z(i+1)-z(i-1))*plant(i)*qsoil(i) * (1 - pox) !* veg
!    end do
    
!*  =============================================================
!*	the methane flux due to ebullition
!*  =================febul(t)====================================

    do j = 1, 2
      do i = 1, ns !ms, ml
        qebul(i,j) = bull(i,j)*(qsoil(i,j) - ch4_crit(i,j)) ! cthresh -> ch4_crit(i)
      end do
    enddo
    
    febul(:) = 0.
    do j = 1, 2
      do i = 2, ns-1
        if (wl(i) .gt. anox_crit(i)) then
          febul(j) = febul(j) + 0.5*(z(i+1)-z(i-1))*qebul(i,j)
        end if
      enddo
      febul(j) = febul(j) + qebul(1,j)*0.5*(z(2)-z(1))
      febul(j) = febul(j) + qebul(ns,j)*0.5*(z(ns)-z(ns-1))
    end do
    
    if (l1 /= 0.) then
      do j = 1, 2
!       A portion of methane ebullition flux (meth_ebul_ice_intercept) 
!       is intercepted by the ice cover
        febul(j) = febul(j)*(1. - meth_ebul_ice_intercept)
        tot_ice_meth_bubbles(j) = tot_ice_meth_bubbles(j) + febul(j)*meth_ebul_ice_intercept*dt
      enddo
    endif
    
    if (l1 == 0.) then
      do j = 1, 2
        if (tot_ice_meth_bubbles(j) /= 0.) then
!         When ice cover disappears, all ice-trapped bubbles are
!         released to the atmosphere instantaneously
          febul(j) = febul(j) + tot_ice_meth_bubbles(j)/dt
          tot_ice_meth_bubbles(j) = 0.
        endif
      enddo
    endif
    
!	===============================================================
!*    the diffusive flux at the soil/water boundary, downwards
!*	=================fdiff(t, dzeta=1)==================================

!    fdiff = - (diff(1) + diff(2))*0.5*(qsoil(1)-qsoil(2))/ & ! ms
!    &         (z(2)-z(1))

!    fdiff = 0.5*h1*ddz(M)/dt*(qwater(M+1,2)-qwater(M+1,1)) - &
!    & lammeth(M)*(qwater(M+1,2)-qwater(M,2))/(h1*ddz(M)) + &
!    & (dhw/dt - dhw0/dt)*0.5*(qwater(M+1,2)-qwater(M,2))
    
    fdiff(1:2) = - lammeth(M)*(qwater(M+1,2,1:2)-qwater(M,2,1:2))/(h1*ddz(M))
    
!    fdiff1 = - diff(ns)*(qsoil(ns)-qsoil(ns-1))/ & ! ml
!    &          (z(ns)-z(ns-1))
    
!	===============================================================
!*	total methane emission ftot
!*  ==============================================================
!   if(w.ge.rns)then

!    ftot = fdiff + febul + fplant
    
    fdiff_lake_surf(:) = - Flux_atm(:) !- lammeth(1)/(h1*ddz(1))*(qwater(2,2)-qwater(1,2))

!    q_sum_new = 0.
!    do i = 2, ns-1 !ms+1, ml-1
!      q_sum_new = q_sum_new + qsoil(i)*0.5*(z(i+1)-z(i-1))
!    end do
!    q_sum_new = q_sum_new + qsoil(ns)*0.5*(z(ns)-z(ns-1)) ! ml
!    q_sum_new = q_sum_new + qsoil(1)*(0.5*(z(2)-z(1)) + 0.5*ddz(M)*h1)

!    do i = 2, ns-1 !ms+1, ml-1
!      rprod_sum = rprod_sum + rprod(i)*dt*0.5*(z(i+1)-z(i-1))
!!      xx = q_old(i)*0.5*(z(i+1)-z(i-1))+ rprod(i)*dt*0.5*(z(i+1)-z(i-1))
!!      oxid_sum = oxid_sum + &
!!      & min(oxid(i)*qsoil(i)*0.5*(z(i+1)-z(i-1)), xx)
!!      bull_sum = bull_sum + &
!!      & min(bull(i)*(qsoil(i)-ch4_crit(i))*dt*0.5*(z(i+1)-z(i-1)), xx) ! cthresh -> ch4_crit(i)
!!      plant_sum = plant_sum + &
!!      & min(plant(i)*qsoil(i)*dt*0.5*(z(i+1)-z(i-1)), xx)
!      oxid_sum = oxid_sum + oxid(i)*qsoil(i)*dt*0.5*(z(i+1)-z(i-1))
!      bull_sum = bull_sum + bull(i)*(qsoil(i)-ch4_crit(i))*dt*0.5*(z(i+1)-z(i-1))
!      plant_sum = plant_sum + plant(i)*qsoil(i)*dt*0.5*(z(i+1)-z(i-1))
!    end do
    
!    rprod_sum = rprod_sum + rprod(ns)*dt*0.5*(z(ns)-z(ns-1)) ! ml
!!    xx = q_old(ns)*0.5*(z(ns)-z(ns-1))+ rprod(ns)*dt*0.5*(z(ns)-z(ns-1))
!!    oxid_sum = oxid_sum + &
!!    & min(oxid(ns)*qsoil(ns)*dt*0.5*(z(ns)-z(ns-1)), xx)
!!    bull_sum = bull_sum + &
!!    & min(bull(ns)*(qsoil(ns)-ch4_crit(ns))*dt*0.5*(z(ns)-z(ns-1)), xx) ! cthresh -> ch4_crit(ns)
!!    plant_sum = plant_sum + &
!!    & min(plant(ns)*qsoil(ns)*dt*0.5*(z(ns)-z(ns-1)), xx)
!    oxid_sum = oxid_sum + oxid(ns)*qsoil(ns)*dt*0.5*(z(ns)-z(ns-1))
!    bull_sum = bull_sum + bull(ns)*(qsoil(ns)-ch4_crit(ns))*dt*0.5*(z(ns)-z(ns-1))
!    plant_sum = plant_sum + plant(ns)*qsoil(ns)*dt*0.5*(z(ns)-z(ns-1))
!    
!    rprod_sum = rprod_sum + rprod(1)*dt*0.5*(z(2)-z(1))
!!    xx = q_old(1)*0.5*(z(2)-z(1))+ rprod(1)*dt*0.5*(z(2)-z(1))
!!    oxid_sum = oxid_sum + &
!!    & min(oxid(1)*qsoil(1)*dt*0.5*(z(2)-z(1)), xx)
!!    bull_sum = bull_sum + &
!!    & min(bull(1)*(qsoil(1)-ch4_crit(1))*dt*0.5*(z(2)-z(1)), xx) ! cthresh -> ch4_crit(1)
!!    plant_sum = plant_sum + &
!!    & min(plant(1)*qsoil(1)*dt*0.5*(z(2)-z(1)), xx)
!    oxid_sum = oxid_sum + oxid(1)*qsoil(1)*dt*0.5*(z(2)-z(1))
!    bull_sum = bull_sum + bull(1)*(qsoil(1)-ch4_crit(1))*dt*0.5*(z(2)-z(1))
!    plant_sum = plant_sum + plant(1)*qsoil(1)*dt*0.5*(z(2)-z(1))
    
!    print*, 'bounds_talik', 0.5*(z(i_talik-1)+z(i_talik)), 0.5*(z(i_talik+1)+z(i_talik))

!    if((q_sum_new-q_sum_old).ne.0.) then 
!!     Check the conservation of the scheme
!      balance = (rprod_sum - oxid_sum - plant_sum - &
!      &  bull_sum + fdiff*dt - & !fdiff1*dt-
!      &  (q_sum_new-q_sum_old))/(q_sum_new-q_sum_old)*100. ! procents
!    endif
!    print*, 'meth_test', qsoil-q_old

!    write(*,*) qwater(:,1,1)
!    write(*,*) qwater(:,1,2)

    deallocate (z)
    deallocate (roots)
    deallocate (a,b,c,f,y)

    return
    END SUBROUTINE METHANE2
    
    SUBROUTINE METHANE_OXIDATION2 &
    & (M, dt, oxyg, qwater1, qwater2, DIC, Tw)
    
    use METH_OXYG_CONSTANTS
    use PHYS_FUNC, only : &
    & REACPOT_ARRHEN 
        
    implicit none
    
!   Input/output variables    
    integer(kind=iintegers), intent(in) :: M
    
    real(kind=ireals), intent(in) :: dt
    
    real(kind=ireals), intent(inout) :: oxyg(1:M+1)
    real(kind=ireals), intent(inout) :: qwater1(1:M+1)
    real(kind=ireals), intent(inout) :: qwater2(1:M+1)
    real(kind=ireals), intent(inout) :: DIC(1:M+1)
    real(kind=ireals), intent(in) :: Tw(1:M+1)
    
!   Local variables and constants

    real(kind=ireals), parameter :: Kelvin0 = 273.15
    real(kind=ireals), parameter :: alpha = 0.1

    integer(kind=iintegers) :: i, j
    real(kind=ireals) :: x(1:2), xx, temp_K ! Work variables
    real(kind=ireals) :: qwater_old(1:2), oxyg_old
    
!    write(*,*) qwater1, qwater2
    
    do i = 2, M
      temp_K = Tw(i) + Kelvin0
      x(1) = dt*REACPOT_ARRHEN(delta_Eq,temp_K,temp0,vq_max) &
      & /(k_o2 + oxyg(i))/(k_ch4 + qwater1(i))
      x(2) = dt*REACPOT_ARRHEN(delta_Eq,temp_K,temp0,vq_max) &
      & /(k_o2 + oxyg(i))/(k_ch4 + qwater2(i))
      qwater_old(1) = qwater1(i)
      qwater_old(2) = qwater2(i)
      oxyg_old = oxyg(i)
      do j = 1, 100
!       Simple iterations to solve the system of three non-linear equations
        xx = oxyg(i)
        oxyg(i) = oxyg(i) - &
        & alpha*(oxyg(i) - (oxyg_old - 2.*oxyg(i)*(x(1)*qwater1(i) + x(2)*qwater2(i) ) ) )
!        qwater1(i) = qwater_old(1) - x(1)*xx*qwater1(i)
        qwater1(i) = qwater1(i) - &
        & alpha*(qwater1(i) - (qwater_old(1) - x(1)*xx*qwater1(i)))
!        qwater2(i) = qwater_old(2) - x(2)*xx*qwater2(i) 
        qwater2(i) = qwater2(i) - &
        & alpha*(qwater2(i) - (qwater_old(2) - x(2)*xx*qwater2(i)))
!        write(*,*) oxyg(i) - (oxyg_old - 2.*oxyg(i)*(x(1)*qwater1(i) + x(2)*qwater2(i))), &
!        & qwater1(i) - (qwater_old(1) - x(1)*xx*qwater1(i)), &
!        & qwater2(i) - (qwater_old(2) - x(2)*xx*qwater2(i)) !, &
!        & oxyg(i), qwater1(i), qwater2(i)
      enddo
      DIC(i) = DIC(i) - 0.5d0*(oxyg(i) - oxyg_old)
    enddo  
    
    END SUBROUTINE METHANE_OXIDATION2    

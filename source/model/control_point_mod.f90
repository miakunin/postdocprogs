MODULE CONTROL_POINT

use LAKE_DATATYPES, only : ireals, iintegers
use COMPARAMS, only : parallel_params

contains
SUBROUTINE CONTROL_POINT_OUT(nx,nx0,nx_max,ny,ny0,ny_max,gs,parparams)

!Writes the control point of LAKE model

use EVOLUTION_VARIABLES, only:  &
!Parallel MPI output from all processes in single file
& l1_2d,h1_2d,hx1_2d,hx2_2d,hy1_2d,hy2_2d,ls1_2d,hs1_2d, &
! Thermocline displacement
& hx1t_2d,hx2t_2d,hy1t_2d,hy2t_2d, &
! Displacements of each layer
& hx1ml_2d,hx2ml_2d,hy1ml_2d,hy2ml_2d, & 
& u_2d,v_2d, &
& E_2d, eps_2d, k_turb_T_flux_2d, &
& Tsoil1_2d,Sals1_2d,wi1_2d,wl1_2d,qsoil_2d, &
& Tw1_2d, Sal1_2d, qwater_2d, oxyg_2d, carbdi_2d,&
& oxygsoil_2d, Tskin_2d, Ti1_2d, Tis1_2d, &
& ueffl_2d, veffl_2d, Tweffl_2d, Saleffl_2d, &
& fl_sn_2d, fl_sn_init_2d, &
& itop_2d, &
& dz_2d, T_2d, wl_2d, dens_2d, &
& snmelt_2d, cdm_2d, &
& time_2d, nstep_2d, i_maxN_2d, &
& dhwfsoil_2d, Elatent_2d, dhw_2d, dhw0_2d, dhi_2d, dhi0_2d, dls0_2d, &
& velfrict_2d, roughness_2d, eflux0_kinem_2d, &
& tot_ice_meth_bubbles_2d, febul0_2d, &
& Eseiches_2d, lamw_2d

use INOUT, only : wrfild_2d_coords, wrfild_3d_coords, &
& WRFILD_2D_LAKE, WRFILD_3D_LAKE

use ARRAYS_GRID, only : gridsize_type

use INOUT_PARAMETERS, only : lake_misc_unit_min, lake_misc_unit_max

implicit none

!Input variables
integer(kind=iintegers), intent(in) :: nx, nx0, ny, ny0, nx_max, ny_max
type(gridsize_type), intent(in) :: gs
type(parallel_params), intent(in) :: parparams

!Output variables

!Local variables
type(wrfild_2d_coords) :: coords
type(wrfild_3d_coords) :: coords3d
integer(kind=iintegers) :: unitout, dnx, dny, recn, nzmax, i, j, k
character(len=40) :: filen

real(kind=ireals), allocatable :: work(:,:,:)

!nzmax = max(gs%M+1,gs%Mice+1,gs%ms,gs%ml,gs%ns*gs%nsoilcols)
!allocate (work(dnx,dny,nzmax))

dnx = nx - nx0 + 1
dny = ny - ny0 + 1

!2D arrays
coords%i0 = nx0
coords%i1 = nx
coords%j0 = ny0
coords%j1 = ny 
coords%i00 = nx0
coords%i11 = nx 
coords%j00 = ny0 
coords%j11 = ny
coords%i0_write = 1
coords%i1_write = nx_max 
coords%j0_write = 1
coords%j1_write = ny_max 
coords%isx = 1
coords%isy = 1 
coords%lrec = 4 ! bytes per real value
coords%recn = 1

filen = './results/control_point/cp2d.bin'
unitout = lake_misc_unit_min

!Parallel MPI output from all processes in single file
call WRFILD_2D_LAKE(l1_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(h1_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hx1_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hx2_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hy1_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hy2_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(ls1_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hs1_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hx1t_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hx2t_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hy1t_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(hy2t_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
if (allocated(hx1ml_2d)) then
  call WRFILD_2D_LAKE(hx1ml_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
  call WRFILD_2D_LAKE(hx2ml_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
  call WRFILD_2D_LAKE(hy1ml_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
  call WRFILD_2D_LAKE(hy2ml_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
endif

call WRFILD_2D_LAKE(Tskin_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(snmelt_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(cdm_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(time_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(dhwfsoil_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(Elatent_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(dhw_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(dhw0_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(dhi_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(dhi0_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(dls0_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1

call WRFILD_2D_LAKE(velfrict_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(roughness_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(eflux0_kinem_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(tot_ice_meth_bubbles_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(Eseiches_2d, coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1


!2D integers arrays
allocate (work(dnx,dny,5))
work(1:dnx,1:dny,1) = real(fl_sn_2d     (1:dnx,1:dny),kind=ireals)
work(1:dnx,1:dny,2) = real(fl_sn_init_2d(1:dnx,1:dny),kind=ireals)
work(1:dnx,1:dny,3) = real(itop_2d      (1:dnx,1:dny),kind=ireals)
work(1:dnx,1:dny,4) = real(nstep_2d     (1:dnx,1:dny),kind=ireals)
work(1:dnx,1:dny,5) = real(i_maxN_2d    (1:dnx,1:dny),kind=ireals)
call WRFILD_2D_LAKE(work(1,1,1), coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(work(1,1,2), coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(work(1,1,3), coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(work(1,1,4), coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
call WRFILD_2D_LAKE(work(1,1,5), coords, unitout, .false., .true., parparams, filen); coords%recn = coords%recn + 1
deallocate(work)

!3D arrays
coords3d%i0 = nx0
coords3d%i1 = nx
coords3d%j0 = ny0
coords3d%j1 = ny 
coords3d%i00 = nx0
coords3d%i11 = nx 
coords3d%j00 = ny0 
coords3d%j11 = ny
coords3d%i0_write = 1
coords3d%i1_write = nx_max 
coords3d%j0_write = 1
coords3d%j1_write = ny_max 
coords3d%isx = 1
coords3d%isy = 1 
coords3d%iss = 1 
coords3d%lrec = 4 ! bytes per real value
recn = 0

filen = './results/control_point/cp3d.bin'
unitout = lake_misc_unit_min

!In-water arrays
coords3d%k0 = 1
coords3d%k1 = gs%M+1
coords3d%k00 = 1
coords3d%k11 = gs%M+1
coords3d%k0_write = 1
coords3d%k1_write = gs%M+1
allocate (work(dnx,dny,gs%M+1))
call SHIFTVI(gs%M+1,dnx,dny,u_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,v_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,E_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,eps_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,k_turb_T_flux_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,Tw1_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,Sal1_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,qwater_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,oxyg_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,carbdi_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,lamw_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,ueffl_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,veffl_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,Tweffl_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%M+1,dnx,dny,Saleffl_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
deallocate(work)

!In-ice arrays
coords3d%k0 = 1
coords3d%k1 = gs%Mice+1
coords3d%k00 = 1
coords3d%k11 = gs%Mice+1
coords3d%k0_write = 1
coords3d%k1_write = gs%Mice+1
allocate (work(dnx,dny,gs%Mice+1))
call SHIFTVI(gs%Mice+1,dnx,dny,Ti1_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%Mice+1,dnx,dny,Tis1_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
deallocate(work)

!In-snow arrays
coords3d%k0 = 1
coords3d%k1 = gs%ms
coords3d%k00 = 1
coords3d%k11 = gs%ms
coords3d%k0_write = 1
coords3d%k1_write = gs%ms
allocate (work(dnx,dny,gs%ms))
call SHIFTVI(gs%ms,dnx,dny,dz_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%ms,dnx,dny,dens_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
deallocate (work)
coords3d%k0 = 1
coords3d%k1 = gs%ml
coords3d%k00 = 1
coords3d%k11 = gs%ml
coords3d%k0_write = 1
coords3d%k1_write = gs%ml
allocate (work(dnx,dny,gs%ml))
call SHIFTVI(gs%ml,dnx,dny,T_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%ml,dnx,dny,wl_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
deallocate (work)

!Sediments surface
coords3d%k0 = 1
coords3d%k1 = gs%nsoilcols
coords3d%k00 = 1
coords3d%k11 = gs%nsoilcols
coords3d%k0_write = 1
coords3d%k1_write = gs%nsoilcols
allocate (work(dnx,dny,gs%nsoilcols))
call SHIFTVI(gs%nsoilcols,dnx,dny,febul0_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
call SHIFTVI(gs%nsoilcols,dnx,dny,oxygsoil_2d,work,1_iintegers)
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
deallocate (work)

!Soil columns
coords3d%k0 = 1
coords3d%k1 = gs%ns*gs%nsoilcols
coords3d%k00 = 1
coords3d%k11 = gs%ns*gs%nsoilcols
coords3d%k0_write = 1
coords3d%k1_write = gs%ns*gs%nsoilcols
allocate (work(dnx,dny,gs%ns*gs%nsoilcols))
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    work(1:dnx,1:dny,k) = Tsoil1_2d(i,j,1:dnx,1:dny)
  enddo
enddo
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    work(1:dnx,1:dny,k) = Sals1_2d(i,j,1:dnx,1:dny)
  enddo
enddo
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    work(1:dnx,1:dny,k) = wi1_2d(i,j,1:dnx,1:dny)
  enddo
enddo
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    work(1:dnx,1:dny,k) = wl1_2d(i,j,1:dnx,1:dny)
  enddo
enddo
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    work(1:dnx,1:dny,k) = qsoil_2d(i,j,1:dnx,1:dny)
  enddo
enddo
call WRFILD_3D_LAKE(work, coords3d, recn, unitout, .false., .true., parparams, filen)
deallocate(work)

END SUBROUTINE CONTROL_POINT_OUT


SUBROUTINE CONTROL_POINT_IN(nx,nx0,nx_max,ny,ny0,ny_max,gs,parparams)

!Writes the control point of LAKE model

use EVOLUTION_VARIABLES, only: ALLOCATE_ARRAYS, &
!Parallel MPI output from all processes in single file
& l1_2d,h1_2d,hx1_2d,hx2_2d,hy1_2d,hy2_2d,ls1_2d,hs1_2d, &
! Thermocline displacement
& hx1t_2d,hx2t_2d,hy1t_2d,hy2t_2d, &
! Displacements of each layer
& hx1ml_2d,hx2ml_2d,hy1ml_2d,hy2ml_2d, & 
& u_2d,v_2d, &
& E_2d, eps_2d, k_turb_T_flux_2d, &
& Tsoil1_2d,Sals1_2d,wi1_2d,wl1_2d,qsoil_2d, &
& Tw1_2d, Sal1_2d, qwater_2d, oxyg_2d, carbdi_2d,&
& oxygsoil_2d, Tskin_2d, Ti1_2d, Tis1_2d, &
& ueffl_2d, veffl_2d, Tweffl_2d, Saleffl_2d, &
& fl_sn_2d, fl_sn_init_2d, &
& itop_2d, &
& dz_2d, T_2d, wl_2d, dens_2d, &
& snmelt_2d, cdm_2d, &
& time_2d, nstep_2d, i_maxN_2d, &
& dhwfsoil_2d, Elatent_2d, dhw_2d, dhw0_2d, dhi_2d, dhi0_2d, dls0_2d, &
& velfrict_2d, roughness_2d, eflux0_kinem_2d, &
& tot_ice_meth_bubbles_2d, febul0_2d, &
& Eseiches_2d, lamw_2d

use INOUT, only : wrfild_3d_coords, RDFILD_3D_LAKE

use ARRAYS_GRID, only : gridsize_type

use INOUT_PARAMETERS, only : lake_misc_unit_min, lake_misc_unit_max

implicit none

!Input variables
integer(kind=iintegers), intent(in) :: nx, nx0, ny, ny0, nx_max, ny_max
type(gridsize_type), intent(in) :: gs
type(parallel_params), intent(in) :: parparams

!Output variables

!Local variables
type(wrfild_3d_coords) :: coords3d
integer(kind=iintegers) :: unitout, dnx, dny, recn, nzmax, i, j, k
character(len=40) :: filen

real(kind=ireals), allocatable :: work(:,:,:)

!nzmax = max(gs%M+1,gs%Mice+1,gs%ms,gs%ml,gs%ns*gs%nsoilcols)
!allocate (work(dnx,dny,nzmax))

dnx = nx - nx0 + 1
dny = ny - ny0 + 1

if (.not. allocated(l1_2d)) &
& call ALLOCATE_ARRAYS(dnx, dny, gs%M, gs%Mice, gs%ns, gs%ms, gs%ml, gs%nsoilcols)

!2D arrays
coords3d%i0 = nx0
coords3d%i1 = nx
coords3d%j0 = ny0
coords3d%j1 = ny 
coords3d%k0 = 1
coords3d%k1 = 1
coords3d%i00 = nx0
coords3d%i11 = nx 
coords3d%j00 = ny0 
coords3d%j11 = ny
coords3d%k00 = 1
coords3d%k11 = 1
coords3d%i0_write = 1
coords3d%i1_write = nx_max 
coords3d%j0_write = 1
coords3d%j1_write = ny_max 
coords3d%k0_write = 1
coords3d%k1_write = 1
coords3d%isx = 1
coords3d%isy = 1 
coords3d%iss = 1
coords3d%lrec = 4 ! bytes per real value

filen = './results/control_point/cp2d.bin'
unitout = lake_misc_unit_min
recn = 0

!Parallel MPI output from all processes in single file

call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, l1_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, h1_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hx1_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hx2_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hy1_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hy2_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, ls1_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hs1_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hx1t_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hx2t_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hy1t_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hy2t_2d)
if (allocated(hx1ml_2d)) then
  call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hx1ml_2d)
  call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hx2ml_2d)
  call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hy1ml_2d)
  call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, hy2ml_2d)
endif

call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, Tskin_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, snmelt_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, cdm_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, time_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, dhwfsoil_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, Elatent_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, dhw_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, dhw0_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, dhi_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, dhi0_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, dls0_2d)

call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, velfrict_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, roughness_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, eflux0_kinem_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, tot_ice_meth_bubbles_2d)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, Eseiches_2d)

!2D integers arrays
allocate(work(dnx,dny,5))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work(1,1,1))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work(1,1,2))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work(1,1,3))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work(1,1,4))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work(1,1,5))
fl_sn_2d     (1:dnx,1:dny) = int(work(1:dnx,1:dny,1))
fl_sn_init_2d(1:dnx,1:dny) = int(work(1:dnx,1:dny,2))
itop_2d      (1:dnx,1:dny) = int(work(1:dnx,1:dny,3))
nstep_2d     (1:dnx,1:dny) = int(work(1:dnx,1:dny,4))
i_maxN_2d    (1:dnx,1:dny) = int(work(1:dnx,1:dny,5))
deallocate(work)

!3D arrays
coords3d%i0 = nx0
coords3d%i1 = nx
coords3d%j0 = ny0
coords3d%j1 = ny 
coords3d%i00 = nx0
coords3d%i11 = nx 
coords3d%j00 = ny0 
coords3d%j11 = ny
coords3d%i0_write = 1
coords3d%i1_write = nx_max 
coords3d%j0_write = 1
coords3d%j1_write = ny_max 
coords3d%isx = 1
coords3d%isy = 1 
coords3d%iss = 1 
coords3d%lrec = 4 ! bytes per real value
recn = 0

filen = './results/control_point/cp3d.bin'
unitout = lake_misc_unit_min

!In-water arrays
coords3d%k0 = 1
coords3d%k1 = gs%M+1
coords3d%k00 = 1
coords3d%k11 = gs%M+1
coords3d%k0_write = 1
coords3d%k1_write = gs%M+1
allocate (work(dnx,dny,gs%M+1))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,u_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,v_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,E_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,eps_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,k_turb_T_flux_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,Tw1_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,Sal1_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,qwater_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,oxyg_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,carbdi_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,lamw_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,ueffl_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,veffl_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,Tweffl_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%M+1,dnx,dny,Saleffl_2d,work,2_iintegers)
deallocate(work)

!In-ice arrays
coords3d%k0 = 1
coords3d%k1 = gs%Mice+1
coords3d%k00 = 1
coords3d%k11 = gs%Mice+1
coords3d%k0_write = 1
coords3d%k1_write = gs%Mice+1
allocate (work(dnx,dny,gs%Mice+1))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%Mice+1,dnx,dny,Ti1_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%Mice+1,dnx,dny,Tis1_2d,work,2_iintegers)
deallocate(work)

!In-snow arrays
coords3d%k0 = 1
coords3d%k1 = gs%ms
coords3d%k00 = 1
coords3d%k11 = gs%ms
coords3d%k0_write = 1
coords3d%k1_write = gs%ms
allocate (work(dnx,dny,gs%ms))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%ms,dnx,dny,dz_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%ms,dnx,dny,dens_2d,work,2_iintegers)
deallocate (work)
coords3d%k0 = 1
coords3d%k1 = gs%ml
coords3d%k00 = 1
coords3d%k11 = gs%ml
coords3d%k0_write = 1
coords3d%k1_write = gs%ml
allocate (work(dnx,dny,gs%ml))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%ml,dnx,dny,T_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%ml,dnx,dny,wl_2d,work,2_iintegers)
deallocate (work)

!Sediments surface
coords3d%k0 = 1
coords3d%k1 = gs%nsoilcols
coords3d%k00 = 1
coords3d%k11 = gs%nsoilcols
coords3d%k0_write = 1
coords3d%k1_write = gs%nsoilcols
allocate (work(dnx,dny,gs%nsoilcols))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%nsoilcols,dnx,dny,febul0_2d,work,2_iintegers)
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
call SHIFTVI(gs%nsoilcols,dnx,dny,oxygsoil_2d,work,2_iintegers)
deallocate (work)

!Soil columns
coords3d%k0 = 1
coords3d%k1 = gs%ns*gs%nsoilcols
coords3d%k00 = 1
coords3d%k11 = gs%ns*gs%nsoilcols
coords3d%k0_write = 1
coords3d%k1_write = gs%ns*gs%nsoilcols
allocate (work(dnx,dny,gs%ns*gs%nsoilcols))
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    Tsoil1_2d(i,j,1:dnx,1:dny) = work(1:dnx,1:dny,k)
  enddo
enddo
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    Sals1_2d(i,j,1:dnx,1:dny) = work(1:dnx,1:dny,k)
  enddo
enddo
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    wi1_2d(i,j,1:dnx,1:dny) = work(1:dnx,1:dny,k)
  enddo
enddo
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    wl1_2d(i,j,1:dnx,1:dny) = work(1:dnx,1:dny,k)
  enddo
enddo
call RDFILD_3D_LAKE(coords3d, recn, unitout, .false., .true., parparams, filen, work)
k = 0
do i = 1, gs%ns
  do j = 1, gs%nsoilcols
    k = k + 1
    qsoil_2d(i,j,1:dnx,1:dny) = work(1:dnx,1:dny,k)
  enddo
enddo
deallocate(work)

END SUBROUTINE CONTROL_POINT_IN


SUBROUTINE SHIFTVI(i1,j1,k1,f1,f2,ind)

implicit none

!Input variables
integer(kind=iintegers), intent(in) :: i1, j1, k1, ind
real(kind=ireals), intent(inout) :: f1(i1,j1,k1)

!Output variables
real(kind=ireals), intent(inout) :: f2(j1,k1,i1)

!Locals
integer(kind=iintegers) :: i

if (ind == 1) then
  do i = 1, i1
    f2(1:j1,1:k1,i) = f1(i,1:j1,1:k1)
  enddo
elseif (ind == 2) then
  do i = 1, i1
    f1(i,1:j1,1:k1) = f2(1:j1,1:k1,i)
  enddo
endif

END SUBROUTINE SHIFTVI


END MODULE CONTROL_POINT

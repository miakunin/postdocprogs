MODULE CONVECTPAR


INTERFACE
SUBROUTINE WaterMix(M,flag, T1, T2, ro1)

use LAKE_DATATYPES, only : ireals, iintegers
! Input variables
integer(kind=iintegers), intent(in)    ::    flag , M
! Output variables
real(kind=ireals)   , intent(inout) ::    T1(:), T2(:), ro1(:) 

END SUBROUTINE WaterMix
END INTERFACE

END MODULE CONVECTPAR

SUBROUTINE WaterMix(M,flag, T1, T2, ro1)

use water_density, only: &
& water_dens_ts, water_dens_ts_kivu

use LAKE_DATATYPES, only : ireals, iintegers
implicit none

! Input variables
integer(kind=iintegers), intent(in)    ::    flag , M
! Output variables
real(kind=ireals)   , intent(inout) ::    T1(:), T2(:), ro1(:)

real(kind=ireals)                   ::    ro2(1:SIZE(ro1))

integer(kind=iintegers) n, i, k, p
logical    firstcall

data firstcall /.true./

!SAVE
                            
do i = 1, M+1
  ro1(i) = WATER_DENS_TS_KIVU(T1(i),0.e0_ireals)
  ro2(i) = ro1(i)
end do    
                       
do i = 1, M+1
  T2(i) = T1(i)
end do        
                   
!if (flag==0) then                
! do n = 1, M-1 
!  k = 1
!  do i = 1, M+1-n
!   if (ro1(i) > ro1(k)) then
!    k = i 
!   end if
!  end do
!  T2(M+1-n) = T1(k)
!  T2(k) = T1(M+1-n)
!  ro2(M+1-n) = ro1(k)
!  ro2(k) = ro1(M+1-n) 
!  do p = 1, M+1-n
!  ro1(p) = ro2(p)
!   T1(p) = T2(p)
!  end do
! end do
!else       

do n = 1, M-2 
  k = 2
    do i = 2, M+1-n
      if (ro1(i) > ro1(k)) then
        k = i 
      end if
    end do
  T2(M+1-n) = T1(k)
  T2(k) = T1(M+1-n)
  ro2(M+1-n) = ro1(k)
  ro2(k) = ro1(M+1-n) 
  do p = 1, M+1-n
    ro1(p) = ro2(p)
    T1(p) = T2(p)
  end do
end do      

!endif
      
if (flag==0) then
  ro1(1) = WATER_DENS_TS_KIVU(0.5*(T1(1)+T1(2)),0.e0_ireals)
!800.969d-7+588.194d-7*(T1(1)+T1(2))/2.- 
!& 811.465d-8*((T1(1)+T1(2))**2)/4.+476.600d-10*
!& ((T1(1)+T1(2))**3)/8.
  ro2(1) = ro1(1)
   c10:do i=M,2,-1
       if (ro2(1)>ro2(i)) then
         T2(i)=(T1(1)+T1(2))/2
         ro2(i)=ro1(1)
         do k=2,i
           T2(k-1)=T1(k)
           ro2(k-1)=ro1(k)
         enddo
         T1=T2
         ro1=ro2
         exit c10
       endif
     enddo c10   
endif

if (firstcall) firstcall=.false.                                       
END SUBROUTINE WaterMix 

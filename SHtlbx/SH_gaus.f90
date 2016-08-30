!!subroutine to calculate gaussian filter coefficients and apply to spherical harmonic coefficients
!!gaussian filter is estimated in the spherical domain as an exponential (according to chambers 2004):
!! Wl = exp (- 0.25 * l^2 * ((Radius/(Rearth))^2.)/ln 2

!!exponentials for a certain radius are only computed once and then stored for reuse at the next call
!!lmax is an optional argument which sets the precomputation of the filter to either a higher degree (for subsequent calls of higher degree SH)
!!coded by Roelof Rietbroek 23 - 5 - 2007

subroutine SH_gaus(clm,radius)
use shtlbx,only:SH_lm,RE,SH_pos
!use modules with constants interfaces
implicit none
double precision,dimension(:), intent(inout) :: clm
double precision, intent(in)::radius
double precision, save:: radold = 0.d0
double precision, save, allocatable, dimension(:) :: wl
double precision :: arg
logical, save :: first = .true.
integer	::l,m,st,nd,lreq
integer,save :: lmax =0

!get maximum degree and order from input vector
call SH_lm(size(clm,1),lreq,m)


!check if it is the first run or whether a new radius is requested or whether the supported lmax is not sufficient
!begin precomputing
if (first .or. abs(radold-radius)> 0.1 .or. (lmax < lreq)) then
!set lmax to degree corresponding to vector length
lmax=lreq

if (first) then
allocate(wl(0:lmax))

else
deallocate(wl)
allocate (wl(0:lmax))
end if


first = .false.
radold = radius
wl =0.d0

!calculate exponentials for smoothing
do,l=0,lreq
  arg = 0.25d0*dble(l*l)*((radius/(RE))**2)/log(2.0)
  wl(l) = exp(-arg)
 ! write(*,*)wl
end do
!write(*,*)'using Gaussian smoothing with halfwidth (km)',int(radius/1000.d0)
end if !end of precomputing


!here the actual multiplication with the SH vector is performed
!write(*,*)'test'
do, l=0,lreq
st=SH_pos(l,0)
nd=SH_pos(l,l)
!write(*,*)st,nd,wl(l)
clm(st:nd)=wl(l)*clm(st:nd)
end do


end subroutine SH_gaus

!!Coded by Roelof Rietbroek, Mon Apr 16 13:41:41 2012
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de
!! added the recursion version ( wahr 1998) to this file

subroutine SH_gaus_recur(clm,radius)
use shtlbx,only:SH_lm,RE,SH_pos,pi
!use modules with constants interfaces
implicit none
double precision,dimension(:), intent(inout) :: clm
double precision, intent(in)::radius
double precision, save:: radold = 0.d0
double precision, save, allocatable, dimension(:) :: wl
double precision :: arg,exparg
logical, save :: first = .true.
integer	::l,m,st,nd,lreq
integer,save :: lmax =0

!get maximum degree and order from input vector
call SH_lm(size(clm,1),lreq,m)


!check if it is the first run or whether a new radius is requested or whether the supported lmax is not sufficient
!begin precomputing
if (first .or. abs(radold-radius)> 0.1 .or. (lmax < lreq)) then
   !set lmax to degree corresponding to vector length
   lmax=lreq

   if (first) then
      allocate(wl(0:lmax))
      
   else
      deallocate(wl)
      allocate (wl(0:lmax))
   end if
   
   
   first = .false.
   radold = radius
   arg=log(2.d0)/(1-cos(radius/RE))
   exparg=exp(-2*arg)
   wl=0
   wl(0)=1!/(2*pi)
   wl(1)=wl(0)*((1+exparg)/((1-exparg)) -1/arg)

   !calculate filter coefficients for smoothing
   do,l=1,lreq-1
      wl(l+1) = -(2*l+1)/arg*wl(l)+wl(l-1)
      if(wl(l+1) < 1d-7)exit !small enough
   end do
   !write(*,*)'using Gaussian smoothing with halfwidth (km)',int(radius/1000.d0)
end if !end of precomputing


!here the actual multiplication with the SH vector is performed
!write(*,*)'test'
do, l=0,lreq
   st=SH_pos(l,0)
   nd=SH_pos(l,l)
   !write(*,*)st,nd,wl(l)
   clm(st:nd)=wl(l)*clm(st:nd)
end do


end subroutine SH_gaus_recur


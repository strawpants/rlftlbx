!!Coded by Roelof Rietbroek, Mon Jul  1 20:17:25 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!module which  performs spherical harmonic synthesis on a geographical grid


module SH_synthesis
implicit none

type SHsynth_struct
integer::lmax=0
integer::nlon=0
integer::nlat=0
double precision,pointer,dimension(:,:,:)::trig=>null()
double precision,pointer,dimension(:,:)::Pnm=>null()
double precision,pointer,dimension(:)::lon=>null()
double precision,pointer,dimension(:)::lat=>null()
logical::save_pnm=.true. ! logical which determines whether Legendre functions are precomputed and stored
end type SHsynth_struct

contains
subroutine SH_synth_init(lmax,lon,lat,SHS,save_pnm)
use shtlbx
use shtools
implicit none
integer,intent(in)::lmax
double precision,intent(in),dimension(:)::lon,lat
type(SHsynth_struct),intent(inout)::SHS
logical,intent(in)::save_pnm

integer::nsh,n,m,i
double precision::d2r

d2r=pi/180.d0

SHS%nlon=size(lon,1)
SHS%nlat=size(lat,1)
allocate(SHS%lon(SHS%nlon),SHS%lat(SHS%nlat))

SHS%lat=lat*d2r
SHS%lon=lon*d2r
SHS%lmax=lmax

nsh=SH_pos(lmax,lmax)

allocate(SHS%trig(SHS%nlon,0:2*lmax,2))
SHS%trig=0.d0
!precompute longitude dependent part
do, m=0,lmax
   SHS%trig(:,m,1)=dcos(m*SHS%lon)
   if(m>0)SHS%trig(:,m,2)=dsin(m*SHS%lon)
end do

SHS%save_pnm=save_pnm
if(save_pnm)then !also precompute legendre functions for all latitudes
   allocate(SHS%Pnm(nsh,SHS%nlat))
   !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(SHS)
   !$OMP DO
   do,i=1,SHS%nlat
      call PlmBar(p=SHS%pnm(:,i),lmax=SHS%lmax,z=dsin(SHS%lat(i)))
   end do
   !$OMP END DO
   !$OMP END PARALLEL

end if


end subroutine SH_synth_init

subroutine SH_dosynth(SHS,clm,slm,GRD)
use shtlbx
use shtools
implicit none
type(SHsynth_struct),intent(in)::SHS
double precision,intent(in)::clm(:),slm(:)
double precision,pointer,intent(inout)::GRD(:,:)

double precision,allocatable,dimension(:)::ctmp,stmp,trigtmp
integer,allocatable::st(:)
integer::i,m,pos
if(.not. associated(GRD))allocate(GRD(SHS%nlat,SHS%nlon))

allocate(st(0:SHS%lmax))

do,m=0,SHS%lmax
   st(m)=SH_pos(m-1,m-1)+1
end do

pos=SH_pos(SHS%lmax,SHS%lmax)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(clm,slm,SHS,GRD,st,pos)
allocate(ctmp(pos),stmp(pos))
!$OMP DO
do,i=1,SHS%nlat
   if(SHS%save_pnm)then
      ctmp=SHS%pnm(:,i)*clm
      stmp=SHS%pnm(:,i)*slm
   else
      call PlmBar(p=ctmp,lmax=SHS%lmax,z=dsin(SHS%lat(i)))
      stmp=ctmp*slm
      ctmp=ctmp*clm
   end if
   do,m=0,SHS%lmax
      GRD(i,:)=SHS%trig(:,m,1)*sum(ctmp(st(m):))+SHS%trig(:,m,2)*sum(stmp(st(m):))
   end do
end do
!$OMP END DO
deallocate(ctmp,stmp)
!$OMP END PARALLEL

deallocate(st)

end subroutine SH_dosynth


end module SH_synthesis

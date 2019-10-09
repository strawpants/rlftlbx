!!assembly of subroutines to construct observation equation associated with GPS deformation data
!!Coded by Roelof Rietbroek, Fri Jan  4 15:43:25 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!partly adapted from

!!Updated by Roelof Rietbroek, Thu Dec 10 18:08:59 2009
!!changed units of the scale parameter to parts per billion
!changed units of the rotation parameter to arcseconds

!!subroutine GPS_obseq_loadsh applies farrels loading theory to station locations
!!returns a design matrix A which relates station residuals in the local spherical frame to mass loading coefficients
!!input:
!!lon,lat: station coordinates in radians (vectors dimension: ndat (amount of data points)
!!lmin,lmax: required minimium and maximum degree of design matrix
!!typ: data type (must contain STAH, STAE or STAN to differentiate between height east and north displacement)

!!Output:
!! A: design matrix dimension ndat x nsh(amount of unknowns of the SH expansion) 
!! tag:unknown parameter names (character array (24 characters long) nsol entries)    
! subroutine GPS_obseq_loadsh(A,tag,typ,lon,lat,lmax,lmin,ltyp)
!!Updated by Roelof Rietbroek, Fri Feb 19 12:38:25 2010
!! made the SH observtion equation more general: also potential changes may be linked to deformations

!!Updated by Roelof Rietbroek, Mon Mar  1 08:58:00 2010
!!fixed bug (scale by Earth Radius) in units of rotation vector ( Helmert parameters and euler plate)

!!Updated by Roelof Rietbroek, Tue Mar 30 12:04:22 2010
!!made sign of the scale in the Helmert parameters positive ( positive height will cause a positive scale)

!!Updated by Roelof Rietbroek, Thu Apr 14 15:43:24 2011
!! added a check such that lmin >=1

!!Updated by Roelof Rietbroek, Tue Aug 21 14:06:01 2012
!! also added the Purcell et al 2011 relationship of uplift to GIA induced potential
!! also add a spherical harmonic expansion of poloidal horizontal deformation



subroutine GPS_obseq_sh(A,tag,typ,lon,lat,lmax,lmin,ltyp,frame,shtyp)
use SHtlbx
use GPStlbx,only:compstr
use shtools
implicit none
double precision, intent(out)::A(:,:)
character(24),intent(out)::tag(:)
character(*),intent(in),dimension(:)::typ
double precision,intent(in)::lon(:),lat(:)
integer,intent(in)::lmax,lmin
integer,intent(in)::ltyp ! load love numbers to be used
integer,intent(in)::shtyp,frame

!private variables
integer::pos,nc,ns,inds,indc,l,m,i,ndat,j
integer,allocatable,dimension(:)::Asvec,Acvec,pcvec,psvec
double precision::latold,lonold
double precision,allocatable,dimension(:)::p,dp,cosmlon,sinmlon,deg,ord,hnm,lnm,knm
character(1)::shtag
integer::stderr

stderr=0
!initializations

select case(shtyp)
case(1) ! surface loading coefficients expressed in equiv water height
   shtag='T'
case(2,5) ! Stokes coefficients ( normalized wrt GM/R)
   shtag='G'
case(3) !Purcell et al 2011 GIA relationship
   shtag='P' ! from 'P' ost Glacial Rebound
case(4)
   shtag='V' ! Poloidal horizontal deformation
case default
   write(stderr,*)"ERROR in GPS_obseq_sh: unknown shtype requested",shtyp
   stop
end select

if(lmin<1)then
   write(stderr,*)"ERROR in GPS_obseq_sh: lmin may be 1 minimum"
   stop
end if



ndat=size(lon,1) !number of data points

latold=-999.d0 !initial values guaranteed to work on first entry
lonold=-999.d0


!construct a spherical harmonic index vector for vectors p and dp (ass. legendre polynomials)
pos=SH_pos(lmax,lmax)

allocate(p(pos),dp(pos))
p=0.d0
dp=0.d0

!load love numbers
allocate(hnm(pos),lnm(pos))
!precompute vectors
allocate(cosmlon(pos),sinmlon(pos))
!degree and order vectors
allocate(deg(pos),ord(pos))

!allocate index vectors
nc=SH_pos(lmax,lmax)-SH_pos(lmin-1,lmin-1) !amount of cosine coefficients
ns=nc-lmax+lmin-1
!cosine vector
allocate(Acvec(nc),pcvec(nc))
!sine vector (no zero order terms)
allocate(Asvec(ns),psvec(ns))!amount of sine coefficients



inds=0
indc=0
do,i=1,pos
   call SH_lm(pos=i,l=l,m=m)
   !make degree and order vectors
   deg(i)=l
   ord(i)=m

!below only  l>=lmin and m .ne. 0 (for Slm coeff)
   if(l<lmin)cycle
   
   indc=indc+1
   Acvec(indc)=SH_tpos(l,m,0,lmax,lmin)
   pcvec(indc)=i
  
   write(tag(Acvec(indc)),'(a3,1x,I3,i3,10x,a4)')shtag//'CN',l,m,'GPSN'
   
   if(m .ne. 0)then

      inds=inds+1
      Asvec(inds)=SH_tpos(l,m,1,lmax,lmin)
      psvec(inds)=i
  
      write(tag(Asvec(inds)),'(a3,1x,I3,i3,10x,a4)')shtag//'SN',l,m,'GPSN'
   end if
end do


!write(*,*)nc,indc
!write(*,*)nc-lmax+lmin-1,inds




select case(shtyp)
case(1)
   !get load love numbers currently Farrels (isostropic)
   !set reference frame to: 'center of figure'
   call SH_loadlove(hnm=hnm,lnm=lnm,typ=ltyp,frame=frame)

   !rescale load love numbers to make them a conversion between station displacement and surface mass loading
   hnm=3*hnm*rho_w/((2*deg+1)*rho_e)
   lnm=3*lnm*rho_w/((2*deg+1)*rho_e)
case(2) ! relate to stokes coefficients
   !get load love numbers currently Farrels (isostropic)
   !set reference frame to: 'center of figure'
   call SH_loadlove(hnm=hnm,lnm=lnm,typ=ltyp,frame=frame)

   hnm=hnm*RE
   lnm=lnm*RE

case(5) !Same as shtyp==2 but also account for solid Earth loading effects
    !(i.e.potential also contains solid Earth contribution)
   allocate(knm(pos))
   call SH_loadlove(hnm=hnm,lnm=lnm,knm=knm,typ=ltyp,frame=frame)

   hnm=hnm*RE/(1+knm)
   lnm=lnm*RE/(1+knm)
case(3) ! get Purcell et al 2011 GIA--Uplift relationship
   call SH_loadGIA_uplift_ratio(rat=hnm,typ=1)
   ! do i=1,pos
   !    write(stderr,*)hnm(i)
   ! end do
   hnm=hnm*RE
   i=SH_pos(1,0)
   hnm(i:i+1)=1.d0 !explicitly set degree 1 components to 1
   !modify tags for degree 1 coefficients (should be plain uplift not potential)
   i=SH_tpos(1,0,0,lmax,lmin)
   tag(i)(1:1)='U'
   i=SH_tpos(1,1,0,lmax,lmin)
   tag(i)(1:1)='U'   
   i=SH_tpos(1,1,1,lmax,lmin)
   tag(i)(1:1)='U'   

case(4) ! horizontal deformation
   lnm=1.d0 !plain factor of 1, no Love number relationship
end select

do,i=1,ndat
   
   !calculate legendre polynomials and cos and sin terms if not different from previous step
   if(latold .ne. lat(i))then 
      latold=lat(i)

      !call SHTOOLS routine to calculate associated legendre polynomials and its derivative
     ! write(*,*)pos,size(p),size(dp),lmax
!      call PlmBar_d1(p=p,dp=dp,lmax=lmax,z=cos(pi/2-latold))

      call PlmBar_d1(p=p,dp1=dp,lmax=lmax,z=sin(latold))
!note that dp is the derivative wrt z and not wrt colatitude

      lonold=lon(i)
      sinmlon=sin(ord*lonold)
      cosmlon=cos(ord*lonold)
      

   end if


   !check which type it is and use appropriate observation equation
  
   if(compstr(typ(i),'*STAH*'))then !height observation equation
      if(shtyp .eq. 4) cycle !no horizontal component seen here
      !cos terms
      A(i,Acvec)=hnm(pcvec)*cosmlon(pcvec)*p(pcvec)
      !sine terms
      A(i,Asvec)=hnm(psvec)*sinmlon(psvec)*p(psvec)

   else if(compstr(typ(i),'*STAE*'))then !east observation equation
      if(shtyp .eq. 3) cycle ! ignore GIA relationship for East component
      !cos terms
      A(i,Acvec)=-ord(pcvec)*lnm(pcvec)*sinmlon(pcvec)*p(pcvec)/(cos(latold))
      !sine terms
      A(i,Asvec)=ord(psvec)*lnm(psvec)*cosmlon(psvec)*p(psvec)/(cos(latold))


   else if(compstr(typ(i),'*STAN*'))then !north observation equation
      if(shtyp .eq. 3) cycle ! ignore GIA relationship for the North component
            !cos terms
      A(i,Acvec)=lnm(pcvec)*cosmlon(pcvec)*dp(pcvec)*cos(latold)
      !sine terms
      A(i,Asvec)=lnm(psvec)*sinmlon(psvec)*dp(psvec)*cos(latold)
!the -cos(latold) originates from the chain rule dsin(colat)/dcolat=dz/dcolat

   end if
   

end do

end subroutine GPS_obseq_sh
! end subroutine GPS_obseq_loadsh


!!subroutine which returns a design matrix relating station deformations to helmert parameters (rigid translation and rotation of an ellipsoid)

Subroutine GPS_obseq_helmert(A,tag,typ,lon,lat)
use GPStlbx, only :compstr
use SHtlbx,only:RE,pi
implicit none
double precision,intent(out)::A(:,:)
character(24),intent(out)::tag(7)

character(*),intent(in),dimension(:)::typ
double precision,intent(in),dimension(:)::lon,lat


integer::i,ndat
double precision,dimension(3)::eh,ee,en
double precision::arcs2rad_RE

A=0.d0
ndat=size(typ,1)
arcs2rad_RE=RE*pi/(180.d0*3600)

!create unknown parameter tags
tag(1)='HTX                 [m]'
tag(2)='HTY                 [m]'
tag(3)='HTZ                 [m]'
tag(4)='HSC               [ppb]'
tag(5)='HRX              [arcs]'
tag(6)='HRY              [arcs]'
tag(7)='HRZ              [arcs]'


do,i=1,ndat
   if(compstr(typ(i),'*STAH*'))then !height displacement
      !calculate unit vector pointing in height direction
      eh(1)=cos(lat(i))*cos(lon(i))
      eh(2)=cos(lat(i))*sin(lon(i))
      eh(3)=sin(lat(i))
    
      !calculate Helmert transformation entries
      A(i,1)=eh(1)
      A(i,2)=eh(2)
      A(i,3)=eh(3)
      A(i,4)=RE*1d-9 ! scale will be dimensionless and in ppb
   !    A(h,7)=-RE

   else if(compstr(typ(i),'*STAE*'))then !East displacement

      !calculate unit vectors pointing in east and north
      ee(1)=-sin(lon(i))
      ee(2)=cos(lon(i))
      ee(3)=0.d0
      
      en(1)=-sin(lat(i))*cos(lon(i))
      en(2)=-sin(lat(i))*sin(lon(i))
      en(3)=cos(lat(i))
      
      !calculate entries for design matrix

      A(i,1)=ee(1)
      A(i,2)=ee(2)
      A(i,3)=ee(3)
      
!       A(i,5)=en(1)
!       A(i,6)=en(2)
!       A(i,7)=en(3)
      
      A(i,5)=en(1)*arcs2rad_RE
      A(i,6)=en(2)*arcs2rad_RE
      A(i,7)=en(3)*arcs2rad_RE

   else if(compstr(typ(i),'*STAN*'))then

      !calculate unit vectors pointing in east and north
      ee(1)=-sin(lon(i))
      ee(2)=cos(lon(i))
      ee(3)=0.d0
      
      en(1)=-sin(lat(i))*cos(lon(i))
      en(2)=-sin(lat(i))*sin(lon(i))
      en(3)=cos(lat(i))
      
      A(i,1)=en(1)
      A(i,2)=en(2)
      A(i,3)=en(3)
      
!       A(i,5)=-ee(1)
!       A(i,6)=-ee(2)
!       A(i,7)=-ee(3)

      A(i,5)=-ee(1)*arcs2rad_RE
      A(i,6)=-ee(2)*arcs2rad_RE
      A(i,7)=-ee(3)*arcs2rad_RE
      
   end if
end do

end Subroutine GPS_obseq_helmert


!subroutine to set up design matrix for GPS stations for the geocenter loading part (alternate way to define degree 1 deformation)
subroutine GPS_obseq_geocent(A,tag,typ,lon,lat,ltyp,frame)
use GPStlbx,only: compstr
use SHtlbx, only : SH_loadlove
implicit none
double precision,intent(out)::A(:,:)
character(24),intent(out)::tag(3)

character(*),intent(in),dimension(:)::typ
double precision,intent(in),dimension(:)::lon,lat
integer,intent(in)::ltyp,frame

double precision,dimension(3)::hnm,lnm,eh,ee,en
integer::i,ndat


ndat=size(typ,1)

!make unknowns tag

tag(1)='SXB'
tag(2)='SYB'
tag(3)='SZB'


!get load degree 1 load numbers (up to degree 1)
call SH_loadlove(hnm=hnm,lnm=lnm,typ=ltyp,frame=frame)


do,i=1,ndat

   if(compstr(typ(i),'*STAH*'))then

!unit vector in up direction
      eh(1)=cos(lat(i))*cos(lon(i))
      eh(2)=cos(lat(i))*sin(lon(i))
      eh(3)=sin(lat(i))
      
      A(i,1)=hnm(2)*eh(1)
      A(i,2)=hnm(2)*eh(2)
      A(i,3)=hnm(2)*eh(3)
      
   else if(compstr(typ(i),'*STAE*'))then
!unit vector in east direction
      ee(1)=-sin(lon(i))
      ee(2)=cos(lon(i))
      ee(3)=0.d0
      

      A(i,1)=lnm(2)*ee(1)
      A(i,2)=lnm(2)*ee(2)
      A(i,3)=lnm(2)*ee(3)
      
   else if(compstr(typ(i),'*STAN*'))then
      en(1)=-sin(lat(i))*cos(lon(i))
      en(2)=-sin(lat(i))*sin(lon(i))
      en(3)=cos(lat(i))


      A(i,1)=lnm(2)*en(1)
      A(i,2)=lnm(2)*en(2)
      A(i,3)=lnm(2)*en(3)
   
   end if

end do




end subroutine GPS_obseq_geocent

!subroutine which constructs a design matrix with a mean for each station
!quite simple separated fgor the sake of clarity
subroutine GPS_obseq_mean(A,tag,typ,domes)
  use GPStlbx, only: compstr
  implicit none
  double precision,intent(out)::A(:,:)
  character(24),intent(out)::tag(:)
  character(*),intent(in)::typ(:)
  character(17),intent(in)::domes(:)

  integer::i

!unknowns tag


do,i=1,size(typ,1)
   if(compstr(typ(i),'*STAH*'))then

      tag(i)='SH  '//domes(i)(1:17)
   else if(compstr(typ(i),'*STAE*'))then
      tag(i)='SE  '//domes(i)(1:17)
   else if(compstr(typ(i),'*STAN*'))then
      tag(i)='SN  '//domes(i)(1:17)
   end if
!design matrix
   A(i,i)=1.d0
end do




end subroutine GPS_obseq_mean


!design matrix for a trend
subroutine GPS_obseq_trend(A,tag,time,typ,domes)
  use GPStlbx, only: compstr
  implicit none
  double precision,intent(out)::A(:,:)
  character(24),intent(out)::tag(:)
  double precision,intent(in)::time !time in years
  character(*),intent(in)::typ(:)
  character(17),intent(in)::domes(:)
  integer::i
!unknowns tag

do,i=1,size(typ,1)
   if(compstr(typ(i),'*STAH*'))then

      tag(i)='SHP '//domes(i)(1:17)
   else if(compstr(typ(i),'*STAE*'))then
      tag(i)='SEP '//domes(i)(1:17)
   else if(compstr(typ(i),'*STAN*'))then
      tag(i)='SNP '//domes(i)(1:17)
   end if
!design matrix
   A(i,i)=time-2000.d0
end do




end subroutine GPS_obseq_trend

!design matrix for harmonic cycle
subroutine GPS_obseq_cycle(A,tag,time,ohm,harmtag,typ,domes)
  use GPStlbx, only: compstr
  implicit none
  double precision,intent(out)::A(:,:)
  character(24),intent(out)::tag(:)
  double precision,intent(in)::time,ohm !time in years
  character(*),intent(in)::typ(:)
  character(2)::harmtag
  character(17),intent(in)::domes(:)

integer::i,ndat,ind
! unknowns tag
ndat=size(typ,1)
ind=0
do,i=1,ndat

   if(compstr(typ(i),'*STAH*'))then
      ind=ind+2
      tag(ind-1)='SHC'//harmtag//domes(i)(2:17)
      tag(ind)='SHS'//harmtag//domes(i)(2:17)
   else if(compstr(typ(i),'*STAE*'))then
      ind=ind+2
      tag(ind-1)='SEC'//harmtag//domes(i)(2:17)
      tag(ind)='SNS'//harmtag//domes(i)(2:17)
   else if(compstr(typ(i),'*STAN*'))then
      ind=ind+2
      tag(ind-1)='SNC'//harmtag//domes(i)(2:17)
      tag(ind)='SNS'//harmtag//domes(i)(2:17)
   end if

   A(i,ind-1)=cos(ohm*time)
   A(i,ind)=sin(ohm*time)
end do


  
end subroutine GPS_obseq_cycle

!subroutine to set up Euler vectors for each plate
!An input list is required which has the plate name
!!Updated by Roelof Rietbroek, Wed Aug 22 16:48:13 2012
!! fixed bug which assigned wrong plates to position

subroutine GPS_euler_plates(A,tag,typ,plate,lon,lat)
  use GPStlbx, only: compstr
  use SHtlbx,only:pi,RE
  implicit none
  double precision,intent(out)::A(:,:)
  character(24),intent(out)::tag(:)
  character(*),intent(in)::typ(:),plate(:)
  double precision,intent(in)::lon(:),lat(:)
  character(24)::plate_accum(200) !reserve space for 200 individual plates
  double precision::arcs2rad_RE
  integer::nplate,ndat,i,j,shft
  logical::match

  plate_accum=''
  A=0.d0
  ndat=size(typ,1)
  arcs2rad_RE=RE*pi/(180.d0*3600)

  nplate=0
  do,i=1,ndat
     match=.false.
     do,j=1,nplate ! search for plate in already processed set
        if(plate(i) .eq. plate_accum(j))then
           match=.true.
           exit
        end if
     end do
     if(.not. match)then ! if not found a new entry will be created
        nplate=nplate+1
        plate_accum(nplate)=trim(plate(i))
        !create new unknown tags ( three for each plate)
        tag((nplate-1)*3+1)='PRX_'//trim(plate(i))
        tag((nplate-1)*3+1)(18:24)='[arcs]'
        tag((nplate-1)*3+2)='PRY_'//trim(plate(i))
        tag((nplate-1)*3+2)(18:24)='[arcs]'
        tag((nplate-1)*3+3)='PRZ_'//trim(plate(i))
        tag((nplate-1)*3+3)(18:24)='[arcs]'
     end if

     !create design matrix
     shft=(j-1)*3
!     write(0,*)shft,nplate,plate_accum(j),plate(i),typ(i)
     if(compstr(typ(i),'*STAH*'))then !height displacement
        A(i,shft+1:shft+3)=0.d0 ! Euler plate rotation has no effect on the up direction
     else if(compstr(typ(i),'*STAE*'))then
        A(i,shft+1)=-arcs2rad_RE*sin(lat(i))*cos(lon(i))
        A(i,shft+2)=-arcs2rad_RE*sin(lat(i))*sin(lon(i))
        A(i,shft+3)=arcs2rad_RE*cos(lat(i))

     else if(compstr(typ(i),'*STAN*'))then
        A(i,shft+1)=arcs2rad_RE*sin(lon(i))
        A(i,shft+2)=-arcs2rad_RE*cos(lon(i))
        A(i,shft+3)=0.d0
     end if
  end do


end subroutine GPS_euler_plates

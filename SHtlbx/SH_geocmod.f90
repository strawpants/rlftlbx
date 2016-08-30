!!subroutine to calculate geocenter potential degree 1 coefficients from am annual geocenter model
!!models are taken from cretaux 2002,which describes the following models:
!!mod=1: cretaux model derived from LAGEOS and TOPEX/POS laser and DORIS measurements
!!mod=2: Bouille model (2000) derived from LAGEOS and DORIS on T/P and SPOT
!!mod=3: Chen (2000) model derived from ocean atmosphere models and LAGEOS laser ranging 
!!Updated by Roelof Rietbroek, Tue Jun 12 10:03:08 2007
!!added new geocenter model:
!!mod=4: Mark-willem Jansen, Juergen Kusche, EJO Schrama, geocenter model
!!model derived from GPS in a combined GRACE GPS inversion (period 2003.0-2006.5)
!!The values hold for the cosine: x=xam*cos(ohm*(t-xph))(similar for y,z)
!!Where ohm is the annual frequency in radians per year (~265.25 days)
!!Coded by Roelof Rietbroek, Wed May 30 16:22:18 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue Nov 10 11:37:25 2009
!!added new (semi)annual geocenter model
!!Updated by Roelof Rietbroek, Wed Apr 27 15:00:57 2011
!!added new annual geocenter model (from special issue paper)
!!Updated by Roelof Rietbroek, Fri May 27 16:44:59 2011
!!fixed bug

!!Updated by Roelof Rietbroek, Sun Jan 20 16:36:58 2013
!! added secondary effect due to the Earth's flattening


subroutine SH_geocmod(c10,c11,s11,c30,c31,s31,time,time1,mod,geocfile)
use shtlbx,only:GRS80_C20,RE,pi
implicit none
integer,intent(in)::mod!parameter specifying which model to use
character(*),intent(in),optional::geocfile
double precision,intent(out)::c10,c11,s11
double precision,intent(in)::time !time in years (eg 2000.5 means somewhere half way of the year 2000) (center period if time1 not provided)
double precision,optional,intent(out)::c30,c31,s31 !affected by the flattening of the Earth
double precision,optional,intent(in)::time1 !end of the observation period also in years 
!!(if this parameter is present param. time is NOT used as the center time of observation but as the begin of the observation period)
double precision::xam,xph,yam,yph,zam,zph !annual phases(in days) and amplitudes(mm) of the geocentermotion
double precision::xsam,ysam,zsam,xsph,ysph,zsph
double precision::ohm,scale,d2y
double precision::x,y,z,xnext,ynext,znext,dumtime,dumtimenext,xt,yt,zt
double precision:: scalen3
integer::last,i,hit,unit
logical::semi
semi=.false.

!setup some constants
! pi=acos(-1.d0)
! RE = 0.6378136460d+07
ohm=2*pi
!scale factor to convert to normalized geopotential coefficients note the conversion for m to mm)
scale=1.d0/(RE*sqrt(3.d0)*1.d3)
scalen3=GRS80_C20*3.d0*sqrt(5.d0/7.d0)/(RE*1d3)

!conversion from days to years
d2y=1.d0/(365.25d0)
!initialize to zero
c10=0.d0
c11=0.d0
s11=0.d0

if(present(c30))c30=0.d0
if(present(c31))c31=0.d0
if(present(s31))s31=0.d0

xt=0.d0
yt=0.d0
zt=0.d0

select case(mod)
case(1)!creteaux model
zam=3.d0
xam=1.1d0
yam=3.7d0
zph=57.d0
xph=16.d0
yph=292.d0

case(2)!Bouille model (LAGEOS only)
zam=3.5d0
xam=2.1d0
yam=2.d0
zph=43.d0
xph=48.d0
yph=327.d0

case(3) !chen model using LAGEOS data and ocean/atmosphere models
zam=2.8d0
xam=2.2d0
yam=3.2d0
zph=46.d0
xph=59.d0
yph=299.d0

case(4) !GPS derive geocenter model (Mark-willem jansen)
zam=3.6d0
xam=2.3d0
yam=2.7d0
zph=18.d0
xph=26.d0
yph=331.d0
case(5)!annual + semi annual motion (2003-2007) from rietbroek 2009
   semi=.true.
   xam=2.1
   xph=75
   yam=4.5
   yph=338
   zam=2.5
   zph=63
   
   xsam=0.3
   xsph=143
   ysam=0.4
   ysph=94
   zsam=1.4
   zsph=107
case(6)!annual + semi annual motion (2003-2008) from Rietbroek et al 2011
   semi=.true.
   xam=2.1
   xph=56
   yam=3.4
   yph=327
   zam=3.0
   zph=18
   
   xsam=0.6
   xsph=162
   ysam=0.2
   ysph=121
   zsam=0.6
   zsph=126
case(7) !custom from file
   unit=17
   open(unit=unit,file=geocfile, form='formatted')
   last=0

   if(present(time1))then !take the average
      hit=0
      do while (last .eq. 0)
         read(unit,*,iostat=last)dumtime,x,y,z
         !      write(0,*)dumtime,x,y,z,last
         if(last .ne. 0)exit
!         write(0,*)dumtime,time1,time,x,y,z,last
         if(dumtime >= time .and. dumtime <=time1)then
            hit=hit+1
            
            ! c11=c11+scale*x*1000.d0
            ! c10=c10+scale*z*1000.d0
            ! s11=s11+scale*y*1000.d0
            
            xt=xt+x*1000.d0
            yt=yt+y*1000.d0
            zt=zt+z*1000.d0
            
            ! write(0,*)hit,x,y,z,c11/(scale*1000),s11/(scale*1000),c10/(scale*1000)
         else if(dumtime > time1)then
            exit
         end if
         
      end do
      
      if(hit .ne. 0)then
         ! c11=c11/dble(hit)
         ! c10=c10/dble(hit)
         ! s11=s11/dble(hit)

         xt=xt/dble(hit)
         yt=yt/dble(hit)
         zt=zt/dble(hit)

      end if
      
   else !linear interpolation
      !read first line
      read(unit,*)dumtime,x,y,z
      if(dumtime > time)return !quick return when nothing was found (values stay zero)
      last=0
      do while (last .eq. 0)
         read(unit,*,iostat=last)dumtimenext,xnext,ynext,znext
         if(last .ne. 0)exit
         if(dumtimenext >= time)then
            xt=(x+(time-dumtime)*(xnext-x)/(dumtimenext-dumtime))*1000.d0
            zt=(z+(time-dumtime)*(znext-z)/(dumtimenext-dumtime))*1000.d0
            yt=(y+(time-dumtime)*(ynext-y)/(dumtimenext-dumtime))*1000.d0
            exit
         end if
         !reset variables
         x=xnext
         y=ynext
         z=znext
         dumtime=dumtimenext
      end do
   end if
   
   close(unit)
   
case default
   return
end select


if(mod <=6)then
   !time1=2003.1
!time=2003.0
!write(*,*)ohm,pi
!now construct the geopontential coefficients
   if (present(time1)) then !use integration over the observation period
      xt=(dsin(ohm*(time1-xph*d2y))-dsin(ohm*(time-xph*d2y)))*xam/((time1-time)*ohm)
      zt=(dsin(ohm*(time1-zph*d2y))-dsin(ohm*(time-zph*d2y)))*zam/((time1-time)*ohm)
      yt=(dsin(ohm*(time1-yph*d2y))-dsin(ohm*(time-yph*d2y)))*yam/((time1-time)*ohm)
      !write(*,*)'c11 c10,s11',c11,c10,s11,time1-xph*d2y,time-xph*d2y

      if(semi)then
         xt=xt+(dsin(2*ohm*(time1-xsph*d2y))-dsin(2*ohm*(time-xsph*d2y)))*xsam/((time1-time)*2*ohm)
         zt=zt+(dsin(2*ohm*(time1-zsph*d2y))-dsin(2*ohm*(time-zsph*d2y)))*zsam/((time1-time)*2*ohm)
         yt=yt+(dsin(2*ohm*(time1-ysph*d2y))-dsin(2*ohm*(time-ysph*d2y)))*ysam/((time1-time)*2*ohm)
      end if

   else !just use the value of the cosine at the center period
      xt=xam*dcos(ohm*(time-xph*d2y))
      zt=zam*dcos(ohm*(time-zph*d2y))
      yt=yam*dcos(ohm*(time-yph*d2y))

      if(semi)then
         xt=xt+xsam*dcos(2*ohm*(time-xsph*d2y))
         zt=zt+zsam*dcos(2*ohm*(time-zsph*d2y))
         yt=yt+ysam*dcos(2*ohm*(time-ysph*d2y))
      end if


   end if
end if
!write(*,*)'c11 c10,s11',c11,c10,s11

!convert geocenter to coefficients

c11=xt*scale
s11=yt*scale
c10=zt*scale



if(present(c30))c30=zt*scalen3
if(present(c31))c31=xt*scalen3*2/sqrt(6.d0)
if(present(s31))s31=yt*scalen3*2/sqrt(6.d0)

!write(0,*)c30,c31,s31


end subroutine SH_geocmod



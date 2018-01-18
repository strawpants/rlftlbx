!module GPStlbx contains interfaces and constants used in SINEX file and GPS proccesing
!!Coded by Roelof Rietbroek, Wed Aug 15 11:56:12 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de



module GPStlbx
implicit none

interface
   subroutine get_sinex(file,nx,ny,x,y,covxy,covxypack,covdiagx,inc1,excl1,inc2,excl2,xlist,gpswk)
     implicit none
     character(*),intent(in)::file !obligatory argument
     !optional arguments
     character(*),intent(in),optional::inc1,excl1,inc2,excl2
     integer,intent(out),optional::nx,ny,gpswk
     double precision,intent(out),optional,dimension(:)::x,y,covdiagx,covxypack
     double precision,intent(out),optional,dimension(:,:)::covxy
     character(*),intent(out),optional,dimension(:)::xlist
   end subroutine get_sinex

   function compstr(str,srch)
     implicit none
     character(*),intent(in)::str,srch
     logical::compstr
   end function compstr

   subroutine GPS_rotmat(type,station,pos,rot)
     implicit none
     character(*),dimension(:),intent(inout)::type
     character(*),dimension(:),intent(in)::station
     double precision,dimension(:),intent(in)::pos
     double precision,dimension(:,:), intent(out)::rot
   end subroutine GPS_rotmat

   subroutine get_sinexblock(file,blockname,nlines,list,inc1,excl1)
     implicit none
     !obligatory arguments
     character(*),intent(in)::file,blockname
     !optional arguments
     integer,intent(out),optional::nlines !amount of lines (without comments)in the block
     character(*),intent(out),optional,dimension(:)::list
     character(*),intent(in),optional::inc1,excl1
   end subroutine get_sinexblock

   function GPS_week(dtstr,rem)
     implicit none
     double precision,optional,intent(out)::rem
     character(*),intent(in)::dtstr
     integer::GPS_week
   end function GPS_week

   function GPS_year(gpswk,dwk)
     implicit none
     integer,intent(in)::gpswk
     double precision,intent(in),optional::dwk !optional GPSweek as double (also fractional week)
     double precision::GPS_year
   end function GPS_year
   
   function sinex_date(date)
     implicit none
     double precision,intent(in)::date
     character(12)::sinex_date
   end function sinex_date

   function GPS_date(dtstr)
     implicit none
     character(*),optional::dtstr
     character(12)::GPS_date  
   end function GPS_date

   function norm_date(sindat)
     implicit none
     character(*),intent(in)::sindat
     character(12)::norm_date
   end function norm_date

   function date_2_jd(dy,mn,yr,hr,min,sec)
     implicit none
     double precision::date_2_jd
     integer::yr,dy,mn,hr,min,sec
   end function date_2_jd

   subroutine jd_2_date(jd,dy,mn,yr,sec)
     implicit none
     double precision,intent(in)::jd
     integer,intent(out)::mn,yr,dy,sec
   end subroutine jd_2_date


    subroutine getYYMMDDSS(decyr,yr,mn,dd,sec)
    implicit none
    integer,intent(out)::yr,mn,dd,sec
    double precision,intent(in)::decyr
    end subroutine getYYMMDDSS  

!sinex write routines
   subroutine init_sinex(file,unit,line)
     implicit none
     integer,intent(out)::unit
     character(*),intent(in),optional::file
     character(52),intent(in)::line
   end subroutine init_sinex
   
   subroutine write_sinexblock(unit,sblock,block,comm)
     implicit none
     integer,intent(in)::unit
     character(*),intent(in)::sblock
     character(*),intent(in),optional::comm
     character(80),dimension(:),intent(in)::block
   end subroutine write_sinexblock

   subroutine write_sinex_MATlow(unit,mat,matpack,type)
     implicit none
     integer,intent(in)::unit
     double precision,intent(in),optional,dimension(:)::matpack
     double precision,intent(in),optional,dimension(:,:)::mat
     character(4),intent(in)::type !either COVA, CORR or INFO
   end subroutine write_sinex_MATlow

   subroutine end_sinex(unit)
     implicit none
     integer,intent(in)::unit
   end subroutine end_sinex


   function get_solu_id(file,station,gpswk,minwks)
     implicit none
     character(*),intent(in)::file
     character(*),intent(in)::station
     integer,intent(in)::gpswk
     character(12)::get_solu_id 
     !optional variables
     integer,optional,intent(in)::minwks !minimum interval for a solution id
   end function get_solu_id

   subroutine get_cumdata(file,inc,pos,vel,reftime,error)
     implicit none
     character(*),intent(in)::file,inc
     integer,intent(out)::error
     double precision, intent(out)::pos,vel,reftime
   end subroutine get_cumdata

   subroutine conv2itrf2000(typ,station,pos,posnew,gpswk,refwk)
     implicit none
     double precision,intent(in),dimension(:)::pos
     double precision,intent(out),dimension(:)::posnew
     character(*),intent(in),dimension(:)::station,typ
     integer,intent(in),optional::gpswk
     integer,intent(in)::refwk

   end subroutine conv2itrf2000

   function get_helmert(gpswk,refwk)
     implicit none
     integer,intent(in)::gpswk,refwk
     double precision::get_helmert(7)
   end function get_helmert

   subroutine conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,dt)
     implicit none
     double precision,intent(in):: tx,ty,tz,rx,ry,rz,d,dt
     double precision,intent(in):: dtx,dty,dtz,drx,dry,drz,dd
     double precision,intent(out):: helm(7)
   end subroutine conv
 
!old style GPS observation equations  
   subroutine GPS_obseq(A,typ,lon,lat,helm,geocent,lmax,lmin)
     implicit none
     double precision,intent(out),dimension(:,:)::A
     double precision,intent(in),dimension(:)::lon,lat
     logical,intent(in)::helm,geocent
     integer,intent(in)::lmax,lmin
     character(*),intent(in),dimension(:)::typ

   end subroutine GPS_obseq
   
!newer style observation equations (separated)
!    subroutine GPS_obseq_loadsh(A,tag,typ,lon,lat,lmax,lmin,ltyp)
   subroutine GPS_obseq_sh(A,tag,typ,lon,lat,lmax,lmin,ltyp,frame,shtyp)
     implicit none
     double precision, intent(out)::A(:,:)
     character(24),intent(out)::tag(:)
     character(*),intent(in),dimension(:)::typ
     double precision,intent(in)::lon(:),lat(:)
     integer,intent(in)::lmax,lmin,ltyp
     integer,intent(in)::shtyp,frame
   end subroutine GPS_obseq_sh


   Subroutine GPS_obseq_helmert(A,tag,typ,lon,lat)
     implicit none
     double precision,intent(out)::A(:,:)
     character(24),intent(out)::tag(7)
     
     character(*),intent(in),dimension(:)::typ
     double precision,intent(in),dimension(:)::lon,lat
   end Subroutine GPS_obseq_helmert

   subroutine GPS_obseq_geocent(A,tag,typ,lon,lat,ltyp,frame)
     implicit none
     double precision,intent(out)::A(:,:)
     character(24),intent(out)::tag(3)
     integer,intent(in)::ltyp,frame

     character(*),intent(in),dimension(:)::typ
     double precision,intent(in),dimension(:)::lon,lat
   end subroutine GPS_obseq_geocent

   subroutine GPS_obseq_mean(A,tag,typ,domes)
     implicit none
     double precision,intent(out)::A(:,:)
     character(24),intent(out)::tag(:)
     character(*),intent(in)::typ(:)
     character(17),intent(in)::domes(:)
   end subroutine GPS_obseq_mean

   subroutine GPS_obseq_trend(A,tag,time,typ,domes)
     implicit none
     double precision,intent(out)::A(:,:)
     character(24),intent(out)::tag(:)
     double precision,intent(in)::time !time in years
     character(*),intent(in)::typ(:)
     character(17),intent(in)::domes(:)
   end subroutine GPS_obseq_trend

   subroutine GPS_obseq_cycle(A,tag,time,ohm,harmtag,typ,domes)
     implicit none
     double precision,intent(out)::A(:,:)
     character(24),intent(out)::tag(:)
     double precision,intent(in)::time,ohm !time in years
     character(*),intent(in)::typ(:)
     character(2)::harmtag
     character(17),intent(in)::domes(:)
   end subroutine GPS_obseq_cycle

   subroutine GPS_euler_plates(A,tag,typ,plate,lon,lat)
     implicit none
     double precision,intent(out)::A(:,:)
     character(24),intent(out)::tag(:)
     character(*),intent(in)::typ(:),plate(:)
     double precision,intent(in)::lon(:),lat(:)
   end subroutine GPS_euler_plates

   subroutine ECEF_2_geodetic(x,y,z,lon,lat,h)
     implicit none
     double precision,intent(in)::x,y,z
     double precision,intent(out)::lon,lat,h
   end subroutine ECEF_2_geodetic

   subroutine geodetic_2_ECEF(lon,lat,h,x,y,z)
     implicit none
     double precision,intent(in)::lon,lat,h
     double precision,intent(out)::x,y,z
   end subroutine geodetic_2_ECEF
end interface

end module GPStlbx

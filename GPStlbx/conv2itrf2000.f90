!fortran subroutine which converts XYZ coordinates given in a certain ITRF frame (dependent on time) to the ITRF2000 frame.
!transformation is made using a seven parameter Helmert transform:

!
!|Xa|   |cx|   |  d -rz  ry ||Xb|
!|Ya| = |cy| + | rz  d  -rx ||Yb|
!|Za|   |cz|   |-ry  rx  d  ||Zb|

subroutine conv2itrf2000(typ,station,pos,posnew,gpswk,refwk)
use GPStlbx,only:get_helmert,compstr
use FORTtlbx
implicit none
double precision,intent(in),dimension(:)::pos
double precision,intent(out),dimension(:)::posnew
character(*),intent(in),dimension(:)::station,typ
integer,intent(in)::refwk
integer,intent(in),optional::gpswk


integer::ndat,xind,yind,zind,arg,i,j,gpswkold,stderr
double precision,allocatable,dimension(:,:)::rot
double precision,allocatable,dimension(:)::trans
double precision::helm(7)
logical, allocatable,dimension(:)::done
logical::multi

!get size of vector
ndat=size(pos,1)
!write(*,*)'size:',ndat
!allocate matrices and initialize
allocate(trans(ndat),rot(ndat,ndat),done(ndat))
trans=0.d0
rot=0.d0
done=.false.
helm=0.d0
gpswkold=0

stderr=0
!get helmert parameters for epoch gpswk
!write(*,*)'test',helm(1)


   helm=get_helmert(gpswk,refwk)


!write(*,*)'test2',helm(1)

!now loop over entries in list vector (containing station and X Y Z entries)and create rotation and translation matrix

do,i=1,ndat
   if(done(i))cycle
   !(re)calculate helmert parameters if a vectro with refernce times is provided

   arg=0
   !search for X Y and Z entries for a particular station
   do,j=1,ndat

      if(.not. compstr(station(j),station(i)))cycle

      if(compstr(typ(j),'*X*'))then
         xind=j
         arg=arg+1
      else if(compstr(typ(j),'*Y*'))then
         yind=j
         arg=arg+1
      else if(compstr(typ(j),'*Z*'))then
         zind=j
         arg=arg+1
      end if
      
      if(arg .eq. 3)then !all components found
         !translation vector
         trans(xind)=helm(1)
         trans(yind)=helm(2)
         trans(zind)=helm(3)
         !rotation matrix
         rot(xind,xind)=helm(4)+1
         rot(yind,yind)=rot(xind,xind)
         rot(zind,zind)=rot(xind,xind)
         
         !off diagonal
         rot(xind,yind)=-helm(7)
         rot(yind,xind)=-rot(xind,yind)

         rot(xind,zind)=helm(6)
         rot(zind,xind)=-rot(xind,zind)

         rot(yind,zind)=-helm(5)
         rot(zind,yind)=-rot(yind,zind)
         

         done(xind)=.true.
         done(yind)=.true.
         done(zind)=.true.
        ! write(*,*)xind,yind,zind,trans(xind),trans(yind),trans(zind)
        ! write(*,*)rot(xind,xind),rot(yind,yind),rot(zind,zind)
       !  write(*,*)station(j),station(i),pos(i)
         exit
      end if
   end do
end do


!now apply transformation


posnew=trans+(matmulblas(rot,pos))

!write(*,*)posnew(1:3) 



end subroutine conv2itrf2000




!Function to retrieve helemert parameters to convert to itrf2000 frame
!adapted from Juergen kusches itrf.f
function get_helmert(gpswk,refwk)
  use GPStlbx,only: GPS_year
  implicit none
  integer,intent(in)::gpswk,refwk
  !the gpswk parameter determine which itrf system is asssumed for the data and the refwk determines the reference week (used in correcting for time drifts)
  double precision::tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,t0,t
  double precision,dimension(7)::helm,helm0,get_helmert
  integer::gpswk_min,gpswk_max,j
!  write(*,*)'entering get_helmert'
!write(*,*)'transforming to itrf2000 week: ',gpswk,GPS_year(gpswk)
 t=GPS_year(refwk)

  gpswk_min = 900
  !     gpswk_max = 1276
  gpswk_max = 1600
      
  !-------------------------------------------------------------------------------
  if ((gpswk.lt.gpswk_min).or.(gpswk.gt.gpswk_max)) then
!     write(*,*)gpswk_min,gpswk,gpswk_max
     stop 'out of range...stop'
  end if
  
  !-------------------------------------------------------------------------------
!  ITR96  (itrf97 -> itrf96) , Weighted transf. of the 47 ITRF(IGS) stations
!1999.58 -0.03 -.05  1.47 -1.43 -.159  .263  .060   cm,ppb,mas
!          .07 -.01   .19  -.19 -.013  .015 -.003   cm,ppb,mas/y

  if (gpswk.lt.1021) then
     write(*,*) 'ITRF96, 1021 > gpswk'
     !        first to ITRF97
     t0= 1999.58d0
     tx = 0.3d0
     ty = 0.5d0
     tz = -14.7d0
     rx = 0.159d0
     ry = -0.263d0
     rz = -0.060d0
     d = 1.430d0
     dtx = -0.7d0
     dty = 0.1d0
     dtz = -1.9d0
     drx = 0.013d0
     dry = -0.015d0
     drz = 0.003d0
     dd = 0.192d0
!ITR97  (itrf00 -> itrf97) , Weighted transf. of the 48 ITRF(IGS) stations
!1998.00  0.60  0.56 -2.01  1.403  0.040 -0.001 0.043   cm,ppb,mas
!        -0.04 -0.08 -0.15  0.012 -0.004  0.001 0.030   cm,ppb,mas/y

     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm0,t-t0)
     !        then to ITRF2000 (IGS00)
     t0= 1998.0d0
     tx = -6.0d0
     ty = -5.6d0
     tz =  20.1d0
     rx = -0.040d0
     ry =  0.001d0
     rz = -0.043d0
     d = -1.403d0
     dtx = 0.4d0
     dty = 0.8d0
     dtz = 1.5d0
     drx = 0.004d0
     dry = -0.001d0
     drz = -0.030d0
     dd = -0.012d0
     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,t-t0)
  
        helm=helm0+helm
  
     !  write(*,'(i8,14g20.10)') gpswk, (helm(j),j=1,7)
  end if
  
  !-------------------------------------------------------------------------------
  
  if ((gpswk.ge.1021).and.(gpswk.lt.1143)) then
     write(*,*) 'ITRF97, 1143 > gpswk >= 1021'
     t0= 1998.0d0
     tx = -6.0d0
     ty = -5.6d0
     tz =  20.1d0
     rx = -0.040d0
     ry =  0.001d0
     rz = -0.043d0
     d = -1.403d0
     dtx = 0.4d0
     dty = 0.8d0
     dtz = 1.5d0
     drx = 0.004d0
     dry = -0.001d0
     drz = -0.030d0
     dd = -0.012d0
     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,t-t0)
    !  write(*,'(i8,14g20.10)') gpswk, (helm(j),j=1,7)
  end if
  
  !-------------------------------------------------------------------------------
  
  if ((gpswk.ge.1143).and.(gpswk.lt.1253)) then
     write(*,*) 'ITRF2000 (IGS00), 1253 > gpswk >= 1143'
     t0= 2000d0
     tx = 0d0
     ty = 0d0
     tz = 0d0
     rx = 0d0
     ry = 0d0
     rz = 0d0
     d = 0d0
     dtx = 0d0
     dty = 0d0
     dtz = 0d0
     drx = 0d0
     dry = 0d0
     drz = 0d0
     dd = 0d0
     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,t-t0)
     !     write(*,'(i8,14g20.10)') gpswk, (helm(j),j=1,7)
  end if
  
  !-------------------------------------------------------------------------------
!    Transformation from IGb00 to IGS00 at epoch 01-JAN-1998:
!      ====================================================================
!               TX     TY     TZ    RX      RY      RZ       D
!               mm     mm     mm    mas     mas     mas     ppb
!      Offset  -0.1    0.2    0.2  0.001   0.004  -0.006   0.116
!      +/-      0.5    0.5    0.5   .017    .018    .020    .07
!      --------------------------------------------------------------------
!              dTX    dTY    dTZ    dRX     dRY     dRZ     dD
!              mm/y   mm/y   mm/y   mas/y   mas/y   mas/y   ppb/y
!      Drift    0.0   -0.2   -0.1  -.004  -0.005   0.000  -0.056
!      +/-      0.5    0.5    0.5   .017    .018    .020    .07
!      --------------------------------------------------------------------
  
  if ((gpswk.ge.1253) .and. (gpswk .lt. 1400)) then
     write(*,*) 'ITRF2000 (IGb00), gpswk >= 1253'
     t0= 1998.0d0
     tx = -0.1d0
     ty =  0.2d0
     tz =  0.2d0
     rx =  0.001d0
     ry =  0.004d0
     rz = -0.006d0
     d = 0.116d0
     dtx = 0.0d0
     dty = -0.2d0
     dtz = -0.1d0
     drx = -0.004d0
     dry = -0.005d0
     drz = 0.0d0
     dd = -0.056d0
     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,t-t0)
     !    write(*,'(i8,14g20.10)') gpswk, (helm(j),j=1,7)
  end if
  !mail 5447
  !4-1)	IGb00 - IGS05 Transformation:
  !-------------------------------------
  
  !	The (14) transformation parameters were estimated using a subset of
  !60 stations from the proposed IGS05 to the current IGS realization (IGb00).
  !The results are below, along with the official ITRF estimates (epoch 2000.0).
  !To ensure consistency between IGb00 and IGS05, the latter did not include
  !the relative to absolute phase center correction for the estimation of
  !the transformation parameters.
  !
  !( Epoch 2000.0 )                IGS          ITRF       IGS-ITRF
  !   R X (mas)   :              -0.0224      0.0000       -0.0224
  !   R Y (mas)   :               0.0341      0.0000        0.0341
  !   R Z (mas)   :              -0.0099      0.0000       -0.0099
  !   T X (m)     :               0.0000      0.0001       -0.0001
  !   T Y (m)     :              -0.0017     -0.0008       -0.0009
  !   T Z (m)     :              -0.0053     -0.0058        0.0005
  !   SCL (ppb)   :               0.8473      0.4000        0.4473
  ! d R X (mas/y) :               0.0033      0.0000        0.0033
  ! d R Y (mas/y) :              -0.0001      0.0000       -0.0001
  ! d R Z (mas/y) :              -0.0161      0.0000       -0.0161
  ! d T X (m/y)   :              -0.0004     -0.0002       -0.0002
  ! d T Y (m/y)   :               0.0007      0.0001        0.0006
  ! d T Z (m/y)   :              -0.0018     -0.0018        0.0000
  ! d SCL (ppb/y) :               0.1748      0.0800        0.0848
  
  
  !mail 5455 
  !Transformations:
  !================
  !
  !	For those of you who may want the transformation from
  !IGS05 (with absolute antenna phase center) and IGb00 (relative
  !antenna phase center) they are (at epoch 2000.0):
  !
  !  R X (mas)   :              -0.0070
  ! R Y (mas)   :               0.0340
  ! R Z (mas)   :              -0.0069
  ! T X (m)     :              -0.0003
  ! T Y (m)     :              -0.0015
  ! T Z (m)     :              -0.0061
  ! SCL (ppb)   :               0.7125
  !d R X (mas/y) :               0.0033
  !d R Y (mas/y) :              -0.0001
  !d R Z (mas/y) :              -0.0161
  !d T X (m/y)   :              -0.0004
  !d T Y (m/y)   :               0.0007
  !d T Z (m/y)   :              -0.0018
  !d SCL (ppb/y) :               0.1748
  
  
  if (gpswk.ge.1400) then
     write(*,*) 'ITRF2005 (IGS05), gpswk >= 1400'
     !first convert to IGS05 to IGb00 (from IGSmail 5447)
   !   t0= 2000.00d0
!      tx = 0.0d0
!      ty = 1.7d0
!      tz = 5.3d0
!      rx = 0.0224d0
!      ry = -0.0341d0
!      rz = 0.0099d0
!      d = -0.8473d0
!      dtx = 0.4d0
!      dty = - 0.7d0
!      dtz = 1.8d0
!      drx = -0.0033d0
!      dry = 0.0001d0
!      drz = -0.0161d0
!      dd = -0.1748d0

 t0= 2000.00d0
     tx = 0.0d0
     ty = -1.7d0
     tz = -5.3d0
     rx = -0.0224d0
     ry = 0.0341d0
     rz = -0.0099d0
     d = 0.8473d0
     dtx = -0.4d0
     dty =  0.7d0
     dtz = -1.8d0
     drx = 0.0033d0
     dry = -0.0001d0
     drz = -0.0161d0
     dd = 0.1748d0
     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm0,t-t0)
     !then convert Igb00 to IGS00
     t0= 1998.0d0
     tx = -0.1d0
     ty =  0.2d0
     tz =  0.2d0
     rx =  0.001d0
     ry =  0.004d0
     rz = -0.006d0
     d = 0.116d0
     dtx = 0.0d0
     dty = -0.2d0
     dtz = -0.1d0
     drx = -0.004d0
     dry = -0.005d0
     drz = 0.0d0
     dd = -0.056d0
     call conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,t-t0)
     !add helmert parameters
    
        helm=helm0+helm
    
     !  write(*,'(i8,14g20.10)') gpswk, (helm(j),j=1,7)
  end if
  
  
  !return values from helm
  get_helmert=helm      
  
  !-------------------------------------------------------------------------------
  
  
end function get_helmert
               

subroutine conv(tx,ty,tz,rx,ry,rz,d,dtx,dty,dtz,drx,dry,drz,dd,helm,dt)
  implicit none
  double precision,intent(in):: tx,ty,tz,rx,ry,rz,d,dt
  double precision,intent(in):: dtx,dty,dtz,drx,dry,drz,dd
  double precision,intent(out):: helm(7)
  double precision:: radius,pi,rho,mas2rad
  radius = 6380000d0
  pi = dacos(-1.0d0)
  rho = 180d0/pi
  mas2rad = (1d0/rho)/3600000.d0


  helm(1)=(tx+dt*dtx)/1000.d0
  helm(2)=(ty+dt*dty)/1000.d0
  helm(3)=(tz+dt*dtz)/1000.d0
  helm(4)=(d+dt*dd)*1e-9
  helm(5)=(rx+dt*drx)*mas2rad
  helm(6)=(ry+dt*dry)*mas2rad
  helm(7)=(rz+dt*drz)*mas2rad
!!factor of -1 below accounts for backtransformation
!   helm(1)=-(tx+dt*dtx)/1000.d0
!   helm(2)=-(ty+dt*dty)/1000.d0
!   helm(3)=-(tz+dt*dtz)/1000.d0
!   helm(4)=-(d+dt*dd)*1e-9
!   helm(5)=-(rx+dt*drx)*mas2rad
!   helm(6)=-(ry+dt*dry)*mas2rad
!   helm(7)=-(rz+dt*drz)*mas2rad
end subroutine conv





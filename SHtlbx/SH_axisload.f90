!program SH_delta calculates the delta function at a certain location exspressed in spherical harmonics
!!adpated from SH_delta.f made by Jkusche
!!Coded by Roelof Rietbroek, Tue Jan 13 10:37:51 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Tue Jul 28 15:14:21 2009
!!added parabolic load
!!added degree wise output (isotropic)

!!Updated by Roelof Rietbroek, Wed Apr 28 15:05:22 2010
!! allowed reading of multiple lon,lat points 

!!Updated by Roelof Rietbroek, Fri Jan 21 17:07:13 2011
! fixed a bug causing out of bound write on degscale

!!Updated by Roelof Rietbroek, Mon Dec  8 15:36:03 2014
!!output degree wise values in terms of unnormalized Legendre functions


program SH_axisload
  use shtools
  use shtlbx
  use forttlbx
  implicit none
  integer::lmax,lmin,itharg,narg,m,l,i,ind,stderr,pos
  double precision,allocatable,dimension(:)::clm,slm,plm,degnorm,degscale,Pl
  double precision::lon,lat,deg2rad,silonm,colonm,psi
  character(200)::dum
  integer::iargc,loadtype,stdout
  logical::iso
  double precision,pointer,dimension(:)::lonv=>null()
  double precision,pointer,dimension(:)::latv=>null()
  integer::np,last,unit,mem,chunk,p

  !defaults
  np=1 ! default is 1 point
  chunk=1000
  unit=13
  stderr=0
  deg2rad=pi/180.d0
  lmin=0
  lmax=200
  lon=0
  lat=0
    
  loadtype=0 !no load type specified
  iso=.false.
  stdout=6

  !!process command line options
  narg=iargc()
  if(narg<1)call help()
  itharg=0
  do,i=1,narg
     itharg=itharg+1
     if(itharg > narg)exit
     call getarg(itharg,dum)
     if(dum(1:1) .eq. '-')then !argument is an option
        select case(dum(2:3))
        case('l=')!limit maximum and minimum degree
           ind=index(dum,',')
           if(ind .ne. 0) then !also a min degree is specified
              read(dum(4:ind-1),*)lmax
              read(dum(ind+1:),*)lmin
           else
              read(dum(4:),*)lmax
           end if
        case('p=') ! get position
           ind=index(dum,',')
           read(dum(4:ind-1),*)lon
           read(dum(ind+1:),*)lat
           call realloc_ptr(lonv,1)
           call realloc_ptr(latv,1)
           lonv(1)=lon
           latv(1)=lat
        case('p:') ! several lonlat points contained in a file
           np=0
           last=0
           mem=0
           open(unit=unit,file=dum(4:),form='formatted')
           
           do while(last .eq. 0)
              read(unit=unit,iostat=last,fmt=*)lon,lat
              if(last .ne. 0) exit
              np=np+1
              if(mem < np)then ! we need to reallocate data ( per chunk)
                 call realloc_ptr(lonv,chunk)
                 call realloc_ptr(latv,chunk)
                 mem=mem+chunk
              end if
              lonv(np)=lon
              latv(np)=lat
           end do

           close(unit)
           
        case('u')!unit load
           loadtype=1 !although it is the default
        case('da')!diskload with radius in angle (degrees)
           loadtype=2
           read(dum(5:),*)psi
           psi=psi*deg2rad
        case('dr') ! diskload width radius in km
           loadtype=2
           read(dum(5:),*)psi
           psi=psi*1000.d0/RE !calculate angle in radians from radius in km
        case('ds')!diskload with surface in km ^2
           loadtype=2
           read(dum(5:),*)psi ! read surface area on the Earth 
           psi=acos(1-psi/(2*pi*(RE/1000.d0)**2)) !calculate angle in radians from radius in km

        case('ca')!diskload with radius in angle (degrees)
           loadtype=3
           read(dum(5:),*)psi
           psi=psi*deg2rad
        case('cr') ! diskload width radius in km
           loadtype=3
           read(dum(5:),*)psi
           psi=psi*1000.d0/RE !calculate angle in radians from radius in km
        case('cs')!diskload with surface in km ^2
           loadtype=3
           read(dum(5:),*)psi ! read surface area on the Earth 
           psi=acos(1-psi/(2*pi*(RE/1000.d0)**2)) !calculate angle in radians from radius in km
        case('i')!isotropic output
           iso=.true.
        case('h')
           call help()
        case default
           write(stderr,*)'ERROR: unknown option selected:',dum(1:3)
           stop
        end select
     else                      !argument is not an option
        
     end if
  end do
  
if(.not. iso .and. .not. associated(lonv))then
   write(stderr,*)"ERROR: -p option must be specified"
   stop
end if

!precompute degree dependent factor for the addition theorem
allocate(degnorm(0:lmax+1),degscale(0:lmax))

do,l=0,lmax+1
!   degnorm(l)=sqrt(2.d0*l+1)
      degnorm(l)=2.d0*l+1
end do


!create degree dependent scale
select case(loadtype)
case(0)
   write(stderr,*)"ERROR: no loading type selected"
   stop
case(1)!unit load
   degscale=degnorm(0:lmax)
case(2)!diskload
   !get legendre polynomials for disk edge
   allocate(Pl(lmax+2))
!    call Plbar(Pl,lmax+1,cos(psi))
   !degree 0
!    degscale(0)=(1-cos(psi)/degnorm(1))/2.d0
!    do,l=1,lmax
!       degscale(l)=(Pl(l)/(degnorm(l)*degnorm(l-1))-&
!            Pl(l+2)/(degnorm(l)*degnorm(l+1)))/2.d0
!    end do

   !unnormalized Legendre polynomials
   call PLegendre(Pl,lmax+1,cos(psi))
   
   degscale(0)=(1-cos(psi))/2.d0
   do,l=1,lmax
      degscale(l)=(Pl(l)-Pl(l+2))/2.d0
   end do
case(3) !parabolic load
   degscale(0)=(1-cos(psi))/3.d0
   do,l=1,lmax
      degscale(l)=((cos((l-1)*psi)-cos(l*psi))/(l-0.5)-&
           (cos((l+1)*psi)-cos((l+2)*psi))/(l+1.5))&
           /((1.d0-cos(psi))*4.d0)
   end do
case default
   write(stderr,*)"ERROR: Unknown loading type selected:",loadtype
   stop
end select

if(iso)then !write degree wise components only (for unnormalized Legendre functions)
   do,l=0,lmax
      write(stdout,*)l,degscale(l)
   end do

   !premature stop of the program no 2D output
   stop
end if
  
pos=SH_pos(lmax,lmax)
allocate(clm(pos),slm(pos),plm(pos))
clm=0.d0
slm=0.d0

do,p=1,np ! loop over all points

   lat=latv(p)
   lon=lonv(p)

   !calculate associated Legendre functions
   call Plmbar(plm,lmax,sin(lat*deg2rad))
   !zero order
   do,l=max(0,lmin),lmax
      pos=SH_pos(l,0)
      clm(pos)=clm(pos)+plm(pos)*degscale(l)/degnorm(l)
   end do
!orders > 0
   do,m=1,lmax
      silonm=sin(m*deg2rad*lon)
      colonm=cos(m*deg2rad*lon)
      !     if(m==90)write(0,'(2F21.18)')pi,atan2(1.d0,1.d0)*4-acos(0.d0)*2
      do,l=max(m,lmin),lmax
         pos=SH_pos(l,m)
         clm(pos)=clm(pos)+colonm*plm(pos)*degscale(l)/degnorm(l)
         slm(pos)=slm(pos)+silonm*plm(pos)*degscale(l)/degnorm(l)
      end do
   end do
end do ! end loop over points

!write to standard output
call SH_write(clm=clm,slm=slm)

end program SH_AXISLOAD

subroutine help()
  implicit none
  integer::unit
  character(4)::frmt
  unit=6
  frmt='(A)'
  write(unit,frmt)"Program SH_axisload calculates spherical harmonic expansions of axis symmetric loads."
  write(unit,frmt)"usage: SH_axisload [OPTIONS] "
  write(unit,frmt)"Where OPTIONS are"
  write(unit,frmt)" -l=lmax[,lmin] : specify maximum and optionally minimum degree (default 200,0)"
  write(unit,frmt)" -p=lon,lat: Specify location of the load function in degrees longitude and latitude"
  write(unit,frmt)" -p:LONLATFILE: Same as above but several lonlat combinations are read from the file and their"
 write(unit,frmt)"   contribution is summed"
  write(unit,frmt)" -u: Unit load"
  write(unit,frmt)" -d[asr]=WIDTH: Uniform disk load (spherical cap) "
  write(unit,frmt)"   -da=WIDTH: WIDTH represents radius angle in degrees"
  write(unit,frmt)"   -dr=WIDTH: WIDTH represents radius angle in km (along the surface of the Earth)"
  write(unit,frmt)"   -ds=WIDTH: WIDTH represents surface area in km^2 (on the surface of the Earth)"
  write(unit,frmt)""
  write(unit,frmt)" -c[asr]=WIDTH: Parabolic cap"
  write(unit,frmt)"   -ca=WIDTH: WIDTH represents radius angle in degrees"
  write(unit,frmt)"   -cr=WIDTH: WIDTH represents radius angle in km (along the surface of the Earth)"
  write(unit,frmt)"   -cs=WIDTH: WIDTH represents surface area in km^2 (on the surface of the Earth)"
  write(unit,frmt)""
  write(unit,frmt)" -i: output isotropic part only (-p option is ignored) degree and value"
  write(unit,frmt)"Program prints to standard output"
!   write(unit,frmt)""
!   write(unit,frmt)""
  stop
end subroutine help

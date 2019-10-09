!!Program to convert Spherical harmonics to the value over a point
!!Coded by Roelof Rietbroek, Fri Aug 22 12:24:28 2008
!!GeoForschunZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!! this program has been adapted from SH_2_nc
!! program prints to standard output (time and value) and takes only one location as input

!!When no SH filename(s) are provided the program reads from the standard input
!!Updated by Roelof Rietbroek, Mon Mar 16 12:20:25 2009
!! added switch to apply a Gaussian filter to the input

!!Updated by Roelof Rietbroek, Fri Feb 19 11:56:24 2010
!! added east and north derivate option

!!Updated by Roelof Rietbroek, Wed Mar 16 11:06:51 2011
!!allow unlimited amount of input files

!!Updated by Roelof Rietbroek, Thu Aug 23 10:34:33 2012
!! allow both east and north componets to be computed at once 
!! (output will have multiple columns in East and North direction)

!!Updated by Roelof Rietbroek, Mon Jan 13 10:02:26 2014
!! also allow an input file with multiple lon,lat points

!!updated By Roelof Rietbroek, Tuesday 21 March 2017, NUllified variable infiles to prevent segfault in various.f90 


 
program SH_2_pnt
use shtools
use shtlbx
use grdtlbx
use FORTtlbx
implicit none
integer::narg,i,ind,nfiles,maxfiles,lmax,lmin,lmaxdum,type,st,nd,pos
integer::l,m,stderr,itharg
parameter(maxfiles=300,stderr=0)
logical::stdin,limdeg,exi,error,smooth
character(200)::dum,lonlatfile
! character(len=200),dimension(maxfiles)::infiles
character(len=200),pointer::infiles(:)=>null()
double precision::lon,lat
double precision::tcentdum,radius
double precision,allocatable,dimension(:)::dp,p,tcent,Yc,Ys
double precision, allocatable,dimension(:,:)::clm,slm,slm_sig,clm_sig
integer::iargc
logical::east,north,grad,plain
double precision,allocatable,dimension(:,:,:)::out
double precision,pointer,dimension(:)::lonv=>null()
double precision,pointer,dimension(:)::latv=>null()
integer::noutcol,shft,nouterr
integer::npoints !!amount of lon lat points 
integer::unit
integer::chunk,memsize
double precision::dumd(2)
double precision::d2r
integer::j
integer:: last

!intialize defaults
unit=13
stdin=.false.
limdeg=.false.
error=.false.
smooth=.false.
nfiles=0
lmin=0
lmax=0
lon=9999.
lat=9999.
east=.false.
north=.false.
plain=.true.
noutcol=1 ! default 1 output column
nouterr=0
lonlatfile=''
npoints=0
chunk=3000
d2r=pi/180.d0
nullify(infiles)
!!!!!!!!!!!!! process command line arguments

!get number of command line arguments
narg=iargc()

!!allocate size for 100 files names
call realloc_ptr(infiles,100)

!check whether the last argument is a filename
call getarg(narg,dum)

if(dum(1:1) .eq. '-' .or. narg .eq. 0) then
   stdin=.true.
   nfiles=1
end if

!Now make a loop over the command line arguments and sort according to filename or options
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg) exit
   !get next argument
   call getarg(itharg,dum)

   if (dum(1:1) .ne. '-')then !add to file name array
      nfiles=nfiles+1
      if(size(infiles,1) < nfiles)call realloc_ptr(infiles,100)
      infiles(nfiles)=trim(dum)

      cycle!cycle loop after setting the file
   end if


!if the argument is an option
   select case(dum(2:2))
   case('l')!limit maximum (and minum degree)
      limdeg=.true.
       ind=index(dum,',')
      if(ind .ne. 0) then !also a min degree is specified
         read(dum(3:ind-1),*)lmax
         read(dum(ind+1:),*)lmin
      else
         read(dum(3:),*)lmax
      end if
   case('e')!also output an error'
      error=.true.
!   case('n')!apply a factor of 1/4pi to the coefficients (when  
     
   case('p')!LONGITUDE and LATITUDE of the considered point
      ind=index(dum,',')
      if(ind > 0)then ! only one point provided
         npoints=1
         allocate(lonv(npoints),latv(npoints))
         read(dum(3:ind-1),*)lon
         read(dum(ind+1:),*)lat
         !convert to radians
         lonv(1)=lon*d2r
         latv(1)=lat*d2r

      else !! read in lon lat pairs from a file
         lonlatfile=trim(dum(3:))
         
         allocate(lonv(chunk),latv(chunk))         
         memsize=chunk
         open(file=trim(lonlatfile),unit=unit,status='old')
         
         npoints=0
         read(unit=unit,iostat=last,fmt=*)dumd(1:2)
         do while ( last .eq. 0 )
            npoints = npoints + 1
               if ( npoints > memsize ) then
                  call realloc_ptr( lonv, chunk )
                  call realloc_ptr( latv, chunk )
                  memsize = memsize + chunk
               end if

            lonv( npoints )=dumd( 1 ) * d2r
            latv( npoints )= dumd( 2 ) * d2r
            read ( unit = unit, iostat = last, fmt = *) dumd(1:2)
         end do
         
         close(unit)
      end if

   case('N')
      north=.true.
   case('E')
      east=.true.      

   case('w')!apply gaussian smoothing
      smooth=.true.
      read(dum(3:),*)radius    
      
   case('h')
      call help()
   case default
      write(stderr,*)'unknown option selected: ',dum(1:2)
      call help()
   end select

end do

!input checks

if(npoints .eq. 0)then
   write(stderr,*)'ERROR: please specify -p option'
   stop
end if

!check for a sensible  range of lon lat points
do, i=1, npoints
   if(abs(lonv(i)) > 360.d0*d2r .or. abs(latv(i)) > 90.d0*d2r )then
      write(stderr,*)'ERROR: lon and/or lat are out of range',lonv(i)/d2r,latv(i)/d2r
      stop
   end if
end do

if(east .or. north)plain=.false.

if(east .and. north)then
   grad=.true.
   noutcol=2 ! add extra outputcolumn
end if

if(error)nouterr=noutcol ! add error columns

!check for possible pole singularity
if(east .and. abs(abs(lat)-pi/2.d0)<1d-1*pi/180.d0)then
   do , i = 1, npoints
      if ( abs( abs( latv( i )) - pi / 2.d0 ) < 1d-1 * pi / 180.d0 ) then
         write(stderr,*)"ERROR in SH_2_pnt: no eastward derivative is allowed in the vicinity of the pole"
         stop
      end if
   end do
end if


!now get some metadata from the first file
if(stdin)then
   call SH_readmeta(type=type,lmax=lmaxdum,tcent=tcentdum)

else
   
   call SH_readmeta(filen=trim(infiles(1)),type=type,lmax=lmaxdum)
end if

!give an error when lmax and lmin are not given but the type is of 0
if(type .eq. 0 .and. .not. limdeg)then
   write(*,*)'For SH files which are of the generic type the -l option must be used'
   call help()
end if

!redefine maximum degree if not explicitly set
if( .not. limdeg ) lmax = lmaxdum


!now allocate memory for the spherical harmonic coefficients and output
pos=SH_pos(lmax,lmax)
allocate(clm(pos,nfiles),slm(pos,nfiles),p(pos),tcent(nfiles))
if(north)allocate(dp(pos))
clm=0.d0
slm=0.d0
p=0.d0
tcent=0.d0
if(north)dp=0.d0
if(error)then
   allocate(clm_sig(pos,nfiles),slm_sig(pos,nfiles))
   clm_sig=0.d0
   slm_sig=0.d0
end if


!read in SH files uin one big data array
if(stdin)then
   tcent(1)=tcentdum
   if(error)then
      call SH_readgrav(clm=clm(:,1),slm=slm(:,1),clm_sig=clm_sig(:,1),slm_sig=slm_sig(:,1),type=type)
   else
      call SH_readgrav(clm=clm(:,1),slm=slm(:,1),type=type)
   end if
else
   do,i=1,nfiles
!      write(stderr,*)'reading ', trim(infiles(i))
      !get time tag from each file
      call SH_readmeta(filen=trim(infiles(i)),tcent=tcent(i))
    !  write(*,*)tcent(i)
      if(error)then
         call SH_readgrav(filen=trim(infiles(i)),clm=clm(:,i),slm=slm(:,i),clm_sig=clm_sig(:,i),slm_sig=slm_sig(:,i),type=type)
      else
         call SH_readgrav(filen=trim(infiles(i)),clm=clm(:,i),slm=slm(:,i),type=type)
      end if
   end do
end if



!!set coefficients below lmin to zero
if(lmin >0)then
   st=1  
   nd=SH_pos(lmin-1,lmin-1)
   clm(st:nd,:)=0.d0
   slm(st:nd,:)=0.d0
   if(error)then
      clm_sig(st:nd,:)=0.d0
      slm_sig(st:nd,:)=0.d0
   end if
end if


! apply Gaussian smoothing if required
if(smooth) then
   do,i=1,nfiles
      call SH_gaus(clm(:,i),radius*1.d3)
      call SH_gaus(slm(:,i),radius*1.d3)
   end do
if(error)then
   do,i=1,nfiles
      call SH_gaus(clm_sig(:,i),radius*1.d3)
      call SH_gaus(slm_sig(:,i),radius*1.d3)
   end do
end if
end if

!take the square of the standard deviations

if(error)then
   slm_sig=slm_sig**2
   clm_sig=clm_sig**2
end if


pos=SH_pos(lmax,lmax)

!!precompute base functions
allocate(Yc(pos),Ys(pos))


!compute legendre polynomials






!! allocate space for output
allocate(out(noutcol+nouterr,nfiles,npoints))
out=0.d0


do, j=1,npoints !loop over all points
   lon=lonv(j)
   lat=latv(j)

   if(north)then
      !calculate derivative of associated Legendre function 
      call PlmBar_d1(p=p,dp1=dp,lmax=lmax,z=sin(lat))
      !replace p with dp and apply Chain rule ( derivative wrt negative colatitude)
      !   p=dp*cos(lat)
   else
      !calculate legendre polynomial
      call PlmBar(p=p, lmax=lmax, z=sin(lat))
      !   if(east)p=p/cos(lat) ! scale with latitude dependent factor
   end if
   
   
   if(east)then
      Yc=0.d0
      Ys=0.d0
      do,l=0,lmax
         do,m=0,l
            ind=SH_pos(l,m)
            Yc(ind)=-m*p(ind)*sin(m*lon)/cos(lat)
            Ys(ind)=m*p(ind)*cos(m*lon)/cos(lat)
         end do
      end do

      shft=1 
      do,i=1,nfiles 
         out(shft,i,j)=dot_productblas(Yc,clm(:,i))+dot_productblas(Ys,slm(:,i))
         if(error)then
            out(shft+noutcol,i,j)=dot_productblas(Yc**2,clm_sig(:,i))+dot_productblas(Ys**2,slm_sig(:,i))
            out(shft+noutcol,i,j)=sqrt(out(shft+noutcol,i,j))
         end if
      end do
      
   end if

   if(north)then
      Yc=0.d0
      Ys=0.d0
      do,l=0,lmax
         do,m=0,l
            ind=SH_pos(l,m)
            Yc(ind)=dp(ind)*cos(m*lon)*cos(lat)
            Ys(ind)=dp(ind)*sin(m*lon)*cos(lat)
         end do
      end do
      
      shft=1
      if(east)shft=2
      do,i=1,nfiles 
         out(shft,i,j)=dot_productblas(Yc,clm(:,i))+dot_productblas(Ys,slm(:,i))
         if(error)then
            out(shft+noutcol,i,j)=dot_productblas(Yc**2,clm_sig(:,i))+dot_productblas(Ys**2,slm_sig(:,i))
            out(shft+noutcol,i,j)=sqrt(out(shft+noutcol,i,j))
         end if
      end do
      
   end if
   
   if(plain)then
      Yc=0.d0
      Ys=0.d0
      do,l=0,lmax
         do,m=0,l
            ind=SH_pos(l,m)
            Yc(ind)=p(ind)*cos(m*lon)
            Ys(ind)=p(ind)*sin(m*lon)
         end do
      end do
      
      shft=1
      do,i=1,nfiles 
         out(shft,i,j)=dot_productblas(Yc,clm(:,i))+dot_productblas(Ys,slm(:,i))
         if(error)then
            out(shft+noutcol,i,j)=dot_productblas(Yc**2,clm_sig(:,i))+dot_productblas(Ys**2,slm_sig(:,i))
            out(shft+noutcol,i,j)=sqrt(out(shft+noutcol,i,j))
         end if
         
      end do
      
   end if
   
end do ! over points

!! print results to standard output
do,i=1,nfiles 
   write(*,*)tcent(i),out(:,i,:)
end do


! if(error)then
!    do,i=1,nfiles 
!       out=dot_productblas(Yc,clm(:,i))+dot_productblas(Ys,slm(:,i))
!       out_sig=dot_productblas(Yc**2,clm_sig(:,i))+dot_productblas(Ys**2,slm_sig(:,i))
   
!    end do
! else
!    do,i=1,nfiles 
!       out=dot_productblas(Yc,clm(:,i))+dot_productblas(Ys,slm(:,i))
!       write(*,*)tcent(i),out
!    end do
! end if

end program SH_2_pnt





!!subroutine to plot help message
subroutine help()
implicit none
character(3)::frmt
integer::stderr
stderr=0

frmt='(A)'


write(stderr,frmt)'Convert spherical harmonic coefficients to the value at a geographical longitude-latitude positions'
write(stderr,frmt)'outputs timetag plus value to standard output'
write(stderr,frmt)'usage SH_2_pnt [OPTIONS] [filename(s)]'
write(stderr,frmt)'Where the filenames represent the spherical harmonic coefficient files'
write(stderr,frmt)'formats supported: '
write(stderr,frmt)'       1:GRACE(GFZ,JPL,CSR) type'
write(stderr,frmt)'       2:GINS type (French GRGS solution)'
write(stderr,frmt)'       3:ICGEM type'
write(stderr,frmt)'       4:(default) short type first line contains the degree supported, start, center end end time'
write(stderr,frmt)'              of observation'
write(stderr,frmt)'         the other lines are SH coefficients'
write(stderr,frmt)'       0:clean no header only SH coefficients'
write(stderr,frmt)'when no input files are provided the program reads the coefficient from standard input (in the format 4)'
write(stderr,frmt)''
write(stderr,frmt)'The OPTIONS can be the following:'
write(stderr,frmt)'  -lLMAX[,LMIN]: Specify maximum degree,LMAX (and optionally minimum degree LMIN) to consider'
write(stderr,frmt)'                 Default takes LMIN=0 and tries to find the maximum degree from the file'
write(stderr,frmt)''
write(stderr,frmt)'  -e: also propagate spherical harmonic coefficient errors (outputs three columns).'
write(stderr,frmt)'      NOTE: simple error propagation (no coefficient correlations assumed).'
write(stderr,frmt)'  -pLON,LAT: Mapping  position: longitude and latitude in degrees (OBLIGATORY argument)'
write(stderr,frmt)'  -pLONLATFILE: map spherical harmonic coefficients to multiple lon,lat points as read from'
write(stderr,frmt)'                the ascii file LONLATFILE.'
write(stderr,frmt)'                Output columns will be time, val_point1, ( err_point1 ), val_point2, ( err_point2 ), etc'
write(stderr,frmt)'  -wRAD: Apply gaussian smoothing with halfwidth radius RAD [km].'
write(stderr,frmt)'  -E: Differentiate the input in the East ward direction'
write(stderr,frmt)'  -N: Differentiate the input in the North ward direction'
write(stderr,frmt)'   When both East and north options are provided, the output will contain East and North columns'
write(stderr,frmt)''
write(stderr,frmt)'  -h: show this help message'
stop
end subroutine help



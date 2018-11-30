!program to create a GPS design matrix and stores it in BIN file format
!reads in an ascii file with 24 character names and longitude an latitude positions


!!Coded by Roelof Rietbroek, Wed Aug  6 09:50:30 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Fri May  8 19:31:59 2009
!!incorporated optional rotation matrix
!!input must now be given per station so station name and longitude and latitude

!!Updated by Roelof Rietbroek, Thu Dec 10 16:51:36 2009
!!input may be cartesian as well

!!Updated by Roelof Rietbroek, Fri Mar 26 15:52:12 2010
!!Frame of reference may be switched (CM,CF,CE,CL,CH)

!!Updated by Roelof Rietbroek, Tue Apr 13 15:15:49 2010
!! Helmert parameters may be applied only to specific observations

!!Updated by Roelof Rietbroek, Wed Aug 22 09:35:45 2012
!! added GIA Stokes Coefficients using the relationship by Purcell et al 2011
!! added horizontal poloidal SH expansion
!!Updated  24 April 2018 add option to compute loading from potential which also
!includes solid Earth loading effect

program GPS_designmat
use binfiletools
use gpstlbx
use shtlbx
use forttlbx
implicit none
integer::unit,i,ind,itharg,narg,lmax,lmin,nunk,nsh,ndat,stderr
integer::last,shift,chunk
parameter(stderr=0) ! maximum amount of station locations
parameter(chunk=500) !size of allocation chunks for input
character(200)::dum,file
type(BINdat)::output,rotmat
logical::helm,geocent,sh,rotate,cartesian
double precision::deg2rad,dumdbl(3)
double precision,allocatable,target::A(:,:)
double precision,allocatable,target::lonlat(:,:)
double precision,pointer,dimension(:)::lon=>null()
double precision,pointer,dimension(:)::lat=>null()
double precision,pointer,dimension(:)::height=>null()
double precision,pointer,dimension(:)::x=>null()
double precision,pointer,dimension(:)::y=>null()
double precision,pointer,dimension(:)::z=>null()
character(24),pointer::dat_d(:)=>null()
character(24),pointer::plate(:)=>null()
character(24),pointer::plate3(:)=>null()
character(24),target,allocatable::unk_d(:)
logical::platecol,eulerplate,potshdirect
integer::shft1,shftp,k,j
double precision,allocatable::tmp(:,:)
integer::iargc,ncol,nplate,ltyp,frame
integer::helm_inc,helm_exc
logical:: giash,horsh, potsh

!! defaults/initializations
potshdirect=.false.
potsh=.false.
platecol=.true.
eulerplate=.false.
unit=5 ! read from standard input by default
deg2rad=pi/180.d0
file=''
frame=3 ! defualt frame has origin in CF
output%file='stdout'
lmax=60
lmin=2
sh=.false.
geocent=.false.
helm=.false.
ltyp=2 !set default love numbers to elastic PREM 
rotate=.false.
helm_inc=0
helm_exc=0

giash=.false.
horsh=.false.

!!process command line options
narg=iargc()
if(narg < 1)call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)

   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('l')!limit maximum and minimum degree
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(3:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(3:),*)lmax
         end if
      case('t')!incorporate Helmert Transformation parameters in the normal equations
         helm=.true.
         if(dum(3:3) .ne. '')then
            select case(dum(3:4))
            case('i=')
               helm_inc=regcompf(trim(dum(5:)))
            case('i:')
               helm_inc=regcompf(trim(dum(5:)),.true.)
            case('e=')
               helm_exc=regcompf(trim(dum(5:)))
            case('e:')
               helm_exc=regcompf(trim(dum(5:)),.true.)
            case default
               write(stderr,*)"ERROR processing the -t argument:",trim(dum)
               stop
            end select
         end if
      case('g')!describe degree 1 variation as geocenter motion
         geocent=.true.
         
      case('s')!create system with Spherical harmonic elastic loading

         select case(dum(3:3))
         case('p')
             if(dum(4:4) == 's')then
                potsh=.true.
            else    
                potshdirect=.true.
            end if
         case('g')
            giash=.true.
         case('v')
            horsh=.true.
         case default
            sh=.true.
         end select

      case('p')!euler plates
         eulerplate=.true.
      case('r')!rotate from cartesian to local frame
         rotate=.true.
      case('F')! specify different filename
         output%file=dum(4:)
      case('E') !Use different load love numbers (different Earth model)
         read(dum(3:),'(i1)')ltyp
      case('O')!set a different origin of the reference frame
         select case(dum(3:4))
         case('CE')
            frame=1
         case('CM')
            frame=2
         case('CF')
            frame=3
         case('CL')
            frame=4
         case('CH')
            frame=5
         case default
            write(stderr,*)'unknown reference frame origin specified, quitting'
            stop
         end select
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
      file=trim(dum)
   end if
end do



!some input checks

if(.not.(helm .or. geocent .or.sh .or. rotate .or. potshdirect .or. potsh .or. eulerplate .or. giash .or. horsh))then
   write(stderr,*)'ERROR: specify either one or more of the options:'
   write(stderr,*)'   -r -s -g -t'
   write(stderr,*)'  no output'
   stop
end if

if(geocent .and. lmin < 2)then
   write(stderr,*)'WARNING: geocenter specified, changing lmin to 2'
   lmin=2
end if

if(potsh .and. potshdirect)then
   write(stderr,*)'ERROR cannot specify both -sp and -sps at the same time'
   stop
end if  

!open file if not provided read from standard input
if(file .ne. '')then
   unit=13
   open(unit=unit,file=file,form='formatted')
   
end if
last=0
ndat=1
allocate(dat_d(chunk))
allocate(plate(chunk))
allocate(x(chunk),y(chunk),z(chunk))
x=0.d0
y=0.d0
z=0.d0
dumdbl=0.d0
!try reading first values and determine the amount of columns
   !read description from the first 24 characters
read(unit,iostat=last,fmt='(A200)')dum
if(last .ne. 0)then ! error reading file
   write(stderr,*)"ERROR reading file:",trim(file)
   stop
end if
dat_d(ndat)=dum(1:24)
!try to read triplets from string
read(dum(25:),fmt=*,iostat=last)x(ndat),y(ndat),z(ndat)

if(last .eq. 0)then 
   ncol=3
   !try reading 3 values and a plate name
   read(dum(25:),fmt=*,iostat=last)x(ndat),y(ndat),z(ndat),plate(1)
   if(last .ne. 0)platecol=.false.
else ! try reading two values and a plate parameter
   ncol=2
   read(dum(25:),fmt=*,iostat=last)x(ndat),y(ndat),plate(1)
   if(last .ne. 0)platecol=.false.
end if

if( eulerplate .and. .not. platecol)then
   write(stderr,*)"ERROR: Euler plate option requires plate names as the final column in the input"
   stop
end if

!determine input type from range 
if(abs(x(1))>360 .or. abs(y(1))>90 .or. abs(z(1))>20000)then
   cartesian=.true.
else
   cartesian=.false.
end if

last=0
do while(last >= 0)!loop over remaining data lines
   ndat=ndat+1
   if(size(x,1)<ndat)then ! reallocate vectors ( not large enough)
      call realloc_ptr(x,chunk)
      call realloc_ptr(y,chunk)
      call realloc_ptr(z,chunk)
      call realloc_ptr(dat_d,chunk)
      call realloc_ptr(plate,chunk)
   end if
   !read description from the first 24 characters
   read(unit=unit,advance='NO',fmt='(a24)',iostat=last)dat_d(ndat)
   if(last < 0)then ! end of file encountered => exit while loop 
      ndat=ndat-1
      exit
   end if
   !read triplets (or twins ) and possibly a plate name
   if(platecol)then
      read(unit=unit,fmt=*,iostat=last)dumdbl(1:ncol),plate(ndat)
   else
      read(unit=unit,fmt=*,iostat=last)dumdbl(1:ncol)
   end if
   if(last .ne. 0)then
      write(stderr,*)"ERROR reading file, missing columns?"
      stop
   end if
   x(ndat)=dumdbl(1)
   y(ndat)=dumdbl(2)
   z(ndat)=dumdbl(3) !will stay zero if ncol==2
end do

if(file .ne. '')then
   close(unit)
end if


if(eulerplate)then ! determine amount of individual plates
   nplate=1
   do,i=2,ndat !loop over
      do,j=1,i-1
         if(plate(i) .eq. plate(j))exit
      end do
      if(j .eq. i)nplate=nplate+1
   end do
end if


!convert lon lat height to XYZ or vice versa
allocate(lon(ndat),lat(ndat),height(ndat))
if(cartesian)then
   do,i=1,ndat
      call ECEF_2_geodetic(x(i),y(i),z(i),lon(i),lat(i),height(i))
   end do
else
   !copy values to lon,lat and height
   lon(1:ndat)=deg2rad*x(1:ndat)
   lat(1:ndat)=deg2rad*y(1:ndat)
   height(1:ndat)=z(1:ndat)
!    do,i=1,ndat
!       call geodetic_2_ECEF(lon(i),lat(i),height(i),x(i),y(i),z(i))
!    end do
end if



!expand location to 3D vectors

allocate(rotmat%side1_d(ndat*3),rotmat%side2_d(ndat*3))
allocate(lonlat(3*ndat,2))

do,i=1,ndat
   rotmat%side1_d(1+(i-1)*3)='STAH   '//dat_d(i)(1:17)
   rotmat%side1_d(2+(i-1)*3)='STAE   '//dat_d(i)(1:17)
   rotmat%side1_d(3+(i-1)*3)='STAN   '//dat_d(i)(1:17)

   rotmat%side2_d(1+(i-1)*3)='STAX   '//dat_d(i)(1:17)
   rotmat%side2_d(2+(i-1)*3)='STAY   '//dat_d(i)(1:17)
   rotmat%side2_d(3+(i-1)*3)='STAZ   '//dat_d(i)(1:17)


   !longitude
   lonlat(1+(i-1)*3,1)=lon(i)
   lonlat(2+(i-1)*3,1)=lon(i)
   lonlat(3+(i-1)*3,1)=lon(i)

   !latitude
   lonlat(1+(i-1)*3,2)=lat(i)
   lonlat(2+(i-1)*3,2)=lat(i)
   lonlat(3+(i-1)*3,2)=lat(i)
   
end do

if(eulerplate)then ! copy plate names for triples (STAH,STAE,STAN)
   allocate(plate3(3*ndat))
   do,i=1,ndat
      plate3((i-1)*3+1)=plate(i)
      plate3((i-1)*3+2)=plate(i)
      plate3((i-1)*3+3)=plate(i)
   end do
end if


!calculate the amount of unknowns
nunk=0
nsh=SH_tpos(lmax,lmax,1,lmax,lmin)
if(sh)nunk=nunk+nsh
if(potshdirect .or. potsh)nunk=nunk+nsh
if(giash)nunk=nunk+nsh
if(horsh)nunk=nunk+nsh

if(geocent)nunk=nunk+3
if(helm)nunk=nunk+7
if(eulerplate)nunk=nunk+nplate*3


!point side description to the right side
output%side1_d=>rotmat%side1_d
! write(0,*)nunk,nsh,rotate
! stop

if(nunk > 0)then ! only create designmatrix if amount of unknowns > 0
   !make design matrix
   allocate(A(3*ndat,nunk)) 
   
  
   A=0.d0
   allocate(unk_d(nunk))
   !get observation equation matrix:
   shift=0
   if(helm)then
      call GPS_obseq_helmert(A=A(1:3*ndat,shift+1:shift+7),&
           tag=unk_d(shift+1:shift+7),typ=output%side1_d&
           ,lon=lonlat(:,1),lat=lonlat(:,2))
      !apply additional criteria
      if(helm_inc >0)then
         do,i=1,3*ndat
            if(.not. regexecf(helm_inc,output%side1_d(i)))A(i,shift+1:shift+7)=0.d0
         end do
      end if

      if(helm_exc >0)then
         do,i=1,3*ndat
            if(regexecf(helm_exc,output%side1_d(i)))A(i,shift+1:shift+7)=0.d0
         end do
      end if

      shift=shift+7
   end if
   
   if(geocent)then
      call GPS_obseq_geocent(A=A(1:3*ndat,shift+1:shift+3),&
           tag=unk_d(shift+1:shift+3),typ=output%side1_d,&
           lon=lonlat(:,1),lat=lonlat(:,2),ltyp=ltyp,frame=frame)
      shift=shift+3
   end if
   
   if(eulerplate)then
      call GPS_euler_plates(A=A(1:3*ndat,shift+1:shift+3*nplate)&
           ,tag=unk_d(shift+1:shift+3*nplate),typ=output%side1_d&
           ,plate=plate3,lon=lonlat(:,1),lat=lonlat(:,2))
      shift=shift+3*nplate
   end if

!   write(*,*)sh,ndat,shift,nsh,lmax,lmin
   if(sh)then
      
      call GPS_obseq_sh(A=A(1:3*ndat,shift+1:shift+nsh),tag=unk_d(shift+1:shift+nsh)&
           ,typ=output%side1_d,lon=lonlat(:,1),lat=lonlat(:,2),lmax=lmax,&
           lmin=lmin,ltyp=ltyp,shtyp=1,frame=frame)     
      shift=shift+nsh
   end if

   if(potshdirect)then
      call GPS_obseq_sh(A=A(1:3*ndat,shift+1:shift+nsh),tag=unk_d(shift+1:shift+nsh)&
           ,typ=output%side1_d,lon=lonlat(:,1),lat=lonlat(:,2),lmax=lmax,&
           lmin=lmin,ltyp=ltyp,shtyp=2,frame=frame)     
      shift=shift+nsh
   end if
   if(potsh)then
      call GPS_obseq_sh(A=A(1:3*ndat,shift+1:shift+nsh),tag=unk_d(shift+1:shift+nsh)&
           ,typ=output%side1_d,lon=lonlat(:,1),lat=lonlat(:,2),lmax=lmax,&
           lmin=lmin,ltyp=ltyp,shtyp=5,frame=frame)     
      shift=shift+nsh
   end if

   if(giash)then
      call GPS_obseq_sh(A=A(1:3*ndat,shift+1:shift+nsh),tag=unk_d(shift+1:shift+nsh)&
           ,typ=output%side1_d,lon=lonlat(:,1),lat=lonlat(:,2),lmax=lmax,&
           lmin=lmin,ltyp=ltyp,shtyp=3,frame=frame)     
      shift=shift+nsh

   end if

   if(horsh)then
      call GPS_obseq_sh(A=A(1:3*ndat,shift+1:shift+nsh),tag=unk_d(shift+1:shift+nsh)&
           ,typ=output%side1_d,lon=lonlat(:,1),lat=lonlat(:,2),lmax=lmax,&
           lmin=lmin,ltyp=ltyp,shtyp=4,frame=frame)     
      shift=shift+nsh
   end if

end if


!convert latitude and longitude back to degrees
lonlat=lonlat/deg2rad


!create rotation matrix if requested
if(rotate)then

   
   allocate(rotmat%pack1(ndat*3*3)) !matrix in packed form
   rotmat%pack1=0.d0
   
   !construct rotation matrix
   
   do,i=1,ndat
      shftp=(i-1)*9
      !packed order is first over rows then over columns

!       !for each station the rotation matrix looks like:
       
!           | cos(lat)cos(lon)     cos(lat)sin(lon)    sin(lat) |
!       R=  |    -sin(lon)              cos(lon)           0     |
!           |-sin(lat)cos(lon)    -sin(lon)sin(lat)    cos(lat) |

      !column 1
      rotmat%pack1(shftp+1)=cos(lat(i))*cos(lon(i))
      rotmat%pack1(shftp+2)=-sin(lon(i))
      rotmat%pack1(shftp+3)=-sin(lat(i))*cos(lon(i))
      !column 2
      rotmat%pack1(shftp+4)=cos(lat(i))*sin(lon(i))
      rotmat%pack1(shftp+5)=cos(lon(i))
      rotmat%pack1(shftp+6)=-sin(lon(i))*sin(lat(i))
      !colummn 3
      rotmat%pack1(shftp+7)=sin(lat(i))
      rotmat%pack1(shftp+8)=0.d0
      rotmat%pack1(shftp+9)=cos(lat(i))
   end do
   
   
   if(nunk .eq. 0)then ! write only rotation matrix to file
      !setup rotatation matrix ( note this matrix will only be written when no unknowns are provided)
      rotmat%file=output%file !use the same file
      rotmat%type='BDFULLVN'
      rotmat%descr='Rotation matrix by GPS_designmat'
      rotmat%mtyp='P' !in packed form
      rotmat%nval1=3*ndat
      rotmat%nval2=3*ndat
      rotmat%pval1=3*ndat*3
      rotmat%pval2=1   
      rotmat%nvec=2
      rotmat%vec=>lonlat
      
      !blockindices
      rotmat%nblocks=ndat
      allocate(rotmat%blockind(ndat))
      do,i=1,ndat
         rotmat%blockind(i)=i*3
      end do
      
      !meta data

      rotmat%nread=4
      allocate(rotmat%readme(rotmat%nread))
      rotmat%readme(1)="Rotation matrix which converts 3D vectors at stations"      
      rotmat%readme(2)="in the local Height, east, North frame"

      rotmat%readme(3)="Vectors contain longitude and latitude of station positions"
      rotmat%readme(4)="Matrix is packed and contains stacked 3 x 3 rotation matrices"

      call write_BINtype(rotmat)
      !stop program
      stop
   else !rotate designmatrix (using blas)
      allocate(tmp(3,nunk))
      do,i=1,ndat

         !copy data in temporary matrix
         shft1=(i-1)*3 ! shift in the station side
         shftp=(i-1)*9 ! shift in the packed matrix

         do,k=1,nunk
            do,j=1,3               
               tmp(j,k)=A(shft1+j,k)
            end do
         end do


!         DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         call dgemm('T','N',3,nunk,3,1.d0,rotmat%pack1(shftp+1),3,&
              tmp,3,0.d0,A(shft1+1,1),3*ndat)
      end do
      
      !adapt pointer of the output matrix
      output%side1_d=>rotmat%side2_d
   end if

end if



!set up output data
output%type='FULL2DVN' !full 2 dimensional matrix with vectors (latitude and longitude)
output%descr='made by GPS_designmat '
output%mtyp='F' ! writes a full matrix to the file
output%nval1=3*ndat
output%nval2=nunk
output%pval1=output%nval1*output%nval2
output%pval2=1
!integer metadata
output%nint=4
allocate(output%ints(output%nint),output%ints_d(output%nint))
output%ints_d(1)='Nobs' 
output%ints_d(2)='Nunknowns' 
if(sh .or. potshdirect )then
   output%ints_d(3)='Lmax'
   output%ints_d(4)='Lmin'
   output%ints(3)=lmax
   output%ints(4)=lmin
else
   output%nint=2
end if

output%ints(1)=ndat*3
output%ints(2)=nunk


!double metadata
output%ndbls=0

!readme meta data
output%nread=7
allocate(output%readme(output%nread))
                 !1                                                                     80!
output%readme=''
if(sh)then
   output%readme(1)='This file holds the design matrix which enables a conversion from mass'
   output%readme(2)='loading coefficients,helmert parameters and geocenter motion to GPS '
else if (potshdirect)then
   output%readme(1)='This file holds the design matrix which enables a conversion from Stokes'
output%readme(2)='coefficients,helmert parameters and geocenter motion to GPS '
end if

if(rotate)then
   output%readme(3)='station residuals (defined in the cartesian frame X,Y,Z).'
else
   output%readme(3)='station residuals (defined in the local Up, North, East frame).'
end if
if(sh)then
   output%readme(4)='Units of the normalized SH coefficients are in meters equivalent water'
else if (potshdirect)then
   output%readme(4)='normalized SH coefficients are dimensionless Stokes'
end if
output%readme(5)='height and geocenter motion and helmert parameters are expressed in'
output%readme(6)='meters. The GPS station residuals are also expressed in meters'
output%readme(7)='Created:'
call date_and_time(output%readme(7)(10:18),output%readme(7)(20:30))

output%nvec=2 ! longitude and latitude vectors are stored as well
output%vec => lonlat


output%side2_d => unk_d(1:nunk)
output%mat1 => A



!write to file
call write_BINtype(output)

end program GPS_designmat

subroutine help()
implicit none
integer::stderr
character(6)::frmt
frmt='(A)'
stderr=0
write(stderr,*)" Program GPS_desigmat creates a design matrix relating GPS station residuals"
write(stderr,*)" to loading parameter,"
write(stderr,*)" Helmert parameters and geocenter motion"
write(stderr,*)" Usage: GPS_desigmat [OPTIONS] [FILE]"
write(stderr,*)" File contains a 17 character description of the station name (STATION NAME LON LAT [HEIGHT]..)"
write(stderr,*)" and two[three] columns with longitude and latitude  an optionally height"
write(stderr,*)" in degrees (starting from column 25) and meters respectively"
write(stderr,*)" Alternatively 3 columns with cartesian coordinates may be provided in meters (X,Y,Z)"
write(stderr,*)" The input type will be automatically determined from the range of values"
write(stderr,*)" If FILE is not provided the program reads from standard input"
write(stderr,*)""
write(stderr,*)"OPTIONS may be:"
write(stderr,*)"-F=OUTPUTFILENAME: specify an output filename for the file (defaults writes to standard output)"
write(stderr,*)"-lLMAX[,LMIN]: specify the maximum and minimum degree for the SH expansion (default :lmax=60 lmin=2)"
write(stderr,*)"-t[i|e][:REGEXFILE|=REGEX]: create designmatrix with helmert parameters. Additional criteria can be"
write(stderr,*)"      supplied with 'i': include or 'e':exclude regular expression REGEX"
write(stderr,*)"      Or regular expressions from a text file REGEXFILE"
write(stderr,*)"      This will only set those entries to non-zero obeying the criteria"
write(stderr,*)"-g: create designmatrix with geocenter motion parameters"
write(stderr,*)"-s: create designmatrix with Spherical harmonic loading coefficients"
write(stderr,*)"-sp: create designmatrix with Spherical harmonic potential(Stokes) coefficients (load only)."
write(stderr,*)"-sps: create designmatrix with Spherical harmonic potential(Stokes) coefficients (load + solid earth)."
write(stderr,*)"-sg: create designmatrix with Spherical harmonic potential(Stokes) coefficients of GIA."
write(stderr,*)"     Using the relationship of Purcell et al 2011 between potential and vertical deformation"
write(stderr,*)"-sv: create designmatrix with Spherical harmonic Poloidal coeffients of the horizontal deformation"

write(stderr,*)"-p: create designmatrix with rotation vector of rigid tectonic Euler plates"
write(stderr,*)"    This option requires an additional column with PLATENAME as input"
write(stderr,*)"-r: rotate cartesian coordinates to local frame system first"
write(stderr,*)"  In case only -r is provided a pure rotation matrix will be constructed"
write(stderr,*)"-ENUM: use love numbers from Earth model NUM:"
write(stderr,*)'    1: Farrells load love numbers (GUTENBERG-BULLEN Earth model)'
write(stderr,*)'    2: Use PREM load love numbers (PREM Earth model) (default)'
write(stderr,*)'    3: Use Guttenberg-Bullen BODY love numbers (degree 2-4 only)'
write(stderr,*)'    4: Use PREM BODY love numbers (degree 2-4 only)'
write(stderr,*)'    5: Use Desai 2002  BODY love numbers (degree 2 only)'
write(stderr,*)""
write(stderr,*)"-OISOFRAME: Set the Origin of reference frame to ISOFRAME"
write(stderr,*)"            Applies to -g and -s ( degree 1)"
write(stderr,*)"            CM: Center of mass of the Earth system"
write(stderr,*)"            CF: Center of surface figure (default)"
write(stderr,*)"            CE: Center of Mass of the solid Earth"
write(stderr,*)"            CL: Center of surface Lateral figure"
write(stderr,*)"            CH: Center of surface Height figure"
stop
end subroutine help

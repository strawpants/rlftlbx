!!Program to convert Spherical harmonics to a netcdf grid file in the coards convention
!!When no SH filename is provided the program reads from the standard input
!!morer fuilenames can be provided to accelarate the computation (optimizes the use of legendre polynomials, but might eat some memory for many and large files)
!!Uses the netcdf library the SHTOOLS library the SHtlbx library and the GRDtlbx library

!!Coded by Roelof Rietbroek, Tue Jun 19 10:36:03 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue Jul  3 15:19:51 2007
!!incorporated error propagation with output to a separate grid

!!Updated by Roelof Rietbroek, Tue Oct  9 10:56:28 2007
!!incorporated translation of the origin (handy for shifting kernels over the earth)

!!Updated by Roelof Rietbroek, Thu Oct 18 11:18:32 2007
!! set default area to: lon =0.5..395.5 lat = -89.5..89.5
!!Updated by Roelof Rietbroek, Mon Aug 17 15:17:46 2009
!!Updated by Roelof Rietbroek, Thu Oct  8 17:09:53 2009
!! added option to switch between pixel and gridline registration ( for ease of use with GMT)
!! changed -d to -I and -r to -R to be consistent with GMT syntax

!!Updated by Roelof Rietbroek, Fri Feb 12 17:58:06 2010
!! allow appending of data to existing file (along the unlimited dimension using the same grid resolution/boundaries)

!!Updated by Roelof Rietbroek, Wed Feb 17 13:26:56 2010
!! allow differentation in north or south direction

!!Updated by Roelof Rietbroek, Tue Apr 13 16:45:51 2010
!! added '=' to command line options

!!Updated by Roelof Rietbroek, Mon Feb  6 14:54:25 2012
!! added OpenMP directives ( latitude loop can now be executed parallel). The linked BLAS library should be a single threaded version
 

!!Updated by Roelof Rietbroek, Thu Mar  8 17:28:46 2012
!! allowed unlimited amount of input files

!!Updated by Roelof Rietbroek, Fri Nov 27 11:48:26 2015
!! nullified infiles at initialization (is not done automatically by newer gfortran compiler)


program SH_2_nc
use shtools
use shtlbx
use netcdf
use grdtlbx
use forttlbx
implicit none
integer::narg,i,j,ind,indo,nfiles,lmax,lmin,lmaxdum,type,nlat,nlon,st,nd,pos !,maxfiles
integer::l,m,ldc,k
!parameter(maxfiles=400)
logical::stdin,limdeg,sissy,exi,error,translate
character(200)::dum,fileout,hist
! character(len=200),dimension(maxfiles)::infiles
character(len=200),pointer::infiles(:)=>null()
double precision::dgrd,lonmin,lonmax,latmax,latmin,mem,colat,time,time2,tcentdum,Tlon,Tlat
double precision,allocatable,dimension(:)::lat,lon,p,tcent,dp
double precision, allocatable,dimension(:,:)::clm,slm,Cvec,Svec,slm_sig,clm_sig,tmp1
double precision,allocatable,dimension(:,:,:)::grd,grd_sig
double precision::offset
integer::iargc,itharg
logical::pixel,append,closef
integer::stderr
logical::north,east
integer::chunk,memsz
!intialize defaults
stderr=0
closef=.true.
pixel=.true.
append=.false.
offset=0.d0
stdin=.false.
limdeg=.false.
sissy=.true.
error=.false.
nfiles=0
lmin=0
lmax=0
!gridresolution
dgrd=1.d0
!grid boundaries
lonmin=0.d0
lonmax=360d0
latmin=-90.d0
latmax=90.d0
translate=.false.
north=.false.
east=.false.
chunk=20
infiles=>null()

call realloc_ptr(infiles,chunk)

memsz=chunk
!!timing for performance measures
call cpu_time(time)
!!!!!!!!!!!!! process command line arguments

!get number of command line arguments
narg=iargc()
itharg=0


!Now make a loop over the command line arguments and sort according to filename or options

do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit! exit loop when last argument has been read in
   !get next argument
   call getarg(itharg,dum)

   if (dum(1:1) .ne. '-')then !add to file name array
      nfiles=nfiles+1

      if(memsz < nfiles)then
         call realloc_ptr(infiles,chunk)
         memsz=memsz+chunk         
      end if
      infiles(nfiles)=trim(dum)
!      write(0,*)trim(dum),memsz,nfiles
      cycle!cycle loop after setting the file
   end if


!if the argument is an option
   select case(dum(2:2))
   case('l')!limit maximum (and minum degree)
      limdeg=.true.
       ind=index(dum,',')
      if(ind .ne. 0) then !also a min degree is specified
         read(dum(4:ind-1),*)lmax
         read(dum(ind+1:),*)lmin
      else
         read(dum(4:),*)lmax
      end if
   case('I')!specify a different grid resolution in degrees (default is 1 degree)
      read(dum(4:),*)dgrd
   case('R')!specify a specific region (default is global)
      ind=index(dum,'/')
      indo=3
      read(dum(indo+1:ind-1),*)lonmin
      
      indo=ind
      ind=ind+index(dum(ind+1:),'/')
      read(dum(indo+1:ind-1),*)lonmax
      
      indo=ind
      ind=ind+index(dum(ind+1:),'/')
      read(dum(indo+1:ind-1),*)latmin

      read(dum(ind+1:),*)latmax
      
      !run some checks straight away
      if((latmax-latmin) <0.d0)then
         write(*,*)'minimum and maximum latitude swopped?: '
         call help()
      else if((lonmax-lonmin) <0.d0)then
         write(*,*)'minimum and maximum longitude swopped?: '
         call help()
      end if

   case('f')!force computation even if it eats up all your resources
      sissy=.false.
   case('e')!also output an error grid'
      error=.true.
!   case('n')!apply a factor of 1/4pi to the coefficients (when  
   case('N')
      north=.true.
   case('E')
      east=.true.
   case('F') !overwrite/append to file
      append=.true.
      closef=.false.
      if(dum(4:4) .eq. ' ')then
         itharg=itharg+1
         call getarg(itharg,fileout)
      else
         fileout=dum(4:)
      end if
      sissy=.false.
   case('t')!translate grid coordinates after mapping
      translate=.true.
      ind=index(dum,',')
      read(dum(4:ind-1),*)Tlon
      read(dum(ind+1:),*)Tlat
   case('g')! use gridline registration of the netcdf grid
      pixel=.false.
   case('h')
      call help()
   case default
      write(*,*)'unknown option selected: ',dum(1:2)
      call help()
   end select

end do

if(east .and. north)then
   write(stderr,*)"ERROR in SH_2_nc: only one of -E and -N is allowed per call"
   stop
end if


!check for possible pole singularity
if(east .and. .not. pixel)then
   if(abs(max(abs(latmin),abs(latmax))-90.d0)<1d-1)then ! to close to the pole!
      write(stderr,*)"ERROR in SH_2_nc: no eastward derivative is allowed in the vicinity of the pole"
      stop
   end if
end if

if(nfiles .eq. 0)then !read from standard input
   if(.not. append)then
      write(stderr,*)"Reading from standard input requires the -F option"
      stop
   end if
   stdin=.true.
   nfiles=1
end if


!check if files already exist
if(.not. append)then
   do,i=1,nfiles
      inquire(FILE=(trim(infiles(i))//'.nc'), EXIST=exi)
      if(exi .and. sissy)then
         write(*,*)'file exists, use -f or -FOUTPUTNAME:',trim(infiles(i))//'.nc'
         call help()
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
if(.not.limdeg)lmax=lmaxdum

!!calculate vector length of latutude and longitude
if(pixel)then
   nlat=int((latmax-latmin)/(dgrd))
   nlon=int((lonmax-lonmin)/(dgrd))
   offset=dgrd/2
else
   nlat=int((latmax-latmin)/(dgrd))+1
   nlon=int((lonmax-lonmin)/(dgrd))+1
end if

!amount of SH coefficients
pos=SH_pos(lmax,lmax)

!! Now run a 'back of the envelope' calculation of how much memory will be used and give a warning if it is more than 400Mb
if (sissy)then
   if(error)then
      mem=(nfiles*4.d0*pos+pos*nlon*2.d0+2*nfiles*nlat*nlon+nlat+nlon+pos)*8.d0/(1024.d0**2)
   else
      mem=(nfiles*2.d0*pos+pos*nlon*2.d0+nfiles*nlat*nlon+nlat+nlon+pos)*8.d0/(1024.d0**2)
   end if
   if(mem> 500.d0)then
      write(*,*)'This program is going to eat around',int(mem),'Mb of your memory, consider smaller batchjobs or larger resolution'
      write(*,*)'Or do not be a sissy and use -f to force the run'
      write(*,*)''
      call help()
   end if

end if

!now allocate memory for the spherical harmonic coefficients

!allocate(clm(pos,nfiles),slm(pos,nfiles),p(pos),tcent(nfiles))
allocate(clm(pos,nfiles),slm(pos,nfiles),tcent(nfiles))
!if(north)allocate(dp(pos)) ! also allocate space for derivative
clm=0.d0
slm=0.d0
tcent=0.d0
if(error)then
   allocate(clm_sig(pos,nfiles),slm_sig(pos,nfiles))
   clm_sig=0.d0
   slm_sig=0.d0
end if


!read in SH files in one big data array
if(stdin)then
   tcent(1)=tcentdum
   if(error)then
      call SH_readgrav(clm=clm(:,1),slm=slm(:,1),clm_sig=clm_sig(:,1),slm_sig=slm_sig(:,1),type=type)
   else
      call SH_readgrav(clm=clm(:,1),slm=slm(:,1),type=type)
   end if
else
   do,i=1,nfiles
      write(*,*)'reading ', trim(infiles(i))
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

write(*,*)'done reading'

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



!take the square of the standard deviations

if(error)then
   do,k=1,nfiles
      do,j=1,pos
         clm_sig(j,k)=clm_sig(j,k)**2
      end do
   end do
   do,k=1,nfiles
      do,j=1,pos
         slm_sig(j,k)=slm_sig(j,k)**2
      end do
   end do

!write(*,*)clm_sig
end if
!!now construct the latitude and longitude vectors 




allocate(lat(nlat),lon(nlon),grd(nlon,nlat,nfiles),Cvec(nlon,pos),Svec(nlon,pos))
lat=0.d0
lon=0.d0
grd=0.d0

if(error)then
   allocate(grd_sig(nlon,nlat,nfiles))
   grd_sig=0.d0
end if
!construct the vectors with respect to their minimum value

do,i=0,nlat-1
   lat(i+1)=latmin+i*dgrd+offset
end do

do,i=0,nlon-1
   lon(i+1)=lonmin+i*dgrd+offset
end do
!write(*,*)lonmin,lonmax,latmin,latmax,lon(5)

!!precompute cos(m longitude) terms and sin ( m longitude) terms
!!this accelerates the computation significantly for large batch jobs
!!(at the cost of some memory though)

if(east)then ! derivate wrt longitude
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,m,i)
   do,l=0,lmax
      do,m=0,l
         ! do,i=1,nlon
         !    Cvec(i,SH_pos(l,m))=-m*sin(m*pi/180.d0*lon(i))
         !    Svec(i,SH_pos(l,m))=m*cos(m*pi/180.d0*lon(i))
         ! end do
         Cvec(:,SH_pos(l,m))=-m*sin(m*pi/180.d0*lon(:))
         Svec(:,SH_pos(l,m))=m*cos(m*pi/180.d0*lon(:))

      end do
   end do
   !$OMP END PARALLEL DO
else
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,m,i)
   do,l=0,lmax
      do,m=0,l
         ! do,i=1,nlon
         !    Cvec(i,SH_pos(l,m))=cos(m*pi/180.d0*lon(i))
         !    Svec(i,SH_pos(l,m))=sin(m*pi/180.d0*lon(i))
         ! end do
         Cvec(:,SH_pos(l,m))=cos(m*pi/180.d0*lon(:))
         Svec(:,SH_pos(l,m))=sin(m*pi/180.d0*lon(:))

      end do
   end do
   !$OMP END PARALLEL DO
end if


!BLAS TRICKING parameter let BLAS routines think that input array grd is 2D with leading dimension ldc (vector containing all  grid points)
ldc=nlat*nlon 

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(clm,slm,clm_sig,slm_sig,lat,lon,grd,grd_sig,Cvec,Svec) &
!$OMP FIRSTPRIVATE(ldc,north,error,east,nlon,nlat,pos,nfiles,lmax,lmin)
!!allocate private arrays (for each thread)
allocate(tmp1(nlon,pos))
allocate(p(pos))
if(north)allocate(dp(pos)) ! also allocate space for derivative
!!now evaluate the spherical harmonics on the grid
!$OMP DO
do,i=1,nlat !first do the loop over latitude(to avoid computing too many legendre polynomials)
   
   if(north)then
      !calculate derivative of associated Legendre function 
      call PlmBar_d1(p=p,dp1=dp,lmax=lmax,z=sin(lat(i)*pi/180.d0))
      !replace p with dp and apply Chain rule ( derivative wrt negative colatitude)
      p=dp*cos(lat(i)*pi/180.d0)
   else
      !calculate legendre polynomial
      call PlmBar(p=p, lmax=lmax, z=sin(lat(i)*pi/180.d0))
      if(east)p=p/cos(lat(i)*pi/180.d0) ! scale with latitude dependent factor
   end if

!create matrix for base functions (COSINE part)
   do,k=1,pos
      ! do,j=1,nlon
      !    tmp1(j,k)=Cvec(j,k)*p(k)
      ! end do
      tmp1(:,k)=Cvec(:,k)*p(k)
   end do

   !call blas multiplication
   call dgemm('N','N',nlon,nfiles,pos,1.d0,tmp1(1,1),nlon,clm(1,1),pos,0.d0,grd(1,i,1),ldc)

   if(error)then ! simple error propagation
      !square matrix tmp1
      do,k=1,pos
         do,j=1,nlon
            tmp1(j,k)=tmp1(j,k)**2
         end do
      end do
      
      !propagate diagonal error
      call dgemm('N','N',nlon,nfiles,pos,1.d0,tmp1(1,1),nlon,clm_sig(1,1),pos,0.d0,grd_sig(1,i,1),ldc)
   end if


!!SINE part !!!!!!!!!!!!
   !reuse tmp1
   do,k=1,pos
      ! do,j=1,nlon
      !    tmp1(j,k)=Svec(j,k)*p(k)
      ! end do
      tmp1(:,k)=Svec(:,k)*p(k)
      
   end do
   !call blas multiplication (update grd)
   call dgemm('N','N',nlon,nfiles,pos,1.d0,tmp1(1,1),nlon,slm(1,1),pos,1.d0,grd(1,i,1),ldc)

   if(error)then ! simple error propagation
      !square matrix tmp1
      do,k=1,pos
         do,j=1,nlon
            tmp1(j,k)=tmp1(j,k)**2
         end do
      end do
      
      !propagate diagonal error (update grd)
      call dgemm('N','N',nlon,nfiles,pos,1.d0,tmp1(1,1),nlon,slm_sig(1,1),pos,1.d0,grd_sig(1,i,1),ldc)
   end if

end do
!$OMP END DO
deallocate(tmp1,p)
if(north)deallocate(dp)
!$OMP END PARALLEL

!!now deallocate some stuff and write the data to netcdf files
deallocate(clm,slm)
if(error)deallocate(clm_sig,slm_sig) 

!!take the square root f the error grdi to obtain standard deviation
if(error)grd_sig=sqrt(grd_sig)


!translate coordinate variables if requested
if(translate)then
   lat=lat+Tlat
   lon=lon+Tlon
end if

if(append)then
   hist='Made by SH_2_nc from spherical harmonics (appended data)'
else
   hist='Made by SH_2_nc from spherical harmonics'
end if

if(east)hist=trim(hist)//"Eastward derivative"
if(north)hist=trim(hist)//"Northward derivative"
write(*,*)'Writing netcdf file(s)'
do,i=1,nfiles ! write data to file(s)
   if(i.eq. nfiles)closef=.true. !make sure last loop closes the file
   if (.not. append)fileout=trim(infiles(i))//'.nc'
   call NC_write(file=fileout,lat=lat,lon=lon,z=grd(:,:,i),&
        history=hist,ctime=tcent(i),&
        pix=pixel,napp=append,nclse = closef)
   if(error)call NC_write(file='err_'//trim(fileout),lat=lat,lon=lon,z=grd_sig(:,:,i)&
        ,history=hist,ctime=tcent(i),&
        pix=pixel,napp=append,nclse=closef)
end do



call cpu_time(time2)
write(*,*)'succesfully finished job in',time2-time,'seconds'

end program SH_2_nc





!!subroutine to plot help message
subroutine help()
implicit none
character(3)::frmt


frmt='(A)'


write(*,frmt)'Convert spherical harmonic coefficients to a geographical longitude-latitude netcdf grid'
write(*,frmt)'output file complies with the netcdf coards convention'
write(*,frmt)'usage SH_2_nc [OPTIONS] [filename(s)]'
write(*,frmt)'Where the filenames represent the spherical harmonic coefficient files'
write(*,frmt)'formats supported: '
write(*,frmt)'       1:GRACE(GFZ,JPL,CSR) type'
write(*,frmt)'       2:GINS type (French GRGS solution)'
write(*,frmt)'       3:ICGEM type'
write(*,frmt)'       4:(default) short type first line contains the degree supported, start, center end end time of observation'
write(*,frmt)'         the other lines are SH coefficients'
write(*,frmt)'       0:clean no header only SH coefficients'
write(*,frmt)'when no input files are provided the program reads the coefficient from standard input (in the format 4)'
write(*,frmt)''
write(*,frmt)'The OPTIONS can be the following:'
write(*,frmt)'  -l=LMAX[,LMIN]: Specify maximum degree,LMAX (and optionally minimum degree LMIN) to consider'
write(*,frmt)'                 Default takes LMIN=0 and tries to find the maximum degree from the file'
write(*,frmt)''
write(*,frmt)'  -I=RES: set the resolution, RES, of the output grid in degrees (default is 1 degree)'
write(*,frmt)''
write(*,frmt)'  -R=LONMIN/LONMAX/LATMIN/LATMAX: specify the region of the output grid in degrees:'
write(*,frmt)'                                 where LONMIN,LONMAX,LATMIN,LATMAX denote the minimum and maximum longitude'
write(*,frmt)'                                 and minimum and maximum latitude respectively'
write(*,frmt)'                                 default (LONMIN=0, LONMAX=360, LATMIN=-90,LATMAX=90)'
write(*,frmt)''
write(*,frmt)'  -f: force the run even if the program is going to use a lot of memory forces also to overwrite netcdffiles'
!write(*,frmt)'  -a: Append all data in the (possibly existing) file ( implies also -f)'
write(*,frmt)'  -E: Differentiate the input in the East ward direction'
write(*,frmt)'  -N: Differentiate the input in the North ward direction'
write(*,frmt)''
write(*,frmt)'  -F=FILEOUT: Specify output netcdf file to append to.'
write(*,frmt)'  -e: also propagate spherical harmonic coefficient errors (outputs in a seperate grid with err_ prepended)'
write(*,frmt)'  -t=TLON,TLAT: translate the origin after mapping by TLON and TLAT degrees.'
write(*,frmt)'                Note the -r option holds for the original coordinate system'
write(*,frmt)'  -g: Give the output file gridline registration ( used in gmt). Default is pixel.'
write(*,frmt)''
write(*,frmt)''
write(*,frmt)'  -h: show this help message'
!$ write(*,frmt)'This version is compiled with OpenMP and allows multi-threading (please set the OMP_NUM_THREADS env.variable)'
stop
end subroutine help



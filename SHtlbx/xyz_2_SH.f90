!!Coded by Roelof Rietbroek, Wed Jun 15 09:11:47 2011
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!! program which calculates Spherical harmonic coefficients from x y z triplets

!!Updated by Roelof Rietbroek, Tue Jan 24 14:28:05 2012
!! added OMP directives which enable multithreading of the outer dataloop
!!Updated by Roelof Rietbroek, Wed Jan 25 10:54:41 2012
!! added weighted area option

!!Updated by Roelof Rietbroek, Tue Nov 20 17:34:18 2012
!! added extra hint in the help file (units steradian)


program xyz_2_SH
use shtools
use shtlbx
use forttlbx
implicit none
integer::narg,iargc,itharg
character(200)::dum,fin
integer::stderr
parameter(stderr=0)
integer::lmax,lmin
integer::unit,last
integer::i,j,l,m
integer::chunk
integer::ndat
integer::memsize
double precision::dumd(4)
double precision,pointer::lon(:),lat(:),z(:)
double precision::d2r
double precision::latold,thres
double precision:: Cosm,Sinm
integer::otyp,ind,npos
integer,allocatable::perm(:)
double precision,allocatable,dimension(:)::clm,slm,p,clmtmp,slmtmp
double precision::lontmp,lattmp,ztmp
logical::verbose,weight
integer::ncol
character(60)::threadmess
!$ integer::OMP_get_thread_num
!$ integer::OMP_get_num_threads
integer::threadload
integer::threadstart

!defaults/initializations
lmax=100
lmin=0
unit=5 !standard input
fin="stdin"
chunk=3000 !allocate with chunks of 3000 long
d2r=pi/180.d0
thres=-1 !set threshold for latitude check ( negative means no sorting)
otyp=4
verbose=.false.
weight=.false.
ncol=3

!!process command line options
narg=iargc()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('h')
         call help()
      case('l')
         ind=index(dum,",")
         if(ind .ne. 0)then
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
      case('t')!specify different output type
         read(dum(3:),*)otyp 
      case('s')!set threshold and sorting option
         read(dum(3:),*)thres
         thres=thres*d2r
      case('v')
         verbose=.true.
      case('w')
         weight=.true.
         ncol=4
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option but a file name
      fin=trim(dum)
      unit=13
   end if
end do


!read data from file/standard input
if(fin .ne. "stdin")open(file=trim(fin),unit=unit,status='old')
allocate(lon(chunk),lat(chunk),z(chunk))
memsize=chunk
last=0
ndat=0

do while(last .eq. 0)
   read(unit=unit,iostat=last,fmt=*)dumd(1:ncol)
   if(last .ne. 0)exit
   ndat=ndat+1
   if(ndat > memsize)then
      call realloc_ptr(lon,chunk)
      call realloc_ptr(lat,chunk)
      call realloc_ptr(z,chunk)
      memsize=memsize+chunk
   end if
   lon(ndat)=d2r*dumd(1)
   lat(ndat)=d2r*dumd(2)
   z(ndat)=dumd(3)
   if(weight)z(ndat)=z(ndat)*dumd(4)! also weight the grid point
end do
if(ndat .eq. 0)then
   write(stderr,*)"ERROR: no data found or reading error"
   if(weight)write(stderr,*)"Did you forget the weights in the 4th column?"
   stop
end if
if(verbose)write(stderr,*)"read",ndat,"triplets"

if(fin .ne. "stdin")close(unit)


!!sort data according to latitude (speeds up calculation when more values with the same latitude are present)
! get a permutation vector to sort along latitude
if(thres >0)then
   if(verbose)write(stderr,*)"calculate permutation vector"
   allocate(perm(ndat))
   call sort_f90(perm,lat(1:ndat))
   !permute data
   call permute(lat(1:ndat),perm)
   call permute(lon(1:ndat),perm)
   call permute(z(1:ndat),perm)
end if

!divide z by 4 pi (for the normalization)
z(1:ndat)=z(1:ndat)/(4*pi)
if(.not. weight)z(1:ndat)=z(1:ndat)/ndat ! normalize with respect to the total amount of data points
!allocate space for SH vector
npos=SH_pos(lmax,lmax)


allocate(clm(npos),slm(npos))
clm=0.d0
slm=0.d0
latold=-9999.d0
threadload=ndat
threadstart=0
if(verbose)write(stderr,*)"calculate coefficients.."
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(clm,slm,lat,lon,z) &
!$OMP FIRSTPRIVATE(latold,lmax,lmin,ndat,thres,npos,verbose,d2r)
allocate(p(npos),clmtmp(npos),slmtmp(npos))
clmtmp=0.d0
slmtmp=0.d0
threadmess="VERBOSE:"
!$ threadload=int(dble(ndat)/OMP_get_num_threads()+0.5)
!$ threadstart=threadload*OMP_get_thread_num()
!$ if(threadload*(OMP_get_thread_num()+1)>ndat) &
!$ threadload=threadload-(threadload*(OMP_get_thread_num()+1)-ndat)
!$ write(threadmess,'(A20,I2)')"VERBOSE: Thread Nr:",OMP_get_thread_num()
!$OMP DO
do,i=1, ndat
   lontmp=lon(i)
   lattmp=lat(i)
   ztmp=z(i)
   if(abs(latold-lattmp)> thres)then
      if(verbose)write(stderr,*)trim(threadmess),i-threadstart,"of",threadload,"(re)calculating Plm for latitude",lattmp*d2r
      call PlmBar(p=p,lmax=lmax,z=sin(lattmp))
   end if
   do,m=0,lmax
      Cosm=ztmp*cos(m*lontmp)
      Sinm=ztmp*sin(m*lontmp)
      do,l=max(m,lmin),lmax
         j=SH_pos(l,m)
         clmtmp(j)=clmtmp(j)+p(j)*Cosm
         slmtmp(j)=slmtmp(j)+p(j)*Sinm
      end do
   end do
   latold=lat(i)!reset latitude
end do
!$OMP END DO
!gather clmtmp and slmtmp in the output vector
!$OMP CRITICAL
clm=clm+clmtmp
slm=slm+slmtmp
!$OMP END CRITICAL
deallocate(p,slmtmp,clmtmp)
!$OMP END PARALLEL

!write goodies to standard output
call SH_write(clm=clm,slm=slm,typ=otyp)





end program xyz_2_SH

subroutine help()
  implicit none
  integer::stderr
  write(stderr,*)"Program xyz_2_SH calculates spherical harmonic coefficients from xyz triplets"
  write(stderr,*)"Usage: xyz_2_SH [OPTIONS] [FILE]"
  write(stderr,*)"Where FILE contains rows with longitude, latitude and value"
  Write(stderr,*) "When FILE is not given the program will read from standard input" 
  write(stderr,*)"The resulting coefficients are printed to standard output"
  write(stderr,*)""
  write(stderr,*)"OPTIONS may be:"
  write(stderr,*)"-l=lmax[,lmin]: use maximum and mininum degree for the output ( default lmax=100, lmin=0)"
  write(stderr,*)'-tNUM: Specify the output type NUM:'
  write(stderr,*)'       1:GRACE(GFZ,JPL,CSR) type'
  write(stderr,*)'       2:GINS type (French GRGS solution)'
  write(stderr,*)'       3:ICGEM type'
  write(stderr,*)'       4:(default) short type first line contains the degree supported, start, center end end time of observation'
  write(stderr,*)'-sTHRES : Sort data according to latitude and set threshold for the latitude change (in degree) for which'
  write(stderr,*)'          a new associated Legendre function should be calculated' 
  write(stderr,*)'-w: Weighing: a 4th column contains the weights of the measurements, weights must be in steradian '
  write(stderr,*)'-v: be verbose'
!$ write(stderr,*)'This version is compiled with OpenMP and allows multi-threading (please set the OMP_NUM_THREADS env.variable)'
write(stderr,*)"See also nc_2_SH for quick FFT methods on equidistant grids"
stop
end subroutine help

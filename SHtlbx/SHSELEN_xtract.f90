!! program to extract normalized real SH (& rates) from binary files of SELEN
!!Coded by Roelof Rietbroek, Wed Oct 15 10:47:59 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Wed Sep 28 15:21:08 2011
!!allowed reading of both double and single precision data


program SHSELEN_xtract
use SHTlbx ! contains SH_pos and SH_write
use FORTtlbx
implicit none
integer::itharg,i,narg,stderr,n,nmax,lmax,sz,l,m
character(200)::dum,file
logical::rate,formatted
double precision::deltat,sqt2
complex*16,allocatable,dimension(:,:)::dat
complex*8,allocatable,dimension(:,:)::dat_s
double precision,allocatable,dimension(:)::clm,slm,vec
integer::datbyte
integer::iargc,fid,opencstream
integer*4::recdum
real::fdum1,fdum2
double precision::ddum1,ddum2
!defaults/initializations
datbyte=16 !amount of bytes in one complex value
stderr=0
rate=.false.
lmax=0
n=-1
nmax=0
formatted=.false.
!get command line arguments
narg=iargc()

if(narg < 4)call help()



itharg=0

do,i=1,narg
   itharg=itharg+1
   if(itharg >narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('r')
         rate=.true.
         read(dum(3:),*)deltat
      case('t')
         read(dum(3:),*)n
      case('l')
         read(dum(3:),*)lmax
      case('n')
         read(dum(3:),*)nmax
      case('f')!file is formatted
         formatted=.true.
      case('s')
         datbyte=8 ! single precision
      case default !unknown option
         write(stderr,*)'Unknown option:',trim(dum)
         stop
      end select
   else !file name
      file=trim(dum)
      
   end if

end do

!input checks
if(lmax==0)then
write(stderr,*)'Lmax must be supplied'
call help()
end if

if(nmax==0)then
write(stderr,*)'maximum amount of time steps NMAX must be supplied'
call help()
end if

if(n<0)then
write(stderr,*)'Time index must be supplied'
call help()
end if

!allocate data array
sz=SH_pos(lmax,lmax)
if( datbyte .eq. 16)then
   allocate(dat(sz,0:nmax-1))
else
   allocate(dat_s(sz,0:nmax-1))
end if

if(formatted)then
!read in data from formatted file
   open(3,file=file,status='unknown')
   read(3,*)dat
   close(3)
else
   !read in data from binary file
   !use C stream access
   ! open(3,file=file,status='unknown',form='unformatted')
   ! read(3)dat
   ! close(3)
   fid=opencstream(0,trim(file)//char(0))
   ! read fortran record marker (integer*4, but platform dependent)
   call cread(fid,4,recdum)
!   write(*,*)" Start",recdum
   if(recdum .ne. nmax*sz*datbyte)then
      write(stderr,*)"Data size from file does not appear to match:"
      write(stderr,*)"Possible causes:"
      write(stderr,*)"nmax and lmax are wrong"
      write(stderr,*)"Endiannes of file is reversed"
      write(stderr,*)"The unit of fortran record marker of the current system does"
      write(stderr,*)"not match that of the system where the file was created ( e.g. recordmarker in bytes versus 4*byte word)"
      write(stderr,*)"requested:",nmax*sz*datbyte,"bytes,available:",recdum
!      write(stderr,*)"with nmax, sz, datbyte, recdum:",nmax,sz,datbyte,recdum
      if(recdum .eq. (nmax+1)*sz*datbyte)then
         write(stderr,*)"HINT: changing NMAX to",nmax+1,"would fit"
      end if
      if(recdum .eq. (nmax)*sz*datbyte/2)then
         write(stderr,*)"HINT: Assuming single precision  would fit, use the -s flag"
      end if

      stop
   end if

!read remaining data:
   ! call cread(fid,8,ddum1)
   ! call cread(fid,8,ddum2)
   ! write(*,*)ddum1,ddum2
   if(datbyte .eq. 16)then
      call cread(fid,nmax*sz*datbyte,dat(1,0))
   else
      call cread(fid,nmax*sz*datbyte,dat_s(1,0))
   end if

   call cread(fid,4,recdum)
!   write(*,*)"after",recdum
   call cclose(fid)
   ! do,l=1,128
   !    write(*,*)l,dat(l,0)
   ! end do
   ! 
   ! write(*,*)recdum

   ! call cswap(4,1,recdum)
   ! write(*,*)"swapped",recdum
 
end if

!allocate real coefficient vectors
allocate(clm(sz),slm(sz))
clm=0
slm=0

if(rate)then
   if(datbyte .eq. 16)then
      dat(:,n)=(dat(:,n)-dat(:,n-1))/deltat ! replace column vector by its difference
   else
      dat_s(:,n)=(dat_s(:,n)-dat_s(:,n-1))/deltat ! replace column vector by its difference
   end if
end if

! !set entries with degree < 2  explicitly to zero
! dat(1:SH_pos(1,1),n)=0.d0

!precompute vector which implements condon Shortley Phase and sqrt(2-delta_om)
allocate(vec(sz))
vec=1
sqt2=sqrt(2.d0)


do,l=1,lmax
   do,m=1,l,2 !odd orders
      vec(SH_pos(l,m))=-sqt2
 !     if(l==4)write(*,*)'odd',l,m,lmax 
   end do
   do,m=2,l,2 !even orders
      vec(SH_pos(l,m))=sqt2 
!      if(l==4)write(*,*)'even',l,m,lmax 
   end do
end do

if(datbyte .eq. 16)then
   clm=dble(dat(:,n))*vec
   clm=clm
   slm=dble(aimag(dat(:,n)))*vec
   slm=-slm
else
   clm=dble(dat_s(:,n))*vec
   clm=clm
   slm=dble(aimag(dat_s(:,n)))*vec
   slm=-slm
end if

!write to standard output
call SH_write(clm=clm,slm=slm)


end program SHSELEN_xtract

subroutine help()
implicit none
integer::unit
character(8)::frmt
unit=6
frmt='(A)'

write(unit,frmt)'Program SHSELEN_xtract extracts real spherical harmonics and rates from SELEN files'
write(unit,frmt)'Input files are the SELEN binary/formatted output files containing complex spherical harmonics'
write(unit,frmt)'Usage: SHSELEN_xtract [OPTIONS] -tTIME -lLMAX -nNMAX FILE'
write(unit,frmt)'Where OPTIONS may be:'
write(unit,frmt)'  -rDELTAT: calculate approximate rate between time and previous time step (with time step DELTAT in years)'
write(unit,frmt)'  -f: Input file is a formatted file (eg shice.dat)'
write(unit,frmt)'  -s: Data are stored in single precision)'
write(unit,frmt)''
write(unit,frmt)'The parameter LMAX and NMAX denote the maximum degree and amount of time steps respectively'
write(unit,frmt)'These must be consistent with the FILE'
write(unit,frmt)'One can consult file data.inc for the NN variable which is NMAX-1'
write(unit,frmt)'TIME denotes the integer index correponding to the desired time (starting at 0)'
write(unit,frmt)'Results are printed to standard output'
stop

end subroutine help

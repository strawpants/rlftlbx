!!Coded by Roelof Rietbroek, Wed Aug 10 16:25:53 2011
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de


program JKcov_2_BIN
use shtlbx
use forttlbx
use binfiletools
implicit none
integer::i,j,narg,iargc,stderr,itharg
character(200)::dum,file
integer::lmax,lmin,ind,q
integer::mrkr,sz,l,m,unit
integer*8::fsize,pcksz,tmpsz
logical::inq,match
type(BINdat)::out !matrix system

!!defaults
lmax=-1
lmin=-1
stderr=0
mrkr=4 ! length of a record marker
inq=.false.
fsize=-1
unit=13

!!process command line options


narg=iargc()
if (narg < 0)call help()

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
         ind=index(dum,',')
         read(dum(4:ind-1),*)lmax
         read(dum(ind+1:),*)lmin
      case('s')!input file size in bytes
         read(dum(4:),*)fsize
      case('i')
         inq=.true.
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is a file
      file=trim(dum)
   end if
end do

!input checks
if(fsize < 0 .and. lmax < 0)then
   write(stderr,*)"ERROR: you must use either -l=lmax,lmin or -s=FILESIZE"
   stop
end if


if(fsize >0)then ! get lmax and lmin from file size
   pcksz=(fsize-2*mrkr)/8
   !try to find out maximum degree
   lmax=0
   lmin=0
   tmpsz=1
   do while(tmpsz <= pcksz)
      lmax=lmax+1
      sz=SH_tpos(lmax,lmax,1,lmax,lmin)
      tmpsz=sz*(sz+1)/2
   end do

   match=.false.
   !now determine lmin
   do,lmin=0,lmax
      sz=SH_tpos(lmax,lmax,1,lmax,lmin)
      tmpsz=sz*(sz+1)/2
      if(tmpsz .eq. pcksz)then
         match=.true.
         exit
      end if
   end do



   if (.not. match)then
      write(stderr,*)"ERROR: no matching lmax and lmin found"
      stop
   else
      write(stderr,*)"VERBOSE: calculated lmax and lmin: ",lmax,lmin
   end if


else ! get file size from lmax and lmin
   sz=SH_tpos(lmax,lmax,1,lmax,lmin)
   pcksz=sz*(sz+1)/2
   write(stderr,*)"VERBOSE: calculated file size: ", pcksz*8+2*mrkr, "bytes"
end if


if(inq)stop !premature exit if required

!set up output matrix
out%nval1=sz
out%nval2=sz
out%pval1=pcksz
out%pval2=1
out%type='SYMVN___'
out%mtyp='P'

out%file='stdout'

out%descr="Converted covariance matrix from JKcov2_BIN"

!integer meta data
out%nint=2
allocate(out%ints(out%nint),out%ints_d(out%nint))
out%ints_d(1)='Lmax'
out%ints(1)=lmax
out%ints_d(2)='Lmin'
out%ints(2)=lmin

!double integer data
!none

! side description
allocate(out%side1_d(out%nval1))
do, m=0,lmax
   do,q=0,min(m,1)
      do,l=max(lmin,m),lmax
         ind=SH_tpos(l,m,q,lmax,lmin)
!         write(stderr,*)ind,l,m,q,out%nval1
         if(q.eq. 0)then
            write(out%side1_d(ind),'(A3,1x,2I3)')'GCN',l,m
         else
            write(out%side1_d(ind),'(A3,1x,2I3)')'GSN',l,m
         end if
      end do
   end do
end do

!load data in matrix

allocate(out%pack1(pcksz))
open(unit=unit,file=trim(file),form='unformatted')
read(unit)out%pack1
close(unit)

!write data to standard output

call write_BINtype(out)




end program JKcov_2_BIN
subroutine help()
implicit none
integer::stderr
stderr=0
write(stderr,*)"Program JKcov_2_BIN converts a binary Fortran matrix"
write(stderr,*)"(used in JK's software to the BIN format used in RR"
write(stderr,*)"Usage JKcov2_BIN [OPTIONS] FILE > OUTPUT"
write(stderr,*)"OPTIONS:"
write(stderr,*)"-s=FILESIZE: input file size in bytes"
write(stderr,*)"-l=LMAX,LMIN : set maximum and minimum degree"
write(stderr,*)"-i: inquire only don't do anything"
stop
end subroutine help

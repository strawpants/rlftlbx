!!Coded by Roelof Rietbroek, Mon Jan 20 14:58:23 2014
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!! adapted from JKcov_2_BIN
!!Updated by Roelof Rietbroek, Wed Feb  5 11:18:04 2014 Changed writing of full format




program BIN_2_JKmat
!use shtlbx
!use forttlbx
use binfiletools
use bin_operations
implicit none
integer::i,j,narg,iargc,stderr,itharg
character(200)::dum,fileout
integer::lmax,lmin,ind,q
integer::mrkr,sz,l,m,unit
integer::flmax,flmin
character(1)::trig
integer*8::fsize,pcksz,tmpsz
logical::expand,match
logical::limdeg
type(BINdat)::indat !matrix system
integer,allocatable::permvec(:)
double precision,pointer:: pack1(:)
integer::st,nd

!!defaults
lmax=-1
lmin=-1
stderr=0
mrkr=4 ! length of a record marker (NOTE: compiler and system dependent!!!)
expand=.false.
limdeg=.false.
fsize=-1
unit=13
fileout='XXXX'
indat%file='stdin'


flmax=-1 !initialize to naughty values 
flmin=-1 !initialize to naughty values 

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
         limdeg=.true.
         ind=index(dum,',')
         read(dum(4:ind-1),*)lmax
         read(dum(ind+1:),*)lmin
      case('F')!output file name
         read(dum(4:),*)fileout
      case('e')
         expand=.true.
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is a file
      indat%file=trim(dum)
   end if
end do

!input checks
if (trim(fileout) .eq. 'XXXX')then
   write(stderr,*)"ERROR: you must provide an output filename with -F="
   stop
end if


!read in meta data (side description, lmax, lmin etc)
call read_BINtype(indat,2)

!find out the maximum and minimum degree supported by the file
do, i=1,indat%nint
   select case(indat%ints_d(i))
   case('Lmax')
      flmax=indat%ints(i)
   case('Lmin')
      flmin=indat%ints(i)
   end select
end do

if( flmax < 0 .or. flmin <0 ) then
   write(stderr,*)"ERROR: no lmax,lmin found in file"
   stop
end if

if ( limdeg ) then
   if( flmax < lmax .or. lmin < flmin )then
      write(stderr,*)"ERROR: requested lmax,lmin is not supported by those of the file:",flmax,flmin
      stop
   end if

   if( flmax == lmax .and. flmin == lmin ) limdeg=.false. !no need to sort data

end if



!read remainder of the data
call read_BINtype(indat)

!expand matrices if needed
if(expand) call BIN_expand(indat)

select case(indat%type)
case('SYMVN___','SYMV1___','SYMV0___')
      write(stderr,*)" Symmetric matrix found"
case('BDSYMV0_','BDSYMVN_')
   write(stderr,*)'Block diagonal symmetric matrix found'
   if(limdeg)then
      write(stderr,*)'ERROR: to restrict lmax, and lmin the matrix needs to be expanded first'
      stop
   end if
case('BDFULLV0','BDFULLVN')
   write(stderr,*)'Block diagonal matrix found'
   if(limdeg)then
      write(stderr,*)'ERROR: to restrict lmax, and lmin the matrix needs to be expanded first'
      stop
   end if
case('FULLSQV0','FULLSQVN')
   write(stderr,*)'Full matrix found'
case('FULL2DVN')
   write(stderr,*)'Full square matrix found'
   if(indat%nval1 .ne. indat%nval2)then
      write(stderr,*)'ERROR: matrix is not square'
      stop
   end if
case default
   write(stderr,*)'ERROR:unsupported matrix type:', indat%type
   stop
end select



!restrict data output based on lmax and lmin
! a permutation vector is needed

if ( limdeg )then
   allocate(permvec(indat%nval1))
   st=0
   nd=indat%nval1
   do,i=1,indat%nval1
      read(indat%side1_d(i),'(1x,a1,2x,i3.3,i3.3)')trig,l,m

      if(l < lmin .or. lmax < l) then 
         permvec(nd)=i
         nd=nd-1
      else
         st=st+1
         permvec(st)=i
      end if

   end do
   
!permute matrix
   call BIN_permute(indat,perm1=permvec,both=.true.)
   
!restrict output matrix 
   select case(indat%type)
   case('SYMVN___','SYMV1___','SYMV0___')
      indat%nval1=st
      indat%pval1=st*(st+1)/2 !upper part of the matrix only
   
   case default
      indat%nval1=st
      indat%nval2=st
      indat%pval1=st*st

   end select
   
end if





select case( indat%mtyp )
   case( 'P' )
      write(stderr,*)"writing ",indat%pval1," 8 byte values"
      write(stderr,*)"with ",2*mrkr," bytes of fortran record markers"

      open(unit=unit,file=trim(fileout),form='unformatted')
      write(unit)indat%pack1(1:indat%pval1)
      close(unit)      
   case( 'F' )
      write(stderr,*)"writing ",indat%nval1*indat%nval1," 8 byte values"
      write(stderr,*)"with ",2*mrkr*indat%nval1*indat%nval1," bytes of fortran record markers!!!"

      open(unit=unit,file=trim(fileout),form='unformatted')
      do,i=1,indat%nval1
         do,j=1,indat%nval1
            write(unit)indat%mat1(j,i)
         end do
      end do
!!      write(unit)indat%mat1(1:indat%nval1,1:indat%nval1)
      close(unit)
   case('U') !additional matrix is needed
      allocate(pack1(indat%pval1))
      sz=0
      do,i=1,indat%nval1
         pack1(sz+1:sz+i)=indat%mat1(1:i,i)
         sz=sz+i
      end do
      write(stderr,*)"writing ",indat%pval1," 8 byte values"
      write(stderr,*)"with ",2*mrkr," bytes of fortran record markers"

      open(unit=unit,file=trim(fileout),form='unformatted')
      write(unit)pack1(1:indat%pval1)
      close(unit)      
end select

!also write the side description to standard output
do,i=1,indat%nval1
   write(*,*)indat%side1_d(i)
end do


end program BIN_2_JKmat


subroutine help()
implicit none
integer::stderr
stderr=0
write(stderr,*)"Program BIN_2_JKmat converts a binary RR matrix to binary Fortran matrix"
write(stderr,*)"(used in JK's software)"
write(stderr,*)"Usage BIN_2_JKmat [OPTIONS] BINRRFILE "
write(stderr,*)"OPTIONS:"
write(stderr,*)"-l=LMAX,LMIN : limit maximum and minimum degree of the output"
write(stderr,*)"-F=OUTPUTFILE: write the binary matrix to OUTPUTFILE"
write(stderr,*)"-e: Expand (block diagonal) matrices into full matrices"
write(stderr,*)"    This also expands symmetric matrices in full format"
write(stderr,*)"The sorting order of the side is written to standard output"
write(stderr,*)"but can be redirected to a file"
write(stderr,*)"see also JKcov_2_BIN"
stop
end subroutine help

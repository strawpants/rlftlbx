!!Program to merge  dense blocks into a blockdiagonal matrix
!!may read and write from/to standard input/output
!!Coded by Roelof Rietbroek, Thu Jan  8 16:39:45 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Fri Apr 15 14:56:30 2011
!! fixed bug: check for work%nvec > 0 before copying


program BIN_BDmerge
use binfiletools
implicit none
integer::i,j,k,blockshift,narg,itharg,stderr,maxin
parameter(maxin=300) ! maximum amount of input systems
type(BINdat)::work(maxin),out
character(200)::dum
integer::sz,szpack,vecstrt,vecnd,packstrt,packnd,nf
logical::pack
integer::iargc
!defaults/initializations
stderr=0

out%file='stdout'
nf=0

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
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option but a filename
      nf=nf+1
      if(nf > maxin)then
         write(stderr,*)"ERROR: maximum amount of inout systems supported is:",maxin
         stop
      end if
      work(nf)%file=trim(dum) !get filename
   end if
end do

!input checks

!if no file names are given try to read the amount of systems (to be read from stdin) from standard input
if( nf .eq. 0)then
   read(5,*)nf
   if(nf > maxin)then
      write(stderr,*)"ERROR: maximum amount of inout systems supported is:",maxin
      stop
   end if
   
   do, i=1,nf
      work(i)%file='stdin'
   end do
end if

if(nf <= 1)then
   write(stderr,*)"ERROR: nothing to merge"
   stop
end if





!read in meta data of first system
call read_BINtype(work(1),2)
select case(work(1)%type)
case('SYMV0___','SYMVN___','BDSYMVN_','BDSYMV0_') !matrix is allowed
   out%type='BDSYMVN_'
   pack=.true. ! symmetric matrices are in packed from
case('FULLSQVN','BDFULLVN','BDFULLV0')
   out%type='BDFULLVN'
   pack=.false.
case default
   write(stderr,*)"ERROR: matrix type is  not supported or not block diagonal:",work(1)%type
   stop
end select

!copy meta info of first system into output system
out%descr='Block diagonal matrix system from BIN_BDmerge'
out%mtyp='P'
out%ndbls=work(1)%ndbls
out%nint=work(1)%nint
if(work(1)%ndbls >0)then
   out%dbls=>work(1)%dbls
   out%dbls_d=>work(1)%dbls_d
end if
if(work(1)%nint >0)then
   out%ints_d=>work(1)%ints_d
   out%ints=>work(1)%ints
end if
out%nvec=work(1)%nvec
out%nread=work(1)%nread
if(out%nread>0)out%readme=>work(1)%readme

!amount of data
out%nval1=work(1)%nval1
out%nval2=out%nval1
out%pval1=work(1)%pval1
out%pval2=1

out%nblocks=work(1)%nblocks

!read remainder of meta data
do,i=2,nf
   call read_BINtype(work(i),2)
   select case(work(1)%type)
   case('SYMV0___','SYMVN___','BDSYMVN_','BDSYMV0_')
      if( .not. pack)then
         write(stderr,*)"ERROR: inconsistent matrix types"
         stop
      end if
   case('BDFULLVN','BDFULLV0','FULLSQVN')
      if(pack)then
         write(stderr,*)"ERROR: inconsistent matrix types"
         stop
      end if
   case default
         write(stderr,*)"ERROR: unsupported matrix type",work(i)%type
         stop
   end select

   !update total amount of data
   out%nval1=out%nval1+work(i)%nval1
   out%nval2=out%nval1
   !packed data amount
   out%pval1=out%pval1+work(i)%pval1
   !update amount of blocks
   out%nblocks=out%nblocks+work(i)%nblocks

end do



!allocate enough data for the output
allocate(out%pack1(out%pval1),out%blockind(out%nblocks))
allocate(out%side1_d(out%nval1))
if(.not. pack)allocate(out%side2_d(out%nval2)) ! allocate second side description

if(out%nvec .ne. 0)allocate(out%vec(out%nval1,out%nvec))

!check side description data and setup blockindex data
vecnd=work(1)%nval1
blockshift=work(1)%nblocks

do,i=1,vecnd
   out%side1_d(i)=work(1)%side1_d(i)
end do

if( .not. pack)then ! fill up second side description
   do,i=1,vecnd
      out%side2_d(i)=work(1)%side2_d(i)
   end do
end if

if(work(1)%nblocks .eq. 1)then
   out%blockind(1)=vecnd
else ! copy block indices from first system1
   do,i=1,work(1)%nblocks
      out%blockind(i)=work(1)%blockind(i)
   end do
end if


!loop over remaining systems
do,i=2,nf

   !check side description data
   do,k=1,work(i)%nval1
      do,j=1,vecnd
         if(out%side1_d(j) .eq. work(i)%side1_d(k))then
            write(stderr,*)"ERROR: blocks overlap"
            stop
         end if
      end do
      vecnd=vecnd+1
      !copy parameter description
      out%side1_d(vecnd)=work(i)%side1_d(k)
      if(.not. pack)out%side2_d(vecnd)=work(i)%side2_d(k)
   end do
   


   !update blockindex
   
   if(work(i)%nblocks .eq. 1)then ! inoput is a dense matrix
      out%blockind(blockshift+1)=out%blockind(blockshift)+work(i)%nval1
      blockshift=blockshift+1
   else ! input is a block diagonal matrix
      do,k=1,work(i)%nblocks
         out%blockind(blockshift+1)=out%blockind(blockshift)+work(i)%blockind(k)
         blockshift=blockshift+1
      end do
   end if
   
end do



!read in data (vector and matrix) from systems
vecstrt=0
vecnd=0
packstrt=0
packnd=0
do,i=1,nf
   vecnd=vecstrt+work(i)%nval1
   packnd=packstrt+work(i)%pval1

   !associate pointers with the correct sections of the output system

   if(work(i)%nvec >0)work(i)%vec=>out%vec(vecstrt+1:vecnd,:)
   work(i)%pack1=>out%pack1(packstrt+1:packnd)
   
   !read in data
   call read_BINtype(work(i))
   
   vecstrt=vecnd
   packstrt=packnd

end do

!write output system

call write_BINtype(out)



end program BIN_BDmerge


subroutine help()
character(5)::frmt
integer::unit
frmt='(A)'
unit=0 ! write to standard error
write(unit,frmt)" Program BIN_BDmerge merges dense diagonal matrix blocks in a block diagonal matrix"
write(unit,frmt)" usage BIN_BDsplit [OPTIONS] [FILES]"
write(unit,frmt)" Where OPTIONS may be:"
write(unit,frmt)" (no options yet)"
write(unit,frmt)" The input blocks may not have common parameters else the result is not necessarily block diagonal"
write(unit,frmt)" NOTE: input filed are read in a segmented way ( so they must be also written in segments when using pipes)"
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
stop
end subroutine help

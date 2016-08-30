!!Program to splitup a block diagonal matrix in separate dense matrix files
!!may read and write from/to standard input/output
!!Coded by Roelof Rietbroek, Thu Jan  8 16:39:45 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Fri Apr 15 14:48:02 2011
!! created better output filename ( with leading zeros, which is easier to sort)


program BIN_BDsplit
use binfiletools
implicit none
integer::i,narg,itharg,stderr,single
type(BINdat)::work,out
character(200)::dum,basename
integer::sz,szpack,vecstrt,vecnd,packstrt,packnd
logical::pack
character(8)::chars
integer::iargc
!defaults/initializations
stderr=0
work%file='stdin'
work%mtyp='P'

out%file='stdout'


basename='stdout'
single=0

!!process command line options
narg=iargc()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('f') ! get basename
         basename=trim(dum(4:))
      case('u')! read single block 
         read(dum(4:),*)single
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
      work%file=trim(dum) !get filename
   end if
end do

!input checks

!write(stderr,*)basename, work%file

!read in data

call read_BINtype(work)

select case(work%type)
case('BDSYMV0_','BDSYMVN_') !matrix is allowed
   out%type='SYMVN___'
   pack=.true. ! symmetric matrices are in packed from
case('BDFULLV0','BDFULLVN')
   out%type='FULLSQVN'
   pack=.false.
case default
   write(stderr,*)"ERROR: matrix type is  not supported or not block diagonal:",work%type
   stop
end select


!general output info
out%descr='dense matrix system from BIN_BDsplit'
out%mtyp='P'
out%ndbls=work%ndbls
out%nint=work%nint
if(work%ndbls >0)then
   out%dbls=>work%dbls
   out%dbls_d=>work%dbls_d
end if
if(work%nint >0)then
   out%ints_d=>work%ints_d
   out%ints=>work%ints
end if
out%nvec=work%nvec
out%nread=work%nread
if(out%nread>0)out%readme=>work%readme


!in the case of standard output (multiple systems) write the amount of systems to standard output
! if( single .eq. 0 .and. basename .eq. 'stdout')then
!    write(6,*)work%nblocks
! end if

sz=0
szpack=0
vecstrt=0
vecnd=0
packstrt=0
packnd=0
do,i=1,work%nblocks   !loop over blocks
   sz=work%blockind(i)-vecnd ! get size of the current block
   !   write(*,*)sz,work%blockind(i)
   if(pack)then ! retrieve size of the matrix
      szpack=sz*(sz+1)/2
   else
      szpack=sz**2
   end if
   vecnd=vecstrt+sz
   packnd=packstrt+szpack
   
   if(( i .eq. single) .or. (single .eq. 0))then
      !construct output system
      
      out%nval1=sz
      out%nval2=out%nval1
      out%pval1=szpack
      out%pval2=1
      out%side1_d=>work%side1_d(vecstrt+1:vecnd) ! point to the right side description
      if(.not. pack)out%side2_d=>work%side2_d(vecstrt+1:vecnd) ! point to the right side description

      if(out%nvec >0)out%vec=>work%vec(vecstrt+1:vecnd,:)
      out%pack1=>work%pack1(packstrt+1:packnd) ! point to right (contiguous) memory section
      
      !write block to file/standard output
      if(basename .ne. 'stdout')then
         select case(work%nblocks)
         case(1:9)
            write(chars,'(i1.1)')i
         case(10:99)
            write(chars,'(i2.2)')i
         case(100:999)
            write(chars,'(i3.3)')i
         case(1000:9999)
            write(chars,'(i4.4)')i
         end select
         out%file=trim(basename)//'_'//trim(adjustl(chars))
         
      end if
      
! write at once
      call write_BINtype(out)
      if( single >0)exit
   end if
      
   !adjust packstrt for the next iteration
   vecstrt=vecnd
   packstrt=packnd
end do





end program BIN_BDsplit


subroutine help()
character(5)::frmt
integer::unit
frmt='(A)'
unit=0 ! write to standard error
write(unit,frmt)" Program BIN_BDsplit outputs (some of) the blocks of a block diagonal matrix"
write(unit,frmt)" The output will be smaller dense matrix system"
write(unit,frmt)" usage BIN_BDsplit [OPTIONS] [FILE]"
write(unit,frmt)" Where OPTIONS may be:"
write(unit,frmt)"  -f=BASENAME: define a BASENAME for the output files"
write(unit,frmt)"    Output file names use BASENAME_BLOCKNUM (default writes all blocks to standard output)."
write(unit,frmt)"  -u=BLOCK: only extract block number BLOCK from file"
write(unit,frmt)""
write(unit,frmt)"NOTE: does not support segmented writes yet"
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
stop
end subroutine help

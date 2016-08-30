!!Coded by Roelof Rietbroek, Wed Jun  6 21:50:01 2012
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!program which converts (some) groops binary data in BIN format
!!fortran module which contains routines to read in groops binary data

module GroopsFile
use forttlbx
implicit none
type Groopsdat
character(200)::filename
character(9)::ftype,type,version
integer*8::nval1=0
integer*8::nval2=0
integer*8::uplo=0
integer::cunit=0
integer*8::pval1=0
character(8)::mattype='XXXXX'
double precision,pointer::pack1(:)
double precision,pointer::mat1(:,:)
end type Groopsdat


integer::kinddbl=8

contains

!routine to read in a groops binary file
subroutine read_GROOPS(dat)
implicit none
integer:: stderr,i
integer::opencstream
character(200)::dum
type(Groopsdat)::dat
integer*8::maxchunk,sz,st
integer::nchunk

!! defaults
stderr=0
maxchunk=1073741824 !chunked read in bytes 

!open file as a stream ( read only)
dat%cunit=opencstream(0,trim(dat%filename)//char(0))


!read the first bytes ( should start with a 'B')
call cread(dat%cunit,9,dat%ftype)
if(dat%ftype(1:1) .ne. 'B')then
   write(stderr,*)"ERROR: not a groops binary file",trim(dat%filename)
   stop
end if

call cread(dat%cunit,9,dat%type)
!write(stderr,'(A)')trim(dat%type)

!read version number
call cread(dat%cunit,9,dat%version)
!write(stderr,'(A)')trim(dat%version)

!read additional padding ( total of 64 bytes must be read now)
call cread(dat%cunit,9,dum)
!write(stderr,'(A)')trim(dum)

select case(dat%type(1:6))
case('matrix')
   !read columns and rows
   call cread(dat%cunit,8,dat%nval1)
   call cread(dat%cunit,8,dat%nval2)
 !  write(stderr,*)dat%nval1,dat%nval2

   select case(dat%mattype)
   case('UPPER')! read a symmetric matrix in packed form in chunks ( avoid 2 Gb limit of read calls)
      dat%nval1=dat%nval2
      dat%pval1=dat%nval2*(dat%nval2+1)/2
      if(.not. associated(dat%pack1))allocate(dat%pack1(dat%pval1))

     nchunk=(kinddbl*dat%pval1)/maxchunk
      do,i=1,nchunk
         st=1+((i-1)*maxchunk)/kinddbl
         sz=maxchunk
         call cread(dat%cunit,sz,dat%pack1(st))
      end do
      !read remainder

      st=1+(nchunk*maxchunk)/kinddbl
      sz=kinddbl*(dat%pval1-st+1)
!     write(0,*)'last',st,sz,nchunk
      if(sz>0)call cread(dat%cunit,sz,dat%pack1(st))
   case('GENERAL')!general square matrix
      dat%pval1=dat%nval2*dat%nval1

      if(.not. associated(dat%mat1))allocate(dat%mat1(dat%nval1,dat%nval2))

      do,i=1,dat%nval2  !loop over columns ( ensures contiguous memory is read in)
         call cread(dat%cunit,kinddbl*dat%nval1,dat%mat1(1:dat%nval1,i))
      end do

   case default
      write(stderr,*)"ERROR: matrix type is not defined:",trim(dat%mattype)
      stop
   end select
   
case default
   write(stderr,*)"ERROR: unknown data type:",dat%type
   stop
end select



!close the file
call cclose(dat%cunit)
end subroutine read_GROOPS


end module GroopsFile




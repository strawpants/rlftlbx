!subroutine to store a specific sinex block into memory and perform search restrictions on it
!on first call the block is saved in a character array
!thecomplet block is saved but the output can be resticted to search(or exclude) for strings in the block lines
subroutine get_sinexblock(file,blockname,nlines,list,inc1,excl1)
use GPStlbx,only:compstr
implicit none
!obligatory arguments
character(*),intent(in)::file,blockname
!optional arguments
integer,intent(out),optional::nlines !amount of lines (without comments)in the block
character(*),intent(out),optional,dimension(:)::list
character(*),intent(in),optional::inc1,excl1
!internal variables
integer::i,j,ndat,funit,stderr,nmax,n1,rdim1
parameter(nmax=100000)!enough
character(80),allocatable,dimension(:)::sblock
character(80)::dum
character(120)::oldblockname,oldfilename
character(12)::date1,date2,date3
logical::fex
integer,allocatable,dimension(:)::ivec1
save oldblockname,oldfilename,sblock, ndat
!defaults initializations
stderr=0
funit=13


!write(*,*)'balh' 

!check whether it is tghe first run for a specific block and file name
if(trim(blockname) .ne. trim(oldblockname) .or. (trim(file) .ne. trim(oldfilename)))then
   oldfilename=trim(file)
   oldblockname=trim(blockname)
 !check if file exists
   inquire(FILE=trim(file),EXIST=fex)
   if(.not. fex) then
      write(stderr,*) 'file '//trim(file)// ' does not exist'
      stop
   end if

!open file and check whether it is a sinex file
   
   open(unit=funit,file=trim(file),status='old')
   write(*,*)'reading block of sinex file',trim(file)
   !read first line with meta data
   read(funit,'(a80)')dum
   if (dum(1:5) .ne. '%=SNX')stop 'does not appear to be a sinex file'

!now calculate amount of data lines in the block

   do,i=1,nmax !loop tio search for the block
      read(funit,'(A80)')dum
      if(index(dum,'+'//trim(blockname)) .ne. 0)then
         ndat=0
         do,j=1,nmax
            read(funit,'(A80)')dum
            if(dum(1:1).eq.'*')cycle !ignore comments
            if(index(dum,'-'//trim(blockname)) .ne. 0)exit
            ndat=ndat+1 !increment
         end do
         exit !exit outer loop
      else if (index(dum,'%ENDSNX') .ne. 0)then!block is not found
         write(stderr,*)'Sinex block ',trim(blockname),' is not found'
         stop
      end if
      
   end do
   !rewind file
   rewind(funit)
   !now we can allocate space for the data array
   if(allocated(sblock))deallocate(sblock)
   allocate(sblock(ndat))
   sblock=''
   
   !now read in data
   do,i=1,nmax !loop tio search for the block
      read(funit,'(A80)')dum
      if(index(dum,'+'//trim(blockname)) .ne. 0)then
         ndat=0
         do,j=1,nmax
            read(funit,'(A80)')dum
            if(dum(1:1).eq.'*')cycle !ignore comments

            if(index(dum,'-'//trim(blockname)) .ne. 0)exit
            ndat=ndat+1 !increment
            sblock(ndat)=dum
         end do
         exit !exit outer loop
      end if
      
   end do
   
   
close(funit)
!write(*,*)'test',inc1
end if !and of intialization part 

allocate(ivec1(ndat))


if(present(inc1) .and. present(excl1))then
   rdim1=3
else if (present(inc1))then
   rdim1=2
else if (present(excl1))then
   rdim1=1
else
   rdim1=0
end if

!now get the number of variables satisfying the restrictions 

!first dimension
select case(rdim1)
case(0)
   n1=ndat
  
   ivec1=(/ (i, i = 1,ndat) /)
  
case(3)!both restrictions
   n1=0
   do,i=1,ndat
      if(.not.compstr(sblock(i),inc1))cycle !cycle when inc1 string is NOT found
      if( compstr(sblock(i),excl1))cycle !cycle when excl1 string is found
      n1=n1+1
      ivec1(n1)=i !set index string
   end do
case(2)
   n1=0
   do,i=1,ndat
      if(.not.compstr(sblock(i),inc1))cycle !cycle when inc1 string is NOT found
      n1=n1+1
      ivec1(n1)=i !set index string
!      write(*,*)list(i)
   end do
case(1)
   n1=0
   do,i=1,ndat
      if( compstr(sblock(i),excl1))cycle !cycle when excl1 string is found
      n1=n1+1
      ivec1(n1)=i !set index string
   end do
end select
   
!now put data in to arguments

if(present(nlines))nlines=n1


if(present(list))then
   if(size(list,1)<n1)then
      write(stderr,*) 'vector list has dimension ',size(list,1),' required ',n1
      stop
   else
      if(len(list(1))<80)stop'list characters width must be at least 80'
    !  write(*,*)sblock(ivec1(1:n1))
      list(1:n1)=sblock(ivec1(1:n1))
   end if
end if



end subroutine get_sinexblock

!!THIS assembly of subroutines should replace those found in BINfiletools3.f90
!!It enables a more modular writing and retrieval of matrices and arrays stored in binary format
!!Coded by Roelof Rietbroek, Wed Feb  6 11:26:07 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!The systems which may be stored in thosse files are general matrix vector systems

!!several types of systems are considered:


!!symmetric matrices

!!block diagonal matrices

!!full matrices 

!!SVD decompositions(planned)



!!A BINARY file is stored as a continuous stream of bytes  (this enables linux piping)
!!The data is stored in 4 segments 

!!FIRST SEGMENT
!!the first segment contains index data in the following order:
!! 2 byte integer: endian test parameter translates on a big endian system as 'IB' if interpreted as ascii and 'BI' on a little endian system
!! 6 bytes character: (remaining part of) the version (8 byte character) the first two characters are derived from the endian parameter
!! 8 byte character: Matrix type:
!!          currently the following types exist: 'SYMV0___','SYMV1___'(old),SYMV2___'(old), 'SYMVN___',BDSYMV0_','BDFULLV0','FULLSQV0'
!!           which contain respectively: symmetric matrix in packed storage for 0 1 2 or an arbitrary amount of vectors, BLock diagonal symmetric matrix with no vectors, Block diagonal matrix with no vectors, Full square matrix
!! 80 byte character: short description of the content
!! 6 * kindint byte integers: nint (amount of integer meta data),ndbls (amount of double meta data),nval1 (full matrix size side 1),nval2(full matrix size side 2),pval1 (stored matrix size side 1),pval2(stored matrix size side 2)
!!2 * kindint byte integers (only in versions > 2.1) nvec(number of stored vectors),number of stored lines in the readme section

!!Type specific index data
!! kindint byte integer: nblocks (amount of diagonal blocks) only for type: 'BDSYMV0_','BDFULLV0'

!!SECOND SEGMENT
!!Contains the meta data
!!nread *80 byte character array readme (only for version > 2.1)
!! nint*24 bytes character array ints_d: description of the integer meta data
!! nint*kindint bytes integer array ints: integer metad data
!! ndbls*24 bytes character array dbls_d: description of the double meta data
!! nint*kinddbl bytes double array dbls: double meta data
!! nval1*24 bytes character array with the description of the side of the array
!!type specific meta data:
!! nblocks*kindint integer array with locations of the end of the blocks in the packed array ( only for types:'BDSYMV0_','BDFULLV0')


!!THIRD SEGMENT
!!vector data
!! kinddbls*nval1 *nvec byte double array: vec(nval1,nvec):

!!FOURTH SEGMENT
!! matrix data
!! Depending on the type a full or packed matrix is written:
!! pval1*pval2* kinddbl byte double array: pack1 or mat1


!USE:
!! use the derived type below to set the way the data is read out or written


!!

!!Updated by Roelof Rietbroek, Mon Jul 21 11:32:50 2008
!! rewrite of the routines. the routines now read and write a derived type array with all the necessary info
!! versions can now hopefully be easier to maintain
!! the module necessary for the type declaration is also included in this file
!! For portability and better stream behavior all the read /write/ open and close routines are now done trhouygh wrapped C routines
!! added a readme text component (80 character width any length) for better self documentation
!! Implemented an endian check and supply byte swapping if necessary
!! The endianness of the file is that which is native on the system, byte swapping is only performed on reads when different endianness is detected
!! the new routines don't use save parameters anymore, all neccessary data is stored in the derived type BINdat
!! the use of the derived type as input makes the routines much shorter

!! module declaring C routines, derived Type and read and write routines 

!!Updated by Roelof Rietbroek, Mon Aug  4 10:33:09 2008
!! routines may now also read and write from standard input/standard output
!! the filename should then be declared 'stdin' and 'stdout' respectively 
!!Updated by Roelof Rietbroek, Thu Feb 12 14:55:40 2009
!! updated sortbintype to use an inplcae symmetric permutation on symmteric matrices

!!Updated by Roelof Rietbroek, Mon May 11 12:11:18 2009
!! added second side description for BDFULL and FULLSQ type ( added compatibility for older version)
!! increased version number
!!Updated by Roelof Rietbroek, Wed Jun 10 17:01:32 2009
!!fixed some bugs with endian swapping issues (appears to work now on a SUN and linux machine, big vs little endian)
!!Updated by Roelof Rietbroek, Wed Aug  5 16:43:28 2009 Modified for using monster files
!! use chunks for reading and writing packed matrices ( may cause trouble for large arrays)
!! use large integers for the packed sizes
!!Updated by Roelof Rietbroek, Fri Aug 14 15:56:59 2009
!! added maxblocksz parameter to speed up packindex routine ( this parameter is derived upon read and will not be written to a file)

!!Updated by Roelof Rietbroek, Tue Nov 17 17:30:40 2009
!!added optional valid return with error

!!Updated by Roelof Rietbroek, Thu Aug  5 10:50:51 2010
!!changed integer meta data to integer*8 to cope with huge number of observations ( updated version number to 2.5)

!!Updated by Roelof Rietbroek, Tue Nov  2 12:25:14 2010
!!added a file integrity check at the beginning of the read routine

module BINfiletools
implicit none


!!! the derived type containg all info (and more) of the files BINdat contains pointers to the actual memory intensive data sections
!!this allows rthe pointers to be allocated within the read subroutines
type BINdat 
!fixed size components
character(200)::file='' !contains file name (MUST be set on write and read)
character(8)::ver,type !version and type strings (type MUST be set on write)
character(80)::descr !short description strings
character(1)::mtyp='P' !'L' or 'U' for reading and writing symmetric matrices if 'L' then mat1 's lower triangle is the valid part
!!'F' when a full matrix write is required or 'P' when a packed matrix write is required
!!this parameter must be set on write/read (it is intialized to read out the packed matrix which will always work)

integer::progress=0 ! progress parameter indicating the progress of read and writes (necessary for segmented read and writes)
!!don't mess with this parameter yourself

integer::cunit=0 ! parameter which keeps track of the C file descriptor (and whether it is opened or not (0))!!don't mess with this parameter yourself

!integer sizes of the matrices and other meta data
integer::nint=0
integer::ndbls=0
integer::nval1=0
integer::nval2=0
integer*8::pval1=0
integer*8::pval2=0
integer*8::ldm=0 ! stores leading dimension of the matrix
integer::nblocks=1 ! default assumes one full block
integer::maxblocksz=0 ! parameter which accelerates indexing function
integer::nvec=0
integer::nread=0 

!!The parameters MUST be set on write (WARNING no check is made for consistency with the actual data provided)
logical::swap=.false. !only use full for read routine ( keeps hold on whether data should be byte swapped on read)(! don't change your self) )
real::rlver=0 ! version number as real ( don't change yourself)
!pointers may be (de)allocated at a later stage (all pointers all nullified on initialization)

character(80),pointer::readme(:)=>null() !readme array with additional info
character(24),dimension(:),pointer::ints_d=>null() !description of integer meta data
character(24),dimension(:),pointer::dbls_d=>null()!description of double meta data
character(24),dimension(:),pointer::side1_d=>null() !description side 1 parameters
character(24),dimension(:),pointer::side2_d=>null() !description side 2 parameters
integer,pointer::blockind(:)=>null() ! for blockdiagonal matrices 
integer*8,pointer::ints(:)=>null() !integer metadata ( changed to integer*8 on 05-08-2010)
double precision,pointer::dbls(:)=>null() !double metadata
double precision,pointer,dimension(:)::pack1=>null() !matrix in packed form
double precision,pointer,dimension(:,:)::vec=>null() !vector data
double precision,pointer,dimension(:,:)::mat1=>null() !full matrix data (or symmetric matrix (only one triangle is used)


!!on read the data arrays (here pointers) may be already allocated and are allowed to be larger than the actual data size
!! if they are unassociated they will be tightly allocated according to the dimensions provided 

end type BINdat

!endian parameter which reads as 'BI' on a little endian machine and 'IB' on a big endian machine
integer*2,parameter::endian=ichar('B')+ichar('I')*256
integer,parameter::kindint=4 ! amount of bytes in an default integer
integer,parameter::kinddbl=8 !amount of bytes in an default double
character(8)::binvers='BINV2.5'  !NOTE the first two letters must always be 'BI', this is used for checking the endianess of the file
double precision::rlbinvers=2.5 ! must be consistent with above
contains

  
  
  !Fortran wrapper for the C function to open a file (returns an index cunit to identify the file descriptor)
  function copen(perm,file)
    integer::copen,opencstream,perm
    character(*),intent(in)::file

    copen=opencstream(perm,trim(file)//char(0))
    
  end function copen




!!!!!!!WRITE ROUTINE!!!!!!!!!!!!!!

Subroutine write_BINtype(dat,stage,err) 
  implicit none
  !!general input
  type(BINdat)::dat
  integer,optional::stage !desired stage to return control to calling unit (for segmented writes)
  integer,optional::err
  
  !!private arguments
  integer::stderr,i,j,istage
  parameter(stderr=0)

  integer*8::st,maxchunk
  integer::nchunk,sz
  logical::verbose
  
  if(present(err))err=0 !no error

  !defaults initializations
  maxchunk=1073741824
  !maxchunk=619760
  verbose=.false.

  dat%ver=binvers !NOTE the first two letters must always be 'BI', this is used for checking the endianess of the file
  
  if(present(stage))then
     istage=stage
  else
     istage=200 !large (everything will be written)
  end if
  
  
  !issue an ERROR when the program attempts to write the file in the wrong order
  if(dat%progress >= istage)then
     write(stderr,*)"ERROR:write_BINtype attempting to write segments in the wrong order"
     if(present(err))then
        err=-1
        return
     else
        stop
     end if
  end if
  
!Check whether the type is defined

select case(dat%type)
case('SYMV0___','SYMVN___','BDSYMV0_','BDSYMVN_','BDFULLV0','BDFULLVN','FULLSQV0','FULLSQVN','FULL2DVN','DIAVN___')
! do nothing
case default
   write(stderr,*)'ERROR write_BINtype: Unknown type selected:',dat%type
   if(present(err))then
        err=-2
        return
     else
        stop
     end if
end select


if(dat%progress .eq. 0)then ! write first segment (also opens file)

   !open file through C routines for writing
   if(trim(dat%file) .eq. 'stdout')then 
      dat%cunit=copen(4,trim(dat%file))
   else
      dat%cunit=copen(1,trim(dat%file))
   end if
   if(verbose)write(stderr,*)"VERBOSE: writing index data"
   !write first data
   call cwrite(dat%cunit,2,endian)!endian parameter
   call cwrite(dat%cunit,6,dat%ver(3:8))!rest of the version name
   call cwrite(dat%cunit,8,dat%type)!type
   call cwrite(dat%cunit,80,dat%descr) !short description
   call cwrite(dat%cunit,kindint,dat%nint)
   call cwrite(dat%cunit,kindint,dat%ndbls)
   call cwrite(dat%cunit,kindint,dat%nval1)
   call cwrite(dat%cunit,kindint,dat%nval2)
   !large integers
   call cwrite(dat%cunit,kindint*2,dat%pval1)
   call cwrite(dat%cunit,kindint*2,dat%pval2)

   call cwrite(dat%cunit,kindint,dat%nvec) !new in this version (version > 2.1)
   call cwrite(dat%cunit,kindint,dat%nread)!new in this version (version > 2.1)



   !type specific records
   select case(dat%type)
   case('BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN')!block diagonal specific data
      if(verbose)write(stderr,*)"VERBOSE: writing type specific index data"
      call cwrite(dat%cunit,kindint,dat%nblocks)
   end select
   
   dat%progress=1 ! reset progress parameter
   
   if(istage .eq. 1)return !leave routine prematurely in case of a segmented write
end if !first stage



if(dat%progress .eq. 1)then !write second segment
   
   if(verbose)write(stderr,*)"VERBOSE: writing info part (readme/integers/doubles)"
   !write readme array to file
   if(dat%nread > 0)then
      call cwrite(dat%cunit,80*dat%nread,dat%readme)
   end if


   if(dat%nint > 0)then
   !      write(*,*)'writing integer data'
      call cwrite(dat%cunit,24*dat%nint,dat%ints_d)
      call cwrite(dat%cunit,2*kindint*dat%nint,dat%ints) ! changed to integer*8 from version 2.5
   end if
   
   !write double meta data to record

   if(dat%ndbls > 0)then
   !      write(*,*)'writing integer data'
      call cwrite(dat%cunit,24*dat%ndbls,dat%dbls_d)
      call cwrite(dat%cunit,kinddbl*dat%ndbls,dat%dbls)
   end if


   !write matrix side description to data
   !first side 
   if(verbose)write(stderr,*)"VERBOSE: writing side description of rows"
   call cwrite(dat%cunit,24*dat%nval1,dat%side1_d)

   
   !!write type specific meta data to file
   !block diagonal indices
   select case(dat%type)
   case('BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN')!block diagonal specific data
      if(verbose)write(stderr,*)"VERBOSE: writing index vector of diagonal blocks"
      call cwrite(dat%cunit,kindint*dat%nblocks,dat%blockind)
   end select

!second side description
   select case(dat%type)
   case('FULL2DVN','BDFULLVN','BDFULLV0','FULLSQVN','FULLSQV0') ! write second side description
      !note 'BDFULLV0' does not have a second side description!!
      ! write second side description
         call cwrite(dat%cunit,24*dat%nval2,dat%side2_d)
   end select
   
   dat%progress=2
   if(istage .eq. 2)return !but proceed in a contiguous write
end if


!write data vectors to file
if(dat%progress .eq. 2)then
   if(verbose)write(stderr,*)"VERBOSE: writing data vectors"
   if(dat%nvec > 0)then
      !write vectors (per vector in order to garantee a contiguous part)
      do,i=1,dat%nvec
         call cwrite(dat%cunit,kinddbl*dat%nval1,dat%vec(1:dat%nval1,i))
      end do
   end if


   dat%progress=3
   if(istage .eq. 3)return
end if




!!write matrices to file
if(dat%progress .eq. 3)then
   if(verbose)write(stderr,*)"VERBOSE: writing data matrix"
   select case(dat%mtyp)
   case('L')! get symmmetric/triangular packed matrix from lower triangle of the full matrix
      do,i=1,dat%nval1
         do,j=1,i !inner loop required for write since dat%mat1(i,1:i) does not refer to a contigious set of data
            !it actually uses a stride argument which cannot be recognized in C !
            !consequently the write of a lower triangular matrix is expected to be somewhat slower than a upper triangular matrix
            call cwrite(dat%cunit,kinddbl,dat%mat1(i,j))
         end do
      end do
   case('U')!get  get symmmetric/triangular packed matrix from upper triangle of the full matrix
      
      do,i=1,dat%nval1 ! this will work since dat%mat1(1:i,i) refers to a contigious area of data

         call cwrite(dat%cunit,kinddbl*i,dat%mat1(1:i,i))
      end do
!    case('B')!write upper triangle of the block diagonal full matrix
!       !write first block from upper triangle
!       shift=0
!       sz=dat%blockind(1)
!       do,i=1,sz
!          call cwrite(dat%cunit,kinddbl*i,dat%mat1(shift+1:shift+i,1))
!          shift=shift+sz ! proceed to the next column
!       end do

!       !loop over remaining blocks
!       do,j=2,dat%nblocks
!          sz=dat%blockind(j)-dat%blockind(j-1) !size of the block
!          do,i=1,sz
!             call cwrite(dat%cunit,kinddbl*i,dat%mat1(shift+1:shift+i,1))
!             shift=shift+sz ! proceed to the next column
!          end do
!       end do

   case('P')!write packed matrix in chunks (new from version 2.4, for monster writes)
      nchunk=(kinddbl*dat%pval1*dat%pval2)/maxchunk
      do,i=1,nchunk
         st=1+((i-1)*maxchunk)/kinddbl
         sz=maxchunk
!         write(0,*)'chunk',i,st,sz
         call cwrite(dat%cunit,sz,dat%pack1(st))
      end do
      !write remainder

      st=1+(nchunk*maxchunk)/kinddbl
      sz=kinddbl*(dat%pval1*dat%pval2-st+1)
!     write(0,*)'last',st,sz,nchunk
      if(sz>0)call cwrite(dat%cunit,sz,dat%pack1(st))

   case('F')!write full matrix directly
      do,i=1,dat%nval2  !loop over columns ( ensures contiguous memory is read in)
         call cwrite(dat%cunit,kinddbl*dat%nval1,dat%mat1(1:dat%nval1,i))
      end do
   end select
end if


!reset pipe progress and close unit


if(trim(dat%file) .ne. 'stdout')call cclose(dat%cunit)

dat%progress=0 ! succesfull closure

end Subroutine write_BINtype



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SUBROUTINE to read binary matrix(from file or linux fifo pipe) into arrays!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Subroutine read_BINtype(dat,stage,err)
implicit none
!!general in/output
type(BINdat)::dat
integer,optional::stage
integer,optional::err

!!private arguments
integer*2::endiantest
integer::stderr,i,j,istage
parameter(stderr=0)
logical::fex


integer*8::st,maxchunk
integer::nchunk,sz
integer,allocatable::itmp(:)  
  
  

!defaults initializations
if(present(err))err=0
istage=200 !read untill the end
maxchunk=1073741824
if(present(stage))istage=stage




!issue an ERROR when the program attempts to read the file in the wrong order
if(dat%progress >= istage)then
   write(stderr,*)"ERROR: read_BINtype: attempting to read segments in the wrong order"
   if(present(err))then
        err=-1
        return
     else
        stop
     end if

end if


if(dat%progress .eq. 0)then !read first stage

   !open file for read
   if(trim(dat%file) .eq. 'stdin')then 
!      write(*,*)'using standard input'
      dat%cunit=copen(3,trim(dat%file))
   else
      !inquire whether the file exists
      inquire(FILE=trim(dat%file),EXIST=fex)
      if(.not. fex)then
         write(stderr,*)'ERROR: read_BINtype: file does not exist: ',trim(dat%file)
         if(present(err))then
            err=-2
            return
         else
            stop
         end if
   end if
      dat%cunit=copen(0,trim(dat%file))
   end if

   
   !read first data
   call cread(dat%cunit,2,endiantest)!endian parameter

   if(endian .eq. endiantest)then
      dat%swap=.false.
   else
      dat%swap=.true. ! file is of opposite endiannes
   end if

   call cread(dat%cunit,6,dat%ver(3:8))!rest of the version name
   dat%ver(1:2)='BI' !replace first two letters

   if(dat%ver(1:4) .ne. 'BINV')then ! give error
      write(stderr,*)'ERROR: read_BINtype: wrong format?'
         if(present(err))then
            err=-2
            return
         else
            stop
         end if
   end if

   read(dat%ver(5:),*)dat%rlver ! read version number in real parameter

   !give error when the version is newer than the version supported by the compiled BINfiletools
   if(dat%rlver > rlbinvers)then
      write(stderr,*)'ERROR: read_BINtype: file has newer version than supported: ',dat%rlver
         if(present(err))then
            err=-2
            return
         else
            stop
         end if
   end if
   
   call cread(dat%cunit,8,dat%type)!type
   call cread(dat%cunit,80,dat%descr) !short description
   call cread(dat%cunit,kindint,dat%nint)
   call cread(dat%cunit,kindint,dat%ndbls)
   call cread(dat%cunit,kindint,dat%nval1)
   call cread(dat%cunit,kindint,dat%nval2)

   if(dat%swap)then
      call cswap(kindint,1,dat%nint)
      call cswap(kindint,1,dat%ndbls)
      call cswap(kindint,1,dat%nval1)
      call cswap(kindint,1,dat%nval2)
   end if

   !compatibility clause ( use of large integers for packed sizes from version 2.4)
   if(dat%rlver<=2.3)then
      call cread(dat%cunit,kindint,i)
      call cread(dat%cunit,kindint,j)
      if(dat%swap)then
         call cswap(kindint,1,i)
         call cswap(kindint,1,j)
      end if
      dat%pval1=i
      dat%pval2=j
   else
      call cread(dat%cunit,kindint*2,dat%pval1)
      call cread(dat%cunit,kindint*2,dat%pval2)
      if(dat%swap)then
         call cswap(kindint*2,1,dat%pval1)
         call cswap(kindint*2,1,dat%pval2)
      end if
   end if



!Compatibility for older versions
   
   if(dat%rlver <= 2.1)then
      dat%nread=0 ! no readme part yet
      
      !get amount of vectors from type description
      select case(dat%type)
      case('SYMV0___','BDSYMV0_','BDFULLV0','BDFULLVN')
         dat%nvec=0
         dat%nval2=dat%nval1
         dat%pval2=1
      case('SYMV1___')
         dat%nvec=1
         dat%nval2=dat%nval1
         dat%pval2=1
      case('SYMV2___')
         dat%nvec=2
         dat%nval2=dat%nval1
         dat%pval2=1
      case('FULLSQV0')
         dat%nvec=0
         dat%nval2=dat%nval1
         dat%pval2=dat%pval1
      end select
      
      !!thus one can now read out any matrix in a vector if requested. how it is read out is controlled by the mtyp parameter
   else
      call cread(dat%cunit,kindint,dat%nvec) !new in this version (version > 2.1)
      call cread(dat%cunit,kindint,dat%nread)!new in this version (version > 2.1)

      if(dat%swap)then
      call cswap(kindint,1,dat%nvec)
      call cswap(kindint,1,dat%nread)
      end if
   end if






   !type specific records
   select case(dat%type)
   case('BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN')!block diagonal specific data
      call cread(dat%cunit,kindint,dat%nblocks)
      if(dat%swap)call cswap(kindint,1,dat%nblocks)
   end select

   dat%progress=1 ! reset progress parameter
   
   if(istage .eq. 1)return !leave routine prematurely in case of a segmented write
end if !first stage


!second stage

if(dat%progress .eq. 1)then !read second segment
   
   !read readme array to file
   if(dat%nread > 0)then
      if(.not. associated(dat%readme))allocate(dat%readme(dat%nread)) ! only allocate when not associated
      call cread(dat%cunit,80*dat%nread,dat%readme)
   end if

   if(dat%nint > 0)then

      if(dat%rlver <= 2.4)then ! compatibility clause ( read normal integers instead of 8byte ones)
         if( .not. associated(dat%ints))allocate(dat%ints(dat%nint)) ! only allocate when not associated
         if( .not. associated(dat%ints_d))allocate(dat%ints_d(dat%nint)) ! only allocate when not associated
         allocate(itmp(dat%nint)) ! temporary integer vector old style
         call cread(dat%cunit,24*dat%nint,dat%ints_d)
         call cread(dat%cunit,kindint*dat%nint,itmp)
         !swap integers if necessary
         if(dat%swap)call cswap(kindint,dat%nint,itmp)
         !copy values in integer*8 version
         forall(i=1:dat%nint)dat%ints(i)=itmp(i)
         deallocate(itmp)
      else
         if( .not. associated(dat%ints))allocate(dat%ints(dat%nint)) ! only allocate when not associated
         if( .not. associated(dat%ints_d))allocate(dat%ints_d(dat%nint)) ! only allocate when not associated
         call cread(dat%cunit,24*dat%nint,dat%ints_d)
         call cread(dat%cunit,2*kindint*dat%nint,dat%ints)
         !swap integers if necessary
         if(dat%swap)call cswap(2*kindint,dat%nint,dat%ints)
      end if
   end if
   
   !read double meta data to record

   if(dat%ndbls > 0)then
      if( .not. associated(dat%dbls))allocate(dat%dbls(dat%ndbls)) ! only allocate when not associated
      if( .not. associated(dat%dbls_d))allocate(dat%dbls_d(dat%ndbls)) ! only allocate when not associated
      call cread(dat%cunit,24*dat%ndbls,dat%dbls_d)
      call cread(dat%cunit,kinddbl*dat%ndbls,dat%dbls)

      !swap doubles 
      if(dat%swap)call cswap(kinddbl,dat%ndbls,dat%dbls)

   end if


   !read matrix side description to data
   !first side 
   if(.not. associated(dat%side1_d))allocate(dat%side1_d(dat%nval1))
   
   call cread(dat%cunit,24*dat%nval1,dat%side1_d)


   !!read type specific meta data to file
   select case(dat%type)
   case('BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN')!block diagonal specific data
      if(.not. associated(dat%blockind))allocate(dat%blockind(dat%nblocks))
      call cread(dat%cunit,kindint*dat%nblocks,dat%blockind)
      if(dat%swap)then
         call cswap(kindint,dat%nblocks,dat%blockind)
      end if
      !calculate maximumsize of the block
      dat%maxblocksz=max(dat%blockind(1),&
           maxval(dat%blockind(2:dat%nblocks)-dat%blockind(1:dat%nblocks-1)))

   end select
   
   !second side description
   select case(dat%type)
   case('BDFULLVN','BDFULLV0','FULLSQVN','FULLSQV0') ! read second side description
      
      if(.not. associated(dat%side2_d))allocate(dat%side2_d(dat%nval2))
      
      ! read second side description if version number >2.2 else copy values from side1
      if(dat%rlver<=2.2)then ! old version
         forall(i=1:dat%nval1)dat%side2_d(i)=dat%side1_d(i)
      else ! newer version ( may have  a different side2 description)
         call cread(dat%cunit,24*dat%nval2,dat%side2_d)
      end if
   case('FULL2DVN')

      if(.not. associated(dat%side2_d))allocate(dat%side2_d(dat%nval2))
      call cread(dat%cunit,24*dat%nval2,dat%side2_d)
   case default !let the second side point to the first side
      dat%side2_d=>dat%side1_d(1:dat%nval1)
   end select
   
   dat%progress=2
   if(istage .eq. 2)return !but proceed in a contiguous read
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!READ VECTOR and MATRIX DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!read vector data 
if(dat%progress .eq. 2)then
   if(dat%nvec > 0)then
      if(.not. associated(dat%vec))allocate(dat%vec(dat%nval1,dat%nvec))
      
      do,i=1,dat%nvec
         call cread(dat%cunit,kinddbl*dat%nval1,dat%vec(1:dat%nval1,i))
         !swap data if required
         if(dat%swap)call cswap(kinddbl,dat%nval1,dat%vec(1:dat%nval1,i))
      end do
   end if

   dat%progress=3
   if(istage .eq. 3)return
end if

!read in third segment data and vectors
if(dat%progress .eq. 3)then
   select case(dat%mtyp)
   case('L')! put symmmetric/triangular packed matrix in lower triangle of the full matrix
      if( .not. associated(dat%mat1))allocate(dat%mat1(dat%nval1,dat%nval2))
      dat%ldm=size(dat%mat1,1)
      do,i=1,dat%nval1
         do,j=1,i !inner loop required for read since dat%mat1(i,1:i) does not refer to a contigious set of data
            call cread(dat%cunit,kinddbl,dat%mat1(i,j))
         end do
      end do
      if(dat%swap)then
         do,i=1,dat%nval2 ! loop over columns (representing contiguous data blocks)
            call cswap(kinddbl,dat%nval2-i+1,dat%mat1(i:,i))
         end do
      end if

   case('U')!put symmmetric/triangular packed matrix from upper triangle of the full matrix
      if( .not. associated(dat%mat1))allocate(dat%mat1(dat%nval1,dat%nval2))
      dat%ldm=size(dat%mat1,1)
      do,i=1,dat%nval2 ! this will work since dat%mat1(1:i,i) refers to a contigious area of data
         call cread(dat%cunit,kinddbl*i,dat%mat1(1:i,i))
      end do
      if(dat%swap)then
         do,i=1,dat%nval2 ! loop over columns (representing contiguous data blocks)
            call cswap(kinddbl,i,dat%mat1(1:i,i))
         end do
      end if
!    case('B')!expand symmetric block diagonal matrix in upper triangles of a full block diagonal matrix
!       if(.not. associated(dat%mat1))then ! allocate enough memory when not done already
!          ind=dat%blockind(1)**2 ! first block size
!          do,i=2,dat%nblocks ! loop over remaining blocks
!             ind=ind+(dat%blockind(i)-dat%blockind(i-1))**2 ! update blcok size
!          end do
!          allocate(dat%mat1(ind,1))
!       end if
      
!       !read first block in upper triangle
!       shift=0
!       sz=dat%blockind(1)
!       do,i=1,sz
!          call cread(dat%cunit,kinddbl*i,dat%mat1(shift+1:shift+i,1))
!          shift=shift+sz ! proceed to the next column
!       end do

!       !loop over remaining blocks
!       do,j=2,dat%nblocks
!          sz=dat%blockind(j)-dat%blockind(j-1) !size of the block
!          do,i=1,sz
!             call cread(dat%cunit,kinddbl*i,dat%mat1(shift+1:shift+i,1))
!             shift=shift+sz ! proceed to the next column
!          end do
         
!       end do

!       !swap data if necessary
!       if(dat%swap)then
!          shift=0
!          sz=dat%blockind(1)
!          do,i=1,sz
!             call cswap(kinddbl,i,dat%mat1(shift+1,1)))
!             shift=shift+sz ! proceed to the next column
!          end do

!          !loop over remaining blocks
!          do,j=2,dat%nblocks
!             sz=dat%blockind(j)-dat%blockind(j-1) !size of the block
!             do,i=1,sz
!                call cswap(kinddbl,i)dat%mat1(shift+1,1))
!                shift=shift+sz ! proceed to the next column
!             end do
         
!          end do
         
!       end if

   case('P')!read packed matrix directly in 1D vector
      if( .not. associated(dat%pack1))allocate(dat%pack1(dat%pval1*dat%pval2))
      dat%ldm=size(dat%pack1,1)
      nchunk=(kinddbl*dat%pval1*dat%pval2)/maxchunk
      do,i=1,nchunk
         st=1+((i-1)*maxchunk)/kinddbl
         sz=maxchunk
         call cread(dat%cunit,sz,dat%pack1(st))
         if(dat%swap)call cswap(kinddbl,sz,dat%pack1(st))
      end do
      !write remainder
      st=1+(nchunk*maxchunk)/kinddbl
      sz=kinddbl*(dat%pval1*dat%pval2-st+1)
      if(sz>0) call cread(dat%cunit,sz,dat%pack1(st))
      if(dat%swap)call cswap(kinddbl,sz/kinddbl,dat%pack1(st))

   case('F')!read full matrix directly in 2D matrix
      if( .not. associated(dat%mat1))allocate(dat%mat1(dat%nval1,dat%nval2))
      dat%ldm=size(dat%mat1,1)
      do,i=1,dat%nval2  !loop over columns ( ensures contiguous memory is read in)
         call cread(dat%cunit,kinddbl*dat%nval1,dat%mat1(1:dat%nval1,i))
      end do
      if(dat%swap)then
         do,i=1,dat%nval2  !loop over columns ( ensures contiguous memory is swapped)
            call cswap(kinddbl,dat%nval1,dat%mat1(1:dat%nval1,i))
         end do
      end if
   end select
end if


if(trim(dat%file) .ne. 'stdin')call cclose(dat%cunit)
dat%progress=0


end Subroutine read_BINtype


!subroutine to sort left and right rows and columns according to index vector and possibly transpose the matrix
!The index vector contains the original indices and the order determines the new scheme
!Zero entries which are ignored. The subroutine compacts the vector (removes zero rows,columns)
!!Updated by Roelof Rietbroek, Thu Feb 12 14:56:32 2009
!!use symmetric permutation
!!Updated by Roelof Rietbroek, Fri Dec 17 15:11:48 2010 (removed bug in using compact where compact2 was needed)
!!Updated by Roelof Rietbroek, Fri Apr 15 14:41:58 2011
!! removed bug: added check for input%nvec > 0 before permuting the vector


subroutine sort_BINtype(input,side1,side2,both)
use forttlbx
implicit none
type(BINdat)::input
integer,dimension(:),optional,target::side1,side2 ! optional sorting indices
logical,optional,intent(in)::both
integer,pointer,dimension(:)::compact=>null() ! may be either allocated by itself or point to a part of side1
integer,pointer,dimension(:)::compact2=>null()
integer::nrows,ncols,stderr,i,j,row,col
double precision,allocatable::tmpvec(:)
double precision::time1,time2
logical::iboth

!defaults
nrows=input%nval1
ncols=input%nval2

stderr=0
iboth=.false.

! if(present(side2) .and. present(side1))then
!    iboth=.true.
! else if(present(both))then
!    iboth=both
! end if


if(present(both))then
   iboth=both ! permute both sides with the same vector
end if
if( iboth .and. .not. present(side1))then
   write(stderr,*)"ERROR: sort_BINtype: symmetric permutation required but no permutation vector given"
      stop
end if


if( iboth .and. present(side1) .and. present(side2))then ! impossible combo
   write(stderr,*)"ERROR: sort_BINtype conflicting arguments"
   write(stderr,*)"       permutation vectors differ for sides"
   write(stderr,*)"       but symmetric permutation is requested"
      stop

end if



!call cpu_time(time1)

!some preliminary input checks




!check matrix system
select case(input%type)
case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')! symmetric
   !copy values in case of one sided permutation or 
   if( .not. iboth)then 
      write(stderr,*)"ERROR: sort_BINtype Symmetric matrices may only have symmetric permutation"
      stop
   end if

   select case(input%mtyp)
   case('U') ! mirror
      do,i=1,input%nval2
         do,j=1,i-1
            input%mat1(i,j)=input%mat1(j,i)
         end do
      end do
   case('L')! matrix is lower (and assumed symmetric) copy values in upper triangle
      do,i=1,input%nval2
         do,j=1,i-1
            input%mat1(j,i)=input%mat1(i,j)
         end do
      end do
   case default
      write(stderr,*)"ERROR: sort_BINtype, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select
   
case('FULLSQV0','FULLSQVN','FULL2DVN')! do nothing just check
   if( iboth .and. nrows .ne. ncols)then
      write(stderr,*)"ERROR: sort_BINtype, Matrix is not square!!"
      stop
   end if
   select case(input%mtyp)
   case('F')! do nothing accept
   case default
      write(stderr,*)"ERROR: sort_BINtype, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select
case('DIAVN___')
   if( .not. iboth)then
      write(stderr,*)"ERROR: sort_BINtype",input%type, "requires both sides to be permuted"
      stop
   end if

   select case(input%mtyp)
   case('P') ! do nothing just accept
   case default
      write(stderr,*)"ERROR: sort_BINtype, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select
! case('BDSYMV0_','BDSYMVN_')
!    blockm=.true.
!    select case(dat%mtyp)
!    case('B')! copy block diagonal upper triangle in lower triangle
!       sz=blockind(1)
!       shift=0
!       do,i=1,sz
!          do,j=1,i-1
!             input%mat1(shift+(j-1)*sz+i,1)=input%mat1(shift+(i-1)*sz+j,1)
!          end do
!       end do
      
!       !loop over remaining blocks
!       do,k=2,dat%nblocks
!          shift=shift+sz*2
!          sz=dat%blockind(k)-dat%blockind(k-1)
!          do,i=1,sz
!             do,j=1,i-1
!                input%mat1(shift+(j-1)*sz+i,1)=input%mat1(shift+(i-1)*sz+j,1)
!             end do
!          end do

!       end do
!    end select
! case('BDFULLV0','BDFULLVN')! do nothing type is accepted   
!    blockm=.true.
case default ! else write error message
   write(stderr,*)"ERROR: sort_BINtype, Matrix system not supported:",input%type
   stop
end select

! call cpu_time(time2)
! write(*,*)'Copying upper<->lower',time2-time1


!create compacted permutation vectors ( remove zeros)
!adapt indices if necessary
if(present(side1))then
   nrows=size(side1)-count(side1 .eq. 0)
   if(nrows .eq. size(side1))then ! we can directly assign, no copying necessary
      compact=>side1
   else ! copying necessary (get rid of zeros)
      allocate(compact(nrows)) ! allocate 
      row=0
      do,i=1,size(side1)
         if(side1(i) .eq. 0)cycle
         row=row+1
         compact(row)=side1(i)
      end do
   end if
end if


!!same for side 2
if(present(side2))then
   ncols=size(side2)-count(side2 .eq. 0)
   if(ncols .eq. size(side2))then ! we can directly assign, no copying necessary
      compact2=>side2
   else ! copying necessary (get rid of zeros)
      allocate(compact2(ncols)) ! allocate
      col=0
      do,i=1,size(side2)
         if(side2(i) .eq. 0)cycle
         col=col+1
         compact2(col)=side2(i)
      end do
   end if
end if


! permute vectors and side description data
if(associated(compact))then
   !vectors and side descriptions
   call permute(A=input%side1_d(1:nrows),perm=compact(1:nrows))


   if(input%nvec >0)call permute_rows(A=input%vec(1:nrows,1:input%nvec),perm=compact(1:nrows))

!permute matrix only when iboth is false
   if( .not. iboth)call permute_rows(A=input%mat1(1:nrows,1:ncols),perm=compact(1:nrows))
end if

!side 2 ( currently only side description)
if( associated(compact2))then
   call permute(A=input%side2_d(1:ncols),perm=compact2(1:ncols))

   !permute matrix only when iboth is false
   if( .not. iboth)call permute_cols(A=input%mat1(1:nrows,1:ncols),perm=compact2(1:ncols))
end if




! permute matrix
if(iboth)then ! permute both sides of the matrix with compact
   ! !matrix
   select case(input%type)
   case('DIAVN___')! permute in 1 step
      call permute(A=input%pack1(1:nrows),perm=compact(1:nrows))
   case('SYMV0___','SYMVN___','SYMV1___','SYMV2___') !symmetric matrix permute in 2 steps
       call permute_rows(A=input%mat1(1:nrows,1:nrows),perm=compact(1:nrows))
       call permute_cols(A=input%mat1(1:nrows,1:nrows),perm=compact(1:nrows))
!     call permute_sym(A=input%mat1(1:nrows,1:nrows),perm=compact(1:nrows)) ! symmetric permutation leaves lower triangle unchanged
   case default ! both sides in two steps
      call permute_rows(A=input%mat1(1:nrows,1:nrows),perm=compact(1:nrows))
      call permute_cols(A=input%mat1(1:nrows,1:nrows),perm=compact(1:nrows))
   end select
end if

end subroutine sort_BINtype

function packindex(dat,row,col)
  implicit none
  type(BINdat)::dat
  integer::row,col
  integer*8::packindex
  integer::mi,ma,st,sz,j

  select case(dat%type)
  case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')
     if(row > col)then
        packindex=0
        return
     end if
        
     packindex=row+(col*(col-1))/2
     return
     
  case('BDSYMV0_','BDSYMVN_')

     !quick return if in lower triangle or outside maximum diagonal band
     if(row > col .or. col-row > dat%maxblocksz)then
        packindex=0
        return
     end if

     ! find the block the points are both belonging to (or not)
     st=0
     packindex=0
     do,j=1,dat%nblocks
        if(row <= dat%blockind(j))then
           if(col <= dat%blockind(j) .and. col > st)then ! point is inside block
              
              packindex=packindex+(row-st)+((col-st)*(col-st-1))/2
              return
           else ! quick return
              packindex=0
              return
           end if
        else ! block still not found
           sz=dat%blockind(j)-st
           packindex=packindex+sz*(sz+1)/2
           !reset st parameter
           st=dat%blockind(j)
        end if
     end do

  case('BDFULLV0','BDFULLVN')
     !quick return point is outside maximum diagonal band
     if(abs(col-row) > dat%maxblocksz)then
        packindex=0
        return
     end if

     ! find the block the points are both belonging to (or not)
     st=0
     packindex=0
     do,j=1,dat%nblocks

        if(row <= dat%blockind(j))then
           if(col <= dat%blockind(j) .and. col > st)then ! point is inside block
              
              packindex=packindex+(col-st-1)*(dat%blockind(j)-st)+row-st ! column major order

              return
           else ! quick return
              packindex=0
              return
           end if
        else ! still not found
           sz=dat%blockind(j)-st
           packindex=packindex+sz**2
           !reset st parameter
           st=dat%blockind(j)
        end if
     end do
  case('FULLSQV0','FULLSQVN','FULL2DVN')
     packindex=(col-1)*dat%nval1+row
  case('DIAVN___')!diagonal matrix
     if(row .eq. col)then
        packindex=row
     else
        packindex=0
     end if
  end select

end function packindex


end module BINfiletools




!!!!!!!!!!!!!!!!!!!BELOW ARE EXTERNAL ROUTINES!!!!!!!!!!!!!!!!!!!




! !!subroutine to swap bytes
! !! input MUST be a contigious patch of memory containing variables of size bytes

!!subroutine is superseded by cswap ( C routine)
! !! the subroutine must be declared with an implict interface else the compiler will complain about non matching types
! !!input may be anything which needs to be swapped
! subroutine swap(input,bytes,sz)
! implicit none
! integer::sz,bytes
! integer*1::input(bytes,sz),dum !the input type is interpreted as a byte array
! integer::i,j,half
! half=bytes/2 !NOTE INTEGER DIVISION (middle byte does not need to be swopped if there is any)
! do,i=1,sz
!    do,j=1,half
!       dum=input(j,i)! copy value in dummy
!       input(j,i)=input(bytes-j+1,i)
!       input(bytes-j+1,i)=dum
!    end do
! end do
! end subroutine swap




!file containing tools to work with binary data files storing normal equations or covariance matrices
!matrix must be square and symmetric
!!Coded by Roelof Rietbroek, Wed Aug 22 12:09:18 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!version 1.0
!binary file contains:
!first record:
!13 parameter meta data:
	!1 parameter (integer) type: 1 (covariance) matrix 2. normal equations
	!2 parameter(integer): size of covariance matrix (here nval)
	!3 parameter (integer):endian test parameter should always be 1
	!4 parameter (integer):lmax (only for spherical harmonics)
	!5 parameter (integer):lmin (only for spherical harmonics)
	!6 parameter(double): reference time
	!7 parameter (double)::unused yet
	!8 parameter(double)::unused yet
	!9 parameter (double): unused yet
	!10 parameter (character(80)): description of data set
	!11 parameter (character(80)): unused yet
	!12 parameter (character(80))::unused yet
	!13 parameter (character(80))::unused yet
!second record
!character array with description of the sides
	!size: len=50, with dimension same as array side:nval
!third record
!Symmetric matrix in packed form
	!1 dimensional array with size nval+(nval*(nval-1))/2
!fourth record (if type .eq. 2)
	!normal equation vector:  A'C b size consistent with covariance array nval



!subroutine to write data to file

subroutine write_sym_mat(file,mat,matpack,vec,list,lmax,lmin,time,d2,d3,d4,descr,char2,char3,char4)
implicit none
character(*),intent(in)::file
double precision,intent(in),dimension(:,:),optional::mat !either provide a full 2 dimesnional matrix 
double precision,intent(in),dimension(:),optional::matpack !or a matrix already in packed form
double precision,intent(in),dimension(:),optional::vec
character(*),intent(in),optional,dimension(:)::list
integer,intent(in),optional::lmax,lmin
double precision,intent(in),optional::time,d2,d3,d4
character(*),optional,intent(in)::descr,char2,char3,char4

!private arguments
integer::endian,unit,typ,ilmax,ilmin,i,j,nval,stderr
double precision::id2,id3,id4,itime
character(80)::idescr,ichar2,ichar3,ichar4
character(50),allocatable,dimension(:)::ilist
double precision,allocatable,dimension(:)::packed

!return error if both packed matrix and normal matrix is provided
if(present(mat) .and. present(matpack))then
   write(stderr,*)' Matrix ambiguity, provide either a packed from matrix or a common matrix'
   stop
else if( (.not. present(mat)) .and. (.not. present(matpack)))then !or if no matrix is provided at all
   write(stderr,*)' No matrix provided, provide either a packed from matrix or a common matrix'
   stop
end if



!get size of data

if(present(mat))then
   nval=size(mat,1)
else !packed matrix is given
   nval=size(matpack,1)
   nval=(int(sqrt(8*nval+1.d0)+0.1)-1)/2
end if



!defaults
unit=2
endian=1
ilmax=0
ilmin=0
itime=0.d0
id2=0.d0
id3=0.d0
id4=0.d0
idescr='Description unknown'
ichar2=''
ichar3=''
ichar4=''



!get data size

nval=size(mat,1)

allocate(ilist(nval))

ilist=''

!process optional arguments
if(present(vec))then
   typ=2 !normal equations stored
else
   typ=1 !matrix stored
end if

if(present(lmax))ilmax=lmax
if(present(lmin))ilmin=lmin
if(present(list))ilist=list

if(present(time))itime=time

if(present(d2))id2=d2
if(present(d3))id3=d3
if(present(d4))id4=d4

if(present(descr))idescr=descr(1:80)

if(present(char2))ichar2=char2
if(present(char3))ichar3=char3
if(present(char4))ichar4=char4

!make packed matrix if square matrix is given
if(present(mat))then
   allocate(packed(nval+(nval*(nval-1))/2))
   do,i=1,nval
      do,j=i,nval
         packed(i+(j*(j-1))/2)=mat(i,j)
      end do
   end do

end if

   



open(unit=unit,file=trim(file),form='unformatted')

!first record meta data
write(unit)typ,nval,endian,ilmax,ilmin,itime,id2,id3,id4,idescr,ichar2,ichar3,ichar4
!second record description of sides of covariance matrix size character array(60)(nval) 
write(unit)ilist
if(present(mat))then
!   write(*,*)'writing mat'
   !third record matrix

   write(unit)packed

else

   write(unit)matpack
end if

!fourth record vector
if(typ .eq. 2)write(unit)vec

close(unit)



end subroutine write_sym_mat

!subroutine which reads the above matrix into arrays and inquires for the meta data
subroutine read_sym_mat(file,mat,matpack,typ,nval,vec,list,lmax,lmin,time,d2,d3,d4,descr,char2,char3,char4)
use GPStlbx
implicit none
character(*),intent(in)::file
double precision,intent(out),optional,dimension(:)::matpack,vec
double precision,intent(out),optional,dimension(:,:)::mat
integer,intent(out),optional::typ,nval,lmax,lmin
character(*),intent(out),optional::descr,char2,char3,char4
double precision,optional,intent(out)::time,d2,d3,d4
character(*),optional,dimension(:),intent(out)::list

character(50),dimension(:),allocatable::ilist
integer::unit,stderr,endian,ilmax,ilmin,i,j,inval,ityp
double precision::id2,id3,id4,itime
character(80)::idescr,ichar2,ichar3,ichar4
double precision,allocatable,dimension(:)::imatpack

!defaults
stderr=0
unit=13
endian=1
ilmax=0
ilmin=0
itime=0.d0
id2=0.d0
id3=0.d0
id4=0.d0
idescr='Description unknown'
ichar2=''
ichar3=''
ichar4=''
!write(*,*)' ',trim(file)
!open file for unformatted acces
open(file=trim(file),unit=unit,form='unformatted')
!read first record (meta data
read(unit)ityp,inval,endian,ilmax,ilmin,itime,id2,id3,id4,idescr,ichar2,ichar3,ichar4

!check whether the endian parameter is converted back to 1 and check whether this could be due to the endiannes
if(endian .ne. 1)then
   !swap bytes to check whether endiannes is the problem
   call i4swap(1,endian)
   if( .not. ltlend() .and. (endian .eq. 1))then
      write(stderr,*)'file, ',trim(file), 'is little endian and system big endian, make sure',&
           ' compilers options are set to read unformatted files as little endian)'
      stop
   else if(endian .eq. 1 .and. ltlend())then
      write(stderr,*)'file, ',trim(file), 'is big endian and system little endian, make sure',&
           ' compilers options are set to read unformatted files as big endian)'
      stop
   else if(endian .ne. 1)then
      write(stderr,*)'file, ',trim(file),'ERROR reading file , (caused by uncommon record header length?)'
      stop
   end if
end if

if(present(nval))nval=inval
if(present(typ))typ=ityp

if(present(lmax))lmax=ilmax
if(present(lmin))lmin=ilmin

if(present(time))time=itime
if(present(d2))d2=id2
if(present(d3))d3=id3
if(present(d4))d4=id4
if(present(descr))descr=idescr
if(present(char2))char2=ichar2
if(present(char3))char3=ichar3
if(present(char4))char4=ichar4

!read second record
allocate(ilist(inval))
read(unit)ilist
if(present(list))then
   if(size(list,1).ne. inval)then
      write(stderr,*)'dimensions of list should be: ',inval
   end if
   list=ilist
end if

!read third record only when vec or mat is requested
if( present(vec) .or. present(mat) .or. present (matpack))then
   
   allocate(imatpack(inval+(inval*(inval-1)/2)))


   read(unit)imatpack

   if(present(matpack))then !if packed matrix is requested
      if(size(matpack,1).ne. inval+(inval*(inval-1))/2)then
         write(stderr,*)'dimensions of matpack should be: ',inval+(inval*(inval-1))/2
         stop
      end if
      matpack=imatpack
   end if

   if(present(mat))then !if matrix is requested
      if(size(mat,1).ne. inval .or. (size(mat,2).ne. inval) )then
         write(stderr,*)'dimensions of mat should be: ',inval,inval
         stop
      end if
      do,i=1,inval
         do,j=i,inval
            mat(i,j)=imatpack(i+(j*(j-1))/2)
            mat(j,i)=mat(i,j)
         end do
      end do
      
   end if

   if(present(vec))then !read fourth record
      if(size(vec,1) .ne. inval)then
         write(stderr,*)'dimension of vec should be ',inval
      end if
      read(unit)vec
   end if


end if

close(unit)

end subroutine read_sym_mat


!function which checks the endiannes of the machine returns true for little endian
!taken from the altim library written by Remko Scharroo
function ltlend()
logical ltlend
integer*2 itest(2)
character*2 ctest(2)
equivalence (ctest,itest)
!data ctest /2h12,2h21/
data ctest /'12','21'/
!write(*,*)itest(1),itest(2)
ltlend=itest(1).gt.itest(2)
end function ltlend
  
! *
! * This function returns .TRUE. if your system uses the LITTLE ENDIAN
! * notation of INTEGERS. For BIG ENDIAN notation, it returns .FALSE.
! *
! * The terms LITTLE ENDIAN and BIG ENDIAN refer to the order of the
! * bytes in an INTEGER word, both for two-byte and four-byte integers.
! * Within a BYTE, the BITS are always ordered most significant to
! * least significant. Within a WORD the BYTES may be ordered:
! * - Most significant to least significant: BIG endian
! *   (IBM RS6000, Sun, SGI, Convex, MacIntosh).
! * - Least significant to most significant: LITTLE endian
! *   (PC, VAX, DEC).
! *
! * Arguments:
! *  (none)
! * Returned value:
! *  LTLEND : .TRUE. for LITTLE endian notation of integers
! *           .FALSE. for BIG endian notation of integers
! *-
! *  5-Jan-1996 - Created by Remko Scharroo
! *  8-Aug-2000 - Changed to check integer*2 in stead of integer*4
! *-----------------------------------------------------------------------
      
  



!**I4SWAP -- Swap one or more 4-byte integers
!*+
      SUBROUTINE I4SWAP (N, I)
      INTEGER*4 N
! *     INTEGER*4 I(N)
! *
! * This routine swaps a 4-byte integer variable or the elements of a
! * 4-byte integer array, i.e. it converts their representation between
! * little and big endian (in either direction).
! *
! * Arguments:
! *  N  (input) : Number of elements to swap
! *  I  (input) : Integer or integer array to swap
! *    (output) : Swapped integer or integer array
! *
! * Examples:
! *     CALL I4SWAP(1,I)    ! Swap integer I
! *     CALL I4SWAP(4,I)    ! Swap integer elements I(1) through I(4)
! *     CALL I4SWAP(3,I(7)) ! Swap integer elements I(7) through I(9)
! *-
! *    Aug-1991 - Created: Remko Scharroo, DUT/SSR&T (c)
! *  2-Jul-1992 - Documented version
! * 12-Nov-1993 - CHARACTER*1 replaced by BYTE
! *-----------------------------------------------------------------------
      BYTE I(4,N),B
      INTEGER*4 K
      DO K=1,N
         B=I(1,K)
         I(1,K)=I(4,K)
         I(4,K)=B
         B=I(2,K)
         I(2,K)=I(3,K)
         I(3,K)=B
      ENDDO
      END



!file with various general small routines and functions

!!Updated by Roelof Rietbroek, Fri Apr 30 11:50:24 2010

function DMtrace(mat)
    double precision,intent(in),dimension(:,:)::mat
    double precision::DMtrace
    integer(8)::i

    if (size(mat,1) .ne. size(mat,2))then
        write(0,*)"ERROR in DMtrace: matrix not square"
        stop 
    end if
    DMtrace=0.d0
    do,i=1,size(mat,1)
        DMtrace=DMtrace+mat(i,i)
    end do

end function

!subroutine to reallocate a 1D pointer (double precision) ( and copy data)
subroutine reallocate_dptr(in,nadd)
implicit none
double precision,pointer::in(:)
integer,intent(in)::nadd

double precision,pointer::tmp(:)=>null()
integer*8::n
n=0


if(associated(in))then
   n=size(in)
   !let tmp point to the same memory as in
   tmp=>in
end if
!nullify in pointer
in=>null()
!allocate in pointer with more data
allocate(in(n+nadd))

in=0.d0 ! initialize to zero
!copy data back
if(n>0)then
   in(1:n)=tmp(1:n)

   !dealocate tmp data
   deallocate(tmp)
end if
end subroutine reallocate_dptr


subroutine reallocate_dptrm(in,xadd,yadd)
implicit none
double precision,pointer::in(:,:)
integer,intent(in)::xadd,yadd

double precision,pointer::tmp(:,:)=> null()
integer*8::nx,ny
nx=0
ny=0
if(associated(in))then
   nx=size(in,1)
   ny=size(in,2)
   !let tmp point to the same memory as in
   tmp=>in
end if

!nullify in pointer
in=>null()

!allocate in pointer with more data
allocate(in(nx+xadd,ny+yadd))
in=0.d0 ! initialize to zero
!copy data back
if(nx > 0 .or. ny > 0)then
   in(1:nx,1:ny)=tmp
   !dealocate tmp data
   deallocate(tmp)
end if
end subroutine reallocate_dptrm

!subroutine to reallocate a 1D pointer (character vector) ( and copy data)
subroutine reallocate_cptr(in,nadd)
implicit none
character(*),pointer,intent(inout)::in(:)
integer,intent(in)::nadd

!character(*),pointer::tmp(:)=> null()
character(len(in(1))),pointer::tmp(:)
integer*8::n

n=0
!!IMPORTANT NOTE the use of associated() requires the pointer to be correctly initiated (allocated or nullified)
!! failing to do so may cause segfaults when the build contains optimization flags
if(associated(in))then
   n=size(in)
   !let tmp point to the same memory as in
   tmp=>in
end if

!nullify in pointer
in=>null()
!allocate in pointer with more data
allocate(in(n+nadd))
in='' ! initialize
if(n>0)then
   !copy data back
   in(1:n)=tmp
   
!dealocate tmp data
   deallocate(tmp)
end if


end subroutine reallocate_cptr

subroutine reallocate_iptr(in,nadd)
implicit none
integer,pointer::in(:)
integer,intent(in)::nadd

integer,pointer::tmp(:) => null()
integer*8::n
n=0
if(associated(in))then
   n=size(in)
   !let tmp point to the same memory as in
   tmp=>in
end if
!nullify in pointer
in=>null()
!allocate in pointer with more data
allocate(in(n+nadd))
in=0 ! initialize to zero
if(n>0)then
   !copy data back
   in(1:n)=tmp

   !dealocate tmp data
   deallocate(tmp)
end if

end subroutine reallocate_iptr

!same but for long integer
subroutine reallocate_liptr(in,nadd)
implicit none
integer*8,pointer::in(:)
integer,intent(in)::nadd

integer*8,pointer::tmp(:) => null()
integer*8::n
n=0
if(associated(in))then
   n=size(in)
   !let tmp point to the same memory as in
   tmp=>in
end if
!nullify in pointer
in=>null()
!allocate in pointer with more data
allocate(in(n+nadd))
in=0 ! initialize to zero
if(n>0)then
   !copy data back
   in(1:n)=tmp

   !dealocate tmp data
   deallocate(tmp)
end if

end subroutine reallocate_liptr

!**I4SWAP -- Swap one or more 4-byte integers
!*+
!       SUBROUTINE I4SWAP (N, I)
!       INTEGER*4 N
! ! *     INTEGER*4 I(N)
! ! *
! ! * This routine swaps a 4-byte integer variable or the elements of a
! ! * 4-byte integer array, i.e. it converts their representation between
! ! * little and big endian (in either direction).
! ! *
! ! * Arguments:
! ! *  N  (input) : Number of elements to swap
! ! *  I  (input) : Integer or integer array to swap
! ! *    (output) : Swapped integer or integer array
! ! *
! ! * Examples:
! ! *     CALL I4SWAP(1,I)    ! Swap integer I
! ! *     CALL I4SWAP(4,I)    ! Swap integer elements I(1) through I(4)
! ! *     CALL I4SWAP(3,I(7)) ! Swap integer elements I(7) through I(9)
! ! *-
! ! *    Aug-1991 - Created: Remko Scharroo, DUT/SSR&T (c)
! ! *  2-Jul-1992 - Documented version
! ! * 12-Nov-1993 - CHARACTER*1 replaced by BYTE
! ! *-----------------------------------------------------------------------
!       BYTE I(4,N),B
!       INTEGER*4 K
!       DO K=1,N
!          B=I(1,K)
!          I(1,K)=I(4,K)
!          I(4,K)=B
!          B=I(2,K)
!          I(2,K)=I(3,K)
!          I(3,K)=B
!       ENDDO
!     END SUBROUTINE I4SWAP


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
  
!*
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


function freeunit()
  implicit none
  integer::freeunit,i
  logical::uopen,uexist
  uopen=.true.
!loop over potential free units
  do,i=15,99
     inquire(unit=i,opened=uopen)

     if(.not. uopen)then
        freeunit=i
        return
     end if
  end do
end function freeunit

!subroutine to rearrange the rows of a double vector according to a given permutation
!algorithm adapted from slatec DDPERM function
Subroutine permute_d(A,perm,inv,err)
implicit none
double precision,intent(inout)::A(:)
integer,intent(inout)::perm(:)
logical,optional,intent(in)::inv
integer,optional,intent(out)::err
double precision::tmp ! automatic array to store a complete row at once
integer::n,i,stderr,indx,indx0,istrt,ncul
logical::iinv
if(present(err))err=0
iinv=.false.
if(present(inv))iinv=inv

stderr=0
ncul=0
n=size(A,1)


if(n .ne. size(perm))then
   write(stderr,*)"ERROR: permutation vector size is not consistent with the rows of A"
   stop
end if

!check whether the permutation is valid
do,i=1,n
   indx=abs(perm(i))
   if(indx >= 1 .and. indx <=n)then
      if(perm(indx)>0)then
         perm(indx)=-perm(indx)
         cycle
      end if
   end if
   if(present(err))then
      err=-1 ! set errorcode and return
      return
   else
      write(stderr,*)"ERROR: invalid permutation"


      stop
   end if
end do
  
!!all perm values are now negative!!!
istrt=0
ncul=0
if(.not. iinv)then
   do while(ncul < n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      
      indx=istrt
      perm(indx)=-perm(indx)
      indx0=perm(indx)
      ncul=ncul+1
      do while(perm(indx0) < 0)
         !swap values
         tmp=A(indx0)
         A(indx0)=A(indx)
         A(indx)=tmp
         perm(indx0)=-perm(indx0)
         indx=indx0
         indx0=perm(indx0)
         ! reset to a positive value (functions as a tag that row has been permuted)
         
         ncul=ncul+1
      end do
      
   end do
else !!! inverse/backward permutation
   do while(ncul <n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      perm(istrt)=-perm(istrt)
      indx=perm(istrt)
      ncul=ncul+1
      do while(indx .ne. istrt)
         !swap values
         tmp=A(indx)
         A(indx)=A(istrt)
         A(istrt)=tmp
         perm(indx)=-perm(indx)
         indx=perm(indx)
         ncul=ncul+1
      end do
   end do

end if
!write(*,*)'ncul',ncul

end Subroutine permute_d

!subroutine to rearrange the rows of a integer vector according to a given permutation
!algorithm adapted from slatec DDPERM function
Subroutine permute_i(A,perm,inv,err)
implicit none
integer,intent(inout)::A(:)
integer,intent(inout)::perm(:)
logical,optional,intent(in)::inv
integer,optional,intent(out)::err
integer::tmp ! automatic array to store a complete row at once
integer::n,i,stderr,indx,indx0,istrt,ncul
logical::iinv
if(present(err))err=0
iinv=.false.
if(present(inv))iinv=inv

stderr=0

n=size(A,1)


if(n .ne. size(perm))then
   write(stderr,*)"ERROR: permutation vector size is not consistent with the rows of A"
   stop
end if

!check whether the permutation is valid
do,i=1,n
   indx=abs(perm(i))
   if(indx >= 1 .and. indx <=n)then
      if(perm(indx)>0)then
         perm(indx)=-perm(indx)
         cycle
      end if
   end if
   if(present(err))then
      err=-1 ! set errorcode and return
      return
   else
      write(stderr,*)"ERROR: invalid permutation"

      stop
   end if
end do
  
!!all perm values are now negative!!!
istrt=0
ncul=0
if(.not. iinv)then
   do while(ncul < n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      
      indx=istrt
      perm(indx)=-perm(indx)
      indx0=perm(indx)
      ncul=ncul+1
      do while(perm(indx0) < 0)
         !swap values
         tmp=A(indx0)
         A(indx0)=A(indx)
         A(indx)=tmp
         perm(indx0)=-perm(indx0)
         indx=indx0
         indx0=perm(indx0)
         ! reset to a positive value (functions as a tag that row has been permuted)
         
         ncul=ncul+1
      end do
      
   end do
else !!! inverse/backward permutation
   do while(ncul <n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      perm(istrt)=-perm(istrt)
      indx=perm(istrt)
      ncul=ncul+1
      do while(indx .ne. istrt)
         !swap values
         tmp=A(indx)
         A(indx)=A(istrt)
         A(istrt)=tmp
         perm(indx)=-perm(indx)
         indx=perm(indx)
         ncul=ncul+1
      end do
   end do

end if

!write(*,*)'ncul',ncul

end Subroutine permute_i

!subroutine to rearrange the rows of a CHARACTER vector according to a given permutation
!algorithm adapted from slatec DDPERM function
Subroutine permute_ch(A,perm,inv,err)
implicit none
Character(*),intent(inout)::A(:)
integer,intent(inout)::perm(:)
logical,optional,intent(in)::inv
integer,optional,intent(out)::err
Character(len(A(1)))::tmp ! automatic array to store a complete string at once
integer::n,i,stderr,indx,indx0,istrt,ncul
logical::iinv
iinv=.false.
if(present(err))err=0
if(present(inv))iinv=inv
stderr=0

n=size(A,1)


if(n .ne. size(perm))then
   write(stderr,*)"ERROR: permutation vector size is not consistent with the rows of A"
   stop
end if

! do,i=1,n
! write(*,*)perm(i)
! end do

!check whether the permutation is valid
do,i=1,n
   indx=abs(perm(i))
   if(indx >= 1 .and. indx <=n)then
      if(perm(indx)>0)then
         perm(indx)=-perm(indx)
         cycle
      end if
   end if
   if(present(err))then
      err=-1 ! set errorcode and return
      return
   else
      write(stderr,*)"ERROR: invalid permutation"
!      write(0,*)n,size(perm,1)
      stop
   end if
end do
  
!!all perm values are now negative!!!
ncul=0
istrt=0
if(.not. iinv)then
   do while(ncul < n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      
      indx=istrt
      perm(indx)=-perm(indx)
      indx0=perm(indx)
      ncul=ncul+1
      do while(perm(indx0) < 0)
         !swap values
         tmp=A(indx0)
         A(indx0)=A(indx)
         A(indx)=tmp
         perm(indx0)=-perm(indx0)
         indx=indx0
         indx0=perm(indx0)
         ! reset to a positive value (functions as a tag that row has been permuted)
         
         ncul=ncul+1
      end do
      
   end do
else !!! inverse/backward permutation
   do while(ncul <n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      perm(istrt)=-perm(istrt)
      indx=perm(istrt)
      ncul=ncul+1
      do while(indx .ne. istrt)
         !swap values
         tmp=A(indx)
         A(indx)=A(istrt)
         A(istrt)=tmp
         perm(indx)=-perm(indx)
         indx=perm(indx)
         ncul=ncul+1
      end do
   end do

end if

!write(*,*)'ncul',ncul

end Subroutine permute_ch

!subroutine to rearrange the rows of a double matrix according to a given permutation
!this is an inplace algorithm
Subroutine permute_rows(A,perm,inv,err)
implicit none
double precision,intent(inout)::A(:,:)
integer,intent(inout)::perm(:)
logical,optional,intent(in)::inv
integer,optional,intent(out)::err
double precision::tmp
integer::n,m,i,j,k,stderr,indx,indx0,istrt
integer::ncul
logical::switch,iinv
stderr=0
if(present(err))err=0
iinv=.false.
if(present(inv))iinv=inv
!   call printaddress(A)
n=size(A,1)
m=size(A,2)


if(n .ne. size(perm,1))then
   write(stderr,*)"ERROR: permutation vector size is not consistent with the rows of A"
   stop
end if

!check whether the permutation is valid
do,i=1,n
   indx=abs(perm(i))
   if(indx >= 1 .and. indx <=n)then
      if(perm(indx)>0)then
         perm(indx)=-perm(indx)
         cycle
      end if
   end if
   if(present(err))then
      err=-1 ! set errorcode and return
      return
   else
      write(stderr,*)"ERROR: invalid permutation"
      stop
   end if
end do

  
! ! !!all perm values are now negative!!!

if( .not. iinv)then ! forward permutation
   switch=.false.
   do,k=1,m !loop over columns

!      write(*,*)A(:,k)
      ncul=0
      istrt=0
      if(switch)then
         do while(ncul < n)
            istrt=istrt+1
            if(perm(istrt) < 0)cycle ! quick cycle if value has been permuted already
            
            indx=istrt
            indx0=perm(indx)
            perm(indx)=-perm(indx)
            ncul=ncul+1
            do while(perm(indx0) > 0)
               !swap values
               tmp=A(indx0,k)
               A(indx0,k)=A(indx,k)
               A(indx,k)=tmp
               indx=indx0
               indx0=perm(indx0)
               perm(indx)=-perm(indx)
               ncul=ncul+1
            end do
            
         end do
      else
         do while(ncul < n)
            istrt=istrt+1
            if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
            
            indx=istrt
            perm(indx)=-perm(indx)
            indx0=perm(indx)
            ncul=ncul+1
            do while(perm(indx0) < 0)
               !swap values
               tmp=A(indx0,k)
               A(indx0,k)=A(indx,k)
               A(indx,k)=tmp
               perm(indx0)=-perm(indx0)
               indx=indx0
               indx0=perm(indx0)
               ! reset to a positive value (functions as a tag that row has been permuted)
               
               ncul=ncul+1
            end do
            
         end do
      end if
      switch= (.not. switch) ! flip logical ( even and odd colum)
   end do ! end loop over columns
!   write(*,*)ncul,n
else ! Inverse/backward permutation
   switch=.false.
   do,k=1,m !loop over columns
      ncul=0
      istrt=0
      if(switch)then
         do while(ncul <n)
            istrt=istrt+1
            if(perm(istrt) < 0)cycle ! quick cycle if value has been permuted already
            
            indx=perm(istrt)
            perm(istrt)=-perm(istrt)
            ncul=ncul+1

            do while(indx .ne. istrt)
               !swap values
               tmp=A(indx,k)
               A(indx,k)=A(istrt,k)
               A(istrt,k)=tmp
               perm(indx)=-perm(indx)
               indx=-perm(indx)
               ncul=ncul+1
            end do
         end do
      else
         do while(ncul <n)
            istrt=istrt+1
            if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
            perm(istrt)=-perm(istrt)
            indx=perm(istrt)
            ncul=ncul+1
            do while(indx .ne. istrt)
               !swap values
               tmp=A(indx,k)
               A(indx,k)=A(istrt,k)
               A(istrt,k)=tmp
               perm(indx)=-perm(indx)
               indx=perm(indx)
               ncul=ncul+1
            end do
         end do
      end if
      switch= (.not. switch) ! flip logical ( even and odd colum)
   end do

end if

if(.not. switch)perm=-perm ! possibly switch back to input values

!write(*,*)'ncul',ncul

end Subroutine permute_rows

!same but permute the columns
!this is also an inplace permutation
Subroutine permute_cols(A,perm,inv,err)
implicit none
double precision,intent(inout)::A(:,:)
integer,intent(inout)::perm(:)
logical,optional,intent(in)::inv
integer,optional,intent(out)::err
double precision::tmp
integer::n,m,i,j,k,stderr,indx,indx0,istrt
integer::ncul
logical::iinv
stderr=0
if(present(err))err=0
iinv=.false. ! default is a forward permutation
if(present(inv))iinv=inv

m=size(A,1)
n=size(A,2)


if(n .ne. size(perm))then
   write(stderr,*)"ERROR: permutation vector size is not consistent with the rows of A"
   stop
end if

!check whether the permutation is valid
do,i=1,n
   indx=abs(perm(i))
   if(indx >= 1 .and. indx <=n)then
      if(perm(indx)>0)then
         perm(indx)=-perm(indx)
         cycle
      end if
   end if
   if(present(err))then
      err=-1 ! set errorcode and return
      return
   else
      write(stderr,*)"ERROR: invalid permutation"
      stop
   end if
end do


if( .not. iinv)then ! forward permutation
   ncul=0
   istrt=0
   !do istrt=1,n
   do while(ncul < n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      
      indx=istrt
      perm(indx)=-perm(indx)
      indx0=perm(indx)
      ncul=ncul+1
      do while(perm(indx0) < 0)
 !        write(*,*)indx,indx0
         do,i=1,m ! inner loop over rows
            !swap values
            tmp=A(i,indx0)
            A(i,indx0)=A(i,indx)
            A(i,indx)=tmp
         end do
         perm(indx0)=-perm(indx0)
         indx=indx0
         indx0=perm(indx0)
         ! reset to a positive value (functions as a tag that row has been permuted)
         
         ncul=ncul+1
      end do
      
   end do
!   write(*,*)ncul,n
else ! Inverse/backward permutation
   ncul=0
   istrt=0
   !do istrt=1,n
   do while(ncul <n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      

      perm(istrt)=-perm(istrt)
      indx=perm(istrt)
      ncul=ncul+1
      do while(indx .ne. istrt)
         
         do,i=1,m ! inner loop over rows
            !swap values
            tmp=A(i,indx)
            A(i,indx)=A(i,istrt)
            A(i,istrt)=tmp
         end do
         perm(indx)=-perm(indx)
         indx=perm(indx)
         ! reset to a positive value (functions as a tag that row has been permuted)
         
         ncul=ncul+1
      end do
      
   end do
end if


!write(*,*)'ncul',ncul
end Subroutine permute_cols


! subroutine to permute both rows and columns of a symmetric matrix without referencing the lower triangle
! principle:
! the new upper triangle is a permutation of the old upper triangle
subroutine permute_sym(A,perm,inv)
implicit none
double precision,intent(inout)::A(:,:)
integer,intent(inout)::perm(:)
logical,optional,intent(in)::inv


!!private variables
double precision::tmp,tmp2(size(perm,1))
integer::n,m,i,j,k,stderr,indx,indx0,istrt
integer::ncul,mn,mx
logical::iinv

stderr=0
write(stderr,*)" this routine does not work up to now :("
stop
iinv=.false. ! default is a forward permutation
if(present(inv))iinv=inv
m=size(A,1)
n=size(A,2)

if(m .ne. n)then
   write(stderr,*)"ERROR permute_sym: matrix A is not square"
   stop
end if

if(n .ne. size(perm))then
   write(stderr,*)"ERROR: permutation vector size is not consistent with the rows of A"
   stop
end if

!check whether the permutation is valid
do,i=1,n
   indx=abs(perm(i))
   if(indx >= 1 .and. indx <=n)then
      if(perm(indx)>0)then
         perm(indx)=-perm(indx)
         cycle
      end if
   end if
   write(stderr,*)"ERROR: invalid permutation"
   stop
end do


!do actual permutation
if( .not. iinv)then ! forward permutation
   ncul=0
   istrt=0
   !do istrt=1,n
   do while(ncul < n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      
      indx=istrt
      perm(indx)=-perm(indx)
      indx0=perm(indx)
      ncul=ncul+1
      do while(perm(indx0) < 0)
         mn=min(indx,indx0)
         mx=max(indx,indx0)

         !column exchange part
         do,i=1,mn-1
            tmp=A(i,mn)
            A(i,mn)=A(i,mx)
            A(i,mx)=tmp
         end do
         
         !row exchange part
         do,i=mx+1,n
            tmp=A(mn,i)
            A(mn,i)=A(mx,i)
            A(mx,i)=tmp
         end do
         
         !cross term part
         do,i=mn+1,mx-1
            tmp=A(mn,i)
            A(mn,i)=A(i,mx)
            A(i,mx)=tmp
         end do
         
         !diagonal exchange
         tmp=A(mn,mx)
         A(mn,mn)=A(mx,mx)
         A(mx,mx)=tmp


         perm(indx0)=-perm(indx0)
         indx=indx0
         indx0=perm(indx0)
         ! reset to a positive value (functions as a tag that row has been permuted)
      
         ncul=ncul+1
      end do
      
   end do
else ! backward permutation
   ncul=0
   istrt=0
   !do istrt=1,n
   do while(ncul <n)
      istrt=istrt+1
      if(perm(istrt) > 0)cycle ! quick cycle if value has been permuted already
      perm(istrt)=-perm(istrt)
      indx=perm(istrt)
      ncul=ncul+1
      do while(indx .ne. istrt)
         mn=min(indx,istrt)
         mx=max(indx,istrt)
         
         !column exchange part
         do,i=1,mn-1
            tmp=A(i,mn)
            A(i,mn)=A(i,mx)
            A(i,mx)=tmp
         end do
         
         !row exchange part
         do,i=mx+1,n
            tmp=A(mn,i)
            A(mn,i)=A(mx,i)
            A(mx,i)=tmp
         end do
         
         !cross term part
         do,i=mn,mx-1
            tmp=A(mn,i)
            A(mn,i)=A(mx+mn-i,mx)
            A(mx+mn-i,mx)=tmp
         end do

         perm(indx)=-perm(indx)
         indx=perm(indx)
         ! reset to a positive value (functions as a tag that row has been permuted)
         
         ncul=ncul+1
      end do
      
   end do
end if
end subroutine permute_sym


!subroutine which gives back a permutation vector for side1, which is sorted according to side2 when permuted with the vector
! entries present in side2 but not in side1 will be ignored
! entries which are present in side1 but not in side2 will be put at the back ( when back=true) or at the front of the permutation vector
! the routine is essentially a fortran wrapper for the C routine get_permvecc(), which uses efficient qsort routines from the Standard C lib 
subroutine get_permvec(side1,side2,back,st,nd,perm,ncom)
implicit none
character(*),intent(in),dimension(:)::side1,side2 ! character description of the sides to be matched
logical,intent(in)::back ! logical which determines to do with entries presetn in side1 but no tin side2
!if back=.true. put those entries at the back of the permutation else in front of the permutation
integer,intent(in)::st,nd ! start and end of the string part to compare 
integer,intent(out)::perm(:) ! permutation vector for side1 to match ( ignoring entries present in side2 but not in side1) side2
integer,intent(out)::ncom ! number of common parameters found

!private variables
integer::n1,n2,strln,get_permvecC
n1=size(side1)
n2=size(side2)

strln=len(side1(1))

if(back)then
   ncom=get_permvecC(n1,n2,strln,st,nd,1,perm(1),side1,side2)
else
   ncom=get_permvecC(n1,n2,strln,st,nd,0,perm(1),side1,side2)
end if

end subroutine get_permvec

!subroutien which guives back a permutation vector for side1, which is sorted according to side2 when permuted with the vector
! entries present in side2 but not in side1 will be ignored
! entries which are present in side1 but not in side2 will be put at the front of the permutation vector
! the routine is essentially a fortran wrapper for the C routine get_permvecc(), which uses efficient qsort routines from the Standard C lib 
subroutine get_permvec2(side1,side2,back,st,nd,perm1,perm2,ncom)
implicit none
character(*),intent(in),dimension(:)::side1,side2 ! character description of the sides to be matched
integer,intent(in)::st,nd ! start and end of the string part to compare 
integer,intent(out)::perm1(:)! permutation vector for side1 to match ( ignoring entries present in side2 but not in side1) side2 
integer,intent(out)::perm2(:) !and for side 2 to match side1 (additional permutation of side1 may be needed)
integer,intent(out)::ncom ! number of common parameters found
logical,intent(in)::back ! logical which determines what to do with unique entries in side1  (stack at front or back)


!private variables
integer::n1,n2,strln,get_permvecC2
n1=size(side1)
n2=size(side2)
if(n1 .ne. size(perm1) .or. n2 .ne. size(perm2))then
   write(0,*)"permutation vectors don't match the character arrays"
   stop
end if


strln=len(side1(1))
if(back)then
   ncom=get_permvecC2(n1,n2,strln,st,nd,1,perm1(1),perm2(1),side1,side2)
else
   ncom=get_permvecC2(n1,n2,strln,st,nd,0,perm1(1),perm2(1),side1,side2)
end if
end subroutine get_permvec2



!old fortran routine ( replaced by get_permvec, which is much quicker)
!subroutine to retrieve the permutation vector
!for side1 to match side2
subroutine get_permvec_old(side1,side2,back,st,nd,perm,ncom)
implicit none
character(*),intent(in),dimension(:)::side1,side2 ! character description of the sides to be matched
logical,intent(in)::back ! logical which determines to do with entries presetn in side1 but no tin side2
!if back=.true. put those entries at the back of the permutation else in front of the permutation
integer,intent(in)::st,nd ! start and end of the string part to compare 
integer,intent(out)::perm(:) ! permutation vector for side1 to match ( ignoring entries present in side2 but not in side1) side2
integer,intent(out)::ncom ! number of common parameters found


!private variables
integer::i,j,n1,n2,ind
integer::nonmatch,st_2,nd_2
logical::found
integer,allocatable::perm2(:)
n1=size(side1)
n2=size(side2)
allocate(perm2(n2))
perm=0
perm2=-1

ncom=0!amount of common parameters found
nonmatch=0
st_2=1 ! initial starting point for search loop
nd_2=n2 ! initial end point for search loop


do,i=1,n1 ! loop over side1 entries
   found=.false.
   do,j=st_2,nd_2  !search for entry in side2
      if(side2(j)(st:nd) .eq. side1(i)(st:nd))then ! found
         found=.true.
         ncom=ncom+1
         perm2(j)=i
         !decrease search loop borders if the match is on one of the boundaries ( this might accelerate some common searches)
         if(st_2 .eq. j)then
            st_2=st_2+1
         else if(nd_2 .eq. j)then
            nd_2=nd_2-1
         end if
!         write(0,*)ncom,j,i,side2(j)(st:nd),side1(i)(st:nd)
         exit !search loop
      end if
   end do
   if( .not. found)then ! put value at the back or front of perm
      nonmatch=nonmatch+1
      if(back)then
         perm(n1-nonmatch+1)=i
      else
         perm(nonmatch)=i
      end if
   end if
end do


!now check for -1 entries in perm2 and fill in the common part of perm

if(back)then
   ind=0
   do,i=1,n2
      if(perm2(i)< 0)cycle ! ignore entry
      ind=ind+1
      perm(ind)=perm2(i)
   end do
else
   ind=n1-ncom
   do,i=1,n2
      if(perm2(i)< 0)cycle ! ignore entry
      ind=ind+1
      perm(ind)=perm2(i)
   end do
end if


deallocate(perm2)

end subroutine get_permvec_old

!!Coded by Roelof Rietbroek, Fri Apr 19 16:29:26 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de
!! function which strips leading directories and suffices from filenames
subroutine basenamef(stripped,path,suf)
implicit none
character(*),intent(out)::stripped
character(*),intent(in)::path
character(*),intent(in),optional::suf

integer::st,nd

st=index(path,'/',back=.true.)
if(present(suf))then
   nd=index(path,trim(suf),back=.true.)
   stripped=trim(path(st+1:nd-1))
else
   stripped=trim(path(st+1:))
end if


end subroutine basenamef


! program test_permute
! use FORTtlbx
! implicit none
! integer,parameter::n=5
! double precision,dimension(n,n)::A
! double precision::vecd(5),time1,time2
!  character(30)::vecch(5)
!  integer::veci(5)
! integer::i,j
! integer::permrow(5)=(/5,1,3,2,4/)
! integer::permcol(5)=(/1,3,2,5,4/)
! !integer::permcol(5)=(/5,1,2,3,4/)

! A=0.d0

! forall(i=1:n,j=1:n,i<=j)A(i,j)=dble(i*j)

! write(*,*)'BEFORE'
! write(*,*)permrow

! do,i=1,5
!    write(*,*)A(i,:)
! end do


! ! call permute_rows(A,permrow)
! ! write(*,*)'AFTER row permutation'
! ! write(*,*)permrow
! ! do,i=1,5
! !    write(*,*)A(i,:)
! ! end do

! ! call permute_rows(A,permrow,.true.)
! ! write(*,*)'AFTER inverse row permutation'
! ! write(*,*)permrow
! ! do,i=1,5
! !    write(*,*)A(i,:)
! ! end do




! ! ! !permute columns

! ! ! write(*,*)'BEFORE'
! ! ! write(*,*)permcol

! ! ! do,i=1,5
! ! !    write(*,*)A(i,:)
! ! ! end do

! ! call permute_cols(A,permcol)
! ! write(*,*)'AFTER col permutation'
! ! write(*,*)permcol
! ! do,i=1,5
! !    write(*,*)A(i,:)
! ! end do

! ! call permute_cols(A,permcol,.true.)
! ! write(*,*)'AFTER inverse col permutation'
! ! write(*,*)permcol
! ! do,i=1,5
! !    write(*,*)A(i,:)
! ! end do

! call permute_sym(A,permrow,.false.)
! write(*,*)'AFTER symmetric permutation'
! write(*,*)permrow
! do,i=1,5
!    write(*,*)A(i,:)
! end do

! call permute_sym(A,permrow,.true.)
! write(*,*)'AFTER inverse symmteric permutation'
! write(*,*)permrow
! do,i=1,5
!    write(*,*)A(i,:)
! end do








! vecch=(/'eins','zwei','drei','vier','funf'/)
! veci=(/1,2,3,4,5/)
! vecd=(/1.d0,2.d0,3.d0,4.d0,5.d0/)
! write(*,*)'vector permutations'
! write(*,*)'BEFORE double'
! write(*,*)permrow
! write(*,*)vecd
! call permute(vecd,permrow)
! write(*,*)'AFTER double'
! !write(*,*)permrow
! write(*,*)vecd

! call permute(vecd,permrow,.true.)
! write(*,*)'AFTER inverse permutation'
! !write(*,*)permrow
! write(*,*)vecd

! write(*,*)'BEFORE integer'
! write(*,*)permrow
! write(*,*)veci
! call permute(veci,permrow)
! write(*,*)'AFTER integer'
! !write(*,*)permrow
! write(*,*)veci
! call permute(veci,permrow,.true.)
! write(*,*)'AFTER inverse permutation '
! !write(*,*)permrow
! write(*,*)veci


! write(*,*)'BEFORE string'
! write(*,*)permrow
! write(*,*)vecch
! call permute(vecch,permrow)
! write(*,*)'AFTER string'
! !write(*,*)permrow
! write(*,*)vecch
! call permute(vecch,permrow,.true.)
! write(*,*)'AFTER inverse permutation'
! !write(*,*)permrow
! write(*,*)vecch

! end program test_permute

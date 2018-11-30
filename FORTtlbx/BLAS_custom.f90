!!this file contains fortran 90 wrappers of the BLAS routines, in ordwer to perform more efficient matrix multiplication then the native matmul (although some compilers already use BLAS internally)
!!For optimal speed it is advised to use the GotoBlAS routines, optimized for your own processor and using parallel methods

!!Coded by Roelof Rietbroek, Mon Oct 29 09:41:29 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!function which left multiplies a 2d matrix with a symmetric matrix (in packed form upper triangle given)
!! packmul2d=mpack*mat


function packmul2d(mpack,mat)
implicit none
double precision,intent(in),dimension(:)::mpack
double precision,intent(in),dimension(:,:)::mat
double precision,dimension(:,:)::packmul2d(size(mat,1),size(mat,2))
double precision::dumsq(size(mat,1),size(mat,1))
integer::row,col,ind,m,n
!get dimension sizes and perform a check
m=size(mat,1)
n=size(mat,2)
packmul2d=0.d0




!create square dummy array (only upper triangle)
dumsq=0.d0
ind=0
do,col=1,m
   do,row=1,col
      ind=ind+1
      dumsq(row,col)=mpack(ind) !upper triangle
!      dumsq(col,row)=mpack(ind)!lower triangle
   end do
end do

!write(*,*)m,n
!call BLAS routine

!!assign pointers
!matp=>mat
!packmulp=>packmul2d
!note the brackets of packmul2d(:,:) must be there (gives !
call DSYMM('L','U',m,n,1.d0,dumsq(:,:),m,mat(:,:),m,0.d0,packmul2d(:,:),m)
!call DSYMM(inside,'U',m,n,1.d0,dumsq,m,matp,m,0.d0,packmulp,m)

!note the brackets of packmul2d(:,:) must be there !!!
!only providing packmul2d gives a segmentation fault ???

end function packmul2d






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function which multiplies a symmetric packed matrix with a vector
! packmul1d=mpack*vec


function packmul1d(mpack,vec)
implicit none
double precision,intent(in),dimension(:)::mpack,vec
double precision,dimension(:)::packmul1d(size(vec,1))
!double precision,pointer,dimension(:)::packmulp,vecp,packp
integer::row,col,ind,m
!get dimension sizes and perform a check
m=size(vec,1)
packmul1d=0.d0

!call BLAS routine
!packp=>mpack
!vecp=>vec
!packmulp=>packmul1d
!note the brackets of packmul2d(:,:) must be there (gives !
!call DSPMV('U',m,1.d0,packp,vecp,1,0.d0,packmulp,1)
call DSPMV('U',m,1.d0,mpack(:),vec(:),1,0.d0,packmul1d(:),1)

end function packmul1d


!!function which replaces the generic matmul with a much quicker BLAS routine
!matmul2d=A*B
function matmul2d(A,B,TRANSA,TRANSB)
implicit none
double precision,intent(in),dimension(:,:)::A,B
integer,intent(in)::TRANSA,TRANSB !specify whether matrix should be transposed 1 is yes 0 is no
double precision::matmul2d(size(A,(1+TRANSA)),size(B,(2-TRANSB))) !automatic array
character(1)::TRA,TRB
integer::k,m,n,lda,ldb,ldc
TRA='N'
TRB='N'
if(TRANSA .eq. 1)TRA='T' !tranpose A
if(TRANSB .eq. 1)TRB='T'

!set indices
m=size(A,1+TRANSA)
n=size(B,2-TRANSB)
k=size(A,2-TRANSA)
lda=size(A,1)
ldb=size(B,1)
ldc=size(matmul2d,1)


!call BLAS routine

call DGEMM(TRA,TRB,m,n,k,1.d0,A,lda,B,ldb,0.d0,matmul2d(:,:),ldc)


end function matmul2d

!!function matmul2d without trasnpose arguments

function matmul2dclean(A,B)
  use forttlbx,only: matmul2d
  implicit none
  double precision, intent(in),dimension(:,:)::A,B
  double precision::matmul2dclean(size(A,1),size(B,2))
!call matmul2d with specific arguments
  matmul2dclean=matmul2d(A,B,0,0)
!  write(*,*)loc(A),loc(B)

end function matmul2dclean


!!function which replaces the generic matmul with a much quicker BLAS routine
!matmul1d=A*vec
function matmul1d(A,vec,TRANSA)
implicit none
double precision,intent(in),dimension(:,:)::A
double precision,intent(in),dimension(:)::vec
integer,intent(in)::TRANSA!specify whether matrix should be transposed 1 is yes 0 is no
double precision::matmul1d(size(A,(1+TRANSA))) !automatic array
character(1)::TRA
integer::m,n
TRA='N'
if(TRANSA .eq. 1)TRA='T' !tranpose A

!set indices
m=size(A,1)
n=size(A,2)

!call BLAS routine
call DGEMV(TRA,m,n,1.d0,A,m,vec,1,0.d0,matmul1d(:),1)

end function matmul1d

function matmul1dclean(A,B)
use forttlbx,only: matmul1d
  implicit none
  double precision, intent(in),dimension(:,:)::A
  double precision, intent(in),dimension(:)::B
  double precision::matmul1dclean(size(A,1))
!call matmul1d with specific arguments
  matmul1dclean=matmul1d(A,B,0)
end function matmul1dclean


!!function which substitues dot_product with a blasversion
function dot_productblas(A,B)
implicit none
double precision,intent(in)::A(:),B(size(A,1)) !B is an automatic array
double precision::dot_productblas,ddot

!write(*,*)shape(A),shape(B),size(A),size(B)
!!call BLAS routine
dot_productblas=ddot(size(A,1),A,1,B,1)

end function dot_productblas

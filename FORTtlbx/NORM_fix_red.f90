!!Updated by Roelof Rietbroek, Tue Jun  9 15:26:50 2009
!! changed arguments to assumed shape arrays and passed only scalars to lapack/BLAS routines



!! subroutines for changing apriori values of normal systems ,(fixing) and reducing normal systems


!! we start with the observation equation:
!! A * x =b + e
!! where x are unknowns
!! A is the design matrix
!! b are the obserrvations
!! and e is the associated observation noise with covaraince matrix P(-1)

!! To solve in a least square approach we can solve the following full NORMAL system:
!!  C*z=d
!! with C =A**(T) * P * A
!! z = x_est - x0
!! d= A**(T) P (b-A * x0)

!!where x0 are the apriori values of the unknowns


!! now we can divide the matrix in 4 blocks

      !full normal matrix
      !	  | C_11  C_12 |   | A_1**(T)* P * A_1     A_1**(T)* P * A_2 |   
      ! C=|            | = |                                         |
      !   | C_21  C_22 |   | A_2**(T)* P * A_1     A_2**(T)* P * A_2 |   


      !full right hand side
      !                               | d_1 |   |  A_1**T P * (b- A_1 * x0_1) -  C_12 * x0_2 !
      ! d= A**T P ( b - A x0) =       |     | = |
      !                               | d_2 |   |  A_2**T P * (b- A_1 * x0_1) -  C_22 * x0_2 !


      !Apriori residual
      ! ltpl= (b - A x0)**T * P * (b-A x0) = (b - A_1 * x0_1)**T * P * (b - A_1 * x0_1)
      !                                      - 2 *x0_2**(T)*A_2**(T) * P * (b - A_1 * x0_1 ) 
      !                                      + x0_2**T C_22 * x0_2






!!!!!!!!!!!!!!!!!!!!!!!!!!!FIXING PARAMETERS TO THEIR APRIORI VALUES!!!!!!!!!!!!!!!!

!!We want to obtain a sub normal system:
!! consider the observation equation subsystem:

!! A_2 z_2 = b

!!the associated normal system becomes

!! Cnew_22 * z_2 = dnew= A_2**(T) P (b -A_2 * x0_2 -A_1 * x0_1) = d_2

!! ltplnew=(b - A_2 x0_2 -A_1 * x0_1)**T * P * (b-A_2 x0_2 -A_1 * x0_1) = ltpl

!! Cnew = C_22 (this is easy )

!! iin other words nothing really changed: we only need to work with subsections of the original 
!! normal equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!CHANGING APRIORI VALUES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!consider the unblocked normal system :

! C x = d 
! C = A**(T) P A
! d= A**(T) P (b- A x0)
! and
! ltpl= (b- A x0) **(T) P (b- A x0)


!imagine the case that we want to update the apriori vector x0 with x1

!x1=x0+dx

!what are the new  Cnew, dnew, ltplnew ?

!Cnew= C (independent of the apriori info and thus unchanged)

!  step 1

!do the update of ltpl in two steps

! first update
! ltpl= ltpl - dx **(T) C dx

! calculate the new right hand side form the old right hand side

! d=d - C * dx


!second update of the ltpl vector
! ltpl= ltpl - dx **(T) d
!where d has changed in the new right hand side now

!done!
 


subroutine setapriori_norm(C,d,dx,ltpl,uplo)
implicit none
character(1),intent(in),optional::uplo
double precision,intent(in)::C(:,:),dx(:)
double precision,intent(inout)::d(:),ltpl
double precision::ddot
integer::n,ildc
character(1)::iuplo
iuplo='U'
if(present(uplo))iuplo=uplo

!write(*,*)iuplo

n=size(dx,1)

!first update of ltpl
ltpl=ltpl-ddot(n,d,1,dx,1)

!calculate new rigth hand side


call dsymv(iuplo,n,-1.d0,C,n,dx,1,1.d0,d,1)


!second update of ltpl (d has changed now)
ltpl=ltpl-ddot(n,d,1,dx,1)

end subroutine setapriori_norm


!!VERSION below uses less memory due to avoidance of implicit array copying by compiler
!!HERE C must be a contiguous set
subroutine setapriori_norm2(C,ldc,d,dx,ltpl,uplo)
implicit none
character(1),intent(in),optional::uplo
integer,intent(in)::ldc
double precision,intent(in)::dx(:),C(:,:)
double precision,intent(inout)::d(:),ltpl
double precision::ddot
integer::n
character(1)::iuplo
iuplo='U'
if(present(uplo))iuplo=uplo

!write(*,*)iuplo

n=size(dx,1)

!first update of ltpl
ltpl=ltpl-ddot(n,d,1,dx,1)

!calculate new rigth hand side


call dsymv(iuplo,n,-1.d0,C(1,1),ldc,dx,1,1.d0,d,1)


!second update of ltpl (d has changed now)
ltpl=ltpl-ddot(n,d,1,dx,1)



end subroutine setapriori_norm2



!!!!!!!!!!!!!!!!!!!!REDUCING SYSTEMS!!!!!!!!!!!!!!!!!!!!!

!!again consider the blocked system
    !full normal matrix
      !	  | C_11  C_12 |   | A_1**(T)* P * A_1     A_1**(T)* P * A_2 |   
      ! C=|            | = |                                         |
      !   | C_21  C_22 |   | A_2**(T)* P * A_1     A_2**(T)* P * A_2 |   


      !full right hand side
      !                               | d_1 |   |  A_1**T P * (b- A_1 * x0_1) -  C_12 * x0_2 !
      ! d= A**T P ( b - A x0) =       |     | = |
      !                               | d_2 |   |  A_2**T P * (b- A_1 * x0_1) -  C_22 * x0_2 !


      !Apriori residual
      ! ltpl= (b - A x0)**T * P * (b-A x0) = (b - A_1 * x0_1)**T * P * (b - A_1 * x0_1)
      !                                      - 2 *x0_2**(T)*A_2**(T) * P * (b - A_1 * x0_1 ) 
      !                                      + x0_2**T C_22 * x0_2


!the inverse of such a block matrix may be written as:

!           |    C_11**(-1)+ C_11**(-1) * C_12 * C_11**(-1) * C_21 * C_11**(-1)               -C_11**(-1) * C_12 * C_22red**(-1)  |   
!C**(-1)=
!           |             - C_22red**(-1)  * C_21 * C_11**(-1)                                          C_22red**(-1)            |   

!where C_22red = C_22 - C21 * C_11**(-1) * C_12

! If we are not interested in the estimates z1 we can still estimate the z2 components by solving the normal system:


! C_22red z2 = d2 - C21 * C_11**(-1) d1

!in other words
! our new C_22 (reduce matrix) becomes C_22 - C21 * C_11**(-1) * C_12
!and our new right hand side:

!d2 becomes d2 -  C21 * C_11**(-1) d1

!and our new apriori estimate becomes 

!ltpl becomes ltpl -   d1**(T)* C_11**(-1) d1
!the steps are as following




subroutine reduce_norm(C11,C12,C22,d1,d2,ltpl)
implicit none
double precision,intent(inout)::C11(:,:),C12(:,:),C22(:,:),d1(:),d2(:),ltpl
double precision::ddot
integer::n1,n2,info,stderr
stderr=0
n1=size(C12,1)
n2=size(C12,2)

!step 1: calculate the cholesky decomposition of C11

!C_11 will be redined as U_11

call DPOTRF('U', n1, C11, n1, info )
if(info .ne. 0)then
   write(stderr,*)' ERROR: MATRIX reduction failed, cannot invert local parameter matrix'
   write(stderr,*)'        Cholesky factorization broke down on parameter:',info
   stop
end if

!thus if this test is not passed the reduction was not possible


!step 2 solve the triangular system U_11**(T) * X = C12 for X
!note C_12 will be redined by X

call DTRSM('L','U','T','N',n1,n2,1.d0,C11,n1,C12,n1)


!step 3 calculate the reduced normal matrix C22red= C22 - X**(T) * X = C22 - C21 * U_11 **(-1) * U_11**(-T) C12

call DSYRK('U','T',n2,n1,-1.d0,C12,n1,1.d0,C22,n2)


!step 4 calculate
!! d1new = U**(-T)* d1 !d1 will be updated
!DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
call dtrsv('U','T','N',n1,C11,n1,d1,1)


!step 5 calculate new right hand side
! d2new=d2 -     C12**(T)  U_11**(-1)  U_11**(-T) d1


!DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
call dgemv('T',n1,n2,-1.d0,C12,n1,d1,1,1.d0,d2,1)


!step 6 calculate the new ltpl

ltpl=ltpl-ddot(n1,d1,1,d1,1)




end subroutine reduce_norm

subroutine reduce_norm2(C11,ldc11,C12,ldc12,C22,ldc22,d1,d2,ltpl)
implicit none
integer,intent(in)::ldc11,ldc12,ldc22
double precision,intent(inout)::C11(:,:),C12(:,:),C22(:,:),d1(:),d2(:),ltpl
double precision::ddot
integer::n1,n2,info,stderr
stderr=0
n1=size(d1)
n2=size(d2)

!step 1: calculate the cholesky decomposition of C11

!C_11 will be redined as U_11

call DPOTRF('U', n1, C11(1,1), ldc11, info )
if(info .ne. 0)then
   write(stderr,*)' ERROR: MATRIX reduction failed, cannot invert local parameter matrix'
   write(stderr,*)'        Cholesky factorization broke down on parameter:',info
   stop
end if

!thus if this test is not passed the reduction was not possible


!step 2 solve the triangular system U_11**(T) * X = C12 for X
!note C_12 will be redined by X

call DTRSM('L','U','T','N',n1,n2,1.d0,C11(1,1),ldc11,C12(1,1),ldc12)


!step 3 calculate the reduced normal matrix C22red= C22 - X**(T) * X = C22 - C21 * U_11 **(-1) * U_11**(-T) C12

call DSYRK('U','T',n2,n1,-1.d0,C12(1,1),ldc12,1.d0,C22(1,1),ldc22)


!step 4 calculate
!! d1new = U**(-T)* d1 !d1 will be updated
!DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
call dtrsv('U','T','N',n1,C11(1,1),ldc11,d1,1)


!step 5 calculate new right hand side
! d2new=d2 -     C12**(T)  U_11**(-1)  U_11**(-T) d1


!DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
call dgemv('T',n1,n2,-1.d0,C12(1,1),ldc12,d1,1,1.d0,d2,1)


!step 6 calculate the new ltpl

ltpl=ltpl-ddot(n1,d1,1,d1,1)


end subroutine reduce_norm2



subroutine solve_norm(C,d,ltpl)
implicit none
double precision,intent(inout)::C(:,:),d(:),ltpl
integer::n1,info,stderr
double precision::ddot
stderr=0
n1=size(d,1)

!Cholesky decomposition C = U_c **(T) U_c
!DPOTRF( UPLO, N, A, LDA, INFO )
call dpotrf('U',n1,C,n1,info)

if(info .ne. 0)then
   write(stderr,*)"ERROR: cannot invert the system"
   write(stderr,*)"ERROR: Cholesky decomposition broke down on parameter",info
   stop
end if

!solve triangular system  (U_c**(T) x= d for x (d is updated) 

!DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

call dtrsv('U','T','N',n1,C,n1,d,1)


!update apriori residual norm

ltpl=ltpl-ddot(n1,d,1,d,1)

!solve triangular system U_c xnew= xold

call dtrsv('U','N','N',n1,C,n1,d,1) !d is updated again



!solve inverse of the matrix
!DPOTRI( UPLO, N, A, LDA, INFO )
call dpotri('U',n1,C,n1,info)

if(info .ne. 0)then
   write(stderr,*)"ERROR: cannot invert the system"
   write(stderr,*)"ERROR: Inversion of Cholesky decomposition broke down on parameter",info
   stop
end if

 
end subroutine solve_norm

subroutine solve_norm2(C,ldc,n,d,ltpl)
implicit none
integer,intent(in)::ldc,n
double precision,intent(inout)::C(:,:),d(n),ltpl
integer::info,stderr
double precision::ddot
stderr=0

!Cholesky decomposition C = U_c **(T) U_c
!DPOTRF( UPLO, N, A, LDA, INFO )
call dpotrf('U',n,C(1,1),ldc,info)


if(info .ne. 0)then
   write(stderr,*)"ERROR: cannot invert the system"
   write(stderr,*)"ERROR: Cholesky decomposition broke down on parameter",info
   stop
end if

!solve triangular system  (U_c**(T) x= d for x (d is updated) 

!DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

call dtrsv('U','T','N',n,C(1,1),ldc,d,1)


!update apriori residual norm
!write(*,*)ltpl
ltpl=ltpl-ddot(n,d,1,d,1)
!write(*,*)ltpl
!solve triangular system U_c xnew= xold

call dtrsv('U','N','N',n,C(1,1),ldc,d,1) !d is updated again



!solve inverse of the matrix
!DPOTRI( UPLO, N, A, LDA, INFO )
call dpotri('U',n,C(1,1),ldc,info)

if(info .ne. 0)then
   write(stderr,*)"ERROR: cannot invert the system"
   write(stderr,*)"ERROR: Inversion of Cholesky decomposition broke down on parameter",info
   stop
end if

 
end subroutine solve_norm2

! !solve normal equation system using a truncated singular value decomposition
subroutine solve_tsvd(C,ldc,n,d,ltpl,ntrunc)
implicit none
integer,intent(in)::ldc,n
double precision,intent(inout)::C(:,:),d(n),ltpl
integer,intent(in)::ntrunc !number of eigenvalues to keep
integer::info,stderr
double precision::ddot,eigvals(n)
double precision,allocatable::workv(:)
integer::lwork,i,j,ilaenv
double precision::tmp(n,ntrunc)
stderr=0

if(ntrunc > n .or. ntrunc <=0)then
   write(stderr,*)"ERROR, requested truncation is erroneous: ",ntrunc
   stop
end if

!get eigen values and eigen vectors
!get optimal blocksize
lwork=ILAENV( 1, 'DSYTRD','U', n, -1, -1, -1 )
lwork = max(1,(lwork+2)*n)
allocate(workv(lwork))

!call lapack 
call dsyev('V','U',n,C(1,1),ldc,eigvals(1),workv(1),lwork,info)
if(info <0)then
   write(stderr,*)"ERROR, illegal value occured at: ",info
   stop
else if (info >0)then
   write(stderr,*)"ERROR, failed to converge"
   stop
end if


!rescale the eigenvectors by the inverse root of the (singular) eigenvalues
do,i=n-ntrunc+1,n !loop over significant eigenvalues
   do,j=1,n
      tmp(j,i-n+ntrunc)=C(j,i)/sqrt(eigvals(i))
   end do
end do
!use eigvals vector to temporary store the right hand side
eigvals(1:n)=d(1:n)



!compute a rank -k update to compute the generalized inverse (C will be overwritten)
call DSYRK('U','N',n,ntrunc,1.d0,tmp(1,1),n,0.d0,C(1,1),ldc)

!apply matrix vector multiplications to obtain the solution vector
call dsymv('U',n,1.d0,C(1,1),ldc,eigvals(1),1,0.d0,d(1),1)

!update apriori residual norm
!write(*,*)ltpl
ltpl=ltpl-ddot(n,eigvals,1,d,1)
 
end subroutine solve_tsvd


subroutine diag_transf(B,A,atb,x0,nunit,nb,add)
implicit none
integer,intent(in)::nb,nunit
double precision,intent(in)::B(:)
double precision,intent(inout),optional::atb(:),x0(:),A(:,:)
logical,intent(in),optional::add
integer::j,k,shift,stderr
logical::iadd
iadd=.false.
stderr=0
if(present(add))iadd=add

if(present(x0) .and. iadd)then
   write(stderr,*)"ERROR (diag_transf): Transformation matrix is not invertible since parameters are added!!"
   stop
end if

if(present(atb))then!simple forward propagation
   if(iadd)then
      shift=nunit+nb
   else
      shift=nunit
   end if
   
   forall(j=1:nb)atb(shift+j)=atb(nunit+j)*B(j)
!   forall(j=nunit+1:nunit+nb)atb(j)=B(j-nunit)*atb(j)
end if

if(present(A))then ! matrix propagation
   if(iadd)then
      shift=nunit+nb
   else
      shift=nunit
   end if

   !!   N_ub B**T
   forall(j=1:shift,k=1:nb)A(j,k+shift)=&
        A(j,k+nunit)*B(k)

   !!            B N_bb B**T
   forall(j=1:nb,k=1:nb, k>=j)A(j+shift,k+shift)=&
        A(j+nunit,k+nunit)*B(j)*B(k)


!    !!   N_ub B**T
!    forall(j=1:nunit,k=nunit+1:nunit+nb)A(j,k)=&
!         A(j,k)*B(k-nunit)

   !!            B N_bb B**T
!    forall(j=nunit+1:nunit+nb,k=nunit+1:nunit+nb, k>=j)A(j,k)=&
!        A(j,k)*B(j-nunit)*B(k-nunit)
end if

if(present(x0))then !inverse propagation B must have no zero entries!
   shift=nunit
   forall(j=1:nb)x0(j+shift)=x0(j+shift)/B(j)
end if

end subroutine diag_transf

!subroutine to transform the error covariance matrix of the underlying solution system
!The apriori cost functional will be adapted to this new situation. The information content of the earlier normal system asembly will be lost however.
!! the input matrix A is assumed to be in the following form:
!! U U U
!! L U U
!! L L L
!! it contains twice the normal matrix
subroutine trans_cov(B,ldb,A,lda,atb,ltpl,nunit,nb,type)
use forttlbx,only:diag_transf,sym_transf,full_transf
implicit none
integer,intent(in)::nunit,nb,lda,ldb,type
double precision,intent(inout)::B(:,:) 
double precision,intent(inout)::A(:,:) !array lda >= (nunit+nb)
double precision,intent(inout)::ltpl,atb(:)
integer::i,j,info
double precision,allocatable::tmp(:)
double precision::ddot

allocate(tmp(nunit+nb))
!Solve for x:= A x= atb ! use the lower part of A ( this will be updated with the Cholesky decomposition)
call dposv('L',nunit+nb,1,A(2,1),lda,atb,nunit+nb,info)
if(info .ne. 0)then
   write(*,*)"ERROR(diag_transerr): Cholesky decomposition failed"
   write(*,*)"ERROR(diag_transerr): Normal system is singular"
   stop
end if

!propagate matrix A  with B
select case(type)
case(1)! B is diagonal 
   call diag_transf(B=B(:,1),A=A,nunit=nunit,nb=nb)
case(2) !B is a symmetric matrix
   call sym_transf(B=B,ldb=ldb,A=A,lda=lda,nunit=nunit,nb=nb)
case(3) !B is a full (square) matrix
   call full_transf(B=B,ldb=ldb,A=A,lda=lda,nunit=nunit,nb=nb,k=nb)
end select


!calculate new atb (in tmp)
call dsymv('U',nunit+nb,1.d0,A(1,1),lda,atb,1,0.d0,tmp,1)

!renew ltpl
ltpl=ddot(nunit+nb,tmp,1,atb,1)

!copy values back in atb
forall(i=1:nunit+nb)atb(i)=tmp(i)
deallocate(tmp)
end subroutine trans_cov


!subroutine to transform normal system matrix and vectors
!values up to nunit are transformed through the unit matrix
! values from nunit+1 until nunit+nb are transformed through B
!Transformation which is performed
!      |  A_uu      A_ub B**T  |   
!  A:= 
!      | ignored   B A_bb B**T |
!
!       |  atb_u  |
! atb:=
!       | B atb_b |
!
!
!       |  x0_u  |
! x0:= 
!       | B**-1 x0_b |
!
!B is a symmetric matrix and contains on output the  Cholesky decomposition of B if x0 is given on input
!A must be large enough to contain the updated matrix lda >= (nb+nunit)

!!Updated by Roelof Rietbroek, Wed Apr 14 15:40:58 2010
!! Allow addition of parameters
!! if the 'add' option is given and true then
! A must have a size lda> max(nb+nunit+nb)
!      |  A_uu      A_ub     A_ub B**T   |   
!  A:= 
!      |  A_bu      A_bb     A_bb B**T   |
!      
!      |  B A_bu   B A_bb    B A_bb B**T |
!
!       |  atb_u  |
! atb:=
!       |  atb_b  |
! 
!       | B atb_b |
!
!
!       
! x0:=  Cannot be determined!!
!       

 

subroutine sym_transf(B,ldb,A,lda,atb,x0,nunit,nb,add)
implicit none
integer,intent(in)::ldb,lda,nb,nunit
double precision,intent(inout)::B(:,:)
double precision,intent(inout),optional::atb(:),x0(:),A(:,:)
logical,intent(in),optional::add
double precision,allocatable::C(:) ! work array only necessary in the case of A
double precision::tmp(nb)
integer::info,stderr,i,j,ldc,shift
logical::iadd
stderr=0
iadd=.false.

if(present(add))iadd=add

!error checks
if(present(x0) .and. iadd)then
   write(stderr,*)"ERROR (sym_transf): Transformation matrix is not invertible since it is rectangular!!"
   write(stderr,*)"ERROR (sym_transf): Not able to propagate back the apriori vector "
   stop
   
end if

if(iadd .and. lda < nunit+2*nb .and. present(A))then
   write(stderr,*)"ERROR (sym_transf): Matrix must be at least",nunit+nb*2
   stop
else if(.not. iadd .and. lda < nunit+nb .and. present(A))then
   write(stderr,*)"ERROR (sym_transf): Matrix must be at least",nunit+nb
   stop
end if


if(present(atb))then!simple propagation (B is unchanged)
   call dsymv('U',nb,1.d0,B(1,1),ldb,atb(nunit+1),0.d0,tmp,1)

   if(iadd)then
      shift=nunit+nb
   else
      shift=nunit
   end if

   atb(shift+1:shift+nb)=tmp         
end if

if(present(A))then ! matrix propagation only if A is present
   if(iadd)then
      allocate(C(nb*(nunit+nb)))!allocate a matrix which is large enough to hold a (nb+nunit) x nb matrix
   else
      allocate(C(nb*max(nunit,nb)))!allocate a matrix which is large enough to hold both a nunit x nb matrix as well as a nb x nb matrix 
   end if
   !calculate A_ub B**T = A_ub B ( use symmetry of B)

   if(nunit>0 .or. iadd)then
      if(iadd)then
         ldc=nunit+nb
      else
         ldc=nunit  !trick dsymm into thinking C is a 2D array: C(nunit,*)
      end if
      
      call dsymm('R','U',ldc,nb,1.d0,B(1,1),ldb,A(1,nunit+1),lda,0.d0,C(1),ldc)
      if(iadd)then
         shift=nunit+nb
      else
         shift=nunit
      end if

      !copy values back in original array
      forall(i=1:ldc,j=1:nb)A(i,j+nunit)=C(i+(j-1)*ldc)

   end if

   !calculate B A_bb B
   !by means of a rank-2 update
   !split up A_bb in two triangular matrices
   !A_bb = L_bb+ U_bb =  L_bb + L_bb**T (where L_bb is the lower triangle of A_bb with a diagonal which is half of the diagonal of A_bb
   forall(i=nunit+1:nunit+nb)A(i,i)=A(i,i)/2.d0 ! divide the diagonal by two

   !triangular matrix matrix multiplication( B is assumed to be full)
   !C=B*U_bb
   !first copy values of B in C
   ldc=nb  !trick dtrmm and dsyr2k into thinking C is a 2D array: C(nb,*)
   do,i=1,nb! loop over columns
      shift=(i-1)*nb ! shift depending on the column number
      do,j=1,nb! loop over rows
         C(shift+j)=B(i,j)
      end do
   end do

   !C=B*U_bb
   call dtrmm('R','U','N','N',nb,nb,1.d0,A(nunit+1,nunit+1),lda,C(1),ldc)
   
   !now calculate symmetric rank 2 operation of the matrix A_bb
   !A_bb:= B A_bb B**T =  C B**T + B C**T

   if(iadd)then
      shift=nunit+nb
      !also restore original diagonal of A
      forall(i=nunit+1:nunit+nb)A(i,i)=A(i,i)*2.d0 !multiply the diagonal by two
   else
      shift=nunit
   end if   

   call dsyr2k('U','N',nb,nb,1.d0,C(1),ldc,B(1,1),ldb,0.d0,A(shift+1,shift+1),lda)

   deallocate(C)!free up memory
end if

if(present(x0))then
   call dposv('U',nb,1,B(1,1),ldb,x0(nunit+1),size(x0,1),info)
!!!B now  contains the Cholesky factorization of B
   if(info .ne. 0)then
      write(*,*)"ERROR(sym_transf): Cholesky decomposition failed"
      write(*,*)"ERROR(sym_transf): B may be positive SEMI definite"
      stop
   end if
end if

end subroutine sym_transf


!! transformation by a non-symmetric matrix B
!! NOTE the sides and lengths of A and x0 and atb must be > max(nb+nunit,k)

!!Updated by Roelof Rietbroek, Wed Apr 14 15:15:31 2010
!! also allowed addign of the parameter space ( introducing rankdefects)

subroutine full_transf(B,ldb,A,lda,atb,x0,nunit,nb,k,add)
implicit none
integer,intent(in)::ldb,lda,nunit,nb,k
double precision,intent(inout)::B(:,:)
double precision,intent(inout),optional::atb(:),x0(:),A(:,:)
logical,optional,intent(in)::add !tag determining whether parameters will be added
double precision,allocatable::C(:)

double precision::tmp(k)
integer::info,ipiv,stderr,i,j,ldc,shift
logical::iadd
stderr=0

iadd=.false.
if(present(add))iadd=add


!error check
if(present(x0))then
   if(nb .ne. k .or. iadd)then
      write(stderr,*)"ERROR (full_transf): Transformation matrix is not invertible since it is rectangular!!"
      write(stderr,*)"ERROR (full_transf): Not able to propagate back the apriori vector "
      stop
   end if
end if


if(lda < nunit+max(nb,k) .and. present(A))then
   write(stderr,*)"ERROR (full_transf): Matrix must be at least",nunit+max(nb,k),lda
   stop
end if

!DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
if(present(atb))then!simple propagation
   ! do,i=1,3
   !    write(0,*)B(i,:),atb(nunit+i)
   ! end do

    call dgemv('N',k,nb,1.d0,B(1,1),ldb,atb(nunit+1),1,0.d0,tmp,1)

    if(iadd)then
       shift=nunit+nb
    else
       shift=nunit
    end if
    
    atb(shift+1:shift+k)=tmp(1:k)         

end if


if(present(A))then ! matrix propagation only if A is present

   allocate(C(k*max(nunit,nb)))! allocate enough memory to store a k x nb or nunit x k matrix (packed in a vector)
   

   !calculate A_ub B**T if needed
   if(nunit>0)then
      !DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      ldc=nunit
      call dgemm('N','T',nunit,k,nb,1.d0,A(1,nunit+1),lda,B(1,1),ldb,0.d0,C(1),ldc)

      if(iadd)then
         shift=nunit+nb
      else
         shift=nunit
      end if
      !copy values back in original array
      forall(i=1:ldc,j=1:k)A(i,j+shift)=C(i+(j-1)*ldc)

   end if


   !calculate A_bb B**T
   !in 2 steps ( per triangle)
   
   !divide up A_bb in two triangles by dividing the diagonal  by 2
   !A_bb = L_bb+ U_bb =  L_bb + L_bb**T (where L_bb is the lower triangle of A_bb with a diagonal which is half of the diagonal of A_bb
   forall(i=nunit+1:nunit+nb)A(i,i)=A(i,i)/2.d0 ! divide the diagonal by two

   !first copy values of B in C
   ldc=k  !trick dtrmm and dsyr2k into thinking C is a 2D array: C(k,*)
   do,i=1,nb! loop over columns
      shift=(i-1)*k ! shift depending on the column number
      do,j=1,k! loop over rows
         C(shift+j)=B(j,i)
      end do
   end do


   !C=B*U_bb (triangular matrix multiplication, C is overwritten)
   call dtrmm('R','U','N','N',k,nb,1.d0,A(nunit+1,nunit+1),lda,C(1),ldc)


   !calculate B A_bb B
   !by means of a rank-2 update

   if(iadd)then
      shift=nunit+nb
   else
      shift=nunit
   end if
  !  WRITE(0,*)"Entering dsyr2k"
  ! write(0,*)'U','N',k,nb,1.d0,C(1),ldc,B(1,1),ldb,0.d0,A(shift+1,shift+1),lda
   call dsyr2k('U','N',k,nb,1.d0,C(1),ldc,B(1,1),ldb,0.d0,A(shift+1,shift+1),lda)


   if(iadd)then ! we need an additional triangular matrix mulitplication to get A_bb B**T
      !copy values of B in reklevant block of A (will be overwritten)

      shift=nunit+nb
      do,i=1,nb! loop over columns of B
!         shift=(i-1)*k ! shift depending on the column number
         do,j=1,k! loop over rows of B
            A(nunit+i,j+shift)=B(j,i)
         end do
      end do
      
      !triangular matrix left multiplication
      
      !A=U_bb B**T ( Note dtrmm operates on 2 different subsections of A)
      call dtrmm('L','U','N','N',nb,k,1.d0,A(nunit+1,nunit+1),lda,A(nunit+1,shift+1),lda)
      
      !add the other triangular matrix multiplication still contained in C**T
      do,i=1,k
         do,j=1,nb
            A(nunit+j,shift+i)=A(nunit+j,shift+i)+C(i+(j-1)*k)
         end do
      end do
      
      !also restore original diagonal of A_bb section
      forall(i=nunit+1:nunit+nb)A(i,i)=A(i,i)*2.d0 !multiply the diagonal by two

   end if

   
   deallocate(C)!free up memory

end if

if(present(x0))then
   ! try solving Bx=x0 for x
   call dgesv(nb,1,B(1,1),ldb,ipiv,x0(nunit+1),size(x0,1),info)
   if(info .ne. 0)then !issue message and stop program if B is singular
      write(stderr,*)"ERROR(full_transf): Transformation matrix is not invertible"
      write(stderr,*)"ERROR (full_transf): Not able to propagate back the apriori vector "
      stop
   end if
   
end if

end subroutine full_transf



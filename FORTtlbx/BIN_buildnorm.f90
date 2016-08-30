!program to built up a normal system from a observation set with (diagonal errors)
!!Coded by Roelof Rietbroek, Thu Oct 29 11:12:24 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Wed Jun  9 22:01:38 2010

!!Updated by Roelof Rietbroek, Thu Aug 26 12:06:10 2010
!! also allow observations not present in the designmatrix to be skipped

!!Updated by Roelof Rietbroek, Wed Jun 29 11:00:01 2011
!!ficed bug in nametag Nunknows -> Nunknowns

!! updated by Roelof Rietbroek Thursday 6 October 2011
!! Allow full error covariance matrices to be propagated

!!Updated by Roelof Rietbroek, Thu Sep  6 13:45:23 2012
!! also allow symmetric matrices as designmatrices

!!Updated by Roelof Rietbroek, Mon Aug 12 09:50:42 2013
!! added option to select different vectors as data



program BIN_buildnorm
use forttlbx
use binfiletools
use bin_operations
implicit none
type(BINdat)::obs,A,out,reusemat
character(200)::dum
integer::i,j,itharg,iargc,narg,stderr
integer::nunk,nobs,s_st,s_end,st,nd,n2
double precision::ltpl,cov_scale
integer,allocatable::perm1(:),perm2(:)
logical::trans,reuse,verbose
character(24),pointer,dimension(:)::side1,side2
integer::info,datcol

!defaults/ initializations
obs%file='stdin'
obs%mtyp='P'
trans=.false.
reuse=.false.
out%file='stdout'
cov_scale=1.d0
s_st=1
s_end=24
stderr=0
verbose=.false.
datcol=1
!!process command line options
narg=iargc()
if(narg <1)call help()

itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('A')!provides designmatrix name
         if(dum(4:4) .eq. '')then
            itharg=itharg+1
            call getarg(itharg,A%file)
         else

            A%file=trim(dum(4:))
         end if
      case('S')!limit global search area of the parameter strings
         st=4
         nd=index(dum(st:),'/')+st-2
         read(dum(st:nd),*)s_st
         read(dum(nd+2:),*)s_end
      case('t')
         trans=.true.
      case('R')!reuse normal matrix from file
         reuse=.true.
         if(dum(4:4) .eq. '')then
            itharg=itharg+1
            call getarg(itharg,reusemat%file)
         else
            reusemat%file=trim(dum(4:))
         end if
      case('V')!vector select
         read(dum(4:),*)datcol
      case('v')!be verbose
         verbose=.true.
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option but a filename
      obs%file=trim(dum)
   end if
end do

!some input checks

!read observation data (meta data only)
if(verbose)write(stderr,*)"Reading observations",obs%file
call read_BINtype(obs,1)

!errror check for no data

if(obs%nval1 .eq. 0)then
   write(stderr,*)'ERROR: No data provided, aborting ..'
   stop
end if

!check its type
select case(obs%type)
case('DIAVN___')!diagonal type (accepted)
   call read_BINtype(obs) ! read remainder of the data
case('SYMVN__')!accept full covariance matrix
   obs%mtyp='U'
   call read_BINtype(obs) ! read remainder of the data
case default
   write(stderr,*)"ERROR: Observation input is not accepted: ",obs%type
   stop
end select

!read design matrix
call read_BINtype(A,1)
if(verbose)write(stderr,*)"Reading design matrix",A%file
!check its type
select case(A%type)
case('FULLSQV0','FULLSQVN','FULL2DVN') !full matrix accepted
   A%mtyp='F' ! read 2d matrix
   call read_BINtype(A) ! read remainder of the data
case('SYMVN___')!allowed but needs to be expanded first
   A%mtyp='U'
   call read_BINtype(A) ! read remainder of the data
   !mirror matrix
   do,j=1,A%nval1
      do,i=1,j-1
         A%mat1(j,i)=A%mat1(i,j)
      end do
   end do
   
   !change type of the matrix
   A%type='FULL2DVN'
   A%mtyp='F'
   A%pval1=A%nval1**2
   
   ! allocate separate memory for column description
   A%side2_d=>null()
   allocate(A%side2_d(A%nval1))
   A%side2_d=A%side1_d !copy values

case default
   write(stderr,*)"ERROR: Transformation input type is not accepted",A%type
   stop
end select

if(trans)then
   nunk=A%nval1
   n2=A%nval2
   side1=>A%side2_d
   side2=>A%side1_d
else
   nunk=A%nval2
   n2=A%nval1
   side1=>A%side1_d
   side2=>A%side2_d
end if

!set up side permutation of the design matrix observations
allocate(perm2(n2),perm1(obs%nval1))
if(verbose)write(stderr,*)"Calculating permutation(s)"
!get permutation vectors
call get_permvec2(obs%side1_d,side1,.true.,s_st,s_end,perm1,perm2,nobs)

!call get_permvec(side1,obs%side1_d,.true.,s_st,s_end,perm,nobs)

! if(nobs .ne. obs%nval1)then
!    write(stderr,*)"ERROR: Design matrix is not complete",nobs,obs%nval1,perm(A%nval2),A%side1_d(perm(A%nval2))
!    stop
! end if

if(nobs .eq. 0 .or. nobs > obs%nval1)then
   write(stderr,*)"ERROR: no matching observations found or sorting error"
   stop
end if


!permute observations matrix (puts non-matching entries at the back)
if(verbose)write(stderr,*)"Permuting observations, ignoring ",obs%nval1-nobs," observations"

call BIN_permute(input=obs,perm1=perm1,both=.true.)


!permute design matrix if necessary ( puts non-matching entries at the back)
if(nobs .ne. n2)then
   if(verbose)write(stderr,*)"Permuting design matrix, ignoring ",n2-nobs," columns"
   if(trans)then
      call BIN_permute(input=A,perm2=perm2)
   else
      call BIN_permute(input=A,perm1=perm2)
   end if
end if

!construct normal matrix and vector

select case(obs%type)
case('DIAVN___')!diagonal type: scale design matrix
   !first take the root of the variances and its inverse
   obs%pack1(1:nobs)=cov_scale/sqrt(obs%pack1(1:nobs))
   
   !scale observation vector
   obs%vec(1:nobs,datcol)=obs%vec(1:nobs,datcol)*obs%pack1(1:nobs)

   !initial ltpl (cost functional)
   ltpl=dot_productblas(obs%vec(1:nobs,datcol),obs%vec(1:nobs,datcol))

   !scale design matrix
   if(trans)then
      do,j=1,nobs
         do,i=1,nunk
            A%mat1(i,j)=A%mat1(i,j)*obs%pack1(j)
         end do
      end do
      !   forall(j=1:nobs,i=1:nunk)A%mat1(i,j)=A%mat1(i,j)*obs%pack1(j)
   else
      do,i=1,nunk
         do,j=1,nobs
            A%mat1(j,i)=A%mat1(j,i)*obs%pack1(j)
         end do
      end do
      !   forall(j=1:nobs,i=1:nunk)A%mat1(j,i)=A%mat1(j,i)*obs%pack1(j)
   end if

case('SYMVN__')!full covariance matrix:apply cholesky decomposition and multiply design matrix
   !calculate Cholesky factorization of the error covariance
   if(verbose)write(stderr,*)"Calculating Cholesky factorization"
   call DPOTRF('U',nobs,obs%mat1(1,1),obs%ldm,info)
   if(info .ne. 0)then
      write(stderr,*)"Cholesky factorization failed on index",info
      stop
   end if

   if(trans)then !calculate A**T U**-T (update A**T)
      call DTRSM('R','U','T','N',nunk,nobs,1.d0,obs%mat1(1,1),obs%ldm,A%mat1(1,1),A%ldm)
   else !calculate U**-1 A (update A)
      call DTRSM('L','U','N','N',nobs,nunk,1.d0,obs%mat1(1,1),obs%ldm,A%mat1(1,1),A%ldm)
   end if

end select
   

!create right hand side vector
if(verbose)write(stderr,*)"Creating right hand side"
allocate(out%vec(nunk,2))
out%vec(:,2)=0.d0 ! set apriori vector to zero

!matrix vector multiplication
!DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

if(trans)then
   call dgemv('N',nunk,nobs,1.d0,A%mat1(1,1),&
        int(A%ldm),obs%vec(1,1),1,0.d0&
        ,out%vec(1,1),1)
else

   call dgemv('T',nobs,nunk,1.d0,A%mat1(1,1),&
        int(A%ldm),obs%vec(1,1),1,0.d0&
        ,out%vec(1,1),1)
end if



if(reuse)then ! resuse normal matrix
   if(verbose)write(stderr,*)"Reusing matrix"   
   
   call read_BINtype(reusemat,1) ! read reuse matrix
   
   select case(reusemat%type)
   case('SYMVN___','SYMV0___','SYMV1___','SYMV2___') ! ok

      if(reusemat%nval1 .ne. nunk)then
         write(stderr,*)"ERROR, Reusing matrix failed, wrong type or size"
         stop
      end if

      reusemat%type='P'  !keep in packed form ( quicker and needs less memory)
      call read_BINtype(reusemat) ! read remainder of the reuse matrix
      out%ldm=reusemat%ldm
      out%pack1=>reusemat%pack1 ! point to packed matrix
   case default
      write(stderr,*)"ERROR in BIN_buildnorm: Cannot reuse this type of matrix",reusemat%type   
      stop
   end select

else ! create normal matrix
if(verbose)write(stderr,*)"Calculating normal matrix"
   !create normal matrix
   allocate(out%mat1(nunk,nunk))
   out%ldm=nunk
   out%mtyp='U' ! type of the matrix will be upper triangular
   !symmetric rank operation
   
   
!DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
   if(trans)then
      call dsyrk('U','N',nunk,nobs,1.d0,A%mat1(1,1),int(A%ldm),0.d0,out%mat1(1,1),int(out%ldm))
   else
      call dsyrk('U','T',nunk,nobs,1.d0,A%mat1(1,1),int(A%ldm),0.d0,out%mat1(1,1),int(out%ldm))
   end if
   
   
end if


!set up output structure
out%nint=2
allocate(out%ints_d(out%nint),out%ints(out%nint))
out%ints_d(1)='Nobs'
out%ints_d(2)='Nunknowns'

out%ints(1)=nobs
out%ints(2)=nunk



out%ndbls=5
allocate(out%dbls_d(out%ndbls),out%dbls(out%ndbls))

out%dbls_d(1)='CTime'
out%dbls_d(2)='STime'
out%dbls_d(3)='ETime'
out%dbls_d(4)='LtPL'
out%dbls_d(5)='Sigma0'


out%dbls=0.d0

!take timetags from observation system
do,i=1,obs%ndbls
   select case(obs%dbls_d(i))
   case('CTime')
      out%dbls(1)=obs%dbls(i)
   case('STime')
      out%dbls(2)=obs%dbls(i)
   case('ETime')
      out%dbls(3)=obs%dbls(i)
   end select
end do

out%dbls(4)=ltpl
out%dbls(5)=1.d0

out%descr='Normal system from BIN_buildnorm'

out%nvec=2
out%type='SYMVN___'
out%nval1=nunk
out%nval2=out%nval1
out%pval1=nunk*(nunk+1)/2
out%pval2=1

!readme part

out%nread=6
allocate(out%readme(out%nread))
out%readme(1)="System is constructed from design matrix:"
out%readme(2)=A%file(1:80)
out%readme(3)=a%descr
out%readme(4)="Used observation set:"
out%readme(5)=obs%file(1:80)
out%readme(6)=obs%descr

!side descriptions
if(trans)then
   out%side1_d=>A%side1_d
else
   out%side1_d=>A%side2_d
end if
if(verbose)write(stderr,*)"Writing matrix"
!output normal equation
call write_BINtype(out) 

end program BIN_buildnorm

subroutine help()
implicit none
integer::stderr
character*8::frmt
stderr=0
frmt='(A)'
write(stderr,frmt)"Program BIN_buildnorm creates a normal system from an observation set"
write(stderr,frmt)"and a given designmatrix"
write(stderr,frmt)"Usage:BIN_buildnorm [OPTIONS] -A=DESIGNMAT OBSFILE"
write(stderr,frmt)"where DESIGNMAT is the file (binary format) containing the design matrix"
write(stderr,frmt)"and OBSFILE is the (binary)file containing the observations and its variances"
write(stderr,frmt)"OPTIONS may be:"
write(stderr,frmt)"-S=ST/END: compare the string only from index ST until END"
write(stderr,frmt)"-t:transpose design matrix"
write(stderr,frmt)"-R=REUSEMAT: just copy the normal matrix from file REUSEMAT"
write(stderr,frmt)"           NO permutation is performed!!"
write(stderr,frmt)"-V=VSELECT: Select the vector column to be used as data vector (default = 1)"
write(stderr,frmt)"-v: be verbose about the process ( prints to standard error)"
! write(stderr,frmt)""
stop
end subroutine help

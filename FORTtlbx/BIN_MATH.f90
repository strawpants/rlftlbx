!! program to perform low level matrix operations from the command line
!!Coded by Roelof Rietbroek, Mon Aug  4 10:35:56 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!todo:
!! support for triangular and diagonal matrices
!! a lot

!!Updated by Roelof Rietbroek, Tue Jan 27 09:50:49 2009
!!Added: eigenvalue decomposition for a symmetric matrix
!!Updated by Roelof Rietbroek, Tue Sep 15 14:29:03 2009
!!added FULL MATRIX --MATRIX multiplication
!!Updated by Roelof Rietbroek, Wed Sep 16 15:59:14 2009
!! added condition number calculation

!!Updated by Roelof Rietbroek, Thu Sep 30 12:42:19 2010
!! added the calculation of the correlation matrix

!!Updated by Roelof Rietbroek, Fri Apr 15 11:10:19 2011
!! added propagation operation A' B A
!! replaced the calculation of the permutation vector for a more efficient method

!!Updated by Roelof Rietbroek, Wed Jun 29 16:36:35 2011
!! ascii output is now compatible with SH_read

!!Updated by Roelof Rietbroek, Mon Jul 25 12:39:26 2011
!! added block diagonal matrix input for matrix vector multiplication
!!fixed bug with v2select (always picked the first column regardless of v2select)

!!Updated by Roelof Rietbroek, Tue Feb 21 17:01:58 2012
!! allowed patching of matrices

!!Updated by Roelof Rietbroek, Thu Aug 23 10:05:13 2012
!!improved checks on the common parameters of matching systems

!!Updated by Roelof Rietbroek, Thu Sep  6 10:13:31 2012
!! also added the option for a symmetric rank k update

!!Updated by Roelof Rietbroek, Thu Mar  7 20:16:27 2013
!added FULL2DVN+FULL2DVN

!!Updated by Roelof Rietbroek, Tue Jun  4 10:31:36 2013
!1 added rank K update for rectangular matrices

!!Updated by Roelof Rietbroek, Mon Aug 17 14:53:37 2015
!! added error propagation of a symemtric matrix through a full2d matrix


program BIN_MATH
use binfiletools
use bin_operations
use forttlbx
implicit none
integer::itharg,narg,istderr,maxf,unit,v2select,nf,stderr,info,lwork
integer::i,j,nrows1,ncols1,nrows2,ncols2,eof,st,nd,stind,ndind,nweights,ilaenv
parameter(stderr=0,maxf=2)!maximum amount of files in one run
character(200)::dum
character(1)::operation
type(BINdat)::work(maxf),out !space for at most maxf data structures (Plus an output system)
logical::trans1,trans2,both2,both1,file2ascii
character(50)::version
double precision::W(maxf)
character(24),pointer::rows1(:)=>null()
character(24),pointer::rows2(:)=>null()
character(24),pointer::cols1(:)=>null()
character(24),pointer::cols2(:)=>null()
character(24),pointer::tmpcols2(:)=>null()

integer,allocatable::sys2_perm1(:),sys2_perm2(:),sys1_perm1(:),sys1_perm2(:)
double precision,allocatable::workv(:),workm(:,:)
character(1)::uplo1,uplo2
integer::iargc
double precision::dumdb,rcond,anorm
integer,allocatable::ipiv(:),iwork(:)
integer::ncom,ncom1
integer::sz,ipack,ivec,shftcol,shftrow
logical::permutesys1,permutesys2
!defaults/initializations
! perm1=>null()
! perm2=>null()
stind=1
ndind=24

!default prints file to standard output
version='1.0 EXPERIMENTAL(not complete)'
unit=13
out%file='stdout'
! incstr=''
! excstr=''
file2ascii=.false.
permutesys1=.false.
permutesys2=.true.
narg=iargc()
if( narg .eq. 0)call help(version)

itharg=0
nf=0
trans1=.false.
trans2=.false.
v2select=0! default vector column is 0 ( matrix is used)
W=1.d0 ! weights
nweights=0
!process command line arguments
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
!   write(0,*)trim(dum)
   if(dum(1:1) .eq. '-')then !option
      select case(dum(2:2))
      case('O')! perform operation
         operation=dum(4:)
!       case('E')!exclude string
!          excstr=trim(excstr)//'|'//trim(dum(4:)) !append string to include string
!       case('I')!include string
!          incstr=trim(incstr)//'|'//trim(dum(4:)) !append string to include string
      case('f') !specify output file name
         if(dum(4:4) .eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,out%file)
         else
            out%file=dum(4:)
         end if
      case('S')!limit global search area of the parameter strings
         st=4
         nd=index(dum(st:),'/')+st-2
         read(dum(st:nd),*)stind
         read(dum(nd+2:),*)ndind
      case('W') !provide weight for each system
         st=4
         do,j=1,maxf
            nd=index(dum(st:),'/')+st-2
            if(nd .eq. st-2)nd=len(dum)!set to end of string if no '/' occurs (only one weight)
            read(dum(st:nd),*)W(j) !read in weight
            if(nd .eq. (len(dum)))exit !last weight has been read
            st=nd+2
         end do
         nweights=j !set number of weights provided
      case('a') ! specify whether file2 is ascii ( contains a vector or diag matrix)
         file2ascii=.true.
         v2select=1 ! set vector selector to 1
      case('t')! transpose matrices
         select case(dum(3:3))
         case('1')
            trans1=.true.
         case('2')
            trans2=.true.
         end select
      case('V')! select vector
         read(dum(4:),*)v2select ! read vector column to be used from argument
      case default ! call help screen
         write(stderr,*)"ERROR: not an option:",trim(dum)
         call help(version)
      end select
   else ! filename
      nf=nf+1
      if(nf > maxf)then
         write(stderr,*)"ERROR: maximum amount of input files is:", maxf
         stop
      end if

      work(nf)%file=trim(dum)
   end if
   
end do

! !possibly strip leading |
! if(incstr(1:1) .eq. '|')incstr=incstr(2:)
! if(excstr(1:1) .eq. '|')excstr=excstr(2:)

!some input checks
if( nf .eq. 0)then
   write(stderr,*)"ERROR: no input files"
   stop
end if

if(nweights .ne. nf .and. nweights .ne. 0)then
   write(stderr,*)"ERROR: amount of weights provided is inconsistent with the amount of files"
   stop
end if

select case(operation)
case('+')
   if(nf .ne. 2)then
      write(stderr,*)"ERROR: operation '",operation,"' requires two input systems"
   stop
   end if

   if(v2select > 0)then
      write(stderr,*)"ERROR: operation '+' is not allowed on vectors"
      stop
   end if

   if(trans1)then
      write(stderr,*)"ERROR: operation '+' does not allow the transpose of the first system"
      stop
   end if

! case('c')
!    if(nf >1)then
!       write(stderr,*)"ERROR: operation '",operation,"' requires only one input system"
!    stop
!    end if
case('x','s','p','P')
   if(nf .ne. 2)then
      write(stderr,*)"ERROR: operation '",operation,"' requires two input systems"
   stop
   end if
case('e','c','C','k') !operations with only one input system
   if(nf .ne. 1)then
      write(stderr,*)"ERROR: operation '",operation,"' requires only one input system"
   stop
   end if

case default ! stop program
   write(stderr,*)"ERROR: Unknown (or unsupported) operation supplied: '",operation,"'"
   stop
end select


if(v2select > 0 .and. trans2)then
   write(stderr,*)"ERROR: matrix transpose of system 2 is illegal, since its vector is used"
   stop
end if



!!!!!!!!!!!!CONFIGURATION!!!!!!!!!!!!

!SORTING PART

!read in meta data of first input
call read_BINtype(work(1),2)


!determine which matrix side will be matched
select case(work(1)%type)
case('SYMV0___','SYMV1___','SYMV2___','SYMVN___')
   nrows1=work(1)%nval1
   ncols1=work(1)%nval1
   rows1=>work(1)%side1_d
   cols1=>work(1)%side1_d
   work(1)%mtyp='U'
   trans1=.false. ! a transpose makes no sense here
   both1=.true.
case('DIAVN___')
   nrows1=work(1)%nval1
   ncols1=work(1)%nval1
   rows1=>work(1)%side1_d
   cols1=>work(1)%side1_d
   work(1)%mtyp='P'
   trans1=.false. ! a transpose makes no sense here
   both1=.true.
case('FULLSQV0','FULLSQVN')
   nrows1=work(1)%nval1
   ncols1=work(1)%nval1
   rows1=>work(1)%side1_d
   cols1=>work(1)%side2_d
   work(1)%mtyp='F'
   both1=.false.
case('FULL2DVN')
   work(1)%mtyp='F'
   if(trans1)then
      nrows1=work(1)%nval2
      ncols1=work(1)%nval1
      rows1=>work(1)%side2_d
      cols1=>work(1)%side1_d
   else   
      nrows1=work(1)%nval1
      ncols1=work(1)%nval2
      rows1=>work(1)%side1_d
      cols1=>work(1)%side2_d
   end if
   both1=.false.
case('BDFULLVN')
   work(1)%mtyp='P'
   if(trans1)then
      nrows1=work(1)%nval2
      ncols1=work(1)%nval1
      rows1=>work(1)%side2_d
      cols1=>work(1)%side1_d
   else   
      nrows1=work(1)%nval1
      ncols1=work(1)%nval2
      rows1=>work(1)%side1_d
      cols1=>work(1)%side2_d
   end if
   both1=.false.
case default 
   write(stderr,*)"ERROR: matrix type is not supported:",work(1)%type
   stop
end select

!read in meta data of second input
if(nf .eq. 2)then
   if(file2ascii)then !read data from ascii file
      if(work(2)%file .eq. 'stdin')then
         write(stderr,*)"ERROR: ascii files may not be read from standard input"
         stop
      end if
      open(unit=unit,file=work(2)%file,form='formatted')
      eof=0
      work(2)%nval1=0
      do while (eof .eq. 0)
         read(unit,iostat=eof,fmt=*)dum(1:1)
         if(eof .ne. 0)exit
         work(2)%nval1=work(2)%nval1+1
      end do
      rewind(unit)
      work(2)%nvec=1
      work(2)%type='DIAVN___' !Create diag systemdummy type to ensure that the system is accepted
      work(2)%nval2=work(2)%nval1
      work(2)%pval1=work(2)%nval1
      work(2)%pval2=1

      allocate(work(2)%vec(work(2)%nval1,1),work(2)%side1_d(work(2)%nval1))
      do,i=1,work(2)%nval1
         read(unit,'(a24)',advance='NO')work(2)%side1_d(i)
         read(unit,*)work(2)%vec(i,1)
      end do
      
      close(unit)
      !copy vector to diagonal
      allocate(work(2)%pack1(work(2)%nval1))
      work(2)%pack1=work(2)%vec(:,1) ! copy vector data to diagonal (data can be used both as  diagonal matix and vector)
   else ! read from binary file
      call read_BINtype(work(2),2)
   end if

   select case(work(2)%type)
   case('SYMV0___','SYMV1___','SYMV2___','SYMVN___')
      nrows2=work(2)%nval1
      ncols2=work(2)%nval1
      rows2=>work(2)%side1_d
      cols2=>work(2)%side1_d
      work(2)%mtyp='U'
      trans2=.false. !makes no sense here
      both2=.true. ! permute both sides to keep symmetry
   case('DIAVN___') ! diagonal matrix
      nrows2=work(2)%nval1
      ncols2=work(2)%nval1
      rows2=>work(2)%side1_d
      cols2=>work(2)%side1_d
      work(2)%mtyp='P' !
      trans2=.false. !makes no sense here
      both2=.true.


   case('FULLSQV0','FULLSQVN')
      nrows2=work(2)%nval1
      ncols2=work(2)%nval1
      rows2=>work(2)%side1_d
      cols2=>work(2)%side2_d
      work(2)%mtyp='F'
      both2=.false.
   case('FULL2DVN')
      work(2)%mtyp='F'
      if(trans2)then
         nrows2=work(2)%nval2
         ncols2=work(2)%nval1
         rows2=>work(2)%side2_d
         cols2=>work(2)%side1_d
      else   
         nrows2=work(2)%nval1
         ncols2=work(2)%nval2
         rows2=>work(2)%side1_d
         cols2=>work(2)%side2_d
      end if
      both2=.false.
   case default 
      write(stderr,*)"ERROR: matrix type is not supported:",work(2)%type
      stop
   end select

   !create permutation vector(s) for the second system

   select case(operation)
   case('+') !both sides may need to be permuted
      if(nrows1 .ne. nrows2 .or. ncols1 .ne. ncols2)then
         write(stderr,*)"ERROR: matrix systems have non matching dimensions"
         stop
      end if
      allocate(sys2_perm1(nrows2))
      call get_permvec(rows2,rows1,.true.,stind,ndind,sys2_perm1,ncom)
      if(ncom .ne. nrows2)then
         write(stderr,*)"ERROR: system 1 has parameters not present in system 2"
         stop
      end if
      ! do,i=1,nrows1
      !    do,j=1,nrows1
      !       if(rows1(i)(stind:ndind) .eq. rows2(j)(stind:ndind))perm1(i)=j
      !    end do
      ! end do
      

      select case(work(2)%type)
      case('FULL2DVN') ! also check second side (make second permutation vector)
         allocate(sys2_perm2(ncols2))
         call get_permvec(cols2,cols1,.true.,stind,ndind,sys2_perm2,ncom)
         ! do,i=1,ncols1
         !    do,j=1,ncols1
         !       if(cols1(i)(stind:ndind) .eq. cols2(j)(stind:ndind))perm2(i)=j
         !    end do
         ! end do
      case default !do nothing
      end select
   case('s')
      if( nrows1 .ne. ncols1)then
         write(stderr,*)"ERROR: matrix of system 1 is rectangular and not invertible"
         stop
      else if(nrows1 .ne. nrows2)then
         write(stderr,*)"ERROR: matrix systems do not match dimensions"
         stop
      end if
      
      allocate(sys2_perm1(nrows2))
      call get_permvec(rows2,rows1,.true.,stind,ndind,sys2_perm1,ncom)

      if(ncom .ne. nrows1)then
         write(stderr,*)"ERROR: Systems share not enough common parameters, ncom:",ncom
         stop
      end if

      ! do,i=1,nrows1
      !    do,j=1,nrows1
      !       if(rows1(i)(stind:ndind) .eq. rows2(j)(stind:ndind))perm1(i)=j
      !    end do
      ! end do
   case('x')
      if(ncols1 .ne. nrows2)then
         write(stderr,*)"ERROR: matrix systems have non matching dimensions"
         stop
      end if
      allocate(sys2_perm1(nrows2))
      call get_permvec(rows2,cols1,.true.,stind,ndind,sys2_perm1,ncom)

      if(ncom .ne. ncols1)then
         write(stderr,*)"ERROR: Systems share not enough common parameters, ncom:",ncom
         stop
      end if

      
!       do,i=1,ncols1
!          do,j=1,ncols1
! !            write(0,*)cols1(i)(stind:ndind),rows2(j)(stind:ndind)
!             if(cols1(i)(stind:ndind) .eq. rows2(j)(stind:ndind))perm1(i)=j
!          end do
!       end do
   case('p')
      select case(work(1)%type)
      case('SYMV0___','SYMV1___','SYMV2___','SYMVN___') ! also check second side (make second permutation vector)
      case default
         write(stderr,*)"ERROR: Matrix 1 must be symmetric"
         stop
      end select

      if(ncols1 .ne. nrows2)then
         write(stderr,*)"ERROR: matrix systems have non matching dimensions"
         stop
      end if
      allocate(sys2_perm1(nrows2))
      call get_permvec(rows2,cols1,.true.,stind,ndind,sys2_perm1,ncom)
   case('P') !both systems need permuting
      select case(work(1)%type)
      case('FULL2DVN','FULLSQVN','SYMV0___','SYMV1___','SYMV2___','SYMVN___') !allowed types
         permutesys1=.true.
      case default
         write(stderr,*)"ERROR: Matrix 1 type is not allowed:",work(1)%type
         stop
      end select

      
      allocate(sys2_perm2(ncols2))
      call get_permvec(cols2,cols1,.true.,stind,ndind,sys2_perm2,ncom)
      
      !create a temporary vector which holds the columns of system 2 after permutation
      allocate(tmpcols2(ncols2))
      tmpcols2=cols2
      call permute(tmpcols2,sys2_perm2)
      !create permutation vector for system 1 (put non-matching entries in front)
      allocate(sys1_perm2(ncols1))
      call get_permvec(cols1,tmpcols2,.false.,stind,ndind,sys1_perm2,ncom1)
      
   end select
   
   
end if


!!DATA READING

!read in remaining part (system 1/2)
call read_BINtype(work(1))
if(nf .eq. 2 .and. .not. file2ascii)call read_BINtype(work(2))
!!END DATA READING

!!PERMUTATION OF SYSTEM 1
if(permutesys1)then
   if(trans1)then
      if(allocated(sys1_perm2) .and. allocated(sys1_perm1))then
         call BIN_permute(input=work(1),perm1=sys1_perm2,perm2=sys1_perm1)
      else if(allocated(sys1_perm1))then
         call BIN_permute(input=work(1),perm2=sys1_perm1,both=both1)
      else if(allocated(sys1_perm2))then
         call BIN_permute(input=work(1),perm1=sys1_perm2,both=both1)
      end if
   else
      if(allocated(sys1_perm2) .and. allocated(sys1_perm1))then
         call BIN_permute(input=work(1),perm1=sys1_perm1,perm2=sys1_perm2)
      else if(allocated(sys1_perm1))then
         call BIN_permute(input=work(1),perm1=sys1_perm1,both=both1)
      else if(allocated(sys1_perm2))then
         call BIN_permute(input=work(1),perm2=sys1_perm2,both=both1)
      end if
   end if
      
end if

!!PERMUTATION OF SYSTEM 2
if(nf .eq. 2 .and. permutesys2)then
   if(trans2)then
      if(allocated(sys2_perm2) .and. allocated(sys2_perm1))then
         call BIN_permute(input=work(2),perm1=sys2_perm2,perm2=sys2_perm1)
      else if(allocated(sys2_perm1))then
         call BIN_permute(input=work(2),perm2=sys2_perm1,both=both2)
      else if(allocated(sys2_perm2))then
         call BIN_permute(input=work(2),perm1=sys2_perm2,both=both2)
      end if
   else
      if(allocated(sys2_perm2) .and. allocated(sys2_perm1))then
         call BIN_permute(input=work(2),perm1=sys2_perm1,perm2=sys2_perm2)
      else if (allocated(sys2_perm1))then
         call BIN_permute(input=work(2),perm1=sys2_perm1,both=both2)
      else if (allocated(sys2_perm2))then
         call BIN_permute(input=work(2),perm2=sys2_perm2,both=both2)
      end if
   end if
      
end if


!system two is a vector system
if(v2select > 0)work(2)%type='VVVVVVVV' !set type of the second system to dummy VVVVVVVV

! or not present (single system)
if(nf .eq. 1)work(2)%type='SSSSSSSS' !set type of the second system to dummy SSSSSSSS


!!!PERFORM OPERATION!!!!

!determine combination of types
select case(work(1)%type) ! type of 1st system
case('SYMV0___','SYMV1___','SYMV2___','SYMVN___')!symmetric matrix
   select case(work(2)%type)
   case('SYMV0___','SYMV1___','SYMV2___','SYMVN___')!symmetric matrix
      !now select operation

      select case(operation)
      case('+')! addition of 2 symmetric matrices
         !addition of 2 symmetric systems ( no transposing necessary)
         do,i=1,work(1)%nval1
            do,j=1,i ! upper triangle only
               work(1)%mat1(j,i)=W(1)*work(1)%mat1(j,i)+W(2)*work(2)%mat1(j,i)
            end do
         end do
         !CONSTRUCT OUTPUT SYSTEM
         out%type='SYMVN___'
         out%descr='BIN_MATH '//work(1)%type//' '//operation//' '//work(2)%type
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1*(out%nval1+1)/2 !symmetric matrix in packed form
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'
         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%mat1=>work(1)%mat1 ! point to right memory section
      case('s')!solve for full square matrix
         
         !copy upper triangle of second system in lower triangle
         do,i=1,work(2)%nval1
            do,j=1,i
               work(2)%mat1(i,j)=work(2)%mat1(j,i)
            end do
         end do
!         write(0,*)size(work(1)%mat1,1),size(work(1)%mat1,2),work(1)%nval1
         ! call LAPACK/BLAS routines
         !Cholesky decomposition
         call dpotrf('U',work(1)%nval1,work(1)%mat1(1,1),work(1)%nval1,info)
         
         if(info .ne. 0)then
            write(stderr,*)"ERROR: Cholesky decomposition failed"
            stop
         end if

         !triangular matrix solve(number 1) & scaling
         call dtrsm('L','U','T','N',work(1)%nval1,work(1)%nval1,W(2)/W(1)&
              ,work(1)%mat1(1,1),work(1)%nval1,work(2)%mat1(1,1),work(1)%nval1)

         !triangular matrix solve(number 2)
         call dtrsm('L','U','N','N',work(1)%nval1,work(1)%nval1,1.d0&
              ,work(1)%mat1(1,1),work(1)%nval1,work(2)%mat1(1,1),work(1)%nval1)


         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULLSQVN' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' **(-1) '//work(2)%type
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1**2 
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='F'
         out%side1_d=>work(2)%side1_d ! point to the right side description
         out%side2_d=>work(2)%side1_d ! Both sides share the same parameters
         out%mat1=>work(2)%mat1 ! point to right memory section

      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('FULL2DVN')!full rectangular matrix or square matrix
      select case(operation)
      case('x') !Symmetric MATRIX -- MATRIX multiplication
         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULL2DVN' ! Matrix system as output
         out%descr='BIN_MATH: '//work(1)%type//operation//work(2)%type

         out%nval1=nrows1 ! Size of the output vector
         out%nval2=ncols2 ! columns of matrix 2
         out%pval1=out%nval1*out%nval2
         out%pval2=1
         out%mtyp='F'

         allocate(out%mat1(out%nval1,out%nval2))
         out%nvec=0 
         out%nread=3
         allocate(out%readme(out%nread))
         out%readme(1)="used input files:"
         out%readme(2)=work(1)%file(1:80)
         out%readme(3)=work(2)%file(1:80)
         out%nint=0
         out%ndbls=0
         out%side1_d=>rows1 ! point to the right side description
         out%side2_d=>cols2

         if(trans2)then
            write(stderr,*)"ERROR, Transpose of second matrix is not allowed: ",work(1)%type,operation,work(2)%type
            write(stderr,*)"ERROR, You may switch the order of the input to obtain the transpose result."
         stop
         end if
         !do calculation (symmetric multiplication and scaling )
         
         call dsymm('L',work(1)%mtyp,nrows1,ncols2,W(1)*W(2),&
              work(1)%mat1(1,1),work(1)%nval1,&
              work(2)%mat1(1,1),work(2)%nval1,&
              0.d0,out%mat1(1,1),out%nval1)

      case('s')
         ! call LAPACK/BLAS routines
         !Cholesky decomposition
         call dpotrf('U',work(1)%nval1,work(1)%mat1(1,1),work(1)%nval1,info)
         
         if(info .ne. 0)then
            write(stderr,*)"ERROR: Cholesky decomposition failed"
            stop
         end if
         
         if(trans2)then
            uplo2='R' !note the transposed matrix is outputted
            write(stderr,*)"MESSAGE: output is transposed!"
            out%descr='BIN_MATH '//work(2)%type//' '//work(1)%type//' **(-1)'

         else
            uplo2='L'
            out%descr='BIN_MATH '//work(1)%type//' **(-1) '//work(2)%type
         end if

         !triangular matrix solve(number 1) & scaling
         call dtrsm(uplo2,'U','T','N',work(1)%nval1,work(2)%nval2,W(2)/W(1)&
              ,work(1)%mat1(1,1),work(1)%nval1,work(2)%mat1(1,1),work(2)%nval1)

         !triangular matrix solve(number 2)
         call dtrsm(uplo2,'U','N','N',work(1)%nval1,work(2)%nval2,1.d0&
              ,work(1)%mat1(1,1),work(1)%nval1,work(2)%mat1(1,1),work(2)%nval1)


         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULL2DVN' ! fully populated rectangular matrix

         out%nval1=work(2)%nval1
         out%nval2=work(2)%nval2
         out%pval1=out%nval1*out%nval2 
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='F'
         out%side1_d=>work(2)%side1_d ! point to the right side description
         out%side2_d=>work(2)%side2_d 
         out%mat1=>work(2)%mat1 ! point to right memory section

      case('P') !patch matrix
         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULL2DVN' ! Matrix system as output
         out%descr='BIN_MATH: '//work(1)%type//operation//work(2)%type

         out%nval1=nrows1+nrows2 ! Size of the output vector
         out%nval2=ncols2+ncols1-ncom ! columns of matrix sum minus the amount of common parameters
         out%pval1=out%nval1*out%nval2
         out%pval2=1
         out%mtyp='F'

         allocate(out%mat1(out%nval1,out%nval2))
         out%mat1=0.d0
         allocate(out%side1_d(out%nval1))
         allocate(out%side2_d(out%nval2))

         out%nvec=0 
         out%nread=3
         allocate(out%readme(out%nread))
         out%readme(1)="used input files:"
         out%readme(2)=work(1)%file(1:80)
         out%readme(3)=work(2)%file(1:80)
         out%nint=0
         out%ndbls=0

         out%side1_d(1:nrows1)=rows1

         out%side1_d(nrows1+1:out%nval1)=rows2 !copy other system below

         !copy first point to the right side description
         out%side2_d(1:ncols1)=cols1
         out%side2_d(ncols1+1:out%nval2)=cols2(ncom+1:ncols2) ! add paramters from system 2 but skip common parameters
         
         !copy values from system 1
         if(trans1)then
            do,j=1,ncols1
               do,i=1,nrows1
                  out%mat1(i,j)=work(1)%mat1(j,i)
                  !additionally mirror the matrix
                  out%mat1(j,i)=work(1)%mat1(j,i)
               end do
            end do            
         else
            do,j=1,ncols1
               do,i=1,nrows1
                  out%mat1(i,j)=work(1)%mat1(i,j)
                  !additionally mirror the matrix
                  out%mat1(j,i)=work(1)%mat1(i,j)
               end do
            end do            
         end if

         !copy values from system 2
         shftcol=ncols1-ncom
         shftrow=nrows1         

         if(trans2)then
            do,j=1,ncols2
               do,i=1,nrows2
                  out%mat1(i+shftrow,j+shftcol)=work(2)%mat1(j,i)
               end do
            end do            
         else
            do,j=1,ncols2
               do,i=1,nrows2
                  out%mat1(i+shftrow,j+shftcol)=work(2)%mat1(i,j)
               end do
            end do            
         end if
      case('p') !propagate symmetric matrix through other matrix

         !OUT=A2' A1 A2 or OUT =A2 A1 A2'

         ! split up A1 in a upper and lower triangular matrix where the diagonal is divided by two
         !A1= L + L'

         do,j=1,nrows1
            work(1)%mat1(j,j)=work(1)%mat1(j,j)/2.d0
         end do

         !copy A2 in a temporary matrix
         allocate(workm(work(2)%nval1,work(2)%nval2))
         workm=work(2)%mat1 !deep copy of data

         if(trans2)then
!            write(0,*)work(2)%nval1,work(2)%nval2,work(1)%nval1
            !compute B = A2 *L
            call dtrmm('R','L','N','N',work(2)%nval1,work(2)%nval2,1.d0,work(1)%mat1(1,1),work(1)%nval1,workm,work(2)%nval1)
            
         
         else

         !compute B = L * A2
            call dtrmm('L','L','N','N',work(2)%nval1,work(2)%nval2,1.d0,work(1)%mat1(1,1),work(1)%nval1,workm,work(2)%nval1)

         end if

         !CONSTRUCT OUTPUT SYSTEM


         out%type='SYMVN___' ! fully populated square matrix
         out%descr='BIN_MATH '//work(2)%type//' * '//work(1)%type//' * '//work(2)%type
         out%nval1=ncols2
         allocate (out%mat1(out%nval1,out%nval1))
         out%nval2=out%nval1
         out%pval1=out%nval1*(out%nval1+1)/2 
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'
        out%side1_d=>cols2 ! point to the right side description 

         if(trans2)then
            !compute OUT= A2 B' + B A2'
            call dsyr2k ('U','N',ncols2,nrows2,1.d0,work(2)%mat1(1,1),work(2)%nval1,workm(1,1),&
                 work(2)%nval1,0.d0,out%mat1(1,1),out%nval1)
         else
            !compute OUT= A2' B + B' A2
            call dsyr2k ('U','T',ncols2,nrows2,1.d0,work(2)%mat1(1,1),work(2)%nval1,workm(1,1),&
                 work(2)%nval1,0.d0,out%mat1(1,1),out%nval1)
         end if
         

         

         
      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('FULLSQV0','FULLSQVN')! full square matrix
      select case(operation)
      case('x') !Symmetric MATRIX -- MATRIX multiplication
         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULLSQVN' ! Matrix system as output
         out%descr='BIN_MATH: '//work(1)%type//operation//work(2)%type

         out%nval1=nrows1 ! Size of the output vector
         out%nval2=ncols2 ! columns of matrix 2
         out%pval1=out%nval1*out%nval2
         out%pval2=1
         out%mtyp='F'

         allocate(out%mat1(out%nval1,out%nval2))
         out%nvec=0 
         out%nread=3
         allocate(out%readme(out%nread))
         out%readme(1)="used input files:"
         out%readme(2)=work(1)%file(1:80)
         out%readme(3)=work(2)%file(1:80)
         out%nint=0
         out%ndbls=0
         out%side1_d=>rows1 ! point to the right side description
         out%side2_d=>cols2

         if(trans2)then
            write(stderr,*)"ERROR, Transpose of second matrix is not allowed: ",work(1)%type,operation,work(2)%type
            write(stderr,*)"ERROR, You may switch the order of the input to obtain the transpose result."
         stop
         end if
         !do calculation (symmetric multiplication and scaling )
         
         call dsymm('L',work(1)%mtyp,nrows1,ncols2,W(1)*W(2),&
              work(1)%mat1(1,1),work(1)%nval1,&
              work(2)%mat1(1,1),work(2)%nval1,&
              0.d0,out%mat1(1,1),out%nval1)

      case('s')
                  ! call LAPACK/BLAS routines
         !Cholesky decomposition
         call dpotrf('U',work(1)%nval1,work(1)%mat1(1,1),work(1)%nval1,info)
         
         if(info .ne. 0)then
            write(stderr,*)"ERROR: Cholesky decomposition failed"
            stop
         end if
         
         if(trans2)then
            uplo2='R' !note the transposed matrix is outputted
            write(stderr,*)"MESSAGE: output is transposed!"
         else
            uplo2='L'
         end if

         !triangular matrix solve(number 1) & scaling
         call dtrsm(uplo2,'U','T','N',work(1)%nval1,work(2)%nval2,W(2)/W(1)&
              ,work(1)%mat1(1,1),work(1)%nval1,work(2)%mat1(1,1),work(2)%nval1)

         !triangular matrix solve(number 2)
         call dtrsm(uplo2,'U','N','N',work(1)%nval1,work(2)%nval2,1.d0&
              ,work(1)%mat1(1,1),work(1)%nval1,work(2)%mat1(1,1),work(2)%nval1)


         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULLSQVN' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' **(-1) '//work(2)%type
         out%nval1=work(2)%nval1
         out%nval2=work(2)%nval2
         out%pval1=out%nval1*out%nval2
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='F'
         out%side1_d=>work(2)%side1_d ! point to the right side description
         out%side2_d=>work(2)%side2_d 
         out%mat1=>work(2)%mat1 ! point to right memory section
      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select

   case('DIAVN___')! diagonal matrix
      select case(operation)
      case('+') ! add diagonal matrix to symmetric matrix
         do,i=1,work(1)%nval1
            do,j=1,i
               work(1)%mat1(j,i)=W(1)*work(1)%mat1(j,i)
            end do
            !add diagonal
            work(1)%mat1(i,i)=work(1)%mat1(i,i)+W(2)*work(2)%pack1(i)
         end do
         
         !CONSTRUCT OUTPUT SYSTEM
         out%type='SYMVN___' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' + '//work(2)%type
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1*(work(1)%nval1+1)/2 
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'
         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%mat1=>work(1)%mat1 ! point to right memory section

      case('p') !propagate symmetric matrix through diagonal matrix
         do,i=1,work(1)%nval1
            do,j=1,i
               work(1)%mat1(j,i)=W(1)*W(2)*work(1)%mat1(j,i)*work(2)%pack1(i)*work(2)%pack1(j)
            end do
         end do

         !CONSTRUCT OUTPUT SYSTEM
         out%type='SYMVN___' ! fully populated square matrix
         out%descr='BIN_MATH '//work(2)%type//' * '//work(1)%type//' * '//work(2)%type
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1*(work(1)%nval1+1)/2 
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'
         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%mat1=>work(1)%mat1 ! point to right memory section
      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('VVVVVVVV')! use vector of second system
      select case(operation)
      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('SSSSSSSS') ! single system operation only
      select case(operation)
      case('e')!eigenvalue decomposition
         
         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULLSQVN' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' eigenvalue decompition'
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1*out%nval1 
         out%pval2=1
         out%nvec=1 ! hold the eigenvalues
         allocate(out%vec(out%nval1,out%nvec))
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='F'

         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%side2_d=>work(1)%side1_d ! second side is within the same parameter space
         out%mat1=>work(1)%mat1 ! point to right memory section


         !get optimal blocksize
         lwork=ILAENV( 1, 'DSYTRD','U', out%nval1, -1, -1, -1 )
         lwork = max(1,(lwork+2)*out%nval1)
         allocate(workv(lwork))

         !call lapack 
         call dsyev('V','U',out%nval1,out%mat1(1,1),out%nval1,out%vec(1,1),workv(1),lwork,info)
         if(info <0)then
            write(stderr,*)"ERROR, illegal value occured at: ",out%side1_d(info)
            stop
         else if (info >0)then
            write(stderr,*)"ERROR, failed to converge"
            stop
         end if
      case('c')! calculate condition number
         !allocate working array
         allocate(workv(work(1)%nval1*3))
         allocate(iwork(work(1)%nval1))
         !calculate matrix norm (use working array)
         anorm=0.d0
         workv=0.d0
         do,i=1,work(1)%nval1
            do,j=1,i-1 !off diagonal elements
               dumdb=abs(work(1)%mat1(j,i))
               workv(j)=workv(j)+dumdb !contribution to column j
               workv(i)=workv(i)+dumdb !contribution to column i
            end do
            workv(i)=workv(i)+abs(work(1)%mat1(i,i))
         end do
!         workv now holds column L1 norms
         anorm=maxval(workv(1:work(1)%nval1))
!         write(*,*)"ANORM",anorm

         !Cholesky decomposition
         write(stderr,*)"Calculating Cholesky factorization"
         call dpotrf('U',work(1)%nval1,work(1)%mat1(1,1),work(1)%nval1,info)
         if(info .ne. 0)then
            write(stderr,*)"ERROR:cholesky factorization failed,info:",info
            stop
         end if

         


         write(stderr,*)"Calculation of the condition number"
         call dpocon('U',work(1)%nval1,work(1)%mat1(1,1),work(1)%nval1&
              ,anorm,rcond,workv(1),iwork,info)
         if(info .ne. 0)then
            write(stderr,*)"ERROR: calculation of condition number failed,info:",info
            stop
         end if
         write(*,*)"Condition number, 1-NORM"
         write(*,*)1/rcond,anorm
         stop !no further output
      case('C') ! calculate correlation matrix
         !CONSTRUCT OUTPUT SYSTEM
         out%type='SYMVN___' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' correlation matrix'
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1*(out%nval1+1)/2
         out%pval2=1
         out%nvec=0 
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'

         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%side2_d=>work(1)%side1_d ! second side is within the same parameter space
         out%mat1=>work(1)%mat1 ! point to right memory section

         !take the quare root of the diagonal
         do,i=1,out%nval1
            out%mat1(i,i)=sqrt(out%mat1(i,i))
         end do

         do,i=1,out%nval1
            do,j=1,i-1 !off diagonal elements first ( keep standard deviations
               out%mat1(j,i)=out%mat1(j,i)/(out%mat1(i,i)*out%mat1(j,j))
            end do
         end do
         !set diagonal to 1
         do,i=1,out%nval1
            out%mat1(i,i)=1
         end do
      case('k') ! rank k operation
         !CONSTRUCT OUTPUT SYSTEM
         out%type='SYMVN___' ! fully symmetric matrix
         out%descr='BIN_MATH '//work(1)%type//' Rank K operation matrix'
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1*(out%nval1+1)/2
         out%pval2=1
         out%nvec=0 
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'

         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%side2_d=>work(1)%side1_d ! second side is within the same parameter space
         allocate(out%mat1(out%nval1,out%nval2))! allocate new memory

         ! mirror upper triangular part to lower part of the input pmatrix
         do,j=1,work(1)%nval1
            do,i=1,j-1 
               work(1)%mat1(j,i)=work(1)%mat1(i,j)
            end do
         end do
         call dsyrk(out%mtyp,'N',out%nval1,out%nval1,W(1),&
              work(1)%mat1(1,1),work(1)%nval1,0.d0,out%mat1(1,1),out%nval1)
         

      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation
         stop
      end select
   case default
      write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
      stop
   end select
case('FULL2DVN','FULLSQV0','FULLSQVN')!full rectangular or square matrix
   select case(work(2)%type)
   case('FULL2DVN','FULLSQV0','FULLSQVN')
      select case(operation)
      case('+') ! add matrices together
         do,i=1,work(1)%nval2
            do,j=1,work(1)%nval1
               work(1)%mat1(j,i)=W(1)*work(1)%mat1(j,i)+W(2)*work(2)%mat1(j,i)
            end do
         end do
         

         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULL2DVN' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' + '//work(2)%type
         out%nval1=work(1)%nval1
         out%nval2=work(1)%nval2
         out%pval1=work(1)%nval1*work(1)%nval2
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='F'
         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%side2_d=>work(1)%side2_d
         out%mat1=>work(1)%mat1 ! point to right memory section (which is now updated)
      case('x') !MATRIX -- MATRIX multiplication
         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULL2DVN' ! Matrix system as output
         out%descr='BIN_MATH: '//work(1)%type//operation//work(2)%type

         out%nval1=nrows1 ! Size of the output vector
         out%nval2=ncols2 ! columns of matrix 2
         out%pval1=out%nval1*out%nval2
         out%pval2=1
         out%mtyp='F'

         allocate(out%mat1(out%nval1,out%nval2))
         out%nvec=0 
         out%nread=3
         allocate(out%readme(out%nread))
         out%readme(1)="used input files:"
         out%readme(2)=work(1)%file(1:80)
         out%readme(3)=work(2)%file(1:80)
         out%nint=0
         out%ndbls=0
         out%side1_d=>rows1 ! point to the right side description
         out%side2_d=>cols2

         if(trans1)then
            uplo1='T'
         else
            uplo1='N'
         end if
         
         if(trans2)then
            uplo2='T'
         else
            uplo2='N'
         end if
         !do calculation ( multiplication and scaling )
         
         call dgemm(uplo1,uplo2,nrows1,ncols2,ncols1, W(1)*W(2),&
              work(1)%mat1(1,1),work(1)%nval1,&
              work(2)%mat1(1,1),work(2)%nval1,&
              0.d0,out%mat1(1,1),out%nval1)
      case('P') !patch matrix
         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULL2DVN' ! Matrix system as output
         out%descr='BIN_MATH: '//work(1)%type//operation//work(2)%type

         out%nval1=nrows1+nrows2 ! Size of the output vector
         out%nval2=ncols2+ncols1-ncom ! columns of matrix sum minus the amount of common parameters
         out%pval1=out%nval1*out%nval2
         out%pval2=1
         out%mtyp='F'

         allocate(out%mat1(out%nval1,out%nval2))
         out%mat1=0.d0
         allocate(out%side1_d(out%nval1))
         allocate(out%side2_d(out%nval2))

         out%nvec=0 
         out%nread=3
         allocate(out%readme(out%nread))
         out%readme(1)="used input files:"
         out%readme(2)=work(1)%file(1:80)
         out%readme(3)=work(2)%file(1:80)
         out%nint=0
         out%ndbls=0

         out%side1_d(1:nrows1)=rows1

         out%side1_d(nrows1+1:out%nval1)=rows2 !copy other system below

         !copy first point to the right side description
         out%side2_d(1:ncols1)=cols1
         out%side2_d(ncols1+1:out%nval2)=cols2(ncom+1:ncols2) ! add paramters from system 2 but skip common parameters
         
         !copy values from system 1
         if(trans1)then
            do,j=1,ncols1
               do,i=1,nrows1
                  out%mat1(i,j)=work(1)%mat1(j,i)
               end do
            end do            
         else
            do,j=1,ncols1
               do,i=1,nrows1
                  out%mat1(i,j)=work(1)%mat1(i,j)
               end do
            end do            
         end if

         !copy values from system 2
         shftcol=ncols1-ncom
         shftrow=nrows1         

         if(trans2)then
            do,j=1,ncols2
               do,i=1,nrows2
                  out%mat1(i+shftrow,j+shftcol)=work(2)%mat1(j,i)
               end do
            end do            
         else
            do,j=1,ncols2
               do,i=1,nrows2
                  out%mat1(i+shftrow,j+shftcol)=work(2)%mat1(i,j)
               end do
            end do            
         end if



      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('DIAVN___')! diagonal matrix
      select case(operation)
      case('+') ! add diagonal matrix to symmetric matrix
         do,i=1,work(1)%nval1
            do,j=1,work(1)%nval1
               work(1)%mat1(j,i)=W(1)*work(1)%mat1(j,i)
            end do
            !add diagonal
            work(1)%mat1(i,i)=work(1)%mat1(i,i)+W(2)*work(2)%pack1(i)
         end do
         

         !CONSTRUCT OUTPUT SYSTEM
         out%type='FULLSQVN' ! fully populated square matrix
         out%descr='BIN_MATH '//work(1)%type//' + '//work(2)%type
         out%nval1=work(1)%nval1
         out%nval2=out%nval1
         out%pval1=out%nval1**2 
         out%pval2=1
         out%nvec=0
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='F'
         out%side1_d=>work(1)%side1_d ! point to the right side description
         out%side2_d=>work(1)%side2_d
         out%mat1=>work(1)%mat1 ! point to right memory section

!         write(stderr,*)"added diagonal",associated(work(2)%pack1),associated(work(1)%mat1)
      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('VVVVVVVV')!VECTOR
      select case(operation)
      case('x')! MATRIX vector multiplication
         !CONSTRUCT OUTPUT SYSTEM
         out%type='VVVVVVVV' ! Vector system as output
         out%descr='BIN_MATH '//work(1)%type//' matrix vector product'

         out%nval1=nrows1 ! Size of the output vector
         out%nval2=out%nval1
!dummy values (unused)
         out%nvec=1 ! Amount of output vectors
         allocate(out%vec(out%nval1,out%nvec))
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%side1_d=>rows1 ! point to the right side description


         if(trans1)then
            uplo1='T'
         else
            uplo1='N'
         end if

         !do calculation
         call dgemv(uplo1,work(1)%nval1,work(1)%nval2,W(1),work(1)%mat1(1,1)&
              ,size(work(1)%mat1,1),work(2)%vec(:,v2select),1,0.d0,out%vec(:,1),1)


      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case('SSSSSSSS')
      select case(operation)
      case('c')!Condition number
         !calculate matrix L1 norm
         anorm=0.d0

         do,i=1,work(1)%nval2 ! loop over columns
            dumdb=0.d0
            do,j=1,work(1)%nval1 !loop over rows
               dumdb=dumdb+abs(work(1)%mat1(j,i))
            end do
            anorm=max(anorm,dumdb)
         end do
         
         !LU decomposition
         info=0
         allocate(ipiv(min(work(1)%nval1,work(1)%nval2)))

         write(stderr,*)"calculating LU factorization"
         call dgetrf(work(1)%nval1,work(1)%nval2,work(1)%mat1&
              ,work(1)%nval1,ipiv,info)
         if(info .ne. 0)then
            write(stderr,*)"ERROR LU decomposition failed,info:",info
            stop
         end if
         



         allocate(workv(4*work(1)%nval2))
         allocate(iwork(work(1)%nval1))
         write(stderr,*)"Calculation of the condition number"
         call dgecon('1',work(1)%nval2,work(1)%mat1(1,1),work(1)%nval1,anorm,&
              rcond,workv,iwork,info)

         write(*,*)"Condition number, 1-NORM"
         write(*,*)1/rcond,anorm
         stop !no further output

      case('k') ! rank k operation
         !CONSTRUCT OUTPUT SYSTEM
         out%type='SYMVN___' ! fully symmetric matrix
         out%descr='BIN_MATH '//work(1)%type//' Rank K operation matrix'
         if(trans1)then
            uplo1='N'
            out%nval1=work(1)%nval1
            out%side1_d=>work(1)%side1_d ! point to the right side description
         else
            uplo1='T'
            out%nval1=work(1)%nval2
            out%side1_d=>work(1)%side2_d ! point to the right side description
         end if

         out%nval2=out%nval1
         out%pval1=out%nval1*(out%nval1+1)/2
         out%pval2=1
         out%nvec=0 
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%mtyp='U'

         
         out%side2_d=>out%side1_d ! second side is within the same parameter space
         allocate(out%mat1(out%nval1,out%nval2))! allocate new memory

         call dsyrk(out%mtyp,uplo1,out%nval1,work(1)%nval1,W(1),&
                          work(1)%mat1(1,1),work(1)%nval1,0.d0,out%mat1(1,1),out%nval1)


      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case default
      write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
      stop
   end select
case('DIAVN___')! diagonal matrix
   write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
   stop

case('BDFULLVN')
   select case(work(2)%type)
      case('VVVVVVVV')!VECTOR
      select case(operation)
      case('x')! Block diagonal matrix vector multiplication
         !CONSTRUCT OUTPUT SYSTEM
         out%type='VVVVVVVV' ! Vector system as output
         out%descr='BIN_MATH '//work(1)%type//' matrix vector product'

         out%nval1=nrows1 ! Size of the output vector
         out%nval2=out%nval1
!dummy values (unused)
         out%nvec=1 ! Amount of output vectors
         allocate(out%vec(out%nval1,out%nvec))
         out%nread=0
         out%nint=0
         out%ndbls=0
         out%side1_d=>rows1 ! point to the right side description


         if(trans1)then
            uplo1='T'
         else
            uplo1='N'
         end if
         !loop over the diagonal blocks
         ipack=0 !index within packed vector
         ivec=0!index within the vector
         do,i=1,work(1)%nblocks
            sz=work(1)%blockind(i)-ivec ! size of the block
            call dgemv(uplo1,sz,sz,W(1),work(1)%pack1(ipack+1)&
              ,sz,work(2)%vec(ivec+1,v2select),1,0.d0,out%vec(ivec+1,1),1)
            !update indices
            ipack=ipack+sz**2
            ivec=ivec+sz
         end do
      case default
         write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
         stop
      end select
   case default
      write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
      stop

   end select
case('TRIVN___')!triangular matrix
   write(stderr,*)"ERROR, not supported (yet): ",work(1)%type,operation,work(2)%type
   stop
end select


!write result to file/standard output
if(out%type .eq. 'VVVVVVVV')then ! write vector output as ascii columns
   if(out%file .ne. 'stdout')then
      open(file=out%file,form='formatted',unit=unit)
   else
      unit=6 ! standard output unit
   end if

   do,i=1,out%nval1
      write(unit=unit,fmt='(A24,1x,G21.14)')out%side1_d(i),out%vec(i,1)

      ! write(unit=unit,fmt='(A24)',advance='NO')out%side1_d(i)
      ! write(unit=unit,fmt='(1x,G21.14)')out%vec(i,1)
   end do
   if(out%file .ne. 'stdout')close(unit)
else ! write to binary file ( everything at once)
   call write_BINtype(out) 

end if



end program BIN_MATH

subroutine help(version)
character(*)::version
character(8)::frmt
integer::stderr
stderr=0! standard error unit
frmt='(A)'
write(stderr,frmt)"Program BIN_MATH version:"//trim(version)
write(stderr,frmt)"Perform low level matrix, matrix-matrix, matrix-vector operations on dense  binary/ascii matrix files"
write(stderr,frmt)"Usage BIN_MATH [OPTIONS] FILE1 [FILE2]"
write(stderr,frmt)" One of  the input files may be read from standard input when its name is 'stdin' "
write(stderr,frmt)"OPTIONS may be:"
write(stderr,frmt)"-O=MATHOP: perform a mathemathical operation on the input:"
write(stderr,frmt)"   MATHOP may be:"
write(stderr,frmt)"   +: Add two matrices OUT=FILE1.mat+FILE2.mat (note subtraction may be achieved by using -W=1.d0/-1.d0)"
write(stderr,frmt)"   s: Solve for matrix or vectors: OUT  FILE1.mat * OUT = FILE2.mat or FILE1.mat * OUT = FILE2.vec"
write(stderr,frmt)"   x: Multiplication OUT=FILE1.MAT*FILE2.mat or FILE1.mat*FILE2.vec (in the case of -V)"
write(stderr,frmt)"   e: Eigenvalue/svd decomposition, calculates eigenvalues and eigenvectors of a matrix"
write(stderr,frmt)"      Or in the case of nonsymmetric matrix the singular value decomposition (outputs 2 files)"
write(stderr,frmt)"   c: Condition number (L1 norm) of a matrix (prints the condition number and L1 norm to standard output)"
write(stderr,frmt)"   C: Correlation matrix of a covariance matrix (only symmetric input matrix allowed)"
write(stderr,frmt)"   p: Propagate (symmetric) matrix FILE1.mat through FILE2.mat: OUT= FILE2.mat' * FILE1.mat * FILE2.mat"
write(stderr,frmt)"      The output is again symmetric"
write(stderr,frmt)"   P: Patch two matrices blockwise together (put the second one underneath the first)."
write(stderr,frmt)"    No parameters will be summed, gaps are filled with zeros"
write(stderr,frmt)"   k: Perform a symmetric rank K operation: C= FILE1.mat'*FILE1.mat (only one input system is allowed)" 
write(stderr,frmt)"-a: FILE2 is a vector in ascii format (2 columns 24 character wide parameter name and a value) "
write(stderr,frmt)"-f=OUTPUTFILENAME saves the result to file OUTPUTFILENAME (default prints to standard output)"
write(stderr,frmt)"-W=weight1/weight2 :scale the input with a constant factor weight1 (file1) weight2 file2"
write(stderr,frmt)"-S=SRCHSTRT/SRCHND: Compare the parameter string only for the substring ranging from SCRHSTRT to SRCHND"
write(stderr,frmt)"-V=V2SELECT: use a vector of system 2. Column of the vector matrix is denoted by V2SELECT."
write(stderr,frmt)"-t1: Transpose matrix of system 1"
write(stderr,frmt)"-t2: Transpose matrix of system 2"
write(stderr,frmt)"A note on the parameter space: The parameter space of the first matrix is used."
write(stderr,frmt)"Files may be read from standard input (specify stdin as file name)"
stop
end subroutine help


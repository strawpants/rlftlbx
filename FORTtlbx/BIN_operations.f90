!collection of routines which can be applied to BINdat matrix structures
!to be added:
!1) Triangular solve
!2) free up memory routine
!3) extract/put  a subset from/in a system

!!Updated by Roelof Rietbroek, Thu Mar 17 13:28:00 2011
!! added better replacement method for side descriptors (replacement can have a different size)
!!Updated by Roelof Rietbroek, Wed Feb 22 13:10:01 2012
!! fixed bug in aepxsion of diagonal matrices
!!Updated by Roelof Rietbroek, Fri Mar 15 10:40:57 2013
!! also expand symmetric matrices into full matrices

module BIN_operations
use binfiletools
contains
  
!subroutine to expand a BINdatsystem to full or symmetric dense matrix
subroutine BIN_expand(input)
implicit none
type(BINdat)::input  
integer::stderr
integer::i,j
integer*8::ind
stderr=0

if(input%mtyp .ne. 'P')then
   write(stderr,*)"WARNING from BIN_expand: matrix is already unpacked"
   return
end if

!allocate dense matrix ( if already associated it will be overwritten)
if(.not. associated(input%mat1))allocate(input%mat1(input%nval1,input%nval2))
input%mat1(1:input%nval1,1:input%nval2)=0.d0 ! initialize to zero

select case(input%type)
case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')

   !unpack in upper part of the dense matrix
   do,i=1,input%nval2 ! loop over columns
      ind=(i*(i-1))/2
      input%mat1(1:i,i)=input%pack1(ind+1:ind+i)
   end do

   !also mirror the matrix
   do,i=1,input%nval2 ! loop over columns
      do,j=1,i-1 !loop untill the first up diagonal
         input%mat1(i,j)=input%mat1(j,i)
      end do
   end do
   
   input%mtyp='F'
   input%type='FULL2DVN' ! change type to FULL2DVN
   input%pval1=input%nval1*input%nval2

case('BDSYMV0_','BDSYMVN_')
   !unpack in upper part of the dense matrix
   do,i=1,input%nval2 ! loop over columns
      do,j=1,i ! loop over rows
         ind=packindex(input,j,i)
         if(ind .ne. 0)input%mat1(j,i)=input%pack1(ind)
      end do
   end do
   input%mtyp='U'
   input%type='SYMVN___' ! change type to SYMMETRIC 
   input%pval1=(input%nval1*(input%nval1+1))/2
case('BDFULLV0','BDFULLVN')
   !unpack in the dense matrix
   do,i=1,input%nval2 ! loop over columns
      do,j=1,input%nval1 ! loop over rows
         ind=packindex(input,j,i)
         if(ind .ne. 0)input%mat1(j,i)=input%pack1(ind)
      end do
   end do
   input%mtyp='F'
   input%type='FULLSQVN' ! change type to FULLSQVN
   input%pval1=input%nval1*input%nval1
case('FULLSQV0','FULLSQVN','FULL2DVN')

   do,i=1,input%nval2 ! loop over columns
      do,j=1,input%nval1 ! loop over rows
         ind=(i-1)*input%nval1+j
         input%mat1(j,i)=input%pack1(ind)
      end do
   end do
   input%mtyp='F'
case('DIAVN___')!diagonal matrix
   do,i=1,input%nval1
      input%mat1(i,i)=input%pack1(i)
   end do
   input%mtyp='U'
   input%type='SYMVN___' ! change type to SYMMETRIC
   input%nval2=input%nval1
   input%pval1=input%nval1*(input%nval1+1)/2
end select

  
!deallocate packed matrix
deallocate(input%pack1)
end subroutine BIN_expand

!subroutine to put meta data in a system 'in'
!updated 13 Aug 2013 ( added optio to left adjust parameter description)
!!Updated by Roelof Rietbroek, Mon Aug 26 21:31:48 2013
!! increased length of dum2 string

subroutine BIN_putmeta(in,dum,iint,idbl)
use forttlbx
implicit none
type(BINdat)::in
character(*)::dum ! input string
logical::app
integer::ind,strlen,strlen1,strlen2,stderr,i
integer*8,optional,intent(in)::iint
double precision,optional,intent(in)::idbl
character(50)::strbefore,strafter,strtmp
character(500)::dum2
integer::last,unit
character(24),pointer::side(:)=>null()
integer::sz
stderr=0
last=0
unit=-1 ! negative if no file is open

!check whether data must be read from a file
if(dum(1:2) .eq. 'F=')then ! read metadata from file
   unit=14
   !open file
   open(unit=unit,file=trim(dum(3:)),form='formatted')
else
   dum2=trim(dum) ! just copy string
end if


do while (last .eq. 0) ! Loop over file ( or at least once)
   if(unit>0)read(unit=unit,iostat=last,fmt='(A500)')dum2 !read next entry from file
   if(last .ne. 0)exit ! exit upon end of file
   
   select case(dum2(1:2))
   case('d=')!update descr
      in%descr=trim(dum2(3:))
   case('I=')!update or append integer data
      ind=index(dum2,'/')
      if(ind .eq. 0)then
         write(stderr,*)"ERROR: error processing input string",trim(dum2)
         stop
      end if
      app=.true.
      do,i=1,in%nint ! seacrh if the parameter is already there
         if(dum2(3:ind-1) .eq. in%ints_d(i))then
            if(present(iint))then
               in%ints(i)=iint
            else
               read(dum2(ind+1:),*)in%ints(i)
            end if
            app=.false.
            exit
         end if
      end do
      
      if(app)then ! append to vector if not found
         !increase data vector by 1
         call realloc_ptr(in%ints_d,1)
         call realloc_ptr(in%ints,1)
         in%nint=in%nint+1
         in%ints_d(in%nint)=dum2(3:ind-1) ! description
         if(present(iint))then
            in%ints(in%nint)=iint
         else
            read(dum2(ind+1:),*)in%ints(in%nint) !value
         end if
      end if
   case('I-')!delete integer meta data
      in%nint=0
      deallocate(in%ints_d,in%ints)
   case('D=')!update or append double precision data
      ind=index(dum2,'/')
      if(ind .eq. 0)then
         write(stderr,*)"ERROR: error processing input string:",trim(dum2)
         stop
      end if
      app=.true.
      do,i=1,in%ndbls ! seacrh if the parameter is already there
         if(dum2(3:ind-1) .eq. in%dbls_d(i))then
            if(present(idbl))then
               in%dbls(i)=idbl
            else
               read(dum2(ind+1:),*)in%dbls(i)
            end if
            app=.false.
            exit
         end if
      end do
      
      if(app)then ! append to vector if not found
         !increase data vector by 1
         call realloc_ptr(in%dbls_d,1)
         call realloc_ptr(in%dbls,1)
         in%ndbls=in%ndbls+1
         in%dbls_d(in%ndbls)=dum2(3:ind-1) ! description
         if(present(idbl))then
            in%dbls(in%ndbls)=idbl
         else
            read(dum2(ind+1:),*)in%dbls(in%ndbls) !value
         end if
      end if
   case('D-')!delete double meta data
      in%ndbls=0
      deallocate(in%dbls_d,in%dbls)
   case('R=')!append Readme
      call realloc_ptr(in%readme,1)
      in%nread=in%nread+1
      in%readme(in%nread)=trim(dum2(3:))
   case('R-')!delete double meta data
      in%nread=0
      if(associated(in%readme))deallocate(in%readme)
   case('P1','P2')!replace strings for parameters
      ind=index(dum2,'/')
      if(ind .eq. 0)then
         write(stderr,*)"ERROR: error processing input string:",trim(dum2)
         stop
      end if
      if(dum2(3:3) .ne. '=')then
         write(stderr,*)"ERROR: error processing input string:",trim(dum2)
         stop
      end if

      strbefore=''
      strafter=''
      strbefore=dum2(4:ind-1)
      strafter=dum2(ind+1:ind+24)
      strlen1=ind-5
      strlen2=len_trim(strafter)
!      write(0,*)strlen1,strlen2,strbefore,strafter
      if(strlen1 > 24 .or. strlen2 >24 )then
         write(stderr,*)"BIN_operations ERROR: strings exceed 24 character limit:",trim(dum2)
         stop
      end if

      !check if side description vector is associated
      if(dum2(1:2) .eq. 'P1')then
         if(.not. associated(in%side1_d))then
            write(stderr,*)"ERROR: Side description 1 is not associated"
            stop
         end if
         !point to the correct side
         side=>in%side1_d
         sz=in%nval1
      else if(dum2(1:2) .eq. 'P2')then
         if(.not. associated(in%side2_d))then
            write(stderr,*)"ERROR: Side description 2 is not associated"
            stop
         end if
         !point to the correct side
         side=>in%side2_d
         sz=in%nval2
      end if

!replacement loop
      do,i=1,sz
         ind=index(side(i),strbefore(1:strlen1))
         if(ind .ne. 0)then
            !construct new string
            strtmp=''
            strtmp=side(i)(1:ind-1)//strafter(1:strlen2)//side(i)(ind+strlen1:24)
            side(i)=strtmp(1:24) !truncate result on the 24th character
         end if
      end do
      
      
   case('T1','T2')!tag parameters
      if(dum2(3:3) .ne. '=')then
         write(stderr,*)"ERROR: error processing input string:",trim(dum2)
         stop
      end if
      strlen=len_trim(dum2(4:))
      ind=24-strlen+1
      if(strlen >24)then
         write(stderr,*)"ERROR: STAMP is larger then descriptor:",trim(dum2)
         stop
      end if

      if(dum2(1:2) .eq. 'T1')then
         if(.not. associated(in%side1_d))then
            write(stderr,*)"ERROR: Side description 1 is not associated"
            stop
         end if
         do,i=1,in%nval1
            in%side1_d(i)(ind:ind+strlen-1)=dum2(4:4+strlen-1) ! replace string part
         end do
         
      else if(dum2(1:2) .eq. 'T2')then
         if(.not. associated(in%side2_d))then
            write(stderr,*)"ERROR: Side description 2 is not associated"
            stop
         end if
         do,i=1,in%nval2
            in%side2_d(i)(ind:ind+strlen-1)=dum2(4:4+strlen-1) ! replace string part
         end do
         
      end if
   case('L1')!left adjust parameters
      do i=1,in%nval1
         in%side1_d(i)=adjustl(in%side1_d(i))
      end do
   case('L2')!left adjust parameters
      if(.not. associated(in%side2_d))then
         write(stderr,*)"ERROR: Left adjust for side 2 requested but side description 2 is not associated"
         stop
      end if
      do i=1,in%nval2
         in%side2_d(i)=adjustl(in%side2_d(i))
      end do
   case default
      write(0,*)"ERROR: error processing input string",trim(dum2)
      stop
   end select
   if(unit<0)exit ! exit loop if not a file
end do ! end do wwhile loop
if(unit>0)close(unit)!close file if opened
end subroutine BIN_putmeta

!updated version of sort_BINtype
!also allow (some) packed structures to be permuted
!!Updated by Roelof Rietbroek, Wed Feb 22 14:18:44 2012
!! also allow symmetric systems to be permuted when perm2 and both are given

subroutine BIN_permute(input,perm1,perm2,both,dryrun)
use forttlbx
implicit none
type(BINdat)::input
integer,dimension(:),optional,target::perm1,perm2 ! optional sorting indices
logical,optional,intent(in)::both,dryrun
!private variables
logical::iboth,idry
integer::nrows,ncols,stderr
integer::i,j,k,sz,st,mxsz,err
integer*8::packind,shft
double precision,pointer::tmp(:,:)=>null()
integer,pointer::iperm1(:)=>null()
nrows=input%nval1
ncols=input%nval2

stderr=0

if(present(dryrun))then
   idry=dryrun
else
   idry=.false.
end if

if(present(both))then
   iboth=both ! permute both sides with the same vector
else
   iboth=.false.
end if


!input checks
if( iboth .and. (.not. present(perm1) .and. .not. present(perm2)))then
   write(stderr,*)"ERROR: BIN_permute: symmetric permutation required but no permutation vector given"
   stop
end if


if( iboth .and. present(perm1) .and. present(perm2))then ! impossible combo
   write(stderr,*)"ERROR: BIN_permute conflicting arguments"
   write(stderr,*)"       permutation vectors differ for sides"
   write(stderr,*)"       but symmetric permutation is requested"
      stop

end if


!check matrix system and 
select case(input%type)
case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')! symmetric
   !copy values in case of one sided permutation or 
   if( .not. iboth)then 
      write(stderr,*)"ERROR: BIN_permute Symmetric matrices may only have symmetric permutation"
      stop
   end if
   if(.not. idry)then
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
      case('P')!needs more memory
         allocate(tmp(input%nval1,input%nval1))
         input%mat1=>tmp
         do,i=1,input%nval2 ! loop over columns
            packind=(i*(i-1))/2
            input%mat1(1:i,i)=input%pack1(packind+1:packind+i)
         end do
         do,i=1,input%nval2
            do,j=1,i-1
               input%mat1(i,j)=input%mat1(j,i)
            end do
         end do
      
      case default
         write(stderr,*)"ERROR: permute_BIN, Matrix storage ",input%mtyp," is not supported for ",input%type
         stop
      end select
   end if
   if(present(perm1))iperm1=>perm1
   if(present(perm2))iperm1=>perm2

   !permute side description
   call permute(A=input%side1_d(1:nrows),perm=iperm1)
   
   

   if( .not. idry)then ! also permute vector and matrix
      if(input%nvec>0)call permute_rows(A=input%vec(1:nrows,1:input%nvec),perm=iperm1)
      call permute_rows(A=input%mat1(1:nrows,1:nrows),perm=iperm1)
      call permute_cols(A=input%mat1(1:nrows,1:nrows),perm=iperm1)
      

      if(input%mtyp .eq. 'P')then ! copy data back in packed array
         do,i=1,input%nval2 ! loop over columns
            packind=(i*(i-1))/2
            input%pack1(packind+1:packind+i)=input%mat1(1:i,i)
         end do
         input%mat1=>null()
         deallocate(tmp)!free up memory
      end if
   end if
   
case('FULLSQV0','FULLSQVN','FULL2DVN')! do nothing just check
   select case(input%mtyp)
   case('F')! do nothing accept
   case default
      write(stderr,*)"ERROR: BIN_permute, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select

   !permute side descriptions
   if(present(perm1))call permute(A=input%side1_d(1:nrows),perm=perm1)
   if(present(perm2))call permute(A=input%side2_d(1:ncols),perm=perm2)
   if( .not. idry)then ! also permute vector and matrix
      if(input%nvec>0 .and. present(perm1))call permute_rows(A=input%vec(1:nrows,1:input%nvec),perm=perm1)
      if(present(perm1))call permute_rows(A=input%mat1(1:nrows,1:ncols),perm=perm1)
      if(present(perm2))call permute_cols(A=input%mat1(1:nrows,1:ncols),perm=perm2)
   end if

case('DIAVN___')
   if( .not. iboth)then
      write(stderr,*)"ERROR: BIN_permute",input%type, "requires both sides to be permuted"
      stop
   end if

   select case(input%mtyp)
   case('P') ! do nothing just accept
   case default
      write(stderr,*)"ERROR: BIN_permute, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select
   if(present(perm1))iperm1=>perm1
   if(present(perm2))iperm1=>perm2
   !permute side descriptions
      call permute(A=input%side1_d(1:nrows),perm=iperm1)

   if( .not. idry)then ! also permute vector and matrix
      if(input%nvec>0)call permute_rows(A=input%vec(1:nrows,1:input%nvec),perm=iperm1)
      call permute(A=input%pack1(1:nrows),perm=iperm1)
   end if

case('BDSYMV0_','BDSYMVN_')
   if( .not. iboth)then
      write(stderr,*)"ERROR: BIN_permute",input%type, "requires both sides to be permuted"
      stop
   end if

   select case(input%mtyp)
   case('P') ! do nothing just accept
   case default
      write(stderr,*)"ERROR: BIN_permute, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select

   if(present(perm1))iperm1=>perm1
   if(present(perm2))iperm1=>perm2
   
   !try to permute the side descrption in the same way as would be done to the matrix
   ! if the permutation fails we know the BLOCK diagonal structure cannot be maintained with the permutation
   
   sz=0
   st=0
   do,i=1,input%nblocks
      sz=input%blockind(i)-st
      call permute(A=input%side1_d(st+1:st+sz),perm=iperm1(st+1:st+sz),err=err)
      if(err .ne. 0)then
         write(stderr,*)"ERROR: permutation violated sparse matrix structure:",input%type
         stop
      end if
      st=input%blockind(i)
      mxsz=max(sz,mxsz)
   end do
   
   !permute vector and matrix data 
   if( .not. idry)then
      !vector data is not sparse
      call permute_rows(A=input%vec(1:nrows,1:input%nvec),perm=iperm1)
!      allocate working array ( must hold largest block)
      allocate(tmp(mxsz,mxsz))
      sz=0
      st=0
      packind=0
      do,i=1,input%nblocks
         sz=input%blockind(i)-st
         
         !copy data in dense blockarray
         do,j=1,sz
            shft=packind+(j-1)*(j)/2
            tmp(1:j,j)=input%pack1(shft+1:shft+j)
         end do
         
         !mirror data in lower triangle
         do,j=1,sz
            do,k=1,j-1
               tmp(j,k)=tmp(k,j)
            end do
         end do

         !permute rows and column of the temporary block array
         
         call permute_rows(A=tmp(1:sz,1:sz),perm=iperm1(st+1:st+sz))
         call permute_cols(A=tmp(1:sz,1:sz),perm=iperm1(st+1:st+sz))
         
         !copy data back in dense blockarray
         do,j=1,sz
            shft=packind+(j-1)*(j)/2
            input%pack1(shft+1:shft+j)=tmp(1:j,j)
         end do

         st=input%blockind(i)
         packind=packind+sz*(sz+1)/2 ! update starting index in packed array
      end do

      deallocate(tmp)
   end if 
case('BDFULLV0','BDFULLVN')! do nothing type is accepted   
   select case(input%mtyp)
   case('P') ! do nothing just accept
   case default
      write(stderr,*)"ERROR: BIN_permute, Matrix storage ",input%mtyp," is not supported for ",input%type
      stop
   end select
   
   !try to permute the side descrption in the same way as would be done to the matrix
   ! if the permutation fails we know the BLOCK diagonal structure cannot be maintained with the permutation
   
   sz=0
   st=0
   do,i=1,input%nblocks
      sz=input%blockind(i)-st
      if(present(perm1))then
         call permute(A=input%side1_d(st+1:st+sz),perm=perm1(st+1:st+sz),err=err)
         if(err .ne. 0)then
            write(stderr,*)"ERROR: permutation of side1 violated sparse matrix structure:",input%type
            stop
         end if
      end if

      if(present(perm2))then
         call permute(A=input%side2_d(st+1:st+sz),perm=perm2(st+1:st+sz),err=err)
         if(err .ne. 0)then
            write(stderr,*)"ERROR: permutation of side2 violated sparse matrix structure:",input%type
            stop
         end if
      end if

      st=input%blockind(i)
      mxsz=max(sz,mxsz)
   end do
   
   !permute vector and matrix data 
   if( .not. idry)then
      !vector data is not sparse
      call permute_rows(A=input%vec(1:nrows,1:input%nvec),perm=perm1)
      !      allocate working array ( must hold largest block)
      allocate(tmp(mxsz,mxsz))
      sz=0
      st=0
      packind=0
      do,i=1,input%nblocks
         sz=input%blockind(i)-st
      

         !copy data in dense blockarray
         do,j=1,sz
            shft=packind+(j-1)*sz
            tmp(1:sz,j)=input%pack1(shft+1:shft+sz)
         end do
         
         !permute rows and column of the temporary block array
         
         if(present(perm1))call permute_rows(A=tmp(1:sz,1:sz),perm=perm1(st+1:st+sz))
         if(present(perm2))call permute_cols(A=tmp(1:sz,1:sz),perm=perm1(st+1:st+sz))
         
         !copy data back in packed array
         do,j=1,sz
            shft=packind+(j-1)*sz
            input%pack1(shft+1:shft+sz)=tmp(1:sz,j)
         end do
         
         st=input%blockind(i)
         packind=packind+sz**2 ! update starting index in packed array
      end do

      deallocate(tmp)
   end if 
case default ! else write error message
   write(stderr,*)"ERROR: BIN_permute, Matrix system not supported:",input%type
   stop
end select

end subroutine BIN_permute


!subroutine which performs a cholesky factorization on the matrix of a bindat structure
!uses lapack
subroutine BIN_cholesky(dat)
implicit none
type(BINdat)::dat
integer::info,stderr,i
integer::sz,st
integer*8::packind
stderr=0

select case(dat%type)
case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')! symmetric 
   select case(dat%mtyp)
   case('P')
      call dpptrf(dat%mtyp,dat%nval1,dat%pack1(1), info)
   case('U','L')
      call dpotrf(dat%mtyp, dat%nval1, dat%mat1(1,1), size(dat%mat1,1), info )
   case default
      write(stderr,*)"ERROR in BIN_cholesky: not supported:",dat%type,dat%mtyp   
      stop
   end select
   
   if(info .ne. 0)then
      write(stderr,*)"ERROR Cholesky decomposition failed"
      write(stderr,*)'        Cholesky factorization broke down on parameter:',info
      stop
   end if
!adapt output type
   dat%type='TRIVN___'
   
case('DIAVN___')
   select case(dat%mtyp)
   case('P')   
      info=0
      do,i=1,dat%nval1
         if(dat%pack1(i) <0)then
            info=-1
            exit
         end if
         dat%pack1(i)=sqrt(dat%pack1(i))
      end do
   case default
      write(stderr,*)"ERROR in BIN_cholesky: not supported:",dat%type,dat%mtyp   
      stop
   end select
case('BDSYMV0_','BDSYMVN_') ! block diagonal
   select case(dat%mtyp)
   case('P')
      sz=0
      st=0
      packind=0
      do,i=1,dat%nblocks
         sz=dat%blockind(i)-st
         
         !perform cholesky factorization on the packed block
         call dpptrf(dat%mtyp,dat%nval1,dat%pack1(packind+1), info)
         if(info .ne. 0)exit

         st=dat%blockind(i)
         packind=packind+sz*(sz+1)/2 ! update starting index in packed array
      end do

      if(info .ne. 0)then
         write(stderr,*)"ERROR Cholesky decomposition failed"
         write(stderr,*)'        Cholesky factorization broke down in block',i,'on parameter:',info
         stop
      end if
   !adapt output type
      dat%type='BDTRIVN_'
   case default
      write(stderr,*)"ERROR in BIN_cholesky: not supported:",dat%type,dat%mtyp   
      stop
   end select
   


case default
   write(stderr,*)"ERROR in BIN_cholesky: not supported:",dat%type,dat%mtyp   
   stop
end select


end subroutine BIN_cholesky


!routine to do a matrix multiplication of two matrix systems
!NOTE the sides must match!!
subroutine BIN_matmul(transa,transb,A,B,C,alpha,beta)
implicit none
type(BINdat),intent(inout)::A,B,C
character(1),intent(in)::transa,transb
double precision,intent(in)::alpha,beta


!private stuff
integer::i,j,stderr,rws,cls,m
character(9)::Atyp,Btyp
logical::bchanged,achanged
achanged=.false.
bchanged=.false.
rws=-1
cls=-1
m=-1

Atyp=A%type//A%mtyp
Btyp=B%type//B%mtyp
stderr=0

!determine the rows of the output
if(transa .eq. 'N')then !rows of A are the rows of the output
   rws=A%nval1
   m=A%nval2
   C%nval1=rws
   if(.not. associated(C%side1_d))allocate(C%side1_d(C%nval1))
   C%side1_d(1:rws)=A%side1_d(1:rws)
else
   rws=A%nval2
   m=A%nval1
   C%nval1=rws
   if(.not. associated(C%side1_d))allocate(C%side1_d(C%nval1))
   C%side1_d(1:rws)=A%side2_d(1:rws)
end if

!determine the ciolumn of the output
if(transb .eq. 'N')then !rows of A are the rows of the output
   cls=B%nval2
   if(m .ne. B%nval1)then
      write(stderr,*)"ERROR in BIN_matmul, dimensions do not match"
      stop
   end if
   
   C%nval2=cls
   if(.not. associated(C%side2_d))allocate(C%side2_d(C%nval2))
   C%side2_d(1:cls)=B%side2_d(1:cls)
else
   cls=B%nval1
   if(m .ne. B%nval2)then
      write(stderr,*)"ERROR in BIN_matmul, dimensions do not match"
      stop
   end if
   
   C%nval2=cls
   if(.not. associated(C%side2_d))allocate(C%side2_d(C%nval2))
   C%side2_d(1:cls)=B%side1_d(1:cls)
end if



!set remainder of the output type
C%pval1=C%nval1*C%nval2
C%pval2=1
C%type='FULL2DVN'
C%mtyp='F'



select case(Atyp)
case('SYMV0___U','SYMVN___U','SYMV1___U','SYMV2___U',&
     'SYMV0___L','SYMVN___L','SYMV1___L','SYMV2___L') !A is symmetric and dense
   select case(Btyp)
   case('DIAVN___P')!other matrix is diagonal (in packed form)
      if( .not. associated(C%mat1))then ! use a pointer to the original B matrix
         C%mat1=>A%mat1
         C%ldm=A%ldm
      end if
      if(associated(C%mat1,A%mat1))Achanged=.true. ! A matrix will be changed on output (C and A share the same memory)

      !mirror input matrix A
      select case(A%mtyp)
      case('U')
         forall(i=1:A%nval1,j=1:A%nval1,j<i)A%mat1(i,j)=A%mat1(j,i)
      case('L')
         forall(i=1:A%nval1,j=1:A%nval1,j<i)A%mat1(j,i)=A%mat1(i,j)
      end select

      !do matrix multiplication (scale columns) and apply scales
      forall(i=1:cls,j=1:rws)C%mat1(j,i)=beta*C%mat1(j,i)+&
           alpha*A%mat1(j,i)*B%pack1(i)

   case('FULLSQV0F','FULLSQVNF','FULL2DVF')!symmetric matrix times a full matrix
      if( transb .eq. 'T')then
         write(stderr,*)"ERROR in BIN_matmul: B may not be transposed in the combination:",A%type,A%mtyp,B%type,B%mtyp
         stop
      end if
      if(.not. associated(C%mat1))then
         allocate(C%mat1(rws,cls))
         C%ldm=rws
      end if

      call  dsymm('L',A%type,rws,cls,alpha,A%mat1(1,1),A%ldm,B%mat1(1,1)&
           ,B%ldm,beta,C%mat1(1,1),C%ldm)
   case('TRIVN___U','TRIVN___L')!B is a triangular matrix
      !mirror input matrix A
      select case(A%mtyp)
      case('U')
         forall(i=1:A%nval1,j=1:A%nval1,j<i)A%mat1(i,j)=A%mat1(j,i)
      case('L')
         forall(i=1:A%nval1,j=1:A%nval1,j<i)A%mat1(j,i)=A%mat1(i,j)
      end select
      
      if(.not. associated(C%mat1))then
         !point to original memory of A
         C%mat1=>A%mat1
         C%ldm=A%ldm
      end if
      
      if(associated(C%mat1,A%mat1))then
         achanged=.true.
         if(beta .ne. 0.d0)then
            write(stderr,*)"ERROR in BIN_matmul:beta!=0 but an update is not allowed"
            stop     
         end if
      else !copy data in new array
         C%mat1=A%mat1
      end if
      
      !triangular matrix multiplication
      call dtrmm('R',transb,'N',rws,cls,alpha,B%mat1(1,1),B%ldm,C%mat1(1,1),C%ldm)
      if(beta .ne. 0.d0)then !add original matrix times a scale
         C%mat1=C%mat1+beta*A%mat1
      end if
   case default
      write(stderr,*)"ERROR in BIN_matmul: unsupported combination:",A%type,A%mtyp,B%type,B%mtyp
      stop
   end select
case('FULLSQV0F','FULLSQVNF','FULL2DVF') !A has a full matrix
   select case(Btyp)
   case('SYMV0___U','SYMVN___U','SYMV1___U','SYMV2___U',&
        'SYMV0___L','SYMVN___L','SYMV1___L','SYMV2___L') !B is symmetric and dense
      if( transa .eq. 'T')then
         write(stderr,*)"ERROR in BIN_matmul: A may not be transposed in the combination:",A%type,A%mtyp,B%type,B%mtyp
         stop
      end if

      if(.not. associated(C%mat1))then
         allocate(C%mat1(rws,cls))
         C%ldm=rws
      end if

      call  dsymm('R',B%type,rws,cls,alpha,B%mat1(1,1),B%ldm,A%mat1(1,1)&
           ,A%ldm,beta,C%mat1(1,1),C%ldm)


   case('DIAVN___P') ! full matrix times a diagonal matrix
      if(transa .eq. 'T')then
         write(stderr,*)"ERROR in BIN_matmul: A may not be transposed in the combination:",A%type,A%mtyp,B%type,B%mtyp
         stop
      end if
      if( .not. associated(C%mat1))then ! use a pointer to the original A matrix
         C%mat1=>A%mat1
         C%ldm=A%ldm
      end if
      if(associated(C%mat1,A%mat1))Achanged=.true. ! A matrix will be changed on output (C and A share the same memory)

      !do matrix multiplication (scale columns) and apply scales
      forall(i=1:cls,j=1:rws)C%mat1(j,i)=beta*C%mat1(j,i)+&
           alpha*A%mat1(j,i)*B%pack1(i)

   case('FULLSQV0F','FULLSQVNF','FULL2DVF') !B has a full matrix
      if( .not. associated(C%mat1))then 
         allocate(C%mat1(rws,cls))
         C%ldm=rws
      end if
      
      call dgemm(transa,transb,rws,cls,m,alpha,A%mat1(1,1),A%ldm,&
           B%mat1(1,1),B%ldm,beta,C%mat1(1,1),C%ldm)
   case('TRIVN___U','TRIVN___L')!B is a triangular matrix
      if(transa .eq. 'T')then
         write(stderr,*)"ERROR in BIN_matmul: A may not be transposed in the combination:",A%type,A%mtyp,B%type,B%mtyp
         stop
      end if
      if(.not. associated(C%mat1))then
         !point to original memory of A
         C%mat1=>A%mat1
         C%ldm=A%ldm
      end if
      if(associated(C%mat1,A%mat1))then
         achanged=.true.
         if(beta .ne. 0.d0)then
            write(stderr,*)"ERROR in BIN_matmul:beta!=0 but an update is not allowed"
            stop     
         end if
      else !copy data in new array
         C%mat1=A%mat1
      end if
      
      !triangular matrix multiplication
      call dtrmm('R',transb,'N',rws,cls,alpha,B%mat1(1,1),B%ldm,C%mat1(1,1),C%ldm)
      if(beta .ne. 0.d0)then !add original matrix times a scale
         C%mat1=C%mat1+beta*A%mat1
      end if
   case default
      write(stderr,*)"ERROR in BIN_matmul: unsupported combination:",A%type,A%mtyp,B%type,B%mtyp
      stop
   end select

case('TRIVN___U','TRIVN___L')!A is a triangular matrix
   select case(Btyp)
      case('SYMV0___U','SYMVN___U','SYMV1___U','SYMV2___U',&
           'SYMV0___L','SYMVN___L','SYMV1___L','SYMV2___L') !B is symmetric and dense
         !mirror input matrix B
         select case(B%mtyp)
         case('U')
            forall(i=1:B%nval1,j=1:B%nval1,j<i)B%mat1(i,j)=B%mat1(j,i)
         case('L')
            forall(i=1:B%nval1,j=1:B%nval1,j<i)B%mat1(j,i)=B%mat1(i,j)
         end select
         
         if(.not. associated(C%mat1))then
            !point to original memory of B
            C%mat1=>B%mat1
            C%ldm=B%ldm
         end if
         if(associated(C%mat1,B%mat1))then
            bchanged=.true.
            if(beta .ne. 0.d0)then
               write(stderr,*)"ERROR in BIN_matmul:beta!=0 but an update is not allowed"
               stop     
            end if
         else !copy data in new array
            C%mat1=B%mat1
         end if

         !triangular matrix multiplication
         call dtrmm('L',transa,'N',rws,cls,alpha,A%mat1(1,1),A%ldm,C%mat1(1,1),C%ldm)
         if(beta .ne. 0.d0)then !add original matrix
            C%mat1=C%mat1+beta*B%mat1
         end if
      case('FULLSQV0F','FULLSQVNF','FULL2DVF') !B has a full matrix
         if(transb .eq. 'T')then
            write(stderr,*)"ERROR in BIN_matmul: B may not be transposed in the combination:",A%type,A%mtyp,B%type,B%mtyp
            stop
         end if
         if(.not. associated(C%mat1))then
            !point to original memory of B
            C%mat1=>B%mat1
            C%ldm=B%ldm
         end if
         if(associated(C%mat1,B%mat1))then
            bchanged=.true.
            if(beta .ne. 0.d0)then
               write(stderr,*)"ERROR in BIN_matmul:beta!=0 but an update is not allowed"
               stop     
            end if
         else !copy data in new array
            C%mat1=B%mat1
         end if

         !triangular matrix multiplication
         call dtrmm('L',transa,'N',rws,cls,alpha,A%mat1(1,1),A%ldm,C%mat1(1,1),C%ldm)
         if(beta .ne. 0.d0)then !add original matrix
            C%mat1=C%mat1+beta*B%mat1
         end if
      case default
         write(stderr,*)"ERROR in BIN_matmul: unsupported combination:",A%type,A%mtyp,B%type,B%mtyp
      stop
   end select
case default
   write(stderr,*)"ERROR in BIN_matmul: unsupported type:",A%type
   stop
end select

if(achanged)write(stderr,*)"WARNING in BIN_matmul: A has changed"
if(bchanged)write(stderr,*)"WARNING in BIN_matmul: B has changed"


end subroutine BIN_matmul

!routine to retrieve or set the diagonal of a matrix system 
subroutine BIN_getputdiag(dat,diag,put)
implicit none
type(BINdat),intent(inout)::dat
double precision::diag(:)
logical,intent(in)::put
integer::i,j,stderr,sz,st,packind
stderr=0

select case(dat%type)
case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')! symmetric 
   select case(dat%mtyp)
   case('P')
      if(put)then
         do,i=1,dat%nval1
            dat%pack1(i*(i+1)/2)=diag(i)
         end do
      else
         do,i=1,dat%nval1
            diag(i)=dat%pack1(i*(i+1)/2)
         end do
      end if
   case('U','L')
      if(put)then
         do,i=1,dat%nval1
            dat%mat1(i,i)=diag(i)
         end do
      else
         do,i=1,dat%nval1
            diag(i)=dat%mat1(i,i)
         end do
      end if
   case default
      write(stderr,*)"ERROR in BIN_getputdiag: not supported:",dat%type,dat%mtyp   
      stop
   end select

case('DIAVN___')
   select case(dat%mtyp)
   case('P')   
      if(put)then
         dat%pack1(1:dat%nval1)=diag
      else
         diag=dat%pack1(1:dat%nval1)
      end if
   case default
      write(stderr,*)"ERROR in BIN_getputdiag: not supported:",dat%type,dat%mtyp   
      stop
   end select

case('BDSYMV0_','BDSYMVN_') ! block diagonal
   select case(dat%mtyp)
   case('P')
      sz=0
      st=0
      packind=0
      do,i=1,dat%nblocks
         sz=dat%blockind(i)-st
         if(put)then
            do,j=1,sz
               dat%pack1(packind+(j*(j+1))/2)=diag(st+j)
            end do
         else
            do,j=1,sz
               diag(st+j)=dat%pack1(packind+(j*(j+1))/2)
            end do
         end if

         st=dat%blockind(i)
         packind=packind+sz*(sz+1)/2 ! update starting index in packed array
      end do
   case default
      write(stderr,*)"ERROR in BIN_getputdiag: not supported:",dat%type,dat%mtyp   
      stop
   end select

case('FULLSQV0','FULLSQVN','FULL2DVN')
   if(dat%nval1 .ne. dat%nval2)then
      write(stderr,*)"ERROR in BIN_getputdiag, matrix is not square"
      stop
   end if
   select case(dat%mtyp)
   case('F')
      if(put)then
         do,i=1,dat%nval1
            dat%mat1(i,i)=diag(i)
         end do
      else
         do,i=1,dat%nval1
            diag(i)=dat%mat1(i,i)
         end do
      end if
   case default
      write(stderr,*)"ERROR in BIN_getputdiag: not supported:",dat%type,dat%mtyp   
      stop
   end select
case default
   write(stderr,*)"ERROR in BIN_getputdiag: not supported:",dat%type,dat%mtyp   
   stop
end select



end subroutine BIN_getputdiag

!subroutine to extract a rectangular subset of a matrix
subroutine BIN_getsub(in,out,x,y,m,n,point,essentials)
implicit none
type(BINdat),intent(in)::in
type(BINdat),intent(out)::out
logical,intent(in)::point ! if true the data points tio the original memory locations
logical,intent(in)::essentials
integer,intent(in)::x,y,m,n
character(9)::intyp
integer::i,j,stderr

stderr=0
intyp=in%type//in%mtyp

out%type='FULL2DVN'
out%mtyp='F'
out%descr=in%descr
out%nval1=m
out%nval2=n
out%pval1=m*n
out%pval2=1
out%nvec=in%nvec
out%nread=in%nread
out%ndbls=in%ndbls
out%nint=in%nint


if(.not. essentials)then
   !double meta data
   if(in%ndbls> 0)then
      if(point)then
         out%dbls=>in%dbls
      else
         if(.not. associated(out%dbls))allocate(out%dbls(out%ndbls))
         out%dbls=in%dbls
      end if
   end if
   !integer meta data
   if(in%nint> 0)then
      if(point)then
         out%ints=>in%ints
      else
         if(.not. associated(out%ints))allocate(out%ints(out%nint))
         out%ints=in%ints
      end if
   end if
   
   !readme data
   if(in%nread> 0)then
      if(point)then
         out%readme=>in%readme
      else
         if(.not. associated(out%readme))allocate(out%readme(out%nread))
         out%ints=in%ints
      end if
   end if
end if

!vector data
if(in%nvec> 0)then
   if(point)then
      out%vec=>in%vec(x:x+out%nval1-1,:)
   else
      if(.not. associated(out%vec))allocate(out%vec(out%nval1,out%nvec))
      out%vec(1:out%nval1,1:out%nval2)=out%vec(x:x+out%nval1-1,:)
   end if
end if

!side description data
if(point)then
   out%side1_d=>in%side1_d(x:x+out%nval1-1)
   out%side2_d=>in%side2_d(x:x+out%nval2-1)
else
   if(.not. associated(out%side1_d))allocate(out%side1_d(out%nval1))
   out%side1_d(1:out%nval1)=in%side1_d(x:x+out%nval1-1)
   if(.not. associated(out%side2_d))allocate(out%side2_d(out%nval2))
   out%side2_d(1:out%nval2)=in%side2_d(y:y+out%nval2-1)
end if

!matrix data

select case(intyp)
   case('SYMV0___U','SYMVN___U','SYMV1___U','SYMV2___U',&
        'SYMV0___L','SYMVN___L','SYMV1___L','SYMV2___L') !input is symmetric and dense
      if(point)then
         out%mat1=>in%mat1(x:x+m-1,y:y+n-1)
         out%ldm=in%ldm !copy leading dimension of original matrix
      else
         if(.not. associated(out%mat1))allocate(out%mat1(m,n))
         out%mat1=in%mat1(x:x+m-1,y:y+n-1)
         out%ldm=size(out%mat1,1)
      end if


      !check whether the subset is square and aligned on the diagonal ( then the output is again symmetric)
      if(m .eq. n .and. x .eq. y)then
!          modify output type
         out%type='SYMVN___'
         out%mtyp=in%mtyp
         out%pval1=m*(m+1)/2
      else
      !possibly mirror part of the matrix falling in lower or upper part of the triangle
         select case(in%mtyp)
         case('U')
            do,i=y,y+n-1
               do,j=i+1,x+m-1 !loop in below diagonal
                  out%mat1(j-x+1,i-y+1)=in%mat1(i,j)
               end do
            end do
         case('L')
            do,i=y,y+n-1
               do,j=x,i-1 !loop above diagonal
                  out%mat1(j-x+1,i-y+1)=in%mat1(i,j)
               end do
            end do
         end select
      end if
   case('FULLSQV0F','FULLSQVNF','FULL2DVNF')! full matrix
      if(point)then
         out%mat1=>in%mat1(x:x+m-1,y:y+n-1)
         out%ldm=in%ldm !copy leading dimension of original matrix
      else
         if(.not. associated(out%mat1))allocate(out%mat1(m,n))
         out%mat1=in%mat1(x:x+m-1,y:y+n-1)
         out%ldm=size(out%mat1,1)
      end if
   case default
      write(stderr,*)"ERROR in BIN_getsub: not supported:",in%type,in%mtyp   
   stop
end select




end subroutine BIN_getsub



end module BIN_operations

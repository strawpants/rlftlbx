!!module file for FORTtlbx library
!!Coded by Roelof Rietbroek, Mon Oct 29 09:42:59 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

module FORTtlbx
implicit none


!!interface which select the appropriate function either:1 (matmul1d:righthand is vector)
!!2. matmul2d: right hand is a matrix
!!3. rithg hand is a vector and no transpose arguments supplied
!!4. right hand is a matrix and no transpose arguments supplied

!IMPORTANT: inputted matrices must be contigious in memory!!! No sections are allowed
interface matmulblas
   !full version
   function matmul2d(A,B,TRANSA,TRANSB)
     implicit none
     double precision,intent(in),dimension(:,:)::A,B
     integer,intent(in)::TRANSA,TRANSB !specify whether matrix should be transposed 1 is yes 0 is no
     double precision::matmul2d(size(A,(1+TRANSA)),size(B,(2-TRANSB))) !automatic array
   end function matmul2d
   
   !!function which replaces the generic matmul with a much quicker BLAS routine
   !matmul1d=A*vec
   function matmul1d(A,vec,TRANSA)
     implicit none
     double precision,intent(in),dimension(:,:)::A
     double precision,intent(in),dimension(:)::vec
     integer,intent(in)::TRANSA!specify whether matrix should be transposed 1 is yes 0 is no
     double precision::matmul1d(size(A,(1+TRANSA))) !automatic array
   end function matmul1d

!quick and dirty versions ( arrays are contiguously stored and no transpose is requested)

   function matmul1dclean(A,B)
     implicit none
     double precision, intent(in),dimension(:,:)::A
     double precision, intent(in),dimension(:)::B
     double precision::matmul1dclean(size(A,1))
   end function matmul1dclean

   function matmul2dclean(A,B)
     implicit none
     double precision, intent(in),dimension(:,:)::A,B
     double precision::matmul2dclean(size(A,1),size(B,2))
    
   end function matmul2dclean
end interface
interface sort_f90
   subroutine isort_f90(perm,ints,inv)
     implicit none
     integer,intent(in)::ints(:)
     integer,intent(out)::perm(:)
     logical,optional,intent(in)::inv
   end subroutine isort_f90

   subroutine dsort_f90(perm,dbls,inv)
     implicit none
     double precision,intent(in)::dbls(:)
     integer,intent(out)::perm(:)
     logical,optional,intent(in)::inv
   end subroutine dsort_f90
! to be added character version
end interface
!!interface which selects the appropriate function for the packed symmetric matrix routines
interface packmul
   function packmul2d(mpack,mat)
     implicit none
     double precision,intent(in),dimension(:)::mpack
     double precision,intent(in),dimension(:,:)::mat
     double precision,dimension(:,:)::packmul2d(size(mat,1),size(mat,2))
   end function packmul2d
   function packmul1d(mpack,vec)
     implicit none
     double precision,intent(in),dimension(:)::mpack,vec
     double precision,dimension(:)::packmul1d(size(vec,1))
   end function packmul1d
end interface

!permutation of various vector types
interface permute
   Subroutine permute_d(A,perm,inv,err)
     implicit none
     double precision,intent(inout)::A(:)
     integer,intent(inout)::perm(:)
     logical,optional,intent(in)::inv
     integer,optional,intent(out)::err
   end Subroutine permute_d

   Subroutine permute_i(A,perm,inv,err)
     implicit none
     integer,intent(inout)::A(:)
     integer,intent(inout)::perm(:)
     logical,optional,intent(in)::inv
     integer,optional,intent(out)::err
   end Subroutine permute_i

   Subroutine permute_ch(A,perm,inv,err)
     implicit none
     Character(*),intent(inout)::A(:)
     integer,intent(inout)::perm(:)
     logical,optional,intent(in)::inv
     integer,optional,intent(out)::err
   end Subroutine permute_ch

end interface

interface realloc_ptr
   
   !character version
   subroutine reallocate_cptr(in,nadd)
     implicit none
     character(*),pointer::in(:)
     integer,intent(in)::nadd
   end subroutine reallocate_cptr

   !double version
   subroutine reallocate_dptr(in,nadd)
     implicit none
     double precision,pointer::in(:)
     integer,intent(in)::nadd
   end subroutine reallocate_dptr

   !integer version
   subroutine reallocate_iptr(in,nadd)
     implicit none
     integer,pointer::in(:)
     integer,intent(in)::nadd
   end subroutine reallocate_iptr

   subroutine reallocate_liptr(in,nadd)
     implicit none
     integer*8,pointer::in(:)
     integer,intent(in)::nadd
   end subroutine reallocate_liptr
   subroutine reallocate_dptrm(in,xadd,yadd)
     implicit none
     double precision,pointer::in(:,:)
     integer,intent(in)::xadd,yadd
   end subroutine reallocate_dptrm
end interface

interface
   function dot_productblas(A,B)
     implicit none
     double precision,intent(in)::A(:),B(size(A,1)) !B and dot_productblas are automatic arrays
     double precision::dot_productblas
   end function dot_productblas

   ! strip directories and possibly the suffix from a path
   subroutine basenamef(stripped,path,suf)
     implicit none
     character(*),intent(out)::stripped
     character(*),intent(in)::path
     character(*),intent(in),optional::suf
   end subroutine basenamef

   !binary file manipulation routines
   function ltlend()
     logical::ltlend
   end function ltlend

   
   function freeunit()
     integer::freeunit
   end function freeunit

   !permutation functions


   ! Permutation of double Matrices
   Subroutine permute_rows(A,perm,inv,err)
     implicit none
     double precision,intent(inout)::A(:,:)
     integer,intent(inout)::perm(:)
     logical,optional,intent(in)::inv
     integer,optional,intent(out)::err
   end Subroutine permute_rows

   Subroutine permute_cols(A,perm,inv,err)
     implicit none
     double precision,intent(inout)::A(:,:)
     integer,intent(inout)::perm(:)
     logical,optional,intent(in)::inv
     integer,optional,intent(out)::err
   end Subroutine permute_cols

   subroutine permute_sym(A,perm,inv)
     implicit none
     double precision,intent(inout)::A(:,:)
     integer,intent(inout)::perm(:)
     logical,optional,intent(in)::inv
   end subroutine permute_sym


   subroutine write_sym_mat(file,mat,matpack,vec,list,lmax,lmin,time,&
        d2,d3,d4,descr,char2,char3,char4)
     implicit none
     character(*),intent(in)::file
     double precision,intent(in),dimension(:,:),optional::mat !either provide a full 2 dimesnional matrix 
     double precision,intent(in),dimension(:),optional::matpack !or a matrix already in packed form
     double precision,intent(in),dimension(:),optional::vec
     character(*),intent(in),optional,dimension(:)::list
     integer,intent(in),optional::lmax,lmin
     double precision,intent(in),optional::time,d2,d3,d4
     character(*),optional,intent(in)::descr,char2,char3,char4
   end subroutine write_sym_mat

   subroutine read_sym_mat(file,mat,matpack,typ,nval,vec,list,lmax,lmin,time,&
        d2,d3,d4,descr,char2,char3,char4)
     implicit none
     character(*),intent(in)::file
     double precision,intent(out),optional,dimension(:)::matpack,vec
     double precision,intent(out),optional,dimension(:,:)::mat
     integer,intent(out),optional::typ,nval,lmax,lmin
     character(*),intent(out),optional::descr,char2,char3,char4
     double precision,optional,intent(out)::time,d2,d3,d4
     character(*),optional,dimension(:),intent(out)::list
   end subroutine read_sym_mat

! !!read and write normal system version 2

!    subroutine read_sym_mat2(file,pipe,V,typ,descr,nint,ndoubles,nval,ints_d,ints,doubles_d,doubles,vec_d,vec,initvec,mat,matpack)
!      implicit none
!      character(*),intent(in)::file
!      double precision,intent(out),optional,dimension(:)::matpack,vec,initvec,doubles
!      double precision,intent(out),optional,dimension(:,:)::mat
!      integer,intent(out),optional,dimension(:)::ints
!      integer,intent(out),optional::typ,nint,ndoubles,nval
!      character(*),intent(out),optional::descr
!      character(24),optional,dimension(:),intent(out)::ints_d,doubles_d,vec_d
!      logical,intent(in),optional::pipe !logical which indicates that the file should only be accessed as a pipe(thus nor rewinding or reopening)
!      character(8),optional, intent(out)::V
!    end subroutine read_sym_mat2

!    subroutine write_sym_mat2(file,descr,ints_d,ints,doubles_d,doubles,vec_d,vec,initvec,mat,matpack)
!      implicit none
!      character(*),intent(in)::file
!      character(24),intent(in),optional::ints_d(:),doubles_d(:),vec_d(:)
!      character(*),intent(in),optional::descr
!      double precision,intent(in),dimension(:,:),optional::mat !either provide a full 2 dimesnional matrix 
!      double precision,intent(in),dimension(:),optional::matpack !or a matrix already in packed form
!      double precision,intent(in),dimension(:),optional::vec,initvec,doubles
!      integer,intent(in),optional::ints(:)
!    end subroutine write_sym_mat2


! !!!BLOCK DIAGONAL MATRIX read and write rouitnes
!    subroutine write_bd_mat2(file,descr,ints_d,ints,doubles_d,doubles,vec_d,bdmat,blockind)
!      implicit none
!      character(*),intent(in)::file
!      character(24),intent(in),optional::ints_d(:),doubles_d(:),vec_d(:)
!      character(*),intent(in),optional::descr
!      double precision,intent(in),dimension(:)::bdmat !block diagonal matrix in vector form 
!      double precision,intent(in),dimension(:),optional::doubles
!      integer,intent(in),optional::ints(:)
!      integer,intent(in)::blockind(:) !index matrix with end indices of the blocks
!    end subroutine write_bd_mat2

!    subroutine read_bd_mat2(file,pipe,V,typ,descr,nint,ndoubles,nval,nblocks,bdsize,ints_d,ints,doubles_d,doubles,vec_d,bdmat,blockind)
!      implicit none
!      character(*),intent(in)::file
!      double precision,intent(out),optional,dimension(:)::doubles
!      double precision,intent(out),optional,dimension(:)::bdmat
!      integer,intent(out),optional,dimension(:)::ints,blockind
!      integer,intent(out),optional::typ,nint,ndoubles,nval,nblocks,bdsize
!      character(*),intent(out),optional::descr
!      character(24),optional,dimension(:),intent(out)::ints_d,doubles_d,vec_d
!      logical,intent(in),optional::pipe !logical which indicates that the file should only be accessed as a pipe(thus nor rewinding or reopening)
!      character(8),optional, intent(out)::V
!    end subroutine read_bd_mat2


!!BIN FILE TOOLS general replaces upper routines


   Subroutine write_BIN(file,type,descr,blockind,ints_d,ints,dbls_d,dbls,side1_d,vec1,vec2,pmat1,mat1,stage,uplo) !(,side2_d,pmat2,mat2)
     implicit none
     !!general input
     character(*),intent(in)::file
     character(1),intent(in),optional::uplo
     character(8),intent(in)::type
     character(80),intent(in),optional::descr
     character(24),intent(in),optional::ints_d(:),dbls_d(:),side1_d(:)!,side2_d(:)
     integer,intent(in),optional::ints(:),stage !progress parameter indicates untill which state needs to be written
     double precision,intent(in),optional::dbls(:)
     

     !!type specific input
     integer,intent(in),optional::blockind(:)
     double precision,intent(in),optional::vec1(:),vec2(:),pmat1(:),mat1(:,:)!,pmat2(:),mat2(:,:)
   end Subroutine write_BIN
   
   Subroutine read_BIN(file,type,descr,nval1,nval2,pval1,pval2,nint,ndbls,nblocks,V&
        ,blockind,ints_d,ints,dbls_d,dbls,side1_d,vec1,vec2,pmat1,mat1,uplo)!side2_d,,pmat2,mat2)
     implicit none
     !!general in/output
     character(*),intent(in)::file
     character(1),intent(in),optional::uplo
     character(8),intent(out),optional::type
     character(80),intent(out),optional::descr
     character(24),intent(out),optional::ints_d(:),dbls_d(:),side1_d(:)!,side2_d(:)
     integer,intent(out),optional::ints(:)
     double precision,intent(out),optional::dbls(:)
     character(8),intent(out),optional::V
     
     !!type specific output
     integer,intent(out),optional::nval1,nval2,pval1,pval2,nint,ndbls,nblocks
     integer,intent(out),optional::blockind(:)
     double precision,intent(out),optional::vec1(:),vec2(:),pmat1(:),mat1(:,:)!,pmat2(:),mat2(:,:)
   end Subroutine read_BIN
   
  

!!!!!!!!!!READUCE TOOLS from NORM_fix_red.f90

   subroutine setapriori_norm(C,d,dx,ltpl,uplo)
     implicit none
     character(1),intent(in),optional::uplo
     double precision,intent(in)::C(:,:),dx(:)
     double precision,intent(inout)::d(:),ltpl
   end subroutine setapriori_norm

   subroutine setapriori_norm2(C,ldc,d,dx,ltpl,uplo)
     character(1),intent(in),optional::uplo
     integer,intent(in)::ldc
     double precision,intent(in)::dx(:),C(:,:)
     double precision,intent(inout)::d(:),ltpl
   end subroutine setapriori_norm2
   
   subroutine reduce_norm(C11,C12,C22,d1,d2,ltpl)
     implicit none
     double precision,intent(inout)::C11(:,:),C12(:,:),C22(:,:),d1(:),d2(:),ltpl
   end subroutine reduce_norm

   
   subroutine reduce_norm2(C11,ldc11,C12,ldc12,C22,ldc22,d1,d2,ltpl)
     integer,intent(in)::ldc11,ldc12,ldc22
     double precision,intent(inout)::C11(:,:),C12(:,:),C22(:,:),d1(:),d2(:),ltpl
   end subroutine reduce_norm2

   subroutine solve_norm(C,d,ltpl)
     implicit none
     double precision,intent(inout)::C(:,:),d(:),ltpl
   end subroutine solve_norm

   subroutine solve_norm2(C,ldc,n,d,ltpl)
     integer,intent(in)::ldc,n
     double precision,intent(inout)::C(:,:),d(n),ltpl
   end subroutine solve_norm2

   subroutine solve_tsvd(C,ldc,n,d,ltpl,ntrunc)
     integer,intent(in)::ldc,n
     double precision,intent(inout)::C(:,:),d(n),ltpl
     integer,intent(in)::ntrunc !number of eigenvalues to keep
   end subroutine solve_tsvd

   subroutine diag_transf(B,A,atb,x0,nunit,nb,add)
     implicit none
     integer,intent(in)::nb,nunit
     double precision,intent(in)::B(:)
     double precision,intent(inout),optional::atb(:),x0(:),A(:,:)
     logical,intent(in),optional::add
   end subroutine diag_transf

   subroutine full_transf(B,ldb,A,lda,atb,x0,nunit,nb,k,add)
     integer,intent(in)::ldb,lda,nunit,nb,k
     double precision,intent(inout)::B(:,:)
     double precision,intent(inout),optional::atb(:),x0(:),A(:,:)
     logical,intent(in),optional::add
   end subroutine full_transf

   subroutine sym_transf(B,ldb,A,lda,atb,x0,nunit,nb,add)
     integer,intent(in)::ldb,lda,nb,nunit
     double precision,intent(inout)::B(:,:)
     double precision,intent(inout),optional::atb(:),x0(:),A(:,:)
     logical,intent(in),optional::add
   end subroutine sym_transf

   subroutine trans_cov(B,ldb,A,lda,atb,ltpl,nunit,nb,type)
     implicit none
     integer,intent(in)::nunit,nb,lda,ldb,type
     double precision,intent(inout)::B(:,:)
     double precision,intent(inout)::A(:,:)
     double precision,intent(inout)::ltpl,atb(:)
   end subroutine trans_cov


!other functions
   pure function power(binvec)
     implicit none
     integer,intent(in)::binvec(:)
     integer::power
   end function power

   function binofpower(power,n)
     implicit none
     integer,intent(in)::power(:),n
     integer::binofpower(size(power))
   end function binofpower

   !subroutine check_normsys(same,typ,nsys,nsysnew,nolocals)
    ! implicit none
    ! integer,intent(in)::typ(:),nsys !type: decimal represenation of binary asssumed columns, nsys:amount of normal systems, ndat(:) amount of unknowns per system
    ! integer,intent(out)::same(:),nsysnew
    ! logical, intent(in)::nolocals !if true then local parameters are ignored in the comparison
   !end subroutine check_normsys

   function locals(typ,col)
     implicit none
     integer,intent(in)::typ(:)
     integer,intent(in),optional::col
     integer::locals(size(typ,1))
   end function locals

   subroutine cp_pack_2_square(packed,square,rowind)
     implicit none
     integer,intent(in)::rowind(:)
     double precision,intent(in)::packed(:)
     double precision,intent(out)::square(:,:)
   end subroutine cp_pack_2_square

   subroutine cp_pack_2_rect(packed,rect,rowind,colind)
     implicit none
     double precision,intent(in)::packed(:)
     double precision,intent(out)::rect(:,:)
     integer,intent(in)::rowind(:),colind(:)
   end subroutine cp_pack_2_rect
   
   subroutine cp_square_2_pack(square,packed,rowind,add)
     implicit none
     integer,intent(in)::rowind(:)
     double precision,intent(inout)::packed(:)!must not alter the other values thus intent(inout)
     double precision,intent(in)::square(:,:)
     logical,intent(in)::add
   end subroutine cp_square_2_pack

   subroutine cp_rect_2_pack(rect,packed,rowind,colind,add)
     implicit none
     integer,intent(in)::colind(:) !index vector with zero entries for undefined variables and 1 entries for defined variables
     integer,intent(in)::rowind(:)
     double precision,intent(inout)::packed(:)!must not alter the other values thus intent(inout)
     double precision,intent(in)::rect(:,:)
     logical,intent(in)::add
   end subroutine cp_rect_2_pack

   function vecind(binvec)
     implicit none
     integer, intent(in)::binvec(:) !assumed shape array 0 or 1 entries
     integer:: vecind(sum(binvec))
   end function vecind

   function subpackind(binvec)
     implicit none
     integer,intent(in)::binvec(:) !assumed shape array
     integer::subpackind((sum(binvec)*(sum(binvec)+1))/2)
   end function subpackind
   

   subroutine next_tree(treeold,treenew,typold,typnew,actionnew,nsysnew,method,overwrite)
     implicit none
     integer,intent(in)::treeold(:),typold(:),method
     integer,intent(out)::treenew(:),typnew(:),actionnew(:),nsysnew
     logical, intent(in)::overwrite
   end subroutine next_tree
   
   subroutine print_tree(tree,action,type,nsys,normfiles)
     implicit none
     integer,intent(in)::tree(:,:),action(:,:),type(:,:),nsys(:)
     character(*),intent(in)::normfiles(:)
   end subroutine print_tree

   subroutine reduce(A_cc,A_cl,A_ll,bc,bl)
     implicit none
     double precision,intent(inout)::A_cc(:,:),A_cl(:,:),A_ll(:,:)
     double precision,intent(inout)::bc(:),bl(:)
   end subroutine reduce

   subroutine solve_para(U,b,A_clUi,bc)
     implicit none
     double precision,intent(in)::U(:,:)
     double precision,intent(inout)::b(:)
     double precision,optional,intent(inout)::bc(:)
     double precision,intent(inout),optional::A_clUi(:,:)!present if local parameters need to be solved
   end subroutine solve_para
   
   subroutine prepare_covlocal(U_cc,A_clUU,iA_ll)
     implicit none
     double precision,intent(in)::U_cc(:,:)
     double precision,intent(inout)::A_clUU(:,:),iA_ll(:,:)!both changed on exit
   end subroutine prepare_covlocal

!!!SUBROUTINES to READ GRACE NORMAL SYSTEMS
   subroutine read_GRACEnormhead(filename,npar,npar_red,nobs,side_d,ltpl,apriori,Etime,Stime,hist,nhist)
     implicit none
     character(*),intent(in)::filename
     character(24),intent(out),optional::side_d(:)
     character(120),intent(out),optional::hist(:)
     double precision,intent(out),optional::apriori(:),ltpl,Etime,Stime
     integer,intent(out),optional::npar,npar_red,nobs,nhist
   end subroutine read_GRACEnormhead
   
   subroutine read_GRACEnormbin(filename,npar,mat,pmat,Atb,ltpl)
     implicit none
     integer,intent(in)::npar
     character(*),intent(in)::filename
     double precision,intent(out),optional::pmat(:),mat(:,:),ltpl,Atb(:)
   end subroutine read_GRACEnormbin

   function regcompf(regex,file)
     implicit none
     integer regcompf
     character(*)::regex
     logical,intent(in),optional::file
   end function regcompf

   function regexecf(nregex,string)
     implicit none
     character(*)::string
     integer::nregex
     logical::regexecf
   end function regexecf



   subroutine get_permvec(side1,side2,back,st,nd,perm,ncom)
     implicit none
     character(*),intent(in),dimension(:)::side1,side2 ! character description of the sides to be matched
     logical,intent(in)::back ! logical which determines to do with entries presetn in side1 but no tin side2
     !if back=.true. put those entries at the back of the permutation else in front of the permutation
     integer,intent(in)::st,nd ! start and end of the string part to compare 
     integer,intent(out)::perm(:) ! permutation vector for side1 to match ( ignoring entries present in side2 but not in side1) side2
     integer,intent(out)::ncom ! number of common parameters found
   end subroutine get_permvec

   subroutine get_permvec2(side1,side2,back,st,nd,perm1,perm2,ncom)
     character(*),intent(in),dimension(:)::side1,side2 ! character description of the sides to be matched
     integer,intent(in)::st,nd ! start and end of the string part to compare 
     integer,intent(out)::perm1(:)! permutation vector for side1 to match ( ignoring entries present in side2 but not in side1) side2 
     integer,intent(out)::perm2(:) !and for side 2 to match side1 (additional permutation of side1 may be needed)
     integer,intent(out)::ncom ! number of common parameters found
     logical,intent(in)::back ! logical which determines what to do with unique entries in side1  (stack at front or back)
   end subroutine get_permvec2

!statistical function
   function pval_stud_t(t,df)
     double precision::pval_stud_t
     double precision,intent(in)::t
     integer,intent(in)::df
   end function pval_stud_t


   subroutine interp1(x,y,xi,yi,dim,method)
     implicit none
     double precision,intent(in)::x(:),xi(:),y(:,:)
     integer,intent(in)::dim !dimension corresponding to the to be interpolated axis
     character(*),intent(in)::method
     double precision,intent(out),pointer::yi(:,:) ! will be (re) allocated within the routine
end subroutine interp1

end interface

end module FORTtlbx

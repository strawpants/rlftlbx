!!-------------------------------------------------------------------------------
!!
!!    SH_createDDK (adapted over the years, orginating from J.Kusche's cov_regul04 program)
!!    Roelof Rietbroek
!!    the program requires full normal matrix N in BINV format (see RLFtlbx), adds
!!    alpha times a regularization matrix R 
!!    and inverts and multiplies with N(Tikhonov-filter)
!!
!!    W = (N + aR)^-1 N
!!
!!    the R entries are defined as l**pow, where l is determined from
!!    mi,ma and assuming the order-degree ordering scheme
!!-------------------------------------------------------------------------------
program SH_createDDK
  use various
  use binfiletools
  use bin_operations
  implicit none
  integer::narg,lmax,lmin,pow,pos,cs,qm,ind,nend,info,nblocks,nint
  integer::i,nval1,pval1
  integer,allocatable::ints(:),blockind(:)
  character(24),allocatable::ints_d(:),side_d(:)
  character(24)::dbls_d(2)
 ! integer, allocatable,dimension(:,:)::indx
  double precision::alpha,mem,deg,dbls(2)
  double precision,allocatable,dimension(:,:)::artmp
  double precision,allocatable,dimension(:)::arpack
  double precision, allocatable,dimension(:)::rdiag,ngg
  character(200)::covfile,outfile,dum
  character(80)::descr
  character(8)::covtype,covtypout
  logical::pwr,blck,altfname

  integer::side,fll,pck,b,sz,j
  pwr=.false.
  altfname=.false.

  !default initialisations
  lmax=0
  lmin=0

  alpha=0.d0
  pow=0
  
  outfile=''
  blck=.false.

  !process command line arguments
  !amount of arguments
  narg=iargc()
  if(narg < 2)call help()
  !get the (approximate error covariance file name
  call getarg(narg,covfile)
  
  !now do a loop over the other arguments
  do,i=1,narg-1
     call getarg(i,dum)
     !call help when some options are no options
     if(dum(1:1).ne. '-')call help()

     select case(dum(2:2))
        case('f')
           outfile=trim(dum(3:))
           altfname=.true.
        case('p')
           pwr=.true.
           ind=index(dum,',')
           read(dum(3:ind-1),*)alpha
           read(dum(ind+1:),*)pow
        case default
           write(*,*)'Unknown option selected'
           call help()
     end select
  end do


  if(.not. pwr)then
     write(0,*)'ERROR: Specify power law parameters'
     call help()
  end if







!!now read in the error covariance data (symmetric in packed form)
  !inquire size and type
  call read_BIN(file=trim(covfile),type=covtype,nval1=nval1,pval1=pval1,nblocks=nblocks,nint=nint)
  

  allocate(ints(nint),ints_d(nint))
  allocate(ngg(pval1))
  allocate(side_d(nval1)) !allocate side description array

  if(covtype(1:8) .eq. 'BDSYMV0_')then !block diagonal matrix only
     blck=.true.
     allocate(blockind(nblocks))

     !get data
     call read_BIN(file=trim(covfile),blockind=blockind,ints_d=ints_d,ints=ints,side1_d=side_d,pmat1=ngg)
     !get lmax and lmin from meta data
     do,i=1,nint
        if(trim(ints_d(i)) .eq. 'Lmax')lmax=ints(i)
        if(trim(ints_d(i)) .eq. 'Lmin')lmin=ints(i)
        
     end do
     !set output type
     covtypout='BDFULLV0'
     
     write(*,*)'Found Block diagonal system'
     

  else if(covtype(1:8) .eq. 'SYMV0___')then !SYMMETRIC MATRIX

     !get data
     call read_BIN(file=trim(covfile),ints_d=ints_d,ints=ints,side1_d=side_d,pmat1=ngg)

     !get lmax and lmin from meta data
     do,i=1,nint
        if(trim(ints_d(i)) .eq. 'Lmax')lmax=ints(i)
        if(trim(ints_d(i)) .eq. 'Lmin')lmin=ints(i)
     end do
     
     !make the matrix one big diagonal block
     allocate(blockind(1))
     blockind(1)=nval1
     nblocks=1

     covtypout='FULLSQV0'
     
     write(*,*)'Found symmetric normal matrix'

  else
     write(0,*)'ERROR expecting a matrix of type BDSYMV0_ or SYMV0___ getting:',covtype        
     stop


  end if
!write(*,*)nval1,pval1,nblocks,nint
!do,i=1,nblocks
!write(*,*)blockind(i)
!end do

!!construct diagonal signal covariance matrix (diagional only)
  allocate(rdiag(nval1))
  rdiag=0.d0
  do,i=1,nval1
!     write(*,*)side_d(i)(5:7)
     read(side_d(i)(5:7),*)deg !read degree from description string
!     write(*,*)deg
     rdiag(i)=alpha*deg**pow

  end do


!   do,l=lmin,lmax
!      where (indx(:,1) .eq. l) rdiag=alpha*dble(l)**pow
!   end do



if(blck)then
   !get size of new block diagonal ar matrix not in packed form anymore
   side=0
   fll=0
   do,i=1,nblocks
      sz=blockind(i)-side !side size
      fll=fll+sz**2       !progression within the full matrix
      side=blockind(i)    !reset
   end do
!   write(*,*)'packed matrix size',fll
   !allocate packedfilter matrix
   allocate(arpack(fll))
   write(*,*)'filter packed matrix size',fll
   !also allocate a temporary square array (maximum size necessary is lmax-lmin+1)
   allocate(artmp(lmax-lmin+1,lmax-lmin+1))
   
else !no block diagonal approach
! allocate a temporary square array 
   allocate(artmp(nval1,nval1))
end if



!loop over diagonal blocks (only one if non block diagonal)
write(*,*)'calculating weight matrix, this can take a while'


side=0
fll=0
pck=0
do,b=1,nblocks
   if(blck)write(*,*)'Doing block: ',b,'of ',nblocks
   sz=blockind(b)-side
  
   !copy data from normal system to square array
   do,i=1,sz
      do,j=1,i
         artmp(i,j)=ngg(pck+j+i*(i-1)/2)
         artmp(j,i)=artmp(i,j)
      end do
!add regularization matrix to packed normal
      ngg(pck+i+i*(i-1)/2)=ngg(pck+i+i*(i-1)/2)+rdiag(side+i)
      
   end do


   !solve per block
   
   !!The Weight matrix can now be computed by solving ngg*W=ar for W
   !!where ar is the original (approximate) error covariance matrix and ngg now additionally holds the signal covariance
   !!the problem is solved using the lapack routine dppsv which uses a cholesky decomposition

   call dppsv( 'U', sz, sz, ngg(pck+1:(pck+sz*(sz+1)/2)), artmp(1:sz,1:sz), sz, info )
  
 !   write(*,*)'Packed index',pck+1,(pck+sz*(sz+1)/2)
!    write(*,*)'temporary index',1,sz
!    write(*,*)'Full BD matrix index',fll+1,fll+sz**2

   if (info.ne.0) then
      write(*,*) 'dppsv solution info =',info
      stop
   end if
   
   !copy filter matrix to final packed array (only neccesary for block diagonal approach (else the matrix will not be packed)

   if(blck)then
      arpack(fll+1:fll+sz**2)=pack(artmp(1:sz,1:sz),.true.)
   end if
      
!reset parameters
   pck=pck+sz*(sz+1)/2
   fll=fll+sz**2
   side=blockind(b)
end do

   
   !!Write weight matrix to unformatted data file

   
   dbls_d(1)='Plaw_power:'
   dbls_d(2)='Plaw_scale:'
   dbls(1)=dble(pow)
   dbls(2)=alpha

   descr='ANISOTROPIC FILTER matrix with power law regularization (alpha*l^pow)'

!construct file name (when a specific one is not requested)
   if(.not. altfname)then
      if(blck)then
         write(outfile,'(a4,i1.1,a1,i3.3,a5,i2.2,a2,i1.1)')'Wbd_',lmin,'-',lmax,'.a_1d',int(log10(alpha)),'p_',pow
      else
         write(outfile,'(a2,i1.1,a1,i3.3,a5,i2.2,a2,i1.1)')'W_',lmin,'-',lmax,'.a_1d',int(log10(alpha)),'p_',pow
      end if
   end if

   write(*,*)'Writing weight matrix to:',outfile

   if(blck)then !block diagonal packed matrix
      call write_BIN(file=trim(outfile),type=covtypout,descr=descr,blockind=blockind,&
           ints_d=ints_d,ints=ints,dbls_d=dbls_d,dbls=dbls,side1_d=side_d,pmat1=arpack)
   else!full matrix
      call write_BIN(file=trim(outfile),type=covtypout,descr=descr,&
           ints_d=ints_d,ints=ints,dbls_d=dbls_d,dbls=dbls,side1_d=side_d,mat1=artmp)
   end if


    
 end program cov_regul05


      subroutine help()
        implicit none
        write(*,*)'Program to construct anisotropic weight matrix for from GRACE (approximate) error'
        write(*,*)'Covariance and signal covariance obeying a power law:SCALE*l^POW'
        write(*,*)'Outputs a weight matrix W=(N+SCALE*R)^-1 N'
        write(*,*)'Where N is the error covariance matrix and R is the signal covariance matrix (diagona im this case)'
        write(*,*)'Usage: Cov_regul05 [OPTIONS] COVFILE'
        write(*,*)'Where the COVFILE is the error covariance file of GRACE'
        write(*,*)'Options:'
       ! write(*,*)'    -lLMAX,LMIN: specify maximum and minimum degree to consider (obligatory)'
        write(*,*)'    -pSCALE,POW: apply the power law as above with specific parameters (SCALE and POW)(obligatory)'
        write(*,*)'    -fOUTPUTFILE: specify other OUTPUTFILE name'
        stop
      end subroutine help

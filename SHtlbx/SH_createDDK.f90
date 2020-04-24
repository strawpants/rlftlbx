!!-------------------------------------------------------------------------------
!!
!!    SH_createDDK (adapted over the years, orginating from J.Kusche's cov_regul04 program)
!!    Roelof Rietbroek
!!    the program requires full or block diagonal normal matrix N in BINV format (see RLFtlbx), adds
!!    alpha times a regularization matrix R    and inverts and multiplies with N(Tikhonov-filter)
!!
!!    W = (N + aR)^-1 N
!!
!!    the R entries are defined as l**pow, where l is determined from
!!    mi,ma and assuming the order-degree ordering scheme
!!-------------------------------------------------------------------------------
program SH_createDDK
  use binfiletools
  use bin_operations
  implicit none
  
  type(BINdat)::nggdat,filtout
  
  
  
  integer::narg,lmax,lmin,pow,ind,info,nblocks
  integer::i
  integer,dimension(:),pointer::blockind=>null()

 ! integer, allocatable,dimension(:,:)::indx
  double precision::alpha,deg,dbls(2)
  double precision,pointer,dimension(:,:)::artmp=>null()
  double precision, allocatable,dimension(:)::rdiag
  character(200)::dum
  logical::pwr,blck,altfname

  integer::side,fll,pck,b,sz,j
  pwr=.false.
  altfname=.false.

  !default initialisations
  lmax=0
  lmin=0

  alpha=0.d0
  pow=0
  
  blck=.false.

  !process command line arguments
  !amount of arguments
  narg=iargc()
  if(narg < 2)call help()
  !get the (approximate error covariance file name
  call getarg(narg,nggdat%file)
  
  !now do a loop over the other arguments
  do,i=1,narg-1
     call getarg(i,dum)
     !call help when some options are no options
     if(dum(1:1).ne. '-')call help()

     select case(dum(2:2))
        case('f')
           filtout%file=trim(dum(3:))
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


  !read in the normal matrix from GRACE in packed form
  call read_BINtype(nggdat)
   do,i=1,nggdat%nint
      if(trim(nggdat%ints_d(i)) .eq. 'Lmax')lmax=nggdat%ints(i)
      if(trim(nggdat%ints_d(i)) .eq. 'Lmin')lmin=nggdat%ints(i)
      
   end do
  



  select case(nggdat%type)
  case ('BDSYMV0_','BDSYMVN_') !block diagonal matrix only
     blck=.true.
     nblocks=nggdat%nblocks
     blockind=>nggdat%blockind
     !set output type
     filtout%type='BDFULLV0'
     filtout%nblocks=nblocks 
     filtout%blockind=>nggdat%blockind
     write(*,*)'Found Block diagonal system'
     

  case('SYMV0___','SYMV1___','SYMV2___','SYMVN___')
     blck=.false.  
     
     !make the matrix one big diagonal block
     allocate(blockind(1))
     blockind(1)=nggdat%nval1
     nblocks=1

     filtout%type='FULLSQV0'
     filtout%mtyp="F"     
     
     write(*,*)'Found symmetric normal matrix'

     case default
       write(0,*)'ERROR expecting a matrix of type BDSYMV or SYMV getting:',nggdat%type        
       stop


   end select

  !construct diagonal signal covariance matrix (diagional only)
  allocate(rdiag(nggdat%nval1))
  rdiag=0.d0
  do,i=1,nggdat%nval1
     read(nggdat%side1_d(i)(5:7),*)deg !read degree from description string
     rdiag(i)=alpha*deg**pow

  end do

if(blck)then
   !get size of new block diagonal matrix which is now not symmetric anymore 
   side=0
   fll=0
   do,i=1,nggdat%nblocks
      sz=nggdat%blockind(i)-side !side size
      fll=fll+sz**2       !progression within the full matrix
      side=nggdat%blockind(i)    !reset
   end do
!   write(*,*)'packed matrix size',fll
   !allocate packedfilter matrix
   allocate(filtout%pack1(fll))
   filtout%pval1=fll
   filtout%pval2=1

   write(*,*)'filter packed matrix size',fll
   !also allocate a temporary square array (maximum size necessary is lmax-lmin+1)
   allocate(artmp(lmax-lmin+1,lmax-lmin+1))
   
else !no block diagonal approach
! allocate a temporary square array 
   allocate(artmp(nggdat%nval1,nggdat%nval1))
   filtout%pval1=nggdat%nval1**2
   filtout%pval2=1
   filtout%mat1=>artmp
end if



!!loop over diagonal blocks (only one if non block diagonal)
write(*,*)'calculating weight matrix'


side=0
fll=0
pck=0
do,b=1,nblocks
   if(blck)write(*,*)'Doing block: ',b,'of ',nblocks
   sz=blockind(b)-side
  
   !copy data from normal system to square array
   do,i=1,sz
      do,j=1,i
         artmp(i,j)=nggdat%pack1(pck+j+i*(i-1)/2)
         !also mirror the data
         artmp(j,i)=artmp(i,j)
      end do
!add regularization matrix to packed normal
      nggdat%pack1(pck+i+i*(i-1)/2)=nggdat%pack1(pck+i+i*(i-1)/2)+rdiag(side+i)
      
   end do


   !solve per block
   
   !!The Weight matrix can now be computed by solving N*W=ar for W
   !!where ar is the original (approximate) error covariance matrix and ngg now additionally holds the signal covariance
   !!the problem is solved using the lapack routine dppsv which uses a cholesky decomposition

   call dppsv( 'U', sz, sz, nggdat%pack1(pck+1:(pck+sz*(sz+1)/2)), artmp(1:sz,1:sz), sz, info )
  
 !   write(*,*)'Packed index',pck+1,(pck+sz*(sz+1)/2)
!    write(*,*)'temporary index',1,sz
!    write(*,*)'Full BD matrix index',fll+1,fll+sz**2

   if (info.ne.0) then
      write(*,*) 'Solving error, dppsv solution info =',info
      stop
   end if
   
   !copy filter matrix to final packed array (only neccesary for block diagonal approach (else the matrix will not be packed)

   if(blck)then
      filtout%pack1(fll+1:fll+sz**2)=pack(artmp(1:sz,1:sz),.true.)
   end if
      
!reset parameters
   pck=pck+sz*(sz+1)/2
   fll=fll+sz**2
   side=blockind(b)
end do

   
   !!!Write weight matrix to unformatted data file

   
   filtout%descr='ANISOTROPIC FILTER matrix with power law regularization (alpha*l^pow)'
   call BIN_alloc_dbls(filtout,2) 
   filtout%dbls_d(1)='Plaw_power:'
   filtout%dbls_d(2)='Plaw_scale:'
   filtout%dbls(1)=dble(pow)
   filtout%dbls(2)=alpha
   filtout%nint=nggdat%nint
   filtout%ints=>nggdat%ints
   filtout%ints_d=>nggdat%ints_d
   filtout%side1_d=>nggdat%side1_d
   filtout%side2_d=>nggdat%side1_d
   filtout%nval1=nggdat%nval1
   filtout%nval2=nggdat%nval1
   
   if(.not. altfname)then
      if(blck)then
         write(filtout%file,'(a4,i1.1,a1,i3.3,a5,i2.2,a2,i1.1)')'Wbd_',lmin,'-',lmax,'.a_1d',int(log10(alpha)),'p_',pow
      else
         write(filtout%file,'(a2,i1.1,a1,i3.3,a5,i2.2,a2,i1.1)')'W_',lmin,'-',lmax,'.a_1d',int(log10(alpha)),'p_',pow
      end if
   end if

  call write_BINtype(filtout) 

 end program SH_createDDK


      subroutine help()
        implicit none
        write(*,*)'Program to construct anisotropic weight matrix for from GRACE (approximate) error'
        write(*,*)'Normal matrix and signal covariance obeying a power law:SCALE*l^POW'
        write(*,*)'Outputs a weight matrix W=(N+SCALE*R)^-1 N'
        write(*,*)'Where N is the normal  matrix and R is the signal covariance matrix (diagona im this case)'
        write(*,*)'Usage: SH_createSSK [OPTIONS] NORMFILE'
        write(*,*)'Where the NORMFILE holds the normal matrix of GRACE'
        write(*,*)'Options:'
       ! write(*,*)'    -lLMAX,LMIN: specify maximum and minimum degree to consider (obligatory)'
        write(*,*)'    -pSCALE,POW: apply the power law as above with specific parameters (SCALE and POW)(obligatory)'
        write(*,*)'    -fOUTPUTFILE: specify other OUTPUTFILE name'
        stop
      end subroutine help

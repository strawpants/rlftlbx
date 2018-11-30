!!program which enables modifications of normal systems
!! possible modifications (in program order)
!! 1) transform normal systems with provided matrix
!! 2) set apriori values of normal system(s) to those values in the provided solution file
!! 3) fix parameters 
!! 4) combine systems
!! 5) reduce normal systems for parameters
!! 6) add regularization matrix
!! 7) solve remaining systems
!! 8) Do 1 iteration of the Foerstner VCE analysis

!!Coded by Roelof Rietbroek, Wed Jan 30 14:43:50 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Wed Mar 26 16:35:09 2008
!!diagonal regularization may also be done with an ascii file ( parameter name .. diagonal of regularization matrix0
!!Updated by Roelof Rietbroek, Fri May  2 13:46:24 2008
!!regularization is now only performed on parameters which are kept (action=3)

!!Updated by Roelof Rietbroek, Tue Dec  2 14:17:40 2008
!!ADDED:
!!1) Transformation is now possible (with diagonal, full and symmetric transformation matrix)
!!2) New binary files are read and outputted (using binfiletools.f90)
!!3) Interface with BLAS is improved (no unneccessary copying of arrays by compiler), BLAS operations are performed on the subsections of the original memory
!!4) output may be written  and read to/from standard output/input
!!5) improved learning curve file. Derived parameters are separated from redefinable ones.

!!Updated by Roelof Rietbroek, Tue Feb 10 22:56:29 2009
!!regularization matrix may also be block diagonal

!!Updated by Roelof Rietbroek, Mon Mar  2 16:10:19 2009
!! USE POSIX extended regular expressions in include and exclude search criteria ( much more flexible search possibilities)  ( updated in version 2.1)

!!Updated by Roelof Rietbroek, Tue Mar  3 15:48:45 2009
!! added 'forget' apriori values ( sets the apriori values to zero without changing ltpl or right hand side)

!!Updated by Roelof Rietbroek, Tue Jun  9 15:48:45 2009
!!Changed input of scalars to norm_fix_red routines to arrays for compatibility ( also changed norm_fix_red routines)

!!Updated by Roelof Rietbroek, Tue Jun 23 10:31:14 2009
!!added forgotten stop statement in error

!!Updated by Roelof Rietbroek, Thu Aug  6 16:20:13 2009
!!removed bug: forgotten close statements for learning curve files 
!! removed bug in regularization conf (do loop was changed in a do while loop but dependnecy on index i remained)

!!Updated by Roelof Rietbroek, Fri Aug 14 12:17:46 2009
!!fixed bug with fixing routine ( removed routine all together)

!!Updated by Roelof Rietbroek, Wed Aug 19 15:27:47 2009
!! fixed bug in reading diagonal ascii regularization matrix

!!Updated by Roelof Rietbroek, Wed Dec  9 17:06:50 2009
!!better handling of reduced/fixed parameters from previous calls ( affects sigma0 calculations)

!!Updated by Roelof Rietbroek, Tue Mar  2 11:22:53 2010
!! added warning when there are no matching parameters ( transformation/regularization)

!!Updated by Roelof Rietbroek, Wed Mar  3 10:39:50 2010
!! extended maximum length of regular expressions

!!Updated by Roelof Rietbroek, Mon Mar 15 11:44:10 2010
!! add a simple apriori replace option (don't change ltpl or the right hand side but do replace the apriori vector)

!!Updated by Roelof Rietbroek, Wed Apr 14 10:23:11 2010
!!Enable adding of parameters to the existing set ( -TA option)

!!Updated by Roelof Rietbroek, Fri Apr 23 12:50:43 2010
!!tried bugfix (need to be checked still)

!!Updated by Roelof Rietbroek, Thu Jul 22 18:23:25 2010
! Allowed reading multiple regular expressions from a file

!!Updated by Roelof Rietbroek, Thu Aug  5 11:22:30 2010
!! allowed giant numbers of observations by changing to integer*8

!!Updated by Roelof Rietbroek, Wed Aug 11 17:49:13 2010
!!start using more efficient sort routines to match up parameters

!!Updated by Roelof Rietbroek, Wed Dec 22 15:17:23 2010
!!In VCE contribution calculation, divide by trace for numerical stability

!!Updated by Roelof Rietbroek, Fri Mar 18 09:33:36 2011
!!make sure extra meta data integers,doubles are copied in the output file

!!Updated by Roelof Rietbroek, Wed Jun 29 10:25:09 2011
!! Added check  for square matrix with -TC option

!!Updated by Roelof Rietbroek, Mon Nov 21 10:09:40 2011
!! added warning when system has too little vectors

!!Updated by Roelof Rietbroek, Thu Nov 24 15:59:53 2011
!!Fixed bug which had eof uninitialized and another one for when maxpara !- indpara

!!Updated by Roelof Rietbroek, Wed Feb 22 17:12:47 2012
!! format the names of the output files such that it can handle up to 1000

!!Updated by Roelof Rietbroek, Thu Feb 23 09:38:11 2012
!! reset the parameter Nreduced when a covariance  transformation is requested
!! adjust and fix VCE estimation such that the original data is returned

!!Updated by Roelof Rietbroek, Wed Feb 20 14:30:52 2013
!! fixed small VCE bug

!!Updated by Roelof Rietbroek, Thu Mar  7 22:54:29 2013
!! added special file _ZERO_ as the apriori values. This will remove the
!! existing a priori information from the system

!!Updated by Roelof Rietbroek, Fri May 31 13:40:11 2013
!! added '+' option to the apriori settings such that thedata in the provided file is added to the existing
!! a priori vector

!!Updated by Roelof Rietbroek, Tue Aug 27 22:38:49 2013
!! added option to solve the system using a truncated SVD approach

!!Updated by Roelof Rietbroek, Wed Nov 20 13:44:35 2013
!! added an option to average the apriori vector when combining multiple systems


!!TODO:
!Use lapack 2.0 rectangular packed storage features?
!use more efficient sorting algorithms from get_permvecc.o and sort_various.o
!Add multiple right hand sides

program NORM_tool
use forttlbx
use binfiletools
use bin_operations
implicit none
integer::nf,nfcmd,narg,itharg,i,j,k,nmax,stderr,eof,ind
parameter(nmax=102,stderr=0)
character(200)::basename,solfile,regfile,conffile,contfile
character(2000)::dum
character(80)::descr
character(200)::normfiles(nmax)
! regular expression indices
integer::apri_inc,red_inc, fix_inc,apri_exc,red_exc,fix_exc,reg_exc,reg_inc,tr_inc,tr_exc
logical::apriori,fix,redu,combine,forcecomb,weight,verbose,warning
logical::autofix,autoredu,regul,solve,learn,learnexist,vce,same,extconfig,vcereloop,diagregul
logical::haltafterlearn,unique,apri_ascii
character(20)::version
integer,dimension(nmax)::comb,ndat,ndat2,nfix,nfixold,nred,validsys
integer::indpara,col,row,clorig,rworig,napri,nreg,shift,nfout
double precision,pointer,dimension(:,:)::upper=>null()
double precision,pointer,dimension(:,:)::lower=>null()
double precision,pointer,dimension(:,:)::P1=>null()
double precision,pointer,dimension(:,:)::P2=>null()
double precision,pointer,dimension(:,:)::P3=>null()
double precision,dimension(nmax)::W
double precision::regscale,trace
character(24),allocatable::unknowns(:)
! character(24),target::outdbls_d(5),outints_d(5)
! double precision,target::outdbls(5)
! integer*8,target::outints(5)
integer::st,nd,rowind,colind,nweights,ns,rnval,snval,unit,maxpara,nremov
integer::msize,nsys,min1,max1,min2,max2
integer*8::ndbig
integer::stind,ndind !start and end parameter of search string
character(4)::chardum
!!sorting arrays
integer,allocatable::sortindex(:,:),aprisort(:),regsort(:),action(:),sorter(:),sortafterB(:,:)

!!temporary arrays
double precision,allocatable::tmpvec1(:),tmpvec2(:)

!define a new derived type for normal systems
 !type normalsystem
type normsys
!!   double precision,pointer::norm(:,:) !removed to prevent dangling pointers
   double precision,pointer::Atb(:)=>null() 
   double precision,pointer::Atbold(:)=>null() 
   double precision,pointer::Apri(:)=>null() 
   double precision,pointer::Apriold(:)=>null() 
   double precision,pointer::contr(:)=>null() 
   double precision::ltpl,sig0,vce,Stime,Ctime,Etime,ltplold
   integer*8::nobs,nobsold,nred,nfix

end type normsys

type(normsys),target::normass(nmax) !assembly of normal systems (maximum of nmax) must be a target as upper, and lower will be pointing to sections of it

!!new derived type variable in version 1.1
type(BINdat)::BINinput(nmax),APRIinput,REGinput,TRANSinput!! input systems Apriori fixing and regularization matrix
logical::diagtrans,covscale,trans,Btrans,symtrans,warned,N_ignore,atb_ignore,x0_ignore,x0_as_b
character(200)::tr_file
double precision,dimension(nmax)::cv_sc
double precision,allocatable::B(:,:)
character(24),pointer::transide1(:)=>null()
character(24),pointer::transide2(:)=>null()
integer,allocatable::transort1(:),transort2(:)
integer::tnval1,tnval2,info
integer::ntimings ! parameters to do timing (maximum of 20 timings)
real::cpu(20),cpu_tot
character(30)::timings(20)
integer::nfmax,maxndat
integer::nb,nk,nunit
logical::covtrans
logical::aprireplace,apriforget,apriadd
integer::iargc
logical::transadd
double precision::trc1,trc2
integer::mxrw 
logical::truncsvd
integer::ntrunc
logical::CombapriAverage
double precision::totalAveweight,Aveweight
!namelists one with derived parameters (should not be changed) and one with changable parameters (may be changed)
!the following may be redefined
namelist/pararedef/weight,descr,basename,regfile,solfile,normfiles,W,cv_sc,covscale,tr_file,&
     regscale,atb_ignore,x0_ignore,N_ignore,covtrans,apriforget,aprireplace,apriadd,truncsvd,&
     ntrunc,CombapriAverage
!parameters below are derived parameters ( don't change)
namelist/paraconf/indpara,nf,rnval,snval,ndat,ndat2,nfix,nred,nremov,apriori,apri_ascii,fix,forcecomb,combine,&
     redu,regul,diagregul,solve,vce,comb,diagtrans,trans,symtrans,transadd,tnval1,tnval2,Btrans,nunit



!defaults initializations

!start of computing
call cpu_time(cpu(1)) ! get start time
ntimings=1

normfiles=''

descr='OUTPUT From NORM_tool'

version='EXPERIMENTALV2.2'
stind=1
ndind=24

comb=0

normass(:)%ltpl=0.d0 
normass(:)%sig0=0.d0
normass(:)%nred=0
normass(:)%nfix=0
normass(:)%vce=0.d0
normass(:)%Stime=0.d0 
normass(:)%Ctime=0.d0 
normass(:)%Etime=0.d0   


covscale=.false. !scalar scaling
cv_sc=1.d0
covtrans=.false.

!!Transformation parameters
diagtrans=.false. !diagonal scaling
trans=.false.
Btrans=.false.
symtrans=.false.
transadd=.false.
tr_file=''
tr_exc=0
tr_inc=0
tnval1=0
tnval2=0
atb_ignore=.false.
x0_ignore=.false.
x0_as_b=.false.
N_ignore=.false.
nunit=0

!apriori parameters
apriforget=.false.
aprireplace=.false.
apriadd=.false.

apriori=.false.
apri_ascii=.false.
apri_exc=0
apri_inc=0
solfile=''
snval=0 
!fixing parameters
fix=.false.
autofix=.false.
fix_exc=0
fix_inc=0
nfix=0

!combination parameters
combine=.false.
forcecomb=.false.
CombapriAverage=.false.
!reduction parameters
redu=.false.
autoredu=.false.
red_exc=0
red_inc=0
nred=0

!regularization parameters
regul=.false.
diagregul=.false.
regfile=''
regscale=1.d0
reg_exc=0
reg_inc=0
rnval=0

!solving parameters
solve=.false.
truncsvd=.false.
ntrunc=0

!vce options 
vce=.false.
vcereloop=.false.

!other parameters
basename='NORMtoolout'
W=1.d0

weight=.false.

verbose=.false.
warning=.true.

learn=.false.
extconfig=.false.

haltafterlearn=.false.

ndat=0
ndat2=0

nf=0
nfcmd=0

! !output meta data description
! !integers
! outints_d(1)='Nobs'
! outints_d(2)='Nunknowns'
! outints_d(3)='Nreduced'
! outints_d(4)='Nfixed'
! outints_d(5)='Modnr'
! !doubles
! outdbls_d(1)='CTime'
! outdbls_d(2)='STime'
! outdbls_d(3)='ETime'
! outdbls_d(4)='LtPL'
! outdbls_d(5)='Sigma0'



!!end initializations


!get command line arguments

narg=iargc() !number of command line arguments
itharg=0
if(narg < 1)call help(trim(version))

do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)

   if(dum(1:1).eq.'-')then !option
      select case(dum(2:2))
      case('T')!transformation parameters
         select case(dum(3:3))
         case('f')
            trans=.true.
            if(dum(5:5) .eq. ' ')then
               itharg=itharg+1
               call getarg(itharg,tr_file)
            else
               tr_file=trim(dum(5:200))
            end if
         case('C')! covariance transformation
            covtrans=.true.
         case('d')
            if(extconfig)cycle 

            diagtrans=.true.
            symtrans=.true.
         case('t')
            if(extconfig)cycle

            Btrans=.true.
         case('i')!specify include terms
            if(extconfig)cycle
            
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  tr_inc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  tr_inc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select

         case('e')!specify exclude terms
            if(extconfig)cycle
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  tr_exc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  tr_exc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
         case('b')
            atb_ignore=.true.
         case('x')
            if(dum(3:3) .eq. 'b')then
               x0_as_b=.true.
            else
               x0_ignore=.true.
            end if
         case('n')! ignore matrix propagation
            N_ignore=.true.
         case('A')
            transadd=.true.
            x0_ignore=.true. ! implies x0_ignore automatically
         end select

      case('a')!apriori options
        
         select case(dum(3:3))
         case('i')!specify include terms
            if(extconfig)cycle
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  apri_inc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  apri_inc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
         case('e')!specify exclude terms
            if(extconfig)cycle
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  apri_exc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  apri_exc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
         case('+')! add the a priori to the vetor already exisisting
            apriadd=.true.
         case('0') ! forget apriori values of the normal system
            apriforget=.true.
         case('f') !specify Solution file name
            if(dum(5:5) .eq. ' ')then
               itharg=itharg+1
               call getarg(itharg,solfile)
            else
               solfile=trim(dum(5:200))
            end if
            apriori=.true. !only when a file name has been provided
         case('r') ! simply replace the apriori values
            aprireplace=.true.
         case('t')! solution file is in plain ascii (textfile)
            if(extconfig)cycle
            apri_ascii=.true.
         case default
            write(stderr,*)'ERROR: apriori option is not complete (try -af/-ae/-ai):'
            stop
         end select


      case('f')!fixing options

         select case(dum(3:3))
         case('i')!specify include terms
            if(extconfig)cycle
!            write(0,*)trim(dum(5:))
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  fix_inc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  fix_inc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
            fix=.true.
         case('e')!specify exclude terms
            if(extconfig)cycle
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  fix_exc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  fix_exc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select

            fix=.true.
         case('a')
            if(extconfig)cycle
            autofix=.true.
            fix=.true.
         case default
            write(stderr,*)'ERROR: fix option is not complete (try -fa/-fe/-fi):'
            stop
         end select

      case('c')!Combination options
         if(extconfig)cycle
         combine=.true.
         select case(dum(3:3))
         case('f')
            forcecomb=.true.
         case('a')
            CombapriAverage=.true.
            forcecomb=.true.
         end select
      case('r')!reducing option
         if(extconfig)cycle

         select case(dum(3:3))
         case('i')!specify include terms
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  red_inc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  red_inc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select

            redu=.true.

         case('e')!specify exclude terms
            ! compile regular expression and get index
                        select case(dum(4:4))
               case(':')
                  red_exc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  red_exc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
            redu=.true.
         case('a')!try an auto reduce
            autoredu=.true.
            redu=.true.
         case default
            write(stderr,*)'ERROR: reduce option is not complete (try -ra/-re/-ri):'
            stop
         end select
       
      case('R')!Regularization options
        
         select case(dum(3:3))
            case('f')
               if(dum(5:5) .eq. ' ')then
                  itharg=itharg+1
                  call getarg(itharg,regfile)
               else
                  regfile=trim(dum(5:200))
               end if

               regul=.true. !can only be regularized when a file has been supplied
            case('d') !diagonal regularization
               if(extconfig)cycle
               regul=.true.
               diagregul=.true.
            case('s')
               read(dum(5:),*)regscale
            case('i')!specify include terms
            if(extconfig)cycle
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  reg_inc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  reg_inc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
               
            case('e')!specify exclude terms
            if(extconfig)cycle
            ! compile regular expression and get index
            select case(dum(4:4))
               case(':')
                  reg_exc=regcompf(trim(dum(5:)),.true.)
               case('=')
                  reg_exc=regcompf(trim(dum(5:)))
               case default 
                  write(stderr,*)"ERROR processing option:",trim(dum)
            end select
         end select
      case('s')!solving options
         solve=.true.
         if(dum(3:4).eq. 't=')then
            truncsvd=.true.
            read(dum(5:),*)ntrunc
         end if
      case('V')!VCE options
            vce=.true.
            if(dum(3:3) .eq. 'p')vcereloop=.true.


!!!!!!!!!!!other parameters
      case('C')! covariance scale of systems
         select case(dum(3:3))
         case('=')
            covscale=.true.
            st=4
            do,j=1,nmax
               nd=index(dum(st:),'/')+st-2
               if(nd .eq. st-2)nd=len(dum)!set to end of string if no '/' occurs (only one weight)
               
               read(dum(st:nd),*)cv_sc(j) !read in covraiance scale weight
               !             write(*,*)st,nd,dum(st:nd),covscale(j)
               if(nd .eq. (len(dum)))exit !last scale has been read
               st=nd+2
            end do
            ns=j !set number of scales provided
         end select
      case('W')!weight of systems
         weight=.true.
            st=4
            
            do,j=1,nmax
               nd=index(dum(st:),'/')+st-2
             if(nd .eq. st-2)nd=len(dum)!set to end of string if no '/' occurs (only one weight)
             
             read(dum(st:nd),*)W(j) !read in weight
             if(nd .eq. (len(dum)))exit !last weight has been read
             st=nd+2
            end do
            nweights=j !set number of weights provided
      case('S')!limit global search area of the parameter strings
         if(extconfig)cycle
         select case(dum(3:3))
         case('=')
            st=4
            nd=index(dum(st:),'/')+st-2
            read(dum(st:nd),*)stind
            read(dum(nd+2:),*)ndind
!          case('1')
!             st=5
!             nd=index(dum(st:),'/')+st-2
!             read(dum(st:nd),*)stind1
!             read(dum(nd+2:),*)ndind1
         end select
      case('D')!use a different description for the output files
         descr=dum(4:84)

      case('L') !use or make a learning curve
         
         select case(dum(3:3))
         case('w')
            learn=.true.
            if(dum(5:5) .eq. ' ')then
               itharg=itharg+1
               call getarg(itharg,conffile)
            else
               conffile=trim(dum(5:200))
            end if
         case('h')
            haltafterlearn=.true.
         case('f')
            extconfig=.true.
            unit=freeunit()
            
            if(dum(5:5) .eq. ' ')then
               itharg=itharg+1
               call getarg(itharg,dum)
               open(unit=unit,file=trim(dum),form='formatted',delim='QUOTE')
            else
               open(unit=unit,file=trim(dum(5:)),form='formatted',delim='QUOTE')
            end if
            
            
            read(unit,*)dum

            read(unit,*)dum

            !read namelists
            read(unit,nml=pararedef) ! write list with redefinable parameters
            read(unit,*)dum
            read(unit,nml=paraconf)
            read(unit,*)dum
            allocate(action(indpara),unknowns(indpara))
            
            if(forcecomb)then
               allocate(sortindex(indpara,nf+1))
               nd=nf+1
            else
               allocate(sortindex(indpara,nf))
               nd=nf
            end if
            if(nunit>0)allocate(sorter(indpara))
            if(apriori)allocate(aprisort(indpara))
            if(regul)allocate(regsort(indpara))
            if(trans)allocate(transort1(indpara),transort2(indpara))
            !read sorting indices
            do,j=1,indpara
               !read unknowns and action (always present)
               read(unit=unit,fmt='(A24,1x,i2)',ADVANCE='NO')unknowns(j),action(j)
               if(nunit>0)read(unit=unit,fmt='(1x,i7)',ADVANCE='NO')sorter(j)
               if(apriori)read(unit=unit,fmt='(1x,i7)',ADVANCE='NO')aprisort(j)
               if(regul)read(unit=unit,fmt='(1x,i7)',ADVANCE='NO')regsort(j)
               if(trans)read(unit=unit,fmt='(1x,i7)',ADVANCE='NO')transort1(j)
               if(trans .and. .not. symtrans)read(unit=unit,fmt='(1x,i7)',ADVANCE='NO')transort2(j)
               !read last sorting indices
               read(unit,*)sortindex(j,1:nd)
            end do
            close(unit)
            if(symtrans)transort2=transort1 ! copy index (sides are the same in the symmetric case)
         end select

      case('v')! verbose
         verbose=.true.
      case('w')!supress warning messages
         warning=.false.
      case('F')! Basename option
         basename=trim(dum(4:200))

      case('h')
         call help(trim(version))
      case default
         write(stderr,*)'ERROR: Unknown option selected:', trim(dum)
         write(stderr,'(A)')" Type -h for help"
         stop
      end select

   else !normal equation file
      
      nfcmd=nfcmd+1
      if(nfcmd > nmax-1)then
         write(stderr,*)'ERROR: Maximum amount of normal systems allowed is ',nmax-1
         stop
      end if
      normfiles(nfcmd)=trim(dum(1:200))
   end if
end do








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!input checks!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(covtrans .and. vce)then
   write(stderr,*)" ERROR: Covariance transformation may not be combined with vce estimation in 1 run"
   stop
end if

if( covtrans .and. (N_ignore .or. x0_as_b .or. x0_ignore .or. atb_ignore))then
   write(stderr,*)" ERROR: Covariance transformation may not be combined with -Tx,-Tb -Txb -Tn"
   stop
end if

if(x0_as_b .and. x0_ignore)then
   write(stderr,*)" ERROR: Apriori vector must be either ignored or treated as a right hand vector, not both"
   stop
end if

if(extconfig .and. (nf .ne. nfcmd) .and. (nfcmd .ne. 0))then
   write(stderr,*)" ERROR: amount of normal files provided does not agree with that from the configuration file"
   stop
else if( .not. extconfig)then
   nf=nfcmd
end if

if(trim(solfile) .eq. '_ZERO_' .and. apriforget )then
   write(stderr,*)" ERROR: using -af='_ZERO_' in combination with -a0 makes no sense"
   stop
end if

if(extconfig)then
   ns=nf
   nweights=nf
end if



!check whether amount of weights provided agrees with the amount of input systems

if(weight .and. nweights .ne. nf )then
   write(stderr,*)"ERROR, amount of weights provided do not match the amount of input systems"
   stop
end if

if(covscale .and. ns .ne. nf )then
   write(stderr,*)"ERROR, amount of covariance scales provided,",ns,", do not match the amount of input systems:",nf
   stop
end if


!check whether -f plain option is provided with the -a option
if(autofix .and. .not. apriori)then
   write(stderr,'(A)')"ERROR: For the plain autofix option -fa an apriori solution file must be provided"
   write(stderr,'(A)')" Type -h for help"
   stop
end if

!check whether autoreduce and force combine is requested at the same time
if(autoredu .and. forcecomb)then
   write(stderr,'(A)')"ERROR: The auto reduce (-ra) and force combine (-cf) option cannot"
   write(stderr,'(A)')"       be requested at the same time"
   stop
end if

if(nf .eq. 0)then
   write(stderr,'(A)')'ERROR: no input files'
   stop
end if

if(vce .and. .not. (combine .and. solve))then
   write(stderr,*)"ERROR: The VCE option can only be used with the combine and solve option"
   stop
end if

if(vce .and. (redu .or. fix .or. apriori .or. trans))then
   write(stderr,*)"ERROR: the VCE option may not be combined with reduction, fixing, or changing apriori values"
   write(stderr,*)"You may want to apply this in a preprocessing step to achieve want you want" 
   stop
end if


if(transadd .and. (covtrans .or. diagtrans .or. symtrans))then
   write(stderr,*)"ERROR: The -TA option may not be used with -TC, -Td or a symmetric matrix"
   stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end input checks!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!some default settings dependent on the input !!





if( .not. extconfig)then ! copy some info in general variables and do sorting (this IF covers all sorting)
!!!!!!!!!!!!!!!!!READ META DATA of input FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(verbose)write(*,'(A)')'VERBOSE: reading parameter descriptions of input systems ..'   
   do,i=1,nf !first two segments per file

      BINinput(i)%file=trim(normfiles(i))
      call read_BINtype(BINinput(i),2)
      !input check
      select case(BINinput(i)%type)
      case('SYMVN___','SYMV1___','SYMV2___') ! symmetric matrices
      case('SYMV0___') ! symmetric matrix but no right hand side
         write(stderr,*)'ERROR: Input system',trim(BINinput(i)%file),'has no right hand side'
         stop 
      case('FULLSQV0','FULLSQVN','FULL2DVN') ! full matrices are not allowed
         write(stderr,*)'ERROR: Input system',trim(BINinput(i)%file),'is not symmetric'
         stop 
      case default
         write(stderr,*)'ERROR: Input system',trim(BINinput(i)%file),'is not supported'
         stop 
      end select

      ndat(i)=BINinput(i)%nval1
   end do
   
   if(trans)then ! read META data of Transformation matrix
      if(verbose)write(*,'(A)')'VERBOSE: reading parameter description of transformation system ..'   
      if(diagtrans)then !diagonal ascii file
         unit=freeunit()
         tnval1=0
         eof=0
         open(unit=unit,file=trim(tr_file),form='formatted',status='old')
         do while (eof .eq. 0)! loop until end of file
            read(unit=unit,fmt='(A1)',iostat=eof)dum(1:1)
            if(eof .ne. 0)exit
            tnval1=tnval1+1

         end do
         rewind(unit)

         tnval2=tnval1
         TRANSinput%nval1=tnval1
         allocate(TRANSinput%side1_d(tnval1))
         do,i=1,tnval1
            read(unit=unit,fmt='(A24)')TRANSinput%side1_d(i)
            !write(*,*)TRANSinput%side1_d(i)
         end do
         
         close(unit)
         !!associate pointers
            transide1=>TRANSinput%side1_d
            transide2=>TRANSinput%side1_d

      else ! full matrix  (may be symmetric)
         TRANSinput%file=trim(tr_file)
         call read_BINtype(TRANSinput,2)
         
         select case(TRANSinput%type)
         case('SYMV0___','SYMVN___','SYMV1___','SYMV2___') ! symmetric matrix
            symtrans=.true.
            if(transadd)then
               write(stderr,*)"ERROR NORM_tool: symmetric matrix is not allowed with -TA"
               stop
            end if
            TRANSinput%mtyp='U'
            tnval1=TRANSinput%nval1
            tnval2=tnval1
            !associate pointers
            transide1=>TRANSinput%side1_d
            transide2=>TRANSinput%side1_d
         case('FULLSQV0','FULLSQVN','FULL2DVN') ! allowed full matrices
            
            !associate pointers
            if(Btrans)then
               tnval1=TRANSinput%nval2
               tnval2=TRANSinput%nval1
               transide1=>TRANSinput%side2_d
               transide2=>TRANSinput%side1_d
            else
               tnval1=TRANSinput%nval1
               tnval2=TRANSinput%nval2
               transide1=>TRANSinput%side1_d
               transide2=>TRANSinput%side2_d
            end if
            
         case default
            write(stderr,*)'ERROR: Transformation matrix system is unsupported:',TRANSinput%type
            stop
         end select
         
      end if

      if(warning .and. (diagtrans .or. symtrans) .and. Btrans)then ! warning
         write(stderr,*)'WARNING: Option -Tt is useless with a diagonal or symmetric transformation matrix'
         stop
      end if

      if(tnval1 .ne. tnval2 .and. covtrans)then
         write(stderr,*)'ERROR: Error covariance transformation is requested but transformation matrix is rectangular'
         stop
      end if
   
   end if

   
   

!!!!!!!!!!!!!!!!!!!!INDEPENDENT PARAMETER CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !get amount of independent parameters
   
   if(verbose)write(*,'(A)')'VERBOSE: searching for independent parameters ..'
   !calculate maximum amount of parameters possible
   maxpara=sum(ndat)+tnval1
   
   allocate(unknowns(maxpara))
   unknowns=''


   if(forcecomb)then
      allocate(sortindex(maxpara,nf+1))
   else
      allocate(sortindex(maxpara,nf))
   end if
   sortindex=0

   !find out system which has the most parameters
   maxndat=maxval(ndat)
   do,i=1,nf
      if(ndat(i).eq. maxndat)then
         nfmax=i
         exit
      end if
   end do

   !now copy variables from this system directly in the unknowns
   
   indpara=maxndat

   forall(i=1:indpara)sortindex(i,nfmax)=i 

   ! copy parameters of this file directly

   unknowns(1:indpara)=BINinput(nfmax)%side1_d(1:indpara)
   !loop over remaining parameters vectors per file
   
   do,i=1,nf
      if(i .eq. nfmax)cycle ! already done
      do,j=1,ndat(i)
         unique=.true.
         
         !check whether parameter is already included
         do,k=1,indpara
            if(BINinput(i)%side1_d(j)(stind:ndind) .eq. unknowns(k)(stind:ndind))then
               sortindex(k,i)=j
               unique=.false.
               exit
            end if
         end do
         
         if(unique)then !not found !add new independent parameter to list
            indpara=indpara+1
            unknowns(indpara)=BINinput(i)%side1_d(j)
            sortindex(indpara,i)=j 
         end if

      end do
   end do


!!!!!!!!!!!!!!!!!!!TRANS CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(trans)then ! append unknowns with row (or column if Btrans) names from transformation matrix

      if(verbose)write(*,'(A)')'VERBOSE: checking which parameters will be changed through transformation'
      allocate(transort1(maxpara))
      transort1=0

      do,j=1,tnval1
         unique=.true.
         !check whether parameter should be excluded based on parameter search strings
         if(tr_exc > 0)then
            if(regexecf(tr_exc,transide1(i)))cycle
         end if

         if(tr_inc > 0)then
            if(.not. regexecf(tr_inc,transide1(i)))cycle
         end if

         
         !check whether parameter is already included
         do,k=1,indpara
            if(transide1(j)(stind:ndind) .eq. unknowns(k)(stind:ndind))then
               transort1(k)=j
               unique=.false.
               exit
            end if
         end do
         

         if(unique)then !not found !add new independent parameter to list
            indpara=indpara+1
            unknowns(indpara)=transide1(j)
            transort1(indpara)=j 
         end if

      end do

      ! now  determine which parameters get removed through a transformation  (set action=1)
      allocate(transort2(indpara))

      if( .not. symtrans)then ! side two must also be checked
         transort2=0
         do,j=1,indpara
            do,k=1,tnval2
               if(unknowns(j)(stind:ndind) .eq. transide2(k)(stind:ndind))then
                  transort2(j)=k
                  exit
               end if
            end do
                       
         end do
      else
         transort2(1:indpara)=transort1(1:indpara)

      end if


      !! CHECK MISSING PARAMETERS!!!!!!!!!!!!!!!!!!
      if(warning)then ! make a check whether some parameters in the transformation matrix will be ignored
         if(verbose)write(*,'(A)')"VERBOSE: checking missing transformation parameters"
         warned=.false.
         do,i=1,nf
            do,j=1,indpara
               if(sortindex(j,i) .eq. 0 .and. transort2(j) > 0)then
                  write(stderr,*)"WARNING: input system",i,"is missing parameter",unknowns(j) 
                  write(stderr,*)"needed for a complete transformation"
                  warned=.true.
                  cycle
               end if
            end do
         end do
         
      !check whether there are no parameters which will be transformed and issue a warning
         if( sum(transort2) .eq. 0)then
            write(stderr,*)"WARNING: -T is requested but no parameters will be transformed"
            warned=.true.
         end if

         if(warned)then
            write(stderr,*)"WARNING: using -w will ignore the above warnings"
            stop
         else
            if(verbose)write(*,'(A)')"VERBOSE: check OK"
         end if
      end if



      
   end if

   
   
!!!!!!!!!!!!!!!!!!! END TRANS CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !allocate action vector
   allocate(action(indpara))
   action=4 !default is keep all parameters (action =4)   

   if(trans .and. .not. symtrans .and. .not. transadd)then !tag parameters which will be removed through transformation
      where((transort1(1:indpara) .eq. 0) .and. (transort2(1:indpara)>0))action(1:indpara)=1
   end if

   
   
!!!!!!!!!!!!!!!!!APRIORI CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !get apriori description vector if provided
   !and check which parameters should be set to the apriori vectors
   
   if(apriori)then
      if(verbose)write(*,'(A)')'VERBOSE: checking which apriori parameters need to be changed'
      
      if(trim(solfile).eq. '_ZERO_')then ! special case (set the  a priori values to 0)
         snval=indpara
         !allocate side descriptor array explicitly
         allocate(APRIinput%side1_d(snval))
         APRIinput%side1_d(1:snval)=unknowns(1:snval) !just copy all the entries
         
      else if(apri_ascii)then!apriori file is an ascii file
         unit=freeunit()
         open(unit=unit,file=trim(solfile),form='formatted')
         !first find out how many lines are in the file
         eof=0
         snval=0
         do while(eof .eq. 0)
            read(unit=unit,fmt=*,iostat=eof)dum(1:1)
            if(eof .eq. 0)snval=snval+1
         end do
         !allocate side descriptor array explicitly
         allocate(APRIinput%side1_d(snval))
         !read parameter names from file
         rewind(unit)
         do,i=1,snval
            read(unit,'(a24)')APRIinput%side1_d(i)
         end do
         
         close(unit)
      else ! input is a solution file
         APRIinput%file=trim(solfile)
         call read_BINtype(APRIinput,2) ! read first two segment of the solutionfile
         snval=APRIinput%nval1
      end if




      !check which apriori values may be changed
      allocate(aprisort(indpara))
      aprisort=0
      do,i=1,snval
         do,j=1,indpara
            if(aprisort(j) .ne. 0 .or. action(j) <= 1)cycle !quick cycle when parameter already found or if parameter will be transformed already
            
            if(unknowns(j)(stind:ndind) .eq. APRIinput%side1_d(i)(stind:ndind))then
               aprisort(j)=i
               exit
            end if
         end do
      end do

!additional include and exclude requirements checks
      if(apri_inc > 0)then
         do,i=1,indpara
            if(aprisort(i) .eq. 0)cycle
            if(.not. regexecf(apri_inc,unknowns(i)))aprisort(i)=0
         end do
      end if

      if(apri_exc > 0)then
         do,i=1,indpara
            if(aprisort(i) .eq. 0)cycle
            if(regexecf(apri_exc,unknowns(i)))aprisort(i)=0
         end do
      end if
      
   end if
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIX PARAMETER CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   
   if(fix)then
      if(verbose)write(*,'(A)')'VERBOSE: checking which parameters need to be fixed'
      !start of with fixing parameters which will not be transformed
      where(action > 1)action=2
      
      if(autofix)then !fix all parameters which were set in the apriori approach
         where(aprisort .eq. 0)action=4 !set back to keep action
      end if
      
      if(fix_inc > 0)then
         do,i=1,indpara
            if(action(i) .ne. 2)cycle
            if(.not. regexecf(fix_inc,unknowns(i)))action(i)=4 ! set back to keep
         end do
      end if

      if(fix_exc > 0)then
         do,i=1,indpara
            if(action(i) .ne. 2)cycle
            if(regexecf(fix_exc,unknowns(i)))action(i)=4 ! set back to keep
         end do
      end if
   end if
   
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COMBINE PARAMETER CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(forcecomb)then
      if(verbose)write(*,'(A)')'VERBOSE: Planning forced combination'
      comb(1:nf+1)=nf+1 !combine all systems to system nf+1
   else if(combine)then
      if(verbose)write(*,'(A)')'VERBOSE: Planning Auto combination'
      !start with the default comb =[1,2,3,4..] !no combination
      forall(i=1:nf)comb(i)=i
      do,i=1,nf-1
         if(comb(i) .ne. i)cycle !cycle when a system has already been found to be similar to an earlier system
         do,k=i+1,nf
            same=.true.!assume the systems are the saem on loop entry 
            !check wheter systems contain the same parameters
            do,j=1,indpara !loop untill proven otherwise (systems are not the same)
               if(action(j) <= 2)cycle !ignore fixed and transformed parameters
               !!NOTE the parameters which will be added through a possible transformation are the same for all systems and need not be checked
               if(.not.((sortindex(j,i)*sortindex(j,k) > 0) .or. (sortindex(j,k)+sortindex(j,i) .eq. 0)))then
                  same=.false.
                  exit
               end if
            end do
            if(same)comb(k)=i
         end do
      end do
   else
      forall(i=1:nf)comb(i)=i
     
   end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REDUCE PARAMETER CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   if(redu)then
      if(verbose)write(*,'(A)')'VERBOSE: checking which parameters need to be reduced'
      !start with setting non-fix parameters to reduce
      where(action > 2)action=3 !sets all keep parameters to reduce first
      
      if(autoredu)then !look for local parameters after the combination
         if(maxval(comb) .eq. 1)then !no autoreduction possible (ignore but give a warning)
            if(warning)then
               write(stderr,*)"WARNING: auto reduce option is ignored, only one system left after combination"
               write(stderr,*)"WARNING: auto reduce option will be set to false on a run with -w"
               stop
            end if
            autoredu=.false.
         else
            !determine the valid systems after the combination
            validsys=1
            do,j=1,nf
               if(comb(j) .ne. j)validsys(j)=0 !tag ambigious system with a zero
            end do
            
            do,i=1,indpara
               if(action(i) <= 2)cycle !cycle straight away when parameter will be transformed or fixed first
               
               !      count((validsys*sortindex(i,1:nf)) > 0) is the number of nonzero entries in the nonambiguous systems
               !if this is one it is a local parameter and will be reduced else it will be kept
               
               if(count((validsys(1:nf)*sortindex(i,1:nf)) > 0) > 1)action(i)=4 !set action to keep since it it not a local parameter
               
            end do
         end if
      end if
      
      !explicit include and exclude terms
      if(red_inc > 0)then
            do,i=1,indpara
               if(action(i) .ne. 3)cycle !no need to check the next condition
               
               if(.not. regexecf(red_inc,unknowns(i)))action(i)=4 !if condition is not satisfied
            end do
      end if

      if(red_exc > 0)then
            do,i=1,indpara
               if(action(i) .ne. 3)cycle !no need to check the next condition
               if(regexecf(red_exc,unknowns(i)))action(i)=4 !if condition is not satisfied
            end do
      end if
         
   end if !reduce
   
   
!!!!!!!!!!!!!REGULARIZATION PARAMETER CHECK!!!!!!!!!!!!!!!!!!!!!!!
   if(regul)then
      if(verbose)write(*,'(A)')'VERBOSE: checking which parameters need to be regularized'

      allocate(regsort(indpara))
      if( .not. diagregul)then ! full matrix
         REGinput%file=trim(regfile)
         call read_BINtype(REGinput,2)
        select case(REGinput%type)
        case('SYMVN___','SYMV0___','SYMV1___','BDSYMVN_','BDSYMV0_')
           !do nothing simply accept
        case default
           write(stderr,*)"ERROR: Regularization matrix is not symmetric or not supported",REGinput%type   
           stop
        end select

         rnval=REGinput%nval1
         !check which parameters may be regularized witgh the input matrix
         regsort=0
         do,i=1,rnval
            do,j=1,indpara
               if(regsort(j) .ne. 0 )cycle !quick cycle when parameter already found
               if(action(j) .ne. 4)cycle !quick cycle when parameter does not need to be regularized (either fixed or reduced already)
               if(unknowns(j)(stind:ndind) .eq. REGinput%side1_d(i)(stind:ndind))then
                  regsort(j)=i
                  exit
               end if
            end do
         end do
      else if(diagregul .and. regfile .ne. '')then !diagonal regularization with file

         regsort=0
         unit=freeunit()
         eof=0
         rnval=0
         open(unit=unit,file=trim(regfile),form='formatted')
   
         do while (eof .eq. 0)! loop until end of file
            read(unit=unit,fmt='(A24)',iostat=eof)dum(1:24)
            if(eof .ne. 0)exit
            rnval=rnval+1
            !loop over unknowns
            do,j=1,indpara
               if(regsort(j) .ne. 0)cycle !quick cycle when parameter already found
               if(action(j) .ne. 4)cycle !quick cycle when parameter does not need to be regularized (either transformed,fixed or reduced already)
               if(unknowns(j)(stind:ndind) .eq. dum(stind:ndind))then
                  regsort(j)=rnval
                  exit
               end if
            end do

         end do
         close(unit)
      else !simple unit diagonal regularization

         forall(i=1:indpara)regsort(i)=i !default assumes all parameters to be regularized
      end if

      !additional include and exclude requirements checks
      if(reg_inc > 0)then
         do,i=1,indpara
            if(regsort(i) .eq. 0)cycle
            if(.not. regexecf(reg_inc,unknowns(i)))regsort(i)=0 !if condition is not satisfied
         end do
      end if

      if(reg_exc > 0)then
         do,i=1,indpara
            if(regsort(i) .eq. 0)cycle
            if(regexecf(reg_exc,unknowns(i)))regsort(i)=0 !if condition is not satisfied
         end do
      end if
      
      !give warning when no parameters will be regularized
      if(warning .and. sum(regsort) .eq. 0)then
         write(stderr,*)"WARNING: -R is requested but no parameters will be regularized"
         stop
      end if

   end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SORTING/PERMUTATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Depending on whether a transformation is requested we need 1 or 2 permutations of the data
   if(verbose)write(*,'(A)')'VERBOSE: calculating necessary permutations'

   if(trans)then
      
      !calculate maximum amount of parameters which will be transformed through the unit matrix
      nunit=count((transort1(1:indpara)+transort2(1:indpara)) .eq. 0)
      if(nunit > 0)then !additional permutation necessary
         allocate(sorter(indpara))
         rowind=0
         shift=0
         forall(i=1:indpara)sorter(i)=i
         do,i=1,indpara
            if(transort1(i)+transort2(i) .eq. 0)then
               rowind=rowind+1
               sorter(rowind)=i
            else
               shift=shift+1
               sorter(nunit+shift)=i
            end if
         end do

!          do,i=1,80
!             write(*,*)i,unknowns(i),sorter(i),transort1(i),transort2(i)
!          end do
      !permute all indices, unknowns etc
         call permute(transort1(1:indpara),sorter)

         if(symtrans)then
            transort2(1:indpara)=transort1(1:indpara)
         else
            call permute(transort2(1:indpara),sorter)
         end if
  
         if(apriori)call permute(aprisort(1:indpara),sorter)
         if(regul)call permute(regsort(1:indpara),sorter)
         do,i=1,nf
            call permute(sortindex(1:indpara,i),sorter)
         end do
         
         call permute(unknowns(1:indpara),sorter)
         call permute(action(1:indpara),sorter)   
         
      end if

   end if ! end pre-permutation


   !NOW create a permutation vector which sorts the action parameter according to value 1 .... 4

   if(.not. allocated(sorter))allocate(sorter(indpara))
   
   !start index
   rowind=0
   
   !1 transform parameters in front

   do,i=1,indpara
      if(action(i) .eq. 1)then
         rowind=rowind+1
         sorter(rowind)=i
      end if
   end do

   !2 then fix parameters
   do,i=1,indpara
      if(action(i) .eq. 2)then
         rowind=rowind+1
         sorter(rowind)=i
      end if
   end do

   !3 then reduce parameters
   do,i=1,indpara
      if(action(i) .eq. 3)then
         rowind=rowind+1
         sorter(rowind)=i
      end if
   end do


   !4 then keep parameters
   do,i=1,indpara
      if(action(i) .eq. 4)then
         rowind=rowind+1
         sorter(rowind)=i
      end if
   end do



   !! We must either permute our indices now or postpone it after the transformation
   if(.not. (trans .and. nunit>0))then ! do permutation of indices now
      if(trans)then
         call permute(transort1(1:indpara),sorter)
         if(symtrans)then
            transort2(1:indpara)=transort1(1:indpara)
         else
            call permute(transort2(1:indpara),sorter)
         end if
      end if
      
      if(apriori)call permute(aprisort(1:indpara),sorter)
      if(regul)call permute(regsort(1:indpara),sorter)
      do,i=1,nf
         call permute(sortindex(1:indpara,i),sorter)
      end do
      
      call permute(unknowns(1:indpara),sorter)
      call permute(action(1:indpara),sorter)

   end if

!amount of parameters which will be fixed per system
   if(trans .and. .not. symtrans)then
      do,i=1,nf
         do,j=1,indpara
            if(action(j) .ne. 2)cycle
            if((sortindex(j,i)+transort1(j))>0)nfix(i)=nfix(i)+1
         end do
      end do
   else
      do,i=1,nf
         do,j=1,indpara
            if(action(j) .ne. 2)cycle
            if(sortindex(j,i)>0)nfix(i)=nfix(i)+1
         end do
      end do
   end if

!amount of parameters which will be reduced per system
   if(trans .and. .not. symtrans)then
      do,i=1,nf
         do,j=1,indpara
            if(action(j) .ne. 3)cycle
            if((sortindex(j,i)+transort1(j))>0)nred(i)=nred(i)+1
         end do
      end do
   else
      do,i=1,nf
         do,j=1,indpara
            if(action(j) .ne. 3)cycle
            if(sortindex(j,i)>0)nred(i)=nred(i)+1
         end do
      end do
   end if

   !calculate the amount of parameters which will be removed

   nremov=count(action <= 2) ! amount of removed parameters before the combination step
   
   !amount of parameters after transformation per system
   
   if(trans .and. .not. symtrans)then
      do,i=1,nf
         ndat2(i)=count(action(1:indpara) > 1 .and. (sortindex(1:indpara,i)+transort1(1:indpara))>0)
         !write(*,*)count(action(1:indpara)>1),count((sortindex(1:indpara,i)+transort1(1:indpara))>0)
      end do
   else
      ndat2=ndat
   end if

   !check the amount of remaining parameters when a forced combination is requested
!!! FORCED COMBINATION PARAMETERS
   if(forcecomb)then
      nfix(nf+1)=0
      ndat(nf+1)=indpara-nremov !amount of independent parameters minus amount of fixed+transformed(removed) parameters
      !set sortindex for the forced combination system
      
      forall(i=1:indpara,action(i)>2)sortindex(i,nf+1)=i-nremov!new index increasing from 1 starting at position nremov+1

      nred(nf+1)=count(action(1:indpara) .eq. 3)
   end if


   
end if ! end if extconfig


if(learn)then
   if(verbose)write(*,'(A,A)')'VERBOSE: Writing configuration file ',trim(conffile)
   unit=freeunit()
   open(unit=unit,file=trim(conffile),form='formatted',delim='QUOTE')
   write(unit,'(A)')"!!!!!!!!!!CONFIGURATION FILE FOR NORM_Tool "//trim(version)//" !!!!!!!!"
   write(unit,'(A)')"!!only change this file if you know what you're doing"
   write(unit,nml=pararedef) ! write list with redefinable parameters
   write(unit,'(A)')"!!Below are the parameters which are derived values, not recommended for changing"
   write(unit,nml=paraconf) !use namelist setup
   !then print the sorting part (cannot be in the namelist because the arrays are allocatable)
   
   if(forcecomb)then
      nd=nf+1
   else
      nd=nf
   end if
   
   !write info line
   write(unit,'(A)',advance='NO')'!!PARAMETER         ACTION '
   if(nunit>0)write(unit,'(A)',advance='NO')'PERMU '
   if(apriori) write(unit,'(A)',advance='NO')'APRISORT '
   if(regul) write(unit,'(A)',advance='NO')'REGSORT '
   if(trans)write(unit,'(A)',advance='NO')'TRANSORT1 '
   if(trans .and. .not. symtrans)write(unit,'(A)',advance='NO')'TRANSORT2 '

   write(unit,'(A)')'SORTINDEX '

   !write indices
   do,j=1,indpara
      write(unit,'(A24,1x,i2)',advance='NO')unknowns(j),action(j)
      if(nunit>0)write(unit,'(1x,i7)',advance='NO')sorter(j)
      if(apriori)write(unit,'(1x,i7)',advance='NO')aprisort(j)
      if(regul)write(unit,'(1x,i7)',advance='NO')regsort(j)
      if(trans)write(unit,'(1x,i7)',advance='NO')transort1(j)
      if(trans .and. .not. symtrans)write(unit,'(1x,i7)',advance='NO')transort2(j)
      write(unit,'(a1)',advance='NO')' '
      write(unit,*)sortindex(j,1:nd)
   end do
   close(unit)
   if(haltafterlearn)stop 
end if

!add timing point
ntimings=ntimings+1
timings(ntimings)='Configuration stage'
call cpu_time(cpu(ntimings))


!!!!!!!!!!!!!!!MEMORY ALLOCATION!!!!!!!!!!!!!!!!!!!!!!!!

!calculate memory used for the largest arrays
msize=0

!normal systems
do,i=1,nf+1
   msize=msize+ndat(i)*(ndat(i)+3)*8/(1024*1024)
end do

!transformation matrices
if(trans .and. .not. symtrans)msize=msize+2*tnval1*tnval2*8/(1024*1024)

!regularization matrix (in packed form)
if(regul .and. .not. diagregul)msize=msize+rnval*(rnval+1)*4/(1024*1024)

if(msize > 1024.d0)then !give a warning message when working array will be larger then 1Gb
   if(warning)then
      write(stderr,*)"WARNING: the memory use of your working arrays will be at least ",msize,' Mb' 
      write(stderr,*)"WARNING: Force run anyway with -w option" 
      stop
   end if
end if


!allocate memory for normal systems
do i=1,nf
   mxrw=max(ndat(i),ndat2(i))
   if(vce .or. covtrans)then
      allocate(BINinput(i)%mat1(mxrw+1,mxrw)) ! add one row to store additional diagonal
   else
      allocate(BINinput(i)%mat1(mxrw,mxrw))
   end if

!allocate memeory for vectors (minimum 2 columns)
   allocate(BINinput(i)%vec(mxrw,max(2,BINinput(i)%nvec)))
   if(warning .and. BINinput(i)%nvec < 2)then
      if(verbose)write(*,'(A)')"WARNING: Using dummy values for the righthand side!"
   end if
   !initialize vectors to 0
   BINinput(i)%vec=0.d0
!allocate space for side description
   allocate(BINinput(i)%side1_d(mxrw)) ! also enough side description data
end do

! mxrw=nda
! if(vce)mxrw=ndat(i)+1
! if(vce .or. covtrans)then ! explicitly allocate more memory for the input matrices to store original matrix in lower triangle
!    do,i=1,nf
!       allocate(BINinput(i)%mat1(ndat(i)+1,ndat(i)))
!   end do
! else if(trans .and. .not. symtrans)then
!    do,i=1,nf
!       allocate(BINinput(i)%mat1(max(ndat(i),ndat2(i)),max(ndat(i),ndat2(i))))
!    end do
! !else the routine read_BIN will allocate enough memory itself
! end if



! ! allocate enough memory for %vec component ( must be at least 2 columns) and matrix (%mat1)
! do,i=1,nf
!       allocate(BINinput(i)%vec(max(ndat(i),ndat2(i)),max(2,BINinput(i)%nvec)))
!       allocate(BINinput(i)%side1_d(max(ndat(i),ndat2(i)))) ! also enough side description data
! end do




if(forcecomb)then ! allocate memory for the combination system (not loaded from a file but created internally)
   allocate(BINinput(nf+1)%vec(ndat(nf+1),2),BINinput(nf+1)%mat1(ndat(nf+1),&
        ndat(nf+1)),BINinput(nf+1)%side1_d(ndat(nf+1)))

   !associate pointers
   normass(nf+1)%Atb=>BINinput(nf+1)%vec(:,1)
   normass(nf+1)%Apri=>BINinput(nf+1)%vec(:,2)
 !  normass(nf+1)%norm=>BINinput(nf+1)%mat1
   
   normass(nf+1)%Atb=0.d0
   normass(nf+1)%Apri=0.d0
!   BINinput(nf+1)%mat1=0.d0

end if

!!!!!!!!!!!!!!!!!!END MEMORY ALLOCATION!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!DATA READING!!!!!!!!!!!!!
if(verbose)write(*,'(A)')'VERBOSE: Reading data from normal systems'

!!! READ META data from input files if not already done!!
if(extconfig)then ! make sure to read the meta data in segments
   do,i=1,nf !first two segments per file
      BINinput(i)%file=trim(normfiles(i))
      call read_BINtype(BINinput(i),2)
   end do
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!RETRIEVE general normal system parameters from META DATA!!!!!!!!!!!!!!!!!!
!set defaults
normass(:)%ctime=9999.d0
normass(:)%etime=-9999.d0 ! very small
normass(:)%stime=9999.d0 ! very large
normass(:)%sig0=1.d0
normass(:)%ltpl=0.d0

!check if values exist already
do,i=1,nf
   do,j=1,BINinput(i)%ndbls
      if(BINinput(i)%dbls_d(j)(1:4) .eq. 'LtPL')then
         normass(i)%ltpl=BINinput(i)%dbls(j)
         
      else if(BINinput(i)%dbls_d(j)(1:6) .eq. 'Sigma0')then

         if( .not. covscale)then
            cv_sc(i)=BINinput(i)%dbls(j)**2
         end if


      else if(BINinput(i)%dbls_d(j)(1:5) .eq. 'STime')then
         normass(i)%stime=BINinput(i)%dbls(j)
      else if(BINinput(i)%dbls_d(j)(1:5) .eq. 'CTime')then
         normass(i)%ctime=BINinput(i)%dbls(j)
      else if(BINinput(i)%dbls_d(j)(1:5) .eq. 'ETime')then
         normass(i)%etime=BINinput(i)%dbls(j)
      end if
   end do

 

   do,j=1,BINinput(i)%nint
      if(BINinput(i)%ints_d(j)(1:4) .eq. 'Nobs')then
         normass(i)%nobs=BINinput(i)%ints(j)

      else if(BINinput(i)%ints_d(j)(1:8) .eq. 'Nreduced')then
        normass(i)%nred=BINinput(i)%ints(j)

     else if (BINinput(i)%ints_d(j)(1:6) .eq. 'Nfixed')then
        normass(i)%nfix=BINinput(i)%ints(j)
     end if
   end do

end do



!!INPUT NORMAL SYSTEMS ( remainder of the data) and first permutation!!!!!!!!!!!!
do,i=1,nf

   !read input system unitll the last part
   BINinput(i)%mtyp='U' ! put symmetric data in upper triangle
   call read_BINtype(BINinput(i))

   !associate pointers

   !right hand side vector points to the first column of vecec
   normass(i)%Atb=>BINinput(i)%vec(:,1)

   !apriori vector points to second column
   normass(i)%Apri=>BINinput(i)%vec(:,2)
   !
   if(apriforget)normass(i)%Apri=0.d0 ! explicitly set the apriori values to 0 ( forget apriori information)

   !resort systems according to sortindex
   ! do,j=1,indpara
   !    write(0,*)j,sortindex(j,i)
   ! end do
   call sort_BINtype(input=BINinput(i),side1=sortindex(1:indpara,i),both=.true.)

  
   if(vce .or. covtrans)then
      !copy upper triangle array to lower triangle
      !mirror plus add a 1 for the row index (that way upper and lower refer to two complete upper/lower symmteric arrays
      !so a symmetric square 3 x 3 matrix A with upper triangle U U U and lower triangle L
      !                                                           U U                    L L
      !                                                             U                    L L L
      !will be stored as:
      !  U U U
      !  L U U
      !  L L U
      !  L L L
      forall(j=1:ndat(i),k=1:ndat(i), j<=k)BINinput(i)%mat1(k+1,j)=BINinput(i)%mat1(j,k)


      if(vce)then
         !copy old parameters
         normass(i)%ltplold=normass(i)%ltpl
         normass(i)%nobsold=normass(i)%nobs
         allocate(normass(i)%Atbold(size(normass(i)%Atb,1)),normass(i)%Apriold(size(normass(i)%Apri,1)))
         forall(j=1:ndat(i))normass(i)%Atbold(j)=normass(i)%Atb(j)
         forall(j=1:ndat(i))normass(i)%Apriold(j)=normass(i)%Apri(j)
      end if
   end if

   !!APPLY sigma0
   !rescale right hand side with covariance scale
   if(cv_sc(i) .ne. 1.d0 .or. W(i) .ne. 0.d0)then
      Normass(i)%Atb=Normass(i)%Atb*(W(i)/cv_sc(i))
   end if
   
   !rescale normal matrix  and LtPL
   if(cv_sc(i) .ne. 1.d0)then
      !rescale ltpl
      normass(i)%ltpl=normass(i)%ltpl/cv_sc(i)
      if(covtrans)then ! rescale complete matrix
         forall(j=1:ndat(i)+1,k=1:ndat(i))BINinput(i)%mat1(j,k)=BINinput(i)%mat1(j,k)/cv_sc(i)
      else !(upper triangle only)
         forall(j=1:ndat(i),k=1:ndat(i), j<=k)BINinput(i)%mat1(j,k)=BINinput(i)%mat1(j,k)/cv_sc(i)
      end if
   end if


   !REDEFINE SORTINDEX (ever increasing with zeros for non-existing parameters)
   rowind=0
   do,j=1,indpara
      if(sortindex(j,i) .eq. 0)cycle
      rowind=rowind+1
      sortindex(j,i)=rowind
   end do

end do

!add timing point
ntimings=ntimings+1
timings(ntimings)='Data reading and sorting'
call cpu_time(cpu(ntimings))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!APPLY TRANSFORMATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(trans)then
   if(verbose)write(*,'(A)')"VERBOSE: transforming parameters..."
   
   if(diagtrans)then !diagonal transformation
      
      !load transformation matrix
      allocate(TRANSinput%pack1(tnval1))
      unit=freeunit()
      open(unit=unit,file=trim(tr_file),form='formatted',status='old')
      do,i=1,tnval1
         read(unit=unit,fmt='(A200)')dum(1:200)
         read(dum(25:),*)TRANSinput%pack1(i)

!          read(unit=unit,fmt='(24x)',advance='NO')
!          read(unit=unit,fmt=*)TRANSinput%pack1(i)
      end do
      close(unit)
      
      allocate(B(tnval1,1)) !small only one column needed
      
      do,i=1,nf!loop over input files
         !count the number of values which are transformed through B
         !fill the transformation matrix B
         nb=0 ! amount of parameters in system i which wil be transformed through B
         do,j=nunit+1,indpara
            if(sortindex(j,i) .eq. 0)cycle ! ignore values which are not present
            nb=nb+1
            B(nb,1)=TRANSinput%pack1(transort2(j))
         end do
         !nb is now the amount of parameters in system i which wil be transformed through B
         if(nb .eq. 0)cycle ! no use to transform this
         shift=ndat(i)-nb

         if(covtrans)then ! error covariance transformation
                        
            call trans_cov(B=B,ldb=size(B,1),A=BINinput(i)%mat1,&
                 lda=size(BINinput(i)%mat1,1),atb=normass(i)%atb,&
                 ltpl=normass(i)%ltpl,nunit=shift,nb=nb,type=1)
            !update the amount of observations (previous levels are forgotten/irretrievable)
            normass(i)%nobs=ndat2(i)
            normass(i)%nred=0 ! added 23 Feb 2012
         else ! parameter transformation
            if(.not. atb_ignore )then
               call diag_transf(B=B(:,1),atb=normass(i)%Atb,nunit=shift,nb=nb,add=transadd)
            end if
            
            if( .not. N_ignore)then
               call diag_transf(B=B(:,1),A=BINinput(i)%mat1,nunit=shift,nb=nb,add=transadd)
            end if
            
            if(x0_as_b)then
               call diag_transf(B=B(:,1),atb=normass(i)%Apri,nunit=shift,nb=nb,add=transadd)
            else if( .not. x0_ignore)then
               call diag_transf(B=B(:,1),x0=normass(i)%Apri,nunit=shift,nb=nb,add=transadd)
            end if
         end if
      end do
   else if(symtrans)then !symmetric matrix transformation
      !load transformation matrix
      TRANSinput%mtyp='U' ! load data in upper triangle
      TRANSinput%file=trim(tr_file)
      call read_BINtype(TRANSinput)

      allocate(B(tnval1,tnval1))

      do,i=1,nf !loop over normal systems
         
         

         nb=count(sortindex(nunit+1:indpara,i)>0) ! amount of parameters which are transformed through B
         if(nb .eq. 0)cycle ! no use to transform this
         shift=ndat(i)-nb ! amount of parameters which are transformed through the unit matrix

         !copy and compact array in transformation matrix
         forall(j=nunit+1:indpara,k=nunit+1:indpara,(sortindex(j,i)*sortindex(k,i))>0.and. j<=k)&
            B(sortindex(j,i),sortindex(k,i))=TRANSinput%mat1(transort1(j),transort2(k))

         !mirror transformation matrix
         forall(j=1:nb,k=1:nb,j<k)B(k,j)=B(j,k)

         if(covtrans)then ! error covariance transformation
            call trans_cov(B=B,ldb=size(B,1),A=BINinput(i)%mat1,&
                 lda=size(BINinput(i)%mat1,1),atb=normass(i)%atb,&
                 ltpl=normass(i)%ltpl,nunit=shift,nb=nb,type=2)
            !update the amount of observations (previous levels are forgotten/irretrievable)
            normass(i)%nobs=ndat2(i)
         else
            if(.not. atb_ignore )then
               
               call sym_transf(B=B,ldb=size(B,1),lda=0,atb=normass(i)%Atb,nunit=shift,nb=nb)
               
            end if
            
            if( .not. N_ignore)then
               call sym_transf(B=B,ldb=size(B,1),A=BINinput(i)%mat1,lda=size(BINinput(i)%mat1,1),nunit=shift,nb=nb)
            end if
            
            if(x0_as_b)then
               call sym_transf(B=B,ldb=size(B,1),lda=0,atb=normass(i)%Apri,nunit=shift,nb=nb)
            else if( .not. x0_ignore)then
               call sym_transf(B=B,ldb=size(B,1),lda=0,x0=normass(i)%Apri,nunit=shift,nb=nb)
            end if
         end if

      end do
   else !full transformation ( possibly eliminating and introducing new variables)


      TRANSinput%mtyp='F' ! full matrix in 2D
      TRANSinput%file=trim(tr_file)
      call read_BINtype(TRANSinput)

      allocate(B(tnval1,tnval2))
      B=0.d0
      !calculate new sorter index for rows after transformation
!      if(.not. allocated(sorter))allocate(sorter(indpara))
      
      do,i=1,nf
 !        sorter=sortindex(:,i)!copy old values straight away
         
         nb=count(sortindex(nunit+1:indpara,i)>0) ! amount of parameters which are transformed through B
         if(nb .eq. 0)cycle ! nothing to transform
         nk=count(transort1(1:indpara)>0)! amount of valid rows in B
         shift=ndat(i)-nb ! amount of parameters transformed with a unit transformation



         !fill up transformation matrix B
         if(Btrans)then ! B is transposed
            do,j=nunit+1,indpara
               if(sortindex(j,i) .eq. 0)cycle
               rowind=0
               do,k=1+nunit,indpara
                  if(transort1(k) .eq. 0)cycle
                  rowind=rowind+1
                  B(rowind,sortindex(j,i)-nunit)=TRANSinput%mat1(transort2(j),transort1(k))
               end do
            end do
         else !B is not transposed
            do,j=nunit+1,indpara
               if(sortindex(j,i) .eq. 0)cycle
               rowind=0
               do,k=1+nunit,indpara
                  if(transort1(k) .eq. 0)cycle
                  rowind=rowind+1
!                   write(0,*)transort1(k),transort2(j),size(TRANSinput%mat1,1),size(TRANSinput%mat1,2)
!                   write(0,*)TRANSinput%mat1(transort1(k),transort2(j))
                  B(rowind,sortindex(j,i)-nunit)=TRANSinput%mat1(transort1(k),transort2(j))
               end do
            end do
         end if

!          write(0,*)'test'
!          stop

         if(covtrans)then ! error covariance transformation
            call trans_cov(B=B,ldb=size(B,1),A=BINinput(i)%mat1,&
                 lda=size(BINinput(i)%mat1,1),atb=normass(i)%atb,&
                 ltpl=normass(i)%ltpl,nunit=nunit,nb=nb,type=3)
            !update the amount of observations (previous levels are forgotten/irretrievable)
            normass(i)%nobs=ndat2(i)
         else
            
            if(.not. atb_ignore )then
               call full_transf(B=B,ldb=size(B,1),lda=0,atb=normass(i)%Atb,nunit=shift,nb=nb,k=nk,add=transadd)
            end if
            
            if( .not. N_ignore)then
               call full_transf(B=B,ldb=size(B,1),A=BINinput(i)%mat1,&
                    lda=size(BINinput(i)%mat1,1),nunit=shift,nb=nb,k=nk,add=transadd)
               
            end if

            if(x0_as_b)then
               call full_transf(B=B,ldb=size(B,1),lda=0,atb=normass(i)%Apri,nunit=shift,nb=nb,k=nk,add=transadd)
            else if( .not. x0_ignore)then
               call full_transf(B=B,ldb=size(B,1),lda=0,x0=normass(i)%Apri,nunit=shift,nb=nb,k=nk,add=transadd)
            end if
         end if
         ndat(i)=ndat2(i)
         BINinput(i)%nval1=ndat2(i) ! also update (BUGfix? 23-04-2010)
         BINinput(i)%nval2=ndat2(i)
            !REDEFINE SORTINDEX (ever increasing with zeros for non-existing parameters)
            !some parameters might have been added/and or removed through transformation

               rowind=0
               do,j=1,indpara
                  if(action(j) .eq. 1)then ! removed through transformation
                     sortindex(j,i)=0
                     cycle !this one is removed now
                  end if
                  if((sortindex(j,i)+ transort1(j)).eq. 0)cycle ! cycle when both indices are zero (not added or kept)
                  rowind=rowind+1
                  sortindex(j,i)=rowind
               end do

      end do ! end loop over files

     !!!SECOND PERMUTATION if necessary
      
      if( nunit>0)then ! we need an additional permutation
         if(verbose)write(*,'(A)')'VERBOSE: Applying second permutation'
         if(.not. symtrans)call permute(transort1(1:indpara),sorter) ! permuted version is still needed
         if(apriori)call permute(aprisort(1:indpara),sorter)
         if(regul)call permute(regsort(1:indpara),sorter)
         call permute(unknowns(1:indpara),sorter)
         call permute(action(1:indpara),sorter)
         do,i=1,nf
            !resort systems according to new sortindex
            call permute(sortindex(1:indpara,i),sorter)
            
            call sort_BINtype(input=BINinput(i),side1=sortindex(:,i),both=.true.)
            
            !REDEFINE SORTINDEX (ever increasing with zeros for non-existing parameters)

               rowind=0
               do,j=1,indpara
                  if(sortindex(j,i) .eq. 0)cycle
                  rowind=rowind+1
                  sortindex(j,i)=rowind
               end do
            
         end do

      end if
      

   end if

   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='Transformation stage'
   call cpu_time(cpu(ntimings))
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END APPLY TRANSFORMATION!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!APPLY APRIORI SOLUTION SUBSTITUTE!!!!!!!!!!!!!!!!


if(apriori)then
   if(verbose)write(*,'(A)')'VERBOSE: Reading new apriori values from file',trim(solfile)
   if(trim(solfile).eq. '_ZERO_')then ! special case (resets the apriori values to zero but change the system take the negated a priori values already present (first system only)
      allocate(APRIinput%vec(snval,1))
      APRIinput%vec(1:snval,1)=0.d0
   else if(apri_ascii)then
      allocate(APRIinput%vec(snval,1))
      unit=freeunit()
      open(unit=unit,file=trim(solfile),form='formatted')
         do i=1,snval

            read(unit=unit,fmt='(A200)')dum(1:200)
            read(dum(25:),*)APRIinput%vec(i,1)

         end do
         close(unit)
   else
      APRIinput%file=trim(solfile)
      call read_BINtype(APRIinput,3) ! read up to the third segment of the solutionfile (vectors)
   end if

   if(.not. allocated(tmpvec2) .and. .not. aprireplace)allocate(tmpvec2(indpara))




   do,i=1,nf !loop over all input systems (loop maybe parallized)
      napri=count((aprisort > 0) .and. (sortindex(1:indpara,i) > 0))
      if(napri .eq. 0)cycle !nothing to update in this system
      
      nd=ndat(i)

      if(aprireplace)then ! just replace apriori values

         if(apriadd)then ! add the already existing a priori vector to the one just read
            do,j=1,indpara
               if((aprisort(j) .eq. 0) .or. (sortindex(j,i) .eq. 0))cycle
               normass(i)%Apri(sortindex(j,i))=normass(i)%Apri(sortindex(j,i))+APRIinput%vec(aprisort(j),1) !reset new apriori value
            end do
         else

            do,j=1,indpara
               if((aprisort(j) .eq. 0) .or. (sortindex(j,i) .eq. 0))cycle
               normass(i)%Apri(sortindex(j,i))=APRIinput%vec(aprisort(j),1) !reset new apriori value
            end do
         end if
      else ! also update ltpl and the right hand side

         tmpvec2=0.d0 !initialize to zero (difference between old and new apriori values)

         if(apriadd)then !
            do,j=1,indpara
               if((aprisort(j) .eq. 0) .or. (sortindex(j,i) .eq. 0))cycle
               tmpvec2(sortindex(j,i))=APRIinput%vec(aprisort(j),1) !tmpvec2 now holds the  a (differences) priori values 
               normass(i)%Apri(sortindex(j,i))=normass(i)%Apri(sortindex(j,i))+APRIinput%vec(aprisort(j),1) !reset new apriori value
            end do
         else
            !make a difference vector and put new apriori values in the apriori estimate
            do,j=1,indpara
               if((aprisort(j) .eq. 0) .or. (sortindex(j,i) .eq. 0))cycle
               tmpvec2(sortindex(j,i))=APRIinput%vec(aprisort(j),1)-normass(i)%Apri(sortindex(j,i)) !tmpvec2 now holds the difference
               normass(i)%Apri(sortindex(j,i))=APRIinput%vec(aprisort(j),1) !reset new apriori value
            end do
         end if
         call setapriori_norm2(C=BINinput(i)%mat1,d=normass(i)%Atb(1:nd),&
              dx=tmpvec2(1:nd),ltpl=normass(i)%ltpl,ldc=ndat(i))

      end if
   end do

   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='Setting apriori values'
   call cpu_time(cpu(ntimings))

end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIXING PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!no real action is performed on fixing parameters, the sorting and consequently selecting the relevant portion of the original normal equations by means of shifting by nfix(i) will cause this.
if(fix .and. verbose)write(*,'(A)')"VERBOSE: fixing parameters"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END FIXING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COMBINING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(combine)then
   if(verbose)write(*,'(A)')"VERBOSE: combining parameters"

   if(CombapriAverage)totalAveweight=0.d0
   do,i=1,nf
      if(comb(i) .eq. i) cycle !no need to add a system to itself
      nd=ndat(i)
      
      !normal system
!       if(i .eq. 1)then
!       do,j=1,indpara
!          write(*,*)sortindex(j,:)
!       end do
!    end if
      !matrix
      forall(j=nremov+1:indpara,k=nremov+1:indpara,(j<=k) .and. (sortindex(j,i)> 0) .and. &
           (sortindex(k,i) >0) .and. (sortindex(j,comb(i))> 0) .and. (sortindex(k,comb(i))>0) )
         BINinput(comb(i))%mat1(sortindex(j,comb(i)),sortindex(k,comb(i)))=&
         BINinput(comb(i))%mat1(sortindex(j,comb(i)),sortindex(k,comb(i)))+&
              BINinput(i)%mat1(sortindex(j,i),sortindex(k,i))
      end forall


      !right hand vector
      forall(j=nremov+1:indpara,(sortindex(j,i)> 0) .and.(sortindex(j,comb(i))> 0))
         normass(comb(i))%Atb(sortindex(j,comb(i)))=normass(comb(i))%Atb(sortindex(j,comb(i)))+&
              normass(i)%Atb(sortindex(j,i))
      end forall
           
      if(forcecomb)then !copy apriori estimate (some may be copied more times in fact)
         if(CombapriAverage)then ! average the apriori vectors
            Aveweight=normass(i)%Etime-normass(i)%Stime
            if(Aveweight <=0.d0)then
               write(stderr,*)"ERROR: zero or negative averaging weight for the apriori vector of system",i
               stop
            end if

            totalAveweight=totalAveweight+Aveweight
!            write(stderr,*)Aveweight,totalAveweight
            if(i==1)then !copy on first element
               forall(j=nremov+1:indpara,(sortindex(j,i)> 0) .and.(sortindex(j,comb(i))> 0))
                  normass(comb(i))%Apri(sortindex(j,comb(i)))=normass(i)%Apri(sortindex(j,i))
               end forall
            else !update the weighted mean of the apriori vector
               forall(j=nremov+1:indpara,(sortindex(j,i)> 0) .and.(sortindex(j,comb(i))> 0))
                  normass(comb(i))%Apri(sortindex(j,comb(i)))=normass(comb(i))%Apri(sortindex(j,comb(i)))+&
                    (Aveweight/totalAveweight)*(normass(i)%Apri(sortindex(j,i))- normass(comb(i))%Apri(sortindex(j,comb(i))))
               end forall
            end if
         else !just copy the values
            forall(j=nremov+1:indpara,(sortindex(j,i)> 0) .and.(sortindex(j,comb(i))> 0))
               normass(comb(i))%Apri(sortindex(j,comb(i)))=normass(i)%Apri(sortindex(j,i))
            end forall
         end if
      end if

      !ltpl

      normass(comb(i))%ltpl=normass(comb(i))%ltpl+normass(i)%ltpl
    
      !start time
      normass(comb(i))%stime=min(normass(comb(i))%stime,normass(i)%stime)
      
      !end time
      normass(comb(i))%etime=max(normass(comb(i))%etime,normass(i)%etime)
      
      !amount of observations
      normass(comb(i))%nobs=normass(comb(i))%nobs+normass(i)%nobs
      
      !total amount of fixed parameters ( from previous operation)
      normass(comb(i))%nfix=normass(comb(i))%nfix+normass(i)%nfix

      !total amount of reduced parameters ( from previous operation)
      normass(comb(i))%nred=normass(comb(i))%nred+normass(i)%nred

      !also copy old meta data (integers)
      do,j=1,BINinput(i)%nint
         call BIN_putmeta(in=BINinput(comb(i)),dum='I='//trim(BINinput(i)%ints_d(j))//'/',iint=BINinput(i)%ints(j))
      end do

      !also copy old meta data (doubles)
      do,j=1,BINinput(i)%ndbls
         call BIN_putmeta(in=BINinput(comb(i)),dum='D='//trim(BINinput(i)%dbls_d(j))//'/',idbl=BINinput(i)%dbls(j))
      end do

   end do

   !set center time (after above loop is finished) (average of start end end time)
   do,i=1,nf+1
      if(comb(i) .ne. i)cycle
      normass(comb(i))%ctime=(normass(comb(i))%stime+normass(comb(i))%etime)/2
   end do

   !!!!!!!!!!!!!!!APRIORI CHECK (CHECK WHETHER APRIORI VALUES OF THE COMBINED SYSTEMS ARE THE SAME)!!!!!!!
   if(warning .and. .not. CombapriAverage)then
      if(verbose)write(*,'(A)')"VERBOSE: Checking whether apriori values of the normal systems agree"
      do,i=1,nf
         if(comb(i) .eq. i) cycle !No need to check itself
         
         do,j=nremov+1,indpara !ignore fixed and transformed parameters
            if(sortindex(j,i) .eq. 0)cycle
            if(abs(normass(i)%Apri(sortindex(j,i))-normass(comb(i))%Apri(sortindex(j,comb(i)))) >1.d-30)then
               !   write(*,*)normass(i)%side_d(sortindex(j,i)),normass(comb(i))%side_d(sortindex(j,comb(i)))
               !   write(*,*)j,normass(i)%Apri(sortindex(j,i)),normass(comb(i))%Apri(sortindex(j,comb(i)))
               
               write(stderr,*)" WARNING: apriori values of combined systems",i,"and",comb(i),"don't agree"
               write(stderr,*)" Use -w to force run"
               stop
               
            end if
         end do
         
      end do
      if(verbose)write(*,'(A)')"VERBOSE: CHECK OK"
   end if



   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='Systems combination'
   call cpu_time(cpu(ntimings))
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END combining!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REDUCING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(redu)then
   if(verbose)write(*,'(A)')"VERBOSE: Reducing normal systems"
   do,i=1,nf+1 
      if(comb(i) .ne. i)cycle !only reduce remaining systems after combination
      !temporarily assign pointers to the right array sections
      P1=>BINinput(i)%mat1(nfix(i)+1:,nfix(i)+1:)
      P2=>BINinput(i)%mat1(nfix(i)+1:,nfix(i)+nred(i)+1:)
      P3=>BINinput(i)%mat1(nfix(i)+nred(i)+1:,nfix(i)+nred(i)+1:)
      call reduce_norm2(C11=P1,ldc11=ndat(i),C12=P2,ldc12=ndat(i),&
           C22=P3,ldc22=ndat(i),d1=normass(i)%Atb(nfix(i)+1:nfix(i)+nred(i)),&
           d2=normass(i)%Atb(nfix(i)+nred(i)+1:ndat(i)),ltpl=normass(i)%ltpl)

   end do
   
   
   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='Reducing parameters'
   call cpu_time(cpu(ntimings))
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END REDUCING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REGULARIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(regul)then
   if(verbose)write(*,*)"VERBOSE: Regularizing normal systems with scale:",regscale
   
   !read regularization matrix

   if( .not. diagregul)then
      

      REGinput%file=trim(regfile)
      call read_BINtype(REGinput)
      
   
      do,i=1,nf+1 !loop over all input systems (loop maybe parallized)
         if(comb(i) .ne. i)cycle !only regularize final systems (after combination)
         
         nreg=count((regsort > 0) .and. (sortindex(1:indpara,i) > 0))
         if(nreg .eq. 0)cycle !nothing to regularize in this system
         
         !regularize matrix
         shift=nfix(i)+nred(i)+1
         
         do,k=shift,indpara
            if(sortindex(k,i) .eq. 0 .or. regsort(k) .eq. 0 )cycle ! parameter not present ( in systesm or regularization matrix)
            do,j=shift,k
               if(sortindex(j,i) .eq. 0 .or. regsort(j) .eq. 0)cycle ! parameter not present ( in system or regularization matrix)
               ind=packindex(REGinput,min(regsort(j),regsort(k)),max(regsort(j),regsort(k)))
               if(ind .eq. 0)cycle ! not in block or triangle
               BINinput(i)%mat1(sortindex(j,i),sortindex(k,i))=&
                    BINinput(i)%mat1(sortindex(j,i),sortindex(k,i))+&
                    regscale*REGinput%pack1(ind)
            end do
         end do

!          forall(j=shift:indpara,k=shift:indpara,( j<=k) .and. (sortindex(j,i) > 0) .and. &
!               (sortindex(k,i) >0) .and. (regsort(j) > 0) .and. (regsort(k) >0))
!             BINinput(i)%mat1(sortindex(j,i),sortindex(k,i))=BINinput(i)%mat1(sortindex(j,i),sortindex(k,i))+&
!                  regscale*REGinput%pack1(min(regsort(j),regsort(k))+max(regsort(j),regsort(k))*(max(regsort(j),regsort(k))-1)/2)
!          end forall
   
      end do
   else if(diagregul .and. regfile .ne. '')then !diagonal regularization with ascii file
      unit=freeunit()

      allocate(REGinput%vec(rnval,1))
      !read diagonal regularization data
      open(unit=unit,file=trim(regfile),form='formatted')
       
      do,i=1,rnval
         read(unit=unit,fmt='(a24)',ADVANCE='NO')dum(1:24)
         read(unit=unit,fmt=*)REGinput%vec(i,1)
      end do
      close(unit)


      do,i=1,nf+1
         if(comb(i) .ne. i) cycle
         forall(j=1:indpara, (regsort(j) .ne. 0) .and. (sortindex(j,i) >0))
            BINinput(i)%mat1(sortindex(j,i),sortindex(j,i))=BINinput(i)%mat1(sortindex(j,i),sortindex(j,i))&
                 +regscale*REGinput%vec(regsort(j),1)
         end forall
      end do

   else !if unit diagonal regularization (no external file)
      do,i=1,nf+1
         if(comb(i) .ne. i) cycle
         forall(j=1:indpara, (regsort(j) .ne. 0) .and. (sortindex(j,i) >0))
            BINinput(i)%mat1(sortindex(j,i),sortindex(j,i))=BINinput(i)%mat1(sortindex(j,i),sortindex(j,i))+regscale
         end forall
      end do
   end if

   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='Regularization'
   call cpu_time(cpu(ntimings))
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END REGULARIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SOLVING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(solve)then
   if(verbose)write(*,'(A)')"VERBOSE: Solving final systems"
   if(verbose .and. truncsvd )write(*,'(A)')"using TSVD"
   do,i=1,nf+1
      if(comb(i) .ne. i)cycle  

      
      shift=nfix(i)+nred(i)
      !temporarily associate pointer with array section
      P1=>BINinput(i)%mat1(shift+1:,shift+1:)
      if(truncsvd)then
         call solve_tsvd(C=P1,ldc=size(BINinput(i)%mat1,1),&
           n=ndat(i)-shift,d=normass(i)%Atb(shift+1:),ltpl=normass(i)%ltpl,ntrunc=ntrunc)
      else
         call solve_norm2(C=P1,ldc=size(BINinput(i)%mat1,1),&
           n=ndat(i)-shift,d=normass(i)%Atb(shift+1:),ltpl=normass(i)%ltpl)
      end if
      !a posteriori covariance scale (also use reduced varaibles from previous operations)
      normass(i)%sig0=sqrt(normass(i)%ltpl/(normass(i)%nobs-(ndat(i)-nfix(i)+normass(i)%nred)))
      
   end do

   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='Solving systems'
   call cpu_time(cpu(ntimings))
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END SOLVING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!VCE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(vce)then
   if(verbose)write(*,'(A)')"VERBOSE: estimating variance components of the root systems"
   if(.not. allocated(tmpvec2))allocate(tmpvec2(indpara))
   do,i=1,nf
      !check whether  destination system is a combination
      nsys=count(comb .eq. comb(i))

      if(nsys .eq. 1)cycle !cycle when there was only one root system (contribution is 1 per definition, vce component was calculated above, during solving step)
     
      
!!!!!!!!!!!CONTRIBUTION VECTOR!!!!!!!!!!

      allocate(normass(i)%contr(ndat(i)))
      normass(i)%contr=0.d0 

      !associate pointer with the original normal system (lower part with rowindex incremented by 1)
      
      upper => BINinput(comb(i))%mat1(1:ndat(comb(i)),1:ndat(comb(i))) !refers to the combined covariance matrix

      lower => BINinput(i)%mat1(2:ndat(i)+1,1:ndat(i))
      
      trc1=0.d0
      
      do,j=1,indpara
         if(sortindex(j,i) .eq. 0)cycle
         trc1=trc1+(lower(sortindex(j,i),sortindex(j,i)))!/ndat(i)
      end do

      trc2=0.d0
      do,j=1,indpara
         if(sortindex(j,comb(i)) .eq. 0)cycle
         trc2=trc2+upper(sortindex(j,comb(i)),sortindex(j,comb(i)))!/ndat(comb(i))
      end do

      

      !calculate contribution vector
!      write(0,*)trc1,trc2,trc1*trc2

      do,j=1,indpara
         if(sortindex(j,i) .eq. 0)cycle
         
         !dotproduct
         do,k=1,indpara
            if(sortindex(k,i) .eq. 0)cycle
            min1=min(sortindex(j,comb(i)),sortindex(k,comb(i)))
            max1=max(sortindex(j,comb(i)),sortindex(k,comb(i)))
            min2=min(sortindex(j,i),sortindex(k,i))
            max2=max(sortindex(j,i),sortindex(k,i))
            
            normass(i)%contr(sortindex(j,i))=normass(i)%contr(sortindex(j,i))+(upper(min1,max1)/trc2)*(lower(max2,min2)/trc1)
           ! normass(i)%contr(sortindex(j,i))=normass(i)%contr(sortindex(j,i))+(upper(min1,max1))*(lower(max2,min2))
            
         end do
         
         normass(i)%contr(sortindex(j,i))=normass(i)%contr(sortindex(j,i))*trc1*trc2/cv_sc(i) !also apply sigma0 to the original values
         !normass(i)%contr(sortindex(j,i))=normass(i)%contr(sortindex(j,i))

     
      end do


  !    write(*,*)'step2'

	!!!!!!!!!!!!!VCE!!!!!!!!!!!!!

      !!update apriori vector with the just estimated solution vector


      do,j=1,indpara
         if(sortindex(j,i) .eq. 0)cycle

         tmpvec2(sortindex(j,i))=normass(comb(i))%Atb(sortindex(j,comb(i)))

       !adjust new apriori estimate
          normass(i)%Apriold(sortindex(j,i))=normass(i)%Apriold(sortindex(j,i))+tmpvec2(sortindex(j,i))
      end do

      P1=>BINinput(i)%mat1(2:,1:)
      call setapriori_norm2(C=P1,ldc=ndat(i)+1,d=normass(i)%ATbold,dx=tmpvec2(1:ndat(i))&
           ,ltpl=normass(i)%ltplold,uplo='L')


      trace=sum(normass(i)%contr)
!!      normass(i)%vce=normass(i)%ltplold/(normass(i)%nobsold-trace/cv_sc(i)) !removed bug again a division by cv_sc??
      normass(i)%vce=normass(i)%ltplold/(normass(i)%nobsold-trace) !removed bug again a division by cv_sc??

      if(verbose)then
         write(*,*)'VERBOSE: VCE, sigma0 of system',i,sqrt(normass(i)%vce),'convergence(--->1):',sqrt(normass(i)%vce/cv_sc(i))

      end if
      !print data to file

      unit=freeunit()

      write(contfile,'(A6,i3.3,a4,i3.3,a4)')'CONTR_',i,'_to_',comb(i),'.txt'
      open(unit=unit,file=trim(contfile),form='formatted')
      write(unit,*)'Contribution of system',i, 'to system',comb(i)
      write(unit,*)'VCE scale:',sqrt(normass(i)%vce)
      do,j=1,indpara
         if(sortindex(j,i) .eq. 0)cycle
         write(unit,'(A24,1x,F12.10)')unknowns(j),normass(i)%contr(sortindex(j,i))
      end do
      close(unit)
      

      
   end do
   
   !add timing point
   ntimings=ntimings+1
   timings(ntimings)='VCE'
   call cpu_time(cpu(ntimings))

end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END VCE!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WRITE to FILE/FIFO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(vcereloop)then !output rescaled normal systems with updated
   if(verbose)write(*,'(A)')"VERBOSE: outputting adjusted normal systems"
   
   nfout=nf

   do,i=1,nf

      !Fill out BIndat type

!      File name
      if(nfout .ne. 1 .and. trim(basename) .ne. 'stdout')then
         write(chardum,'(a1,i3.3)')'_',i
      else
         chardum=''
      end if      
      !file name
      BINinput(i)%file=trim(basename)//chardum
      !system type
      BINinput(i)%type='SYMVN___'
      !matrix type
      BINinput(i)%mtyp='L'
      !short description 
      BINinput(i)%descr=descr

      !adapt original integer meta data
      do,j=1,BINinput(i)%nint
         if(BINinput(i)%ints_d(j)(1:5) .eq. 'Modnr')BINinput(i)%ints(i)=BINinput(i)%ints(i)+1
      end do

      !adapt original double meta data
      do,j=1,BINinput(i)%ndbls
         if(BINinput(i)%dbls_d(j)(1:4) .eq. 'LtPL')then
            BINinput(i)%dbls(j)=normass(i)%ltplold
         else if(BINinput(i)%dbls_d(j)(1:6) .eq. 'Sigma0')then
            BINinput(i)%dbls(j)=sqrt(normass(i)%vce)
         end if
      end do

      !side description data
      forall(j=1:indpara,sortindex(j,i)>0)BINinput(i)%side1_d(sortindex(j,i))=unknowns(j)


      !do a segmented writein order to facilitate "simultaneaous" piping of more files
      !first two segments segment index data plus meta data 

      call write_BINtype(BINinput(i),2)

   end do

   !third segment vector and meta data
   do,i=1,nf

      !copy adapted vector data
      do,j=1,ndat(i)
         BINinput(i)%vec(j,1)=normass(i)%Atbold(j)
         BINinput(i)%vec(j,2)=normass(i)%Apriold(j)
      end do

!      adapt matrix data ( tolerate the associate memory leak)
      !make a temporary pointer to the full matrix

      BINinput(i)%mat1=>BINinput(i)%mat1(2:ndat(i)+1,1:ndat(i)) !points to lower section only!!!
      !information on the first row is lost ( tolerated memory leak)

      call write_BINtype(BINinput(i))

   end do
   
   stop !don't output final files any more

end if




if(verbose)write(*,'(A)')"VERBOSE: outputting final systems"
nfout=0
do,i=1,nf+1
   if(comb(i) .ne. i)cycle
   nfout=nfout+1
end do

do,i=1,nf+1
   if(comb(i) .ne. i) cycle !only output final systems
   shift=nfix(i)+nred(i)
   ndbig=ndat(i)-shift ! amount of remaining data
  ! write(*,*)shift,nfix(i),nred(i),nd

   !construct output system filename

   if(nfout .ne. 1 .and. trim(basename) .ne. 'stdout')then
      write(chardum,'(a1,i3.3)')'_',i
   else
      chardum=''
   end if

   !file name
   BINinput(i)%file=trim(basename)//chardum
   if(verbose)write(*,'(A,A)')"VERBOSE: writing system: ",trim(BINinput(i)%file)
   !system type
   BINinput(i)%type='SYMVN___'
   !matrix type
   BINinput(i)%mtyp='U'
   !short description 
   BINinput(i)%descr=descr
   
   !sizes
   BINinput(i)%nval1=ndbig
   BINinput(i)%nval2=BINinput(i)%nval1
   BINinput(i)%pval1=(ndbig*(ndbig+1))/2
   BINinput(i)%pval2=1
   BINinput(i)%nvec=2
   BINinput(i)%nread=0
   
  !update Integer meta data
   call BIN_putmeta(in=BINinput(i),dum='I=Nobs/',iint=normass(i)%nobs)
   call BIN_putmeta(in=BINinput(i),dum='I=Nreduced/',iint=normass(i)%nred+nred(i))
   call BIN_putmeta(in=BINinput(i),dum='I=Nfixed/',iint=normass(i)%nfix+nfix(i))
   call BIN_putmeta(in=BINinput(i),dum='I=Nunknowns/',iint=ndbig)


   !also update double data
   call BIN_putmeta(in=BINinput(i),dum='D=CTime/',idbl=normass(i)%Ctime)
   call BIN_putmeta(in=BINinput(i),dum='D=STime/',idbl=normass(i)%Stime)
   call BIN_putmeta(in=BINinput(i),dum='D=ETime/',idbl=normass(i)%Etime)
   call BIN_putmeta(in=BINinput(i),dum='D=LtPL/',idbl=normass(i)%ltpl)
   call BIN_putmeta(in=BINinput(i),dum='D=Sigma0/',idbl=normass(i)%sig0)

   ! outints(1)=normass(i)%nobs
   ! outints(2)=nd
   ! outints(3)=nred(i)+normass(i)%nred ! this call and previous calls 
   ! outints(4)=nfix(i)+normass(i)%nfix ! this call and previous calls 
   ! outints(5)=1
   ! BINinput(i)%nint=5
   ! BINinput(i)%ints=>outints
   ! BINinput(i)%ints_d=>outints_d
   ! !double meta data
   ! outdbls(1)=normass(i)%Ctime
   ! outdbls(2)=normass(i)%Stime
   ! outdbls(3)=normass(i)%Etime
   ! outdbls(4)=normass(i)%ltpl
   ! outdbls(5)=normass(i)%sig0
   ! BINinput(i)%ndbls=5
   ! BINinput(i)%dbls=>outdbls
   ! BINinput(i)%dbls_d=>outdbls_d

   ! copy side description data
   forall(j=1:indpara,sortindex(j,i)>0 .and. action(j).eq. 4 )BINinput(i)%side1_d(sortindex(j,i)-shift)=&
           unknowns(j)

   !write first two segments
   call write_BINtype(BINinput(i),2)

  if(verbose)write(*,'(A)')'VERBOSE:first two segments'
end do

!write third and fourth segments
do,i=1,nf+1
   if( comb(i) .ne. i)cycle

   shift=nfix(i)+nred(i)
   nd=ndat(i)-shift

   !vector data
   BINinput(i)%vec=>BINinput(i)%vec(shift+1:ndat(i),1:2)

   !matrix data
   
   BINinput(i)%mat1=>BINinput(i)%mat1(shift+1:ndat(i),shift+1:ndat(i))

   
   call write_BINtype(BINinput(i))

   
 if(verbose)write(*,'(A)')'VERBOSE:third and fourth segments'
end do
  

!add timing point
ntimings=ntimings+1
timings(ntimings)='File Output'
call cpu_time(cpu(ntimings))
if(verbose)then
   cpu_tot=cpu(ntimings)-cpu(1)
   write(*,'(A)')'VERBOSE: Timing report:'
   write(*,'(A,F10.2)')'VERBOSE: Total computing time (seconds):',cpu_tot
   write(*,'(A,i2,A)')'VERBOSE: Division in',ntimings,' stages:'
   do,i=2,ntimings
      write(*,'(A,F10.2,A,F10.2,A)')"VERBOSE: "//timings(i)//':',cpu(i)-cpu(i-1),'sec',(cpu(i)-cpu(i-1))*100/cpu_tot,"%"
   end do

end if

if(verbose)write(*,'(A)')'VERBOSE: finished!'


end program NORM_tool

subroutine help(version)
implicit none
character(*)::version
character(4)::frmt
frmt='(A)'
write(*,frmt)'Program NORM_tool modifies and/or solves normal systems provided '
write(*,frmt)'in the format readable by read_BINtype'
write(*,frmt)'Version '//version
write(*,frmt)' The program may perform all or some of the actions below in the'
write(*,frmt)'  following order:'
write(*,frmt)' 1) Transformation of the systems'
write(*,frmt)' 2) Set apriori values of normal systems to those provided in a '
write(*,frmt)'    solution file or normal file'
write(*,frmt)' 3) Fix parameters to their a priori values and eliminate them '
write(*,frmt)'    from the normal systems'
write(*,frmt)' 4) Combine systems into new normal systems'
write(*,frmt)' 5) Reduce normal systems for specific parameters'
write(*,frmt)' 6) Add a regularization matrix to the normal systems'
write(*,frmt)' 7) Solve systems'
write(*,frmt)' 8) Estimate variance components of the root systems (VCE)'
write(*,frmt)''
write(*,frmt)'Usage NORM_tool [OPTIONS] [NORMFILE(S)]'
write(*,frmt)' Where NORMFILE(S) are one or more input files containing normal'
write(*,frmt)' systems with various parameters'
write(*,frmt)" Setting NORMFILE(S) to 'stdin' will read from standard input "
write(*,frmt)' [OPTIONS] may be one or more of the following operations and are'
write(*,frmt)' performed in the following order:'
write(*,frmt)''
write(*,frmt)' 1) Transformation of the systems.'
write(*,frmt)'   PARAMETER transformation'
write(*,frmt)'   The default transformation of system N,b,x0,ltpl by '
write(*,frmt)'   transformation B transforms the parameter space:'
write(*,frmt)'   N* = B N B^T '
write(*,frmt)'   b* = B b'
write(*,frmt)'   Solve for x0*: x0 = B^T x0* (only possible when B is '
write(*,frmt)'   invertible!!) '
write(*,frmt)'   ltpl*=ltpl (unchanged) '
write(*,frmt)'   Each parameter is transformed either through B or through'
write(*,frmt)'   a unit transformation.'
write(*,frmt)'   ERROR COVARIANCE transformation:'
write(*,frmt)'   Furthermore a covariance transformation may also be requested.'
write(*,frmt)'   This essentially scales the error covariance by B**-1'
write(*,frmt)'   in the new parameter space:'
write(*,frmt)'   N* = B N B^T : (no need to invert the original normal matrix)'
write(*,frmt)'   b*=N* N^-1b: Note that the normal system must be invertible!!'
write(*,frmt)'   x0*=x0'
write(*,frmt)'   ltpl*= (N^-1b)^T  N*  (N^-1b): NOTE Information of the previous'
write(*,frmt)'   ltpl (fit with original data) will be lost!'
write(*,frmt)'   OPTIONS:'
write(*,frmt)'  -Tf=TRANSFILE: TRANSFILE is the file name of the full (or '
write(*,frmt)'                 diagonal) transformation matrix B'
write(*,frmt)'  -Td: Diagonal transformation, TRANSFILE is assumed to be'
write(*,frmt)'       diagonal and an ascii file(2 columns : Parameter,value)'
write(*,frmt)'  -TC: Apply error covariance transformation (discussed above)'
write(*,frmt)'  -TA: ADD pararameter space to existing parameters'
write(*,frmt)'       The default removes parameters which have been transformed'
write(*,frmt)'  -Ti[=|:]INCREGEX: Only include rows of the transformation matrix'
write(*,frmt)'                which satisfy the POSIX regular expression'
write(*,frmt)"                INCREGEX (when -Ti=). In the case of -Ti: INCREGEX"
write(*,frmt)"                is assumed to be a file containing regular expressions"
write(*,frmt)"  -Te[=|:]EXCREGEX  Don't include rows of the transformation matrix"
write(*,frmt)"                with names which match EXCREGEX (case -Te=)"
write(*,frmt)"                case -Te:EXREGEX is assumed to be a file with reg. expr."
write(*,frmt)"   The -Ti and -Te options may not be used with symmetric "
write(*,frmt)"   transformations matrices"
write(*,frmt)'  -Tt: Transpose B -> B^T'
write(*,frmt)'  -Tb: Ignore right hand side b ( keep original values)'
write(*,frmt)'  -Tx: Ignore apriori vector x0 ( keep original values)'
write(*,frmt)'  -Txb: Treat x0 as b. Thus: x0*=B x0'
write(*,frmt)'  -Tn: Ignore transformation of the normal matrix. (keep the old'
write(*,frmt)'       part)'
write(*,frmt)'  The options -Tx -Txb -Tx and -Tn must be used with caution as'
write(*,frmt)'  following steps may be producing utter nonsense. '
write(*,frmt)'  They are merely there to facilitate propagation of solution'
write(*,frmt)'  systems and to avoid unneccessary computations.'
write(*,frmt)' '
write(*,frmt)' 2) Setting APRIORI values:'
write(*,frmt)'  -af=SOLFILE : set apriori values of the normal files to those'
write(*,frmt)'                solution vectors provided in SOLFILE.'
write(*,frmt)'                This will also change the value of the apriori'
write(*,frmt)'                residual norm LtPl and the right hand side of the'
write(*,frmt)'                system'
write(*,frmt)'                The difference between the current apriori values'
write(*,frmt)'                and the new values are influencing the new ltpl '
write(*,frmt)'                and right hand sides only.'
write(*,frmt)"                A special case is when SOLFILE='_ZERO_'"
write(*,frmt)"                In that case, the new apriori value will be zero"
write(*,frmt)"                (can be used with the restrictions -ai and -ae)"
write(*,frmt)'  -a0: Set the apriori values to zero without changing ltpl and'
write(*,frmt)'       the right hand side ( forget apriori information)'
write(*,frmt)'       This option may be used in conjunction with -af (the '
write(*,frmt)'       absolute new apriori values will be propagated)' 
write(*,frmt)'  -ar: Replace the apriori values with those as provided'
write(*,frmt)'       but do not change the LtPL or the right hand side vector'
write(*,frmt)"  -a+: Add the apriori vector of SOLFILE to the one already existing."
write(*,frmt)'  -ai[=|:]INCREGEX: Only use parameters from SOLFILE which match the'
write(*,frmt)'                INCREGEX regular expression (case -ai=)'
write(*,frmt)"                case -ai: Use regular expr from file INCREGEX"
write(*,frmt)"  -ae[=|:]EXCREGEX: Don't use parameters from SOLFILE which match the"
write(*,frmt)"                regular expression EXCREGEX"
write(*,frmt)"  -at: Solution file is a text file in ascii (1st column Parameter"
write(*,frmt)"       name (24 characters), 2nd column apriori value "
write(*,frmt)''
write(*,frmt)' 3) FIXING variables to their A PRIORI values:'
write(*,frmt)'  -fa: Auto Fix parameters with just updated apriori values (from '
write(*,frmt)'       the above -a option), the plain -fa option may only be'
write(*,frmt)'       used in conjunction with -af[ie=]'
write(*,frmt)'  -fi[=|:]INCREGEX : Only fix parameters which match the regular '
write(*,frmt)'                 expression INCREGEX'
write(*,frmt)"  -fe[=|:]EXCREGEX: DON'T fix parameters which match the regular "
write(*,frmt)"                expression EXCREGEX"
write(*,frmt)''
write(*,frmt)' 4) COMBINING identical normal systems'
write(*,frmt)'  -c: Auto Combine (add) all input systems which have exactly the'
write(*,frmt)'      same variables into a new normal system '
write(*,frmt)'  -cf: Force a combination into 1 final system'
write(*,frmt)'  -ca: Create a new a priori vector for the combination system, '
write(*,frmt)'       which is obtained by averaging the a priori vectors. The '
write(*,frmt)'       time span of the input system determines the averaging weight.'
write(*,frmt)'       This option automatically implies -cf'
write(*,frmt)''
write(*,frmt)' 5) REDUCING VARIABLES:'
write(*,frmt)'  -ra: Auto reduce local variables in the systems (cannot be '
write(*,frmt)'       combined with -cf)'
write(*,frmt)'       Thus variables which are only present in one of the input'
write(*,frmt)'       systems are reduced. If only one system is left after the'
write(*,frmt)'       combination step no action is taken.'
write(*,frmt)'  -ri[=|:]INCREGEX: Reduce all systems (after possible combination)'
write(*,frmt)'      for variables which match the regular expression INCREGEX'
write(*,frmt)"  -re[=|:]EXCREGEX: DON'T Reduce all systems (aftre possible"
write(*,frmt)"      combination) for variables which match the regular"
write(*,frmt)"      expression EXCREGEX"
write(*,frmt)'  NOTE: The reduction will only affect variables which have not'
write(*,frmt)'        been removed by transformation/fixing'
write(*,frmt)''
write(*,frmt)' 6) REGULARIZATION'
write(*,frmt)'  -Rf=REGFILE: Regularize the normal system(s) with the '
write(*,frmt)'               regularization matrix provided in REGFILE'
write(*,frmt)'  -Rd: Use a diagonal regularization'
write(*,frmt)'   If the options -Rf and -Rd are both provided the REGFILE must'
write(*,frmt)'   be an ascii file with rows consisting of 1)'
write(*,frmt)'   parameter names (24 characters) and 2) their diagonal'
write(*,frmt)'   regularization entries.'
write(*,frmt)'   If no REGFILE is provided the regularization matrix is assumed'
write(*,frmt)'   to be unit diagonal'
write(*,frmt)'  -Rs=SCALE: Scale the regularization matrix by REGSCALE before'
write(*,frmt)'   adding it.'
write(*,frmt)'  -Ri[=|:]INCREGEX: Only incorporate the regularization for variables'
write(*,frmt)'   which match the regular expression INCREGEX '
write(*,frmt)"  -Re[=|:]EXCREGEX: DON'T incorporate the regularization for variables"
write(*,frmt)"   which match the regular expression EXCREGEX"
write(*,frmt)'  NOTE: The reduction will only affect variables which have not'
write(*,frmt)'  been removed by transformation/fixing/reducing'
write(*,frmt)''
write(*,frmt)' 7) SOLVE system(s): perform a cholesky decomposition and solve'
write(*,frmt)'    the final system(s).'
write(*,frmt)'  -s: In combination with -c[f] only the amount of systems after'
write(*,frmt)'      the combination step are solved else all (modified) input'
write(*,frmt)'      systems will be solved.'
write(*,frmt)'  -st=N: solve the systems using the N most significant singular'
write(*,frmt)'         values in a truncated SVD approach'
write(*,frmt)' 8) VCE estimation:'
write(*,frmt)'  -V: Do one iteration of the Foerstner VCE estimation. The solve'
write(*,frmt)'      option must be on.'
write(*,frmt)'  -Vp: Prepare the next iteration. This writes the'
write(*,frmt)'       rescaled/adjusted original normal systems to files/pipes' 
write(*,frmt)''
write(*,frmt)' FURTHER options:'
write(*,frmt)'  -W=WEIGHT1/WEIGHT2/...  : Apply a weight to each of the input'
write(*,frmt)'     normal systems. This only affects the right hand side'
write(*,frmt)'  -C=COVSCALE1/COVSCALE2/...  : Apply a scale to each of the input'
write(*,frmt)'     normal systems. This affects both sides.'
write(*,frmt)'     Default uses the value of Sigma0^2 if provided in the normal'
write(*,frmt)'     system file.'
write(*,frmt)'  -F=BASENAME: Set the output basename, default is NORMtoolout'
write(*,frmt)'     (_01, _02, etc will be appended when more files are outputted)'
write(*,frmt)"   Setting BASENAME to 'stdout' will redirect ALL output to "
write(*,frmt)"   standard output "
write(*,frmt)'  -S=SRCHSTRT/SRCHND: Compare only the substring in the parameter'
write(*,frmt)'     names ranging from SRCHSTR untill SRCHND'
write(*,frmt)'     (1<= SRCHSTR <= SRCHND <= 24)'
write(*,frmt)'  -v: print additional information to screen during process'
write(*,frmt)'  -w: ignore all warning messages'
write(*,frmt)'  -D=DESCRIPTION: put a different description in the header part '
write(*,frmt)'     of the output files. The description can have a maximum'
write(*,frmt)'     length of 80 characters.'
write(*,frmt)''   
write(*,frmt)' LEARNING options:'
write(*,frmt)'  Reads or writes a learning curve file. This options may use to'
write(*,frmt)'  skip the sorting part of the program, which might be interesting'
write(*,frmt)'  when doing multiple calls with the same setup or to check'
write(*,frmt)'  whether the program will do what you requested.'
write(*,frmt)'  NOTE: Options give after the -L option will overule learning '
write(*,frmt)'        curve settings.'
write(*,frmt)'  NOTE: Options influencing the sorting may not be redefined and'
write(*,frmt)'        will be ignored!!'
write(*,frmt)'  -Lf=CONFFILE: reads a learning curve from a file CONFFILE'
write(*,frmt)'  -Lw=CONFFILE: Create a Learning curve file and write this to'
write(*,frmt)'                CONFFILE.' 
write(*,frmt)'  -Lh: Halt after writing the learning curve file' 
write(*,frmt)'PIPING: The normal systems may be read from and written to linux'
write(*,frmt)'        fifo pipes and standard output ,enabling fast linking of'
write(*,frmt)'        several calls to NORM_tool.'
stop
end subroutine help

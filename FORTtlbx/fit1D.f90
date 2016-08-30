!!Command line program to fit 1 dimensional (time)series with several types of (non)-linear models
!!reads from standard input and from a file a 2 or 3 three column series
!writes fitted coefficients, reconstructed series ( fromt estimated coefficients) and residuals to standard output.

!todo (05-01-08)
!basically everything:

!non-linear fits
!discontinuities
!exponential fits

!! uses BLAS routines
!!Updated by Roelof Rietbroek, Thu Sep 24 16:07:39 2009
!!removed the length limitation of the input series ( data is now dynamically reallocated)
!!added parameter for errorscaling with posteriori

!!Updated by Roelof Rietbroek, Mon Mar  1 10:33:46 2010
!!fixed reading reference time from command line

!!Updated by Roelof Rietbroek, Mon Aug 23 09:16:14 2010
!! explicitly open standard output with longer record length (avoid unwanted line breaks)

!!Updated by Roelof Rietbroek, Fri Aug 27 11:17:29 2010
!!allow a general polynomial fit
!!make parameters more general in design 

!!Updated by Roelof Rietbroek, Sun Sep 12 16:41:29 2010
!! also show (weighted) RMS of series before and after

!!Updated by Roelof Rietbroek, Sun Sep 19 21:11:25 2010
!! added null hypothesis testing print out Pvalues 

!!Updated by Roelof Rietbroek, Tue Mar  8 10:02:38 2011
!! commented out explicit open statement on standard output (yields a runtime error on gfortran)

!!Updated by Roelof Rietbroek, Fri Aug 31 15:04:43 2012
!!also added the fitting of a power law (using log10)

!!Updated by Roelof Rietbroek, Wed Mar 20 09:30:28 2013
!! also allowed the output of the (rescaled) errors in the data

!!Updated by Roelof Rietbroek, Mon Oct 27 12:48:25 2014
!!add an option to automatically remove intercept
!! also allow a minimum polynomial degree 
!default center time centers the data in the middle

!!Updated by Roelof Rietbroek, Wed Mar 18 12:56:30 2015
!!also allow the use of an empirical autocovariance function of order n

!!Updated by Roelof Rietbroek, Thu Jul 16 11:31:08 2015
!! fixed bug (nullify para_des and para_typ pointers) which caused segfault with gfortran 5

!!Updated by Roelof Rietbroek, Tue Aug 18 17:32:32 2015
!!added option to keep solution vector from the first iteration while modifying the error with an AR model


program fit1D
use forttlbx
use autoregressive
implicit none
character(200)::dum,file,parastring,dum2
character(20)::version,date,time
integer::stderr,itharg,narg,unit,nmax,i,j
logical::stdin
double precision::Ohm,t0,ltpl,ltplold,ddot,pi,am,ph,Ohmtmp
double precision,pointer::epoch(:)=>null()
double precision,pointer::input(:)=>null()
double precision,pointer::sigma(:)=>null()
double precision,pointer::cpepoch(:),cpinput(:),cpb(:)
double precision::covnew(2,2),prop(2,2),work(3)
double precision,allocatable::A(:,:),cpA(:,:),b(:),N(:,:),prob(:)
double precision,allocatable::residuals(:),bdcov(:,:)
double precision,allocatable::Solbackup(:)
double precision::sigma0,intercept,gam
integer::npara,ndat,stdout,nshift,ncol,last
logical::weight,phase,header,colout,coloutfull,colouterr
logical::errscale,pvalue
logical::t0provided
integer::iargc,chunk,nalloc
integer::len_sysdep
character(10),pointer::para_des(:)=> null() !pointer holding parameter descriptions
integer,pointer::para_typ(:)=>null() ! pointer holding type description
integer::poly_ord,ind,min_poly_ord
logical::removeintercept,autoregr
integer::npass,info,k,np
double precision::thresratio
type(ArModel):: Ar
double precision::rmsapri,rmspost
logical::arErroronly
! 1:mean
! 2:trend
! 3:annual harmonic (sine and cosine)
! 4:semiannual harmonic (sine and cosine)
! 5:general polynomial



!defaults
chunk=1000
pi=acos(-1.d0)
Ohm=2*pi
t0=0.d0
version='Experimental3'
stdin=.true. ! assume input file is read from standard input
!unit for output to standard error
stderr=0
stdout=6
unit=5

weight=.false.
header=.true.
phase=.false.
npara=0
ncol=2
colout=.false.
coloutfull=.false.
parastring=''
sigma0=1.d0
errscale=.false.
poly_ord=2
min_poly_ord=0
pvalue=.false.
colouterr=.false.
removeintercept=.false.
t0provided=.false.
autoregr=.false.
npass=1
Ar%ord=1
arErroronly=.false.
!!process command line options
narg=iargc()

! write(0,*)pval_stud_t(2.d0,60)
! stop
if (narg < 1)call help(version) !no arguments make no sense

itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:4))
      case('fs')! fit semi annual
         call realloc_ptr(para_des,2)
         call realloc_ptr(para_typ,2)

         para_des(npara+1)='SA_COS'
         para_des(npara+2)='SA_SIN'
         
         para_typ(npara+1)=4
         para_typ(npara+2)=-4
         npara=npara+2
         
      case('fa','fa=')!fit annual term
         if(dum(4:4) .eq. '=')read(dum(5:),*)Ohm ! possibly read different annual frequency

         call realloc_ptr(para_des,2)
         call realloc_ptr(para_typ,2)

         para_des(npara+1)='AN_COS'
         para_des(npara+2)='AN_SIN'
         
         para_typ(npara+1)=3
         para_typ(npara+2)=-3
         npara=npara+2

      case('fm')! fit a mean 
         call realloc_ptr(para_des,1)
         call realloc_ptr(para_typ,1)

         para_des(npara+1)='MEAN'
         para_typ(npara+1)=1
         npara=npara+1
      case('ft','ft=') !fit a trend
         if(dum(4:4) .eq. '=')then
            read(dum(5:),*)t0
            t0provided=.true.
         end if

         call realloc_ptr(para_des,1)
         call realloc_ptr(para_typ,1)

         para_des(npara+1)='TREND'
         para_typ(npara+1)=2
         npara=npara+1
      case('fp','fp=') !fit a polynomial
         if(dum(4:4) .eq. '=')then
            ind=index(dum,',')
            if(ind .eq. 0)then
               dum2=dum(5:)
             !  read(dum(5:),*)poly_ord
            else
               dum2=dum(5:ind-1)
            !   read(dum(5:ind-1),*)poly_ord
               read(dum(ind+1:),*)t0
               t0provided=.true.
            end if
            !also check for a minimum polynomial order
            ind=index(dum2,'-')
            if(ind .eq. 0)then
               read(dum2,*)poly_ord
            else
               read(dum2(1:ind-1),*)min_poly_ord
               read(dum2(ind+1:),*)poly_ord
            end if
            
         end if

         call realloc_ptr(para_des,poly_ord-min_poly_ord+1)
         call realloc_ptr(para_typ,poly_ord-min_poly_ord+1)
         do,j=min_poly_ord,poly_ord
            write(para_des(npara+1),'(a5,i1)')'POLY_',j
            para_typ(npara+1)=5
            npara=npara+1
         end do

!         npara=npara+poly_ord-min_poly_ord+1
         ! write(*,*)npara,min_poly_ord,poly_ord,t0
         ! do,j=1,npara
         !    write(*,*)para_des(j),para_typ(j)
         ! end do
         ! stop
      case('fP') !fit power law a x^p
         call realloc_ptr(para_des,2)
         call realloc_ptr(para_typ,2)
         para_des(npara+1)='PSCALE'
         para_typ(npara+1)=6
         para_typ(npara+2)=6
         para_des(npara+2)='POWER'
         npara=npara+2 ! scale plus exponent
      case('s')!scale errors
         errscale=.true.
      case('n')! don't print header
         header=.false.
      case('o')! output columns with t,dat,fit
         colout=.true.
      case('os')! output columns with t,dat,fit
         colout=.true.
         coloutfull=.true. ! also output separate fits fro each estimated parameter
      case('ose')!output separate fit and scaled errors
         colout=.true.
         coloutfull=.true.
         colouterr=.true.
      case('oe')!output fit and scaled errors
         colout=.true.
         colouterr=.true.
      case('p')!express harnmonics in phases
         phase=.true.
      case('t') !retrieve probability that the true parameter values are larger than the estimated parameters given the null hypothesis from student t-distribution
         pvalue=.true.
      case('w')! weighted least squares
         weight=.true.
         ncol=3
      case('i')
         removeintercept=.true.
      case('ar','ar=')
         autoregr=.true.
         npass=2
         if(dum(4:4) .eq. '=')then
            read(dum(5:),*)Ar%ord
         end if

      case('are')
         autoregr=.true.
         npass=2
         if(dum(5:5) .eq. '=')then
            read(dum(6:),*)Ar%ord
         end if
         arErroronly=.true.
      case('h') ! ask for help

         call help(version)
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:3)
         stop
      end select
   else!argument is not an option ( but a file name)
      unit=14 ! change unit number
      stdin=.false.
      file=dum
   end if
end do



!some input checks
if( npara .eq. 0 )then
   write(stderr,*)"ERROR: fit1D, no fit selected"
   stop
end if

if(phase .and. colout)then
   write(stderr,*)"ERROR: -p and -o option are not allowed simultaneously"
   stop
end if

if(phase .and. pvalue)then
   write(stderr,*)"WARNING: t-test from phase parameters are based on the wrong pdf"
end if


if(npara .ne. 2 .and. minval(para_typ).eq. 6 & 
.and. maxval(para_typ) .eq. 6)then
   write(stderr,*)"ERROR: -fP may  not be combined with other -f options"
   stop
end if

if(para_typ(1) .eq. 6 .and. coloutfull  )then
   write(stderr,*)"ERROR: -fP may not be combined with the -os option"
   stop
end if

if(para_typ(1) .eq. 6 .and. autoregr  )then
   write(stderr,*)"ERROR: -fP may not be combined with the -a option"
   stop
end if


if(para_typ(1) .eq. 6 .and. weight  )then
   write(stderr,*)"ERROR: -fP may not be combined with the -w option"
   stop
end if

if(weight .and. autoregr)then
   write(stderr,*)"WARNING: supplied sigmas will be replaced by empirical autoregrariances from an AR model on the second pass"
end if

!read input from file or standard input

if(.not. stdin)open(unit=unit,file=trim(file),form='formatted')
last=0
ndat=1
!allocate pointers
allocate(epoch(chunk),input(chunk))
if(weight .or. autoregr)allocate(sigma(chunk))
nalloc=chunk
!try reading first line to find out whether data has enough columns
read(unit=unit,fmt=*,iostat=last)work(1:ncol)

if( last < 0)then
   write(stderr,*)"ERROR: not enough columns"
   stop
else if(last > 0)then
   write(stderr,*)"ERROR: format error"
   stop
end if
epoch(ndat)=work(1)
input(ndat)=work(2)
if(weight)sigma(ndat)=work(3)

do while (last .eq. 0)
   ndat=ndat+1
   if(ndat > nalloc)then !reallocate arrays ( add a chunk) (& copy data)
      call realloc_ptr(epoch,chunk)
      call realloc_ptr(input,chunk)
      if(weight)call realloc_ptr(sigma,chunk)
      nalloc=nalloc+chunk

   end if

   read(unit=unit,fmt=*,iostat=last)work(1:ncol)
   if( last .ne. 0)then
      ndat=ndat-1
      exit ! exit on end of file
   end if
   !put data in arrays
   epoch(ndat)=work(1)
   input(ndat)=work(2)
   if(weight)sigma(ndat)=work(3)
end do

if(.not. stdin)close(unit)

!end data reading

!check for impossible autocovariance requests
if(autoregr .and. Ar%ord > ndat)then
   write(stderr,*)"ERROR: maximum lag of autocovariance exceeds the amount of time points"
   stop
end if

!remove intercept if requested
if(removeintercept)then
   intercept=sum(input(1:ndat))/ndat
   do,i=1,ndat
      input(i)=input(i)-intercept
   end do
end if


!compute the apriori root mean square
rmsapri=0.d0
do,i=1,ndat
   rmsapri=rmsapri+(input(i)**2)/ndat
end do
rmsapri=sqrt(rmsapri)

!use the center time of the input as the reference
if(.not. t0provided)then
   t0=(minval(epoch(1:ndat))+maxval(epoch(1:ndat)))/2.d0
end if


!copy input data
allocate(cpepoch(ndat),cpinput(ndat),cpb(npara))
cpinput=input(1:ndat)
cpepoch=epoch(1:ndat)

!allocate memory for the design matrix and normal matrix
allocate(A(ndat,npara),N(npara,npara),b(npara))



!!!!!!
!construct design matrix A
nshift=0
do while (nshift < npara)
   select case(para_typ(nshift+1))
   case(1)! mean
      do,i=1,ndat
         A(i,nshift+1)=1.d0
      end do
      nshift=nshift+1

   case(2)!trend
      do,i=1,ndat
         A(i,nshift+1)=(epoch(i)-t0)
      end do
      nshift=nshift+1

   case(3,4)!annual or semiannual cosine and sine
      if(para_typ(nshift+1) .eq. 3)then
         ohmtmp=Ohm
      else
         ohmtmp=2*Ohm
      end if

      do,i=1,ndat
         A(i,nshift+1)=cos(ohmtmp*epoch(i))
         A(i,nshift+2)=sin(ohmtmp*epoch(i))
      end do
      nshift=nshift+2
   case(5) ! polynomial fit
      do,j=min_poly_ord,poly_ord
         do,i=1,ndat
            A(i,nshift+1)=(epoch(i)-t0)**(j)
         end do
         nshift=nshift+1
      end do
!      nshift=nshift+poly_ord-min_poly_ord+1
   case(6) ! fit power law using log10
      ! copy the log10 of the original data 
     


      !modify epoch data and sigma
      
      do,i=1,ndat
         input(i)=LOG10(input(i))
         epoch(i)=LOG10(epoch(i))
         A(i,nshift+1)=1.d0
         A(i,nshift+2)=epoch(i)
      end do

      nshift=2
      exit ! exit loop since only 2 parameters are allowed
   end select

end do

!copy design matrix
allocate(cpA(ndat,npara))
cpA=A

allocate(residuals(ndat))

!!!!!!construct normal matrix N

do, np=1,npass !possibly two passes needed

!in the case of a weighted least squares: temporary scale the observations and design matrix wioth the standard deviations
   if(weight)then

      do,i=1,npara
         do,j=1,ndat
            !scale rows of the design matrix 
            A(j,i)=cpA(j,i)/sigma(j)
         end do
      end do
      !scale observation vector
      do,j=1,ndat
         input(j)=input(j)/sigma(j)
      end do
end if


!normal matrix N
call dsyrk('U','T',npara,ndat,1.d0,A,ndat,0.d0,N,npara)

!right hand side vector (b)
call dgemv('T',ndat,npara,1.d0,A,ndat,input,1,0.d0,b,1)

!apriori residual fit
ltpl=ddot(ndat,input,1,input,1)
ltplold=ltpl
!!solve normal system by means of a Cholesky decomposition
call solve_norm(N,b,ltpl)

!possibly create a backup of the solution vector
if(np==1 .and. arErroronly)then
   allocate(Solbackup(npara))
   Solbackup=b(1:npara)
end if

!!compute the post-fit residuals
residuals=cpinput(1:ndat)
!subtract estiamted fit from inoput o create residuals(residuals is updated)
call dgemv('N',ndat,npara,-1.d0,cpA(1,1),ndat,b(1),1,1.d0,residuals(1),1)

!compute posteriori rms
!compute the apriori root mean square
rmspost=0.d0
do,i=1,ndat
   rmspost=rmspost+(residuals(i)**2)/ndat
end do
rmspost=sqrt(rmspost)

if(autoregr .and. np .eq. 1)then
   weight=.false. !disable weighting by sigma on next pass

   !compute AR model from residuals
   call ComputeArYW(residuals,Ar)

   !compute banded autocovariance matrix from AR model
   call makeautoCovBD(Ar)

   !decorrelate right hand side and design matrix

   !perform a check on the strength of the first offdiagional. Possibly increase the strength of the diagonal
   thresratio=2.1
   if(abs(Ar%autocovBD(1,1)/Ar%autocovBD(2,1)) <= thresratio)then
      write(stderr,*)"WARNING: Empirical autocovariance may result in an unstable positive semi-definite matrix"
      write(stderr,*)"         Artificially adding noise on the diagonal of the covariance matrix.."
      Ar%autocovBD(1,:)=abs(Ar%autocovBD(2,:))*thresratio; !make the diagonal at least 210% larger then the off diagonal
   end if


   !compute a cholesky decomposition of the banded matrix
   call dpbtrf('L',ndat,Ar%ord,Ar%autocovBD(1,1),size(Ar%autocovBD,1),info)
   if(info .ne. 0)then
      write(stderr,*)"ERROR: Cholesky decomposition failed on empirical banded covariance matrix",info
   
      stop
   end if

   !dpbtrs (UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO)
   !decorrelate data vector
   input(1:ndat)=cpinput(1:ndat)
   call dpbtrs ('L', ndat,Ar%ord, 1, Ar%autocovBD(1,1),size(Ar%autocovBD,1), input(1),ndat,info)

   if(info .ne. 0)then
      write(stderr,*)"ERROR: decorrelation failed",info
      stop
   end if
   
   !decorrelate design matrix
   A=cpA;
   call dpbtrs ('L', ndat,Ar%ord, npara, Ar%autocovBD(1,1),size(Ar%autocovBD,1), A(1,1),ndat,info)

   if(info .ne. 0)then
      write(stderr,*)"ERROR: decorrelation failed",info
      stop
   end if
   
   
end if
   
end do !end passes


!possibly replace solution vector with initial estimates
if(ArErroronly)then
   b(1:npara)=Solbackup(1:npara)
end if

!convert harmonic to magnitude and phase if desired
if(phase)then
   nshift=0
   do while (nshift < npara)
      select case(para_typ(nshift+1))
      case(3,4) ! annual or semiannual
         if(para_typ(nshift+1) .eq. 3)then
            ohmtmp=Ohm
         else
            ohmtmp=2*Ohm
         end if
         
         am=sqrt(b(nshift+1)**2+b(nshift+2)**2)
         ph=atan2(b(nshift+2),b(nshift+1))/ohmtmp;
         !error propagation

         !make propagator matrix
         prop(1,1)=b(nshift+1)/am
         prop(1,2)=b(nshift+2)/am

         prop(2,1)=b(nshift+2)/(b(nshift+1)**2*(1+(b(nshift+2)/b(nshift+1))**2)*ohmtmp)
         prop(2,2)=1/(b(nshift+1)*(1+(b(nshift+2)/b(nshift+1))**2)*ohmtmp)

         covnew=N(nshift+1:nshift+2,nshift+1:nshift+2)
         covnew(2,1)=covnew(1,2) ! mirror

         covnew=matmul(prop,matmul(covnew,transpose(prop)))
         
         para_des(nshift+1)(4:7)='AMPL'
         para_des(nshift+2)(4:8)='PHASE'
         b(nshift+1)=am
         b(nshift+2)=ph
         N(nshift+1,nshift+1)=covnew(1,1)
         N(nshift+2,nshift+2)=covnew(2,2)

         nshift=nshift+2
      case default ! just shift one item
         nshift=nshift+1
      end select
         
   end do
end if

!!!OUTPUT!!!!!!!!!
!!explicitly open standard output with long enough record length ( needed for some compilers but not for gfortran)
len_sysdep=(npara+2)*25
!open(unit=stdout,form='formatted',recl=len_sysdep)

if(para_typ(1) .eq. 6)then !power law
   cpb=b ! copy values
   b(1)=10**b(1) ! reverse the log10 action of the intersect -> scale
end if


sigma0=sqrt(ltpl/(ndat-npara)) !compute posteriori sigma
!print information on the estimated parameters in the header


if(header)then
   call date_and_time(date,time)
   write(stdout,fmt=*)'!!!output from fit1D ',&
        date(7:8)//"-"//date(5:6)//"-"//date(1:4),' ',time(1:2)//":"//time(3:4),'!!!'   
   write(stdout,fmt=*)'NOBS',ndat
   write(stdout,fmt=*)'NPARA',npara
   write(stdout,fmt=*)'LTPL_APRI',ltplold
   write(stdout,fmt=*)'LTPL_POST',ltpl

   write(stdout,fmt=*)'SIGMA0',sigma0
   write(stdout,fmt=*)'RMS_APRI',rmsapri
   write(stdout,fmt=*)'RMS_POST',rmspost


   !report fitted AR model parameters if fitted
   if(autoregr)then
      do,i=1,Ar%ord
         write(stdout,fmt='(1x,A7,I1,A1,1x,G12.6)')'AR_MOD(',i,')',Ar%para(i)
      end do
         write(stdout,fmt=*)'AR_SIG^2',Ar%sig2
   end if
   

   if(pvalue)then
      allocate(prob(npara))
      do,i=1,npara
         prob(i)=pval_stud_t(b(i)/(sigma0*sqrt(N(i,i))),ndat-npara)
      end do
   end if

   !write estimated parameters to file
   if(.not. errscale)sigma0=1 !reset to 1 if rescaling is not desired
   if(pvalue)then
      do,i=1,npara
         write(stdout,fmt=*)para_des(i),b(i),sigma0*sqrt(N(i,i)),prob(i)
      end do
   else
      do,i=1,npara
         write(stdout,fmt=*)para_des(i),b(i),sigma0*sqrt(N(i,i))
      end do
   end if
   !output string with order of the parameters
   if(colout)write(unit=stdout,fmt='(A)',ADVANCE='NO')'!! TIME DAT'
   if(coloutfull)then
      write(unit=stdout,fmt=*)para_des(1:npara)
   else if(colout)then
      write(unit=stdout,fmt=*)' FIT'
   end if

end if





if( .not. colout)stop ! quit when no more output is requested


if(para_typ(1) .eq. 6)then !power law
   b(1)=cpb(1) ! copy fitted log10 value back in solution
end if

!!adapt design matrix ( multiply each column with the estimated parameter)
do,i=1,npara
   do,j=1,ndat
      A(j,i)=cpA(j,i)*b(i)
   end do
end do

!! write column data 


if(para_typ(1) .eq. 6)then ! select orginal data
   !deallocate modified data and point to original data
   deallocate(input)
   deallocate(epoch)
   input=>cpinput
   epoch =>cpepoch
   
   do,i=1,ndat
      A(i,1)=10**(sum(cpA(i,:)))
      A(i,2)=0.d0
   end do
   
end if

!possibly allocate sigma vector
if(.not. (weight .or.autoregr) .and. colouterr)then
   allocate(sigma(ndat))
   sigma=1.d0 !vector needs to exist set to 1
end if

!possibly rescale errors
if(errscale .and. colouterr)then
!   write(0,*)'rescaling error with: ', sigma0
   sigma=sigma0*sigma
end if

if(coloutfull)then
   if(colouterr)then
      do,i=1,ndat
         write(stdout,fmt=*)cpepoch(i),cpinput(i),A(i,:),sigma(i)
      end do
   else
      do,i=1,ndat
         write(stdout,fmt=*)cpepoch(i),cpinput(i),A(i,:)
      end do
   end if
else
   if(colouterr)then
      do,i=1,ndat
         write(stdout,fmt=*)cpepoch(i),cpinput(i),sum(A(i,:)),sigma(i)
      end do
   else
      do,i=1,ndat
         write(stdout,fmt=*)cpepoch(i),cpinput(i),sum(A(i,:))
      end do
   end if
end if




end program fit1D

subroutine help(version)
character(*):: version
integer::unit=0
character(4)::frmt
frmt='(A)'

write(unit,frmt)"Program fit1D version"//trim(version)
write(unit,frmt)"Fits curves to 1D (time) series"
write(unit,frmt)"usage: fit1D [OPTIONS] [FILE]"
write(unit,frmt)"Reads from a file or from standard input"
write(unit,frmt)"input has 2 or 3 column2 ( time, value, standard dev) "
write(unit,frmt)"denoted T,DAT,DEV below"
write(unit,frmt)"writes to standard output"
write(unit,frmt)""
write(unit,frmt)"OPTIONS may be:"
write(unit,frmt)"fits"
write(unit,frmt)"-fm: fit a mean"
write(unit,frmt)"-ft[=T0]: fit a trend [and optionally adapt reference time]"
write(unit,frmt)"-fa[=OMEGA]: fit an annual harmonic [and optionally adapt annual period]"
write(unit,frmt)"-fs: fit a semiannual harmonic"
write(unit,frmt)"-fp[=[MINORD-]ORD,T0]: fit a polynomial of order ORD (default=2) through the data"
write(unit,frmt)"   optionally a minimum order MINORD can be provided to restrict the minimum fitted polynomial" 
write(unit,frmt)"        Note: cannot be combined with -ft and -fm since it causes a rank defect"
write(unit,frmt)"-fP: fit a (Non-linear) power law d = a * t^p, (unknows, a and p)"
write(unit,frmt)"     The data will be transformed with log10 to enable a linear fit"
write(unit,frmt)"     Of course this may violate the Gaussiannity of the input data"
write(unit,frmt)"Other fits or not yet supported but planned ( also non-linear)"
write(unit,frmt)" The observation equation is composed out of 1 or more of the following terms"
write(unit,frmt)" DAT = MEAN"
write(unit,frmt)"     + TREND(T-T0) "
write(unit,frmt)"     + sum(i=0,ORD) { POLY_i (T-T0)^i }"
write(unit,frmt)"     + A_an cos(OMEGA T) + B_an sin (OMEGA T)"
write(unit,frmt)"     + A_sa cos(2 OMEGA T)+B_sa sin (2 OMEGA T)"
write(unit,frmt)"Defaults: OMEGA = 2 PI, the reference time T0 defaults to the center time of the input series"
write(unit,frmt)""
write(unit,frmt)"Output options:"
write(unit,frmt)"-n: Don't print header information (eg. estimated parameter values)"
write(unit,frmt)"-s: scale errors by the posteriori sigma0"
write(unit,frmt)"-o[s][e]: Print data in columns ( t, dat, [ total fit | fit per parameter])"
write(unit,frmt)"  The optional parameter s causes the total fit to be separated per parameter"
write(unit,frmt)"  Columns obey the same sequence as the order of the -f options on the command line"
write(unit,frmt)"  append 'e' to also print the (rescaled) errors in the last column"
write(unit,frmt)"-p: Output harmonic information in terms of Amplitude and Phase"
write(unit,frmt)"-w: use weighted least squares ( third column with standard deviations"
write(unit,frmt)"    required)."
write(unit,frmt)"-ar[=n]: Estimate an Autoregressive model of order n  (default=1), for the residuals, and use"
write(unit,frmt)"         this as an empirical error-covariance"
write(unit,frmt)"         Note: this will replace the -w option in the second run"
write(unit,frmt)"-ear[=n]: only use the AR model to correct the error estimate, but don't reestimate the fit "
write(unit,frmt)"-t: print out pvalues (in third column of estimated parameters) based on the "
write(unit,frmt)"     appropriate student-t distribution."
write(unit,frmt)"    the p-values represent the probability to obtain a test statistic that is more extreme"
write(unit,frmt)"    than observed given the null hypothesis (the value is 0)."
write(unit,frmt)"-i: Remove the intercept from the data before fitting"
write(unit,frmt)"-h This help message"
write(unit,frmt)""
stop


end subroutine help

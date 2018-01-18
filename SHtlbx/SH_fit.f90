!! Program SH_fit fits functionals through spherical harmonic coefficients

!!Coded by Roelof Rietbroek, Thu Mar 14 14:03:12 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!!Updated by Roelof Rietbroek, Fri Apr 19 16:20:53 2013
!!also provide the possibility to write the SH residuals



program SH_fit
use SHtlbx
use FORTtlbx
use FORTtime
implicit none
integer::narg,i,j,lmax,lmin,lmaxf,gtyp,ind,pos,itharg,l,m,stderr,nf,posm
character(200)::dum
character(200),pointer::filen(:)=>null()
logical::limdeg,weight,sigscale,rmsoutput,sigscaleresi
double precision::mean,tcent,tstart,tend,t0,ltpl,Ohm,ohmtmp,ddot
double precision,allocatable, dimension(:)::time,rhs,input,clmrms
double precision,allocatable,dimension(:)::stime,etime,sigma0
double precision,allocatable,dimension(:,:)::clm,clm_sig,clm_f,clm_sigf
double precision,allocatable, dimension(:,:)::A,Atmp,N,Nreg
integer,allocatable::perm(:)
integer::nmonths
parameter(nmonths=12)
integer::iargc,chunk,npara,nshift
character(10),pointer::para_des(:)=>null() !pointer holding parameter descriptions
integer,pointer::para_typ(:)=>null() ! pointer holding type description
integer::poly_ord,otyp,sts,nds,stc,ndc
character(200)::basen,fileout
logical::resSH
logical::usereg !uses an additional regularization
logical::withclim
double precision:: regscale
type(time_t)::epoch
integer:: fmonth
!get command line arguments
narg=iargc()
lmin=0
limdeg=.false.
weight=.false.
stderr=0
chunk=50
npara=0
basen='SHFIT_'
sigscale=.false.
rmsoutput=.false.
resSH=.false.
otyp=4
Ohm=2*pi
t0=2003.0
sigscaleresi=.false.
allocate(filen(chunk))
usereg=.false. 
withclim=.false.
itharg=0
nf=0
lmax=-10!gives an error if not redefined


if(narg < 1)call help()!call help straight away


do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .ne. '-')then
      nf=nf+1
      if(size(filen,1) < nf)then
         call realloc_ptr(filen,chunk)
      end if
      filen(nf)=dum
   else !if argument is an option
      select case(dum(2:3))
      case('fs')! fit semi annual
         call realloc_ptr(para_des,2)
         call realloc_ptr(para_typ,2)
         
         para_des(npara+1)='SA_COS'
         para_des(npara+2)='SA_SIN'
         
         para_typ(npara+1)=4
         para_typ(npara+2)=-4
         npara=npara+2
         
      case('fa')!fit annual term
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
      case('ft') !fit a trend
         if(dum(4:4) .eq. '=')then
            read(dum(5:),*)t0
         end if

         call realloc_ptr(para_des,1)
         call realloc_ptr(para_typ,1)

         para_des(npara+1)='TREND'
         para_typ(npara+1)=2
         npara=npara+1
      case('fp') !fit a polynomial
         if(dum(4:4) .eq. '=')then
            ind=index(dum,',')
            if(ind .eq. 0)then
               read(dum(5:),*)poly_ord
            else
               read(dum(5:ind-1),*)poly_ord
               read(dum(ind+1:),*)t0
            end if
         end if

         call realloc_ptr(para_des,poly_ord+1)
         call realloc_ptr(para_typ,poly_ord+1)
         do,j=0,poly_ord
            write(para_des(npara+j+1),'(a5,i1)')'POLY_',j
            para_typ(npara+j+1)=5
         end do

         npara=npara+poly_ord+1
      case('fc')!monthly climatology
         withclim=.true.         
         call realloc_ptr(para_des,nmonths)
         call realloc_ptr(para_typ,nmonths)

         para_des(npara+1)='JAN'
         para_des(npara+2)='FEB'
         para_des(npara+3)='MAR'
         para_des(npara+4)='APR'
         para_des(npara+5)='MAY'
         para_des(npara+6)='JUN'
         para_des(npara+7)='JUL'
         para_des(npara+8)='AUG'
         para_des(npara+9)='SEP'
         para_des(npara+10)='OCT'
         para_des(npara+11)='NOV'
         para_des(npara+12)='DEC'
         
         para_typ(npara+1:npara+nmonths)=6
         
         npara=npara+nmonths
      
      case('l=')!limit maximum (and possibly minum degree)
         limdeg=.true.
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
      case('w')!weighting
         weight=.true.
      case('s')
         sigscale=.true.
      case('r')!output the variance after the fit (or if no fit)
         rmsoutput=.true.
      case('Rs')!also write the residuals to seperate files (scaled with posteriori sigma0)
         sigscaleresi=.true.
         resSH=.true.
         sigscale=.true.
      case('R')!also write the residuals to seperate files
         resSH=.true.

      case('F=')! use a different basename
         basen=trim(dum(4:))
      case default
         write(stderr,*)'unknown option selected, quitting'
         call help()
      end select
   end if
end do

if(nf<1)then
   write(stderr,*)'No files specified, quitting'
   stop
end if


!possibly when a polynomial is estimates with monthly climatology we need to
!apply an additional regularization
if(withclim .and. npara .ne. 12)then
    usereg=.true.
    allocate(Nreg(npara,npara))
    Nreg=0.d0
    !we're going to create a regularization matrix, Nreg, which constrains the sum of
    !all the monthly values to be 0
    do,i=1,npara
        do, j=1,i
            if(para_typ(i)==6 .and. para_typ(j)==6)then
                Nreg(i,j)=1
                Nreg(j,i)=1
            end if
        end do
    end do

    !do,i=1,npara
       !write(0,'(14F2.0)')Nreg(i,:)
    !end do

end if




!read in input data
do,i=1,nf
   write(stderr,*)"reading: ",trim(filen(i))
   
!now get time tags and maximum degree supported from main files
   call SH_readmeta(filen=trim(filen(i)),lmax=lmaxf,tcent=tcent,tend=tend,tstart=tstart,type=gtyp)
   
   if(i==1)then !! allocate and determine lmax if not given explixitly
      if(lmax <0)lmax=lmaxf
      pos=SH_pos(lmax,lmax)
      allocate(time(nf),stime(nf),etime(nf),clm(2*pos,nf))
      
      time=0.d0
      clm=0.d0
      if(weight)then
         allocate(clm_sig(2*pos,nf))
         clm_sig=0.d0
      end if
      stc=1
      ndc=pos
      sts=pos+1
      nds=2*pos
   end if

   if(lmaxf < lmax )then
      write(stderr,'(A)')"WARNING: requested maximum degree",lmax,"is not supported by file:",trim(filen(i))
      
   end if

   !read coefficients
   if(weight)then
      call SH_readgrav(filen=trim(filen(i)),clm=clm(stc:ndc,i),&
          &slm=clm(sts:nds,i),clm_sig=clm_sig(stc:ndc,i),slm_sig=clm_sig(sts:nds,i),type=gtyp)
   else
      call SH_readgrav(filen=trim(filen(i)),clm=clm(stc:ndc,i),slm=clm(sts:nds,i),type=gtyp)
   end if
   time(i)=tcent
   stime(i)=tstart
   etime(i)=tend
end do

posm=SH_pos(lmin-1,lmin-1)
!perform the fit
if(npara >0)then
   allocate(A(nf,npara))
   ! make design matrix  (unweighted!)   
   A=0.d0

   nshift=0
   do while (nshift < npara)
      select case(para_typ(nshift+1))
      case(1)! mean
         do,i=1,nf
            A(i,nshift+1)=1.d0
         end do
         nshift=nshift+1
         
      case(2)!trend
         do,i=1,nf
            A(i,nshift+1)=(time(i)-t0)
         end do
         nshift=nshift+1
         
      case(3,4)!annual or semiannual cosine and sine
         if(para_typ(nshift+1) .eq. 3)then
            ohmtmp=Ohm
         else
            ohmtmp=2*Ohm
         end if
         
         do,i=1,nf
            A(i,nshift+1)=cos(ohmtmp*time(i))
            A(i,nshift+2)=sin(ohmtmp*time(i))
         end do
         nshift=nshift+2
      case(5) ! polynomial fit
         do,j=0,poly_ord
            do,i=1,nf
               A(i,nshift+j+1)=(time(i)-t0)**(j)
            end do
         end do
         nshift=nshift+poly_ord+1
     case(6)!monthly climatology
        do,i=1,nf
            epoch=FTime_fromDecYr(time(i))
            !write(0,*)fmonth(epoch)
            A(i,nshift+fmonth(epoch))=1.d0
        end do
        nshift=nshift+nmonths
     end select

      end do


   

   allocate(clm_f(2*pos,npara),clm_sigf(2*pos,npara))

   if(sigscale)then
      allocate(sigma0(2*pos))
      sigma0=1.d0 ! default posteriori scale of error-variance 
   end if

   clm_f=0.d0
   clm_sigf=0.d0

   !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(clm,clm_sig,clm_sigf,clm_f,A,sigma0) &
   !$OMP FIRSTPRIVATE(pos,posm,lmax,nf,npara,sigscale,weight,sigscaleresi)
   !allocate goodies in each thread
   allocate(Atmp(nf,npara),N(npara,npara),rhs(npara),input(nf))

   !$OMP DO
   do,i=1,2*pos    !loop over the coefficients (may be parallized)
      if(i>pos)then !sine order zero coefficient shall be skipped
         call SH_lm(i-pos,l,m)
         if(m.eq.0)cycle
      end if
      if(i>=1 .and. i<=posm)cycle ! restriction by minimum degree
      if(i>=pos+1 .and. i<=pos+posm)cycle ! restriction by minimum degree
      
      if(weight)then !decorrelate design matrix and coefficients
         do,j=1,nf
            Atmp(j,:)=A(j,:)/clm_sig(i,j)
            input(j)=clm(i,j)/clm_sig(i,j)
         end do
      else ! just copy
         Atmp=A
         input=clm(i,:)
      end if

      !build normal matrix
      call dsyrk('U','T',npara,nf,1.d0,Atmp,nf,0.d0,N,npara)   
      !the right hand side
      call dgemv('T',nf,npara,1.d0,Atmp,nf,input,1,0.d0,rhs,1)
      ! apriori residual fit
      ltpl=ddot(nf,input,1,input,1)
    
      !possibly apply a regularization
      if (usereg)then
            regscale=Mtrace(N)/Mtrace(Nreg)
            N=N+regscale*Nreg 

      end if
      !solve the goodies
      call solve_norm(N,rhs,ltpl)

      !copy solution
      clm_f(i,:)=rhs
      ! and propagated error
      if(sigscale)then
         sigma0(i)=sqrt(ltpl/(nf-npara))
         do,j=1,npara
            clm_sigf(i,j)=sqrt(N(j,j))*sigma0(i)
         end do
         if(sigscaleresi)then !possibly rescale the original input errors
            do,j=1,nf
               clm_sig(i,j)=clm_sig(i,j)*sigma0(i)
            end do
         end if
      else
         do,j=1,npara
            clm_sigf(i,j)=sqrt(N(j,j))
         end do
      end if

      !overwrite the original coefficients with the residual
      call dgemv('N',nf,npara,1.d0,A,nf,rhs,1,0.d0,input,1)
      do,j=1,nf
         clm(i,j)=clm(i,j)-input(j)
      end do
   end do
   !$OMP END DO
   deallocate(Atmp,N,rhs,input)
   !$OMP END PARALLEL

! write the results to files   
   do,i=1,npara
      fileout=trim(basen)//trim(para_des(i))//'.sh'
      call SH_write(filen=trim(fileout),clm=clm_f(stc:ndc,i),&
           &slm=clm_f(sts:nds,i),clm_sig=clm_sigf(stc:ndc,i),&
           &slm_sig=clm_sigf(sts:nds,i),typ=otyp)
   end do

end if

if(resSH)then !also output residual SH files
   
   do,i=1,nf !loop over files
      call basenamef(stripped=fileout,path=trim(filen(i)))
      fileout=trim(basen)//trim(fileout)
      !restrict minimum degree when needed
      if(posm >= 1)then
         clm(1:posm,:)=0.d0
         clm(pos+1:pos+posm,:)=0.d0
      end if
      
      write(stderr,*)"writing residual file: ",trim(fileout)
      !write residuals to file
      if(weight)then
               !write to file (as errors with zero mean!)
         call SH_write(filen=trim(fileout),clm=clm(stc:ndc,i),&
              &slm=clm(sts:nds,i),clm_sig=clm_sig(stc:ndc,i),slm_sig=clm_sig(sts:nds,i)&
              &,tcent=time(i),tstart=stime(i),tend=etime(i),typ=otyp)
      else
         call SH_write(filen=trim(fileout),clm=clm(stc:ndc,i),&
              &slm=clm(sts:nds,i),tcent=time(i),tstart=stime(i),tend=etime(i),typ=otyp)

      end if

   end do

end if

if(sigscale)then ! output a SH file with the sigma0's
   fileout=trim(basen)//'POSTSIGMA0.sh'
   call SH_write(filen=trim(fileout),clm=sigma0(stc:ndc),slm=sigma0(sts:nds),typ=otyp)
end if

if(rmsoutput)then
   ! compute the variance of the residual
   allocate(clmrms(2*pos))
   clmrms=0.d0
   do,i=1,2*pos    !loop over the coefficients ( may be parallized)
      if(i>pos)then !sine order zero coefficient shall be skipped
         call SH_lm(i-pos,l,m)
         if(m.eq.0)cycle
      end if
      if(i>=1 .and. i<=posm)cycle ! restriction by minimum degree
      if(i>=pos+1 .and. i<=pos+posm)cycle ! restriction by minimum degree
      
      do,j=1,nf
         clmrms(i)=clmrms(i)+(clm(i,j)**2)/nf
      end do
      !take the square root
      clmrms(i)=sqrt(clmrms(i))
   end do
   allocate(input(pos))
   input=0.d0 !dummy vector containing zeros

   fileout=trim(basen)//'POSTRMS.sh'
   !write to file (as errors with zero mean!)
   call SH_write(filen=trim(fileout),clm=input,slm=input,clm_sig=clmrms(stc:ndc),slm_sig=clmrms(sts:nds),tcent=t0,typ=otyp)

end if




end program SH_fit


!!help function
subroutine help()
integer::stderr
character(6)::frmt
stderr=0
frmt='(A)'
write(stderr,frmt)'SH_fit computes fits through sets of Spherical harmonic coefficients'
write(stderr,frmt)'usage SH_fit [options] SHFILES'
write(stderr,frmt)"Where SHFILES are the files to be fitted ( each has a different epoch)"
write(stderr,frmt)"The fits and residual variance are written to files starting with the"
write(stderr,frmt)"base SHFIT_ ( may be changed below)"
write(stderr,frmt)'Where OPTIONS may be:'
write(stderr,frmt)' -l=LMAX,LMIN: limit solution to coefficients between degree LMIN and LMAX'
write(stderr,frmt)' -w: use a weighted fit (standard error deviations must be present)'
write(stderr,frmt)" -fm: Fit a mean through the coefficients"
write(stderr,frmt)" -ft[=t0]: Fit a trend through the coefficients"
write(stderr,frmt)"  optionally provide a different center time t0 (defaults to 2003)"
write(stderr,frmt)" -fa: Fit an annual harmonic through the coefficients"
write(stderr,frmt)" -fs: Fit an (semi)annual harmonic through the coefficients"
write(stderr,frmt)" -fp[=ORD,T0]: fit a polynomial of order ORD (default=2) through the data"
write(stderr,frmt)"  Note: cannot be combined with -ft and -fm since it causes a rank defect"
write(stderr,frmt)" -fc: Fit a monthly climatology through the data (mean jan,feb,etc..)"
write(stderr,frmt)"  Multiple -f options are allowed"
write(stderr,frmt)"  -F=BASENAME: use a different base for the output (default is SHFIT_)"
write(stderr,frmt)"  -s: scale the errors by the posteriori sigma computed"
write(stderr,frmt)" In addition, a file with the posteriori sigmas is written"
write(stderr,frmt)"  -r Output the rms of the post-fit residual (also works without an -f option)"
write(stderr,frmt)"  -R[s] Also out put the residuals in separate SH files"
write(stderr,frmt)"   The 's' option rescale the input errors with the posteriori sigma, and automatically implies the -s option"
write(stderr,frmt)"  See also: fit1D"
!$ write(*,frmt)'This version is compiled with OpenMP and allows multi-threading (please set the OMP_NUM_THREADS env.variable)'
stop
end subroutine help

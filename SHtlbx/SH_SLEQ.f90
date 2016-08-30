!!Program which performs the sea level equations on timvarying viscous loads
!! Features"
!!*the Program uses a Load Love number formalism (time dependent Love numbers are computed externally)
!!*Solves the Sea Level equation in the Spectral domain
!! that's all for now
!!Coded by Roelof Rietbroek, Wed Mar 27 20:07:52 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de


program SH_SLEQ
use binfiletools
use shtlbx
use shtools
use netcdf
use grdtlbx
use forttlbx
use sh_synthesis
use gia_func
implicit none
integer::itharg,narg,i,lmax,lmaxice,lmaxl,stderr,ltyp,j
character(2)::FR
character(200)::dum,dir,Lovefelas
type(BINdat)::out,load,love
type(geogrid)::topo
logical::rotfeedback,verbose
integer::iargc,ncid,nt,l,m,numFR
double precision::dt,tmin,tmax,h1dum,k1dum,l1dum
double precision,pointer::lon(:),lat(:)
double precision,allocatable,dimension(:,:)::Ice,QSea,dIce,dSea,dQsea,dU,dV,dN
double precision,allocatable,dimension(:,:)::hn,ln,kn
double precision,pointer,dimension(:,:)::Oce
double precision,allocatable,dimension(:)::timeload,time,dtime
double precision,pointer::icef(:,:),lovef(:,:)
double precision,allocatable,dimension(:)::clmdum,slmdum
integer,allocatable,dimension(:)::deg,ord,trig
integer::nsh,nshoce,nlat,nlon,trigdum,ind
integer::glaciter,iglac,it
type(SHsynth_struct)::SHS

!!defaults
glaciter=1 !maximum number of iterations over the Glacial cycles
stderr=0
ltyp=-1 !use Love number from the viscous file with smallest time
load%file=''
love%file=''
FR='CM' ! center of common mass frame
lmax=0
rotfeedback=.false.
verbose=.false.
call getenv('WORK_DIR',dir)
topo%file=trim(dir)//'/data/etopo1/etopo1_0125deg.nc'
dt=1
load%mtyp='F'
love%mtyp='F'
Lovefelas=''

!!process command line options
narg=iargc()
if(narg <1) call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:3))
      case('L=')
         if(dum(4:4) .eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,dum)
            love%file=trim(dum)
         else
            love%file=trim(dum(4:))
         end if
      case('E=') ! specify purely elastic load love number model for the elastic limit
         read(dum(4:),*)ltyp
         if(ltyp .eq. 0)then
            Lovefelas=trim(dum(5:))
         end if
      case('F=') ! specify purely elastic load love number model for the elastic limit
         FR=dum(4:5)
      case('r')
         rotfeedback=.true.
      case('T=')
         topo%file=trim(dum(4:))
      case('t=')
         read(dum(4:),*)dt
      case('l=')!possibly restrict maximum degree
         read(dum(4:),*)lmax
      case('v')!be verbose
         verbose=.true.
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
      load%file=trim(dum)
   end if
end do

!input checks
if(love%file .eq. '')then
   write(stderr,*)'ERROR: No loading Love numbers provided'
   stop
end if

if(load%file .eq. '')then
   write(stderr,*)'ERROR: No load history provided'
   stop
end if

if(rotfeedback)then
   write(stderr,*)'WARNING: The -r option is currently a dummy'
   stop
end if

select case(FR)
case('CE')
   numFR=1
case('CM')
   numFR=2
case('CF')
   numFR=3
case('CL')
   numFR=4
case('CH')
   numFR=5
case default
   write(stderr,*)'ERROR: unknown reference frame origin specified, quitting'
   stop
end select

!read in loading history
call read_BINtype(load)

!get the maximum degree from the integer data
lmaxice=0
do,i=1,load%nint
   if(load%ints_d(i) .eq. 'Lmax')lmaxice=load%ints(i)
end do

if(lmax >0 .and. lmax >lmaxice)then
   write(stderr,*)'ERROR: the maximum degree requested is not supported by the load'
   stop
end if

if(lmax.eq. 0)lmax=lmaxice !if not explicitly set, take the maximum degree of the file


!extract the epochs from the loading history
allocate(timeload(load%nval2))
do, i=1,load%nval2
   read(load%side2_d(i)(5:),*)timeload(i)
end do

!make a time axis containing the internal discretization
tmin=minval(timeload)
tmax=maxval(timeload)
nt=int((tmax-tmin)/dt)+1

!optionally adjust tmin to fit the boundaries
tmin=tmax-dt*(nt-1)
allocate(time(nt),dtime(nt))
do,i=1,nt
   time(i)=tmin+dt*(i-1)
   dtime(i)=(i-1)*dt !needed time differences
end do


!load the Visco elastic Love numbers
call read_BINtype(love)

!get the maximum degree from the integer data
lmaxl=0
do,i=1,love%nint
   if(love%ints_d(i) .eq. 'Lmax')lmaxl=love%ints(i)
end do

if(lmax >lmaxl)then
   write(stderr,*)'ERROR: the maximum degree requested is not supported by the Love numbers'
   stop
end if

!get amount of SH coefficients (stacked form)
nsh=SH_tpos(l=lmax,m=lmax,q=1,lmax=lmax,lmin=0)
!for full spectral resolution one needs to resolve the OCean function up to degree 2*lmax
nshoce=SH_tpos(l=lmax*2,m=lmax*2,q=1,lmax=lmax*2,lmin=0)
allocate(Oce(nshoce,nt))





!interpolate ice history
call interp1(x=timeload,y=load%mat1,xi=time,yi=icef,dim=2,method='linear')

!set up state vectors
allocate(Ice(nsh,nt),QSea(nsh,nt)) !absolute values
allocate(dIce(nsh,nt),dSea(nsh,nt),dQSea(nsh,nt),dU(nsh,nt),dV(nsh,nt),dN(nsh,nt)) !incremental values
Ice=0.d0 !absolute ice load
QSea=0.d0 !absolute (quasi spectral sea level)
dSea=0.d0 !incremental relative sea level (masked over the Ocean)
dQsea=0.d0 !incremental Quasi spectral sea level
dU=0.d0 !incremental uplif
dV=0.d0 !incremtal horizontal deformation
dN=0.d0 !incremental geoid

allocate(deg(nsh),ord(nsh),trig(nsh))

do,i=1,load%nval1 !loop over SH coefficients
   read(load%side1_d(i),'(4X,I3,I3)')l,m
   if(l>lmax)cycle
   if(load%side1_d(i)(2:2) .eq. 'C')then
      trigdum=0
   else
      trigdum=1
   end if
   ind=SH_tpos(l=l,m=m,q=trigdum,lmax=lmax,lmin=0)
   deg(ind)=l
   ord(ind)=m
   trig(ind)=trigdum
   Ice(ind,:)=icef(i,:)
end do

!set first time difference of the Love numbers to the elastic limit
love%vec(1,1)=0.d0 !this will prevent the interpolation routine from winching


!interpolate Love numbers to appropriate time differences
call interp1(x=love%vec(:,1),y=love%mat1,xi=dtime,yi=lovef,dim=1,method='linear')

!make separate matrices with load Love numbers
allocate(hn(lmax+1,nt),ln(lmax+1,nt),kn(lmax+1,nt))
hn=0.d0
ln=0.d0
kn=0.d0
do,i=1,size(lovef,2)
   read(love%side2_d(i),'(4x,I4)')l
   select case(love%side2_d(i)(1:2))
   case('hn')
      hn(l+1,:)=lovef(:,i)
   case('ln')
      ln(l+1,:)=lovef(:,i)
   case('kn')
      kn(l+1,:)=lovef(:,i)
   case default
      write(stderr,*)"ERROR: error processing Love number description: ",love%side2_d(i)
   end select
end do
deallocate(lovef) !since it is not needed anymore (although small)

!optionally replace the Love numbers in the elastic limit by compressible ones
if(ltyp >= 0)then
   call SH_loadlove(hnm=hn(:,1),knm=kn(:,1),lnm=ln(:,1),typ=ltyp,frame=numFR,iso=.true.,fileop=Lovefelas)
end if


!adopt appropriate reference system for every time step
do, i=1,nt
   call SH_Love_d1_trans(h1in=hn(i,2),l1in=ln(i,2),k1in=kn(i,2),h1out=h1dum,l1out=l1dum,k1out=k1dum,FR=FR)
   hn(2,i)=h1dum
   ln(2,i)=l1dum
   kn(2,i)=k1dum
end do


do, iglac=1,glaciter

!make the Ocean functions
   call GIA_Ocean_func(topo=topo,lmax=lmax,Tice=Ice(:,1:1),Sea=Qsea(:,1:1)&
        &,Oce=Oce(:,1:1),SHS=SHS)
   

ind=SH_pos(2*lmax,2*lmax)
allocate(clmdum(ind),slmdum(ind))
slmdum=0.d0
clmdum=0.d0
do, i=1,1
   do,l=0,2*lmax
      do,m=0,l
         ind=SH_tpos(l,m,0,lmax*2,0)
         clmdum(SH_pos(l,m))=Oce(ind,1)
         if(m>0)then
            ind=SH_tpos(l,m,1,lmax*2,0)
            slmdum(SH_pos(l,m))=Oce(ind,1)
         end if
      end do
   end do

   write(dum,'(A3,I2.2,A3)')'Oce',i,'.sh'
   call SH_write(clm=clmdum,slm=slmdum,filen=trim(dum))
end do
!contruct the Ocean Product-2-sum matrices
!call SH_Prod2Sum_mats(mats=Omats,lmax=lmax,SHin=Oce)

!construct the ice load increments ( and remove floating ice)
do,i=1,nt
   
end do


   do,it=1,nt !loop over the time points
      do,i=1,it !apply convolution of the load
         
         
      end do

      !solve the Sea Level Equation 


   end do ! end loop over the time points


!check for convergence

end do ! end loop over glacial cycles

!construct the observables and output to file(s)



!TO BE COntinued


! read in 

! !!setup oce matrix
! oce%mtyp='U' ! symmetric matrix


! !!read in meta data

! call read_BINtype(oce,2)

! !! check if the matrxi type is correct
! select case(oce%type)
! case('SYMVN___','SYMV0___','SYMV1___','SYMV2___')
!    !read the remainder
!    call read_BINtype(oce)
! case default
!    write(stderr,*)'ERROR: Ocean matrix is not symmetric??'
!    stop 
! end select

! if(oce%side1_d(1)(2:10) .ne. 'CN   0  0')then
!    write(stderr,*)'ERROR: Ocean matrix side must start with 00 component'
!    stop 
! end if

! !!determine lmax
! !!and create an index vector

! ! do,i=1,oce%nint
! !    if(trim(oce%ints_d(i)) .eq. 'Lmax')lmax=oce%ints(i)
! ! end do

! !if not found try to retrieve from side description
! allocate(deg(oce%nval1),ord(oce%nval1),trig(oce%nval1))
! do,i=1,oce%nval1
!    read(oce%side1_d(i),'(4x,i3,i3)')deg(i),ord(i)
!    if(oce%side1_d(i)(2:2).eq. 'C')then
!       trig(i)=0
!    else
!       trig(i)=1
!    end if
!    lmax=max(deg(i),lmax)
! end do

! if(regmat)then
!    allocate(out%mat1(oce%nval1,oce%nval2)) ! needs a working copy of the oceanmat
!    out%mat1=oce%mat1 !copy values
! else
!    out%mat1=>oce%mat1 ! point to the right matrix ( no addtional memeory necessary)
! end if

! !!mirror matrix
! do,i=1,oce%nval1
!    do,j=1,i-1
!       out%mat1(i,j)=out%mat1(j,i)
!    end do
! end do

! !get ocean ratio (RR not needed anymore)
! !A_00=oce%mat1(1,1)

! !!load load love numbers
! pos=SH_pos(lmax,lmax)
! allocate(hnm(pos),lnm(pos),knm(pos))
! call SH_loadlove(hnm=hnm,knm=knm,lnm=lnm,typ=ltyp,frame=FR)

! if(rotfeedback)then ! load body love numbers if required
!    !load tidal love numbers (k and h) (maximum of degree 2 only)
!    pos=SH_pos(2,2)
!    allocate(knm_T(pos),hnm_T(pos))
!    call SH_loadlove(hnm=hnm_T,knm=knm_T,typ=4,frame=1) !note frame is irrelevant
!    gamma2=-(1.d0+knm_T(pos)-hnm_T(pos))/g
! end if


! !!output systems
! out%file='stdout'
! out%nval1=oce%nval1
! out%nval2=out%nval1
! out%pval2=1

! if(regmat)then
!    out%nread=4
!    if(plaw)out%nread=out%nread+2
!    if(regusef)out%nread=out%nread+2
!    out%nint=4
!    out%ndbls=5
!    out%nvec=3 ! store used love numbers in a vectors
! else
!    out%nread=13
!    out%nint=2
!    out%ndbls=3
!    out%nvec=3 ! store used love numbers in a vectors
! end if
   
! if(rotfeedback)then
!    out%ndbls=out%ndbls+2
!    out%nread=out%nread+1
! end if




! allocate(out%vec(out%nval1,out%nvec))
! allocate(out%readme(out%nread))
! allocate(out%ints_d(out%nint),out%ints(out%nint))
! allocate(out%dbls_d(out%ndbls),out%dbls(out%ndbls))
! out%readme=''


! if(regmat)then

!    out%type='SYMVN___'
!    out%mtyp='U'
!    out%pval1=out%nval1*(out%nval1+1)/2
!    out%descr="Regularization matrix regularizing surface load minus self consistent sea level"
!    if(plaw)then
!       inr=inr+1
!       out%readme(inr)="Power law used for error-covariance:"
!       inr=inr+1
!       write(out%readme(inr),*)sc,"n^",pn
!    end if
!    if(regusef)then
!       inr=inr+1
!       out%readme(inr)="File used for diagonal error-covariance:"
!       inr=inr+1
!       out%readme(inr)=trim(regfile)
!    end if

!    inr=inr+1
!    out%readme(inr)="Ocean file used:"
!    inr=inr+1
!    out%readme(inr)=oce%file(1:80)
!    inr=inr+1
!    out%readme(inr)="Vector data contain dummy (zero) value a apriori and right hand side"
   

!    if(rotfeedback)then
!       inr=inr+1
!       out%readme(inr)="Constructed using linearized rotational feedback of the sea level"
!    end if
!    ii=ii+1
!    out%ints_d(ii)='Lmax'
!    out%ints(ii)=lmax
!    ii=ii+1
!    out%ints_d(ii)='Lmin'
!    out%ints(ii)=0
!    ii=ii+1
!    out%ints_d(ii)='Nobs'
!    out%ints(ii)=out%nval1
!    ii=ii+1
!    out%ints_d(ii)='Nunknowns'
!    out%ints(ii)=out%nval1

!    idb=idb+1
!    out%dbls_d(idb)='Rho_w'
!    out%dbls(idb)=rho_w
!    idb=idb+1
!    out%dbls_d(idb)='Rho_e'
!    out%dbls(idb)=rho_e
!    idb=idb+1
!    out%dbls_d(idb)='LtPL'
!    out%dbls(idb)=0.d0
!    idb=idb+1
!    out%dbls_d(idb)='Sigma0'
!    out%dbls(idb)=1.d0

   
!    if(rotfeedback)then
!       idb=idb+1
!       out%dbls_d(idb)="k2_body"
!       out%dbls(idb)=knm_T(SH_pos(2,0))
!       idb=idb+1
!       out%dbls_d(idb)="h2_body"
!       out%dbls(idb)=hnm_T(SH_pos(2,0))
!    end if

! else
!    out%type='FULLSQVN'
!    out%mtyp='F'
!    out%pval1=out%nval1*out%nval2
!    out%descr="sea level greens function from SH_sealevmat"
!    inr=inr+1
!    out%readme(inr)='Spectral Greens function for the sea level equation'
!    inr=inr+1
!    out%readme(inr)='which is mass consistent and self gravitating'
!    inr=inr+1
!    out%readme(inr)='File contains the matrix  M'

!    inr=inr+1
!    out%readme(inr)="Which relates quasi spectral sea level S' to a forcing F"
!    inr=inr+1
!    out%readme(inr)='in the following way:'
!    inr=inr+1
!    if(invert)then
!       out%readme(inr)=" M * F = S'"
!    else
!       out%readme(inr)=" M * S' = F"
!    end if
!    inr=inr+1
!    out%readme(inr)="The Forcing function expresses the sea level caused by the"

!    inr=inr+1
!    out%readme(inr)="external forcing only but not by the sea level itself"
!    inr=inr+1   
! !   out%readme(inr)="F_00 denotes the eustatic sea level change (mass taken from ocean)"
!    out%readme(inr)="F_00 denotes MASS TAKEN FROM THE OCEAN per steradian"
!    inr=inr+1
!    out%readme(inr)="Ocean file used:"
!    inr=inr+1
!    out%readme(inr)=oce%file(1:80)
!    inr=inr+1
!    out%readme(inr)="load love set from the model are put in vectordata (hnm,lnm,knm)"
   

!    if(rotfeedback)then
!       inr=inr+1
!       out%readme(inr)=&
!      "Constructed using linearized rotational feedback of the sea level"
!    end if
!    ii=ii+1
!    out%ints_d(ii)='Lmax'
!    out%ints(ii)=lmax
!    ii=ii+1
!    out%ints_d(ii)='Lmin'
!    out%ints(ii)=0
   
!    idb=idb+1
!    out%dbls_d(idb)='Rho_w'
!    out%dbls(idb)=rho_w
!    idb=idb+1
!    out%dbls_d(idb)='Rho_e'
!    out%dbls(idb)=rho_e

      
!    if(rotfeedback)then
!       idb=idb+1
!       out%dbls_d(idb)="k2_body"
!       out%dbls(idb)=knm_T(SH_pos(2,0))
!       idb=idb+1
!       out%dbls_d(idb)="h2_body"
!       out%dbls(idb)=hnm_T(SH_pos(2,0))
!    end if

! end if








! !second side description
! !adapt original side description data

! allocate(oce%side2_d(oce%nval2)) ! allocate the second dside description vector

! if(regmat)then
!    do,i=1,oce%nval1
!       oce%side1_d(i)(1:1)='T' ! change G in T (Land surface load)
!       oce%side2_d(i)=oce%side1_d(i)
!    end do

! else

!    do,i=1,oce%nval1
!       oce%side1_d(i)(1:1)='T' ! change G in T (equivalent water height)
!       oce%side1_d(i)(17:)='SEALEVEL'
!       oce%side2_d(i)=oce%side1_d(i)
!       oce%side2_d(i)(16:24)='SEAQUASI'
!    end do

! end if
! out%side1_d=>oce%side1_d

! if(invert)then
!     out%side2_d=>oce%side1_d
!     out%side1_d=>oce%side2_d
! else
!     out%side2_d=>oce%side2_d
!     out%side1_d=>oce%side1_d
! end if

! !sort and put love numbers in vector
! out%vec(1,:)=1/1.d-99 !( Put  for the degree zero love numbers)
! do,i=2,oce%nval1
!    pos=SH_pos(deg(i),ord(i))
! !   write(0,*)pos,deg(i),ord(i)
!    out%vec(i,1)=hnm(pos)
!    out%vec(i,2)=lnm(pos)
!    out%vec(i,3)=knm(pos)
! end do


! !calculate changes due to rotational feedback
! if(rotfeedback)then
!    !construct matrix relating surface loading changes to 'apparent' potential changes
!    !based on linearized conservation of momentum
!    Ji3tmp=Ji3_2_m(ltyp)
!    RotMat=matmul(m_2_rotpot(),matmul(Ji3tmp,SurfL_2_Ji3()))
!    allocate(Rotrows(4,out%nval1)) ! array to temporary store the effect of rotational feedback


!    !create index vector (length 4, c00,c20,c21,s21)
!    pos=0
!    i=0
!    do while(pos< 4 .and. i<oce%nval1)
!       i=i+1
!       select case(deg(i))
!       case(0)!degree and order 0
!          rotind(1)=i
!          pos=pos+1
!       case(2)
!          select case(ord(i))
!          case(0)!degree 2 order 0
!             rotind(2)=i
!             pos=pos+1
!          case(1) !degree 2 order 1
!             if(trig(i) .eq. 0)then ! cosine
!                rotind(3)=i
!                pos=pos+1
!             else
!                rotind(4)=i
!                pos=pos+1
!             end if
!          end select
!       end select
!    end do

!    !Copy ocean matrix values in RotRows
!    do,j=1,out%nval1
!       do,i=1,4
!          Rotrows(i,j)=oce%mat1(rotind(i),j)
!       end do
!    end do
   
!    !calculate matrix multiplication
!    Rotrows=matmul(Rotmat,Rotrows)
!    Rotrows=Rotrows*gamma2

! end if


! !adapt oce%mat1 step 1 EARTH MODEL dependent part
! !convolve non zero degree part with G_N-U greens function ( geoid minus uplift)
! do,i=1,out%nval1 !loop over columns
!    do,j=2,out%nval1 ! loop over rows ( skip degree 0)
!       !calculate row dependent factor
!       fact=-(3.d0*rho_w/rho_e)*(1.d0+out%vec(j,3)-out%vec(j,1))/(2*deg(j)+1)    
!       out%mat1(j,i)=fact*out%mat1(j,i)
!    end do
!    !simply negate the first row
!    out%mat1(1,i)=-out%mat1(1,i)
! end do


! if(rotfeedback)then !add the matrix from rotational feedback
!    !update matrix with rotational feedback part
!    do,j=1,out%nval1
!       do,i=2,4 !skip degree 0 to ensure mass consistency ( is small anyway)
!          out%mat1(rotind(i),j)=out%mat1(rotind(i),j)+Rotrows(i,j)
!       end do
!    end do
   

! end if


! !degree 0 row has a special treatment ( divide by O_00 to retrieve eustatic sea level)
! ! do, i=1,oce%nval1
! !    out%mat1(1,i)=out%mat1(1,i)/A_00
! ! end do

! !add unit diagonal
! do,i=2,oce%nval1
!    out%mat1(i,i)=1+out%mat1(i,i)
! end do



! !!invert if requested ( LU decomposition)
! if(invert)then
!    call cpu_time(time1)
!    allocate(ipiv(out%nval1))

!    !calculate matrix L1 norm
!    anorm=0.d0
!    do,i=1,out%nval2 ! loop over columns
!       dumdb=0.d0
!       do,j=1,out%nval1 !loop over rows
!          dumdb=dumdb+abs(out%mat1(j,i))
!       end do
!       anorm=max(anorm,dumdb)
!    end do
   
!    allocate(workv(4*out%nval2))
!    allocate(iwork(out%nval2))

!    !first LU factorization
!    call dgetrf(out%nval1,out%nval1,out%mat1,out%nval1,ipiv,info)

!    if( info .ne. 0)then
!       write(stderr,*)"ERROR: factorization failed"
!       stop
!    end if
!    rcond=0.d0


! !    ! calculate condition number
!     call dgecon('1',out%nval2,out%mat1(1,1),out%nval1,anorm,&
!                rcond,workv,iwork,info)

!    if( info .ne. 0)then
!       write(stderr,*)"ERROR: calcuating condition number"
!       stop
!    end if

!    idb=idb+1
!    out%dbls_d(idb)='Cond_No_L1'
!    out%dbls(idb)=1/rcond

!    !invert 
!    !optimal lwork:

!    lwork=out%nval1*ilaenv( 1, 'DGETRI', ' ', out%nval1, -1, -1, -1 )

!    allocate(work(lwork))
!    call dgetri(out%nval1,out%mat1,out%nval1,ipiv,work,lwork,info)
!    if( info .ne. 0)then
!       write(stderr,*)"ERROR: inversion failed"
!       stop
!    end if
!    call cpu_time(time2)
!    inr=inr+1
!    write(out%readme(inr),*)"Matrix inversion took",time2-time1,"seconds" 


! end if

! if(regmat)then !apply post processing to obtain a regularization matrix
!    !When the above matrix is called G apply
!    ! M =  O (I - G * D* (I-O) )  
!    ! and perform a rank k update :
!    !Mout= M' * M


!    if(rotfeedback)then
!       deallocate(Rotrows)
!       allocate(Rotrows(out%nval1,4))
      
!       !Copy matrix values in RotRows
!       do,i=1,4
!          do,j=1,out%nval1
!             Rotrows(j,i)=out%mat1(j,rotind(i))
!          end do
!       end do
      
!       !special treatment of degree 0
!       Rotmat(:,1)=Rotmat(:,1) !-Rotmat(:,1)/A_00
!       !calculate right side matrix multiplication
      
!       Rotrows=matmul(Rotrows,Rotmat)
!       Rotrows=Rotrows*gamma2

!    end if
   
!    ! right multiply with a diagonal matrix
!    !convolve non zero degree part with G_N-U greens function ( geoid minus uplift)
!    do,i=1,out%nval1 !loop over columns
!       if(deg(i)==0)then !special case for degree 0
!          fact=1.d0 !-1/A_00
!       else
!          fact=(3.d0*rho_w/rho_e)*(1.d0+out%vec(i,3)-out%vec(i,1))/(2*deg(i)+1)    
!       end if

!       do,j=1,out%nval1 ! loop over rows 
!          out%mat1(j,i)=fact*out%mat1(j,i)
!       end do

!    end do

!    if(rotfeedback)then !add the matrix from rotational feedback
!       !update matrix with rotational feedback part
!       do,i=1,4
!          do,j=1,out%nval1
!             out%mat1(j,rotind(i))=out%mat1(j,rotind(i))+Rotrows(j,i)
!          end do
!       end do
   

!    end if


!    ! right multiplication with -(I-O) = O-I
   
!    ! allocate matrix and copy values
!    allocate(work2d(out%nval1,out%nval1))
!    work2d=out%mat1 !copy values
   
!    call dsymm('R',oce%mtyp,out%nval1,out%nval2,1.d0,oce%mat1(1,1),oce%nval1,out%mat1(1,1),out%nval1,-1.d0,work2d(1,1),out%nval1)
   
!    ! add unit diagonal
!    do ,i=1, out%nval1
!       work2d(i,i)=work2d(i,i)+1.d0
!    end do

!    !left multiply by the ocean matrix
!    call dsymm('L',oce%mtyp,out%nval1,out%nval2,1.d0,oce%mat1(1,1),oce%nval1,work2d(1,1),out%nval1,0.d0,out%mat1,out%nval1)
   
!    if(plaw)then !apply degree dependent error-covariance
!       do,i=1,out%nval1
!          out%mat1(i,1:out%nval1)=out%mat1(i,1:out%nval1)/sqrt(sc*deg(i)**pn)
!       end do
!    end if

!    if(regusef)then
!       stc=1
!       ndc=pos
!       sts=pos+1
!       nds=2*pos
!       allocate(clmtmp(2*pos),clm_sig(2*pos))
!       call SH_readgrav(filen=trim(regfile),clm=clmtmp(stc:ndc),&
!            &slm=clmtmp(sts:nds),clm_sig=clm_sig(stc:ndc),&
!            slm_sig=clm_sig(sts:nds))
!       do,i=1,out%nval1
!          ind=SH_pos(deg(i),ord(i))
!          out%mat1(i,1:out%nval1)=out%mat1(i,1:out%nval1)/(clm_sig(ind+trig(i)*pos))
!       end do

!    end if


! !   DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
!    ! rank k update
!     call dsyrk('U','T',out%nval1,out%nval1,1.d0,out%mat1(1,1),out%nval1,0.d0,work2d(1,1),out%nval1)

!     deallocate(out%mat1)
!     out%mat1=>work2d !point to updated matrix ( warning may yield a dangling pointer, when out%mat1 is not deallocated)

!     !reset vec part to zero
!     out%vec=0.d0
!     out%nvec=2

!  end if



! !write to file

! call write_BINtype(out)


end program SH_SLEQ


subroutine help()
implicit none
integer::unit=6
character(8)::frmt='(A)'
character(8)::version
version="ALPHA 0"
write(unit,frmt)"Program SH_SLEQ computes the time varying variations of sea level due to a (glacial) load"
write(unit,frmt)" Usage: SH_SLEQ -L=LOVENUMBERFILE  LOADINGFILE > OUTPUTFILE"
write(unit,frmt)"Where LOVENUMBERFILE contains the time dependent Love numbers"
write(unit,frmt)"and LOADINGFILE contains the (glacial load) as spherical harmonic coefficients"
write(unit,frmt)" The OUTPUTFILE is a binary file containing SH coefficients of QS,RS,N,U,V"
write(unit,frmt)" for every time step, with additional info."
write(unit,frmt)"OPTIONS may be:"
write(unit,frmt)"-r: incorporate rotational feedback"
write(unit,frmt)"-T=TOPOFILE: use a different topography file (netcdf format)"
write(unit,frmt)" The default uses a 0.125 degree equidistant sampled version of etopo1 (bedrock)"
write(unit,frmt)"-t=DKYR:Discretize in steps of DKYR kiloyears (Default uses 1 KYRS)"
write(unit,frmt)" The load an load Love numbers will be interpolated accordingly"
write(unit,frmt)"-F=FR: output the variations in a different frame (default is CM)"
write(unit,frmt)"       alternatively choose CF, CE, CL, CH"
write(unit,frmt)"-ENUM: optionally replace the elastic limit by a standard Love number set"
write(unit,frmt)" (see SH_dump -h for the NUM codes)"
write(unit,frmt)"-l=LMAX: restrict the input to LMAX degrees"
write(unit,frmt)"-h: display this help message"
write(unit,frmt)" VERSION: "//trim(version)//" (pretty much a random number generator)"
!$ write(*,frmt)'This version is compiled with OpenMP and allows multi-threading (please set the OMP_NUM_THREADS env.variable)'
stop
end subroutine help

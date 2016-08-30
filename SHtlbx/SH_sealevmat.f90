
!!Program which creates a matrix which solves the sea level equation for a given load distribution
!!the program requires a product 2 sum conversion matrix made by the program SH_prod2sum
!!Input:
!!Product 2 sum matrix for the ocean function O() (in binary form readable by read_BINtype)
!!optionally a different load love number set ( default uses PREM elastic load love numbers)
!!output is a full matrix which is either
!! M1:   M1 S = L
!!or 
!!M2:    M2 L = S               ( M2 = M1^-1)
!!where L is the direct sea level change due to loading/tides and S' is the quasi spectral sea level
!! S=O(S')
!! S' is smoother than S and is suitable for spatial interpretation, While S should be used for loading.
!!
!!
!!Coded by Roelof Rietbroek, Wed Jan 28 10:14:17 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Wed Sep 16 14:16:50 2009
!! removed bug ( use row scaling instead of column scaling)
!! L_00 input is expected to be eustatic sea level rise

!!Updated by Roelof Rietbroek, Mon Oct 12 12:03:14 2009
!! added calculation of condition number

!!Updated by Roelof Rietbroek, Wed Feb  3 15:21:17 2010
!!started adding rotational feedback mechanism

!!Updated by Roelof Rietbroek, Thu Aug 30 14:54:23 2012
!! added option to produce a regularization matrix, which regularizes the difference
!! between the global surface load and the self consistent sea level over the ocean
!! fixed bug in calculation of the rotational feedback (still needs to be checked properly)
!!Updated by Roelof Rietbroek, Sun Mar  3 13:18:41 2013
!!ensured mass conservation when rotational feedback is used

!!Updated by Roelof Rietbroek, Wed Mar  6 21:40:51 2013
!! allowed degree dependent power law in the regularization matrix
!!Updated by Roelof Rietbroek, Thu Mar  7 20:48:00 2013
!! output the regularization matrix as a normal equation system ( allows VCE estimation of the regularization parameter)

!!Updated by Roelof Rietbroek, Fri Mar 15 09:24:18 2013
!! also allows degree and order dependent weights to be read in
!! and fixed bug (forgot square root in error-covariance of the power law)

!!Updated by Roelof Rietbroek, Mon Mar 18 15:00:29 2013
!!the degree 0 rows is NOT scaled by the Ocean surface anymore!!
!!So the degree 0 coefficient is now not the eustatic sea level anymore but the workld average
!! new files can be recognized by a tag in the readme data

program SH_sealevmat
use binfiletools
use shtlbx 
implicit none
integer::itharg,narg,i,ltyp,j,pos,FR,lmax,stderr
character(200)::dum,regfile
type(BINdat)::out,oce
logical::invert,rotfeedback
integer::lwork,info,ilaenv
integer,allocatable,dimension(:)::ipiv,deg,ord,trig,iwork
double precision::fact,time1,time2
double precision, allocatable,dimension(:)::work,hnm,knm,lnm,knm_T,hnm_T
integer::iargc
double precision::anorm,dumdb,rcond !,A_00
double precision,allocatable::workv(:)
double precision,allocatable,target::work2d(:,:) 
double precision::Rotmat(4,4)
integer::rotind(4)
double precision::gamma2
logical::regmat,plaw,regusef
double precision,allocatable::Rotrows(:,:),clm_sig(:),clmtmp(:)
double precision::Ji3tmp(3,3)
integer::ind,stc,ndc,sts,nds
double precision::sc,pn
integer::inr,idb,ii !keeps track of the current position in the readme and double and integer meta data


!!dfefaults
stderr=0
ltyp=2
invert=.true.
oce%file=''
FR=3 ! center of figure frame
lmax=0
rotfeedback=.false.
regmat=.false.
plaw=.false.
inr=0
idb=0
ii=0
regfile=''
regusef=.false.

!!process command line options
narg=iargc()
if(narg <1) call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('E') ! specify load love number model
         read(dum(3:),*)ltyp
      case('n')
         invert=.false.
      case('r')
         rotfeedback=.true.
      case('R')
         regmat=.true.
         if(dum(3:4) .eq. 'f=')then
            regfile=trim(dum(5:))
            regusef=.true.
         else
            ind=index(dum,',')
            if(ind .ne. 0)then
               plaw=.true.
               read(dum(3:ind-1),*)sc
               read(dum(ind+1:),*)pn
            end if
         end if
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
      oce%file=trim(dum)
   end if
end do

!input checks
if(oce%file .eq. '')then
   write(stderr,*)'ERROR: No ocean matrix provided'
   stop
end if

if(plaw .and. regusef)then
   write(stderr,*)'ERROR: power law and errors from file may not be combined'
   stop
end if

if(regmat)invert=.true. ! change to true for sure


!!setup oce matrix
oce%mtyp='U' ! symmetric matrix


!!read in meta data

call read_BINtype(oce,2)

!! check if the matrxi type is correct
select case(oce%type)
case('SYMVN___','SYMV0___','SYMV1___','SYMV2___')
   !read the remainder
   call read_BINtype(oce)
case default
   write(stderr,*)'ERROR: Ocean matrix is not symmetric??'
   stop 
end select

if(oce%side1_d(1)(2:10) .ne. 'CN   0  0')then
   write(stderr,*)'ERROR: Ocean matrix side must start with 00 component'
   stop 
end if

!!determine lmax
!!and create an index vector

! do,i=1,oce%nint
!    if(trim(oce%ints_d(i)) .eq. 'Lmax')lmax=oce%ints(i)
! end do

!if not found try to retrieve from side description
allocate(deg(oce%nval1),ord(oce%nval1),trig(oce%nval1))
do,i=1,oce%nval1
   read(oce%side1_d(i),'(4x,i3,i3)')deg(i),ord(i)
   if(oce%side1_d(i)(2:2).eq. 'C')then
      trig(i)=0
   else
      trig(i)=1
   end if
   lmax=max(deg(i),lmax)
end do

if(regmat)then
   allocate(out%mat1(oce%nval1,oce%nval2)) ! needs a working copy of the oceanmat
   out%mat1=oce%mat1 !copy values
else
   out%mat1=>oce%mat1 ! point to the right matrix ( no addtional memeory necessary)
end if

!!mirror matrix
do,i=1,oce%nval1
   do,j=1,i-1
      out%mat1(i,j)=out%mat1(j,i)
   end do
end do

!get ocean ratio (RR not needed anymore)
!A_00=oce%mat1(1,1)

!!load load love numbers
pos=SH_pos(lmax,lmax)
allocate(hnm(pos),lnm(pos),knm(pos))
call SH_loadlove(hnm=hnm,knm=knm,lnm=lnm,typ=ltyp,frame=FR)

if(rotfeedback)then ! load body love numbers if required
   !load tidal love numbers (k and h) (maximum of degree 2 only)
   pos=SH_pos(2,2)
   allocate(knm_T(pos),hnm_T(pos))
   call SH_loadlove(hnm=hnm_T,knm=knm_T,typ=4,frame=1) !note frame is irrelevant
   gamma2=-(1.d0+knm_T(pos)-hnm_T(pos))/g
end if


!!output systems
out%file='stdout'
out%nval1=oce%nval1
out%nval2=out%nval1
out%pval2=1

if(regmat)then
   out%nread=4
   if(plaw)out%nread=out%nread+2
   if(regusef)out%nread=out%nread+2
   out%nint=4
   out%ndbls=5
   out%nvec=3 ! store used love numbers in a vectors
else
   out%nread=13
   out%nint=2
   out%ndbls=3
   out%nvec=3 ! store used love numbers in a vectors
end if
   
if(rotfeedback)then
   out%ndbls=out%ndbls+2
   out%nread=out%nread+1
end if




allocate(out%vec(out%nval1,out%nvec))
allocate(out%readme(out%nread))
allocate(out%ints_d(out%nint),out%ints(out%nint))
allocate(out%dbls_d(out%ndbls),out%dbls(out%ndbls))
out%readme=''


if(regmat)then

   out%type='SYMVN___'
   out%mtyp='U'
   out%pval1=out%nval1*(out%nval1+1)/2
   out%descr="Regularization matrix regularizing surface load minus self consistent sea level"
   if(plaw)then
      inr=inr+1
      out%readme(inr)="Power law used for error-covariance:"
      inr=inr+1
      write(out%readme(inr),*)sc,"n^",pn
   end if
   if(regusef)then
      inr=inr+1
      out%readme(inr)="File used for diagonal error-covariance:"
      inr=inr+1
      out%readme(inr)=trim(regfile)
   end if

   inr=inr+1
   out%readme(inr)="Ocean file used:"
   inr=inr+1
   out%readme(inr)=oce%file(1:80)
   inr=inr+1
   out%readme(inr)="Vector data contain dummy (zero) value a apriori and right hand side"
   

   if(rotfeedback)then
      inr=inr+1
      out%readme(inr)="Constructed using linearized rotational feedback of the sea level"
   end if
   ii=ii+1
   out%ints_d(ii)='Lmax'
   out%ints(ii)=lmax
   ii=ii+1
   out%ints_d(ii)='Lmin'
   out%ints(ii)=0
   ii=ii+1
   out%ints_d(ii)='Nobs'
   out%ints(ii)=out%nval1
   ii=ii+1
   out%ints_d(ii)='Nunknowns'
   out%ints(ii)=out%nval1

   idb=idb+1
   out%dbls_d(idb)='Rho_w'
   out%dbls(idb)=rho_w
   idb=idb+1
   out%dbls_d(idb)='Rho_e'
   out%dbls(idb)=rho_e
   idb=idb+1
   out%dbls_d(idb)='LtPL'
   out%dbls(idb)=0.d0
   idb=idb+1
   out%dbls_d(idb)='Sigma0'
   out%dbls(idb)=1.d0

   
   if(rotfeedback)then
      idb=idb+1
      out%dbls_d(idb)="k2_body"
      out%dbls(idb)=knm_T(SH_pos(2,0))
      idb=idb+1
      out%dbls_d(idb)="h2_body"
      out%dbls(idb)=hnm_T(SH_pos(2,0))
   end if

else
   out%type='FULLSQVN'
   out%mtyp='F'
   out%pval1=out%nval1*out%nval2
   out%descr="sea level greens function from SH_sealevmat"
   inr=inr+1
   out%readme(inr)='Spectral Greens function for the sea level equation'
   inr=inr+1
   out%readme(inr)='which is mass consistent and self gravitating'
   inr=inr+1
   out%readme(inr)='File contains the matrix  M'

   inr=inr+1
   out%readme(inr)="Which relates quasi spectral sea level S' to a forcing F"
   inr=inr+1
   out%readme(inr)='in the following way:'
   inr=inr+1
   if(invert)then
      out%readme(inr)=" M * F = S'"
   else
      out%readme(inr)=" M * S' = F"
   end if
   inr=inr+1
   out%readme(inr)="The Forcing function expresses the sea level caused by the"

   inr=inr+1
   out%readme(inr)="external forcing only but not by the sea level itself"
   inr=inr+1   
!   out%readme(inr)="F_00 denotes the eustatic sea level change (mass taken from ocean)"
   out%readme(inr)="F_00 denotes MASS TAKEN FROM THE OCEAN per steradian"
   inr=inr+1
   out%readme(inr)="Ocean file used:"
   inr=inr+1
   out%readme(inr)=oce%file(1:80)
   inr=inr+1
   out%readme(inr)="load love set from the model are put in vectordata (hnm,lnm,knm)"
   

   if(rotfeedback)then
      inr=inr+1
      out%readme(inr)=&
     "Constructed using linearized rotational feedback of the sea level"
   end if
   ii=ii+1
   out%ints_d(ii)='Lmax'
   out%ints(ii)=lmax
   ii=ii+1
   out%ints_d(ii)='Lmin'
   out%ints(ii)=0
   
   idb=idb+1
   out%dbls_d(idb)='Rho_w'
   out%dbls(idb)=rho_w
   idb=idb+1
   out%dbls_d(idb)='Rho_e'
   out%dbls(idb)=rho_e

      
   if(rotfeedback)then
      idb=idb+1
      out%dbls_d(idb)="k2_body"
      out%dbls(idb)=knm_T(SH_pos(2,0))
      idb=idb+1
      out%dbls_d(idb)="h2_body"
      out%dbls(idb)=hnm_T(SH_pos(2,0))
   end if

end if








!second side description
!adapt original side description data

allocate(oce%side2_d(oce%nval2)) ! allocate the second dside description vector

if(regmat)then
   do,i=1,oce%nval1
      oce%side1_d(i)(1:1)='T' ! change G in T (Land surface load)
      oce%side2_d(i)=oce%side1_d(i)
   end do

else

   do,i=1,oce%nval1
      oce%side1_d(i)(1:1)='T' ! change G in T (equivalent water height)
      oce%side1_d(i)(17:)='SEALEVEL'
      oce%side2_d(i)=oce%side1_d(i)
      oce%side2_d(i)(16:24)='SEAQUASI'
   end do

end if
out%side1_d=>oce%side1_d

if(invert)then
    out%side2_d=>oce%side1_d
    out%side1_d=>oce%side2_d
else
    out%side2_d=>oce%side2_d
    out%side1_d=>oce%side1_d
end if

!sort and put love numbers in vector
out%vec(1,:)=1/1.d-99 !( Put  for the degree zero love numbers)
do,i=2,oce%nval1
   pos=SH_pos(deg(i),ord(i))
!   write(0,*)pos,deg(i),ord(i)
   out%vec(i,1)=hnm(pos)
   out%vec(i,2)=lnm(pos)
   out%vec(i,3)=knm(pos)
end do


!calculate changes due to rotational feedback
if(rotfeedback)then
   !construct matrix relating surface loading changes to 'apparent' potential changes
   !based on linearized conservation of momentum
   Ji3tmp=Ji3_2_m(ltyp)
   RotMat=matmul(m_2_rotpot(),matmul(Ji3tmp,SurfL_2_Ji3()))
   allocate(Rotrows(4,out%nval1)) ! array to temporary store the effect of rotational feedback


   !create index vector (length 4, c00,c20,c21,s21)
   pos=0
   i=0
   do while(pos< 4 .and. i<oce%nval1)
      i=i+1
      select case(deg(i))
      case(0)!degree and order 0
         rotind(1)=i
         pos=pos+1
      case(2)
         select case(ord(i))
         case(0)!degree 2 order 0
            rotind(2)=i
            pos=pos+1
         case(1) !degree 2 order 1
            if(trig(i) .eq. 0)then ! cosine
               rotind(3)=i
               pos=pos+1
            else
               rotind(4)=i
               pos=pos+1
            end if
         end select
      end select
   end do

   !Copy ocean matrix values in RotRows
   do,j=1,out%nval1
      do,i=1,4
         Rotrows(i,j)=oce%mat1(rotind(i),j)
      end do
   end do
   
   !calculate matrix multiplication
   Rotrows=matmul(Rotmat,Rotrows)
   Rotrows=Rotrows*gamma2

end if


!adapt oce%mat1 step 1 EARTH MODEL dependent part
!convolve non zero degree part with G_N-U greens function ( geoid minus uplift)
do,i=1,out%nval1 !loop over columns
   do,j=2,out%nval1 ! loop over rows ( skip degree 0)
      !calculate row dependent factor
      fact=-(3.d0*rho_w/rho_e)*(1.d0+out%vec(j,3)-out%vec(j,1))/(2*deg(j)+1)    
      out%mat1(j,i)=fact*out%mat1(j,i)
   end do
   !simply negate the first row
   out%mat1(1,i)=-out%mat1(1,i)
end do


if(rotfeedback)then !add the matrix from rotational feedback
   !update matrix with rotational feedback part
   do,j=1,out%nval1
      do,i=2,4 !skip degree 0 to ensure mass consistency ( is small anyway)
         out%mat1(rotind(i),j)=out%mat1(rotind(i),j)+Rotrows(i,j)
      end do
   end do
   

end if


!degree 0 row has a special treatment ( divide by O_00 to retrieve eustatic sea level)
! do, i=1,oce%nval1
!    out%mat1(1,i)=out%mat1(1,i)/A_00
! end do

!add unit diagonal
do,i=2,oce%nval1
   out%mat1(i,i)=1+out%mat1(i,i)
end do



!!invert if requested ( LU decomposition)
if(invert)then
   call cpu_time(time1)
   allocate(ipiv(out%nval1))

   !calculate matrix L1 norm
   anorm=0.d0
   do,i=1,out%nval2 ! loop over columns
      dumdb=0.d0
      do,j=1,out%nval1 !loop over rows
         dumdb=dumdb+abs(out%mat1(j,i))
      end do
      anorm=max(anorm,dumdb)
   end do
   
   allocate(workv(4*out%nval2))
   allocate(iwork(out%nval2))

   !first LU factorization
   call dgetrf(out%nval1,out%nval1,out%mat1,out%nval1,ipiv,info)

   if( info .ne. 0)then
      write(stderr,*)"ERROR: factorization failed"
      stop
   end if
   rcond=0.d0


!    ! calculate condition number
    call dgecon('1',out%nval2,out%mat1(1,1),out%nval1,anorm,&
               rcond,workv,iwork,info)

   if( info .ne. 0)then
      write(stderr,*)"ERROR: calcuating condition number"
      stop
   end if

   idb=idb+1
   out%dbls_d(idb)='Cond_No_L1'
   out%dbls(idb)=1/rcond

   !invert 
   !optimal lwork:

   lwork=out%nval1*ilaenv( 1, 'DGETRI', ' ', out%nval1, -1, -1, -1 )

   allocate(work(lwork))
   call dgetri(out%nval1,out%mat1,out%nval1,ipiv,work,lwork,info)
   if( info .ne. 0)then
      write(stderr,*)"ERROR: inversion failed"
      stop
   end if
   call cpu_time(time2)
   inr=inr+1
   write(out%readme(inr),*)"Matrix inversion took",time2-time1,"seconds" 


end if

if(regmat)then !apply post processing to obtain a regularization matrix
   !When the above matrix is called G apply
   ! M =  O (I - G * D* (I-O) )  
   ! and perform a rank k update :
   !Mout= M' * M


   if(rotfeedback)then
      deallocate(Rotrows)
      allocate(Rotrows(out%nval1,4))
      
      !Copy matrix values in RotRows
      do,i=1,4
         do,j=1,out%nval1
            Rotrows(j,i)=out%mat1(j,rotind(i))
         end do
      end do
      
      !special treatment of degree 0
      Rotmat(:,1)=Rotmat(:,1) !-Rotmat(:,1)/A_00
      !calculate right side matrix multiplication
      
      Rotrows=matmul(Rotrows,Rotmat)
      Rotrows=Rotrows*gamma2

   end if
   
   ! right multiply with a diagonal matrix
   !convolve non zero degree part with G_N-U greens function ( geoid minus uplift)
   do,i=1,out%nval1 !loop over columns
      if(deg(i)==0)then !special case for degree 0
         fact=1.d0 !-1/A_00
      else
         fact=(3.d0*rho_w/rho_e)*(1.d0+out%vec(i,3)-out%vec(i,1))/(2*deg(i)+1)    
      end if

      do,j=1,out%nval1 ! loop over rows 
         out%mat1(j,i)=fact*out%mat1(j,i)
      end do

   end do

   if(rotfeedback)then !add the matrix from rotational feedback
      !update matrix with rotational feedback part
      do,i=1,4
         do,j=1,out%nval1
            out%mat1(j,rotind(i))=out%mat1(j,rotind(i))+Rotrows(j,i)
         end do
      end do
   

   end if


   ! right multiplication with -(I-O) = O-I
   
   ! allocate matrix and copy values
   allocate(work2d(out%nval1,out%nval1))
   work2d=out%mat1 !copy values
   
   call dsymm('R',oce%mtyp,out%nval1,out%nval2,1.d0,oce%mat1(1,1),oce%nval1,out%mat1(1,1),out%nval1,-1.d0,work2d(1,1),out%nval1)
   
   ! add unit diagonal
   do ,i=1, out%nval1
      work2d(i,i)=work2d(i,i)+1.d0
   end do

   !left multiply by the ocean matrix
   call dsymm('L',oce%mtyp,out%nval1,out%nval2,1.d0,oce%mat1(1,1),oce%nval1,work2d(1,1),out%nval1,0.d0,out%mat1,out%nval1)
   
   if(plaw)then !apply degree dependent error-covariance
      do,i=1,out%nval1
         out%mat1(i,1:out%nval1)=out%mat1(i,1:out%nval1)/sqrt(sc*deg(i)**pn)
      end do
   end if

   if(regusef)then
      stc=1
      ndc=pos
      sts=pos+1
      nds=2*pos
      allocate(clmtmp(2*pos),clm_sig(2*pos))
      call SH_readgrav(filen=trim(regfile),clm=clmtmp(stc:ndc),&
           &slm=clmtmp(sts:nds),clm_sig=clm_sig(stc:ndc),&
           slm_sig=clm_sig(sts:nds))
      do,i=1,out%nval1
         ind=SH_pos(deg(i),ord(i))
         out%mat1(i,1:out%nval1)=out%mat1(i,1:out%nval1)/(clm_sig(ind+trig(i)*pos))
      end do

   end if


!   DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
   ! rank k update
    call dsyrk('U','T',out%nval1,out%nval1,1.d0,out%mat1(1,1),out%nval1,0.d0,work2d(1,1),out%nval1)

    deallocate(out%mat1)
    out%mat1=>work2d !point to updated matrix ( warning may yield a dangling pointer, when out%mat1 is not deallocated)

    !reset vec part to zero
    out%vec=0.d0
    out%nvec=2

 end if



!write to file

call write_BINtype(out)


end program SH_sealevmat


subroutine help()
implicit none
integer::unit=6
character(8)::frmt='(A)'

write(unit,frmt)"Program SH_sealevmat constructs a matrix used for solving the sea level equation"
write(unit,frmt)"The matrix M1 or M2 is written to standard output, and holds for"
write(unit,frmt)" M1*S'= L"
write(unit,frmt)" M2*L= S'  (M2=M1^(-1))"
write(unit,frmt)" Where M1 is constructed from the sea level equation:"
write(unit,frmt)" S'= N-U + DV/g (over the ocean only)"
write(unit,frmt)" N is the geoid, U is the ocean floor uplift"
write(unit,frmt)" The term S'_00 = DV/g is solved to be mass consistent"
write(unit,frmt)" The actual sea level used for loading can be retrieved by:"
write(unit,frmt)" S=O*S', where O is the ocean convolution matrix (input for the program)"
write(unit,frmt)" L is the loading function expressed in relative sea level."
write(unit,frmt)" It describes the effect of the external loading only on the sea level"
write(unit,frmt)" L_00 is the mass per steradian taken out of the Ocean"
!write(unit,frmt)" L_00/O_00 is the eustatic sea level change (mass added to the ocean by the loading)"
write(unit,frmt)" For surface loading this would involve the load love numbers, hn',kn',ln' while, for tidal "
write(unit,frmt)" loading body love numbers, hn,kn,ln, would be required"
write(unit,frmt)""
write(unit,frmt)"Usage SH_sealevmat [OPTIONS] OCEANMAT"
write(unit,frmt)" Where OCEANMAT is the file which holds the ocean convoluton matrix"
write(unit,frmt)"Where OPTIONS might be:"
write(unit,frmt)" -Ef=LOADLOVEFILE"
write(unit,frmt)" -E[NUM]: predefined load love numbers;"
write(unit,frmt)"   NUM=2:PREM(default) , NUM=1:Gutenberg-Bullen(farrel 1972)"
write(unit,frmt)" -n: Output matrix M1 ( don't invert)"
write(unit,frmt)" -r: Rotational feedback. Allow sea level to change the instantaneous"
write(unit,frmt)"     rotation axis of the Earth. (apply linearized conservation of angular momentum)"
write(unit,frmt)"-R[SCALE,POWER]: construct a regularization Matrix., Which reduces the difference  between"
write(unit,frmt)"   The surface load and the mass consistent equilibrium approach to the continental load"
write(unit,frmt)"   A power law (with SCALE n^POWER) can be optionally assumed for the covariance of this difference."
write(unit,frmt)" -Rf=SHFILE: SAme as -R but uses the diagonal errors of SHFILE as weights."

write(unit,frmt)"   NOTE: the -R option ignores the -n option"


stop
end subroutine help

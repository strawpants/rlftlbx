!!Coded by Roelof Rietbroek, Tue Jul 28 15:33:31 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de


!!program converts isotropic spherical harmonics (degree dependent only) to spatial domain (colatitude)

!!Updated by Roelof Rietbroek, Thu Nov 29 11:19:25 2012
!! added possibility for calculating the derivative of the Legendre polynomial
!! also added the computation of surface loading Greens functions
!!Updated by Roelof Rietbroek, Wed Dec  5 17:29:16 2012
!!Set default to unnormalized Legendre polynomials and add option for adding a normalization
!!Updated by Roelof Rietbroek, Thu Dec  6 14:43:35 2012
!! added geoid greensfunction
!!Updated by Roelof Rietbroek, Mon Dec  8 16:42:06 2014
!!fixed bug in dfault output domain (outputted to 0,pi degrees instread of 0 180 degrees) 
!! Updated 18 April 2017 fixd bug (rho_w should have bveen rho_e)

program SHiso_2_colat
use shtools
use shtlbx
implicit none
integer::narg,iargc,l,i,itharg,lmax,nmax
integer::lmax_f,lmin
integer::unit,stderr,last,ind
character(200)::dum,file,lovefile
logical::stdin
parameter(nmax=100000)!maximum amount of the degree supported
double precision::Cl(nmax+1),res,rmin,rmax,d2r
double precision,allocatable::Pl(:),dPl(:),loadnm(:),colat(:)
double precision,allocatable::Fout(:,:)
character(5)::type
logical::asymp,limdeg,sinmult
integer::ltyp,FR,ncolat,ncol
double precision::loadinfty,fact,sintmp,tmplogsqrt
logical::norm,logx

!defaults initializations
rmin=0.d0
rmax=180
d2r=pi/180.d0
stdin=.true.
stderr=0
unit=5
lmax=100
ncol=1
lmax_f=0
limdeg=.false.
lmin=0
ltyp=2 !load lovenumber typ (2 is PREM Earth model)
FR=3 ! center of figure is default
lovefile=''
loadinfty=0.d0
sinmult=.false.
type='PLAIN'
asymp=.false.
!default output resolution
res=0.1*d2r
Cl=0.d0
norm=.false.
logx=.false.


!!process command line options
narg=iargc()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('d')!set resolution
         if(dum(3:3) .eq. 'l')then
            read(dum(5:),*)res
            logx=.true.
         else
            read(dum(4:),*)res
         end if
      case('r')!limit output region
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min range is specified
            read(dum(4:ind-1),*)rmin
            read(dum(ind+1:),*)rmax
         else
            read(dum(4:),*)rmax
         end if
      case('l')!limit maximum (and minum degree)
         limdeg=.true.
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
      case('D')! Compute the derivative
         type='DERIV'
      case('G')
         select case(dum(3:3))
         case('u')
            type='GRN_U'
         case('U')
            type='GRN_U'
            sinmult=.true.
         case('v')
            type='GRN_V'
         case('V')
            type='GRN_V'
            sinmult=.true.
         case('n')
            type='GRN_N'
         case('N')
            type='GRN_N'
            sinmult=.true.
         case default
            write(stderr,*)'ERROR: processing option:',trim(dum)
            stop
         end select
         if(dum(4:4) .eq. 'a') asymp=.true. ! include asymtpotic beahvior of the Love numbers
         if(asymp)ncol=2 !output two columns
      case('E') !Use different load love numbers (different Earth model)
         read(dum(3:),'(i1)')ltyp
         if(ltyp .eq. 0)then
            if(dum(4:).eq. ' ')then
               itharg=itharg+1
               call getarg(itharg,lovefile)
            else
               lovefile=dum(4:)
            end if
         end if
      case('O')!set a different origin of the reference frame
         select case(dum(3:4))
         case('CE')
            FR=1
         case('CM')
            FR=2
         case('CF')
            FR=3
         case('CL')
            FR=4
         case('CH')
            FR=5
         case default
            write(stderr,*)'unknown reference frame origin specified, quitting'
            stop
         end select
      case('N') ! use normalized polynomials
         norm=.true.
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is a file
      stdin=.false.
      unit=13
      file=trim(dum)
   end if
end do

! set up colatitude variations
if(logx)then
   if(rmin < tiny(1.d0))then
      write(stderr,*)'ERROR: a logarithmic X -axis requires rmin > 0 '
      stop
   end if
   ncolat=int((log10(rmax)-log10(rmin))/res)+1
   allocate(colat(ncolat))
   do,i=1,ncolat
      colat(i)=d2r*10**(log10(rmin)+res*(i-1))
   end do
else
   rmax=rmax*d2r
   rmin=rmin*d2r
   res=res*d2r
   ncolat=int((rmax-rmin)/res)+1
   allocate(colat(ncolat))
   do,i=1,ncolat
      colat(i)=rmin+res*(i-1)
   end do
end if

!read in the data
select case(type)
case('DERIV','PLAIN')
   if(.not. stdin)open(file=trim(file),unit=unit,form='formatted')
   last=0
   do while (last .eq. 0)
      read(unit=unit,iostat=last,fmt=*)l,Cl(l+1)
!      write(0,*)l,Cl(l+1)
      lmax_f=max(l,lmax_f)
   end do
   if(.not. stdin)close(unit)
   if(.not. limdeg)lmax=lmax_f

   if(lmin>0)Cl(1:lmin)=0.d0 ! set lower coefficients to zero

!   write(0,*)lmin,lmax
case('GRN_U') ! needs to read in load love numbers
   if(.not. limdeg)then
      write(stderr,*)'ERROR: -G option requires -l option'
      stop
   end if
   allocate(loadnm(lmax+1))
   call SH_loadlove(hnm=loadnm,typ=ltyp,frame=FR,fileop=trim(lovefile),iso=.true.)

   if(asymp)loadinfty=loadnm(lmax+1) !set 'asymptotic' load love number 

   fact=3.d0/(4*pi*rho_e*(RE**2)) ! = a/me
   do, l=lmin,lmax
      Cl(l+1)=fact*(loadnm(l+1)-loadinfty)
   end do

   if(norm)then
      do, l=lmin,lmax
         Cl(l+1)=Cl(l+1)/sqrt(2.d0*l+1.d0)
      end do
   end if


case('GRN_N') ! needs to read in load love numbers
   if(.not. limdeg)then
      write(stderr,*)'ERROR: -G option requires -l option'
      stop
   end if
   allocate(loadnm(lmax+1))
   call SH_loadlove(knm=loadnm,typ=ltyp,frame=FR,fileop=trim(lovefile),iso=.true.)
   

   if(asymp)then
      loadinfty=lmax*loadnm(lmax+1) !set 'asymptotic' load love number 
      lmin=max(lmin,1) !possibly reset lmin ( may not be zero)
   end if

   fact=3.d0/(4*pi*rho_e*(RE**2)) ! = a/me
   do, l=lmin,lmax

      
      if(asymp)then
         Cl(l+1)=fact*(loadnm(l+1)-loadinfty/dble(l))
      else
         Cl(l+1)=fact*(loadnm(l+1)+1.d0)
      end if
   end do

   if(norm)then
      do, l=lmin,lmax
         Cl(l+1)=Cl(l+1)/sqrt(2.d0*l+1.d0)
      end do
   end if

   
case('GRN_V')
   if(.not. limdeg)then
      write(stderr,*)'ERROR: -G option requires -l option'
      stop
   end if
   allocate(loadnm(lmax+1))
   call SH_loadlove(lnm=loadnm,typ=ltyp,frame=FR,fileop=trim(lovefile),iso=.true.)

   if(asymp)loadinfty=lmax*loadnm(lmax+1) !set 'asymptotic' load love number 

   fact=3.d0/(4*pi*rho_e*(RE**2)) ! = a/me
   do, l=max(lmin,1),lmax
      Cl(l+1)=fact*(dble(l)*loadnm(l+1)-loadinfty)/dble(l)
   end do

   if(norm)then
      do, l=lmin,lmax
         Cl(l+1)=Cl(l+1)/sqrt(2.d0*l+1.d0)
      end do
   end if

end select

!convert to colatitude 

allocate(Pl(lmax+1))
allocate(Fout(ncolat,ncol))
Fout=0.d0

select case(type)
case('PLAIN')
   do,i=1,ncolat
      if(norm)then
         call PlBar(Pl,lmax,cos(colat(i)))
      else
         call PLegendre(Pl,lmax,cos(colat(i)))
      end if
      Fout(i,1)=dot_product(Pl(1:lmax+1),Cl(1:lmax+1))
   end do
case('DERIV')
   allocate(dPl(lmax+1))
   do,i=1,ncolat
      if(norm)then
         call PlBar_d1(Pl,dPl,lmax,cos(colat(i)))
      else
         call PLegendre_d1(Pl,dPl,lmax,cos(colat(i)))
      end if
      ! apply chain rule to get derivative 
      dPl=-dPl*sin(colat(i))
      Fout(i,1)=dot_product(dPl(2:lmax+1),Cl(2:lmax+1))
   end do
case('GRN_U')
   sintmp=1.d0
   do,i=1,ncolat
      if(norm)then
         call PlBar(Pl,lmax,cos(colat(i)))
      else
         call PLegendre(Pl,lmax,cos(colat(i)))
      end if

      if(sinmult)sintmp=sin(colat(i)/2.d0)
      Fout(i,1)=sintmp*dot_product(Pl(1:lmax+1),Cl(1:lmax+1))
   end do

   if(asymp)then
      sintmp=1.d0
      do,i=1,ncolat
         if( .not. sinmult)sintmp=sin(colat(i)/2.d0)
         Fout(i,2)=fact*loadinfty/(2*sintmp)
      end do
   end if

case('GRN_N')
   sintmp=1.d0
   do,i=1,ncolat
      if(norm)then
         call PlBar(Pl,lmax,cos(colat(i)))
      else
         call PLegendre(Pl,lmax,cos(colat(i)))
      end if

      if(sinmult)sintmp=sin(colat(i)/2.d0)
      Fout(i,1)=sintmp*dot_product(Pl(1:lmax+1),Cl(1:lmax+1))
   end do

   if(asymp)then
      sintmp=1.d0
      do,i=1,ncolat
         sintmp=sin(colat(i)/2.d0)
         tmplogsqrt=sqrt(1-cos(colat(i)))
         tmplogsqrt=log(tmplogsqrt*(tmplogsqrt+sqrt(2.d0))/2.d0)
         
         if(colat(i) <= tiny(colat(i)) .and. sinmult)tmplogsqrt=0.d0 ! to fix limit theta -> 0 ( prevents nans)

         if(sinmult)then
            Fout(i,2)=fact*(1.d0/(2.d0)-sintmp*loadinfty*tmplogsqrt)
         else
            Fout(i,2)=fact*(1.d0/(2.d0*sintmp)-loadinfty*tmplogsqrt)
         end if
      end do
   end if

case('GRN_V')
   allocate(dPl(lmax+1))
   sintmp=1.d0
   do,i=1,ncolat
      if(norm)then
         call PlBar_d1(Pl,dPl,lmax,cos(colat(i)))
      else
         call PLegendre_d1(Pl,dPl,lmax,cos(colat(i)))
      end if

      ! apply chain rule to get derivative 
      dPl=-dPl*sin(colat(i))

      if(sinmult)sintmp=sin(colat(i)/2.d0)
      Fout(i,1)=sintmp*dot_product(dPl(2:lmax+1),Cl(2:lmax+1))
   end do

   if(asymp)then
      do,i=1,ncolat
         sintmp=sin(colat(i)/2.d0)
         if(sinmult)then
            Fout(i,2)=-fact*loadinfty*(cos(colat(i)/2.d0)/(2.d0)*((1.d0+2*sintmp)/(1.d0+sintmp)))
         else
            Fout(i,2)=-fact*loadinfty*(cos(colat(i)/2.d0)/(2.d0*sintmp)*((1.d0+2*sintmp)/(1.d0+sintmp)))
         end if
      end do
   end if


end select


! write output
unit=6 ! set to standard output
do,i=1,ncolat
   write(unit,*)colat(i)/d2r,Fout(i,1:ncol)
end do

end program SHiso_2_colat


subroutine help()
implicit none
integer::stderr
character(8)::frmt
stderr=0
frmt='(A)'
write(stderr,frmt)"SHiso_2_colat converts degree wise SH coefficients"
write(stderr,frmt)" into the spatial domain(colatitude)"
write(stderr,frmt)"usage: SHiso_2_colat [OPTIONS] [FILE]"
write(stderr,frmt)"Where FILE is a 2 column ascii file with degree and value"
write(stderr,frmt)"If FILE is not given it will read from standard input"
write(stderr,frmt)"And OPTIONS are:"
write(stderr,frmt)"-d[l]=RES: set resolution of output in degrees (default is 0.1 degree)"
write(stderr,frmt)"  If 'l' is provided nodes on a logarithmic (base 10) X-axis are outputted."
write(stderr,frmt)"-r=RMIN,RMAX: limit output region in degrees: default is 0- 180"
write(stderr,frmt)"-l=LMAX[,LMIN]: limit degrees to LMAX and optionally LMIN"
write(stderr,frmt)"-D: compute the derivative in colatitude direction"
write(stderr,frmt)"-G[uvnUVN][a]: Compute the surface loading Greens function"
write(stderr,frmt)"           specify 'u' for the uplift, 'v' for the horizontal deformation"
write(stderr,frmt)"           or 'n' for geoid change."
write(stderr,frmt)"           when 'a' is appended an asymptotic behavior is assumed for l>lmax"
write(stderr,frmt)"           (see Farrell 1972)"
write(stderr,frmt)"           The capital U,V and N multiply the Greensfunction with sin(theta/2)"
write(stderr,frmt)"           This prevents asymptotic beahvior for theta -> 0)"
write(stderr,frmt)"-ENUM:  Use the following choices of Earth models for the Love numbers NUM is:"
write(stderr,frmt)"    1: Farrells load love numbers (GUTENBERG-BULLEN Earth model)"
write(stderr,frmt)"    2: Use PREM load love numbers (PREM Earth model) (default)"
write(stderr,frmt)"    3: Use Guttenberg-Bullen BODY love numbers (degree 2-4 only)"
write(stderr,frmt)"    4: Use PREM BODY love numbers (degree 2-4 only)"
write(stderr,frmt)"    5: Use Desai 2002  BODY love numbers (degree 2 only)"
write(stderr,frmt)"    6: Use PREM load Love umbers from Wang et al 2012"
write(stderr,frmt)"    7: Use ak135 load Love umbers from Wang et al 2012"
write(stderr,frmt)"    8: Use iasp91 load Love umbers from Wang et al 2012"
write(stderr,frmt)"    0CUSTOMFILE: use a custom file (same format as the above models"
write(stderr,frmt)"-OISOFRAME: Set the Origin of reference frame to ISOFRAME"
write(stderr,frmt)"            Applies to -g and -s ( degree 1)"
write(stderr,frmt)"            CM: Center of mass of the Earth system"
write(stderr,frmt)"            CF: Center of surface figure (default)"
write(stderr,frmt)"            CE: Center of Mass of the solid Earth"
write(stderr,frmt)"            CL: Center of surface Lateral figure"
write(stderr,frmt)"            CH: Center of surface Height figure"
write(stderr,frmt)"-N: use normalized Legendre polynomials"
write(stderr,frmt)""
write(stderr,frmt)"Output is a 2 column ascii file with colatitude (degrees) and value"
write(stderr,frmt)"The options -D and -G cannot be combined"
stop
end subroutine help

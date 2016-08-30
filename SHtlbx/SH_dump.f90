!!Program to read in SH files and dump to screen in various formats (ICGEM)
!performs on the fly conversions and adaptions if requested
!!Coded by Roelof Rietbroek, Thu May 31 13:36:11 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue Jul  3 15:20:28 2007
!!incorporated error propagation and -s and -a options now also accept the next argument (after a space) as their filenames
!!so both -sFILE as -s FILE are allowed

!!Updated by Roelof Rietbroek, Tue Sep 18 18:10:03 2007
!! separate (calibrated) error files are now allowed SH_read also accepts GRACE calibrated error files (keyword is CALSDV)

!!Updated by Roelof Rietbroek, Fri May 23 10:41:36 2008
!!One can choose between two different Earth models now (PREM and GUTENBERG-BULLEN)

!!Updated by Roelof Rietbroek, Tue Jul  8 15:37:51 2008
!! one can also supply a custom degree wise scale file to be applied

!!Updated by Roelof Rietbroek, Wed Oct 15 17:12:35 2008
!! fixed bug geoid scale can now also be used with inv

!!Updated by Roelof Rietbroek, Fri Jan 30 11:40:27 2009
!! added static sea level in the functional
!! avoid use of g but rather use GM/RE^2
!! change name of some logicals
!! changed help function
!! added forcing potential mode ( sets kn love number to zero)
!! added Earth models ( also tidal forcing)

!!Updated by Roelof Rietbroek, Thu Apr  9 11:19:46 2009
!!incorporated adding of DOT terms from CSR/GFZ files

!!Updated by Roelof Rietbroek, Wed Sep 16 10:25:45 2009
!!added optional negative sign for constructing direct sealevel
!!Updated by Roelof Rietbroek, Mon Sep 21 11:47:30 2009
!! added unit input: this enables pure scaling factors to be outputted

!!Updated by Roelof Rietbroek, Fri Oct  9 16:30:16 2009
!!added support for conversion between loading potential and loading+solid Earth potential

!!added support for new annual geocenter model ( from JGR paper)

!!Updated by Roelof Rietbroek, Fri Dec 18 10:56:56 2009
!! added support for rotating Spherical harmonics under three Euler angles (zy'z')

!!Updated by Roelof Rietbroek, Thu Jan 21 11:39:20 2010
!!added option for standard centrifugal potential

!!Updated by Roelof Rietbroek, Fri Feb  5 17:57:10 2010
!! added tag for conservation of angular momentum (forced by loading potential)

!!Updated by Roelof Rietbroek, Thu Feb 11 12:09:16 2010
!! added time tag forcing

!!Updated by Roelof Rietbroek, Wed Feb 17 11:09:36 2010
!! Added horizontal deformation part (not yet differentiated to direction)

!!Updated by Roelof Rietbroek, Wed Apr 21 08:36:17 2010
!!Added conversion to microgal

!!Updated by Roelof Rietbroek, Wed Sep  8 17:30:31 2010
!! allow choice of reference frame

!!Updated by Roelof Rietbroek, Tue Apr 26 10:50:35 2011
!! allow a c20 replacement from the Ries or JIGOG series
!!Updated by Roelof Rietbroek, Wed Apr 27 14:31:38 2011
!! changed the application of c20 replacement before -s -a option
!!Updated by Roelof Rietbroek, Wed Apr 27 15:06:12 2011
!!added new annual geocenter motion model

!!Updated by Roelof Rietbroek, Fri May 13 15:21:30 2011
!!added gravity anomalies

!!Updated by Roelof Rietbroek, Mon May 16 15:39:57 2011
!fixed normalization (also apply to slm and fixed factorial)

!!Updated by Roelof Rietbroek, Thu Sep  1 15:58:43 2011
!!fixed bug in processing -m option

!!Updated by Roelof Rietbroek, Mon Apr 16 11:10:11 2012
!! also allow restoring of the timevariable coefficients from icgem files

!!Updated by Roelof Rietbroek, Wed Aug 29 17:06:31 2012
!! set order zero coefficients to zero when -u is selected

!!Updated by Roelof Rietbroek, Fri Feb 22 15:09:01 2013
!! added automatic retrieval of EOP parameters

!!Updated by Roelof Rietbroek, Tue Mar  5 21:48:35 2013
!!enabled integration over time of the harmonic

!!Updated by Roelof Rietbroek, Mon Mar 18 17:18:55 2013
!! added option to the -m5 option to avoid the division by the Ocean/Land ratio

!!Updated by Roelof Rietbroek, Wed Nov 20 11:15:11 2013
!!modified help output of the -a and -s option to be consistent

!!Updated by Roelof Rietbroek, Fri Dec 20 15:16:30 2013
!! added option to keep or sibtract the iers 2010 mean pole when computing the 
! effect of the Earth rotation




program SH_dump
use SHtlbx
use SHTOOLS
use EOP_tools
implicit none
integer::narg,i,lmax,lmin,pos,typ,maxf,add,subt,cmod,otyp,ind,st,nd,ftyp,gtyp,FR,etyp
integer::ltyp,itharg,l,m,unit,ms,stderr
parameter(maxf=50)
character(200)::filen
character(200)::dum,errfile,geocfile,scalefile,lovefile
character(len=200),dimension(maxf)::addfile,subtfile
double precision,allocatable,dimension(:)::clm,slm,vlm,hnm,knm,lnm,clm_add,slm_add,clm_subt,slm_subt
double precision,allocatable,dimension(:)::clm_sig,slm_sig
double precision::radius,tcent,mu,tstart,tend,fact,tmpfac,c10,c11,s11,pfac,norm,tmp
logical:: geoid,water,geowater,smooth,smoothrecur,geoc,filter,inv,limdeg,uplift,stdin,ptide,error,compl,custom
logical::normalize,customdegree,rigid,sealev,dot,harm
logical::loadNdirect
integer::lmaxdot,posdot
double precision::dotref,c00sc
double precision,allocatable,dimension(:,:)::clm_tv,slm_tv
integer::iargc
logical::unitin,rotate,rotpot,polarmotion
logical::horiz,cons_angmom,cons_anginfo,forcetime
integer::ind2
double precision::angles(3),sq2pi,sq4pi,mpolar(3)
double precision,allocatable::dj(:,:,:),Csh(:,:),Csh_rot(:,:)
double precision:: OHM_R
double precision::Rmat(4,4)
double precision::d2coef(4) !temporary vector to hold rotational feedback coef
double precision::cent_pot(2)!degree 00 and 20 coefficients of perfectly aligned rotation axis
double precision::itcent,itstart,itend
double precision::tmat(5)
logical::radialacc
logical::gravanom
integer::c20type
integer::idum,ierr
logical::geocJ2
double precision::c30,c31,s31
double precision::mjd(1:1)
double precision::Ji3tmp(3,3) ! needed for buggy gfortan?

!initialize parameters
unitin=.false.
c00sc=1.d0
pfac=0.d0
fact=1.d0
add=0
subt=0
cmod=0
ftyp=0
otyp=4
lmax=0
lmin=0
ltyp=2 !load lovenumber typ (2 is PREM
FR=3 !frame type for degree 1 load love numbers (3 means center of figure, 1 center of mass)
itharg=0
ms=0 !integer to check for ambiguous command line arguments
stderr=0 !standard error unit
polarmotion=.false.
mpolar=0.d0 ! set polar motion angles to 0
c20type=0
geocJ2=.false.

!set all logicals to false
error=.false.

geoid=.false.
water=.false.
geowater=.false.
sealev=.false.
uplift=.false.
LoadNdirect=.false.
smooth=.false.
smoothrecur=.false.
geoc=.false.
filter=.false.
inv=.false.
limdeg=.false.
cons_angmom=.false.
cons_anginfo=.false.
horiz=.false.
radialacc=.false.
gravanom=.false.
stdin=.true.
ptide=.false.
errfile=''
compl=.false.
custom=.false.
customdegree=.false.
normalize=.false.
tstart=0.d0
tend=0.d0
filen=''
geocfile=''
scalefile=''
lovefile=''
gtyp=-9999
rigid=.false.
rotate=.false.
rotpot=.false.
forcetime=.false.
!add rates defaults
dot=.false. ! add dot terms
harm=.false. ! add harmonic terms (annual semi annual)
lmaxdot=4
dotref=2000. !refernce in year

!!set up processing strategy from command line
!get number of command line arguments
narg=iargc()
if(narg .eq. 0)call help()



!then make a loop over other arguments
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit! exit loop when last argument has been read in
   call getarg(itharg,dum)
   if(dum(1:1) .ne. '-')then !argument is the filename
      filen=dum
      stdin=.false.
   else !when argument is an option
      select case(dum(2:2))
      case('m')!multiplication by a predefined factor (possibly degree dependent) to apply to the data
         select case(dum(3:3))
         case('s')!apply custom scale factor
            ms=ms+1
            custom=.true.
            read(dum(4:),*)fact
         case('d')!degree wise scaling from file
            ms=ms+1
            customdegree=.true.
            if(dum(4:).eq. ' ')then
               itharg=itharg+1
               call getarg(itharg,scalefile)
            else
               scalefile=trim(dum(4:))
            end if
!            write(*,*)scalefile
         case default
            ind=len_trim(dum)
            if(dum(ind:ind) .ne. '-' .and. dum(ind:ind) .ne. 'i' .and. dum(ind:ind) .ne. 'o')then
               read(unit=dum(3:ind),fmt=*,iostat=ierr)idum
            else
               read(unit=dum(3:ind-1),fmt=*,iostat=ierr)idum
            end if

            if(ierr .ne. 0)then
               write(stderr,*)'ERROR: SH_dump:unknown option: ',dum(1:3)
               stop
            end if
            select case(idum)   
            case(1)!geoid height
               ms=ms+1
               geoid=.true.
               fact=RE !(= GM/(RE*gamma))
            case(2)!equivalent water height
               ms=ms+1
               water=.true.
               fact=RE*rho_e/(rho_w*3)
            case(3)!sea level height (no ocean loading effect, more consistent with altimetry)
               ms=ms+1
               geowater=.true.
               fact=RE*rho_e/(rho_w*3)
            case(4)!ocean floor deformation due to loading effect only
               ms=ms+1
               uplift=.true.
               fact=RE
            case(5)
               ms=ms+1
               sealev=.true.
               fact=RE
               if(dum(4:4) .eq. '-')then ! apply a negative scale and convert to eustatic sealevel rise
                  c00sc=-1.d0/A_oce
               else if(dum(4:4) .eq. 'o')then !just assume one for the degree 0 scale
                  c00sc=1.d0
               else
                  c00sc=1.d0/A_oce
               end if
            case(6)!convert direct to direct+loading potential
               ms=ms+1
               fact=1.d0
               loadNdirect=.true.
            case(7)!convert to changes in centrifugal potential
               ms=ms+1
               cons_angmom=.true.
               fact=RE*rho_e/(rho_w*3)
               if(dum(4:4).eq. 'i')cons_anginfo=.true. ! supply extra info on derived polar motion
            case(8)!horizontal deformation due to loading potential
               ms=ms+1
               fact=RE
               horiz=.true.
            case(9)!convert stokes coefficients to radial acceleration in microgal
               radialacc=.true.
               fact=GM*(1d8/(RE**2)) !conversion to microgal
            case(10)!convert stokes coefficients to gravity anomalies in microgal
               gravanom=.true.
               fact=GM*(1d8/(RE**2)) !conversion to microgal
               
            case default
               write(stderr,*)'ERROR: SH_dump:unknown option: ',dum(1:3)
               stop
            end select

         end select
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
      case('d') !add dot terms
         ind=index(dum,',')
         
         if(index(dum(3:3),'a') > 0)then
            harm=.true.
            if(ind >0)then
               read(dum(4:ind-1),*)lmaxdot
               read(dum(ind+1:),*)dotref
            end if
         else
            dot=.true.
            if(ind >0)then
               read(dum(3:ind-1),*)lmaxdot
               read(dum(ind+1:),*)dotref
            end if
         end if
    case('j')!replace C20 values
         read(dum(3:),*)c20type
      case('R') ! assume rigid Earth
         rigid=.true.
      case('w')!apply gaussian smoothing
         smooth=.true.
         if(dum(3:3) .eq. 'r')then
            smoothrecur=.true.
            read(dum(4:),*)radius
         else
            read(dum(3:),*)radius
         end if
!      case('f')!filter coefficients
!          write(*,*)'Filter options are not yet supported, check out the program SH_filter'
!          call help()
!          filter=.true.
!          read(dum(3:),*)ftyp !filter type (swenson filter, chambers filter, kusche filter)
      case('i')!inverse multiplication (for example to convert equivalent water height to geopotential coefficients  
         inv=.true.
      case('a')!add Spherical harmonic field (more fields are allowed max 5)
         add=add+1
         if(add > maxf)then
            write(stderr,*)"ERROR:maximum number of files which may be added",maxf
            stop
         end if
         if(dum(3:3) .eq. ' ')then !get next argument
            itharg=itharg+1
            call getarg(itharg,dum)
            addfile(add)=dum
         else!or when the filename is glued to the option tag
            addfile(add)=dum(3:)
         end if
      case('s')!subtract SH field (more fields are allowed max 5)
         subt=subt+1
         if(subt > maxf)then
            write(stderr,*)"ERROR:maximum number of files which may be subtracted",maxf
            stop
         end if
         if(dum(3:3) .eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,dum)
            subtfile(subt)=dum
         else
            subtfile(subt)=dum(3:)
         end if

      case('g','G')!applygeocenter motion
         geoc=.true.
         if(dum(2:2).eq. 'G')geocJ2=.true.
         read(dum(3:),'(i1)')cmod
         
         if(cmod .eq. 7)then !read geocenterfilename from command line
            if(dum(4:4) .eq. ' ')then !get next argument
               itharg=itharg+1
               call getarg(itharg,geocfile)
            else!or when the filename is glued to the option tag
               geocfile=dum(4:)
            end if
         end if
      case('t')!Specify output type of spherical harmonic set (standard GRACE,ICGEM,GINS,stripped)
         select case(dum(3:3))
         case('o')
            read(dum(4:),*)otyp  
         case('i')
            read(dum(4:),*)gtyp  
         end select

      case('l')!limit maximum and minimum degree
         limdeg=.true.
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(3:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(3:),*)lmax
         end if
      case('p')!add or subtract a permanent tide from the c2,0 geopotential coefficient (in other words convert between tide free and zero tide models)
         ptide=.true.
         if(dum(3:3).eq. 'a')then
            pfac=1.d0
         else if(dum(3:3) .eq. 's')then
            pfac=-1.d0
         else
            write(stderr,*)'ERROR SH_dump:specify -pa or ps to add or subtract the permanent tide'
            call help()
         end if
      case('e')!also calculate simple error propagation
         error=.true.
         !now check whether a file is appended (separate error file)

         if(dum(3:3) .eq. 'f')then !option -ef is given , thus filename expected
            
            if(dum(4:4).ne. ' ')then
               errfile=dum(4:)
            else
               
               itharg=itharg+1
               call getarg(itharg,errfile)
               
            end if
         end if
      case('r')!euler rotation of the underlyuing reference frame
         rotate=.true.
         if(polarmotion)then
            write(stderr,*)"ERROR: Specifying -Zp and -r together is not allowed"
            stop
         end if
         !find commas
         ind=index(dum,',')
         ind2=index(dum(ind+1:),',')+ind

         read(dum(3:ind-1),*)angles(1)
         read(dum(ind+1:ind2-1),*)angles(2)
         read(dum(ind2+1:),*)angles(3)
         !convert degrees to radians
         angles=angles*pi/180.d0
      case('u')!input is a unit vector
         unitin=.true.
      case('N') !apply normalization factor with i the inverse normalization
         normalize=.true.
      case('c')!take the complement of the vector
         compl=.true.
      case('Z')!default rotational potential (bulge is aligned to z axis)
         rotpot=.true.
         OHM_R=OHM_E ! initial magnitude of rotation vector

         if(dum(3:3) .eq. 'p' .or. dum(3:3) .eq. 'P' )then ! polar motion is supplied
            if(rotate)then
               write(stderr,*)"ERROR: Specifying -Zp and -r together is not allowed"
               stop
            end if
            polarmotion=.true.
            !find commas
            ind=index(dum,',')
            if(ind .eq. 0)then ! look for a date in MJD
               read(dum(4:),*)mjd(1:1)
               if(dum(3:3) .eq. 'p')then
                  call get_polarmotion(mjd(1:1),mpolar(1:1),mpolar(2:2),mpolar(3:3),.true.)
               else
                  call get_polarmotion(mjd(1:1),mpolar(1:1),mpolar(2:2),mpolar(3:3),.false.)
               end if
               !convert to arcseconds
               mpolar(1)=mpolar(1)*3600*180/pi
               mpolar(2)=mpolar(2)*3600*180/pi
            else ! read EOP data from command line option
               
               ind2=index(dum(ind+1:),',')+ind
               
               read(dum(4:ind-1),*)mpolar(1)
               read(dum(ind+1:ind2-1),*)mpolar(2)
               read(dum(ind2+1:),*)mpolar(3)

            end if

            !Set and calculatte rotation parameters
            rotate=.true.
            angles(1)=0 ! irrelevant (rotation about axis of symmetry)
            angles(2)=sqrt(mpolar(1)**2+mpolar(2)**2)*pi/(180*3600)
            angles(3)=-atan2(mpolar(1),mpolar(2))-pi/2
            
            !Calculate magnitude of the perturbed rotation vector

!            OHM_R=OHM_E-0.843994804d-9*mpolar(3) ! iers 2003 conventions mpolar(3) is excess length of day in seconds
            OHM_R=OHM_E*(1-mpolar(3)/86400.d0) ! iers 2003 conventions mpolar(3) is excess length of day in seconds
            !set up centrifugal potential scales for a perfectly aligned rotation axis

            !degree and order 0
            cent_pot(1)=RE**3/(GM*3.d0)
            !zonal degree 2 order 0
            cent_pot(2)=-cent_pot(1)*sqrt(5.d0)/5.d0
!             !convert arcsec to dimensionless units and excess LOD (in seconds) to dimensionless units
!             mpolar(1)=mpolar(1)*pi/(OHM_E*3600*180)
!             mpolar(2)=mpolar(2)*pi/(OHM_E*3600*180)
!             mpolar(3)=OHM_E*mpolar(3)/(2*pi)
!             mpolar(3)=-mpolar(3)/(1+mpolar(3))

            !now the Earths rotation vector is represented ( in radians per second)
!                   /  mpolar(1)  \ 
!            omega=|   mpolar(2)   | * OHM_E
!                   \ 1+mpolar(3) /
         else if(dum(3:3) .ne. '')then
            write(stderr,*)"ERROR, -Z option not understood: ",trim(dum)
            stop
         end if
      case('T')!specify times explicitly
         forcetime=.true.
         ind=index(dum,',')
         ind2=index(dum(ind+1:),',')+ind
         read(dum(3:ind-1),*)itstart
         read(dum(ind+1:ind2-1),*)itcent
         read(dum(ind2+1:),*)itend

      case default
         write(stderr,*)'ERROR SH_dump:unknown option: ',dum(2:2)
         call help()
      end select
   end if
end do



!input checks
if(ms >1)then !ambiguous options were requested
   write(stderr,*)"ERROR SH_dump: Ambiguous -m options applied"
   stop
end if

if(gtyp .ne. -9999 .and. .not. limdeg)then
   write(stderr,*)'ERROR SH_dump: Using -ti requires the use of -lLMAX'
   call help()
end if

if( unitin .and. .not. limdeg)then
   write(stderr,*)'ERROR SH_dump: Having a unit input requires the use of -lLMAX'
   stop
end if

if(unitin .and. rotpot)then
   write(stderr,*)'ERROR SH_dump: Cannot choose -u and -Z at the same time'
   stop
end if

if(rotate .and. error)then
   write(stderr,*)"ERROR: error propagation is not allowed with rotations"
   stop
end if

if(.not. limdeg .and. (rotpot .or. cons_angmom))then !set maximum degree to 2 if not specified explicitly
   limdeg=.true.
   lmax=2
   lmin=0
end if

if(cons_angmom .and. inv)then
   write(stderr,*)"ERROR: -m7 does not support inverse (yet)"
   stop
end if


if(c20type >0 .and. inv)then
   write(stderr,*)"ERROR: Option -j does not support inverse (-i)"
   stop
end if


! write(*,*)'subtfiles ',trim(subtfile(1)),' ',trim(subtfile(2))
! write(*,*)'addfiles ',trim(addfile(1)),' ',trim(addfile(2))
! write(*,*)lmax,lmin
! write(*,*)'stdin',stdin,' ',trim(filen)
! stop
!Write(*,*)'options',geoid,water,geowater,smooth,geoc,filter,inv,limdeg,uplift,subt,add
!WRITE(*,*)'TEST',trim(filen)
!now check some metadata of the main file
if(gtyp .eq. -9999 .and. stdin .and. .not. (unitin .or. rotpot))then
   call SH_readmeta(type=gtyp,lmax=ind,tcent=tcent,tstart=tstart,tend=tend)
else if(.not. stdin .and. gtyp .eq. -9999 .and. .not. (unitin .or. rotpot))then
   call SH_readmeta(filen=trim(filen),type=gtyp,lmax=ind,tcent=tcent,tstart=tstart,tend=tend)
end if

!explicitly use different times
if(forcetime)then
   tcent=itcent
   tstart=itstart
   tend=itend
end if


!adapt lmax if not specifically limited

!give an error when lmax and lmin are not given but the type is of 0
if(gtyp .eq. 0 .and. .not. limdeg)then
   write(stderr,*)'ERROR SH_dump:For SH files which are of the generic type the -l option must be used'
   call help()
end if
if(.not.limdeg)lmax=ind

!now allocate required vectors
pos=SH_pos(lmax,lmax)
allocate(clm(pos),slm(pos))
if(unitin)then
   clm=1.d0
   slm=1.d0
   do,l=lmin,lmax
      slm(SH_pos(l,0))=0.d0
   end do
else if(rotpot)then ! rotational potential
   slm=0.d0
   clm=0.d0
   !degree and order 0
   clm(1)=(OHM_R)**2*cent_pot(1)
   !zonal degree 2 order 0
   clm(SH_pos(2,0))=(OHM_R)**2*cent_pot(2)
else !initialize to zero
   clm=0.d0
   slm=0.d0
end if



!allocate additional multiplication vector when required
if(water .or. geowater .or. normalize .or. cons_angmom .or. radialacc .or. gravanom)allocate(vlm(pos))

!allocate load love number vectors if required
if(water .or. geowater .or. uplift .or. sealev .or. loadNdirect .or. horiz)allocate(hnm(pos),knm(pos),lnm(pos))

!allocate addition files if required
if(add > 0)allocate(clm_add(pos),slm_add(pos))

!allocate subtraction files if required
if(subt > 0)allocate(clm_subt(pos),slm_subt(pos))

!allocate error vector if required
if(error)then
allocate(clm_sig(pos),slm_sig(pos))
clm_sig=0.d0
slm_sig=0.d0
end if

if(dot .or. harm)then
   posdot=SH_pos(lmaxdot,lmaxdot)
   allocate(clm_tv(posdot,5),slm_tv(posdot,5))
   clm_tv=0.d0
   slm_tv=0.d0
end if

!now read in main file
if(stdin .and. .not. (unitin .or. rotpot))then
   if(error .and. ( dot .or. harm) )then
      call SH_readgrav(clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=gtyp&
           ,clm_tv=clm_tv,slm_tv=slm_tv)!use formal errors and dot terms contained in file
   else if(dot .or. harm)then
      call SH_readgrav(clm=clm,slm=slm,type=gtyp,clm_tv=clm_tv,slm_tv=slm_tv)!use dot terms contained in file
   else if(error)then
      call SH_readgrav(clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=gtyp)!use formal errors contained in file
   else
      call SH_readgrav(clm=clm,slm=slm,type=gtyp)
   end if
else if( .not. (unitin .or. rotpot))then
   if(error .and. ( dot .or. harm))then
      call SH_readgrav(filen=trim(filen),clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,&
           type=gtyp,clm_tv=clm_tv,slm_tv=slm_tv)
   else if(dot .or. harm)then
      call SH_readgrav(filen=trim(filen),clm=clm,slm=slm,type=gtyp,clm_tv=clm_tv,slm_tv=slm_tv)
   else if(error)then
      call SH_readgrav(filen=trim(filen),clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=gtyp)
   else
      call SH_readgrav(filen=trim(filen),clm=clm,slm=slm,type=gtyp)
   end if
end if




!!ADD rates back ( at center time)
if(dot .or. harm)then
   tmat=0.d0
   if(dot)tmat(1)=(tcent-dotref)
   if(harm)then
      if(tend>tstart)then !integrate over time, when a valid time period exists
         !annual
         tmat(2)=-(cos(2*pi*(tend-dotref))-cos(2*pi*(tstart-dotref)))/(2*pi*(tend-tstart))
         tmat(3)=(sin(2*pi*(tend-dotref))-sin(2*pi*(tstart-dotref)))/(2*pi*(tend-tstart)) 

         !semi-annual
         tmat(4)=-(cos(4*pi*(tend-dotref))-cos(4*pi*(tstart-dotref)))/(4*pi*(tend-tstart))
         tmat(5)=(sin(4*pi*(tend-dotref))-sin(4*pi*(tstart-dotref)))/(4*pi*(tend-tstart)) 

      else ! just take the center time point
         !annual
         tmat(2)=sin(2*pi*(tcent-dotref))
         tmat(3)=cos(2*pi*(tcent-dotref))
         !semi-annual
         tmat(4)=sin(4*pi*(tcent-dotref))
         tmat(5)=cos(4*pi*(tcent-dotref))
      end if
   end if
!   write(0,*)tend,tstart,dotref,tmat
!   write(0,*)clm_tv(SH_pos(4,4),:),clm(SH_pos(4,4))
   !update clm vector
   call dgemm('N','N',posdot,1,5,1.d0,clm_tv(1,1),posdot,tmat,5,1.d0,clm(1),pos)
   !update slm vector
   call dgemm('N','N',posdot,1,5,1.d0,slm_tv(1,1),posdot,tmat,5,1.d0,slm(1),pos)
end if


!apply normalization factor

if(normalize)then
      do,l=0,lmax
         !zero order term
         vlm(SH_pos(l,0))=sqrt(2.d0*dble(l)+1.d0)
         do,m=1,l
            vlm(SH_pos(l,m))=sqrt(2*(2.d0*dble(l)+1.d0))*qfact(l-m,l+m)
         end do
      end do

      if(inv)then
         clm=clm*vlm
         slm=slm*vlm
      else
         clm=clm/vlm
         slm=slm/vlm
      end if
      
      if(error)then
         if(inv)then
            clm_sig=clm_sig/vlm
            slm_sig=slm_sig/vlm
         else
            clm_sig=clm_sig*vlm
            slm_sig=slm_sig*vlm
         end if
      end if
      vlm=0.d0
end if


!When a separate error file was provided replace the vules in c(s)lm_sig with those of error file

if(errfile .ne. '')then !only occurs with -e option
   clm_sig=0.d0
   slm_sig=0.d0
   call SH_readmeta(filen=trim(errfile),type=etyp)
   call SH_readgrav(filen=trim(errfile),clm=clm_sig,slm=slm_sig,type=etyp)
end if


!replace c20 values
if(c20type>0 .and. lmin <=2 .and. lmax >=2) then
   clm(SH_pos(2,0))=c20_replace(time=tstart,time1=tend,type=c20type)
end if


!write(*,*)'read', filen
!now read in subtract files and subtract from main file
do,i=1,subt!
!write(*,*)'subtracting ',trim(subtfile(i))
call SH_readgrav(filen=trim(subtfile(i)),clm=clm_subt,slm=slm_subt)

!check whether the mainfile has degree 1 and/or zero coefficients
!and if not don't subtract those (set them to zero in slm_subt)


if(abs(clm(1)) < 1.d-100)clm_subt(1)=0.d0
if(abs(clm(2)) < 1.d-100)clm_subt(2)=0.d0
if(abs(clm(3)) < 1.d-100)clm_subt(3)=0.d0
if(abs(slm(3)) < 1.d-100)slm_subt(3)=0.d0
!subtract
clm=clm-clm_subt
slm=slm-slm_subt
end do


!apply filter if requested (after subtraction of a possible reference field)
! if(filter) then
!    !call SH_filter(clm,slm,ftyp)!not yet supported Thu May 31 16:40:25 2007
! end if


!same for add files
do,i=1,add
!write(*,*)'adding ',trim(addfile(i))
call SH_readgrav(filen=trim(addfile(i)),clm=clm_add,slm=slm_add)
!add
clm=clm+clm_add
slm=slm+slm_add
end do

!add or subtract permanent tide if requested
if(ptide) then
!permanent tide defined as
clm(4)=clm(4)+pfac*(-0.31460d0)*0.30190d0/(RE*sqrt(4*pi))
end if

!apply geocenter motion if requested

if(geoc) then
   if(geocJ2)then
      call SH_geocmod(c10=c10,c11=c11,s11=s11,c30=c30,c31=c31,s31=s31,time=tstart,time1=tend,mod=cmod,geocfile=geocfile)
      ! clm(SH_pos(3,0))=c30
      ! clm(SH_pos(3,1))=c31
      ! slm(SH_pos(3,1))=s31

       clm(SH_pos(3,0))=clm(SH_pos(3,0))+c30
       clm(SH_pos(3,1))=clm(SH_pos(3,1))+c31
       slm(SH_pos(3,1))=slm(SH_pos(3,1))+s31
   else
      call SH_geocmod(c10=c10,c11=c11,s11=s11,time=tstart,time1=tend,mod=cmod,geocfile=geocfile)
      
   end if

   !call SH_geocmod(c10=c10,c11=c11,s11=s11,time=tcent,mod=cmod,geocfile)
      clm(2)=c10
      clm(3)=c11
      slm(3)=s11


end if



!now apply smoothing if requested
if(smoothrecur) then
   call SH_gaus_recur(clm,radius*1.d3)
   call SH_gaus_recur(slm,radius*1.d3)

   if(error)then
      call SH_gaus_recur(clm_sig,radius*1.d3)
      call SH_gaus_recur(slm_sig,radius*1.d3)
   end if
else if(smooth) then
   call SH_gaus(clm,radius*1.d3)
   call SH_gaus(slm,radius*1.d3)

   if(error)then
      call SH_gaus(clm_sig,radius*1.d3)
      call SH_gaus(slm_sig,radius*1.d3)
   end if
end if


if(compl)then !take the complement of the vector

   !construct complementary vector of the basin function
!construct normal vector (normalization vector)
   clm(1)=1-clm(1)
   clm(2:)=-clm(2:)
   slm(2:)=-slm(2:)
end if


!now apply type specific multiplication
if(water .or. geowater .or. uplift .or. sealev .or. loadNdirect  .or. horiz)then
!get the load love numbers
   if(rigid)then
      lnm=0.d0
      hnm=0.d0
      knm=1.d0
   else
      call SH_loadlove(hnm=hnm,knm=knm,lnm=lnm,typ=ltyp,frame=FR,fileop=trim(lovefile))
      knm=knm+1.d0
   end if

end if



if(water .or. geowater .or. cons_angmom)then!construct degree dependent vector
 
   do,i=0,lmax
      tmpfac=dble(2*i+1)
      st=SH_pos(i,0)
      nd=SH_pos(i,i)
      vlm(st:nd)=tmpfac
   end do
else if(radialacc)then !multiply with l+1

   do,i=0,lmax
      tmpfac=dble(i+1)
      st=SH_pos(i,0)
      nd=SH_pos(i,i)
      vlm(st:nd)=tmpfac
   end do
else if(gravanom)then
   do,i=0,lmax
      tmpfac=dble(i-1)
      st=SH_pos(i,0)
      nd=SH_pos(i,i)
      vlm(st:nd)=tmpfac
   end do
end if

if(water) then
!apply factors for equivalent water height
!write(0,*)clm(1)
   if(inv)then
     clm=clm*knm/(vlm*fact)
     slm=slm*knm/(vlm*fact)
     if(error)then
        clm_sig=clm_sig*knm/(vlm*fact)
        slm_sig=slm_sig*knm/(vlm*fact)
     end if
   else
      clm=clm*vlm*fact/knm
      slm=slm*vlm*fact/knm
      if(error)then
         clm_sig=clm_sig*vlm*fact/knm
         slm_sig=slm_sig*vlm*fact/knm
      end if
   end if



else if(geowater)then
!apply factors for sea level (wrt to a mean solid Earth)


   if(inv)then
     clm=clm*knm/(vlm*fact-RE*hnm)
     slm=slm*knm/(vlm*fact-RE*hnm)
     if(error)then
         clm_sig=clm_sig*knm/(vlm*fact-RE*hnm)
         slm_sig=slm_sig*knm/(vlm*fact-RE*hnm)
     end if
   else
      clm=clm*vlm*fact/knm-RE*clm*hnm/knm
      slm=slm*vlm*fact/knm-RE*slm*hnm/knm
      if(error)then
          clm_sig=clm_sig*(vlm*fact/knm-RE*hnm/knm)
          slm_sig=slm_sig*(vlm*fact/knm-RE*hnm/knm)
      end if
   end if


else if(uplift)then
!apply factors for sea level (wrt to a solid rigid Earth)


   if(inv)then
     clm=clm/(hnm*fact)
     slm=slm/(hnm*fact)
     if(error)then
     clm_sig=clm_sig/(hnm*fact)
     slm_sig=slm_sig/(hnm*fact)
     end if
   else
      clm=fact*clm*hnm
      slm=fact*slm*hnm
      if(error)then
      clm_sig=fact*clm_sig*hnm
      slm_sig=fact*slm_sig*hnm
      end if
   end if
else if (horiz)then ! horizontal deformation

   if(inv)then
     clm=clm/(lnm*fact)
     slm=slm/(lnm*fact)
     if(error)then
     clm_sig=clm_sig/(lnm*fact)
     slm_sig=slm_sig/(lnm*fact)
     end if
   else
      clm=fact*clm*lnm
      slm=fact*slm*lnm
      if(error)then
      clm_sig=fact*clm_sig*lnm
      slm_sig=fact*slm_sig*lnm
      end if
   end if

else if(geoid)then
   if(inv)then
      clm=clm/fact
      slm=slm/fact
      if(error)then
         clm_sig=clm_sig/fact
         slm_sig=slm_sig/fact
      end if
   else
      clm=fact*clm
      slm=fact*slm
      if(error)then
         clm_sig=fact*clm_sig
         slm_sig=fact*slm_sig
      end if
   end if
else if(sealev)then

   if(inv)then
      !degree zero:Eustatic sea level
       clm(1)=clm(1)/(c00sc*RE*rho_e/(rho_w*3))

      do,i=2,pos
         clm(i)=clm(i)/(fact*(knm(i)-hnm(i)))
         slm(i)=slm(i)/(fact*(knm(i)-hnm(i)))
      end do

      if(error)then
         clm_sig(1)=clm_sig(1)/(c00sc*RE*rho_e/(rho_w*3))
         do,i=2,pos
            clm_sig(i)=clm_sig(i)/(fact*(knm(i)-hnm(i)))
            slm_sig(i)=slm_sig(i)/(fact*(knm(i)-hnm(i)))
         end do
     end if
   else
      !degree zero: Eustatic sea level (spread over the ocean surface)
       clm(1)=clm(1)*c00sc*(RE*rho_e/(rho_w*3))

      do,i=2,pos
         clm(i)=fact*clm(i)*(knm(i)-hnm(i))
         slm(i)=fact*slm(i)*(knm(i)-hnm(i))
      end do
      if(error)then
         clm_sig(1)=clm_sig(1)*c00sc*(RE*rho_e/(rho_w*3))
         do,i=2,pos
            clm_sig(i)=fact*clm_sig(i)*(knm(i)-hnm(i))
            slm_sig(i)=fact*slm_sig(i)*(knm(i)-hnm(i))
       end do
      end if
   end if
else if (loadNdirect)then
   if(inv)then
     clm=clm/(knm*fact)
     slm=slm/(knm*fact)
     if(error)then
     clm_sig=clm_sig/(knm*fact)
     slm_sig=slm_sig/(knm*fact)
     end if
   else
      clm=fact*clm*knm
      slm=fact*slm*knm
      if(error)then
      clm_sig=fact*clm_sig*knm
      slm_sig=fact*slm_sig*knm
      end if
   end if
   

else if(cons_angmom)then
   !convert loading potential to surface loading
   clm=fact*vlm*clm
   slm=fact*vlm*slm
   
! xtract relevant coefficients from input
   d2coef(1)=clm(SH_pos(0,0))
   d2coef(2)=clm(SH_pos(2,0))
   d2coef(3)=clm(SH_pos(2,1))
   d2coef(4)=slm(SH_pos(2,1))

   
   !create transformation matrix
   Rmat=0.d0


   if(cons_anginfo)then
      !create intermediate transformation matrix
      Ji3tmp=Ji3_2_m(ltyp)
      Rmat(1:3,:)=matmul(Ji3tmp,SurfL_2_Ji3())
      !calculate polar motion
      mpolar=matmul(Rmat(1:3,:),d2coef)
!      mpolar(3)=-(OHM_E*sqrt(mpolar(1)**2+mpolar(2)**2+(mpolar(3)+1)**2)-OHM_E)/0.843994804d-9
      mpolar(3)=86400*(1.d0-sqrt(mpolar(1)**2+mpolar(2)**2+(mpolar(3)+1)**2))
      mpolar(1)=mpolar(1)*180*3600/(pi)
      mpolar(2)=-mpolar(2)*180*3600/(pi)
      write(stderr,*)"Derived polar motion, X[''],Y[''],LOD[s]",mpolar
   end if
   
   ! write(0,*)'S2J',SurfL_2_Ji3()
   ! write(0,*)"Ltyp",ltyp
   Ji3tmp=Ji3_2_m(ltyp)
   ! write(0,*)'J2m',Ji3tmp
   ! write(0,*)'m2L',m_2_rotpot()
   ! stop

   
   Rmat=matmul(m_2_rotpot(),matmul(Ji3tmp,SurfL_2_Ji3()))
   d2coef=(RE/GM)*matmul(Rmat,d2coef)
   
   !put coefficients back in the output vector
   clm=0.d0
   slm=0.d0

   clm(SH_pos(0,0))=d2coef(1)
   clm(SH_pos(2,0))=d2coef(2)
   clm(SH_pos(2,1))=d2coef(3)
   slm(SH_pos(2,1))=d2coef(4)
   
else if(radialacc .or. gravanom)then
   if(inv)then
     clm=clm/(vlm*fact)
     slm=slm/(vlm*fact)
     if(error)then
        clm_sig=clm_sig/(vlm*fact)
        slm_sig=slm_sig/(vlm*fact)
     end if
   else
      clm=clm*vlm*fact
      slm=slm*vlm*fact
      if(error)then
         clm_sig=clm_sig*vlm*fact
         slm_sig=slm_sig*vlm*fact
      end if
   end if
else if(custom)then
   clm=fact*clm
   slm=fact*slm
   if(error)then
      clm_sig=fact*clm_sig
      slm_sig=fact*slm_sig
   end if

else if(customdegree)then
   unit=17

!   write(0,*)trim(scalefile)
   open(unit=unit,file=trim(scalefile))
   ierr=0
   do,i=0,lmax
      read(unit=unit,fmt=*,iostat=ierr)l,fact
      if(ierr .ne. 0) exit
      if(l>lmax)exit
      st=SH_pos(l,0)
      nd=SH_pos(l,l)
      
      clm(st:nd)=clm(st:nd)*fact
      slm(st:nd)=slm(st:nd)*fact
      if(error)then
         clm_sig(st:nd)=clm_sig(st:nd)*fact
         slm_sig(st:nd)=slm_sig(st:nd)*fact
      end if
   end do
   close(unit)

end if

if(rotate)then
   allocate(dj(lmax+1,lmax+1,lmax+1))
   !construct part of the wigner D rotation matrix
   call djpi2(dj, lmax)

   nd=SH_pos(lmax,lmax)

   allocate(Csh(2,nd),Csh_rot(2,nd)) ! allocate complex Spehrical harmonic vector

   sq4pi=sqrt(pi)*2.d0
   sq2pi=sqrt(2*pi)
   !copy data in temporary array and convert to complex SH in Varshalovich norm
   do,l=0,lmax
      !order zero part
      i=SH_pos(l,0)
      Csh(1,i)=sq4pi*clm(i)
      Csh(2,i)=0.d0
      do,m=1,l
         i=SH_pos(l,m)
         !real part
         Csh(1,i)=sq2pi*clm(i)*(1.d0-2.0*mod(m,2))
         !imaginary part
         Csh(2,i)=-sq2pi*slm(i)*(1.d0-2.0*mod(m,2))
      end do
   end do

   
   !rotate the SH coefficients
   call SHRotateCoef(angles, Csh, Csh_rot, dj, lmax)
   
   !copy back in original array and convert back to geodesy 4pi normalized 
   do,l=0,lmax
      !order zero part
      i=SH_pos(l,0)
      clm(i)=Csh_rot(1,i)/sq4pi
      do,m=1,l
         i=SH_pos(l,m)
         !real part
         clm(i)=Csh_rot(1,i)*(1.d0-2.0*mod(m,2))/sq2pi
         !imaginary part
         slm(i)=-Csh_rot(2,i)*(1.d0-2.0*mod(m,2))/sq2pi
      end do
   end do
   
   if(polarmotion)then ! calcuate residuals wrt to mean centrifugal bulge
      clm(1)=clm(1)-(OHM_E)**2*cent_pot(1)
      clm(SH_pos(2,0))=clm(SH_pos(2,0))-(OHM_E)**2*cent_pot(2)
   end if
end if

!now write spherical harmonic coefficients to standard output in a specific format
!write(*,*)otyp,tstart,tcent,tend


!set coefficients below lmin to zero
if (limdeg .and. lmin > 0)then
st=1
nd=SH_pos(lmin-1,lmin-1)
clm(st:nd)=0.d0
slm(st:nd)=0.d0
if(error)then
   clm_sig(st:nd)=0.d0
   slm_sig(st:nd)=0.d0
end if
end if


if(error)then
call SH_write(clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,tstart=tstart,tcent=tcent,tend=tend,typ=otyp)
else
call SH_write(clm=clm,slm=slm,tstart=tstart,tcent=tcent,tend=tend,typ=otyp)
end if
end program SH_dump

!!subroutine which displays the help function when wrong syntax is tried
subroutine help()
use shtlbx, only :A_oce
implicit none
character(3)::frmt
integer::stderr
stderr=0
frmt='(A)'
write(stderr,frmt)'usage SH_dump [options] filename'
write(stderr,frmt)'SH_dump reads in various GRACE formats and dumps them to standard output (or a file)'
write(stderr,frmt)'Further processing can be performed with the following options:'
write(stderr,frmt)''
write(stderr,frmt)'-mNUM: output predefined parameters from geopotential files, where NUM is one of the numbers:'
write(stderr,frmt)'   1: geoid [m] (but make sure an ellipsoid is subtracted) (applies Bruns formula Dv/g)'
write(stderr,frmt)'   2: equivalent water height [m] (appropriate static field must be extracted) '
write(stderr,frmt)'      Multiplication with: (radius_earth * density_earth * (2l+1)/(3*rho_water*(1+kn))'
write(stderr,frmt)'   3: geocentric equivalent water height (no ocean loading, to compare with altimetry) [m]'
Write(stderr,frmt)'      Multiplication with:(radius_earth*density_earth*(2l+1)/(3*rho_water*(1+kn))'
Write(stderr,frmt)'      -hn*radius_earth)'
write(stderr,frmt)'   4: Uplift (hn Phi_load/g) INPUT IS LOADING POTENTIAL!'
write(stderr,frmt)'   5[-][o]: Static( no circulation assumed) sea level caused by loading potential: (1+kn-hn)Phi_load/g'
write(stderr,frmt)'      Note: if a - is supplied the degree 0 term denotes the eustatic sea level change'
write(stderr,frmt)'      If the degree 0 term ( in terms of mass) is taken OUT of the ocean. '
write(stderr,frmt)"      The option 'o' instead of '-' implies a scaling (1+k0-h0)=1"
write(stderr,frmt)'      Default degree 0 term is the eustatic sea level when the degree 0 is put IN the ocean'
write(stderr,*)'      Ocean surface to global surface ratio used:',A_oce
write(stderr,frmt)'      INPUT IS LOADING POTENTIAL!!'
write(stderr,frmt)'   6: Convert LOADING (surface or tidal) potential to LOADING+SOLID EARTH potential (1+kn)Phi_load'
write(stderr,frmt)'   7[i]: Convert surface LOADING potential to centrifugal potential pertubations, using conservation'
write(stderr,frmt)'      of angular momentum. This effectively takes only degree 20, 21 and 00 coefficients and applies'
write(stderr,frmt)"      linearized elastic perturbation theory. If 'i' is supplied the derived polar motion"
write(stderr,frmt)"      is printed to standard error"
write(stderr,frmt)"   8: Convert Loading Potential to Horizontal deformation. To obtain North and East components, "
write(stderr,frmt)"      differentiate the output in North and East direction in a postprocessing step"
write(stderr,frmt)"   9: Convert Stokes coefficients to radial acceleration in microgal"
write(stderr,frmt)"   10: Convert Stokes coefficients to Gravity anomalies in microgal"
write(stderr,frmt)'   sSCALE: custom scale factor of size SCALE'
write(stderr,frmt)'   dSCALEFILE: use a degree wise scale from file (format file: 1st column:degree, 2nd column:scale)'
write(stderr,frmt)'-ENUM:  Use the following choices of Earth models for the Love numbers NUM is:'
write(stderr,frmt)'    1: Farrells load love numbers (GUTENBERG-BULLEN Earth model)'
write(stderr,frmt)'    2: Use PREM load love numbers (PREM Earth model) (default)'
write(stderr,frmt)'    3: Use Guttenberg-Bullen BODY love numbers (degree 2-4 only)'
write(stderr,frmt)'    4: Use PREM BODY love numbers (degree 2-4 only)'
write(stderr,frmt)'    5: Use Desai 2002  BODY love numbers (degree 2 only)'
write(stderr,frmt)'    0CUSTOMFILE: use a custom file (same format as the above models'
write(stderr,frmt)'-R: Assume rigid Earth (sets love numbers to zero)( may be used to retrieve loading potential'
write(stderr,frmt)'    from equiv. wat)'
write(stderr,frmt)' '
write(stderr,frmt)'-i: take the inverse of the factors above (handy for converting back to geopotential coefficients)'
write(stderr,frmt)''
write(stderr,frmt)'-w[r]RAD: Apply gaussian smoothing with halfwidth radius RAD [km].'
write(stderr,frmt)'  use -wrRAD to use the recursive method from Wahr 1998, else use the method from Chambers 2004'
write(stderr,frmt)''
! write(stderr,frmt)'-fNUM: Filter coefficients with NUM representing the following techniques:'
! write(stderr,frmt)'       1: Swenson wahr filter'
! write(stderr,frmt)'       2: Chambers filter (similar to swenson wahr but somewhat smoother for oceanographic applications)'
! write(stderr,frmt)'       3: Kusche anisotropic filter' 
write(stderr,frmt)''
write(stderr,frmt)'-aFIELD: Add Spherical harmonic field from file FIELD to the main SH field (up to 50 -a are allowed)'
write(stderr,frmt)'-sFIELD: Subtract Spherical harmonic field from file FIELD to the main SH field (up to 50 -s are allowed)'
write(stderr,frmt)'         The -w option (above) is applied after the subtraction of the files, but before adding of the files'
write(stderr,frmt)''
write(stderr,frmt)'-[gG]NUM:  Apply (semi)annual geocenter model according to the following models NUM:'
write(stderr,frmt)'        1: Cretaux 2002 model(LAGEOS,TOPEX/POS laser and DORIS measurements'
write(stderr,frmt)'        2: Bouille model 2000, LAGEOS, T/P and SPOT DORIS measurements'
write(stderr,frmt)'        3: Chen(1999)  LAGEOS'
write(stderr,frmt)'        4: Mark-Willem Jansen,Kusche,Schrama model, GPS from combined GPS-GRACE inversion (period 2003-2006.5)'
write(stderr,frmt)'        5: Roelof Rietbroek from combined FESOM+GPS+GRACE (annual + semi annual) (2003-2007)'
write(stderr,frmt)'        6: Rietbroek et al 2011 from combined FESOM+GPS+GRACE (annual + semi annual) (2003-2008)'
write(stderr,frmt)'        7 CUSTOMGEOCENTERFILE : Apply a geocenter model from a custom file.'
write(stderr,frmt)'          columns must be in chronological order and consist of time (years), x(m), y(m) z(m).'
write(stderr,frmt)'          Values will be averaged over the period' 
write(stderr,frmt)'        note that the time must be extractable from the input files (GRACE,ICGEM,GINS or default output format)'
write(stderr,frmt)"        Use the large 'G' to also compute the effect of the Earth's flattening in the shifted reference frame "
write(stderr,frmt)"        This effects only the coefficients C30, C31 and S31"
write(stderr,frmt)''
write(stderr,frmt)'-jTYPE: Replace c20 values ( in potential CHANGES)'
write(stderr,frmt)'        TYPE may be:'
write(stderr,frmt)'        1: Ries et al series ( from SLR)'
write(stderr,frmt)'        2: JIGOG Rietbroek et al 2011 series'
write(stderr,frmt)'        Value will be replaced BEFORE -s and -a options'

write(stderr,frmt)''


write(stderr,frmt)'-toNUM: Specify the output type NUM:'
write(stderr,frmt)'       1:GRACE(GFZ,JPL,CSR) type'
write(stderr,frmt)'       2:GINS type (French GRGS solution)'
write(stderr,frmt)'       3:ICGEM type'
write(stderr,frmt)'       4:(default) short type first line contains the degree supported, start, center end end time'
write(stderr,frmt)'          of observation'
write(stderr,frmt)'         the other lines are SH coefficients'
write(stderr,frmt)'       6:Linear format each lines contains G[S/C]N  DEGORD VALUE STD'
write(stderr,frmt)'       0:clean no header only SH coefficients'
write(stderr,frmt)'-tiNUM: Specify the input type NUM (see above):'
write(stderr,frmt)'       This option requires -lLMAX to be specified as it will skip the read of meta data. '
write(stderr,frmt)'       This option can help in reading non-default data from standard input'
write(stderr,frmt)'-lLMAX[,LMIN]: Limit the solution to LMAX (and optionally LMIN): this causes the program to consider'
write(stderr,frmt)'               only coefficients up to degree and order LMAX (and ignore coefficients with degree < LMIN)'
write(stderr,frmt)'               Default uses the maximum degree of the main file'
write(stderr,frmt)''
write(stderr,frmt)'-pN: add (N=a) or subtract (N=s) the permanent (deformation) tide from the c20 coefficient '
write(stderr,frmt)'     (tidal love number is set to 0.3019)'
Write(stderr,frmt)'     this is useful for converting between tide free and zero tide systems'
write(stderr,frmt)''
write(stderr,frmt)'-e[f ERRFILE]: Also calculate simple error propagation (no correlation assumed),'
write(stderr,frmt)'    possible error estimates in the add or subtract files are NOT propagated only that of the main file'
write(stderr,frmt)'    Furthermore errors originating from geocenter models are not incorporated either'
write(stderr,frmt)'    if the option is provided as -ef ERRFILE then error coefficients from the file ERRFILE will be used.'
write(stderr,frmt)''
write(stderr,frmt)'-c:  Take the complementary part of the input data.(only really useful for basin expansions) '
write(stderr,frmt)''
write(stderr,frmt)'-d[a][LMAXDOT,REFYEAR]: Add (conventional) rates or harmonics back to coefficients (only supported formats)'
write(stderr,frmt)'  select a to restore the (semi)annual harmonics (a combination of -d and -da restores both) '
write(stderr,frmt)'                     Defaults: t, LMAXDOT=4, REFYEAR=2000, rates are assumed to per year'

write(stderr,frmt)''
! write(stderr,frmt)'-G: File is generic this will allow a generic file to be read from standard input'
! write(stderr,frmt)'     (checks for filetype and time extraction will be skipped)'
! write(stderr,frmt)''
write(stderr,frmt)'-N: Normalize (4pi normalization) the coefficients. NOTE: only the main shfile will be normalized not'
write(stderr,frmt)'    the subtract/add files '
write(stderr,frmt)'    The -N option with -i applies the inverse normalization'
write(stderr,frmt)'-u: Set the input to 1. This option is useful to find out what scales are applied to the coefficients'
write(stderr,frmt)'-Z[pPX,Y,LOD]: Set the input to the (mean) centrifugal potential'
write(stderr,frmt)"  If the polar motion (X,Y in arcsec and LOD excess length of day in seconds) are supplied the"
write(stderr,frmt)'  PERTURBATIONS due to the centrifugal are calculated by rotating the perturbed centrifugal bulge.'
write(stderr,frmt)'  This automatically implies the option -r-GAMMA, -BETA, -ALPHA with:'
write(stderr,frmt)'  ALPHA=atan2(X,Y)+pi/2, BETA=sqrt(X^2+Y^2), GAMMA=0 (rotation of the bulge around its axis of symmetry'
write(stderr,frmt)'  may be anything.'
write(stderr,frmt)'  Alternatively, one may write -ZpTIME or zPTIME, where TIME is the time in MJD.(EOP values are obtained by'
write(stderr,frmt)"  interpolation of the standard EOP time series. Use the small 'p' to use the polar motion relative to the"
write(stderr,frmt)"   mean iers (2010) pole or large 'P' to keep the contribution of the mean pole."
write(stderr,frmt)'  if -l is not given the lmax is automatically limited to degree 2'
write(stderr,frmt)''
write(stderr,frmt)"-rALPHA,BETA,GAMMA: Rotate spherical harmonics with three Euler Angles, around the Z axis , then the"
write(stderr,frmt)"           new Y axis and then the new Z axis. Angles ALPHA, BETA, GAMMA in degrees. The set of"
write(stderr,frmt)"           coeffients will hold for the new rotated reference system. "
write(stderr,frmt)"           Use -GAMMA, -BETA, -ALPHA to rotate the body itself in the current reference frame."
write(stderr,frmt)"           Equivalently, the rotation can be interpreted as an initial rotation "
write(stderr,frmt)"           about the z axis by GAMMA, then a rotation of BETA about the INITIAL X axis, and then a "
write(stderr,frmt)"           rotation of ALPHA about the INITIAL Z axis."
write(stderr,frmt)'           The rotation is performed after all other options'
write(stderr,frmt)''
write(stderr,frmt)"-OISOFRAME: Set the Origin of reference frame to ISOFRAME"
write(stderr,frmt)"            Applies to -g and -s ( degree 1)"
write(stderr,frmt)"            CM: Center of mass of the Earth system"
write(stderr,frmt)"            CF: Center of surface figure (default)"
write(stderr,frmt)"            CE: Center of Mass of the solid Earth"
write(stderr,frmt)"            CL: Center of surface Lateral figure"
write(stderr,frmt)"            CH: Center of surface Height figure"
write(stderr,frmt)'-TTSTART,TCENTER,TEND: explicitly specify timetags for file output, TSTART,TCENTER,TEND'
write(stderr,frmt)'                       in decimal years'
write(stderr,frmt)'example: SH_dump -m1 -g2 -w400 -l50 -f1 SH_file.coef'
write(stderr,frmt)'         Calculates equivalent water height, applies a geocenter model, uses filtering and smoothing'
write(stderr,frmt)'         And it only considers degrees smaller or equal than 50'
write(stderr,frmt)'Another note: the program SH_dump can read from standard input as well so you can use more SH_dumps'
write(stderr,frmt)'behind each other:'
write(stderr,frmt)'example:  SH_dump options filename | SH_dump otheroptions >output'
stop
end subroutine help

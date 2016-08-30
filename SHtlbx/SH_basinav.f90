!!file to make a basin average from a SH mask file and a input file
!!files can be read from standard input so that they can be preprocessed and piped with SH_dump
!!Coded by Roelof Rietbroek, Tue Jun  5 09:52:58 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue Jul  3 15:31:50 2007
!!incorporated error simple error propagation

!!Updated We 23 April 2008, More input files are now allowed (much quicker since the basin file must only be read in once)
!!Updated by Roelof Rietbroek, Fri Nov  5 16:07:51 2010
!!ffixed small bug in format for an error message

!!Updated by Roelof Rietbroek, Mon Jun 25 10:16:54 2012
!! also allow output in gigaton (providing input is in equivalent water)

!!Updated by Roelof Rietbroek, Tue Mar 19 16:59:37 2013
!!optionally apply a automatic rescaling


program SH_basinav
use SHtlbx
use FORTtlbx
implicit none
integer::narg,i,lmax,lmin,lmaxf,lmaxb,gtyp,ind,stpos,ndpos,pos,itharg,l,m,maxf,stderr,nf
parameter(maxf=500)
character(200)::dum,filen(maxf),basinf
logical::limdeg,stdin,error,smooth, gigaton,rescale
double precision::mean,tcent,mean_sig,tmp,norm,landav,rad,fact,rescalefac
double precision,allocatable, dimension(:)::clm,slm,bclm,bslm,clm_sig,slm_sig,combclm,combslm
double precision,allocatable, dimension(:)::clmo,slmo
double precision,allocatable,dimension(:)::bclm_sq,bslm_sq,bclm_orig,bslm_orig
integer::iargc

!get command line arguments
narg=iargc()
basinf=''
lmin=0
limdeg=.false.
stdin=.true. !the default
error=.false.
stderr=0

smooth=.false.
gigaton=.false.
rescale=.false.
itharg=0
nf=0
lmax=-10!gives an error if not redefined
rescalefac=1.d0

if(narg < 1)call help()!call help straight away


do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .ne. '-')then
      nf=nf+1
      filen(nf)=dum
      stdin=.false.
    
   else !if argument is an option
      select case(dum(2:2))
      case('l')!limit maximum (and possibly minum degree)
         limdeg=.true.
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
         read(dum(3:ind-1),*)lmax
         read(dum(ind+1:),*)lmin
         else
         read(dum(3:),*)lmax
         end if
      case('b')!basinfile
         if(dum(3:3) .eq. ' ')then !also accept a space after the -b option
            itharg=itharg+1
            call getarg(itharg,dum)
            basinf=dum
         else
            basinf=dum(3:)
         end if
      case('e')!error propagation
         error=.true.
      case('w')!apply gaussian smoothing
         read(dum(3:),*)rad
         rad=rad*1.d3
         smooth=.true.
      case('G')
         gigaton=.true.
      case('R')
         rescale=.true.
      case default
         write(stderr,*)'unknown option selected, quitting'
         call help()
      end select
   end if
end do


if(stdin)nf=1 !reset amount of files to 1 if reading from standard input



!Check if basinfile has been set
if(basinf .eq. '')then
write(stderr,*)'basin file must be specified, quitting'
call help()
end if 




!write(*,*)lmax,lmin


!!load basin file

call SH_readmeta(filen=trim(basinf),lmax=lmaxb,type=gtyp)


!check if lmax and lmin are specified for generic SH files
if(gtyp .eq. 0 .and. .not. limdeg)then
   write(stderr,*)'ERROR:For SH files which are of the generic type the -l option must be used'
   stop
end if

!limit lmax if requested

if(limdeg .and. lmaxb < lmax)then
   write(stderr,*)"ERROR: Basin coefficient file",trim(basinf),"supports only up to degree",lmaxb 
   stop
else if( limdeg)then
   lmaxb=lmax
end if



pos=SH_pos(lmaxb,lmaxb)
!allocate basin vectors
allocate(bclm(pos),bslm(pos))

call SH_readgrav(filen=trim(basinf),clm=bclm,slm=bslm)


!check whether bc00 is present
if(abs(bclm(1)) < 1.d-100)then
   write(stderr,*)'basin Spherical harmonic set must have a 00 entry, quitting'
   call help()
end if

if(rescale)then !copy the (unfiltered) basin coefficients
   allocate(bclm_orig(pos),bslm_orig(pos))
   bclm_orig=bclm
   bslm_orig=bslm
end if



!smooth basin function if necessary
if(smooth)then
   call SH_gaus(bclm,rad)
   call SH_gaus(bslm,rad)
end if

if(error) then !take also the squared of the basin vector ( needed for error propagation)
   allocate(bclm_sq(pos),bslm_sq(pos))
   bclm_sq=bclm**2
   bslm_sq=bslm**2
end if

if(gigaton)then
   fact=rho_w*1d-12*bclm(1)*(4*pi*RE**2)
else
   fact=1.d0
end if

!loop over input files

do,i=1,nf



!now get time tags and maximum degree supported from main file
if(stdin)then

   call SH_readmeta(lmax=lmaxf,tcent=tcent,type=gtyp)
   
else
   call SH_readmeta(filen=trim(filen(i)),lmax=lmaxf,tcent=tcent,type=gtyp)
   
end if



!check if lmax and lmin are specified for generic SH files
if(gtyp .eq. 0 .and. .not. limdeg)then
   write(stderr,'(A)')'ERROR:For SH files which are of the generic type the -l option must be used'
   stop
end if

!limit lmax if requested

if(limdeg .and. lmaxb < lmax)then
   write(stderr,*)"ERROR: coefficient file",trim(filen(i)),"supports only up to degree",lmaxf 
   stop
else if( limdeg)then
   lmaxf=lmax
end if

if(lmaxb < lmaxf)then
   write(stderr,*)"ERROR: basin coefficients file",trim(basinf),"supports only up to degree",lmaxb, "requested",lmaxf
   stop
end if



!allocate coefficient vectors (only when necessary)
pos=SH_pos(lmaxf,lmaxf)

if( .not. allocated(clm))then
   allocate(clm(pos),slm(pos))
   if(error)allocate(clm_sig(pos),slm_sig(pos))
else
   if(size(clm,1) < pos)then
      deallocate(clm,slm)
      allocate(clm(pos),slm(pos))
      if(error)then
      deallocate(clm_sig,slm_sig)
      allocate(clm_sig(pos),slm_sig(pos))
      end if
   end if
end if



clm=0.d0
slm=0.d0

if(error)then
   clm_sig=0.d0
   slm_sig=0.d0
end if



!load main file coefficients into files
if(stdin) then
   if(error)then
      call SH_readgrav(clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=gtyp)
   else
      call SH_readgrav(clm=clm,slm=slm,type=gtyp)
   end if
else
   if(error)then
      call SH_readgrav(filen=trim(filen(i)),clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=gtyp)
   else
      call SH_readgrav(filen=trim(filen(i)),clm=clm,slm=slm,type=gtyp)
   end if
end if






!now convolve the two sets which each other (and limit the minimum degree)

stpos=SH_pos(lmin,0)
ndpos=SH_pos(lmaxf,lmaxf)

! landav=dot_product(clm(stpos:ndpos),combclm(stpos:ndpos))+dot_product(slm(stpos:ndpos),combslm(stpos:ndpos))
   
! write(*,*)'landav',landav/combclm(1)

if(i==1 .and. rescale)then !compute a rescaling factor
   rescalefac=dot_productblas(bclm_orig(stpos:ndpos),bclm(stpos:ndpos))+dot_productblas(bslm_orig(stpos:ndpos),bslm(stpos:ndpos))
   rescalefac=bclm(1)/(rescalefac)
   write(stderr,*)"APPLYING RESCALE FACTOR:",rescalefac
end if





mean=dot_productblas(clm(stpos:ndpos),bclm(stpos:ndpos))+dot_productblas(slm(stpos:ndpos),bslm(stpos:ndpos))



if(error)then
!take squares
clm_sig(stpos:ndpos)=clm_sig(stpos:ndpos)**2
slm_sig(stpos:ndpos)=slm_sig(stpos:ndpos)**2

mean_sig=dot_productblas(clm_sig(stpos:ndpos),bclm_sq(stpos:ndpos))&
     +dot_productblas(slm_sig(stpos:ndpos),bslm_sq(stpos:ndpos))
mean_sig=sqrt(mean_sig)*fact*rescalefac/(bclm(1))
end if

!divide by bc00 term ( and possibly convert to gigaton)
mean=mean*fact*rescalefac/(bclm(1))

!write results to standard output
if(error)then
   write(*,*)tcent,mean,mean_sig
else
   write(*,*)tcent,mean
end if

end do ! end loop over files


end program SH_basinav


!!help function
subroutine help()
integer::stderr
character(6)::frmt
stderr=0
frmt='(A)'
write(stderr,frmt)'Calculate the basinaverage of multiple SH files'
write(stderr,frmt)'The averages are printed to standard output'
write(stderr,frmt)'usage SH_basinav [options] -bBASINFILE MAINFILES'
write(stderr,frmt)'Where options are:'
write(stderr,frmt)' -lLMAX,LMIN: limit solution to coefficients between degree LMIN and LMAX'
write(stderr,frmt)' -e: also calculate the error propagation from coefficient errors in the SH files'
write(stderr,frmt)' -wRADIUS, apply gaussian smoothing with halfwidth radius RADIUS [km]'
write(stderr,frmt)' -G: Output equivalent water height in Gigaton (only useful when input is in equivalent water height)'
write(stderr,frmt)' -R: apply a simple rescaling to the output to account for the Gaussian filter and truncation'
write(stderr,frmt)'    NOTE: this assumes that the BASINFILE is UNFILTERED'
stop
end subroutine help

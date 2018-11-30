!!program to create mass loading normal systems from OBP data given in netcdf grids
!data in netdcf grids must be provided as triplets in units of Db

!!Coded by Roelof Rietbroek, Fri Feb 15 14:18:43 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Tue May 26 09:41:11 2009
!!Added cosine latitude scaling for errrors
!!Updated by Roelof Rietbroek, Wed Jun 10 15:42:50 2009
!! writes binary data using the newer routine (platform indenpendent)

!!Updated by Roelof Rietbroek, Mon Aug 10 16:48:05 2009
!!big rewrite:
!! one file per run ( normal matrix may be 'stolen' from other file)
!! use of diagonal error covariance


program OBP_normeq
use forttlbx,only : dot_productblas,freeunit
use shtlbx
use netcdf
use gpstlbx
use grdtlbx
use binfiletools
implicit none
integer::stderr,i,narg,itharg,lmax,lmin,j,ind,ncid,nciderr,ncerr,chunkstart
integer::nsh,nunk,nchunks,ndat,datsz,datst,datnd,shift,strlen,unit,eof
parameter(stderr=0)
logical::geocent,globalmean,sh,reject,equiref
character(20)::chargpswk
character(24),pointer::Unlist(:)=>null()
character(200)::dum,outfile,resfile,rejectfile,errfile
integer::count
double precision::db2eqw,sigma0,alpha,ltpl,latreject,lonreject
double precision,allocatable::A(:,:)
double precision,pointer::AtA(:,:)=>null()
double precision,pointer::ab(:)=>null()
double precision,pointer::apriori(:)=>null()
double precision,allocatable::lat(:),lon(:),obp(:)
integer,allocatable::selector(:)
logical::cosscale,diagcov,usenorm
double precision,allocatable::diagerr(:) ! vector with precomputed cos(lat)
integer::k,l
type(BINdat)::out,useN ! output derived type
integer::iargc,gpswk

!defaults/initializations
sigma0=0.01 !assume standard deviation of 1cm equivalent water height for the pressure data
alpha=1.d0/(sigma0**2) ! only when unit or cos(lat scale) error covraiance is used
cosscale=.false.
diagcov=.false.
!default writes to standard output
out%file='stdout'

!read in upper matrix from use system
useN%mtyp='U'

!conversion factor to convert from Db to equivalent water height in m (~1)
db2eqw=1.d4/(g*rho_w)


lmax=60 
lmin=0

geocent=.false.
globalmean=.false.
sh=.false.
reject=.false.
equiref=.false.
usenorm=.false.

!get command line arguments
itharg=0
narg=iargc()
gpswk=0

!process command line arguments
if(narg < 1) call help()

do,i=1,narg
   itharg=itharg+1
   if(itharg> narg)exit
   call getarg(itharg,dum)
!  write(*,*)itharg,dum
   if(dum(1:1) .eq. '-')then !option
      select case( dum(2:2))
      case('l')!limit maximum and minimum degree
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
      case('R') !reject pixels
         reject=.true.
         if(dum(4:4) .eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,rejectfile)
         else
            rejectfile=dum(4:)
         end if
      case('g')!describe degree 1 variation as geocenter motion
         geocent=.true.
         
      case('s')!create system with Spherical harmonic elastic loading
         if(dum(3:3) .eq. 'e')equiref=.true.
         sh=.true.

      case('b')!estimate global mean mass bias
         globalmean=.true.

      case('e')! use covariance from file
         select case(dum(3:3))
         case('d')!diagonal covariance
            diagcov=.true.
            alpha=1.d0 ! reset alpha scale to 1
            if(dum(4:5) .eq. '= ')then
               itharg=itharg+1
               call getarg(itharg,errfile)
            else
               errfile=dum(5:)
            end if
         case('c')!cosine latitude scaled
            cosscale=.true.   
         case default
            write(0,*)"ERROR:not supported yet:",dum(1:3)
         end select
         
      case('u')!use NORMAL matrix from file
         usenorm=.true.
         useN%file=trim(dum(4:))
      case('t')!explicitly set time
         read(dum(4:),*)gpswk
      case default
         write(stderr,*)'Unknown option selected, quitting'
         call help()
      end select


   else
      resfile=trim(dum)
   end if
end do


if(geocent .and. lmin <=1)then
   lmin=2 !reset minimum degree
end if


if(.not. (globalmean .or. geocent .or. sh))then
   write(stderr,*)'At least one of the options g/s/b must be specified to create a normal system'
   stop
end if

!make meta data headers for files
out%nint=7
allocate(out%ints_d(out%nint),out%ints(out%nint))

out%ints=0
out%ints(6)=1
out%ints(3)=0
if(sh)then
   out%ints(4)=lmax
   out%ints(5)=lmin
end if

out%ints_d(1)='Nobs'
out%ints_d(2)='Nunknows'
out%ints_d(3)='Nreduced'
out%ints_d(4)='Lmax'
out%ints_d(5)='Lmin'
out%ints_d(6)='Modnr'
out%ints_d(7)='GPSweek'

out%ndbls=6
allocate(out%dbls_d(out%ndbls),out%dbls(out%ndbls))

out%dbls=0.d0
out%dbls(6)=1.d0
out%dbls(4)=1.d0
out%dbls(5)=0.d0 !no apriori estimate yet


out%dbls_d(1)='CTime'
out%dbls_d(2)='STime'
out%dbls_d(3)='ETime'
out%dbls_d(4)='Weight'
out%dbls_d(5)='LtPL'
out%dbls_d(6)='Sigma0'

if( equiref)then
   out%descr='Normal system from OBP data (relative to static equilibrium sea level)'
else
   out%descr='Normal system from OBP data'
end if
out%nvec=2
out%type='SYMVN___'
out%mtyp='U'
out%nread=5
allocate(out%readme(out%nread))
out%readme=''
out%readme(1)='Normal equation derived from FESOM pseudo observations'
out%readme(2)='Used error-covariance:'
if(diagcov)then
   out%readme(3)='diagonal error covariance from file:'
   out%readme(4)=errfile(1:80)
else
   out%readme(3)='unit diagonal'
end if

if(cosscale)then
   out%readme(5)='Diagonal error scaled by cosine latitude'
end if

!!!!!!!!!!!!!!!!!!!!!NETCDF data reading!!!!!!!!!!!
!read in data from netcdf data file
ncerr=nf90_open(trim(resfile), mode=NF90_NOWRITE, ncid=ncid)

if(ncerr .ne. NF90_NOERR )then
   write(stderr,*)"error opening NETCDF file",trim(resfile)
   stop
end if

!get length of triples (use the first dimension id (latitude) to set it)
ncerr=nf90_inquire_dimension(ncid=ncid, dimid=1, len=ndat)

allocate(lat(ndat),lon(ndat),obp(ndat))



if(gpswk .eq. 0)then ! try to get the gps week from the file
   !get GPS week
   ncerr=nf90_inquire_attribute(ncid=ncid, varid=NF90_GLOBAL, name='GPS_week', len=strlen)
   
   ncerr=nf90_get_att(ncid=ncid, varid=NF90_GLOBAL, name='GPS_week', values=chargpswk(1:strlen))
   !read gps week into an integer
   read(chargpswk(1:strlen),*)gpswk
end if

out%ints(7)=gpswk



write(stderr,*)'Reading file:',trim(resfile)



out%dbls(1)=GPS_year(int(out%ints(7)))
out%dbls(2)=GPS_year(int(out%ints(7)),dble(out%ints(7)))
out%dbls(3)=GPS_year(int(out%ints(7)),dble(out%ints(7)+1))
   
!check whether pixels are rejected from REJECTIONFILE

!get latitude (assumed to be the first variable)
      
ncerr=nf90_get_var(ncid=ncid, varid=1, values=lat)
      
      
!get longitude (assumed to be the second variable)
ncerr=nf90_get_var(ncid=ncid, varid=2, values=lon)

allocate(selector(ndat))
forall(j=1:ndat)selector(j)=j !assume first that all pixels are valid (value is there index )
      
if(reject)then
         
   eof=0 !end of file detector
   unit=freeunit()
   open(unit=unit,file=rejectfile,form='formatted')
   
   do while (eof .eq. 0)! loop untill end of file
      read(unit=unit,fmt=*,iostat=eof)lonreject,latreject
      
      !loop over all data points
      do,j=1,ndat
         if((lon(j) .eq. lonreject) .and. (lat(j) .eq. latreject))selector(j)=0 !set selector to zero for rejected pixel
      end do
      
   end do
   close (unit)
   
         
   !now eliminate zeros from the data list (shift all zeros to the end and reset ndat)
   count=0
   do,j=1,ndat
      if(selector(j) .ne. 0)then
         count=count+1
         selector(count)=j
      end if
      
   end do
   write(stderr,*)'Rejection procedure removed ',ndat-count,'points'            
   !reset ndat
   ndat=count
end if

!resort longitude and latitude
do,j=1,ndat
   lat(j)=lat(selector(j))*pi/180.d0
   lon(j)=lon(selector(j))*pi/180.d0
end do
!precompute cosine lat if necessary
if(cosscale .or. diagcov)then
   allocate(diagerr(ndat))
   diagerr=1.d0 !initialize to 1
end if

if (diagcov)then !read external file
   
   ncerr=nf90_open(trim(errfile), mode=NF90_NOWRITE, ncid=nciderr)
      !get err (assumed to be the third variable)
   ncerr=nf90_get_var(ncid=nciderr, varid=3, values=diagerr)
      !close the data set
   ncerr=nf90_close(nciderr)
   
   !scale to equivalent wwater height
   
   diagerr=db2eqw*diagerr
   
end if

if(cosscale)then !additionally scale the error by cosine(latitude) and invert
   do,j=1,ndat
      diagerr(j)=cos(lat(j))/diagerr(j)
   end do
end if




!get OBP (assumed to be the third variable)
ncerr=nf90_get_var(ncid=ncid, varid=3, values=obp)

!close the data set
ncerr=nf90_close(ncid)

!convert longitude and latitude to radiansand obp to equiv water height and apply selection
!the loop is allowed since selector(j) is always larger or equal than j
if(cosscale .or. diagcov)then
   do,j=1,ndat
      obp(j)=obp(selector(j))*db2eqw*diagerr(j)
   end do
else
   do,j=1,ndat
      obp(j)=obp(selector(j))*db2eqw
   end do
end if




!!!!!!!!!!END NETCDF data READING!!!!!!!!!!!


!calculate the amount of unknowns
nunk=0
nsh=SH_tpos(lmax,lmax,1,lmax,lmin)
if(sh)nunk=nunk+nsh
if(geocent)nunk=nunk+3
if(globalmean)nunk=nunk+1
out%nval1=nunk
out%nval2=nunk
out%pval1=out%nval1*(out%nval1+1)/2
out%pval2=1
  
!calculate a reasable amount of chunks (each observation equation matrix will be ~100Mb)
nchunks=8*nunk*ndat/(1024*1024*100)

if(nchunks .eq. 0)nchunks=1
!allocate space for normal equations 

   




allocate(out%vec(nunk,out%nvec),out%mat1(nunk,nunk))
AtA=>out%mat1
Ab=>out%vec(:,1)
Apriori=>out%vec(:,2)

!make a character array with description of the unknowns
allocate(out%side1_d(nunk))
Unlist=>out%side1_d
Unlist=''
         
apriori=0.d0 !currently there is no apriori info

   
!get observation equation matrix:

datsz=ndat/nchunks
chunkstart=1 !start index for the remaining do loop !defaults to 1
allocate(A(datsz,nunk)) !(maximum size needed)


!first part calculate (possibly) inconsistent part (data size not equal to ndat/nchunks)
datsz=mod(ndat,nchunks)

if(datsz .eq. 0)then 
   datsz=ndat/nchunks
   chunkstart=2 !set the start index of the remaining do loop to 2
end if

   
shift=0
      
if(geocent)then
   call OBP_obseq_geocent(A=A(1:datsz,shift+1:shift+3),tag=Unlist(shift+1:shift+3)&
        ,lon=lon(1:datsz),lat=lat(1:datsz))
   shift=shift+3
end if

   
if(globalmean)then
   call OBP_obseq_globalmean(A=A(1:datsz,shift+1:shift+1),tag=Unlist(shift+1:shift+1))
   shift=shift+1
end if

if(sh)then
   if(equiref)then
      call OBP_obseq_loadsh_equiref(A=A(1:datsz,shift+1:shift+nsh),tag=Unlist(shift+1:shift+nsh),&
           lon=lon(1:datsz),lat=lat(1:datsz),lmax=lmax,lmin=lmin)     
   else
      call OBP_obseq_loadsh(A=A(1:datsz,shift+1:shift+nsh),tag=Unlist(shift+1:shift+nsh),&
           lon=lon(1:datsz),lat=lat(1:datsz),lmax=lmax,lmin=lmin)     
   end if
   shift=shift+nsh
end if

!optionally rescale A with diagerr
if(cosscale .or. diagcov)then
   do,k=1,nunk
      do,l=1,datsz
         A(l,k)=A(l,k)*diagerr(l)
      end do
   end do
end if
   
!create normal matrix
if(usenorm)then !steal the matrix from an external file
   call read_BINtype(useN)
   out%mat1=>useN%mat1 !point the output matrix to the matrix just read
else ! calculate matrix
   call dsyrk('U','T',nunk,datsz,alpha,A,size(A,1),0.d0,AtA,nunk)
end if

!right hand side
call dgemv('T',datsz,nunk,alpha,A,size(A,1),obp,1,0.d0,Ab,1)


!initial ltpl
ltpl=alpha*dot_productblas(obp(1:datsz),obp(1:datsz))

!make a loop over remaining data chunks
!update start index
datst=datsz

!new datasize (chunksize)
datsz=ndat/nchunks


!update new end index
datnd=datst+datsz


do,j=chunkstart,nchunks   !loop over remaining data chunks
   write(stderr,*)'Doing chunk',j,'of',nchunks
   write(stderr,*)'Used datapoints',datnd,'of',ndat
   shift=0
   if(geocent)then
      call OBP_obseq_geocent(A=A(1:datsz,shift+1:shift+3),lon=lon(datst+1:datnd),lat=lat(datst+1:datnd))
      shift=shift+3
   end if
   
   if(globalmean)then
      call OBP_obseq_globalmean(A=A(1:datsz,shift+1:shift+1))
      ! write(*,*)'A bias',A(1,shift+1:shift+1),A(datsz,shift+1:shift+1)
      shift=shift+1
      
   end if
   
   if(sh)then
      if(equiref)then
         call OBP_obseq_loadsh_equiref(A=A(1:datsz,shift+1:shift+nsh),lon=lon(datst+1:datnd),&
              lat=lat(datst+1:datnd),lmax=lmax,lmin=lmin)     
      else
         call OBP_obseq_loadsh(A=A(1:datsz,shift+1:shift+nsh),lon=lon(datst+1:datnd),&
              lat=lat(datst+1:datnd),lmax=lmax,lmin=lmin)     
      end if
      shift=shift+nsh
   end if
   !optionally rescale A with diagerr
   if(cosscale .or. diagcov)then
      do,k=1,nunk
         do,l=1,datsz
            A(l,k)=A(l,k)*diagerr(l)
         end do
      end do
   end if

   !UPDATE normal matrix
   if( .not. usenorm)then ! only do this when it actually needs to be calculated
      call dsyrk('U','T',nunk,datsz,alpha,A,size(A,1),1.d0,AtA,nunk)
   end if
   !AND right hand side
   call dgemv('T',datsz,nunk,alpha,A,size(A,1),obp(datst+1),1,1.d0,Ab,1)

   !update ltpl

   ltpl=ltpl+alpha*dot_productblas(obp(datst+1:datnd),obp(datst+1:datnd))
      
   !update index parameters
   datst=datst+datsz
   datnd=datnd+datsz
end do
   
   
out%dbls(5)=ltpl
out%ints(1)=ndat
out%ints(2)=nunk


!write data to standard output
call write_BINtype(out) 

write(stderr,*)'Finished!'









end program OBP_normeq


subroutine help()
implicit none
character(4)::frmt
integer::stderr
stderr=0
frmt='(A)'

write(stderr,frmt)' Program OBP_normeq creates Ocean bottom pressure normal equation for Spherical harmonic'
write(stderr,frmt)' mass loading coefficients'
write(stderr,frmt)' File is outputted to standard output as unformatted binary file readable by'
write(stderr,frmt)' the FORTRAN read_BINtype routine'
write(stderr,frmt)' Observation equation is based on the surface mass loading of an elastic Earth'
write(stderr,frmt)' Usage: OBP_normeq [OPTIONS] RESFILE'
write(stderr,frmt)' Where RESFILE is the residual file in netcdf format containing triples lon,lat (in degrees),obp (in db)'
write(stderr,frmt)''
write(stderr,frmt)' OPTIONS can be the following:'
write(stderr,frmt)'  -l=LMAX[,LMIN]: limit maximum (and possibly minimum) degree of the output to LMAX and LMIN'
write(stderr,frmt)'                defaults assumes lmax =60 and lmin=1 (or lmin=2 when -g is specified)'
write(stderr,frmt)''
write(stderr,frmt)'  -R=REJECTIONFILE: Reject pixels set in ascii file REJECTIONFILE (two columns longitude,latitude in degrees) '
write(stderr,frmt)'  -g: Express degree 1 loading  variations as geocenter motion'
write(stderr,frmt)'  -s[e]: Construct system with spherical harmonic loading coefficients.'
write(stderr,frmt)'         The option [e] causes the bottom pressure to be reffered to the dynamic part of the sea level only.'
write(stderr,frmt)'  -b: Estimate an unknown mass bias parameter'
write(stderr,frmt)'  -ec: add optional cos(latitude) scaling to the diagonal error covariance.'
write(stderr,frmt)'       We will have an error of errcm/cos(lat).'
write(stderr,frmt)'  -ed=ERRFILE: use a diagonal error covariance from netcdf file ERRFILE.'
write(stderr,frmt)'     the error is assumed to be in the 3rd netcdf variable having the same length as the mainfile variable'
write(stderr,frmt)'  -u=USENORMFILE: avoid computation of the normal matrix by using one from the file USENORMFILE'
write(stderr,frmt)'  -t=GPSWEEK: explicitly set the GPS week of the normal equation'
write(stderr,frmt)'  The routines currently assume a standard deviation of 1 cm'
write(stderr,frmt)'  NOTE:The latitude is assumed to be the first Netcdf variable, the longitude the second and the bottom'
write(stderr,frmt)'  pressure the third '
stop
end subroutine help

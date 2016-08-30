!!program to create stokes coefficients from a netcdf file
!!It uses subroutines from the SHTOOLS library to calculate legendre polynomials
!!reference: swenson wahr Methods for inferring regional surface-mass anomalies from Gravity Recovery and Climate Experiment (GRACE)
!!Measurement of time variable gravity, Journal of geophysical research, Vol. 107, 2002
!!NOTE:The program uses a very simple pointwise integration, which is nice and fast but not always optimal
!!Things to do: apply a 2D FFT and then convert to SH coefficients (source: Global Spherical harmonic computation by two-dimensional fourier methods)
!!Possibly this method can be combined with the non-equidistant FFT method of :
!A. Dutt, V. Rokhlin, Fast fourier transforms for nonequispaced data, SIAM J. Sci Comput.,1993
!this effectively yield resampling on an appropriate grid.

!!adapted from an earlier program by Roelof Rietbroek, 11-5-2007
!! input is a equispaced netcdf file a possible maximum resolution degree and an outputfilename

!!Updated by Roelof Rietbroek, Thu Mar 19 14:13:00 2009
!! added support for a fft approach by the Wieczorek SHTOOLS routine SHExpandDH
!!note only grids with N x 2N points are allowed ( N must be even)

!!Updated by Roelof Rietbroek, Sun Jun 28 22:08:52 2009
!!added greenwich center detection

!!Updated by Roelof Rietbroek, Tue Aug  2 17:16:19 2011
!!added better checks for the netcdf files
!!3D arrays with different names are now also allowed
!!Updated by Roelof Rietbroek, Wed Aug  3 09:46:43 2011
!!added time tag recognition
!!Updated by Roelof Rietbroek, Wed Mar  4 14:49:11 2015
!! also search for 'Latitude' and Longitude' for dimension/variable names

!!! TODO(27 March 2013): use the NC_read subroutine to avoid duplicate code

!!Updated by Roelof Rietbroek, Tue Apr 21 15:24:54 2015
!!fixed bug: fft approach shifted results by half a grid cell!
!!fixed bug: only count values on the pole once, and avoid double counting of eastern boundary values (lon=+180 or lon=360)
!!Updated by Roelof Rietbroek, Wed May 20 11:27:26 2015
!! fized a small bug which prevented greenwhich grids to be treated properly



program nc_2_SH
  use shtlbx
  use shtools
  use grdtlbx
  use netcdf
  implicit none
  integer::i,j,lmax,n,l,m,h,rectind,ndot,nind,start,loopend,ncid,otyp,lonvar,latvar
  integer:: ndims, nvars, natts, unlimdimid, nformat,varlen,dimlen,nlon,nlat,zvar,pos,status
  integer::dimids(3)
  double precision::step,arg,deltalat,deltalon,step2,time,time2,quant
  double precision, allocatable,dimension(:)::lon,lat,dOhm
  double precision, allocatable, dimension(:)::p,q,slm,clm
  double precision,allocatable,dimension(:,:)::smat,cmat
  character(5)::a,b
  character(150)::fileout,dum,filein,varname,dimname,Zunk
  logical::prexit,stdout,usefft,greenwich
  integer::stderr,sample,lmaxdum,lonpos
  double precision,allocatable::cilm(:,:,:),grd(:,:)
  integer::iargc
  integer::varndim,vartype
  logical::verbose
  integer::zslice
  integer::ind,xdim,ydim
  integer::tdim,tvar,nt
  integer,dimension(3)::pvec,cnt
  double precision,allocatable::times(:)
  double precision::ctime,stime,etime,dt
  logical::nanreplace,isnan
  integer:: nzvardim
  logical::pixel
  double precision:: lonleft,tolerance
  call cpu_time(time)
  !initialise defaults
  h=0
  zvar=-1 !will yield an error when no appropriate netcdf variable was found
  lonvar=-1!same for longitude
  latvar=-1!and latitude
  tvar=-1
  xdim=-1
  ydim=-1
  tdim=-1
  ctime=0.d0
  etime=0.d0
  stime=0.d0
  prexit=.false. !when set to true the program does run soome tests on the grid but exits prematurely
  lmax=250
  fileout='####'
  otyp=4
  stdout=.false.

  usefft=.true.
  stderr=0
  greenwich=.false.
  verbose=.false.
  Zunk='z' !check for 'z' netdf variables
  !get command line arguments
  n=iargc()
  zslice=1
  nanreplace=.false.
  pixel=.false. !pixel or gridline registration
  
  if (n < 1) then
     call help()
  else
     !get inputfile name
     call getarg(n,filein)
     if(filein(1:1).eq. '-')call help()
     do,i=1,n-1
        call getarg(i, dum)
        if(dum(1:1).eq. '-')then
           select case(dum(2:2))
           case('l')!limit maximum degree 
              read(dum(3:),'(I3)')lmax

            case('f') !specify output filename
               fileout=dum(3:)
               stdout=.false.
            case('s')!write to standard output
               stdout=.true.
            case('t')!specify different output type
              read(dum(3:),*)otyp 
           case('n') ! don't use the fft approach of Wieczorek
              usefft=.false.
           case('V')
              verbose=.true.
           case('N')
              nanreplace=.true.
           case('Z')! different Z value name
              ind=index(dum,':')
              if(ind .ne. 0)then
                 Zunk=dum(3:ind-1)
                 read(dum(ind+1:),*)zslice
              else
                 Zunk=dum(3:)
                 zslice=1
              end if
           case default
              call help()
           end select
        end if
     end do

     if(fileout(1:4) .eq. '####' .and. .not. stdout) then !only redefine when no new output file is provided
        fileout=trim(filein)//'.sh'
     end if


  end if




  !Open netcdf file and retrieve step size and nodes
  call handle_err(nf90_open(trim(filein),0,ncid))
  if(verbose)write(stderr,*)'VERBOSE: opening: ',trim(filein)
  !inquire some grid properties
  call handle_err(nf90_inquire(ncid, nDimensions = ndims, nVariables=nvars,&
       nAttributes=natts, unlimitedDimId=unlimdimid))

  !get the lon/lat dimensions
  do,i=1,ndims
     call handle_err(nf90_inquire_dimension(ncid=ncid,dimid=i,len = ind,name=dimname))
     select case(dimname)
        case('x','nx','nlon','lon','Longitude','longitude')
           xdim=i
           nlon=ind
        case('y','ny','nlat','lat','Latitude','latitude')
           ydim=i
           nlat=ind
        case('t','time')
           tdim=i
           nt=ind
           allocate(times(nt))
     end select
  end do


! try to figure out the ID of the latitude and longitude and Z variables
  do,i=1,nvars
     call handle_err(nf90_inquire_variable(ncid=ncid, varid=i, name=varname, xtype=vartype, ndims=varndim))
     select case(varname)
     case('x','lon','longitude','Longitude') !accept as longitude
        if(verbose)write(stderr,*)"VERBOSE: Taking ",trim(varname), " as longitude"
        if(varndim .ne. 1)then
           write(stderr,*)"nc_2_SH ERROR: longitude has more than 1 dimension"
           prexit=.true.
        end if
        lonvar=i
     case('y','lat','latitude','Latitude') !accept as latitude
        if(verbose)write(stderr,*)"VERBOSE: Taking ",trim(varname), " as latitude"
        if(varndim .ne. 1)then
           write(stderr,*)"nc_2_SH ERROR: latitude has more than 1 dimension"
           prexit=.true.
        end if
        latvar=i
     case('t','time')
        if(verbose)write(stderr,*)"VERBOSE: Taking ",trim(varname), " as time"
        tvar=i
     case default ! check for unknown Z
        if(varname .eq. Zunk)then
           if(verbose)write(stderr,*)"VERBOSE: Taking ",trim(varname), "(Variable ",i,") as Z value. slice number:",zslice
           zvar=i
           nzvardim=varndim
           if(varndim >3)then
              write(stderr,*)"nc_2_SH ERROR: ",trim(varname)," has more than 3 dimensions"
              prexit=.true.
           end if
        end if
     end select
  end do


  if( zvar < 0)then
     write(stderr,*)"nc_2_SH ERROR: Z value not found"
     prexit=.true.
  end if

  if( latvar < 0)then
     write(stderr,*)"nc_2_SH ERROR: latitude value not found"
     prexit=.true.
  end if

  if( lonvar < 0)then
     write(stderr,*)"nc_2_SH ERROR: longitude value not found"
     prexit=.true.
  end if

  

!  write(*,*)'fileid',ncid,'ndimensions',ndims,'nvariables'&
!  ,nvars,'nattributes',natts,'id unlim dimension',unlimdimid




  !allocate the dimension vectors
  allocate(lon(nlon),lat(nlat))
  
  !get the longitude and latitude vectors
  call handle_err(nf90_get_var(ncid,lonvar, lon))
  call handle_err(nf90_get_var(ncid,latvar, lat))

  !!check grid: coverage, registration, intervals

  if(minval(lon) < -90)greenwich=.true.

!now convert latitude to colatitude in radians
  lat=(90.d0-lat)*pi/(180.d0)
 !convert longitude to radians
  lon=lon*pi/180.d0

!do some checks and determine step size
  deltalat=abs(lat(nlat)-lat(1))
  deltalon=abs(lon(nlon)-lon(1))

  step=(deltalat)/(nlat-1.d0)
  step2=(deltalon)/(nlon-1.d0)

  if(abs((step - step2)/step) > 1.d-5)then
     write(stderr,*)'grid is not equidistant in latitude and longitude'
     prexit=.true.
  end if
  
  !now check for global coverage
  tolerance=1e3*epsilon(0.d0)
  if((pi-deltalat-tolerance) > step)then
     write(0,*)deltalat, step*180.0/pi, lat(1), lat(nlat)
     write(stderr,*)'this grid does not have global coverage in latitude direction'
     prexit=.true.
  end if

  if((2*pi-deltalon-tolerance)> step)then
     write(stderr,*)'this grid does not have global coverage in longitude direction'
     prexit=.true.
  end if

  if(deltalat > 1.5d0*pi)then
     write(stderr,*)'it looks like this grid has transposed latitude and longitude vectors'
     prexit=.true.
  end if



  !check for pixel or gridline registration (grid must be global)
  if(greenwich)then
     lonleft=-180.d0*pi/(180.d0)
  else
     lonleft=0.d0
  end if
  

  if(lonleft-tolerance < lon(1) .and.  lon(1) < lonleft+tolerance)then
     pixel=.false.
  else if (lonleft+step/2.d0-tolerance < lon(1) .and.  lon(1) < lonleft+step/2.d0+tolerance)then
     pixel=.true.
     write(stderr,*)'Grid has pixel registration',step,step2,nlon,nlat
  else
     write(stderr,*)'Could not determine registration mode of grid'
     write(stderr,*)lonleft,tolerance,step, lon(1)
     prexit=.true.
  end if



  
  !check if input grid satisfies conditions
  if(usefft)then
     if(pixel)then
        write(stderr,*)"Grid needs global gridline registration for FFT approach"
        prexit=.true.
     end if

     if( nlon .eq. nlat)then
        sample=1
     else if( nlon-1 .eq. 2*(nlat-1))then
        sample=2
     else
        write(stderr,*)"Grid is not suited for the FFT approach"
        prexit=.true.
     end if

     if( mod(nlat-1,2) .ne. 0)then
        write(stderr,*)"Nlat-1 must be even for the FFT approach"
        prexit=.true.
     end if

     lmaxdum=(nlat-1)/2-1
     if(lmaxdum < lmax)then
        write(stderr,*)"Grid spacing only allows a maximum degree",lmaxdum ,"for the FFT approach"
        prexit=.true.
     end if
  end if

  
  !get timetags (if present)
  if(tdim>0)call handle_err(nf90_get_var(ncid,tvar, times))
  
  !check which dimension corresponds to which dimension of the

  call handle_err(nf90_inquire_variable(ncid=ncid, varid=zvar, dimids=dimids))
  
  !create a permutation vector for the dimension
  if(dimids(1) .eq. ydim .and. dimids(2) .eq. xdim)then
     pvec(1)=2
     pvec(2)=1
     pvec(3)=3
  else if (dimids(1) .eq. xdim .and. dimids(2) .eq. ydim)then
     pvec(1)=1
     pvec(2)=2
     pvec(3)=3
  else
      write(stderr,*)'nc_2_SH ERROR: dimensions of the gridvariable do not appear to fit'
      prexit=.true.
  end if

  if(dimids(3).eq. tdim)then !slice represents a epoch: extract time
     if(zslice .eq. 1 .and. nt>1)then
        dt=times(zslice+1)-times(zslice)
     else if(nt>1)then
        dt=times(zslice)-times(zslice-1)
     end if
     ctime=times(zslice)
     stime=ctime-dt/2
     etime=ctime+dt/2
  else !try looking for a global variable
          
     status=nf90_get_att(ncid=ncid,varid=NF90_GLOBAL,name="Time_year", values=ctime)
     if(status .eq. NF90_NOERR)then !only when found
        etime=ctime
        stime=ctime
     end if
  end if
 




!exit prematurely when checking errors occurred
  if(prexit)call help()

!allocate space for the spherical harmonic coefficients
   pos=SH_pos(lmax,lmax)
   allocate(clm(pos),slm(pos))
   clm=0.d0
   slm=0.d0
   cnt(3)=zslice
   if(usefft)then

      !read data from netcdf grid
      allocate(grd(nlat-1,nlon-1))
      do,i=1,nlat-1!loop over latitude
         cnt(pvec(2))=i
         do,j=1,nlon-1 !loop over longitude
            !get gridvalue from netcdf file
            cnt(pvec(1))=j            
            if(greenwich)then
               if(j>nlon/2)then
                  lonpos=j-nlon/2
               else
                  lonpos=j+nlon/2
               end if
            else
               lonpos=j
            end if


            call handle_err(nf90_get_var(ncid=ncid,varid=zvar,values=grd(nlat-i+1,lonpos),start=cnt))
            ! if(swop)then
            !    status=nf90_get_var(ncid=ncid,varid=zvar,values=grd(nlat-i+1,lonpos),start=(/i,j,zslice/))
            ! else
            !    status=nf90_get_var(ncid=ncid,varid=zvar,values=grd(nlat-i+1,lonpos),start=(/j,i,zslice/))
            ! end if
         end do
      end do
      
      if(nanreplace)then
         do,i=1,nlat-1
            do,j=1,nlon-1
               if(isnan(grd(i,j)))grd(i,j)=0.d0
            end do
         end do
      end if

           !allocate SH matrix
      allocate(cilm(2,lmax+1,lmax+1))
      call SHExpandDH(grid=grd, n=nlat-1, cilm=cilm, lmax=lmaxdum, sampling=sample,lmax_calc=lmax)
      
      !put coefficients in the right vectors
      do,l=0,lmax
         i=SH_pos(l,0)
         !degree 0 term ( cosine only)
         clm(i)=cilm(1,l+1,1)
         do,m=1,l
            i=SH_pos(l,m)
            clm(i)=cilm(1,l+1,m+1)
            slm(i)=cilm(2,l+1,m+1)
         end do
      end do
      deallocate(grd,cilm)
   else ! old approach ( not preffered)
      !allocate spherical harmonic variables and initialize
      allocate(p(pos))
      
      clm=0.d0
      slm=0.d0
      p=0.d0
      
      !precompute cosine and sine terms to accelerate computation
      
      allocate(dOhm(nlat),cmat(pos,nlon),smat(pos,nlon))
      do,l=0,lmax
         do,m=0,l
            !pos=SH_pos(l,m)
            cmat(SH_pos(l,m),:)=cos(m*lon(1:nlon))
            smat(SH_pos(l,m),:)=sin(m*lon(1:nlon))
         end do
      end do
      
      !precompute dOmega
      dOhm=sin(lat)*step**2
      cnt(3)=zslice
      
      !now apply the integral equation to convert to spherical harmonic degrees
      
   
      !step=4*step
      do,i=1,nlat!loop over latitude
         !write(*,*)'Processing latitude',90.d0-lat(i)*180./pi
         cnt(pvec(2))=i
         
         !get the 4pi normalized legendre polynomial for the colatitude
         
         call PlmBar(p=p, lmax=lmax, z=cos(lat(i)))
         do,j=1,nlon !loop over longitude
            if(.not. pixel .and. j==nlon)exit !prevent counting the right most longitude value twice for a gridline registred grid
            cnt(pvec(1))=j            
            !get gridvalue from netcdf file
            !if(i .eq. 1) write(*,*)lon(j)*180/pi
            call handle_err(nf90_get_var(ncid=ncid,varid=zvar,values=quant,start=cnt))
            if(nanreplace .and. isnan(quant))cycle !assume it's zero
            ! if(swop)then
            !    status=nf90_get_var(ncid=ncid,varid=zvar,values=quant,start=(/i,j/))
            ! else
            !    status=nf90_get_var(ncid=ncid,varid=zvar,values=quant,start=(/j,i/))
            ! end if
            !divide the quantity by 4pi
            quant=quant/(4*pi)
            ! if(lon(j)*180.d0/pi <35 .and. lon(j)*180.d0/pi>30.d0 .and. &
            !                &               lat(i)*180.d0/pi > 55.d0 .and. lat(i)*180.d0/pi < 58.d0  )write(*,*)lon(j)*180/pi,90-lat(i)*180/pi,quant
            !cycle longitude loop if gridpoint is zero
            if(abs(quant) < 1.d-100)cycle
            ! 			do,l=0,lmax
            ! 				do,m=0,l
            !        				pos=SH_pos(l,m)
            
            ! 				 if (m .ne. 0) then !only calculate when m not equal to zero
            ! 				 slm(pos)=slm(pos)+quant*p(pos)*sin(m*lon(j))*sin(lat(i))*step**2
            ! 				 end if
            ! 				 clm(pos)=clm(pos)+quant*p(pos)*cos(m*lon(j))*sin(lat(i))*step**2
            
            ! 				end do
            ! 			end do
            !fortran 90 style
            clm=clm+p*quant*cmat(:,j)*dOhm(i)
            slm=slm+p*quant*smat(:,j)*dOhm(i)
            if((i==1 .or. i==nlat) .and. .not. pixel)cycle !only add pole values once in the integration for a gridline registration!!!
         end do
         
         
         
      end do
      
   end if
   
!close netcdf data
call handle_err(nf90_close(ncid))

!write(*,*)'test'
!now write the spherical harmonics to a file
 if (stdout) then
    call SH_write(clm=clm,slm=slm,typ=otyp,tstart=stime,tcent=ctime,tend=etime)
 else
   call SH_write(clm=clm,slm=slm,typ=otyp,filen=trim(fileout))
 end if

call cpu_time(time2)
if(verbose)write(stderr,*)'VERBOSE: computation took',time2-time,'seconds'


end program nc_2_SH
	
subroutine help()
  integer::stderr
  stderr=0
  write(stderr,*)'Convert a netcdf file to 4 pi normalized spherical harmonics (writes to standard output)'
  write(stderr,*)'Usage: nc_2_SH [OPTIONS] inputfile'
  write(stderr,*)'Where inputfile is the name of the netcdf grid'
  write(stderr,*)'The following OPTIONS are supported: ' 
  write(stderr,*)' -fOUTPUTFILE: Specify a different name for the output file'
  write(stderr,*)'               Where OUTPUTFILE is the possible outputfilename (default appends .sh to input file name)'
  write(stderr,*)''
  write(stderr,*)'-s: write to standard output'
  write(stderr,*)'-tNUM: Specify the output type NUM:'
  write(stderr,*)'       1:GRACE(GFZ,JPL,CSR) type'
  write(stderr,*)'       2:GINS type (French GRGS solution)'
  write(stderr,*)'       3:ICGEM type'
  write(stderr,*)'       4:(default) short type first line contains the degree supported, start, center end end time of observation'
  write(stderr,*)'         the other lines are SH coefficients'
  write(stderr,*)'       0:clean no header only SH coefficients'
  write(stderr,*)''
  write(stderr,*)'-ZVARNAME:ZSLICE: takes a different netcdf variable (named VARNAME) for the Z value.'
  write(stderr,*)'                   ZSLICE is the number of the slice (starting at 1).'
  write(stderr,*)'-V verbose'
  write(stderr,*)'-N: replace NaN by zeros'
  write(stderr,*)'-lLMAX: specify a different maximum degree LMAX (default is 250)'
  write(stderr,*)'-n: Use a simple integration without FFT'
  write(stderr,*)'NOTE: The netcdf file must be a equidistant and defined globally'
     	
    	
  stop
end subroutine help

!!subroutine to write a latitude longitude grid to a netcdf file
!!the netcdf file will be in COARDS convention
!!latitude and longitude must be in degrees
!!Coded by Roelof Rietbroek, Wed Jun 20 16:38:41 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Wed Jun 27 15:10:51 2007
!!incorporated also optional x and y attributes (name units)

!!Updated by Roelof Rietbroek, Tue Oct  2 10:38:50 2007
!!added time tag to global atttributes

!!Updated by Roelof Rietbroek, Thu Oct  8 17:02:12 2009
!! added the optional pixel/gridline registration 

!!Updated by Roelof Rietbroek, Fri Feb 12 14:33:36 2010
!! allow appending data to an existing dataset ( along the unlimited dimension)

!!Updated by Roelof Rietbroek, Tue Feb 16 11:50:50 2010
!also allow gmt switch (creates a 2d grid without unlimited dimension)

!!Updated by Roelof Rietbroek, Wed Feb 16 17:26:40 2011
!!fixed bug with non-initialized grd2d

!!Updated by Roelof Rietbroek, Wed Mar 27 21:07:27 2013
!! added a reading routine for 'standard' geographical grids


subroutine NC_write(file,lat,lon,z,history,zunits,zname,xunits,&
     xname,yunits,yname,ctime,pix,napp,nclse)
use netcdf
implicit none
character(*),intent(in)::file
character(*),optional,intent(in)::history,zunits,zname,xunits,xname,yunits,yname
double precision,dimension(:),intent(in)::lat,lon
double precision,optional,intent(in)::ctime
logical,intent(in),optional::pix ! use pixel registration
double precision,dimension(:,:),intent(in)::z
logical,intent(in),optional::nclse,napp !optional arguments allowing appending grids to an existing dataset
double precision::lonmax,lonmin,latmax,latmin,zmax,zmin,dlon,dlat
double precision::zrange(2)
integer::nlat,nlon,status,i,j,nx,ny
integer::dimidx,dimidy,varidx,varidy,varidz
logical::swop
logical::iclose,iappend ! internal values of nclse and napp
logical::exi,grd2d
integer::node_offset
integer,save::ncid=-1
integer::dimidunl,varidtime
integer,save::nunl=0
integer::nlon2,nlat2,stderr,zdims
double precision::ctime2

stderr=0
swop=.false.
node_offset=0
iclose=.true.
iappend=.false.
grd2d=.false.

if(present(nclse))iclose=nclse
if(present(napp))iappend=napp

!if iappend is false and iclose is true make a simple 2d grid
if(iclose .and. .not. iappend)grd2d=.true.

!get the length of the latitude and longitude
if(present(pix))then
   if(pix)node_offset=1 ! use pixel registration ( used primarily in GMT)
end if

nlat=size(lat,1)
nlon=size(lon,1)

nx=size(z,1)
ny=size(z,2)

if(nlat .eq. nx .and. nlon .eq. ny) swop=.true. 
!if(nx==ny)swop=.true.
!get maximum and minimum values assumes monotonically (de)increase
latmin=min(lat(1),lat(nlat))
latmax=max(lat(1),lat(nlat))

lonmin=min(lon(1),lon(nlon))
lonmax=max(lon(1),lon(nlon))

if(node_offset .eq. 1)then ! adapt lon and lat ranges
   dlon=(lonmax-lonmin)/(nlon-1)
   lonmin=lonmin-dlon/2
   lonmax=lonmax+dlon/2

   dlat=(latmax-latmin)/(nlat-1)
   latmin=latmin-dlat/2
   latmax=latmax+dlat/2
end if

zmin=z(1,1)
zmax=z(1,1)

if(swop) then
   do,i=1,nlat
      do,j=1,nlon
         zmin=min(zmin,z(i,j))
         zmax=max(zmax,z(i,j))
      end do
   end do
else

   do,i=1,nlat
      do,j=1,nlon
         zmin=min(zmin,z(j,i))
         zmax=max(zmax,z(j,i))
      end do
   end do

end if



if(ncid .eq. -1)then ! (re)open existing netcdf file in define mode
   !inquire whther file exists
   inquire(FILE=file,EXIST=exi)
   if(exi .and. iappend)then ! reopen file
      call handle_err(nf90_open(path=file, mode=NF90_WRITE, ncid=ncid))
   else ! create new file ( possibly overwriting old one)
      !!open new netcdf grid in define mode
      call handle_err(nf90_create(path=file,cmode=NF90_CLOBBER,ncid=ncid))
      if(iappend)iappend=.false. ! make sure that netcdf is completely initialized not appended to
   end if
end if

!!set global attributes
!!set history
if(iappend)then !check for redefinitions and whether dimensions agree
   !get id of necessary variables

   !time
   call handle_err(nf90_inq_varid(ncid, 'time', varidtime))
   
   !z
   call handle_err(nf90_inq_varid(ncid, 'z', varidz))

   
   !retrieve dimid from lon,lat and time (unlimited)
   !x
   call handle_err( nf90_inq_dimid(ncid, "x", dimidx))
   
   !y
   call handle_err(nf90_inq_dimid(ncid, "y", dimidy))

   !time
   call handle_err(nf90_inquire(ncid=ncid,unlimitedDimId=dimidunl))
   
   !get dimension values
   !x
   call handle_err(nf90_inquire_dimension(ncid, dimidx,len=nlon2))
   !y
   call handle_err(nf90_inquire_dimension(ncid, dimidy,len=nlat2))

   !time (unlimited but get current value)
   call handle_err(nf90_inquire_dimension(ncid=ncid, dimid=dimidunl,len = nunl))
      
   !check wheter new dimensions agree with old one 
   if(nlon2 .ne. nlon .or. nlat2 .ne. nlat)then
      write(stderr,*)"ERROR in NC_write: number of lon/lat nodes has changed"
      stop
   end if

   !check whether grid has 3 dimensions ( required by append)

   call handle_err(nf90_inquire_variable(ncid=ncid, varid=varidz,ndims=zdims))
   if(zdims .ne. 3)then
      write(0,*)"ERROR in NC_WRITE: cannot append grid to this netcdf file: ",trim(file)
      stop
   end if
   !retrieve current zmin and zmax
   call handle_err(nf90_get_att(ncid, varidz, 'actual_range', zrange))

   !update minimum and maximum of z
   zrange(1)=min(zmin,zrange(1))
   zrange(2)=max(zmax,zrange(2))
   call handle_err(nf90_put_att(ncid, varidz,  'actual_range', zrange))

   !retrieve current Ctime
   call handle_err(nf90_get_att(ncid=ncid,varid=NF90_GLOBAL,name='Time_year', values=ctime2))

   !calculate new ctime (average)
   ctime2=(ctime+nunl*ctime2)/(nunl+1)

   !update Ctime
   call handle_err(nf90_put_att(ncid=ncid,varid=NF90_GLOBAL,name='Time_year', values=ctime2))

else ! set up new  netcdf file
   if(present(history))then
      call handle_err(nf90_put_att(ncid=ncid,varid=NF90_GLOBAL,name='history', values=trim(history)))
   else
      call handle_err(nf90_put_att(ncid=ncid,varid=NF90_GLOBAL,name='history', values='unknown'))
   end if
   !!set convention tag
   call handle_err(nf90_put_att(ncid=ncid,varid=NF90_GLOBAL,name='Conventions', values='COARDS'))

   !!set node_registration (grid line or pixel registration)
   call handle_err(nf90_put_att(ncid=ncid,varid=NF90_GLOBAL,name='node_offset', values=node_offset))

   if(present(ctime))then
      ! write(*,*)time
      call handle_err(nf90_put_att(ncid=ncid,varid=NF90_GLOBAL,name='Time_year', values=ctime))
   end if
   
   
   !!setup dimensions
   
   !longitude
   call handle_err(nf90_def_dim(ncid=ncid, name='x', len=nlon, dimid=dimidx))

   
   !latitude
   call handle_err(nf90_def_dim(ncid=ncid, name='y', len=nlat, dimid=dimidy))

   if( .not. grd2d)then
      !unlimited dimension(time)
      call handle_err(nf90_def_dim(ncid, "time", nf90_unlimited, dimidunl))
      nunl=0 ! 0 records yet
   end if



   !define variables

   !longitude
   call handle_err(nf90_def_var(ncid=ncid, name='x', xtype=NF90_DOUBLE, dimids=dimidx, varid=varidx))


   !latitude
   call handle_err(nf90_def_var(ncid=ncid, name='y', xtype=NF90_DOUBLE, dimids=dimidy, varid=varidy))


   !time
   if(present(ctime) .and. .not. grd2d)then
      call handle_err(nf90_def_var(ncid=ncid, name='time', xtype=NF90_DOUBLE, dimids=dimidunl, varid=varidtime))
   end if

   !z-variable
   if(grd2d)then
      if(swop)then
         call handle_err(nf90_def_var(ncid=ncid, name='z', xtype=NF90_DOUBLE, dimids=(/dimidy,dimidx/), varid=varidz))
      else
         call handle_err(nf90_def_var(ncid=ncid, name='z', xtype=NF90_DOUBLE, dimids=(/dimidx,dimidy/), varid=varidz))
      end if
   else
      if(swop)then
         call handle_err(nf90_def_var(ncid=ncid, name='z', xtype=NF90_DOUBLE, dimids=(/dimidy,dimidx,dimidunl/), varid=varidz))
      else
         call handle_err(nf90_def_var(ncid=ncid, name='z', xtype=NF90_DOUBLE, dimids=(/dimidx,dimidy,dimidunl/), varid=varidz))
      end if
   end if

   !write(*,*)status,NF90_NOERR,nx,ny,nlon,nlat
   
   
   !set attributes to variables
   
   
   !longitude
   if(present(xunits) .and. present(xname))then
      call handle_err(nf90_put_att(ncid=ncid,varid=varidx,name='long_name', values=trim(xname)))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidx,name='units', values=trim(xunits)))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidx,name='actual_range', values=(/lonmin,lonmax/)))
   else
      call handle_err(nf90_put_att(ncid=ncid,varid=varidx,name='long_name', values='Longitude'))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidx,name='units', values='degrees_east'))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidx,name='actual_range', values=(/lonmin,lonmax/)))
   end if
   
   !latitude
   if(present(yunits) .and. present(yname))then
      call handle_err(nf90_put_att(ncid=ncid,varid=varidy,name='long_name', values=trim(yname)))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidy,name='units', values=trim(yunits)))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidy,name='actual_range', values=(/latmin,latmax/)))
   else
      call handle_err(nf90_put_att(ncid=ncid,varid=varidy,name='long_name', values='Latitude'))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidy,name='units', values='degrees_north'))
      call handle_err(nf90_put_att(ncid=ncid,varid=varidy,name='actual_range', values=(/latmin,latmax/)))
   end if

   !z variable
   if(present(zunits))then
      call handle_err(nf90_put_att(ncid=ncid,varid=varidz,name='units', values=zunits))
   else
      call handle_err(nf90_put_att(ncid=ncid,varid=varidz,name='units', values='unknown'))
   end if
   
   if(present(zname))then
      call handle_err(nf90_put_att(ncid=ncid,varid=varidz,name='long_name', values=zname))
   else
      call handle_err(nf90_put_att(ncid=ncid,varid=varidz,name='long_name', values='unknown'))
   end if

   call handle_err(nf90_put_att(ncid=ncid,varid=varidz,name='actual_range', values=(/zmin,zmax/)))

   !!now put the actual data in the variables
   !!end the define mode
   call handle_err(nf90_enddef(ncid))

   !longitude
   call handle_err(nf90_put_var(ncid=ncid, varid=varidx, values=lon))

   !latitude
   call handle_err(nf90_put_var(ncid=ncid, varid=varidy, values=lat))
   
end if !end netcdf definition

if(grd2d)then
   !write grid values
   call handle_err(nf90_put_var(ncid=ncid,varid=varidz,values=z,start=(/1,1/)))
else
   nunl=nunl+1 ! increase record number before writing
   if(present(ctime))then ! append time
      call handle_err(nf90_put_var(ncid=ncid, varid=varidtime, values=ctime, start=(/nunl/)))
   end if
   
   !write (or append) grid values
   call handle_err(nf90_put_var(ncid=ncid,varid=varidz,values=z,start=(/1,1,nunl/)))
end if
if(iclose)then
   !!close netcdf file
   call handle_err(nf90_close(ncid))
   ncid=-1 ! set ncid to nonsense value (indicates file is closed)
end if

end subroutine NC_write

!! subroutine NC_read reads in a geographical grid slice from a netcdf file
!! output grid has dimension z2d(nlon,nlat)
!! input arrays must be pointers and may be unassociated
subroutine NC_read(file,lon,lat,z2d,ctime,etime,stime,zname,zslice)
use netcdf
!use grdtlbx
implicit none
character(*),intent(in)::file
double precision,pointer,dimension(:),intent(inout)::lon,lat
double precision,pointer,dimension(:,:),intent(inout)::z2d
double precision,optional,intent(out)::ctime,etime,stime
character(*),intent(in),optional::zname
integer,intent(in),optional::zslice
!local variables
integer::ncid,ndims,nvars,natts,unlimdimid,status
integer::i,j,nlon,nlat,nt,ind,xdim,ydim,varndim,vartype,tvar
integer::lonvar,latvar,izslice,stderr,tdim,zvar,znd
character(150)::dimname,varname,Zunk
double precision,allocatable::times(:),work(:,:)
logical::prexit,trans
integer,dimension(3)::cnt,dimids,strt
double precision::dt
stderr=0
dimids=0


if(present(zslice))then
   izslice=zslice
else
   izslice=1
end if

prexit=.false.

if(present(zname))then
   Zunk=trim(zname)
else
   Zunk='z'
end if

!Open netcdf file and retrieve step size and nodes
call handle_err(nf90_open(trim(file),0,ncid))

!inquire some grid properties
call handle_err(nf90_inquire(ncid, nDimensions = ndims, nVariables=nvars,&
     nAttributes=natts, unlimitedDimId=unlimdimid))

  !get the lon/lat dimensions
  do,i=1,ndims
     call handle_err(nf90_inquire_dimension(ncid=ncid,dimid=i,len = ind,name=dimname))
     select case(dimname)
        case('x','nx','nlon','lon')
           xdim=i
           nlon=ind
        case('y','ny','nlat','lat')
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
     case('x','lon','longitude') !accept as longitude
        if(varndim .ne. 1)then
           write(stderr,*)"nc_2_SH ERROR: longitude has more than 1 dimension?"
           prexit=.true.
        end if
        lonvar=i
     case('y','lat','latitude') !accept as latitude
        if(varndim .ne. 1)then
           write(stderr,*)"nc_2_SH ERROR: latitude has more than 1 dimension"
           prexit=.true.
        end if
        latvar=i
     case('t','time')
        tvar=i
     case default ! check for unknown Z
        if(varname .eq. Zunk)then
           zvar=i
           if(varndim .eq.2)znd=varndim
           if(varndim >3  )then
              write(stderr,*)"NC_READ: ",trim(varname)," has too many dimensions"
              prexit=.true.
           end if
           if(varndim < 2  )then
              write(stderr,*)"NC_READ: ",trim(varname)," has too few dimensions"
              prexit=.true.
           end if
        end if
     end select
  end do


  if( zvar < 0)then
     write(stderr,*)"NC_read ERROR: Z value not found"
     prexit=.true.
  end if

  if( latvar < 0)then
     write(stderr,*)"NC_read ERROR: latitude value not found"
     prexit=.true.
  end if

  if( lonvar < 0)then
     write(stderr,*)"NC_read ERROR: longitude value not found"
     prexit=.true.
  end if



  !allocate the dimension vectors
  if(.not. associated(lon))allocate(lon(nlon))
  if(.not. associated(lat))allocate(lat(nlat))
  
  !get the longitude and latitude vectors
  call handle_err(nf90_get_var(ncid,lonvar, lon))
  call handle_err(nf90_get_var(ncid,latvar, lat))
  
  !get timetags (if present)
  if(tdim>0)call handle_err(nf90_get_var(ncid,tvar, times))
  
  !check which dimension corresponds to which dimension of the
  call handle_err(nf90_inquire_variable(ncid=ncid, varid=zvar, dimids=dimids(1:znd)))

  !check whether the grid's tranpose is stored
  if(dimids(1) .eq. ydim .and. dimids(2) .eq. xdim)then
     trans=.true.
  else if (dimids(1) .eq. xdim .and. dimids(2) .eq. ydim)then
     trans=.false.
  else
      write(stderr,*)'NC_read ERROR: dimensions of the gridvariable do not appear to fit'
      prexit=.true.
  end if

  if(znd >2 .and. dimids(3).eq. tdim)then !slice represents a epoch: extract time
     if(izslice .eq. 1 .and. nt>1)then
        dt=times(izslice+1)-times(izslice)
     else if(nt>1)then
        dt=times(izslice)-times(izslice-1)
     end if
     if(present(ctime))ctime=times(izslice)
     if(present(stime))stime=times(izslice)-dt/2
     if(present(etime))etime=times(izslice)+dt/2
  else !try looking for a global variable
     if(present(ctime))then
        status=nf90_get_att(ncid=ncid,varid=NF90_GLOBAL,name="Time_year", values=ctime)

        if(status .eq. NF90_NOERR)then !only when found
           if(present(etime))etime=ctime
           if(present(stime))stime=ctime
        end if
     end if
  end if

 if(.not. associated(z2d))allocate(z2d(nlon,nlat))

 if(trans)then
    strt(3)=izslice
    strt(1)=1
    strt(2)=1
    cnt(1)=nlat
    cnt(2)=nlon
    cnt(3)=1
    allocate(work(nlat,nlon))
    call handle_err(nf90_get_var(ncid=ncid,varid=zvar,values=work,start=strt(1:znd),count=cnt(1:znd)))
    z2d=transpose(work)
    deallocate(work)
 else
    strt(3)=izslice
    strt(1)=1
    strt(2)=1
    cnt(1)=nlon
    cnt(2)=nlat
    cnt(3)=1
    call handle_err(nf90_get_var(ncid=ncid,varid=zvar,values=z2d,start=strt(1:znd),count=cnt(1:znd)))
 end if
 
  if(prexit)then
     write(stderr,*)"NC_read too many ERRORS: aborting"
     stop
  end if


  !close netcdf data
  call handle_err(nf90_close(ncid))
  if(allocated(times))deallocate(times)
  
end subroutine NC_read

subroutine handle_err(status)
  use netcdf
  integer, intent (in) :: status
  integer::stderr
  stderr=0
     
  if(status /= nf90_noerr) then
     write(stderr,*)trim(nf90_strerror(status))
     stop "Stopped"
  end if
end subroutine handle_err

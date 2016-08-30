!subroutine which returns a rotation matrix which converts between cartesian coordinates and local north east up coordinates
!the type vector contains on entry strings with X, Y, Z entries and on exit it will contain the entries of the output vector in the order H,E,N
!example: 
! input vector type and station correspond to the sides of the rotation matrix
!index:   type:  station (ID): pos:
! 1       STAX    STAT1    x(STAT1)
! 2       STAY    STAT1    y(STAT1)
! 3	  STAZ    STAT1    z(STAT1)  
! ..      STAX    STAT2    x(STAT2)
!
!the order of the above input rows maybe changed (the output matrix will change accordingly)

!in the typpe parameter only the CAPTITALS X,Y and Z are checked, avoid ambiguities, When neither X,Y or Z is present an error message is issued 

!The output is contained in the matrix rotmat and the corresponding entries(for 1 side) in type
!so X becomes H
! Y becomes E
! Z becomes N
! stations remain the same


!!Coded by Roelof Rietbroek, Thu Aug 16 13:49:22 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

subroutine GPS_rotmat(type,station,pos,rot)
implicit none
character(*),dimension(:),intent(inout)::type
character(*),dimension(:),intent(in)::station
double precision,dimension(:),intent(in)::pos
double precision,dimension(:,:), intent(out)::rot

!no optional arguments
integer::i,j,ix,iy,iz,xind,yind,zind,arg,ndat,stderr
logical,allocatable,dimension(:)::done
character(len(type(1))),allocatable,dimension(:)::newtype
double precision::cir,mag,rx,ry,rz

!initialize
rot=0.d0
stderr=0

!first do a few dimension checks
ndat=size(pos,1)
allocate(done(ndat),newtype(ndat))
done=.false.


do,i=1,(ndat) !loop to search for stations
   if(done(i))cycle !cycle when a row is encountered which is already done (from a previous station)
   
   !loop to gather necessary x,y,z coordinates for the particular station
   arg=0
   do,j=1,ndat
      if(index(station(j),trim(station(i))) .eq. 0) cycle !cycke when a different station is encountered
      ix=index(type(j),'X')!index of the letter in the type string
      iy=index(type(j),'Y')
      iz=index(type(j),'Z')
      if(ix > 0)then
         arg=arg+1
         rx=pos(j)
         xind=j
         
         newtype(j)=type(j)(1:ix-1)//'H'//type(j)(ix+1:)
         
      else if(iy > 0)then
         arg=arg+1
         ry=pos(j)
         yind=j
         newtype(j)=type(j)(1:iy-1)//'E'//type(j)(iy+1:)
      else if(iz > 0)then
         arg=arg+1
         rz=pos(j)
         zind=j
         newtype(j)=type(j)(1:iz-1)//'N'//type(j)(iz+1:)
      else
         write(stderr,*)'Type must contain X,Y or Z'
         stop
      end if
      if(arg .eq. 3)exit !all three components collected: exit loop
   end do
   
   !construct rotation matrix for station
   cir=sqrt(rx*rx+ry*ry)
   mag=sqrt(rx*rx+ry*ry+rz*rz)
   !first row
   rot(xind,xind)=rx/mag
   rot(xind,yind)=ry/mag
   rot(xind,zind)=rz/mag
   !second row
   rot(yind,xind)=-ry/cir
   rot(yind,yind)=rx/cir
   rot(yind,zind)=0.d0
   !third row
   rot(zind,xind)=-rot(xind,zind)*rot(yind,yind)
   rot(zind,yind)=rot(xind,zind)*rot(yind,xind)
   rot(zind,zind)=cir/mag

   !tag entries to be done
   done(xind)=.true.
   done(yind)=.true.
   done(zind)=.true.

   
end do

!replace type list with new one
type=newtype

end subroutine GPS_rotmat

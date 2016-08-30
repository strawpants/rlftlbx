!!Coded by Roelof Rietbroek, Thu Dec  3 15:28:58 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!two subroutine to convert ECEF coordinates to geodetic lon,lat height coordinates and vice versa
!!uses the WGS84 ellipsoid

subroutine ECEF_2_geodetic(x,y,z,lon,lat,h)
implicit none
double precision,intent(in)::x,y,z
double precision,intent(out)::lon,lat,h
double precision,save::ell_a,ell_b,ell_e,ell_et,pi
logical,save :: first=.true.
double precision::p,th,Rn
if(first)then
   ! WGS84 ellipsoid
   ell_a = 6378137; !semi major axis (in meters)
   ell_e = 6.69437999014d-3 ! eccentricity squared
   
   !derived values
   ell_b=sqrt(ell_a**2*(1-ell_e))!semi minor axis
   ell_et=(ell_a**2-ell_b**2)/ell_b**2 !second eccentrencity squared
   pi=acos(-1.d0)   
   first=.false.
end if

!longitude is exactly know
lon=atan2(y,x)
if(lon < 0)lon=2*pi+lon ! shift longitude to 0..2pi

p=sqrt(x**2+y**2)
th=atan2(ell_a*z,ell_b*p) !geocentric latitude
lat = atan2((z+ell_et*ell_b*sin(th)**3),(p-ell_e*ell_a*cos(th)**3)); !geodetic latitude
Rn=ell_a/sqrt(1-ell_e*(sin(lat))**2)

if(p<1)then ! very close to the pole ( use approximation Rn~=semi-minor axis)
   h=abs(z)-ell_b
else
   h=p/cos(lat)-Rn
end if

end subroutine ECEF_2_geodetic


subroutine geodetic_2_ECEF(lon,lat,h,x,y,z)
implicit none
double precision,intent(in)::lon,lat,h
double precision,intent(out)::x,y,z
double precision::ell_a,ell_e,Rn
! WGS84 ellipsoid
ell_a = 6378137; !semi major axis (in meters)
ell_e = 6.69437999014d-3 ! eccentricity squared

Rn=ell_a/sqrt(1-ell_e*(sin(lat))**2)

x=(Rn+h)*cos(lat)*cos(lon)
y=(Rn+h)*cos(lat)*sin(lon)
z=(Rn*(1-ell_e)+h)*sin(lat)

end subroutine geodetic_2_ECEF

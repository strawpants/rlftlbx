!!Coded by Roelof Rietbroek, Tue Jun 25 19:41:14 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

! subroutine which computes the ocean function in from the present day topography and optionally the Ice thickness and a sea level change (change from the present day topography)
!the topography is given on a grid while the Ice and Sea level change are provided as Spherical harmonic coefficients

!input
! topop: present day topography in meter Nx2N
! lon(2N), lat(N): corresponding longitudes and latitudes 
! lmax: maximum degree and order of the input ice and Sea 
! Tice: thickness of the ice sheet ( in spherical harmonic coefficients)
! Sea: sea level change relative to the present day value (in SH)

!output
! Oce: Ocean function in spherical harmonics up to degree and order 2*lmax

!the SH ordering follows the function SH_tpos
!conversion from the grid to SH is done using an FFT approach

module GIA_func
implicit none
contains
subroutine GIA_Ocean_func(topo,lmax,Tice,Sea,Oce,SHS)
use shtlbx
use shtools
use netcdf
use grdtlbx
use forttlbx
use sh_synthesis
implicit none
integer,intent(in)::lmax
double precision,intent(in),optional::Tice(:,:),Sea(:,:)
type(geogrid),intent(inout)::topo
double precision,intent(inout)::Oce(:,:)
type(SHsynth_struct),intent(inout)::SHS

!private variables
integer::nsh,i,j,l,m,k,pos
double precision,allocatable::mask(:,:),ocetmp(:,:,:)
double precision,allocatable,dimension(:)::clmtmp,slmtmp
double precision,pointer,dimension(:,:)::tmpgrd=>null()
logical::fex
integer::stderr,ind,lmaxdum,nt
integer,allocatable,dimension(:)::deg,ord,tri,iclm,islm
double precision::coscolat,grounded
stderr=0
nt=0
nsh=SH_tpos(2*lmax,2*lmax,1,2*lmax,0)
if(size(Oce,1) .ne. nsh)then
   write(stderr,*)"ERROR: input vector Oce has the wrong dimensions"
   stop
end if
   

! read in topography file when not done already
if( .not. associated(topo%z))then

   inquire(FILE=trim(topo%file),EXIST=fex)
   if(fex)then
      call NC_read(file=topo%file,lon=topo%lon,lat=topo%lat,z2d=topo%z)
      topo%nlat=size(topo%lat,1)
      topo%nlon=size(topo%lon,1)
   else
      write(stderr,*)"ERROR: file not found:",trim(topo%file)
      stop
   end if

   !dimension check
   if(topo%nlon .ne. size(topo%z,1) .or. topo%nlat .ne. size(topo%z,2))then
      write(stderr,*)"ERROR: topography grid does not fit the longitude and latitude vectors"
      stop
   end if

   !check whether the grid is usable for the FFT approach
   if(topo%nlon .ne. 2*topo%nlat)then
      write(stderr,*)"ERROR: topography grid is not suited for the FFT approach"
      write(stderr,*)"ERROR: expecting 2NxN points"
      stop
   end if
   
   if(lmax > topo%nlat/2-1)then
      write(stderr,*)"ERROR: topography grid requires more than",2*(lmax+1)," latitude points"
      stop
   end if

   write(stderr,*)topo%lon(1),topo%lat(1),topo%z(1,1)
   write(stderr,*)topo%lon(topo%nlon),topo%lat(topo%nlat),topo%z(topo%nlon,topo%nlat)
end if




!make a mask from the topography
allocate(mask(topo%nlat,topo%nlon))

!find out number of time tags
if(present(Sea))nt=size(Sea,2)
if(present(Tice) .and. nt >0)then
   if(size(Tice,2) .ne. nt)then !start crying about inconsistent size
      write(stderr,*)"ERROR GIA_Ocean_func: Ice thickness and sea level change differ in size"
      stop
   end if
else if (present(Tice))then
   nt=size(Tice,2)
end if

if(nt == 0)nt=1 !reset to 1 if no optional goodies are provided

!stack the Ocean function in the output vector

allocate(ocetmp(2,lmax*2+1,lmax*2+1)) !allocate temporary array for the Spherical harmonics

!possibly apply restrictions based on the additional sea level and/or ice thickness
if(present(Sea) .or. present(Tice))then
   
   allocate(tmpgrd(topo%nlat,topo%nlon))
   if(SHS%nlon .eq. 0)then !initialize synthesis structure
      call SH_synth_init(lmax=lmax,lon=topo%lon,lat=topo%lat,SHS=SHS,save_pnm=.true.)
   end if
   
   !make index vectors to map SH representations
   pos=SH_pos(lmax,lmax)
   allocate(iclm(pos),islm(pos),clmtmp(pos),slmtmp(pos))
   
   do,l=0,lmax
      do,m=0,l
         ind=SH_tpos(l,m,0,lmax,0)
         iclm(SH_pos(l,m))=ind
         if(m>0)then
            ind=SH_tpos(l,m,1,lmax,0)
         else
            ind=0
         end if
         islm(SH_pos(l,m))=ind
      end do
   end do
end if

do,i=1,nt !loop over all SHsets (may be parallelized)
   !copy topography (transpose and and flip latitude order)
   forall(j=1:topo%nlon,k=1:topo%nlat)mask(topo%nlat-k+1,j)=topo%z(j,k)
   
   if(present(Sea))then
      
      !copy data in other SH vectors (different format)
      forall(j=1:pos)clmtmp(j)=Sea(iclm(j),i)
      forall(j=1:pos,islm(j)>0)slmtmp(j)=Sea(islm(j),i)
      !do the synthesis on a grid
      call SH_dosynth(SHS=SHS,clm=clmtmp,slm=slmtmp,GRD=tmpgrd)
      
      !update the mask grid
      mask=mask-tmpgrd
   end if
   
   
   if(present(Tice))then ! check for marine grounded ice
      !copy data in other SH vectors (different format)
      forall(j=1:pos)clmtmp(j)=Tice(iclm(j),i)
      forall(j=1:pos,islm(j)>0)slmtmp(j)=Tice(islm(j),i)
      !do the synthesis on a grid
      call SH_dosynth(SHS=SHS,clm=clmtmp,slm=slmtmp,GRD=tmpgrd)
      
      !check for marine grounded ice
      where(tmpgrd > -mask*rho_w/rho_ice)mask=mask+tmpgrd*rho_ice/rho_w !artificially add the ice thickness to the topography where it is grounded
   end if
   
   where(mask >= 0)mask=0 ! make the Ocean function on the grid   
   where(mask < 0)mask=1 ! make the Ocean function on the grid


!compute the spherical harmonic coefficients up to degree 2*lmax
   !compute SH coefficeints using FFT
   call SHExpandDH(grid=mask, n=topo%nlat, cilm=ocetmp, lmax=lmaxdum, sampling=2,lmax_calc=lmax*2)   
   
   !copy values at the appropriate place
   do,l=0,2*lmax
      do,m=0,l
         ind=SH_tpos(l,m,0,lmax*2,0)
         Oce(ind,i)=ocetmp(1,l+1,m+1)
         if(m>0)then
            ind=SH_tpos(l,m,1,lmax*2,0)
            Oce(ind,i)=ocetmp(2,l+1,m+1)
         end if
      end do
   end do

end do





!deallocate goodies
deallocate(ocetmp,mask)
if(present(Sea) .or. present(Tice))deallocate(iclm,islm,tmpgrd,clmtmp,slmtmp)

end subroutine GIA_Ocean_func
end module GIA_func

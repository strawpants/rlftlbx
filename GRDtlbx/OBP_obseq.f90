!!assembly of subroutines to construct observation equation associated with ocean bottom pressure

!!Coded by Roelof Rietbroek, Thu Feb 14 13:35:25 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!adapted from GPS_obseq2

!!Updated by Roelof Rietbroek, Tue May 13 15:09:44 2008
!!fixed a bug in  OBP_obseq_globalmean
!!The design matrix was only set to 1 for the first measurement point 



!creates an design matrix which relates point measurements of bottom pressure to spherical harmonic coefficients of equivalent water height
subroutine OBP_obseq_loadsh(A,tag,lon,lat,lmax,lmin)
use SHtlbx
use shtools
implicit none
double precision, intent(out)::A(:,:)
character(24),intent(out),optional::tag(:)
double precision,intent(in)::lon(:),lat(:)
integer,intent(in)::lmax,lmin

!private variables
integer::pos,nc,inds,indc,l,m,i,ndat
integer,allocatable,dimension(:)::Asvec,Acvec,pcvec,psvec
double precision::latold,lonold
double precision,allocatable,dimension(:)::p,cosmlon,sinmlon,deg,ord

!initializations




ndat=size(lon,1) !number of data points

latold=-999.d0 !initial values guaranteed to work on first entry
lonold=-999.d0


!construct a spherical harmonic index vector for vectors p and dp (ass. legendre polynomials)
pos=SH_pos(lmax,lmax)

allocate(p(pos))
p=0.d0



!precompute vectors
allocate(cosmlon(pos),sinmlon(pos))
!degree and order vectors
allocate(deg(pos),ord(pos))

!allocate index vectors
nc=SH_pos(lmax,lmax)-SH_pos(lmin-1,lmin-1) !amount of cosine coefficients
!cosine vector
allocate(Acvec(nc),pcvec(nc))
!sine vector (no zero order terms)
allocate(Asvec(nc-lmax+lmin-1),psvec(nc-lmax+lmin-1))!amount of sine coefficients

inds=0
indc=0
do,i=1,pos
   call SH_lm(pos=i,l=l,m=m)
   !make degree and order vectors
   deg(i)=l
   ord(i)=m

!below only  l>=lmin and m .ne. 0 (for Slm coeff)
   if(l<lmin)cycle
   
   indc=indc+1
   Acvec(indc)=SH_tpos(l,m,0,lmax,lmin)
   pcvec(indc)=i
  
   if(present(tag))write(tag(Acvec(indc)),'(a3,1x,I3,i3,10x,a4)')'TCN',l,m,'OBPM'
   
   if(m .ne. 0)then

      inds=inds+1
      Asvec(inds)=SH_tpos(l,m,1,lmax,lmin)
      psvec(inds)=i
      if(present(tag))write(tag(Asvec(inds)),'(a3,1x,I3,i3,10x,a4)')'TSN',l,m,'OBPM'
   end if
end do







do,i=1,ndat
   
   !calculate legendre polynomials and cos and sin terms if not different from previous step (can be speed advantage when provided smartly on input)
   if(latold .ne. lat(i))then 
      latold=lat(i)

      !call SHTOOLS routine to calculate associated legendre polynomials and its derivative
     ! write(*,*)pos,size(p),size(dp),lmax
!      call PlmBar_d1(p=p,dp=dp,lmax=lmax,z=cos(pi/2-latold))

      call PlmBar(p=p,lmax=lmax,z=sin(latold))



   end if

   if(lonold .ne. lon(i))then !only calculate if lon is different from the previous step
      lonold=lon(i)
      sinmlon=sin(ord*lonold)
      cosmlon=cos(ord*lonold)
   end if

   A(i,Acvec)=cosmlon(pcvec)*p(pcvec)
   A(i,Asvec)=sinmlon(psvec)*p(psvec)
end do



end subroutine OBP_obseq_loadsh



!subroutine to set up design matrix for GPS stations for the geocenter loading part (alternate way to define degree 1 deformation)
subroutine OBP_obseq_geocent(A,tag,lon,lat)
use SHtlbx
implicit none
double precision,intent(out)::A(:,:)
character(24),intent(out),optional::tag(3)
double precision,intent(in),dimension(:)::lon,lat

double precision::geo2deg1
double precision,dimension(3)::eh
integer::i,ndat

!geocenter to mass load conversion

geo2deg1=rho_e/(rho_w*sqrt(3.d0)) 

!geo2deg1=(rho_e/(rho_w*sqrt(3.d0)))*(1-(1+0.021+0.268))*rho_w/rho_e 




ndat=size(lat,1)

!make unknowns tag
if(present(tag))then
   tag(1)='SXB'
   tag(2)='SYB'
   tag(3)='SZB'
end if



do,i=1,ndat



!unit vector in up direction
      eh(1)=cos(lat(i))*cos(lon(i))
      eh(2)=cos(lat(i))*sin(lon(i))
      eh(3)=sin(lat(i))
      
      A(i,1)=geo2deg1*eh(1)
      A(i,2)=geo2deg1*eh(2)
      A(i,3)=geo2deg1*eh(3)
      
end do




end subroutine OBP_obseq_geocent

subroutine OBP_obseq_globalmean(A,tag)
  implicit none
  double precision,intent(out)::A(:,:)
  character(24),intent(out),optional::tag(:)

  if(present(tag))tag='T_GLB_BIAS          OBPM'
  
  A(:,1:1)=1.d0

end subroutine OBP_obseq_globalmean


!observation equation which relates the pressure to the dynamic part of the sea level
!thus static equilibrium sea level is subtracted
subroutine OBP_obseq_loadsh_equiref(A,tag,lon,lat,lmax,lmin)
use SHtlbx
use shtools
implicit none
double precision, intent(out)::A(:,:)
character(24),intent(out),optional::tag(:)
double precision,intent(in)::lon(:),lat(:)
integer,intent(in)::lmax,lmin

!private variables
integer::pos,nc,inds,indc,l,m,i,ndat,FR,ltyp
integer,allocatable,dimension(:)::Asvec,Acvec,pcvec,psvec
double precision::latold,lonold
double precision,allocatable,dimension(:)::p,cosmlon,sinmlon,deg,ord,hnm,knm,fact

!initializations
FR=3
ltyp=2



ndat=size(lon,1) !number of data points

latold=-999.d0 !initial values guaranteed to work on first entry
lonold=-999.d0


!construct a spherical harmonic index vector for vectors p and dp (ass. legendre polynomials)
pos=SH_pos(lmax,lmax)

allocate(p(pos))
p=0.d0



!precompute vectors
allocate(cosmlon(pos),sinmlon(pos))
!degree and order vectors
allocate(deg(pos),ord(pos))

!allocate index vectors
nc=SH_pos(lmax,lmax)-SH_pos(lmin-1,lmin-1) !amount of cosine coefficients
!cosine vector
allocate(Acvec(nc),pcvec(nc))
!sine vector (no zero order terms)
allocate(Asvec(nc-lmax+lmin-1),psvec(nc-lmax+lmin-1))!amount of sine coefficients


inds=0
indc=0
do,i=1,pos
   call SH_lm(pos=i,l=l,m=m)
   !make degree and order vectors
   deg(i)=l
   ord(i)=m

!below only  l>=lmin and m .ne. 0 (for Slm coeff)
   if(l<lmin)cycle
   
   indc=indc+1
   Acvec(indc)=SH_tpos(l,m,0,lmax,lmin)
   pcvec(indc)=i
  
   if(present(tag))write(tag(Acvec(indc)),'(a3,1x,I3,i3,7x,a7)')'TCN',l,m,'DYNOBPM'
   
   if(m .ne. 0)then

      inds=inds+1
      Asvec(inds)=SH_tpos(l,m,1,lmax,lmin)
      psvec(inds)=i
      if(present(tag))write(tag(Asvec(inds)),'(a3,1x,I3,i3,7x,a7)')'TSN',l,m,'DYNOBPM'
   end if
end do


!load love numbers
allocate(hnm(pos),knm(pos),fact(pos))
call SH_loadlove(hnm=hnm,knm=knm,typ=ltyp,frame=FR)
!construct factors
fact=0.d0
do,i=1,pos
   fact(i)=(3.d0*rho_w/rho_e)*(1.d0+knm(i)-hnm(i))/(2*deg(i)+1)
   fact(i)=1.d0-fact(i)
   write(*,*)deg(i),fact(i)
end do




do,i=1,ndat
   
   !calculate legendre polynomials and cos and sin terms if not different from previous step (can be speed advantage when provided smartly on input)
   if(latold .ne. lat(i))then 
      latold=lat(i)

      !call SHTOOLS routine to calculate associated legendre polynomials and its derivative
     ! write(*,*)pos,size(p),size(dp),lmax
!      call PlmBar_d1(p=p,dp=dp,lmax=lmax,z=cos(pi/2-latold))

      call PlmBar(p=p,lmax=lmax,z=sin(latold))



   end if

   if(lonold .ne. lon(i))then !only calculate if lon is different from the previous step
      lonold=lon(i)
      sinmlon=sin(ord*lonold)
      cosmlon=cos(ord*lonold)
   end if
   
   A(i,Acvec)=cosmlon(pcvec)*p(pcvec)*fact(pcvec)
   A(i,Asvec)=sinmlon(psvec)*p(psvec)*fact(psvec)
end do


end subroutine OBP_obseq_loadsh_equiref

!!module SHtlbx contains interfaces and constants used in GRACE processing and Spherical armonic analysis
!!Coded by Roelof Rietbroek, Thu May 23 14:59:08 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Thu May 31 14:59:25 2007
!!Updated by Roelof Rietbroek, Thu Jan 21 11:58:16 2010
!Added Mean Earth rotation
!!Updated by Roelof Rietbroek, Fri Feb  5 18:13:15 2010
!added moments of inertia and rotfeedback interface
!added c20_replace function

!!Updated by Roelof Rietbroek, Sun Jan 20 16:39:29 2013
!! added C20 Stokes coefficient as a constant from GRS80
!!Updated by Roelof Rietbroek, Mon Jun 23 16:43:13 2014
!! fixed Initialization of the principal moments of inertia


module SHtlbx
implicit none
save
!mean Earth radius in m y
double precision::RE = 0.6378136460d+07
!gravitational constant of the earth 
double precision::GM = 0.3986004415d+15

!gravitational universal constant
double precision:: GU = 6.6738480d-11

!density of sea water in kg/m^3
double precision, parameter:: rho_w =1.025d+3
!mean density of the Earth
double precision, parameter :: rho_e = 5517.0d+0

!density of (compacted) ice
double precision,parameter ::rho_ice = 931.d0 !kg/m^3 taken from G. Spada and friends

!mean acceleration of the Earth
double precision,parameter::g=9.80665d0

!ratio of the oceans surface to the total surface (taken from ocean.sh)
!double precision,parameter::A_oce=0.70797655106184
double precision,parameter::A_oce=7.0802102076271d-1

!mean earth rotation (in radians per second) relative to the stars 
double precision,parameter::OHM_E=7.292115467064d-5
double precision,parameter::GRS80_C20=-4.84166854896120d-04

!Mean moments of interita from the earth ( kg m^2)
!Updated 23 June 2014 : Based on a comment of B. Uebbing, changed E+37 to d+37 to properly initiate double values 
double precision,parameter::Ixx_A=8.0102d+37  !~= (I_xx+Iyy)/2
double precision,parameter::Izz_C=8.0365d+37

!Chandler frequency ( radians per second) from iers constants
double precision,parameter::Ohm_chandler=1.679064144d-7

!pi
!double precision,parameter::pi=2*acos(0.d0)
double precision,parameter::pi=3.1415926535897932384&
&6264338327950288419716939937510d0

interface
    subroutine SH_gaus(clm,radius)
     double precision,dimension(:), intent(inout) :: clm
     double precision, intent(in)::radius
   end subroutine SH_gaus

   subroutine SH_gaus_recur(clm,radius)
     double precision,dimension(:), intent(inout) :: clm
     double precision, intent(in)::radius
   end subroutine SH_gaus_recur
   
   subroutine SH_lm(pos,l,m)
     integer, intent(in)::pos
     integer, intent(out)::l,m
   end subroutine SH_lm
   
   function SH_pos(l,m)
     integer, intent(in)::l,m
     integer:: SH_pos
   end function SH_pos
   
   function SH_opos(l,m,lmax)
     implicit none
     integer, intent(in)::l,m,lmax
     integer::SH_opos
   end function SH_opos
   

    function SHgDWpos(l,m,q,lmax,lmin)
        integer*8 SHgDWpos
        integer,intent(in)::l,m,q,lmax,lmin
    end function 

   function SH_tpos(l,m,q,lmax,lmin)
     implicit none
     integer,intent(in)::l,m,lmax,lmin,q
     integer::SH_tpos
   end function SH_tpos

   subroutine SH_tlm(pos,l,m,q,lmax,lmin)
     implicit none
     integer,intent(in)::lmax,lmin,pos
     integer,intent(out)::l,m,q
   end subroutine SH_tlm

   subroutine SH_readgrav(filen,clm,slm,clm_sig,slm_sig,skip,frmt,type,str,clm_tv,slm_tv)
     double precision,dimension(:),optional,intent(out):: clm,slm,clm_sig,slm_sig
     double precision,optional,intent(out),dimension(:,:)::clm_tv,slm_tv ! optional time variable coefficients in file
     integer, intent(in),optional::skip,type
     character(*),optional,intent(in)::frmt,str
     character(*), intent(in),optional::filen
   end subroutine SH_readgrav

!!$   function SH_ftype(filen)
!!$     character(*),intent(in)::filen
!!$     integer::SH_ftype
!!$   end function SH_ftype
!!$   
   subroutine getctime(string,typ,tst,tc,tnd)
     character(*),intent(in)::string
     integer, intent(in)::typ
     double precision,intent(out)::tst,tc,tnd
   end subroutine getctime
   
   subroutine SH_readmeta(filen,type,lmax,mu,rad,tcent,tstart,tend)
     character(*),intent(in),optional::filen
     integer,intent(out),optional::type,lmax
     double precision,optional,intent(out)::mu,rad,tcent,tstart,tend
   end subroutine SH_readmeta
                       
   subroutine SH_write(clm,slm,clm_sig,slm_sig,tstart,tcent,tend,comm1,comm2,comm3,filen,typ)
   double precision,intent(in),optional,dimension(:)::clm,slm,clm_sig,slm_sig
   double precision,intent(in),optional::tstart,tcent,tend
   integer,intent(in),optional::typ
   character(*),intent(in),optional::comm1,comm2,comm3,filen
   end subroutine SH_write

   subroutine SH_geocmod(c10,c11,s11,c30,c31,s31,time,time1,mod,geocfile)
     integer,intent(in)::mod!parameter specifying which model to use
     character(*),intent(in),optional::geocfile
     double precision,intent(out)::c10,c11,s11
     double precision,intent(in)::time !time in years (eg 2000.5 means somewhere half way of the year 2000) (center period if time1 not provided)
     double precision,optional,intent(out)::c30,c31,s31 !affected by the flattening of the Earth
     double precision,optional,intent(in)::time1 !end of the observation period also in years 
   end subroutine SH_geocmod

   subroutine SH_loadlove(hnm,knm,lnm,typ,frame,fileop,iso)
     character(*),intent(in),optional::fileop
     double precision,intent(out),optional,dimension(:)::knm,hnm,lnm
     integer,intent(in)::typ,frame
     logical,intent(in),optional::iso
   end subroutine SH_loadlove

   subroutine SH_loadGIA_uplift_ratio(rat,typ)
     double precision,intent(out),dimension(:)::rat
     integer,intent(in)::typ
   end subroutine SH_loadGIA_uplift_ratio

   subroutine SH_readW(ar,file,indx,lmax,lmin,oldtype,packed)
     double precision,dimension(:,:),intent(out),optional::ar
     integer,intent(inout),optional::lmax,lmin
     integer,intent(out),optional::indx(:,:)
     logical,intent(in),optional::oldtype,packed
     character(*),intent(in)::file
   end subroutine SH_readW

   subroutine SH_tripint_lm(l1,m1,q1,l2,m2,q2,qout,triplm)
     implicit none
     integer,intent(in)::l1,l2,m1,m2,q1,q2,qout
     double precision,intent(out),dimension(:)::triplm
   end subroutine SH_tripint_lm

   function SH_overlapint(l1,m,l2)
     implicit none
     double precision::SH_overlapint
     integer,intent(in)::l1,m,l2
   end function SH_overlapint



   function qfact(top,bot)
     implicit none
     integer,intent(in)::top,bot
     double precision::qfact
   end function qfact

   !rotational feedback functions
   function SurfL_2_Ji3()
     implicit none
     double precision:: SurfL_2_Ji3(3,4)
   end function SurfL_2_Ji3

   function m_2_rotpot()
     implicit none
     double precision::m_2_rotpot(4,3)
   end function m_2_rotpot

   function Ji3_2_m(ltyp)
     implicit none
     double precision::Ji3_2_m(3,3)
     integer,intent(in)::ltyp
   end function Ji3_2_m

   function c20_replace(time,time1,type)
     implicit none
     double precision,intent(in)::time
     double precision,optional,intent(in)::time1
     integer,intent(in)::type
     double precision:: c20_replace
   end function c20_replace

   subroutine SH_Love_d1_trans(h1in,l1in,k1in,h1out,l1out,k1out,FR)
     double precision,intent(in)::h1in,l1in,k1in
     double precision,intent(out)::h1out,l1out,k1out
     character(*),intent(in)::FR
   end subroutine SH_Love_d1_trans

   ! subroutine GIA_Ocean_func(topo,lmax,Tice,Sea,Oce,SHS)
   !   use grdtlbx
   !   use sh_synthesis
   !   integer,intent(in)::lmax
   !   double precision,intent(in),optional::Tice(:,:),Sea(:,:)
   !   type(geogrid),intent(inout)::topo
   !   double precision,intent(inout)::Oce(:,:)
   !   type(SHsynth_struct),intent(inout)::SHS
   ! end subroutine GIA_Ocean_func

end interface

end module SHtlbx

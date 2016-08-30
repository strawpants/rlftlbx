!module SHtlbx contains interfaces and constants used in GRACE processing and Spherical armonic analysis
!!Coded by Roelof Rietbroek, Thu May 23 14:59:08 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Thu May 31 14:59:25 2007
module SHtlbx
implicit none
save
!mean Earth radius in m y
double precision::RE = 0.6378136460d+07
!gravitational constant of the earth 
double precision::GM = 0.3986004415d+15
!density of sea water in kg/m^3
double precision, parameter:: rho_w =1.025d+3
!mean density of the Earth
double precision, parameter :: rho_e = 5517.0d+0

!mean acceleration of the Earth
double precision,parameter::g=9.80665d0

!ratio of the oceans surface to the total surface (taken from ocean.sh)
double precision,parameter::A_oce=7.0802102076271d-1
!pi
double precision,parameter::pi=3.1415926535897932384&
6264338327950288419716939937510d0

!mean earth rotation (in radians per second) relative to the stars 
double precision,parameter::OHM_E=7.292115467064d-5

!Mean moments of interita from the earth ( kg m^2)
double precision,parameter::Ixx_A=8.0101E+37  !~=Iyy
double precision,parameter::Izz_C=8.0365E+37

!Chandler frequency ( radians per second)
double precision,parameter::Ohm_chandler=1.679064144d-7


interface
   subroutine SH_gaus(clm,radius)
     double precision,dimension(:), intent(inout) :: clm
     double precision, intent(in)::radius
   end subroutine SH_gaus
   
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

   subroutine SH_readgrav(filen,clm,slm,clm_sig,slm_sig,skip,frmt,type,str,clm_dot,slm_dot)
     double precision,dimension(:),optional,intent(out):: clm,slm,clm_sig,slm_sig
     double precision,optional,intent(out),dimension(:)::clm_dot,slm_dot ! optional conventional rates in file
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

   subroutine SH_geocmod(c10,c11,s11,time,time1,mod,geocfile)
     integer,intent(in)::mod!parameter specifying which model to use
     character(*),intent(in),optional::geocfile
     double precision,intent(out)::c10,c11,s11
     double precision,intent(in)::time !time in years (eg 2000.5 means somewhere half way of the year 2000) (center period if time1 not provided)
     double precision,optional,intent(in)::time1 !end of the observation period also in years 
   end subroutine SH_geocmod

   subroutine SH_loadlove(hnm,knm,lnm,typ,frame,fileop)
     character(*),intent(in),optional::fileop
     double precision,intent(out),optional,dimension(:)::knm,hnm,lnm
     integer,intent(in)::typ,frame
   end subroutine SH_loadlove

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

end interface

end module SHtlbx

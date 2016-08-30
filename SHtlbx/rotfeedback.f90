!!Coded by Roelof Rietbroek, Wed Feb  3 15:46:27 2010
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

! subroutines to be used for calculating (linearized) rotational feedback

!references:
! Peltier & Luthcke 2009 JGR "on the origins of Earth rotation anomalies"
! Milne & Mitrovica 1997, Geoph J int: Post glacial sea-level change on a rotation Earth"
! Wu & peltier 1984, Geop J. R. Astr Soc., Pleistocene deglaciation and the Earth's rotation a new analysis

! Calculate matrix relating surface loading to Ji3 moments of inertia
! the matrix M gives a linear relation between the surface loading coefficients ( in eqh) and changes in moment of inertia
!                 / C00 \
! / J13  \       |  C20  |
!|  J23   |  = M |  C21  |  
! \ J33  /        \ S21 /
!equations from Wu and peltier 1984 converted to real 4pi normalized (geodesy) spherical harmonics
! Note the equations ( 43) from Milne 1997 appear to be erronuous (depending on whether the Condon-Shortley phase is applied)
!!Updated by Roelof Rietbroek, Wed Sep 22 20:53:44 2010
!! corrected for proper normalization (was inversely applied)
 

function SurfL_2_Ji3()
use shtlbx,only : pi,RE,rho_w
implicit none
double precision:: SurfL_2_Ji3(3,4)

!initialize
SurfL_2_Ji3=0.d0

! !J13
! SurfL_2_Ji3(1,3)=-4.d0/5.d0*pi*(RE**4)*sqrt(6.d0/5.d0)*rho_w

! !J23
! SurfL_2_Ji3(2,4)=-4.d0/5.d0*pi*(RE**4)*sqrt(6.d0/5.d0)*rho_w


! !J33
! SurfL_2_Ji3(3,1)=8.d0/3.d0*pi*(RE**4)*rho_w
! SurfL_2_Ji3(3,2)=-8.d0/3.d0*pi*(RE**4)*rho_w/(5*sqrt(5.d0))

!J13
SurfL_2_Ji3(1,3)=-4.d0/5.d0 *sqrt(10.d0/6.d0)*pi*(RE**4)*rho_w

!J23
SurfL_2_Ji3(2,4)=-4.d0/5.d0 *sqrt(10.d0/6.d0)*pi*(RE**4)*rho_w


!J33
SurfL_2_Ji3(3,1)=8.d0/3.d0*pi*(RE**4)*rho_w
SurfL_2_Ji3(3,2)=-8.d0/3.d0*pi*(RE**4)*rho_w/(sqrt(5.d0))

end function SurfL_2_Ji3

!subroutine providing a matrix to relate the change in polar motion to the change in rotational potential (linearized)
!from Peltier 2009 ( and corrected fro the sign -m -> m for coefficient 21

      
!  / C_r00 \          / m1 \
! |  C_r20  |  =  A  |  m2  |
! |  C_r21  |         \ m3 /
!  \ S_r21 / 

!the Earths rotation vector:

!   / ohm1 \            /  m1  \
!  |  ohm2  | = OHM_E  |   m2   |
!   \ ohm3 /            \ 1+m3 /


function m_2_rotpot()
use SHtlbx,only:OHM_E,RE
implicit none
double precision::m_2_rotpot(4,3)

m_2_rotpot=0.d0

m_2_rotpot(1,3)=((OHM_E*RE)**2)*2.d0/3.d0

m_2_rotpot(2,3)=-2.d0/3.d0*((OHM_E*RE)**2)/sqrt(5.d0)

m_2_rotpot(3,1)=-((OHM_E*RE)**2)/sqrt(15.d0)


m_2_rotpot(4,2)=-((OHM_E*RE)**2)/sqrt(15.d0)




end function m_2_rotpot


!subroutine to create a matrix which relates the rigid pertubations of the moments of inertia to polar wander
!references Peltier 2009, Milne 1997
!Notes:
! only the elastic part is used!
! Chandler wobble is empirical (making fluid body love number empirical as well) 
         
!    / m1 \           / J13 \ 
!   |  m2  |  =  B   |  J23  |
!    \ m3 /           \ J33 /
             

function Ji3_2_m(ltyp)
use SHtlbx, only:OHM_E,Ixx_A,Izz_C,Ohm_chandler,SH_loadlove
implicit none
double precision::Ji3_2_m(3,3)
integer,intent(in)::ltyp
double precision::knm(6),k2
!load elastic k degree 2 load love numbers
!write(0,*)knm,ltyp
call SH_loadlove(knm=knm,typ=ltyp,frame=1) !not: the frame of reference is irrelevant
k2=1.d0+knm(4)

Ji3_2_m=0.d0

Ji3_2_m(1,1)=OHM_E*k2/(Ixx_A*Ohm_chandler)

Ji3_2_m(2,2)=Ji3_2_m(1,1)

Ji3_2_m(3,3)=-k2/Izz_c


end function Ji3_2_m

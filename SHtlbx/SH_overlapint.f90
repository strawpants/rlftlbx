! function to calculate the special overlap integral of two associated (4pi normalized) Legendre functions:
!                     / +1
! K(l1,m,l2) = int   |      P      (x) P      (x)  dx
!                   / -1     l1,m-1     l2,m+1




!using 1 sum of wigner 3j symbols
!references:
! @ARTICLE{mavromatis:1999,
!    author = {{Mavromatis}, H.~A.},
!     title = "{A single-sum expression for the overlap integral of two associated Legendre polynomials  }",
!   journal = {Journal of Physics A Mathematical General},
!      year = 1999,
!     month = apr,
!    volume = 32,
!     pages = {2601-2603},
!    adsurl = {http://adsabs.harvard.edu/abs/1999JPhA...32.2601M},
!   adsnote = {Provided by the SAO/NASA Astrophysics Data System}
! }
! @ARTICLE{Dong:2002,
!    author = {{Dong}, S.-H.},
!     title = "{COMMENT:  Comment on `A single-sum expression for the overlap integral of two associated Legendre polynomials' }",
!   journal = {Journal of Physics A Mathematical General},
!      year = 2002,
!     month = may,
!    volume = 35,
!     pages = {4187-4188},
!    adsurl = {http://adsabs.harvard.edu/abs/2002JPhA...35.4187D},
!   adsnote = {Provided by the SAO/NASA Astrophysics Data System}
! }

!!Coded by Roelof Rietbroek, Thu May  8 10:41:47 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

function SH_overlapint(l1,m,l2)
use SHTOOLS
!use SHtlbx, only: qfact
implicit none
double precision::SH_overlapint
integer,intent(in)::l1,m,l2

!private variables 
double precision::wign3j(l1+l2+1),wign3j0(l1+l2+1)
double precision::scale,tmp
integer::k,jmin,jmax,jmin0,top1,top2,bot1,bot2

!initialization
SH_overlapint=0.d0


!selection rules


!quick return if l1 and l2 have different parity
if(mod(l1,2) .ne. mod(l2,2))return

!or if order is zero
if(m .eq. 0)return

!or if l1 < m-1 or l2< m+1 ! non existent Legendre function
if((m-1 > l1) .or. (m+1 > l2))return

!Unormalized
!scale=(4*((-1)**(m-1)))*sqrt(qfact(l1+m-1,l1-m+1)*qfact(l2+m+1,l2-m-1))

!normalized scale
if(m .eq. 1)then
!   scale=(4*((-1)**(m-1)))*sqrt((2*l1+1.d0))*sqrt(4.d0*l2+2.d0)
   scale=(4*((-1)**(m-1)))*sqrt(2.d0)
else
   scale=(8.d0*((-1)**(m-1)))
!   scale=(4*((-1)**(m-1)))*2.d0*sqrt((2*l1+1.d0))*sqrt(2.d0*l2+1.d0)
end if



!get Wigner 3-j symbols

!call Wigner3j(w3j=wign3j0, jmin=jmin, jmax=jmax, j2=l1, j3=l2, m1=0, m2=0, m3=0)
!call Wigner3j(w3j=wign3j, jmin=jmin, jmax=jmax, j2=l1, j3=l2, m1=-2, m2=-m+1, m3=m+1)

call Wigner3j(w3j=wign3j0, jmin=jmin0, jmax=jmax, j2=l1, j3=l2, m1=0, m2=0, m3=0)
call Wigner3j(w3j=wign3j, jmin=jmin, jmax=jmax, j2=l1, j3=l2, m1=-2, m2=-m+1, m3=m+1)


!calculate the overlap integral


do,k=jmin,jmax
   if(mod(k,2) .eq. 0)then
      ! use integers which behave as O(l^2) ( to prevent too large integers)
      top1=(2*k+1)*(2*l1+1)
      bot1=(k-1)*(k+2)
      top2=(2*k+1)*(2*l2+1)
      bot2=k*(k+1)
      
      tmp=sqrt((top1/dble(bot1))*(top2/dble(bot2)))

!      tmp=sqrt(dble((k-1)*(k+2)))*sqrt(dble((k+1)*k))
!      write(*,*)tmp
      SH_overlapint=SH_overlapint+wign3j0(k-jmin0+1)*&
           wign3j(k-jmin+1)*tmp
   end if
end do

SH_overlapint=SH_overlapint*scale




end function SH_overlapint





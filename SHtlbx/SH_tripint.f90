
!subroutine which calculates the triple integral of three real spherical harmonic base functions (fully normalized)
!The degree order and sine/cosine parameter of the first two base functions are fixed
!the third base function has free degree and orders for which all values are outputted
!the vector goes in triplm
!ordering scheme as from SH_pos(l,m)
!note 
subroutine SH_tripint_lm(l1,m1,q1,l2,m2,q2,qout,triplm)
use shtlbx,only: SH_pos
use shtools
implicit none
integer,intent(in)::l1,l2,m1,m2,q1,q2,qout
double precision,intent(out),dimension(:)::triplm
double precision::wign3(l1+l2+1),wign3a(l1+l2+1),wign0a(l1+l2+1)
integer::m3,m3a,lmin3,lmin3a,lmax3,lmax3a,l
!save parameters
integer,save::l1sv=-99
integer,save::l2sv=-99
integer,save::lmin0,lmax0
double precision,save,allocatable::wign0(:),degfac(:)
double precision::lonint1,lonint2,normfac,scale1,scale2,norm
double precision::sq2
sq2=sqrt(2.d0)
!check index


!initialize
triplm=0.d0

!first selection criteria: determine if qout,q1,q2 contain an odd number of cosines (q=0) else return zero straight away
!write(*,*)'test',mod((qout+q1+q2),2),qout,q1,q2
if(mod((qout+q1+q2),2) .ne. 0)return
 

!second criteria: check whether input contains a zero order sine coefficient request
!this would yield zero over the longitude integral

if((q1 .eq. 1) .and. m1 .eq. 0)return
if((q2 .eq. 1) .and. m2 .eq. 0)return



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!now there are two possibilities for the order
m3=m1+m2
m3a=abs(m2-m1)

!we can quickly return if m3 is zero order and a sine component is requested (qout ==1)
if(qout .eq. 1 .and. m3 .eq. 0)return


!calculate the integral over longitude for each of those two possibilities (m3 and m3a)

!lonint1 corresponding to m3 (always existent)
if( (q1 .eq. 1) .and. (q2 .eq. 1))then
!  lonint1=-4*pi/2.d0
   lonint1=-1.d0/2.d0
else
!   lonint1=4*pi/2.d0
 lonint1=1.d0/2.d0
end if
if(m1 .ne. 0)lonint1=lonint1/2.d0
if(m2 .ne. 0)lonint1=lonint1/2.d0


!now make a scalar (longitude integral) for the second possibility of order |m2-m1|

if(m1 > m2)then
   if(qout .eq. 1 .and. q2 .eq. 1)then
!      lonint2=-4*pi/2.d0
      lonint2=-1.d0/2.d0
   else
!      lonint2=4*pi/2.d0
      lonint2=1.d0/2.d0
   end if
   if(m2 .ne. 0)lonint2=lonint2/2.d0
   if(m3a .ne. 0)lonint2=lonint2/2.d0
else
   if(qout .eq. 1 .and. q1 .eq. 1)then
!      lonint2=-4*pi/2.d0
     lonint2=-1.d0/2.d0
   else
!     lonint2=4*pi/2.d0
      lonint2=1.d0/2.d0
   end if
   if(m1 .ne. 0)lonint2=lonint2/2.d0
   if(m3a .ne. 0)lonint2=lonint2/2.d0
end if



!calculate wigner 3j symbols
! subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3)
!corresponds to
!            /  j  j2 j3 \
!            \  m1 m2 m3 /

!now calculate the zero order wigner-3j symbols

!wlen=l1+l2+1
!allocate(wign0(wlen),wign3(wlen),wign3a(wlen))

!            /  l3  l1 l2 \  ==  / l1  l2  l3 \ 
!            \  0   0   0 /      \ 0   0   0  /

!only calculate wign0 when necessary

!!precomputation of degree dependent part
if(l1sv .ne. l1 .or. l2sv .ne. l2)then
!   write(0,*)'calculating wign0',l1,l2
   if(allocated(wign0))deallocate(wign0)
   allocate(wign0(l1+l2+1))
   call Wigner3j(wign0, lmin0, lmax0, l1, l2, 0, 0,0)

   !precalculate square roots
   if(allocated(degfac))deallocate(degfac)
   allocate(degfac(0:l1+l2))
   do,l=abs(l1-l2),l1+l2,2 ! with step 2 
  !    if(mod(l1+l2+l,2) .eq. 1)write(0,*)'ERROR1'
      degfac(l)=sqrt((2*l+1.d0)*(2*l1+1.d0)*(2*l2+1.d0))
   end do
   l1sv=l1
   l2sv=l2

end if

!call wigner 3j symbol for first valid order



!use symmetry property
      !/ l1   l2  l3 \  ==    /  l3 l1 l2 \ (which is inputted, yields more values for l3)
      !\ m1   m2 -m3 /        \ -m3 m1 m2 /

call Wigner3j(wign3, lmin3, lmax3, l1, l2, -m3, m1, m2)



!call wigner 3j symbol for second valid order
if(m3a .ne. m3 )then


   if(m1>m2)then
      !/ l3     l2   l1 \  == / l2  l1  l3   \
      !\ m1-m2  m2  -m1 /     \ m2 -m1  m1-m2/
    
      call Wigner3j(wign3a, lmin3a, lmax3a, l2, l1, m3a, m2, -m1)

   else if(m1<m2)then
      !/ l3     l1   l2 \ == / l1  l2  l3   \ 
      !\ m2-m1  m1  -m2 /    \ m1 -m2 m2-m1 /
    
      call Wigner3j(wign3a, lmin3a, lmax3a, l1, l2, m3a, m1, -m2)
  !    !apply symmetry factor
   !   if(mod(m2,2) .ne. 0)wign3a=-1*wign3a
   else if( qout .ne. 1)then !don't calculate when m3a is zero and sinevector is requested
    
      call Wigner3j(wign3a, lmin3a, lmax3a, l1, l2, m3a, m1, -m2)
   end if
 
end if

! write(*,*)0,'0',wign0
! write(*,*)m3,'3',wign3
! write(*,*)m3a,'3a',wign3a
! write(*,*)lmin0,lmax0,lmin3,lmax3,lmin3a,lmax3a
!now do a loop over the degrees (most restrictive)
!first order

! !test: write Clebsch Gordon coefficients to screen
! write(*,*)'clebsh-Gordon for l1,0,l2,0,0',l1,0,l2,0,0
! do,l=lmin0,lmax0
!    write(*,*)l1,l2,l,0,0,0,(-1)**(l1-l2)*sqrt(2*l+1.)*wign0(l-lmin0+1)
! end do

! write(*,*)'clebsh-Gordon for l1,m1,l2,m2,m3',l1,m1,l2,m2,m3
! do,l=lmin3,lmax3
!    write(*,*)l1,l2,l,m1,m2,m3,(-1)**(m3+l1-l2)*sqrt(2*l+1.)*wign3(l-lmin3+1)
! end do


! write(*,*)'clebsh-Gordon for l1,m1,l2,m2,m3a',l1,m1,l2,m2,m3a
! do,l=lmin3a,lmax3a
!    if(m1>m2)then
!       write(*,*)l,l2,l1,m3a,m2,m3a+m2,(-1)**(m3a+m2+l-l2)*sqrt(2*l1+1.)*wign3a(l-lmin3a+1)
!    else
!       write(*,*)l,l1,l2,m3a,m1,m3a+m1,(-1)**(m3a+m1+l-l1)*sqrt(2*l2+1.)*wign3a(l-lmin3a+1)
!    end if
!end do


!calculate the scale factor to apply to wigner3j symbols (conversion from wigner to clebsch gordan)


scale1=2.d0*(1-2*mod(m3,2))


if(m1 >= m2)then
   scale2=2.d0*(1-2*mod(m1,2))!mirror if m1 is odd
!    scale2=2.d0*(-1)**m1 !mirror if m1 is odd
else
   scale2=2.d0*(1-2*mod(m2,2)) !mirror if m2 is odd
!    scale2=2.d0*(-1)**m2 !mirror if m2 is odd
end if





norm=1.d0
if(m1 .ne. 0)norm=norm*sq2
if(m2 .ne. 0)norm=norm*sq2
if(m3 .ne. 0)norm=norm*sq2
! if(m1 .ne. 0)norm=norm*2.d0
! if(m2 .ne. 0)norm=norm*2.d0
! if(m3 .ne. 0)norm=norm*2.d0

do,l=lmin3+mod(l1+l2+lmin3,2),lmax3,2

   !calculate normalization factor
!   normfac=sqrt(norm*(2*l+1.d0)*(2*l1+1.d0)*(2*l2+1.d0))
 !  if(mod(l1+l2+l,2) .eq. 1)write(0,*)'ERROR2'
   triplm(SH_pos(l,m3))=scale1*norm*degfac(l)*lonint1*wign0(l-lmin0+1)*wign3(l-lmin3+1)
end do


!do the same for the second order

if(m3a .eq. 0 .and. qout .eq. 1)return

if(m3a .ne. m3)then
   norm=1.d0
   if(m1 .ne. 0)norm=norm*sq2
   if(m2 .ne. 0)norm=norm*sq2
   if(m3a .ne. 0)norm=norm*sq2
!    if(m1 .ne. 0)norm=norm*2.d0
!    if(m2 .ne. 0)norm=norm*2.d0
!    if(m3a .ne. 0)norm=norm*2.d0
   do,l=lmin3a+mod(l1+l2+lmin3a,2),lmax3a,2

      !calculate normalization factor
!      normfac=sqrt(norm*(2*l+1.d0)*(2*l1+1.d0)*(2*l2+1.d0))
!      if(mod(l1+l2+l,2) .eq. 1)write(0,*)'ERROR3'
      triplm(SH_pos(l,m3a))=scale2*norm*degfac(l)*lonint2*wign0(l-lmin0+1)*wign3a(l-lmin3a+1)
   end do
end if


end subroutine SH_tripint_lm

!!test program for tripint subroutine
!program output can be compared to analytical values from blewitt 2003
!note that a normalization transformation is applied
! program test_tripint
! use SHtlbx
! implicit none
! character(3)::q1s,q2s
! integer::q1,q2,l1,l2,m1,m2,pos2,i,l3,m3
! double precision,allocatable,dimension(:)::tripYC,tripYS,tripYCi,tripYSi
! double precision::conv,conv2
! q1=0
! l1=1 
! m1=1

! q2=0
! l2=1
! m2=0




! if(q1 .eq. 0)then
!    q1s='COS'
! else
!    q1s='SIN'
! end if

! if(q2 .eq. 0)then
!    q2s='COS'
! else
!    q2s='SIN'
! end if





! pos2=SH_pos((l1+l2),(l1+l2))
! write(*,*)pos2,mod(0,2)
! allocate(tripYC(pos2),tripYS(pos2),tripYCi(pos2),tripYSi(pos2))


! call SH_tripint_lm(l1,m1,q1,l2,m2,q2,0,tripYC)
! call SH_tripint_lm(l1,m1,q1,l2,m2,q2,1,tripYS)

! call SH_tripint_lm(l2,m2,q2,l1,m1,q1,0,tripYCi)
! call SH_tripint_lm(l2,m2,q2,l1,m1,q1,1,tripYSi)
! do,i=1,pos2
!    call SH_lm(i,l3,m3)
!    !undo normalization to check with output from Blewitt 2003
!    conv=sqrt(qfact(l3-m3,l3+m3)*(2*l3+1)*qfact(l2-m2,l2+m2)*(2*l2+1))
! !   conv=sqrt(qfact(l1-m1,l1+m1)*qfact(l2-m2,l2-m2)*qfact(l3-m3,l3+m3)*(2*l1+1)*(2*l2+1)*(2*l3+1))
!   ! if(m1 .ne. 0)conv=conv*sqrt(2.d0)
!    if(m2 .ne. 0)conv=conv*sqrt(2.d0)
!    if(m3 .ne. 0)conv=conv*sqrt(2.d0)
!    conv2=sqrt(qfact(l1-m1,l1+m1)*(2*l1+1))
!    if(m1 .ne. 0)conv2=conv2*sqrt(2.d0)
! !   write(*,*)conv2,conv,qfact(0,4)
!   write(*,*)q1s,l1,m1,q2s,l2,m2,l3,m3,'COS',tripYC(i)*conv2/conv,'SIN',tripYS(i)*conv2/conv
! !   write(*,*)q1s,l1,m1,q2s,l2,m2,l3,m3,'SIN',tripYC(i),tripYCi(i),'COS',tripYS(i),tripYSi(i)
! !   write(*,*)q1s,l1,m1,q2s,l2,m2,l3,m3,'COS',tripYC(i)*sqrt(qfact(l1+m1,l1-m1))*conv2/conv,&
! !'SIN',tripYS(i)*sqrt(qfact(l1+m1,l1-m1))*conv2/conv
! end do


! end program test_tripint


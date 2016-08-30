!!Coded by Roelof Rietbroek, Mon Sep  7 16:14:09 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!! the module!!
!faster but more memory intensive calculation of triple SH basefunction integrals
module SH_trpnt
  !derived type to hold a vector ( varying size per entry)
  type degvec
     double precision,pointer::scale(:)=>null() !nullify pointer
  end type degvec
  
  type trpnt_memstruct
     !arrays with precomputed order dependent scales
     double precision,pointer::scale3a(:,:)=>null()
     double precision,pointer::scale3b(:,:)=>null()
     !array with precomputed degree dependent scales
!     type(degvec),allocatable::degarr(:)
     type(degvec),pointer::degarr(:)
     
  end type trpnt_memstruct
  
contains
  !return valid order a
  function trp_m3a(m1,m2)
    implicit none
    integer::trp_m3a
    integer,intent(in)::m1,m2
    trp_m3a=m1+m2
  end function trp_m3a

  !return valid order b
  function trp_m3b(m1,m2)
    implicit none
    integer::trp_m3b
    integer,intent(in)::m1,m2
    trp_m3b=abs(m1-m2)
  end function trp_m3b


  function trp_qfaca(q1,q2,q3,m1,m2)
    implicit none
    double precision::trp_qfaca
    integer,intent(in)::q1,q2,q3,m1,m2

    integer::m3

    m3=trp_m3a(m1,m2)
    
    if(mod((q1+q2+q3),2) .ne. 0)then !quick return with 0 value for longitude integral
       trp_qfaca=0.d0
       return
    else if((q1==1 .and. m1==0) .or. (q2==1 .and. m2==0) .or. (q3==1 .and. m3==0))then
       !zero order sine terms present => return 0
       trp_qfaca=0.d0
       return
    else if(q1+q2 .eq. 2)then !two sine terms will result in a negative longitude integral
       trp_qfaca=-1.d0
       return
    else !positive longitude integral
       trp_qfaca=1.d0
       return
    end if
    
  end function trp_qfaca


  function trp_qfacb(q1,q2,q3,m1,m2)
    implicit none
    double precision::trp_qfacb
    integer,intent(in)::q1,q2,q3,m1,m2

    integer::m3
    m3=trp_m3b(m1,m2)
    
    if(mod((q1+q2+q3),2) .ne. 0)then !quick return with 0 value for longitude integral
       trp_qfacb=0.d0
       return
    else if((q1==1 .and. m1==0) .or. (q2==1 .and. m2==0) .or. (q3==1 .and. m3==0))then
       !zero order sine terms present => return 0
       trp_qfacb=0.d0
       return
    else if(q1+q2+q3 .ne. 2)then ! + when no two Sine terms are present (always positive)
       trp_qfacb=1.d0
       return
    else if(m1>m2)then ! two sine terms are presetn and m1>m2
       if(q2+q3 .eq. 2)then !two sine terms will result in a negative longitude integral
          trp_qfacb=-1.d0
       else
          trp_qfacb=1.d0
       end if
       return
    else !two sine terms present and m2>=m1
       if(q1+q3 .eq. 2)then !two sine terms will result in a negative longitude integral
          trp_qfacb=-1.d0
       else
          trp_qfacb=1.d0
       end if
       return
    end if
    
  end function trp_qfacb


  function trp_degind(l1,l2)
    implicit none
    integer::trp_degind
    integer,intent(in)::l1,l2
    integer::row,col
    row=min(l1,l2)+1
    col=max(l1,l2)+1
    trp_degind=row+col*(col-1)/2
  end function trp_degind

  subroutine trpnt_vec_lm(l1,l2,m1,m2,m3a,m3b,lmin3a,lmin3b,lmax3a,lmax3b,veca,vecb,trpnt)
    implicit none
    type(trpnt_memstruct),intent(in),optional::trpnt
    double precision,intent(out)::veca(:),vecb(:)
    integer,intent(in)::l1,l2,m1,m2
    integer,intent(out)::m3a,m3b,lmin3a,lmin3b,lmax3a,lmax3b
    double precision,pointer::degvec(:)
    double precision::ordsca,ordscb
!    double precision::wign3a(l1+l2+1),wign3b(l1+l2+1)
    integer::l,ind
    integer::st
    integer::dum
    !initialize vectors
    veca=0.d0
    vecb=0.d0

    !quick return if possible
    if(m1>l1)return
    if(m2>l2)return
    
    ind=trp_degind(l1,l2)
    
    if(present(trpnt))then ! use precomputed values
       degvec=>trpnt%degarr(ind)%scale ! point to the right vector
       ordsca=trpnt%scale3a(m2,m1)
       ordscb=trpnt%scale3b(m2,m1)
    else ! calculate values on the fly
       allocate(degvec(0:l1+l2))
       call trpnt_degscale(l1,l2,degvec)
       ordscb=trpnt_mscaleb(m1,m2)
       ordsca=trpnt_mscalea(m1,m2)
    end if


    !retrieve first valid order
    m3a=trp_m3a(m1,m2)
    
    lmin3a=max(abs(l1-l2),abs(m3a))
    
    !call wigner 3j symbol for first valid order m3a

    !use symmetry property
    !/ l1   l2  l3 \  ==    /  l3 l1 l2 \ (which is inputted, yields varying values for l3)
    !\ m1   m2 -m3 /        \ -m3 m1 m2 /
  !  write(0,*)"wign3a",size(wign3a),l1,l2    
!    call Wigner3j(wign3a, lmin3a, lmax3a, l1, l2, -m3a, m1, m2)


    call Wigner3j_light(veca(lmin3a+1:), dum, lmax3a, l1, l2, -m3a, m1, m2)



    st=lmin3a+mod(l1+l2+lmin3a,2)
    do,l=st,lmax3a,2 ! loop with step 2 ( only over valid parity)
!       veca(l+1)=trpnt%degarr(ind)%scale(l)*trpnt%scale3a(m1,m2)*wign3a(l-lmin3a+1)
!       veca(l+1)=trpnt%degarr(ind)%scale(l)*trpnt%scale3a(m2,m1)*veca(l+1)

!       veca(l+1)=degsc(l)*trpnt_mscalea(m1,m2)*veca(l+1)
       veca(l+1)=degvec(l)*ordsca*veca(l+1)
!        veca(l+1)=trpnt%degarr(ind)%scale(l)*trpnt%scale3a(m2,m1)*veca(l+1)

     end do


    !second valid order
    m3b=trp_m3b(m1,m2)
    
    !check for quick return: m3a==m3b
    if(m3a .eq. m3b )return !nothing to do

    lmin3b=max(abs(l1-l2),abs(m3b))
    !call wigner 3j symbol for second valid order m3b ( if not the same)
    if(m1>=m2)then
       !/ l3     l2   l1 \  == / l2  l1  l3   \
       !\ m1-m2  m2  -m1 /     \ m2 -m1  m1-m2/
!       write(0,*)"wign3b_1",size(wign3a),l1,l2    
!       call Wigner3j(wign3b,lmin3b, lmax3b, l2, l1, m3b, m2, -m1)
       call Wigner3j_light(vecb(lmin3b+1:), dum, lmax3b, l2, l1, m3b,&
            m2, -m1)
       
    else if(m1<m2)then
       !/ l3     l1   l2 \ == / l1  l2  l3   \ 
       !\ m2-m1  m1  -m2 /    \ m1 -m2 m2-m1 /
 !      write(0,*)"wign3b_2",size(wign3a),l1,l2    
!        call Wigner3j(wign3b, lmin3b, lmax3b, l1, l2, m3b, m1, -m2)
       call Wigner3j_light(vecb(lmin3b+1:), dum, lmax3b, l1, l2, m3b,&
            m1, -m2)
    end if

    st=lmin3b+mod(l1+l2+lmin3b,2)
    do,l=st,lmax3b,2 ! loop with step 2 ( only over valid parity)
!        vecb(l+1)=trpnt%degarr(ind)%scale(l)*trpnt%scale3b(m1,m2)*wign3b(l-lmin3b+1)

!       vecb(l+1)=degsc(l)*trpnt_mscaleb(m1,m2)*vecb(l+1)
        vecb(l+1)=degvec(l)*ordscb*vecb(l+1)
!          vecb(l+1)=trpnt%degarr(ind)%scale(l)*trpnt%scale3b(m2,m1)*vecb(l+1)
    end do
    if(.not. present(trpnt))deallocate(degvec) ! deallocate vector if it is not precomputed
     
  end subroutine trpnt_vec_lm

  subroutine trpnt_degscale(l1,l2,degsc)
    implicit none
    integer,intent(in)::l1,l2
    double precision,intent(inout)::degsc(0:)
    integer::lmin0,lmax0,l
    double precision::wign0(l1+l2+1)
 
    call Wigner3j_light(wign0, lmin0, lmax0, l1, l2, 0, 0,0)

    degsc=0.d0
    do,l=lmin0,l1+l2,2 ! with step 2 
       degsc(l)=sqrt(2*l+1.d0)*sqrt((2*l1+1.d0)*(2*l2+1.d0))*wign0(l-lmin0+1)
    end do

  end subroutine trpnt_degscale

  !subroutine to precompute degree,order and q dependent scales
  subroutine trpnt_init(trpnt,lmax,memout)
    use shtools
    implicit none
    type(trpnt_memstruct),intent(out)::trpnt
    integer,intent(in)::lmax
    integer,optional,intent(out)::memout

    !local parameters
    double precision::wign0(2*lmax+1)
    integer::l1,l2,lmin0,lmax0,ind,l
    integer*8::mem
    double precision::sq2,scale
    integer::m1,m2,m3
    sq2=sqrt(2.d0)


    mem=0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !degree dependent part
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ind=trp_degind(lmax,lmax)
    allocate(trpnt%degarr(ind)) ! allocate enough pointers to hold a symmetric (packed) array of vectors
    
    do,l2=0,lmax
       do,l1=0,l2 ! only upper part
          ind=trp_degind(l1,l2) ! index in trpnt%degarr
          allocate(trpnt%degarr(ind)%scale(0:l1+l2)) ! allocate vector at entry ind
          mem=mem+l1+l2+1
          trpnt%degarr(ind)%scale=0.d0 !initialize to zero
          call trpnt_degscale(l1,l2,trpnt%degarr(ind)%scale)
       end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !order dependent part
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !m3a part
    allocate(trpnt%scale3a(0:lmax,0:lmax))
    mem=mem+(lmax+1)**2
    do,m2=0,lmax
       do,m1=0,lmax ! symmetric
         trpnt%scale3a(m1,m2)=trpnt_mscalea(m1,m2)
       end do
    end do

    !m3b part
    allocate(trpnt%scale3b(0:lmax,0:lmax))
    mem=mem+(lmax+1)**2
    do,m1=0,lmax
       do,m2=0,lmax
         trpnt%scale3b(m2,m1)=trpnt_mscaleb(m1,m2)
       end do
    end do


    if(present(memout))memout=mem/(128) !return working array memory in Kb
  end subroutine trpnt_init

  function trpnt_mscalea(m1,m2)
    implicit none
    double precision::trpnt_mscalea
    integer,intent(in)::m1,m2
    integer::m3
    double precision::sq2
    sq2=sqrt(2.d0)
    m3=trp_m3a(m1,m2)
    trpnt_mscalea=(1-2*mod(m3,2)) ! (-1)^m3
    
    if(m1 .ne. 0)trpnt_mscalea=trpnt_mscalea*sq2/2.d0
    if(m2 .ne. 0)trpnt_mscalea=trpnt_mscalea*sq2/2.d0
    if(m3 .ne. 0)trpnt_mscalea=trpnt_mscalea*sq2
  end function trpnt_mscalea

  function trpnt_mscaleb(m1,m2)
    implicit none
    double precision::trpnt_mscaleb
    integer,intent(in)::m1,m2
    integer::m3
    double precision::sq2
    sq2=sqrt(2.d0)
    m3=trp_m3b(m1,m2)
    if(m1<=m2)then
          trpnt_mscaleb=(1-2*mod(m2,2)) ! (-1)^m2
       if(m1 .ne. 0)trpnt_mscaleb=trpnt_mscaleb*sq2/2.d0
       if(m2 .ne. 0)trpnt_mscaleb=trpnt_mscaleb*sq2
       if(m3 .ne. 0)trpnt_mscaleb=trpnt_mscaleb*sq2/2.d0
    else
          trpnt_mscaleb=(1-2*mod(m1,2)) ! (-1)^m1
       if(m1 .ne. 0)trpnt_mscaleb=trpnt_mscaleb*sq2
       if(m2 .ne. 0)trpnt_mscaleb=trpnt_mscaleb*sq2/2.d0
       if(m3 .ne. 0)trpnt_mscaleb=trpnt_mscaleb*sq2/2.d0
    end if
  end function trpnt_mscaleb

subroutine Wigner3j_light(w3j, jmin, jmax, j2, j3, m1, m2, m3)
!	This subroutine will calculate the Wigner 3j symbols
!
!		j  j2 j3
!		m1 m2 m3
!
!	for all allowable values of j. The returned values in the array j are 
!	calculated only for the limits
!
!		jmin = max(|j2-j3|, |m1|)
!		jmax = j2 + j3
!
!	To be non-zero, m1 + m2 + m3 = 0. In addition, it is assumed that all j and m are 
!	integers.
!
!	This routine is based upon the stable non-linear recurence relations of Luscombe and 
!	Luban (1998) for the "non classical" regions near jmin and jmax. For the classical 
!	region, the standard three term recursion relationship is used (Schulten and Gordon 1975). 
!
!       This is a rewrite of the routine Wigner3j from Mark Wieczorek
!
!       Differences:
!       * We do downard and upward recursion on the fly, eliminating the need for (implicit) allocation of vectors
!       * fixed bug with index 0 adressing
!       * Use forward and backward recursion (from Luscombe & Luban 1998)
!         untill values encounter the classical region, not necessarily the midpoint.
!       * Calculation of ratios x_over_y etc instead of x and y seperately
!       * The sign doesn't need to be fixed afterwards
!       * The common (unknown) scale will be not associated with the midpoint but with 
!         The end point at j=jmax. The values from jmin:jmeet-1 will need to be scaled by the ratio:
!         w3j(jmax)/w3j(jmin) which is determined by overlapping the upward and downward recursion
!         by 1 point and scaling the bottom and top to match the midpoint
!
!
!	Calling parameters
!		IN	
!			j2, j3, m1, m2, m3 	Integer values.
!		OUT	
!			w3j			Array of length jmax - jmin + 1.
!			jmin, jmax		Minimum and maximum values
!						out output array.
!	Dependencies: None
!	
!	Originally Written by Mark Wieczorek August (2004)
!
!	2009-2010: rewrite by Roelof Rietbroek
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer, intent(in) ::	j2, j3, m1, m2, m3
  integer, intent(out) ::	jmin, jmax
  real*8, intent(out) ::	w3j(:)
  
  real*8::ratio,ratio_old,dum
  real*8::tmp1,tmp2,maxtop,tol,meet
  integer :: j,jnum,stderr,jmeet,jup,jlo
  logical::classic

  
  stderr=0
  w3j = 0.0d0
  jmin = max(abs(j2-j3), abs(m2+m3))
  jmax = j2 + j3
  jnum = jmax - jmin + 1

  if (size(w3j) < jnum) then
     write(stderr,*) "Error --- Wigner3j"
     write(stderr,*) "W3J must be dimensioned (J2+J3+1-JMIN) where J2, J3  and JMIN are ", j2, j3,jmin
     write(stderr,*) "Input array is dimensioned ", size(w3j)
     stop
  endif
  
  !check for quick (zero) returns
  if (abs(m2) > j2 .or. abs(m3) > j3) then
     return
  elseif (m1 + m2 + m3 /= 0) then
     return
  elseif (jmax < jmin) then
     return
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Only one term is present (quicky)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (jnum == 1) then
     w3j(1) = (1-2*mod(abs(j2-j3-m1),2))/ sqrt(2.0d0*jmin+1.0d0) ! (-1)^(j2-j3+m2+m3)*sqrt(2*jmax+1) 
     return	
  endif



!starting degrees for iteration
jup=jmax+1
jlo=jmin-1

!set up parameters for 2 term downward recursion in non-classical region
!initial top value (make sure it is the right sign)
w3j(jindex(jmax))=(1.d0-2*mod(abs(j2-j3-m1),2))*1.d-3 !values at the boundaries are generally small
!the scale 1d-3 is just an initial order of magnitude

!initial bottom value
w3j(jindex(jmin))=w3j(jindex(jmax))

!initial value for the maximum encountered
maxtop=abs(w3j(jindex(jmax))) ! maxtop keeps track of the maximum 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! two term downward recursion ( non-classical region, in the direction of increasing w3j)
!!!     
!!! |-------Non-classic----------|---- Classic --- | --<<<--Non-Classic--<<<<--|  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!assume we're not in the classic region
classic=.false.

ratio_old=1.d0 !may be anything but 0 as initial value will be multiplied by zero in the first step


!iterate from top downward untill classic region is reached
do while(.not. classic .and. jup - 1 > jlo+1) !don't calculate the values at jmin
    jup=jup-1

   !get ratios
   call x_y_over_z(jup,tmp1,tmp2)
   !calculate new ratio from old ratio ! note x_over_z(jmax)=0, so initial value will be set ok
   ratio=-tmp1/ratio_old-tmp2
   

   !next value of wigner 3j symbol

   w3j(jindex(jup-1))=ratio*w3j(jindex(jup))
   if(abs(w3j(jindex(jup-1))/w3j(jindex(jup)))<1.d0)classic=.true.  ! also holds in the case ratio=0
   ratio_old=ratio
   !update maximum
   maxtop=max(maxtop,abs(w3j(jindex(jup-1))))
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! two term upward recursion ( in non-classical region in the direction of increasing w3j)
!!! 
!!! |-->>>--Non-classic--->>>---|---- Classic --- | +++++++++Non-Classic+++++++ |  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!assume we're not in the classic region
classic=.false.

!starting degree for iteration
ratio_old=1.d0 !may be anything from 0 as initial value ( will be multiplied by zero in the first step)

do while(.not. classic .and. jlo + 1 < jup - 1) !calculate up to the maximum of jup-2
   jlo=jlo+1
   call y_z_over_x(jlo,tmp1,tmp2)
   ratio=-tmp1-tmp2/ratio_old
   !calculate new ratio from old ratio
!   ratio=-y_over_x(jlo)-z_over_x(jlo)/ratio_old
   !next value of wigner 3j symbol

   w3j(jindex(jlo+1))=ratio*w3j(jindex(jlo))


!write(0,*)'up2',jlo+1,tmp1,tmp2,y_over_x(jlo),z_over_x(jlo)

   if(abs(w3j(jindex(jlo+1))/w3j(jindex(jlo)))<1.d0)classic=.true. ! also holds in the case ratio=0
   ratio_old=ratio


end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! choose an initial meeting point ( note the value of w3j should be non-zero at that point)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! start with the (approximate) center of the classical region
jmeet=(jlo+jup)/2+1 !always add a 1 to avoid jmeet==jmin, jmeet=jmax is allowed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Three term recursion downward in classical region (untill meeting point)
!!
!!! |++++++++Non-classic++++++|------ meetpoint -<<<<<<- | +++++++++Non-Classic+++++++ |  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



do,j=jup-1,jmeet+1,-1 ! Value at jmeet is also calculated
   call x_y_over_z(j,tmp1,tmp2)
   w3j(jindex(j-1))=-tmp2*w3j(jindex(j))-tmp1*w3j(jindex(j+1))
   !update maximum
   maxtop=max(maxtop,abs(w3j(jindex(j-1))))
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Three term recursion upward in classical region  (untill meeting point)
!!
!!! |++++++++Non-classic++++++|->>>>>- meetpoint +++++++++ | +++++++++Non-Classic+++++++ |  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do,j=jlo+1,jmeet-2 ! do not calculate value at jmeet
   call y_z_over_x(j,tmp1,tmp2)
   w3j(jindex(j+1))=-tmp2*w3j(jindex(j-1))-tmp1*w3j(jindex(j))
end do


!scale upper and lower parts to a common midpoint 
! we use the upward three term recursion to search for a non-zero (large enough) meetingpoint
tol=0.05 ! set tolerance ( 0< tol<1) point will be accepted as a meeting point if it is larger then 5% of the maximum value encountered ( in the downward recursion) 
dum=w3j(jindex(jmeet-1)) !copy initial value
do,j=jmeet-1,jmax-1 ! loop to search for a non-zero wigner value to apply scale ( should be a quick end)
   call y_z_over_x(j,tmp1,tmp2)
   !apply three term upward recursion
   if(j == jmin)then !w3j(jindex(jmin-1)) is zero ( only one term)
      dum=-tmp1*w3j(jindex(j))
   else
      dum=-tmp1*w3j(jindex(j))&
           -tmp2*w3j(jindex(j-1))
   end if
   !check for validity of the meeting point
   if(abs(w3j(jindex(j+1))) > tol*maxtop)exit ! exit loop when a large enough value is found
   w3j(jindex(j+1))=dum ! overwrite value with that of the upward recursion (meeting point will be shifted)
end do

jmeet=j+1 ! update meeting point
!!! The meeting point will be now be shifted upwards untill a large enough common value has been found
!!! |                Factor1                |->                  Factor2                  |  
!!! |++++++++Non-classic++++++| ++++++++ meetpoint ++++++++ | +++++++++Non-Classic+++++++ |  
!!! At the meetinpoint we have one value from the downward recursion and a overlap value(dum) from the upward recursion
!! now we need to scale the upward and downward recursion such that the value at the meetingpoint will be the same


!determine a reasonable common value ( maximum top will be close to 1) at the meeting point and make sure sign is correct
if(w3j(jindex(jmeet))*w3j(jindex(jmax)) >= 0.d0)then !meeting point has the same sign as maximum
   
   meet=(1.d0-2*mod(abs(j2-j3-m1),2))*tol
else
   meet=-(1.d0-2*mod(abs(j2-j3-m1),2))*tol
end if
   




ratio=meet/dum ! scale factor for the values from upward recursion


!scale lower part and apply correct sign
do,j=1,jindex(jmeet-1)
   w3j(j)=w3j(j)*ratio
end do

!write(0,*)"meet,dum",jmeet,meet,dum

ratio=meet/w3j(jindex(jmeet)) !scale factor for the values from downward recursion

!scale upper part
do,j=jindex(jmeet),jindex(jmax)
   w3j(j)=w3j(j)*ratio
end do
!write(0,*)"meet,jmeet",jmeet,meet,w3j(jindex(jmeet))

!write(0,*)"jmeet",jmeet,jlo,jup,jmin,jmax
!normalize w3j
call normw3j()

!finished

contains 
  integer function jindex(j)
    integer :: j
    jindex = j-jmin+1
  end function jindex

  subroutine normw3j
    real*8:: norm
    integer j
    
    norm = 0.0d0
    do j = jmin, jmax
       norm = norm + dble(2*j+1) * w3j(jindex(j))**2
    enddo
     w3j(1:jnum) = w3j(1:jnum) / sqrt(norm)

  end subroutine normw3j


! !function y over z
! !                               /                                                                      \
! !  Y(j)                 (2j+1) |  ( m2 + m3 ){ j2( j2 + 1 ) - j3( j3 + 1 ) } - ( m2 - m3 ) j ( j + 1 )  |      
! !                               \                                                                      /
! ! ------ =     ---------------------------------------------------------------------------------------------------------
! !                                 /        (j2 - j3)^2        (j2 + j3 +1)^2                (m2 - m3)^2      \
! !  Z(j)           (j+1)  j^3 sqrt|  ( 1 -  ------------   ) ( --------------  - 1 ) ( 1 -  -------------   )  |
! !                                 \           j^2                  j^2                         j^2           /
  

! ! valid only for  jmin < j <=jmax
! !

! !separating j(j+1)

! !                               /              j2( j2 + 1 )   j3( j3 + 1 )                \
! !  Y(j)                 (2j+1) |  ( m2 + m3 ){ ------------ - ------------ } - ( m2 - m3 ) |      
! !                               \              j ( j + 1 )    j ( j + 1 )                 /
! ! ------ =     -------------------------------------------------------------------------------------------------
! !                             /        (j2 - j3)^2        (j2 + j3 +1)^2                (m2 - m3)^2      \
! !  Z(j)              j^2 sqrt|  ( 1 -  ------------   ) ( --------------  - 1 ) ( 1 -  -------------   )  |
! !                             \           j^2                  j^2                         j^2           /
! !




!   function y_over_z(j)
!     real*8::y_over_z,top,bot,ratio
!     integer::j
! !     if(j<=jmin .or. j > jmax)then
! !        write(0,*)"error,out of bound, y_over_z"
! !        stop
! !     end if

!     top=((m2+m3)*(dble(j2*(j2+1)-j3*(j3+1))/dble(j*(j+1))) - ( m2 - m3 ) )

!     bot=sqrt((1.d0-(dble((j2-j3))/dble(j))**2)*&
!          ((dble((j2+j3+1))/dble(j))**2 - 1.d0)*&
!          (1.d0-(dble((m2+m3))/dble(j))**2))

!     y_over_z=(dble(2*j+1)/dble(j**2))*top/bot

!   end function y_over_z


! !                               /              j2( j2 + 1 )   j3( j3 + 1 )                \
! !  Y(j)     j ( j + 1 ) (2j+1) |  ( m2 + m3 ){ ------------ - ------------ } - ( m2 - m3 ) |      
! !                               \              j ( j + 1 )    j ( j + 1 )                 /
! ! ------ =     -------------------------------------------------------------------------------------------------
! !                               /          (j2 - j3)^2        (j2 + j3 +1 )^2              (m2 - m3)^2      \
! !  X(j)     j ( j + 1 )^3  sqrt|  ( 1 -  ------------   ) (  --------------  - 1 ) ( 1 -  -------------   )  |
! !                               \           (j+1)^2             (j+1)^2                     (j+1)^2         /
! !

! ! for jmin <= j < jmax

! !in the case of jmin = j = 0 we will have a special case because then also m2 + m3 = 0 = j causing ( m2 + m3 ) / j = 1 
! ! this case will also work below but a two term recursion may not be used as it will cause a division by zero
!   function y_over_x(j)
!     real*8::y_over_x,top,bot,ratio
!     integer::j
! !     if(j >= jmax .or. j < jmin)then
! !        write(0,*)"error,y_over_x out of bound"
! !        stop
! !     end if
!     if(j .eq. 0)then
!        ratio=1.d0
!     else
!        ratio=dble(m2+m3)/dble(j)
!     end if

!     top=ratio*dble(j2*(j2+1)-j3*(j3+1))/dble(j+1) - ( m2 - m3 )

!     bot=sqrt((1.d0-(dble((j2-j3))/dble(j+1))**2)*&
!          ((dble((j2+j3+1))/dble(j+1))**2 - 1.d0)*&
!          (1-(dble((m2+m3))/dble(j+1))**2))
    
!     y_over_x=(dble(2*j+1)/dble((j+1)**2))*top/bot

!   end function y_over_x

!                                /          (j2 - j3)^2        (j2 + j3 +1 )^2              (m2 - m3)^2      \
!  X(j)     j  ( j + 1 ) ^3 sqrt|  ( 1 -  ------------   ) (  --------------  - 1 ) ( 1 -  -------------   )  |
!                                \           (j + 1)^2             (j + 1)^2                  (j + 1)^2      / 
! ------ =     -------------------------------------------------------------------------------------------------
!                               /          (j2 - j3)^2        (j2 + j3 +1 )^2              (m2 - m3)^2      \
!  Z(j)     ( j + 1 ) j^3  sqrt|  ( 1 -  ------------   ) (  --------------  - 1 ) ( 1 -  -------------   )  |
!                               \           (j)^2                (j)^2                        (j)^2         / 
!



!                         
!  X(j)             j            / ( j + 1 )^2 - ( j2 - j3 )^2    
!                                       
! ------ =     ----------- sqrt |  ----------------------------  *  
!                         
!  Z(j)         ( j + 1 )        \   ( j )^2 - ( j2 - j3 )^2      
!                                
!

!                    ( j2 + j3 + 1 )^2 - ( j + 1 )^2      ( j + 1 )^2 - ( m2 + m3 )^2   \ 
                                                                      
!                    --------------------------------- * -----------------------------  |
!                     
!                     ( j2 + j3 + 1 )^2  -  ( j )^2          ( j )^2 - ( m2 + m3 )^2      / 


! for jmin < j <= jmax

!   function x_over_z(j)
!     real*8::x_over_z
!     integer::j
    
! !     if(j<=jmin .or. j > jmax)then
! !        write(0,*)"error,x_over_z out of bound"
! !        stop
! !     end if
!     x_over_z= sqrt( (dble((j+1)**2-(j2-j3)**2)/dble((j)**2-(j2-j3)**2))*&
!          (dble((j2+j3+1)**2-(j+1)**2)/dble((j2+j3+1)**2-(j)**2))*&
!          (dble((j+1)**2-(m2+m3)**2)/dble((j)**2-(m2+m3)**2)))
    
!     x_over_z=x_over_z*(dble((j))/dble(j+1))
    
!   end function x_over_z

!   !reciprocal of z_over_x valid for jmin <= j < jmax
!   function z_over_x(j)
!     real*8::z_over_x
!     integer::j
    
! !     if(j<jmin .or. j >=jmax)then
! !        write(0,*)"error,z_over_x out of bound"
! !        stop
! !     end if

!     if( j.eq. 0)then ! special case 
!        z_over_x=0.d0
!        return
!     end if
       
!     z_over_x= sqrt( (dble((j)**2-(j2-j3)**2)/dble((j+1)**2-(j2-j3)**2))*&
!          (dble((j2+j3+1)**2-(j)**2)/dble((j2+j3+1)**2-(j+1)**2))*&
!          (dble((j)**2-(m2+m3)**2)/dble((j+1)**2-(m2+m3)**2)))
    
!     z_over_x=z_over_x*(dble((j+1))/dble(j))
    
!   end function z_over_x

  subroutine x_y_over_z(j,x_o_z,y_o_z)
    real*8,intent(out)::x_o_z,y_o_z
    integer,intent(in)::j
!   real*16::tmp1,tmp2 ! use for improved accuracy
   real*8::tmp1,tmp2


    !common denominator part
    tmp1=(j**2-(j2-j3)**2)
    tmp1=tmp1*((j2+j3+1)**2-j**2)
    tmp1=tmp1*(j**2-(m1)**2)

    !xratio
    tmp2=((j+1)**2-(j2-j3)**2)
    tmp2=tmp2*((j2+j3+1)**2-(j+1)**2)
    tmp2=tmp2*((j+1)**2-(m1)**2)

    x_o_z=(dble(j)/(j+1))*sqrt(tmp2/tmp1)

    !yratio
    
    tmp2=((-m1)*(j2*(j2+1)-j3*(j3+1))-&
         (m2-m3)*(j*(j+1)))/sqrt(tmp1)

    y_o_z=(1.d0+dble(j)/(j+1))*tmp2

  end subroutine x_y_over_z

  subroutine y_z_over_x(j,y_o_x,z_o_x)
    real*8,intent(out)::y_o_x,z_o_x
    integer,intent(in)::j
!    real*16::tmp1,tmp2 ! use for improved accuracy
    real*8::tmp1,tmp2
    if(j .eq. 0)then ! quick return
       z_o_x=0.d0
       y_o_x=-(m2-m3)*(j+1)/sqrt((j2+j3+1)**2-1.d0)
       return
    end if

    !common denominator part
    tmp1=((j+1)**2 - (j2-j3)**2)
    tmp1=tmp1*((j2+j3+1)**2 - (j+1)**2)
    tmp1=tmp1*((j+1)**2 - m1**2)

    ! zratio
    tmp2=((j)**2-(j2-j3)**2)
    tmp2=tmp2*((j2+j3+1)**2-(j)**2)
    tmp2=tmp2*((j)**2-(m1)**2)

    z_o_x=(dble(j+1)/(j))*sqrt(tmp2/tmp1)
    
    !yratio
    tmp2=((-m1)*(j2*(j2+1)-j3*(j3+1))-&
         (m2-m3)*(j*(j+1)))/sqrt(tmp1)

    y_o_x=(2.d0+1.d0/j)*tmp2

  end subroutine y_z_over_x

end subroutine Wigner3j_light



end module SH_trpnt

!!File with two methods to retrieve a Spherical harmonic vector index from inputted degree and order and vice versa
!!The SH vector must be defined accordingly
!!Where the index of a certain component with degree l and order m is defined as
!!index	= (l*(l+1))/2+m+1
!!
!! coded by Roelof Rietbroek 10-5-2007

!!Function SH_pos(l,m) outputs SH vector index as a function of degree l and order m
!!integer:: l,m,SH_pos
function SH_pos(l,m)
implicit none
integer,intent(in):: l,m
integer::SH_pos
!return negative index if order is larger than degree
if (m > l) then 
SH_pos=-1
return
end if

SH_pos=(l*(l+1))/2+m+1

end function SH_pos

!!function which gives the indices corresponding to a order degree sorting

function SH_opos(l,m,lmax)
implicit none
integer, intent(in)::l,m,lmax
integer::SH_opos
if (m > l) then 
SH_opos=-1
return
end if


SH_opos=m*(lmax+1)-(m*(m-1))/2+1+l-m


end function SH_opos

!!function which gives position parameter as a function of degree and order
!!this index is used in programs from Juergen Kusche
function SH_tpos(l,m,q,lmax,lmin)
implicit none
integer,intent(in)::l,m,lmax,lmin,q
integer::SH_tpos
integer,parameter::nend=(360+1)*(360+2)/2
integer::or,de,SH_pos,pos(nend,0:1),qm,cs,lma,lmi,ind

save lmi,lma,pos

if (lma .ne. lmax .or. lmi .ne. lmin)then !precompute an index using the degree sorted scheme index
!write(*,*)'blah die blah'
lma=lmax
lmi=lmin
ind=0
pos=0
      do or=0,lma 
         qm=0
         if (or.gt.0) qm=1 
         do cs=0,qm
            do de=max(lmi,or),lma
            ind=ind+1
            pos(SH_pos(de,or),cs)=ind
            end do
         end do
      end do

end if 

!now retrieve the right index
SH_tpos=pos(SH_pos(l,m),q)
end function SH_tpos

!inverse of the above returns the degree, order and sign of a coefficient for a matrix in tight form
subroutine SH_tlm(pos,l,m,q,lmax,lmin)
implicit none
integer,intent(in)::lmax,lmin,pos
integer,intent(out)::l,m,q
integer,parameter::nend=(360*(360+1))/2
integer::degords(nend,3),cs,de,or,qm,lmi,lma,ind
save degords,lma,lmi !only calculate the index vector if neccessary

if (lma .ne. lmax .or. lmi .ne. lmin)then !precompute an index using the degree sorted scheme index
!write(*,*)'blah die blah'
lma=lmax
lmi=lmin
ind=0
degords=0
      do or=0,lma 
         qm=0
         if (or.gt.0) qm=1 
         do cs=0,qm
            do de=max(lmi,or),lma
            ind=ind+1
            degords(ind,1)=de!store degree
            degords(ind,2)=or!store order
            degords(ind,3)=cs!store sign (cosine 0 sine 1)
            end do
         end do
      end do

end if 

!return degree order and sign

l=degords(pos,1)
m=degords(pos,2)
q=degords(pos,3)


end subroutine SH_tlm






!!subroutine SH_lm calculates the degree and order from a given SH vector index
!!integer::pos,l,m
subroutine SH_lm(pos,l,m)
implicit none
integer, intent(in)::pos
integer, intent(out)::l,m
double precision :: dum

dum=(sqrt(1.d0+8*(pos-1))-1.d0)/2.d0
l=int(dum)
m=int((dum-l)*(l+1))
end subroutine SH_lm


!subroutine which calculates the quotient of two factorials:
!            top!
! qfact   = -----
!            bot!
function qfact(top,bot)

implicit none
integer,intent(in)::top,bot

integer::mini,maxi,i
double precision::qfact

mini=min(top,bot)+1
maxi=max(top,bot)
qfact=1
if(top > bot)then

   do,i=mini,maxi
      qfact=qfact*i
   end do
   return
else

   do,i=mini,maxi
      qfact=qfact/dble(i)
   end do
   return

end if

end function qfact

!test program for the routines
!uncomment the part below to produce an executable
!program testind
!implicit none
!integer ::l,m,SH_pos,pos,lt,mt

!do, m=0,30
!   do, l=m,30
!      pos=SH_pos(l,m)
!      call SH_lm(pos,lt,mt)
!      write(*,*)l,m,lt,mt,pos
!   end do
!
!end do
!
!end program testind

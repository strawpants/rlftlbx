!!Coded by Roelof Rietbroek, Mon Jun 24 21:17:31 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de
!! this subroutine transforms degree 1 Love numbers in another (isomorphic reference frame)
!! See the paper by Blewitt 2003, and the Phd dissertation of R. Rietbroek

subroutine SH_Love_d1_trans(h1in,l1in,k1in,h1out,l1out,k1out,FR)
implicit none
double precision,intent(in)::h1in,l1in,k1in
double precision,intent(out)::h1out,l1out,k1out
character(*),intent(in)::FR

!local variables
double precision::isomorph !isomorphic frame parameter
integer::stderr
stderr=0

select case(FR)
case('CM')
   isomorph=1+k1in
case('CF')
   isomorph=(h1in+2.d0*l1in)/3.d0
case('CE')
   isomorph=k1in
case('CL')
   isomorph=l1in
case('CH')
   isomorph=h1in
case default
   write(stderr,*)"ERROR: unknown reference origin requested:", FR
   stop
end select

h1out=h1in-isomorph
l1out=l1in-isomorph
k1out=k1in-isomorph


end subroutine SH_Love_d1_trans

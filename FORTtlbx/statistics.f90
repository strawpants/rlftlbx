!!Coded by Roelof Rietbroek, Sun Sep 19 14:38:05 2010
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!set of statistical functions (taken from numerical recipes in fortran 90 and adapted in double precision)

!returns the student-t probability of P( T > t , T < -t ) for df degrees of freedom
function pval_stud_t(t,df)
  implicit none
  double precision::pval_stud_t
  double precision,intent(in)::t
  integer,intent(in)::df
  double precision::betai_s
  pval_stud_t=betai_s(0.5d0*df,0.5d0,df/(df+t**2))
end function pval_stud_t

FUNCTION betai_s(a,b,x)
IMPLICIT NONE
double precision, INTENT(IN) :: a,b,x
double precision :: betai_s
!Returns the incomplete beta function Ix (a, b).

double precision :: betacf_s,gammln_s 
double precision :: bt
integer,parameter::stderr=0
!write(stderr,*)"betai",a,b,x

if(x < 0.0 .or. x > 1.0)then
   write(stderr,*)"betai_s: x doesn't satisfy 0 <= x <= 1",x
   stop
end if

if (x == 0.0 .or. x == 1.0) then
   bt=0.0
else
   !Factors in front of the continued fraction.
   bt=exp(gammln_s(a+b)-gammln_s(a)-gammln_s(b)&
        +a*log(x)+b*log(1d0-x))
end if
if (x < (a+1.0)/(a+b+2d0)) then
   !Use continued fraction directly.
   betai_s=bt*betacf_s(a,b,x)/a
else
   !Use continued fraction after making the symmetry transformation.
   betai_s=1d0-bt*betacf_s(b,a,1d0-x)/b

end if
END FUNCTION betai_s


FUNCTION betacf_s(a,b,x)
IMPLICIT NONE
double precision, INTENT(IN) :: a,b,x
double precision :: betacf_s
INTEGER, PARAMETER :: MAXIT=100
double precision, PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
!Used by betai: Evaluates continued fraction for incomplete beta function by modified
!Lentz’s method (§5.2).
double precision :: aa,c,d,del,h,qab,qam,qap
INTEGER :: m,m2
integer,parameter::stderr=0

! write(stderr,*)"betacf"

qab=a+b
!These q’s will be used in factors that occur in the coefficients (6.4.6).
qap=a+1d0
qam=a-1d0
c=1d0
!First step of Lentz’s method.
d=1d0-qab*x/qap
if (abs(d) < FPMIN) d=FPMIN
d=1.0/d
h=d
do m=1,MAXIT
   m2=2*m
   aa=m*(b-m)*x/((qam+m2)*(a+m2))
   d=1d0+aa*d
!One step (the even one) of the recurrence.
if (abs(d) < FPMIN) d=FPMIN
c=1d0+aa/c
if (abs(c) < FPMIN) c=FPMIN
d=1d0/d
h=h*d*c
aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
d=1d0+aa*d
!Next step of the recurrence (the odd one).
if (abs(d) < FPMIN) d=FPMIN
c=1d0+aa/c
if (abs(c) < FPMIN) c=FPMIN
d=1.0/d
del=d*c
h=h*del
if (abs(del-1d0) <= EPS) exit !Are we done?
end do
if (m > MAXIT)then
   write(stderr,*)'a or b too big, or MAXIT too small in betacf_s'
   stop
end if

betacf_s=h
END FUNCTION betacf_s


FUNCTION gammln_s(xx)
IMPLICIT NONE
interface
   FUNCTION arth_d(first,increment,n)
     double precision, INTENT(IN) :: first,increment
     INTEGER, INTENT(IN) :: n
     double precision, DIMENSION(n) :: arth_d
   end FUNCTION arth_d
end interface
double precision, INTENT(IN) :: xx
double precision :: gammln_s
!Returns the value ln[Γ(xx)] for xx > 0.
double precision :: tmp,x
!Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
!accuracy is good enough.
double precision :: stp = 2.5066282746310005d0
double precision, DIMENSION(6) :: coef = (/76.18009172947146d0,&
-86.50532032941677d0,24.01409824083091d0,&
-1.231739572450155d0,0.1208650973866179d-2,&
-0.5395239384953d-5/)
integer,parameter::stderr=0

!write(stderr,*)"gammaln",xx

if(xx <= 0.0)then
   write(stderr,*)"gammln_s xx must be > 0"
   stop
end if

x=xx
tmp=x+5.5
tmp=(x+0.5d0)*log(tmp)-tmp
gammln_s=tmp+log(stp*(1.000000000190015d0+&
sum(coef(:)/arth_d(x+1.0d0,1.0d0,size(coef))))/x)
END FUNCTION gammln_s

FUNCTION arth_d(first,increment,n)
implicit none
double precision, INTENT(IN) :: first,increment
INTEGER, INTENT(IN) :: n
double precision, DIMENSION(n) :: arth_d
INTEGER :: k,k2
double precision :: temp
INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
integer,parameter::stderr=0
! write(stderr,*)"arth_d"

if (n > 0) arth_d(1)=first
if (n <= NPAR_ARTH) then
   do k=2,n
      arth_d(k)=arth_d(k-1)+increment
   end do
else
   do k=2,NPAR2_ARTH
      arth_d(k)=arth_d(k-1)+increment
   end do
   temp=increment*NPAR2_ARTH
   k=NPAR2_ARTH
   do while(k<n)
!   if (k >= n) exit
   k2=k+k
   arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
   temp=temp+temp
   k=k2
end do
end if
END FUNCTION arth_d

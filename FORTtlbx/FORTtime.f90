!! fortran module to work with datetime functionality
!! wrapper around the C ctime.h functionality
!! See also ctimeF.c


module FORTtime
implicit none

type time_t
    integer(8):: TData!holds the equivalent of the time_t in C
end type time_t

contains

function FTime_fromDecyr(decyr)
implicit none
double precision,intent(in)::decyr      
type(time_t)::FTime_fromDecyr,fmktime
integer::mndy(12)
integer::i,dpy,yr,doy,sec
data mndy/31,28,31,30,31,30,31,31,30,31,30,31/


yr=int(decyr)
if(mod(yr,4) .eq. 0)then
    mndy(2)=29
    dpy=366
else
    dpy=365
end if
!day of year
doy=int((decyr-yr)*dpy)
!seconds
sec=((decyr-yr)*dpy-doy)*86400

FTime_fromDecyr=fmktime(yr,0,doy+1,0,0,sec)

end function



end module    


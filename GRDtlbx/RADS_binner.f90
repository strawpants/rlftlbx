!!Coded by Roelof Rietbroek, Wed Feb  5 15:05:42 2014
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de
!!program RADS_binner bins altimetry data along-track in space and time

program RADS_binner
implicit none
integer::itharg,iargc,narg
integer::i
character(300)::dum
character
logical::listsat

!defaults
listsat=.false.

!!process command line options
narg=iargc()

if(narg <1)call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('a')
         satcode=dum(3:)
      case('l')
         listsat=.true.
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
   end if
end do



end program RADS_binner

subroutine help()
implicit none
character(4)::frmt
integer::stderr
stderr=0
frmt='(A)'

write(stderr,frmt)' Program RADS_binner bins rads altimetry sea level anomalies in space and time '
write(stderr,frmt)" Usage RADS_binner [OPTIONS] > OUTPUTBINFILE"
write(stderr,frmt)" where OPTIONS can be:"
write(stderr,frmt)" -aSATCODE: Specify the satellite of interest"
write(stderr,frmt)" -l: List available satellites"
stop
end subroutine help



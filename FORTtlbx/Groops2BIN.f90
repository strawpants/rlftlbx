!!Coded by Roelof Rietbroek, Tue Jun  5 11:21:28 2012
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!!program which converts (some) groops binary data in BIN format


program Groops2BIN
use GroopsFile
use forttlbx
implicit none
integer:: stderr,i,narg,iargc,itharg
character(200)::dum
type(Groopsdat)::dat,rhs
logical::normal


!! defaults
stderr=0
normal=.false.


!!process command line options
narg=iargc()

if (narg <1)call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('m')! specify matrix type
         select case(dum(2:2))
            case('1')!general matrix
               dat%mattype='GENERAL'
            case('2')
               dat%mattype='UPPER'
            case('3')
               normal=.true.
               dat%mattype='UPPER'
               rhs%mattype='GENERAL'
         end select
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option but a file
      dat%filename=trim(dum)
   end if
end do

dat%mattype='UPPER'
call read_groops(dat)

end program Groops2BIN


subroutine help()
implicit none
character(4)::frmt
integer::stderr
stderr=0
frmt='(A)'


write(stderr,frmt)'Program Groops2BIN converts groops binary data to BINV format'
write(stderr,frmt)'Usage Groops2BIN [OPTIONS] GROOPSFILE'
write(stderr,frmt)'Where OPTIONS may be:'
stop
end subroutine help

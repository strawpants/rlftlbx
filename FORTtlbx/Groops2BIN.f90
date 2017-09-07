!!Coded by Roelof Rietbroek, Tue Jun  5 11:21:28 2012
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!!program which converts (some) groops binary data in BIN format


program Groops2BIN
use GroopsFile
use forttlbx
use binfiletools
implicit none
integer:: i,narg,iargc,itharg
character(200)::dum,fout
type(BINdat)::dat
logical::normal
integer::ind,lmax,lmin,l,m,q
integer::s1scheme,s2scheme
integer::iinsert,nmeta

!! defaults
!stderr=0
normal=.false.
lmax=-1
lmin=0
s1scheme=1 !default to degreeWise sorting scheme
s2scheme=0
nmeta=0


!!process command line options
fout='stdout'
narg=iargc()

if (narg <1)call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
     select case(dum(2:2))
     case('l')!set maximum and minimum degree for SH vectors
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
         !select case(dum(2:2))
         !case('1')!normal equation system
            !normal=.true.
            !dat%mattype='UPPER'
           !!case('1')!general matrix
               !!dat%mattype='GENERAL'
            !!case('2')
               !!dat%mattype='UPPER'
            !!case('3')
               !!normal=.true.
               !!dat%mattype='UPPER'
               !!rhs%mattype='GENERAL'
         !end select
     case('s')!description of side1 (obligatory argument)
        ind=index(dum,'DegreeWise')
        if(ind .ne. 0)then
            s1scheme=1
        end if  
    case('F')
        fout=trim(dum(4:))
     !case('s2')!description of side2 (not needed when -n is selected)
        !write(stderr,*)"ERROR:not implemented yet"
        !stop(1)
     !case('n')
         !normal=.true.
    ! case('M')!append meta data
        
     case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:3)
         stop
      end select
   else!argument is not an option but a file
      dat%file=trim(dum)
   end if
end do

!if(normal)then
    call read_GROOPSNEQ(dat)    
!else
!    call read_GROOPSBIN(dat)
!end if  

!!setup side description
select case(s1scheme)
case(1)
    if(lmax<0)then
         write(stderr,*)'ERROR: missing -l option: lmax,lmin must be set'
         stop(1)
    end if
    call getDegreeWise_desc(dat%side1_d,lmax,lmin)

    !also adds lmax and lmin to meta data
    iinsert=size(dat%ints)+1
    dat%nint=dat%nint+2
    call realloc_ptr(dat%ints_d,2)
    call realloc_ptr(dat%ints,2)
    dat%ints_d(iinsert)='Lmax'
    dat%ints(iinsert)=lmax
    iinsert=iinsert+1
    dat%ints_d(iinsert)='Lmin'
    dat%ints(iinsert)=lmin

case default
         write(stderr,*)'ERROR: -s1 option is not understood'
         stop(1)
end select



!check whether side description has the expected amount of entries
if(size(dat%side1_d) .ne. dat%nval1)then
         write(stderr,*)'ERROR: provided side description is inconsistent in number'
         stop(1)
end if


!write data to file or standard input
dat%file=fout
call write_BINtype(dat)

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
write(stderr,frmt)"  -l=lmax[,lmin]: Assume data describes Stokes coefficients"
write(stderr,frmt)"        and Sets a maximum and minimum degree"
write(stderr,frmt)" -s1=SORTINGSCHEME Set sorting scheme of the rows"
write(stderr,frmt)"    defaults to -s1=DegreeWise"
write(stderr,frmt)" -F=Fileout: write binary data to Fileout(default stdout)"
write(stderr,frmt)'NOTE: this program can currently only convert NEQS'
stop
end subroutine help

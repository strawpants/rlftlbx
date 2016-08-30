!small program to make a diagonal regularization matrix

!!Coded by Roelof Rietbroek, Fri Apr  4 15:28:35 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de


program SH_diagreg
use SHtlbx
implicit none
integer::l,m,lmax,lmin,ind,itharg,narg,stderr,i
double precision::dumdbl,scale1,scale2,power
character(150)::dum
logical::powlaw,grad,inverse,equi
integer::iargc


scale1=1.d0
scale2=1.d0
powlaw=.false.
inverse=.false.
grad=.false.
equi=.false.


stderr=0
lmax=100
lmin=0
itharg=0

narg=iargc()

if(narg < 1)call help()

do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit! exit loop when last argument has been read in
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then
      
      select case(dum(2:2))
         case('t')
            select case(dum(3:3))
            case('p')
               powlaw=.true.
               read(dum(5:),*)power
            case('k')
               powlaw=.true.
               scale1=10d10
               power=4.d0
            case('s') ! surface gradient
               scale1=4*pi*(GM/RE)**2
               grad=.true.
               lmin=1
            end select
         case('s')
            read(dum(4:),*)scale2

         case('l')
            ind=index(dum,',')
            if(ind .ne. 0) then !also a min degree is specified
               read(dum(4:ind-1),*)lmax
               read(dum(ind+1:),*)lmin
            else
               read(dum(4:),*)lmax
         end if
      case('e')!express in equivalent water height (each entry will be additionally scaled by RE*rho_e/(rho_w*3)*(2l+1)/(1+kl))
         equi=.true.
      case('i')
         inverse=.true.
         
         case default
            write(stderr,*)'Unknown option selected:',dum(1:2)
            call help()
      end select



   else
      write(stderr,*)'Unknown argument: ',trim(dum)
      call help()
   end if
   
   
end do

if(.not. (powlaw .or. grad))then
   write(stderr,*)'No output type specified'
   call help()
end if



!reset scale
scale1=scale1*scale2


if(powlaw)then
   if(inverse)then
      scale1=1.d0/scale1
      power=-power
   end if
   do,l=lmin,lmax
      dumdbl=scale1*(l**power)
      do,m=0,l
      write(*,'(A3,1x,i3,i3,16x,G12.9)')'TCN',l,m,dumdbl
      write(*,'(A3,1x,i3,i3,16x,G12.9)')'TSN',l,m,dumdbl
      end do
   end do

else if(grad)then
   if(inverse)then

      do,l=lmin,lmax
         !m .eq. 1 entries
         dumdbl=scale1*2
         write(*,'(A3,1x,i3,i3,16x,G12.9)')'TCN',l,1,1.d0/dumdbl
         write(*,'(A3,1x,i3,i3,16x,G12.9)')'TSN',l,1,1.d0/dumdbl         
         do,m=2,l
            dumdbl=scale1*(l+1)*l
            write(*,'(A3,1x,i3,i3,16x,G12.9)')'TCN',l,m,1.d0/dumdbl
            write(*,'(A3,1x,i3,i3,16x,G12.9)')'TSN',l,m,1.d0/dumdbl
         end do
      end do
   else
      do,l=lmin,lmax
         !m .eq. 1 entries
         dumdbl=scale1*2
         write(*,'(A3,1x,i3,i3,16x,G12.9)')'TCN',l,1,dumdbl
         write(*,'(A3,1x,i3,i3,16x,G12.9)')'TSN',l,1,dumdbl         
         do,m=2,l
            dumdbl=scale1*l*(l+1)
            write(*,'(A3,1x,i3,i3,16x,G12.9)')'TCN',l,m,dumdbl
            write(*,'(A3,1x,i3,i3,16x,G12.9)')'TSN',l,m,dumdbl
         end do
      end do
   end if

end if






end program SH_diagreg

subroutine help()
character(5)::frmt
integer::unit
frmt='(A)'
unit=6
write(unit,frmt)"Program SH_diagreg"
write(unit,frmt)"Prints the diagonal of various geopotential regularization matrices to standard output" 
write(unit,frmt)"Usage SH_diagreg [OPTIONS]"
write(unit,frmt)" Where options may be:"
write(unit,frmt)" -ts: minimize surface gradient" 
write(unit,frmt)" -tk: Kaula's rule (10e10*l^4)"
write(unit,frmt)" -tp=POWER: general degree dependent power law, l^POWER"
write(unit,frmt)" -l=lmax,lmin"
write(unit,frmt)" -s=SCALE: Apply a constant scale to the matrix"
stop



end subroutine help

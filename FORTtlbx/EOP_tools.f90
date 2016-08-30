!!Coded by Roelof Rietbroek, Fri Feb 22 11:15:31 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!! set of subroutines to work with the EOP parameters
!!Updated by Roelof Rietbroek, Fri Dec 20 15:12:00 2013
!! added option to subtract the mean iers pole or not to obtain the the polar motion

module EOP_tools
implicit none
type EOPdat
   character(200):: file='EOPcurrent' ! default eop file name 
   double precision,pointer, dimension(:,:)::dat=>null()
   double precision,pointer,dimension(:,:)::time=>null()
   integer::nt=0
   character(3)::type='C04'
end type EOPdat

!contains all the goodies:
type(EOPdat)::EOP

contains

subroutine read_EOP()
  character(200)::dir,dum
  character(100)::frmt
  integer::un,last,nmax,stderr,i,sz
  stderr=0
  nmax=90000 !! allows a maximum of 90000 days to be read in (will do untill 2100)
  

  un=13
  last=0
  allocate(EOP%dat(nmax,12))
  allocate(EOP%time(nmax,4))
  EOP%dat=0.d0
  EOP%time=0.d0

  ! get the environment variable WORK_DIR
  call getenv('WORK_DIR',dir)
  dir=trim(dir)//'/data/EOP/'

  !open for reading
  open(unit=un,file=trim(dir)//trim(EOP%file),status='old')
  !check the first the first 13 header lines
  do i=1,13
     read(unit=un,fmt='(A200)')dum
     !test for C01 file
     if(index(dum,': EOP(IERS)') .ne. 0)EOP%type='EOP'  

  end do

  select case(EOP%type)
  case('EOP') ! for example an C01 file
     frmt='(2X,F10.3,1X,F9.6,1X,F9.6,1X,F10.7,3X,F9.6,4X,F9.6,3X,F9.6,1X,F9,6,1X,F10.7,1X,F9.6,2X,F9.6)'
     frmt='(11F)'
     sz=1
     do while(last >=0) !until dead us part (EOF)

        read(unit=un,fmt=*,iostat=last)EOP%time(sz,4),EOP%dat(sz,1:3),EOP%dat(sz,5:6),EOP%dat(sz,7:9),EOP%dat(sz,11:12)
!        write(0,*)last,EOP%time(sz,:),EOP%dat(sz,:)
        if(last<0)exit ! acceptable dead
        if(last>0)then ! reading error
           write(stderr,*)"ERROR in read_EOP(), truncated file?"
           stop
        end if
        !     write(0,*)last,sz,EOP%time(sz,:),EOP%dat(sz,:)
        sz=sz+1
     end do
     EOP%nt=sz
  case('C04') ! C04 file
     frmt='(3(F4.0),F7.0,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2F12.6)'
     read(unit=un,fmt=*)! read an additional empty header line

     sz=1
     do while(last >=0) !until dead us part (EOF)

        read(unit=un,fmt=trim(frmt),iostat=last)EOP%time(sz,:),EOP%dat(sz,:)
!        write(0,*)EOP%time(sz,:),EOP%dat(sz,:)
        if(last<0)exit ! acceptable dead
        if(last>0)then ! reading error
           write(stderr,*)"ERROR in read_EOP(), truncated file?"
           stop
        end if
        !     write(0,*)last,sz,EOP%time(sz,:),EOP%dat(sz,:)
        sz=sz+1
     end do
     EOP%nt=sz
     
  end select

  close(un)


end subroutine read_EOP

!compute for a list of timetags the residual polar motion ( relative to the IERS2010 mean polar motion
subroutine get_polarmotion(time,m1,m2,lod,subtractmeaniers)
double precision,intent(in)::time(:) ! time in mjd
double precision,intent(out)::m1(:),m2(:),lod(:) !output in radians
double precision,target::xpoly0(4),xpoly1(2)
double precision,target::ypoly0(4),ypoly1(2)
double precision, pointer::xpoly(:),ypoly(:)
logical,optional::subtractmeaniers
double precision::t0,ts,xm,ym,tt,mm1,mm2,wg,d2r
integer::n0,n1,n,stderr,i,j,k
logical::isubtractmeaniers !internal logical value whether to subtract the mean pole
isubtractmeaniers=.true.

if(present(subtractmeaniers))isubtractmeaniers=subtractmeaniers

d2r=2*acos(0.d0)/180.d0
stderr=0

if(size(time,1) .ne. size(m1,1) .or. size(time,1) .ne. size(m2,1) .or. size(time,1) .ne. size(lod,1))then
   write(stderr,*)"ERROR in get_polarmotion.f90: array size inconistent"
   stop
end if

if(EOP%nt==0)call read_EOP() ! read file if not already done so


t0=51544 !2000 ( reference time)
ts=55197 !2010 (split time)
n0=3 ! polynomial degree before 2010
n1=1 ! polynomial degree after 2010

!coefficients ( from IERS2010)
xpoly0(1)=55.974
xpoly0(2)=1.8243
xpoly0(3)=0.18413
xpoly0(4)=0.007024

ypoly0(1)=346.346
ypoly0(2)=1.7896
ypoly0(3)=-0.10729
ypoly0(4)=-0.000908

xpoly1(1)=23.513
xpoly1(2)=7.6141

ypoly1(1)=358.891
ypoly1(2)=-0.6287


!loop over all time tags
do,i=1,size(time,1)

   if(time(i)<=ts)then
      n=n0
      xpoly=>xpoly0
      ypoly=>ypoly0
   else
      n=n1
      xpoly=>xpoly1
      ypoly=>ypoly1
   end if
   
   !compute the mean pole a this time point
   xm=0.d0
   ym=0.d0
   do,j=0,n !loop over polynomial coefficients
      tt=((time(i)-t0)/365.25)**j
      xm=xm+xpoly(j+1)*tt
      ym=ym+ypoly(j+1)*tt
   end do

   !convert mas to radians
   xm=xm*d2r*1e-3/3600
   ym=ym*d2r*1e-3/3600
   
   !interpolate the EOP values from the file
   ! get the right timepoint
   do,j=1,EOP%nt
      if(EOP%time(j,4)>time(i))exit
   end do

   if(j==(EOP%nt+1) .or. j==1)then
      write(stderr,*)"ERROR in get_polarmotion.f90: time out of range"
      stop
   end if
   wg=(dble(time(i))-EOP%time(j-1,4))/(EOP%time(j,4)-EOP%time(j-1,4))

   mm1=EOP%dat(j-1,1)+wg*(EOP%dat(j,1)-EOP%dat(j-1,1))
   mm2=EOP%dat(j-1,2)+wg*(EOP%dat(j,2)-EOP%dat(j-1,2))
   lod(i)=EOP%dat(j-1,4)+wg*(EOP%dat(j,4)-EOP%dat(j-1,4))
   !convert mas to radians
   mm1=mm1*d2r/3600
   mm2=mm2*d2r/3600
   
   ! m1(i)=mm1
   ! m2(i)=mm2
   if(isubtractmeaniers)then
      
      m1(i)=mm1-xm
      m2(i)=mm2-ym
   else
      m1(i)=mm1
      m2(i)=mm2
   end if

end do


end subroutine get_polarmotion

end module EOP_tools

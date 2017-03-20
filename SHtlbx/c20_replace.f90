!!Coded by Roelof Rietbroek, Tue Apr 26 10:00:36 2011
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!! fortran routine which gives a C20 replacement value from a time series of C20


function c20_replace(time,time1,type)
implicit none
double precision,intent(in)::time
double precision,optional,intent(in)::time1
integer,intent(in)::type
double precision:: c20_replace
character(200)::dir,file
integer::un,stderr,hit,last
logical::fex
double precision::dumtime,c20,c20_sig
double precision::dumtime_next,c20_next,c20_sig_next
!defaults/initializations
un=13
stderr=0
c20_replace=0.d0

!get environment variable $RLFTLBX_DATA and adapt to the right path
call getenv('RLFTLBX_DATA',dir)
dir=trim(dir)//'/C20_replace/'

!load data from file
select case(type)
case(1)!Ries C20 replacement
   file=trim(dir)//'C20_TNries.txt'
case(2)! JIGOG C20 replacement
   file=trim(dir)//'C20_JIGOGRL02.txt'
end select

!check whether file exists
inquire(FILE=trim(file),EXIST=fex)
if(.not. fex)then
   write(stderr,*)"ERROR in:c20_replace.f90, file does not exist:"
   write(stderr,*)trim(file)
   write(stderr,*)" did you set the RLFTLBX_DATA environment variable?"
   stop
end if

!!read data
last=0
open(unit=un,file=file,status='old')

if(present(time1))then !take the average
   hit=0
   do while (last .eq. 0)
      read(unit=un,fmt=*,iostat=last)dumtime,c20,c20_sig
      if(last .ne. 0)exit
      if(dumtime >= time .and. dumtime <=time1)then
         hit=hit+1
         c20_replace=c20_replace+c20
      else if(dumtime > time1)then
         exit
      end if

   end do
   
   if(hit .ne. 0)then
      c20_replace=c20_replace/dble(hit)
   end if

else !linear interpolation
   !read first line
   read(unit=un,fmt=*,iostat=last)dumtime,c20,c20_sig
   if(dumtime > time)return !quick return when nothing was found (values stay zero)
   last=0
   do while (last .eq. 0)
      read(unit=un,fmt=*,iostat=last)dumtime_next,c20_next,c20_sig_next
      if(last .ne. 0)exit
      if(dumtime_next >= time)then
         c20_replace=(c20+(time-dumtime)*(c20_next-c20)/(dumtime_next-dumtime))
         exit
      end if
      !reset variables
      c20=c20_next
      dumtime=dumtime_next
   end do
end if

close(un)



end function c20_replace

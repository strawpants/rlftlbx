!!subroutine to read in load love numbers (possibly order dependent)
!!currently supports only farrels load love numbers
!!Coded by Roelof Rietbroek, Fri Jun  1 17:16:36 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Thu Jan  3 11:35:10 2008
!!Updated by Roelof Rietbroek, Fri May 23 10:27:35 2008
!!also incorporated PREM load love numbers
!!Updated by Roelof Rietbroek, Fri Jan 30 16:47:42 2009
!!added body love numbers as well
!! changed method of reading ( now degree 1 terms need not be present in the file)


!!Updated by Roelof Rietbroek, Mon Dec  3 14:53:01 2012
!! added loadlove numbers from Wang et al 2012
!!also allowed reading of degree only vectors (no orders)


subroutine SH_loadlove(hnm,knm,lnm,typ,frame,fileop,iso)
use SHtlbx,only:SH_lm,SH_pos
implicit none
double precision,intent(out),optional,dimension(:)::knm,hnm,lnm
double precision::h,l,k
character(*),intent(in),optional::fileop !optional file with different love numbers (typ must be set to zero)
logical,intent(in),optional::iso ! true when hnm, knm and lnm are degree dependent only
integer,intent(in)::typ,frame
integer::un,lmax,i,m,n,st,nd,last
character(200)::dir,dum,file
character(20)::srch
logical::iiso !internal isotropic request variable

iiso=.false.

if(present(iso))iiso=iso

un=13

!get the maximum degree required and initialize vectors (set to zero)
if(present(hnm))then
   hnm=0.d0
   if(iiso)then
      lmax=size(hnm,1)-1
   else
      call SH_lm(size(hnm,1),lmax,m)
   end if
end if

if(present(knm))then
   knm=0.d0
   if(iiso)then
      lmax=size(knm,1)-1
   else
      call SH_lm(size(knm,1),lmax,m)
   end if
end if

if(present(lnm))then
   lnm=0.d0
   if(iiso)then
      lmax=size(lnm,1)-1
   else
      call SH_lm(size(lnm,1),lmax,m)
   end if
end if


!get environment variable $RLFTLBX_DATA
call getenv('RLFTLBX_DATA',dir)
dir=trim(dir)//'/love/'

select case(typ)
case(1)!Farrels (1972) load love numbers
   file=trim(dir)//'farrel.love'
case(2) !PREM load love numbers
   file=trim(dir)//'PREM.love'
case(3)
   file=trim(dir)//'G-Bbody.love'
case(4)
   file=trim(dir)//'PREMbody.love'
case(5)
   file=trim(dir)//'Desaibody.love'
case(6) ! PREM from Wang et al 2012 
   file=trim(dir)//'PREM_Wang.love'
case(7) ! from Wang et al 2012
   file=trim(dir)//'ak135_Wang.love'
case(8) ! from Wang et al 2012
   file=trim(dir)//'iasp91_Wang.love'
case(0)
   if(.not. present(fileop))then
      write(*,*)'SH_loadlove: custom love number file must be specified'
      stop
   end if
   file=trim(fileop)
case default
   write(*,*)'unknown love number models specified'
   stop
end select



   !select search string for the frame type
   select case(frame)
   case(1)
      srch='CE'
   case(2)
      srch='CM'
   case(3)
      srch='CF'
   case(4)
      srch='CL'
   case(5)
      srch='CH'
   case default
      write(*,*)'unknown reference frame specified, quitting'
      stop
   end select
!write(*,*)trim(dir)//'farrel.love'
open(unit=un,file=file,status='old')
! !skip the comment lines and the degree 0 term and read in degree 1 terms
! do,i=1,40
! read(un,'(A200)')dum
! if(index(dum,trim(srch)).ne. 0)read(dum,*)n,h,l,k
! if(index(dum,'CH').ne. 0)exit !loop start reading other data
! end do

! !
! if(present(hnm))hnm(2:3)=h
! if(present(lnm))lnm(2:3)=l
! if(present(knm))knm(2:3)=k



!now read the love numbers
last=0
n=0

do while ( last >= 0 .and. n<=lmax) !exit when all necessary love numbers are read
   n=0
   read(unit=un,fmt='(A)',iostat=last)dum
   if(last .ne. 0)cycle ! end of file
   
   if(index(dum,'#') .ne. 0)cycle ! comment
   read(dum,fmt=*,iostat=last)n,h,l,k ! read data 
   !write(*,*)dum,last
   if (n>lmax)then
      exit ! quick exit
   else if(n .eq. 1)then !also read TAG at the end
      if(index(dum,trim(srch)) .eq. 0)cycle ! degree 1 found but not in the requested frame
      ! else if( n .ne. 0)then
      !    if(last .ne. 0)cycle
   end if

   !set the interval for the vectors which all have degree n

   if(iiso)then
      if(present(hnm))hnm(n+1)=h
      if(present(lnm))lnm(n+1)=l
      if(present(knm))knm(n+1)=k
   else
      st=SH_pos(n,0)
      nd=SH_pos(n,n)
      !write(*,*)n,h,l,k
      if(present(hnm))hnm(st:nd)=h
      if(present(lnm))lnm(st:nd)=l
      if(present(knm))knm(st:nd)=k
   end if

end do

close(un)


end subroutine SH_loadlove

!!Coded by Roelof Rietbroek, Tue Aug 21 14:29:08 2012
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de
!!calculates GIA potential versus Uplift relations
!!Purcell et al 2011
!!TO DO: add Wahr et al relationship?

subroutine SH_loadGIA_uplift_ratio(rat,typ)
use SHtlbx,only:SH_lm,SH_pos
implicit none
double precision,intent(out),dimension(:)::rat
integer,intent(in)::typ
double precision::ratio
integer::un,lmax,m,n,st,nd,last
character(200)::dir,dum,file
un=13

!get the maximum degree required and initialize vectors (set to zero)
rat=0.d0
call SH_lm(size(rat,1),lmax,m)



!get environment variable $RLFTLBX_DATA
call getenv('RLFTLBX_DATA',dir)
dir=trim(dir)//'/love/'

select case(typ)
case(1)! Purcell et al 2011
   file=trim(dir)//'Purcell2011GIA.ratio'
case default
   write(*,*)'unknown GIA uplift relationship specified'
   stop
end select


!write(*,*)trim(dir)//'farrel.love'
open(unit=un,file=file,status='old')

!now read the  ratios
last=0
n=0

do while ( last >= 0 .and. n<=lmax) !exit when all necessary love numbers are read
   n=0
   read(unit=un,fmt='(A)',iostat=last)dum
   if(last .ne. 0)cycle ! end of file
   
   if(index(dum,'#') .ne. 0)cycle ! comment
   read(dum,fmt=*,iostat=last)n,ratio ! read data 
   !write(*,*)dum,last
   if (n>lmax)exit
   !set the interval for the vectors which all have degree n
   st=SH_pos(n,0)
   nd=SH_pos(n,n)
!write(*,*)n,h,l,k
   rat(st:nd)=ratio
end do

close(un)


end subroutine SH_loadGIA_uplift_ratio

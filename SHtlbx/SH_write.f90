!!subroutine to write Spherical harmonic data to standard output (or a file optionally) in more formats
!!Coded by Roelof Rietbroek, Fri Jun  1 12:25:47 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Thu Mar  4 12:16:23 2010
!! allow writing of the linear format ( saves a lot of time awking)



subroutine SH_write(clm,slm,clm_sig,slm_sig,tstart,tcent,tend,comm1,comm2,comm3,filen,typ)
use SHtlbx,only:SH_lm,SH_pos
implicit none
double precision,intent(in),optional,dimension(:)::clm,slm,clm_sig,slm_sig
double precision,intent(in),optional::tstart,tcent,tend
double precision::t0,t1,t2
integer,intent(in),optional::typ
character(*),intent(in),optional::comm1,comm2,comm3,filen
character(len=150),dimension(40)::header
character(6)::str,str2!keyword
character(60)::frmt,dum
integer::un,hdr,out,ftyp,pos,lmax,l,m,i
logical:: linear
t0=0.d0
t1=0.d0
t2=0.d0
out=0
str='      '
linear=.false.

!write(*,*)'entered SH_write'  

!check if enough input vectors are provided (only two combinations are accepted)
!1: clm,slm,clm_sig,slm_sig
!2: clm,slm 
!all other combinations are wrong and yield an error message

if (present(clm_sig) .and. present(slm_sig) .and. present(slm_sig) .and. present(clm_sig))then
out=1 !full set outputted

else if (present(clm) .and. present(slm))then
out=2
else
write(*,*)'both slm_sig and clm_sig must be specified in routine SH_write'
stop
end if

!get maximum degree from vector

call SH_lm(size(clm,1),lmax,m)
!write(*,*)lmax,m
!open file if output goes to file
if(present(filen))then
    if(filen == '-')then
        un=6 !default output
    else !open file
        un=13
        open(unit=un,file=filen,status='new')
    end if
else!standard output (default)
    un=6
end if

!setup time tags
if(present(tstart))t0=tstart
if(present(tcent))t1=tcent 
if(present(tend))t2=tend


!set typ specific information
if (present(typ))then
ftyp=typ
else
ftyp=4!default typ, clean output only coefficients and degree and order
end if

select case (ftyp)
case(1)!GRACE GFZ,CSR,JPL type
write(*,*)'output type not yet supported'
stop
case(2)
write(*,*)'output type not yet supported'
stop
case(3) !ICGEM format
write(*,*)'output type not yet supported'
stop
case(4)
hdr=1 !one header line with time tag
!construct header line (double are converted to real otherwise it becomes slow for some reason)y
write(header(1)(1:),'(A6,I5,3f12.6)')'META',lmax,t0,t1,t2
if(out.eq. 1)then
   !construct output format
   frmt='(A6,Tl5,2I5,1x,e21.14,1x,e21.14,1x,e21.14,1x,e21.14)'
else
   frmt='(A6,Tl5,2I5,1x,e21.14,1x,e21.14)'
end if
case(6)
linear=.true. !one coefficient per line ( and optionally a sigma)
hdr=0 ! don't print any header lines
if(out.eq. 1)then
   frmt='(A3,1X,I3,I3,14X,e21.14,1x,e21.14)'
else
   frmt='(A3,1X,I3,I3,14X,e21.14)'
end if
str='GCN'
str2='GSN'

case default
hdr=0
if(out.eq. 1)then
   frmt='(A6,Tl5,2I5,1x,e21.14,1x,e21.14,1x,e21.14,1x,e21.14)' 
else
   frmt='(A6,Tl5,2I5,1x,e21.14,1x,e21.14)' 
end if
end select


!write header to standard output (or file)
do,i=1,hdr
write(un,*)trim(adjustl(header(i)(:)))
end do

!write(*,*)'blah'
if(out .eq. 1)then !output also
   if(linear)then
      !write coefficients to standard output(or file)
      do,m=0,lmax
         do,l=m,lmax
            pos=SH_pos(l,m)
            write(un,trim(frmt))str,l,m,clm(pos),clm_sig(pos)
            if(m>0)write(un,trim(frmt))str2,l,m,slm(pos),slm_sig(pos)
         end do
      end do
   else
      !write coefficients to standard output(or file)
      do,m=0,lmax
         do,l=m,lmax
            pos=SH_pos(l,m)
            write(un,trim(frmt))str,l,m,clm(pos),slm(pos),clm_sig(pos),slm_sig(pos)
         end do
      end do
   end if
else
   if(linear)then
      !write coefficients to standard output(or file)
      do,m=0,lmax
         do,l=m,lmax
            pos=SH_pos(l,m)
            write(un,trim(frmt))str,l,m,clm(pos)
            if(m>0)write(un,trim(frmt))str2,l,m,slm(pos)
         end do
      end do
   else
      do,m=0,lmax
         do,l=m,lmax
            pos=SH_pos(l,m)
            write(un,trim(frmt))str,l,m,clm(pos),slm(pos)
            !write(*,*)str,l,m,clm(pos),slm(pos),0.d0,0.d0
         end do
      end do
   end if
end if
!write(*,*)'blah',lmax
if(un .ne. 6)close(un)
end subroutine SH_write

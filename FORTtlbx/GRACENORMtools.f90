!!subroutines to read in GRACE normals in binary form
!!planned is also a write routine

!!Coded by Roelof Rietbroek, Tue Feb 12 14:10:51 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue Feb 10 22:17:54 2009
!!history of the STATIS entries may be retrieved

!!Updated by Roelof Rietbroek, Wed Feb  9 11:25:04 2011
!! remove fortran importable 'open' code and replace with stream access
!! add automatic byte swapping

!!Updated by Roelof Rietbroek, Wed Feb 15 17:17:33 2012
!! added better time retrieval
!! updated by Roelof Rietbroek wed 13 Sept 2017
!! only update time when valid date was read in

!!subroutine read_GRACEnorm reads in GRACE normal file header part
subroutine read_GRACEnormhead(filename,npar,npar_red,nobs,side_d,ltpl,apriori,Etime,Stime,hist,nhist)
use forttlbx,only:freeunit
use gpstlbx
implicit none
character(*),intent(in)::filename
character(24),intent(out),optional::side_d(:)
character(120),intent(out),optional::hist(:)
double precision,intent(out),optional::apriori(:),ltpl,Etime,Stime
integer,intent(out),optional::npar,npar_red,nobs,nhist

!private variables
character(120)::dum
integer::unit,i,inpar,inpar_red,inobs,reqlevel,st,inhist,maxhist,stderr
double precision::iltpl,rem,iEtime,iStime
character(12)::frmt,date1
integer::yyst,mmst,ddst,yynd,mmnd,ddnd,sec
Double precision::jd,jd0,diy
logical::nd
!defaults initializations
stderr=0
inpar=0
inpar_red=0
inobs=0

iStime=9e6
iEtime=0
inhist=0
unit=freeunit()
nd=.false. ! end of header
if(present(hist))maxhist=size(hist,1)

!!!!!!!!header part!!!!!!!!!!!
open(unit=unit,file=trim(filename),form='formatted')

!!!loop over header lines

do while( .not. nd)
   read(unit,'(A)')dum
   
   select case(dum(1:6))
      case('STATIG')
         read(dum(7:),*)inobs,inpar,inpar_red,iltpl
      case('STATIS')
         inhist=inhist+1
         if(present(hist))then
            if(maxhist < inhist)then
               write(stderr,*)"ERROR: history array is too small",inhist
               stop
            end if
            hist(inhist)=dum
         end if

         read(dum,'(33x,i2,i2,i2,1x,i2,i2,i2)')yyst,mmst,ddst,yynd,mmnd,ddnd
         !only update time if we have read a valid date from the STATIS line   
         if( (mmst*ddst*mmnd*ddnd) > 0 )then
            !retrieve time info
            !start time
            if(yyst >50)then
                jd=date_2_jd(ddst,mmst,yyst+1900,0,0,0)
            else
                jd=date_2_jd(ddst,mmst,yyst+2000,0,0,0)
            end if

            iStime=min(iStime,jd)


            !end time
             if(yynd >50)then
                jd=date_2_jd(ddnd,mmnd,yynd+1900,0,0,0)
            else
                jd=date_2_jd(ddnd,mmnd,yynd+2000,0,0,0)
            end if

            iEtime=max(iEtime,jd)
        end if  



      case('NORMAL')
         read(unit,'(A)')dum
         nd=.true. !end of header encountered ( loop stops)
   end select

end do

!!convert julian day to decimal year
call jd_2_date(iStime,ddst,mmst,yyst,sec)
        ! write(0,*)ddst,mmst,yyst
jd0=date_2_jd(1,1,yyst,0,0,0)
diy=date_2_jd(1,1,yyst+1,0,0,0)-jd0
iStime=yyst+(iStime-jd0)/diy

!!convert julian day to decimal year
call jd_2_date(iEtime,ddnd,mmnd,yynd,sec)
!write(0,*)ddnd,mmnd,yynd
jd0=date_2_jd(1,1,yynd,0,0,0)
diy=date_2_jd(1,1,yynd+1,0,0,0)-jd0
iEtime=yynd+(iEtime-jd0)/diy


!put values into variables if requested
if(present(npar))npar=inpar
if(present(npar))npar_red=inpar_red
if(present(ltpl))ltpl=iltpl
if(present(nobs))nobs=inobs
if(present(Stime))Stime=iStime
if(present(Etime))Etime=iEtime
if(present(nhist))nhist=inhist

!force premature exit if only indices are requested
reqlevel=0
!data from the header file
if(present(side_d))reqlevel=1
if(present(apriori))reqlevel=1

if(reqlevel > 0)then !only go through this bit if requested
   !get side description
   if(present(side_d))then
      st=0
      do,i=1,floor(inpar/5.d0)
         
         read(unit,'(5A)')side_d(st+1:st+5)

         st=st+5
      end do
      !read last part
      if(st < inpar)then
         write(frmt,'(a1,i1,a2)')'(',inpar-st,'A)'
         
         read(unit,frmt)side_d(st+1:inpar)
      end if
      
   else !read in dummies
      do,i=1,ceiling(inpar/5.d0)
         read(unit,*)dum(1:3)
      end do
   end if
   
   
   
   !get apriori vector
   if(present(apriori))then
      st=0
      do,i=1,floor(inpar/6.d0)
         read(unit,'(6F20.14)')apriori(st+1:st+6)
         st=st+6
      end do
      !read last part
      if(st < inpar)then
         write(frmt,'(a1,i1,A7)')'(',inpar-st,'F20.14)'
         read(unit,frmt)apriori(st+1:inpar)
      end if
   else !read in dummies
      do,i=1,ceiling(inpar/6.d0)
         read(unit,*)dum(1:3)
      end do
   end if
end if

close(unit)

end subroutine read_GRACEnormhead


!!!subroutine to read in the binary part of the data

subroutine read_GRACEnormbin(filename,npar,mat,pmat,Atb,ltpl)
use forttlbx,only:ltlend
implicit none
integer,intent(in)::npar
character(*),intent(in)::filename
double precision,intent(out),optional::pmat(:),mat(:,:),ltpl,Atb(:)


!private arguments
integer::unit,i,st,len_sysdep
double precision::iAtb(npar) !automatic array
integer,parameter::kinddbl=8 ! amount of bytes for a double
integer::opencstream

logical::needswap

needswap=ltlend() !byte swapping is necessary when the local machine is little endian

!open direct access file as a stream for reading
unit=opencstream(0,trim(filename)//char(0))

!!unit=freeunit()
!!!!note the record length (unit) is sytem dependent!!!!!!!!!!
!!we can find the required record length with the following test statement
!!thus find out what the system dependent record length should be for writing a double vector of length npar+2
!inquire(IOLENGTH=len_sysdep)iAtb(:),iAtb(1:2)

!example the binary files were written on sun machines were the record length is described in the amount of bytes (npar+2)*8
!on linux systems however the record length must be specified in the amount of 4byte words thus (npar+2)*2

!open(unit=unit,file=trim(filename),form='unformatted',convert='BIG_ENDIAN',access='DIRECT',recl=len_sysdep)

    

if(present(mat))then
   do,i=1,npar
      call cread(unit,kinddbl,iAtb(i)) ! read one double from the stream
      call cread(unit,kinddbl*npar,mat(i,1)) ! read npar doubles in row i of mat
      call cskip(unit,kinddbl) !skip unused value
!      read(unit=unit,rec=i)iAtb(i),mat(i,1:npar)
   end do
end if


if(present(pmat))then
   st=0
   do,i=1,npar
      call cread(unit,kinddbl,iAtb(i)) ! read one double from the stream
      if(needswap)call cswap(kinddbl,1,iAtb(i))
      
     call cread(unit,kinddbl*i,pmat(st+1)) ! read i doubles in pmat
     if(needswap)call cswap(kinddbl,i,pmat(st+1))
     !skip symmetric values and final unused value until next record
     call cskip(unit,kinddbl*(npar-i+1))

      !read directly into packed matrix (use symmetry properties)
!      read(unit=unit,rec=i)iAtb(i),pmat(st+1:st+i)
     
      st=st+i
   end do
end if



if(present(Atb))Atb=iAtb

if(present(ltpl))then
   call cread(unit,kinddbl,ltpl)
   if(needswap)call cswap(kinddbl,1,ltpl)
end if

call cclose(unit)

! if(present(ltpl))read(unit=unit,rec=npar+1)ltpl
! close(unit)


end subroutine read_GRACEnormbin



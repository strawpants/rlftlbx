!set of subroutines to create sinex files
!important note: the current version is by no means compatible with the official sinex version because many mandatory block are still missing
!!Coded by Roelof Rietbroek, Tue Aug 21 10:34:27 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de


!!Updated by Roelof Rietbroek, Thu May 22 16:01:10 2008
!! also writes to standard output if no file name is specified


subroutine init_sinex(file,unit,line)
use GPStlbx,only: GPS_date
implicit none
integer,intent(out)::unit
character(*),intent(in),optional::file
character(52),intent(in)::line !to be consistent with the sinex format line must be 52 characters long
!and should contain the following format (example):
!NRC 06:301:00000 06:310:00000 P 00875 2 X E
!columns consistent with
!1. Institute providing 2 data data start time 3. data end time 4. observation code 5.amount of unknowns (as characters) 6. constraint code 7. provided data types

logical::opennd
integer::i
character(12)::sindate
character(80)::frstline,dum
!open file unit

!defaults
unit=6 !standard output is default

if(present(file))then
do,i=7,99 !try unit numbers
   inquire(UNIT=i,OPENED=opennd)
   if(.not. opennd)then
      unit=i
      open(unit=unit,file=trim(file),status='new')
      exit!exit loop when an unused unit number was found
   end if
end do
end if


write(frstline,'(a15,a12,1x,a52)')'%=SNX 2.02 GFZ ',GPS_date(),line

!write to file
write(unit,'(a80)')frstline
dum='*Experimental sinex file still not compliant with the full sinex description'
write(unit,'(a80)')adjustl(dum)
dum='*-----------------------------------------------------------------------'
write(unit,'(a80)')adjustl(dum)
end subroutine init_sinex

!subroutine which writes a sinex block named sblock to an open unit
subroutine write_sinexblock(unit,sblock,block,comm)
  implicit none
  integer,intent(in)::unit
  character(*),intent(in)::sblock
  character(*),intent(in),optional::comm
  character(80),dimension(:),intent(in)::block
  character(80)::dum
  integer::stderr,i,nlines
  logical::opennd

!defaults
stderr=0
!inquire whether unit is opened already  
inquire(UNIT=unit,OPENED=opennd)
if(.not. opennd)then
   write(stderr,*)'sinex file not initialized yet'
   stop
end if
!get amount of data lines
nlines=size(block,1)
dum='*-----------------------------------------'
write(unit,'(a80)')adjustl(dum)
dum='+'//trim(sblock)
write(unit,'(a80)')adjustl(dum)

if(present(comm))then
   dum=comm
   write(unit,'(a80)')adjustl(dum)
end if

do,i=1,nlines
   write(unit,'(a80)')block(i)
end do
dum='-'//trim(sblock)
write(unit,'(a80)')adjustl(dum)
end subroutine write_sinexblock


!!subroutine which writes a lower triangular SOLUTION/MATRIX_ESTIMATE SINEX block
subroutine write_sinex_MATlow(unit,mat,matpack,type)
implicit none
integer,intent(in)::unit
double precision,intent(in),optional,dimension(:)::matpack
double precision,intent(in),optional,dimension(:,:)::mat
character(4),intent(in)::type !either COVA, CORR or INFO
integer::row,col,nmat,stderr
logical::pack,opennd
character(80)::frmt1,frmt2,frmt3,dum


!initialize
stderr=0
pack=.false.

if(present(matpack))pack=.true.
if(present(mat).and. pack)then !do an input check
   write(stderr,*)'Only one of mat or matpack can be specified, quitting'
   stop
end if

if(.not.present(mat) .and. .not.present(matpack))then
   write(stderr,*)'Either one of mat or matpack must be specified, quitting'
   stop
end if

!get size of matrix
if(pack)then
   nmat=size(matpack,1)
   nmat=(int(sqrt(8*nmat+1.d0)+0.1)-1)/2
else
   nmat=size(mat,1)
end if


!define sinex formats
frmt1='(1x,i5,1x,i5,1x,e21.14)'
frmt2='(1x,i5,1x,i5,1x,e21.14,1x,e21.14)'
frmt3='(1x,i5,1x,i5,1x,e21.14,1x,e21.14,1x,e21.14)'


!check if file is opened
inquire(UNIT=unit,OPENED=opennd)
if(.not. opennd)then
   write(stderr,*)'sinex file not initialized yet'
   stop
end if

dum='*-----------------------------------------'
write(unit,'(a80)')adjustl(dum)
!open block
dum='+SOLUTION/MATRIX_ESTIMATE L '//trim(type)
write(unit,'(a80)')adjustl(dum)
dum='*Indx1 Indx2 ____Indx2+0__________ ____Indx2+1__________ ____Indx2+2__________'
write(unit,'(a80)')adjustl(dum)
!loop over entries


if(pack)then !if packed matrix


   do,row=1,nmat
      col=1
      do while(col < row-1)
         write(unit,frmt3)row,col,matpack(row+(col*(col-1))/2),&
              matpack(row+((col+1)*(col))/2),matpack(row+((col+2)*(col+1))/2)
         col=col+3 !increment column
      end do
   
      !write incomplete line
      if(col .eq. row)then
         write(unit,frmt1)row,col,matpack(row+(col*(col-1))/2) !only one entry left
      else if(row-col .eq. 1)then
         write(unit,frmt2)row,col,matpack(row+(col*(col-1))/2),&
              matpack(row+((col+1)*(col))/2) !two entries left
         !else do nothing
      end if
   end do


else
   do,row=1,nmat

      col=1
      do while(col < row-1)
         write(unit,frmt3)row,col,mat(row,col),mat(row,col+1),mat(row,col+2)
         col=col+3 !increment column
      end do
      
      !write incomplete line
      if(col .eq. row)then
         write(unit,frmt1)row,col,mat(row,col) !only one entry left
      else if(row-col .eq. 1)then
         write(unit,frmt2)row,col,mat(row,col),mat(row,col+1) !two entries left
         !else do nothing
      end if
   end do
end if






!close block
dum='-SOLUTION/MATRIX_ESTIMATE L '//trim(type)
write(unit,'(a80)')adjustl(dum)
end subroutine write_sinex_MATlow











subroutine end_sinex(unit)
implicit none
integer,intent(in)::unit
character(80)::dum
dum='*----------------------------'
write(unit,'(a80)')adjustl(dum)
dum='%ENDSNX'
write(unit,'(a80)')adjustl(dum)

close(unit)
end subroutine end_sinex

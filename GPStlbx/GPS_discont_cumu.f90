!file which contains subroutines and functions to handle igs discontinuity files and cumulative files
!general approach is an initialization on first call (reading in of files and then saves workspace variables to memory ready for additional calls)


!function which finds out the solution id of the requested week and station 
!returns a search string with the appropriate solution id incorporated
function get_solu_id(file,station,gpswk,minwks)
use GPStlbx,only :compstr, get_sinexblock, GPS_week, GPS_date
implicit none
character(*),intent(in)::file
character(*),intent(in)::station
integer,intent(in)::gpswk
character(12)::get_solu_id
!optional variables
integer,optional,intent(in)::minwks !minimum interval for a solution id

!local variables
integer::ndat,i,date1,date2,stderr
character(120)::oldfile
character(80),dimension(:),allocatable::discon
save oldfile,ndat,discon

!defaults:
stderr=0

!on first run
if(trim(oldfile) .ne. trim(file))then
   oldfile=trim(file)!reset file name
!write(*,*)trim(file),trim(oldfile)
   call get_sinexblock(file=trim(file),blockname='SOLUTION/DISCONTINUITY',nlines=ndat)
 ! write(*,*)ndat
   if(allocated(discon))then
      deallocate(discon)
   end if
   allocate(discon(ndat))
   call get_sinexblock(File=trim(file),Blockname='SOLUTION/DISCONTINUITY',List=discon)

   

end if !end if block of first run

!after initialization

!do loop over stations and search for input station

!set default (remains if not found)
get_solu_id='############'

do,i=1,ndat
  ! write(*,*)discon(i),station
   if(.not.compstr(discon(i),'*'//station//'*'))cycle !station not found go to next entry
   
   !exclude those entries which have a velocity tag
   if(compstr(discon(i),'* V -*'))cycle 
!retrieve interval of entry
   date1=GPS_week(discon(i)(17:28))
   date2=GPS_week(discon(i)(30:41))

   if(date2 .eq. 0)then !set date to current date
      date2=GPS_week(GPS_date())
   end if
   

  

   !check whether requested week is outside the time interval
!also exclude when discontinuity occurs in this week
   if((date1 >= gpswk) .or. (date2 <= gpswk))cycle

   !check if interval is longer than minimum
   if(present(minwks))then
      if((date2-date1) < minwks)then
         write(stderr,*)'Solution period shorter than ',minwks,' weeks : ',discon(i)(1:13)
         cycle
      end if
   end if
!write(*,*)discon(i)(17:28),date1,gpswk,discon(i)(30:41),date2,trim(station)
!also exclude when discontinuity occurs in this week
 

   !when arrived here we have the right entry
   get_solu_id=discon(i)(2:13)
   exit !
end do


end function get_solu_id


!subroutine which loads a cumulative SOLUTION/ESTIMATE block into memory and returns position nand velocity depending on search string
subroutine get_cumdata(file,inc,pos,vel,reftime,error)
use GPStlbx,only:get_sinex,compstr,GPS_year,GPS_week
implicit none
character(*),intent(in)::file,inc
integer,intent(out)::error
double precision, intent(out)::pos,vel,reftime


!local variables
integer::ndat,i,ntest
character(120)::oldfile
character(80),dimension(:),allocatable::cumdat
double precision,allocatable,dimension(:)::posdat
integer,allocatable,dimension(:)::tag

save oldfile,ndat,cumdat,tag,posdat

!on first run
if(trim(oldfile) .ne. trim(file))then
   oldfile=trim(file)!reset file name
   call get_sinex(file=trim(file),nx=ndat,inc1='* STA*|* VEL*')
   allocate(cumdat(ndat),posdat(ndat),tag(ndat))
   tag=0
   call get_sinex(file=trim(file),nx=ndat,x=posdat,inc1='* STA*|* VEL*',xlist=cumdat)
   

!analyze block and tag
   do,i=1,ndat
      call get_sinex(file=trim(file),nx=ntest,inc1='*STA*'//cumdat(i)(15:26)//'*|*VEL*'//cumdat(i)(15:26)//'*')
   if(ntest .eq. 6)tag(i)=1!all parameters present
   end do
   
end if !end if block of first run

!after initialization
!dummy values
pos=1.d20
vel=1.d20
error=-2 !default not found at all
!do loop over stations and search for input string



do,i=1,ndat
   if(.not.compstr(cumdat(i),inc))cycle !station not found go to next entry

!if found get data
   if(tag(i).eq. 0)then
      error=-1
      exit !not all 3D data is present in the file (position and velocity in XYZ)
   end if
   !retrieve components
   if(compstr(cumdat(i),'* STA*'))pos=posdat(i)
   if(compstr(cumdat(i),'* VEL*'))vel=posdat(i)
   
   if((vel < (1.d20)) .and. (pos < (1.d20)))then !both values have been found and set
      reftime=GPS_year(GPS_week(cumdat(i)(27:40)))
      !write(*,*)reftime,cumdat(i)(27:40),GPS_week(cumdat(i)(27:40))
      error=0
      exit
   end if
      
end do




end subroutine get_cumdata

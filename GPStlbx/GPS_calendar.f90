!command line tool to conert GPS weeks to dates and vice versa
!!Coded by Roelof Rietbroek, Thu Aug 30 13:52:03 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Thu Oct 18 12:13:27 2007
!!Updated by Roelof Rietbroek, Fri Nov  2 11:55:51 2007
!!added -c option (current date)
!!round output to seconds
!!fixed bug on command line option (not all digits where read in, causing round off errors in seconds)
!!fixed bug in sinex string format ( doy starts with :001: instead of :000:)
!!Updated by Roelof Rietbroek, Mon Dec 10 16:07:48 2007
!!The program can now aslo read from standard input (useful for piping)
!! example: echo '2006.874' | GPS_calendar now also works


!!GPStlbxlibrary needed

!!Updated by Roelof Rietbroek, Thu Sep 24 16:47:15 2009
!!convert evevrything to Julian day first
!!added Julian day

!!Updated by Roelof Rietbroek, Tue Mar 23 14:06:11 2010
!!allow multiple inputs and added cycles for JASON,TOPEX 

!!Updated by Roelof Rietbroek, Fri Apr 29 11:12:24 2011
!! when reading from standard input also output further columns without conversion

!!Updated by Roelof Rietbroek, Mon May 16 10:34:18 2011
!! for single column output remove the printed space at the beginning (compiler specific)

!!Updated by Roelof Rietbroek, Fri Jun  1 13:15:53 2012
!! added JASON 2 to satellites


program GPS_calendar
use GPStlbx
use FORTtlbx
implicit none
integer::i,narg,input,output,wk
double precision::dummy,dwk
character(12)::date,time
character(200)::dum
integer::iargc,itharg,stderr
integer::ind,last,n,j
integer::yr,dy,sec,hr,min,y0,y1,mn
double precision::jd,diy,jd0,jd1
integer::nsat,satid
parameter(nsat=4)
double precision,dimension(nsat)::cyc_len,cyc_start,cyc_end
character(10),dimension(nsat)::satn
character(200)::datein
character(40)::datetmp
integer::st,nd
!defaults initializations
input=0
output=3 !output format 3 stands for truncated GPSweek
stderr=0
sec=0
n=0
st=0
nd=0

!altimeter cycle lengths and starts ( in days)
!ENVISAT
cyc_len(1)=35.0
cyc_start(1)=2452368.39881944d0
cyc_end(1)=cyc_start(1)+365*100 ! plus 100 years
satn(1)='ENVISAT'

!TOPEX
cyc_len(2)=9.91564189814815d0
cyc_start(2)=2448889.46068287d0
cyc_end(2)=2453471.16363426d0
satn(2)='TOPEX'

!JASON1
cyc_len(3)=9.91564189814815d0
cyc_start(3)=2452289.71670139d0
cyc_end(3)=cyc_start(3)+365*100 ! plus 100 years
satn(3)='JASON1'

!JASON2
cyc_len(4)=9.91564189814815d0
cyc_start(4)=2454659.5555555555d0
cyc_end(4)=cyc_start(3)+365*100 ! plus 100 years
satn(4)='JASON2'



!process command line argument

narg=iargc()


itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg>narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq.'-')then !option
      select case(dum(2:2))

      case('s')!output A SINEX DATE STRING
         output=1
      case('d')!output a dd-mm-yyyy string
         output=2
      case('y')!output in year format
         output=4
      case('c')!use current date
         n=n+1
         call date_and_time(date=date,time=time)
         read(date,'(i4,i2,i2)')yr,mn,dy
         read(time,'(i2,i2,i2)')hr,min,sec
         jd=date_2_jd(dy,mn,yr,hr,min,sec)
         write(datein,*)"JD",jd
      case('f')! output in fractional GPS week
         output=5
      case('j')
         output=6
      case('m') !modified julian day
         output=7
      case('a')!output as altimeter mission cycle
         output=8
         satid=0
         do,j=1,nsat
            if(index(dum,trim(satn(j)))>0)then
               satid=j
               exit
            end if
         end do

         if(satid .eq. 0)then
            write(stderr,*)"ERROR, processing -a option, mispelled satelitte?"
            stop
         end if
      case('D')!for multicolumn input describe the start and end index of the datestring within the line
         ind=index(dum,"/")
         read(dum(4:ind-1),*)st
         read(dum(ind+1:),*)nd

      case default
         write(*,*)'unknown option selected'
         call help(satn,nsat)
      end select
   else !no option but date
      datein=trim(dum) !copy line to string
      n=1 !means input is set
   end if
end do


!if no dates are given read list from standard input (one per line)
! if(n .eq. 0)then
!    last=0
!    do while (last .eq. 0)
!       read(unit=*,fmt='(a40)',iostat=last)dum
!       if(last .ne. 0)exit
!       n=n+1
!       call realloc_ptr(datein,1)
!       datein(n)=trim(dum)
!    end do
! end if


!loop over input dates /lines of a file
last=0
do while(last .eq. 0)
   
   !read next line if info is coming from standard input
   if(n.eq. 0)then
      read(unit=*,fmt='(a200)',iostat=last)datein
      if(last .ne. 0)exit ! exit loop: end of file reached
   end if

   ! extract relevant date part from string
   if(st.ne.nd)then
      datetmp=datein(st:nd)
   else
      datetmp=trim(datein)
   end if

   !get input type from string
   if(index(datetmp,'-') .ne. 0)then
      input=2  
   else if(index(datetmp,':') .ne. 0)then !sinex format
      input=1

   else if (index(datetmp,'JD').ne. 0)then ! input is a Julian day
      input=6
      if(index(datetmp,'MJD') .ne. 0)input=7 ! set input to modified Julian day
      ind=index(datetmp,'JD')
      
      datetmp=adjustl(datetmp(ind+2:))
   else if(index(datetmp,'.').ne. 0)then
      input=4
   else if(len(trim(datetmp))>5)then ! Altimeter cycle
      input=8
   else ! input is a GPS week
      input=3
      if(output .eq. 3)output=2 !set output to normal date format
   end if


   !now convert everything in a julian day
   select case(input)
   case(1) !convert sinex string to JD
      read(datetmp,'(i2.2,1x,i3.3,1x,i5.5)')yr,dy,sec
      if(yr<50)then
         yr=yr+2000
      else
         yr=yr+1900
      end if
      jd=date_2_jd(1,1,yr,0,0,0)+dy-1+dble(sec)/86400.d0
   case(2) ! convert date to julian day
      ind=index(datetmp,'-')
      read(datetmp(ind-2:),'(i2.2,1x,i2.2,1x,i4.4)')dy,mn,yr
      
      ind=index(datetmp,':')!check if a time is supplied
      if(ind .ne. 0)then
         
         read(datetmp(ind-2:),'(i2.2,1x,i2.2,1x,i2.2)')hr,min,sec
      else
         hr=0
         min=0
         sec=0
      end if
      jd=date_2_jd(dy,mn,yr,hr,min,sec)
         
   case(3)! input is in GPS weeks
      read(datetmp,*)dummy
      jd0=date_2_jd(6,1,1980,0,0,0) ! start of the GPS week system
      jd=jd0+(dummy)*7+3.5 !GPS week also contained a 0 week)
   case(4) ! input is in years
      read(datetmp,*)dummy
      y0=int(dummy)
      !calculate the amount of valid days in the year
      jd0=date_2_jd(1,1,y0,0,0,0)
      jd1=date_2_jd(1,1,y0+1,0,0,0)
      diy=jd1-jd0
      jd=jd0+(dummy-y0)*diy
   case(6)! input is already a julian day
      read(datetmp,*)jd
   case(7)
      read(datetmp,*)jd
      jd=jd+2400000.5
   case(8) ! altimeter cycle ( use center of the cycle as the reference point)
      do,i=1,nsat
         ind=index(datetmp,trim(satn(i)))
         if(ind>0)then
            read(datetmp(ind+len(trim(satn(i))):),*)dummy
            jd=cyc_start(i)+cyc_len(i)*(dummy-1.0)+cyc_len(i)/2.d0
            if(jd < cyc_start(i) .or. jd> cyc_end(i))then
               write(stderr,*)"WARNING:cycle out of range for altimeter",satn(i)
               cycle 
            end if
            exit
         end if
      end do


   end select


   !output to screen
   select case(output)
   case(1) ! output as sinex string
      call jd_2_date(jd,dy,mn,yr,sec)
      jd0=date_2_jd(1,1,yr,0,0,0)
      dy=jd-jd0+1
      if(yr<=1950 .or. yr >2050)then
         write(stderr,*)"WARNING:out of range of sinex date"
         cycle
      end if

      !strip century from year
      yr=yr-int(yr/100)*100
      write(dum,'(i2.2,1a,i3.3,1a,i5.5)')yr,':',dy,':',sec
      
   case(2) ! output as a gregorian date
      call jd_2_date(jd,dy,mn,yr,sec)
      hr=sec/3600
      min=sec/60-hr*60
      sec=sec-min*60-hr*3600
      write(dum,'(i2.2,1a,i2.2,1a,i4.4,1x,i2.2,1a,i2.2,1a,i2.2)')&
           dy,'-',mn,'-',yr,hr,':',min,':',sec
   case(3)!output as a GPS week
      jd0=date_2_jd(6,1,1980,0,0,0) ! start of the GPS week system
      wk=(jd-jd0)/7
      if(wk < 0 .or. wk >9999)then
         write(stderr,*)"WARNING:out of range of GPSweeks"
         cycle
      end if
      write(dum,'(I4.4)')wk
   case(4)!output as decimal year
      call jd_2_date(jd,dy,mn,yr,sec)
      jd0=date_2_jd(1,1,yr,0,0,0)
      diy=date_2_jd(1,1,yr+1,0,0,0)-jd0
      write(dum,'(F14.9)')yr+(jd-jd0)/diy
   case(5)
      jd0=date_2_jd(6,1,1980,0,0,0) ! start of the GPS week system
      dwk=(jd-jd0)/7.0
      write(dum,'(F14.9)')dwk
   case(6)
      write(dum,*)'JD',jd
   case(7)
      write(dum,*)'MJD',jd-2400000.5
   case(8) ! output altimeter cycle
      if(jd < cyc_start(satid) .or. jd > cyc_end(satid))then
         write(stderr,*)"WARNING:cycle out of range for altimeter",satn(satid)
         cycle 
      end if
      write(dum,*)trim(satn(satid)),int((jd-cyc_start(satid))/cyc_len(satid))+1
   end select
   
   !output dum string with original data
   if(st .ne. nd)then
      write(*,*)datein(1:st-1),trim(dum),trim(datein(nd+1:))
   else
      write(*,'(A)')trim(dum)
   end if

   if(n .eq. 1)exit !quick exit from loop when only 1 date is read from the command line

end do


end program GPS_calendar

!subroutine to plot help message
subroutine help(satn,n)
implicit none
integer::n,i
character(10)::satn(n)
write(*,*)'Program GPS_calendar converts dates'
write(*,*)'usage GPS_calendar OPTIONS [DATE]'
write(*,*)'OPTIONS can be:'
write(*,*)'  -s: Output as a sinex date string yy:doy:sssss'
write(*,*)'  -d: Output as a normal date string dd-mm-yyyy'
write(*,*)'  -y: Output as fractional year value eg 2000.87'
write(*,*)'  -c: Use current date and time (no DATE required)'
write(*,*)'  -f: Output in fractional GPS week'
write(*,*)'  -j: Output in Julian days (starts with the code JD)'
write(*,*)'  -m: Output in Modified Julian days (starts with the code MJD)'
write(*,*)'  -aSAT: Output as altimeter mission cycle. SAT can be one of:'
write(*,*)'   ',satn
write(*,*)'  -D=ST/ND: For multiple column data. Only convert the date '
write(*,*)'   which starts at character index ST and ends at ND'
write(*,*)'   Other columns remain unaffected and are simply outputted'
write(*,*)''
write(*,*)'The above formats are also autmatically recognized from the command line'
write(*,*)' Also a GPSweek can be provided on the command line (without decimal point)'
write(*,*)'Default conversions:'
write(*,*)' sinex string,date string and frac. year are converted to GPSweek'
write(*,*)' GPSweek is converted to normal date string (center of week)'
write(*,*)'DATES is an argument list containing dates in the formats described'
write(*,*)'in the options above. Spaces need to be protected. If no DATES are provided'
write(*,*)' The program will from standard input ( one date per line'
stop
end subroutine help

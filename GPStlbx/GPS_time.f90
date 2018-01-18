!function to convert a date string found in sinex file to a GPSweek number
function GPS_week(dtstr,rem)
implicit none
double precision,optional,intent(out)::rem
character(*),intent(in)::dtstr
integer::GPS_week
integer::col1,col2,yr,dy,scnds,leapyrs,comyrs,totdays

!get the indices of the double colons in the string
col1=index(dtstr,':')
col2=index(dtstr(col1+1:),':')+col1
!write(*,*)dtstr,col1,col2
!read in year day and seconds

read(dtstr(col1-2:col1-1),'(i2.2)')yr
read(dtstr(col1+1:col2-1),'(I3.3)')dy
read(dtstr(col2+1:),'(i5.5)')scnds
!write(*,*)dtstr(1:col1-1),'&',dtstr(col1+1:col2-1),'&',dtstr(col2+1:),yr,dy,scnds
!quick return when 00:000:00000 string is encountered
if((yr .eq. 0)  .and. (dy .eq. 0) .and. (scnds .eq. 0))then
GPS_week=0
return
end if

!correct year for the right century
if(yr <= 50)then
   yr=yr+2000
else
   yr=yr+1900
end if


!convert to GPS_week

!get amount of completed leap years since 1980
!integer division
leapyrs=(yr-1980+3)/4


!complete common years 
comyrs=yr-1980-leapyrs

!calculate amount of days since 6 january 1980

totdays=leapyrs*366+comyrs*365+dy-6 !last term accounts for first few days of 1980

!integer division to obtain GPSweek number
GPS_week=totdays/7

if(present(rem))then
rem=dble(totdays+dble(scnds)/86400.d0)/7.d0
end if

end function GPS_week

!function which returns a decimal year from a GPSweek number
function GPS_year(gpswk,dwk)
implicit none
integer,intent(in)::gpswk
double precision,intent(in),optional::dwk !optional GPSweek as double 
double precision::GPS_year,dys,dys2,dysinyr
integer::yr

if(present(dwk))then
   dys=7.d0*dwk !calculate days with doubkle precision GPSweek
else
   dys=dble(7*gpswk)+3.5 !center time of GPSweek
end if


!find out year
dys2=361.d0
if(dys<dys2)then
GPS_year=1980.d0+(dys+5)/366.d0
return
end if


do,yr=1981,2050
   if(mod(yr,4).eq. 0)then !leap year
      dysinyr=366.d0
   else
      dysinyr=365.d0
   end if

   if(dys <(dys2+dysinyr))then
      GPS_year=yr+(dys-dys2)/dysinyr
!      write(*,*)yr,dys-dys2
      return
   end if
      dys2=dys2+dysinyr
end do

end function GPS_year

!returns a date in sinex format: yy:doy:sssss
!input time in years eg. 2000.769
function sinex_date(date)
implicit none
double precision,intent(in)::date
double precision::dy
character(12)::sinex_date
integer::yr,dpy,doy,sec
dpy=365
yr=int(date)
if(mod(yr,4) .eq. 0)dpy=366

dy=(date-yr)*dpy

doy=int(dy)
sec=int((dy-doy)*86400.d0+0.5d0) 
if(sec .eq. 86400)then
   doy=doy+1
   sec=0
end if

if(yr >= 2000)then
   yr=yr-2000
else
   yr=yr-1900
end if

write(sinex_date,'(i2.2,a1,i3.3,a1,i5.5)')yr,':',doy+1,':',sec

end function sinex_date

!get year,momth and day and seconds from decimal year 
subroutine getYYMMDDSS(decyr,yr,mn,dd,sec)
implicit none
integer,intent(out)::yr,mn,dd,sec
double precision,intent(in)::decyr
integer::mndy(12)
integer::i,dpy
data mndy/31,28,31,30,31,30,31,31,30,31,30,31/


yr=int(decyr)
if(mod(yr,4) .eq. 0)then
    mndy(2)=29
    dpy=366
else
    dpy=365
end if

dd=int((decyr-yr)*dpy)

sec=((decyr-yr)*dpy-dd)*86400

mn=0
!now find month 
do,i=1,12
   if(dd==mndy(i))then
      mn=i
      exit
   end if
   !subtract the amount of days in the month considered
   dd=dd-mndy(i)
end do


end subroutine



function getMonth(decyr)
implicit none
double precision,intent(in)::decyr
integer:: getMonth
integer::yr,dd,sec !dummy vars

write(0,*)decyr
!call getYYMMDDSS(decyr,yr,getMonth,dd,sec)


end function getMonth

!!parse a string with time tags
!function parseDecYrStr(frmt,decyr)
!implicit none
!character(400)::parseDecYrStr
!character(*),inent(in)::frmt
!double precision,intent(in)::decyr
!integer::i,ind,strln
!integer::iocur,iocurprev
!integer::yr,mm,dd,sec

!strln=len_trim(frmt)
!call getYYMMDDSS(decyr,yr,mm,dd,sec)


!!find and substitute occurences of %yy %mm %dd %ss


!iocur=1
!iocurprev=1

!do while(ind < strlen)
    !iocur=ind+index(frmt(ind:),'%')
    !parseDecYrStr=trim(parseDecYrStr)//frmt(:ind-1    

!end do



!end function

!function which returns the current time or argument ( dd-mm-yyyy)  as a sinex date string
function GPS_date(dtstr)
implicit none
character(*),optional::dtstr
character(12)::GPS_date,date,time
integer::yr,mnth,dy,hr,min,sec,doy

!get current time
if(present(dtstr))then
   dtstr=adjustl(dtstr)
   !write(*,'(a9)')dtstr
   read(dtstr,'(i2,1x,i2,1x,i4)')dy,mnth,yr
   !write(*,*)dy,mnth,yr
   hr=0
   min=0
   sec=0
   if(yr >= 2000)then
      yr=yr-2000
   else
      yr=yr-1900
   end if

else
   call date_and_time(date=date,time=time)
   read(date,'(2x,i2,i2,i2)')yr,mnth,dy
   read(time,'(i2,i2,i2)')hr,min,sec
end if

   ! !subtract day since 0 january does not exist
!    dy=dy-1

!construct sinex date string

select case(mnth)
case(1)
   doy=dy
case(2)
   doy=31+dy
case(3) 
   doy=31+28+dy
case(4) 
   doy=31+28+31+dy
case(5) 
   doy=31+28+31+30+dy
case(6) 
   doy=31+28+31+30+31+dy
case(7) 
   doy=31+28+31+30+31+30+dy
case(8) 
   doy=31+28+31+30+31+30+31+dy
case(9) 
   doy=31+28+31+30+31+30+31+31+dy
case(10) 
   doy=31+28+31+30+31+30+31+31+30+dy
case(11)
   doy=31+28+31+30+31+30+31+31+30+31+dy
case(12)
   doy=31+28+31+30+31+30+31+31+30+31+30+dy
end select

if(mod(yr,4).eq. 0 .and. mnth > 2)doy=doy+1 !leap year and month is larger than february



write(GPS_date,'(i2.2,a1,i3.3,a1,i5.5)')yr,':',doy,':',hr*3600+min*60+sec




end function GPS_date

!function which returns sinex date as dd-mm-yyyy string
function norm_date(sindat)
implicit none
character(*),intent(in)::sindat
character(12)::norm_date
integer::mndy(12)
integer::i,mnth,col1,col2,dy,yr,sec
data mndy/31,28,31,30,31,30,31,31,30,31,30,31/

!get the indices of the double colons in the string
col1=index(sindat,':')
col2=index(sindat(col1+1:),':')+col1
!write(*,*)dtstr,col1,col2
!read in year day and seconds
read(sindat(col1-2:col1-1),'(I2)')yr
read(sindat(col1+1:col2-1),'(I3)')dy
read(sindat(col2+1:),*)sec


if(mod(yr,4) .eq. 0)mndy(2)=29 !leap year

if(yr<50)then
   yr=yr+2000
else
   yr=yr+1900
end if
!now find month and day of month
do,i=1,12
   if(dy<=mndy(i))then
      
      mnth=i
      exit
   end if
   dy=dy-mndy(i)
end do

!add one day to create right date
!if(sec >=0)dy=dy+1


!construct date string

write(norm_date,'(i2.2,a1,i2.2,a1,i4.4)')dy,'-',mnth,'-',yr


end function norm_date



function date_2_jd(dy,mn,yr,hr,min,sec)
implicit none
double precision::date_2_jd
integer::yr,dy,mn,hr,min,sec
integer::c
integer::ind

ind = dy-32075+1461*(yr+4800+(mn-14)/12)/4+367*&
     (mn-2-(mn-14)/12*12)/12-3*((yr+4900+(mn-14)/12)/100)/4

date_2_jd=ind+(hr-12)/24.d0+min/1440.d0+dble(sec)/86400.d0

end function date_2_jd


subroutine jd_2_date(jd,dy,mn,yr,sec)
  implicit none
  double precision,intent(in)::jd
  integer,intent(out)::mn,yr,dy,sec
  integer::p,s1,s2,s3,s4,n,i,e,q

  
  p=int(jd+0.5d0)
  s1=p+68569
  n=(4*s1)/146097
  s2=s1-(146097*n+3)/4
  i=(4000*(s2 + 1))/1461001
  s3=s2 - (1461*i)/4 + 31
  q=(80*s3)/2447
  e = s3 - (2447*q)/80
  s4 = q/11
  mn = q + 2 - 12*s4 
  yr = 100*(n - 49) + i + s4 
  dy= int(e + jd - p + 0.5) 
  sec=((e + jd - p + 0.5)-dy)*86400
  
end subroutine jd_2_date

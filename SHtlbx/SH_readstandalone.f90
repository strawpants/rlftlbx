!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Subroutines which tries to read potential coefficients from SH files
!!It tries to determine whether the file is in a known format (can be forced)
!!A skip can also be optionally provided and when really knowing the file by heart one can specify a fortran format string
!!The length of the input vectors determine the maximum degree which is read in.
!!Coded by Roelof Rietbroek, Tue May 25 14:20:55 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue May 29 14:17:53 2007

subroutine SH_readgrav(filen,clm,slm,clm_sig,slm_sig,skip,frmt,type,str)
  use SHtlbx
  implicit none
  double precision,dimension(:),optional,intent(out):: clm,slm,clm_sig,slm_sig
  integer, parameter::nvar =6
  integer, intent(in),optional::skip,type
  character(*),optional,intent(in)::frmt,str
  character(*), intent(in),optional::filen
  character(200)::dum,frm
  character(20)::srch
  integer:: typ,i,j,skp,uni,lmax,st,pos,l,m,ncol,last,srchnd,err,stderr
  double precision,allocatable,dimension(:)::dumload
  logical::stdin,fex
  stdin=.false.
  stderr=0
  lmax=0
  skp=0 !default no skip of header lines
  frm=''!initialize to empty format
  typ=0
 !default unit for open file
     uni=13
 
 !check if filename is present
  if(.not. present(filen))then
    stdin=.true. !read from standard input
    uni=5
  end if

  ! check filetype
  if(present(type)) then
     typ=type
  else
!     write(*,*)'test',trim(filen),type,present(type)
    if(stdin)then
       call SH_readmeta(type=typ)
    else 
       call SH_readmeta(filen=trim(filen),type=typ)
    end if
  end if

!default type specifix parameters (will be changed later when custom parameters are inputted in the routine)
  select case(typ)
  case(1)!GSM,GFZ JPL format
        srch='GRCOF2' !default search string
        srchnd=6 !end of default search string
        !start the start index to read from
        st=srchnd+1!last non blank character
  case(2)!GINS format
        srch='   ' !default search string three spaces
        srchnd=3
        !start the start index to read from
        st=srchnd+1!last non blank character
     !explicitly specify format for GINS type(to guarantee the read in case the GINS file contain DOT,S1A,etc entries) 
        !(2i3,a3,2e21.14,2e13.6,1x,i2)
        frm='(2i3,3x,2e21.14,1x,2e13.6)'
        skp=6
  case(3)!ICGEM format
     	srch='gfc'!default icgem search string
     	srchnd=3 !Include a trailing space to prevent ambiguity with gfct!!
      !start the start index to read from
      st=srchnd+2!last non blank character
  case(4)!standard interchange format
     srch=''
     srchnd=0
     st=srchnd+1
     if(stdin) then
        skp=0
     else
        skp=1
     end if
  case(5)!GRACE calibrated error format
     srch='CALSDV  '
     srchnd=8
     st=srchnd+1
     skp=0
  case default!generic
       srch=''
       srchnd=0
       st=srchnd+1
       skp=0
    end select

!!!!!optional parameters areincorporated below

  !set format if requested
  if(present(frmt))then
     frm=trim(frmt)
  end if

  !specify search string if required
  if (present(str)) then ! if the search string 
     srch=trim(str)
     srchnd=len(str)!last entry of string (including trailing spaces)
     st=srchnd+1
  end if

!!!set up variable specific items

  !get lmax requested and the maximum column (and initioalize vectors)

  if(present(clm))then
     call SH_lm(size(clm,1),lmax,m)
     ncol=1
     clm=0.d0
  end if

  if (present(slm))then
     call SH_lm(size(slm,1),lmax,m)
     ncol=2
     slm=0.d0
  end if

  if (present(clm_sig))then
     call SH_lm(size(clm_sig,1),lmax,m)
     ncol=3
     clm_sig=0.d0
  end if

  if(present(slm_sig))then
     call SH_lm(size(slm_sig,1),lmax,m)
     ncol=4
     slm_sig=0.d0
  end if


  !open file if not reading from standard input
 if (.not.stdin) then
    !inquire filename
    inquire(FILE=trim(filen),EXIST=fex)
    if( .not. fex)then
       write(stderr,*)'File ',trim(filen),' does not exist, quitting'
       stop
    end if
    open(unit=uni,file=filen,status='old')
 end if
 !skip header if requested
  if(present(skip)) then
     skp=skip
  end if

  do,i=1,skp
     read(uni,'(A200)')dum
!     write(*,*)dum
  end do


  

  !now get the data line by line as a string
  allocate(dumload(ncol))
  do,i=1,100000 !many lines for a maximum (stops when reaching end of file)
     read(unit=uni,fmt='(A200)',iostat=last)dum
     !write(*,*)dum
     if (last .ne. 0 .or. (Len_Trim(dum).eq. 0 .and. i >100)) exit !end of file encountered or empty line read (sometimes end of file shifted))
!      write(*,*)i,trim(dum)
     !cycle loop when search string is not found
     !write(*,*)index(trim(dum),srch(1:srchnd)),trim(dum)

     if (index(trim(dum),srch(1:srchnd)) .eq. 0)cycle !search for a string with a
     dumload=0.d0
     !write(*,*)dum,st
     if(trim(frm) .eq. '') then
       read(dum(st:),*,iostat=err) l,m,dumload
     !write(*,*)l,m,dumload
       if(err .ne. 0)then
         
         read(dum(st:),*) l,m,dumload(1)!try reading only one value (in case zonal coefficients are left as blank)
         !write(*,*)l,m,dumload
      end if


     else
        
       read(dum,trim(frm)) l,m,dumload
     
    end if
    !if (m .eq. 1)    write(*,*)trim(dum),dumload
     if(l>lmax)cycle !cycle when degree requested is smaller than the degree read (cycle accelerates the reading)
      !    write(*,*)st,trim(dum(st:))
     ! write(*,*)'test'

     !now get the values and put them in the rights vector

     pos=SH_pos(l,m)
     if(present(clm))clm(pos)=dumload(1)
     if(present(slm))slm(pos)=dumload(2)
     if(present(clm_sig))clm_sig(pos)=dumload(3)
     if(present(slm_sig))slm_sig(pos)=dumload(4)


    ! write(*,*)l,m,clm(pos),slm(pos),clm_sig(pos),slm_sig(pos)
  end do

  if(.not. stdin)close(uni)

end subroutine SH_readgrav


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!function which tries to check what kind of file is inputted (must be fullfilename)
!!outputs 0: generic (not recognized)
!!        1: GRACE (GFZ,JPL,CSR,GRGS)
!!        2  GRGS (GINS)
!!	  3: ICGEM format
!!coded by Roelof Rietbroek 23-5-2007

!!$function SH_ftype(filen)
!!$implicit none
!!$character(*),intent(in)::filen
!!$character(200)::dum
!!$integer::SH_ftype,i,un
!!$
!!$SH_ftype=0
!!$!open file and reads first line
!!$!write(*,*)filen
!!$un=13
!!$open(unit=un,file = filen)
!!$read(un,'(A200)')dum
!!$
!!$
!!$if (index(dum,'FIRST') .ne. 0)then 
!!$   SH_ftype=1 !conventional GRACE type
!!$else if ((index(dum,'FIELD - EIGEN').ne. 0) )then
!!$   SH_ftype=2 !GINS GRGS type
!!$else !now try a somewhat more rigourous approach to see if the format is ICGEM
!!$   do,i=1,80 !check a maximum  of 80 header lines for a key word
!!$      read(un,'(A200)')dum
!!$!write(*,*)trim(dum)
!!$      if (index(dum,'product_type') .ne. 0) then
!!$         SH_ftype=3
!!$         exit
!!$      end if
!!$   end do
!!$end if
!!$
!!$close(13)
!!$
!!$
!!$end function SH_ftype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!subroutine to extract metadata from SH file header (or standard input)
!!extracts: tcent(if possible) center time of observation (in year.fractionofyear)
!!file type:
!!GRACE(GFZ,CSR,JPL):type=1
!!GINS(french from GRGS):type=2
!!ICGEM: type=3
!!type=4 suitable for piping (first line contains lmax and times)
!!generic: type=0
!!and extracts maximum degree,lmax, of the file (if supported)
!!all arguments are optional(except the filename
subroutine SH_readmeta(filen,type,lmax,mu,rad,tcent,tstart,tend)
use SHtlbx
implicit none
integer,optional,intent(out)::type,lmax
double precision,optional,intent(out)::mu,rad,tcent,tstart,tend!strat en end of obsrecation period
double precision::dum1,dum2,dum3,dum4,t0,t1,t2
character(*),intent(in),optional::filen
character(200)::dum
integer::uni,i,ftyp,lm,stderr
logical::stdin,fex
stderr=0

stdin=.false.
uni=13
!defaults

ftyp=0 !default generic type
!open file
!write(*,*)trim(filen),Len_Trim(filen)

!check if filen is present otherwise read from standard input
if(present(filen))then
!write(*,*)trim(filen)
!inquire filename
   inquire(FILE=trim(filen),EXIST=fex)
   if( .not. fex)then
      write(stderr,*)'File ',trim(filen),' does not exist, quitting'
      stop
   end if
   open(unit=uni,file=trim(filen))
else
   uni=5
   stdin=.true.
end if
!first check which  file type it is

!read first line
!write(*,*)dum

read(uni,'(A200)')dum
if (index(dum,'FIRST') .ne. 0)then 
  ftyp=1 !conventional GRACE type
  
else if ((index(dum(1:6),'FIELD').ne. 0) )then
   ftyp=2 !GINS GRGS type
else if (index(dum,'META') .ne.0)then
ftyp=4
else !now try a somewhat more rigourous approach to see if the format is ICGEM
   do,i=1,80 !check a maximum  of 80 header lines for a key word
      read(uni,'(A200)')dum
      if (index(dum,'product_type') .ne. 0) then
         ftyp=3
         exit
      else if(index(dum,'CALSDV') .ne. 0)then !GRACE calibrated error file
         ftyp=5
         exit
      end if
   end do
end if

!set typ to ftyp if required
if(present(type))type=ftyp

!now rewind file if other parameters are requested
if(present(lmax) .or. present(mu) .or. present(rad) .or. present(tcent) .or. present(tstart) .or. present(tend)) then
  

      if(stdin)then
         if(ftyp .ne. 4)then
            write(*,*)'cannot rewind this standard input,try the default output-type(4) of SH_dump'
            stop
         end if
      else
         rewind(uni)
      end if
  


else !return to subroutine caller
 if(.not. stdin) close(uni)
   return
end if


!run through the header again but now extract more required stuff
select case(ftyp)

case(1)!GRACE GFZ,CSR,JPL type
read(uni,'(A200)')dum
!get observation times
if(present(tcent) .or. present(tstart) .or. present(tend)) then
call getctime(dum,ftyp,t0,t1,t2)
if(present(tstart))tstart=t0
if(present(tcent))tcent=t1
if(present(tend))tend=t2
end if

do,i=1,20
   read(uni,'(A200)')dum
   if(index(dum,'EARTH') .ne. 0)then
      read(dum(6:),*)dum1,dum2
      if(present(mu))mu=dum1
      if(present(rad))rad=dum2
   else if((index(dum,'SHM ').ne. 0).and.(present(lmax)))then
      read(dum(4:),*)lmax !read lmax supported by file
   end if
end do

case(2)!GINS typ
   read(uni,'(A200)')dum
   if(present(tcent) .or. present(tstart) .or. present(tend)) then
      call getctime(dum,ftyp,t0,t1,t2)
      if(present(tstart))tstart=t0
      if(present(tcent))tcent=t1
      if(present(tend))tend=t2
   end if
   
   read(uni,'(/,4e20.14)')dum1,dum2,dum3,dum4
   if(present(mu))mu=dum3
   if(present(rad))rad=dum1
   read(uni,'(A200)')dum
   read(uni,'(A200)')dum
   if((index(dum,'MAXIMAL DEGREE').ne. 0).and.(present(lmax)))read(dum(17:),*)lmax !read lmax supported by file

case(3)!ICGEM type
do,i=1,80
   read(uni,'(A200)')dum
   if(present(lmax) .and. (index(dum,'max_degree') .ne. 0))read(dum(11:),*)lmax
   if(present(mu) .and. (index(dum,'earth_gravity_constant') .ne. 0))read(dum(23:),*)mu
   if(present(rad) .and. (index(dum,'radius') .ne. 0))read(dum(7:),*)rad
 !Get observation times 
 if(index(dum,'time_period_of_data') .ne. 0) then

    if(present(tcent) .or. present(tstart) .or. present(tend)) then
      call getctime(dum,ftyp,t0,t1,t2)
      if(present(tstart))tstart=t0
      if(present(tcent))tcent=t1
      if(present(tend))tend=t2
   end if

   end if
   if((index(dum,'end_of_head') .ne. 0))exit !exit loop when last header line is reached
end do

case(4)!!suitable for piping uses the previous dum doesn't call getctime
   read(dum(6:),*)lm,t0,t1,t2
   if(present(lmax))lmax=lm
   if(present(tstart))tstart=t0
   if(present(tcent))tcent=t1
   if(present(tend))tend=t2
   if(present(mu))mu=GM
   if(present(rad))rad=RE

case(5)!GRACE calibrated error file (contains as much info as a generic file)

   if(present(lmax)) lmax=-1!give impossible value
   if(present(mu))mu=GM!default Earths gravitational constant
   if(present(rad))rad=RE!default mean Earth radius
   if(present(tstart))tstart=0.d0
   if(present(tcent))tcent=0.d0
   if(present(tend))tend=0.d0

case default !generic case assumes no header and thus no metadata

   if(present(lmax)) lmax=-1!give impossible value
   if(present(mu))mu=GM!default Earths gravitational constant
   if(present(rad))rad=RE!default mean Earth radius
   if(present(tstart))tstart=0.d0
   if(present(tcent))tcent=0.d0
   if(present(tend))tend=0.d0

end select

if(.not. stdin)close(uni)! close unit
end subroutine SH_readmeta



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Function which gets the center times for specific strings belonging to 
!!ICGEM,GRACE ort GINS formatted SH files
!!center time is in years
subroutine getctime(string,type,tst,tc,tnd)
implicit none
character(*),intent(in)::string
integer,intent(in)::type
double precision, intent(out)::tst,tc,tnd
double precision::year1,year2,day1,day2,yeardays,yeardaysold
integer::i,leap,yr1,mnth1,dy1,err,mid,jd,yr2,mnth2,dy2
logical::first

!(re)initialize centertime
tst=0.d0
tnd=0.d0
tc=0.d0
select case(type)
case(1)!GRACE GFZ,JPL,CSR format
   if ((index(string,'CSR') .ne. 0) .and. (index(string,'0000_0001') .ne. 0)) then
            Read(string(18:21),*) Year1
            Read(string(26:29),*) Year2
            Read(string(22:24),*) Day1
            Read(string(30:32),*) Day2
            else !If Release 3 or higher
            Read(string(13:16),*) Year1
            Read(string(21:24),*) Year2
            Read(string(17:19),*) Day1
            Read(string(25:27),*) Day2
    end If

    !construct start time
    if(mod(int(year1),4) .eq. 0)then
                leap=366
    else
                leap=365
    end if

    tst=year1+day1/leap

    !construct end time
    if(mod(int(year2),4) .eq. 0)then
                leap=366
    else
                leap=365
    end if

    tnd=year2+day2/leap

    !construct centertime
    tc=(tst+tnd)/2.d0
    
case(2)!GINS format
   !search for the day strings (days since 1-1-1950)
	!start searching for a "_" from the back
   mid=index(string,'_',.true.)
  	read(string(mid-5:mid-1),'(F5.0)',iostat=err)day1
	read(string(mid+1:mid+5),'(F5.0)',iostat=err)day2
 if(err .ne. 0)return!return to calling function
! write(*,*)day1,day2
 !get time points
 tst=day1
 tnd=day2
 
!now get the reference years
 yeardays=10957.d0 !days from 1-1-1950 for the the year 1980
 yeardaysold=yeardays
first=.true.
 do,i=0,99 !hundred years from 1980 will probably do
    if(mod(i,4).eq. 0)then !account for leap years
       leap=366
    else
       leap=365
    end if
    yeardays=yeardays+leap
!write(*,*)leap,yeardays
    if(yeardays>= tst .and. first ) then
        !now construct reference date
       tst=dble(1980+i)+(tst-yeardaysold)/dble(leap)
       first=.false.
    end if
    
    if (yeardays >= tnd ) then
       tnd=dble(1980+i)+(tnd-yeardaysold)/dble(leap)
       !exit loop when end time has been set

       exit
    else
       yeardaysold=yeardays
    end if
 end do
 !finally construct the reference date
 tc=(tnd+tst)/2.d0

case(3)!ICGEM format
   mid=index(string,' - ')
   read(string(mid-8:),'(I4,I2,I2,3x,I4,I2,I2)')yr1,mnth1,dy1,yr2,mnth2,dy2
!   write(*,*)yr1,mnth1,dy1,yr2,mnth2,dy2
 !get start time
   if(mod(yr1,4).eq. 0)then !account for leap years
       leap=366
    else
       leap=365
    end if
   ! write(*,*)yr,mnth,dy,leap,(jd(yr,mnth,dy)-jd(yr,1,1))/dble(leap)
    tst=dble(yr1)+(jd(yr1,mnth1,dy1)-jd(yr1,1,1))/dble(leap)
 !get end time
    if(mod(yr2,4).eq. 0)then !account for leap years
       leap=366
    else
       leap=365
    end if
   ! write(*,*)yr,mnth,dy,leap,(jd(yr,mnth,dy)-jd(yr,1,1))/dble(leap)
    tnd=dble(yr2)+(jd(yr2,mnth2,dy2)-jd(yr2,1,1))/dble(leap)
     
!get centertime
    tc=(tst+tnd)/2.d0
case default!generic format (set to zero)
tst=0.d0
tc=0.d0
tnd=0.d0
end select

end subroutine getctime

!yields the julian day for a given gregorian date
    FUNCTION JD (YEAR,MONTH,DAY)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).
!
      INTEGER,intent(in)::YEAR,MONTH,DAY
      INTEGER::I,J,K,JD
!
      I= YEAR
      J= MONTH
      K= DAY
!
      JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4
!
      RETURN
      END FUNCTION JD


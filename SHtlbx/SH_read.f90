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

!!Updated by Roelof Rietbroek, Fri Jan 25 13:52:26 2008
!!added comment serach string to explicitly ignore comment lines
!!Updated by Roelof Rietbroek, Fri Jun 27 23:52:21 2008
!!fixed bug found by Indridi Einarsson: ICGEM format was not recoginzed when "product_type" keyword occurred in the first line
!!also changed the do i=1,100000 to a do while loop such that large expansion can be read

!!Updated by Roelof Rietbroek, Wed Aug 12 15:33:00 2009
!!fixed a bug in the linear format

!!Updated by Roelof Rietbroek, Wed Sep 30 10:31:52 2009
!!skip emtpy lines
!!Updated by Roelof Rietbroek, Tue Nov  2 09:51:49 2010
!! make linear format string consistent with SH_write
!!Updated by Roelof Rietbroek, Thu Feb 24 13:38:09 2011
!!added: maxdum which fixes a bug which occurred in gfortran with linear format (format inconsistent with number of doubles)

!!Updated by Roelof Rietbroek, Tue Apr 26 12:37:54 2011
!!allow grace files to be read from standard input

!!Updated by Roelof Rietbroek, Mon May 16 10:14:39 2011
!!make a better check for the SHM flag ( only check first 5 characters of a line)

!!Updated by Roelof Rietbroek, Mon Apr 16 09:54:36 2012
!! allow icgem timevariable components to be read

!!Updated by Roelof Rietbroek, Fri Apr 20 14:23:45 2012
!! fixed bug ( unitialized nsrch for the linear format


!!Updated by Roelof Rietbroek, Sun Mar 10 12:02:16 2013
!! made linear format consistent with BIN_swiss again

!!Updated by Roelof Rietbroek, Thu Dec  3 09:14:08 2015
!!modified SH_readmeta such that it also allows spaces before the keywords in the icgem header


subroutine SH_readgrav(filen,clm,slm,clm_sig,slm_sig,skip,frmt,type,str,clm_tv,slm_tv)
  use SHtlbx,only:SH_readmeta,SH_pos
  implicit none
  double precision,dimension(:),optional,intent(out):: clm,slm,clm_sig,slm_sig
  integer, intent(in),optional::skip,type
  double precision,optional,intent(out),dimension(:,:)::clm_tv,slm_tv! optional time variable components harmonics (trend, yearly semi yearly)
  character(*),optional,intent(in)::frmt,str
  character(*), intent(in),optional::filen
  character(200)::dum,frm
  character(20)::srch,cmtsrch,srch2
  character(20)::srchs(10) ! array of search strings with special meaning (format dependent)
  integer::nsrch,frst
  integer:: typ,i,j,skp,uni,lmax,st,pos,l,m,last,srchnd,err,stderr,cmtsrchnd,lmaxtv
  integer::maxdum
  integer::isrch
  double precision::dumload(6)

  logical::stdin,fex
  cmtsrch='vxjhcne' !initialize comment search string to bullocks
  last=0
  cmtsrchnd=7
  stdin=.false.
  stderr=0
  lmax=0
  skp=0 !default no skip of header lines
  frm=''!initialize to empty format
  typ=0

  maxdum=4 !maximum amount of doubles readable
  nsrch=0
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

!check for requested dot terms and see if it matches the type
  if(present(clm_tv) .or. present(slm_tv))then
     select case(typ)
     case(1,3)!format supports time variable coefficients
        lmaxtv=0
        l=0
        !get maximum degree for lmaxtv
        if(present(clm_tv))then
           call SH_lm(size(clm_tv,1),lmaxtv,m)
           if(size(clm_tv,2) < 5)then
              write(stderr,*)"ERROR: SH_readgrav, input array for timevariable coefficients too small"
              stop
           end if
        end if
        if(lmaxtv ==0)then
           write(stderr,*)"ERROR: SH_readgrav, Timevariable coefficients requested but lmax=0? "
           stop
        end if

     case default!not supported
        write(stderr,*)"Simultaneous reading of Timevariable fields are not supported for this file type"
        stop
     end select
  end if


!default type specifix parameters (will be changed later when custom parameters are inputted in the routine)
  select case(typ)
  case(1)!GSM,GFZ JPL format
     srchs(1)='GRCOF2' !default search string
     srchs(2)='GRDOTA' ! search string for DOT terms
     nsrch=2 
     cmtsrch='CMMNT'
     cmtsrchnd=5
     srchnd=6 !end of default search string
     !start the start index to read from
     st=srchnd+1!last non blank character

     frst=10 ! search only the first 10 characters of a line
  case(2)!GINS format
        srchs(1)='   ' !default search string three spaces
        srchnd=3
        nsrch=1 
        !start the start index to read from
        st=srchnd+1!last non blank character
     !explicitly specify format for GINS type(to guarantee the read in case the GINS file contain DOT,S1A,etc entries) 
        !(2i3,a3,2e21.14,2e13.6,1x,i2)
        frm='(2i3,3x,2e21.14,1x,2e13.6)'
        skp=6
        frst=10 ! search only the first 10 characters of a line
  case(3)!ICGEM format
     srchs(1)='gfc '!default icgem search string ( with trailing space!)
     srchs(2)='gfct' !bias
     srchs(3)='trnd' !trend
     srchs(4)='dot'  ! also a trend (conventional)
     srchs(5)='asin' !sin harmonic amplitude
     srchs(6)='acos' !cos harmonic amplitude
     nsrch=6
     srchnd=4
     !start the start index to read from
     st=srchnd+1!last non blank character
     frst=5 ! search only the first 5 characters of a line
     maxdum=5 ! also allows for phases to be read in
  case(4)!standard interchange format
     srchs(1)=''
     nsrch=1
     srchnd=0
     st=srchnd+1
     if(stdin) then
        skp=0
     else
        skp=1
     end if
     frst=1 ! search only the first 1 characters of a line
  case(5)!GRACE calibrated error format
     srchs(1)='CALSDV  '
     nsrch=1
     srchnd=8
     st=srchnd+1
     skp=0
     frst=10 ! search only the first 10 characters of a line
  case(6,7) !linear ascii format 24 character name then columns 1 coefficient per row
     frst=24 ! search only the first 24 characters of a line
     nsrch=2
     srchs(1)='CN '
     srchs(2)='SN '
     srchnd=3
     st=srchnd+1
     skp=0
     maxdum=2
!      frm='(4x,i3,i3,15x,G18.10,1x,G18.10)'
     frm='(4x,i3,i3,15x,G22.14,1x,G22.14)' ! updated 10-03-2013
  case default!generic
       srchs(1)=''
       nsrch=1
       srchnd=0
       st=srchnd+1
       skp=0
       frst=1 ! search only the first 1 characters of a line
    end select

!!!!!optional parameters areincorporated below

  !set format if requested
  if(present(frmt))then
     frm=trim(frmt)
  end if

  !specify search string if required
  if (present(str)) then ! if the search string 
     srchs(1)=trim(str)
     nsrch=1
     srchnd=len(str)!last entry of string (including trailing spaces)
     st=srchnd+1
  end if

!!!set up variable specific items

  !get lmax requested and the maximum column (and initioalize vectors)

  if(present(clm))then
     call SH_lm(size(clm,1),lmax,m)
     clm=0.d0
  end if

  if (present(slm))then
     call SH_lm(size(slm,1),lmax,m)
     slm=0.d0
  end if

  if (present(clm_sig))then
     call SH_lm(size(clm_sig,1),lmax,m)
     clm_sig=0.d0
  end if

  if(present(slm_sig))then
     call SH_lm(size(slm_sig,1),lmax,m)
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
     read(uni,'(A)')dum
!     write(*,*)dum
  end do


  

  !now get the data line by line as a string
  last=0
  do while(last .eq. 0) !untill the end of file is reached (or I/O error)
     read(unit=uni,fmt='(A)',iostat=last)dum
     if (last .ne. 0 )exit !end of file encountered or empty line read (sometimes end of file shifted))
     if(dum .eq. '')cycle ! skip empty lines
     !cycle loop when search string is not found
     !write(*,*)index(trim(dum),srch(1:srchnd)),trim(dum)
     if(index(trim(dum),cmtsrch(1:cmtsrchnd)) > 0)cycle !cycle when a comment string was found

     !only proceed the loop when the line contains a valid search string
     isrch=0
     do,i=1,nsrch
        if(index(dum(1:frst),srchs(i)(1:srchnd)) .ne. 0)then
           isrch=i
           exit
        end if
     end do
     if(isrch == 0)cycle ! no valid line was found

     dumload=0.d0
     
     if(trim(frm) .eq. '') then
        read(dum(st:),*,iostat=err) l,m,dumload(1:maxdum)
!        read(dum(st:),*) l,m,dumload
        !write(*,*)l,m,dumload
        !         if(err .ne. 0)then
        
        !            read(dum(st:),*,iostat=err) l,m,dumload(1)!try reading only one value (in case zonal coefficients are left as blank)
        if(err > 0)then ! if the error is smaller than zeroa end of record was detected ( no problem values will be zero)
           write(stderr,*)'WARNING: formaterror',err,', skipping line:',trim(dum)
           cycle
        end if
        !         end if
        
     else ! with format

        read(dum,trim(frm),iostat=err) l,m,dumload(1:maxdum)
       !read(dum,trim(frm)) l,m,dumload
!        if(l>=99)write(0,*)'test',err,trim(dum),l,m,dumload(1:2)         
        if(err > 0)then
           write(stderr,*)'SH_read WARNING: format error',err,', skipping line:',trim(dum)
           cycle
        end if
     end if

      !now get the values and put them in the right vector
      if(l>lmax)cycle !cycle when degree requested is smaller than the degree read (cycle accelerates the reading)
         

      pos=SH_pos(l,m)
      select case(typ)
      case(1)!GRCOF2 format
         select case(isrch) 
         case(1) !Cosine Sine coefficients
            if(present(clm))clm(pos)=dumload(1)
            if(present(slm))slm(pos)=dumload(2)
            if(present(clm_sig))clm_sig(pos)=dumload(3)
            if(present(slm_sig))slm_sig(pos)=dumload(4)
         case(2) !Conventional rates
            if(l>lmaxtv)cycle
            i=1
            if(present(clm_tv))clm_tv(pos,i)=dumload(1)
            if(present(slm_tv))slm_tv(pos,i)=dumload(2)
         end select
      case(3) !ICGEM format
         select case(isrch) 
         case(1,2) !Cosine Sine coefficients (gfc or gfct)
            if(present(clm))clm(pos)=dumload(1)
            if(present(slm))slm(pos)=dumload(2)
            if(present(clm_sig))clm_sig(pos)=dumload(3)
            if(present(slm_sig))slm_sig(pos)=dumload(4)
         case(3,4) !trends
            if(l>lmaxtv)cycle
            i=1
            if(present(clm_tv))clm_tv(pos,i)=dumload(1)
            if(present(slm_tv))slm_tv(pos,i)=dumload(2)
         case(5)!sin harmonic
            if(l>lmaxtv)cycle
            !write(0,*)dumload
            if(dumload(5) .eq. 1.d0)then
               i=2
            else
               i=4
            end if
            if(present(clm_tv))clm_tv(pos,i)=dumload(1)
            if(present(slm_tv))slm_tv(pos,i)=dumload(2)

         case(6)!cos harmonic
            if(l>lmaxtv)cycle
            if(dumload(5) .eq. 1.d0)then
               i=3
            else
               i=5
            end if
            if(present(clm_tv))clm_tv(pos,i)=dumload(1)
            if(present(slm_tv))slm_tv(pos,i)=dumload(2)
         end select
      case(6,7)!linear format
         select case(isrch) 
         case(1) !Cosine coefficients
            if(present(clm))clm(pos)=dumload(1)
            if(present(clm_sig))clm_sig(pos)=dumload(2)
         case(2)
            if(present(slm))slm(pos)=dumload(1)
            if(present(slm_sig))slm_sig(pos)=dumload(2)
         end select
      case default ! allo other formats
         if(present(clm))clm(pos)=dumload(1)
         if(present(slm))slm(pos)=dumload(2)
         if(present(clm_sig))clm_sig(pos)=dumload(3)
         if(present(slm_sig))slm_sig(pos)=dumload(4)
      end select
      
      ! select case(isrch) 
      ! case(1) !default method
      !    if(separ)then
      !       if(present(clm))clm(pos)=dumload(1)
      !       if(present(clm_sig))clm_sig(pos)=dumload(2)
      !    else
      !       if(present(clm))clm(pos)=dumload(1)
      !       if(present(slm))slm(pos)=dumload(2)
      !       if(present(clm_sig))clm_sig(pos)=dumload(3)
      !       if(present(slm_sig))slm_sig(pos)=dumload(4)
      !    end if

      ! case(2) 
      !    if(separ)then
      !       if(present(slm))slm(pos)=dumload(1)
      !       if(present(slm_sig))slm_sig(pos)=dumload(2)
      !    else if(dot)then
      !       if(present(clm))clm(pos)=dumload(1)
      !       if(present(slm))slm(pos)=dumload(2)
      !       if(present(clm_sig))clm_sig(pos)=dumload(3)
      !       if(present(slm_sig))slm_sig(pos)=dumload(4)
      !    end if
      ! end select
      ! if(separ)then
      !    if(index(trim(dum),srch(1:srchnd)) .ne. 0)then
      !       if(present(clm))clm(pos)=dumload(1)
      !       if(present(clm_sig))clm_sig(pos)=dumload(2)
      !    else 
      !       if(present(slm))slm(pos)=dumload(1)
      !       if(present(slm_sig))slm_sig(pos)=dumload(2)
      !    end if
      ! else if (dot)then
      !    if(index(trim(dum),srch(1:srchnd)) .ne. 0)then
      !       if(present(clm))clm(pos)=dumload(1)
      !       if(present(slm))slm(pos)=dumload(2)
      !       if(present(clm_sig))clm_sig(pos)=dumload(3)
      !       if(present(slm_sig))slm_sig(pos)=dumload(4)
      !    else! read dot terms
      !       if(l>lmaxtv)cycle
      !       if(present(clm_dot))clm_dot(pos)=dumload(1)
      !       if(present(slm_dot))slm_dot(pos)=dumload(2)
      !    end if
      ! else
      !    if(present(clm))clm(pos)=dumload(1)
      !    if(present(slm))slm(pos)=dumload(2)
      !    if(present(clm_sig))clm_sig(pos)=dumload(3)
      !    if(present(slm_sig))slm_sig(pos)=dumload(4)
      ! end if
    !if (m .eq. 1)    write(*,*)trim(dum),dumload
    
      !    write(*,*)st,trim(dum(st:))
     ! write(*,*)'test'



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
!!Updated by Roelof Rietbroek, Tue Nov 27 17:47:45 2007
!! added end of file detection (previously gave end of file runtime errors for shortfiles (< 80 lines))
!!Updated by Roelof Rietbroek, Thu Apr 14 16:58:33 2011
!! added two possibilities to determine linear type

subroutine SH_readmeta(filen,type,lmax,mu,rad,tcent,tstart,tend)
use SHtlbx,only:getctime,GM,RE
implicit none
integer,optional,intent(out)::type,lmax
double precision,optional,intent(out)::mu,rad,tcent,tstart,tend!strat en end of obsrecation period
double precision::dum1,dum2,dum3,dum4,t0,t1,t2
character(*),intent(in),optional::filen
character(200)::dum
integer::uni,i,ftyp,lm,stderr,last,indd,ind
logical::stdin,fex
integer::nh
nh=200 !default amount of headers to be read
stderr=0

stdin=.false.
uni=13
!defaults

ftyp=0 !default generic type
!open file
!write(*,*)trim(filen),Len_Trim(filen)y

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

read(uni,'(A)')dum
if (index(dum,'FIRST') .ne. 0)then 
  ftyp=1 !conventional GRACE type
  
else if ((index(dum(1:6),'FIELD').ne. 0) )then
   ftyp=2 !GINS GRGS type
else if (index(dum,'META') .ne.0)then
ftyp=4
else !now try a somewhat more rigourous approach to see if the format is ICGEM
   if(.not. stdin)rewind(uni)!bug fix
   do,i=1,nh !check a maximum  of nh header lines for a key word
      read(uni,'(A200)',iostat=last)dum
      if(last .ne. 0)exit !exit loop when end of file is reached (generic type)
      if (index(dum,'product_type') .ne. 0) then
         ftyp=3
         exit
      else if(index(dum,'CALSDV') .ne. 0)then !GRACE calibrated error file
         ftyp=5
         exit
      else if( (index(dum,'BINV') .ne. 0) )then !linear SH file
         ftyp=6
         exit
      else if((index(dum(2:4),'CN ') .ne. 0) .or. (index(dum(2:4),'SN ') .ne. 0))then
         ftyp=7 !linear type but without header
      end if
   end do
end if

!set typ to ftyp if required
if(present(type))type=ftyp

!now rewind file if other parameters are requested
if(present(lmax) .or. present(mu) .or. present(rad) .or. present(tcent) .or. present(tstart) .or. present(tend)) then

      if(stdin)then
         if(ftyp .ne. 4 .and. ftyp .ne. 1 .and. ftyp .ne. 6)then
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
!read(uni,'(A)')dum
!get observation times
if(present(tcent) .or. present(tstart) .or. present(tend)) then
call getctime(dum,ftyp,t0,t1,t2)
if(present(tstart))tstart=t0
if(present(tcent))tcent=t1
if(present(tend))tend=t2
end if

do,i=1,20
   read(uni,'(A)')dum
   if(index(dum,'EARTH') .ne. 0)then
      read(dum(6:),*)dum1,dum2
      if(present(mu))mu=dum1
      if(present(rad))rad=dum2
   else if(index(dum(1:5),'SHM ').ne. 0)then
      if(present(lmax))read(dum(4:),*)lmax !read lmax supported by file
      exit !exit loop, no more meta info is following
   end if
end do

case(2)!GINS typ
   read(uni,'(A)')dum
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
do,i=1,300 !check the first 300 lines
   read(uni,'(A)')dum
   ind=index(dum,'max_degree')
   if(present(lmax) .and. ( ind .ne. 0))read(dum(ind+11:),*)lmax
   ind=index(dum,'earth_gravity_constant')
   if(present(mu) .and. (ind .ne. 0))read(dum(ind+23:),*)mu
   ind=index(dum,'radius')
   if(present(rad) .and. (ind .ne. 0))read(dum(ind+7:),*)rad
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
case(6) !linear SH type with header
   do,i=1,nh !check the first nh lines for META data
      read(uni,'(A)')dum
      indd=index(dum,'Lmax ')
      if(indd .ne. 0 .and. present(lmax))read(dum(indd+5:),*)lmax
!       indd=index(dum,'Lmin ')
!       if(indd .ne. 0 .and. present(lmin))read(dum(indd:),*)lmin

      indd=index(dum,'CTime ')
      if(indd .ne. 0 .and. present(tcent))read(dum(indd+6:),*)tcent
      indd=index(dum,'ETime ')
      if(indd .ne. 0 .and. present(tend))read(dum(indd+6:),*)tend
      indd=index(dum,'STime ')
      if(indd .ne. 0 .and. present(tstart))read(dum(indd+6:),*)tstart
      indd=index(dum,'FILE DATA: ')
      if(indd .ne. 0)exit !exit loop data is following
   end do
case default !generic case assumes no header and thus no metadata

   if(present(lmax)) lmax=-1!give impossible value
   if(present(mu))mu=GM!default Earths gravitational constant
   if(present(rad))rad=RE!default mean Earth radius
   if(present(tstart))tstart=0.d0
   if(present(tcent))tcent=0.d0
   if(present(tend))tend=0.d0

end select

if(.not. stdin)then
close(uni)! close unit

end if
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


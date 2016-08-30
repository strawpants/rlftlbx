!!Coded by Roelof Rietbroek, Tue Sep 22 12:47:56 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Thu Mar 17 15:09:42 2011
!! fixed bug in initializing trend/mean/seas/semi (buggy with gfortran)
!!Updated by Roelof Rietbroek, Wed May 29 15:43:04 2013
!1 added a polynomial expansion to the options
!!Updated by Roelof Rietbroek, Mon Aug 26 21:32:17 2013
!!increased length of dum string
!! added option to squeeze out zero valued columns in the output

!!program which constructs a transformation matrix with added time dependencies



program BIN_TIME_MAT
use binfiletools
use forttlbx
implicit none
integer::iargc,narg,itharg,i,j,ind
character(500)::dum
type(BINdat)::in,out
double precision::t,t0,ohm
integer::chunk,unit
logical::trend,seas,semis,mean,ascii,poly
integer::status,stderr
integer::treg,mreg,areg,sreg,preg,npoly
logical::acc,squeeze
integer,allocatable::ssort(:)

!initialize
ascii=.false.
treg=-1
sreg=-1
areg=-1
mreg=-1
preg=-1
ohm=acos(-1.d0)*2.d0 !seasonal frequency
stderr=0
unit=5
t=-1d4
t0=2005.d0
narg=iargc()
if(narg < 1)call help()
in%file='stdin'
chunk=1000 !initial size of in%side1_d for ascii reads 
trend=.false.
seas=.false.
semis=.false.
mean=.false.
poly=.false.
npoly=-1
squeeze=.false.
!!process command line options


itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('t')!get time from command line
         if(dum(3:3).eq. '0')then
            read(dum(5:),*)t0
         else
            read(dum(4:),*)t
         end if
      case('T')! add trend
         trend=.true.
         if(dum(3:3).eq.'=')treg=regcompf(trim(dum(4:)))

      case('A')!add seasonal
         seas=.true.
         if(dum(3:3).eq.'=')areg=regcompf(trim(dum(4:)))
      case('S')!add semiseasonal
         semis=.true.
         if(dum(3:3).eq.'=')sreg=regcompf(trim(dum(4:)))
      case('M')
         mean=.true.
         if(dum(3:3).eq.'=')mreg=regcompf(trim(dum(4:)))
      case('p')
         poly=.true.
         ind=index(dum,'=')
         if(ind==0)then
            read(dum(3:),*)npoly
         else
            read(dum(3:ind-1),*)npoly
            preg=regcompf(trim(dum(ind+1:)))
         end if
      case('a')!file is an ascii file
         ascii=.true.
      case('h')
         call help()
      case('s')
         squeeze=.true.
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
      in%file=trim(dum)
      unit=13
   end if
end do

!input checks
if(ascii .and. t<-9999)then
   write(stderr,*)"ERROR: ascii input requires the -t flag to be set"
   stop
end if

if( .not. (mean .or. trend .or. seas .or. semis .or. poly))then
   write(stderr,*)"ERROR: no output use at least 1 on -A -T -S -M -p"
   stop
end if

if(poly .and. npoly .eq. -1)then
   write(stderr,*)"ERROR: polynomial not properly initialized"
   stop
end if

!read the input description

if(ascii)then
   allocate(in%side1_d(chunk))
   in%nval1=0
   status=0
   if(in%file .ne. 'stdin')open(file=in%file,unit=unit,form='formatted')
   do while (status .eq. 0)
      in%nval1=in%nval1+1
      if(in%nval1>size(in%side1_d))call realloc_ptr(in%side1_d,chunk)
      read(unit=unit,fmt='(a24)',iostat=status)in%side1_d(in%nval1)
      if(status .ne. 0)then
         in%nval1=in%nval1-1
         exit
      end if
   end do
   if(in%file .ne. 'stdin')close(13)
else
   call read_BINtype(in,2) !only up to stage 2 no data is needed
end if


!try to find the timetag ( if not forced with -t)
if(t<-9999)then
   do,i=1,in%ndbls
      if(in%dbls_d(i) .eq. 'CTime')t=in%dbls(i)
   end do
end if

!check if a time has been found
if(t< -9999)then
      write(stderr,*)"ERROR: no CTime found in file, use -t"
   stop
end if

!now setup output matrix
out%file='stdout'
out%type='FULL2DVN'
out%mtyp='F'
out%ndbls=1
allocate(out%dbls(out%ndbls),out%dbls_d(out%ndbls))
out%dbls_d(1)='CTime'
out%dbls(1)=t

out%nvec=0
out%nint=0
out%nread=0
out%descr="Transformation matrix with added time dependencies from BIN_MAT_TIME"

out%nval2=in%nval1 !amount of columns stays the same
out%side2_d=>in%side1_d

!determine the amount of output rows
out%nval1=0
if(mean)then
   if(mreg < 0)then !don't check for regular expressions
      out%nval1=out%nval1+in%nval1
   else
      do,i=1,in%nval1
         if(.not. regexecf(mreg,in%side1_d(i)))cycle
         out%nval1=out%nval1+1
      end do
   end if
end if

!trend

if(trend)then
   if(treg < 0)then ! don't check for regular expressions
      out%nval1=out%nval1+in%nval1
   else
      do,i=1,in%nval1
         if(.not. regexecf(treg,in%side1_d(i)))cycle
         out%nval1=out%nval1+1
      end do
   end if
end if

!seasonal harmonic ( two added parameters per match)
if(seas)then
   if(areg < 0)then !don't check for regular expressions
      out%nval1=out%nval1+2*in%nval1
   else
      do,i=1,in%nval1
         if(.not. regexecf(areg,in%side1_d(i)))cycle
         out%nval1=out%nval1+2
      end do
   end if
end if

!semi seasonal harmonic ( two added parameters per match)
if(semis)then
   if(sreg < 0)then !don't check for regular expressions
      out%nval1=out%nval1+2*in%nval1
   else
      do,i=1,in%nval1
         if(.not. regexecf(sreg,in%side1_d(i)))cycle
         out%nval1=out%nval1+2
      end do
   end if
end if


!polynomial
if(poly)then
   if(preg < 0)then !don't check for regular expressions
      out%nval1=out%nval1+(npoly+1)*in%nval1
   else
      do,i=1,in%nval1
         if(.not. regexecf(preg,in%side1_d(i)))cycle
         out%nval1=out%nval1+(npoly+1)
      end do
   end if
end if


out%pval1=out%nval1*out%nval2
out%pval2=1
!allocate memory
allocate(out%mat1(out%nval1,out%nval2))
out%mat1=0.d0 !initialize to zero ( matrix will be sparse)
allocate(out%side1_d(out%nval1))

allocate(ssort(out%nval2))
ssort=0

ind=0
do,i=1,in%nval1 ! loop over columns
   if(mean)then
      acc=.true.
      if(mreg >= 0)acc=regexecf(mreg,in%side1_d(i))
      if(acc)then
         ind=ind+1
         out%side1_d(ind)=in%side1_d(i)
         out%mat1(ind,i)=1.d0
         ssort(i)=1
      end if
   end if

   if(trend)then
      acc=.true.
      if(treg >= 0)acc=regexecf(treg,in%side1_d(i))
      if(acc)then
         ind=ind+1
         out%side1_d(ind)='TR_'//in%side1_d(i)(1:21)
         out%mat1(ind,i)=t-t0
         ssort(i)=1
      end if
   end if

   if(seas)then
      acc=.true.
      if(areg >= 0)acc=regexecf(areg,in%side1_d(i))
      if(acc)then
         ind=ind+1
         out%side1_d(ind)='AC_'//in%side1_d(i)(1:21)
         out%mat1(ind,i)=cos(ohm*(t-t0))
         ind=ind+1
         out%side1_d(ind)='AS_'//in%side1_d(i)(1:21)
         out%mat1(ind,i)=sin(ohm*(t-t0))
         ssort(i)=1
      end if
   end if

   if(semis)then
      acc=.true.
      if(sreg >= 0)acc=regexecf(sreg,in%side1_d(i))
      if(acc)then
         ind=ind+1
         out%side1_d(ind)='SC_'//in%side1_d(i)(1:21)
         out%mat1(ind,i)=cos(2*ohm*(t-t0))
         ind=ind+1
         out%side1_d(ind)='SS_'//in%side1_d(i)(1:21)
         out%mat1(ind,i)=sin(2*ohm*(t-t0))
         ssort(i)=1
      end if
   end if

   if(poly)then
      acc=.true.
      if(preg >= 0)acc=regexecf(preg,in%side1_d(i))
      if(acc)then
         do,j=0,npoly !loop over polynomial degree
            ind=ind+1
            write(out%side1_d(ind),'(A1,I1,A1,A21)')'P',j,'_',in%side1_d(i)(1:21)
            out%mat1(ind,i)=(t-t0)**j

         end do
         ssort(i)=1
      end if
   end if
end do

!
if(squeeze)then
   ind=0

   do,i=1,out%nval2
      if(ssort(i) .ne. 0)then
         ind=ind+1
         if(ind .ne. i)then !shift to the left if necessary
            out%mat1(:,ind)=out%mat1(:,i)
            out%side2_d(ind)=out%side2_d(i)
         end if
      end if
   end do
   out%nval2=ind !restrict output
   out%pval1=out%nval1*out%nval2   
end if
!write matrix to file/stdout
call write_BINtype(out)

end program BIN_TIME_MAT

subroutine help()
implicit none
integer::stderr
character(8)::frmt
stderr=0
frmt='(a)'
write(stderr,frmt)"Program BIN_TIME_MAT creates a transformation matrix with added"
write(stderr,frmt)" time dependencies ( mean, trend, (semi)seasonal harmonic)"
write(stderr,frmt)"Usage BIN_TIME_MAT [OPTIONS][INPUTFILE] > OUTMAT"
write(stderr,frmt)"The INPUTFILE is a binary matrix from which the row description is taken."
write(stderr,frmt)" If none is given the program will read from standard input"
write(stderr,frmt)"OPTIONS may be:"
write(stderr,frmt)"-T[=REGEX]: Add a secular time dependency (trend), for the parameters"
write(stderr,frmt)"   which obey the optional regular expression REGEX"
write(stderr,frmt)"-A[=REGEX]: add a seasonal time dependency (annual)"
write(stderr,frmt)"-S[=REGEX]: add a semiseasonal time dependency (semi anual)"
write(stderr,frmt)"-M[=REGEX]: Add a mean (just recycles input parameters)"
write(stderr,frmt)"-pN[=REGEX]: Add a polynomial of degree N"
write(stderr,frmt)"-t=TIME Set time in years for the considered output matrix"
write(stderr,frmt)"-t0=REFTIME Set the reference time in years for the considered output matrix"
write(stderr,frmt)"            Default is 2005"
write(stderr,frmt)"   Default tries to retrieve the time point from the input file"
write(stderr,frmt)"-a The input is an ascii file ( 24 characters per row)"
write(stderr,frmt)"   The -t option must be given in the case of -a"
write(stderr,frmt)"-s: squeeze out unused columns in the output matrix"
stop
end subroutine help

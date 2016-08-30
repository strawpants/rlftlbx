!program to convert SINEX file to BINARY BINFiletools files
!no parameter transformation is performed
!!Coded by Roelof Rietbroek, Fri Apr 10 18:02:31 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

program SNX2BIN
use sinextools
use binfiletools
use gpstlbx
implicit none
integer::itharg,i,narg,stderr
logical::norm,apriori
character(200)::dum
character(10)::date,time
double precision::ctime,etime,stime,rem
type(SNXdat)::snx
type(BINdat)::binout
integer::iargc,shft
!defaults/initializations
norm=.false.
apriori=.false.
stderr=0

!!process command line options
narg=iargc()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('n')!look for normal equation
         norm=.true.
      case('a')!look for an apriori matrix only
         apriori=.true.
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option ( but a file name)
      snx%file=trim(dum)
   end if
end do


!input checks
if(norm .and. apriori)then
   write(stderr,*)'ERROR: normal matrix and apriori matrix cannot be extracted at the same time'
   stop
end if
!read data from sinex file

call read_SNX(snx)

!extract time in gpsweeks

i=GPS_week(snx%tstart,stime)

i=GPS_week(snx%tend,etime)


ctime=(etime+stime)/2


!setup BINdat structure
binout%file='stdout'
binout%type='SYMVN___'
binout%descr='Converted from SINEX by SNX2BIN'

!parameter space
binout%nval1=snx%nval
binout%nval2=binout%nval1
binout%pval1=binout%nval1*(binout%nval1+1)/2
binout%pval2=1
allocate(binout%side1_d(binout%nval1))



!header info

!integers
binout%nint=5
allocate(binout%ints(binout%nint),binout%ints_d(binout%nint))
binout%ints=0

binout%ints_d(1)='Nobs'
binout%ints(1)=snx%nobs
binout%ints_d(2)='Nunknowns'
binout%ints(2)=snx%nval
binout%ints_d(3)='Nunk_orig'
binout%ints(3)=snx%nunk
binout%ints_d(4)='GPSweek'
binout%ints(4)=int(ctime)
binout%ints_d(5)='Nreduced'
binout%ints(5)=snx%nunk-snx%nval


!doubles

binout%ndbls=5+3*snx%nstat !also put lon,lat,height of the stations in the file
allocate(binout%dbls(binout%ndbls),binout%dbls_d(binout%ndbls))
binout%dbls=0.d0
binout%dbls_d(1)="LtPL"
binout%dbls(1)=snx%ltpl
binout%dbls_d(2)="Sigma0"
binout%dbls(2)=snx%sig0_c

binout%dbls_d(3)="CTime"
binout%dbls(3)=GPS_year(9999,ctime)
binout%dbls_d(4)="STime"
binout%dbls(4)=GPS_year(9999,stime)
binout%dbls_d(5)="ETime"
binout%dbls(5)=GPS_year(9999,etime)

!readme part
binout%nread=1
allocate(binout%readme(binout%nread))

call date_and_time(date,time)
binout%readme(1)="created:"//date(7:8)//"-"//date(5:6)//"-"//date(1:4)//' '//time(1:2)//":"//time(3:4)


!construct side description
do, i=1,binout%nval1
   binout%side1_d(i)=snx%para(i)//' '//snx%site(i)//' '//snx%sitec(i)//' '//snx%soluid(i)
   write(binout%side1_d(i)(20:),'(i5)')i ! add index number to parameter string
end do

!allocate vector data
if(apriori)then !no vector data
   binout%nvec=0
else ! add two vectors ( apriori vector and solution/right and side)
   binout%nvec=2
   allocate(binout%vec(binout%nval1,binout%nvec))
   binout%vec=0.d0
end if

if(norm)then
   if(snx%Nmflag .ne. 1)then
      write(stderr,*)"ERROR: SINEX file does not have a normal matrix block"
      stop
   end if 

   if(snx%Nvflag .ne. 1)then
      write(stderr,*)"ERROR: SINEX file does not have a normal vector block"
      stop
   end if

   binout%mtyp=snx%Nmform
   !associate pointer to matrix data
   binout%mat1=>snx%Nm

   !vector data
   

   !right hand side of the normal equation
   do,i=1,binout%nval1
      binout%vec(i,1)=snx%Nv(i)
   end do

   if(snx%Avflag .ne. 1)then
      write(stderr,*)"WARNING: SINEX file does not have a apriori vector"
      write(stderr,*)"apriori estimate set to zero"

   else
      !apriori vector
      do,i=1,binout%nval1
!         write(*,*)binout%side1_d(i),snx%Sv(i),snx%Av(i)
         binout%vec(i,2)=snx%Av(i)
      end do
   end if
else if (apriori)then !extract apriori matrix
   if(snx%Amflag .ne. 1)then
      write(stderr,*)"ERROR: SINEX file does not have a apriori matrix block"
      stop
   end if 

   binout%mtyp=snx%Amform
   !associate pointer to matrix data
   binout%mat1=>snx%Am

else ! extract solution 
   
   binout%mtyp=snx%Smform
   !associate pointer to matrix data 
   binout%mat1=>snx%Sm

   !right hand side of the normal equation
   do,i=1,binout%nval1
      binout%vec(i,1)=snx%Sv(i)
   end do

   
end if

!put station info in meta info doubles
do,i=1,snx%nstat
   shft=(i-1)*3+5
   binout%dbls_d(shft+1)='LON '//snx%domes(i)
   binout%dbls(shft+1)=snx%lon(i)
   binout%dbls_d(shft+2)='LAT '//snx%domes(i)
   binout%dbls(shft+2)=snx%lat(i)
   binout%dbls_d(shft+3)='HEI '//snx%domes(i)
   binout%dbls(shft+3)=snx%height(i)

end do


!write to BINdat type

call write_BINtype(binout)


end program SNX2BIN

subroutine help()
implicit none
integer::unit
character(8)::frmt
unit=0 ! write help to standard error
frmt='(A)'
write(unit,frmt)"Program SNX2BIN: converts a SINEX file to"
write(unit,frmt)" the binary format readable by BINfiletools"
write(unit,frmt)"usage: SNX2BIN [OPTIONS] [SINEXFILE]"
write(unit,frmt)"The program reads from standard input when "
write(unit,frmt)"no SINEXFILE is provided"
write(unit,frmt)""
write(unit,frmt)"OPTIONS may be:"
write(unit,frmt)" -n : Extract the normal equation from the sinex file"
write(unit,frmt)"     This option will look for blocks with SOLUTION/NORMAL_..."
write(unit,frmt)" -a : Extract the apriori matrix ( no vectors extracted) from the sinex file"
write(unit,frmt)"     This option will look for blocks with SOLUTION/MATRIX_APRIORI..."
write(unit,frmt)""
write(unit,frmt)"Output is written to standard output"
write(unit,frmt)""
stop

end subroutine help

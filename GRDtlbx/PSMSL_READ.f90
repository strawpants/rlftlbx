!!Coded by Roelof Rietbroek, Wed Oct 14 10:54:22 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!modified from psmslrd.f ( commented out at the bottom)

!!Updated by Roelof Rietbroek, Tue Oct 26 16:21:34 2010
!! updated for new data format
!! removed transpose option
!! write in binary form
!!Updated by Roelof Rietbroek, Fri Oct 29 10:11:57 2010
!!Added slice option (extract slice at one specific time point)

!!Updated by Roelof Rietbroek, Thu Sep 22 08:09:57 2011
!! added correct behavior when search area contains the meridian
!!Updated by Roelof Rietbroek, Fri Oct  5 15:28:38 2012
!! fixed bug in time restricted search



program PSMSL_READ
use forttlbx
use binfiletools
use psmsl_tool
implicit none
character(200)::dum
integer::iargc,narg,itharg
integer::regind,ind,indo
type(PSMSL_struc),target::dat
type(PSMSL_stat),pointer::stat ! used for shortcutting
integer::stderr,stdout,nstat,i,j,k,iyr,rec
double precision::lonmin,lonmax,latmin,latmax
double precision::lons
double precision::etime,stime
double precision:: slicet,slicet_a,slicet_b,w_a,w_b
integer::gaplen,ngaps
integer::st,nd,nmatch
character(30)::frmt
integer::quitafter
integer,pointer,dimension(:)::tmp,tmp2
character(200)::rootdir
integer,allocatable::valid(:)
logical::t_restr,g_restr
integer::output,step
double precision::sigma0
integer::mis_a,mis_b
type(BINdat)::out
logical::contains_greenwhich
!defaults/initializations

output=2 ! 0 for info 1 for more info 2 for data
stderr=0
stdout=6
quitafter=9999999 ! quit after we have found  quitafter valid stations
stime=0 ! extreme values ( will be adjusted)
etime=5000 ! extreme values ( will be adjusted)
sigma0=1d-3 ! default assumes 1 mm error

lonmin=0
lonmax=360
latmin=-90
latmax=90
gaplen=99999 ! initially very large tolerance for gaps
regind=0
t_restr=.false.
g_restr=.false.

contains_greenwhich=.false.

!get default PSMSL root
call getenv('WORK_DIR',dum)
rootdir=trim(dum)//'/data/PSMSL_new/'


!default filelist metric monthly data

!!process command line options
narg=iargc()
if( narg < 1)call help(dat%dir)
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('i')! print station info (name + geographical position)
         output=0
         if(dum(3:3) .eq. '+')output=1
      case('d')!also ouptput missing days
         output=3
      case('q')! quit reading the file after first match (speeds uop search for single stations)
         if(dum(3:3).eq. '=')then
            read(dum(4:),*)quitafter
         else
            quitafter=1
         end if
      case('S')!extract time slice
         if(dum(3:3).eq. '=')then
            ind=index(dum,'/')
            if(ind .eq. 0)then
               read(dum(4:),*)slicet
            else
               read(dum(4:ind-1),*)slicet
               read(dum(ind+1:),*)sigma0
            end if
         else
            write(stderr,*)"ERROR processing -S option:",trim(dum)
            stop
         end if
         output=4
      case('T') ! restrict in time
         ind=index(dum,'/')
         read(dum(4:ind-1),*)stime
         read(dum(ind+1:),*)etime
         t_restr=.true.
      case('R')! restrict in space
         ind=index(dum,'/')
         indo=3
         read(dum(indo+1:ind-1),*)lonmin
         
         indo=ind
         ind=ind+index(dum(ind+1:),'/')
         read(dum(indo+1:ind-1),*)lonmax

         if(lonmin .eq. -180 .and. lonmax .eq. 180)then !exception rule indicating searching the complete domain
            lonmin=0
            lonmax=360
         end if
         if(lonmax <0)lonmax=lonmax+360
         if(lonmin <0)lonmin=lonmin+360

         indo=ind
         ind=ind+index(dum(ind+1:),'/')
         read(dum(indo+1:ind-1),*)latmin
         
         read(dum(ind+1:),*)latmax
         !start crying when boundaries are not ok
         if(latmin>=latmax)then
            write(stderr,*)'ERROR: Minimum latitude >= Maximum latitude?'
            stop
         end if


         if(lonmin>lonmax)then
            contains_greenwhich=.true.
            write(stderr,*)'Message: Speciaal voor Laura met Greenwhich Meridiaan! Groetjes, Roelof'
         end if

         if(contains_greenwhich)then !transform left boundary back to minus number again
            lonmin=lonmin-360
         end if
!         write(0,*)lonmin,lonmax,latmin,latmax
      case('N')!restrict station name/code
         !compile regular expression
         regind=regcompf(trim(dum(4:)))
      case('C')! continuous observation period
         if(dum(3:3) .eq. '=')then
            read(dum(4:),*)gaplen
         else
            gaplen=0
         end if
         g_restr=.true.
      case('h')
         call help(dat%dir)
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is the data directory or shortcut
      if(dum .eq. 'RLR')then
         dat%dir=trim(rootdir)//'rlr_monthly/'
      else if (dum .eq. 'RLA')then
         dat%dir=trim(rootdir)//'rlr_annual/'
      else if(dum .eq. 'MET')then
         dat%dir=trim(rootdir)//'met_monthly/'
      else
         dat%dir=trim(dum)
      end if
   end if
end do


if(dat%dir .eq. '')then
   write(stderr,*)" ERROR: no database specified"
   stop
end if

if(stime >etime)then
   write(stderr,*)"ERROR: requested start time is later than end time"
   stop
end if


!initialize database (read file inventory and set coastal codes, minimum and maximum times,..)
call psmsl_init(dat)
allocate(valid(dat%nstat))
valid=0

!snap etime/stime to nearest datapoint(center) and move to boundary 
if(t_restr)then
   stime=int(stime/dat%step)*dat%step ! snap to data point below
   etime=int(etime/dat%step)*dat%step+dat%step ! snap to data point above
end if

!determine boundaries for interpolation of time slice
if(output .eq. 4)then
   slicet_a=int((slicet+dat%step/2)/dat%step)*dat%step-dat%step/2 ! snap to data point below
   slicet_b=slicet_a+dat%step
end if


!Apply possible search criteria
do,i=1,dat%nstat
   stat=>dat%stat(i) ! shortcut
!first criteria which don't require a data/info read
   !region
   
   if(stat%lon < 0)then
      lons=stat%lon+360 ! make longitude consistent with 0-360 bounds
   else
      lons=stat%lon
   end if
   if(contains_greenwhich .and. lons > 180 )lons=lons-360 !center data on Greenwhich meridian

   if(lons > lonmax .or. lons < lonmin)cycle
   
   if(stat%lat > latmax .or. stat%lat < latmin)cycle

   !regular expression in name or coastline
   if(regind>0)then
      if(.not. &
           regexecf(regind,stat%sid//stat%name//stat%coast))cycle
   end if
   
   !rectrict on time constraints
!   if(t_restr .and. (stat%t0 > stime .or. stat%tn < etime))cycle !( no data in requested period)
   if(t_restr .and. (stat%t0 > etime .or. stat%tn < stime))cycle !( no data in requested period)

  

   !load station doc and data (if needed)

   if(g_restr .or. output >= 2)call psmsl_rdstat(dat,i,.false.) ! data is needed
   if(output .eq. 1)call psmsl_rdstat(dat,i,.true.) !documentation needed


   
   if(g_restr)then
      !restrict on (amount of) data gaps
      if(t_restr)then
         st=psmsl_tindex(dat,i,stime+dat%step/2)
         nd=psmsl_tindex(dat,i,etime-dat%step/2)
      else if(output .eq. 4)then ! check for missing values in the interpolation values
         st=psmsl_tindex(dat,i,slicet_a)
         nd=st+1
      else
         st=1
         nd=stat%ndat
      end if

      !count amount of datagaps
      ngaps=0
      do,j=st,nd
         if(st < 1 .or. nd > stat%ndat)then ! out of time range (count as gap)
            ngaps=ngaps+1
         else if(stat%H(j) .eq. -99999)then
            ngaps=ngaps+1
         end if
      end do
      if(ngaps > gaplen)cycle 
   end if

   valid(i)=1 ! tag as a valid station
   if(sum(valid) .eq. quitafter)exit !quick exit
end do
nmatch=sum(valid)
if(nmatch .eq. 0)then
   write(stderr,*)"NO tide gauges satisfying criteria found"
   stop
else
   write(stderr,*)"Found", nmatch," matching Tide gauges"
end if

select case(output)
case(0) !print plain TG station info
    do,i=1,dat%nstat
       stat=>dat%stat(i)
       if(valid(i) .eq. 1)write(stdout,'(A24,1x,2F7.2)')stat%sid//' '//stat%name(1:18),stat%lon,stat%lat
    end do
case(1)!more documentation will be printed
   do,i=1,dat%nstat
      stat=>dat%stat(i)
      if(valid(i) .eq. 1)then
         write(stdout,*)">>>"
         write(stdout,'(A46,1x,2F7.2)')stat%sid//' '//stat%name,stat%lon,stat%lat
         write(stdout,'(A60)')stat%coast
         write(stdout,'(a12,f10.4,1x,a1,1x,f10.4)')"Time period:",stat%t0,"-",stat%tn
         write(stdout,'(a19,1x,a1)')"Flag for attention:",stat%qcflag
         do,j=1,stat%ncom
            write(stdout,'(A80)')stat%comm(j)
         end do
         write(stdout,*)"<<<"
      end if
   end do
case(2,3) ! write data to binary file
   !set up output structure
   out%file='stdout'
   out%descr='Results from TG query of PSMSL database'
   out%ndbls=4
   allocate(out%dbls(out%ndbls),out%dbls_d(out%ndbls))
   out%dbls_d(1)='STime'
   out%dbls_d(2)='CTime'
   out%dbls_d(3)='ETime'
   out%dbls_d(4)='TStep'

   if(.not. t_restr)then !automatically retrieve from restricted dataset
      
      !retrieve minimum and maximum time of restricted set
      etime=0 ! to be replcaed 
      stime=5000 !to be replaced
      do,i=1,dat%nstat
         if(valid(i) .eq. 0)cycle
         stat=>dat%stat(i)
         etime=max(etime,stat%tn)
         stime=min(stime,stat%t0)
      end do
   end if
!   write(0,*)dat%stat(1)%t0,stime,etime,psmsl_tindex(dat,1,etime),psmsl_tindex(dat,1,stime)
   out%nval2=psmsl_tindex(dat,1,etime-dat%step/2)-psmsl_tindex(dat,1,stime+dat%step/2)+1
   if(output .eq. 3) then
      step=2
   else
      step=1
   end if

   out%nval2=out%nval2*step ! also include missing days
 



   out%dbls(1)=stime
   out%dbls(3)=etime
   out%dbls(2)=(stime+etime)/2.d0
   out%dbls(4)=dat%step

   out%nval1=sum(valid) ! amount of valid stations
   out%pval1=out%nval1*out%nval2
   out%pval2=1
   out%type='FULL2DVN'
   out%mtyp='F'

! put data in array
   allocate(out%mat1(out%nval1,out%nval2))
   out%mat1=-99999 ! defualt to mi9ssing value
   out%nvec=2 !holds longitude and latitude
   allocate(out%vec(out%nval1,out%nvec))
   
   
   allocate(out%side1_d(out%nval1))
   allocate(out%side2_d(out%nval2))
   st=0
   do,i=1,dat%nstat
      if(valid(i) .eq. 0)cycle
      stat=>dat%stat(i)
      st=st+1
      out%side1_d(st)=stat%sid//' '//stat%name(1:18)
      out%vec(st,1)=stat%lon
      out%vec(st,2)=stat%lat
   end do
   
   !set up time descriptor of column ( time in year) tags are center time of observation
   if(output .eq. 3)then ! double output 
      do,i=1,out%nval2,step
         write(out%side2_d(i),'(A4,F10.5)')'TYR_',stime+(i-1)/step*dat%step+dat%step/2
         write(out%side2_d(i+1),'(A4,F10.5)')'MIS_',stime+(i-1)/step*dat%step+dat%step/2 ! missing day tag 
      end do
   else
      do,i=1,out%nval2
         write(out%side2_d(i),'(A4,F10.5)')'TYR_',stime+(i-1)*dat%step+dat%step/2
      end do
   end if
   !fill up matrix data

   st=0
   do,i=1,dat%nstat
      if(valid(i) .eq. 0)cycle
      stat=>dat%stat(i)
      st=st+1
!       if(output .eq. 2)then
!          tmp=>stat%H
!       else
!          tmp=>stat%misd
!       end if

      do,j=1,out%nval2,step
         nd=psmsl_tindex(dat,i,stime+(j-1)/step*dat%step+dat%step/2)
!         write(0,*)nd,stat%ndat,dat%yearly
         if(nd < 1 .or. nd>stat%ndat)cycle !out of range
         out%mat1(st,j)=1d-3*stat%H(nd)
         if(output .eq. 3)out%mat1(st,j+1)=stat%misd(nd)
      end do
   end do
   call write_BINtype(out)

case(4) !output slice
   !set up output structure
   out%file='stdout'
   out%descr='Interpolated time slice from the PSMSL database'
   out%ndbls=4
   allocate(out%dbls(out%ndbls),out%dbls_d(out%ndbls))
   out%dbls_d(1)='STime'
   out%dbls_d(2)='CTime'
   out%dbls_d(3)='ETime'
   out%dbls_d(4)='TStep'


   out%nval1=sum(valid) ! amount of valid stations
   out%nval2=out%nval1

   out%dbls(1)=slicet-dat%step/2
   out%dbls(2)=slicet
   out%dbls(3)=slicet+dat%step/2
   out%dbls(4)=dat%step

   out%pval1=out%nval1
   out%pval2=1
   out%type='DIAVN___'
   out%mtyp='P'

   out%nvec=3

   !allocate arrays
   allocate(out%pack1(out%pval1))
   allocate(out%vec(out%nval1,out%nvec))
   out%vec=-99999 ! default to missing value
   out%pack1=sigma0**2 ! start with a scaled unit diagonal

   allocate(out%side1_d(out%nval1))
   out%side2_d=>out%side1_d ! point to the same description



   !calculate interpolation weights
   w_a=(slicet_b-slicet)/dat%step
   w_b=(slicet-slicet_a)/dat%step

   st=0
   do,i=1,dat%nstat
      if(valid(i) .eq. 0)cycle
      stat=>dat%stat(i)
      st=st+1
      out%side1_d(st)=stat%sid//' '//stat%name(1:18)
      out%vec(st,2)=stat%lon
      out%vec(st,3)=stat%lat

      !determine index of datapoint 'A'
      nd=psmsl_tindex(dat,i,slicet_a)
      !selection criteria
      if(nd < 1 .or. nd+1>stat%ndat)cycle !iether one or all data points is out of range
      if(stat%H(nd) .eq. -99999 .or. stat%H(nd+1) .eq. -99999)cycle ! one of the datapoints is missing
      
      !put goodies in arrays ( and scale to meter)
      out%vec(st,1)=(stat%H(nd)*w_a+stat%H(nd+1)*w_b)*1d-3     
      !also propagate error based on missing days
      mis_a=stat%misd(nd)
      mis_b=stat%misd(nd+1)
      if(abs(mis_a) > 31)mis_a=30 ! set to 30 missing days for interpolated values
      if(abs(mis_b) > 31)mis_b=30 ! set to 30 missing days for interpolated values
      out%pack1(st)=((31/(31-mis_a))*w_a+(31/(31-mis_b))*w_b)*out%pack1(st)
   end do
   
   call write_BINtype(out)
end select


end program PSMSL_READ

subroutine help(dir)
implicit none
character(8)::frmt
character(*)::dir
integer::unit
unit=0
frmt='(A)'

write(unit,frmt)"PSMSL_READ reads the PSMSL data base and applies selection criteria"
write(unit,frmt)"Usage PSMSL_READ [OPTIONS] DATABASE"
write(unit,frmt)"DATABASE is the PSMSL root directory. specify either the full dir or"
write(unit,frmt)"the shortcuts RLR:"//trim(dir)//'rlr_monthly (Revised Local Reference)'
write(unit,frmt)"              RLA:"//trim(dir)//'rlr_annual  (Revised Local Reference, yearly)'
write(unit,frmt)"              MET:"//trim(dir)//'met_monthly (metric dataset)'
write(unit,frmt)"OPTIONS may be:" 
write(unit,frmt)"-i[+]: print brief station name and location" 
write(unit,frmt)"  Append the + sign to get more detailed info on stations"
write(unit,frmt)"-d: Also output the amount of missing days per station."
write(unit,frmt)"    note: -i and -d may not be combined"
write(unit,frmt)"    matrix data will have extra time columns containing the missing days"
write(unit,frmt)"-S=CTIME[/SIGMA0]: extract a single time slice (over the TG network), centered on"
write(unit,frmt)"          CTIME, from the database in the form of an observation matrix"
write(unit,frmt)"          Data will be (linearly) interpolated from yearly/monthly values."
write(unit,frmt)"          Output is a diagonal matrix in binary form, representing errors derived"
write(unit,frmt)"          from the amount of missing days. 3 associated vectors containing"
write(unit,frmt)"          the tide gauge observations, lon and lat are also supplied"
write(unit,frmt)"          If provided scale the diagonal by SIGMA0^2 (default SIGMA0=1e-3: 1mm)"
write(unit,frmt)"          In combination with -C: this will ignore data with missing values"
write(unit,frmt)"          in the boundaries of the time slice or with a longer period ( with -T)"
!write(unit,frmt)"-t: transpose output"
!write(unit,frmt)"-n=VALUE: set the NODATA nodes to VALUE ( default is 99999)"
write(unit,frmt)"-q[=NVALID]: quick exit after a matching station has been found"
write(unit,frmt)"            When NVALID is given the program quits after NVALID stations"
write(unit,frmt)"The following search criteria may be used to restrict the output:"
write(unit,frmt)"-N=REGEX: only output stations which name and/or country (code and name)"
write(unit,frmt)"          match the regular expression REGEX"
write(unit,frmt)"-R=LONMIN/LONMAX/LATMIN/LATMAX : only output stations which are contained within"
write(unit,frmt)"    the geographical region bounded by LONMIN/LONMAX/LATMIN/LATMAX ( in degrees)"
write(unit,frmt)"    Negative longitudes are allowed"
write(unit,frmt)"-T=YSTART/YEND: Analyze the data only in a restricted time span"
write(unit,frmt)"                to the years between YSTART/YEND"
write(unit,frmt)"                Stations with no data in this period will be skipped."
! write(unit,frmt)"-L=LENGTH: select stations which have at least LENGTH nodes of data (years or months)"
! write(unit,frmt)"             The period considered maybe restricted by using the -T option"
write(unit,frmt)"-C[=GAP]: Reject stations which are not continuous in the requested period"
write(unit,frmt)"             If GAP is provided up to GAP missing nodes are tolerated."
write(unit,frmt)"Results are printed to standard output (units meter) in binary form (readable by BIN_swiss)" 

stop
end subroutine help


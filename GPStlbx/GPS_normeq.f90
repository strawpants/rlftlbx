!program to create normal equations for  SH-surface load coefficients from sinex files and a reference (in time) 
!uses a namelist setup to exclude or include data entries with specific string parameters. For example one can chose to include only Height data points for example or eclude specifi stations.

!!Coded by Roelof Rietbroek, Fri Oct 19 15:35:13 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Mon Mar 17 15:08:17 2008
!!added the subtraction of a mean


program GPS_normeq
use GPStlbx
use SHtlbx
use FORTtlbx
implicit none
integer::narg,itharg,i,funit,nmax,nfiles,nstat,row,col,nunk,ind,gpswk,ndat
integer::l,m,q,lmax,lmin,f,lond,lonm,latd,latm,info,stderr,last,j,shift,mmax
parameter(mmax=1000)
parameter(nmax=800)!maximum amount of input residual files
double precision::lons,lats,lontmp,lattmp,time1,time2
double precision, allocatable,dimension(:,:)::A,cov,AtA
double precision, allocatable,dimension(:)::lat,lon,Ab,residu ,apriori
character(24),allocatable,dimension(:)::Unlist
character(80),allocatable,dimension(:)::list
character(17),allocatable,dimension(:)::domes
character(80)::outfile,descr
character(80),dimension(4)::dum80
character(4)::charwk
character(120)::dum,incfile,excfile,meanfile
character(120),dimension(nmax)::resfile
character(1000)::inclstr,exclstr !long strings holding selection criteria
logical::diag,incl,excl,helm,geocent,ancyc,mean,sancyc,sh,trend,submean
integer::ints(7),nsh,nmean
double precision::doubles(6),LtPL,time,meanstat(mmax)
character(24)::ints_d(7),doubles_d(6),meanlist(mmax)
integer::iargc
!defaults initializations
stderr=0
funit=15
narg=0
itharg=0
diag=.false.!diagonal system (default takes full covariance)
nfiles=0
outfile='GPSnormeq2'
descr='GPS NORMAL system from GPS_normeq'
sh=.false.
helm=.false.
geocent=.false.
ancyc=.false.
sancyc=.false.
mean=.false.
trend=.false.

incl=.false. !if include terms are present (default includes all entries)
excl=.false. !if exclude terms are present (default excludes nothing)
inclstr=''
exclstr=''
lmin=1!default lmin 
lmax=7!and lmax values
incfile=''
excfile=''
shift=0

!get command line arguments
narg=iargc()


!process command line arguments
if(narg < 1) call help()

do,i=1,narg
   itharg=itharg+1
   if(itharg> narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !option
      select case( dum(2:2))
      case('i')! include terms from file
         incl=.true.
         if(dum(3:3) .ne. ' ')then
            incfile=dum(3:)
         else
            itharg=itharg+1
            call getarg(itharg,incfile)
         end if

      case('e')! exclude terms from file
         excl=.true.
         if(dum(3:3) .ne. ' ')then
            excfile=dum(3:)
         else
            itharg=itharg+1
            call getarg(itharg,excfile)
         end if
      case('I')!append term to include string
         if(inclstr .eq. '')then
            inclstr='*'//trim(dum(3:))
         else
            inclstr=trim(inclstr)//'*|*'//trim(dum(3:))
         end if
         incl=.true.
      case('E')!append term to exclude string
         if(exclstr .eq. '')then
            exclstr='*'//trim(dum(3:))
         else
            exclstr=trim(exclstr)//'*|*'//trim(dum(3:))
         end if
         excl=.true.
      case('l')!limit maximum and minimum degree
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(3:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(3:),*)lmax
         end if
      case('t')!incorporate Helmert Transformation parameters in the normal equations
         helm=.true.
         
      case('g')!describe degree 1 variation as geocenter motion
         geocent=.true.
         
      case('s')!create system with Spherical harmonic elastic loading
         sh=.true.
         
      case('m')!calculate a mean per station
         mean=.true.
         if(dum(3:3) .eq. 't')then
           
            trend=.true.
         end if
         
      case('a')!fit annual and semi annual cycle
         ancyc=.true.         
         if(dum(3:3) .eq. 's')then!semi annual cycle
            sancyc=.true.
         end if
      case('M') !subtract mean file
         submean=.true.
         if(dum(3:3) .ne. ' ')then
            meanfile=dum(3:)
         else
            itharg=itharg+1
            call getarg(itharg,meanfile)
         end if
      case('d')!apply simple least squares (diagonal only)
         diag=.true.
      case('F')!use a different basename
         outfile=dum(3:)
      case default
         write(stderr,*)'Unknown option selected, quitting'
         call help()
      end select
   else !append residual file to list
      nfiles=nfiles+1
     resfile(nfiles)=dum
   end if
end do

!input checks
if(nfiles <1)then
   write(stderr,*)'No residual file supplied, quitting'
   stop
end if

if(geocent .and. lmin <=1)then
   lmin=2 !reset minimum degree
end if

if(.not.(sh .or. geocent .or. helm .or. mean .or. ancyc .or. sancyc .or.trend))then
   write(stderr,*)'At least one of -t, -s, -a, -m, must be specified'
   write(stderr,*)'type GPS_normeq -h for more info'
   stop
end if



!load include terms from file
last=0
if(incfile .ne. '')then
   open(unit=funit,file=trim(incfile))
   do, i=1,nmax
      
      read(unit=funit,fmt=*,iostat=last)dum
      if(last .ne. 0) exit
      !append to include string
      if(inclstr .eq. '')then
         inclstr='*'//trim(dum)
      else
         inclstr=trim(inclstr)//'*|*'//trim(dum)
      end if
      
   end do

   
   close(funit)
end if


!load exclude terms from file
last=0
if(excfile .ne. '')then
   open(unit=funit,file=trim(excfile))
   do, i=1,nmax
     
      read(unit=funit,fmt=*,iostat=last)dum
      if(last .ne. 0) exit
      !append to exclude string
      if(exclstr .eq. '')then
         exclstr='*'//trim(dum)
      else
         exclstr=trim(exclstr)//'*|*'//trim(dum)
      end if

      
   end do

   
   close(funit)
end if

!add wild character at the back of incl and excl string
if(inclstr .ne. '')inclstr=trim(inclstr)//'*'
if(exclstr .ne. '')exclstr=trim(exclstr)//'*'

!write(*,*)trim(inclstr)
!write(*,*)trim(exclstr)
!stop


!make meta data headers for files
ints=0
ints(6)=1
if(sh)then
   ints(4)=lmax
   ints(5)=lmin
end if

ints_d(1)='Nobs'
ints_d(2)='Nunknows'
ints_d(3)='Nreduced'
ints_d(4)='Lmax'
ints_d(5)='Lmin'
ints_d(6)='Modnr'
ints_d(7)='GPSweek'

doubles=0.d0
doubles(6)=1.d0
doubles(4)=1.d0

doubles_d(1)='CTime'
doubles_d(2)='STime'
doubles_d(3)='ETime'
doubles_d(4)='Weight'
doubles_d(5)='LtPL'
doubles_d(6)='Sigma0'


!initial allocation (should be enough for most scenarios ansd if not this will be readjusted autmatically)
ndat=500
allocate(residu(ndat),list(ndat),lat(ndat),lon(ndat),cov(ndat,ndat),domes(ndat))


nunk=4000
allocate(AtA(nunk,nunk),Ab(nunk),apriori(nunk), Unlist(nunk))
apriori=0.d0
Unlist=''      
allocate(A(ndat,nunk))
A=0.d0 !(initialize)


!load mean to subtract
if(submean)then
   open(unit=funit,iostat=last,file=trim(meanfile),form='formatted')
   meanlist=''
   meanstat=0.d0
   last=0
   do,i=1,mmax
      read(unit=funit,fmt='(A24)',advance='NO',iostat=last)meanlist(i)
      if(last .ne. 0)then
         nmean=i-1
         exit
      end if
      read(unit=funit,fmt=*,iostat=last)meanstat(i)
   end do
   close(funit)
end if


!loop over SINEX residual files

do,f=1,nfiles
  
   !retrieve amount of data
   if(incl .and. excl)then !both include and exclude criteria are present
      call get_sinex(file=trim(resfile(f)),nx=ndat,inc1=trim(inclstr),excl1=trim(exclstr),gpswk=gpswk)
   else if(incl)then !only include terms present
      call get_sinex(file=trim(resfile(f)),nx=ndat,inc1=trim(inclstr),gpswk=gpswk)
   else if(excl)then !only exclude terms are present
      call get_sinex(file=trim(resfile(f)),nx=ndat,excl1=trim(exclstr),gpswk=gpswk)
   else
      call get_sinex(file=trim(resfile(f)),nx=ndat,gpswk=gpswk)
   end if

!   convert gps week to years
   time=GPS_year(gpswk)

   !reallocate arrays only if necessary
   if(size(residu,1)< ndat)then
      deallocate(residu,list,lat,lon,domes,cov)
      allocate(residu(ndat),list(ndat),lat(ndat),lon(ndat),cov(ndat,ndat),domes(ndat))
   end if
   
   
   
   
   !read in residuals and covariance

   if(incl .and. excl)then !both include and exclude criteria are present
      call get_sinex(file=trim(resfile(f)),x=residu(1:ndat),xlist=list(1:ndat),covxy=cov(1:ndat,1:ndat),&
           inc1=trim(inclstr),excl1=trim(exclstr),inc2=trim(inclstr),excl2=trim(exclstr))
   else if(incl)then !only include terms present
      call get_sinex(file=trim(resfile(f)),x=residu(1:ndat),xlist=list(1:ndat),covxy=cov(1:ndat,1:ndat)&
           ,inc1=trim(inclstr),inc2=trim(inclstr))
   else if(excl)then !only exclude terms are present
      call get_sinex(file=trim(resfile(f)),x=residu(1:ndat),xlist=list(1:ndat),covxy=cov(1:ndat,1:ndat)&
           ,excl1=trim(exclstr),excl2=trim(exclstr))
   else
      call get_sinex(file=trim(resfile(f)),x=residu(1:ndat),xlist=list(1:ndat),covxy=cov(1:ndat,1:ndat))
   end if

  
  !  do,i=1,ndat
!       do,j=i,ndat
!          write(*,*)i,j,covpack(i+((j-1)*j)/2)
!       end do
!    end do

   !!construct latitude and longitude tags for each entry (station positions are retrieved from the SITE/ID SINEX block)
   do, i=1,ndat
!      call cpu_time(time1)
      call get_sinexblock(file=trim(resfile(f)),blockname='SITE/ID',&
           list=dum80(1:4)(1:80),nlines=nstat,inc1='*'//list(i)(15:21)//'*')
 !     call cpu_time(time2)
  !    write(*,*)'check',nstat,i,time2-time1
      if(nstat .eq. 0)then
         !try removing the B or A
         call get_sinexblock(file=trim(resfile(f)),blockname='SITE/ID',&
              list=dum80(1:4)(1:80),nlines=nstat,inc1='*'//list(i)(15:19)//'*')
         if(nstat .eq. 0)then!still not found
            write(stderr,*)'ERROR: Station, ',list(i)(15:21),' Not found in SITE/ID block'
            stop
         end if
      else if(nstat > 1)then
         write(stderr,*)'WARNING: Station, ',list(i)(15:21),' has ambiguous entries in SITE/ID block'
         write(stderr,*)'Taking.. ',dum80(1)
!         stop
      end if
      
      !add domes number to character array
      domes(i)=dum80(1)(1:8)//dum80(1)(10:18)
      read(dum80(1)(44:67),'(i4,1x,i2,1x,f4.1,1x,i3,1x,i2,1x,f4.1)')lond,lonm,lons,latd,latm,lats
!      write(*,*)'test'//dum80(1)(44:67)
      lontmp=dble(lond)+sign(lonm,lond)/60.d0
      lon(i)=(lontmp+sign(lons,lontmp)/3600.d0)*pi/180.d0
!      lon(i)=(dble(lond)+lonm/60.d0+lons/3600.d0)*pi/180.d0
      lattmp=dble(latd)+sign(latm,latd)/60.d0
      lat(i)=(lattmp+sign(lats,lattmp)/3600.d0)*pi/180.d0
      !lat(i)=(dble(latd)+latm/60.d0+lats/3600.d0)*pi/180.d0

   end do

   !subtract mean if requested
   if(submean)then
      do,i=1,ndat
         if(compstr(list(i),'*STAH*'))then
            do,j=1,nmean
               if(compstr(meanlist(j),'SH*'//list(i)(15:21)//'*'))then
                  residu(i)=residu(i)-meanstat(j)
                  exit
               end if
            end do

         else if(compstr(list(i),'*STAE*'))then 
            do,j=1,nmean
               if(compstr(meanlist(j),'SE*'//list(i)(15:21)//'*'))then
                  residu(i)=residu(i)-meanstat(j)
                  exit
               end if
            end do
         else if(compstr(list(i),'*STAN*'))then 
            do,j=1,nmean
               if(compstr(meanlist(j),'SN*'//list(i)(15:21)//'*'))then
                  residu(i)=residu(i)-meanstat(j)
                  exit
               end if
            end do
         end if
      end do
   end if



!!!!!!!!!!!!below we construct the observation equation matrix

!calculate amount of unknowns

   
   nunk=0
   nsh=SH_tpos(lmax,lmax,1,lmax,lmin)
   if(sh)nunk=nunk+nsh
   if(geocent)nunk=nunk+3
   if(helm)nunk=nunk+7
   if(mean)nunk=nunk+ndat
   if(ancyc)nunk=nunk+2*ndat
   if(sancyc)nunk=nunk+2*ndat
   if(trend)nunk=nunk+ndat

!allocate space for normal equations 

   if(size(Ab,1) < nunk)then
      deallocate(AtA,Unlist,apriori,Ab)
      allocate(AtA(nunk,nunk),Ab(nunk),apriori(nunk), Unlist(nunk))
      apriori=0.d0
      Unlist=''      
   end if

   if(size(A,1) < ndat .or. size(A,2) < nunk)then
      deallocate(A)
      allocate(A(ndat,nunk))

      A=0.d0 !(initialize)
   end if


   !get observation equation matrix:
   shift=0
   if(helm)then
      call GPS_obseq_helmert(A=A(1:ndat,shift+1:shift+7),tag=Unlist(shift+1:shift+7),typ=list(1:ndat)&
           ,lon=lon(1:ndat),lat=lat(1:ndat))
      shift=shift+7
   end if
   
   if(geocent)then
      call GPS_obseq_geocent(A=A(1:ndat,shift+1:shift+3),tag=Unlist(shift+1:shift+3),typ=list(1:ndat),&
           lon=lon(1:ndat),lat=lat(1:ndat))
      shift=shift+3
   end if
   
   if(mean)then
      call GPS_obseq_mean(A=A(1:ndat,shift+1:shift+ndat),tag=Unlist(shift+1:shift+ndat),&
           typ=list(1:ndat),domes=domes(1:ndat))
      shift=shift+ndat 
   end if

   if(trend)then
      call GPS_obseq_trend(A=A(1:ndat,shift+1:shift+ndat),tag=Unlist(shift+1:shift+ndat),&
           time=time,typ=list(1:ndat),domes=domes(1:ndat))
      shift=shift+ndat 
   end if

   if(ancyc)then
      !annual cycle
      call GPS_obseq_cycle(A=A(1:ndat,shift+1:shift+2*ndat),tag=Unlist(shift+1:shift+2*ndat)&
           ,time=time,ohm=2*pi,harmtag='AN',typ=list(1:ndat),domes=domes(1:ndat))
      shift=shift+2*ndat
   end if

   if(sancyc)then
      !semiannual cycle
      call GPS_obseq_cycle(A=A(1:ndat,shift+1:shift+2*ndat),tag=Unlist(shift+1:shift+2*ndat)&
           ,time=time,ohm=4*pi,harmtag='SA',typ=list(1:ndat),domes=domes(1:ndat))
      shift=shift+2*ndat
   end if   

!   write(*,*)sh,ndat,shift,nsh,lmax,lmin
   if(sh)then

      call GPS_obseq_loadsh(A=A(1:ndat,shift+1:shift+nsh),tag=Unlist(shift+1:shift+nsh),&
           typ=list(1:ndat),lon=lon(1:ndat),lat=lat(1:ndat),lmax=lmax,lmin=lmin)     
      shift=shift+nsh
   end if


   
   
!   call GPS_obseq(A=A,typ=list,lon=lon,lat=lat,helm=helm,geocent=geocent,lmax=lmax,lmin=lmin)
   

!    do,i=1,ndat
!       write(*,*)list(i)(8:19),A(i,11:21)
!    end do
!    write(*,*)Unlist(11:21)(1:10)   
   !contruct normal equations
   Ab(1:ndat)=0.d0
   AtA(1:ndat,1:ndat)=0.d0


   !invert covariance matrix using a cholesky decomposition

   if(diag)then !only use diagonal

      do,col=1,ndat
         do,row=1,col
            if(row .ne. col)then
               cov(row,col)=0.d0
               cov(col,row)=0.d0
            else !diagonal
               cov(row,col)=1/sqrt(cov(row,col))
            end if
         end do
      end do

   else !full covariance

      write(stderr,*)'Calculating Cholesky factorization of observation covariance matrix..'


      call dpotrf('U',ndat,cov(1:ndat,1:ndat),ndat,info )
      if(info .ne. 0)call xerbla('DPOTRF',info)


      write(stderr,*)'done'
   end if

   


!solve triangular system U**T  B= A for B where (U**(T) is the cholesky fact of cov
!DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
  call DTRSM('L','U','T','N',ndat,nunk,1.d0,cov(1:ndat,1:ndat),ndat,A(1:ndat,1:nunk),ndat)
!A is overwritten!

   write(stderr,*)'Making left hand of normal equations (AtA)..'
!   reuse A (is now U**-T A 
!   DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
   call DSYRK('U','T',nunk,ndat,1.d0,A(1:ndat,1:nunk),ndat,0.d0,AtA(1:nunk,1:nunk),nunk)
   !AtA is overwritten


   write(stderr,*)'Making right hand of normal equations (Ab) and norm btPb..'

   !modify apriori residual
   !solve U**(T)X=residu for X
   !DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
   call DTRSV('U','T','N',ndat,cov(1:ndat,1:ndat),ndat,residu(1:ndat),1)


   !calculate apriori norm (b-Ax0)'P(b-Ax0) in thise case b'P b= (residu' U**(-1))' ( U**(-T) residu)
   LtPL=dot_productblas(residu(1:ndat),residu(1:ndat))

   
   Ab(1:nunk)=matmulblas(A(1:ndat,1:nunk),residu(1:ndat),1)




!write data to binary file


   write(charwk,'(i4.4)')gpswk
   doubles(1)=time
   doubles(2)=time-3.5/365.25
   doubles(3)=time+3.5/365.25
   doubles(5)=LtPL
   ints(1)=ndat
   ints(2)=nunk
   ints(7)=gpswk

   

   call write_BIN(file=trim(outfile)//'.'//charwk,type='SYMV2___',descr=descr,ints_d=ints_d,ints=ints,dbls_d=doubles_d&
        ,dbls=doubles,side1_d=Unlist(1:nunk),vec1=Ab(1:nunk),vec2=apriori(1:nunk),mat1=AtA(1:nunk,1:nunk))

end do  !end loop over files


write(stderr,*)'finished!'


end program GPS_normeq

subroutine help()
implicit none
character(4)::frmt

frmt='(A)'

write(*,frmt)' Program GPS_normeq creates GPS normal equations for Spherical harmonic coefficients'
write(*,frmt)' Files are outputted as unformatted binary file readable by the FORTRAN read_sym_mat2 routine'
write(*,frmt)' default outputfile name is GPSnormeq2.GPSWK where GPSWK is the gpsweek number' 
write(*,frmt)' Observation equation is based on the surface mass loading of on an elastic Earth'
write(*,frmt)' Usage: GPS_normeq [OPTIONS] RESFILES'
write(*,frmt)' Where RESFILES are the residual files outputted from GPS_residual'
write(*,frmt)''
write(*,frmt)' OPTIONS can be the following:'
write(*,frmt)'  -IINCTERM: only include entries which contain the ascii string INCTERM in their SOLUTION/ESTIMATE'
write(*,frmt)'             block (see residual file), more -I options enlarge the valid set'
write(*,frmt)'             (at least one of the search strings is valid)'
write(*,frmt)''
write(*,frmt)'  -EEXCLTERM: exclude entries which contain the ascii term EXCLTERM, more -E options will decrease'
write(*,frmt)'              the valid set'
write(*,frmt)'             (entries which match any of the EXCL terms are discarded)'
write(*,frmt)'  -i INCLUDEFILE: use include entries from file INCLUDEFILE (one entry per line)'
write(*,frmt)''
write(*,frmt)'  -e EXCLUDEFILE: use exclude entries from file EXCLUDEFILE (one entry per line)'
write(*,frmt)''
write(*,frmt)'  -lLMAX[,LMIN]: limit maximum (and possibly minimum) degree of the output to LMAX and LMIN'
write(*,frmt)'                defaults assumes lmax =7 and lmin=1 (or lmin=2 when -g is specified)'
write(*,frmt)''
write(*,frmt)'  -t: Also incorporate a 7 parameter HELMERT transformation in the Observation equation'
write(*,frmt)'  -g: Express degree 1 loading  variations as geocenter motion'
write(*,frmt)'  -s: Construct system with spherical harmonic loading coefficients'
write(*,frmt)'  -m: Construct a mean unknown for each station and direction '
write(*,frmt)'  -mt: Construct mean and linear drift unknowns for each station and direction '
write(*,frmt)'  -a: Construct an annual fit for each station and direction'
write(*,frmt)'  -as: Construct an annual and semiannual fit for each station and direction'
write(*,frmt)'  -M MEANFILE: subtract a mean from the data'
write(*,frmt)'  -d: Only use the diagional of the SINEX files'
write(*,frmt)'  -FBASENAME: use BASENAME as the output file basename (default is GPSnormeq1)'
stop
end subroutine help

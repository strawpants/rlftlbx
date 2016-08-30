!program to convert new style residual file to old style residual files used by the program invert02 from J. Kusche
!Also calculates a multiyearmean for the data

!!Coded by Roelof Rietbroek, Thu Sep 27 13:08:08 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Tue Nov  6 15:53:45 2007
!!now converts newer sinex files (containing position Covariance etc ) outputted from GPS_residual 

!!Updated by Roelof Rietbroek, Thu Mar  6 09:56:59 2008
!! no need for a discontinuity file anymore


program GPS_conv2old
use FORTtlbx
use GPStlbx
implicit none
integer::i,j,itharg,narg,ind,nfile,ndat,stderr,ndis,wk,k
integer,dimension(:),allocatable::aveind
double precision,allocatable, dimension(:)::cv,lon,lat,h,e,n,h_sig,n_sig,e_sig
double precision,allocatable, dimension(:,:)::ave
double precision::year,lats,lons,lontmp,lattmp,pi
integer::tmin,tmax,latd,lond,lonm,latm,ufunit,num,maxstat
logical::mean,new
!logical::cov
character(4)::gpswk
character(13)::mfname
character(120)::dum,infiles(1000)
character(80),dimension(:),allocatable::stattotlist,statlist,hlist,elist,nlist
integer::nstat,nstattot
integer::iargc
!initializations
maxstat=1000
pi=acos(-1.d0)
stderr=0
ufunit=13
itharg=0
!cov=.false.
mean=.false.
tmin=99999!set to large value
tmax=0!set to small value
nfile=0

!get command line arguments

narg=iargc()


if(narg <1)call help()

do,i=1,narg
   if(itharg >= narg)exit !exit loop after last argument
   itharg=itharg+1
   call getarg(itharg,dum)

   if(dum(1:1).eq.'-')then !argument is an option

      select case(dum(2:2))
!         case('c')!input file are covariance matrices
 !           cov=.true.
         case('m') !calculate mean over certain period (discontinuity file must be given!)
            mean=.true.
           !  ind=index(dum(3:),',')
!             if(ind .eq. 0)then
!                tmin=0
!                tmax=2000
!             else
!                read(dum(3:ind-1),*)tmin
!                read(dum(ind+1:),*)tmax
!             end if

!          case('d')!specify discontinuity file
!             if(dum(3:3) .eq. ' ')then !also accept a space after the -d option
!                itharg=itharg+1
!                call getarg(itharg,dum)
!                disfile=trim(dum)
!             else
!                disfile=trim(dum(3:))
!             end if


         case default
            write(stderr,*)'Unknown option ',dum(2:2)
            call help
      end select
   else!argument is an input file
      nfile=nfile+1
      infiles(nfile)=trim(dum)
   end if



end do

!do some input checks
!if(mean .and. cov)then !impossible combination
!   write(stderr,*)' Converting covariance matrices and making a time mean cannot be done at the same time'
!   stop
!end if

if (nfile .eq. 0)then
   write(stderr,*)'No input files given'
   stop
end if

! if(. (index(disfile,' ') .eq. 1))then
   
!    write(stderr,*)'Discontinuity file must be provided'
!    stop
! end if

! !write(*,*)trim(disfile),index(trim(disfile),'')

! !if(.not. cov)then
!    !read in discontinuity file
!    call get_sinexblock(file=trim(disfile),blockname='SITE/ID',nlines=ndis)
!    allocate(discon(ndis))
   
!    call get_sinexblock(file=trim(disfile),blockname='SITE/ID',list=discon)
! !end if

if(mean)then !allocate arrays for  mean values and amount of valid data points
   allocate(ave(maxstat,3),aveind(maxstat), stattotlist(maxstat))
   ave=0.d0
   aveind=0
end if

allocate(statlist(maxstat))
allocate(h(maxstat),e(maxstat),n(maxstat),hlist(maxstat),elist(maxstat),nlist(maxstat))
allocate(h_sig(maxstat),e_sig(maxstat),n_sig(maxstat))
allocate(lon(maxstat),lat(maxstat))
allocate(cv(maxstat*3*(3*maxstat+1)/2))
!loop over input files
nstattot=0

do,i=1,nfile
!if covariance file readin covariance file
!    if(cov)then
!       write(stderr,*)'Converting covariance file.. ',trim(infiles(i))
!       call read_sym_mat(file=trim(infiles(i)),nval=ndat)
!       if(allocated(cv))then
!          deallocate(cv)
!       end if

!       allocate(cv(ndat*(ndat+1)/2))
!       call read_sym_mat(file=trim(infiles(i)),matpack=cv,time=year)

!       !now write entries to unformaTted data file
!       write(gpswk,'(I4.4)')GPS_week(sinex_date(year))
!       open(unit=ufunit,file='igc.'//gpswk,form='unformatted')
!       write(ufunit)cv
!       close(ufunit)

!    else if(mean)then!calculate mean of residual files
   
   if(mean)then!calculate mean of residual files
      !get meta data
      call get_sinex(file=trim(infiles(i)),nx=ndat,inc1='*STAH *',gpswk=wk)
      !get 
      tmin=min(wk,tmin)
      tmax=max(wk,tmax)


      h=0.d0
      e=0.d0
      n=0.d0
     elist=''
     hlist=''
     nlist=''

      !read data in arrays
      !height
      call get_sinex(file=trim(infiles(i)),x=h(1:ndat),inc1='*STAH *',xlist=hlist(1:ndat))
      call get_sinex(file=trim(infiles(i)),x=e(1:ndat),inc1='*STAE *',xlist=elist(1:ndat))
      call get_sinex(file=trim(infiles(i)),x=n(1:ndat),inc1='*STAN *',xlist=nlist(1:ndat))


      !read in station positions
      call get_sinexblock(file=trim(infiles(i)),blockname='SITE/ID',nlines=nstat)

   
    call get_sinexblock(file=trim(infiles(i)),blockname='SITE/ID',list=statlist(1:nstat))


    !add stations to existing list if not found
    
    do,j=1,nstat
       new=.true. !assume a new station untill proven otherwise
       do,k=1,nstattot
          if(statlist(j)(2:8) .eq. stattotlist(k)(2:8))then
             new=.false.
             exit
          end if
       end do
       !make new entry
       if(new)then
          nstattot=nstattot+1
          stattotlist(nstattot)=statlist(j)
       end if
    end do
    
    

    
    

      !now do a loop over all stations in the discontinuity file and retrieve station position
      
      do,j=1,nstattot
         num=0
         do,k=1,ndat
!            if(compstr(hlist(k),'ZWEN'))write(*,*)trim(hlist(k)),trim(elist(k)),trim(nlist(k))
            if(compstr(hlist(k),'*'//stattotlist(j)(2:8)//'*'))then
              
            num=num+1
            ave(j,1)=ave(j,1)+h(k)
            aveind(j)=aveind(j)+1 !increment amount of data points
            end if

            if(compstr(elist(k),'*'//stattotlist(j)(2:8)//'*'))then
            num=num+1
            ave(j,2)=ave(j,2)+e(k)

            end if

            if(compstr(nlist(k),'*'//stattotlist(j)(2:8)//'*'))then
            num=num+1
            ave(j,3)=ave(j,3)+n(k)

            end if

            if(num >= 3)exit !exit loop when all components have been collected

         end do
      end do
      

   else !convert residual files and write position files
      write(stderr,*)'Converting residual file ',trim(infiles(i))

      !get meta data
      call get_sinex(file=trim(infiles(i)),nx=ndat,inc1='*STAH *',gpswk=wk)
            
    
    

      h=0.d0
      e=0.d0
      n=0.d0
      h_sig=0.d0
      e_sig=0.d0
      n_sig=0.d0
      !read data in arrays
      !height
      call get_sinex(file=trim(infiles(i)),x=h(1:ndat),covdiagx=h_sig(1:ndat),inc1='*STAH *',xlist=hlist(1:ndat))
      call get_sinex(file=trim(infiles(i)),x=e(1:ndat),covdiagx=e_sig(1:ndat),inc1='*STAE *')
      call get_sinex(file=trim(infiles(i)),x=n(1:ndat),covdiagx=n_sig(1:ndat),inc1='*STAN *')

      !get covariance matrix

      call get_sinex(file=trim(infiles(i)),covxypack=cv(1:ndat*3*(3*ndat+1)/2),inc1='*STA*',inc2='*STA*')

      !write covariance matrix
      write(mfname,'(A4,I4.4)')'igc.',wk
      open(unit=ufunit,file=trim(mfname),form='unformatted')
      write(ufunit)cv(1:ndat*3*(3*ndat+1)/2)
      close(ufunit)


      !retrieve longitude and latitude from SITE/ID block

      !read in station positions
      call get_sinexblock(file=trim(infiles(i)),blockname='SITE/ID',nlines=nstat)

   
      call get_sinexblock(file=trim(infiles(i)),blockname='SITE/ID',list=statlist(1:nstat))

!      write(*,*)nstat,ndat, size(lon,1),size(statlist,1)
      
      do,j=1,nstat
         do,k=1,ndat
            if(compstr(hlist(k),'*'//statlist(j)(2:8)//'*'))then

               read(statlist(j)(44:67),'(i4,1x,i2,1x,f4.1,1x,i3,1x,i2,1x,f4.1)')lond,lonm,lons,latd,latm,lats
 !              write(*,*)'start'//statlist(j)(45:67)
!               write(*,*)statlist(j)(45:47),statlist(j)(49:51),statlist(j)(53:56),statlist(j)(58:60)

               lontmp=dble(lond)+sign(lonm,lond)/60.d0
               lon(k)=(lontmp+sign(lons,lontmp)/3600.d0)
               !      lon(i)=(dble(lond)+lonm/60.d0+lons/3600.d0)*pi/180.d0
               lattmp=dble(latd)+sign(latm,latd)/60.d0
               lat(k)=(lattmp+sign(lats,lattmp)/3600.d0)
!               write(*,*)lon(k),lat(k)
               exit
            end if
         end do
      end do


!       do,j=1,ndiscon
!          do,k=1,ndat
!             write(*,*)discon(j)(2:8),hlist(k)
!             if(compstr(hlist(k),'*'//discon(j)(2:8)//'*'))then
!                read(discon(j)(45:67),'(i4,1x,i2,1x,f4.1,1x,i3,1x,i2,1x,f4.1)')lond,lonm,lons,latd,latm,lats
!                lontmp=dble(lond)+sign(lonm,lond)/60.d0
!                lon(k)=(lontmp+sign(lons,lontmp)/3600.d0)*pi/180.d0
!                !      lon(i)=(dble(lond)+lonm/60.d0+lons/3600.d0)*pi/180.d0
!                lattmp=dble(latd)+sign(latm,latd)/60.d0
!                lat(k)=(lattmp+sign(lats,lattmp)/3600.d0)*pi/180.d0
               
!                !        lon(k)=(dble(lond)+lonm/60.d0+lons/3600.d0)
!                !              lat(k)=(dble(latd)+latm/60.d0+lats/3600.d0)
!                !               write(*,*)lond,lonm,lons,latd,latm,lats, lon(k),lat(k)
               
!                exit
!             end if

!          end do
!       end do

      !write data to file
      !residual file
      write(mfname,'(A4,I4.4)')'igs.',wk
      open(unit=ufunit,file=trim(mfname))
      
      do,j=1,ndat
         write(ufunit,'(1x,A4,1x,6E13.5)')hlist(j)(15:18),h(j),e(j),n(j),h_sig(j),e_sig(j),n_sig(j)
       !  write(*,'(1x,A4,1x,6E13.5)')hlist(j)(15:18),h(j),e(j),n(j),h_sig(j),e_sig(j),n_sig(j)
      end do
      close(ufunit)

!position file (just a copy of site/ID block data)
      write(mfname,'(A4,I4.4)')'igp.',wk
      open(unit=ufunit,file=trim(mfname))
      
      do,j=1,ndat
         write(ufunit,'(2F12.5,1x,A4)')lon(j),lat(j),hlist(j)(15:18)
      end do
      close(ufunit)

   end if
end do

!write mean file
if(mean)then
   !average the sum of the residuals
   
   do,i=1,nstattot

   !   write(*,*)i,aveind(i),stattotlist(i)(2:8)
      ave(i,:)=ave(i,:)/dble(aveind(i))
   end do
!write(*,*)nstattot,size(stattotlist)

!write mean data to file
   write(mfname,'(A4,I4.4,A1,I4.4)')'igs.',tmin,'-',tmax
   open(unit=ufunit,file=trim(mfname))
   
   do,i=1,nstattot
      if(aveind(i) .eq. 0)cycle !ignore stations which had no entries
   !   write(*,*)discon(i)(2:5),ave(i,1),ave(i,2),ave(i,3),aveind(i)
      write(ufunit,*)stattotlist(i)(2:5),ave(i,1),ave(i,2),ave(i,3)
   end do
   
   close(ufunit)
end if



end program GPS_conv2old

subroutine help()
character(3)::frmt
frmt='(A)'
write(*,frmt)'Program GPS_conv2old converts residuals in new format to old format (required by Kusches'
write(*,frmt)' inversion software)'
write(*,frmt)'Usage: GPS_conv2old [OPTIONS] INPUTFILES'
write(*,frmt)'Where INPUTFILES are residual files'
write(*,frmt)'OPTIONS:'
!write(*,frmt)'  -c: inputfiles are covriance matrices' !depreciated
write(*,frmt)'  -m: calculate the mean of the files'
!write(*,frmt)'  -d DISCON: specify discontinuity file'
stop
end subroutine help

!fortran program which makes residuals from sinex files and refers them against a cumulative file or reference file
!!Coded by Roelof Rietbroek, Wed Aug 22 09:04:06 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Mon Oct 22 17:10:39 2007
!!now writes the covariance information also in the sinex file instead of in a seperate binary file

!!Updated by Roelof Rietbroek, Fri Nov  9 16:27:23 2007
!!fixed bug: writes variance instead of standard deviations to solution/estimate block

program GPS_residual
use GPStlbx
use FORTtlbx
implicit none
integer::stderr,narg,itharg,i,j,k,maxf,maxdat,minwks
parameter(maxf=500)!maximum amount of file which program can handle
parameter(maxdat=2000)!maximum amount of data points in a weekly igs file
character(120)::outfile1,outfile2,cumfile,disfile,dum
character(120),dimension(:)::sinexfile(maxf)
logical::itrf,lim
double precision,dimension(:)::sindat(maxdat),refdat(maxdat)
double precision,dimension(:,:)::cov(maxdat,maxdat)
double precision,allocatable,dimension(:,:)::rot
double precision,allocatable,dimension(:)::postmp
integer::validstat(maxdat),funit,f,error,nfile,ind,ndat,nval,wk,refwk1,nsiteid
integer,allocatable,dimension(:)::ivec
character(120),dimension(:)::list(maxdat)
character(7),allocatable,dimension(:)::typelist
character(80)::desc,dum2
character(52)::hdrline
character(5)::a
character(80),allocatable,dimension(:)::ref,siteid
character(12)::srchstrng
double precision::pos,vel,reftime
logical::velocity,rotate
integer::iargc
!defaults/initializations
cumfile=' '
disfile=' '
sinexfile=' '
stderr=0
velocity=.true.
nfile=0
itharg=0
desc='Covariance matrix of igs file rotated in local frame'
dum2=''
outfile1='residual'
outfile2='cov'
minwks=0
nfile=0
refwk1=0
rotate=.true.
itrf=.false.
!first get command line arguments

narg=iargc()
if(narg <1)call help()

!process command line arguments


do,i=1,narg
   itharg=itharg+1
  ! write(*,*)itharg,narg
   if(itharg > narg)exit !exit loop when last argument has been read in
!write(*,*)'before',itharg,narg,trim(dum)
   call getarg(itharg,dum)
!write(*,*)'after',itharg,narg,trim(dum)
   if(dum(1:1).eq. '-') then !process when argunment is an option
      select case(dum(2:2))
      case('d')!specify discontinuity file
         if(dum(3:3) .eq. ' ')then !also accept a space after the -d option
            itharg=itharg+1
            call getarg(itharg,dum)
            disfile=dum
         else
            disfile=dum(3:)
         end if
      case('c')!specify cumulative file
         if(dum(3:3) .eq. ' ')then !also accept a space after the -d option
            itharg=itharg+1
            call getarg(itharg,cumfile)
         else
            cumfile=dum(3:)
         end if
     case('r')!specify reference GPS week of cumulative file
        read(dum(3:),*)refwk1
      case('m')!consider mean only no velocity correction
         velocity=.false.

      case('i')!refer all measurements to ITRF2000
         itrf=.true.
      case('l') !Limit solution epoch to minimum amount of weeks'
         lim=.true.
         read(dum(3:),*)minwks
      case('x')!don't rotate measurements to local frame
         rotate=.false.
      case default
         write(stderr,*)'unknown option: ',dum(2:2)
         call help()
      end select
   else !argument is an input sinex file
      nfile=nfile+1
      
      sinexfile(nfile)=trim(dum)
     
   end if
end do

!check if cumulative GPS reference week is given
if(refwk1 .eq. 0 .and. itrf )then
    write(stderr,*)'reference week of cumulative file must be specified'
    write(stderr,*)'use -rWK'
    stop
 end if


do,f=1,nfile!loop over files
   !load weekly sinex file in array
   sindat=0.d0
   refdat=0.d0
   cov=0.d0
  
   list=' '
   
   call get_sinex(trim(sinexfile(f)),nx=ndat,x=sindat,inc1='* STA*',inc2='* STA*',covxy=cov,xlist=list,gpswk=wk)
   validstat=0
   

   !loop over station positions contained in the sinexfile to retrieve reference positions from cumulative file
   
   do,i=1,ndat
      
      !skip stations which have more than one solution in one week
      call get_sinex(trim(sinexfile(f)),nx=error,inc1='*'//list(i)(15:22)//'*')


      select case(error)
         case(:2)
            write(stderr,*)'WARNING: missing X,Y,Z components in the weekly sinex'
            cycle
         case(6:)
            write(stderr,*)'WARNING: ignoring station with multiple entries',list(i)(15:22)
            cycle
      end select


      !get solution id string to search for in cumulative file (if solution id is given as ----
      if((list(i)(23:26) .ne. '----'))then
         
         srchstrng=list(i)(15:22)//'   '//list(i)(26:26)

      else
         srchstrng=get_solu_id(file=trim(disfile),station=list(i)(15:21),gpswk=wk,minwks=minwks)
      end if
 
     
      if(srchstrng .eq.'############')cycle !not a valid station or too short time interval
      !get x,y or z position and x y or z velocity from cumulative file

      if(compstr(list(i),'*STAX*'))call get_cumdata(file=trim(cumfile),inc='*X*'//srchstrng//'*',pos=pos,&
           vel=vel,reftime=reftime,error=error)
      if(compstr(list(i),'*STAY*'))call get_cumdata(file=trim(cumfile),inc='*Y*'//srchstrng//'*',pos=pos,&
           vel=vel,reftime=reftime,error=error)
      if(compstr(list(i),'*STAZ*'))call get_cumdata(file=trim(cumfile),inc='*Z*'//srchstrng//'*',pos=pos,&
           vel=vel,reftime=reftime,error=error)
      !error handling
      if(error .eq. -1)then
         write(stderr,*)'WARNING: ',srchstrng,': has not all three X,Y,Z components'
         cycle !go to next station when entry was not found and ignore
      else if(error .eq. -2)then
          write(stderr,*)'WARNING: ',srchstrng,': is not found in cumulative file'
          !!!
          !now try to search for a ---- solution id (older cumulative files)
          !get x,y or z position and x y or z velocity from cumulative file

          if(compstr(list(i),'*STAX*'))call get_cumdata(file=trim(cumfile),inc='*X*'//srchstrng(1:8)//'----*',pos=pos,&
               vel=vel,reftime=reftime,error=error)
          if(compstr(list(i),'*STAY*'))call get_cumdata(file=trim(cumfile),inc='*Y*'//srchstrng(1:8)//'----*',pos=pos,&
               vel=vel,reftime=reftime,error=error)
          if(compstr(list(i),'*STAZ*'))call get_cumdata(file=trim(cumfile),inc='*Z*'//srchstrng(1:8)//'----*',pos=pos,&
               vel=vel,reftime=reftime,error=error)

         !  if(error .eq. -2)then !third option just use what we got 
!              if(compstr(list(i),'*STAX*'))call get_cumdata(file=trim(cumfile),inc='*X*'//srchstrng(1:8)//'*',pos=pos,&
!                   vel=vel,reftime=reftime,error=error)
!              if(compstr(list(i),'*STAY*'))call get_cumdata(file=trim(cumfile),inc='*Y*'//srchstrng(1:8)//'*',pos=pos,&
!                   vel=vel,reftime=reftime,error=error)
!              if(compstr(list(i),'*STAZ*'))call get_cumdata(file=trim(cumfile),inc='*Z*'//srchstrng(1:8)//'*',pos=pos,&
!                   vel=vel,reftime=reftime,error=error)
!           end if

          if(error .ne. 0)cycle !still not found => cycle loop
          !give warning message
          write(stderr,*)'using ',srchstrng(1:8),' instead'
        
      end if
      validstat(i)=1 !set flag to 1 for valid station
     
      if(velocity)then
         refdat(i)=pos+vel*(GPS_year(wk)-reftime)
      else
         refdat(i)=pos
      end if
     
!     if(compstr(list(i),'*ARTU*'))write(*,*)pos,vel
   end do
   
   !now make an index vector for the valid stations
  ! write(*,*)'test'
   if(allocated(ivec))deallocate(ivec)
   !amount of valid stations
   nval=sum(validstat)
   allocate(ivec(nval))
   ind=0
   do,i=1,ndat
      if(validstat(i).eq.0)cycle
      ind=ind+1
      ivec(ind)=i
   end do

!adjust to common itrf2000 frame if required 
   if(itrf)then
      allocate(postmp(nval))

     !call conv2itrf2000(typ=list(ivec)(7:13),station=list(ivec)(15:21),pos=refdat(ivec),gpswk=wk)
!      write(*,*)'before',list(ivec(20)),refdat(ivec(20)),20,ivec(20),nval
      call conv2itrf2000(typ=list(ivec)(7:13),station=list(ivec)(15:21),pos=refdat(ivec),posnew=postmp,gpswk=refwk1,refwk=wk)
    
     refdat(ivec)=postmp
 !    write(*,*)'after',list(ivec(20)),refdat(ivec(20)),20,ivec(20),nval

!         call  conv2itrf2000(typ=list(ivec)(7:13),station=list(ivec)(15:21),pos=refdat(ivec),posnew=postmp,gpswk=wk)
     call  conv2itrf2000(typ=list(ivec)(7:13),station=list(ivec)(15:21),pos=sindat(ivec),posnew=postmp,gpswk=wk,refwk=wk)
     !refdat(ivec)=postmp
     sindat(ivec)=postmp
     deallocate(postmp)
  !    write(*,*)'after',sindat(ivec(1:3)),refdat(ivec(1:3))
   end if

   
   !rotate data to local frame 
   if(rotate)then
      !get rotation matrix
      if(allocated(rot))deallocate(rot)
      allocate(rot(nval,nval))
      
      !use a temporary variable to correctly represent an intent(inout) attribute 
      if(allocated(typelist))deallocate(typelist)
      allocate(typelist(nval))
      typelist=list(ivec)(7:13)
      call GPS_rotmat(type=typelist,station=list(ivec)(15:21),pos=refdat(ivec),rot=rot)
      !copy back
      list(ivec)(7:13)=typelist

!the construction below would have an unchanged list(ivec)(7:13)
!      call GPS_rotmat(type=list(ivec)(7:13),station=list(ivec)(15:21),pos=refdat(ivec),rot=rot)
      
      !rotate data
      sindat(ivec)=matmulblas(rot,(sindat(ivec)-refdat(ivec)))
       
      !propagate covariance matrix through rotation matrix
      cov(ivec,ivec)=matmulblas(rot,matmulblas(cov(ivec,ivec),rot,0,1))
      !   sindat(ivec)=sindat(ivec)-refdat(ivec)
   else !don't rotate data just subtract
      sindat(ivec)=sindat(ivec)-refdat(ivec)
   end if


   !make character array to write to file
   if(allocated(ref))deallocate(ref)
   allocate(ref(nval))
   do,i=1,nval
       write(ref(i),'(1x,i5,a41,E21.15,1x,e11.6)')i,list(ivec(i))(7:47),sindat(ivec(i)),sqrt(cov(ivec(i),ivec(i)))
      !write(*,'(1x,i5,a41,E21.15,1x,e11.6)')i,list(ivec(i))(7:47),sindat(ivec(i)),cov(ivec(i),ivec(i))
   end do
   !write data to files

   !copy the site data block of the sinex file
   call get_sinexblock(file=trim(sinexfile(f)),blockname='SITE/ID',nlines=nsiteid)
   if(allocated(siteid))deallocate(siteid)
   allocate(siteid(nsiteid))
   call get_sinexblock(file=trim(sinexfile(f)),blockname='SITE/ID',list=siteid)
   
   !make character array from gpsweek number
   write(a,'(a1,i4.4)')'.',wk
   !make header line for file
   write(hdrline,'(a32,i5,a7)')'GFZ '//sinex_date(GPS_year(wk)-(3.5/365.25))//' '//sinex_date(GPS_year(wk)+(3.5/365.25))//' P '&
        ,nval,' 2  X  '
   !initialize file
   call init_sinex(file=trim(outfile1)//a,unit=funit,line=hdrline)
   !write data to file
   !siteid block
   call write_sinexblock(unit=funit,sblock='SITE/ID',block=siteid,&
        comm='*Code Pt __Domes__ T _Station Description__ _Longitude_ _Latitude__ _Height')

   !write solution estimate block
   call write_sinexblock(unit=funit,sblock='SOLUTION/ESTIMATE',block=ref,&
        comm='*Index _Type_ Code Pt Soln _Ref_Epoch__ Unit S __Estimated Value____ _Std_Dev___')
   !write covariance matrix
   call write_sinex_MATlow(unit=funit,mat=cov(ivec,ivec),type='COVA')
   call end_sinex(funit)
   
   !write covariance data to unformatted file
   !call write_sym_mat(file=trim(outfile2)//a,mat=cov(ivec,ivec),list=list(ivec)(1:50)&
    !    ,time=GPS_year(wk),descr='rotated IGS-GPS COVARIANCE matrix')
   

end do !end loop over files


end program GPS_residual

!prints help message
subroutine help()
implicit none
character(4)::frmt

frmt='(A)'
write(*,frmt)'Usage: GPS_residual [OPTIONS]  SINEXFILE(S)'
write(*,frmt)'Calculates residuals with respect to GPS cumulative solution'
write(*,frmt)'Residuals are rotated in height east and north frame'
!write(*,frmt)'Outputs for every sinexfile an SINEX block with residuals and a unformatted covariance matrix'
write(*,frmt)'Outputs for every sinexfile a modified SINEX file with residuals and their covariance '
write(*,frmt)'(and possibly rotated to H,E,N)'
write(*,frmt)''
write(*,frmt)'OPTIONS are:'
write(*,frmt)'-d IGSDISCONTINUITTYFILE: specify the name of the IGS discontinuity file (required)'
write(*,frmt)'-c IGSCUMULATIVEFILE: specify the name of the IGS cumulative file (required)'
write(*,frmt)'-m : Do not use velocities in the cumulative file'
write(*,frmt)'-rWK:Use WK as a ITRF reference week for the cumulative file'
write(*,frmt)'-i : Convert residuals to the ITRF2000 frame'
write(*,frmt)'-lWKS: Only use data which have a continuous solution over amount at least the amount of weeks: WKS'
write(*,frmt)'-x: do not rotate data to local frame'

end subroutine help


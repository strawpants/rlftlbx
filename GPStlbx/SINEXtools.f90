!!Coded by Roelof Rietbroek, Fri Apr 10 13:56:53 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!! A collection of routines to read in SINEX files
!!design philosophy:
!!1) data must be read chronologically (NO REWIND or REOPENING of the file is allowed)
!! this must obeyed in order to read from pipes/standard input
!!2) The SINEX data is read in a derived type specifically made for SINEX
!!   Fields are pointers and may therefore be automatically 
!!   allocated within this routine ( or outside if preferred) 
!!  This allows great flexibility (read: easy updating) in this module and the programs which are calling it

!!TODO: basically a lot
!! writing of SINEX files?
!! extend SINEX derived structure to conform the full SINEX format


Module SINEXtools
  implicit none
  !version number of this tool
  character(10)::version='BETA0.1'

  !!define the SINEX derived type
  type SNXdat
     character(200)::file='stdin' ! SINEX filename ( defaults read from standard input)
     double precision::vers ! sinex version of the file
     character(12)::created,tstart,tend ! creation, Observation start end end time in sinex format
     character(1)::obsc,constr ! single character holding global observation code and constraint parameter
     character(1)::solcont(6) ! solution content

     !parameter space, stations and solution ID
     character(6),pointer,dimension(:):: para =>null()
     character(4),pointer,dimension(:):: site =>null()
     character(2),pointer,dimension(:):: sitec =>null()
     character(4),pointer,dimension(:):: soluid =>null()

     !siteid part
     integer::nstat=0
     character(19),pointer,dimension(:)::domes=>null()
     character(22),pointer,dimension(:)::stat_d=>null()
     double precision,pointer,dimension(:)::lon=>null()
     double precision,pointer,dimension(:)::lat=>null()
     double precision,pointer,dimension(:)::height=>null()
     
     !solution part
     double precision,pointer,dimension(:)::Sv=>null() ! solution vector
     double precision,pointer,dimension(:)::Svsig=>null() ! sigma solution vector
     character(1),pointer,dimension(:)::Scon=>null() ! constaint code
     double precision,pointer,dimension(:,:)::Sm=>null() ! Matrix 
     character(1)::Smform
     character(4)::Smtype
     integer::Smflag=0
     integer::Svflag=0

     !normal matrix part
     double precision,pointer,dimension(:,:)::Nm=>null() ! Matrix
     double precision,pointer,dimension(:)::Nv=>null() ! right hand side vector of normal system
     character(1),pointer,dimension(:)::Ncon=>null() ! constaint code
     integer::Nmflag=0    
     integer::Nvflag=0
     character(1)::Nmform

     !solution/statistics part
     integer::nval ! amount of unknowns in the SINEX file
     integer::nunk ! original amount of unknowns
     integer::nobs ! amount of used observations
     integer::dof ! degrees of freedom (= nobs-nval)
     double precision::sig0_c=1.d0 ! variance factor ( defualts to 1)
     double precision::vtpv,ltpl ! quadratic cost functional ( posteriori/apriori)
     integer::SSflag=0

     !apriori part
     double precision,pointer,dimension(:,:)::Am=>null() ! Matrix
     double precision,pointer,dimension(:)::Av=>null() ! right hand side vector of normal system
     integer::Amflag=0    
     integer::Avflag=0
     character(4)::Amtype
     character(1)::Amform
     
  end type SNXdat

  contains !! below are the associated routines
    subroutine read_SNX(dat)
      use forttlbx
      implicit none
      type(SNXdat)::dat
      integer::unit,stderr,row,col,last,lastb
      logical::fex,warn
      character(80)::dum! covers exactly 1 full line
      double precision::dumvec(3)
      character(120)::form
      integer::chunk,szalloc,lon,lonm,lat,latm
      double precision::lons,lats

      !defaults
      unit=5 ! read from standard input by default
      stderr=0
      warn=.true.
      chunk=200

      !open file
      if(dat%file .ne. 'stdin')then
         unit=13
         !inquire filename
         inquire(FILE=trim(dat%file),EXIST=fex)
         if( .not. fex)then
            write(stderr,*)'ERROR:read_SNX File ',trim(dat%file),' does not exist, quitting'
            stop
         end if
         open(unit=unit,file=dat%file,form='formatted',status='old')
      end if
      
      !read first line
      form='(a5,1x,f4.2,5x,a12,5x,a12,1x,a12,1x,a1,1x,I5.5,1x,a1,1x,a6,6(1x,a1))'
      read(unit=unit,fmt=form,iostat=last)dum(1:5),dat%vers,dat%created,dat%tstart&
           ,dat%tend,dat%obsc,dat%nval,dat%constr,dat%solcont
      if (dum(1:5) .ne. '%=SNX')then
         write(stderr,*)"ERROR:read_SNX:",trim(dat%file),"does not appear to be a sinex file"
         stop
      end if

      ! read until end of file occurrs
      last=0
      do while(last .eq. 0) ! loop to search for blocks
         read(unit=unit,fmt='(A1)',iostat=last,advance='NO')dum(1:1)!read only one byte to speed up scanning for blocks

         if(dum(1:1) .ne. '+')then ! only pass when the beginning of a block is encountered
            read(unit=unit,fmt='(1x)',iostat=last)
            cycle
         end if
         
         !read remainder of the line
         read(unit=unit,fmt='(A79)',iostat=last)dum(2:)
         !select the appropriate block ( unsupported ones will be ignored)

        ! write(0,*)dum(2:)
         select case(dum(2:))
         case('SITE/ID')
            last=0
            dat%nstat=0
            szalloc=chunk
            if(.not. associated(dat%stat_d))allocate(dat%stat_d(chunk))
            if(.not. associated(dat%domes))allocate(dat%domes(chunk))
            if(.not. associated(dat%lon))allocate(dat%lon(chunk))
            if(.not. associated(dat%lat))allocate(dat%lat(chunk))
            if(.not. associated(dat%height))allocate(dat%height(chunk))
            
            do while(last .eq. 0)
               read(unit=unit,fmt='(A80)',iostat=last)dum
               if(last .ne.  0 .or. dum(1:1) .eq. '-') exit !end of file/block is reached
               if(dum(1:1).eq. '*')cycle !skip comments
               dat%nstat=dat%nstat+1
               if(dat%nstat > szalloc)then ! reallocate data
                  call realloc_ptr(dat%domes,chunk)
                  call realloc_ptr(dat%stat_d,chunk)
                  call realloc_ptr(dat%lon,chunk)
                  call realloc_ptr(dat%lat,chunk)
                  call realloc_ptr(dat%height,chunk)
                  szalloc=szalloc+chunk
               end if
               
               !extract info from string
               dat%domes(dat%nstat)=dum(2:20)
               dat%stat_d(dat%nstat)=dum(22:43)
! AIS1  A 49998S001 P AIS1 49998S001         228 24  1.7  55  4  8.7    32.2
               read(dum(44:),'(2(1x,i3,1x,i2,1x,F4.1),1x,F7.1)')&
                    lon,lonm,lons,lat,latm,lats,dat%height(dat%nstat)
               dat%lon(dat%nstat)=lon+sign(lonm,lon)/60.d0+sign(lons,dble(lon))/3600.d0
               dat%lat(dat%nstat)=lat+sign(latm,lat)/60.d0+sign(lats,dble(lat))/3600.d0
               
            end do
         case('SOLUTION/ESTIMATE') ! solution estimate block
!               write(0,*)'SOLUTION/ESTIMATE'
               !allocate vectors
               if(.not. associated(dat%Sv))allocate(dat%Sv(dat%nval))
               if(.not. associated(dat%Svsig))allocate(dat%Svsig(dat%nval))
               if(.not. associated(dat%para))allocate(dat%para(dat%nval))
               if(.not. associated(dat%Scon))allocate(dat%Scon(dat%nval))
               if(.not. associated(dat%soluid))allocate(dat%soluid(dat%nval))
               if(.not. associated(dat%sitec))allocate(dat%sitec(dat%nval))
               if(.not. associated(dat%site))allocate(dat%site(dat%nval))
               form='(i5,1x,a6,1x,a4,1x,a2,1x,a4,1x,a12,1x,a4,1x,a1,1x,E21.15,1x,E11.6)'
!                form='(i5,1x,a6,1x,a4,1x,a2,1x,a4,1x,a12,1x,a4,1x,a1,1x,E21.15,1x,E11.6)'
               lastb=0
               do while (lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%Svflag=1 ! set flag to read in and exit loop
                     exit
                  end if

                  !read data
                     read(unit=unit,fmt=form,iostat=lastb)row,dat%para(row),dat%site(row),&
                       dat%sitec(row),dat%soluid(row),dum(1:12),&
                       dum(14:17),dat%Scon(row),dat%Sv(row),dat%Svsig(row)     
!                      read(unit=unit,fmt=form,iostat=lastb)row,dat%para(row),dat%site(row),&
!                        dat%sitec(row),dat%soluid(row),dum(1:12),&
!                        dum(14:17),dat%Scon(row),dat%Sv(row),dat%Svsig(row)     
!                     write(*,*)dat%para(row),dat%Sv(row)
               end do
            case('SOLUTION/NORMAL_EQUATION_VECTOR') ! normal equation vector
!               write(0,*)'SOLUTION/NORMAL_EQUATION_VECTOR'
               !allocate vectors
               if(.not. associated(dat%Nv))allocate(dat%Nv(dat%nval))
               if(.not. associated(dat%para))allocate(dat%para(dat%nval))
               if(.not. associated(dat%Ncon))allocate(dat%Ncon(dat%nval))
               if(.not. associated(dat%soluid))allocate(dat%soluid(dat%nval))
               if(.not. associated(dat%sitec))allocate(dat%sitec(dat%nval))
               if(.not. associated(dat%site))allocate(dat%site(dat%nval))
               form='(i5,1x,a6,1x,a4,1x,a2,1x,a4,1x,a12,1x,a4,1x,a1,1x,E21.15)'
               lastb=0
               do while (lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%Nvflag=1 ! set flag to read in and exit loop
                     exit
                  end if
                  !read data
                  read(unit,form)row,dat%para(row),dat%site(row),&
                       dat%sitec(row),dat%soluid(row),dum(1:12),&
                       dum(14:17),dat%Ncon(row),dat%Nv(row)
               end do
            case('SOLUTION/APRIORI')
!                write(0,*)'SOLUTION/APRIORI'
               !allocate vectors
               if(.not. associated(dat%Av))allocate(dat%Av(dat%nval))

               form='(i5,41x,E21.15)'
!                form='(i5,1x,a6,1x,a4,1x,a2,1x,a4,1x,a12,1x,a4,1x,a1,1x,E21.15)'
               lastb=0
               do while (lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%Avflag=1 ! set flag to read in and exit loop
                     exit
                  end if
                  !read data
!                   read(unit,form)row,dat%para(row),dat%site(row),&
!                        dat%sitec(row),dat%soluid(row),dum(1:12),&
!                        dum(14:17),dat%Acon,dat%Av(row)
                  
                   read(unit,form)row,dat%Av(row)
!                   write(0,*)row,dat%Av(row)
               end do
            case('SOLUTION/STATISTICS')! some important statistical values
  !             write(0,*)'SOLUTION/STATISTICS'
               lastb=0
               form='(a79)'
               do while(lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%SSflag=1 ! set flag to read in and exit loop
                     exit
                  end if
                  !read description
                  read(unit=unit,fmt=form,iostat=lastb)dum
 !                 write(0,*)dum(1:30),dumvec(1)
                  select case(dum(1:30))
                     case('NUMBER OF OBSERVATIONS')
                        read(dum(31:),*)dat%nobs
                     case('NUMBER OF UNKNOWNS')
                        read(dum(31:),*)dat%nunk
                     case('SQUARE SUM OF RESIDUALS (VTPV)')
                        read(dum(31:),*)dat%vtpv
                     case('NUMBER OF DEGREES OF FREEDOM')
                        read(dum(31:),*)dat%dof
                     case('VARIANCE FACTOR')
                        read(dum(31:),*)dat%sig0_c

                     case('WEIGHTED SQUARE SUM OF O-C')
                        read(dum(31:),*)dat%ltpl
!                      case('PHASE MEASUREMENTS SIGMA')
!                      case('CODE MEASUREMENTS SIGMA')
!                      case('SAMPLING INTERVAL (SECONDS)')
                  end select
                  
               end do
            case('SOLUTION/MATRIX_ESTIMATE L CORR',&
                 'SOLUTION/MATRIX_ESTIMATE U CORR',&
                 'SOLUTION/MATRIX_ESTIMATE L COVA',&
                 'SOLUTION/MATRIX_ESTIMATE U COVA',&
                 'SOLUTION/MATRIX_ESTIMATE L INFO',&
                 'SOLUTION/MATRIX_ESTIMATE U INFO') ! CORR INFO or COVA solution matrix
 !              write(0,*)'SOLUTION/MATRIX_ESTIMATE'
               dat%Smform=dum(27:27) ! set matrix form ( Upper or lower)
               dat%Smtype=dum(29:32) ! set content COVA/CORR or INFO

               !allocate memory in packed matrix
               if(.not. associated(dat%Sm))allocate(dat%Sm(dat%nval,dat%nval))
               
               
               
               lastb=0
               form='(i5,1x,i5,1x,E21.14,1x,E21.14,1x,E21.14)'
               do while (lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%Smflag=1 ! set flag to read in and exit loop
                     exit
                  end if
                  !read data
                  read(unit,form)row,col,dumvec
                  dat%Sm(row,col)=dumvec(1)
                  if(col<=dat%nval-2)then
                     dat%Sm(row,col+1)=dumvec(2)
                     dat%Sm(row,col+2)=dumvec(3)
                  else if(col<=dat%nval-1)then
                     dat%Sm(row,col+1)=dumvec(2)
                  end if


               end do
               
            case('SOLUTION/NORMAL_EQUATION_MATRIX L'&
                 ,'SOLUTION/NORMAL_EQUATION_MATRIX U'&
                 ,'SOLUTION/NORMAL_EQUATION_MATRIX L INFO'&
                 ,'SOLUTION/NORMAL_EQUATION_MATRIX U INFO') ! unconstrained normal matrix
                  
               dat%Nmform=dum(34:34) ! upper or lower
               
               !allocate memory in packed matrix
               if(.not. associated(dat%Nm))allocate(dat%Nm(dat%nval,dat%nval))
               
               !read data ( 1 byte) to check for comments or end of block
               
               lastb=0
               form='(i5,1x,i5,1x,E21.14,1x,E21.14,1x,E21.14)'
               do while (lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%Nmflag=1 ! set flag to read in and exit loop
                     exit
                  end if

                  read(unit,form)row,col,dumvec
                  dat%Nm(row,col)=dumvec(1)
                  if(col<=dat%nval-2)then
                     dat%Nm(row,col+1)=dumvec(2)
                     dat%Nm(row,col+2)=dumvec(3)
                  else if(col<=dat%nval-1)then
                     dat%Nm(row,col+1)=dumvec(2)
                  end if

                  
!                   if(dat%Nmform .eq. 'U')then !upper triangular ( put in upper triangular packed matrix)
!                      dat%Nm(row+((col-1)*(col))/2)=dumvec(1)
!                      if(col <= (row-2))then
!                         dat%Nm(row+((col)*(col+1))/2)=dumvec(2)
!                         dat%Nm(row+((col+1)*(col+2))/2)=dumvec(3)
!                      else if(col <= (row-1))then
!                         dat%Nm(row+((col)*(col+1))/2)=dumvec(2)
!                      end if
!                   else !lower triangular
!                      dat%Nm(col+((row-1)*(row))/2)=dumvec(1)
!                      if(col <= (row-2))then
!                         dat%Nm(col+1+((row-1)*(row))/2)=dumvec(2)
!                         dat%Nm(col+2+((row-1)*(row))/2)=dumvec(3)
!                      else if(col <= (row-1))then
!                         dat%Nm(col+1+((row-1)*(row))/2)=dumvec(2)
!                      end if
!                   end if
               end do
            case('SOLUTION/MATRIX_APRIORI L CORR',&
                 'SOLUTION/MATRIX_APRIORI U CORR',&
                 'SOLUTION/MATRIX_APRIORI L COVA',&
                 'SOLUTION/MATRIX_APRIORI U COVA',&
                 'SOLUTION/MATRIX_APRIORI L INFO',&
                 'SOLUTION/MATRIX_APRIORI U INFO') ! CORR INFO or COVA solution matrix
               dat%Amform=dum(26:26) ! set matrix form ( Upper or lower)
               dat%Amtype=dum(28:31) ! set content COVA/CORR or INFO
               !allocate memory in packed matrix
               if(.not. associated(dat%Am))allocate(dat%Am(dat%nval,dat%nval))
               dat%Am=0.d0
               
               
               lastb=0
               form='(i5,1x,i5,1x,E21.14,1x,E21.14,1x,E21.14)'
               do while (lastb .eq. 0)
                  !read data ( 1 byte) to check for comments or end of block
                  read(unit=unit,fmt='(A1)',iostat=lastb,advance='NO')dum(1:1)
                  if(dum(1:1) .eq. '*')then
                     read(unit=unit,fmt='(1x)')
                     cycle ! ignore comments
                  else if(dum(1:1).eq. '-')then ! end of block encountered
                     read(unit=unit,fmt='(1x)')
                     dat%Amflag=1 ! set flag to read in and exit loop
                     exit
                  end if
                  !read data
                  read(unit,form)row,col,dumvec
                  dat%Am(row,col)=dumvec(1)
                  if(col<=dat%nval-2)then
                     dat%Am(row,col+1)=dumvec(2)
                     dat%Am(row,col+2)=dumvec(3)
                  else if(col<=dat%nval-1)then
                     dat%Am(row,col+1)=dumvec(2)
                  end if
               end do

            case default ! issue a warning
               if(warn)write(stderr,*)"WARNING: read_SNX, Block ignored:"&
                    ,trim(dum(2:))
         end select
         
      end do

      if(dat%file .ne. 'stdin')close(unit) ! close file if not standard input
    end subroutine read_SNX
end Module SINEXtools

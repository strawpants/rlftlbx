!!program to convert GRACE normal systems to BINtype normal system
!!adapted from GRACE_normeq2
!!Coded by Roelof Rietbroek, Tue Feb 10 21:05:45 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!Updated by Roelof Rietbroek, Mon Oct 19 16:07:50 2009
!!Write out data in phases ( this might be quicker for configuration puposes of large files)

!!Updated by Roelof Rietbroek, Wed Jun 29 11:00:45 2011
!!fixed bug in Nunknowns nametag
!! Updated 13 sept 2017
!! also allow automatic replacement of G[CS]NV tags with G[CS]N+space

program GRACE2BIN
use binfiletools
use gpstlbx
use forttlbx
implicit none
integer::i,j,itharg,narg,stderr,lmax,lmin,l,m,ind
parameter(stderr=0)
character(120)::dum,basename
character(120),allocatable::hist(:)
double precision::ltpl,Stime,Etime,Ctime

integer::npara,nobs,npar_red,nhist
type(BINdat)::normout
logical::ignoapri,vnreplace
integer::iargc
!defaults initializations
itharg=0
lmax=0
lmin=9999
normout%file='stdout'
ignoapri=.false.
vnreplace=.false.
basename=''

!!!!!!!!!!!!!command line arguments processing!!!!!!!!!!!!!!!!!!!!!!!!!!1
!read in command line arguments
narg=iargc()

if(narg < 1)call help() !call help without any arguments


do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit !exit loop after last argument
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
         case('n')! ignore apriori values ( set  to zero)
            ignoapri=.true.
        case('r')!replace input which has VN in it with V+space
        vnreplace=.true.
        case default
         call help()
      end select

   else !argument is a geopotential header file
      ind=index(dum,'.h',.true.)
      if(ind .eq. 0)then
         write(stderr,*)"ERROR: input file does not appear to be a valid header file"
         stop
      end if
      basename=dum(1:ind-1)
   end if

   
end do

!!input checks
if(basename .eq. '')then
   write(0,*)"ERROR:no file specified"
   stop
end if


!!!!!!!end command line processing


!setup some meta data descriptions for the output system
normout%descr='GRACE normal system from K-band and GPSphase measurements'
normout%type='SYMVN___'
normout%mtyp='P'
normout%nvec=2

!integer meta data
normout%nint=7
allocate(normout%ints_d(normout%nint),normout%ints(normout%nint))
normout%ints_d(1)='Nobs'
normout%ints_d(2)='Nunknowns'
normout%ints_d(3)='Nreduced'
normout%ints_d(4)='Lmax'
normout%ints_d(5)='Lmin'
normout%ints_d(6)='Modnr'
normout%ints_d(7)='GPSweek'

!double meta data
normout%ndbls=5
allocate(normout%dbls_d(normout%ndbls),normout%dbls(normout%ndbls))
normout%dbls_d(1)='CTime'
normout%dbls_d(2)='STime'
normout%dbls_d(3)='ETime'
normout%dbls_d(4)='LtPL'
normout%dbls_d(5)='Sigma0'

!get info from header data

call read_GRACEnormhead(filename=trim(basename)//'.h',npar=npara,npar_red=npar_red,nobs=nobs,nhist=nhist)

allocate(hist(nhist))

!allocate matrix and vector data
normout%nval1=npara
normout%nval2=npara
normout%pval1=npara*(npara+1)/2
normout%pval2=1
allocate(normout%side1_d(npara),normout%vec(npara,normout%nvec))
allocate(normout%pack1(normout%pval1))

!readme part
normout%nread=nhist*2+1
allocate(normout%readme(normout%nread))
 

!get side description,apriori vector and time tags
if(ignoapri)then
   call read_GRACEnormhead(filename=trim(basename)//'.h',side_d=normout%side1_d,&
        Stime=Stime,Etime=Etime,hist=hist,ltpl=ltpl)
   normout%vec(:,2)=0.d0
else
   call read_GRACEnormhead(filename=trim(basename)//'.h',side_d=normout%side1_d,apriori=normout%vec(:,2)&
        ,Stime=Stime,Etime=Etime,hist=hist,ltpl=ltpl)
end if
Ctime=(Etime+Stime)/2.d0






!retrieve maximum and minimum degree
do,j=1,npara
   if(normout%side1_d(j)(1:3) .eq. 'GCN' .or. normout%side1_d(j)(1:3) .eq. 'GSN')then
      read(normout%side1_d(j),'(4x,i3,i3)')l,m
      lmax=max(lmax,l)
      lmin=min(lmin,l)
   end if
end do


if(vnreplace)then
    do,i=1,normout%nval1
        ind=index(normout%side1_d(i),'NV')
        if(ind .eq. 3)then
            normout%side1_d(i)(3:4)='N '
        end if
    end do    
    
end if

!put in the meta data
normout%ints(1)=nobs
normout%ints(2)=npara
normout%ints(3)=npar_red
normout%ints(4)=lmax
normout%ints(5)=lmin
normout%ints(6)=1
normout%ints(7)=GPS_week(sinex_date(Ctime))
 

normout%dbls(1)=Ctime
normout%dbls(2)=STime
normout%dbls(3)=ETime
normout%dbls(4)=ltpl
normout%dbls(5)=1.d0

!put in readme data
normout%readme(1)="Data converted from original file with history:"
j=1
do,i=1,nhist
   j=j+1
   normout%readme(j)=hist(i)(1:50)
   j=j+1
   normout%readme(j)=hist(i)(51:)
end do


!now write  first stages to file

call write_BINtype(normout,2)


!get binary part
call read_GRACEnormbin(filename=trim(basename)//'.b',npar=npara,pmat=normout%pack1,&
     Atb=normout%vec(:,1),ltpl=ltpl)

!write remaining (possibly big) part to file
call write_BINtype(normout) 


end program GRACE2BIN

!help routine prints usage information to screen
subroutine help()
character(4)::frmt
frmt='(A)'
write(*,frmt)'Program which converts GRACE GFZ normal systems in binary form to BINtype'
write(*,frmt)' Usage: GRACE2BIN [OPTIONS] GRACEFILE'
write(*,frmt)' Where GRACEFILE denotes the header file (ends with .h) of the normal system'
write(*,frmt)' OPTIONS may be:'
write(*,frmt)' -n: Ignore apriori values (set to zero)'
write(*,frmt)' -r: Replace the V in the side description with a space for G[CS]NV'
stop
end subroutine help

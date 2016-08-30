!!Coded by Roelof Rietbroek, Wed Sep  2 15:56:47 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

! Program to calculate a transformation matrix linking Spherical harmonic coefficients to pattern scales

program SH_transmat
use shtlbx
use binfiletools
implicit none
integer::i,j,itharg,narg,stderr,iargc
integer::maxf,lmax,lmin,ind,unit
parameter(maxf=1000) !maximum maount of input files
character(120)::dum,shfiles(maxf),tagfile,date,time
type(BINdat)::out
integer::nrow,ncol,l,m,q,nf,flmax,ftyp,last,pos
double precision,allocatable::clm(:)
double precision::W(maxf)
integer,allocatable::sort(:)
integer::nw,st,nd

!!defaults initializations
stderr=0
lmax=0 !maximum degree of output rows
lmin=0 ! minimimum degrees of output rows
nf=0 !amount of files
tagfile=''
unit=13
nw=0
W=1.d0

!!process command line options
narg=iargc()
if (narg <1)call help()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)

   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('n')!argument is the file name of column tags
         if(dum(4:4).eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,tagfile)
         else
            tagfile=dum(4:)
         end if
      case('l')!limit maximum and minimum degree
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
      case('W')!scale each input file with an individual scale
            st=4
            do,j=1,maxf
               nd=index(dum(st:),'/')+st-2
               if(nd .eq. st-2)nd=len(dum)!set to end of string if no '/' occurs (only one weight)
             
               read(dum(st:nd),*)W(j) !read in weight
             if(nd .eq. (len(dum)))exit !last weight has been read
             st=nd+2
            end do
            nw=j !set number of weights provided
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option ( but a file)
      nf=nf+1
      shfiles(nf)=trim(dum)
   end if
end do


!some input checks
if(nf==0)then
   write(stderr,*)"ERROR:SH_transmat, no input files supplied"
   stop
end if

if(nw .ne. nf .and. nw .ne. 0)then 
   write(stderr,*)"ERROR:SH_transmat, # of weights do not agree with # files"
   stop
end if

!try to determine the value of lmax if it is not provided on the command line

if(lmax==0)then
!    Retrieve lmax from the first file
   call SH_readmeta(filen=trim(shfiles(1)),type=ftyp,lmax=lmax)
end if

!Column description
ncol=nf
allocate(out%side2_d(ncol))
out%side2_d=''
if(tagfile .eq. '')then
   do,i=1,ncol

      write(out%side2_d(i),'(a3,i3)')'PAT',i
   end do

else
   open(unit=unit,file=tagfile)
   do,i=1,nf
      read(unit=unit,fmt='(A24)',iostat=last)out%side2_d(i)
      if(last <0)then
         write(stderr,*)"ERROR:SH_transmat, end of file reached for tagfile"
         stop
      end if
   end do
   close(unit)
end if


!Row description
nrow=SH_tpos(lmax,lmax,1,lmax,lmin)
allocate(out%side1_d(nrow))

out%side1_d=''
do,l=lmin,lmax
   !zero order part
   ind=SH_tpos(l,0,0,lmax,lmin)
   
   write(out%side1_d(ind),'(a4,i3,i3)')'GCN ',l,0
   do,m=1,l
      ind=SH_tpos(l,m,0,lmax,lmin)
      write(out%side1_d(ind),'(a4,i3,i3)')'GCN ',l,m
      ind=SH_tpos(l,m,1,lmax,lmin)
      write(out%side1_d(ind),'(a4,i3,i3)')'GSN ',l,m
   end do
end do

!setup output structure
out%nval1=nrow
out%nval2=ncol
out%file='stdout'
out%type='FULL2DVN'
out%mtyp='F'
out%descr='Fingerprint transformation matrix from SH_transmat'
out%pval1=out%nval1*out%nval2
out%pval2=1
out%nvec=0
! integers
out%nint=2
allocate(out%ints(out%nint),out%ints_d(out%nint))
out%ints_d(1)='Lmax'
out%ints(1)=lmax
out%ints_d(2)='Lmin'
out%ints(2)=lmin

!doubles
out%ndbls=0

!readme part
out%nread=4+nf
allocate(out%readme(out%nread))
out%readme=""
call date_and_time(date,time)
out%readme(1)='Created: '//date(7:8)//"-"//date(5:6)//"-"//date(1:4)//' '//time(1:2)//":"//time(3:4)
write(out%readme(2),*)'Scales and used input files:(',nf,')'
do,i=1,nf
   write(out%readme(2+i),'(E10.4,1x,A69)')W(i),shfiles(i)(1:69)
end do


allocate(out%mat1(out%nval1,out%nval2))


!read in SHfiles and put right column
pos=SH_pos(lmax,lmax)
allocate(clm(2*pos)) !clm contains clm (1:pos) AND slm (pos+1:2*pos)

!precompute index matrix
allocate(sort(nrow))
do,i=1,nrow
   call SH_tlm(i,l,m,q,lmax,lmin)
   sort(i)=(q*pos)+SH_pos(l,m)
end do

do,i=1,nf ! Loop over files
   call SH_readmeta(trim(shfiles(i)),ftyp,flmax)
   if(flmax<lmax)then
      write(stderr,*)"WARNING: file",trim(shfiles(i)),"only supports up to degree lmax:",flmax
      !stop
   end if
   call SH_readgrav(filen=trim(shfiles(i)),clm=clm(1:pos),slm=clm(pos+1:),type=ftyp)
   !put the values in the matrix
   out%mat1(:,i)=W(i)*clm(sort)
end do


!output matrix (at once)
call write_BINtype(out)


end program SH_transmat


subroutine help()
  implicit none
character(8)::frmt
integer::unit
unit=0
frmt='(A)'
write(unit,frmt)"Program SH_transmat calculates a transformation matrix in BINtype format"
write(unit,frmt)"Linking SH coefficients to predefined patterns"
write(unit,frmt)"Usage SH_transmat [OPTIONS] PATTERNSHFILES"
write(unit,frmt)"Where the PATTERNSHFILES contain the predefined patterns"
write(unit,frmt)"The output ( to standard output) is a matrix A"
write(unit,frmt)"| c_00 |     |  pattern_1_scale  |"
write(unit,frmt)"| c_10 |     |  pattern_2_scale  |"
write(unit,frmt)"|   .  | = A |        .          |"
write(unit,frmt)"|   .  |     |        .          |"
write(unit,frmt)"| c_nn |     |  pattern_k_scale  |"
write(unit,frmt)" Matrix A has size ~nmax^2 by k"
write(unit,frmt)" The names of the rows are based on potential Coefficients G[SC]N DEGORD"
write(unit,frmt)" while the column names default to PATT_1 PATT_2 etc."
write(unit,frmt)" "
write(unit,frmt)"The following options are allowed:"
write(unit,frmt)"-l=LMAX,LMIN: restrict the degree range of the rows of matrix A"
write(unit,frmt)"-n=NAMEFILE: use different names from the FILE NAMEFILE (1 name per row) "
write(unit,frmt)"-W=A1/A2/../An : apply scales to input files: number of weights must"
write(unit,frmt)"                 match the amopunt of input files"
write(unit,frmt)"   The file must contain at least as many rows as provided patternfiles"
write(unit,frmt)"   The order must correspond to the order of the input file"
stop
end subroutine help

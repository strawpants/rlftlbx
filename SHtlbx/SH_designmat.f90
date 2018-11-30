!! program which calculates a design matrix which propagates Spherical harmonics on a grid or point or on several regional base functions


!!TODO: a lot

!!Coded by Roelof Rietbroek, Tue Apr 28 12:42:43 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Wed Mar 10 09:07:58 2010
!! Additional multiplication with a diagonal spherical harmonic coefficient file is possible

!!Updated by Roelof Rietbroek, Fri Apr 30 16:51:56 2010
!! additional right hand multiplication with given matrix is possible (This allows for larger SH to point conversions when patterns are given)



program SH_designmat
use binfiletools
use forttlbx
use shtools
use shtlbx
implicit none
character(20)::version
character(200)::dum,pointfile
integer::itharg,narg,i,j,stderr,lmin,lmax,k
integer::l,m,qm,q,no,ind,indo,nsh,npara,npara2
integer::nlon,nlat,shft,gtyp
double precision::deg2rad
double precision::dres,lonmin,lonmax,latmin,latmax
logical::equigrid
double precision,allocatable::trivec(:,:,:)
double precision,allocatable::p(:)
double precision,pointer::tmp(:,:)=>null()
double precision,pointer::lon(:)=>null()
double precision,pointer::lat(:)=>null()
integer,allocatable::deg(:),ord(:),tri(:)
type(BINdat)::out,trans
integer::iargc,unit,chunk,mem_sz,last,rows,ipos
logical::stdin,rightmult
double precision::latold,lonold,tol

!!defaults
tol=0.d0
stderr=0
rightmult=.false.
stdin=.false.
deg2rad=pi/180.d0
version='ALPHA1.0'
lmax=30
lmin=0
unit=13
chunk=10000
!grid parameters
equigrid=.false. ! default assumes apoint file as input
dres=3
lonmin=0.d0
lonmax=360.d0
latmin=-90.d0
latmax=90.d0
pointfile=''
!!process command line options
narg=iargc()

itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('d')! equidistant grid sampling
         equigrid=.true.
         if(dum(3:3) .ne. '' )then! read in different grid interval
            read(unit=dum(4:),fmt=*,iostat=last)dres
            if(last .ne. 0 .or. dum(3:3) .ne. '=')then
               write(stderr,*)"Error processing option",trim(dum)
               stop
            end if
          end if
      case('l') ! adapt maximum and minimum degree
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(unit=dum(4:ind-1),fmt=*,iostat=last)lmax
            read(unit=dum(ind+1:),fmt=*,iostat=last)lmin
         else
            read(unit=dum(4:),fmt=*,iostat=last)lmax
         end if
         if(last .ne. 0 .or. dum(3:3) .ne. '=')then
            write(stderr,*)"Error processing option",trim(dum)
            stop
         end if
      case('r')
         equigrid=.true.
         
         ind=index(dum,'/')
         indo=3
         read(unit=dum(indo+1:ind-1),fmt=*,iostat=last)lonmin
         
         indo=ind
         ind=ind+index(dum(ind+1:),'/')
         read(unit=dum(indo+1:ind-1),fmt=*,iostat=last)lonmax
      
         indo=ind
         ind=ind+index(dum(ind+1:),'/')
         read(unit=dum(indo+1:ind-1),fmt=*,iostat=last)latmin

         read(unit=dum(ind+1:),fmt=*,iostat=last)latmax
         if(last .ne. 0 .or. dum(3:3) .ne. '=')then
            write(stderr,*)"Error processing option",trim(dum)
            stop
         end if

      !run some checks straight away
         if((latmax-latmin) <0.d0)then
            write(*,*)'minimum and maximum latitude swopped?: '
            call help(version)
         else if((lonmax-lonmin) <0.d0)then
            write(*,*)'minimum and maximum longitude swopped?: '
            call help(version)
         end if
      case('R') !additional right multiplication
         rightmult=.true.
         if(dum(4:4) .eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,trans%file)
         else
            trans%file=dum(4:)
         end if
      case('e') !use different tolerance for latitude comparison

         read(unit=dum(4:),fmt=*,iostat=last)tol
         if(last .ne. 0 .or. dum(3:3) .ne. '=')then
            write(stderr,*)"Error processing option",trim(dum)
            stop
         end if
      case('h')
         call help(version)
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
      pointfile=trim(dum)
      equigrid=.false.
   end if
end do

!input checks
if(pointfile .eq. '' .and. .not. equigrid)then ! read point file from standard input
   stdin=.true.
   unit=5
end if

if(rightmult)then ! load transformation matrix in memory

   call read_BINtype(trans,2)

   !check if supported and read remaining part
   select case(trans%type)
   case('DIAVN___') ! diagonal matrix ok, do nothing
   case('FULL2DVN') ! full matrix is ok read in 2d matrix
      trans%mtyp='F'
   case default
      write(stderr,*)"ERROR: type not supported"
      stop
   end select

   !possibly override lmax,lmin settings
   do, i=1,trans%nint
      select case(trans%ints_d(i))
      case('Lmax','lmax','nmax','Nmax')
         lmax=trans%ints(i)
      case('Lmin','lmin','nmin','Nmin')
         lmin=trans%ints(i)
      end select
   end do
   !read remainder
   call read_BINtype(trans)
end if

!SH parameter space
nsh=SH_tpos(lmax,lmax,1,lmax,lmin)

if(rightmult)then !take order of coefficients from matrix
   if(nsh .ne. trans%nval1)then
      write(stderr,*)"ERROR: side1 of transformation matrix is not consistent"
      write(stderr,*)"with amount of SH coefficients requested"
      stop
   end if

   !get degrees trigonometric sign and order from parameter
   allocate(deg(nsh),ord(nsh),tri(nsh))
   do,i=1,trans%nval1
      read(trans%side1_d(i),'(1x,A1,2x,i3,i3)')dum(1:1),deg(i),ord(i)

      select case(dum(1:1))
      case('C')
         tri(i)=0
      case('S')
         tri(i)=1
      case default
         write(stderr,*)"ERROR procesing SH parameter",trans%side1_d(i)
         stop
      end select
   end do
   !let the new side description point to the first side of trans
   out%side2_d => trans%side2_d
   npara2=trans%nval2
else !set up side desciption and index
   allocate(out%side2_d(nsh))
   !create index vectors and column side description
   !index vector
   allocate(deg(nsh),ord(nsh),tri(nsh))
   tri=0
   
   no=0
   do m=0,lmax
      qm=0
      if (m.gt.0) qm=1 
      do q=0,qm
         do l=max(lmin,m),lmax
            no=no+1
            deg(no)=l
            ord(no)=m
            tri(no)=q
            
            if(q .eq. 0)then    !cosine terms
               write(out%side2_d(no),'(a3,1x,i3,i3)')'GCN',l,m
            else                !sine term
               write(out%side2_d(no),'(a3,1x,i3,i3)')'GSN',l,m
            end if
         end do
      end do
   end do
   npara2=nsh

end if

if(equigrid)then ! do some precomputation and set up side 1
   !use pixel registration
   nlon=int((lonmax-lonmin))/dres
   nlat=int((latmax-latmin))/dres
   npara=nlon*nlat
   out%descr="created by SH_desigmat for mapping SH to an equidistant grid" 
   
   !create side description and fill up longitude and latitude vectors
   
   allocate(out%side1_d(npara))
   allocate(lon(nlon),lat(nlat))
   
   do,i=1,nlat
      lat(i)=latmin+(i-1)*dres+dres/2.d0         
   end do

   do,i=1,nlon
      lon(i)=lonmin+(i-1)*dres+dres/2.d0         
   end do

   !precompute longitude order vector
   allocate(trivec(1:nlon,0:1,0:lmax))
   trivec=0.d0
   !order zero vector
   trivec(:,0,0)=1.d0
   do,m=1,lmax
      do,i=1,nlon
         trivec(i,0,m)=cos(m*lon(i)*deg2rad)
      end do
      do,i=1,nlon
         trivec(i,1,m)=sin(m*lon(i)*deg2rad)
      end do
   end do
   
   !setup side description
   do,i=1,nlat
      ind=(i-1)*nlon
      do,j=1,nlon
         write(out%side1_d(ind+j),'(A3,1x,F8.3,1x,F7.3)')"GRD",lon(j),lat(i)
      end do
   end do
   
    allocate(out%mat1(npara,npara2),p(SH_pos(lmax,lmax)))
    if(rightmult)allocate(tmp(nlon,nsh)) ! allocate temporary matrix

    !create design matrix
    do,i=1,nlat
       !calculate Legendre function for given latitude
       call PlmBar(p=p, lmax=lmax, z=sin(lat(i)*deg2rad))
       shft=(i-1)*nlon
       !point straight away to right memory section in output array
       if(.not. rightmult)tmp=>out%mat1(shft+1:shft+nlon,1:nsh)
       !fill up temporary matrix
       do,k=1,nsh
          tmp(1:nlon,k)=p(SH_pos(deg(k),ord(k)))*trivec(1:nlon,tri(k),ord(k))
       end do
       if(rightmult)then !Additionally rightmultiply with (possibly diagonal) matrix
          select case(trans%type)
          case('FULL2DVN') !call blas matrix matrix multiplication
             call dgemm('N','N',nlon,npara2,nsh,1.d0,tmp(1,1),nlon,trans%mat1(1,1),nsh,0.d0,out%mat1(shft+1,1),npara)
          case('DIAVN__ ')!right multiply with diagonal matrix
             do,k=1,nsh
                out%mat1(shft+1:shft+nlon,k)=tmp(1:nlon,k)*trans%pack1(k)
             end do
          end select
       end if
    end do
else ! read lonlat values from pointfile
   out%descr="created by SH_desigmat for mapping SH to vectorized lon,lat " 
   allocate(out%side1_d(chunk))
   allocate(lon(chunk),lat(chunk))
   if(.not. stdin)open(unit=unit,file=pointfile,form='formatted')
   mem_sz=chunk
   npara=0
   last=0
   do while (last .eq. 0)! loop until end of file

      read(unit=unit,fmt='(A)',iostat=last)dum
      if(last .ne. 0)exit
      npara=npara+1

      if(npara > mem_sz)then ! reallocate memory if needed
         call realloc_ptr(out%side1_d,chunk)
         call realloc_ptr(lon,chunk)
         call realloc_ptr(lat,chunk)
         mem_sz=size(lon)
      end if
      out%side1_d(npara)=dum(1:24)
      read(dum(25:),*)lon(npara),lat(npara)
!      write(stderr,*)dum(1:24),lon(npara),lat(npara)
   end do

   if(.not. stdin)close(unit)


!create design matrix
   !define an appropriate amount of rows to perform dgemm routines

   
   allocate(out%mat1(npara,npara2),p(SH_pos(lmax,lmax)))

   
   if(rightmult)then    !allow a temporary matrix of 300 mb maximum
      chunk=(300*1024*1024/8)/nsh
      allocate(tmp(chunk,nsh))
   else
      chunk=npara
      tmp=>out%mat1 ! point to target memory straigth away (there is no use to do it in chunks)
   end if

   
   

   !allocate longitude order vector (hold only one longitude but all orders)
   allocate(trivec(1:1,0:1,0:lmax))
   trivec=0.d0
   trivec(1,0,0)=1.d0 ! order zero factor is independent from longitude
   latold=-9999.0
    !create design matrix
    do,i=1,(npara/chunk+1)
       shft=(i-1)*chunk
       rows=min(chunk,npara-shft)

       if(rows .eq. 0)exit !safety valve; everything is done already

       !loop over points within chunk
       do,j=1,rows
          !calculate index in absolute position vector
          ipos=shft+j

          !calculate Legendre function for given latitude (if recomputation is needed)
          if(abs(lat(ipos)-latold)>tol)then
             call PlmBar(p=p, lmax=lmax, z=sin(lat(ipos)*deg2rad))
             latold=lat(ipos) !update valid latitude used in Legendre function
          end if

          !precompute trigometric factors ( ir recomputation is needed)
          if(abs(lon(ipos)-lonold)>tol)then
             do,m=1,lmax
                trivec(1,0,m)=cos(m*lon(ipos)*deg2rad)
                trivec(1,1,m)=sin(m*lon(ipos)*deg2rad)
             end do
             lonold=lon(ipos)
          end if


       !fill up temporary matrix
       do,k=1,nsh !loop over columns of tmp
          ind=SH_pos(deg(k),ord(k))
          tmp(j,k)=p(ind)*trivec(1,tri(k),ord(k))
       end do
    end do !end j loop

    if(rightmult)then !Additionally rightmultiply with (possibly diagonal) matrix
       select case(trans%type)
       case('FULL2DVN') !call blas matrix matrix multiplication
          !DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
          call dgemm('N','N',rows,npara2,nsh,1.d0,tmp(1,1),chunk,trans%mat1(1,1),nsh,0.d0,out%mat1(shft+1,1),npara)
       case('DIAVN___')!right multiply with diagonal matrix
          do,k=1,nsh
             out%mat1(shft+1:shft+rows,k)=tmp(1:rows,k)*trans%pack1(k)
          end do
       end select
    end if
 end do !end i loop


end if








!setup remaining parameters in the derived type BINdat
out%file='stdout'
out%mtyp='F'
out%ldm=npara
out%nval2=npara2
out%nval1=npara
out%pval1=out%nval1*out%nval2
out%pval2=1
out%type='FULL2DVN'
out%nvec=0
out%ndbls=0
if(tol .ne. 0)then
   out%ndbls=1
   allocate(out%dbls_d(out%ndbls),out%dbls(out%ndbls))
   out%dbls_d(1)="Tolerance(degree)"
   out%dbls(1)=tol
end if

if(rightmult)then ! provide  some info in the file
   out%nint=0
   out%nread=4
   allocate(out%readme(out%nread))
   out%readme(1)="Right multiplied with matrix from file:"
   out%readme(2)=trans%file(1:80)
   out%readme(3)="in the  spectral domain with"
   write(out%readme(4),*)"Lmax=",lmax," and Lmin=",lmin
else
   out%nread=0
   out%nint=2
   allocate(out%ints_d(out%nint),out%ints(out%nint))
   out%ints_d(1)='Lmax'
   out%ints(1)=lmax
   out%ints_d(2)='Lmin'
   out%ints(2)=lmin
end if








!write to file/standard output
call write_BINtype(out)



end program SH_designmat

subroutine help(version)
implicit none
character(*)::version
character(8)::frmt
integer::unit
unit=6
frmt='(A)'
write(unit,frmt)"SH_designmat constructs a linear conversion matrix" 
write(unit,frmt)"Version: "//trim(version) 
write(unit,frmt)"Converting Spherical harmonic coefficients to (grid) points"
write(unit,frmt)""
write(unit,frmt)"Usage SH_designmat [OPTIONS] [POINTFILE]"
write(unit,frmt)"Where POINTFILE contains a 24 character parameter code and"
write(unit,frmt)" a longitude, latitude  (in degrees) per row. If not provided SH_desigmat reads"
write(unit,frmt)" from standard input"
write(unit,frmt)"and OPTIONS may be one or more of the following:"
write(unit,frmt)"-d[=RES]: conversion to equidistant grid with optional sampling interval."
write(unit,frmt)"         default takes RES=3 (degrees)" 
write(unit,frmt)"-r=lonmin/lonmax/latmin/latmax: Create a synthesis on a equidistant "
write(unit,frmt)"       grid bounded by minimum and maximum longitude and latitude"
write(unit,frmt)"       default: 0/360/-90/90 the values are centered on cell midpoints" 
write(unit,frmt)" The -d and -r option are ignored when a POINTFILE is provided"
write(unit,frmt)"-l=LMAX[,LMIN]: limit spherical harmonic maximum and minimum degree"
write(Unit,frmt)"               default: lmin=0, lmax=30" 
write(unit,frmt)"-R=MATFILE: Right multiply the matrix with a matrix (diagonal or full) from file"
write(unit,frmt)"   Matrix must be in binary format (readable) by BIN_swiss"
write(unit,frmt)"   matrix side1 must be purely SH coefficients matching the -l option"
write(unit,frmt)"-e=tol: Specify different tolerance for comparing next latitude (default=0)"
write(unit,frmt)" This can speed up computation when input is sorted according to latitude"
write(unit,frmt)"Program outputs a matrix in BINtype format to standard output"
write(unit,frmt)"" 
write(unit,frmt)"" 
write(unit,frmt)"" 
write(unit,frmt)"" 

stop


end subroutine help

!program which creates a Spherical harmonics product to sum conversion matrix for fully normalized coefficients
!matrix is written in unformatted form readable by the read_sym_mat routine from the FORTtlbx library

!!Coded by Roelof Rietbroek, Mon Sep 10 14:13:26 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!the matrix can be use to either mask(read convolve) spherical harmonic coefficients with another SH set
!Or the matrix can be used to regularize Spherical harmonic coefficients in an inversion

!required input parameters are 1.the maximum degree of the coefficient to solve for
!2.A spherical harmonic coefficient set (eg land or ocean mask) (fully normalized)
!note that the inputted spherical harmonic set must have at least a maximum solution degree of lmax*2
!else a warning message is given

!matrix is symmetric and stored in packed upper triangle format
! (same format as use by lapack routines) 


!Check the routines of SH_readgrav to see which SH files can be read in

!!Updated by Roelof Rietbroek, Mon Mar 10 09:56:37 2008
!!Output is now in new BIN format 
!!this program together with the routine SH_tripint can be made a lot quicker. Many calculations are redundant.
!!for example the longitude integral can be precomputed, there is a degree dependent

!!Updated by Roelof Rietbroek, Tue Jun 23 09:55:34 2009
!!added option to create a matrix of the complementary part of the SH set.


program SH_prod2sum2
use SHtlbx
use FORTtlbx
use binfiletools
implicit none
character(200)::dum,shfile
character(24),allocatable,dimension(:),target::list
integer::lmax,lmin,i,j,l1,l2,m1,m2,itharg,narg,nresv,pos
integer::indc1,inds1,indc2,inds2,ind,nsh,lmaxf,stderr,type
integer::deg1(3),done
double precision,allocatable,dimension(:)::triplm,clm,slm
double precision,target,allocatable::packmat(:)
double precision::fact
logical::load,compl
type(BINdat)::out ! output type
integer::iargc
!defaults/initializations
nresv=0
lmax=90
lmin=0
load=.false.
fact=rho_e/(sqrt(3.d0)*rho_w) !factor to descale degree 1 coefficients to geocenter motion
out%file='stdout' ! write to standard output
stderr=0
compl=.false.



!command line options

!get number of arguments

narg=iargc()

If(narg <1)call help()
!last argument is the SH file

call getarg(narg,shfile)
if(shfile(1:1) .eq. '-')call help()

!loop over remaining arguments
itharg=0
narg=narg-1 !don't count last option

do,i=1,narg
   if(itharg > narg)exit
   itharg=itharg+1
   call getarg(itharg,dum)
!   write(*,*)itharg,dum
   if(dum(1:1).ne. '-')then
      write(stderr,*)'argument ',itharg,' appears to be a file'
      write(stderr,*)'only one is allowed'
      call help()
   end if
   
   select case(dum(2:2))
!       case('f')!set output file name
         
!          if(dum(3:3) .eq. ' ')then!allow a space behind argument (easier command line typing)
!             itharg=itharg+1
!             call getarg(itharg,out%file)
!          else
!             out%file=dum(3:)
!          end if
      case('r')!make a regularization matrix for mass loading coefficients
         load=.true.
      case('l')!get maximum and minimum degree
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(3:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(3:),*)lmax
         end if
      case('c')!calculate complementary part
         compl=.true.
      case default
         write(stderr,*)'Unknown option ',dum(1:2),' quitting'
         call help()
   end select


end do


!get spherical harmonic coefficients from file
!first meta data

call SH_readmeta(filen=trim(shfile),type=type,lmax=lmaxf)
! write(*,*)type,lmaxf,trim(shfile)
!check if maximum degree supported by the file is sufficient
if(lmaxf < 2*lmax)then
   write(stderr,*)'The requested degree lmax ',lmax, ' is not supported by the file: ',trim(shfile)
   write(stderr,*)'lmax can be at most ',lmaxf/2,' for this file'
   call help()
end if

!if so allocate memory and read in data
pos=SH_pos(lmax*2,lmax*2)
allocate(clm(pos),slm(pos))
!allocate triple integral vector
allocate(triplm(pos))



clm=0.d0
slm=0.d0
triplm=0.d0

call SH_readgrav(filen=trim(shfile),clm=clm,slm=slm,type=type)

!do,i=1,pos
 !  write(*,*)i,clm(i),slm(i)
!end do
!now allocate memory for the packed array

nsh=SH_tpos(l=lmax,m=lmax,q=1,lmax=lmax,lmin=lmin)


allocate(packmat(((nsh+nresv)*(nsh+nresv+1))/2))



!initialize to zero
packmat=0.d0

!!! construct matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!loop over all degree and orders(riows then columns) 
done=0
do,l1=lmin,lmax
   write(stderr,*)'calculating degree:',l1
   do,l2=lmin,lmax
      do,m1=0,l1
      
     
      !get index of row 

      indc1=SH_tpos(l=l1,m=m1,q=0,lmax=lmax,lmin=lmin)+nresv      !cosine
      inds1=SH_tpos(l=l1,m=m1,q=1,lmax=lmax,lmin=lmin)+nresv      !sine

!loop over column degree and orders
         do,m2=0,l2
            !get index of column
            indc2=SH_tpos(l=l2,m=m2,q=0,lmax=lmax,lmin=lmin)+nresv      !cosine
            inds2=SH_tpos(l=l2,m=m2,q=1,lmax=lmax,lmin=lmin)+nresv      !sine
            
             !  write(*,*)allocated(triplm),allocated(packmat)
            
           !now four combinations are possible:
            !1: indc1, indc2 : row=cosine, column=cosine
            !2: indc1 ,inds2 : row=cosine, column=sine
            !3: inds1 ,indc2 : row=sine,   column=cosine
            !4: inds1 ,inds2 : row=sine,   column=sine

            !get triple integral vectors 
!           if(( indc1+(indc2*(indc2-1))/2) <1)write(*,*)l1,m1,l2,m2,indc1,indc2

            if(indc1 <= indc2)then !only calculate upper part of the matrix
               !only cosine coefficients have non-zero influence
               call SH_tripint_lm(l1=l1,m1=m1,q1=0,l2=l2,m2=m2,q2=0,qout=0,triplm=triplm)
!               write(*,*)'1st ',indc1+(indc2*(indc2-1))/2,(nsh*(nsh+1))/2,dot_product(triplm,clm)
               packmat(indc1+(indc2*(indc2-1))/2)=dot_productblas(triplm,clm)
            end if

            if((indc1 <= inds2) .and. (m2 .ne. 0))then !only calculate upper part of the matrix
               !only sine coefficients have non-zero influence
               call SH_tripint_lm(l1=l1,m1=m1,q1=0,l2=l2,m2=m2,q2=1,qout=1,triplm=triplm)
 !             write(*,*)'2nd ',indc1+(inds2*(inds2-1))/2,(nsh*(nsh+1))/2,dot_product(triplm,slm)
               packmat(indc1+(inds2*(inds2-1))/2)=dot_productblas(triplm,slm)
            end if
            
            if((inds1 <= indc2) .and. (m1 .ne. 0))then !only calculate upper part of the matrix
               !only sine coefficients have non-zero influence
               call SH_tripint_lm(l1=l1,m1=m1,q1=1,l2=l2,m2=m2,q2=0,qout=1,triplm=triplm)
  !             write(*,*)'3rd ',inds1+(indc2*(indc2-1))/2,(nsh*(nsh+1))/2,dot_product(triplm,slm)
               packmat(inds1+(indc2*(indc2-1))/2)=dot_productblas(triplm,slm)
            end if

            if((inds1 <= inds2) .and. (m1 .ne. 0) .and. (m2 .ne. 0))then !only calculate upper part of the matrix
               !only cosine coefficients have non-zero influence
               call SH_tripint_lm(l1=l1,m1=m1,q1=1,l2=l2,m2=m2,q2=1,qout=0,triplm=triplm)
   !            write(*,*)'4th ',inds1+(inds2*(inds2-1))/2,(nsh*(nsh+1))/2,dot_product(triplm,clm)
               packmat(inds1+(inds2*(inds2-1))/2)=dot_productblas(triplm,clm)
            end if

         end do !m2
      end do !l2
      
   end do !m1  
end do !l1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!matrix is now constructed!!!!!!!!!!!!!!

!write matrix to unformatted file
!write(*,*)lmax,lmin,trim(shfile),trim(outfile)
!write(*,*)packmat
!make list of parameter
allocate(list(nsh+nresv))
if(load)then
   
   if (lmin <=1 .and. lmax >= 1)then !convert degree 1 values to geocenter motion values
      deg1(1)=SH_tpos(l=1,m=0,q=0,lmax=lmax,lmin=lmin)+nresv
      deg1(2)=SH_tpos(l=1,m=1,q=0,lmax=lmax,lmin=lmin)+nresv
      deg1(3)=SH_tpos(l=1,m=1,q=1,lmax=lmax,lmin=lmin)+nresv

      write(list(deg1(1)),'(a3)')'SZB'
      write(list(deg1(2)),'(a3)')'SXB'
      write(list(deg1(3)),'(a3)')'SYB'

      !rescale matrix 
      
      do,i=1,nsh+nresv
         do,j=i,nsh+nresv
            select case((count(j .eq. deg1) +count(i .eq. deg1)))
            case(2)

                  packmat(i+j*(j-1)/2)=packmat(i+j*(j-1)/2)*fact*fact
               
            case(1)

               packmat(i+j*(j-1)/2)=packmat(i+j*(j-1)/2)*fact
            end select
         end do
         
      end do


   end if

   do,l1=lmin,lmax
      if(l1 .eq. 1)cycle !degree 1 is already accounted for
      write(list(SH_tpos(l=l1,m=0,q=0,lmax=lmax,lmin=lmin)+nresv),'(a3,1x,i3,i3)')'TCN',l1,0
      do,m1=1,l1 !order 1 and higher
         write(list(SH_tpos(l=l1,m=m1,q=0,lmax=lmax,lmin=lmin)+nresv),'(a3,1x,i3,i3)')'TCN',l1,m1
         write(list(SH_tpos(l=l1,m=m1,q=1,lmax=lmax,lmin=lmin)+nresv),'(a3,1x,i3,i3)')'TSN',l1,m1
      end do
   end do
else
do,l1=lmin,lmax
      write(list(SH_tpos(l=l1,m=0,q=0,lmax=lmax,lmin=lmin)+nresv),'(a3,1x,i3,i3)')'GCN',l1,0
      do,m1=1,l1 !order 1 and higher

         write(list(SH_tpos(l=l1,m=m1,q=0,lmax=lmax,lmin=lmin)+nresv),'(a3,1x,i3,i3)')'GCN',l1,m1
         write(list(SH_tpos(l=l1,m=m1,q=1,lmax=lmax,lmin=lmin)+nresv),'(a3,1x,i3,i3)')'GSN',l1,m1
      end do

   end do
end if


if(compl)then !take the complementary part of the matrix
   packmat=-packmat !first take the negative of the matrix
   
   !then add a 1 to the diagonal ( add unit matrix)
   
   do,i=1,nsh+nresv
      packmat(i*(i+1)/2)=packmat(i*(i+1)/2)+1.d0
   end do

end if



!make meta data for output file
ind=len_trim(shfile)
out%descr='Product 2 sum conversion matrix from:'//shfile(max(1,ind-40):ind)
out%nint=2
allocate(out%ints_d(out%nint),out%ints(out%nint))
out%ints_d(1)='Lmax'
out%ints_d(2)='Lmin'
out%ints(1)=lmax
out%ints(2)=lmin
out%mtyp='P'
out%type='SYMVN___'
out%nvec=0
out%ndbls=0
out%pack1=>packmat ! point to the right memory
out%side1_d=>list ! side description
out%nval1=nsh+nresv
out%nval2=out%nval1
out%pval1=out%nval1*(out%nval1+1)/2
out%pval2=1
out%nread=7
allocate(out%readme(out%nread))
out%readme(1)="Product 2 sum matrix, M, is obtained from wigner3j symbols"
out%readme(2)=" The matrix holds for fully normalized Spherical harmonics "
out%readme(3)=" applications: Regularization over specific areas"
out%readme(4)=" calculation of regionally specific harmonic base functions"
out%readme(5)=" This is not a projection matrix on a geographically limited area !!"
out%readme(6)=" which is impossible for a bandlimited spectrum. "
if(compl)then
   out%readme(7)="COMPLEMENTARY matrix, Mc, calculated: Mc=I-M"
else
   out%readme(7)=""
end if

!write to file/fifo
call write_BINtype(out)

!call write_BIN(file=trim(outfile),type='SYMV0___',descr=descr,ints_d=ints_d,ints=ints,side1_d=list,pmat1=packmat)

 
end program SH_prod2sum2

!!subroutine which prints a help message according to above program
subroutine help()
implicit none

write(*,'(A)')'Program SH_prod2sum_mat'
write(*,'(A)')'Constructs a spherical harmonic product-to-sum conversion matrix in binary format readable by read_BINtype'
write(*,'(A)')'Usage SH_prod2sum [OPTIONS] SHFILE'
write(*,'(A)')'Where the SHFILE is a fully normalized spherical harmonic coefficient file'
write(*,'(A)')'OPTIONS can be:'
write(*,'(A)')''
write(*,'(A)')'-r: Make a regularization matrix for mass loading coefficients'
!write(*,'(A)')'        (This can be handy for reserving space for helmert parameters later)'
write(*,'(A)')''
write(*,'(A)')'-lLMAX[,LMIN]: Output values for SH coefficients between LMIN and LMAX only (default:90,0)'
write(*,'(A)')'               note that the LMAX must be at most twice as large as the maximum degree supported in the file' 
write(*,'(A)')'-c: Take the complementary part of the matrix: Mc=I-M'
write(*,'(A)')' Program writes to standard output'
!write(*,'(A)')'-fOUTFILE: use the name OUTFILE for the output file (default: prod2sum.mat)'

stop
end subroutine help

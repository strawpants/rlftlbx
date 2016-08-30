!!Coded by Roelof Rietbroek, Thu Jun 17 11:57:22 2010
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!program which creates a block diagonal filter matrix
program SH_decormatSW
use forttlbx
use binfiletools
implicit none
character(200)::dum
integer::narg,iargc,i,j,itharg
integer::pdeg,lmax,lmin
integer::mstart,mstop
integer::stderr,chunk,mem
type(BINdat)::W !filter system
integer*8::pind,pindold
integer::sz,l,m,tri,lst,pfit,col,ind
integer::memsize,nsh,nblck,lda,info
double precision, dimension(:,:),allocatable::tmpA,AtA,ACAT,A
!!defaults
stderr=0
lmax=120
lmin=2
pdeg=-1
mstart=-1
mstop=-1
chunk=10000


!!process command line options
narg=iargc()
if (narg <2)call help() ! go to help straight away ( at least 2 options are obligatory)


itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:3))
      case('P=')!retrieve degree of the polynomial
         read(dum(4:),*)pdeg
      case('l=')!retrieve degree of the polynomial
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(4:),*)lmax
         end if
      case('M=')!retrieve degree of the polynomial
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(4:ind-1),*)mstart
            read(dum(ind+1:),*)mstop
         else
            read(dum(4:),*)mstart
         end if
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
   end if
end do



!input checks
if(mstart<0)then
   write(stderr,*)"ERROR:must supply -M option"
   stop
end if

if(pdeg<0)then
   write(stderr,*)"ERROR:must supply -P option"
   stop
end if



pfit=2*(pdeg+1) ! amounf of unkonw polynomials coefficients per order block ( twice polynomial coefficients for odd and even parity)

if(mstop<0)mstop=lmax-pfit-1 ! set default mstop so that unknowns can still be fitted in the largest block


if(mstart > min(mstop,lmax-pfit-1))then
   write(stderr,*)"ERROR:mstart > min(mstop,lmax-2*(pdeg+1)-1)"
   stop
end if

if(mstop > lmax-pfit-1)then
   write(stderr,*)"ERROR:mstop > lmax-2*(pdeg+1)-1 (polynomial not solvable)"
   stop
end if




!set up Filter system
!amount of blocks
if(mstart .eq. 0)then
   W%nblocks=(mstop-mstart+1)*2-1 ! don't calculate Sine order 0 block
else
   W%nblocks=(mstop-mstart+1)*2
end if


nsh=0
nblck=0
pind=0
allocate(W%side1_d(chunk))
allocate(W%blockind(W%nblocks))

memsize=chunk
do,m=mstart,mstop ! loop over order blocks
   do,tri=0,min(m,1) !Cosine/Sine loop
      do,l=max(lmin,m),lmax ! loop over degrees
         nsh=nsh+1
         if(nsh>memsize)then
            call realloc_ptr(W%side1_d,chunk)
            memsize=memsize+chunk
         end if

         if(tri .eq. 0)then !Cosine
            write(W%side1_d(nsh),'(A3,1x,2I3)')'GCN',l,m
         else !Sine
            write(W%side1_d(nsh),'(A3,1x,2I3)')'GSN',l,m
         end if
      end do
      nblck=nblck+1
      W%blockind(nblck)=nsh ! end of block in absolute side index
      sz=lmax-max(lmin,m)+1
      pind=pind+sz*(sz+1)/2 ! add packed block to matrix
   end do
end do

W%nval1=nsh
W%nval2=W%nval1
W%pval1=pind
W%pval2=1

! Allocate space for maximum design matrix
lda=lmax+1
allocate(A(0:lmax,pfit))
A=0.d0
!precompute design matrix (will be a chess board)
do,i=0,pdeg ! loop over polynomial expansion coefficients
   col=2*i+1
   do,l=0,lmax,2 ! loop over rows of A (even degrees)
      A(l,col)=(l)**i
   end do
   
   col=2*i+2
   do,l=1,lmax,2 ! loop over rows of A (odd degrees)
      A(l,col)=(l)**i
   end do
end do

!allocate Packed matrix
allocate(W%pack1(W%pval1))
allocate(AtA(pfit,pfit))
allocate(tmpA(0:lmax,pfit)) ! Same shape as A ( will hold a copy)

allocate(ACAT(lmax+1,lmax+1)) !allocate space for a dense block matrix ( only upper half will be used)

pind=0
do,m=mstart,mstop ! loop over order blocks
   lst=max(lmin,m) ! starting degree of block
   sz=lmax-lst+1 ! size of the block
   
   !copy data in tmpA
   tmpA(lst:lmax,:)=A(lst:lmax,:)
   !calculate AtA = tmpA **T * tmpA ( relevant degrees only)
   call dsyrk('U','T',pfit,sz,1.d0,tmpA(lst,1),lda,0.d0,AtA(1,1),pfit)
   
   !calculate cholesky decomposition U of AtA = (U**T U)
!   DPOTRF( UPLO, N, A, LDA, INFO )
   call dpotrf('U',pfit,AtA(1,1),pfit,info)
   if(info .ne. 0)then
      write(stderr,*)"ERROR:Cholesky decomposition failed for block with order:",m
   end if

!calculate tmpA = tmpA U**-1 ( tmpA is overwritten with the result)
!DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

   call dtrsm('R','U','N','N',sz,pfit,1.d0,AtA(1,1),pfit,tmpA(lst,1),lda)

!Initialize ACAT and set unit diagonal to 1
   ACAT=0.d0
   do,i=1,lmax+1
      ACAT(i,i)=1.d0
   end do

!Calculate ACAT = ACAT - Atmp * Atmp**T
   
   call dsyrk('U','N',sz,pfit,-1.d0,tmpA(lst,1),lda,1.d0,ACAT(1,1),lda)



! now copy values from the dense matrix in the packed matrix
   do,i=1,sz
      do,j=1,i
         W%pack1(pind+j+i*(i-1)/2)=ACAT(j,i) ! only upper triangle
      end do
   end do

!   write(0,*)'COS order',m,lst,sz,pind,ACAT(1,1),ACAT(2,1),ACAT(2,2),W%pack1(pind+1:pind+3)
   pindold=pind   

   pind=pind+sz*(sz+1)/2
   !copy complete Cosine block in the sine block
   if(m .ne. 0)then
      W%pack1(pind+1:pind+sz*(sz+1)/2)=W%pack1(pindold+1:pindold+sz*(sz+1)/2)
 !     write(0,*)'SIN order',m,lst,sz,pind,ACAT(1,1),ACAT(2,1),ACAT(2,2),W%pack1(pind+1:pind+3)
      pind=pind+sz*(sz+1)/2
   end if
end do


!set up remaining filter system parameters
W%type='BDSYMVN_' ! symmetric blockdiagonal matrix
W%mtyp='P'
W%descr="Swenson and Wahr 2006 decorrelation matrix"
W%file='stdout'
W%nint=5
allocate(W%ints_d(W%nint),W%ints(W%nint))
W%ints_d(1)='Lmax'
W%ints(1)=lmax
W%ints_d(2)='Lmin'
W%ints(2)=lmin
W%ints_d(3)='max Poly degree'
W%ints(3)=pdeg
W%ints_d(4)='Mstart'
W%ints(4)=mstart
W%ints_d(5)='Mstop'
W%ints(5)=mstop
W%nread=0


call write_BINtype(W)


end program SH_decormatSW




subroutine help()
integer::unit
character(8)::frmt
frmt='(A)'
unit=0 ! print to standard error
write(unit,frmt)"Program SH_decormatSW creates a block diagonal decorrelation matrix" 
write(unit,frmt)" according to Swenson and Wahr 2006 (but with a generalized polynomial)" 
write(unit,frmt)"Usage SH_decormatSW [OPTIONS]" 
write(unit,frmt)"Where options are:" 
write(unit,frmt)"-P=N: set the maximum degree of the polynomial fit" 
write(unit,frmt)"-M=START[,STOP]: start (and optionally stop) the filtering at "
write(unit,frmt)" order START (STOP: default LMAX-(N+2))" 
write(unit,frmt)"-l=LMAX[,LMIN]: set the maximum (default 120) (and optionally) minimum degree (default is 2) of the filter matrix" 
write(unit,frmt)"A block diagonal matrix in binary form (readable by BIN_swiss) will be printed"
write(unit,frmt)" to standard output" 
write(unit,frmt)" Currently no defaults are set for N and START but a reeasonable value is N=4 and START=6"
write(unit,frmt)" used in (Chen et al 2008) "
stop

end subroutine help

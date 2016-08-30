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

!!Updated by Roelof Rietbroek, Fri Sep  4 14:36:06 2009
!!created openmp version with as much precomputation as possible (at the cost of some memory)

!!Updated by Roelof Rietbroek, Wed Oct 10 13:15:24 2012
!! increased maximum supported degree to 180


program SH_prod2sum2
use sh_trpnt
use SHtlbx
use FORTtlbx
use binfiletools
implicit none
character(200)::dum,shfile
character(24),allocatable,dimension(:),target::list
integer::lmax,lmin,i,j,l1,l2,m1,m2,itharg,narg,pos
integer::ind,lmaxf,stderr,type
integer*8::indp,indccol,indscol,indcrow,indsrow,nsh ! large integer to refer to index in packed array ( neccessary for big arrays > 2 Gb)
integer::deg1(3)
integer,allocatable::sort(:)
! double precision,allocatable,dimension(:)::clm,slm
double precision,target,allocatable::packmat(:)
double precision::fact
logical::load,compl
type(BINdat)::out ! output type
integer::iargc
integer::lmaxsup
parameter(lmaxsup=180) !unfortunately we need to allocate some extra memory in order to avoid allocatable private arrays in openmp (doesn't work in buggy ifort)
double precision,dimension(2*lmaxsup+1)::veca,vecb
double precision,dimension((lmaxsup+1)*(lmaxsup+2))::cworkcol,sworkcol
double precision,dimension((2*lmaxsup+2)*(2*lmaxsup+1)/2)::clm,slm
double precision::cddota,cddotb,sddota,sddotb,ddot
integer::starta,startb,m3a,m3b,q1,q2,q3
integer::sza,szb,shifta,shiftb
integer::lmin3a,lmin3b,lmax3a,lmax3b
integer::memout
character(120)::message,date,time
logical::verbose
character(40)::timetag(10)
double precision::t0,t3,timing(10),t1,t2
integer::nt
integer::meth! timing function method
logical::precomp
type(trpnt_memstruct)::tripint
!$ integer::OMP_get_thread_num


!defaults/initializations
meth=0
!$ meth=1 !reset to omp wall time method

lmax=90
lmin=0
load=.false.
fact=rho_e/(sqrt(3.d0)*rho_w) !factor to descale degree 1 coefficients to geocenter motion
out%file='stdout' ! write to standard output
stderr=0
compl=.false.
verbose=.false.
message="finished degree"
precomp=.true.
t0=time_func(meth)

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
      case('v')
         verbose=.true.
      case('n')! don't precompute common values
         precomp=.false.
      case default
         write(stderr,*)'Unknown option ',dum(1:2),' quitting'
         call help()
   end select


end do


if(lmax > lmaxsup)then
   write(stderr,*)"The requested maximum degree",lmax," corresponds to"
   write(stderr,*)" building a monster and is therefore larger then supported",lmaxsup
   stop
end if

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
!allocate(clm(pos),slm(pos))



clm=0.d0
slm=0.d0
call SH_readgrav(filen=trim(shfile),clm=clm(1:pos),slm=slm(1:pos),type=type)

!sort clm and slm according to order instead of degree
allocate(sort(pos))
ind=0
do,m1=0,lmax*2
   do,l1=m1,lmax*2
      ind=ind+1

      sort(ind)=SH_pos(l1,m1)
!      write(0,*)l1,m1,clm(SH_pos(l1,m1)),slm(SH_pos(l1,m1))
   end do
end do


!resort according to order
call permute(clm(1:pos),sort,.false.)
call permute(slm(1:pos),sort,.false.)


!Initialize/precompute varaibles for tripleint routines
if(precomp)then
   call trpnt_init(tripint,lmax,memout)
   
   if(verbose)write(stderr,*)"Precomputation used",memout,"Kb of memory"
end if

!now allocate memory for the packed array
nsh=SH_tpos(l=lmax,m=lmax,q=1,lmax=lmax,lmin=lmin)
allocate(packmat(((nsh)*(nsh+1))/2))
!initialize to zero
packmat=0.d0

!!! construct matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nt=1

timetag(1)='degree loop'
timetag(2)='trpnt_vec_lm'
timetag(3)='dot product'
timetag(4)='filling workcols'
timetag(5)='updating shared memory'

!OMP directive to parallize degree loop
!Note since the workload is not the same per loop it is adviced to use a dynamic thread asssignment (use the schedule derivative)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(packmat,clm,slm,verbose,stderr,meth,timetag,tripint) &
!$OMP FIRSTPRIVATE(nsh,lmax,lmin)
!$OMP DO SCHEDULE(GUIDED,1) 
!do,l1=lmin,lmax
do,m1=0,lmax
   !$ if(verbose)write(stderr,*)"Assigning order",m1,"To Thread",OMP_get_thread_num()
   if(verbose)then
      timing=0.d0
      t1=time_func(meth)
   end if
!   do,m1=0,l1
   do,l1=max(m1,lmin),lmax 

            
      !get index of col 
      indccol=SH_tpos(l=l1,m=m1,q=0,lmax=lmax,lmin=lmin)      !cosine
      if(m1>0)then
         indscol=SH_tpos(l=l1,m=m1,q=1,lmax=lmax,lmin=lmin)      !sine
         sworkcol(1:indscol)=0.d0
      else !set to zero ( always in lower triangle)
         indscol=0 
      end if
      !Initialize working column to zero ( only necessary part)

      cworkcol(1:indccol)=0.d0

      !loop over row degree and orders
      do,l2=lmax,lmin,-1
         do,m2=0,l2
            !get index of column
            indcrow=SH_tpos(l=l2,m=m2,q=0,lmax=lmax,lmin=lmin)      !cosine
            if(m2>0)then
               indsrow=SH_tpos(l=l2,m=m2,q=1,lmax=lmax,lmin=lmin)      !sine
            else
               indsrow=nsh+1 !always in lower triangle
            end if

            !quick cycle if all valid values are in the lower part of the triangle
            if((indccol<indcrow) .and. (indccol<indsrow) .and. (indscol< indcrow) .and. (indscol< indsrow))then
               cycle
            end if

!             calculate l1,m1,l2,m2 dependent part of the triple integral (only changes veca and vecb)
!            if(verbose)t2=time_func(meth)
            if(precomp)then
               call trpnt_vec_lm(l1,l2,m1,m2,m3a,m3b,lmin3a,lmin3b,&
                    lmax3a,lmax3b,veca,vecb,tripint)
            else
               call trpnt_vec_lm(l1,l2,m1,m2,m3a,m3b,lmin3a,lmin3b,lmax3a,&
                 lmax3b,veca,vecb)
            end if
!            if(verbose)timing(2)=timing(2)+time_func(meth)-t2 !add time spend in trpnt_vec_lm
            !get shift parameters (0 or 1) to align with valid degree parity
            shifta=mod(l1+l2+lmin3a,2) !sum of the valid degrees must be even
            shiftb=mod(l1+l2+lmin3b,2) !sum of the valid degrees must be even

            !start indices in clm and slm

            starta=SH_opos(lmin3a,m3a,lmax*2)
            startb=SH_opos(lmin3b,m3b,lmax*2)
            !lenght of the to be computed dot products
            sza=lmax3a-lmin3a+1-shifta
            szb=lmax3b-lmin3b+1-shiftb

!DOT PRODUCT USING BLAS ( NOT WORKING ON SUN PLATFORM :(  )
            !calculate dotproduct (with stride 2 to ignore zero entries) of clm times valid veca/vecb
!             cddota=ddot(sza,veca(lmin3a+shifta+1),2,clm(starta+shifta),2)
!             if(m3a .ne. m3b)cddotb=ddot(szb,vecb(lmin3b+shiftb+1),2,clm(startb+shiftb),2)
!             !calculate dotproduct (with stride 2 to ignore zero entries) of slm times valid veca/vecb
!             if(m3a .ne. 0)sddota=ddot(sza,veca(lmin3a+shifta+1),2,slm(starta+shifta),2)
!             if(m3b .ne. 0 .and. (m3a .ne. m3b))sddotb=ddot(szb,vecb(lmin3b+shiftb+1),2,slm(startb+shiftb),2)

!Straightforward dot product ( using stride 2)

            cddota=0.d0
            do,i=1,sza,2
               cddota=cddota+veca(lmin3a+shifta+i)*clm(starta+shifta+i-1)
            end do

            cddotb=0.d0
            if(m3a .ne. m3b)then
               do,i=1,szb,2
                  cddotb=cddotb+vecb(lmin3b+shiftb+i)*clm(startb+shiftb+i-1)
               end do
            end if
               
            
            sddota=0.d0
            if(m3a .ne. 0)then
               do,i=1,sza,2
                  sddota=sddota+veca(lmin3a+shifta+i)*slm(starta+shifta+i-1)
               end do
            end if
            sddotb=0.d0
            if(m3b .ne. 0 .and. (m3a .ne. m3b))then
               do,i=1,szb,2
                  sddotb=sddotb+vecb(lmin3b+shiftb+i)*slm(startb+shiftb+i-1)
               end do
            end if
            
           !now four non zero combinations are possible:
            !1: indccol, indcrow :column=cosine, row=cosine
            !2: indccol ,indsrow :column=cosine, row=sine
            !3: indscol ,indcrow :column=sine,   row=cosine
            !4: indscol ,indsrow :column=sine,   row=sine
   !         if(verbose)t2=time_func(meth)
            if(indccol >= indcrow)then !only calculate upper part of the matrix
               !only cosine coefficients have non-zero influence
               !ind=indccol+(indcrow*(indcrow-1))/2
               !packmat(ind)=trp_qfaca(0,0,0,m1,m2)*cddota
               cworkcol(indcrow)=trp_qfaca(0,0,0,m1,m2)*cddota
               if(m3a .ne. m3b)then
!                  packmat(ind)=packmat(ind)+trp_qfacb(0,0,0,m1,m2)*cddotb
                  cworkcol(indcrow)=cworkcol(indcrow)+trp_qfacb(0,0,0,m1,m2)*cddotb
               end if
            end if

            
            if((indccol >= indsrow) .and. (m2 .ne. 0).and. (m3a .ne. 0))then !only calculate upper part of the matrix
               !only sine coefficients have non-zero influence
               !ind=indccol+(indsrow*(indsrow-1))/2
               !packmat(ind)=trp_qfaca(0,1,1,m1,m2)*sddota
               cworkcol(indsrow)=trp_qfaca(0,1,1,m1,m2)*sddota
               if(m3a .ne. m3b)then
                  !packmat(ind)=packmat(ind)+trp_qfacb(0,1,1,m1,m2)*sddotb
                  cworkcol(indsrow)=cworkcol(indsrow)+trp_qfacb(0,1,1,m1,m2)*sddotb
               end if
            end if
            
            if((indscol >= indcrow) .and. (m1 .ne. 0) .and. (m3a .ne. 0))then !only calculate upper part of the matrix
               !only sine coefficients have non-zero influence
               !ind=indscol+(indcrow*(indcrow-1))/2
               !packmat(ind)=trp_qfaca(1,0,1,m1,m2)*sddota
               sworkcol(indcrow)=trp_qfaca(1,0,1,m1,m2)*sddota
               if(m3a .ne. m3b)then
                  !packmat(ind)=packmat(ind)+trp_qfacb(1,0,1,m1,m2)*sddotb
                  sworkcol(indcrow)=sworkcol(indcrow)+trp_qfacb(1,0,1,m1,m2)*sddotb
               end if
            end if

            if((indscol >= indsrow) .and. (m1 .ne. 0) .and. (m2 .ne. 0))then !only calculate upper part of the matrix
               !only cosine coefficients have non-zero influence
               !ind=indscol+(indsrow*(indsrow-1))/2
               !packmat(ind)=trp_qfaca(1,1,0,m1,m2)*cddota
               sworkcol(indsrow)=trp_qfaca(1,1,0,m1,m2)*cddota
               if(m3a .ne. m3b)then
                  !packmat(ind)=packmat(ind)+trp_qfacb(1,1,0,m1,m2)*cddotb
                  sworkcol(indsrow)=sworkcol(indsrow)+trp_qfacb(1,1,0,m1,m2)*cddotb
               end if
            end if
     !       if(verbose)timing(4)=timing(4)+time_func(meth)-t2 !add time spend on workcolfill
         end do !m2
         
      end do !l2
      ! copy data in packed (shared) array in bulk ( getting it out of l2,m2 loop hopefully avoids some cache misses)
    !  if(verbose)t2=time_func(meth)
      !sine column
      if(m1 > 0)then
         indp=(indscol)*(indscol-1)/2 !start of column in the packed matrix
      
         packmat(indp+1:indp+indscol)=sworkcol(1:indscol)
!          write(0,*)indscol,"Thread ",OMP_get_thread_num()&
!               ,m1,maxval(sworkcol(1:indscol)),minval(sworkcol(1:indscol))
      end if
      
      !cosine column
      indp=(indccol)*(indccol-1)/2
      packmat(indp+1:indp+indccol)=cworkcol(1:indccol)
!       write(0,*)indccol,"Thread ",OMP_get_thread_num(),&
!            m1,maxval(sworkcol(1:indccol)),minval(cworkcol(1:indccol))
!      if(verbose)timing(5)=timing(5)+time_func(meth)-t2 !add time spend on updating shared memory
!       if(indccol .eq. 81)then
!          do,i=1,indccol
!             indp=i +(indccol)*(indccol-1)/2
!             write(0,*)i,cworkcol(i),packmat(indp)
!          end do
!          stop
!       end if
   end do !m1  
  if(verbose)then
     timing(1)=timing(1)+time_func(meth)-t1
     !$ write(message,*)"Thread ",OMP_get_thread_num()," finished order"
     write(stderr,*)trim(message),m1,"in",timing(1),"sec"
    
  end if
end do !l1
!$OMP END DO
!$OMP END PARALLEL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!matrix is now constructed!!!!!!!!!!!!!!

!write matrix to unformatted file
!write(*,*)packmat
!make list of parameter
allocate(list(nsh))
if(load)then
   if(verbose)write(stderr,*)"Adapting matrix to loading regularization"
   if (lmin <=1 .and. lmax >= 1)then !convert degree 1 values to geocenter motion values
      deg1(1)=SH_tpos(l=1,m=0,q=0,lmax=lmax,lmin=lmin)
      deg1(2)=SH_tpos(l=1,m=1,q=0,lmax=lmax,lmin=lmin)
      deg1(3)=SH_tpos(l=1,m=1,q=1,lmax=lmax,lmin=lmin)

      write(list(deg1(1)),'(a3)')'SZB'
      write(list(deg1(2)),'(a3)')'SXB'
      write(list(deg1(3)),'(a3)')'SYB'

      !rescale matrix 
      
      do,i=1,nsh
         do,j=i,nsh
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
      write(list(SH_tpos(l=l1,m=0,q=0,lmax=lmax,lmin=lmin)),'(a3,1x,i3,i3)')'TCN',l1,0
      do,m1=1,l1 !order 1 and higher
         write(list(SH_tpos(l=l1,m=m1,q=0,lmax=lmax,lmin=lmin)),'(a3,1x,i3,i3)')'TCN',l1,m1
         write(list(SH_tpos(l=l1,m=m1,q=1,lmax=lmax,lmin=lmin)),'(a3,1x,i3,i3)')'TSN',l1,m1
      end do
   end do
else
do,l1=lmin,lmax
      write(list(SH_tpos(l=l1,m=0,q=0,lmax=lmax,lmin=lmin)),'(a3,1x,i3,i3)')'GCN',l1,0
      do,m1=1,l1 !order 1 and higher

         write(list(SH_tpos(l=l1,m=m1,q=0,lmax=lmax,lmin=lmin)),'(a3,1x,i3,i3)')'GCN',l1,m1
         write(list(SH_tpos(l=l1,m=m1,q=1,lmax=lmax,lmin=lmin)),'(a3,1x,i3,i3)')'GSN',l1,m1
      end do

   end do
end if


if(compl)then !take the complementary part of the matrix
   if(verbose)write(stderr,*)"Taking complementary part of matrix"
   packmat=-packmat !first take the negative of the matrix
   
   !then add a 1 to the diagonal ( add unit matrix)
   
   do,i=1,nsh
      packmat(i*(i+1)/2)=packmat(i*(i+1)/2)+1.d0
   end do

end if



!make meta data for output file
ind=len_trim(shfile)
out%descr='Produc_2_sum conversion matrix from:'//shfile(max(1,ind-40):ind)
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
out%pack1=>packmat ! point to the right memory section
out%side1_d=>list ! "           " side description
out%nval1=nsh
out%nval2=out%nval1
out%pval1=out%nval1*(out%nval1+1)/2
out%pval2=1
out%nread=9
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
call date_and_time(date,time)
t3=time_func(meth)
write(out%readme(8),*)"Created: "//date(7:8)//"-"//date(5:6)//"-"//date(1:4),' ',time(1:2)//":"//time(3:4)&
     //", computation took:",t3-t0,"sec"


!$ out%readme(9)="Constructed with openMP enabled"

!write to file/fifo
if(verbose)write(stderr,*)"Outputting matrix"
t3=time_func(meth)
call write_BINtype(out)

if(verbose)then
   write(stderr,*)"Total excution time",t3-t0,"seconds"
   t0=time_func(meth)
   write(stderr,*)"Matrix output in",t0-t3,"seconds"

end if

contains
double precision function time_func(meth)
  implicit none
  integer,intent(in)::meth
  integer::count,rate
!$  double precision::omp_get_wtime

  select case(meth)
  case(0)
     call system_clock(count=count,count_rate=rate)
     time_func=dble(count)/rate
  case(1)
     !$    time_func=omp_get_wtime()
  case default
     time_func=0.0
  end select

end function time_func

end program SH_prod2sum2

!!subroutine which prints a help message according to above program
subroutine help()
implicit none
character(120)::parames1,parames2
parames1="This program is compiled with openMP disabled"
parames2="For large degrees you are adviced to enable openMP compiler options"
!possibly redefine parallelization messages
!$ parames1="This program is compiled with openMP enabled."
!$ parames2="You may set the OMP_NUM_THREADS environment variable for parallelization"


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
write(*,'(A)')"-n: Don't precompute common elements"
write(*,'(A)')'-v: Verbose: display additional messages'
write(*,'(A)')' Program writes to standard output'
write(*,'(A)')parames1
write(*,'(A)')parames2



stop
end subroutine help




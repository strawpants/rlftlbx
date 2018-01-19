!!Coded by Roelof Rietbroek, Fri Feb  8 16:10:31 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!adapted from SH_filter

!!program SH_filter2 applies anisotropic filter(kusche 2007) to Spherical harmonic coefficinets sets
!!both full filter matrcies and blockdiagonal matrices are supported
!!matrix format is that readable by read_BIN routines

!!Updated by Roelof Rietbroek, Fri Oct 10 16:16:22 2008
!!copies entries with degrees lower than the filter lmin from the original file to the output file 
!!Updated by Roelof Rietbroek, Tue Jan 13 12:23:29 2009
!! new matrix type may be read in

!!Updated by Roelof Rietbroek, Thu Jan 15 15:29:51 2009
!!also symmetric matrix may be read in (for eg ocean masking)

!!Updated by Roelof Rietbroek, Fri Feb 13 10:18:10 2009
!! fixed bug in error propagation ( forgot to square results to get standard deviations)
!! permuted columns at once 

!!Updated by Roelof Rietbroek, Fri Jun 18 12:06:07 2010
!! added multiplication for block diagonal symmetric matrix

!!Updated by Roelof Rietbroek, Tue Jun 22 16:22:25 2010
!! added transposing of filter matrix

!!Updated by Roelof Rietbroek, Wed Jul  4 14:56:26 2012
!!optionally apply a remove-restore field

!!Updated by Roelof Rietbroek, Mon Mar 18 13:29:58 2013
!!The complement of the filter matrx can now also be given (I-W)
!! this is convenient for the application of masks and their complement
!! also removed restrictions on the maount of input files

!!Updated by Roelof Rietbroek, Sat Nov 21 11:09:12 2015
!! default now strips the leading directories from the file name


!!Updated by Roelof Rietbroek, Thu Dec  3 11:49:52 2015
!! better treatment of remove-restore operation when files have different maximum degrees

!!Updated by Roelof Rietbroek, Tue Mar 29 12:40:43 2016
!! output directory can now be set


program SH_filter2
use binfiletools
use forttlbx
use shtlbx
implicit none
integer::narg,i,nfiles,filt,lmax,lmin,lmaxW,lminW,nd,p1,typ,lmaxdum,itharg
integer::l,m,stderr,side,full,sz,j,ind,k
parameter(stderr=0)
integer,allocatable,dimension(:)::pos,ints
character(400)::dum,ref_sh
character(20)::app
character(1)::trig
Character(400),pointer,dimension(:)::SHfiles
logical::stdin,stdout,limdeg,error,blck,sym,replsuf
double precision,allocatable,dimension(:,:)::tmp,tmp_sig,SH,time,blocktmp
double precision,allocatable,dimension(:)::ref_field
type(BINdat)::Wsys
integer::iargc,ldsh
character(1)::trans
logical::rem_rest,compl
character(400)::outdir

stdin=.false.
stdout=.false.
limdeg=.false.
error=.false.
compl=.false.
nfiles=0
lmax=0
lmin=0
app='.filt'
trans='N' 
ref_sh=''
itharg=0
blck=.false.
lmaxW=0
lminW=1d6
replsuf=.false.
rem_rest=.false.
outdir=""
allocate(SHfiles(200))

!process command line options
narg=iargc()
if(narg < 1) call help() !call help when too litle arguments are provided
do,i=1,narg
   itharg=itharg+1
   if(itharg> narg) exit

   call getarg(itharg,dum)
   if(dum(1:1) .ne. '-')then
      nfiles=nfiles+1
      if(nfiles > size(SHfiles,1))then !realllocate more space
         call realloc_ptr(SHfiles,200)
      end if
      SHfiles(nfiles)=trim(dum)
      cycle
   end if

!command line options
   select case(dum(2:2))
   case('W')
      if(dum(3:3) .eq. ' ')then
         itharg=itharg+1
         call getarg(itharg,Wsys%file)
      else
         Wsys%file=trim(dum(3:))
         
      end if
   case('T') ! transpose matrix
      trans='T'
   case('l') !read lmax and lmin
      p1=index(dum,',')
      read(dum(3:p1-1),*)lmax
      read(dum(p1+1:),*)lmin   
      limdeg=.true.
   case('e')!propagate diagonal errors
      error=.true.
      
   case('s')
      stdout=.true.
   case('a')
      app='.'//trim(dum(3:))
   case('r')!replace suffix
      replsuf=.true.
      app='.'//trim(dum(3:))
   case ('R')! apply a remove and restore action
      if(dum(3:3) .eq. ' ')then
         itharg=itharg+1
         call getarg(itharg,ref_sh)
      else
         ref_sh=trim(dum(3:))
      end if
      rem_rest=.true.
   case('c') !take the complement
      compl=.true.
   case('o') !Output directory
      outdir=trim(dum(3:))//"/"
      
      if(dum(3:3) .eq. ' ')then
         itharg=itharg+1
         call getarg(itharg,outdir)
      else
         outdir=trim(dum(3:))
      end if
      !possibly add a slash add the end
      if ( outdir(len_trim(outdir):len_trim(outdir)) .ne. '/')then
        outdir=trim(outdir)//"/"
      end if

   case default
      write(stderr,*)'unkown option selected' 
      call help()
   end select

end do

!check if output must be read from standard input and written to stdout
if(nfiles .eq. 0)then
   stdin=.true.
   stdout=.true.
   nfiles=1 !set to one to go through the loop anyway
end if   




!!!anisotropic kusche filter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!now if kusche type filter load unformatted weightfile

!!read index data
if(.not. stdout)write(stderr,*)'reading Weightfile ',trim(Wsys%file)
call read_BINtype(Wsys,2)

select case(Wsys%type)
case('BDFULLV0','BDFULLVN')
   Wsys%mtyp='P' ! keep packed form
   blck=.true.
   if(.not. stdout)write(stderr,*)'Block diagonal filter  matrix found'
case('BDSYMV0_','BDSYMVN_')
   Wsys%mtyp='P' ! keep packed form
   blck=.true.
   if(.not. stdout)write(stderr,*)'Block diagonal symmetric filter  matrix found'
case('FULLSQV0','FULLSQVN')
   Wsys%mtyp='F' ! unpack matrix
   allocate(Wsys%blockind(1))
   Wsys%blockind(1)=Wsys%nval1
   Wsys%nblocks=1
   if(.not. stdout)write(stderr,*)'full matrix filter found'
case('SYMVN___','SYMV1___','SYMV0___')
   Wsys%mtyp='U' ! fill upper triangle
   allocate(Wsys%blockind(1))
   Wsys%blockind(1)=Wsys%nval1
   Wsys%nblocks=1
   sym=.true.
   if(.not. stdout)write(stderr,*)'Symmetric matrix filter found'
   
case default
   write(stderr,*)'Unsupported filter type provided: ',Wsys%type
   stop
end select

!check whether matrix is square (obligatory)
if(Wsys%nval1 .ne. Wsys%nval2)then
   write(stderr,*)'ERROR: Filter matrix is not square '
   stop
end if

!read remaining data
call read_BINtype(Wsys)


!retrieve lmax and lmin

do,i=1,Wsys%nint
   select case(Wsys%ints_d(i))
   case('Lmax','lmax','Nmax','nmax')
      
      lmaxW=Wsys%ints(i)
   case('Lmin','lmin','nmin','Nmin')
      lminW=Wsys%ints(i)
   end select
end do




!if not found try to retrieve from side description
if(lmaxW .eq. 0 .or. lminW .eq. 1d6)then
   do,i=1,Wsys%nval1
      read(Wsys%side1_d(i),'(1x,a1,2x,i3.3,i3.3)')trig,l,m
      lmaxW=max(l,lmaxW)
      lminW=min(l,lminW)
   end do
end if

!check with input request

if(lmaxW < lmax)then
   write(stderr,*)'Requested maximum degree may be maximum:', lmaxW
   write(stderr,*)'Output is only up to lmax=:', lmaxW
   
end if

if(lminW > lmin)then
   write(stderr,*)'Requested minimum degree not supported by filter(lminfilter=', lminW,')'
   write(stderr,*)'Filtered only from lmin=', lminW, 'Lower coefficients unchanged'
end if


if(trans .eq. 'T')then
   select case(Wsys%type)
   case('SYMVN___','SYMV0___','SYMV1___','BDSYMV0_','BDSYMVN_')
      write(stderr,*)'WARNING: Transpose option is ignored for the symmetric matrix'
   end select
end if

select case(Wsys%type)
case('BDFULLV0','BDFULLVN','BDSYMV0_','BDSYMVN_')
   !retriev maximum blocksize
   ldsh=0
   side=0
   do,i=1,Wsys%nblocks
      ldsh=max(Wsys%blockind(i)-side,ldsh)
      side=Wsys%blockind(i)
   end do
   allocate(SH(ldsh,nfiles))!maximum size needed for a block
                                     !diagonal approach
case('FULLSQV0','FULLSQVN','SYMVN___','SYMV0___','SYMV1___')
   allocate(SH(Wsys%nval1,nfiles))!maximum size needed for a full matrix approach
end select
      
   

   
   nd=SH_pos(lmaxW,lmaxW)

!note nd > nval1 (also m=0 sine coefficients are stored)
!temporary storage vector which has the sine part stacked on top of the cosine part
   allocate(tmp(nd*2,nfiles))
   tmp=0.d0
   allocate(time(nfiles,3)) !time array holding start center and end times of the datasets

   if(error)then
      allocate(tmp_sig(nd*2,nfiles))
      tmp_sig=0.d0
   end if


!make permutation vector
   allocate(pos(2*nd))
   pos=0
   do,i=1,Wsys%nval1
      read(Wsys%side1_d(i),'(1x,a1,2x,i3.3,i3.3)')trig,l,m
!      write(*,*)i,trig,l,m
      ind=SH_pos(l,m)

      if(trig .eq. 'S')then
         pos(ind+nd)=i
      else
         pos(ind)=i
      end if
   end do  

!fill up the remainder of the permutation vector
   do,i=Wsys%nval1+1,2*nd
      do,j=1,2*nd ! find  a zero spot
         if(pos(j) .eq. 0)then
            pos(j)=i
            exit
         end if
      end do
   end do
   

! read in the reference field
if(rem_rest)then

   call SH_readmeta(filen=trim(ref_sh),type=typ,lmax=lmaxdum)

   allocate(ref_field(2*nd))
   ref_field=0.d0
   write(stderr,*)'Reading reference Field:',trim(ref_sh)
   call SH_readgrav(filen=trim(ref_sh),clm=ref_field(1:nd),slm=ref_field(nd+1:nd*2),type=typ)
   
end if

!now do a loop over the input files to read them in matrix
   do,i=1,nfiles
      if(.not. stdout)write(*,*)'reading ',trim(SHfiles(i))


!read in meta data and spherical harmonic coefficients
      if(stdin)then
         call SH_readmeta(type=typ,lmax=lmaxdum,tcent=time(i,1),tstart=time(i,2),tend=time(i,3))
      else
         call SH_readmeta(filen=trim(SHfiles(i)),type=typ,lmax=lmaxdum,tcent=time(i,1),tstart=time(i,2),tend=time(i,3))
      end if


      !now read in the real file (directly stack the vector with cosine and sine terms)
      if(stdin)then
         if(error)then
            call SH_readgrav(clm=tmp(1:nd,i),slm=tmp(nd+1:nd*2,i),clm_sig=tmp_sig(1:nd,i),slm_sig=tmp_sig(nd+1:nd*2,i),type=typ)
         else
            call SH_readgrav(clm=tmp(1:nd,i),slm=tmp(nd+1:nd*2,i),type=typ)
         end if
      else
         if(error)then
            call SH_readgrav(filen=trim(SHfiles(i)),clm=tmp(1:nd,i),slm=tmp(nd+1:nd*2,i),&
                 clm_sig=tmp_sig(1:nd,i),slm_sig=tmp_sig(nd+1:nd*2,i),type=typ)
         else
            call SH_readgrav(filen=trim(SHfiles(i)),clm=tmp(1:nd,i),slm=tmp(nd+1:nd*2,i),type=typ)
         end if
      end if

      !possibly constrain maximum degree if those coefficients are not available
      if(.not. limdeg .and. lmaxdum < lmaxW)then
         limdeg=.true.
         lmax=lmaxdum
      end if
      
      if(limdeg)then !set coefficients lower and higher than lmin,lmax to zero
         if(lmin > 0)then
            p1=SH_pos(lmin-1,lmin-1)
            write(stderr,*)'Limiting degrees,lmin,lmax'
            tmp(1:p1,i)=0.d0
            tmp(nd+1:nd+p1,i)=0.d0

            if(error)then
               tmp_sig(1:p1,i)=0.d0
               tmp_sig(nd+1:nd+p1,i)=0.d0
            end if
            if(rem_rest)then
               ref_field(1:p1)=0.d0
               ref_field(nd+1:nd+p1)=0.d0
            end if

            !limit maximum degree
            p1=SH_pos(lmax,lmax)+1
            if(p1<=nd)then
               tmp(p1:nd,i)=0.d0
               tmp(nd+p1:2*nd,i)=0.d0
               if(error)then
                  tmp_sig(p1:nd,i)=0.d0
                  tmp_sig(nd+p1:2*nd,i)=0.d0
               end if

               if(rem_rest)then
                  ref_field(p1:nd)=0.d0
                  ref_field(nd+p1:2*nd)=0.d0
               end if
               
            end if
         end if
      end if


      if(rem_rest)tmp(:,i)=tmp(:,i)-ref_field

      

   end do !end do loop over files


   !permute tmp and tmp_sig  tmp->x x(pos)=tmp and tmp_sig->x x(pos)=tmp_sig
   !thus backward permutation
   ! the tmp arrays now will hold the packed parameters(present in matrix system) first and then the other parameters
   call permute_rows(tmp(:,:),pos,.true.)
   if(error)call permute_rows(tmp_sig(:,:),pos,.true.)


!!!!!!!take matrix complement if requested!!
if(compl)then
   select case(Wsys%mtyp)
   case('P') !packed matrix
         Wsys%pack1=-Wsys%pack1 !negate
         !add one to the diagonal
         do,i=1,Wsys%nval1
            ind=packindex(Wsys,i,i)
            Wsys%pack1(ind)=1.d0+Wsys%pack1(ind)
         end do
   case default
         Wsys%mat1=-Wsys%mat1 !negate
         !add one to the diagonal
         do,i=1,Wsys%nval1
            Wsys%mat1(i,i)=1.d0+Wsys%mat1(i,i)
         end do
   end select

end if

!!!!!!!!!!!!!!!!!!!!!!WEIGHTING!!!!!!!!!!!!!!!!!!!
!loop over blocks (actual weighting)
select case(Wsys%type)
case('BDFULLV0','BDFULLVN') ! blockdiagonal matrix
   side=0! keep track of position in the side of the matrix
   full=0! kepp track of position in packed block diagonal matrix
   do,i=1,Wsys%nblocks !in case of the full matrix this loop will have been done only once
      sz=Wsys%blockind(i)-side
      
      !copy data in SH array
      do,k=1,nfiles
         do,j=1,sz
            SH(j,k)=tmp(side+j,k)
         end do
      end do

      !CALL BLAS ( trick BLAS into thinking it gets a sz x sz matrix)
      call dgemm(trans,'N',sz,nfiles,sz,1.d0,Wsys%pack1(full+1)&
           ,sz,SH(1,1),size(SH,1),0.d0,tmp(side+1,1),nd*2)
      
      if(error)then! simple error propagation
         !square entries of tmp_sig and copy in SH array
         do,k=1,nfiles
            do,j=1,sz
               SH(j,k)=tmp_sig(side+j,k)**2
            end do
         end do
         
         !square entries of filter matrix
         do,k=1,sz**2
            Wsys%pack1(full+k)=Wsys%pack1(full+k)**2
         end do
         !call BLAS
         call dgemm(trans,'N',sz,nfiles,sz,1.d0,Wsys%pack1(full+1)&
              ,sz,SH(1,1),size(SH,1),0.d0,tmp_sig(side+1,1),nd*2)
         
         ! square root of the entries
         do,k=1,nfiles
            do,j=1,sz
              tmp_sig(side+j,k)=sqrt(tmp_sig(side+j,k))
            end do
         end do
      end if
      !reset indices
      side=Wsys%blockind(i)
      full=full+sz**2 ! end of the position within the packed matrix
   end do


case('BDSYMV0_','BDSYMVN_') ! blockdiagonal symmetric matrix


   !allocate temporary matrix to hold upper triangle of block
   allocate(blocktmp(ldsh,ldsh))
   side=0! keep track of position in the side of the matrix
   full=0! kepp track of position in packed block diagonal matrix
   do,i=1,Wsys%nblocks !in case of the full matrix this loop will have been done only once

      sz=Wsys%blockind(i)-side

      !copy data in SH array
      do,k=1,nfiles
         do,j=1,sz
            SH(j,k)=tmp(side+j,k)
         end do
      end do

      !copy packed triangle in upper triangle of blocktmp
      do,j=1,sz
         do,k=1,j
            blocktmp(k,j)=Wsys%pack1(full+k+j*(j-1)/2)
         end do
      end do

      !DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      !CALL BLAS (Symmetric matrix multiplication)
      call dsymm('L','U',sz,nfiles,1.d0,blocktmp(1,1),ldsh,SH(1,1),&
           ldsh,0.d0,tmp(side+1,1),nd*2)
      if(error)then! simple error propagation
         !square entries of tmp_sig and copy in SH array
         do,k=1,nfiles
            do,j=1,sz
               SH(j,k)=tmp_sig(side+j,k)**2
            end do
         end do
         
         !square entries of filter matrix
         do,k=1,sz*(sz+1)/2
            Wsys%pack1(full+k)=Wsys%pack1(full+k)**2
         end do
         !call BLAS
         call dsymm('L','U',sz,nfiles,1.d0,blocktmp(1,1),ldsh,SH(1,1),&
              ldsh,0.d0,tmp_sig(side+1,1),nd*2)

         
         ! square root of the entries
         do,k=1,nfiles
            do,j=1,sz
              tmp_sig(side+j,k)=sqrt(tmp_sig(side+j,k))
            end do
         end do
      end if
      !reset indices
      side=Wsys%blockind(i)
      full=full+sz*(sz+1)/2 ! end of the position within the packed matrix
   end do

case('FULLSQV0','FULLSQVN') !full square matrix
   sz=Wsys%nval1
   !copy data in SH array
   do,k=1,nfiles
      do,j=1,sz
         SH(j,k)=tmp(j,k)
      end do
   end do
   
   !point to the right
   call dgemm(trans,'N',sz,nfiles,sz,1.d0,Wsys%mat1(1,1)&
        ,Wsys%nval1,SH(1,1),size(SH,1),0.d0,tmp(1,1),nd*2)
   
   if(error)then! simple error propagation
      !square entries of tmp_sig and copy in SH array
      do,k=1,nfiles
         do,j=1,sz
            SH(j,k)=tmp_sig(j,k)**2
         end do
      end do
      
      !square entries of filter matrix
      do,k=1,sz
         do,j=1,sz
            Wsys%mat1(j,k)=Wsys%mat1(j,k)**2
         end do
      end do
      
      !call BLAS
      call dgemm(trans,'N',sz,nfiles,sz,1.d0,Wsys%mat1(1,1)&
           ,Wsys%nval1,SH(1,1),size(SH,1),0.d0,tmp_sig(1,1),nd*2)
      
      !square root of the results
         do,k=1,nfiles
            do,j=1,sz
              tmp_sig(j,k)=sqrt(tmp_sig(j,k))
            end do
         end do

   end if

case('SYMVN___','SYMV0___','SYMV1___')
   sz=Wsys%nval1
   !copy data in SH array
   do,k=1,nfiles
      do,j=1,sz
         SH(j,k)=tmp(j,k)
      end do
   end do
   !call blas

   call dsymm('L','U',sz,nfiles,1.d0,Wsys%mat1(1,1)&
        ,Wsys%nval1,SH(1,1),size(SH,1),0.d0,tmp(1,1),nd*2)
   
   if(error)then! simple error propagation
      !square entries of tmp_sig and copy in SH array
      do,k=1,nfiles
         do,j=1,sz
            SH(j,k)=tmp_sig(j,k)**2
         end do
      end do
      
      !square entries of filter matrix (upper triangle only)
      do,k=1,sz
         do,j=1,k
            Wsys%mat1(j,k)=Wsys%mat1(j,k)**2
         end do
      end do
      
      !call BLAS
      call dsymm('L','U',sz,nfiles,1.d0,Wsys%mat1(1,1)&
           ,Wsys%nval1,SH(1,1),size(SH,1),0.d0,tmp_sig(1,1),nd*2)
      
      !square root of the results
      do,k=1,nfiles
         do,j=1,sz
            tmp_sig(j,k)=sqrt(tmp_sig(j,k))
         end do
      end do

   end if
end select

!!!!!!!!!!!!!!END WEIGHTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!loop over output files
!   tmp=0.d0
!   if(error)tmp_sig=0.d0

! permute back values of tmp and tmp_sig to original scheme
call permute_rows(tmp(:,:),pos,.false.)
if(error)call permute_rows(tmp_sig(:,:),pos,.false.)

   do,i=1,nfiles

      if(rem_rest)tmp(:,i)=tmp(:,i)+ref_field ! restore reference field

      if(limdeg)then
         p1=SH_pos(lmax,lmax)
      else
         p1=nd
      end if
!Write the coefficients back to a file or standard output
     
      if(stdout)then
         if(error)then
            call SH_write(clm=tmp(1:p1,i),slm=tmp(nd+1:nd+p1,i),clm_sig=tmp_sig(1:p1,i),&
                 slm_sig=tmp_sig(nd+1:nd+p1,i),tcent=time(i,1),tstart=time(i,2),tend=time(i,3),typ=4)
         else
            call SH_write(clm=tmp(1:p1,i),slm=tmp(nd+1:nd+p1,i),tcent=time(i,1),tstart=time(i,2),tend=time(i,3),typ=4)
         end if
      else
         if(replsuf)then
            ind=index(SHfiles(i),'.',.true.)
            dum=SHfiles(i)(1:ind-1)//trim(app)
         else
            dum=trim(SHfiles(i))//trim(app)
         end if
         !trim leading directories
         ind=index(dum,'/',.true.)
         dum=trim(outdir)//dum(ind+1:);
         if(error)then
            call SH_write(clm=tmp(1:p1,i),slm=tmp(nd+1:nd+p1,i),clm_sig=tmp_sig(1:p1,i),&
                 slm_sig=tmp_sig(nd+1:nd+p1,i),tcent=time(i,1),tstart=time(i,2),tend=time(i,3)&
                 ,filen=trim(dum),typ=4)
         else
            call SH_write(clm=tmp(1:p1,i),slm=tmp(nd+1:nd+p1,i),tcent=time(i,1),tstart=time(i,2),tend=time(i,3)&
             ,filen=trim(dum),typ=4)
         end if
      end if
      

   end do
   


   

 end program SH_filter2


subroutine help()
  implicit none
  integer::unit
  unit=0
  write(unit,*)'Program which filters spherical harmonic coefficients through a matrix SH_filt= WMAT*SHold'
  write(unit,*)'Filter matrix may be symmetric, fully populated or block diagonal'
  write(unit,*)'The output are the files with .filt appended.'
  write(unit,*)''
  write(unit,*)'usage SH_filter [OPTIONS] [SHFILES]'
  write(unit,*)'The outputfiles are the input files with a suffix appended (filt) and leading directories stripped'
  write(unit,*)''
  write(unit,*)'OPTIONS can be the following:'
  write(unit,*)'   -W FILTERFILE: use the weights in the filter file'
  write(unit,*)'	  The filter file may contain either a block diagonal filter matrix, a full or symmetric matrix'
  write(unit,*)'   -T : transpose filter matrix'
  write(unit,*)'   -s: write results to standard output'
  write(unit,*)'   -lLMAX,LMIN: limit the input to these degrees (sets coefficients with l<lmin and l>lmax to zero)'
  write(unit,*)'   -e: also propagate errors (diagonals from the original coefficient covariance matrix)'
  write(unit,*)'   -a|rAPP: Append|replace extension (.APP) to the output (APP may be maximum 19 characters)'
  write(unit,*)'   -R REFFIELD: apply a remove and restore action in the filtering (Using reference field REFFIELD'
  write(unit,*)'   -c Use the complement (I-W) of the filter matrix (handy for mask operations)'
  write(unit,*)'   -o OUTDIR  put output files in the directory OUTDIR (default is current directory)'
  
  stop
end subroutine help
 

!!subroutine to retrieve data from a sinex file
!!most of the arguments are optional so onloy necessary ifo is retrieved
!!For efficiency the retrieval methods are programmed according to the order they occur in the sinex file.
!!Coded by Roelof Rietbroek, Tue Aug 14 12:02:36 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany 
!!email: roelof@gfz-potsdam.de
!!
!! first version retrieves SOLUTION/ESTIMATE block, solution SOLUTION/MATRIX_ESTIMATE
!!upon first call the blocks are read in and stored
!!upon subsequent(or first call) calls search queries can be performed to retrieve all or a selected portion of the data

!!Updated by Roelof Rietbroek, Wed Oct 24 11:51:28 2007
!!covariance matrix is now stored in packed form (half the size) and can be retrieved in packed format
!!Updated by Roelof Rietbroek, Tue March 21 2017
!! addedd space in between stop statement to comply with new fortran compiler



subroutine get_sinex(file,nx,ny,x,y,covxy,covxypack,covdiagx,inc1,excl1,inc2,excl2,xlist,gpswk)
use GPStlbx,only: compstr,GPS_week,sinex_date,GPS_year
implicit none
character(*),intent(in)::file !obligatory argument
!optional arguments
character(*),intent(in),optional::inc1,excl1,inc2,excl2 
integer,intent(out),optional::nx,ny,gpswk
double precision,intent(out),optional,dimension(:)::x,y,covdiagx,covxypack
double precision,intent(out),optional,dimension(:,:)::covxy
character(*),intent(out),optional,dimension(:)::xlist

integer::i,j,k,n1,n2,stderr,funit,nmax,nobs,ind,ind2,npos,rdim1,rdim2
parameter(nmax=1000000)!enough(mximum amount of lines in a SINEX file)
character(120)::fileold
character(80)::dum,mattyp
double precision,allocatable,dimension(:)::solmat
double precision, allocatable, dimension(:)::solest,sigsolest
double precision::version,dumvec(3),rem1,rem2
character(47),allocatable,dimension(:)::list
character(12)::date1,date2,date3
logical::fex,upper,cov=.false.
integer,allocatable,dimension(:)::ivec1,ivec2
save fileold,solest,sigsolest,solmat,list,nobs,cov

!defaults
stderr=0
funit=13
mattyp=''
upper=.false.


if(trim(fileold) .ne. trim(file))then !check whether file has been called in previous call else read in data
   fileold=trim(file)

   !check if file exists
   inquire(FILE=trim(file),EXIST=fex)
   if(.not. fex) then
      write(stderr,*) 'file '//trim(file)// ' does not exist'
      stop
   end if
!open sinex file
   open(unit=funit,file=trim(file),status='old')
   write(*,*)'reading sinex file ',trim(file)
!read first line with meta data
   read(funit,'(a5,1x,f4.2,5x,a12,5x,a12,1x,a12,3x,I5.5)')dum(1:5),version,date1,date2,date3,nobs
   if (dum(1:5) .ne. '%=SNX')stop 'does not appear to be a sinex file'

   !deallocate if still allocated from the previous call
   if(allocated(solest))deallocate(solest,sigsolest,list)
!allocate arrays
   allocate(solest(nobs),sigsolest(nobs),list(nobs))
   
   !initialize variables
   solest=0.d0
   sigsolest=0.d0
   list=''



   !now read in solest array from the SOLUTION/ESTIMATE block
   do,i=1,nmax !loop tio search for solution/estimate block
      read(funit,'(A80)')dum
      !write(*,*)dum
      if(index(dum,'+SOLUTION/ESTIMATE') .ne. 0)then !start solution estimate block
         write(*,*)'reading SOLUTION/ESTIMATE block'
         ind=0
         do,j=1,nmax !loop within block
            read(funit,'(a80)')dum
            if(dum(1:1).eq.'*')cycle !when a comment is encountered
            if(index(dum,'-SOLUTION/ESTIMATE') .ne. 0)exit !exit block when block is closed
           
            ind=ind+1
            read(dum,'(a47,E21.15,1x,e11.6)')list(ind),solest(ind),sigsolest(ind)
            !write(*,*)ind,trim(list(ind)),solest(ind),sigsolest(ind),nobs

         end do
         exit !search loop after block has been found and  read in
      else if(index(dum,'%ENDSNX') .ne. 0)then !exit loop because end of file
      write(stderr,*)'Block SOLUTION/ESTIMATE is not found'
         exit
      end if
      
   end do
   
   !now look for a covariance matrix in the SOLUTION/MATRIX_ESTIMATE block
   
   do,i=1,nmax
      read(funit,'(a80)')dum
      if(index(dum,'+SOLUTION/MATRIX_ESTIMATE') .ne. 0)then !start solution estimate block
         
         cov=.true.
         !retrieve matrix type (eg L COVA)
         mattyp=dum((index(dum,'MATE')+4):)
         ind=0
         ind2=0
        write(*,*)'reading SOLUTION/MATRIX_ESTIMATE block ',trim(mattyp)
        if(index(mattyp,'U') .ne. 0)upper=.true. !upper triangle is given
         !reallocate solmat and initialize

         if(allocated(solmat))deallocate(solmat)
         allocate(solmat(nobs+((nobs-1)*(nobs))/2))
         solmat=0.d0

         do,j=1,nmax !loop within block
            read(funit,'(a80)')dum
            if(dum(1:1).eq.'*')cycle
            if(index(dum,'-SOLUTION/MATRIX_ESTIMATE').ne. 0)exit !stop solution block


               read(dum,'(1x,i5,1x,i5,1x,E21.14,1x,E21.14,1x,E21.14)')ind,ind2,dumvec
              
               if(upper)then !upper triangular matrix
                  solmat(ind+((ind2-1)*(ind2))/2)=dumvec(1)
                  if(ind2 <= (ind-2))then
                     solmat(ind+((ind2)*(ind2+1))/2)=dumvec(2)
                     solmat(ind+((ind2+1)*(ind2+2))/2)=dumvec(3)
                  else if(ind2 <= (ind-1))then
                     solmat(ind+((ind2)*(ind2+1))/2)=dumvec(2)
                  
                  end if
               else !transpose lower triangular matrix in upper triangular matrix
                  solmat(ind2+((ind-1)*(ind))/2)=dumvec(1)
                  if(ind2 <= (ind-2))then
                     solmat(ind2+1+((ind-1)*(ind))/2)=dumvec(2)
                     solmat(ind2+2+((ind-1)*(ind))/2)=dumvec(3)
                  else if(ind2 <= (ind-1))then
                     solmat(ind2+1+((ind-1)*(ind))/2)=dumvec(2)
                  
                  end if
!                  write(stderr,*)ind, ind2,size(solmat,1),solmat(ind2+((ind-1)*(ind))/2)&
!                       ,solmat(ind2+((ind)*(ind+1))/2),solmat(ind2+((ind+1)*(ind+2))/2)
               end if

            
            !write(*,*)ind,ind2,solmat(ind,ind2),dumvec,j
         end do
         exit !exit loop after block has been read in
      else if(index(dum,'%ENDSNX') .ne. 0)then !exit loop because end of file
      write(stderr,*)'SOLUTION/MATRIX_ESTIMATE block not found'
         exit
      end if
   end do

!old loops concerned with non-packed matrix (obsolete) 
    !now put in standard deviations on the diagonal and mirror matrix (at least if a covariance block is found)
   ! if(mattyp .ne. '')then
!       if(index(mattyp,'U') .ne. 0)then !upper triangkle is given
!          do,j=1,nobs
!             do,k=1,j-1
!                solmat(j,k)=solmat(k,j)
!             end do
!             !solmat(j,j)=sigsolest(j)**2!plug in standard dev
!          end do
!       else !lower triangle is given
!          do,j=1,nobs
!             do,k=1,j-1
!                solmat(k,j)=solmat(j,k)
!             end do
!            ! solmat(j,j)=sigsolest(j)**2!plug in standard dev
!          end do
         
!       end if
!    end if
   
   close(funit)

end if !reading in of the file for the first time


!now if the subrotuine has been called for the first time all the necessary data has been read in
!below the optional arguments are elaborated and the required data is extracted
!optional arguments in the routine search in the 



!Now check if there are search restrictions for the first dimension and set the restriction type rdim1=
!3 for both including and excluding search
!2 for only including search
!1 for only excluding search
!0 no restriction

!write(*,*)'test'
!allocate index arrays 
allocate(ivec1(nobs),ivec2(nobs))


if(present(inc1) .and. present(excl1))then
   rdim1=3
else if (present(inc1))then
   rdim1=2
else if (present(excl1))then
   rdim1=1
else
   rdim1=0
end if

!do the same for the second dimension
if(present(inc2) .and. present(excl2))then
   rdim2=3
else if (present(inc2))then
   rdim2=2
else if (present(excl2))then
   rdim2=1
else
   rdim2=0
end if


!now get the number of variables satisfying the restrictions 

!first dimension
select case(rdim1)
case(0)
   n1=nobs
  
   ivec1=(/ (i, i = 1,nobs) /)
  
case(3)!both restrictions
   n1=0
   do,i=1,nobs
      if(.not.compstr(list(i),inc1))cycle !cycle when inc1 string is NOT found
      if( compstr(list(i),excl1))cycle !cycle when excl1 string is found
      n1=n1+1
      ivec1(n1)=i !set index string
   end do
case(2)
   n1=0
   do,i=1,nobs
      if(.not.compstr(list(i),inc1))cycle !cycle when inc1 string is NOT found
      n1=n1+1
      ivec1(n1)=i !set index string
!      write(*,*)list(i)
   end do
case(1)
   n1=0
   do,i=1,nobs
      if( compstr(list(i),excl1))cycle !cycle when excl1 string is found
      n1=n1+1
      ivec1(n1)=i !set index string
   end do
end select
   

!second dimension

select case(rdim2)
case(0)
   n2=nobs
   ivec2=(/ (i, i = 1,nobs) /)
case(3)!both restrictions
   n2=0
   do,i=1,nobs
      if(.not.compstr(list(i),inc2))cycle !cycle when inc2 string is NOT found
      if( compstr(list(i),excl2))cycle !cycle when excl2 string is found
      n2=n2+1
      ivec2(n2)=i !set index string
   end do
case(2)
   n2=0
   do,i=1,nobs
      if(.not.compstr(list(i),inc2))cycle !cycle when inc2 string is NOT found
      n2=n2+1
      ivec2(n2)=i !set index string
   end do
case(1)
   n2=0
   do,i=1,nobs
      if( compstr(list(i),excl2))cycle !cycle when excl2 string is found
      n2=n2+1
      ivec2(n2)=i !set index string
   end do
end select

!now we're ready to substitue the values in the required input arguments

!indices
if(present(nx))nx=n1
if(present(ny))ny=n2


!input check (no packed matrix can be requested if n1 .ne. n2)
if ( (n1 .ne. n2) .and. present(covxypack))then
   write(stderr,*)' Packed matrix can not be provided for nx .ne. ny'
   stop
end if

!1dimensional data vectors

if(present(x))then
   if(size(x,1)<n1)then
      write(stderr,*) 'vector x has dimension ',size(x,1),' required ',n1
      stop
   else

      x(1:n1)=solest(ivec1(1:n1))
   end if
end if

if(present(y))then
   if(size(y,1)<n2)then
      write(stderr,*) 'vector y has dimension ',size(y,1),' required ',n2
      stop
   else
      y(1:n2)=solest(ivec2(1:n2))
   end if
end if

!covariance matrix (square form)

if((.not.cov) .and. present(covxy))then
!input check
   if(size(covxy,1) < n1 .or. size(covxy,2) < n2)then
      write(stderr,*) 'vector covxypack has dimension ',size(covxy,1),'x',size(covxy,2),' required ',n1,'x',n2
      stop
   end if
 !  write(*,*)cov,present(covxy)
   write(stderr,*)'WARNING: no full covariance matrix specified providing diagonal only'
   covxy=0.d0
   do,i=1,n1
      do,j=1,n2
         if(ivec1(i) .eq. ivec2(j))then
            covxy(i,j)=sigsolest(ivec1(i))**2
         end if
      end do
   end do
else if(present(covxy))then
   !loop over entries
   do,i=1,n1
      do,j=i,n2
!         covxy(1:n1,1:n2)=solmat(ivec1(1:n1),ivec2(1:n2))
         covxy(i,j)=solmat(ivec1(i)+(ivec2(j)-1)*(ivec2(j))/2)
         covxy(j,i)=covxy(i,j)
      end do
   end do
end if

!covariance matrix (packed form upper triangular)

if((.not.cov) .and. present(covxypack))then
!input check
   if(size(covxypack,1) < (n1+((n1-1)*n1)/2))then
      write(stderr,*) 'vector covxypack has dimension ',size(covxypack,1),' required ',(n1+((n1-1)*n1)/2)
      stop
   end if

   write(stderr,*)'WARNING: no full covariance matrix specified providing diagonal only'
   covxypack=0.d0
   do,i=1,n1
      do,j=1,n2
         if(ivec1(i) .eq. ivec2(j))then
            covxypack(i+((j-1)*j)/2)=sigsolest(ivec1(i))**2
         end if
      end do
   end do
else if(present(covxypack))then
   !loop over entries
   do,i=1,n1
      do,j=i,n2
!         covxy(1:n1,1:n2)=solmat(ivec1(1:n1),ivec2(1:n2))
         covxypack(i+((j-1)*j)/2)=solmat(ivec1(i)+(ivec2(j)-1)*(ivec2(j))/2)

      end do
   end do
end if


!diagoanl of covarince matrix only
if(present(covdiagx))then
if(size(covdiagx,1) < n1)then
      write(stderr,*) 'vector covdiagx has dimension ',size(covdiagx,1),' required ',n1
      stop
   else
      covdiagx(1:n1)=sigsolest(ivec1(1:n1))
   end if   

end if

!return list of data points corresponding to the x restictions

if(present(xlist))then
   if(size(xlist,1)<n1)then
      write(stderr,*) 'vector xlist has dimension ',size(xlist,1),' required ',n1
      stop
   else
      if(len(xlist(1))<80)stop 'xlist characters width must be at least 80'

      xlist(1:n1)=list(ivec1(1:n1))
   end if
end if

!return reference gpsweek number if requested (average of the dates in header line)
if(present(gpswk))then
   gpswk=GPS_week(sinex_date((GPS_year(GPS_week(date2,rem1))+GPS_year(GPS_week(date3,rem2)))/2.d0))
   gpswk=int((rem1+rem2)/2.d0)
end if

end subroutine get_sinex

!function which checks whether a string is present in another
!wild character * allowed and the | (or) sign checks the substrings regardless of there position
!returns true when string is found
function compstr(str,srch)
implicit none
integer::nmax
parameter(nmax=5000)
logical::compstr,first,last
character(*),intent(in)::str,srch
integer::st,st2,nd,nd2,test,test2,test3,start,i,j,ln,wild1,wild2
compstr=.false.



st=1
nd=len(srch)
st2=1
ln=len(str)

do,i=1,nmax
   if(compstr)exit !keep trying when no match occured in previous run
   !search for a | character
   test=index(srch(st:),'|')
   if(test .ne. 0)then
      nd=test-2+st
   else
      nd=len(srch)
   end if
   
   !now search for a * character within the substring srch(st:nd)
   st2=st
   nd2=nd
   start=1
!  write(*,*)'first search string ',srch(st:nd),i

!check whether first and last entry are a * or not
    wild1=index(srch(st2:nd),'*')
    wild2=index(srch(1:nd),'*',.true.)
   
   !four scenarios: *term*, *term, term*, term
   if((wild1 .eq. 1) .and. (wild2 .eq. nd))then !take of the wild characters of both sides and check for a fit of first and last term
      st2=st2+1
      nd=nd-1
      first=.false.
      last=.false.
   else if( wild1 .eq. 1)then !strip * of the front and check whether the last term fits
      st2=st2+1
      first=.false.
      last=.true.
   else if (wild2 .eq. nd)then !strip of the * at the back and check if the first term fits
      nd=nd-1
      last=.false.
      first=.true.
   else
      first=.true.
      last =.true.
   end if

   
   
   do,j=1,nmax

      test2=index(srch(st2:nd),"*")
!      write(*,*)st2,nd,' ',test2
      
      if(test2 .ne. 0)then
!         write(*,*)srch(st2:nd2)
         nd2=test2-2+st2
         
         
      else !last search part encountered
         nd2=nd
      end if
      
      !write(*,*)'test arguments:',start,st2,nd2
      !write(*,*)'arg1:',trim(str(start:)),' arg2:',trim(srch(st2:nd2))
      test3=index(str(start:),srch(st2:nd2))
      if(test3 .ne. 0)then
         compstr=.true.
      else
         compstr=.false.
         exit
      end if

      if(first)then !check whether first entry without
         if(test3 .ne. 1)then
            compstr=.false.
            exit
         end if
         first=.false.
      end if
      
      if(nd2 .eq. nd)then
   !      if(last)write(*,*)test3,ln-start-nd2+st2+1
         if(last .and. (test3 .ne. (ln-start-nd2+st2+1)))then
            compstr=.false.
            
         end if
         exit
      end if

      start=test3+start+nd2-st2+1
      st2=test2+st2




      
   end do
  



   st=test+st !reset start index
   if(test .eq. 0)exit !exit loop after last | found
end do


end function compstr


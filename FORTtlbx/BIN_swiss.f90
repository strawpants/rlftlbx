!!program which outputs information on a binary file made by write_sym_mat
!!

!!Coded by Roelof Rietbroek, Fri Nov  9 10:19:18 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de
!!TODO:
!!1: do some statistics on the data (eg trace, condition number etc)
!!2: Allow some data tweaking (changing header data)

!!Updated by Roelof Rietbroek, Wed Jan  2 14:16:08 2008
!!version 2 uses the version 2 matrices  (different file structure)

!!Updated by Roelof Rietbroek, Thu Jul 24 12:24:04 2008
!! new version now reads in the matrices using version 2.2 routines (using derived type)

!!Updated by Roelof Rietbroek, Mon Aug  4 10:10:13 2008
!! now also reads from standard input
!!Updated by Roelof Rietbroek, Wed Nov 12 10:38:57 2008
!! now prints maximal 6 values per line for matrix values
!!Updated by Roelof Rietbroek, Tue Jan 27 10:34:28 2009
!!also print a specific column of a matrix
!!Updated by Roelof Rietbroek, Thu Aug  6 10:22:35 2009
!! added support for large files (large index integers)

!!Updated by Roelof Rietbroek, Thu Oct 29 14:13:23 2009
!!added diagonal matrix

!!Updated by Roelof Rietbroek, Tue Dec  8 11:40:18 2009
!! converted BINfileinfo program into BIN_swiss ( swiss army tool)

!!Updated by Roelof Rietbroek, Mon Mar  1 15:07:10 2010
!! added support for longer regular expressions

!!Updated by Roelof Rietbroek, Fri Apr 16 14:09:47 2010
!!added support for reading regular expressions from file

!!Updated by Roelof Rietbroek, Wed Jun 16 16:34:42 2010
!! allow plain scaling of matrix

!!Updated by Roelof Rietbroek, Wed Jul 21 08:53:52 2010
!!possibly restict output based on numerical vector criteria

!!Updated by Roelof Rietbroek, Wed Jul 28 13:44:13 2010
!!allowed tagging of parameters at the end

!!Updated by Roelof Rietbroek, Sun Aug 29 15:16:25 2010
!!explicitly open standard output with a longer record length

!!Updated by Roelof Rietbroek, Tue Oct  5 13:52:09 2010
!! allow sides to be sorted/permuted with a predefined parameter sequence

!!Updated by Roelof Rietbroek, Wed Apr 27 12:41:39 2011
!!also print an end of header tag when vector/matrix data is printed ( this aids the reading of SH coefficients by using pipes)
!! only prints this when also a header is requested

!!Updated by Roelof Rietbroek, Thu Jun 16 09:38:15 2011
!! allow restrictions on matrix columns

!!Updated by Roelof Rietbroek, Wed Aug 22 11:32:57 2012
!! added a check for illegal diagonal output of a non-square matrix

!!Updated by Roelof Rietbroek, Mon Aug 27 09:35:54 2012
!! allow to change the type of the matrix ( under certain conditions)

!!Updated by Roelof Rietbroek, Thu Sep 13 13:58:55 2012
!! fixed small bug ( output format of floats differed when using transpose option)

!!Updated by Roelof Rietbroek, Fri Sep 21 10:37:26 2012
!! added option to export vector part as 2d matrix

!!Updated by Roelof Rietbroek, Sat Mar  9 16:42:39 2013
!! fixed bug in parse_n_print (no spaces when negative numbers where present)

!!Updated by Roelof Rietbroek, Mon Apr 29 17:27:32 2013
!! produce an error message when one of the sides turns out to have zero length

!!Updated by Roelof Rietbroek, Tue Aug 13 15:08:36 2013
!! added option to left adjust side description

!!Updated by Roelof Rietbroek, Mon Aug 26 21:31:18 2013
!! increased length of dum parameter



program BIN_swiss
use binfiletools
use bin_operations
use forttlbx
implicit none
integer::stderr,i,j,itharg,narg,get,st
type(BINdat)::dat
integer::skip,sz,stage,ind
logical::header,root,diag,transp
logical::binary,ascii,dryrun
!defaults initializations
character(200)::parse ! parse string defining the info which is printed to screen
character(504)::dum
character(500),pointer::metastr(:)
character(6)::dum2
integer::iargc
integer::st1,st2,nd1,nd2,mn
integer::reg_inc1,reg_inc2 !regular expression ids
integer,allocatable::perm1(:),perm2(:)
logical::expand,matpermute,excl1,excl2,match
double precision::matscale,vecscale
character(24),pointer::restrict(:)
integer:: nrestr
double precision::tmp
integer,allocatable::keep1(:),keep2(:)
integer::ind1,ind2,col
double precision::mnval,mxval
integer::stdout
integer::unit,nsort,last,memsz,ncom,chunk
parameter(chunk=100)
logical::sort1,sort2
character(24),pointer::sort_d(:),tmpvec_d(:),newside_d(:)
character(200)::sort1f,sort2f
integer,allocatable,dimension(:)::tmpperm1,tmpperm2
character(8)::outtype
logical::v2mat
integer::indold,sz2
!default initializations

!standard error unit
unit=13
sort1=.false.
sort2=.false.
memsz=0
sort1f=''
sort2f=''
matscale=1.d0
vecscale=1.d0
stderr=0
stdout=6
parse=''
header=.true.
root=.false.
diag=.false.
transp=.false.
excl1=.false.
excl2=.false.
dat%file='stdin' !default reads file from standard input
reg_inc1=0
reg_inc2=0
binary=.false. 
ascii=.false.
expand=.false.
matpermute=.false.
stage=2 ! defsult reads only meta info from the file
mn=0 ! amount of metadatat update commands
!get command line arguments
nrestr=0 ! amount of restrict commands
narg=iargc()
st1=0
st2=0
outtype=''
v2mat=.false.
!if (narg < 1)call help()


!loop over remaining arguments (options)

itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit !exit loop when last argument has been encountered
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('s')!write solution to screen (plus square of diagonal of matrix)
         root=.true.
         stage=4
         ascii=.true.
      case('n') ! don't print header
         header=.false.
      case('d')!print diagonal of matrix (no square root taken)
         diag=.true.
         stage=4
         ascii=.true.
      case('p') !custom parse string
         parse=trim(dum(4:204))
         parse=trim(parse)//',' !! append , for emptying buffer
         ! check from parse string what data needs to be read out
         if(index(parse,'v') .ne. 0)stage=3 ! vector data
         if(index(parse,'m') .ne. 0)stage=4 !full matrix data
         if(index(parse,'d') .ne. 0)stage=4 ! diagonal of the matrix
         if(index(parse,'r') .ne. 0)stage=4 !root of the diagonal
         if(index(parse,'c') .ne. 0)stage=4 !column of the matrix
         ascii=.true.
      case('t')!transpose matrix output no vector printing is allowed!
         transp=.true.
      case('S')! apply additional scales
         select case(dum(3:4))
         case('m=')
            read(dum(5:),*)matscale
         case('v=')
            read(dum(5:),*)vecscale
         case default
            write(stderr,*)'Error processing option:',trim(dum)
            stop
         end select
      case('R')!restrict based on numerical range
         nrestr=nrestr+1
         call realloc_ptr(restrict,1)
         restrict(nrestr)=trim(dum(4:27))
         matpermute=.true.
      case('M')!append metadata operations
         mn=mn+1
         call realloc_ptr(metastr,1)
         metastr(mn)=trim(dum(3:))
      case('P') ! permute sides 1 and/or 2
         matpermute=.true.
         select case(dum(3:4))
         case('1:')
            sort1=.true.
            sort1f=dum(5:)
         case('2:')
            sort2=.true.
            sort2f=dum(5:)
         case default
            write(stderr,*)'Error processing option:',trim(dum)
            stop
         end select
         
      case('I','E')!include/exclude regular expression
         select case(dum(3:4))
         case('1=')
            if(reg_inc1 .ne. 0)then
               write(stderr,*)'Only one occasion of -I1 -E1 is allowed'
               stop
            end if
            reg_inc1=regcompf(trim(dum(5:)))
         case('1:')
            if(reg_inc1 .ne. 0)then
               write(stderr,*)'Only one occasion of -I1 -E1 is allowed'
               stop
            end if
            reg_inc1=regcompf(trim(dum(5:)),.true.)
         case('2=')
            if(reg_inc2 .ne. 0)then
               write(stderr,*)'Only one occasion of -I2 -E2 is allowed'
               stop
            end if
            reg_inc2=regcompf(trim(dum(5:)))
         case('2:')
            if(reg_inc2 .ne. 0)then
               write(stderr,*)'Only one occasion of -I2 -E2 is allowed'
               stop
            end if
            reg_inc2=regcompf(trim(dum(5:)),.true.)
         case default
            write(stderr,*)'Option is incorrect ',trim(dum)
            stop
         end select

         select case(dum(2:3))
         case('E1')
            excl1=.true.
         case('E2')
            excl2=.true.
         end select

         matpermute=.true.
      case('e')!expand sparse matrix
         expand=.true.
      case('c') ! change the matrix type
         select case(dum(3:3))
         case('s')
            outtype='SYMVN___'
            case default 
               write(stderr,*)'ERROR: processing option:',trim(dum)
         end select
      case('b')!binary output
         binary=.true.
         header=.false. ! set header output to false
         stage=4
      case('V')! export vectors as matrix
         v2mat=.true.
         select case(dum(3:3))
         case('=')!read tags from command line
            indold=3
            do while(indold < len_trim(dum))
               !get index of backslash
               ind=index(dum(indold+1:),'/')
               if(ind .eq. 0)then
                  ind=len_trim(dum)+1 ! or end of command line argument
               else
                  ind=ind+indold !index from start of the string
               end if
               sz2=sz2+1
               call realloc_ptr(newside_d,1) ! add element to description
               newside_d(sz2)=dum(indold+1:ind-1)
               indold=ind!reset starting index
            end do
         case(':')!read tags from file
            open(unit=unit,file=trim(dum(4:)),form='formatted',status='old')
            last=0
            do while(last .eq. 0)
               read(unit=unit,iostat=last,fmt='(A24)')dum
               if(last .ne. 0)exit
               sz2=sz2+1
               call realloc_ptr(newside_d,1)
               newside_d(sz2)=dum(1:24)
              
            end do
            close(unit)
         case default
            write(stderr,*)'ERROR: processing -V option:',dum(1:3)
            stop     
         end select
      case('h')
         call help()
      case default
         write(stderr,*)'unknown option selected:',trim(dum)
         write(stderr,*)''
         call help()
      end select
      else !not an option but file name
         
         dat%file=trim(dum(1:))
      end if
end do

!check for bad requests
if(root .and. diag)then
   write(stderr,*)"ERROR:cannot request -d and -s at the same time"
   stop
end if

if(transp .and. binary)then
   write(stderr,*)"ERROR: -t (transpose) option cannot be combined with binary output"
   stop
end if

if(binary .and. ascii)then
   write(stderr,*)"ERROR:cannot request binary and ascii output at the same time"
   stop
end if

if((ascii .and. stage >3) .and. matpermute)then
   write(stderr,*)"ERROR:cannot permute matrix and write its values in ascii at the same time"
   write(stderr,*)"      You can pipe the permuted output (-b) in a new instance of BIN_swiss"
   stop
end if


if(v2mat .and. .not. binary)then
   write(stderr,*)"ERROR: -V option is only allowed in combinatiojn with the -b option"
   stop
end if

!read neccessary segments

call read_BINtype(dat,2)

if(matpermute .or. outtype .ne. '')then !read some matrices in full format
   select case(dat%type)
   case('SYMV0___','SYMVN___','SYMV1___','SYMV2___')
      dat%mtyp='U'
   case('FULLSQV0','FULLSQVN','FULL2DVN')
      dat%mtyp='F'
   case default
      dat%mtyp='P'
   end select
end if

if(stage >2)call read_BINtype(dat,stage) ! read remaining data

!update/append meta info if requested
if(mn > 0)then
   do,i=1,mn
      call BIN_putmeta(dat,metastr(i))
   end do
end if

!expand sparse matrices if requested
if(expand .and. stage > 3 )then
   call BIN_expand(dat)
end if

if(reg_inc1 >0 .or. nrestr >0)then
   allocate(keep1(dat%nval1))
   keep1=0 ! default is to keep all ( 0: keep 1 :remove)
end if

if(nrestr >0)then ! search for rows to exclude based on numerical bounds
   do,j=1,nrestr ! loop over all criteria
      ind1=index(restrict(j),',')
      ind2=index(restrict(j),',',.true.)
      if(ind1 >= ind2)then
         write(stderr,*)"ERROR processing restriction string:",restrict(j)
         stop
      end if
      read(restrict(j)(1:ind1-1),*)mnval
      read(restrict(j)(ind2+1:),*)mxval
      select case(restrict(j)(ind1+1:ind1+1))
      case('v')!rectrict based on vector criteria
         read(restrict(j)(ind1+2:ind1+2),*)col
         if(col > dat%nvec)then
            write(stderr,*)"ERROR:requested vector column does not exist"
            stop
         end if

         do,i=1,dat%nval1
            if(dat%vec(i,col)<mnval .or.dat%vec(i,col)>mxval)keep1(i)=1
         end do
      case('c')
         read(restrict(j)(ind1+2:ind1+2),*)col
         if(col > dat%nval2)then
            write(stderr,*)"ERROR:requested matrix column does not exist"
            stop
         end if
         if(dat%mtyp .eq. 'P')then
            do,i=1,dat%nval1
               tmp=dat%pack1(packindex(dat,i,col))
               if(tmp<mnval .or. tmp>mxval)keep1(i)=1
            end do
         else
            do,i=1,dat%nval1
               tmp=dat%mat1(i,col)
               if(tmp<mnval .or. tmp>mxval)keep1(i)=1
            end do
         end if
      case('d')!restrict based on diagonal info or matrix column
         if(dat%mtyp .eq. 'P')then
            do,i=1,dat%nval1
               tmp=dat%pack1(packindex(dat,i,i))
               if(tmp <mnval .or. tmp >mxval)keep1(i)=1
            end do
         else
            do,i=1,dat%nval1
               tmp=dat%mat1(i,i)
               if(tmp <mnval .or. tmp >mxval)keep1(i)=1
            end do
         end if
         case default
            write(stderr,*)"ERROR:processing restrict string: ",restrict(j)
            stop
      end select
   end do
end if

if(reg_inc1>0)then ! make permutation vector for side1
   do,i=1,dat%nval1
      if(keep1(i) .eq. 1) cycle ! no need to check, already excluded
      match=regexecf(reg_inc1,dat%side1_d(i))
      if(excl1)match=.not. match
      if(.not. match)keep1(i)=1
   end do
end if




!create permutation vector
if(reg_inc1 + nrestr >0)then
   !check for zero length output
   if(sum(keep1) .eq. dat%nval1)then
      write(stderr,*)"ERROR: side 1 has zero length"
      stop
   end if

   allocate(perm1(dat%nval1))
   call sort_f90(perm1,keep1)
   st1=dat%nval1-sum(keep1)
else if(matpermute)then
   st1=dat%nval1
   allocate(perm1(st1))
   forall(i=1:st1)perm1(i)=i
end if

!optionally sort according to predefined vector
if(sort1)then
   !read data from file
   open(unit=unit,file=sort1f,form='formatted')
   last=0
   nsort=0
   allocate(sort_d(chunk))
   sort_d=''
   memsz=chunk
   do while(last .eq. 0)

      if(nsort+1 > memsz)then
         call realloc_ptr(sort_d,chunk)
         memsz=memsz+chunk
      end if
      read(unit=unit,fmt='(A24)',iostat=last)sort_d(nsort+1)
      if(last .ne. 0)exit
      nsort=nsort+1
   end do
   close(unit)

   allocate(tmpvec_d(st1),tmpperm1(st1),tmpperm2(nsort))
   tmpvec_d=dat%side1_d(perm1(1:st1))
   
!calculate the permutation vector of the non-excluded part
   call get_permvec2(tmpvec_d,sort_d(1:nsort),.true.,1,24,tmpperm1,tmpperm2,ncom)

   !modify the already existing permutation vector
   perm1(1:st1)=perm1(tmpperm1(1:st1))
   deallocate(tmpvec_d,tmpperm1,tmpperm2,sort_d)
end if


if(reg_inc2>0)then ! make permutation vector for side2
   allocate(keep2(dat%nval2))
   keep2=0
   do,i=1,dat%nval2
      match=regexecf(reg_inc2,dat%side2_d(i))
      if(excl2)match=.not. match ! negate in the case of exclusion
      if(.not. match)keep2(i)=1
   end do
end if




if(reg_inc2 >0)then
   if(sum(keep2) .eq. dat%nval2)then
      write(stderr,*)"ERROR: side 2 has zero length"
      stop
   end if
   allocate(perm2(dat%nval2))
   call sort_f90(perm2,keep2)
   st2=dat%nval2-sum(keep2)
else if(sort2f .ne. '')then
   st2=dat%nval2
   allocate(perm2(st2))
   forall(i=1:st2)perm2(i)=i
end if

!optionally sort according to predefined vector
if(sort2)then
   !read data from file
   open(unit=unit,file=sort2f,form='formatted')
   last=0
   nsort=0
   allocate(sort_d(chunk))
   sort_d=''
   memsz=chunk
   do while(last .eq. 0)

      if(nsort+1 > memsz)then
         call realloc_ptr(sort_d,chunk)
         memsz=memsz+chunk
      end if
      read(unit=unit,fmt='(A24)',iostat=last)sort_d(nsort+1)
      if(last .ne. 0)exit
      nsort=nsort+1
   end do
   close(unit)

   allocate(tmpvec_d(st2),tmpperm1(st2),tmpperm2(nsort))
   tmpvec_d=dat%side2_d(perm2(1:st2))
   
!calculate the permutation vector of the non-excluded part
   call get_permvec2(tmpvec_d,sort_d(1:nsort),.true.,1,24,tmpperm1,tmpperm2,ncom)

   !modify the already existing permutation vector
   perm2(1:st2)=perm2(tmpperm1(1:st2))
   deallocate(tmpvec_d,tmpperm1,tmpperm2,sort_d)
end if





! if(reg_inc1 + reg_inc2 +nrestr > 0)then ! permute systems
if(matpermute)then ! permute systems
   if(stage > 2)then
      dryrun=.false.
   else
      dryrun=.true.
   end if

   select case(dat%type)
   case('SYMV0___','SYMV1___','SYMV2___','SYMVN___','BDSYMV0_','BDSYMVN_','DIAVN___')!symmetric matrices (permute both sides)
      if(allocated(perm2))then
         write(stderr,*)"WARNING: permutation, requested for side 2 of symmetric matrix, is ignored (expand matrix first)"
      end if
      call BIN_permute(dat,perm1=perm1,both=.true.,dryrun=dryrun)
      dat%nval1=st1
      dat%nval2=dat%nval1
   case default
      if(st1 .ne. 0 .and. st2 .ne. 0)then
          call BIN_permute(dat,perm1=perm1,perm2=perm2,dryrun=dryrun)
          dat%nval1=st1
          dat%nval2=st2
       elseif(st1 .ne. 0)then
          call BIN_permute(dat,perm1=perm1,dryrun=dryrun)
          dat%nval1=st1
       elseif(st2 .ne. 0)then
          call BIN_permute(dat,perm2=perm2,dryrun=dryrun)
          dat%nval2=st2
       end if
       !change output type if output matrix is rectangular
       if(dat%type .eq. 'FULLSQVN' .and. dat%nval1 .ne. dat%nval2)dat%type='FULL2DVN'

   end select
   dat%pval1=packindex(dat,dat%nval1,dat%nval2)
   dat%pval2=1
      

end if

!Possibly matrix scaling
if(matscale .ne. 1.d0)then
   select case(dat%mtyp)   
   case ('P')
      dat%pack1=dat%pack1*matscale
   case('F','U','L')
      dat%mat1=dat%mat1*matscale
   end select
end if

!Possibly matrix scaling
if(vecscale .ne. 1.d0)then
   dat%vec=dat%vec*vecscale
end if

if(v2mat)then ! export vectors as new matrix
   if(associated(dat%pack1))deallocate(dat%pack1) ! throw away some goodies      
   if(associated(dat%mat1))deallocate(dat%mat1) ! throw away some goodies      
   ! modify data
   dat%mat1=> dat%vec ! point matrix to vector now
   dat%type='FULL2DVN' ! new matrix type
   dat%mtyp='F'
   dat%nval2=sz2
   dat%pval1=dat%nval1*dat%nval2
   dat%pval2=1
   dat%nvec=0 ! no vectors anymore
   dat%side2_d=>newside_d
end if

!!END MODIFYING DATA
!this part does not work with gfortran but does with ifort
!!open up standard output with long enough record length
! if(ascii)then
!    if(index(parse,'m') .ne. 0)then
!       if(transp)then
!           open(unit=stdout,recl=24+20*(dat%nval1+dat%nvec),form='formatted')
!       else
!          open(unit=stdout,recl=24+20*(dat%nval2+dat%nvec),form='formatted')
!       end if
!    end if
! end if


!PRINT INDEX INFORMATION TO SCREEN
if(header)then
   write(stdout,*)'MATRIX type: ',dat%type
   write(stdout,*)'Version: ',dat%ver
   !info on endianess
   if(dat%swap)then
      write(stdout,*)'Endianness of file has been swapped'
   end if

   select case(dat%type)
      case('SYMV0___','SYMV1___','SYMV2___','SYMVN___')
         dat%side2_d=>dat%side1_d
         write(stdout,*)'Symmetric matrix in packed storage,',dat%nvec,' associated vectors'
      case('BDSYMV0_','BDSYMVN_')
         write(stdout,*)'Block diagonal symmetric matrix in packed storage,',dat%nvec,' associated vectors'
         write(stdout,*)'Amount of diagonal blocks: ',dat%nblocks
         dat%side2_d=>dat%side1_d
      case('BDFULLV0','BDFULLVN')
         write(stdout,*)'Block diagonal full matrix in packed storage,',dat%nvec,' associated vectors'
         write(stdout,*)'Amount of diagonal blocks: ',dat%nblocks
      case('FULLSQV0','FULLSQVN')
         write(stdout,*)'Full square matrix,',dat%nvec,' associated vectors'
      case('FULL2DVN')
         write(stdout,*)'Full 2D matrix,',dat%nvec,' associated vectors'
      case('DIAVN___')!diagonal matrix
         write(stdout,*)'Diagonal matrix (packed),',dat%nvec,' associated vectors'
      case default
         write(stderr,*)'Unknown matrix type: ',dat%type
         stop
   end select
   write(stdout,*)'Full matrix dimensions:',dat%nval1, ' x ',dat%nval2
   write(stdout,*)'stored matrix dimensions:',dat%pval1, ' x ',dat%pval2

   write(stdout,*)'File description:'
   write(stdout,*)dat%descr

   write(stdout,*)'META data:',dat%nint,'integers',dat%ndbls,'doubles'

   if(dat%nint >0)then
      write(stdout,*)'INTEGER META data:'
      do,i=1,dat%nint
         write(stdout,*)i,dat%ints_d(i),dat%ints(i)
      end do
   end if

   if(dat%ndbls >0)then
      write(stdout,*)'DOUBLE META data:'
      do,i=1,dat%ndbls
         write(stdout,*)i,dat%dbls_d(i),dat%dbls(i)
      end do
   end if

   if(dat%nread > 0)then
      write(stdout,*)'README part:'
      do,i=1,dat%nread
         write(stdout,'(a80)')dat%readme(i)
      end do
   end if

!write a string indicating data is coming
   if(ascii)write(stdout,*)"FILE DATA:"
   
end if

if(ascii)then

   !construct automatic parse string if needed;
   if(root .or. diag)then
      parse='s'
      do,i=1,dat%nvec
         write(dum2,'(i6)')i
         
         parse=trim(parse)//",v"//adjustl(dum2)
      end do
      
      if(root)parse=trim(parse)//',r,'
      if(diag)parse=trim(parse)//',d,'
   end if
   


   !check if the transpose condition does not have a v in the parse string\
   if(transp .and. index(parse,'v') .ne. 0)then
      write(stderr,*)"ERROR: using the transpose condition does not allow vectors to be printed"
      stop
   end if
   

   if(transp)then
      do,i=1,dat%nval2
         !do,i=1,10
         call parse_n_print(dat,parse,i,transp)
      end do
   else
      do,i=1,dat%nval1
         !do,i=1,10
         call parse_n_print(dat,parse,i,transp)
      end do
   end if
else if(binary)then
   select case (outtype) !optionally change output type
   case('SYMVN___')
      if(dat%nval1 .ne. dat%nval2)then
         write(stderr,*)"ERROR: changing matrix type to symmetric is only allowed on square matrices"
         stop
      end if
      dat%type=outtype
!      write(stderr,*)dat%mtyp,dat%nval1,dat%nval2,dat%nread
      dat%mtyp='U'
      dat%pval1=dat%nval1*(dat%nval1+1)/2
      dat%pval2=1
   end select

   dat%file='stdout'
   call write_BINtype(dat)
end if


end program BIN_swiss


!subroutine to parse format string and print to standard output
subroutine parse_n_print(dat,parse,n,transp)
use binfiletools
implicit none
type(BINdat)::dat
character(*)::parse
character(1)::op
integer::n,m
!private arguments
integer::st,nd,i,stdout,nl,stderr
integer*8::ind
double precision::dum,tmp,dumcol
logical::transp

stdout=6
stderr=0

if(transp)then
   st=0
   op='n' ! set operator to read a new value
   
   dum=0.d0
   !write(unit=stdout,advance='NO',fmt='(i8)')packindex(dat,n,n)
   do while(st < len_trim(parse))
      !!find the occurence of the first operator
      nd=scan(parse(st+1:),'*/,+- ')+st
      !   nd=index(parse(st+1:),"\")+st ! end of column string "\"
      
      select case(parse(st+1:st+1))
      case('s')! parameter character string
         !write out character and a space
!         write(stdout,*)'test',n
         write(unit=stdout,advance='NO',fmt='(a24,a1)')dat%side2_d(n),' '
         !ignore the remaining part
         st=nd
!       case('v') !vector argument
!          read(parse(st+2:nd-1),*)m ! extract integer from parse string
         
!          call oper(dum,dat%vec(n,m),op)
         
!          st=nd-1
      case('d')! diagonal
         if(dat%nval1 .ne. dat%nval2)then
            write(stderr,*)'ERROR:diagonal output of a rectangular matrix is not allowed'
            stop
         end if
         ind=packindex(dat,n,n)
         call oper(dum,dat%pack1(ind),op)
         !      write(stdout,*)dat%pack1(ind)
         st=nd-1
      case('r')! root of the diagonal
         if(dat%nval1 .ne. dat%nval2)then
            write(stderr,*)'ERROR:diagonal output of a rectangular matrix is not allowed'
            stop
         end if
         ind=packindex(dat,n,n)
         call oper(dum,sqrt(dat%pack1(ind)),op)
         st=nd-1
      case('c')! column of the full matrix
         read(parse(st+2:nd-1),*)m ! extract integer from parse string
         ind=packindex(dat,m,n)
         dumcol=0.d0
         if(ind .ne. 0)dumcol=dat%pack1(ind)
         call oper(dum,dumcol,op)
         st=nd-1
      case('i')!index of the row or column
         write(unit=stdout,advance='NO',fmt='(i6,1x)')n
         !ignore the remaining part
         st=nd
      case('m')! full matrix, only prints stored values (so upper triangle or blocks)
         nl=0
         do,i=1,dat%nval1
            ind=packindex(dat,i,n)
            if( ind .ne. 0)then
               nl=nl+1
!                if(mod(nl,6) .eq. 0)then !print with line end
!                   write(unit=stdout,fmt='(G18.10)')dat%pack1(ind)
!                else
                  write(unit=stdout,advance='NO',fmt='(G22.14)')dat%pack1(ind)
!                end if
            end if
         end do
         st=nd !no operations allowed on the 
      case('*','/','+','-')
         op=parse(st+1:st+1)
         st=nd
      case(',') ! empty buffer and add a space
         write(unit=stdout,advance='NO',fmt='(G22.14,a1)')dum,' '
         op='n' !reset operator to rerad in a new value
         st=nd
      case default ! assume it is a scalar number
         read(parse(st+1:nd-1),*)tmp
         call oper(dum,tmp,op)
         st=nd-1
      end select
   end do
   
else

   st=0
   op='n' ! set operator to read a new value
   
   dum=0.d0
   !write(unit=stdout,advance='NO',fmt='(i8)')packindex(dat,n,n)
   do while(st < len_trim(parse))
      !!find the occurence of the first operator
      nd=scan(parse(st+1:),'*/,+- ')+st
      !   nd=index(parse(st+1:),"\")+st ! end of column string "\"
      
      select case(parse(st+1:st+1))
      case('s')! parameter character string
         !write out character and a space

         write(unit=stdout,advance='NO',fmt='(a24,a1)')dat%side1_d(n),' '
         !ignore the remaining part
         st=nd
      case('v') !vector argument
         read(parse(st+2:nd-1),*)m ! extract integer from parse string
         
         call oper(dum,dat%vec(n,m),op)
         
         st=nd-1
      case('d')! diagonal
         if(dat%nval1 .ne. dat%nval2)then
            write(stderr,*)'ERROR:diagonal output of a rectangular matrix is not allowed'
            stop
         end if
         ind=packindex(dat,n,n)
         call oper(dum,dat%pack1(ind),op)
         !      write(stdout,*)dat%pack1(ind)
         st=nd-1
      case('r')! root of the diagonal
         if(dat%nval1 .ne. dat%nval2)then
            write(stderr,*)'ERROR:diagonal output of a rectangular matrix is not allowed'
            stop
         end if
         ind=packindex(dat,n,n)
         call oper(dum,sqrt(dat%pack1(ind)),op)
         st=nd-1
      case('c')! column of the full matrix
         read(parse(st+2:nd-1),*)m ! extract integer from parse string
         ind=packindex(dat,n,m)
         dumcol=0.d0
         if(ind .ne. 0)dumcol=dat%pack1(ind)
         call oper(dum,dumcol,op)
         st=nd-1
      case('i')!index of the row or column
         write(unit=stdout,advance='NO',fmt='(i6,1x)')n
         !ignore the remaining part
         st=nd
      case('m')! full matrix, only prints stored values (so upper triangle or blocks)
         nl=0
         
         do,i=1,dat%nval2
            ind=packindex(dat,n,i)
!            write(0,*)"ind",ind,n,i
            if( ind .ne. 0)then
               nl=nl+1
!                if(mod(nl,6) .eq. 0)then
!                   write(unit=stdout,fmt='(G18.10)')dat%pack1(ind)
!                else
                  write(unit=stdout,advance='NO',fmt='(G22.14)')dat%pack1(ind)
!               end if
            end if
         end do
         st=nd !no operations allowed on the 
      case('*','/','+','-')
         op=parse(st+1:st+1)
         st=nd
      case(',') ! empty buffer and add a space
         write(unit=stdout,advance='NO',fmt='(G22.14,a1)')dum,' '
         op='n' !reset operator to rerad in a new value
         st=nd
      case default ! assume it is a scalar number
         read(parse(st+1:nd-1),*)tmp
         call oper(dum,tmp,op)
         st=nd-1
      end select
       
      
      
   end do
end if

!write a newline
write(unit=stdout,advance='NO',fmt='(/)')


end subroutine parse_n_print

subroutine oper(dat1,dat2,op)
implicit none
character(1)::op
double precision::dat1,dat2

select case(op)
case('*') ! multiply
!   write(0,*)dat1,dat2
   dat1=dat1*dat2
case('/')!divide
   dat1=dat1/dat2
case('+')!add
   dat1=dat1+dat2
case('-')! subtract
   dat1=dat1-dat2
case('n')! take the second value

   dat1=dat2 

end select

end subroutine oper
      

!help function prints help info to screen
subroutine help()
implicit none
character(4)::frmt
integer::stderr

stderr=0
frmt='(A)'


write(stderr,frmt)'Program BIN_swiss prints out file information and/or modifies input data/descriptions..'
write(stderr,frmt)'from a binary data file containing matrices'
write(stderr,frmt)'Usage: BIN_swiss [OPTIONS] [BINFILE]'
write(stderr,frmt)'Where OPTIONS can be:'
write(stderr,frmt)'  -p=PARSE: Prints out data in user defined columns. PARSE is a string which contains commas'
write(stderr,frmt)'            to separate the columns and letter codes to determine what should be printed:'
write(stderr,frmt)'            for example -p=s,v1,v2,r constructs output which column are the 24character'
write(stderr,frmt)'            parameter description, vector 1, vector 2 and the root of the diagonal of the matrix'
write(stderr,frmt)'            the following parameters are recognized: s (parameter description), vn(vector n),'
write(stderr,frmt)'            d(diagonal of the matrix), r(square root of the diagonal of matrix),'
write(stderr,frmt)'            cn(nth column of the matrix),m( matrix ) and i (index) the index number. '
write(stderr,frmt)'            The last code will produce only output for the defined part of the matrix'
write(stderr,frmt)'            Simple arithmetic operators (+-/*) are allowed on vectors and the column/diagonal'
write(stderr,frmt)'            of the matrix.'
write(stderr,frmt)'            For example: -p=v1+v2,100*r Prints out the element by element'
write(stderr,frmt)'            addition of vectors 1 and 2 and scales the root of the diagonal by 100'
write(stderr,frmt)'            The order of operator execution is from left to right so v1+c2/3 is executed as (v1+c2)/3'
write(stderr,frmt)'            For scalars provided on the command line no scientific notation is allowed (interferes with -)'
write(stderr,frmt)'  -d: print all vectors to screen plus diagonal of the matrix (parse string s,v1,v2...vn,d)'
write(stderr,frmt)'  -s: same as -d except that the diagonal entries are the square root (parse string s,v1,v2...vn,r)'
write(stderr,frmt)'  -t: Transpose system: prints values of the other side (for example s refers to description'
write(stderr,frmt)'      of the second side)'
write(stderr,frmt)'  -n: do not print header to screen'
write(stderr,frmt)''
write(stderr,frmt)'  -Sm=MATRIXSCALE Apply scale to matrix'
write(stderr,frmt)'  -Sv=VECSCALE Apply scale to vectors'
write(stderr,frmt)''
write(stderr,frmt)'  -[I|E]1=REGEX: only include(I) or exclude (E) rowdata which obeys the'
write(stderr,frmt)'                 regular expression REGEX(ascii output of'
write(stderr,frmt)'                 the matrix is not allowed)'
write(stderr,frmt)'  -[I|E]1:REGEXFILE: Same as above but reads regular expressions from a file REGEXFILE'
write(stderr,frmt)'  -[I|E]2[=|:]REGEX: only include(I) or exclude (E) columndata which obeys the'
write(stderr,frmt)'                 regular expression REGEX(ascii output'
write(stderr,frmt)'                 of the matrix is not allowed)'
write(stderr,frmt)'  -[I|E]2:REGEXFILE: Same as above but reads regular expressions from a file REGEXFILE'
write(stderr,frmt)'  -P[12]:PERMFILE Permute the side (1 and/or 2) by a supplying in PERMFILE a desired parameter'
write(stderr,frmt)'   order (one 24 char parameter  per line)'
write(stderr,frmt)''
write(stderr,frmt)'  -R=MIN,PARA,MAX: restrict the output of PARA based on numerical range MIN, MAX'
write(stderr,frmt)"   PARA is a string similar to the parse string mentioned above. 'v2' means vector "
write(stderr,frmt)"   column 2, 'd' means diagonal of the matrix (needs to be square) Multiple -R options are allowed"
write(stderr,frmt)''
write(stderr,frmt)'  -MR=READMESTRING: Append string to readme data'
write(stderr,frmt)'  -MD=DNAME/DVALUE: Update/Append meta double value'
write(stderr,frmt)'  -MI=INAME/IVALUE: Update/Append meta integer value'
write(stderr,frmt)'  -Md=DESCR: Update description string'
write(stderr,frmt)'  -MI-: delete integer meta data'
write(stderr,frmt)'  -MD-: delete double meta data'
write(stderr,frmt)'  -MR-: delete readme meta data'
write(stderr,frmt)'  -MP1=STRBEFORE/STRAFTER: Adapt parameter description of side 1. Replace the'
write(stderr,frmt)'       part matching STRBEFORE with STRAFTER. The result will be truncated to'
write(stderr,frmt)'       24 characters'
write(stderr,frmt)'  -MP2=STRBEFORE/STRAFTER: Same as -MP1 but then for side description 2'
write(stderr,frmt)'  -MT1=STAMP: Add stamp to the end of row parameter description'
write(stderr,frmt)'  -MT2=STAMP: Add stamp to the end of column parameter description'
write(stderr,frmt)'  -ML1 : Left adjust parameter description of side 1'
write(stderr,frmt)'  -ML2 : Left adjust parameter description of side 2'
write(stderr,frmt)'  -MF=INFOFILE: multiple meta data strings are provided in the INFOFILE'
write(stderr,frmt)'      For example a line with D=DNAME/DVALUE ' 
write(stderr,frmt)'  -b: output in binary form (readable again by BIN_swiss), automatically implies '
write(stderr,frmt)'      -n (and no -d, -p -s allowed)'
write(stderr,frmt)'  -e: expand sparse matrix (WARNING: some sparse systems might not fit in your memory when'
write(stderr,frmt)'      expanded)'
write(stderr,frmt)'  -cs: Make the matrix symmetric (ignores the lower triangle of the matrix).'
write(stderr,frmt)'       This option has only effect in combination with the -b option'
write(stderr,frmt)'  -V=TAG1/TAG2/...: Export the vector part of the matrix as a new matrix. New names are denoted by tags'
write(stderr,frmt)'  -V:TAGFILE:  Same as above but side description is read from file TAGFILE'
write(stderr,frmt)'       This option has only effect in combination with the -b option'

!write(stderr,frmt)

stop
end subroutine help


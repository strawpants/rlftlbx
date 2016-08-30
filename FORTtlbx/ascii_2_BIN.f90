!program to create a diagonal binary file from an ascii input
!!Updated by Roelof Rietbroek, Wed Mar 10 14:13:57 2010
!! made metadata adding consistent with BIN_swiss

!!Updated by Roelof Rietbroek, Fri Apr 23 10:20:54 2010
!! also allow vectors to be put in a matrix
!!Updated by Roelof Rietbroek, Wed Feb 22 11:21:56 2012
!! fixed bug (use a different unit for the TAGSfile)




program ascii_2_BIN
use binfiletools
use bin_operations
use forttlbx
implicit none
integer::i,j,itharg,narg,iargc
type(BINdat)::out
integer::nmxcol,ncol,chunk,mem,stderr
parameter(nmxcol=2000,chunk=100000,stderr=0)
double precision::dummies(nmxcol)
double precision,pointer::work(:) ! to read data into
integer::unit,last,ndat
character(200)::infile
character(30000)::dum
logical::mat2d
integer::sz2,ind,indold,unit2

!defaults/initializations
out%file='stdout'
out%descr='Constructed by ascii_2_BIN'

mat2d=.false.
sz2=0

unit=5 ! default is standard input
unit2=14
infile=''

!!process command line options
narg=iargc()
itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:2))
      case('M')!insert meta data
         call BIN_putmeta(out,trim(dum(3:))) ! append meta data
      case('T')! put columns in a 2d matrix
         mat2d=.true.
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
               call realloc_ptr(out%side2_d,1) ! add element to character side2 description
               out%side2_d(sz2)=dum(indold+1:ind-1)
               indold=ind!reset starting index
            end do
         case(':')!read tags from file
            open(unit=unit2,file=trim(dum(4:)),form='formatted',status='old')
            last=0
            do while(last .eq. 0)
               read(unit=unit2,iostat=last,fmt='(A24)')dum
               if(last .ne. 0)exit
               sz2=sz2+1
               call realloc_ptr(out%side2_d,1)
               out%side2_d(sz2)=dum(1:24)
              
            end do
            close(unit2)
         case default
            write(stderr,*)'ERROR: processing -T option:',dum(1:3)
            stop     
         end select
      case('h')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:2)
         stop
      end select
   else!argument is not an option but the name of the asciifile
      unit=13
      infile=trim(dum)
   end if
end do

!open file
if(unit .ne. 5)open(unit=unit,file=trim(infile),form='formatted',status='old')

!try getting the amount of valid columns by scanning the first line
read(unit=unit,fmt='(A30000)',iostat=last)dum
ncol=0
do,i=1,nmxcol

   read(dum(25:),*,iostat=last)dummies(1:i)
!   write(0,*)last,i,dum(25:100)
   if(last .ne. 0)then
      ncol=i-1
      exit
   end if
end do

if(ncol <1)then
   write(stderr,*)"ERROR: no data found"
   stop
end if
if(mat2d .and. ncol < sz2)then
   write(stderr,*)"ERROR: amount of columns found are less then needed:",ncol,sz2
   stop
end if

!allocate work array
allocate(work(ncol*chunk),out%side1_d(chunk))

mem=chunk

!put first data line in the output structure and work array
ndat=1

out%side1_d(ndat)=dum(1:24)
work(1:ncol)=dummies(1:ncol)
!read all data in work array

last=0
do while (last .eq. 0)! loop until end of file
   read(unit=unit,fmt='(A30000)',iostat=last)dum
   if(last .ne. 0)exit
   ndat=ndat+1
   
   if(ndat > mem)then ! reallocate more memory if needed
      call realloc_ptr(work,chunk*ncol)
      mem=mem+chunk
      call realloc_ptr(out%side1_d,chunk)
   end if


   ind=(ndat-1)*ncol


   read(dum(25:),*)work(ind+1:ind+ncol)
   out%side1_d(ndat)=dum(1:24)
end do

if(unit .ne. 5)close(unit)



!set up output structures
if(mat2d)then
   out%nval1=ndat
   out%nval2=sz2
   out%pval1=ndat*sz2
   out%pval2=1
   out%type='FULL2DVN' !full matrix
   out%mtyp='P'
   out%nvec=ncol-sz2
   allocate(out%pack1(out%pval1))
   do,i=1,sz2
      ind=(i-1)*ndat

      out%pack1(ind+1:ind+ndat)=work(i+out%nvec:ncol*ndat:ncol)
   end do
   if(out%nvec>0)then
      allocate(out%vec(ndat,out%nvec))
      do,i=1,out%nvec
         out%vec(:,i)=work(i:ncol*ndat:ncol) ! copy data with stride ncol
      end do

   end if
   

else
   out%nval1=ndat
   out%nval2=ndat
   out%pval1=ndat
   out%pval2=1
   out%type='DIAVN___'
   out%mtyp='P'
   out%nvec=ncol-1

   allocate(out%pack1(out%pval1))
   out%pack1(1:ndat)=work(1+out%nvec:ncol*ndat:ncol)


   if(out%nvec>0)then
      allocate(out%vec(ndat,out%nvec))
      do,i=1,out%nvec
         out%vec(:,i)=work(i:ncol*ndat:ncol) ! copy data with stride ncol
      end do

   end if
end if

!write to output
call write_BINtype(out)



end program ascii_2_BIN


subroutine help()
integer::stderr
character(8)::frmt
stderr=0
frmt='(A)'
write(stderr,frmt)"Program ascii_2_BIN converts ascii input to a diagonal binary matrix"
write(stderr,frmt)"Usage ascii_2_BIN [OPTIONS] [ASCIIFILE]"
write(stderr,frmt)"Where ASCIIFILE contains a 24 character parameter code and several data columns"
write(stderr,frmt)"When no file is provided the program will read from standard input"
write(stderr,frmt)"The default is that the last data column is assumed to be the diagonal of the matrix"
write(stderr,frmt)" while the columns before are put in the vector"
write(stderr,frmt)""
write(stderr,frmt)"OPTIONS may be:"
write(stderr,frmt)"  -T=TAG1/TAG2/TAG3/..: Put the data columns in a 2d matrix instead of in the vector"
write(stderr,frmt)'  TAG1/TAG2/.. etc denote the new names of the columns and must be consistent'
write(stderr,frmt)"  with the amount of provided columns."
write(stderr,frmt)"  -T:TAGSFILE: Same as above but names are read from a file TAGSFILE (1 parameter per line)"
write(stderr,frmt)'  -MR=READMESTRING: Append string to readme data'
write(stderr,frmt)'  -MD=DNAME/DVALUE: Update/Append meta double value'
write(stderr,frmt)'  -MI=INAME/IVALUE: Update/Append meta integer value'
write(stderr,frmt)'  -Md=DESCR: Update description string'
write(stderr,frmt)'  -MI-: delete integer meta data'
write(stderr,frmt)'  -MD-: delete double meta data'
write(stderr,frmt)'  -MR-: delete readme meta data'
write(stderr,frmt)'  -MP1=STRBEFORE/STRAFTER: Adapt parameter description of side 1. Replace the'
write(stderr,frmt)'       part matching STRBEFORE with STRAFTER. Lengths of STRBEFORE and STRAFTER '
write(stderr,frmt)'       must be equal'
write(stderr,frmt)'  -MP2=STRBEFORE/STRAFTER: Same as -MP1 but then for side description 2'
write(stderr,frmt)'  -MF=INFOFILE: multiple meta data strings are provided in the INFOFILE'
write(stderr,frmt)'      For example a line with D=DNAME/DVALUE ' 
write(stderr,frmt)""
write(stderr,frmt)"The matrix system will be printed to standard output"
stop
end subroutine help

!!Coded by Roelof Rietbroek, Wed Jun  6 21:50:01 2012
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!program which converts (some) groops binary data in BIN format
!!fortran module which contains routines to read in groops binary data

module GroopsFile
use forttlbx
use binfiletools
use shtlbx,only:SHgDWpos
implicit none
integer,parameter::stderr=0

type GroopsMAT
character(200)::file=''
character(1)::ftype
character(12)::type=''
character(12)::version=''
integer::mtype=0
integer::nrow=0
integer::ncol=0
integer*8::npack=0
integer::uplo=0
integer::cunit=0
integer::kinddbl=8
double precision,pointer::pack1(:)=>null()
double precision,pointer::mat1(:,:)=>null()
end type GroopsMAT

type GroopsXML
integer::unit=15
character(200),pointer::buffer(:)
character(200)::file
end type GroopsXML

contains
!!reads a groops normal equation system froma binary file
subroutine read_GROOPSNEQ(dat)
type(BINdat), intent(inout)::dat
type(GroopsMAT)::Gnorm,Grhs
type(GroopsXML)::Ginfo
integer:: i
character(200)::dum,fbase
integer*8::maxchunk,sz,st
integer::nchunk
integer::cunit

!defaults

maxchunk=1073741824 !chunked read in bytes 

!extract the base name of the files
i=index(dat%file,'.',.true.)
if(i .ne. 0)then

    fbase=dat%file(1:i-1)
    !possibly strip another suffix (.info or .righthandside)
    i=index(fbase,'.',.true.)
    if(i .ne. 0)then
        fbase=fbase(1:i-1)
    end if

end if


!read normal equation matrix
Gnorm%file=trim(fbase)//'.dat'//char(0)
call read_GROOPSBIN(Gnorm)

!read the righthandside
Grhs%file=trim(fbase)//'.rightHandSide.dat'//char(0)
call read_GROOPSBIN(Grhs)
if(Grhs%ncol .ne. 1)then
    write(stderr,*)"ERROR: cannot not handle systems with multiple right hand sides"
    stop  
end if
Ginfo%file=trim(fbase)//'.info.xml'//char(0)
call read_Groopsxml(Ginfo)

!write(stderr,*)Grhs%nrow,Gnorm%nrow
!!set up and gather stuff for the  normal matrix
dat%pack1=>Gnorm%pack1

dat%nvec=2
allocate(dat%vec(Grhs%nrow,dat%nvec))

dat%vec=0
dat%vec(:,1)=Grhs%pack1(:)
dat%descr="Converted Groopsnormal system"

dat%nint=2
allocate(dat%ints_d(dat%nint),dat%ints(dat%nint))
dat%ints=0

dat%ndbls=5
allocate(dat%dbls_d(dat%ndbls),dat%dbls(dat%ndbls))
dat%dbls=0.0

dat%mtyp='P'
!select case(Gnorm%uplo)
!case(0)
    !dat%mtyp='L'
!case(1)
    !dat%mtyp='U'
!end select


dat%ints_d(1)='Nobs'
dat%ints(1)=getXMLint(Ginfo,'observationCount')
dat%ints_d(2)='Nunknowns'
dat%ints(2)=getXMLint(Ginfo,'parameterCount')


dat%dbls_d(1)='CTime'
dat%dbls_d(2)='STime'
dat%dbls_d(3)='ETime'
dat%dbls_d(4)='LtPL'
dat%dbls(4)=getXMLdbl(Ginfo,'lPl')
dat%dbls_d(5)='Sigma0'
dat%dbls(5)=1.0

dat%type='SYMVN___'
dat%nval1=Gnorm%nrow
dat%nval2=dat%nval1
dat%pval1=Gnorm%npack
dat%pval2=1



end subroutine read_GROOPSNEQ




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function readUint(cunit)
integer::readUint
integer,intent(in)::cunit
character(8)::Uint

!read an unisgned integer
call cread(cunit,8,Uint)
!interpret as normal integere
readUint=transfer(Uint,readUint)

end function readUint



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readString(cunit,string)
integer,intent(in)::cunit
character(*),intent(out)::string
integer*8::ssize

ssize=readUint(cunit)

if (ssize > len(string))then
write(0,*)"ERROR:GroopsFile.f90 cannot read string (too few chareacters reserved)"
stop 
end if
string=''
call cread(cunit,ssize,string)
!append a trailing null char
string(ssize+1:ssize+1)=char(0)
string=trim(string)!just to be sure it's fortran compatible
end subroutine readString 








!! read matrix from a file
subroutine read_GROOPSBIN(Gfile)
type(GroopsMAT),intent(inout)::Gfile
integer::opencstream,i
integer*8::maxchunk,st
integer*4::sz
parameter(maxchunk=1073741824)
integer::nchunk

character(200)::dummy
!open matrix file as stream ( read only)
Gfile%cunit=opencstream(0,Gfile%file)


!read the first bytes ( should start with a 'B')
call cread(Gfile%cunit,1,Gfile%ftype)
if(Gfile%ftype .ne. 'B')then
   write(stderr,*)"ERROR: not a groops binary file",trim(Gfile%file)
   stop  
end if

!read i version and type strings 

call readString(Gfile%cunit,Gfile%type)
call readString(Gfile%cunit,Gfile%version)

if(Gfile%version(1:3) .ne. '1.1')then
    write(stderr,*)"ERROR no support for version 1.1"
    stop  
end if

!read remainder up to 64 bits
call cread(Gfile%cunit,modulo(1+len_trim(Gfile%version)+len_trim(Gfile%type)-2,8),dummy)


!read the matrix type
Gfile%mtype=readUint(Gfile%cunit)


select case(Gfile%mtype)
case(0)!GENERAL
!read dimensions
    Gfile%nrow=readUint(Gfile%cunit)
    Gfile%ncol=readUint(Gfile%cunit)
    Gfile%npack=Gfile%nrow*Gfile%ncol

      !if(.not. associated(gfile%pack1))allocate(gfile%pack1(gfile%npack))
      !call creadmany(gfile%cunit,Gfile%npack*kinddbl,gfile%pack1(1))
case(1,2) !SYMMETRIC or TRIANGULAR

    !read dimensions and storage type
    Gfile%uplo=readUint(Gfile%cunit)
    Gfile%nrow=readUint(Gfile%cunit)
    Gfile%ncol=Gfile%nrow
    Gfile%npack=Gfile%nrow*(Gfile%nrow+1)/2
      !if(.not. associated(gfile%pack1))allocate(gfile%pack1(gfile%npack))
      !call creadmany(gfile%cunit,Gfile%npack*kinddbl,gfile%pack1(1))
     !nchunk=(kinddbl*gfile%npack)/maxchunk
      !do,i=1,nchunk
         !st=1+((i-1)*maxchunk)/kinddbl
         !sz=maxchunk
         !call cread(gfile%cunit,sz,gfile%pack1(st))
      !end do
      !!read remainder

      !st=1+(nchunk*maxchunk)/kinddbl
      !sz=kinddbl*(gfile%npack-st+1)
!!     write(0,*)'last',st,sz,nchunk
      !if(sz>0)call cread(gfile%cunit,sz,gfile%pack1(st))
end select

if(.not. associated(gfile%pack1))allocate(gfile%pack1(gfile%npack))
call cread(gfile%cunit,Gfile%npack*kinddbl,gfile%pack1(1))
!write(0,*)gfile%version,gfile%type,gfile%mtype,gfile%nrow,gfile%ncol



call cclose(gfile%cunit)

end subroutine read_groopsbin

!!reads a groops xml file in a buffer  
subroutine read_groopsXML(GF)
    type(GroopsXML)::GF
    integer:: stat=0
    integer:: nline=0
    
    allocate(GF%buffer(4))

    open(unit=GF%unit,file=GF%file,form='formatted',status='old')
    
    do while(stat .eq. 0)
        nline=nline+1
        if(size(GF%buffer)< nline)then
            call reallocate_cptr(GF%buffer,100)
        end if
        read(GF%unit,iostat=stat,fmt='(A200)')GF%buffer(nline)
        if(stat .ne. 0)then
            nline=nline-1
            exit
        end if
        nline=nline+1
    end do
    close(GF%unit)
end subroutine read_groopsXML

!function which seaches in a xml buffer for a certain integer parameter
function getXMLint(GF,search)
    integer*8::getXMLint
    type(GroopsXML)::GF
    character(*),intent(in)::search
    integer::i,istart,iend
    

    do,i=1,size(GF%buffer)
            istart=index(GF%buffer(i),'<'//trim(search)//'>')
            if(istart .eq. 0)then
                cycle
            end if
            iend=index(GF%buffer(i),'</'//trim(search)//'>')

            read(GF%buffer(i)(istart+len_trim(search)+2:iend-1),*)getXMLint
    end do


end function getXMLint

!function which seaches in a xml buffer for a certain double parameter
function getXMLdbl(GF,search)
    double precision::getXMLdbl
    type(GroopsXML)::GF
    character(*),intent(in)::search
    integer::i,istart,iend
    

    do,i=1,size(GF%buffer)
            istart=index(GF%buffer(i),'<'//trim(search)//'>')
            if(istart .eq. 0)then
                cycle
            end if
            iend=index(GF%buffer(i),'</'//trim(search)//'>')

            read(GF%buffer(i)(istart+len_trim(search)+2:iend-1),*)getXMLdbl
    end do


end function getXMLdbl

subroutine getDegreeWise_desc(descr,lmax,lmin)
    integer,intent(in)::lmax,lmin
    character(24),intent(inout),pointer::descr(:)

    integer::nsh,indx,l,m
    nsh=SHgDWpos(lmax,lmax,1,lmax,lmin)

    if(.not. associated(descr))then
        allocate(descr(nsh))
    end if

    indx=0
    do,l=lmin,lmax
        do,m=0,l
               indx=indx+1
               write(descr(indx),'(a3,1x,i3,i3)')'GCN',l,m
               if(m .ne. 0)then
                   indx=indx+1
                    write(descr(indx),'(a3,1x,i3,i3)')'GSN',l,m
               end if 
        end do
    end do

end subroutine getDegreeWise_desc



end module GroopsFile




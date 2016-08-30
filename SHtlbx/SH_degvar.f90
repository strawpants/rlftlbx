!!program to calculate (error) degree variances of a spherical harmonic set
!!prints to standard output
!! may read from standard input

!!Coded by Roelof Rietbroek, Wed Apr 30 14:19:04 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Wed Dec  7 14:44:30 2011
!! increased the length of the allowable filename to 200
!!Updated by Roelof Rietbroek, Thu Jan  5 11:19:59 2012
!! also allow logarithmic output


program SH_degvar
use SHtlbx
use FORTtlbx
implicit none
character(200)::dum,shfile
integer::l,m,itharg,narg,st,nd,pos,lmax,lmin,ind,ftyp,i,stderr,stdout
integer::lmaxdum
logical::stdin,error,generic,limdeg,averaged,l10
double precision,allocatable,dimension(:)::clm,slm,clm_sig,slm_sig
double precision::sigvar,sigerrvar
integer::iargc
!defaults
limdeg=.false.
stdin=.false.
error=.false.
generic=.false.
averaged=.false.
lmax=0
lmin=0
lmaxdum=0
itharg=0
stderr=0
stdout=6
shfile=''
l10=.false.



!get command line arguments

narg=iargc()

do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1).eq. '-')then
      select case(dum(2:2))
      case('l')
         limdeg=.true.
         ind=index(dum,',')
         if(ind .ne. 0) then !also a min degree is specified
            read(dum(3:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         else
            read(dum(3:),*)lmax
         end if
      case('e')
         error=.true.
      case('G')
         generic=.true.
         ftyp=0
      case('a')
         averaged=.true.
      case('L')
         l10=.true.
      case('h')
         call help()
      case default
         write(stderr,*)"ERROR SH_degvar: Unrecognized option",dum(1:2)
         call help()

      end select
   else !filename
      shfile=dum
   end if

end do

if(shfile .eq. '')then
   stdin=.true.
end if


if(generic .and. .not. limdeg)then
   write(stderr,*)"ERROR SH_degvar: The -l option must be used with generic files"
   stop
end if


!get meta data of the file (or standard output)
!we do this only when limdeg is not specified (this enables generic files to be also read from standard input)

if(.not. generic)then

   if(stdin)then
      call SH_readmeta(type=ftyp,lmax=lmaxdum)
   else
      call SH_readmeta(filen=shfile,type=ftyp,lmax=lmaxdum)
   end if
end if


if(.not. limdeg)then
   if(lmaxdum < lmax)then
      write(stderr,*)"ERROR SH_degvar: maximum degree supported by file is less than requested"
      stop
   end if

   lmax=lmaxdum

end if





pos=SH_pos(lmax,lmax)

if(error)then
   allocate(clm_sig(pos),slm_sig(pos),clm(pos),slm(pos))
   clm_sig=0.d0
   slm_sig=0.d0
   clm=0.d0
   slm=0.d0
else
   allocate(clm(pos),slm(pos))
   clm=0.d0
   slm=0.d0
end if

!read in data


if(error)then
   if(stdin)then
      call SH_readgrav(clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=ftyp)
   else
      call SH_readgrav(filen=shfile,clm=clm,slm=slm,clm_sig=clm_sig,slm_sig=slm_sig,type=ftyp)
   end if

else
   if(stdin)then
      call SH_readgrav(clm=clm,slm=slm,type=ftyp)
   else
      call SH_readgrav(filen=shfile,clm=clm,slm=slm,type=ftyp)
   end if
end if



!loop over degrees

do,l=lmin,lmax
   
   st=SH_pos(l,0)
   nd=SH_pos(l,l)


      sigvar=dot_productblas(clm(st:nd),clm(st:nd))+dot_productblas(slm(st:nd),slm(st:nd))

   if(averaged)sigvar=sigvar/(2.d0*l+1.d0)
   if(l10)sigvar=log10(sigvar)
   if(error)then
      sigerrvar=dot_productblas(clm_sig(st:nd),clm_sig(st:nd))+dot_productblas(slm_sig(st:nd),slm_sig(st:nd))
      
      if(averaged)sigerrvar=sigerrvar/(2.d0*l+1.d0)

      if(l10)sigerrvar=log10(sigerrvar)
      write(stdout,*)l,sigvar,sigerrvar
   else
      write(stdout,*)l,sigvar
   end if
end do
   

 
end program SH_degvar


subroutine help()
implicit none
character(8)::frmt
integer::unit
unit=6
frmt='(A)'

write(unit,frmt)"Program SH_degvar prints (error) degree variances to screen"
write(unit,frmt)"Usage SH_degvar [OPTIONS] [SHFILE]"
write(unit,frmt)"Where SHFILE is a file with spherical harmonic coefficients"
write(unit,frmt)"And [OPTIONS] may be:"
write(unit,frmt)""
write(unit,frmt)"-lLMAX[,LMIN]: Limit output to maximum degree LMAX (and possibly minimum degree LMIN)"
write(unit,frmt)"-e: also calculate error degree variances of the file (output will have three columns)"
write(unit,frmt)"-G: Read generic file from standard input "
write(unit,frmt)"-a: Averaged signal degree variance, Divide the signal degree variance by 2l+1, to "
write(unit,frmt)"    obtain an averaged variance per coefficient"
write(unit,frmt)"-L: Output the log10 values of the signal (and error)"
write(unit,frmt)"-h: prints this help message"

stop
end subroutine help




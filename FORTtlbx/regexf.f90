! this file contains two C wrappers in fortran 
! this file cooperates with REGEX.c
 

!!Coded by Roelof Rietbroek, Mon Mar  2 14:41:30 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Tue Apr 13 14:14:12 2010
!! added the possibility to read a file name  and compile from its contents

!!Updated by Roelof Rietbroek, Wed Aug  4 16:59:28 2010
!! also allow trailing spaces in the search string when calling regexecf 


!function to compile a regular expression
function regcompf(regex,file)
implicit none
  integer regcompf,compregex
  character(*)::regex
  logical,intent(in),optional::file
  
  logical::ifile
  integer::maxlen
  parameter(maxlen=500000)
  character(maxlen)::dumchar !maximum regex length
  character(80)::dum
  integer::unit,stat,stderr
  stderr=0
  ifile=.false.

  if(present(file))ifile=file

  if(ifile)then ! input is a file name  with one regex per line
       unit=13
       stat=0
     !load file
     open(unit=unit,file=trim(regex),form='formatted')
     read(unit=unit,iostat=stat,fmt='(A80)')dum

     dumchar='('//trim(dum)//')'
     do while(stat .eq. 0) ! append additional regular expressions
        read(unit=unit,iostat=stat,fmt='(A80)')dum
        if(stat <0)exit ! end of file found
        if(len_trim(dumchar)+3+len_trim(dum) >maxlen)then
           write(stderr,*)"ERROR in regcompf: maximum length of regular expression reached"
           stop
        end if
        dumchar=trim(dumchar)//'|('//trim(dum)//')'

     end do

     close(unit)
     !compile regular expression
     regcompf=compregex(trim(dumchar)//char(0))
  else
     !put in a null character at the end of the string and compile regular expression
     regcompf=compregex(trim(regex)//char(0))
     !  write(0,*)'testing',regcompf  
  end if
end function regcompf


!function to match a string  with a compiled (index provided with nregex) regular expression

function regexecf(nregex,string)
implicit none
  character(*)::string
  integer::nregex,execregex
  logical::regexecf

  if(execregex(nregex,string//char(0)) .eq. 0)then
     regexecf=.true. ! set to true in case of a match
  else
     regexecf=.false. ! no match
  end if
!  write(0,*)regexecf
end function regexecf




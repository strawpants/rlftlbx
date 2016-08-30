!!Coded by Roelof Rietbroek, Wed Jul 21 16:40:47 2010
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de



!fortran 90 wrappers for C sorting routines
!!Updated by Roelof Rietbroek, Wed Jun 15 09:59:27 2011
!! added double precision

!!TODO: alphabetically

subroutine isort_f90(perm,ints,inv)
implicit none
integer,intent(in)::ints(:)
integer,intent(out)::perm(:)
logical,optional,intent(in)::inv

integer::iinv
iinv=0
if(present(inv))then
   if(inv)iinv=1
end if
   !call C routine
call isort_c(perm(1),size(perm),iinv,ints(1))

end subroutine isort_f90

subroutine dsort_f90(perm,dbls,inv)
implicit none
double precision,intent(in)::dbls(:)
integer,intent(out)::perm(:)
logical,optional,intent(in)::inv

integer::iinv
iinv=0
if(present(inv))then
   if(inv)iinv=1
end if
   !call C routine
call dsort_c(perm(1),size(perm),iinv,dbls(1))

end subroutine dsort_f90

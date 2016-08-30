!!module file which descibes the interfaces for the grdtlbx library
!!Coded by Roelof Rietbroek, Wed Jun 20 16:40:18 2007
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

module grdtlbx
  implicit none

  !added a grid type holding a longitude and latitude grid
  !and a file name
  type geogrid
     double precision,pointer,dimension(:)::lon => null()
     double precision,pointer,dimension(:)::lat => null()
     double precision,pointer,dimension(:,:)::z => null()
     integer::nlon,nlat
     character(300)::file
  end type geogrid

  interface
     subroutine NC_write(file,lat,lon,z,history,zunits,zname,xunits,&
          xname,yunits,yname,ctime,pix,napp,nclse)
       implicit none
       character(*),intent(in)::file
       character(*),optional,intent(in)::history,zunits,zname,xunits,xname,yunits,yname
       double precision,dimension(:),intent(in)::lat,lon
       double precision,optional,intent(in)::ctime
       logical,intent(in),optional::pix ! use pixel registration
       double precision,dimension(:,:),intent(in)::z
       logical,intent(in),optional::nclse,napp !optional arguments allowing appending grids to an existing dataset
     end subroutine NC_write

     ! added 27 March 2013
     subroutine NC_read(file,lon,lat,z2d,ctime,etime,stime,zname,zslice)
       implicit none
       character(*),intent(in)::file
       double precision,pointer,dimension(:),intent(inout)::lon,lat
       double precision,pointer,dimension(:,:),intent(inout)::z2d
       double precision,optional,intent(out)::ctime,etime,stime
       character(*),intent(in),optional::zname
       integer,intent(in),optional::zslice
     end subroutine NC_read

!      subroutine NC_write(file,lat,lon,z,history,zunits,zname,xunits,xname,yunits,yname,ctime,pix)
!        implicit none
!        character(*),intent(in)::file
!        character(*),optional,intent(in)::history,zunits,zname,xunits,xname,yunits,yname
!        logical,intent(in),optional::pix ! use pixel registration
!        double precision,dimension(:),intent(in)::lat,lon
!        double precision,optional,intent(in)::ctime
!        double precision,dimension(:,:),intent(in)::z
!      end subroutine NC_write


     subroutine OBP_obseq_loadsh(A,tag,lon,lat,lmax,lmin)
       implicit none
       double precision, intent(out)::A(:,:)
       character(24),intent(out),optional::tag(:)
       double precision,intent(in)::lon(:),lat(:)
       integer,intent(in)::lmax,lmin
     end subroutine OBP_obseq_loadsh

     subroutine OBP_obseq_loadsh_equiref(A,tag,lon,lat,lmax,lmin)
       implicit none
       double precision, intent(out)::A(:,:)
       character(24),intent(out),optional::tag(:)
       double precision,intent(in)::lon(:),lat(:)
       integer,intent(in)::lmax,lmin
     end subroutine OBP_obseq_loadsh_equiref

     subroutine OBP_obseq_geocent(A,tag,lon,lat)
       implicit none
       double precision,intent(out)::A(:,:)
       character(24),intent(out),optional::tag(3)
       double precision,intent(in),dimension(:)::lon,lat
     end subroutine OBP_obseq_geocent


     subroutine OBP_obseq_globalmean(A,tag)
       implicit none
       double precision,intent(out)::A(:,:)
       character(24),intent(out),optional::tag(:)
     end subroutine OBP_obseq_globalmean

  end interface

end module grdtlbx

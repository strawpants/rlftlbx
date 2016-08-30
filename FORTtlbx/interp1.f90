!!Coded by Roelof Rietbroek, Wed May 22 19:39:17 2013
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

!fortran subroutine which interpolates data (1d)
!the input/output may have multiple columns/rows which will be interpolated in the same way

!x and xi should be monotonically increasing!!!

subroutine interp1(x,y,xi,yi,dim,method)
implicit none
double precision,intent(in)::x(:),xi(:),y(:,:)
integer,intent(in)::dim !dimension corresponding to the to be interpolated axis
character(*),intent(in)::method
double precision,intent(inout),pointer::yi(:,:) ! will be (re) allocated within the routine

!local variables
integer::nx,nxi,nd,dimd,stderr
integer::i,j,nw,ii
double precision,allocatable::wght(:)
double precision::dx
stderr=0

nx=size(x,1)

nxi=size(xi,1)
dimd=2
if(dim .eq. 2)dimd=1 !switch dimensions

nd=size(y,dimd)

if(.not. associated(yi))then

   if(dim .eq. 1)then
      allocate(yi(nxi,nd))
   else
      allocate(yi(nd,nxi))
   end if
else
   !check dimensions
   if(size(yi,dim) .ne. nxi .or. size(yi,dimd) .ne. nd)then
      write(stderr,*)"ERROR in interp1: yi has wrong dimensions:",size(yi,1),size(yi,2)
      stop
   end if
end if


select case(trim(method))
case('linear')! linear interpolation
   nw=2 !width of the interpolation window
   allocate(wght(nw))
   ii=0
   do,i=1,nxi
      !find the point directly under the interpolated point

      do while(xi(i) > x(ii+1))
         ii=ii+1
         if(ii+1>nx)exit
      end do
      if(ii .eq. 0)then
         if(xi(ii+1) .eq. x(i))ii=1 !special case when the start is on the boundary
      end if
      if(ii .eq. 0 .or. ii .eq. nx+1)then
         write(stderr,*)"ERROR in interp1: requested point is outside the data domain"
         stop
      end if
      !compute the weight
      dx=x(ii+1)-x(ii)

      wght(1)=(x(ii+1)-xi(i))/dx
      wght(2)=(xi(i)-x(ii))/dx
      !apply the weights to the data
      if(dim .eq. 1)then
         call dgemv('T',2,nd,1.d0,y(ii,1),nx,wght(1),1,0.d0,yi(i,1),nxi)
      else
         call dgemv('N',nd,2,1.d0,y(1,ii),nd,wght(1),1,0.d0,yi(1,i),1)

      end if
   end do
case default
   write(stderr,*)"ERROR in interp1: unknwon interpolation method",trim(method)
   stop
end select



end subroutine interp1

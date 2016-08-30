!program which constructs an east-west or north-south gradient regularization matrix expressed in spherical harmonics

!!Coded by Roelof Rietbroek, Mon Dec 15 11:20:36 2008
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!!Updated by Roelof Rietbroek, Mon May  4 16:43:42 2009
!!Added option for applying power law increase to the sides of the matrix ( geophysically justified?)

program SH_gradMAT
use shtlbx
use forttlbx
use binfiletools
implicit none
integer::i,j,m,l,l2,ind,shift,lmax,lmin,sz,itharg,startl
integer::stderr,radius,nbl,narg,colind,rowind,start,shiftold
type(BINdat)::output
logical::east_west,powerlaw
character(4)::tmp
character(200)::dum
double precision::ratio,scale1,scale2,scale3,scale4,scale5,deg1,deg2,pow
integer::iargc



!!defaults
lmax=0
lmin=1
output%file='stdout'
east_west=.true.
stderr=0
ratio=1.d0
!power law
powerlaw=.false.
pow=0.d0

!!process command line options
narg=iargc()
if(narg <1)call help()

itharg=0
do,i=1,narg
   itharg=itharg+1
   if(itharg > narg)exit
   call getarg(itharg,dum)
   if(dum(1:1) .eq. '-')then !argument is an option
      select case(dum(2:3))
      case('l=') !specify maximum and minimum degree
         ind=index(dum,',')
         if(ind .eq. 0)then
            read(dum(4:),*)lmax
         else
            read(dum(4:ind-1),*)lmax
            read(dum(ind+1:),*)lmin
         end if
      case('f=')!specify output file name
         if(dum(4:4) .eq. ' ')then
            itharg=itharg+1
            call getarg(itharg,output%file)
         else 
            output%file=dum(4:)
         end if
      case('s ')!north-south version
         east_west=.false.
      case('r=')!radius height of regularization (default takes the Earths surface)
         read(dum(4:),*)ratio
         ratio=RE/(1000.d0*ratio) ! factor
      case('p=')!additionally apply power law
         powerlaw=.true.
         read(dum(4:),*)pow
         pow=pow/2.d0
      case('h ','he')
         call help()
      case default
         write(stderr,*)'ERROR: unknown option selected:',dum(1:3)
         stop
      end select
   else!argument is not an option
   end if
end do



!input checks
if(lmax < lmin)then
   write(stderr,*)"ERROR: maximum degree must be larger or equal than minimum degree"
   stop
end if

If(ratio > 1.d0)then
   write(stderr,*)"ERROR: radius of regularization must be external of the Earth"
   stop
end If

!construct output system and META data
output%type='BDSYMVN_' !block diagonal system with symmetric blocks
output%mtyp='P'
output%nvec=0
output%nval1=SH_tpos(lmax,lmax,1,lmax,lmin)
output%nval2=output%nval1
output%nread=11
allocate(output%readme(output%nread))
output%nint=2
output%ndbls=2
allocate(output%ints(output%nint),output%dbls(output%ndbls),&
     output%ints_d(output%nint),output%dbls_d(output%ndbls))
allocate(output%side1_d(output%nval1))


output%descr="Regularization matrix made by SH_gradmat"
output%ints_d(1)="Lmax"
output%ints(1)=lmax
output%ints_d(2)="Lmin"
output%ints(2)=lmin
output%dbls_d(1)='Scale'
output%dbls(1)=pi*(GM/RE)*(GM/RE)
output%dbls_d(2)='Radius'
output%dbls(2)=(RE/ratio)/1000.d0



output%readme(1)="Regularization matrix derived from minimizing the following integral: "
output%readme(2)="   // lambda=2pi, theta=pi"
output%readme(3)="  ||"
if(east_west)then
   output%readme(4)="  ||   /         d V          \ 2    "
   output%readme(5)="  ||  | ---------------------  |    r^2 sin (theta) dlambda dtheta "
   output%readme(6)="  ||   \ r sin(theta) dlambda /      "
else
   output%readme(4)="  ||   /   d V    \ 2    "
   output%readme(5)="  ||  | ---------  |    r^2 sin (theta) dlambda dtheta "
   output%readme(6)="  ||   \  dtheta  /      "
end if
output%readme(7)="  ||"
output%readme(8)=" // lambda=0, theta=0"
output%readme(9)=" r is the radius of the regularization"
output%readme(10)=" V is the potential at r, theta (colatitude),lambda (longitude)"

if(powerlaw)then
   write(dum,*)pow
   output%readme(11)="The potential is scaled with a degree dependent power law: l^"//dum(1:5)
else
   output%readme(11)=""
end if


!packed block diagonal structure

output%pval2=1
output%nblocks=2*lmax+1
allocate(output%blockind(output%nblocks))


!zero order block (Cosine only)
sz=lmax-lmin+1
output%blockind(1)=sz
output%pval1=sz*(sz+1)/2


nbl=1
do,m=1,lmax ! loop over remaining order blocks
   sz=lmax-max(lmin,m)+1 ! side of the block
   nbl=nbl+1

   !Cosine part
   output%blockind(nbl)=output%blockind(nbl-1)+sz ! end location in the parameter space per block
   output%pval1=output%pval1+sz*(sz+1)/2 ! packed matrix size update

   !Sine part
   nbl=nbl+1
   output%blockind(nbl)=output%blockind(nbl-1)+sz ! end location in the parameter space per block
   output%pval1=output%pval1+sz*(sz+1)/2 ! packed matrix size update
end do

!construct side description data
ind=0
do,m=0,lmax ! loop over orders
   !COSINE part
   do,l=max(lmin,m),lmax ! loop over non-zero degrees
      ind=ind+1
      write(output%side1_d(ind),"(a4,i3,i3)")'GCN ',l,m
   end do
   
   !SINE part
   if(m > 0)then
      do,l=max(lmin,m),lmax ! loop over non-zero degrees
         ind=ind+1
         write(output%side1_d(ind),"(a4,i3,i3)")'GSN ',l,m
      end do
   end if

end do


!construct east-west matrix
allocate(output%pack1(output%pval1))
output%pack1=0.d0 !initialize all values to zero


if(east_west)then
   !fill matrix
   
   !zero order block is zero for East-west gradient (but store anyway for compatibility)
   nbl=1
   sz=output%blockind(nbl) ! size of the current block
   shift=sz*(sz+1)/2 ! update progression in the packed matrix
   
   do,m=1,lmax ! loop over non-zero orders
      nbl=nbl+1
      sz=output%blockind(nbl)-output%blockind(nbl-1) ! size of the current block
      !   write(stderr,*)"shift,block number,size",shift,nbl,sz
      if(m .eq. 1)then
         scale1=1/sqrt(8.d0)
      else
         scale1=1/4.d0
      end if
      
      start=max(lmin,m)
      do,l=start,lmax ! loop over column degrees
         colind=l-start+1 !index of the column
         !loop over degrees of the row degrees (with step 2 and only the upper triangle of the matrix)
         
         if(mod(l,2) .ne. mod(start,2))then ! skip rows which are even
            startl=start+1 ! align start index to the first parameter with the same parity
         else
            startl=start
         end if
         
         deg1=ratio**(l+1)
         scale2=((2*l+1)/dble((2*l-1)))*(l+m)*(l+m-1)
         scale3=((2*l+1)/dble((2*l-1)))*(l-m)*(l-m-1)
         
         do,l2=startl,l,2 !loop over non-vanishing overlap integral components
            rowind=l2-start+1
            
            deg2=ratio**(l2+1)*deg1
            
            scale4=sqrt(scale2*((2*l2+1)/dble((2*l2-1)))*(l2-m)*(l2-m-1))
            scale5=sqrt(scale3*((2*l2+1)/dble((2*l2-1)))*(l2+m)*(l2+m-1))
            !         write(stderr,*)'l m index',l,m,shift+rowind+colind*(colind-1)/2
            output%pack1(shift+rowind+colind*(colind-1)/2)=scale1*deg2*scale4*SH_overlapint(l-1,m,l2-1)+&
                 scale1*deg2*scale5*SH_overlapint(l2-1,m,l-1)
            !         if(l .eq. l2)write(stderr,*)'First',m,l,l2,SH_overlapint(l-1,m,l2-1),SH_overlapint(l2-1,m,l-1),output%pack1(shift+rowind+colind*(colind-1)/2)
         end do
         
         !update diagonal part
         ind=shift+colind*(colind+1)/2 ! index of digonal element in packed matrix
         output%pack1(ind)=output%pack1(ind)+deg1*deg1*(scale2+scale3)
         !     write(stderr,*)'Added',m,l,l2,output%pack1(ind)
      end do
      
      !update shift parameter
      shiftold=shift
      shift=shift+sz*(sz+1)/2
      
      !copy values in the SINE part of the matrix
      
      ind=0
      do,i=1,sz*(sz+1)/2
         ind=ind+1
         output%pack1(shift+ind)=output%pack1(shiftold+ind)
      end do
      
      !again update progression in matrix:
      shift=shift+sz*(sz+1)/2
      nbl=nbl+1 ! shift block number to SINE block   
   end do

else !Construct North_south matrix if requested
   !fill matrix
   
   

   !zero order block is zero for East-west gradient and diagonal for North-South gradient)
   !It is treated as a special case
   shift=0
   nbl=1
   sz=output%blockind(nbl)! size of the current block
   start=max(lmin,0)
   do,l=start,lmax
      colind=l-start+1
      !diagonal part only
      ind=shift+colind*(colind+1)/2 ! index of diagonal element in packed matrix
      deg1=ratio**(l+1)
      scale2=((l)*(l+1))*4
      output%pack1(ind)=deg1*deg1*scale2
   end do

   shift=shift+sz*(sz+1)/2 ! update progression in the packed matrix
  
   
   
   do,m=1,lmax ! loop over non-zero orders
      nbl=nbl+1
      sz=output%blockind(nbl)-output%blockind(nbl-1) ! size of the current block
      !   write(stderr,*)"shift,block number,size",shift,nbl,sz
      if(m .eq. 1)then
         scale1=1/sqrt(8.d0)
      else
         scale1=1/4.d0
      end if
      
      start=max(lmin,m)
      do,l=start,lmax ! loop over column degrees
         colind=l-start+1 !index of the column
         !loop over degrees of the row degrees (with step 2 and only the upper triangle of the matrix)
         
         if(mod(l,2) .ne. mod(start,2))then ! skip rows which are even
            startl=start+1 ! align start index to the first parameter with the same parity
         else
            startl=start
         end if
         
         deg1=ratio**(l+1)
         scale2=dble((l+m)*(l-m+1))
         scale3=dble((l-m)*(l+m+1))
         
         do,l2=startl,l,2 !loop over non-vanishing overlap integral components
            rowind=l2-start+1
            
            deg2=ratio**(l2+1)*deg1
            
            scale4=sqrt(scale2*(l2-m)*(l2+m+1))
            scale5=sqrt(scale3*(l2+m)*(l2-m+1))
            !         write(stderr,*)'l m index',l,m,shift+rowind+colind*(colind-1)/2
            output%pack1(shift+rowind+colind*(colind-1)/2)=-scale1*deg2*scale4*SH_overlapint(l,m,l2)-&
                 scale1*deg2*scale5*SH_overlapint(l2,m,l)
            !         if(l .eq. l2)write(stderr,*)'First',m,l,l2,SH_overlapint(l-1,m,l2-1),SH_overlapint(l2-1,m,l-1),output%pack1(shift+rowind+colind*(colind-1)/2)
         end do
         
         !update diagonal part
         ind=shift+colind*(colind+1)/2 ! index of digonal element in packed matrix
         output%pack1(ind)=output%pack1(ind)+deg1*deg1*(scale2+scale3)
         !     write(stderr,*)'Added',m,l,l2,output%pack1(ind)
      end do
      
      !update shift parameter
      shiftold=shift
      shift=shift+sz*(sz+1)/2
      
      !copy values in the SINE part of the matrix
      
      ind=0
      do,i=1,sz*(sz+1)/2
         ind=ind+1
         output%pack1(shift+ind)=output%pack1(shiftold+ind)
      end do
      
      !again update progression in matrix:
      shift=shift+sz*(sz+1)/2
      nbl=nbl+1 ! shift block number to SINE block   
   end do

end if


!!!!!OPTIONAL power law normalization

if(powerlaw)then
   nbl=1
   sz=output%blockind(nbl)! size of the current block
   shift=0
   if (.not. east_west)then ! also scale zero order block ( diagonal only) in case of north_south gradient

      start=max(lmin,0)
      do,l=start,lmax
         colind=l-start+1
         !diagonal part only
         ind=shift+colind*(colind+1)/2 ! index of diagonal element in packed matrix
         deg1=l**(pow) ! power law scale
         output%pack1(ind)=deg1*deg1*output%pack1(ind)
      end do
   end if
   
   shift=shift+sz*(sz+1)/2 ! update progression in the packed matrix


   do,m=1,lmax ! loop over non-zero orders
      nbl=nbl+1
      sz=output%blockind(nbl)-output%blockind(nbl-1) ! size of the current block
      

      start=max(lmin,m)
      do,l=start,lmax ! loop over column degrees
         colind=l-start+1 !index of the column
         !loop over degrees of the row degrees (with step 2 and only the upper triangle of the matrix)
         
         if(mod(l,2) .ne. mod(start,2))then ! skip rows which are even
            startl=start+1 ! align start index to the first parameter with the same parity
         else
            startl=start
         end if
         
         deg1=l**(pow)
         
         do,l2=startl,l,2 !loop over non-vanishing overlap integral components
            rowind=l2-start+1
            deg2=l2**(pow)
            ind=shift+rowind+colind*(colind-1)/2
            output%pack1(ind)=deg1*deg2*output%pack1(ind)
         end do
         
      end do
      
      !update shift parameter
      shiftold=shift
      shift=shift+sz*(sz+1)/2
      
      !copy values in the SINE part of the matrix
      
      ind=0
      do,i=1,sz*(sz+1)/2
         ind=ind+1
         output%pack1(shift+ind)=output%pack1(shiftold+ind)
      end do
      
      !again update progression in matrix:
      shift=shift+sz*(sz+1)/2
      nbl=nbl+1 ! shift block number to SINE block   
   end do

end if




! !MAT_N-S = N_surf - MAT_E-W 
! !where N_surf is the diagonal regularization matrix which regularizes the surface gradient of the potential, with on the diagonal:
! !N_surf = 4 pi (GM/R)^2 l(l+1)


! if(.not. east_west)then
   
! !change sign of the matrix
!    do,i=1,output%pval1
!       output%pack1(i)=-output%pack1(i)
!    end do
   
!    !Add diagonal  degree dependent matrix
!    !zero order part

!    do,l=lmin,lmax
!       deg1=ratio**(2*l+2)*4.d0
!       rowind=l-lmin+1
!       ind=rowind*(rowind+1)/2
!       output%pack1(ind)=output%pack1(ind)+deg1*l*(l+1)
!    end do
!    nbl=1
!    sz=output%blockind(nbl)-output%blockind(nbl-1)
!    shift=sz*(sz+1)/2

!    do,m=1,lmax !loop over remaining orders

!       start=max(lmin,m)

!       !Cosine block
!       nbl=nbl+1 ! next block
!       !size
!       sz=output%blockind(nbl)-output%blockind(nbl-1)
      
!       do,l=start,lmax
!          deg1=ratio**(2*l+2)*4.d0
!          rowind=l-start+1
!          ind=rowind*(rowind+1)/2
!          output%pack1(ind+shift)=output%pack1(ind+shift)+deg1*l*(l+1)
!       end do

!       !Copy to Sine Block
!       shiftold=shift
!       shift=shift+sz*(sz+1)/2
!       nbl=nbl+1
      
!       do,l=start,lmax
!          rowind=l-start+1
!          ind=rowind*(rowind+1)/2
!          output%pack1(ind+shift)=output%pack1(ind+shiftold)
!       end do
!       shift=shift+sz*(sz+1)/2
   

!    end do

! end if





!Write matrix to file/standard output


call write_BINtype(output) 


end program SH_gradMAT

subroutine help()
implicit none
character(8)::frmt
integer::unit
frmt='(A)'
unit=6

write(unit,frmt)"Program SH_gradmat constructs a East-West or North-South Potential Gradient regularization matrix"
write(unit,frmt)"Usage: SH_gradmat -l=LMAX[,LMIN] [ OPTIONS ]"
write(unit,frmt)""
write(unit,frmt)"Command line options:"
write(unit,frmt)"  -l=LMAX[,LMIN]: Specify maximum (LMAX) and optionally minimum (LMIN) degree for the output matrix"
write(unit,frmt)"   LMIN has a default of 1"
write(unit,frmt)"  -f=FILEOUT: set the output file name, default writes to standard output"
write(unit,frmt)"  -s: construct North-South gradient matrix (default is an East-West gradient matrix)"
write(unit,frmt)"  -r=RADIUS: Takes a different radius where the gradient minimization is applied"
write(unit,frmt)"   The radius is measured in kms from the Earths center. Default takes the Earths radius"
write(unit,frmt)"  -p=POW: Additionally apply power law to the gradient matrix"
write(unit,frmt)"   The matrix PHI will be scaled to the new matrix PHI_pow according to: "
write(unit,frmt)"   PHI_pow = D PHI D. Where D is a diagonal matrix containing entries with l^(POW/2)"


! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""
! write(unit,frmt)""

stop
end subroutine help

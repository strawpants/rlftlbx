!!Coded by Roelof Rietbroek, Wed Mar 18 22:00:39 2015
!!UniBonn IGG,
!!Nussallee17,Bonn ,
!!53115,Bonn,Germany
!!email: roelof@geod.uni-bonn.de

module AutoRegressive
  implicit none

  type ARmodel
     integer::ord,ndat
     double precision,pointer::Para(:)=>null()
     double precision::sig2
     !autocovariance matrix in banded form
     double precision,pointer::autocovBD(:,:)=>null()
  end type ARmodel
contains
  !compute an autoregressive model from a series
  subroutine ComputeArYW(input,Ar)
    double precision,dimension(:),intent(in)::input
    type(ArModel),intent(inout):: Ar
    double precision::EmpAutoCov(0:Ar%ord)
    double precision::mat(Ar%ord,Ar%ord)
    double precision::gam(0:Ar%ord)
    double precision::mean
    integer::lag,j,i,ord
    integer::info
    
    ord=Ar%ord
    if(ord .eq. 0)then
       write(0,*)"ERROR: AR order must be larger than 0"
       stop
    end if
    Ar%ndat=size(input,1)
    allocate(Ar%Para(ord))
    gam=0.d0

    mean=0.d0
    do,i=1,Ar%ndat
       mean=mean+input(i)/Ar%ndat
    end do
    
    !compute the necessary empirical autocovariances
    do,lag=0,ord
       do,j=1+lag,Ar%ndat
          gam(lag)=gam(lag)+((input(j)-mean)*(input(j-lag)-mean))/(Ar%ndat-1)
       end do
 !      write(0,*)lag,gam(lag)
    end do

    !Set up Yule Walker Equations
    do,i=1,ord
       do,j=1,i
          mat(j,i)=gam(i-j)
       end do
        Ar%para(i)=gam(i)
    end do

    !solve for the AR parameters
    call dpotrf ('U', ord, mat(1,1),ord, info)
    if (info .ne. 0)then
       write(0,*)"ERROR: Cholesky inversion failed in computeARYW()",info
       stop
    end if
    call dpotrs('U',ord,1,mat(1,1),ord,Ar%para(1),ord,info)
    if (info .ne. 0)then
       write(0,*)"ERROR: computeARYW() failed",info
       stop
    end if

    
    !solve the  sig^2
    Ar%sig2=gam(0)
    do,i=1,ord
       Ar%sig2=Ar%sig2-Ar%para(i)*gam(i)
    end do
    
       
  end subroutine ComputeArYW
  
  !function returning an banded symetric lower triangular toeplitz matrix in lapack/blas format
  subroutine makeautoCovBD(Ar)
    type(Armodel),intent(inout)::Ar
    double precision,pointer::autocorr(:)=>null()
    double precision,pointer::autocorrUpdate(:)=>null()
    integer::i,j,iter, maxiter
    double precision::denom,gam,corrdif,thres

    
    allocate(Ar%autocovBD(Ar%ord+1,Ar%ndat))
    Ar%autocovBD=0.d0
    
    allocate(autocorr(0:Ar%ord),autocorrUpdate(0:Ar%ord))
    autocorr=0.d0
    autocorr(0)=1.d0
    autocorrUpdate(0)=1.d0
    corrdif=1.d0
    thres=1e-4
    maxiter=20
    iter=0
    
    !estimate correlation function recursively
    do while (corrdif > thres .and. iter <maxiter)
       corrdif=0.d0
       do,i=1,Ar%ord
          autocorrUpdate(i)=0.d0
          do,j=1,Ar%ord
             autocorrUpdate(i)=autocorrUpdate(i)+Ar%para(j)*autocorr(abs(j-i))
          end do
          corrdif=max(corrdif,abs(autocorrUpdate(i)-autocorr(i)))
          
       end do
       iter=iter+1
       autocorr=autocorrUpdate
    end do
    
    if(iter == maxiter)then
       write(0,*)"WARING:autocorrelation has reached maximum iteration",maxiter
    end if

    !estimate gam(0)
    denom=1.d0
    do,i=1,Ar%ord

       denom=denom-Ar%para(i)*autocorr(i)
    end do
    
    gam=Ar%sig2/denom

    !fill up banded matrix
    do,i=0,Ar%ord
       Ar%autocovBd(i+1,:)=gam*autocorr(i)
!       write(0,*)Ar%autocovBd(i+1,1),autocorr(i),iter
    end do
    
  end subroutine makeautoCovBD
  
end module AutoRegressive

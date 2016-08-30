program freenet

! -------------------------------------------------------------------------
! Purpose:    Set up free network restriction to be added to the
!             normal equation system
!
! Author:     M. Fritsche
!
! -------------------------------------------------------------------------
  implicit none

! Local parameters
! ----------------
  integer, parameter :: lfnsnx = 10
  integer, parameter :: maxsta = 1000
  real,    parameter :: aell = 6371.d3
  real,    parameter :: sig0 = 1.d0
  real,    parameter :: zero = 0.d0
  real,    parameter :: one  = 1.d0

! Local typ
! ---------
  type t_xyz
    real,    dimension(3) :: x
    integer, dimension(3) :: i
  end type t_xyz

! Local variables
! ---------------
  character(len=206) :: line
  character(len=4)   :: partyp

  integer :: i
  integer :: k
  integer :: ik
  integer :: ios
  integer :: npar
  integer :: nsta
  integer :: ipar
  integer :: ista1
  integer :: ista2
  integer :: nwrite

  real :: value
  real :: X
  real :: Y
  real :: Z
  real, dimension(7,7) :: BTB
  real, dimension(3,7) :: Bmat1
  real, dimension(3,7) :: Bmat2
  real, dimension(7,7) :: Pmat
  real, dimension(3,3) :: regMat
  real, dimension(7)   :: sigma
  real, dimension(:), allocatable :: neq

  type(t_xyz), dimension(maxsta) :: snxcrd

  logical :: crdsav = .false.
  logical :: wrtflg

! Read sinex file
! ---------------
  open (unit=lfnsnx,file='./input.snx',status='OLD',form='FORMATTED',iostat=ios)
  npar = 0
  nsta = 0
  do
    read (lfnsnx,'(A)',iostat=ios) line
    if     (line(1:18) .eq. '-SOLUTION/ESTIMATE' .or. ios .ne. 0) then
      exit
    elseif (line(1:18) .eq. '+SOLUTION/ESTIMATE') then
      crdsav = .true.
      cycle
    end if
    if (line(1:1) .eq. '*' .or.  .not.crdsav) cycle
    read (line,'(1X,I5,1X,A4,36X,e21.15)') ipar, partyp, value
    npar = max(ipar,npar)
    if (partyp .eq. 'STAX') nsta = nsta + 1
    if     (partyp .eq. 'STAX') then
      snxcrd(nsta)%x(1) = value
      snxcrd(nsta)%i(1) = ipar
    elseif (partyp .eq. 'STAY') then
      snxcrd(nsta)%x(2) = value
      snxcrd(nsta)%i(2) = ipar
    elseif (partyp .eq. 'STAZ') then
      snxcrd(nsta)%x(3) = value
      snxcrd(nsta)%i(3) = ipar
    end if
  end do
  close (unit=lfnsnx)

! Create the (B'B) Matrix
! -----------------------
  BTB = 0.d0
  do ista1 = 1, nsta
    X  = snxcrd(ista1)%x(1) / aell
    Y  = snxcrd(ista1)%x(2) / aell
    Z  = snxcrd(ista1)%x(3) / aell

    Bmat1 = RESHAPE(SOURCE=(/  one, zero, zero,     &
                              zero,  one, zero,     &
                              zero, zero,  one,     &
                              zero,    Z,   -Y,     &
                                -Z, zero,    X,     &
                                 Y,   -X, zero,     &
                                 X,    Y,    Z  /), SHAPE=(/ 3, 7 /) )

    BTB = BTB + MATMUL( TRANSPOSE(Bmat1), Bmat1 )
  end do !ista1

! Set sigma for each pseudoobservation
! ------------------------------------
  sigma = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)

! Compute the weight matrix (B'B)^-1 * P * (B'B)^-1
! -------------------------------------------------
  Pmat = 0.d0
  do i = 1, 7
    if (sigma(i) /= 0.0) then
      Pmat(i,i) = ( sig0 / abs(sigma(i)) )**2
    else
      do k = 1, 7
        BTB(i,k) = 0.d0
        BTB(k,i) = 0.d0
      end do
      BTB(i,i) = 1.d0
    end if
  end do

  if ( matinv(BTB, 7) == 0 ) then
    DO i = 1, 7
      if (sigma(i) == 0.0) then
        BTB(i,i) = 0.0    
      end if
    end do
    Pmat = MATMUL( BTB, MATMUL( Pmat, BTB ) )
  end if

! Allocate upper triangular NEQ
! -----------------------------
  allocate(neq(npar*(npar+1)/2))
  neq = 0.d0

! Create the regularization matrix and add the weights
! ----------------------------------------------------
  do ista1 = 1, nsta
    X  = snxcrd(ista1)%x(1) / aell
    Y  = snxcrd(ista1)%x(2) / aell
    Z  = snxcrd(ista1)%x(3) / aell

    Bmat1 = RESHAPE(SOURCE=(/  one, zero, zero,     &
                              zero,  one, zero,     &
                              zero, zero,  one,     &
                              zero,    Z,   -Y,     &
                                -Z, zero,    X,     &
                                 Y,   -X, zero,     &
                                 X,    Y,    Z  /), SHAPE=(/ 3, 7 /) )

    do ista2 = 1, nsta
      X  = snxcrd(ista2)%x(1) / aell
      Y  = snxcrd(ista2)%x(2) / aell
      Z  = snxcrd(ista2)%x(3) / aell

      Bmat2 = RESHAPE(SOURCE=(/  one, zero, zero,     &
                                zero,  one, zero,     &
                                zero, zero,  one,     &
                                zero,    Z,   -Y,     &
                                  -Z, zero,    X,     &
                                   Y,   -X, zero,     &
                                   X,    Y,    Z  /), SHAPE=(/ 3, 7 /) )

      regMat = MATMUL( Bmat1, MATMUL( Pmat, TRANSPOSE(Bmat2) ) )

      do i = 1,3
        do k = 1,3
          IF (snxcrd(ista2)%i(k) < snxcrd(ista1)%i(i)) cycle
          ik = ikf(snxcrd(ista1)%i(i),snxcrd(ista2)%i(k))
          neq(ik) = neq(ik) + regMat(i,k)
        end do
      end do

    end do
  end do

! Write Solution/Matrix_Apriori Block
! -----------------------------------
  open (unit=lfnsnx,file='./freenet.snx',status='UNKNOWN',form='FORMATTED',iostat=ios)

  write (lfnsnx,'(A)') '+SOLUTION/MATRIX_APRIORI L COVA'
  write (lfnsnx,'(A)') &
    '*PARA1 PARA2 ____PARA2+0__________ ____PARA2+1__________ ____PARA2+2__________'

  nWrite = 0
  line   = ''
  wrtflg = .false.
  do i = 1, npar
    do k = 1, i
      ik = ikf(i,k)
      if ( neq(ik) .ne. 0.d0 ) wrtflg = .true.

      if     (nWrite .eq. 0) then
        write (line( 1:34), '(2(1X,I5), 1X,E21.14)') i, k, neq(ik)
      elseif (nWrite .eq. 1) then
        write (line(35:56), '(1X,E21.14)') neq(ik)
      else
        write (line(57:78), '(1X,E21.14)') neq(ik)
      end if

      nWrite = nWrite + 1

      if (nWrite == 3 .or. i == k) then
        if ( wrtflg ) write (lfnsnx,'(A)') line(1:LEN_TRIM(line))
        nWrite = 0
        line   =  ''
        wrtflg = .false.
      end if

    end do
  end do
  write (lfnsnx,'(A)') '-SOLUTION/MATRIX_APRIORI L COVA'
  close (unit=lfnsnx)

! Free memory
! -----------
  deallocate(neq)

! Exit program
! ------------
  call exit(0)

! Definition of additional subroutines / functions ------------------------

contains

function matinv(aa,nn)

! -------------------------------------------------------------------------
!
! Purpose:    Compute matrix inversion using Gauss-Jordan elimination.
!             Returns 1 if matrix is singular and 0 otherwise.
!
! -------------------------------------------------------------------------
  implicit none

! Return value
! ------------
  integer                :: matinv

! List of Parameters
! ------------------
  integer                :: nn
  real, dimension(nn,nn) :: aa

! Local Variables
! ---------------
  real                   ::  big,dum,pivinv

  integer                :: ii,icol,irow,jj,kk,ll,mm
  integer, dimension(nn) :: indxc
  integer, dimension(nn) :: indxr
  integer, dimension(nn) :: ipiv

  ipiv(:) = 0
  matinv  = 0

  do ii = 1, nn
    big = 0.d0

    do jj = 1, nn
      if (ipiv(jj) /= 1) then
        do kk = 1, nn
          if (ipiv(kk) == 0) then
            if ( ABS(aa(jj,kk)) >= big) then
              big  = ABS(aa(jj,kk))
              irow = jj
              icol = kk
            end if
          else if (ipiv(kk) > 1) then
            matinv = 1
            return
          end if
        end do
      end if
    end do

    ipiv(icol) = ipiv(icol)+1
    if (irow /= icol) then
      do ll = 1, nn
        dum = aa(irow,ll)
        aa(irow,ll) = aa(icol,ll)
        aa(icol,ll) = dum
      end do
    end if
    indxr(ii) = irow
    indxc(ii) = icol
    if (ABS(aa(icol,icol)) < 1.d-10) then
      matinv = 1
      return
    end if

    pivinv        = 1.d0 / aa(icol,icol)
    aa(icol,icol) = 1.d0
    aa(icol,:)    = aa(icol,:)*pivinv

    do mm = 1, nn
      if (mm /= icol) then
        dum = aa(mm,icol)
        aa(mm,icol) = 0.d0
        aa(mm,:) = aa(mm,:) - aa(icol,:) * dum
      end if
    end do
  end do

  do ll = nn, 1, -1
    if (indxr(ll) /= indxc(ll)) then
      do kk = 1,nn
        dum = aa(kk,indxr(ll))
        aa(kk,indxr(ll)) = aa(kk,indxc(ll))
        aa(kk,indxc(ll)) = dum
      end do
    end if
  end do

end function matinv

! -------------------------------------------------------------------------

function ikf(i,k)
  implicit none
  integer :: i, k, ikf

  if (k .gt. i) then
    ikf = k*(k-1)/2+i
  else
    ikf = i*(i-1)/2+k
  end if
end function ikf

! -------------------------------------------------------------------------

end program freenet

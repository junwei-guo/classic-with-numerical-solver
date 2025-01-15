!Module contains various useful numerical solvers, J. Guo, University of Calgary, Aug-21, 2024
!Define real kind and integer kind
module dataType
  implicit none
  
  integer, parameter :: sp    = selected_real_kind(6)  ! Single precision
  integer, parameter :: dp    = selected_real_kind(15) ! Double precision
  integer, parameter :: int8  = selected_int_kind(2)     ! 8-bit integer (range -128 to 127)
  integer, parameter :: int16 = selected_int_kind(4)     ! 16-bit integer (range -32768 to 32767)
  integer, parameter :: int32 = selected_int_kind(9)     ! 32-bit integer (range about -2.1E9 to 2.1E9)
  integer, parameter :: int64 = selected_int_kind(18)    ! 64-bit integer (range about -9.2E18 to 9.2E18)


  !Default kinds for real and integer
  integer, parameter :: rk    = selected_real_kind(15)  ! equivalent to real(8)
  integer, parameter :: ik    = selected_int_kind(9)
contains
  
end module  dataType


module numSolver
  use, intrinsic :: ieee_arithmetic
  use dataType
  type :: GenOdeSolver
    ! integer   :: neqn   = 2  
    integer   :: iflag  = 1
    integer ( kind = 8 ) :: iwork(5)
    real(rk),dimension(:),allocatable :: work
    real(rk) :: abserr   = 0.00001D+00
    real(rk) :: relerr   = 0.00001D+00
    contains
    !procedure :: Solve
  end type GenOdeSolver
    
  type :: rkSolver
    integer   :: neqn   = 2

  
    contains

  end type rkSolver

  type :: scalarSolver
    real(rk),dimension(2)   :: IntervalBounds
    real(rk)                :: InitialGuess
    real(rk)                :: DDfDDx_max               !Upper bound of the second derivative
    real(rk)                :: MachineError = 2.0D-16   !Machine error found on a x64 intel PC.
    real(rk)                :: errE         = 1.0D-2    !Default: 1.0D-8
    real(rk)                :: errT         = 1    !Default: 1.0D-8
    real(rk)                :: x,fx                     
    contains
    procedure :: CheckMachineError
    procedure :: FindMinimum
    procedure :: FindRoot
  end type scalarSolver


  

  contains
 
  subroutine FindMinimum(this,func)
    class(scalarSolver)    :: this
    real(rk), external :: func
    real(rk), external :: glomin
    this%fx = glomin(    this%IntervalBounds(1), this%IntervalBounds(2), this%InitialGuess, this%DDfDDx_max, this%MachineError, this%errE,this%errT,func,this%x)


  end subroutine FindMinimum

  subroutine FindRoot(this,func)
    class(scalarSolver)    :: this
    real(rk), external :: func
    real(rk), external :: zero


    this%x = zero(    this%IntervalBounds(1), this%IntervalBounds(2), this%MachineError, this%errT,   func)

    this%fx = func(this%x)

    if (ieee_is_nan(this%x)) then
      print *, 'x is NaN'
    end if
  end subroutine FindRoot
  


  subroutine CheckMachineError(this)
    class(scalarSolver)    :: this
    real(rk) :: epsilon, one
    epsilon = 1.0
    one = 1.0
    do
        epsilon = epsilon / 2.0
        if (one + epsilon == one) exit
    end do
    this%MachineError = 2.0*epsilon
  end subroutine CheckMachineError


end module numSolver

module numDiag 
  use dataType
  implicit none
  public  :: initializeEnergBalDiag
  public  :: printEnergBalDiag

  real(rk)    :: timeTotal       !Total program excution time
  real(rk)    :: energBalTol = 0.01
  real(rk)    :: nItrEnergBalVeg, nItrEnergBalNoVeg
  real(rk)    :: timeEnergBalVeg, timeEnergBalNoVeg
  real(rk)    :: nSolEnergBalVeg, nSolEnergBalNoVeg 
  real(rk)    :: resiEnergBalVeg, resiEnergBalNoVeg
  contains
  

  subroutine initializeEnergBalDiag()
    implicit none 

    nItrEnergBalNoVeg = 0.0     !Total number of iterations
    nItrEnergBalVeg   = 0.0
    timeEnergBalNoVeg = 0.0     !Total time in solving
    timeEnergBalVeg   = 0.0

    nSolEnergBalVeg   = 0.0     !Number of solver called
    nSolEnergBalNoVeg = 0.0     

    resiEnergBalVeg   = 0.0     !Summation of the residual energy squared
    resiEnergBalNoVeg = 0.0

    timeTotal         = 0.0

  end subroutine initializeEnergBalDiag


  subroutine printEnergBalDiag()
    implicit none 

    print *,  "energBalTol           = "//trim(dou2str(energBalTol))//";"
    print *,  "timeTotal (s)         = "//trim(dou2str(timeTotal))//";"
    print *,  "nItrEnergBalVeg (M)   = "//trim(dou2str(nItrEnergBalVeg/1000000.0))//"; nItrEnergBalNoVeg (M) = "//trim(dou2str(nItrEnergBalNoVeg/1000000.0))//";"
    print *,  "timeEnergBalVeg (s)   = "//trim(dou2str(timeEnergBalVeg))//"; timeEnergBalNoVeg (s) = "//trim(dou2str(timeEnergBalNoVeg))//";"
    print *,  "nSolEnergBalVeg (M)   = "//trim(dou2str(nSolEnergBalVeg/1000000.0))//"; nSolEnergBalNoVeg (M) = "//trim(dou2str(nSolEnergBalNoVeg/1000000.0))//";"
    print *,  "resiEnergBalVeg       = "//trim(dou2str(resiEnergBalVeg))//"; resiEnergBalNoVeg     = "//trim(dou2str(resiEnergBalNoVeg))//";"

  end subroutine printEnergBalDiag

  function dou2str(number) result(temp_str)
  implicit none
  real(rk), intent(in) :: number
  character(len=256) :: temp_str,str
  ! Write the number to a string

  write(temp_str , *) number
  
  end function dou2str
end module numDiag 





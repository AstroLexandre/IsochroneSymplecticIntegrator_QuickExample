!==================================================================
!  File name : run_IsochroneSymplecticIntegrator_QuickExample.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     A simple example demonstrating the use of isochrone
!     splitting for a test particle evolving in a 
!     Plummer potential. The output file gives, at each time 
!     step, the position, velocity, and the energy 
!     of the test particle.
!     This program is provided purely for educational purposes.
!     You are welcome to reuse, incorporate, or modify it in
!     your own work, provided that we are cited.
!
!  Author              : Alexandre Bougakov
!  Email               : Alexandre.Bougakov@obspm.fr
!  Creation date       : 2026-02-20
!  Last modification   : 2026-02-27
!  Version             : 1.0
!
!  License :
!     This code is distributed under the MIT license.
!     See the LICENSE file for more details.
!
!  Attribution       :
!     If you use this code, please cite the corresponding 
!     publication.
!  Reference article : Bougakov et al., A&A, 697, A106 (2025)
!  DOI               : 10.1051/0004-6361/202553886
!
!  Copyright (c) 2026 Alexandre Bougakov
!==================================================================

program run_IsochroneSplitting_QuickExample
  
  use MathConstant
  use Isochrone
  use Plummer
  use HamiltonianSplitting
  implicit none

  character(256):: paramFile ! namelist containing the initialisation of parameters
  character(256):: outFile   ! output main file contain position, velocity & energy at each time
  integer:: Nout             ! register number of the output file
  integer(kind=8):: Nstep    ! number of time step realised
  integer:: Nsave            ! number of timesteps between each output
  real(prc):: X(6), t        ! state vector and time
  real(prc):: ti, tf, Dt     ! time start, time end, time step
  real(prc):: eta, kappa     ! potential parameters
  real(prc):: mu, b          ! isochrone splitting parameters
  character(20):: hdrfmt     ! header format
  character(22):: strfmt     ! output format

  ! initialisation
  eta = 0.0_prc
  kappa = 0.0_prc
  X(:) = 0.0_prc
  ti = 0.0_prc
  tf = 0.0_prc
  Dt = 0.0_prc
  Nsave = 0
  mu = 0.0_prc
  b = 0.0_prc

  ! setting the formats
  write(hdrfmt,'(A,I0,A,I0,A)') "(A,(X,A", precision(1.0_prc)+6, "),7(X,A", precision(1.0_prc)+8, "))"
  write(strfmt,'(5(A,I0))') "(ES", precision(1.0_prc)+8, ".", precision(1.0_prc)+1, ",7(X,ES", precision(1.0_prc)+8, ".", precision(1.0_prc)+1, "))"

  if(COMMAND_ARGUMENT_COUNT() == 2) then
     call GET_COMMAND_ARGUMENT(1,paramFile)
     call GET_COMMAND_ARGUMENT(2,outFile)
  else
     write(0,*) "Use: <program name> <paramFile.txt> <outFile.dat>"
     STOP
  endIf

  call SetParam(paramFile)

  !----- isochrone check ----
  ! outFile = "out/IsochroneTest.dat"
  ! eta = 1.0_prc
  ! kappa = 0.5_prc
  ! mu = 1.0_prc
  ! b = 0.5_prc
  !----------------------

  if(abs(mu) < epsilon(1.0_prc)) then
    mu = 0.0_prc
    write(0,*) "WARNING: mu = 0 assumed."
  end if

  ! opening the output file
  open(newunit=Nout,file=trim(outFile),status='REPLACE',action='WRITE')
  write(Nout,hdrfmt) "#", "t", "x", "y", "z", "vx", "vy", "vz", "H"
  ! start the integration
  t = ti
  Nstep = 0_8
  do while(abs(t - ti) <= tf)
    ! output
    if(Nsave >= 0 .and. mod(Nstep,Nsave) == 0) write(Nout,strfmt) t, X(:), Hamiltonian(X, t, potWrap)
    ! performs one time step
    call integstep_symplectic(X, t, Dt, mu, b, DpotWrap)

    Nstep = Nstep + 1_8
  end do
  ! stop the integration
  write(Nout,strfmt) t, X(:), Hamiltonian(X, t, potWrap)

  ! closing the output file
  close(Nout)

contains

  function potWrap(X,t)
    real(prc):: potWrap
    real(prc),intent(in):: X(6), t

    ! set the potential
    !----- test            ----
    potWrap = PotentialPlummer(X,t,eta,kappa)
    !----- isochrone check ----
    ! potWrap = PotentialIsochrone(X,t,eta,kappa)
    !----------------------

  end function potWrap
  !
  !----------------------------------------------------
  !
  function DpotWrap(X,t)
    real(prc):: DpotWrap(6)
    real(prc),intent(in):: X(6), t

    ! set the Dpotential
    !----- test            ----
    DpotWrap(:) = DPotentialPlummer(X,t,eta,kappa)
    !----- isochrone check ----
    ! DpotWrap(:) = DPotentialIsochrone(X,t,eta,kappa)
    !----------------------

  end function DpotWrap

  !---------------------------------------------------------------------------------------------
  ! Reading and setting the parameters
  !
  subroutine SetParam(paramFile)
    character(*),intent(in):: paramFile

    integer:: Nin               ! register number of the input files
    integer:: ios
    character(256):: line       ! current line read from the input file
    character(256):: key, value ! parameter name extracted from the line and the corresponding parameter value as a string

    open(newunit=Nin, file=trim(paramFile), status='old', action='read')

    do
      read(Nin, '(A)', iostat=ios) line
      if(ios /= 0) exit
      line = adjustl(line)

      if(line(1:1) == '#' .or. trim(line) == '') cycle

      read(line, '(A)', iostat=ios) line
      if(index(line,'=') == 0) cycle

      key = adjustl(line(1:index(line,'=')-1))
      value = adjustl(line(index(line,'=')+1:))

      select case(trim(key))
        ! set the potential parameters
        case('eta')
          read(value,*) eta
        case('kappa')
          read(value,*) kappa

        ! set the state vector (x, y, z, vx, vy, vz): position and velocity in Cartesian coordinates
        case('x')
          read(value,*) X(1)
        case('y')
          read(value,*) X(2)
        case('z')
          read(value,*) X(3)
        case('vx')
          read(value,*) X(4)
        case('vy')
          read(value,*) X(5)
        case('vz')
          read(value,*) X(6)

        ! set the integration parameters
        case('ti')
          read(value,*) ti
        case('tf')
          read(value,*) tf
        case('Dt')
          read(value,*) Dt

        ! save every Nsave points
        case('Nsave')
          read(value,*) Nsave

        ! set the isochrone splitting parameters
        case('mu')
          read(value,*) mu
        case('b')
          read(value,*) b

        case default
          write(0,*) 'FATAL ERROR: unknown parameter -> ', trim(key)
          write(0,*) "Program stopped."
          STOP
      end select
    end do

    close(Nin)

  end subroutine SetParam
  !---------------------------------------------------------------------------------------------
end program run_IsochroneSplitting_QuickExample
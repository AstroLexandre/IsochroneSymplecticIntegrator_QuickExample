!==================================================================
!  File name : HamiltonianSplitting.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     The Hamiltonian is split into two parts: H = A + 𝜀B,
!     where 𝜀 must be small to take advantage 
!     of hierarchical integrators.
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

module HamiltonianSplitting
  
  use MathConstant
  use SymplecticIntegrator
  use IsochronePropagator
  use Isochrone
  implicit none

  public:: Hamiltonian, integstep_symplectic

contains
!***********************************************************************************************

  ! Value of the Hamiltonian function.
  ! X, t: state vector and time
  ! 
  function Hamiltonian(X, t, pot)
    real(prc):: Hamiltonian
    real(prc),intent(in):: X(6), t
    interface
     function pot(XX,tt)
       use MathConstant, only: prc
       real(prc):: pot
       real(prc),intent(in):: XX(6), tt
     end function pot
    end interface
     
    ! kinetic part
    Hamiltonian = 0.5_prc*dot_product(X(4:6),X(4:6))

    ! potential part
    Hamiltonian = Hamiltonian + pot(X, t)

  end function Hamiltonian

!***********************************************************************************************

  ! Performs one time step on (X,t) using:
  ! 
  !   - The Leapfrog scheme
  !   - The time step Dt
  !   - If mu = 0:
  !       A = kinetic energy and B = potential ;
  !     else if mu /= 0 && b = 0:
  !       A = Kepler Hamiltonian and B = potential ;
  !     else:
  !       A = isochrone Hamiltonian and B = potential;
  !
  !     where A and B refer to what is written above.
  ! 
  subroutine integstep_symplectic(X, t, Dt, mu, b, Dpot)
    real(prc),intent(inOut):: X(6)
    real(prc),intent(inOut):: t
    real(prc),intent(in):: Dt, mu, b
    interface
     function Dpot(XX,tt)
       use MathConstant, only: prc
       real(prc):: Dpot(6)
       real(prc),intent(in):: XX(6), tt
     end function Dpot
    end interface

    call LeapFrog(phiA, phiB, X, t, Dt)

  contains
    !----------------------------------------------------------
    ! Subroutine performing one step hh with Hamiltonian A
    ! 
    subroutine phiA(XX, tt, hh)
      real(prc),intent(inOut):: XX(:), tt
      real(prc),intent(in):: hh

      if(mu == 0.0_prc) then
        ! Performs a time step with only the kinetic part of the Hamiltonian function
        ! (i.e. using a void potential).     
        XX(1:3) = XX(1:3) + XX(4:6)*hh ! linear drift with constant velocity

      else
        ! isochrone drift
        call step_isochrone(mu, b, XX(1:3), XX(4:6), hh)
      end if

      tt = tt + hh ! updating the time variable

    end subroutine phiA
    !----------------------------------------------------------
    ! Performing one time step hh with Hamiltonian B
    ! 
    subroutine phiB(XX, tt, hh)
      real(prc),intent(inOut):: XX(:), tt
      real(prc),intent(in):: hh

      real(prc):: Dphi(6)

      ! velocity kick
      Dphi(:) = Dpot(XX, tt) - DPotentialIsochrone(XX, tt, mu, b)
      XX(4:6) = XX(4:6) - Dphi(1:3)*hh ! evolution given by Hamilton's equations
                                       ! XX(1:3) remains unchanged.

    end subroutine phiB
    !----------------------------------------------------------
  end subroutine integstep_symplectic

!***********************************************************************************************
end module HamiltonianSplitting
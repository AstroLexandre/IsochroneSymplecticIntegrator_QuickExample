!==================================================================
!  File name : SymplecticIntegrator.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     Defines the Leapfrog integration scheme.
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

module SymplecticIntegrator
  
  use MathConstant
  implicit none

  public:: LeapFrog

contains
!***********************************************************************************************

  ! The Leapfrog scheme is implemented, but you are encouraged to implement
  ! higher-order hierarchical integrators (if needed).
  ! See e.g.: CeMDA, Vol. 80, Issue 1, pp. 39-62 (2001)
  !           (SABA1 & SBAB1 corresponds to the Leapfrog scheme in the paper)
  !
  subroutine LeapFrog(phiA, phiB, X, t, Dt)
    interface
     subroutine phiA(XX,tt,hh)
       use MathConstant, only: prc
       real(prc),intent(inOut):: XX(:), tt
       real(prc),intent(in):: hh
     end subroutine phiA
     subroutine phiB(XX,tt,hh)
       use MathConstant, only: prc
       real(prc),intent(inOut):: XX(:), tt
       real(prc),intent(in):: hh
     end subroutine phiB
    end interface
    real(prc),intent(inOut):: X(:), t
    real(prc),intent(in):: Dt

    call phiA( X, t, 0.5_prc*Dt )
    call phiB( X, t,         Dt )
    call phiA( X, t, 0.5_prc*Dt )

  end subroutine LeapFrog

!***********************************************************************************************
end module SymplecticIntegrator
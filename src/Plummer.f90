!==================================================================
!  File name : Plummer.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     Computes the Plummer potential and its gradient.
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

module Plummer
  
  use MathConstant
  implicit none

  public:: PotentialPlummer, DPotentialPlummer

contains
!***********************************************************************************************

  ! See Monthly Notices of the Royal Astronomical Society, Vol. 71, p.460-470
  !
  function PotentialPlummer(X,t,mu,b)
    real(prc):: PotentialPlummer
    real(prc),intent(in):: X(6), t
    real(prc),intent(in):: mu, b

    PotentialPlummer = -mu/sqrt(dot_product(X(1:3), X(1:3))+b**2)

  end function PotentialPlummer

!***********************************************************************************************

  ! Gradient of the Plummer potential
  ! Vector Out : dPhi_dx, dPhi_dy, dPhi_dz, 0, 0, 0
  !
  function DPotentialPlummer(X,t,mu,b)
    real(prc):: DPotentialPlummer(6)
    real(prc),intent(in):: X(6), t
    real(prc),intent(in):: mu, b

    DPotentialPlummer(:) = 0.0_prc
    DPotentialPlummer(1:3) = mu/sqrt(dot_product(X(1:3), X(1:3))+b**2)**3*X(1:3)

  end function DPotentialPlummer
!***********************************************************************************************
end module Plummer
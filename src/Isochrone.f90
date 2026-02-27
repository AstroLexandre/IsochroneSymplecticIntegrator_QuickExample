!==================================================================
!  File name : Isochrone.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     Computes the isochrone potential and its gradient.
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

module Isochrone
  
  use MathConstant
  implicit none

  public:: PotentialIsochrone, DPotentialIsochrone

contains
!***********************************************************************************************

  ! See Annales d'Astrophysique, Vol. 22, p.126
  !
  function PotentialIsochrone(X,t,mu,b)
    real(prc):: PotentialIsochrone
    real(prc),intent(in):: X(6), t
    real(prc),intent(in):: mu, b

    PotentialIsochrone = -mu/(b+sqrt(b**2+dot_product(X(1:3), X(1:3))))

  end function PotentialIsochrone

!***********************************************************************************************

  ! Gradient of the isochrone potential
  ! Vector Out : dPhi_dx, dPhi_dy, dPhi_dz, 0, 0, 0
  !
  function DPotentialIsochrone(X,t,mu,b)
    real(prc):: DPotentialIsochrone(6)
    real(prc),intent(in):: X(6), t
    real(prc),intent(in):: mu, b
    
    real(prc):: rb

    DPotentialIsochrone(:) = 0.0_prc

    rb = sqrt(b**2+dot_product(X(1:3), X(1:3)))
    DPotentialIsochrone(1:3) = mu*X(1:3)/(b+rb)**2/rb

  end function DPotentialIsochrone
!***********************************************************************************************
end module Isochrone
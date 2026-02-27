!==================================================================
!  File name : MathConstant.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     This module defines various mathematical elements used in 
!     the program.
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

module MathConstant
  implicit none

  integer,parameter,public:: prc = selected_real_kind(15) ! Double (kind=8) : 15 ; Long Double (kind=10) : 16 (not available on all systems) ; Quad (kind=16) : 20
  real(prc),parameter,public:: PI = acos(-1.0_prc)
  real(prc),parameter,public:: twoPI = 2.0_prc*PI

  public:: sqrt3, cross_product

contains
!***********************************************************************************************

  ! Cubic root without error with the argument is negative
  ! 
  function sqrt3(x)
    real(prc):: sqrt3
    real(prc),intent(in):: x
    
    real(prc),parameter:: pw = 1.0_prc/3.0_prc

    sqrt3 = sign(abs(x)**pw, x)
    
  end function sqrt3

!***********************************************************************************************

  ! cross product of two 3-dimensional vectors
  ! 
  function cross_product(U,V)
    real(prc):: cross_product(3)
    real(prc),intent(in):: U(3), V(3)

    cross_product(1) = U(2)*V(3) - U(3)*V(2)
    cross_product(2) = U(3)*V(1) - U(1)*V(3)
    cross_product(3) = U(1)*V(2) - U(2)*V(1)

  end function cross_product

!***********************************************************************************************
end module MathConstant
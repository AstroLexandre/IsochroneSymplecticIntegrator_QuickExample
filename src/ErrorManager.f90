!==================================================================
!  File name : ErrorManager.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     Handles errors in the program.
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

module ErrorManager
  
  implicit none

  public:: EndProgPot

contains
!***********************************************************************************************

  ! Subroutine used to stop the program 
  !
  subroutine EndProgPot(phi)
    character(*),intent(in):: phi

    write(0,*) phi
    write(0,*) "This is a programming error and should not happen!"
    write(0,*) "Program stopped."
    STOP

  end subroutine EndProgPot

!***********************************************************************************************
end module ErrorManager
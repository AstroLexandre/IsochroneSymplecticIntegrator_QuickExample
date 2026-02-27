!==================================================================
!  File name : KeplerEquation.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     This program numerically solves Kepler's equation.
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

module KeplerEquation

  use MathConstant
  implicit none

  public:: invKep_par, invKep_hyp, invKepvar

contains
!***********************************************************************************************

  ! Inversion of Kepler equation in the elliptic case in variational form (Danby, 1988).
  ! 
  ! The input parameters are:
  !    DM = M - M0 = mean anomaly variation
  !    k = e*cos(u0)
  !    h = e*sin(u0)
  ! where u0 is the eccentric anomaly at initial point.
  ! 
  ! The function returns the variation of eccentric anomaly Du = u - u0.
  ! 
  ! /!\ The value of the eccentricity must be in [0:1] (but no check is realised to speed up computations).
  ! 
  function invKepvar(DM, k, h)
     real(prc):: invKepvar
     real(prc),intent(in):: DM, k, h

     ! NB: Due to round-off errors, the method may oscillate up and down around the root
     ! for some isolated parameter values. We set the 5*epsilon limit
     ! in order to avoid digging into the round-off noise too often, but
     ! the vast majority of points converge below this value in the same number
     ! of iterations (in double precision: mostly 2, somewhere 1 and 3).
     real(prc),parameter:: eps = 5.0_prc*epsilon(1.0_prc)

     ! maximum number of iterations (set to the number of digits of the real kind)
     integer,parameter:: IMAX = precision(1.0_prc)

     real(prc):: e2, h2, f, df, ddfS2
     real(prc):: DMs, Du, dDu, ksin, hcos, hsin, kcos
     real(prc):: Du_old1, Du_old2
     integer:: i

     DMs = modulo(DM,twoPI) ! DM in [-2*pi,0] or [0,2pi]
     if(abs(DMs) > PI) DMs = DMs - sign(twoPI,DMs) ! DM in [-pi,pi]

     !-------------------------------------------
     ! initialisations

     h2 = h**2
     e2 = h2 + k**2 ! should be <= 1

     if((h2 < 0.09_prc*e2) .and. (abs(DMs) < 0.2_prc)) then
        Du = DMs + (sqrt3(6.0_prc*DMs)-DMs)*e2
     else
        Du = DMs - h + sign(0.85_prc, h*cos(DMs-h)+k*sin(DMs-h) )*sqrt(e2)
     endIf

     Du_old1 = Du ! previous 1 value of Du
     Du_old2 = Du ! previous 2 value of Du

     dDu = MAX(abs(Du),1.0_prc) ! needed to have dDu/Du very different from zero at initialisation

     !-------------------------------------------
     ! iterations

     do i = 0,IMAX

        ksin = k*sin(Du)
        hcos = h*cos(Du)
        f = Du - ksin - hcos + h - DMs

        if(abs(f) <  eps .or. abs(dDu) < epsilon(1.0_prc)*MAX(abs(Du),1.0_prc)) then ! convergence reached
           invKepvar = Du
           !write(*,'(ES22.15,2(X,ES22.15),X,I3)') atan2(h,k), sqrt(e2), DM, i ! for convergence tests
           return
        endIf

        kcos = k*cos(Du)
        hsin = h*sin(Du)
        df = 1.0_prc - kcos + hsin

        dDu = -f/df ! Newton method
        ddfS2 = 0.5_prc*(ksin + hcos)
        dDu = -f/(df + dDu*ddfS2) ! Halley method
        dDu = -f/(df + dDu*ddfS2 + dDu**2*(kcos-hsin)/6.0_prc) ! higher order

        Du = Du + dDu

        ! in case Du oscillates up and down around the root
        if(Du == Du_old2) then
           dDu = 0.5_prc*dDu
           Du = Du - dDu ! we try to break the infinite loop by taking the central value
        endIf

        Du_old2 = Du_old1 ! updating the previous values of Du
        Du_old1 = Du

     endDo

     invKepvar = Du
     !write(*,'(ES22.15,2(X,ES22.15),X,I3)') atan2(h,k), sqrt(e2), DM, i ! for convergence tests

     write(0,'(A)') " WARNING: maximum number of iterations reached in invKepvar."
     write(0,'(A,ES22.15,2(X,ES22.15))') " Input arguments are: (DM,k,h) = ", DM, k, h
     write(0,'(A,ES22.15)') " Convergence state is: ", abs(f)

  end function invKepvar

!***********************************************************************************************

  ! Inversion of Kepler equation in the hyperbolic case (Serafin, 1998). Returns the eccentric
  ! anomaly.
  ! 
  ! /!\ The eccentricity must be >= 1 (but no test is realised to speed up computations).
  ! 
  function invKep_hyp(M,e)
     real(prc):: invKep_hyp
     real(prc),intent(in):: M, e

     ! precision reachable
     real(prc),parameter:: eps = epsilon(1.0_prc)

     ! maximum number of iterations (ajusted according to the number of digits of the real kind)
     integer,parameter:: IMAX = 2.0_prc*precision(1.0_prc)

     real(prc):: Ms
     real(prc):: esinh, ecosh
     real(prc):: f, df, ddfS2
     real(prc):: u, du
     real(prc):: u_old1, u_old2
     integer:: i

     Ms = abs(M) ! M in [0,+inf]

     !-------------------------------------------
     ! initialisations

     if(e == 1.0_prc) then
        u = sqrt3(6.0_prc*Ms)
     else
        u = MIN( asinh(Ms/(e-1.0_prc)) , sqrt3(6.0_prc*Ms/e) )
     endIf

     u_old1 = u ! previous 1 value of u
     u_old2 = u ! previous 2 value of u

     du = MAX(abs(u),1.0_prc) ! needed to have du/u very different from zero at initialisation

     !-------------------------------------------
     ! iterations

     do i = 0,IMAX

        esinh = e*sinh(u)
        f = esinh - u - Ms

        if(abs(f) < eps .or. abs(du) < epsilon(1.0_prc)*MAX(abs(u),1.0_prc) ) then ! convergence reached
           invKep_hyp = sign(1.0_prc,M)*u
           !write(*,'(2(ES22.15,1X),I0)') M, e, i ! for convergence tests
           return
        endIf

        ecosh = e*cosh(u)
        df = ecosh - 1.0_prc

        du = -f/df ! Newton method
        ddfS2 = esinh*0.5_prc
        du = -f/(df + du*ddfS2) ! Halley method
        du = -f/(df + du*ddfS2 + du**2*ecosh/6.0_prc) ! higher order

        u = u + du

        ! in case u oscillates up and down around the root
        if(u == u_old2) then
           du = 0.5_prc*du
           u = u - du ! we try to break the infinite loop by taking the central value
        endIf

        u_old2 = u_old1 ! updating the previous values of u
        u_old1 = u

     endDo

     invKep_hyp = sign(1.0_prc,M)*u
     !write(*,'(2(ES22.15,1X),I0)') M, e, i ! for convergence tests

     write(0,'(A)') " WARNING: maximum number of iterations reached in invKep_hyp."
     write(0,'(A,ES22.15,X,ES22.15)') " Input arguments are: (M,e) = ", M, e
     write(0,'(A,ES22.15)') " Convergence state is: ", abs(f)

  end function invKep_hyp


!***********************************************************************************************

  ! Inversion of Kepler equation in the exact parabolic case. The solution is in closed form
  ! and no iterations are needed.
  ! 
  function invKep_par(M)
     real(prc):: invKep_par
     real(prc),intent(in):: M

     real(prc):: Q, R

     R = 3.0_prc*M
     Q = sqrt(1.0_prc + R**2)
     invKep_par = sqrt3(R+Q) + sqrt3(R-Q)
     
  end function invKep_par

!***********************************************************************************************
end module KeplerEquation
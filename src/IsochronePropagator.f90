!==================================================================
!  File name : IsochronePropagator.f90
!  Project   : IsochroneSymplecticIntegrator_QuickExample
!  Program   : run_IsochroneSymplecticIntegrator_QuickExample.x
!
!  Description :
!     This program analytically propagates the isochrone
!     trajectory of the particle from t0 to t0+Dt.
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

module IsochronePropagator

  use MathConstant
  use ErrorManager
  use KeplerEquation
  implicit none

  public:: step_isochrone
  private:: StepIsoP

  !----------------------------------------------------
  ! Parameters of the isochrone step.
  ! They are all variables needed to calculate the 
  ! analytical trajectories in isoStep.
  ! 
  type StepIsoP
    real(prc):: r02, v02, rv0, u0, u, r2, rdr_dt, r2dnu_dt, e, ang, ang0, Iang, Iang0, beta2, z, na2, zbeta, muz, coef1_2
    real(prc):: L(3)
    character(3):: Ener
  end type StepIsoP
  !----------------------------------------------------

contains
  !***********************************************************************************************

  ! Two-body propagator. Performs a step Dt in the Cartesian position and velocity
  ! (X,V) in the exact Hénon-isochrone two-body motion with constants mu and beta.
  ! 
  ! All cases are supported!
  ! 
  ! The unit length is supposed to be chosen so that characteristic distances are of order
  ! unity:
  !      10*epsilon can be considered as extremely small as pericentre distance
  !      1/(10*epsilon) can be considered as extremely large as semi-major axis
  ! 
  !
  subroutine step_isochrone(mu, beta, X, V, Dt)
    real(prc), intent(in):: mu, beta
    real(prc), intent(inOut):: X(3), V(3)
    real(prc), intent(in):: Dt

    real(prc), parameter:: eps = 1.0e1_prc*epsilon(1.0_prc)

    type(StepIsoP):: particule

    real(prc):: DM, Du, k0, h0, z2, n, r0beta, S
    real(prc):: r2dnu_dt2, r2dnu_dt_bmumu, tem

    if(abs(mu) <= epsilon(1.0_prc)) call EndProgPot("FATAL ERROR in isostep: mu -> 0 !!!")

    particule%r02 = dot_product(X, X)
    particule%v02 = dot_product(V, V)
    particule%rv0 = dot_product(X, V)

    particule%L(:) = cross_product(X, V)                   ! vecteur moment cinétique
    particule%r2dnu_dt = norm2(particule%L)                ! moment cinétique

    particule%beta2 = beta**2                              ! b^2
    r0beta = sqrt(particule%r02+particule%beta2)           
    particule%z = particule%v02/mu - 2.0_prc/(beta+r0beta) ! 2h/mu
  
    !--------------------- rosette case ---------------------
    if(particule%z < -eps) then ! z = -1/a

      z2 = particule%z**2
      n = sqrt(-mu*z2*particule%z)      ! moyen mouvement
      particule%na2 = n/z2
      k0 = 1.0_prc+r0beta*particule%z   ! e*cos(u0)
      h0 = particule%rv0/particule%na2  ! e*sin(u0)
      particule%e = sqrt(k0**2+h0**2)

      if(particule%e < epsilon(1.0_prc)) then
        particule%u0 = 0.0_prc
      else
        particule%u0 = atan2(h0, k0)    ! anomalie excentrique initiale
      endif

      DM = n*Dt                                      ! variation anomalie moyenne
      Du = invKepvar(DM, k0, h0) + ceiling_force(DM) ! variation in eccentric anomaly
      particule%u = particule%u0 + Du                ! final eccentric anomaly

      ! We cannot use the formula: r^2 = a^2*(1 - e*cos(u))^2 - b^2
      ! because it generates dramatic cancellation errors when r << b. We use the equivalent formula below:
      ! particule%r2 = (1.0_prc - particule%e*cos(particule%u))**2/z2 - particule%beta2
      S = (2.0_prc*(1.0_prc + particule%z*r0beta)*sin(0.5_prc*Du)**2 + h0*sin(Du))/particule%z
      particule%r2 = particule%r02 + S**2 - 2.0_prc*r0beta*S

      particule%rdr_dt = particule%na2*particule%e*sin(particule%u) ! produit scalaire entre le vecteur position et vitesse

      particule%Ener = "inf"

    !-------------------- hyperbolic-like case --------------
    else if(particule%z > eps) then ! z = 1/a

      z2 = particule%z**2
      n = sqrt(mu*z2*particule%z)
      particule%na2 = n/z2
      
      k0 = 1.0_prc+r0beta*particule%z   ! e*cosh(u0)
      h0 = particule%rv0/particule%na2  ! e*sinh(u0)

      particule%zbeta = beta*particule%z
      particule%muz = mu/particule%z

      particule%e = sqrt((1.0_prc+particule%zbeta)**2+particule%r2dnu_dt**2/particule%muz)
      particule%u0 = asinh(h0/particule%e)

      particule%u = invKep_hyp(h0-particule%u0+n*Dt, particule%e) ! final eccentric anomaly

      particule%r2 = (particule%e*cosh(particule%u)-1.0_prc)**2/z2-particule%beta2
      particule%rdr_dt = particule%na2*particule%e*sinh(particule%u) ! produit scalaire entre le vecteur position et vitesse

      particule%Ener = "sup"

    !--------------------- parabolic-like case --------------
    else

      r2dnu_dt2 = particule%r2dnu_dt**2
      r2dnu_dt_bmumu = r2dnu_dt2+2.0_prc*mu*beta

      tem = particule%rv0*(particule%rv0**2+3.0_prc*r2dnu_dt_bmumu)/(6.0_prc*mu**2)+Dt

      particule%rdr_dt = sqrt3(3.0_prc*mu**2*tem+sqrt(9.0_prc*mu**4*tem**2+r2dnu_dt_bmumu**3))
      particule%rdr_dt = particule%rdr_dt-r2dnu_dt_bmumu/particule%rdr_dt ! produit scalaire entre le vecteur position et vitesse

      particule%r2 = ((r2dnu_dt2+particule%rdr_dt**2)/(2.0_prc*mu)+beta)**2-particule%beta2

      particule%coef1_2 = sqrt(r2dnu_dt2+4.0_prc*mu*beta)

      particule%Ener = "egl"

    endIf

    call calMotion(X(:), V(:), particule)

  !---------------------------------------------------------------------------------------
  contains
    !---------------------------------------------------------------------------------------
    ! Subroutine mettant à execution le calcul de trajectoire analytique.
    !
    !----------------------------------|
    ! Dans le cas générale             |
    !__________________________________|
    !  Mouvement classique suivant la valeur de l'énergie : cercle, ellipse & rosette pour 
    !  le cas lié ; courbe s'apparentant à une parabole ou une hyperbole pour le cas non-lié.
    !
    !----------------------------------|
    ! Dans le cas moment cinétique nul |
    !__________________________________|
    !  Le mouvement est rectiligne.
    !
    !  /!\ Le vecteur X(:) est normalisé en fonction du contexte considéré et il sert de 
    !      point de référence pour déterminer la direction vers laquelle se propage la 
    !      particule. Le vecteur normalisé est alors appelé vecteur XVn ; le vecteur X(:) 
    !      initial assumera cette fonction afin d'éviter l'allocation de mémoire supplémentaire.
    !
    !  Mouvement lié : Le vecteur XVn indique le mouvement initial de la particule : il est 
    !  --------------- défini à partir du vecteur position ou de la vitesse via une 
    !                  normalisation. Un seuil très large est choisi pour déterminer quand 
    !                  utiliser le vecteur X(:) ou V(:) dans la définition de XVn afin 
    !                  d'éviter la division par zéro.
    !
    !  Mouvement non-lié : Le vecteur XVn est défini à partir du vecteur vitesse normalisé.
    !  -------------------                                           
    !
    subroutine calMotion(X, V, particule)
      real(prc), intent(inOut):: X(3), V(3)
      type(StepIsoP), intent(inOut):: particule

      real(prc), parameter:: PI_two = PI/2.0_prc

      real(prc):: Sx, Sv

      if(particule%r2dnu_dt > epsilon(1.0_prc)) then ! On vérifie si l'orbite est bâton ou non

        select case(particule%Ener)
          case('inf')
            particule%zbeta = beta*particule%z
            particule%muz = sqrt(-mu/particule%z)
            particule%coef1_2 = sqrt(particule%r2dnu_dt**2+4.0_prc*beta*mu)

            particule%ang = tan(particule%u/2.0_prc)
            particule%ang0 = tan(particule%u0/2.0_prc)
            particule%Iang = ceiling_floor(particule%u)
            particule%Iang0 = ceiling_floor(particule%u0)

          case('sup')
            particule%muz = sqrt(particule%muz)
            particule%coef1_2 = sqrt(particule%r2dnu_dt**2+4.0_prc*beta*mu)

            particule%ang = tanh(particule%u/2.0_prc)
            particule%ang0 = tanh(particule%u0/2.0_prc)
            particule%Iang = 0.0_prc
            particule%Iang0 = 0.0_prc

          case('egl')
            particule%zbeta = 0.0_prc
            particule%muz = 1.0_prc
            particule%e = 0.0_prc ! en réalité, e = 1, mais l'architecture du propagateur analytique impose e = 0

            particule%ang = particule%rdr_dt
            particule%ang0 = particule%rv0
            particule%Iang = 0.0_prc
            particule%Iang0 = 0.0_prc
          
          case default
            call EndProgPot("Programming Error in isostep: cannot set the classical motion!")
        end select

        particule%L(:) = particule%L(:)/particule%r2dnu_dt ! normalize angular momentum
        call classic_motion(X(:), V(:), particule)

      else

        if(particule%r02 == 0.0_prc .and. particule%v02 == 0.0_prc) return ! la particule ne bouge pas

        select case(particule%Ener)
          case('inf')
            call rectilinear_osc(Sx, Sv, modulo(particule%u-particule%rdr_dt/particule%na2, 4.0_prc*PI))

            if(particule%u0 < -PI_two) then
              X(:) = -X(:)/sqrt(particule%r02)
            else if(abs(particule%u0) <= PI_two) then
              X(:) = V(:)/sqrt(particule%v02)
            else
              X(:) = X(:)/sqrt(particule%r02)
            endif

          case('sup')
            call rectilinear_exp(Sx, Sv, particule%e*sinh(particule%u)-particule%u)
            X(:) = V(:)/sqrt(particule%v02)

          case('egl')
            call rectilinear_exp(Sx, Sv, particule%rdr_dt)
            X(:) = V(:)/sqrt(particule%v02)
          
          case default
            call EndProgPot("Programming Error in isostep: cannot set the rectilinear motion!")
        end select

        call rectilinear_motion(X(:), V(:), particule, Sx, Sv)

      endif

    end subroutine calMotion
    !
    !---------------------------------------------------------------------------------------
    ! 
    subroutine classic_motion(X, V, particule)
      real(prc), intent(inOut):: X(3), V(3)
      type(StepIsoP), intent(in):: particule

      real(prc):: Xe(2), LX(3), exc(3)
      real(prc):: XtoV(2, 2), matr(3, 2)
      real(prc):: nu, nu0

      call polar_routine(nu, nu0, particule)

      XtoV(1, :) = (/particule%rdr_dt, -particule%r2dnu_dt/)
      XtoV(2, :) = (/particule%r2dnu_dt, particule%rdr_dt/)

      LX(:) = cross_product(particule%L, X)
      exc(:) = (cos(nu0)*cross_product(LX, particule%L)-sin(nu0)*LX(:))/sqrt(particule%r02) ! Vecteur excentricité au précédent péricentre
      
                                                         ! exc                       L x exc                   |  L
      matr(:, 1) = exc(:)                                ! cosw*cosO-sinw*sinO*cosI, -sinw*cosO-cosw*sinO*cosI |, sinO*sinI
      matr(:, 2) = cross_product(particule%L, exc)       ! cosw*sinO+sinw*cosO*cosI, -sinw*sinO+cosw*cosO*cosI |, -cosO*sinI
                                                         ! -sinw*sinI              , -cosw*sinI                |, cosI

      Xe(:) = sqrt(particule%r2)*(/cos(nu), sin(nu)/)

      X(:) = matmul(matr, Xe)
      V(:) = matmul(matr, matmul(XtoV, Xe))/particule%r2 ! Ve(:) = matmul(XtoV, Xe)/r2

    end subroutine classic_motion
    !
    !---------------------------------------------------------------------------------------
    ! 
    subroutine rectilinear_motion(X, V, particule, Sx, Sv)
      real(prc), intent(inOut):: X(3), V(3)
      real(prc), intent(in):: Sx, Sv
      type(StepIsoP), intent(in):: particule

      V(:) = sign(sqrt(mu*(particule%z+2.0_prc/(beta+sqrt(particule%r2+particule%beta2)))), Sv)*X(:)
      X(:) = sign(sqrt(particule%r2), Sx)*X(:)

    end subroutine rectilinear_motion

    !------------------------------------------------|
    ! Les deux subroutines suivantes cherchent la    |
    ! direction vers laquelle la particule se dirige |
    ! à l'instant t+Dt dans le cas d'une chute libre.|
    !________________________________________________|
    !---------------------------------------------------------------------------------------
    ! Chute libre de la particule dans le cas lié : la particule oscille sur une orbite 
    ! rectiligne en passant par l'origine.
    !
    subroutine rectilinear_osc(Sx, Sv, modM)
      real(prc), intent(Out):: Sx, Sv
      real(prc), intent(in):: modM

      if(modM <= PI) then
        Sx = 1.0_prc
        Sv = 1.0_prc
      else if(modM > PI .and. modM <= twoPI) then
        Sx = 1.0_prc
        Sv = -1.0_prc
      else if(modM > twoPI .and. modM <= 3.0_prc*PI) then
        Sx = -1.0_prc
        Sv = -1.0_prc
      else ! if(modM > threePI .and. modM <= fourPI) then
        Sx = -1.0_prc
        Sv = 1.0_prc
      endif

    end subroutine rectilinear_osc
    !
    !---------------------------------------------------------------------------------------
    ! Chute libre de la particule dans le cas non-lié : une fois passé le péricentre, 
    ! la particule ne reviendra plus jamais.
    !
    subroutine rectilinear_exp(Sx, Sv, modM)
      real(prc), intent(Out):: Sx, Sv
      real(prc), intent(in):: modM

      if(modM <= 0.0_prc) then
        Sx = -1.0_prc
        Sv = 1.0_prc
      else
        Sx = 1.0_prc
        Sv = 1.0_prc
      endif

    end subroutine rectilinear_exp

    !---------------------------------------------------------------------------------------
    ! On calcule l'angle polaire de la particule au nouvel instant de temps. 
    ! La formule originale a été modifiée afin d'éviter les erreurs d'arrondi, car une 
    ! certaine quantité pouvait tendre vers zéro. En conséquence, la division par cette 
    ! quantité était indisponible dans le calcul de l'angle polaire. 
    ! En effet, cette erreur apparaît lorsqu'on est à la limite de l'orbite bâton, càd : 
    !       * beta/a = 1-e (dans le cas de la rosette) 
    !       * beta/a = e-1 (dans le cas hyperbolique).
    ! De plus, cette reformulation à l'avantage d'être commune aux trois types de trajectoires.
    !
    subroutine polar_routine(nu, nu0, particule)
      real(prc), intent(Out):: nu, nu0
      type(StepIsoP), intent(in):: particule

      real(prc):: coef1, coef2, coef3

      coef1 = particule%r2dnu_dt/particule%coef1_2
      coef2 = particule%muz*(particule%e+1.0_prc-particule%zbeta)/particule%coef1_2
      coef3 = particule%muz*(particule%e+1.0_prc+particule%zbeta)/particule%r2dnu_dt

      nu = nu_angle(particule%ang, particule%Iang, coef1, coef2, coef3)
      nu0 = nu_angle(particule%ang0, particule%Iang0, coef1, coef2, coef3)

    end subroutine polar_routine
    !
    !---------------------------------------------------------------------------------------
    ! Calcul de l'anomalie vraie à un instant t.
    !
    function nu_angle(angle, int_angle, coef1, coef2, coef3)
      real(prc):: nu_angle
      real(prc), intent(in):: angle, int_angle, coef1, coef2, coef3

      nu_angle = modulo(coef1*(atan(coef2*angle)+int_angle)+atan(coef3*angle)+int_angle, twoPI)

    end function nu_angle
    !
    !---------------------------------------------------------------------------------------
    ! On vérifie le signe de l'anomalie excentrique, puis on applique la fonction appropriée 
    ! pour récupérer la partie entière de l'angle (par défaut ou par excès). Elle est ensuite
    ! multiplié par pi. Ce mécanisme permet de rendre les fonctions arctan continues dans la 
    ! fonction nu_angle().
    !
    function ceiling_floor(angle)
      real(prc):: ceiling_floor
      real(prc), intent(in):: angle

      if(angle >= 0.0_prc) then
        ceiling_floor = ceiling((angle-PI)/twoPI)*PI
      else
        ceiling_floor = floor((angle+PI)/twoPI)*PI
      endif

    end function ceiling_floor

    !---------------------------------------------------------------------------------------
    ! On récupère l'information sur le nombre de tours que la particule a effectués après 
    ! un pas de temps Dt (pour le cas lié).
    !
    ! REMARQUE : on évite l'erreur numérique avec ce test.
    !
    function ceiling_force(angle)
      real(prc):: ceiling_force
      real(prc), intent(in):: angle

      if(angle >= 0.0_prc) then
        ceiling_force = ceiling((angle-PI)/twoPI)*twoPI
      else
        ceiling_force = (ceiling((angle+PI)/twoPI)-1.0_prc)*twoPI
      endif

    end function ceiling_force

  end subroutine step_isochrone

!***********************************************************************************************
end module IsochronePropagator
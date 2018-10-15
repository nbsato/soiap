module ewald_module
  implicit none
  real(8), parameter :: M_PI       = 3.141592653589793238462d0 !< Pi, circular constant
  real(8), parameter :: M_2_SQRTPI = 1.12837916709551257390d0 !< 2/sqrt(Pi)

contains

  subroutine calcPotentialEwald( P, R, uv, bv, Omega )
    real(8), intent(out) :: P
    real(8), intent(in) :: R(3)
    real(8), intent(in) :: uv(3,3)
    real(8), intent(in) :: bv(3,3)
    real(8), intent(in) :: Omega

    real(8) :: alpha
    integer :: Nr, Ng
    integer :: iGa, iGb, iGc
    integer :: ira, irb, irc
    real(8) :: G(3), GG, GR, RL(3), dl
    real(8) :: derfc ! prototype

    alpha = sqrt(0.84d0*M_PI/sum(uv(:,1)*uv(:,1)))
    Ng = 3
    Nr = 3

    P=0.0d0

    do iGa=-Ng, Ng
       do iGb=-Ng, +Ng
          do iGc=-Ng, +Ng
             G(:) = bv(:,1) * dble(iGa) &
                  + bv(:,2) * dble(iGb) &
                  + bv(:,3) * dble(iGc)
             GG = sum(G(:)*G(:))
             GR = sum(G(:)*R(:))

             if( GG == 0.0d0 ) then
                P = P - M_PI/(alpha*alpha)
             else
                P = P + cos(GR) &
                     * 4.0d0*M_PI/GG * exp(-GG/(4.0d0*alpha*alpha))
             end if
          end do
       end do
    end do
    P = P * (1.d0/Omega)

    do ira=-Nr, +Nr
       do irb=-Nr,+Nr
          do irc=-Nr, +Nr
             RL(:) = R(:) &
                  + uv(:,1) * dble(ira) &
                  + uv(:,2) * dble(irb) &
                  + uv(:,3) * dble(irc)
             dl = sqrt(sum(RL(:)*RL(:)))

             if( dl == 0.0d0 ) then
                P = P - alpha*M_2_SQRTPI
             else
                P = P + derfc(alpha*dl)/dl
             end if
          end do
       end do
    end do

  end subroutine calcPotentialEwald

  subroutine calcForceEwald( F, R, uv, bv, Omega )
    real(8), intent(out) :: F(3)
    real(8), intent(in) :: R(3)
    real(8), intent(in) :: uv(3,3)
    real(8), intent(in) :: bv(3,3)
    real(8), intent(in) :: Omega

    real(8) :: alpha
    integer :: Nr, Ng
    integer :: iGa, iGb, iGc
    integer :: ira, irb, irc
    real(8) :: G(3), GG, GR, df, RL(3), dl
    real(8) :: derfc ! prototype

    alpha = sqrt(0.84d0*M_PI/sum(uv(:,1)*uv(:,1)))
    Ng = 3
    Nr = 3

    F(:) = [ 0.0d0, 0.0d0, 0.0d0 ]
    do iGa=-Ng, Ng
       do iGb=-Ng, +Ng
          do iGc=-Ng, +Ng
             G(:) = bv(:,1) * dble(iGa) &
                  + bv(:,2) * dble(iGb) &
                  + bv(:,3) * dble(iGc)
             GG = sum(G(:)*G(:))
             GR = sum(G(:)*R(:))

             if( GG == 0.d0 ) then
                df = 0.0d0
             else
                df = (-1.0d0)*sin(GR) &
                     * 4.0d0*M_PI/GG*exp(-GG/(4.0d0*alpha*alpha))
             end if
             F(:) = F(:) + df*G(:)
          end do
       end do
    end do
    F(:) = F(:) * (1.d0/Omega)

    do ira=-Nr, +Nr
       do irb=-Nr,+Nr
          do irc=-Nr, +Nr
             RL(:) = R(:) &
                  + uv(:,1) * dble(ira) &
                  + uv(:,2) * dble(irb) &
                  + uv(:,3) * dble(irc)
             dl = sqrt(sum(RL(:)*RL(:)))

             if( dl == 0.0d0 ) then
             else
                df = ( - (alpha*dl)*gauss(alpha*dl) - derfc(alpha*dl) )/dl**2
                F(:) = F(:) + df*RL(:)/dl
             end if
          end do
       end do
    end do
  end subroutine calcForceEwald

  subroutine calcStressEwald( S, R, uv, bv, Omega )
    real(8), intent(out) :: S(3,3)
    real(8), intent(in) :: R(3)
    real(8), intent(in) :: uv(3,3)
    real(8), intent(in) :: bv(3,3)
    real(8), intent(in) :: Omega

    real(8) :: alpha
    integer :: Nr, Ng
    integer :: iGa, iGb, iGc
    integer :: ira, irb, irc
    real(8) :: G(3), GG, GR, RL(3), dl
    real(8) :: derfc ! prototype
    real(8), parameter :: delta(3,3) &
         =  reshape( (/ &
         1.0d0, 0.0d0, 0.0d0, &
         0.0d0, 1.0d0, 0.0d0, &
         0.0d0, 0.0d0, 1.0d0 /), (/3,3/) )
    integer :: la, mu

    alpha = sqrt(0.84d0*M_PI/sum(uv(:,1)*uv(:,1)))
    Ng = 3
    Nr = 3

    S(:,:)=0.0d0

    do iGa=-Ng, Ng
       do iGb=-Ng, +Ng
          do iGc=-Ng, +Ng
             G(:) = bv(:,1) * dble(iGa) &
                  + bv(:,2) * dble(iGb) &
                  + bv(:,3) * dble(iGc)
             GG = sum(G(:)*G(:))
             GR = sum(G(:)*R(:))

             if( GG == 0.0d0 ) then
                S(:,:) = S(:,:) - delta(:,:)*M_PI/(alpha*alpha)
             else
                S(:,:) = S(:,:) &
                     + delta(:,:)*cos(GR) &
                     * 4.0d0*M_PI/GG * exp(-GG/(4.0d0*alpha*alpha))

                do la=1, 3
                   do mu=1, 3
                      S(la,mu) = S(la,mu) &
                           - cos(GR) &
                           * 4.0d0*M_PI/GG * exp(-GG/(4.0d0*alpha*alpha)) &
                           * 2.0d0*G(la)*G(mu)/GG &
                           * (1.0d0+GG/(4.0d0*alpha*alpha))
                   end do
                end do
             end if
          end do
       end do
    end do
    S = S * (1.d0/Omega)

    do ira=-Nr, +Nr
       do irb=-Nr,+Nr
          do irc=-Nr, +Nr
             RL(:) = R(:) &
                  + uv(:,1) * dble(ira) &
                  + uv(:,2) * dble(irb) &
                  + uv(:,3) * dble(irc)
             dl = sqrt(sum(RL(:)*RL(:)))

             if( dl == 0.0d0 ) then
             else
                do la=1, 3
                   do mu=1, 3
                      S(la,mu) = S(la,mu) &
                           + alpha*RL(la)*RL(mu)/dl**2 &
                           * ( gauss(alpha*dl) + derfc(alpha*dl)/(alpha*dl) )
                   end do
                end do
             end if
          end do
       end do
    end do

  end subroutine calcStressEwald

  real(8) function gauss( x )
    real(8), intent(in) :: x

    gauss = M_2_SQRTPI*exp(-x*x)
  end function gauss

end module ewald_module

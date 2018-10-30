! ------------------------------------------------------------------------
! Copyright (C) 2017 Nobuya Sato, Hiori Kino, and Takashi Miyake
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
! ------------------------------------------------------------------------

! Angular-dependent potential (ADP) for the Nd-Fe-B system
! A. Kubo, J. Wang, and Y. Umeno, Modelling Simul. Mater. Sci. Eng. 22, 065014 (2014).
module adp_kwu14_module

  use, non_intrinsic :: adp_base_module, only: adp_base_type, kr

  implicit none

  private

  type, public :: adp_kwu14_parameter_phi_type

    private

    real(kr), public :: c1
    real(kr), public :: eta1
    real(kr), public :: c2
    real(kr), public :: eta2
    real(kr), public :: kphi
    real(kr), public :: thetaphi
    real(kr), public :: rc
    real(kr), public :: h

  end type adp_kwu14_parameter_phi_type

  type, public :: adp_kwu14_parameter_rho_type

    private

    real(kr), public :: a0
    real(kr), public :: b0
    real(kr), public :: c0
    real(kr), public :: r0
    real(kr), public :: y
    real(kr), public :: gamma
    real(kr), public :: rc
    real(kr), public :: h

  end type adp_kwu14_parameter_rho_type

  type, public :: adp_kwu14_parameter_f_type

    private

    real(kr), public :: f0
    real(kr), public :: f2
    real(kr), public :: p1
    real(kr), public :: p2
    real(kr), public :: p3

  end type adp_kwu14_parameter_f_type

  type, public :: adp_kwu14_parameter_u_type

    private

    real(kr), public :: dd1 ! D_1
    real(kr), public :: d1 ! d_1
    real(kr), public :: dd2 ! D_2
    real(kr), public :: d2 ! d_2
    real(kr), public :: ku
    real(kr), public :: thetau
    real(kr), public :: rc
    real(kr), public :: h

  end type adp_kwu14_parameter_u_type

  type, public :: adp_kwu14_parameter_w_type

    private

    real(kr), public :: qq1 ! Q_1
    real(kr), public :: q1 ! q_1
    real(kr), public :: qq2 ! Q_2
    real(kr), public :: q2 ! q_2
    real(kr), public :: kw
    real(kr), public :: thetaw
    real(kr), public :: rc
    real(kr), public :: h

  end type adp_kwu14_parameter_w_type

  type, public :: adp_kwu14_parameter_type

    private

    type(adp_kwu14_parameter_phi_type), public :: parameter_phi_Nd_Nd
    type(adp_kwu14_parameter_phi_type), public :: parameter_phi_Nd_Fe
    type(adp_kwu14_parameter_phi_type), public :: parameter_phi_Nd_B
    type(adp_kwu14_parameter_phi_type), public :: parameter_phi_Fe_Fe
    type(adp_kwu14_parameter_phi_type), public :: parameter_phi_Fe_B
    type(adp_kwu14_parameter_phi_type), public :: parameter_phi_B_B
    type(adp_kwu14_parameter_rho_type), public :: parameter_rho_Nd
    type(adp_kwu14_parameter_rho_type), public :: parameter_rho_Fe
    type(adp_kwu14_parameter_rho_type), public :: parameter_rho_B
    type(adp_kwu14_parameter_f_type), public :: parameter_f_Nd
    type(adp_kwu14_parameter_f_type), public :: parameter_f_Fe
    type(adp_kwu14_parameter_f_type), public :: parameter_f_B
    type(adp_kwu14_parameter_u_type), public :: parameter_u_Nd_Nd
    type(adp_kwu14_parameter_u_type), public :: parameter_u_Nd_Fe
    type(adp_kwu14_parameter_u_type), public :: parameter_u_Nd_B
    type(adp_kwu14_parameter_u_type), public :: parameter_u_Fe_Fe
    type(adp_kwu14_parameter_u_type), public :: parameter_u_Fe_B
    type(adp_kwu14_parameter_u_type), public :: parameter_u_B_B
    type(adp_kwu14_parameter_w_type), public :: parameter_w_Nd_Nd
    type(adp_kwu14_parameter_w_type), public :: parameter_w_Nd_Fe
    type(adp_kwu14_parameter_w_type), public :: parameter_w_Nd_B
    type(adp_kwu14_parameter_w_type), public :: parameter_w_Fe_Fe
    type(adp_kwu14_parameter_w_type), public :: parameter_w_Fe_B
    type(adp_kwu14_parameter_w_type), public :: parameter_w_B_B

  end type adp_kwu14_parameter_type

  type, extends(adp_base_type), public :: adp_kwu14_type

    private

    type(adp_kwu14_parameter_type), public :: parameter

  contains

    private

    procedure, public :: cutoff
    procedure, public :: phi
    procedure, public :: phi_derivative
    procedure, public :: rho
    procedure, public :: rho_derivative
    procedure, public :: f
    procedure, public :: f_derivative
    procedure, public :: u
    procedure, public :: u_derivative
    procedure, public :: w
    procedure, public :: w_derivative

    procedure, nopass, public :: psi
    procedure, nopass, public :: psi_derivative

  end type adp_kwu14_type

  ! Tables 2 and 3
  ! energy in eV, length in Ang
  type(adp_kwu14_parameter_type), parameter, public :: adp_kwu14_parameter = adp_kwu14_parameter_type( &
      & parameter_phi_Nd_Nd=adp_kwu14_parameter_phi_type( &
      &     c1=+1.8141e+2_kr, &
      &     eta1=+2.2097_kr, &
      &     c2=-2.4199_kr, &
      &     eta2=+2.2447_kr, &
      &     kphi=+2.4729_kr, &
      &     thetaphi=+3.7907_kr, &
      &     rc=7._kr, &
      &     h=+2.4990_kr), &
      & parameter_phi_Nd_Fe=adp_kwu14_parameter_phi_type( &
      &     c1=-1.0920_kr, &
      &     eta1=+2.0563e-1_kr, &
      &     c2=+1.0813e+2_kr, &
      &     eta2=+5.2096_kr, &
      &     kphi=+9.9557e-1_kr, &
      &     thetaphi=-1.0117_kr, &
      &     rc=7._kr, &
      &     h=+4.6332_kr), &
      & parameter_phi_Nd_B=adp_kwu14_parameter_phi_type( &
      &     c1=+5.9121e-1_kr, &
      &     eta1=-3.6373e-1_kr, &
      &     c2=-2.8066_kr, &
      &     eta2=-4.6774e-1_kr, &
      &     kphi=+8.2479e-2_kr, &
      &     thetaphi=-1.3876_kr, &
      &     rc=7._kr, &
      &     h=+2.1270_kr), &
      & parameter_phi_Fe_Fe=adp_kwu14_parameter_phi_type( &
      &     c1=-7.1938e-2_kr, &
      &     eta1=+2.9722e-1_kr, &
      &     c2=-1.4808e+1_kr, &
      &     eta2=+2.7382_kr, &
      &     kphi=-2.5918_kr, &
      &     thetaphi=+7.5596e-1_kr, &
      &     rc=7._kr, &
      &     h=+5.7408_kr), &
      & parameter_phi_Fe_B=adp_kwu14_parameter_phi_type( &
      &     c1=-4.1688e-4_kr, &
      &     eta1=-6.7311e-1_kr, &
      &     c2=-2.3465e+1_kr, &
      &     eta2=+5.1131_kr, &
      &     kphi=+2.6352_kr, &
      &     thetaphi=-6.1955_kr, &
      &     rc=7._kr, &
      &     h=+1.4540_kr), &
      & parameter_phi_B_B=adp_kwu14_parameter_phi_type( &
      &     c1=-8.3127e-2_kr, &
      &     eta1=-4.3385e-1_kr, &
      &     c2=+1.0180_kr, &
      &     eta2=-6.5785e-1_kr, &
      &     kphi=-5.7598e-2_kr, &
      &     thetaphi=+1.5221_kr, &
      &     rc=7._kr, &
      &     h=+7.1841e-1_kr), &
      & parameter_rho_Nd=adp_kwu14_parameter_rho_type( &
      &     a0=+1.3674e-4_kr, & ! value in Table 2 is wrong [A. Kubo, private communication.]
      &     b0=-1.0263e-10_kr, &
      &     c0=+1.1448e-1_kr, &
      &     r0=+6.9298e-1_kr, &
      &     y=+2.4459e+1_kr, &
      &     gamma=+6.6080_kr, &
      &     rc=7._kr, &
      &     h=+4.0716_kr), &
      & parameter_rho_Fe=adp_kwu14_parameter_rho_type( &
      &     a0=+8.2167e-3_kr, &
      &     b0=+2.4869e+2_kr, &
      &     c0=-6.6259e-1_kr, &
      &     r0=-9.6300e-1_kr, &
      &     y=+4.9180_kr, &
      &     gamma=+9.6095e-1_kr, &
      &     rc=7._kr, &
      &     h=+6.9137_kr), &
      & parameter_rho_B=adp_kwu14_parameter_rho_type( &
      &     a0=+1.0263_kr, &
      &     b0=+2.2677e+2_kr, &
      &     c0=+1.1457e-2_kr, &
      &     r0=-2.7254e-1_kr, &
      &     y=+3.4945_kr, &
      &     gamma=+3.0052_kr, &
      &     rc=7._kr, &
      &     h=+2.3852e-1_kr), &
      & parameter_f_Nd=adp_kwu14_parameter_f_type( &
      &     f0=-4.4823_kr, &
      &     f2=+1.4308e+1_kr, &
      &     p1=+1.4102e+1_kr, &
      &     p2=+1.0452e+1_kr, &
      &     p3=+2.2000e-7_kr), &
      & parameter_f_Fe=adp_kwu14_parameter_f_type( &
      &     f0=-6.6388_kr, &
      &     f2=+1.7102_kr, &
      &     p1=-1.8407_kr, &
      &     p2=+1.5533_kr, &
      &     p3=+1.2000_kr), &
      & parameter_f_B=adp_kwu14_parameter_f_type( &
      &     f0=-6.2818_kr, &
      &     f2=+1.1276e+1_kr, &
      &     p1=-5.6487_kr, &
      &     p2=-2.2536_kr, &
      &     p3=+4.5776_kr), &
      & parameter_u_Nd_Nd=adp_kwu14_parameter_u_type( &
      &     dd1=+1.5940e-3_kr, &
      &     d1=-1.4728_kr, &
      &     dd2=-8.2390e-5_kr, &
      &     d2=-9.2954_kr, &
      &     ku=+3.1605e-3_kr, &
      &     thetau=+4.7005_kr, &
      &     rc=7._kr, &
      &     h=+5.0456_kr), &
      & parameter_u_Nd_Fe=adp_kwu14_parameter_u_type( &
      &     dd1=-9.0774e+1_kr, &
      &     d1=+1.9393_kr, &
      &     dd2=-1.3551e+1_kr, &
      &     d2=+2.5824_kr, &
      &     ku=+1.3600_kr, &
      &     thetau=+6.2832_kr, &
      &     rc=7._kr, &
      &     h=+5.6382_kr), &
      & parameter_u_Nd_B=adp_kwu14_parameter_u_type( &
      &     dd1=-5.3249e+1_kr, &
      &     d1=+2.0463_kr, &
      &     dd2=+3.2819e-1_kr, &
      &     d2=-3.8249_kr, &
      &     ku=-2.5310e-3_kr, &
      &     thetau=-1.5528_kr, &
      &     rc=7._kr, &
      &     h=+4.9208_kr), &
      & parameter_u_Fe_Fe=adp_kwu14_parameter_u_type( &
      &     dd1=+6.9480e-4_kr, &
      &     d1=-6.7207e-1_kr, &
      &     dd2=-6.2276_kr, &
      &     d2=+4.7280_kr, &
      &     ku=-1.7343_kr, &
      &     thetau=+5.2968_kr, &
      &     rc=7._kr, &
      &     h=+2.0710_kr), &
      & parameter_u_Fe_B=adp_kwu14_parameter_u_type( &
      &     dd1=+2.8149e-1_kr, &
      &     d1=+3.1964e-1_kr, &
      &     dd2=+3.6899e-2_kr, &
      &     d2=-5.6051e-1_kr, &
      &     ku=+5.5220_kr, &
      &     thetau=+1.5756_kr, &
      &     rc=7._kr, &
      &     h=+6.3083_kr), &
      & parameter_u_B_B=adp_kwu14_parameter_u_type( &
      &     dd1=+2.4279e+1_kr, &
      &     d1=+2.1148_kr, &
      &     dd2=-3.6605e-4_kr, &
      &     d2=-5.4051_kr, &
      &     ku=-3.3995e-2_kr, &
      &     thetau=+5.8210_kr, &
      &     rc=7._kr, &
      &     h=+5.4862_kr), &
      & parameter_w_Nd_Nd=adp_kwu14_parameter_w_type( &
      &     qq1=+6.3335e-4_kr, &
      &     q1=-2.0938_kr, &
      &     qq2=+5.9700e-6_kr, &
      &     q2=-1.0649e+1_kr, &
      &     kw=-8.9474e-3_kr, &
      &     thetaw=-1.6691_kr, &
      &     rc=7._kr, &
      &     h=+8.9804_kr), &
      & parameter_w_Nd_Fe=adp_kwu14_parameter_w_type( &
      &     qq1=+5.5943e+1_kr, &
      &     q1=+2.6823_kr, &
      &     qq2=-3.6012_kr, &
      &     q2=+3.7140_kr, &
      &     kw=-1.9276_kr, &
      &     thetaw=+4.4443_kr, &
      &     rc=7._kr, &
      &     h=+2.3101_kr), &
      & parameter_w_Nd_B=adp_kwu14_parameter_w_type( &
      &     qq1=+1.3374e+2_kr, &
      &     q1=+2.5516_kr, &
      &     qq2=+1.4397e-3_kr, &
      &     q2=-9.9980_kr, &
      &     kw=+3.4000e-5_kr, &
      &     thetaw=-1.5710_kr, &
      &     rc=7._kr, &
      &     h=+5.6555_kr), &
      & parameter_w_Fe_Fe=adp_kwu14_parameter_w_type( &
      &     qq1=-1.0149e+2_kr, &
      &     q1=+3.0399_kr, &
      &     qq2=-1.5912e+1_kr, &
      &     q2=+3.3707_kr, &
      &     kw=-1.1232_kr, &
      &     thetaw=+2.6762e-1_kr, &
      &     rc=7._kr, &
      &     h=+7.7167_kr), &
      & parameter_w_Fe_B=adp_kwu14_parameter_w_type( &
      &     qq1=+1.0705e+1_kr, &
      &     q1=+2.6656_kr, &
      &     qq2=+0.0000_kr, &
      &     q2=+0.0000_kr, &
      &     kw=+0.0000_kr, &
      &     thetaw=+0.0000_kr, &
      &     rc=7._kr, &
      &     h=+4.0185_kr), &
      & parameter_w_B_B=adp_kwu14_parameter_w_type( &
      &     qq1=+4.6504e-1_kr, &
      &     q1=+1.2802_kr, &
      &     qq2=+3.6899e+1_kr, &
      &     q2=+7.8769_kr, &
      &     kw=+9.0266e-1_kr, &
      &     thetaw=-6.8788e-1_kr, &
      &     rc=7._kr, &
      &     h=+7.3414_kr))

contains

  pure real(kr) function cutoff(this)

    implicit none

    class(adp_kwu14_type), intent(in) :: this

    cutoff = max( &
        & this%parameter%parameter_phi_Nd_Nd%rc, &
        & this%parameter%parameter_phi_Nd_Fe%rc, &
        & this%parameter%parameter_phi_Nd_B%rc, &
        & this%parameter%parameter_phi_Fe_Fe%rc, &
        & this%parameter%parameter_phi_Fe_B%rc, &
        & this%parameter%parameter_phi_B_B%rc, &
        & this%parameter%parameter_rho_Nd%rc, &
        & this%parameter%parameter_rho_Fe%rc, &
        & this%parameter%parameter_rho_B%rc, &
        & this%parameter%parameter_u_Nd_Nd%rc, &
        & this%parameter%parameter_u_Nd_Fe%rc, &
        & this%parameter%parameter_u_Nd_B%rc, &
        & this%parameter%parameter_u_Fe_Fe%rc, &
        & this%parameter%parameter_u_Fe_B%rc, &
        & this%parameter%parameter_u_B_B%rc, &
        & this%parameter%parameter_w_Nd_Nd%rc, &
        & this%parameter%parameter_w_Nd_Fe%rc, &
        & this%parameter%parameter_w_Nd_B%rc, &
        & this%parameter%parameter_w_Fe_Fe%rc, &
        & this%parameter%parameter_w_Fe_B%rc, &
        & this%parameter%parameter_w_B_B%rc)

  end function cutoff

  ! Eq. (6)
  pure real(kr) function phi(this, z1, z2, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z1
    integer, intent(in) :: z2
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_phi_type) :: param

    if (is_Nd(z1) .and. is_Nd(z2)) then
      param = this%parameter%parameter_phi_Nd_Nd
    else if ((is_Nd(z1) .and. is_Fe(z2)) .or. (is_Fe(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_phi_Nd_Fe
    else if ((is_Nd(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_phi_Nd_B
    else if (is_Fe(z1) .and. is_Fe(z2)) then
      param = this%parameter%parameter_phi_Fe_Fe
    else if ((is_Fe(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Fe(z2))) then
      param = this%parameter%parameter_phi_Fe_B
    else if (is_B(z1) .and. is_B(z2)) then
      param = this%parameter%parameter_phi_B_B
    else
      phi = 0._kr
      return
    end if

    phi = this%psi((r - param%rc) / param%h) &
        & * (param%c1 * exp(-param%eta1 * r) &
        &    + param%c2 * cos(param%kphi * r + param%thetaphi) / r ** param%eta2)

  end function phi

  pure real(kr) function phi_derivative(this, z1, z2, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z1
    integer, intent(in) :: z2
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_phi_type) :: param
    real(kr) :: x
    real(kr) :: e
    real(kr) :: a

    if (is_Nd(z1) .and. is_Nd(z2)) then
      param = this%parameter%parameter_phi_Nd_Nd
    else if ((is_Nd(z1) .and. is_Fe(z2)) .or. (is_Fe(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_phi_Nd_Fe
    else if ((is_Nd(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_phi_Nd_B
    else if (is_Fe(z1) .and. is_Fe(z2)) then
      param = this%parameter%parameter_phi_Fe_Fe
    else if ((is_Fe(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Fe(z2))) then
      param = this%parameter%parameter_phi_Fe_B
    else if (is_B(z1) .and. is_B(z2)) then
      param = this%parameter%parameter_phi_B_B
    else
      phi_derivative = 0._kr
      return
    end if

    x = (r - param%rc) / param%h
    e = exp(-param%eta1 * r)
    a = param%kphi * r + param%thetaphi

    phi_derivative = &
        & this%psi_derivative(x) / param%h * (param%c1 * e + param%c2 * cos(a) / r ** param%eta2) &
        & + this%psi(x) * (param%c1 * (-param%eta1) * e &
        &                  + param%c2 * (-sin(a)) * param%kphi / r ** param%eta2 &
        &                  + param%c2 * cos(a) * (-param%eta2) / r ** (param%eta2 + 1._kr))

  end function phi_derivative

  ! Eq. (7)
  pure real(kr) function rho(this, z, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_rho_type) :: param
    real(kr) :: d
    real(kr) :: e

    if (is_Nd(z)) then
      param = this%parameter%parameter_rho_Nd
    else if (is_Fe(z)) then
      param = this%parameter%parameter_rho_Fe
    else if (is_B(z)) then
      param = this%parameter%parameter_rho_B
    else
      rho = 0._kr
      return
    end if

    if (r > param%r0) then
      d = r - param%r0
    else
      d = 0._kr ! rho(r) is originally not specified for the case r < r0 while it cannot be defined.
    end if
    e = exp(-param%gamma * d)

    rho = this%psi((r - param%rc) / param%h) &
        & * (param%a0 * d ** param%y * e * (1._kr + param%b0 * e) + param%c0)

  end function rho

  pure real(kr) function rho_derivative(this, z, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_rho_type) :: param
    real(kr) :: x
    real(kr) :: d
    real(kr) :: e

    if (is_Nd(z)) then
      param = this%parameter%parameter_rho_Nd
    else if (is_Fe(z)) then
      param = this%parameter%parameter_rho_Fe
    else if (is_B(z)) then
      param = this%parameter%parameter_rho_B
    else
      rho_derivative = 0._kr
      return
    end if

    x = (r - param%rc) / param%h
    if (r > param%r0) then
      d = r - param%r0
    else
      d = 0._kr ! rho(r) is originally not specified for the case r < r0 while it cannot be defined.
    end if
    e = exp(-param%gamma * d)

    rho_derivative = &
        & this%psi_derivative(x) / param%h * (param%a0 * d ** param%y * e * (1._kr + param%b0 * e) + param%c0) &
        & + this%psi(x) * (param%a0 * param%y * d ** (param%y - 1._kr) * e * (1._kr + param%b0 * e) &
        &                  + param%a0 * d ** param%y * (-param%gamma) * e * (1._kr + 2._kr * param%b0 * e))

  end function rho_derivative

  ! Eq. (8)
  pure real(kr) function f(this, z, rhobar)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z
    real(kr), intent(in) :: rhobar

    type(adp_kwu14_parameter_f_type) :: param
    real(kr) :: x

    if (is_Nd(z)) then
      param = this%parameter%parameter_f_Nd
    else if (is_Fe(z)) then
      param = this%parameter%parameter_f_Fe
    else if (is_B(z)) then
      param = this%parameter%parameter_f_B
    else
      f = 0._kr
      return
    end if

    x = rhobar - 1._kr

    f = param%f0 &
        & + 0.5_kr * param%f2 * x ** 2 &
        & + param%p1 * x ** 3 &
        & + param%p2 * x ** 4 &
        & + param%p3 * x ** 5

  end function f

  pure real(kr) function f_derivative(this, z, rhobar)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z
    real(kr), intent(in) :: rhobar

    type(adp_kwu14_parameter_f_type) :: param
    real(kr) :: x

    if (is_Nd(z)) then
      param = this%parameter%parameter_f_Nd
    else if (is_Fe(z)) then
      param = this%parameter%parameter_f_Fe
    else if (is_B(z)) then
      param = this%parameter%parameter_f_B
    else
      f_derivative = 0._kr
      return
    end if

    x = rhobar - 1._kr

    f_derivative = &
        & param%f2 * x &
        & + 3._kr * param%p1 * x ** 2 &
        & + 4._kr * param%p2 * x ** 3 &
        & + 5._kr * param%p3 * x ** 4

  end function f_derivative

  ! Eq. (9)
  pure real(kr) function u(this, z1, z2, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z1
    integer, intent(in) :: z2
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_u_type) :: param

    if (is_Nd(z1) .and. is_Nd(z2)) then
      param = this%parameter%parameter_u_Nd_Nd
    else if ((is_Nd(z1) .and. is_Fe(z2)) .or. (is_Fe(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_u_Nd_Fe
    else if ((is_Nd(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_u_Nd_B
    else if (is_Fe(z1) .and. is_Fe(z2)) then
      param = this%parameter%parameter_u_Fe_Fe
    else if ((is_Fe(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Fe(z2))) then
      param = this%parameter%parameter_u_Fe_B
    else if (is_B(z1) .and. is_B(z2)) then
      param = this%parameter%parameter_u_B_B
    else
      u = 0._kr
      return
    end if

    u = this%psi((r - param%rc) / param%h) &
        & * (param%dd1 * exp(-param%d1 * r) &
        &    + param%dd2 * cos(param%ku * r + param%thetau) / r ** param%d2)

  end function u

  pure real(kr) function u_derivative(this, z1, z2, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z1
    integer, intent(in) :: z2
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_u_type) :: param
    real(kr) :: x
    real(kr) :: e
    real(kr) :: a

    if (is_Nd(z1) .and. is_Nd(z2)) then
      param = this%parameter%parameter_u_Nd_Nd
    else if ((is_Nd(z1) .and. is_Fe(z2)) .or. (is_Fe(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_u_Nd_Fe
    else if ((is_Nd(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_u_Nd_B
    else if (is_Fe(z1) .and. is_Fe(z2)) then
      param = this%parameter%parameter_u_Fe_Fe
    else if ((is_Fe(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Fe(z2))) then
      param = this%parameter%parameter_u_Fe_B
    else if (is_B(z1) .and. is_B(z2)) then
      param = this%parameter%parameter_u_B_B
    else
      u_derivative = 0._kr
      return
    end if

    x = (r - param%rc) / param%h
    e = exp(-param%d1 * r)
    a = param%ku * r + param%thetau

    u_derivative = & 
        & this%psi_derivative(x) / param%h * (param%dd1 * e + param%dd2 * cos(a) / r ** param%d2) &
        & + this%psi(x) * (param%dd1 * (-param%d1) * e &
        &                  + param%dd2 * (-sin(a)) * param%ku / r ** param%d2 &
        &                  + param%dd2 * cos(a) * (-param%d2) / r ** (param%d2 + 1._kr))

  end function u_derivative

  ! Eq. (10)
  pure real(kr) function w(this, z1, z2, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z1
    integer, intent(in) :: z2
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_w_type) :: param

    if (is_Nd(z1) .and. is_Nd(z2)) then
      param = this%parameter%parameter_w_Nd_Nd
    else if ((is_Nd(z1) .and. is_Fe(z2)) .or. (is_Fe(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_w_Nd_Fe
    else if ((is_Nd(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_w_Nd_B
    else if (is_Fe(z1) .and. is_Fe(z2)) then
      param = this%parameter%parameter_w_Fe_Fe
    else if ((is_Fe(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Fe(z2))) then
      param = this%parameter%parameter_w_Fe_B
    else if (is_B(z1) .and. is_B(z2)) then
      param = this%parameter%parameter_w_B_B
    else
      w = 0._kr
      return
    end if

    w = this%psi((r - param%rc) / param%h) &
        & * (param%qq1 * exp(-param%q1 * r) &
        &    + param%qq2 * cos(param%kw * r + param%thetaw) / r ** param%q2)

  end function w

  pure real(kr) function w_derivative(this, z1, z2, r)

    implicit none

    class(adp_kwu14_type), intent(in) :: this
    integer, intent(in) :: z1
    integer, intent(in) :: z2
    real(kr), intent(in) :: r

    type(adp_kwu14_parameter_w_type) :: param
    real(kr) :: x
    real(kr) :: e
    real(kr) :: a

    if (is_Nd(z1) .and. is_Nd(z2)) then
      param = this%parameter%parameter_w_Nd_Nd
    else if ((is_Nd(z1) .and. is_Fe(z2)) .or. (is_Fe(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_w_Nd_Fe
    else if ((is_Nd(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Nd(z2))) then
      param = this%parameter%parameter_w_Nd_B
    else if (is_Fe(z1) .and. is_Fe(z2)) then
      param = this%parameter%parameter_w_Fe_Fe
    else if ((is_Fe(z1) .and. is_B(z2)) .or. (is_B(z1) .and. is_Fe(z2))) then
      param = this%parameter%parameter_w_Fe_B
    else if (is_B(z1) .and. is_B(z2)) then
      param = this%parameter%parameter_w_B_B
    else
      w_derivative = 0._kr
      return
    end if

    x = (r - param%rc) / param%h
    e = exp(-param%q1 * r)
    a = param%kw * r + param%thetaw

    w_derivative = &
        & this%psi_derivative(x) / param%h * (param%qq1 * e + param%qq2 * cos(a) / r ** param%q2) &
        & + this%psi(x) * (param%qq1 * (-param%q1) * e &
        &                  + param%qq2 * (-sin(a)) * param%kw / r ** param%q2 &
        &                  + param%qq2 * cos(a) * (-param%q2) / r ** (param%q2 + 1._kr))

  end function w_derivative

  ! Eq. (11)
  pure real(kr) function psi(x)

    implicit none

    real(kr), intent(in) :: x

    if (x < 0) then
      psi = x ** 4 / (1._kr + x ** 4)
    else
      psi = 0._kr
    end if

  end function psi

  pure real(kr) function psi_derivative(x)

    implicit none

    real(kr), intent(in) :: x

    if (x < 0) then
      psi_derivative = 4._kr * x ** 3 / (1._kr + x ** 4) ** 2
    else
      psi_derivative = 0._kr
    end if

  end function psi_derivative

  pure logical function is_Nd(z)

    implicit none

    integer, intent(in) :: z

    integer, parameter :: z_Nd = 60

    is_Nd = z == z_Nd

  end function is_Nd

  pure logical function is_Fe(z)

    implicit none

    integer, intent(in) :: z

    integer, parameter :: z_Fe = 26

    is_Fe = z == z_Fe

  end function is_Fe

  pure logical function is_B(z)

    implicit none

    integer, intent(in) :: z

    integer, parameter :: z_B = 5

    is_B = z == z_B

  end function is_B

end module adp_kwu14_module

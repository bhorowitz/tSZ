double precision function integrand3(lnM,z)
  USE cosmo
  USE sigma
  USE growth
  USE angular_distance
  USE nfw
  USE multipole
! The SZ power spectrum is given by C_l = int dz dV/dz int dlnM dn/dlnM*Tsz^2.
! This routine returns the integrand, dV/dz*dn/dlnM*Tsz^2. Note that Tsz is 
! in units of uK, Rayleigh-Jeans limit. This version uses Arnaud et al.'s 
! pressure profile.
! December 14, 2010: E.Komatsu
  IMPLICIT none
  double precision, intent(IN) :: lnM,z
  double precision :: xin=1d-5,xout=6d0 ! pressure profile is cut at r=xout*r500
  double precision :: omega,tol=1d-6,zin=30d0,pi=3.1415926535d0
  double precision :: mvir,rvir,rs,delc,m500,r500,l500
  double precision :: rhoc,Ez,da,dvdz,transform2d,transform2dp,tsz,tszp,pgnfw,rombint1
  double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh,lnnu,dlnnudlnRh,Rh
  double precision :: g,scalefactor,m180d,m200d
  double precision :: mgas

  real :: chebev,lnsigma2,lnRh,dlnsigma2dlnRh
  integer :: i
  external pgnfw,rombint1
  mvir=dexp(lnM)
! compute omega's, rhoc, and delc
  omega=om0*(1d0+z)**3d0/Ez(z)**2d0 ! Omega(z); E(z)=H(z)/H0
  rhoc=2.775d11*Ez(z)**2d0 ! critical density in units of h^2 M_sun/Mpc^3
  ! Eq.(6) of Bryan & Norman (1998), valid only for a flat universe!
  delc=18d0*pi*pi+82d0*(omega-1d0)-39d0*(omega-1d0)**2d0
  ! One may also use Eq.(C19) of Nakamura & Suto (1997):
  ! delc=omega*18d0*pi*pi*(1d0+0.4093d0*(1d0/omega-1d0)**0.9052d0)
! compute virial and NFW parameters
  rvir=(3d0*mvir/4d0/pi/delc/rhoc)**(1d0/3d0) ! virial radius, h^-1 Mpc
  ! Concentration parameter from Duffy et al. (2008)
  cvir=7.85d0*(mvir/2d12)**(-0.081)/(1d0+z)**0.71
  ! One may also use the concentration parameter from Seljak (2000) 
  ! cvir=10d0*(mvir/3.42d12)**(-0.2)/(1d0+z)
  rs=rvir/cvir ! NFW scale radius, h^-1 Mpc
! convert to m500
  CALL mvir2mdel(mvir,rs,cvir,500d0*rhoc,m500)
  r500=(3d0*m500/4d0/pi/500d0/rhoc)**(1d0/3d0) ! h^-1 Mpc
  l500=da(z)/r500
! calculate the Fourier transform of Tsz
  ! See, e.g., Appendix D1 of Komatsu et al., arXiv:1001.4538
  !            Eq.(2) of Komatsu&Seljak (2002)
  transform2d=1.65d0*(h0/0.7d0)**2d0*Ez(z)**(8d0/3d0) &
       *(m500/3d14/0.7d0)**(2d0/3d0+0.12d0) &
       *rombint1(pgnfw,xin,xout,tol,ell/l500) &
       /0.5176d0*(4d0*3.14159d0)/l500**2d0*(r500/h0)
  transform2dp=1.65d0*(h0/0.7d0)**2d0*Ez(z)**(8d0/3d0) &
       *(m500/3d14/0.7d0)**(2d0/3d0+0.12d0) &
       *rombint1(pgnfw,xin,xout,tol,ellp/l500) &
       /0.5176d0*(4d0*3.14159d0)/l500**2d0*(r500/h0)
  Tsz=-2d0*283d0*(transform2d/50d0) ! uK, Rayleigh-Jeans limit
  Tszp=-2d0*283d0*(transform2dp/50d0) ! uK, Rayleigh-Jeans limit

! calculate mass function
  ! Sheth&Tormen and Tinker et al.'s mass functions are given for the 
  ! overdensity mass M200d (with respect to the mean mass density rather than 
  ! the critical density).
  CALL mvir2mdel(mvir,rs,cvir,200d0*omega*rhoc,m200d)
  Rh=(3d0*m200d/4d0/pi/om0/2.775d11)**(1d0/3d0) ! h^-1 Mpc
  ! Alternatively, one may wish to use Jenkins et al.'s mass function,
  ! which is given for the overdensity mass M180d (with respect to the mean 
  ! mass density rather than the critical density).
  ! CALL mvir2mdel(mvir,rs,cvir,180d0*omega*rhoc,m180d)
  ! Rh=(3d0*m180d/4d0/pi/om0/2.775d11)**(1d0/3d0) ! h^-1 Mpc
  lnRh=real(dlog(Rh))
  lnsigma2=CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
  dlnsigma2dlnRh=CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
  scalefactor=g(z)/g(zin)*(1d0+zin)/(1d0+z)
  lnnu=2d0*dlog(deltac)-dble(lnsigma2)-2d0*dlog(scalefactor) ! ln(nu)
  dlnnudlnRh=-dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
  dndlnRh=(3d0/4d0/pi)*dlnnudlnRh*mf(lnnu,z)/Rh**3d0
  ! dndlnRh=(3d0/4d0/pi)*dlnnudlnRh*mf(lnnu)/Rh**3d0 ! if using ST or Jenkins
  dndlnMh=dndlnRh/3d0 ! in units of h^3 Mpc^-3
! volume element
  dvdz=(1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
! This is the integrand: dV/dz*dn/dlnM*Tsz^2, to be integrated over lnM and z
  integrand3=dvdz*dndlnMh*Tsz**2d0*mgas(lnM,mcrit)*Tszp**2d0
  return
end function integrand3
!----------------------------------------------------------------------------
double precision function integrand2(lnM,z)
  USE cosmo
  USE sigma
  USE growth
  USE angular_distance
  USE nfw
  USE multipole
! The SZ power spectrum is given by C_l = int dz dV/dz int dlnM dn/dlnM*Tsz^2.
! This routine returns the integrand, dV/dz*dn/dlnM*Tsz^2. Note that Tsz is 
! in units of uK, Rayleigh-Jeans limit. This version uses Arnaud et al.'s 
! pressure profile.
! December 14, 2010: E.Komatsu
  IMPLICIT none
  double precision, intent(IN) :: lnM,z
  double precision :: xin=1d-5,xout=6d0 ! pressure profile is cut at r=xout*r500
  double precision :: omega,tol=1d-6,zin=30d0,pi=3.1415926535d0
  double precision :: mvir,rvir,rs,delc,m500,r500,l500
  double precision :: rhoc,Ez,da,dvdz,transform2d,tsz,pgnfw,rombint1
  double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh,lnnu,dlnnudlnRh,Rh
  double precision :: g,scalefactor,m180d,m200d
  double precision :: mgas,bias,lnnu2,bias_fix
  real :: chebev,lnsigma2,lnRh,dlnsigma2dlnRh
  integer :: i
  external pgnfw,rombint1
  mvir=dexp(lnM)
! compute omega's, rhoc, and delc
  omega=om0*(1d0+z)**3d0/Ez(z)**2d0 ! Omega(z); E(z)=H(z)/H0
  rhoc=2.775d11*Ez(z)**2d0 ! critical density in units of h^2 M_sun/Mpc^3
  ! Eq.(6) of Bryan & Norman (1998), valid only for a flat universe!
  delc=18d0*pi*pi+82d0*(omega-1d0)-39d0*(omega-1d0)**2d0
  ! One may also use Eq.(C19) of Nakamura & Suto (1997):
  ! delc=omega*18d0*pi*pi*(1d0+0.4093d0*(1d0/omega-1d0)**0.9052d0)
! compute virial and NFW parameters
  rvir=(3d0*mvir/4d0/pi/delc/rhoc)**(1d0/3d0) ! virial radius, h^-1 Mpc
  ! Concentration parameter from Duffy et al. (2008)
  cvir=7.85d0*(mvir/2d12)**(-0.081)/(1d0+z)**0.71
  ! One may also use the concentration parameter from Seljak (2000) 
  ! cvir=10d0*(mvir/3.42d12)**(-0.2)/(1d0+z)
  rs=rvir/cvir ! NFW scale radius, h^-1 Mpc
! convert to m500
  CALL mvir2mdel(mvir,rs,cvir,500d0*rhoc,m500)
  r500=(3d0*m500/4d0/pi/500d0/rhoc)**(1d0/3d0) ! h^-1 Mpc
  l500=da(z)/r500
! calculate the Fourier transform of Tsz
  ! See, e.g., Appendix D1 of Komatsu et al., arXiv:1001.4538
  !            Eq.(2) of Komatsu&Seljak (2002)
  transform2d=1.65d0*(h0/0.7d0)**2d0*Ez(z)**(8d0/3d0) &
       *(m500/3d14/0.7d0)**(2d0/3d0+0.12d0) &
       *rombint1(pgnfw,xin,xout,tol,ell/l500) &
       /0.5176d0*(4d0*3.14159d0)/l500**2d0*(r500/h0)
  Tsz=-2d0*283d0*(transform2d/50d0) ! uK, Rayleigh-Jeans limit
! calculate mass function
  ! Sheth&Tormen and Tinker et al.'s mass functions are given for the 
  ! overdensity mass M200d (with respect to the mean mass density rather than 
  ! the critical density).
  CALL mvir2mdel(mvir,rs,cvir,200d0*omega*rhoc,m200d)
  Rh=(3d0*m200d/4d0/pi/om0/2.775d11)**(1d0/3d0) ! h^-1 Mpc
  ! Alternatively, one may wish to use Jenkins et al.'s mass function,
  ! which is given for the overdensity mass M180d (with respect to the mean 
  ! mass density rather than the critical density).
  ! CALL mvir2mdel(mvir,rs,cvir,180d0*omega*rhoc,m180d)
  ! Rh=(3d0*m180d/4d0/pi/om0/2.775d11)**(1d0/3d0) ! h^-1 Mpc
  lnRh=real(dlog(Rh))
  lnsigma2=CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
  dlnsigma2dlnRh=CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
  scalefactor=g(z)/g(zin)*(1d0+zin)/(1d0+z)
  lnnu=1.0d0*(2d0*dlog(deltac)-dble(lnsigma2)-2d0*dlog(scalefactor)) ! ln(nu)
  dlnnudlnRh=-dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
  dndlnRh=(3d0/4d0/pi)*dlnnudlnRh*mf(lnnu,z)/Rh**3d0
  ! dndlnRh=(3d0/4d0/pi)*dlnnudlnRh*mf(lnnu)/Rh**3d0 ! if using ST or Jenkins
  dndlnMh=dndlnRh/3d0 ! in units of h^3 Mpc^-3
! volume element
  dvdz=(1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
! This is the integrand2: dV/dz*dn/dlnM*Tsz^2, to be integrated over lnM and z
  !print*,dvdz,dndlnMh,Tsz**2d0,mgas(lnM,mcrit)
  integrand2=dndlnMh*Tsz*mgas(lnM,mcrit)*bias(lnnu)
   !print*,log10(EXP(lnM)),exp(0.5*lnnu),bias(lnnu),integrand2

  return
end function integrand2
!-----------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION Ez(z)
  USE cosmo
  ! Ez = H(z)/H0 (dimensionless)
  ! x=a
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: z
  double precision :: x,ok0,or0=0d0 ! radiation density is ignored!
  ok0=1d0-om0-ode0
  x=1d0/(1d0+z)
  Ez=dsqrt(om0/x**3d0+ok0/x**2d0+or0/x**4d0+ode0/x**(3d0+3d0*w))
  return
END FUNCTION Ez
!-----------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION pgnfw(x,y) ! x=r/r500 & y=l/l500
  USE cosmo
  IMPLICIT none
  double precision :: x,y
  double precision :: c=1.177d0,g=0.3081d0,a=1.0510d0,b=5.4905d0,P0=8.403d0
  pgnfw=P0*(0.7d0/h0)**1.5d0/(c*x)**g/(1d0+(c*x)**a)**((b-g)/a) &
       *x**2d0*dsin(y*x)/(y*x)
  return
END FUNCTION pgnfw
!-------------
DOUBLE PRECISION FUNCTION mgas(mass,mcrit)
  double precision :: mass,mcrit
  mgas=1/(1+EXP(mcrit-mass))
  return
END FUNCTION mgas
!———————
DOUBLE PRECISION FUNCTION bias(lnnu)
    double precision :: lnnu
    double precision:: delta_c, delta, y, A, a1,B,b1,C,c1,v
    v = EXP(0.5d0*lnnu) ! silly conventions…
    delta_c = 1.686d0
    delta = 200.0d0 !definition of overdensity size
    y = LOG10(delta)
    A= 1.0d0 + 0.24d0*y*EXP(-1.0d0*(4.0d0/y)**4.0d0)
    a1= 0.44d0*y - 0.88d0
    B= 0.183d0
    b1= 1.5d0
    C= 0.019d0 + 0.107d0*y + 0.19d0*EXP(-1*(4.0d0/y)**4.0d0)
    c1= 2.4d0
    !print*,A,a1,B,b1,C,c1
    bias=1.0d0-(A*v**a1)/(v**a1+delta_c**a1)+B*v**b1 + C*v**c1
  return
END FUNCTION bias

DOUBLE PRECISION FUNCTION bias_fix(lnM)
    double precision :: lnM
    bias_fix=2.0d0*lnM - 26.0d0
    return
END FUNCTION bias_fix

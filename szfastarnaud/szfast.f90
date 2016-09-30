PROGRAM szfast
  USE cosmo
  USE linearpk
  USE sigma
  USE growth
  USE angular_distance
  USE multipole
! This code calculates the power spectrum of the Sunyaev-Zel'dovich effect.
! The methodology is described in Komatsu & Seljak, 336, 1256 (2002).
! The output, "multipole_szpower.txt", contains:
! 1st: multipole
! 2nd: l(l+1)Cl/twopi [uK^2] in the Rayleigh-Jeans limit [i.e., dT/T=-2y]
! December 14, 2010: E.Komatsu
  IMPLICIT none
  double precision :: z1=1d-5,z2=5d0 ! redshift limits
  character(len=128) :: filename
  integer :: npk,i,j,n=1001 ! 201 is good enough for l>100
  double precision :: cl,dcldlnx,lnx,dlnx,lnx1,lnx2
  external dcldlnx
! Specify three cosmological parameters
! The data type has been defined in MODULE cosmo.
  obh2=0.02262d0
  om0=0.277d0
  ode0=0.723d0
  h0=0.702d0
  w=-1d0
! read in and tabulate P(k)
  filename='wmap5baosn_max_likelihood_matterpower_at_z=30.dat'
  npk=896 ! # of lines in the file
  CALL open_linearpk(filename,npk)
! fit sigma^2(R) to Chebyshev polynomials
  CALL compute_sigma2
  CALL close_linearpk
! compute the growth factor
  CALL setup_growth
! compute and tabulate da(z) 
  CALL setup_da
! compute the SZ power spectrum, l(l+1)Cl/twopi [uK^2] 
! in the RJ limit, dT/T=-2y
  open(1,file='multipole_szpower.txt')
  lnx1=dlog(1d0+z1)
  lnx2=dlog(1d0+z2)
  dlnx=(lnx2-lnx1)/dble(n-1)
  do j=1,16
     ell=10d0**(1d0+dble(j-1)*0.2)
     ! integrate by the trapezoidal rule
     cl=0d0
     cl=dcldlnx(lnx1)*(0.5d0*dlnx)
     cl=cl+dcldlnx(lnx2)*(0.5d0*dlnx)
     do i=2,n-1
        lnx=lnx1+dble(i-1)*dlnx
        cl=cl+dcldlnx(lnx)*dlnx
     enddo
     print*,ell,cl*ell*(ell+1d0)/2d0/3.14159d0
     write(1,'(2E15.5)')ell,cl*ell*(ell+1d0)/2d0/3.14159d0
  enddo
  close(1)
END PROGRAM szfast
!--------------------------------------------------------------
double precision function dcldlnx(lnx)
  ! dcldlnx = (1+z) (dV/dz) int dlnM dn/dlnM Tsz^2
  ! lnx=ln(1+z)
  IMPLICIT none
  double precision, intent(IN) :: lnx
  double precision :: M1=5d11,M2=5d15 ! mass limits, h^-1 Msun
  double precision :: integrand,lnM1,lnM2,lnM,z
  integer :: i
  external integrand
  z=exp(lnx)-1d0
  lnM1=dlog(M1) ! minimum mass, h^-1 Msun
  lnM2=dlog(M2) ! maximum mass, h^-1 Msun
  CALL qgaus1(integrand,lnM1,lnM2,dcldlnx,z) ! Gaussian quadrature
  dcldlnx=dcldlnx*(1d0+z)
  return
end function dcldlnx

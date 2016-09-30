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
  character(len=128) :: filename,file2
  character(len=32) :: arg
  integer :: npk,i,j,n=1001 ! 201 is good enough for l>100
  double precision :: cl,dcldlnx,lnx,dlnx,lnx1,lnx2
  double precision :: cl2, dcl2dlnx
  !double precision, dimension(15) :: ell_array
  !ell_array = (/ 3000.0d0 ,  1247.0d0 ,   959.0d0 ,   738.0d0 ,   567.5d0,   436.5d0,   335.50d0, 257.50d0,   198.0d0 ,   152.50d0,   117.0d0 ,    89.50d0,    68.50d0,    52.50d0, 40.0d0 /)
  external dcldlnx

  REAL, DIMENSION(15) :: ell_list
  ell_list = (/ 3000.0, 1247.0, 959.0, 738.0, 567.5, 436.5, 335.50, 257.50, 198.0, 152.50, 117.0, 89.50, 68.50, 52.50, 40.0 /)
! Specify three cosmological parameters
! The data type has been defined in MODULE cosmo.
  obh2=0.02289d0
  om0=0.260d0
  ode0=0.723d0
  h0=0.717d0
  w=-1d0
  CALL getarg(1, arg)
  read (arg,*) mcrit 
! read in and tabulate P(k)
  CALL getarg(2, arg)
  read (arg,*) filename
  !filename='test_out'!’wmap5baosn_max_likelihood_matterpower_at_z=30.dat'
  npk=705 ! # of lines in the file
  CALL open_linearpk(filename,npk)
! fit sigma^2(R) to Chebyshev polynomials
  CALL compute_sigma2
  !CALL close_linearpk
! compute the growth factor
  CALL setup_growth
! compute and tabulate da(z) 
  CALL setup_da
! compute the SZ power spectrum, l(l+1)Cl/twopi [uK^2] 
! in the RJ limit, dT/T=-2y
  CALL getarg(3, arg)
  read (arg,*) file2
  open(1,file=file2)
  lnx1=dlog(1d0+z1)
  lnx2=dlog(1d0+z2)
  dlnx=(lnx2-lnx1)/dble(n-1)
  do j=1,15
     ell=dble(ell_list(j))
     ! integrate by the trapezoidal rule
     cl=0d0
     cl=dcldlnx(lnx1)*(0.5d0*dlnx)
     cl=cl+dcldlnx(lnx2)*(0.5d0*dlnx)

     cl2=0d0
     cl2=dcl2dlnx(lnx1,ell)*(0.5d0*dlnx)
     cl2=cl2+dcl2dlnx(lnx2,ell)*(0.5d0*dlnx)

     do i=2,n-1
        lnx=lnx1+dble(i-1)*dlnx
        cl=cl+dcldlnx(lnx)*dlnx
        cl2=cl2+dcl2dlnx(lnx,ell)*dlnx
     enddo
     print*,ell,cl*ell*(ell+1d0)/2d0/3.14159d0,cl2*ell*(ell+1d0)/2d0/3.14159d0
     write(1,*)ell,cl,cl2
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

!——
double precision function dcl2dlnx(lnx,ell)
  ! dcl2dlnx = (1+z) (dV/dz) (int dlnM dn/dlnM Tsz*f_gas*bias)**2 p(l+1/2,z)                                         
  ! lnx=ln(1+z)                                                                            
  !IMPLICIT none
  double precision, intent(IN) :: lnx
  double precision :: M1=5d11,M2=5d15 ! mass limits, h^-1 Msun                             
  double precision :: integrand2,lnM1,lnM2,lnM
  double precision :: z
  double precision :: Ez,da,dvdz,comoving
  double precision :: linear_pk,g
  integer :: i
  external integrand2
  external linear_pk
  z=exp(lnx)-1d0
  !print*,z,(1d0+z)**2d0*da(z)**2d0,Ez(z)
  dvdz=(1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
  lnM1=dlog(M1) ! minimum mass, h^-1 Msun                                                  
  lnM2=dlog(M2) ! maximum mass, h^-1 Msun                                                  
  CALL qgaus1(integrand2,lnM1,lnM2,dcl2dlnx,z) ! Gaussian quadrature
  !print*,dcl2dlnx**2.0d0,dvdz,(g(z)/(1.0d0+z))**2.0d0,linear_pk(DBLE((ell+1)/da(z)))
  comoving = da(z)*(1+z)
  !print*,linear_pk(DBLE((ell+0.5d0)/comoving))
  dcl2dlnx=dcl2dlnx**2.0d0*(1d0+z)**1.0d0*dvdz*linear_pk(DBLE((ell+0.5d0)/comoving))*(g(z)/(1.0d0+z)*(31.0d0)/g(30.0d0))**2.0d0
  return
end function dcl2dlnx


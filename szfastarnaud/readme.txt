*** Power Spectrum of the Sunyaev-Zel'dovich Effect ***
This version uses Arnaud et al.'s pressure profile.
December 14, 2010: E.Komatsu

Ref: Komatsu & Seljak, MNRAS, 336, 1256 (2002)
     Arnaud et al., A&A, 517, 92 (2010)

This code calculates the power spectrum of the Sunyaev-Zel'dovich (SZ)
effect:

Cl = int dz dV/dz int dlnM dn/dlnM Tsz^2

where dV/dz is the comoving volume element, dn/dlnM is the halo mass 
function, and Tsz is the 2d Fourier transform of the SZ effect in the 
Rayleigh-Jeans limit, i.e., Tsz=-2y*2.725e6 [uK]. To calculate Cl at
different frequencies, you will need to multiply it by [gnu/(-2)]^2, 
where gnu is the appropriate spectrum function of the SZ effect.

** This code is different from the original version used for 
Komatsu & Seljak (2002) in three respects:
1. Mass function: this code uses Tinker et al. (2008); KS02 used 
Jenkins et al. (2001). 
2. Concentration: this code uses Duffy et al. (2008); KS02 used
Seljak (2000).
3. Pressure profile: this code uses Arnaud et al. (2010); KS02 used
Komatsu&Seljak (2001).


To compute the SZ power spectrum, it is necessary to use the linear 
matter power spectrum. We provide the sample data, 
"wmap5baosn_max_likelihood_matterpower.dat," which was generated using CAMB 
code for the maximum likelihood parameters given in Table I of Komatsu et al.
(2008) [WMAP 5-year interpretation paper] with "WMAP5+BAO+SN". The input file 
for CAMB is also provided (wmap5baosn_max_likelihood_params.ini). NOTE THAT 
THIS POWER SPECTRUM IS COMPUTED AT Z=0.

In the code, we use the matter power spectrum evolved back to z=30. This
isprovided as "wmap5baosn_max_likelihood_matterpower_at_z=30.dat".
(The above file at z=0 is not used.)

<< NOTE ON NFW CONCENTRATION PARAMETER >>

The default concentration parameter is Duffy et al. (2008):
cvir=7.85d0*(mvir/2d12)**(-0.081)/(1d0+z)**0.71

You may also use the concentration parameter of Seljak (2000):
cvir=10d0*(mvir/3.42d12)**(-0.2)/(1d0+z)

which was originally used by Komatsu&Seljak (2002).

<< NOTE ON MASS FUNCTION >>

The default mass function is Eq.(3) of Tinker et al. (2008) with the redshift
evolution factors given by Eq.(5)-(8).

Two more mass functions are provided: 
- Sheth&Tormen (mf_shethtormen.f90) and 
- Jenkins et al. (mf_jenkins.f90; which was originally used by Komatsu&Seljak 2002).

To use these mass functions, change the following:

1. Change "mf_tinker_redshift.o" to "mf_jenkins.o" or "mf_shethtormen.o"
in Makefile, 

2. Change "mf(lnnu,z)" to "mf(lnnu)" in integrand.f90 [because Jenkins et al.'s 
and Sheth&Tormen's functions do not have explicit redshift dependence].
Specifically: 
dndlnRh=(3d0/4d0/pi)*dlnnudlnRh*mf(lnnu,z)/Rh**3d0

to 

dndlnRh=(3d0/4d0/pi)*dlnnudlnRh*mf(lnnu)/Rh**3d0

in integrand.f90.

3.  Also, since Jenkins et al. use a different mass definition, you need to 
comment out

  ! Sheth&Tormen and Tinker et al.'s mass functions are given for the 
  ! overdensity mass M200d (with respect to the mean mass density rather than 
  ! the critical density).
  CALL mvir2mdel(mvir,rs,cvir,200d0*omega*rhoc,m200d)
  Rh=(3d0*m200d/4d0/pi/om0/2.775d11)**(1d0/3d0) ! h^-1 Mpc

and use

  ! Alternatively, one may wish to use Jenkins et al.'s mass function,
  ! which is given for the overdensity mass M180d (with respect to the mean 
  ! mass density rather than the critical density).
  CALL mvir2mdel(mvir,rs,cvir,180d0*omega*rhoc,m180d)

in "integrand.f90".



- To compile and use the program, edit Makefile and simply "./make"
- It will generate an executable called "szfast"
- Running szfast will generate the data file called "multipole_szpower.txt", 
which contains:

1st: multipole
2nd: l(l+1)Cl/twopi [uK^2]

For your convenience, the result from the default options (Duffy et al.'s 
concentration parameter and Tinker et al.'s mass function) is provided as
"multipole_szpower_default.txt".

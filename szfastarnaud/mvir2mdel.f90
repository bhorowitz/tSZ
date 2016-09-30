!-----------------------------------------------------------------------------
!
! This subroutine converts the virial mass, Mvir, to an over-density mass,
! Mdel, for a given over-density, delrho. The dark matter density profile
! is given by a Navarro-Frenk-White profile, which is parametrized by 
! a scale radius (rs) and a concentration parameter (c). 
!
! INPUTS: Mvir, rs, c, rhodel [double precision]
! OUTPUT: Mdel [double precision]
!
! For details, see, e.g., Eq (D15) of Komatsu et al. (WMAP 7-year paper),
! arXiv:1001.4538.
!
! December 10, 2010
! E. Komatsu
! (originally written in February 14, 2002)
!
!-----------------------------------------------------------------------------
SUBROUTINE mvir2mdel(Mvir,rs,c,delrho,Mdel)
  IMPLICIT none
  double precision, intent(IN) :: Mvir,rs,c,delrho
  double precision, intent(OUT):: Mdel
  integer :: i,j,iter,imax=50
  double precision, allocatable :: mtest(:)
  double precision :: ltest,m1,m2,ml,mu,lnMdel,var(4)
  double precision :: zbrent,mv2md
  external zbrent,mv2md
  var(1:4)=(/Mvir,rs,c,delrho/)
  allocate(mtest(0:imax))
! Find initial mass range for zbrent (non-linear equation solver),
! within which Mdel will be found.
  mtest(0)=Mvir ! initial guess for the overdensity mass
  ltest=mv2md(dlog(mtest(0)),var)
  if(ltest.le.0d0)then
     do i=1,imax
        mtest(i)=mtest(i-1)*2.d0
        ltest=mv2md(dlog(mtest(i)),var)
        if (ltest.gt.0d0)then
           m1=dlog(mtest(i))
           m2=dlog(mtest(i-1))
           exit
        endif
     enddo
  else
     do i=1,imax
        mtest(i)=mtest(i-1)/2.d0
        ltest=mv2md(dlog(mtest(i)),var)
        if(ltest.lt.0d0)then
           m1=dlog(mtest(i))
           m2=dlog(mtest(i-1))
           exit
        endif
     enddo
  endif
! Lower & upper limits for M, between which Mdel will be found.
  ml=min(m1,m2) 
  mu=max(m1,m2)
! Find Mdel using zbrent.
  lnMdel=zbrent(mv2md,ml,mu,1d-4,iter,var)
  Mdel=dexp(lnMdel)
  return
END SUBROUTINE mvir2mdel
!-----------------------------------------------------------------------------
double precision function mv2md(lnM,var)
! zbrent will solve the equation, mv2md(lnMdel)=0, for lnMdel.
  IMPLICIT none
  double precision, intent(IN) :: lnM,var(4)
  double precision :: Mdel,rdel,func,fourpi
  double precision :: Mvir,rs,c,delrho
  external func
  fourpi=4d0*3.141592653589793238d0
  Mdel=dexp(lnM)
  Mvir=var(1)
  rs=var(2)
  c=var(3)
  delrho=var(4)
  rdel=(Mdel/(fourpi/3d0*delrho))**(1d0/3d0)
  mv2md=Mdel/Mvir-func(rdel/rs)/func(c) 
        ! See, e.g., Eq.(D15) of Komatsu et al., arXiv:1001.4538.
  return
END function mv2md
!-----------------------------------------------------------------------------
double precision function func(x)
  IMPLICIT none
  double precision, intent(IN) :: x
  func=dlog(1d0+x)-x/(1d0+x)
  return
END function func

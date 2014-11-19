!>
!
! flux at the wall
!
subroutine wallFlux(nxyz,ql,flx,specRadius,gamma)
!
implicit none
!
real*8, intent(in) :: nxyz(2)
real*8, intent(in) :: ql(4)
real*8, intent(inout) :: flx(4)
real*8, intent(inout) :: specRadius
real*8, intent(in) :: gamma
!
real*8 :: p,u,v,a,ri
!
p=(gamma-1)*(ql(4)-0.5*(ql(3)**2+ql(2)**2)/ql(1))
ri=1./ql(1)
u=ql(2)*ri
v=ql(3)*ri
a=sqrt(gamma*p*ri);
!
specRadius=(sqrt(u*u+v*v)+a)*sqrt(nxyz(1)**2+nxyz(2)**2)
!
flx(1)=0
flx(2)=p*nxyz(1)
flx(3)=p*nxyz(2)
flx(4)=0
!
return
end subroutine wallFlux
!>
!
! Jacobian of wall flux
!
subroutine wallFluxJacobian(nxyz,ql,lmat,rmat,gamma)
!
implicit none
!
real*8, intent(in) :: nxyz(2)
real*8, intent(in) :: ql(4)
real*8, intent(inout) :: lmat(4,4)
real*8, intent(inout) :: rmat(4,4)
real*8, intent(in) :: gamma
!
real*8 :: u,v,umag2,gm1
!
u=ql(2)/ql(1)
v=ql(3)/ql(1)
umag2=u*u+v*v
gm1=gamma-1
!
lmat(1,:)=0.
lmat(4,:)=0.
!
lmat(2,1)=gm1*0.5*umag2*nxyz(1)
lmat(2,2)=-gm1*u*nxyz(1)
lmat(2,3)=-gm1*v*nxyz(1)
lmat(2,4)=gm1*nxyz(1)
!
lmat(3,1)=gm1*0.5*umag2*nxyz(2)
lmat(3,2)=-gm1*u*nxyz(2)
lmat(3,3)=-gm1*v*nxyz(2)
lmat(3,4)=gm1*nxyz(2)
!
rmat=0.
!
return
end subroutine wallFluxJacobian

subroutine flux_visc_2d(faceNormal2,vol,ql2,qr2,flux,&
     rey,pr,prtr,gamma,c2b,rgas)
!
implicit none
!
real*8, intent(in) :: faceNormal2(2)
real*8, intent(in) :: vol
real*8, intent(in) :: ql2(4)
real*8, intent(in) :: qr2(4)
real*8, intent(inout) :: flux(4)
real*8, intent(in) :: rey,pr,prtr,gamma,c2b,rgas
!
real*8 :: faceNormal(3)
real*8 :: faceSpeed=0.
real*8 :: ql(6),qr(6)
real*8 :: fluxv(6)
!
faceNormal(1:2)=faceNormal2
faceNormal(3)=0.
!
ql(1:3)=ql2(1:3)
ql(4)=0.
ql(5)=ql2(4)
ql(6)=0.
!
qr(1:3)=qr2(1:3)
qr(4)=0.
qr(5)=qr2(4)
qr(6)=0.
!
call flux_visc(faceNormal,faceSpeed,vol,ql,qr,fluxv,&
     rey,pr,prtr,gamma,c2b,Rgas)
!
flux(1)=fluxv(1)
flux(2)=fluxv(2)
flux(3)=fluxv(3)
flux(4)=fluxv(5)
!
return
end subroutine flux_visc_2d
!
! compute the 2nd order viscous fluxes in a given direction
! no cross terms here.
!
subroutine flux_visc(faceNormal,faceSpeed,vol,ql,qr,flux,&
                    rey,pr,prtr,gamma,c2b,Rgas)
!
implicit none
!
real*8, intent(in)  :: ql(6)
real*8, intent(in)  :: qr(6)
real*8, intent(in)  :: faceNormal(3)
real*8, intent(in) ::  faceSpeed(1)
real*8, intent(in) :: vol
real*8, intent(out) :: flux(6)
real*8, intent(in) :: rey,pr,prtr,gamma,c2b,Rgas
!
real*8 :: ur(3),ul(3)
integer :: i,ip
real*8  :: gm1
real*8  :: rj1,rj2,rj3
real*8  :: a1,a2,a3,a4,a5,a6,a7
real*8  :: dre,sigma,rgasi,gm1Pr
real*8  :: t1,mu,t2,mut,rr,fac
real*8  :: rhoi,rhoi1,vi
real*8  :: ux,vx,wx,uux,vvx,wwx,uvx,vwx,wux,tx
real*8  :: plft,prht
!
real*8, save  :: third=1./3
!
gm1=gamma-1.
dre=1./rey
sigma=pr/prtr
rgasi=1./Rgas
gm1Pr=1./gm1/Pr
!
!
rhoi=1./ql(1)
rhoi1=1./qr(1)
!
! sutherland's law for laminar viscosity
!
! pressure (left)
plft=(gamma-1)*(ql(5)-0.5*(ql(2)**2+ql(3)**2+ql(4)**2)/ql(1))
! temperature (left)
t1=rgasi*plft*rhoi
! sutherland
mu=(c2b+1.)*t1*sqrt(t1)/(c2b+t1)
! pressure (right)
prht=(gamma-1)*(qr(5)-0.5*(qr(2)**2+qr(3)**2+qr(4)**2)/qr(1))
! temperature (right)
t2=rgasi*prht*rhoi1
mu=0.5*(mu+(c2b+1)*t2*sqrt(t2)/(c2b+t2))
!
! turbulent viscosity
!
mut=0.5*(ql(6)+qr(6))
!
! face normals (scaled by face area), i.e not unit normals
!
vi=1./vol
rj1=faceNormal(1)
rj2=faceNormal(2)
rj3=faceNormal(3)
!
a4=rj1**2+rj2**2+rj3**2
rr=sqrt(a4)
a1=a4+third*rj1**2
a2=a4+third*rj2**2
a3=a4+third*rj3**2
a5=third*rj1*rj2
a6=third*rj2*rj3
a7=third*rj3*rj1
!
! 2nd order gradients at the interface
! need to change here to make the routine higher order
!
ur=qr(2:4)/qr(1)
ul=ql(2:4)/ql(1)
!
ux=ur(1)-ul(1)
vx=ur(2)-ul(2)
wx=ur(3)-ul(3)
uux=ux*(ur(1)+ul(1))
vvx=vx*(ur(2)+ul(2))
wwx=wx*(ur(3)+ul(3))
uvx=(ur(1)*ur(2)-ul(1)*ul(2))
vwx=(ur(2)*ur(3)-ul(2)*ul(3))
wux=(ur(3)*ur(1)-ul(3)*ul(1))
!
tx = (t2-t1)*gm1Pr
!
! no density contribution
!
flux(1)=0.
!
! momentum equations
!
flux(2) = a1*ux+a5*vx+a7*wx
flux(3) = a5*ux+a2*vx+a6*wx
flux(4) = a7*ux+a6*vx+a3*wx
!
fac=dre*vi !vol(i) !*rr
flux(2:4) = (mu+mut)*fac*flux(2:4)
!
! energy equation
!
flux(5) = 0.5*a1*uux + 0.5*a2*vvx + 0.5*a3*wwx + &
              a5*uvx +     a6*vwx +     a7*wux + &
              a4*tx
flux(5) = (mu+mut*sigma)*fac*flux(5)
flux(6) = 0.
!
return
end subroutine flux_visc
   
   
  
   
   


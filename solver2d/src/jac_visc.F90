#ifdef UNIT_TEST
program testjac
  !
  implicit none
  real*8 :: ql(6),qr(6)
  real*8 :: dql(6),dqr(6)
  real*8 :: ql1(6),qr1(6)
  real*8 :: df1(6)
  real*8 :: flux(6),flux1(6)
  real*8 :: dfluxl(6,6)
  real*8 :: dfluxr(6,6)
  real*8 :: dfluxl1(6,6)
  real*8 :: dfluxr1(6,6)
  real*8 :: faceNormal(3)
  real*8 :: faceSpeed
  real*8 :: vol = 1.4e-4
  real*8 :: rey = 1e5
  real*8 :: pr=0.72
  real*8 :: prtr=0.3333
  real*8 :: gamma=1.4
  real*8 :: c2b=0.3678
  real*8 :: Rgas=1./1.4
  real*8 :: pinf
  real*8 :: drand
  !
  pinf=1./gamma
  !
  ql(1)=1.
  ql(2)=0.3
  ql(3)=0.02
  ql(4)=0.01
  ql(5)=pinf/(gamma-1)+0.5*(ql(2)**2+ql(3)**2+ql(4)**2)/ql(1)
  !
  qr(1)=1.
  qr(2)=0.2
  qr(3)=0.01
  qr(4)=0.005
  qr(5)=pinf/(gamma-1)+0.5*(qr(2)**2+qr(3)**2+qr(4)**2)/qr(1)
  !
  faceNormal(1)=-2.
  faceNormal(2)=-3.
  faceNormal(3)=1.
  faceSpeed=0.
  !
  call jac_visc(faceNormal,faceSpeed,vol,ql,qr,flux,dfluxl,dfluxr,&
       rey,pr,prtr,gamma,c2b,Rgas)
  write(6,*) 'flux=',flux
  call flux_visc(faceNormal,faceSpeed,vol,ql,qr,flux,&
       rey,pr,prtr,gamma,c2b,Rgas)
  !
  write(6,*) 'flux=',flux
  !
  dql=0
  dqr=0
  dql=1e-6
  dqr=1e-6
  !
  ql1=ql+dql
  qr1=qr
  !
  call jac_visc(faceNormal,faceSpeed,vol,ql1,qr1,flux1,dfluxl1,dfluxr1,&
       rey,pr,prtr,gamma,c2b,Rgas)
  !
  write(6,*) 'flux1=',flux1
  !
  write(6,*) 'flux1-flux=',flux1-flux
  write(6,*) 'dfluxl*dql=',matmul(dfluxl,dql)
  df1=matmul(dfluxl,dql)
  write(6,*) 'error=',df1-(flux1-flux)
  write(6,*) '-----------------------------------------------------------'
  !
  ql1=ql
  qr1=qr+dqr
  !
  call jac_visc(faceNormal,faceSpeed,vol,ql1,qr1,flux1,dfluxl1,dfluxr1,&
       rey,pr,prtr,gamma,c2b,Rgas)
  !
  write(6,*) 'flux2=',flux1
  write(6,*) 'flux2-flux=',flux1-flux
  write(6,*) 'dfluxr*dqr=',matmul(dfluxr,dqr)
  write(6,*) 'error=', matmul(dfluxr,dqr)-(flux1-flux)
  write(6,*) '-----------------------------------------------------------'
  !
end program testjac
#endif
!
!> 2-D viscous flux and jacobian
!>
subroutine jac_visc_2d(faceNormal2,vol,ql2,qr2,flux,&
     lmat,rmat,rey,pr,prtr,gamma,c2b,rgas)
!
implicit none
!
real*8, intent(in) :: faceNormal2(2)
real*8, intent(in) :: vol
real*8, intent(in) :: ql2(4)
real*8, intent(in) :: qr2(4)
real*8, intent(inout) :: flux(4)
real*8, intent(inout) :: lmat(4,4),rmat(4,4)
real*8, intent(in) :: rey,pr,prtr,gamma,c2b,rgas
!
real*8 :: faceNormal(3)
real*8 :: faceSpeed=0.
real*8 :: ql(6),qr(6)
real*8 :: fluxv(6)
real*8 :: dfluxl(6,6),dfluxr(6,6)
integer :: imap(4)
integer :: i,j
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
call jac_visc(faceNormal,faceSpeed,vol,ql,qr,fluxv,dfluxl,dfluxr,&
     rey,pr,prtr,gamma,c2b,Rgas)

!
flux(1)=fluxv(1)
flux(2)=fluxv(2)
flux(3)=fluxv(3)
flux(4)=fluxv(5)
!
imap=(/1,2,3,5/)
!
do i=1,4
   do j=1,4
      lmat(i,j)=dfluxl(imap(i),imap(j))
      rmat(i,j)=dfluxr(imap(i),imap(j))
   enddo
enddo
!
return
end subroutine jac_visc_2d

!
! compute Jacobian of the 2nd order viscous fluxes in a given direction
! no cross terms here.
!
!
! Jay Sitaraman
! 12/05/13
!
subroutine jac_visc(faceNormal,faceSpeed,vol,ql,qr,flux,dfluxl,dfluxr,&
                    rey,pr,prtr,gamma,c2b,Rgas)
!
implicit none
!
real*8, intent(in)  :: ql(6)
real*8, intent(in)  :: qr(6)
real*8, intent(in)  :: faceNormal(3)
real*8, intent(in) ::  faceSpeed(1)
real*8, intent(in) :: vol
real*8, intent(inout) :: flux(6)
real*8, intent(inout) :: dfluxl(6,6),dfluxr(6,6)
real*8, intent(in) :: rey,pr,prtr,gamma,c2b,Rgas
!
real*8 :: ur(3),ul(3)
integer :: i,ip,n,m
real*8  :: gm1
real*8  :: rj1,rj2,rj3
real*8  :: a1,a2,a3,a4,a5,a6,a7
real*8  :: dre,sigma,rgasi,gm1Pr
real*8  :: t1,mu,t2,mut,rr,fac
real*8  :: rhoi,rhoi1,vi,rhoi2,rhoi12
real*8  :: ux,vx,wx,uux,vvx,wwx,uvx,vwx,wux,tx
real*8  :: plft,prht
!
! variables for linearization
!
real*8 :: ulmag,urmag
real*8 :: dmul(5)
real*8 :: dmur(5)
real*8 :: dmutl,dmutr
real*8 dplft(5),dprht(5)
real*8 dt1l(5),dt2r(5)
real*8 :: dull(3,5),durr(3,5)
real*8 duxl(5),duxr(5)
real*8 dvxl(5),dvxr(5)
real*8 dwxl(5),dwxr(5)
real*8 duuxl(5),duuxr(5)
real*8 dvvxl(5),dvvxr(5)
real*8 dwwxl(5),dwwxr(5)
real*8 duvxl(5),duvxr(5)
real*8 dvwxl(5),dvwxr(5)
real*8 dwuxl(5),dwuxr(5)
real*8 dtxl(5),dtxr(5)
!
real*8, save  :: third=1./3
!
gm1=gamma-1.
dre=1./rey
sigma=pr/prtr
rgasi=1./Rgas
gm1Pr=1./gm1/Pr
!
rhoi=1./ql(1)
rhoi1=1./qr(1)
rhoi2=rhoi*rhoi
rhoi12=rhoi1*rhoi1
!
ul=ql(2:4)/ql(1)
ulmag=ul(1)*ul(1)+ul(2)*ul(2)+ul(3)*ul(3)
!
ur=qr(2:4)/qr(1)
urmag=ur(1)*ur(1)+ur(2)*ur(2)+ur(3)*ur(3)
!
! velocity derivatives
!
dull(1,1)=-ul(1)*rhoi
dull(1,2)=rhoi
dull(1,3:5)=0.
!
dull(2,1)=-ul(2)*rhoi
dull(2,2)=0.
dull(2,3)=rhoi
dull(2,4:5)=0.
!
dull(3,1)=-ul(3)*rhoi
dull(3,2:3)=0
dull(3,4)=rhoi
dull(3,5:5)=0
!
durr(1,1)=-ur(1)*rhoi1
durr(1,2)=rhoi1
durr(1,3:5)=0.
!
durr(2,1)=-ur(2)*rhoi1
durr(2,2)=0.
durr(2,3)=rhoi1
durr(2,4:5)=0.
!
durr(3,1)=-ur(3)*rhoi1
durr(3,2:3)=0
durr(3,4)=rhoi1
durr(3,5)=0
!
! pressure (left)
!
plft=gm1*(ql(5)-0.5*ql(1)*ulmag)
!
dplft(1)=0.5*gm1*ulmag
dplft(2)=-gm1*ul(1)
dplft(3)=-gm1*ul(2)
dplft(4)=-gm1*ul(3)
dplft(5)=gm1
!
! temperature (left)
!
t1=rgasi*plft*rhoi
!
do n=1,5
   dt1l(n)=rgasi*dplft(n)*rhoi
enddo
dt1l(1)=dt1l(1)-rgasi*plft*rhoi2
!
mu=(c2b+1.)*t1*sqrt(t1)/(c2b+t1)
do n=1,5
   dmul(n)= (c2b+1)*1.5*dt1l(n)*sqrt(t1)/(c2b+t1)+&
           -(c2b+1)*t1*sqrt(t1)/((c2b+t1)**2)*dt1l(n)
enddo
!
! pressure (right)
!
prht=gm1*(qr(5)-0.5*qr(1)*urmag)
!
dprht(1)=gm1*0.5*urmag
dprht(2)=-gm1*ur(1)
dprht(3)=-gm1*ur(2)
dprht(4)=-gm1*ur(3)
dprht(5)=gm1
!
! temperature (right)
!
t2=rgasi*prht*rhoi1
!
do n=1,5
   dt2r(n)=rgasi*dprht(n)*rhoi1
enddo
dt2r(1)=dt2r(1)-rgasi*prht*rhoi12
!
mu=0.5*(mu+(c2b+1)*t2*sqrt(t2)/(c2b+t2))
!
do n=1,5
   dmul(n)=0.5*dmul(n)
   dmur(n)=0.5*((c2b+1)*1.5*dt2r(n)*sqrt(t2)/(c2b+t2)+&
               -(c2b+1)*t2*sqrt(t2)/((c2b+t2)**2)*dt2r(n))
enddo
!mu=1.
!dmul=0
!dmur=0
!
! turbulent viscosity
!
mut=0.5*(ql(6)+qr(6))
dmutl=0.5
dmutr=0.5
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
! pure 2nd order gradients
! of velocity
!
ux=ur(1)-ul(1)
vx=ur(2)-ul(2)
wx=ur(3)-ul(3)
!
do n=1,5
   duxl(n)=-dull(1,n)
   duxr(n)=durr(1,n)
   !
   dvxl(n)=-dull(2,n)
   dvxr(n)=durr(2,n)
   !
   dwxl(n)=-dull(3,n)
   dwxr(n)=durr(3,n)
enddo
!
! velocity tensor product
! terms (2nd order difference)
!
uux=ux*(ur(1)+ul(1))
vvx=vx*(ur(2)+ul(2))
wwx=wx*(ur(3)+ul(3))
uvx=(ur(1)*ur(2)-ul(1)*ul(2))
vwx=(ur(2)*ur(3)-ul(2)*ul(3))
wux=(ur(3)*ur(1)-ul(3)*ul(1))
!
! derivatives of this term w.r.t ql and qr
!
do n=1,5
   !
   duuxl(n)=duxl(n)*(ur(1)+ul(1))+ux*dull(1,n)
   duuxr(n)=duxr(n)*(ur(1)+ul(1))+ux*durr(1,n)
   !
   dvvxl(n)=dvxl(n)*(ur(2)+ul(2))+vx*dull(2,n)
   dvvxr(n)=dvxr(n)*(ur(2)+ul(2))+vx*durr(2,n)
   !
   dwwxl(n)=dwxl(n)*(ur(3)+ul(3))+wx*dull(3,n)
   dwwxr(n)=dwxr(n)*(ur(3)+ul(3))+wx*durr(3,n)
   !
   duvxl(n)=-dull(1,n)*ul(2)-ul(1)*dull(2,n)
   !
   duvxr(n)=durr(1,n)*ur(2)+ur(1)*durr(2,n)
   !
   dvwxl(n)= -dull(2,n)*ul(3)-ul(2)*dull(3,n)            
   !
   dvwxr(n)=durr(2,n)*ur(3)+ur(2)*durr(3,n)
   !
   dwuxl(n)=-dull(3,n)*ul(1)-ul(3)*dull(1,n)            
   !
   dwuxr(n)=durr(3,n)*ur(1)+ur(3)*durr(1,n) 
   !
enddo
!
tx = (t2-t1)*gm1Pr
do n=1,5
   dtxl(n)=-dt1l(n)*gm1Pr
   dtxr(n)=dt2r(n)*gm1Pr
enddo
!
! no density contribution
!
flux(1)=0.
do n=1,6
   dfluxl(1,n)=0
   dfluxr(1,n)=0
enddo
!
! momentum equations
!
flux(2) = a1*ux+a5*vx+a7*wx
flux(3) = a5*ux+a2*vx+a6*wx
flux(4) = a7*ux+a6*vx+a3*wx
!
do n=1,5
   !
   dfluxl(2,n)=a1*duxl(n)+a5*dvxl(n)+a7*dwxl(n)
   dfluxl(3,n)=a5*duxl(n)+a2*dvxl(n)+a6*dwxl(n)
   dfluxl(4,n)=a7*duxl(n)+a6*dvxl(n)+a3*dwxl(n)
   !
   dfluxr(2,n)=a1*duxr(n)+a5*dvxr(n)+a7*dwxr(n)
   dfluxr(3,n)=a5*duxr(n)+a2*dvxr(n)+a6*dwxr(n)
   dfluxr(4,n)=a7*duxr(n)+a6*dvxr(n)+a3*dwxr(n)
   !
enddo
!
fac=dre*vi 
!
do m=2,4
   do n=1,5
      dfluxl(m,n)=fac*((mu+mut)*dfluxl(m,n)+dmul(n)*flux(m))
      dfluxr(m,n)=fac*((mu+mut)*dfluxr(m,n)+dmur(n)*flux(m))
   enddo
   n=6
   dfluxl(m,n)=fac*dmutl*flux(m)
   dfluxr(m,n)=fac*dmutr*flux(m)
enddo
!
flux(2:4) = (mu+mut)*fac*flux(2:4)
!
! energy equation
!
!
do n=1,5
   dfluxl(5,n)=0.5*a1*duuxl(n)+0.5*a2*dvvxl(n)+0.5*a3*dwwxl(n)+&
               a5*duvxl(n) + a6*dvwxl(n) + a7*dwuxl(n)+ &
               a4*dtxl(n)

   dfluxr(5,n)=0.5*a1*duuxr(n)+0.5*a2*dvvxr(n)+0.5*a3*dwwxr(n)+&
               a5*duvxr(n) + a6*dvwxr(n) + a7*dwuxr(n)+ &
               a4*dtxr(n)
enddo
!
flux(5) = 0.5*a1*uux + 0.5*a2*vvx + 0.5*a3*wwx + &
              a5*uvx +     a6*vwx +     a7*wux + &
              a4*tx
!
do n=1,5
   dfluxl(5,n)=fac*(dmul(n)*flux(5)+(mu+mut*sigma)*dfluxl(5,n))
   dfluxr(5,n)=fac*(dmur(n)*flux(5)+(mu+mut*sigma)*dfluxr(5,n))
enddo
!
n=6
dfluxl(5,n)=dmutl*sigma*fac*flux(5)
dfluxr(5,n)=dmutr*sigma*fac*flux(5)
!
flux(5) = (mu+mut*sigma)*fac*flux(5)
!
!
do n=1,6
   dfluxl(6,n)=0
   dfluxr(6,n)=0
enddo
!
return
end subroutine jac_visc
   
   
  
   
   

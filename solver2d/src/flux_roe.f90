

!------------------------------------------------------------------------------
!#     Computes total convective flux across a face using Roe's approximate
!#                 Riemann solver given left and right conserved states
!#
!#     nxyz(3)     - Three components of face normal vector. These can be
!#                   dimensional (i.e including face area) or non-dimensional.
!#                   Returned flux does not include face area.
!#     ql(5),qr(5) - Vector of conserved variables (ro, ro*u, ro*v, ro*w, Et)
!#     flux(5)     - Vector of flux across face. Not multiplied by face area.
!#     specRadius  - spectral radius (largest eigen value)
!#     gam         - Ratio of specific heats
!#
!------------------------------------------------------------------------------
!
!################# (C) Copyright Karthik Mani 2005 ############################
!
!------------------------------------------------------------------------------
subroutine flux_roe(nxyz,ql,qr,flux,specRadius,gam)
!------------------------------------------------------------------------------
implicit none
!------------------------------------------------------------------------------
real*8  :: nxyz(3)
real*8  :: ql(5),qr(5)
real*8  :: flux(5)
real*8  :: gam
real*8  :: specRadius

real*8  :: gm1
real*8  :: nxd,nyd,nzd,area
real*8  :: rol,ul,vl,wl,pl,hl,cl
real*8  :: ror,ur,vr,wr,pr,hr,cr
real*8  :: uconl,uconr
real*8  :: nx,ny,nz
real*8  :: fm(5)
real*8  :: robar,srl,srr,denom
real*8  :: ubar,vbar,wbar,hbar,uconbar,uvwbar,cbar
real*8  :: dp,dro,du,dv,dw,ducon,duvw,dc
real*8  :: eig1,eig2,eig3
real*8  :: term
real*8  :: lamL,lamR,lamB

real*8  :: fact,A,B,term1,term2,del1,del2

real*8  :: eps1,eps2,eps3
real*8  :: eig1L,eig2L,eig3L
real*8  :: eig1R,eig2R,eig3R

real*8  :: fl(5),fr(5)

!------------------------------------------------------------------------------

       gm1 = gam - 1.0

       nxd = nxyz(1)
       nyd = nxyz(2)
       nzd = nxyz(3)
        
       area = sqrt(nxd*nxd + nyd*nyd + nzd*nzd)

       nx = nxd/area
       ny = nyd/area
       nz = nzd/area

       !----> back calculate left and right primitive states

       rol = ql(1)
       ul  = ql(2)/ql(1)
       vl  = ql(3)/ql(1)
       wl  = ql(4)/ql(1)
       pl  = gm1*( ql(5) - 0.5 * rol * (ul*ul + vl*vl + wl*wl) )
       hl  = (ql(5) + pl)/rol
       cl  = sqrt(gam*pl/rol)
        
       ror = qr(1)
       ur  = qr(2)/qr(1)
       vr  = qr(3)/qr(1)
       wr  = qr(4)/qr(1)
       pr  = gm1*( qr(5) - 0.5 * ror * (ur*ur + vr*vr + wr*wr) )
       hr  = (qr(5) + pr)/ror
       cr  = sqrt(gam*pr/ror)

       !----> primitive state differences
       
       dro = ror - rol
       du  =  ur - ul
       dv  =  vr - vl
       dw  =  wr - wl
       dp  =  pr - pl
       dc  =  cr - cl

       !----> face normal velocities

       uconr = ur*nx + vr*ny + wr*nz
       uconl = ul*nx + vl*ny + wl*nz

       !----> Roe average state

       fact = sqrt(ror/rol)
       A    = 1.0 /(1.0 + fact)
       B    = fact/(1.0 + fact)

       robar = rol*fact
       ubar  = ul*A + ur*B
       vbar  = vl*A + vr*B
       wbar  = wl*A + wr*B
       hbar  = hl*A + hr*B
       cbar = gm1*(hbar - 0.5*(ubar*ubar + vbar*vbar + wbar*wbar))
       cbar = sqrt(cbar)
       uconbar = ubar*nx + vbar*ny + wbar*nz

       !----> Eigenvalues

       eig1 = dabs(uconbar)
       eig2 = dabs(uconbar + cbar)
       eig3 = dabs(uconbar - cbar)

       specRadius=eig1+cbar

       !--------------------------------------------------!
       !---------------> Roe flux <-----------------------!
       !--------------------------------------------------!

       term1 = -eig1 + 0.5*(eig2 + eig3)
       term2 = 0.5*(eig2 - eig3)
       del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar
       del2  = term1*(uconr - uconl)*robar + term2*dp/cbar

       fm(:) = eig1*( qr(:) - ql(:) )
  
       fm(1) =  fm(1) + del1
       fm(2) =  fm(2) + del1*ubar + del2*nx
       fm(3) =  fm(3) + del1*vbar + del2*ny
       fm(4) =  fm(4) + del1*wbar + del2*nz
       fm(5) =  fm(5) + del1*hbar + del2*uconbar

       !----------------------------------------------------!
       !------> native flux for left and right states <-----!
       !----------------------------------------------------!

       fl(1) = rol*uconl
       fr(1) = ror*uconr

       fl(2) = rol*uconl*ul + nx*pl
       fr(2) = ror*uconr*ur + nx*pr

       fl(3) = rol*uconl*vl + ny*pl
       fr(3) = ror*uconr*vr + ny*pr

       fl(4) = rol*uconl*wl + nz*pl
       fr(4) = ror*uconr*wr + nz*pr

       fl(5) = ( ql(5) + pl ) * uconl
       fr(5) = ( qr(5) + pr ) * uconr


       !----> Total flux       

       flux(:) = 0.5*( fl(:) + fr(:) - fm(:) )

!------------------------------------------------------------------------------
end subroutine flux_roe

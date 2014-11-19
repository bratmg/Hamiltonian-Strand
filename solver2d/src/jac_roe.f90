subroutine jac_roe2d(nxyz2,ql2,qr2,lmat2,rmat2,gam,imode)
implicit none
real*8  :: nxyz2(2)
real*8  :: ql2(4),qr2(4)
real*8  :: lmat2(4,4),rmat2(4,4)
real*8  :: gam
integer :: imode
!
real*8  :: nxyz(3)
real*8  :: ql(5),qr(5)
real*8  :: lmat(5,5),rmat(5,5)
real*8  :: area
integer  :: idiv(4)
integer :: i,j
!
idiv(1)=1
idiv(2)=2
idiv(3)=3
idiv(4)=5
!
nxyz(1)=nxyz2(1)
nxyz(2)=nxyz2(2)
nxyz(3)=0;
area=sqrt(nxyz2(1)**2+nxyz2(2)**2)
!
ql(1:3)=ql2(1:3)
ql(4)=0;
ql(5)=ql2(4);
!
qr(1:3)=qr2(1:3)
qr(4)=0;
qr(5)=qr2(4);
!
call jac_roe(nxyz,ql,qr,lmat,rmat,gam,imode)
!
! reduce and multiply by area
!
do i=1,4
 do j=1,4
  lmat2(i,j)=lmat(idiv(i),idiv(j))*area
  rmat2(i,j)=rmat(idiv(i),idiv(j))*area
 enddo
enddo
!
return
end subroutine jac_roe2d
!
!------------------------------------------------------------------------------
!#     Computes flux Jacobian matrix for subroutine "flux_roe" that can be
!#               used in Newton solvers and adjoint codes.
!#
!#     nxyz(3)             - Three components of face normal vector.
!#                           These can be dimensional (i.e including face area)
!#                           or non-dimensional.
!#                           Returned flux does not include face area.
!#     ql(5),qr(5)         - Conserved variables (ro, ro*u, ro*v, ro*w, Et)
!#     lmat(5,5),rmat(5,5) - Left and right state flux Jacobian matrices
!#     gam                 - Ratio of specific heats
!#     imode               - 0 = approximate linearization where the eigenvalues
!#                               are treated as constants
!#                           1 = exact linearization
!#
!------------------------------------------------------------------------------
!
!################# (C) Copyright Karthik Mani 2006 ############################
!
!------------------------------------------------------------------------------

subroutine jac_roe(nxyz,ql,qr,lmat,rmat,gam,imode)
!------------------------------------------------------------------------------

implicit none
!------------------------------------------------------------------------------

real*8  :: nxyz(3)
real*8 :: ql(5),qr(5)
real*8  :: lmat(5,5),rmat(5,5)
real*8  :: gam
integer :: imode


real*8  :: gm1
real*8  :: nxd,nyd,nzd,area,nx,ny,nz
real*8  :: rol,ul,vl,wl,pl,hl
real*8  :: ror,ur,vr,wr,pr,hr
real*8  :: uconl,uconr
real*8  :: ubar,vbar,wbar,hbar,uconbar,cbar,robar
real*8  :: dp,dro,du,dv,dw
real*8  :: eig1,eig2,eig3

real*8  :: fact,A,B,term1,term2,del1,del2

integer :: k,i,j
real*8  :: dro_dql(5),dro_dqr(5)
real*8  :: du_dql(5),du_dqr(5)
real*8  :: dv_dql(5),dv_dqr(5)
real*8  :: dw_dql(5),dw_dqr(5)
real*8  :: dp_dql(5),dp_dqr(5)
real*8  :: ducon_dql(5),ducon_dqr(5)
real*8  :: ddel1_dql(5),ddel1_dqr(5)
real*8  :: ddel2_dql(5),ddel2_dqr(5)

real*8  :: dq5_dql(5),dq5_dqr(5)
real*8  :: dh_dql(5),dh_dqr(5)
real*8  :: dfact_dql(5),dfact_dqr(5)
real*8  :: dA_dql(5),dA_dqr(5)
real*8  :: dB_dql(5),dB_dqr(5)
real*8  :: drobar_dql(5),dubar_dql(5),dvbar_dql(5),dwbar_dql(5)
real*8  :: drobar_dqr(5),dubar_dqr(5),dvbar_dqr(5),dwbar_dqr(5)
real*8  :: dhbar_dql(5),duconbar_dql(5),dcbar_dql(5)
real*8  :: dhbar_dqr(5),duconbar_dqr(5),dcbar_dqr(5)

real*8  :: deig1_dql(5),deig2_dql(5),deig3_dql(5)
real*8  :: deig1_dqr(5),deig2_dqr(5),deig3_dqr(5)
real*8  :: dterm1_dql(5),dterm1_dqr(5)
real*8  :: dterm2_dql(5),dterm2_dqr(5)
real*8  :: imat(5,5)
real*8  :: cl,cr,dc_dql(5),dc_dqr(5)
real*8  :: t1a,t1b,t2a,t2b,t3a,t3b
real*8  :: eps1,eps2,eps3

real*8  :: lmat1(5,5),rmat1(5,5)
!------------------------------------------------------------------------------

      gm1 = gam - 1.0

      nxd = nxyz(1)
      nyd = nxyz(2)
      nzd = nxyz(3)

      area = sqrt(nxd*nxd + nyd*nyd + nzd*nzd)

      nx = nxd/area
      ny = nyd/area
      nz = nzd/area

      !------> back calculate primitive state

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

      !-----> primitive state differences

      dro = ror - rol
      du  =  ur - ul
      dv  =  vr - vl
      dw  =  wr - wl
      dp  =  pr - pl

      !----> face normal velocities

      uconr = ur*nx + vr*ny + wr*nz
      uconl = ul*nx + vl*ny + wl*nz
!------------------------------------------------------------------------------!
!-------> linearization of left and right primitive states <-------------------!
!------------------------------------------------------------------------------!
      !---> left state

      dro_dql    = 0.0
      dro_dql(1) = 1.0

      du_dql(:) =  0.0
      du_dql(1) = -ul /rol
      du_dql(2) =  1.0/rol

      dv_dql(:) =  0.0
      dv_dql(1) = -vl /rol
      dv_dql(3) =  1.0/rol

      dw_dql(:) =  0.0
      dw_dql(1) = -wl /rol
      dw_dql(4) =  1.0/rol

      dp_dql(1) =  0.5*gm1*( ul*ul + vl*vl + wl*wl )
      dp_dql(2) = -gm1*ul
      dp_dql(3) = -gm1*vl
      dp_dql(4) = -gm1*wl
      dp_dql(5) =  gm1

      dq5_dql(:) = 0.0
      dq5_dql(5) = 1.0

      dh_dql(:) = -(ql(5) + pl)*dro_dql(:)/rol/rol + (1.0/rol)*(dq5_dql(:) + dp_dql(:))
 
      dc_dql(:) = (0.5*gam/cl)*( (1.0/rol)*dp_dql(:) - (pl/rol/rol)*dro_dql(:) )

      ducon_dql(1) = -uconl/rol
      ducon_dql(2) =  nx   /rol
      ducon_dql(3) =  ny   /rol
      ducon_dql(4) =  nz   /rol
      ducon_dql(5) =  0.0

!------------------------------------------------------------------------------!
      !---> right state

      dro_dqr    = 0.0
      dro_dqr(1) = 1.0

      du_dqr(:) =  0.0
      du_dqr(1) = -ur /ror
      du_dqr(2) =  1.0/ror

      dv_dqr(:) =  0.0
      dv_dqr(1) = -vr /ror
      dv_dqr(3) =  1.0/ror

      dw_dqr(:) =  0.0
      dw_dqr(1) = -wr /ror
      dw_dqr(4) =  1.0/ror

      dp_dqr(1) =  0.5*gm1*( ur*ur + vr*vr + wr*wr)
      dp_dqr(2) = -gm1*ur
      dp_dqr(3) = -gm1*vr
      dp_dqr(4) = -gm1*wr
      dp_dqr(5) =  gm1

      dq5_dqr(:) = 0.0
      dq5_dqr(5) = 1.0

      dh_dqr(:) = -(qr(5) + pr)*dro_dqr(:)/ror/ror + (1.0/ror)*(dq5_dqr(:) + dp_dqr(:))

      dc_dqr(:) = (0.5*gam/cr)*( (1.0/ror)*dp_dqr(:) - (pr/ror/ror)*dro_dqr(:) )

      ducon_dqr(1) = -uconr/ror
      ducon_dqr(2) =  nx   /ror
      ducon_dqr(3) =  ny   /ror
      ducon_dqr(4) =  nz   /ror
      ducon_dqr(5) =  0.0

!------------------------------------------------------------------------------!
!----------------------------> Roe average state <-----------------------------!
!------------------------------------------------------------------------------!

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

!------------------------------------------------------------------------------!
!--------------------------> Eigenvalues <-------------------------------------!
!------------------------------------------------------------------------------!

      eig1 = abs(uconbar)
      eig2 = abs(uconbar + cbar)
      eig3 = abs(uconbar - cbar)

!------------------------------------------------------------------------------!
!--------------> approximate linearization section <---------------------------!
!------------------------------------------------------------------------------!

if( imode==1 ) then
      term1 = -eig1 + 0.5*(eig2 + eig3)
      term2 = 0.5*(eig2 - eig3)
      del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar
      del2  = term1*(uconr - uconl)*robar + term2*dp/cbar

      ddel1_dql(:) = - term1*dp_dql(:)/cbar/cbar - term2*robar*ducon_dql(:)/cbar
      ddel1_dqr(:) = + term1*dp_dqr(:)/cbar/cbar + term2*robar*ducon_dqr(:)/cbar

      ddel2_dql(:) = - term1*ducon_dql(:)*robar - term2*dp_dql(:)/cbar
      ddel2_dqr(:) = + term1*ducon_dqr(:)*robar + term2*dp_dqr(:)/cbar

      goto 111
     
endif        

!------------------------------------------------------------------------------!
!-----------> linearization of Roe averaged state <----------------------------!
!------------------------------------------------------------------------------!

      dfact_dql(:) = (0.5/fact)*(-ror/rol/rol)*dro_dql(:)
      dfact_dqr(:) = (0.5/fact)*(1.0/rol)*dro_dqr(:)

      dA_dql(:) = -dfact_dql(:)/(1.0+fact)/(1.0+fact)
      dA_dqr(:) = -dfact_dqr(:)/(1.0+fact)/(1.0+fact)

      dB_dql(:) = dfact_dql(:)/(1.0 + fact)/(1.0 + fact)
      dB_dqr(:) = dfact_dqr(:)/(1.0 + fact)/(1.0 + fact)

      drobar_dql(:) = dro_dql(:)*fact + rol*dfact_dql(:)
      drobar_dqr(:) =                   rol*dfact_dqr(:)

      dubar_dql(:) = du_dql(:)*A + ul*dA_dql(:)               + ur*dB_dql(:)
      dubar_dqr(:) =               ul*dA_dqr(:) + du_dqr(:)*B + ur*dB_dqr(:)

      dvbar_dql(:) = dv_dql(:)*A + vl*dA_dql(:)               + vr*dB_dql(:)
      dvbar_dqr(:) =               vl*dA_dqr(:) + dv_dqr(:)*B + vr*dB_dqr(:)

      dwbar_dql(:) = dw_dql(:)*A + wl*dA_dql(:)               + wr*dB_dql(:)
      dwbar_dqr(:) =               wl*dA_dqr(:) + dw_dqr(:)*B + wr*dB_dqr(:)

      dhbar_dql(:) = dh_dql(:)*A + hl*dA_dql(:)               + hr*dB_dql(:)
      dhbar_dqr(:) =               hl*dA_dqr(:) + dh_dqr(:)*B + hr*dB_dqr(:)

      dcbar_dql(:) = gm1*( dhbar_dql(:) - ubar*dubar_dql(:)     &
                                        - vbar*dvbar_dql(:)     &
                                        - wbar*dwbar_dql(:) )
      dcbar_dql(:) = dcbar_dql(:)*0.5/cbar

      dcbar_dqr(:) = gm1*( dhbar_dqr(:) - ubar*dubar_dqr(:)     &
                                        - vbar*dvbar_dqr(:)     &
                                        - wbar*dwbar_dqr(:) )
      dcbar_dqr(:) = dcbar_dqr(:)*0.5/cbar 

      duconbar_dql(:) = dubar_dql(:)*nx + dvbar_dql(:)*ny + dwbar_dql(:)*nz
      duconbar_dqr(:) = dubar_dqr(:)*nx + dvbar_dqr(:)*ny + dwbar_dqr(:)*nz

!------------------------------------------------------------------------------!
!------------------> linearization of Eigenvalues <----------------------------!
!------------------------------------------------------------------------------!

      if(uconbar>=0.0) then
                deig1_dql(:) = duconbar_dql(:)
                deig1_dqr(:) = duconbar_dqr(:)
      elseif(uconbar< 0.0) then
                deig1_dql(:) = -duconbar_dql(:)
                deig1_dqr(:) = -duconbar_dqr(:)
      endif

      if( (uconbar + cbar) >= 0.0 ) then
                deig2_dql(:) = ( duconbar_dql(:) + dcbar_dql(:) )
                deig2_dqr(:) = ( duconbar_dqr(:) + dcbar_dqr(:) )
      elseif( (uconbar + cbar) < 0.0 ) then
                deig2_dql(:) = -( duconbar_dql(:) + dcbar_dql(:) )
                deig2_dqr(:) = -( duconbar_dqr(:) + dcbar_dqr(:) )
      endif

      if( (uconbar - cbar) >= 0.0 ) then
                deig3_dql(:) = ( duconbar_dql(:) - dcbar_dql(:) )
                deig3_dqr(:) = ( duconbar_dqr(:) - dcbar_dqr(:) )
      elseif( (uconbar - cbar) < 0.0 ) then
                deig3_dql(:) = -( duconbar_dql(:) - dcbar_dql(:) )
                deig3_dqr(:) = -( duconbar_dqr(:) - dcbar_dqr(:) )
      endif

!------------------------------------------------------------------------------!
      term1 = -eig1 + 0.5*(eig2 + eig3)
      term2 = 0.5*(eig2 - eig3)
      del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar
      del2  = term1*(uconr - uconl)*robar + term2*dp/cbar

      dterm1_dql(:) = -deig1_dql(:) + 0.5*( deig2_dql(:) + deig3_dql(:) )
      dterm1_dqr(:) = -deig1_dqr(:) + 0.5*( deig2_dqr(:) + deig3_dqr(:) )

      dterm2_dql(:) = 0.5*( deig2_dql(:) - deig3_dql(:) )
      dterm2_dqr(:) = 0.5*( deig2_dqr(:) - deig3_dqr(:) )

      ddel1_dql(:) = dterm1_dql(:)*dp/cbar/cbar - term1*dp_dql(:)/cbar/cbar - 2.0*term1*dp*dcbar_dql(:)/cbar/cbar/cbar
      ddel1_dql(:) = ddel1_dql(:) + dterm2_dql(:)*robar*( uconr-uconl )/cbar + term2*drobar_dql(:)*(uconr-uconl)/cbar &
                                  - term2*robar*ducon_dql(:)/cbar - dcbar_dql(:)*term2*robar*(uconr-uconl)/cbar/cbar

      ddel1_dqr(:) = dterm1_dqr(:)*dp/cbar/cbar + term1*dp_dqr(:)/cbar/cbar - 2.0*term1*dp*dcbar_dqr(:)/cbar/cbar/cbar
      ddel1_dqr(:) = ddel1_dqr(:) + dterm2_dqr(:)*robar*( uconr-uconl )/cbar + term2*drobar_dqr(:)*(uconr-uconl)/cbar &
                                  + term2*robar*ducon_dqr(:)/cbar - dcbar_dqr(:)*term2*robar*(uconr-uconl)/cbar/cbar

      ddel2_dql(:) = dterm1_dql(:)*(uconr-uconl)*robar - term1*ducon_dql(:)*robar + term1*(uconr-uconl)*drobar_dql(:)
      ddel2_dql(:) = ddel2_dql(:) + dterm2_dql(:)*dp/cbar - term2*dp_dql(:)/cbar - dcbar_dql(:)*term2*dp/cbar/cbar

      ddel2_dqr(:) = dterm1_dqr(:)*(uconr-uconl)*robar + term1*ducon_dqr(:)*robar + term1*(uconr-uconl)*drobar_dqr(:)
      ddel2_dqr(:) = ddel2_dqr(:) + dterm2_dqr(:)*dp/cbar + term2*dp_dqr(:)/cbar - dcbar_dqr(:)*term2*dp/cbar/cbar

!------------------------------------------------------------------------------!
111 continue
!------------------------------------------------------------------------------!
!-----------------------> Roe flux Jacobian <----------------------------------!
!------------------------------------------------------------------------------!


      !------------> common linearization terms

      lmat(:,:) = 0.0
      rmat(:,:) = 0.0
     
      imat(:,:) = 0.0
      do k = 1, 5
           imat(k,k) = 1.0
      enddo

      lmat(:,:) = lmat(:,:) - eig1*imat(:,:)
      rmat(:,:) = rmat(:,:) + eig1*imat(:,:) 

      lmat(1,:) = lmat(1,:) + ddel1_dql(:)
      rmat(1,:) = rmat(1,:) + ddel1_dqr(:)

      lmat(2,:) = lmat(2,:) + ddel1_dql(:)*ubar + ddel2_dql(:)*nx
      rmat(2,:) = rmat(2,:) + ddel1_dqr(:)*ubar + ddel2_dqr(:)*nx

      lmat(3,:) = lmat(3,:) + ddel1_dql(:)*vbar + ddel2_dql(:)*ny
      rmat(3,:) = rmat(3,:) + ddel1_dqr(:)*vbar + ddel2_dqr(:)*ny

      lmat(4,:) = lmat(4,:) + ddel1_dql(:)*wbar + ddel2_dql(:)*nz
      rmat(4,:) = rmat(4,:) + ddel1_dqr(:)*wbar + ddel2_dqr(:)*nz

      lmat(5,:) = lmat(5,:) + ddel1_dql(:)*hbar + ddel2_dql(:)*uconbar
      rmat(5,:) = rmat(5,:) + ddel1_dqr(:)*hbar + ddel2_dqr(:)*uconbar



      !------> additional terms for exact linearization

      if(imode/=1) then
        do i = 1, 5
          do j = 1, 5
              lmat(i,j) = lmat(i,j) + ( qr(i) - ql(i) )* deig1_dql(j)
              rmat(i,j) = rmat(i,j) + ( qr(i) - ql(i) )* deig1_dqr(j)
          enddo
        enddo


        lmat(2,:) = lmat(2,:) + del1*dubar_dql(:)
        rmat(2,:) = rmat(2,:) + del1*dubar_dqr(:)

        lmat(3,:) = lmat(3,:) + del1*dvbar_dql(:)
        rmat(3,:) = rmat(3,:) + del1*dvbar_dqr(:)

        lmat(4,:) = lmat(4,:) + del1*dwbar_dql(:)
        rmat(4,:) = rmat(4,:) + del1*dwbar_dqr(:)

        lmat(5,:) = lmat(5,:) + del1*dhbar_dql(:) + del2*duconbar_dql(:)
        rmat(5,:) = rmat(5,:) + del1*dhbar_dqr(:) + del2*duconbar_dqr(:)
      endif

!------------------------------------------------------------------------------!
      !-------------------------------------------------------!
      !-----------> Compute native flux Jacobian <------------!
      !-------------------------------------------------------!
      !-----> left state

      lmat1(1,:) = dro_dql(:)*uconl + rol*ducon_dql(:)

      lmat1(2,:) = dro_dql(:)*uconl*ul + rol*ducon_dql(:)*ul + rol*uconl*du_dql(:)
      lmat1(3,:) = dro_dql(:)*uconl*vl + rol*ducon_dql(:)*vl + rol*uconl*dv_dql(:)
      lmat1(4,:) = dro_dql(:)*uconl*wl + rol*ducon_dql(:)*wl + rol*uconl*dw_dql(:)

      lmat1(3,:) = lmat1(3,:) + ny*dp_dql(:)
      lmat1(2,:) = lmat1(2,:) + nx*dp_dql(:)
      lmat1(4,:) = lmat1(4,:) + nz*dp_dql(:)

      lmat1(5,:) = ( dq5_dql(:) + dp_dql(:) )*uconl + (ql(5) + pl)*ducon_dql(:)

      !======================================================================
      !-----> right state

      rmat1(1,:) = dro_dqr(:)*uconr + ror*ducon_dqr(:)

      rmat1(2,:) = dro_dqr(:)*uconr*ur + ror*ducon_dqr(:)*ur + ror*uconr*du_dqr(:)
      rmat1(3,:) = dro_dqr(:)*uconr*vr + ror*ducon_dqr(:)*vr + ror*uconr*dv_dqr(:)
      rmat1(4,:) = dro_dqr(:)*uconr*wr + ror*ducon_dqr(:)*wr + ror*uconr*dw_dqr(:)

      rmat1(3,:) = rmat1(3,:) + ny*dp_dqr(:)
      rmat1(2,:) = rmat1(2,:) + nx*dp_dqr(:)
      rmat1(4,:) = rmat1(4,:) + nz*dp_dqr(:)

      rmat1(5,:) = ( dq5_dqr(:) + dp_dqr(:) )*uconr + (qr(5) + pr)*ducon_dqr(:)

      !======================================================================

      lmat(:,:) = 0.5*( lmat1(:,:) - lmat(:,:) )
      rmat(:,:) = 0.5*( rmat1(:,:) - rmat(:,:) )

!      write(6,*) '---------- lmat f90 ----------------'
!      do i=1,5
!         write(6,"(5(1X,E14.8))") (lmat(i,j),j=1,5)
!      enddo
!      write(6,*) '---------- rmat f90 ----------------'
!      do i=1,5
!         write(6,"(5(1X,E14.8))") (rmat(i,j),j=1,5)
!      enddo
!------------------------------------------------------------------------------

end subroutine jac_roe

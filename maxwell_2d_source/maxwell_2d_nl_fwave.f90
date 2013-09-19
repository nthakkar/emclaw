! ===============================================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq,num_aux)
! ===============================================================================
!
!   # This version outputs f-waves.

!   # Riemann solver for the time dependent nonlinear Maxwell em equations in 1d,
!   # in this case eps and mu dependend on time and position,
!   # variable coefficients
!   #   kappa1(q,t,r)*(q1)_t + (q2)_x           = -(eps)_t*(q1)
!   #   kappa2(q,t,r)*(q2)_t + (q1)_y           = -(eps)_t*(q2)
!   #   kappa3(q,t,r)*(q3)_t + (q2)_y - (q1)_x  = -(mu)_t*(q3) 
!   #
!   # where q1=Eq, q2=E, q3=H, eps=f(x,t), and mu=g(x,t)
!   # and kappa1 = eps + 2*chi2_e_r*E1 + 3*chi3_e*E1^2
!   # and kappa2 = eps + 2*chi2_e_r*E2 + 3*chi3_e*E2^2
!   # and kappa3 = mu  + 2*chi2_m_r*H3 + 3*chi3_m*H3^2
!   #
!   # aux(1,i) = eps(i)
!   # aux(2,i) = mu(i)
!   # aux(3,i) = eps_t(i)
!   # aux(4,i) = mu_t(i)

!   # function f(x_i,t_i) gives the permittivity value at the ith cell at step t_i
!   # function g(x_i,t_i) gives the permeability value at the ith cell at step t_i
!   #    For RIP:   f(x_i,t_i)=f(x_i-v*t_i), and g(x_i,t_i)=g(x_i-v*t_i)
!   # the system assumes the em functions to be some constant value + transient part

!   # On input, ql contains the state vector at the left edge of each cell
!   #           qr contains the state vector at the right edge of each cell

!   # On output, fwave contains the waves as jumps in f,
!   #            s the speeds,
!   #
!   #            amdq = A^- Delta q,
!   #            apdq = A^+ Delta q,
!   #                   the decomposition of the flux difference minus the source term
!   #                       f(qr(i-1)) - f(ql(i)) - \psi(q,x,t)
!   #                   into leftgoing and rightgoing parts respectively.
!   #

!   # Note that the ith Riemann problem has left state qr(:,i-1)
!   #                                    and right state ql(:,i)
!   # From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    double precision :: auxl(num_aux,1-mbc:maxm+mbc)
    double precision :: auxr(num_aux,1-mbc:maxm+mbc)
    double precision :: fwave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)
    double precision :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision :: chi2(3)
    double precision :: chi3(3)
    integer :: i, mx, mbc, maxm, num_aux, meqn, mwaves, m, ixy

    double precision :: eta1i, eta1im, eta2i, eta2im, eta3i, eta3im
    double precision :: etat1i, etat1im, etat2i, etat2im, etat3i, etat3im
    double precision :: q1i, q1im, q2i, q2im, q3i, q3im
    double precision :: ci, cim, zi, zim
    double precision :: kappa1, kappa2, kappa3, eo, mo, zo, co
    double precision :: df1, df2, df3, psi1, psi2, psi3, dx, dy
    double precision :: dq1, dq2, dq3, beta1, beta2, beta3
    double precision :: vac1, vac2, vac3
    common /cparam/  dx, dy, chi2, chi3, eo, mo, co, zo, vac1, vac2, vac3

!     # split the jump in q at each interface into waves

    do 20 i = 2-mbc, mx+mbc
        eta1i   = auxl(1,i  )
        eta1im  = auxr(1,i-1)
        eta2i   = auxl(2,i  )
        eta2im  = auxr(2,i-1)
        eta3i   = auxl(3,i  )
        eta3im  = auxr(3,i-1)
        etat1i  = auxl(4,i  )
        etat1im = auxr(4,i-1)
        etat2i  = auxl(5,i  )
        etat2im = auxr(5,i-1)
        etat3i  = auxl(6,i  )
        etat3im = auxr(6,i-1)

        q1i     = ql(1,i)
        q1im    = qr(1,i-1)
        q2i     = ql(2,i)
        q2im    = qr(2,i-1)
        q3i     = ql(3,i)
        q3im    = qr(3,i-1)

!     # calculate velocity c = 1/sqrt(eps*mu) and impedance Z = sqrt(eps/mu)
        ci      = co 
        cim     = co
        zi      = zo
        zim     = zo

        psi1 = -0.5d0*(etat1i*q1i + etat1im*q1im)
        psi2 = -0.5d0*(etat2i*q2i + etat2im*q2im)
        psi3 = -0.5d0*(etat3i*q3i + etat3im*q3im)

!     # flux difference
        
        dq1 = (q1i - q1im)
        dq2 = (q2i - q2im)
        dq3 = (q3i - q3im)

        kappa1 = 0.5d0*(eta1i + eta1im + 2.d0*chi2(1)*(q1i + q1im) + 3.d0*chi3(1)*((q1i + q1im)**2))
        kappa2 = 0.5d0*(eta2i + eta2im + 2.d0*chi2(2)*(q2i + q2im) + 3.d0*chi3(2)*((q2i + q2im)**2))
        kappa3 = 0.5d0*(eta3i + eta3im + 2.d0*chi2(3)*(q3i + q3im) + 3.d0*chi3(3)*((q3i + q3im)**2))
        !   Normal & perpendicular waves
!   ------------
        if (ixy==1) then
            df2 = dq3/vac2 - dx*psi2
            df3 = dq2/vac3 - dx*psi3
            beta1 = (df3*zi - df2)/(zi+zim)
            beta2 = 0.d0
            beta3 = (df3*zim + df2)/(zi+zim)
            fwave(1,1,i) = 0.d0
            fwave(2,1,i) = beta1*(-zim)/kappa2
            fwave(3,1,i) = beta1/kappa3
            fwave(1,2,i) = 0.d0
            fwave(2,2,i) = beta3*(zi)/kappa2
            fwave(3,2,i) = beta3/kappa3
            s(1,i) = -cim
            s(2,i) = ci 
        else
            df1 = dq3/vac1 - dy*psi1
            df3 = dq1/vac3 - dy*psi3
            beta1 = -1.d0*(df1 + df3*zi)/(zi+zim)
            beta2 = 0.d0
            beta3 = (df1 - df3*zim)/(zi+zim)
            fwave(1,1,i) = beta1*(zim)/kappa1
            fwave(2,1,i) = 0.d0
            fwave(3,1,i) = beta1/kappa3
            fwave(1,2,i) = beta3*(-zi)/kappa1
            fwave(2,2,i) = 0.d0
            fwave(3,2,i) = beta3/kappa3
            s(1,i) = -cim
            s(2,i) = ci 
        endif
        
    20 END DO


!     # compute the leftgoing and rightgoing fluctuations:
!     # Note s(1,i) < 0   and   s(2,i) > 0.

    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = fwave(m,1,i)
            apdq(m,i) = fwave(m,2,i)
    220 END DO

    return
    end subroutine rpn2

    subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
    fwave,s,amdq,apdq,num_aux)
!     =====================================================

!     # This version outputs f-waves.

!     # Riemann solver for the time dependent nonlinear em equations in 1d,
!     # in this case eps and mu dependend on time and position,
!     #  variable coefficients
!     #   kappa1(q,t,x) (q1)_t +  (q2)_x  = -eps_t E
!     #   kappa2(q,t,x) (q2)_t +  (q1)_x = -mu_t H
!     # where q1=E, q2=H, eps=f(x,t), and mu=g(x,t)
!     # and kappa1 = eps + 2*chi2_e*E + 3*chi3_e*E^2
!     # and kappa2 = mu  + 2*chi2_m*H + 3*chi3_m*H^2

!     # aux(1,i) = eps(i)
!     # aux(2,i) = mu(i)
!     # aux(3,i) = eps_t(i)
!     # aux(4,i) = mu_t(i)

!     # function f(x_i,t_i) gives the permittivity value at the ith cell at step t_i
!     # function g(x_i,t_i) gives the permeability value at the ith cell at step t_i
!     #    For RIP:   f(x_i,t_i)=f(x_i-v*t_i), and g(x_i,t_i)=g(x_i-v*t_i)
!     # the system assumes the em functions to be some constant value + transient part

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell

!     # On output, fwave contains the waves as jumps in f,
!     #            s the speeds,
!     #
!     #            amdq = A^- Delta q,
!     #            apdq = A^+ Delta q,
!     #                   the decomposition of the flux difference minus the source term
!     #                       f(qr(i-1)) - f(ql(i)) - \psi(q,x,t)
!     #                   into leftgoing and rightgoing parts respectively.
!     #

!     # Note that the ith Riemann problem has left state qr(:,i-1)
!     #                                    and right state ql(:,i)
!     # From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    double precision :: auxl(num_aux,1-mbc:maxm+mbc)
    double precision :: auxr(num_aux,1-mbc:maxm+mbc)
    double precision :: fwave(mwaves,meqn,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)

    integer :: i, mx, mbc, maxm, num_aux, meqn, mwaves, m

    double precision :: epsi, epsim, mui, muim, Ei, Eim, Hi, Him, ci, cim, zi, zim
    double precision :: epsti, epstim, muti, mutim, kappa1, kappa2
    double precision :: df1, df2, psi1, psi2, b1, b2, dx, chi2_e, chi2_m, chi3_e, chi3_m

    common /cparam/  dx, chi2_e, chi2_m, chi3_e, chi3_m

!     # split the jump in q at each interface into waves

    do 20 i = 2-mbc, mx+mbc
        epsi   = auxl(1,i  )
        epsim  = auxr(1,i-1)
        mui   = auxl(2,i  )
        muim  = auxr(2,i-1)
        epsti   = auxl(3,i  )
        epstim  = auxr(3,i-1)
        muti   = auxl(4,i  )
        mutim  = auxr(4,i-1)

        Ei     = ql(1,i  )
        Eim    = qr(1,i-1)
        Hi     = ql(2,i  )
        Him    = qr(2,i-1)

!     # calculate velocity c = 1/sqrt(eps*mu) and impedance Z = sqrt(eps/mu)
        ci  = 1.d0 
        cim = 1.d0
        zi  = 1.d0
        zim = 1.d0

        psi1 = - 0.5d0*(epsti*Ei + epstim*Eim)
        psi2 = - 0.5d0*(muti*Hi  + mutim*Him )

!     # flux difference minus source term
        df1 = Hi - Him - dx*psi1
        df2 = Ei - Eim - dx*psi2

        b1 = (zi * df2 - df1) / (zim + zi)
        b2 = (zim * df2 + df1) / (zim + zi)

        kappa1 = 0.5d0*(epsi+epsim + 2*chi2_e*(Ei + Eim))
        kappa2 = 0.5d0*(mui+muim   + 2*chi2_m*(mui+muim))

!        a1 = -b1 / cim
!        a2 = b2 / ci
  
!     # Compute the waves.
    
        fwave(1,1,i) = b1 *(-zim) / kappa1
        fwave(2,1,i) = b1 / kappa2
        s(1,i) = -cim
    
        fwave(1,2,i) = b2*(zi) / kappa1
        fwave(2,2,i) = b2 / kappa2
        s(2,i) = ci


    20 END DO


!     # compute the leftgoing and rightgoing fluctuations:
!     # Note s(1,i) < 0   and   s(2,i) > 0.

    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = fwave(m,1,i)
            apdq(m,i) = fwave(m,2,i)
    220 END DO

    return
    end subroutine rp1

!     ==================================================================
    subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,num_aux,wave,s,amdq,apdq)
!     ==================================================================

!   # Riemann solver for Maxwell's  equations in 3d, with time-space 
!   # varying material properties.

!   # in this case eps and mu dependend on time and position,
!   # variable coefficients

!   # Note that although there are 6 eigenvectors per direction, two 
!   # eigenvalues are always zero and so we only need to compute 4 waves.

!   # Solve Riemann problems along one slice of data.
!   # This data is along a slice in the x-direction if ixyz=1
!   #                               the y-direction if ixyz=2.
!   #                               the z-direction if ixyz=3.

!   # On output, fwave contains the waves as jumps in f, and
!   #            s the speeds,
!   #
!   #            amdq = A^- Delta q,
!   #            apdq = A^+ Delta q,
!   #                   the decomposition of the flux difference minus the source term
!   #                       f(qr(i-1)) - f(ql(i)) - psi(q,x,t)
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

    do i = 2-mbc, mx+mbc
        eta1i   = auxl(1,i  )
        eta1im  = auxr(1,i-1)
        eta2i   = auxl(2,i  )
        eta2im  = auxr(2,i-1)
        eta3i   = auxl(3,i  )
        eta3im  = auxr(3,i-1)
        eta4i   = auxl(4,i  )
        eta4im  = auxr(4,i-1)
        eta5i   = auxl(5,i  )
        eta5im  = auxr(5,i-1)
        eta6i   = auxl(6,i  )
        eta6im  = auxr(6,i-1)
        etat1i  = auxl(7,i  )
        etat1im = auxr(7,i-1)
        etat2i  = auxl(8,i  )
        etat2im = auxr(8,i-1)
        etat3i  = auxl(9,i  )
        etat3im = auxr(9,i-1)
        etat4i  = auxl(10,i  )
        etat4im = auxr(10,i-1)
        etat5i  = auxl(11,i  )
        etat5im = auxr(11,i-1)
        etat6i  = auxl(12,i  )
        etat6im = auxr(12,i-1)

        q1i     = ql(1,i)
        q1im    = qr(1,i-1)
        q2i     = ql(2,i)
        q2im    = qr(2,i-1)
        q3i     = ql(3,i)
        q3im    = qr(3,i-1)
        q4i     = ql(4,i)
        q4im    = qr(4,i-1)
        q5i     = ql(5,i)
        q5im    = qr(5,i-1)
        q6i     = ql(6,i)
        q6im    = qr(6,i-1)

!     # calculate velocity c = 1/sqrt(eps*mu) and impedance Z = sqrt(eps/mu)
        ci      = co 
        cim     = co
        zi      = zo
        zim     = zo

        psi1 = -0.5d0*(etat1i*q1i + etat1im*q1im)
        psi2 = -0.5d0*(etat2i*q2i + etat2im*q2im)
        psi3 = -0.5d0*(etat3i*q3i + etat3im*q3im)
        psi4 = -0.5d0*(etat4i*q4i + etat4im*q4im)
        psi5 = -0.5d0*(etat5i*q5i + etat5im*q5im)
        psi6 = -0.5d0*(etat6i*q6i + etat6im*q6im)        

!     # flux difference
        
        dq1 = (q1i - q1im)
        dq2 = (q2i - q2im)
        dq3 = (q3i - q3im)
        dq4 = (q4i - q4im)
        dq5 = (q5i - q5im)
        dq6 = (q6i - q6im)

        kappa1 = 0.5d0*(eta1i + eta1im + 2.d0*chi2(1)*(q1i + q1im) + 3.d0*chi3(1)*((q1i + q1im)**2))
        kappa2 = 0.5d0*(eta2i + eta2im + 2.d0*chi2(2)*(q2i + q2im) + 3.d0*chi3(2)*((q2i + q2im)**2))
        kappa3 = 0.5d0*(eta3i + eta3im + 2.d0*chi2(3)*(q3i + q3im) + 3.d0*chi3(3)*((q3i + q3im)**2))
        kappa4 = 0.5d0*(eta4i + eta4im + 2.d0*chi2(4)*(q4i + q4im) + 3.d0*chi3(4)*((q4i + q4im)**2))
        kappa5 = 0.5d0*(eta5i + eta5im + 2.d0*chi2(5)*(q5i + q5im) + 3.d0*chi3(5)*((q5i + q5im)**2))
        kappa6 = 0.5d0*(eta6i + eta6im + 2.d0*chi2(6)*(q6i + q6im) + 3.d0*chi3(6)*((q6i + q6im)**2))


        if (ixyz == 1) then
            df2 = dq6/vac2 - dx*psi2
            df3 = dq5/vac3 - dx*psi3
            df5 = dq3/vac5 - dx*psi5
            df6 = dq2/vac6 - dx*psi6
            
            beta1 = ( df2 + df6*zi )/(zi + zim)
            beta2 = (-df3 + df5*zi )/(zi + zim)
            beta3 = (-df2 + df6*zim)/(zi + zim)
            beta4 = ( df3 + df5*zim)/(zi + zim)
            
            fwave(1,1,i)   = 0.0d0
            fwave(2,1,i)   = beta1*(zim)/kappa2
            fwave(3:5,1,i) = 0.0d0
            fwave(6,1,i)   = beta1/kappa6

            fwave(1:2,2,i) = 0.0d0
            fwave(3,2,i)   = beta2*(-zim)/kappa2
            fwave(4,2,i)   = 0.0d0
            fwave(5,2,i)   = beta2/kappa5
            fwave(6,2,i)   = 0.0d0

            fwave(1,3,i)   = 0.0d0
            fwave(2,3,i)   = beta3*(-zi)/kappa2
            fwave(3:5,3,i) = 0.0d0
            fwave(6,3,i)   = beta3/kappa6

            fwave(1:2,4,i) = 0.0d0
            fwave(3,4,i)   = beta4*(zi)/kappa2
            fwave(4,4,i)   = 0.0d0
            fwave(5,4,i)   = beta4/kappa5
            fwave(6,4,i)   = 0.0d0

            s(1:2,i) = -cim
            s(3:4,i) = ci

        else if (ixyz == 2) then
            df1 = dq6/vac1 - dy*psi1
            df3 = dq4/vac3 - dy*psi3
            df4 = dq3/vac4 - dy*psi4
            df6 = dq1/vac6 - dy*psi6

            beta1 = (-df1 + df6*zi )/(zi + zim)
            beta2 = ( df3 + df4*zi )/(zi + zim)
            beta3 = ( df1 + df6*zim)/(zi + zim)
            beta4 = (-df3 + df4*zim)/(zi + zim)

            fwave(1,1,i)   = beta1*(-zim)/kappa1
            fwave(2:5,1,i) = 0.0d0
            fwave(6,1,i)   = beta1/kappa6

            fwave(1:2,2,i) = 0.0d0
            fwave(3,2,i)   = beta2*(zim)/kappa3
            fwave(4,2,i)   = beta2/kappa4
            fwave(5:6,2,i) = 0.0d0

            fwave(1,3,i)   = beta3*(zi)/kappa1
            fwave(2:5,3,i) = 0.0d0
            fwave(6,3,i)   = beta3/kappa6

            fwave(1:2,4,i) = 0.0d0
            fwave(3,4,i)   = beta4*(-zi)/kappa3
            fwave(4,4,i)   = beta4/kappa4
            fwave(5:6,4,i) = 0.0d0

            s(1:2,i) = -cim
            s(3:4,i) = ci

        else if (ixyz == 3) then
            df1 = dq5/vac1 - dz*psi1
            df2 = dq4/vac2 - dz*psi2
            df4 = dq2/vac4 - dz*psi4
            df5 = dq1/vac5 - dz*psi5

            beta1 = ( df1 + df5*zi )/(zi + zim)
            beta2 = (-df2 + df4*zi )/(zi + zim)
            beta3 = (-df1 + df5*zim)/(zi + zim)
            beta4 = ( df2 + df4*zim)/(zi + zim)

            fwave(1,1,i)   = beta1*(zim)/kappa1
            fwave(2:4,1,i) = 0.0d0
            fwave(5,1,i)   = beta1/kappa5
            fwave(6,1,i)   = 0.0d0

            fwave(1,2,i)   = 0.0d0
            fwave(2,2,i)   = beta2*(-zim)/kappa2
            fwave(3,2,i)   = 0.0d0
            fwave(4,2,i)   = beta2/kappa4
            fwave(5:6,2,i) = 0.0d0

            fwave(1,3,i)   = beta1*(-zi)/kappa1
            fwave(2:4,3,i) = 0.0d0
            fwave(5,3,i)   = beta1/kappa5
            fwave(6,3,i)   = 0.0d0

            fwave(1,4,i)   = 0.0d0
            fwave(2,4,i)   = beta2*(zi)/kappa2
            fwave(3,4,i)   = 0.0d0
            fwave(4,4,i)   = beta2/kappa4
            fwave(5:6,4,i) = 0.0d0

            s(1:2,i) = -cim
            s(3:4,i) = ci

        endif
    
    enddo


!     # compute the leftgoing and rightgoing flux differences:
!     # Note s(i,1) < 0   and   s(i,2) > 0.

    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = s(1,i)*wave(m,1,i)
            apdq(m,i) = s(2,i)*wave(m,2,i)
    220 end do


    return

    end subroutine rpn3


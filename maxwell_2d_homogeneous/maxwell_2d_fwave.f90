!   =====================================================
    subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
    wave,s,amdq,apdq,num_aux)
!   =====================================================
!
!   # Riemann solver for the TM-EM wave equations in 2D, with time-
!   # space varying, anisotropic, material properties eps and mu
!   #
!   # Note that although there are 3 eigenvectors, the second 
!   # eigenvalue is always zero and so we only need to compute 2
!   # waves.
!   #
!   # Solves Riemann problems along one slice of data.
!   # 
!   # This version uses interleaved arrays.
!   # This version outputs f-waves.
!   #
!   # in this case eps and mu dependend on time and position,
!   #
!   #           (q1)_t - (q3/mu3)_y               = 0
!   #           (q2)_t + (q3/mu3)_x               = 0
!   #           (q3)_t + (q2/eps2)_x - (q1/eps1)_y = 0
!   #
!   # The definition of eps and mu comes from the pyclaw routine
!   #
!   # On input, 
!   # ql contains the state vector at the left edge of each cell
!   # qr contains the state vector at the right edge of each cell
!   #
!   # On output, wave contains the waves as jumps,
!   #            s the speeds,
!   #
!   #           amdq = A^- Delta q,
!   #           apdq = A^+ Delta q,
!   #               the decomposition of the flux difference
!   #                       f(qr(i-1)) - f(ql(i))
!   #               into leftgoing and rightgoing parts respectively.
!   #
!   # Note that the ith Riemann problem has left state qr(:,i-1)
!   #                                  and right state ql(:,i)
!   # From the basic clawpack routines, this routine is called with
!   # ql = qr

!   ------------
    implicit double precision (a-h,o-z)
!   ------------
    dimension auxl(num_aux,1-mbc:maxm+mbc)
    dimension auxr(num_aux,1-mbc:maxm+mbc)
    dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
    dimension    s(mwaves,1-mbc:maxm+mbc)
    dimension   ql(meqn,1-mbc:maxm+mbc)
    dimension   qr(meqn,1-mbc:maxm+mbc)
    dimension apdq(meqn,1-mbc:maxm+mbc)
    dimension amdq(meqn,1-mbc:maxm+mbc)

!   local arrays
!   ------------
    dimension dq(3)
    dimension beta(3)

    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
!   ------------


    do 20 i = 2-mbc, mx+mbc
        eta1i   = auxl(1,i)
        eta1im  = auxr(1,i-1)
        eta2i   = auxl(2,i)
        eta2im  = auxr(1,i-1)
        q1i    = ql(1,i)
        q1im   = qr(1,i-1)
        q2i    = ql(2,i)
        q2im   = qr(2,i-1)
        q3i    = ql(3,i)
        q3im   = qr(3,i-1)

!   # calculate velocity    v = 1/sqrt(eps*mu) 
!   # and impedance         Z = sqrt(eps/mu)
        vi      = 1.0 / dsqrt(eta1i * eta2i)
        vim     = 1.0 / dsqrt(eta1im * eta2im)
        zi      = dsqrt(eta1i / eta2i)
        zim     = dsqrt(eta1im / eta2im)
!   # flux difference
        
        dq(1) = q1i/eta1i - q1im/eta1im
        dq(2) = q2i/eta1i - q2im/eta1im
        dq(3) = q3i/eta2i - q3im/eta2im

!   Normal & perpendicular waves
!   ------------
        if (ixy==1) then   
            beta(1) = (-dq(3)+dq(2)*zi)/(zi+zim)
            beta(2) = 0.
            beta(3) = (dq(3)+dq(2)*zim)/(zi+zim)
            wave(1,1,i) = 0.
            wave(2,1,i) = beta(1) * (-zim)
            wave(3,1,i) = beta(1)
            wave(1,2,i) = 0.
            wave(2,2,i) = beta(3) * (zi)
            wave(3,2,i) = beta(3)
            s(1,i) = -vim
            s(2,i) = vi 
        else
            beta(1) = -(dq(3)+dq(1)*zi)/(zi+zim)
            beta(2) = 0
            beta(3) = (dq(3)-dq(1)*zim)/(zi+zim)
            wave(1,1,i) = beta(1) * (zim)
            wave(2,1,i) = 0.
            wave(3,1,i) = beta(1)
            wave(1,2,i) = beta(3) * (-zi)
            wave(2,2,i) = 0.
            wave(3,2,i) = beta(3)
            s(1,i) = -vim
            s(2,i) = vi 
            !if (beta(1).gt.1.e-10) then
            !    write(*,*) s(1,i), beta(1)
            !endif 
        endif

    20 END DO


!   # compute the leftgoing and rightgoing fluctuations:
!   # Note s(1,i) < 0   and   s(2,i) > 0.

    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = wave(m,1,i)
            apdq(m,i) = wave(m,2,i)
    220 END DO

    return
    end subroutine rpn2

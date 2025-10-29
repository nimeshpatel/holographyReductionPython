	subroutine check_converge(c,cold,ndim,maxiter,niter,
     &             itermi,ncent,ndish,nsub,err)
C
C       CHECK CONVERGENCE OF MISELL
C       !! CURENTLY USING ONLY MAX ITERATION NUMBER CHECK
C       ORIGINAL VERSION: X.Z. 11/2/95
C
        include 'holis.inc'
C
        double complex c(ndim,ndim)
        double complex cold(ndim,ndim)
C
        nsub2 = nsub*nsub
        ndish2 = ndish*ndish
C
        err = 0.
        npoint = 0
C
        do i = 1, ndim
        do j = 1, ndim
           n2 = (i-ncent)**2 + (j-ncent)**2
           am = dsqrt(dreal(c(i,j)*dconjg(c(i,j))))
           am1 = dsqrt(dreal(cold(i,j)*dconjg(cold(i,j))))
C
           if (n2.ge.nsub2.and.n2.le.ndish2.
     &     and.am.ne.0.and.am1.ne.0.) then
              phase = datan2(dimag(c(i,j)), dreal(c(i,j)))
              phase1 = datan2(dimag(cold(i,j)), dreal(cold(i,j)))
              err = err + (phase-phase1)**2
              npoint = npoint + 1
           endif
        enddo
        enddo
C
        err = dsqrt(err/npoint)
C
        if (niter.ge.maxiter) then
           itermi = 1
        endif
C
        return
        end

	subroutine change_amplitude(c,am,ndim)
C
C       SET THE C ARRAY VALUE TO ZERO FOR ZERO AMPLITUDE POINTS
C       ORIGINAL VERSION: 11/20/95 X.Z.
C
        include 'holis.inc'
C
        double complex c(ndim,ndim)
        dimension am(ndim, ndim)
C
        do i = 1, ndim
        do j = 1, ndim
           amp = dsqrt(dreal(c(i,j)*dconjg(c(i,j))))
           if (amp.ne.0.) then
              phase = datan2(dimag(c(i,j)), dreal(c(i,j)))
              c(i,j) = dcmplx(am(i,j)*dcos(phase),
     &                 am(i,j)*dsin(phase))
           else
              c(i,j) = dcmplx(0.0,0.0)
           endif
        enddo
        enddo
C
        return
        end

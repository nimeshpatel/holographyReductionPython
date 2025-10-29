	subroutine fftshift(c,ndim)
C
C       SHIFT ARRAY VALUES AFTER FFT TO CONFORM TO DFT
C       CONVENTION
C       ORIGINAL VERSION: 4/22/95 X.Z.
C
        include 'holis.inc'
C
        double complex c(ndim,ndim),cc
C
        do j = 1, ndim
        do i = 1, ndim/2
          cc = c(i,j)
          c(i,j) = c(i+ndim/2,j) 
          c(i+ndim/2,j) = cc
        enddo
        enddo
c
        do i = 1, ndim
        do j = 1, ndim/2
          cc = c(i,j)
          c(i,j) = c(i,j+ndim/2)
          c(i,j+ndim/2) = cc
        enddo
        enddo
C
        return
        end

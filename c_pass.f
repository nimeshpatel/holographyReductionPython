	subroutine c_pass(c,cold,ndim)
C
C       PASSING DATA BETWEEN TWO COMPLEX ARRAIES
C       ORIGINAL VERSION: 11/20/95 X.Z.
C
        include 'holis.inc'
C
        double complex c(ndim,ndim), cold(ndim,ndim)
C
        do i = 1, ndim
        do j = 1, ndim
           cold(i,j) = c(i,j)
        enddo
        enddo
C
        return
        end

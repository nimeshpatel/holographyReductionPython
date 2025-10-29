	subroutine initialize(ccum, ndim)
C
        include 'holis.inc'
C
C       INITIALIZING COMPLEX ARRAY
C       ORIGINAL VERSION: 11/20/95 X.Z.
C
        double complex ccum(ndim, ndim)
C
        do i = 1, ndim
        do j = 1, ndim
           ccum(i,j) = 0.
        enddo
        enddo
C
	return
	end
        

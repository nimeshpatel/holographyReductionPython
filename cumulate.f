	subroutine cumulate(c, ccum, ndim)
C
C       SUM TO ARRAY CCUM
C       ORIGINAL VERSION: 11/25/95 X.Z.
C
        include 'holis.inc'
C
        double complex c(ndim, ndim)
        double complex ccum(ndim, ndim)
C
        do i = 1, ndim
        do j = 1, ndim
           ccum(i,j) = ccum(i,j) + c(i,j)
        enddo
        enddo
C
	return
	end
        

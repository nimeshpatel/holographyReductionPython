	subroutine divide(c, ccum, ndim, nmap)
C
C       SIMPLE DIVISION
C       ORIGINAL VERSION: 11/25/95 X.Z. 
C
        include 'holis.inc'
C
        double complex c(ndim, ndim)
        double complex ccum(ndim, ndim)
C
        do i = 1, ndim
        do j = 1, ndim
           c(i,j) = ccum(i,j)/nmap
        enddo
        enddo
C
	return
	end
        

	subroutine write_resiph(am,ph,c,ndim)
C
        include 'holis.inc'
C
C       WRITE RESIDULE PHASE INTO DATA FILE
C       ORIGINAL VERSION: 3/12/93 X.Z.
C
        dimension am(ndim, ndim), ph(ndim, ndim)
        double complex c(ndim, ndim)
C
        do i=1, ndim
        do j=1, ndim
1          write(10,*) ph(i,j)
 	enddo
	enddo
C
	return
	end
        

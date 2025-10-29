	subroutine write_amph(am, ph, ndim)
C
C       WRITE AMPLITUDE AND PHASE ARRAYS INTO
C       DATA FILES
C       ORIGINAL VERSION: 4/10/93
C
        include 'holis.inc'
C
        dimension am(ndim, ndim), ph(ndim, ndim)
C
        do i=1, ndim
        do j=1, ndim
           write(7,*) am(i,j)
1          write(8,*) ph(i,j)
 	enddo
	enddo
C
	return
	end
        

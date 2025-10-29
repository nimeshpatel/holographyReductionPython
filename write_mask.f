	subroutine write_mask(ndim)
C
        include 'holis.inc'
C
c	writes out an integer mask file
c	original version 11/26/97, TK
C
        integer mask(nmax, nmax)
	common /mask/ mask
C
        do i=1, ndim
        do j=1, ndim
1          write(7,*) mask(i,j)
 	enddo
	enddo
C
	return
	end
        

	subroutine apply_mask5(ph,ndim)
C
C	MASKING THE real APERTURE DATA FILE WITH THE DATA POINTS
C       OUTSIDE THE MAIN REFLECTOR AND INSIDE THE SUBREFLECTOR
C       PUT TO ZERO
c	Also includes a simple cross pattern, tilted 7.5 deg from
c	45 deg; no shadows included. 
C
c	uses the mask pattern in maskfile
C
        include 'holis.inc'
C
        dimension ph(ndim, ndim)
C
        do i = 1, ndim
        do j = 1, ndim
	 read(10,*) mask
	 if (mask.eq.0) then
          ph(i,j) = -9999.0
	 endif
        enddo
        enddo
C
        return
        end 

	subroutine separate_ap(c, am, ph, ndim)
C
        include 'holis.inc'
C
C       SEPARTE THE AMPLITUDE AND PHASE FROM
C       THE COMPLEX APERTURE ARRAY
C       ORIGINAL VERSION: 3/12/93 X.Z.
c	modified to include -9999 blanking 11/26/97 TK
C
        dimension am(ndim, ndim), ph(ndim, ndim)
        double complex c(ndim, ndim)
	integer mask(nmax,nmax)
	common /mask/ mask
C
        do i=1, ndim
        do j=1, ndim
           am(i,j)=dsqrt(dreal(c(i,j)*dconjg(c(i,j))))
           if (am(i,j).eq.0.) then
              ph(i,j) = 0.
           elseif (c(i,j).eq.dcmplx(-9999.0,-9999.0)) then
		ph(i,j) = -9999.0
	   else	
              ph(i,j)=datan2(dimag(c(i,j)),dreal(c(i,j)))
           endif
 	enddo
	enddo
C
	return
	end
        

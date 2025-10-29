	subroutine read_amph(am, ph, c, ndim)
C
        include 'holis.inc'
C
C       READ IN AMPLITUDE AND PHASE
C       ORIGINAL VERSION: 5/18/93 X.Z.
C
        dimension am(ndim, ndim), ph(ndim, ndim)
        double complex c(ndim, ndim)
C
        read(2,*) ((am(i,j),j=1,ndim),i=1,ndim)
        read(3,*) ((ph(i,j),j=1,ndim),i=1,ndim)
C
        do i=1, ndim
        do j=1, ndim
           c(i,j)=dcmplx((am(i,j)*dcos(ph(i,j))),
     &            (am(i,j)*dsin(ph(i,j))))
 	enddo
	enddo
C
	return
	end
        

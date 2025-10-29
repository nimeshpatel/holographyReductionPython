	subroutine undefocus(c,ndim,ncent,
     &             fprim,fm,pi,alambda,df,dprim,rate)
C
C       UNDO THE DEFOCUS ON APERTURE PATTERN
C       ORIGINAL VERSION: 11/20/95
C
        include 'holis.inc'
C
        double complex c(ndim, ndim)
        double complex cj
C
        cj=dcmplx(0.,1.)
        dx = dprim/rate/ndim
C
        do i = 1, ndim
        do j = 1, ndim
           radius1= dsqrt(dfloat((i-ncent)**2+(j-ncent)**2))*dx 
           phi = -defocus1(radius1,alambda,fprim,fm,df,pi)
           c(i,j) = c(i,j)*(dcos(phi)+cj*dsin(phi))
        enddo
        enddo
C
        return
        end

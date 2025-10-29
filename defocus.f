	subroutine defocus(c,ndim,ncent,
     &             fprim,fm,pi,alambda,df,dprim,rate)
C
        include 'holis.inc'
C
C       CREATE DEFOCUSED COMPLEX APERATURE PATTERN
C       ORIGINAL VERSION: 6/20/93 X.Z.
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
           phi = defocus1(radius1,alambda,fprim,fm,df,pi)
           c(i,j) = c(i,j)*(dcos(phi)+cj*dsin(phi))
        enddo
        enddo
C
        return
        end

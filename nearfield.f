       subroutine nearfield(c,ndim,ncent,ndish,nsub,
     &            dprim,distan,rate,freq,pi,clight)
C
       include 'holis.inc'
C
C      PERFORM NEAR FILED CORRECTION
C      original version: 5/23/93
C
       double complex c(ndim,ndim)
       double complex cj
C
       freq1=freq*1.e9
       alambda=clight/freq1
       resol=dprim/rate/ndim
       wk=2*pi/alambda
       cj = dcmplx(0.,1.)
C
       do m=1,ndim
       do n=1,ndim
          x=(dfloat(n-ncent)*resol)
          y=(dfloat(m-ncent)*resol)
          radius=dsqrt(x**2 + y**2)
          phi = wk*radius**2/2/distan
          c(m,n)=c(m,n)*(dcos(phi) + cj*dsin(phi))
       enddo
       enddo
C
       return
       end

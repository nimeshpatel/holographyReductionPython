        subroutine rms1d(ph1,nmasked,rms1)
C	finds rms of the first nmasked points of the 
c	1-d array ph1
c	original version 26/11/97 TK.
C
        include 'holis.inc'
C
        dimension ph1(nmasked)
C	write(6,*)nmasked
        rms1=0.D0
        sum=0.D0
        sumsq=0.D0
        ncent=ndim/2 + 1
C
        do i=1,nmasked
	   sum=sum+ph1(i)
	   sumsq=sumsq+ph1(i)**2
20	   continue
        enddo
C
	ave=sum/nmasked
	avesq=sumsq/nmasked
	rms1=dsqrt(avesq-ave**2)
C	write(6,*)rms1
c
	return
	end

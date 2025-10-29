        subroutine phasefit(ph,ph1,ndim,ndimsq,
     &  nfit,ndish,nsub,sig,u,v,w,coef,diffp,
     &  dprim,rate,freq,fprim,fmag,pi,clight,nmasked)
C
C       LEAST SQUARE FITTING ON APERTURE PHASE LARGE
C       SCALE ERROR TERMS
C       ORIGINAL VERSION: 4:26/93 X.Z.
c----------
C	rewritten to do masking and fitting correctly, 11/22/97,TK.
c----------
        include 'holis.inc'
        external func2
        dimension ph(ndim, ndim)
        dimension diffp(ndim, ndim)
        dimension ph1(nmasked), sig(nmasked), coef(nfit)
        dimension v(nfit,nfit), u(nmasked,nfit), w(nfit)
	integer mask(nmax,nmax)
	common /mask/ mask
c----------------------------------------------------------------------
C
        ndata=nmasked
        mp=nmasked
        ma=nfit
        np=nfit
        ncent=ndim/2 +1
C
        do i = 1,nmasked
           sig(i)=1.
        enddo
C
        call svdfit1(ph1,sig,ndata,coef,ma,u,v,w,mp,np,chisq,
     &              ndim,ncent,ndish,nsub,fprim,fmag,diffp,func2,
     &              dprim,rate,freq,pi,clight) 
C
        write (6,*) 'chisq=', chisq
        write (66,*) 'chisq=', chisq
C
C MAKE PHASE CORRECTION
C
        freq1=freq*1.e9
        alambda=clight/freq1
        fm=fprim*fmag
C
        do m=1,ndim
        do n=1,ndim
           x=dreal(n-ncent)
           y=dreal(m-ncent)
           radius=dsqrt(x**2 + y**2)
           radius1=radius*dprim/rate/ndim
           if (mask(m,n).eq.1) then
              ph(m,n)=ph(m,n)-coef(1)-coef(2)*x-coef(3)*y
     &               -coef(5)* diffp(m,n)
     &               -coef(4)*4*pi/alambda*(
     &              (((radius1/2./fprim)**2)/(1+(radius1/2./fprim)**2))+
     &              (((radius1/2./fm)**2)/(1+(radius1/2./fm)**2)))
     &               -coef(6)*(x**2)
     &               -coef(7)*(y**2)
     &               -coef(8)*(xy)
           endif
C
100     continue
        enddo
        enddo
C
	return
	end
        

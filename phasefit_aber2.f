        subroutine phasefit_aber2(ph,ph_um,ph1,ndim,ndimsq,
     &  nfit,ndish,nsub,sig,u,v,w,coef,diffp,
     &  dprim,rate,freq,fprim,fmag,pi,clight,nmasked)
C
C       LEAST SQUARE FITTING ON APERTURE PHASE LARGE
C       SCALE ERROR TERMS
C       ORIGINAL VERSION: 4:26/93 X.Z.
c----------
C       rewritten to do masking and fitting correctly, 11/22/97,TK.
c----------
        include 'holis.inc'
        external func2_aber2
        dimension ph(ndim, ndim)
        dimension ph_um(ndim, ndim)
        dimension diffp(ndim, ndim)
        dimension ph1(nmasked), sig(nmasked), coef(nfit)
        dimension v(nfit,nfit), u(nmasked,nfit), w(nfit)
        integer mask(nmax,nmax)
        integer dofit(nfitmax)
        common /mask/ mask
        common /fits/ dofit
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
        call svdfit1_aber2(ph1,sig,ndata,coef,ma,u,v,w,mp,np,chisq,
     &              ndim,ncent,ndish,nsub,fprim,fmag,diffp,func2_aber2,
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
        open(20,file='phasefit.Ep',status='old')
C
        do m=1,ndim
        do n=1,ndim
           x=dfloat(n-ncent)
           y=dfloat(m-ncent)
           radius=dsqrt(x**2 + y**2)
           radius1=radius*dprim/rate/ndim
              correction=-coef(1)-coef(2)*x-coef(3)*y
C     &               -coef(5)* diffp(m,n)
C    1 line modified, to fit aztigmatism in y in place of diffraction
C    5 Jan 07; TK
     &               -coef(5)*(radius1**2)*(1.414*(x+y)/radius)**2
     &               -coef(4)*4*pi/alambda*(
     &              (((radius1/2./fprim)**2)/(1+(radius1/2./fprim)**2))+
     &              (((radius1/2./fm)**2)/(1+(radius1/2./fm)**2)))
     &               -coef(6)*(radius1**2)*(x/radius)**2
     &               -coef(7)*(radius1**3)*(x/radius)
     &               -coef(8)*(radius1**3)*(y/radius)
C    3 lines modified for astigmatism and coma, TK 19 Jan 2005
C     &               -coef(6)*(x**2)
C     &               -coef(7)*(y**2)
C     &               -coef(8)*(xy)
              ph_um(m,n)=ph_um(m,n)+correction
           
           if (mask(m,n).eq.1.or.ph(m,n).ne.-9999.0) then
              ph(m,n)=ph(m,n)+correction
           else 
           ph(m,n)=-9999.0
           endif
C
c       type *,correction
        write(20,*) correction
100     continue
        enddo
        enddo
C
        return
        end


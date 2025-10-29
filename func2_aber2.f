           subroutine func2_aber2(i,afunc,ma,ndim,ncent,ndish,nsub,
     &             fprim,fmag,diffp,dprim,rate,freq,pi,clight)
C
C          ORIGINAL VERSION: 3/15/93 X.Z.
c       modified to work with new masking and the 1-d index
c       arrays in common /indices/
C
           include 'holis.inc'
C
           dimension  afunc(ma),diffp(ndim,ndim)
        integer ii(nmaxsq),ij(nmaxsq)
        integer dofit(nfitmax)
        common /indices/ ii,ij
        common /fits/ dofit
C
           freq1=freq*1.e9
           alambda=clight/freq1      ! in m
           fm=fprim*fmag
           n=ii(i)
           m=ij(i)
           x=dfloat(m-ncent)
           y=dfloat(n-ncent)
           radius=dsqrt(x**2 + y**2)
           radius1=radius*dprim/rate/ndim
           afunc(1)=dofit(1)*1.
           afunc(2)=dofit(2)*x
           afunc(3)=dofit(3)*y
C        
C        UNIT AMOUNT OF DEFOCUS CALCULATED USING RUZE FORMULA
C
           afunc(4)=
     &              dofit(4)*(4*pi/alambda*(
     &              ( ((radius1/2/fprim)**2)/(1+
     &              (radius1/2/fprim)**2) ) +
     &              ( ((radius1/2/fm)**2)/(1+
     &              (radius1/2/fm)**2) ) 
     &              ))
c    &             x**2+y**2
C        UNIT DIFFRACTION PHASE 
C
C           afunc(5)=dofit(5)*diffp(m,n)
C           afunc(6)=dofit(6)*x**2
C           afunc(7)=dofit(7)*y**2
C           afunc(8)=dofit(8)*x*y
C       Fit for astigmatism and coma, TK 19 jan 2005
C       Fit also for astigmatism in y in place of diffraction fit, which 
C       we never do; 05 jan 2007; TK.
C           afunc(5)=dofit(5)*diffp(m,n)
           afunc(5)=dofit(5)*(radius1**2)*(1.414*(x+y)/radius)**2
           afunc(6)=dofit(6)*(radius1**2)*(x/radius)**2
           afunc(7)=dofit(7)*(radius1**3)*(x/radius)
           afunc(8)=dofit(8)*(radius1**3)*(y/radius)

C
1000    return
        end

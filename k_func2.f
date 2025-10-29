	   subroutine func2(i,afunc,ma,ndim,ncent,ndish,nsub,
     &             fprim,fmag,diffp,dprim,rate,freq,pi,clight)
C
C          ORIGINAL VERSION: 3/15/93 X.Z.
c	modified to work with new masking and the 1-d index
c	arrays in common /indices/
C
           include 'holis.inc'
C
           dimension  afunc(ma),diffp(ndim,ndim)
	integer ii(nmaxsq),ij(nmaxsq)
	common /indices/ ii,ij
C
           freq1=freq*1.e9
           alambda=clight/freq1      ! in m
	   fm=fprim*fmag
           n=ii(i)
           m=ij(i)
           x=dreal(m-ncent)
           y=dreal(n-ncent)
	   radius=dsqrt(x**2 + y**2)
	   radius1=radius*dprim/rate/ndim
           afunc(1)=1.
           afunc(2)=x
           afunc(3)=y
C        
C        UNIT AMOUNT OF DEFOCUS CALCULATED USING RUZE FORMULA
C
           afunc(4)=
     &              4*pi/alambda*(
     &              ( ((radius1/2/fprim)**2)/(1+
     &              (radius1/2/fprim)**2) ) +
     &              ( ((radius1/2/fm)**2)/(1+
     &              (radius1/2/fm)**2) ) 
     &              )
c    &             x**2+y**2
C        UNIT DIFFRACTION PHASE 
C
           afunc(5)=diffp(m,n)
           afunc(6)=x**2
           afunc(7)=y**2
           afunc(8)=x*y
C
1000    return
        end

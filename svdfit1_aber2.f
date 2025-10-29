      SUBROUTINE svdfit1_aber2(ph1,sig,ndata,coef,
     *ma,u,v,w,mp,np,chisq,ndim,ncent,ndish,
     *nsub,fprim,fmag,diffp,funcs,dprim,rate,freq,pi,clight)
C
      include 'holis.inc'
C
      dimension coef(ma),sig(ndata),ph1(ndata),u(mp,np),v(np,np),w(np)
      PARAMETER (TOL=1.e-10)
C     USES svbksb,svdcmp
C     dimension afunc(nmaxsq),b(nfitmax) This is a bug M.Saito 1999.7.1
      dimension afunc(nmaxsq),b(ndata)
      dimension diffp(ndim,ndim)
c

      do 12 i=1,ndata
       call funcs(i,afunc,ma,ndim,ncent,ndish,nsub,
     *fprim,fmag,diffp,dprim,rate,freq,pi,clight)
        tmp=1./sig(i)
        do 11 j=1,ma
          u(i,j)=afunc(j)*tmp
11      continue
        b(i)=ph1(i)*tmp
12    continue
c
      call svdcmp1(u,ndata,ma,mp,np,w,v)
      wmax=0.
      do 13 j=1,ma
        if(w(j).gt.wmax)wmax=w(j)
13    continue
      thresh=TOL*wmax
      do 14 j=1,ma
        if(w(j).lt.thresh) w(j)=0.
14    continue
      call svbksb1(u,w,v,ndata,ma,mp,np,b,coef)

      write(6,*) 'coef(1)=', coef(1), '   ! constant offset'
      write(66,*) 'coef(1)=', coef(1), '   ! constant offset'
      write(6,*) 'coef(2)=', coef(2), '   ! tilt in x direction'
      write(66,*) 'coef(2)=', coef(2), '   ! tilt in x direction'
      write(6,*) 'coef(3)=', coef(3), '   ! tilt in y direction'
      write(66,*) 'coef(3)=', coef(3), '   ! tilt in y direction'
      write(6,*) 'coef(4)=', coef(4)*1000, '   ! defocus in mm'
      write(66,*) 'coef(4)=', coef(4)*1000, '   ! defocus in mm'
      write(6,*) 'coef(5)=', coef(5), '   ! astigmatism-45 term' 
      write(66,*) 'coef(5)=', coef(5), '   ! astigmatism-45 term'
      write(6,*) 'coef(6)=', coef(6), '   ! astigmatism term' 
      write(66,*) 'coef(6)=', coef(6), '   ! astigmatism term'
      write(6,*) 'coef(7)=', coef(7), '   ! coma-x term' 
      write(66,*) 'coef(7)=', coef(7), '   ! coma-x term'
      write(6,*) 'coef(8)=', coef(8), '   ! coma-y term' 
      write(66,*) 'coef(8)=', coef(8), '   ! coma-y term'
      chisq=0.
c
      do 16 i=1,ndata
       call funcs(i,afunc,ma,ndim,ncent,ndish,nsub,
     *fprim,fmag,diffp,dprim,rate,freq,pi,clight)
        sum=0.
        do 15 j=1,ma
          sum=sum+coef(j)*afunc(j)
15      continue
        chisq=chisq+((ph1(i)-sum)/sig(i))**2
16    continue
      return
      END

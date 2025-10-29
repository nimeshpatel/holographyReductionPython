      function ran1(idum)
C
C     NUMERICAL RECIPES ROUTINE
C     CHANGED TO DOUBLE PRECISION
C
      include 'holis.inc'
C
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     *ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-EPS)
C
      integer iv(NTAB)
C
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*IQ)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
C
      return
      end

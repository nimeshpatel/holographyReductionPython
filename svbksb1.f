      SUBROUTINE svbksb1(u,w,v,m,n,mp,np,b,x)
C
C     NUMERICAL RECIPES ROUTINE
C     CONVERTED TO DOUBEL PRECISION
C
      include 'holis.inc'
C
      dimension b(mp),u(mp,np),v(np,np),w(np),x(np)
      INTEGER i,j,jj
      dimension tmp(nmaxsq)
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
c       write(6,*) u(i,j), b(j)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END

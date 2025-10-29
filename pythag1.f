      FUNCTION pythag1(a,b)
C
      include 'holis.inc'
C
C     NUMERICAL RECIPES ROUTINE
C     CHANGED TO DOUBLE PRECISION
C
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag1=absa*dsqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag1=0.
        else
          pythag1=absb*dsqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END

C
C	TWO-DIMENSIONAL,RADIX 2,DECIMATION IN FREQ. FFT ALGORITHM
C
C	BERNARD ARAMBEPOLA,ENGINEERING DEPT.,TRUMPINGTON ST.
C	CAMBRIDGE, CB2 1PZ ,ENGLAND
C
C	A       -TWO DIMENSIONAL COMPLEX INPUT SEQUENCE
C	N       -NO. OF POINTS ALONG EACH DIMENSION
C	INVERS  -.FALSE. FOR FORWARD TRANSFORM, .TRUE. FOR INVERSE
C	X       -ONE DIMENSIONAL COMPLEX ARRAY OF N ELEMENTS
C	IBIT    -ONE DIMENSIONAL INTEGER ARRAY OF N ELEMENTS
C
C
C	ON OUTPUT...
C	THE 2-D COMPLEX ARRAY 'A' CONTAINS THE 2-D DFT(OR INVERSE)
C	IN NATURAL ORDER
C
C
      subroutine fft(a, n, invers, x, ibit, pi)
C
      include 'holis.inc'
C
      double complex a(n, n), x(n) 
      dimension ibit(n)
      double complex  aa, bb
C
C
C	STEP 1 : FORMING THE COMPLEX EXPONENTIALS
C
      logical l1, invers
      n1 = n / 2
      qi = pi / n1
      n2 = n1 / 2
      n3 = (n1 + n2) + 1
      n5 = n1 + 2
      n6 = n
      n7 = n1
      x(1) = (1.0,0.0)
      x(n2 + 1) = (0.0,-1.0)
      if (invers) x(n2 + 1) = (0.0,1.0)
      x(n1 + 1) = (-1.0,0.0)
      x(n3) = (0.0,1.0)
      if (invers) x(n3) = (0.0,-1.0)
      do 9 i = 2, n2
      qqi = qi * (i - 1)
      xr = dcos(qqi)
      xi = - dsin(qqi)
      if (invers) xi = - xi
      x(i) = dcmplx(xr,xi)
      x(n7) = dcmplx(- xr,xi)
      x(n5) = dcmplx(- xr,- xi)
      x(n6) = dcmplx(xr,- xi)
      n7 = n7 - 1
      n5 = n5 + 1
c
c
c	STEP 2 : 2-D,RADIX 2,D.I.F.,F.F.T. STAGES
c
c
    9 n6 = n6 - 1
      itw = n1
      itw1 = n
      itv = 1
   30 ic = 1
      id = itw
      do 20 i = 1, itv
      ie = 1
      ig = itw
      do 23 ii = 1, itv
      k1 = 1
      do 21 i1 = ic, id
      k2 = 1
      do 22 j1 = ie, ig
      i2 = i1 + itw
c
c	2-D BUTTERFLY ADDITIONS AND SUBTRACTIONS
c
      j2 = j1 + itw
      aa = a(i1,j1) + a(i1,j2)
      a(i1,j2) = a(i1,j1) - a(i1,j2)
      bb = a(i2,j1) + a(i2,j2)
      a(i2,j2) = a(i2,j1) - a(i2,j2)
      a(i1,j1) = aa + bb
      a(i2,j1) = aa - bb
      aa = a(i1,j2) + a(i2,j2)
      a(i2,j2) = a(i1,j2) - a(i2,j2)
c
c	TWIDDLE FACTOR MULTIPLICATIONS
c
      a(i1,j2) = aa
      k3 = (k1 + k2) - 1
      if (k3 .eq. 1) goto 22
      a(i1,j2) = a(i1,j2) * x(k2)
      a(i2,j1) = a(i2,j1) * x(k1)
      a(i2,j2) = a(i2,j2) * x(k3)
   22 k2 = k2 + itv
   21 k1 = k1 + itv
      ie = ie + itw1
   23 ig = ig + itw1
      ic = ic + itw1
   20 id = id + itw1
      if (itw .eq. 1) goto 100
      itv = itv + itv
      itw1 = itw
      itw = itw / 2
c
c
c	STEP 3 : 2-D BIT REVERSAL
c
c
c	FORMING THR BIT REVERSED ARRAY 'IBIT'
c
      goto 30
  100 j = 1
      ibit(n) = n
      n3 = n - 1
      do 10 i = 1, n3
      ibit(i) = j
      n2 = n1
   11 if (n2 .ge. j) goto 10
      j = j - n2
      n2 = n2 / 2
      goto 11
c
c	BIT-REVERSING THE 2-D ARRAY
c
   10 j = j + n2
      do i = 1, n
         if (ibit(i).lt.i) goto 13
         i1 = ibit(i)
         l1 = .false.
         if (i1 .eq. i) l1 = .true.
C
         do j = 1, n
            if (l1 .and. (ibit(j) .le. j)) goto 12
            j1 = ibit(j) 
            aa = a(i,j)
            a(i,j) = a(i1,j1)
            a(i1,j1) = aa   
c
c
c       STEP 4 : SCALING OPERATIONS
c
c
12          continue 
         enddo 
13       continue
      enddo

      do 50 i = 1, n
      do 50 j = 1, n
   50 a(i,j) = a(i,j) / n
c
      return 
      end

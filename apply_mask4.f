c	subroutine apply_mask4(c,ndim,ncent,pixbylen,out_mask,ain_mask,quad_hw)
	subroutine apply_mask4(c,ndim,ncent,pixbylen,out_mask,
     &  ain_mask,quad_hw)
C
C	MASKING THE APERTURE DATA FILE WITH THE DATA POINTS
C       OUTSIDE THE MAIN REFLECTOR AND INSIDE THE SUBREFLECTOR
C       PUT TO ZERO
c	Also includes a simple cross pattern, tilted 7.5 deg from
c	45 deg; no shadows included. 
C
c	creating mask and the 1-d phase and index arrays
C	Original Version 1998, TK.
C
        include 'holis.inc'
C
        double complex c(ndim, ndim)
        integer mask(nmax,nmax)
	common /mask/  mask
C
C	write(6,*)ndim,ncent,pixbylen,out_mask,ain_mask,quad_hw

C
        do i=1, ndim
        do j=1, ndim
           mask(i,j)=1
        enddo
        enddo
       
	pi = 3.1415926535897932384626
c       	theta = (45.0-7.5)*pi/180.
       	theta = 90.0*pi/180. 
	slope1 = dtan(theta)
	slope2 = -1./slope1
	n_out = dnint(out_mask*pixbylen/2.0)
	n_in   = dnint(ain_mask*pixbylen/2.0)
	n_quad_hw = dnint(quad_hw*pixbylen)
C	write(6,*)n_out,n_in,n_quad_hw
C
        do i=1, ndim
        do j=1, ndim
	x = dfloat(j-ncent)
	y = dfloat(i-ncent)
C
c	dish and subreflector
C
           radius=dsqrt(dfloat((i-ncent)**2 + (j-ncent)**2))
           if (radius.ge.n_out.or.radius.le.n_in) then
              mask(i,j)=0
           endif
C	quadrupod legs 
	y1 = x*slope1
	y2 = x*slope2
	if ( (dabs(x).le.n_quad_hw) .or.
     &       (dabs(x-0).le.n_quad_hw) .or.
     &       (dabs(dabs(x)-0).le.n_quad_hw) .or.
     &       (dabs(dabs(y)-0).le.n_quad_hw) .or.
     &	     (dabs(y).le.n_quad_hw) ) then
		mask(i,j) = 0
	endif
c
        enddo
        enddo
C
	k=0
        do i = 1, ndim
        do j = 1, ndim
	 if (mask(i,j).eq.0) then
          c(i,j) = dcmplx(-9999.0,-9999.0)
	  k=k+1
	 endif
c	write(6,*) mask(i,j),k
        enddo
        enddo
C
        return
        end 

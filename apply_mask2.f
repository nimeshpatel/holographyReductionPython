        subroutine apply_mask2(c,ndim,ncent,ndish,nsub)
C
C	MASKING THE APERTURE DATA FILE WITH THE DATA POINTS
C       OUTSIDE THE MAIN REFLECTOR AND INSIDE THE SUBREFLECTOR
C       PUT TO ZERO
C       MAY ADD STRUTS MODEL LATER
C       ORIGINAL VERSION: 11/23/95 X.Z.
C
        include 'holis.inc'
C
        double complex c(ndim, ndim)
C
        do i=1, ndim
        do j=1, ndim
C           radius=dsqrt(dreal((i-ncent)**2 + (j-ncent)**2) )
           radius=dsqrt(dfloat((i-ncent)**2 + (j-ncent)**2) )
C           radius=dsqrt(dreal((i-ncent)**2 + (j-ncent)**2) )
           radius=dsqrt(dfloat((i-ncent)**2 + (j-ncent)**2) )
           if (radius.gt.ndish.or.radius.lt.nsub) then
              c(i,j)=dcmplx(0.,0.)
           endif
        enddo
        enddo
C
        return
        end 

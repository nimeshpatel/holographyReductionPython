	subroutine create_guess(c,cold,ndim,taper,radphase,
     &             ncent,ndish,nsub,idum)
C
C       CREATES THE INITIAL GUESS OF THE COMPLEX APERTURE
C       PATTERN FOR THE MISELL ALGORITHM, USING THE SPECIFIED
C       EDGE TAPER AND RAMDOM PHASE FLUCTUATION
C
C       HISTORY:
C       ORIGINAL VERSION: 11/28/95 X.Z.
C
        include 'holis.inc'
C
        double complex c(ndim,ndim), cold(ndim,ndim),cphi
C
        aloge = 0.4342944
        ndish2 = ndish*ndish
        nsub2 = nsub*nsub
        if (taper.ne.0.) then
            sig2 = 10.*ndish2*aloge/taper
        endif
C
        do i = 1, ndim
        do j = 1, ndim
C
           c(i,j) = dcmplx(0.0, 0.0)
C
           nx = j - ncent
           ny = i - ncent
           n2 = nx**2 + ny**2
C
           if (n2.lt.ndish2.and.n2.gt.nsub2) then
              ph = radphase*gasdev(idum)
              cphi = dcmplx(dcos(phi), dsin(phi))
              if (taper.ne.0.) then
                 c(i,j) = dexp(- n2/2./sig2)*cphi
              else
                 c(i,j) = cphi
              endif
           endif
C
           cold(i,j) = c(i,j)
C 
        enddo
        enddo
C
        return
        end

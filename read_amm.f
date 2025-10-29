	subroutine read_amm(am1,am2,am3,am4,ndim,namp)
C
C       READ FOUR DEFOCUSED AMPLITUDE FILES
C       ORIGINAL VERSION: 11/22/95 X.Z.
C
        include 'holis.inc'
C
        dimension am1(ndim, ndim), am2(ndim, ndim)
        dimension am3(ndim, ndim), am4(ndim, ndim)
C
        read(2,*) ((am1(i,j),j=1,ndim),i=1,ndim)
        read(3,*) ((am2(i,j),j=1,ndim),i=1,ndim)
C
        if (namp.eq.3.or.namp.eq.4) then
           read(72,*)((am3(i,j),j=1,ndim),i=1,ndim)
        endif
C
        if (namp.eq.4) then
           read(73,*)((am3(i,j),j=1,ndim),i=1,ndim)
        endif
C
	return
	end
        

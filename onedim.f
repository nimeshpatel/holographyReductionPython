	subroutine onedim(ph, ph1, ndim, nmasked)
C
        include 'holis.inc'
C	construct a 1-d array ph1 from ph including only the unmasked pixels.
c	also generate the arraus ii and ij containing the indices i and j in ph
c	of the elements in ph1 
c	original version 11/26/97 TK
C
        dimension ph(ndim, ndim),ph1(nmasked)
	integer mask(nmax,nmax), ii(nmaxsq), ij(nmaxsq)
	common /mask/ mask
	common /indices/ ii, ij
C
	k = 0
        do i=1,ndim
          do j=1,ndim
	   if (ph(i,j).eq.-9999.0) mask(i,j)=0 
            if (mask(i,j).eq.1) then
              k=k+1
              ph1(k) = ph(i,j)
              ii(k) = i
              ij(k) = j
            endif
          enddo
        enddo
        nmasked=k
C
	return
	end
        

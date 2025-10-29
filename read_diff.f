	subroutine read_diff(diffp,ndim)
C
C       READING IN DIFFRACATION PATTERN
C       ORIGINAL VERSION: 4/25/95 X.Z.
C
        include 'holis.inc'
C
        dimension diffp(ndim,ndim)
C
        do i=1, ndim
        do j=1, ndim
           read(9,*) diffp(i,j)
	enddo
	enddo
c	offset the data to be centered on 33,33 instead of
c	32,32  and write it into a file. 12/10/97 TK, Masao
	do i=ndim-1,1,-1
	do j=ndim-1,1,-1
		diffp(i+1,j+1)=diffp(i,j)
 	enddo
	enddo
	do i = 1,ndim
	diffp(i,1)=0.0
	diffp(1,i)=0.0
	enddo
        do i=1, ndim
        do j=1, ndim
	write(72,*) diffp(i,j)
	enddo
	enddo
C
	return
	end
        

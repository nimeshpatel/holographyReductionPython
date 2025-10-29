	subroutine read_tk(ph, ninp)
C
        include 'holis.inc'
C 
C       READING IN RAW AMP AND PHASE
C       ORIGINAL VERSION: 7/16/93 X.Z.
C
        dimension ph(ninp, ninp)
C
        read(99,*) ((ph(i,j),j=1,ninp),i=1,ninp)
C
	return
	end
        

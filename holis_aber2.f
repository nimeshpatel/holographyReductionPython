	program holis
C
C       CALCULATES THE APERTURE FIELD OF AN ANTENNA FROM THE
C       MEASURED RADIATION PATTERNS, USING DIFFERENT TYPES OF
C       "HOLOGRAPHY" TECHNIQUES.
C       
C       THREE TYPES OF ALGORITHMS ARE IMPLEMENTED. THESE ARE:
C
C       1) STANDARD PHASE-COHERENT HOLOGRAPHY, WHICH REQUIRES 1 FAR-FIELD
C          AMPLITUDE MAP AND 1 FAR-FIELD PHASE MAP.
C       2) MISELL'S PHASE-RETRIEVAL ALGORITHM, WHICH REQUIRES 2          
C          FAR-FIELD AMPLITUDE MAPS OBTAINED AT DIFFERENT FOCUS
C          SETTINGS.
C       3) GLOBAL CHI-SQUARED FITTING APPROACH. IT CAN TAKE N<5 FAR
C          FIELD AMPLITUDE MAPS, AS WELL AS N PHASE MAPS, OBTAINED
C          AT DIFFERENT FOCUS SETTINGS. WEIGHT CAN BE ADJUSTED
C          BETWEEN THE AMPLITUDE AND PHASE DATA USED. IN PARTICULAR,
C          ZERO WEIGHT TO THE PHASE CORRESPONDS TO USING THE
C          DEFOCUSED AMPLITUDE MAPS ONLY FOR THE FITTING.
C
C       RADIATION PATTERNS OBTAINED AT THE NEAR FIELD CAN BE ACCOMODATED
C       
C       APERTURE MASKING AND DIFFRACTION PATTERNS ARE INCORPORATED
C       INTO THE INITIAL GUESS FOR THE PHASE-RETRIEVAL AND FITTING
C       APPROACHES.
C
C       DIFFRACTION IS SUBTRACTED FROM THE FINAL APERTURE FIELD 
C       DISTRIBUTIONS.
C
C       MASKING IS APPLIED WHEN ESTIMATING THE SURFACE RMS
C       
C       SIGN CONVENSION: THE BEAM DATA IS READ IN AS INCREASING
C       AZIMUTH VALUES AT A SERIES OF (INCREASING) ELEVATIONS. THE DATA 
C       FILES ARE OF DIMENSION (NDIM,NDIM), WHERE NDIM IS OF POWERS OF 2. 
C       THE MAXIMUM NMAX IS SPECIFIED IN A SEPARATE INCLUDE FILE 'HOLIS.INC'. 
C       THE CENTER OF THE MAP IS ALWAYS AT (NDIM/2+1, NDIM/2+1).
C
C       POSITIVE DEFOCUS DENOTES THE MOVEMENT OF THE SECONDARY AWAY
C       FROM THE MAIN REFLECTOR
C      
C       DOUBLE PRECISION CALCULATION THROUGHOUT
C
C       HISTORY:
C	ORIGINAL VERSION: XIAOLEI ZHANG 6/25/93
C       LAST REVISION: 11/20/96 X.Z.
c
c	Modified to include nearfield and defocus correction
c	flags/parameters and to read in an unwrapped file. 01/Mar/97 tk
C
C-------------------------------------------------------------------------------
C
        include 'holis.inc'
C     
C       COPY OF THE CONTENTS OF THE INCLUDE FILE GIVEN BELOW
C 
C	PARAMETER (NMAX=256)
C	PARAMETER (NMAXSQ=NMAX*NMAX)
C	PARAMETER (NFITMAX=5)
C
C	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C
        double complex c(nmax,nmax)                  !FFT DATA ARRAY 
        double complex ccum(nmax,nmax)                
        double complex cmask(nmax,nmax)              !MASKING ARRAY 
        double complex cold(nmax,nmax)               !OLD C DATA ARRAY 
        double complex x1(nmax)                      !FFT WORKING ARRAY
        dimension am(nmax,nmax), ph(nmax,nmax)       !AMPLITUDE AND PHASE ARRAYS
	dimension am_um(nmax,nmax),ph_um(nmax,nmax)  !unmasked amp & phase
        dimension am1(nmax,nmax), am2(nmax,nmax)     !MISELL AMPLITUDE ARRAYS
        dimension am3(nmax,nmax), am4(nmax,nmax)     !MISELL AMPLITUDE ARRAYS
        dimension diffp(nmax,nmax)                   !UNIT DIFFRACTION PHASE MAP
        dimension ph1(nmaxsq), sig(nmaxsq)           !1d phase and error arrays
        dimension u(nmaxsq,nfitmax)                  !U, V, W: WORKING ARRAY FOR
        dimension v(nfitmax,nfitmax)                 ! SVD
        dimension w(nfitmax)
        dimension coef(nfitmax)                      !COEFFICIENTS OF PHASEFIT
        integer ibit(nmax)                           !INTEGER ARRAY FOR FFT
c-------11/22/97, TK.
	integer mask(nmax,nmax)		             ! mask pattern
	integer ii(nmaxsq), ij(nmaxsq)		     ! valid indices after masking
	integer dofit(nfitmax)	                     ! which functions to fit
c-------
        character inamp*15, inphase*15               !INPUT AMP&PHASE FILENAMES
        character inamp1*15, inamp2*15               !INPUT AMP FILES FOR MISELL
        character inamp3*15, inamp4*15               ! AND MAXIMUM ENTROPY
        character outamp*15, outphase*15             !OUTPUT AMP & PHASE FILES
        character resiphase*15                       !RESIDUAL PHASE AFTER FIT
        character resiph_um*15                       !umRESIDUAL PHASE AFTER FIT
        character fph_um*15                       !umRESIDUAL PHASE AFTER FIT
        character fam_um*15                       !umRESIDUAL PHASE AFTER FIT
        character difffile*15                        !UNIT DIFFRACTION FILE NAME
	character maskfile*15			     !maskfile name
        logical   inverse                            !LOGICAL VARIABLE FOR FFT
c
	common /mask/ mask			     ! integer mask array
	common /indices/ ii,ij			     ! 1-d indices arrays
	common /fits/ dofit
C
        pi=3.1415926535897932384626
        clight=2.99792e8
C
C------------------------------------------------------------------------------
C
C	DEFINITION OF SOME OF THE PARAMETERS USED IN THE PROGRAM: 
C
C       X1          COMPLEX MDIM BY MDIM WORKING ARRAY FOR FFT
C       IBIT        INTEGER MDIM WORKING ARRAY FOR FFT
C	C           COMPLEX MDIM BY MDIM ARRAY FOR FAR-FIELD AND APERTURE
C                   FIELD
c
C------------------------------------------------------------------------------
C
C       THE SIGN CONVENSION
C
C       WHICH DISTINGUISHES THE FORWARD AND INVERSE FFT IS OPPOSITE TO
C       WITH THAT USED IN THE NUMERICAL RECIPIES, PP 381, AS WELL AS
C       WITH D. MORRIS, IEEE/AP-33 1985, PP 749., BUT AGREES WITH
C       THAT USED IN PACKAGE "MATLAB". UNDER THE CURRENT CONVENSION, 
C       FROM APERTUE TO FAR FIELD IS INVERSE FFT, FROM FAR FIELD TO 
C       APERTURE IS FORWARD FFT
C
C------------------------------------------------------------------------------
C    
C       OPEN LOG FILE
C
        open(66, status='unknown', file='holis.log')
C
C       CHOOSE HOLOGRAPHY SOLUTION ALGORITHMS
C
 4      write(6,*) 'Which algorithm you choose to use(1 for phase coherent, 
     &  2 for misell phase recovery, 3 for global fitting) ?'  
        read(5,*) i
C
        if (i.eq.1) then
           write (6,*) 'Use with-phase holography algorithm'
           write (6,*) '(You must first fill out "withphase.prm")'
           write (66,*) 'Use with-phase holography algorithm'
           write (66,*) 
           goto 1
	elseif (i.eq.2) then
           write (6,*) 'Use Misell phase retrieval algorithm' 
           write (6,*) '(You must first fill out "Misell.prm")'
           write (66,*) 'Use Misell phase retrieval algorithm' 
           write (66,*) 
           goto 2
        elseif (i.eq.3) then
           write (6,*) 'Use global fitting holography algorithm'
           write (6,*) '(You must first fill out "fitting.prm")'
           write (66,*) 'Use global fitting holography algorithm'
           write (66,*) 
           goto 3
        else
           goto 4
        endif
C
C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C
C       IMPLEMENTATION OF PHASE COHERENT HOLOGRAPHY 
C
 1      open (51, status='old', file='withphase_aber.prm')
C
C       THIS SECTION PERFORMS THE INVERSE FAST FOURIER TRANSFORM (FFT)
C       TO OBTAIN THE COMPLEX APERTURE FIELD DISTRIBUTION FROM THE
C       INPUT COMPLEX FAR-FIELD DISTRIBUTION. 
C
C       READING IN THE PRRAMETERS FOR THE TELESCOPE AND THE ALGORITHM 
C
	call read_prm1(inamp,inphase,outamp,outphase,resiphase,
     &                ndim,rate,distan,freq,samp_itvl,
     &                dprim,dsec,fprim,fmag,difffile,
     &                nfit,ncent,ndish,nsub,ndimsq,near_f,defoc,
     &                out_mask,ain_mask,quad_hw,maskfile,resiph_um,
     &                fam_um,fph_um) 
C       
        freq1 = freq*1.e9
        alambda = clight/freq1
	pixbylen = 2*ndish/dprim
C
C       READING IN FAR FIELD AMPLITUDE AND PHASE, TRANFER BACK 
C       FROM THE SUBROUTINE AS COMPLEX ARRAY C.
C       THE DATA WERE READ IN AS INCREASING AZIMUTH (I) AT A
C       SERIES OF (INCREASING) ELEVATIONS (J).
C
        open(2,file=inamp,status='old')
        open(3,file=inphase,status='old')
C
        call read_amph(am, ph, c, ndim)
C
C       DO ORIGIN SHIFT BEFORE FFT
C
        call fftshift(c,ndim)       !IF THE FAR FIELD PATTERN IS PEAKED IN
C                                   !THE CENTER OF THE MAP, I.E. MEASURED
C
C       FORWARD FOURIER TRANSFORM FROM FAR FIELD TO APERTURE 
C
	inverse=.false.              !forward fft to get aperture field
C
 	call fft(c, ndim, inverse, x1, ibit, pi)
C
C       DO ORIGIN SHIFT AFTER FFT
C
        call fftshift(c,ndim)       !ONLY IF THE FAR FIELD PATTERN IS MEASURED
C                                   !DO NOT DO THIS IF FAR FIELD IS FROM FFT
C                                   !OF APERTURE FIELD !!!
C
C       PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT 
	if(near_f.eq.1) then
C
        call nearfield(c,ndim,ncent,ndish,nsub,
     &                 dprim,distan,rate,freq,pi,clight)
	endif
C
C       TAKE OUT A ROUGH DEFOCUS (USUALLY DO NOT NEED TO DO THIS) 
C
	if (defoc.ne.0.0) then
         df = defoc 
         fm = fprim*fmag
         call undefocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,df,dprim,rate)
	endif
C
        goto 10000
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C
C       IMPLEMENTATION OF MISELL'S PHASE RECOVERY HOLOGRAPHY
C
 2      open (52, status='old', file='misell.prm')
C
C       READING IN THE PRRAMETERS FOR THE TELESCOPE AND THE ALGORITHM     
C
        call read_prm2(namp,inamp1,inamp2,inamp3,inamp4,
     &                defo1,defo2,defo3,defo4,
     &                outamp,outphase,resiphase,
     &                ndim,rate,distan,freq,samp_itvl,
     &                dprim,dsec,fprim,fmag,difffile,
     &                nfit,maxiter,
     &                ncent,ndish,nsub,ndimsq,
     &                taper,radphase,nmap,idum)
C
        freq1 = freq*1.e9
        alambda = clight/freq1
C
C       READING IN THE TWO FAR FIELD AMPLITUDE MAPS
C       THE DATA WERE READ IN AS INCREASING AZIMUTH (I) AT A            
C       SERIES OF (INCREASING) ELEVATIONS (J).
C
        open(2,file=inamp1,status='old')
        open(3,file=inamp2,status='old')
C
        if (namp.eq.3.or.namp.eq.4) then
           open(72,file=inamp3,status='old')
        endif
C
        if (namp.eq.4) then
           open(73,file=inamp4,status='old')
        endif
C
        call read_amm(am1,am2,am3,am4,ndim,namp)
C
        niter = 0                   !NUMBER OF ITERATIONS
        itermi = 0                  !INDICATOR FOR TERMINATION OF ITERATIONS
        imap = 0                    !INDICATOR OF NUMBER IF MAPS MADE WITH
C                                   ! DIFFERENT INITIAL SEEDS
        call initialize(ccum, ndim)
C
C       ITERATION ON MAPS WITH DIFFERENT INITIAL SEEDS
C
8888    imap = imap + 1
C
C       CREATE THE INITIAL GUESS OF THE COMPLEX APERTURE PATTERN
C
        call create_guess(c,cold,ndim,taper,radphase,
     &       ncent,ndish,nsub,idum)
C
        fm = fprim*fmag
C
C       MISELL ITERATION FOR EACH INDIVIDUAL MAP STARTS HERE
C
C       DEFOCUS THE FIRST APERTURE MAP
C
999     call defocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,defo1,dprim,rate)
C
C       PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT
C
        call nearfield(c,ndim,ncent,ndish,nsub,
     &                 dprim,-distan,rate,freq,pi,clight)
C
C       PROPAGATE INTO FAR-FIELD: WITH THE FIRST DECOCUS VALUE
C
        call fftshift(c,ndim)       
        inverse=.true.              !INVERSE FFT TO GET FAR FIELD
        call fft(c, ndim, inverse, x1, ibit, pi)
        call fftshift(c,ndim)     
C
C       CHANGE AMPLITUDE TO THE FIRST MEASURED AMPLITAUDE
C
        call change_amplitude(c,am1,ndim)
C
C       PROPAGATE BACK ONTO APERTURE: WITH THE FIRST REPLACED AMPLITUDE 
C
        call fftshift(c,ndim)       
        inverse=.false.              !FORWARD FFT TO GET APERTURE FIELD
        call fft(c, ndim, inverse, x1, ibit, pi)
        call fftshift(c,ndim)
C
C       UNDEFOCUS THE FIRST APERTURE PATTERN
C
        call undefocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,defo1,dprim,rate)
C
C       PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT 
C 
        call nearfield(c,ndim,ncent,ndish,nsub, 
     &                 dprim,distan,rate,freq,pi,clight) 
C 
C       DEFOCUS THE APERTURE FIELD TO THE SECOND DEFOCUS VALUE
C
        call defocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,defo2,dprim,rate)
C
C       PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT  
C  
        call nearfield(c,ndim,ncent,ndish,nsub,  
     &                 dprim,-distan,rate,freq,pi,clight)  
C  
C       PROPAGATE INTO FAR-FIELD: WITH THE SECOND DEFOCUS VALUE 
C
        call fftshift(c,ndim)
        inverse=.true.              !INVERSE FFT TO GET FAR FIELD
        call fft(c, ndim, inverse, x1, ibit, pi)
        call fftshift(c,ndim)
C
C       CHANGE AMPLITUDE TO THE SECOND MEASURED AMPLITAUDE 
C
        call change_amplitude(c,am2,ndim)
C
C       PROPAGATE BACK ONTO APERTURE: WITH THE SECOND REPLACED AMPLITUDE
C
        call fftshift(c,ndim)
        inverse=.false.              !FORWARD FFT TO GET APERTURE FIELD
        call fft(c, ndim, inverse, x1, ibit, pi)
        call fftshift(c,ndim)
C
C       UNDEFOCUS THE SECOND APERTURE FIELD 
C
        call undefocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,defo2,dprim,rate)
C
C       PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT   
C   
        call nearfield(c,ndim,ncent,ndish,nsub,   
     &                 dprim,distan,rate,freq,pi,clight)   
C 
        if (namp.eq.3.or.namp.eq.4) then
C
C          DEFOCUS THE APERTURE FIELD TO THE THIRD DEFOCUS VALUE
C
           call defocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,defo3,dprim,rate)
C
C          PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT   
C    
           call nearfield(c,ndim,ncent,ndish,nsub,    
     &                    dprim,-distan,rate,freq,pi,clight)    
C 
C          PROPAGATE INTO FAR-FIELD: WITH THE THIRD DEFOCUS VALUE
C
           call fftshift(c,ndim)
           inverse=.true.              !INVERSE FFT TO GET FAR FIELD
           call fft(c, ndim, inverse, x1, ibit, pi)
           call fftshift(c,ndim)
C
C          CHANGE AMPLITUDE TO THE THIRD MEASURED AMPLITUDE
C
           call change_amplitude(c,am3,ndim)
C
C          PROPAGATE BACK ONTO APERTURE: WITH THE THIRD REPLACED AMPLITUDE
C
           call fftshift(c,ndim)
           inverse=.false.              !FORWARD FFT TO GET APERTURE FIELD
           call fft(c, ndim, inverse, x1, ibit, pi)
           call fftshift(c,ndim)
C
C          UNDEFOCUS THE THIRD APERTURE FIELD
C
           call undefocus(c,ndim,ncent,fprim,
     &                 fm,pi,alambda,defo3,dprim,rate)
C
C          PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT    
C       
           call nearfield(c,ndim,ncent,ndish,nsub,     
     &                    dprim,distan,rate,freq,pi,clight)     
C 
           if (namp.eq.4) then 
C 
C             DEFOCUS THE APERTURE FIELD TO THE FOURTH DEFOCUS VALUE 
C 
              call defocus(c,ndim,ncent,fprim, 
     &                 fm,pi,alambda,defo4,dprim,rate) 
C 
C             PERFORM SECOND ORDER CORRECTION FOR NEAR FIELD EFFECT     
C        
              call nearfield(c,ndim,ncent,ndish,nsub,      
     &                       dprim,-distan,rate,freq,pi,clight)      
C  
C             PROPAGATE INTO FAR-FIELD: WITH THE SECOND DEFOCUS VALUE 
C 
              call fftshift(c,ndim) 
              inverse=.true.              !INVERSE FFT TO GET FAR FIELD
              call fft(c, ndim, inverse, x1, ibit, pi) 
              call fftshift(c,ndim) 
C 
C             CHANGE AMPLITUDE TO THE FOURTH MEASURED AMPLITAUDE 
C 
              call change_amplitude(c,am4,ndim) 
C 
C             PROPAGATE BACK ONTO APERTURE: WITH THE FOURTH REPLACED AMPLITUDE 
C 
              call fftshift(c,ndim) 
              inverse=.false.              !FORWARD FFT TO GET APERTURE FIELD 
              call fft(c, ndim, inverse, x1, ibit, pi) 
              call fftshift(c,ndim) 
C 
C             UNDEFOCUS THE FOURTH APERTURE FIELD 
C 
              call undefocus(c,ndim,ncent,fprim, 
     &                 fm,pi,alambda,defo4,dprim,rate) 
C 
              call nearfield(c,ndim,ncent,ndish,nsub, 
     &                       dprim,distan,rate,freq,pi,clight)
C
           endif
        endif
C
C       COMPARE AND SEE IF THE MAXIMUM
C       NUMBER OF ITERARTIONS HAVE BEEN REACHED
C
        call apply_mask2(c,ndim,ncent,ndish,nsub)           !ONLY DISH AND SUB
c       call apply_mask3(c,cmask,ndim,ncent,ndish,nsub,
c    &       dprim,rate)
C
        niter = niter + 1
C
        call check_converge(c,cold,ndim,maxiter,
     &       niter,itermi,ncent,ndish,nsub,err)
C
        write(6,*) ' '
        write(66,*) ' '
        write(6,*) 'imap =', imap
        write(66,*) 'imap =', imap
        write(6,*) 'Misell algorithm performed',niter,' iterations'
        write(6,*) ' rms convergence error is', err,' radian'
        write(66,*) 'Misell algorithm performed',niter,' iterations'
        write(66,*) ' rms convergence error is', err,' radian'
C
        call c_pass(c,cold,ndim)
C
        if ( itermi.ne.1) goto 999
C
        call cumulate(c, ccum, ndim)
C
        idum = idum + 95
        niter = 0
        itermi = 0
        if (imap.lt.nmap) goto 8888
C
        call divide(c, ccum, ndim, nmap)
        goto 10000
C------------------------------------------------------------------------
C
C       IMPLEMENTATION OF GLOBAL FITTING HOLOGRAPHY
C
 3      open (53, status='old', file='fitting.prm')
        write(6,*) 'this algorithm not implemented yet'
        stop
C
C------------------------------------------------------------------------
C
C
10000  continue    ! THE FOLLOWING IS COMMON TO ALL THREE ALRORITHMS
C
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C
        call separate_ap(c,am,ph,ndim)
	do i = 1,nmax
	do j = 1,nmax
	ph_um(i,j)=ph(i,j)
	am_um(i,j)=am(i,j)
	enddo
	enddo
C
	open(10,file=fam_um,status='unknown')
	call write_resiph(am,am_um,c,ndim)
	close(10)
	open(10,file=fph_um,status='unknown')
	call write_resiph(am,ph_um,c,ndim)
	close(10)
C       SEPARATE THE COMPLEX MATRIX INTO THE REAL AMP AND PHASE MATRICES
C
C       APPLYING THE MASK ON THR APERTURE FIELD
C
c
c	a new array mask containing the mask pattern included
c	11/2297, TK. apply_mask4 does also masking for the quadrupod
c	legs in a simple way, by masking an appropriately tilted
c	cross. The array mask containing the mask pattern of 1s and 0s
c	is in common /mask/; the masked pixels are replaced by
c	(-9999.0,-9999.0) 
c 
        call apply_mask4(c,ndim,ncent,pixbylen,out_mask,
     &                   ain_mask,quad_hw)
	open(7, file=maskfile,status='unknown')
	call write_mask(ndim)
	close(7)
c       call apply_mask3(c,cmask,ndim,ncent,ndish,nsub, !INCLUDE STRATS
c    &       dprim,rate)
C
c	the masked pixels have -9999 in both phase and amp.
c
        call separate_ap(c,am,ph,ndim)
C
C       WRITE OUT SOLVED APERTURE AMPLITUDE AND PHASE
C
        open(7, file=outamp, status='unknown')
        open(8, file=outphase, status='unknown')
C
C
C       APPLYING THE MASK ON THE APERTURE PHASE FILE
C
c       call apply_mask1(ph,ph1,ndim,ncent,ndish,nsub)    !DISH AND SUB
C
        call write_amph(am, ph, ndim)
	close(7)
	close(8)
C
c modified to read an unwrapped phase file 10/Jan/97, tks
c
	write (6,*) 'do you want to read in an unwrapped
     &  file(tk.dat) (1/0)?'  
	read(5,*) in_unwrapped 
	if (in_unwrapped .eq. 1) then
        open (99, file='tk.dat', status='unknown')
        call read_tk(ph, ndim)
	endif
C
C       CALCULATE RMS BEFORE PHASEFIT
C
c       the call to onedim included to do masking and fitting properly; 
c	convert the 2-D phase data array into a 1-D array and genberate 
c	the 1-D arrays ii and ij in common /indices/ containing the 
c	corresponding i and j values. nmasked is the no of valid data 
c	points after masking.
c       11/2297, TK
c
	call onedim(ph,ph1,ndim,nmasked)
c
        call rms1d(ph1,nmasked,rms1)
C
        write(6,*) ' '
        write(66,*) ' '
        write(6,*) 'rms before the phasefit =', rms1, '  radian'
        write(66,*) 'rms before the phasefit =', rms1, '  radian'
C
C       READING IN UNIT DIFFRACTION MAP BEFORE DOING PHASE FIT
C
        open(9, file=difffile, status='old')
        call read_diff(diffp,ndim)
C
C       PERFORMING LINEAR LEAST SQUARE FIT FOR THE CONSTANT OFFSET,
C       THE TILT TERMS, DEFOCUS AND DIFFRACTION
C
c        old call commented; replaced by new one for correct masking & fitting 
c	 11/22/97, TK.
c
c        call phasefit(ph,ndim,ndimsq,
c     &       nfit,ndish,nsub,ph1,sig,u,v,w,coef,diffp,
c     &       dprim,rate,freq,fprim,fmag,pi,clight)
c	 new call
        write(6,*) 'phasefit_aber2 call'
        call phasefit_aber2(ph,ph_um,ph1,ndim,ndimsq,
     &  nfit,ndish,nsub,sig,u,v,w,coef,diffp,
     &  dprim,rate,freq,fprim,fmag,pi,clight,nmasked)
C
C       CALCULATE RMS AFTER PHASEFIT
	call onedim(ph,ph1,ndim,nmasked)
c
        call rms1d(ph1,nmasked,rms1)
C
c        call rms(ph,ndim,ndish,nsub,rms1)
C
        write(6,*) 'rms after the phasefit =', rms1, '  radian'
        write(66,*) 'rms after the phasefit =', rms1, '  radian'
C
C       APPLYING THE MASK ON THE APERTURE PHASE FILE
	open(10,file=maskfile,status='old')
C
        call apply_mask5(ph,ndim)    !DISH AND SUB
	close(10)
C
C       WRITE OUT APERTURE PHASE AFTER THE FIT
C
        open(10, file=resiphase, status='unknown')
C
        call write_resiph(am, ph, c, ndim)
	close(10)
        open(10, file=resiph_um, status='unknown')
        call write_resiph(am, ph_um, c, ndim)
	close(10)
c
        stop
	end

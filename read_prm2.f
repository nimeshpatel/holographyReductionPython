        subroutine read_prm2(namp,inamp1,inamp2,inamp3,inamp4,
     &                defo1,defo2,defo3,defo4,
     &                outamp,outphase,resiphase,
     &                ndim,rate,distan,freq,samp_itvl,
     &                dprim,dsec,fprim,fmag,difffile,
     &                nfit,maxiter,
     &                ncent,ndish,nsub,ndimsq,taper,radphase,
     &                nmap,idum)
C
C       SUBROUTINE FOR READING PARAMETERS IN MISELL CALCULATION
C
C       HISTORY: 
C       ORIGINAL VERSION 11/27/95 X.Z.
C
        include 'holis.inc'
C
        character inamp1*15, inamp2*15            !INPUT AMPLITUDE FILENAMES 
        character inamp3*15, inamp4*15            !INPUT AMPLITUDE FILENAMES 
        character outamp*15, outphase*15          !OUTPUT AMP & PHASE FILENAMES
        character difffile*15                     !DIFFRACTION PHASE FILE 
        character resiphase*15                    !RESIDUAL PHASE FILE 
        character x*80                            !AUX VAR FOR PARAMETER INPUT
C
        write (6,*) ' '
        write (66,*) ' '
        write (6,*) '-- Reading parameters from "misell.prm" --'
        write (66,*) '-- Reading parameters from "misell.prm" --'
        write (6,*) ' '
        write (66,*) ' '
C
80      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 80
        else
           read (x,1003) a, namp
           write (6,*) 'Number of input amplitude maps =', namp
           write (66,*) 'Number of input amplitude maps =', namp
        endif
C
1       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 1
        else
           read (x,1002) a, inamp1
           write (6,*) 'First far-field amplitude file name =', inamp1
           write (66,*) 'First far-field amplitude file name =', inamp1
        endif
C
2       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 2
        else
           read (x,1002) a, inamp2
           write (6,*) 'Second far field amplitude name =', inamp2
           write (66,*) 'Second far field amplitude name =', inamp2
        endif
C
81      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 81
        else    
           read (x,1002) a, inamp3
           write (6,*) 'Third far-field amplitude file name =', inamp3
           write (66,*) 'Third far-field amplitude file name =', inamp3
        endif
C
82       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 82
        else
           read (x,1002) a, inamp4
           write (6,*) 'Fourth far field amplitude name =', inamp4
           write (66,*) 'Fourth far field amplitude name =', inamp4
        endif
C
C
91      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 91
        else
           read (x,1004) a, defo1
           write (6,*) 'The amount of the first defocus =', defo1
           write (66,*) 'The amount of the first defocus =', defo1
        endif
C
92      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 92
        else
           read (x,1004) a, defo2
           write (6,*) 'The amount of the second defocus =', defo2
           write (66,*) 'The amount of the second defocus =', defo2
        endif
C
93      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 93
        else
           read (x,1004) a, defo3
           write (6,*) 'The amount of the third defocus =', defo3
           write (66,*) 'The amount of the third defocus =', defo3
        endif
C
94      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 94
        else
           read (x,1004) a, defo4
           write (6,*) 'The amount of the fourth defocus =', defo4
           write (66,*) 'The amount of the fourth defocus =', defo4
        endif
C
3       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 3
        else
           read (x,1002) a, outamp
           write (6,*) 'Aperture field amplitude file name =', outamp
           write (66,*) 'Aperture field amplitude file name =', outamp
        endif
C
4       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 4
        else
           read (x,1002) a, outphase
           write (6,*) 'Aperture field phase file name =', outphase
           write (66,*) 'Aperture field phase file name =', outphase
        endif
C
5       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 5
        else
           read (x,1002) a, resiphase
           write (6,*) 'Residual aperture phase file name =', resiphase
           write (66,*) 'Residual aperture phase file name =', resiphase
        endif
C
C
6       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 6
        else
           read (x,1003) a, ndim
           write (6,*) 'Size N of the N by N array =', ndim
           write (66,*) 'Size N of the N by N array =', ndim
        endif
C
7       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 7
        else
           read (x,1004) a, rate
           write (6,*) 'Nyquist sampling rate =', rate
           write (66,*) 'Nyquist sampling rate =', rate
        endif
C
8       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 8
        else
           read (x,1004) a, distan
           write (6,*) 'Distance of the transmitter (meters)=',
     &                  distan
           write (66,*) 'Distance of the transmitter (meters)=',
     &                  distan
        endif
C
9       read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 9
        else
           read (x,1004) a, freq
           write (6,*) 'Observing Frequency (GHz) =', freq
           write (66,*) 'Observing Frequency (GHz) =', freq
        endif
C
10     read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 10
        else
           read (x,1004) a, samp_itvl
           write (6,*) 'Samping intervel in the far field (arcsec) =', 
     &                  samp_itvl
           write (66,*) 'Samping intervel in the far field (arcsec) =', 
     &                  samp_itvl
        endif
C
11      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 11
        else
           read (x,1004) a, dprim
           write (6,*) 'Diameter of the primary mirror (m) =', dprim
           write (66,*) 'Diameter of the primary mirror (m) =', dprim
        endif
C
12      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 12
        else
           read (x,1004) a, dsec
           write (6,*) 'Diameter of the secondary mirror (m) =', dsec
           write (66,*) 'Diameter of the secondary mirror (m) =', dsec
        endif
C
13      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 13
        else
           read (x,1004) a, fprim
           write (6,*) 'Focal length of the primaray (m) =', fprim
           write (66,*) 'Focal length of the primaray (m) =', fprim
        endif
C
14      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 14
        else
           read (x,1004) a, fmag
           write (6,*) 'Magnification of the Cassegrain system =', fmag
           write (66,*) 'Magnification of the Cassegrain system =', fmag
        endif
C
15      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 15
        else
           read (x,1002) a, difffile
           write (6,*) 'Name of phase diffraction file =', difffile
           write (66,*) 'Name of the phase diffraction file =', difffile
        endif
C
16      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 16
        else
           read (x,1003) a, nfit
           write (6,*) 'Number of the large scale fitting parameter =', 
     &                  nfit
           write (66,*) 'Number of the large scale fitting parameter =', 
     &                  nfit
        endif
C
86      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 86
        else
           read (x,1003) a, maxiter
           write (6,*) 'Maximum number of iterations =', 
     &                  maxiter
           write (66,*) 'Maximum number of iterations =', 
     &                  maxiter
        endif
C
17      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 17
        else
           read (x,1004) a, taper
           write (6,*) 'Edge taper of the initial guess aperture =', 
     &                  taper
           write (66,*) 'Edge taper of the initial guess aperture =', 
     &                  taper
        endif
C
18      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 18
        else
           read (x,1004) a, radphase
           write (6,*) 'Random phase of the initial guess aperture =', 
     &                  radphase
           write (66,*) 'Random phase of the initial guess aperture =', 
     &                  radphase
        endif
C
19      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 19
        else
           read (x,1003) a, nmap
           write (6,*) 'Number of maps with different initial seeds =',
     &                  nmap
           write (66,*) 'Number of maps with different initial seeds =',
     &                  nmap
        endif
C
20      read (52, 1001) x
        if (x(1:1).eq.'!') then
           goto 20
        else
           read (x,1003) a, idum
           write (6,*) 'First random initial seed =',
     &                  idum
           write (66,*) 'First random initial seed =',
     &                  idum
        endif
C
        ncent=ndim/2 + 1      ! center of the dish in mesh
        ndish=ndim/2*rate     ! radius of the dish in pixels
        ratio=dsec/dprim      ! ratio of the secondary to primary dia
        nsub=ndish*ratio      ! radius of the subreflector in pixels
        ndimsq=ndim*ndim      ! dimension of the 1D array
C
 1001   format (A)
 1002   format (A40, T50, A15)
 1003   format (A40, T50, I5)
 1004   format (A40, T50, f15.7)
C
	return
	end
C

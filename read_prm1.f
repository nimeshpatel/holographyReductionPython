        subroutine read_prm1(inamp,inphase,outamp,outphase,
     &                resiphase,ndim,rate,distan,freq,samp_itvl,
     &                dprim,dsec,fprim,fmag,difffile,
     &                nfit,ncent,ndish,nsub,ndimsq,near_f,defoc,
     &                out_mask,ain_mask,quad_hw,maskfile,resiph_um,
     &		      fam_um,fph_um)
C
C       SUBROUTINE FOR READING PARAMETERS IN WITHPHASE CALCULATION
C
C       ORIGINAL VERSION: 3/15/93, X.Z.
C       LAST REVISION: 11/27/95, X.Z.
c	modified to include masking parameters: 11/26/97, TK
C
        include 'holis.inc'
C
	integer dofit(nfitmax)                       ! which functions to fit
        character inamp*15, inphase*15               !INPUT AMP & PHASE FILENAME VARIABLE
        character outamp*15, outphase*15             !OUTPUT AMP & PHASE FILENAME VARIABLE
        character difffile*15                        !DIFFRACTION PHASE FILE FOR FITTING
        character resiphase*15                       !RESIDUAL PHASE FILE AFTER FITTING
        character resiph_um*15                       !umRESIDUAL PHASE FILE AFTER FITTING
        character fam_um*15                           !um aperture amp FILE 
        character fph_um*15                           !um aperture PHASE FILE 
	character maskfile*15		             !maskfile name	
        character x*80                               !AUX VARIABLE FOR PARAMETER INPUT

	common /fits/ dofit
   
   
C
        write (6,*) ' '
        write (66,*) ' '
        write (6,*) '-- Reading parameters from "withphase.prm" --'
        write (66,*) '-- Reading parameters from "withphase.prm" --'
        write (6,*) ' '
        write (66,*) ' '
C
1       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 1
        else
           read (x,1002) a, inamp
           write (6,*) 'Far field amplitude file name =', inamp
           write (66,*) 'Far field amplitude file name =', inamp
        endif
C
2       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 2
        else
           read (x,1002) a, inphase
           write (6,*) 'Far field phase file name =', inphase
           write (66,*) 'Far field phase file name =', inphase
        endif
C
3       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 3
        else
           read (x,1002) a, outamp
           write (6,*) 'Aperture field amplitude file name =', outamp
           write (66,*) 'Aperture field amplitude file name =', outamp
        endif
C
4       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 4
        else
           read (x,1002) a, outphase
           write (6,*) 'Aperture field phase file name =', outphase
           write (66,*) 'Aperture field phase file name =', outphase
        endif
C
5       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 5
        else
           read (x,1002) a, resiphase
           write (6,*) 'Residual aperture phase file name =', resiphase
           write (66,*) 'Residual aperture phase file name =', resiphase
        endif
105       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 105
        else
           read (x,1002) a, resiph_um
           write (6,*)'umResidual aperture phase file name:',resiph_um
           write (66,*)'umResidual aperture phase file name:',resiph_um
        endif
C
106       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 106
        else
           read (x,1002) a, fam_um
           write (6,*) 'um aperture amp file name =', fam_um
           write (66,*) 'um aperture amp file name =', fam_um
        endif
C
107       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 107
        else
           read (x,1002) a, fph_um
           write (6,*) 'um aperture phase file name =', fph_um
           write (66,*) 'um aperture phase file name =', fph_um
        endif
C
C
6       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 6
        else
           read (x,1003) a, ndim
           write (6,*) 'Size N of the N by N array =', ndim
           write (66,*) 'Size N of the N by N array =', ndim
        endif
C
7       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 7
        else
           read (x,1004) a, rate
           write (6,*) 'Nyquist sampling rate =', rate
           write (66,*) 'Nyquist sampling rate =', rate
        endif
C
8       read (51, 1001) x
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
9       read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 9
        else
           read (x,1004) a, freq
           write (6,*) 'Observing Frequency (GHz) =', freq
           write (66,*) 'Observing Frequency (GHz) =', freq
        endif
C
10     read (51, 1001) x
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
11      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 11
        else
           read (x,1004) a, dprim
           write (6,*) 'Diameter of the primary mirror (m) =', dprim
           write (66,*) 'Diameter of the primary mirror (m) =', dprim
        endif
C
12      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 12
        else
           read (x,1004) a, dsec
           write (6,*) 'Diameter of the secondary mirror (m) =', dsec
           write (66,*) 'Diameter of the secondary mirror (m) =', dsec
        endif
C
13      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 13
        else
           read (x,1004) a, fprim
           write (6,*) 'Focal length of the primaray (m) =', fprim
           write (66,*) 'Focal length of the primaray (m) =', fprim
        endif
C
14      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 14
        else
           read (x,1004) a, fmag
           write (6,*) 'Magnification of the Cassegrain system =', fmag
           write (66,*) 'Magnification of the Cassegrain system =', fmag
        endif
C
15      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 15
        else
           read (x,1002) a, difffile
           write (6,*) 'Name of phase diffraction file =', difffile
           write (66,*) 'Name of the phase diffraction file =', difffile
        endif
C
16      read (51, 1001) x
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
17      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 17
        else
           read (x,1003) a, near_f
           write (6,*) 'Near field correction = ', 
     &                  near_f
           write (66,*) 'Near field correction =', 
     &                  near_f
        endif
c
18      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 18
        else
           read (x,1004) a, defoc
           write (6,*) 'Defocus correction = ', 
     &                  defoc
           write (66,*) 'Defocus correction =', 
     &                  defoc
        endif
c
        ncent=ndim/2+1        ! center of the dish in mesh
        ndish=ndim/2*rate     ! radius of the dish in pixels
        ratio=dsec/dprim      ! ratio of the secondary to primary dia
        nsub=ndish*ratio      ! radius of the subreflector in pixels
        ndimsq=ndim*ndim      ! dimension of the 1D array
C
19      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 19
        else
           read (x,1004) a, out_mask 
           write (6,*) 'Outer dia for masking = ', 
     &			out_mask
           write (66,*) 'Outer dia for masking =', 
     &                  out_mask
        endif
c
20      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 20 
        else
           read (x,1004) a, ain_mask
           write (6,*) 'Inner dia for masking = ', 
     &                  ain_mask
           write (66,*) 'Inner dia for masking =', 
     &                  ain_mask
        endif
c
21      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 21 
        else
           read  (x,1004) a, quad_hw
           write (6,*) 'half-width for quadrupod masking = ', 
     &                  quad_hw
           write (66,*) 'half-width for quadrupod masking =', 
     &                  quad_hw
        endif
c
22      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 22
        else
           read (x,1002) a, maskfile
           write (6,*) 'maskfile name = ', 
     &                 maskfile 
           write (66,*) 'maskfile name=', 
     &                  maskfile
        endif
C
	do i= 1,nfit
24      read (51, 1001) x
        if (x(1:1).eq.'!') then
           goto 24
        else
           read (x,1003) a, dofit(i) 
           write (6,*) 'fitting function i=', i,' ON/OFF',dofit(i)
           write (66,*) 'fitting function i=', i,' ON/OFF',dofit(i)
C          write (6,*) a, dofit(i)
C          write (66,*) a, dofit(i) 
        endif
	enddo
c
 1001   format (A)
 1002   format (A40, T50, A15)
 1003   format (A40, T50, I5)
 1004   format (A40, T50, f15.7)
C
	return
	end
C

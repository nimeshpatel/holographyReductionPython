	program holograf3New
C-----------------------------------------------------------------------
C IMAGING ROUTINE FOR PLOTING THE GRAY SCALE AND CONTOUR MAPS OF A 2D
C ARRAY. IT IS BUILT BASED ON THE PGPLOT SUBROUTINES PGGRAY AND PGCONT
C
C MADE TO ACCOMPANY THE SMA HOLOGRAPHY PACKAGE HOLIS
C
C HISTORY:
C MAY 1993 BY X. ZHANG
C 1997 sometime 1997: plot panel boundaries, handle magic number masking: TK
C DEC 12 1997 BY MS/NP
C APR 15 1998 BY Masao Saito
C JAN 06 2000 BY Masao Saito (specify max and min separately)
C May 2004 TK to handle APEX 
C 2018,19,20 TK, GLT version - inverted elevation data.
c compile as gfortran holografGLT.f readprm_plt.f -L/usr/X11R6/lib -lX11 -lpgplot -o holografGLT
C-----------------------------------------------------------------------
C
c	include 'holis.inc'
	parameter (mxi=512, mxj=512, pi=3.1415)
        real      f(mxi, mxj)
        real      f2(mxi, mxj)
        real      tr(6),xpnt(4),ypnt(4)
        character filename1*80, filename2*80, dummy*35
        character label1*80, label2*80
	character x*80, awhite*80, ablack*80, awhite2*80, ablack2*80
        character afmin*80, afmax*80
        character postscript*1
c
        character nodefile*20
C
C OPEN DEVICE FOR GRAPHICS. IT SPECIFIES A PLOTTING AREA OF 2 BY 1 WINDOWS.
C IT ALSO PROMPTS THE USER FOR PLOTING DEVICE TYPE.
C
 	call pgbeg(0, '?', 2, 1) 
C
C INPUT 2D DATA ARRAY, WITH THE USER SPECIFIED WHITE/BLACK CORRESPONDING
C VALUES
C
	print *, 'Name for the input data array?'
c       accept *, filename1
	read(5,*), filename1
	print *, 'Name for the input data array?'
c       accept *, filename2
	read(5,*), filename2
	print *, 'Dimension N of the input N by N data array (less than or
     & equal to 512) ?'
c       accept *, n_1
	read(5,*), n_1
	print *, 'Dimension N of the input N by N data array (less than or
     & equal to 512) ?'
c       accept *, n_2
	read(5,*), n_2
C
	open(unit=41, file=filename1, form='formatted')
	open(unit=42, file=filename2, form='formatted')
C
c        do j=1,n_1
         do j=(n_1+1),2,-1
        do i=1,n_1
	   read(41,*) f(i,j)
	   if(f(i,j).eq.-9999.0) f(i,j)=10000.0
	enddo
	enddo
c        do j=1,n_2
        do j=n_2+1,2,-1
        do i=1,n_2
	   read(42,*) f2(i,j)
	   if(f2(i,j).eq.-9999.0) f2(i,j)=10000.0
	enddo
	enddo
C
C FIND THE MAXIMUM AND MINIMUM VALUES OF THE INPUT DATA ARRAY
C
        fmin=f(1,1)
        fmax=f(1,1) 
        fmin2=f(1,1)
        fmax2=f(1,1) 
	do i=1,n_1
	do j=1,n_1
	   fmin=min(f(i,j), fmin)
           fmax=max(f(i,j), fmax)
	   fmin2=min(f2(i,j), fmin2)
           fmax2=max(f2(i,j), fmax2)
	enddo
	enddo
C
C PRINT OUT THE MAXIMUM AND MINIMUM VALUES ON THE SCREEN
C
 	print *, 'The maximum value in the data file is', fmax 
 	print *, 'The minimum value in the data file is', fmin 
C
C PASS THE MAXIMUM AND MINIMUM VALUES TO CHARACTER VARIABLES
C
	print *, 'OK1'
        write(afmax,*) fmax
        write(afmin,*) fmin
	print *, 'OK2'

C ASK USER TO SPECIFY THE LIMITS THAT CORRESPONDING TO BLACK AND WHITE.
C DEFAULT BLACK=FMIN, WHITE=FMAX.
C
 	print *, 'Please specify value for WHITE (default=maximum):' 
        read (5,'(A)') x 
        if (x.eq.'  ') then
           white=fmax
        else
           read(x,*) white
        endif
C
 	print *, 'Please specify value for BLACK (default=minimum):' 
        read (5,'(A)') x 
        if (x.eq.'  ') then
           black=fmin
	else
           read(x,*) black
        endif
        write (6,*) 'white=', white
        write (6,*) 'black=', black
C
C PASS THE WHITE AND BLACK VALUES TO CHARACTER VARIABLES
C
        write(awhite,*) white
        write(ablack,*) black
C
C NOW DO THE GRAY SCALE PLOT
C
 	print *, 'NOTE THAT THE BLACK AND WHITE ON THE SCREEN' 
 	print *, 'MAY BE REVERSED FROM THE SPECIFIED VALUE (BUT'
        print *, 'WILL BE CORRECT FOR A POSTSCRIPT PRINTOUT) !!!' 
C
C CLEAR THE SCREEN. SET UP WINDOW AND VIEWPORT
C
	call pgpage
	call pgsvp(0.05,0.95,0.05,0.95)
        x1=0.5   
        y1=0.5
        x2=real(n_1)+0.5
        y2=real(n_1)+0.5
 	call pgwnad(x1,x2,y1,y2)
       	call pgsci(5)
        call pgsch(2.0)
 	call pgmtxt('t',3.0,0.55,0.5, filename1)
C 	call pgmtxt('t',3.0,1.05,0.5, filename1)
        call pgsch(1.0)
 	call pgmtxt('t',2.0,0.2,0.0,'White/Black: ')
        dummy=awhite//'/'//ablack
 	call pgmtxt('t',2.0,0.6,0.5,dummy)
	print *, 'Wedge label1?'
c       accept *, label1
	read(5,*), label1
 	call pgmtxt('B',8.0,0.48,0.3,label1)
	call pgsci(1)
C
C DRAW THE MAP
C
        tr(1)=0
        tr(2)=1.0
        tr(3)=0.
        tr(4)=0.
        tr(5)=0.
        tr(6)=1.0
	call pggray(f,mxi,mxj,1,n_1,1,n_1,black,white,tr)
        call pgsci(5)
        call pgbox('bcnts',0.0,0,'bcnts',0.0,0)
C
C DRAW A WEDGE
C
	call pgwedg('B',3.0,4.0, black,white,' ')
c
c	draw the panel boundaries
c
	print *, 'Want to draw panel boundaries(1/0)?'
c	accept *, pan_ans
	read(5,*), pan_ans
	if (pan_ans.eq.1) then
	print *,'panels!'
       open (51, status='old', file='panelplt.prm')
	offset=0.0
       call readprm_plt(nodefile,ndim,
     &             Ra1,Ra2,Ra3,Ra4,Ra5,Ra6,Ra7,Ra8,Ra9,
     &             n1,n2,n3,n4,n5,n6,n7,n8,
     &             dists,fprim,dprim,
     &             dsec,rate,dr,dtheta,offset)
        print *,ndim,Ra1,Ra2,Ra3,Ra4,Ra5,Ra6,Ra7,Ra8,Ra9,
     &          n1,n2,n3,n4,n5,n6,n7,n8,dists,fprim,dprim,
     &          dsec,rate,dr,dtheta,offset

        open(3,file=nodefile,status='old',form='formatted')
c
        ncent=ndim/2 + 1
        ndish=ndim/2*rate
        nsub=ndish*dsec/dprim
        dx=dprim/2./ndish
        npanels=n1+n2+n3+n4+n5+n6+n7+n8
c
c	draw ring 1
	do i = 1,n1
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(2) = vr2*cos(phi)/dx+ncent
	ypnt(2) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(3) = vr3*cos(phi)/dx+ncent
	ypnt(3) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(4) = vr4*cos(phi)/dx+ncent
	ypnt(4) = vr4*sin(phi)/dx+ncent
	phi = (pi-vp1-offset)
	xpnt(1) = vr1*cos(phi)/dx+ncent
	ypnt(1) = vr1*sin(phi)/dx+ncent
	
	call pgline(4,xpnt,ypnt)
	enddo

c	draw ring 2
	
	do i = 1,n2
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp3-offset)
	xpnt(1) = vr3*cos(phi)/dx+ncent
	ypnt(1) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(2) = vr4*cos(phi)/dx+ncent
	ypnt(2) = vr4*sin(phi)/dx+ncent
	call pgline(2,xpnt,ypnt)
	enddo
	
c	draw ring 3	

	do i = 1,n3
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	phi = (pi-vp1-offset)
	xpnt(4) = vr1*cos(phi)/dx+ncent
	ypnt(4) = vr1*sin(phi)/dx+ncent
	call pgline(4,xpnt,ypnt)
	enddo

c	draw ring 4
	do i = 1,n4
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 5
	do i = 1,n5
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 6
	do i = 1,n6
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 7
	do i = 1,n7
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 8
	do i = 1,n8
	print *,i
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
	endif
c	done drawing panels!
	close(3)
C
C PRINT OUT THE MAXIMUM AND MINIMUM VALUES ON THE SCREEN
C
        print *, 'The maximum value in the data file is', fmax2
        print *, 'The minimum value in the data file is', fmin2
C
C PASS THE MAXIMUM AND MINIMUM VALUES TO CHARACTER VARIABLES
C
        write(afmax,*) fmax2
        write(afmin,*) fmin2

C ASK USER TO SPECIFY THE LIMITS THAT CORRESPONDING TO BLACK AND WHITE.
C DEFAULT BLACK=FMIN, WHITE=FMAX.
C
        print *, 'Please specify value for WHITE (default=maximum):'
        read (5,'(A)') x
        if (x.eq.'  ') then
           white2=fmax2
        else
           read(x,*) white2
        endif
C
        print *, 'Please specify value for BLACK (default=minimum):'
        read (5,'(A)') x
        if (x.eq.'  ') then
           black2=fmin2
        else
           read(x,*) black2
        endif
        write (6,*) 'white=', white2
        write (6,*) 'black=', black2
C
C PASS THE WHITE AND BLACK VALUES TO CHARACTER VARIABLES
C
        write(awhite2,*) white2
        write(ablack2,*) black2
C
C CLEAR THE SCREEN. SET UP WINDOW AND VIEWPORT
C
        call pgpage
        call pgsvp(0.05,0.95,0.05,0.95)
        x1=0.5
        y1=0.5
        x2=real(n_2)+0.5
        y2=real(n_2)+0.5
        call pgwnad(x1,x2,y1,y2)
        call pgsci(5)
        call pgsch(2.0)
C       call pgmtxt('t',3.0,1.05,0.5, filename2)
        call pgmtxt('t',3.0,0.55,0.5, filename2)
        call pgsch(1.0)
        call pgmtxt('t',2.0,0.2,0.0,'White/Black: ')
        dummy=awhite2//'/'//ablack2
        call pgmtxt('t',2.0,0.6,0.5,dummy)
	print *, 'Wedge label2?'
c       accept *, label2
	read(5,*), label2
 	call pgmtxt('B',8.0,0.48,0.3,label2)
        call pgsci(1)
C
C DRAW THE MAP
C
        tr(1)=0
        tr(2)=1.0
        tr(3)=0.
        tr(4)=0.
        tr(5)=0.
        tr(6)=1.0
        call pggray(f2,mxi,mxj,1,n_2,1,n_2,black2,white2,tr)
        call pgsci(5)
        call pgbox('bcnts',0.0,0,'bcnts',0.0,0)
C
C
C DRAW A WEDGE
C
	call pgwedg('B',3.0,4.0, black2,white2,' ')
C-----------------------------------------------------------
C
c       draw the panel boundaries
c
c        print *, 'Want to draw panel boundaries(1/0)?'
c       accept *, pan_ans
c        read(5,*), pan_ans
        if (pan_ans.eq.1) then
        print *,'panels!'
c       open (51, status='old', file='panelplt.prm')
c       call readprm_plt(nodefile,ndim,Ra1,Ra2,Ra3,
c     &             Ra4,Ra5,n1,n2,n3,n4,dists,fprim,dprim,
c     &             dsec,rate,dr,dtheta,offset)
c        print *,ndim,Ra1,Ra2,Ra3,
c     &             Ra4,Ra5,n1,n2,n3,n4,dists,fprim,dprim,
c     &             dsec,rate,dr,dtheta,offset

        open(3,file=nodefile,status='old',form='formatted')
c
        ncent=ndim/2 + 1
        ndish=ndim/2*rate
        nsub=ndish*dsec/dprim
        dx=dprim/2./ndish
        npanels=n1+n2+n3+n4
c
c       draw ring 1
        do i = 1,n1
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(2) = vr2*cos(phi)/dx+ncent
        ypnt(2) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(3) = vr3*cos(phi)/dx+ncent
        ypnt(3) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(4) = vr4*cos(phi)/dx+ncent
        ypnt(4) = vr4*sin(phi)/dx+ncent
        phi = (pi-vp1-offset)
        xpnt(1) = vr1*cos(phi)/dx+ncent
        ypnt(1) = vr1*sin(phi)/dx+ncent

        call pgline(4,xpnt,ypnt)
        enddo

c       draw ring 2

        do i = 1,n2
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp3-offset)
        xpnt(1) = vr3*cos(phi)/dx+ncent
        ypnt(1) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(2) = vr4*cos(phi)/dx+ncent
        ypnt(2) = vr4*sin(phi)/dx+ncent
        call pgline(2,xpnt,ypnt)
        enddo

c       draw ring 3

        do i = 1,n3
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(1) = vr2*cos(phi)/dx+ncent
        ypnt(1) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(2) = vr3*cos(phi)/dx+ncent
        ypnt(2) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(3) = vr4*cos(phi)/dx+ncent
        ypnt(3) = vr4*sin(phi)/dx+ncent
        phi = (pi-vp1-offset)
        xpnt(4) = vr1*cos(phi)/dx+ncent
        ypnt(4) = vr1*sin(phi)/dx+ncent
        call pgline(4,xpnt,ypnt)
        enddo

c       draw ring 4
        do i = 1,n4
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(1) = vr2*cos(phi)/dx+ncent
        ypnt(1) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(2) = vr3*cos(phi)/dx+ncent
        ypnt(2) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(3) = vr4*cos(phi)/dx+ncent
        ypnt(3) = vr4*sin(phi)/dx+ncent
        call pgline(3,xpnt,ypnt)
        enddo
c	draw ring 5
	do i = 1,n5
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 6
	do i = 1,n6
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 7
	do i = 1,n7
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 8
	do i = 1,n8
	print *,i
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
	endif
c       done drawing panels!
        close(3)
        print *, 'Produce postscript output (y/n, default=n)?'
        read (5, '(A)') postscript
C
        if (postscript.eq.'y') then
C
           call pgbeg(0, '/ps', 2, 1)
           call pgscf(2)
C
C CLEAR THE SCREEN. SET UP WINDOW AND VIEWPORT
C
           call pgpage
           call pgsvp(0.05,0.95,0.05,0.95)
           x1=0.5
           y1=0.5
           x2=real(n_1)+0.5
           y2=real(n_1)+0.5
           call pgwnad(x1,x2,y1,y2)
           call pgsci(5)
           call pgsch(2.0)
           call pgmtxt('t',3.0,0.55,0.5, filename1)
           call pgsch(1.0)
           call pgmtxt('t',2.0,0.2,0.0,'White/Black: ')
           dummy=awhite//'/'//ablack
           call pgmtxt('t',2.0,0.6,0.5,dummy)
 	   call pgmtxt('B',8.0,0.48,0.3,label1)
           call pgsci(1)
C
C DRAW THE MAP
C
           tr(1)=0
           tr(2)=1.0
           tr(3)=0.
           tr(4)=0.
           tr(5)=0.
           tr(6)=1.0
           call pggray(f,mxi,mxj,1,n_1,1,n_1,black,white,tr)
           call pgsci(5)
           call pgbox('bcnts',0.0,0,'bcnts',0.0,0)
C
C DRAW A WEDGE
C
           call pgwedg('B',3.0,4.0, black,white,' ')
c	draw panels
       if (pan_ans.eq.1) then
        
	open(3,file=nodefile,status='old',form='formatted')
c
        ncent=ndim/2 + 1
        ndish=ndim/2*rate
        nsub=ndish*dsec/dprim
        dx=dprim/2./ndish
        npanels=n1+n2+n3+n4
c
c       draw ring 1
        do i = 1,n1
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(2) = vr2*cos(phi)/dx+ncent
        ypnt(2) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(3) = vr3*cos(phi)/dx+ncent
        ypnt(3) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(4) = vr4*cos(phi)/dx+ncent
        ypnt(4) = vr4*sin(phi)/dx+ncent
        phi = (pi-vp1-offset)
        xpnt(1) = vr1*cos(phi)/dx+ncent
        ypnt(1) = vr1*sin(phi)/dx+ncent

        call pgline(4,xpnt,ypnt)
        enddo

c       draw ring 2

        do i = 1,n2
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp3-offset)
        xpnt(1) = vr3*cos(phi)/dx+ncent
        ypnt(1) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(2) = vr4*cos(phi)/dx+ncent
        ypnt(2) = vr4*sin(phi)/dx+ncent
        call pgline(2,xpnt,ypnt)
        enddo

c       draw ring 3

        do i = 1,n3
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(1) = vr2*cos(phi)/dx+ncent
        ypnt(1) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(2) = vr3*cos(phi)/dx+ncent
        ypnt(2) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(3) = vr4*cos(phi)/dx+ncent
        ypnt(3) = vr4*sin(phi)/dx+ncent
        phi = (pi-vp1-offset)
        xpnt(4) = vr1*cos(phi)/dx+ncent
        ypnt(4) = vr1*sin(phi)/dx+ncent
        call pgline(4,xpnt,ypnt)
        enddo

c       draw ring 4
        do i = 1,n4
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(1) = vr2*cos(phi)/dx+ncent
        ypnt(1) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(2) = vr3*cos(phi)/dx+ncent
        ypnt(2) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(3) = vr4*cos(phi)/dx+ncent
        ypnt(3) = vr4*sin(phi)/dx+ncent
        call pgline(3,xpnt,ypnt)
        enddo
        endif
c       done drawing panels!
	close(3)
C
C CLEAR THE SCREEN. SET UP WINDOW AND VIEWPORT
C
        call pgpage
        call pgsvp(0.05,0.95,0.05,0.95)
        x1=0.5
        y1=0.5
        x2=real(n_2)+0.5
        y2=real(n_2)+0.5
        call pgwnad(x1,x2,y1,y2)
        call pgsci(5)
        call pgsch(2.0)
        call pgmtxt('t',3.0,0.55,0.5, filename2)
        call pgsch(1.0)
        call pgmtxt('t',2.0,0.2,0.0,'White/Black: ')
        dummy=awhite2//'/'//ablack2
        call pgmtxt('t',2.0,0.6,0.5,dummy)
        call pgsci(1)
C
C DRAW THE MAP
C
        tr(1)=0
        tr(2)=1.0
        tr(3)=0.
        tr(4)=0.
        tr(5)=0.
        tr(6)=1.0
        call pggray(f2,mxi,mxj,1,n_2,1,n_2,black2,white2,tr)
        call pgsci(5)
        call pgbox('bcnts',0.0,0,'bcnts',0.0,0)
C
C DRAW A WEDGE
C
	call pgwedg('B',3.0,4.0, black,white,' ')
        endif
C
C
C CLOSE THE DEVICE AND EXIT.
C
c       draw the panel boundaries
c
c        print *, 'Want to draw panel boundaries(1/0)?'
c       accept *, pan_ans
c        read(5,*), pan_ans
        if (pan_ans.eq.1) then
        print *,'panels!'
c       open (51, status='old', file='panelplt.prm')
c       call readprm_plt(nodefile,ndim,Ra1,Ra2,Ra3,
c     &             Ra4,Ra5,n1,n2,n3,n4,dists,fprim,dprim,
c     &             dsec,rate,dr,dtheta,offset)
        print *,ndim,Ra1,Ra2,Ra3,
     &             Ra4,Ra5,n1,n2,n3,n4,dists,fprim,dprim,
     &             dsec,rate,dr,dtheta,offset

        open(3,file=nodefile,status='old',form='formatted')
c
        ncent=ndim/2 + 1
        ndish=ndim/2*rate
        nsub=ndish*dsec/dprim
        dx=dprim/2./ndish
        npanels=n1+n2+n3+n4+n5+n6+n7+n8
c
c       draw ring 1
        do i = 1,n1
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(2) = vr2*cos(phi)/dx+ncent
        ypnt(2) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(3) = vr3*cos(phi)/dx+ncent
        ypnt(3) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(4) = vr4*cos(phi)/dx+ncent
        ypnt(4) = vr4*sin(phi)/dx+ncent
        phi = (pi-vp1-offset)
        xpnt(1) = vr1*cos(phi)/dx+ncent
        ypnt(1) = vr1*sin(phi)/dx+ncent

        call pgline(4,xpnt,ypnt)
        enddo

c       draw ring 2

        do i = 1,n2
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp3-offset)
        xpnt(1) = vr3*cos(phi)/dx+ncent
        ypnt(1) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(2) = vr4*cos(phi)/dx+ncent
        ypnt(2) = vr4*sin(phi)/dx+ncent
        call pgline(2,xpnt,ypnt)
        enddo

c       draw ring 3

        do i = 1,n3
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(1) = vr2*cos(phi)/dx+ncent
        ypnt(1) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(2) = vr3*cos(phi)/dx+ncent
        ypnt(2) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(3) = vr4*cos(phi)/dx+ncent
        ypnt(3) = vr4*sin(phi)/dx+ncent
        phi = (pi-vp1-offset)
        xpnt(4) = vr1*cos(phi)/dx+ncent
        ypnt(4) = vr1*sin(phi)/dx+ncent
        call pgline(4,xpnt,ypnt)
        enddo

c       draw ring 4
        do i = 1,n4
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

        phi = (pi-vp2-offset)
        xpnt(1) = vr2*cos(phi)/dx+ncent
        ypnt(1) = vr2*sin(phi)/dx+ncent
        phi = (pi-vp3-offset)
        xpnt(2) = vr3*cos(phi)/dx+ncent
        ypnt(2) = vr3*sin(phi)/dx+ncent
        phi = (pi-vp4-offset)
        xpnt(3) = vr4*cos(phi)/dx+ncent
        ypnt(3) = vr4*sin(phi)/dx+ncent
        call pgline(3,xpnt,ypnt)
        enddo
c	draw ring 5
	do i = 1,n5
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 6
	do i = 1,n6
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 7
	do i = 1,n7
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
c	draw ring 8
	do i = 1,n8
	print *,i
C
C          READING IN VERTEX COORDINATES FROM TABLE
C
           read (3,'(A)') panelname
           read (3,*) vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
           print *,vr1,vp1,vr2,vp2,vr3,vp3,vr4,vp4
C          DETERMINE THE MAX AND MIN OF R, THETA ON THE PANEL
C          ALSO DETERMINE THE MAX MIN OF R< THETA FOR FITTING
C          GIVING SAFeTY MARGIN OF DR, DTHETA
C
           rmaxfit=vr2-dr
           rminfit=vr1+dr
           thetamaxfit=vp4-dtheta
           thetaminfit=vp2+dtheta

	phi = (pi-vp2-offset)
	xpnt(1) = vr2*cos(phi)/dx+ncent
	ypnt(1) = vr2*sin(phi)/dx+ncent
	phi = (pi-vp3-offset)
	xpnt(2) = vr3*cos(phi)/dx+ncent
	ypnt(2) = vr3*sin(phi)/dx+ncent
	phi = (pi-vp4-offset)
	xpnt(3) = vr4*cos(phi)/dx+ncent
	ypnt(3) = vr4*sin(phi)/dx+ncent
	call pgline(3,xpnt,ypnt)
	enddo
	endif
c       done drawing panels!
        close(3)
        call pgend
	print *, dprim,ndim,ndish,rate,dx
        end

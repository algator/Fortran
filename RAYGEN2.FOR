        program raygen2
c       last change/update/compile: Mon. July 25,1994 [6:52 p.m.]
c
c		RAYGEN2.FOR 
c		version 1
c
c 
c	Author: Alexandre Y. Fong, Cirrus Technology Consulting
c				   Copyright 1994
c
c 	A program to generate a simulated Rayleigh extinction
c	profile with clouds.
c..............................................................
c
c
c	arrays and constants
c
c
c	rez: lidar resolution = 40 m
c	rng: range
c	npts: number of total range points = 100
c	pr:pressure
c	tmp:temperature
c	rh:relative humidity
c	ntot: number of points in sonde data
c..............................................................
c
                real rng(100)
		DIMENSION RAYN(100),SIGR(100),BR(100)
		real sig(100),k(100),hcl(10),excl(10),wid(10) 
		real rez,R,g,larry,tkel
                real pr(100),tmp(100),dew(20),windvel(20)
                integer windir(20),rh(20),ht(20),ih
                integer isw,npts,ihcl(10),iwid(10),inc
                character*60 outfile,sonde
c
c	open parameter, input and output files
c
		open (unit=8,file='master.dat',status='old')
c		write(*,'(a\)') '\name of file to write to?'
		read(8,'(a)') outfile
		open (unit=20,file=outfile,status='new')
c		write(*,*) '\do you want to use standard lapse rate (0) or'
c		write(*,*) '\load a sonde (1)?'
		read(8,*) isw
c		if (isw.eq.1) then
c		write(*,*) '\number of points in sonde data?'
		read(8,*) ntot
c		write(*,'(a\)') '\name of sonde file to use for Rayleigh?'
		read(8,'(a)') sonde
               if (isw.eq.1) then
		open (unit=9,file=sonde,status='old')
		else
		end if
c                else
c               write(*,*) '\number of points?'
		read(8,*) npts
c		write(*,*) '\resolution (metres)?'
                read(8,*) rez
c		end if
c..............................................................
c
c	Prepare arrays:
c
	p0=1013.2
	g=9.8
	R=8.31432
	M=28.9644
	do 1 i=1,npts
	rng(i)=i*rez
	RAYN(i)=0.0
	SIGR(i)=0.0
	BR(i)=0.0
	sig(i)=0.0
	k(i)=0.0
1	continue
c..............................................................
c
c
c	NOTE: 
c	load sonde data (look for tabulated data, always starts at pr=1000
c			 and is preceeded by 3 empty lines)
c
c
	if (isw.eq.1) then
	npts=ntot
	rewind 9
56      read(9,*,end=58) larry
	if(larry.eq.1000.) then
	backspace (9)
	do 55 kl=1,ntot
c	         ^^^^
c	how many points in sonde data?
c
        read(9,57,end=58) pr(kl),ht(kl),tmp(kl),
     +rh(kl),dew(kl),windir(kl),windvel(kl)
57      format(1x,f6.1,1x,I5,1x,f5.1,2x,I2,2x,f5.1,2x,I3,3x,f4.1)
c	if sonde are in feet...
c	ht(kl)=0.3*ht(kl)
55	continue
        else
        goto 56
        end if
	do 88 i=1,ntot
        write(*,57) pr(i),ht(i),tmp(i),
     +rh(i),dew(i),windir(i),windvel(i)
88	continue
        else 
c..............................................................
c
c	lapse rate version
c	
	do 2 i=1,npts
        tmp(i)=15.0-6.5*(rng(i)/1000.)
	tkel=273-tmp(i)
        pr(i)=p0*exp((-1.0*g*M*(rng(i)/1000.))/(R*tkel))
c        write(*,*) rng(i),pr(i),tmp(i)
2	continue
         END IF
c
c..............................................................
c
c	calculate rayliegh value
c
c	** Call rayliegh profile subroutine - modified for 1064 nm **
c
58		call rayshiv(npts,tmp,pr,sig,br,rayn)
c
c	add aerosol correction?
c
c	load model and add to above profile
c..............................................................
c
c	time to insert clouds
c
c	write(*,*) '\number of cloud layers'
	read(8,*) inc
	do 77 i=1,inc
c	write(*,*) '\in increasing base height, 
c     + cloud layer height with maximum extinction value (1/km) and
c     +	corresponding geometric thickness'
        read(8,*) hcl(i),excl(i),wid(i)
        excl(i)=excl(i)/1000.
	ihcl(i)= ifix(hcl(i)/rez)
	iwid(i)= ifix(wid(i)/(2.*rez))
	do 75 l=(ihcl(i)-iwid(i)),(ihcl(i)+iwid(i))
        ih=iwid(i)-abs(l-ihcl(i))
	if (ih.eq.0.) then
        goto 75
        else
        end if
	sig(l)=(excl(i)/iwid(i))*ih
c	write(*,*) ihcl(i),iwid(i),ih,sig(l)
75	continue
77	continue

c..............................................................
c
c	add noise
c
c..............................................................
c
c	write optical extinction data to file 19 in Andre's format  
c
	do 3 i=1,npts
	if (sig(i).eq.0.) then
	k(i)=0.
	else
	k(i)=br(i)/sig(i)
	end if
	sig(i)=sig(i)*1000.
	br(i)=br(i)*1000.
	if (isw.eq.1) then
        write(20,7773) ht(i),sig(i),br(i),k(i)
	else
        write(20,7772) rng(i),sig(i),br(i),k(i)
	end if
7772    format(1x,f7.0,1x,e10.3,1x,e10.3,1x,e10.3)
7773    format(1x,I5,1x,e10.3,1x,e10.3,1x,e10.3)
3	continue
c..............................................................
c
c	close all files...
c
		close(8)
		close(9)
		close(20)
c
                stop
                end
c
c..................................................................
c
	subroutine rayshiv(ntot,temp,prez,sigr,br,rayn)
c		rayshiv.for	
c	A variation S.R. Pal's program for generating a number density 
c 	reference extinction and backscatter coefficient 
c	from temperature and pressure known as raysig.for
c
c	MODIFIED FOR 1064nm ************************************
c	
c-------RAYSIG.FOR FOR CALCULATING RAYLEIGH NUMBER DENSITY FROM
c       PRESSURE AND TEMPERATURE  AND RAYLEIGH ATTENUATION 
c	X1(I)=ALOG10(X(I))
c	Y1(I)=ALOG10(Y(I))
c 
c  RAYLEIGH SCAT CROSSSECTION FOR  532 nm IS OMEGA=5.2170E-31 m.sqred
c                             FOR 1064 nm IS OMEGA=3.2597E-32 m2

	DIMENSION RAYN(100),SIGR(100),BR(100)
	REAL PREZ(100),TEMP(100)
	REAL WTM,RGAS,AVOG,OMEGAR
	INTEGER NTOT

	WTM=28.9644
	RGAS=8.31432
	AVOG=6.022E+23
 	F=AVOG*100./RGAS
	OMEGAR=3.2597E-32
	PHASER=3.0/(8.0*3.14159)
	DO 70 I=1,NTOT
		TEMP(I)=273.+TEMP(I)
		RAYN(I)=F*PREZ(I)/TEMP(I)
  		SIGR(I)=OMEGAR*RAYN(I) 
		BR(I)=PHASER*SIGR(I)
c	write(*,*) RAYN(i),SIGR(i),BR(i)
		TEMP(I)=TEMP(I)-273
70	CONTINUE
c	write(*,'(a\)') 'Rayleigh profile done'
	return 
	end

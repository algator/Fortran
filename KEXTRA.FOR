        program kextra
c       last change/update/compile: Thurs. Sept 8,1994 [12:08 noon]
c
c		Kextra.FOR 
c		version 1
c
c 
c	Author: Alexandre Y. Fong, Cirrus Technology Consulting
c				   Copyright 1994
c
c 	A program to  extract individual traces from
c	*riv.dat or extinction files and write them to column form
c..............................................................
c
c
c	arrays and constants
c
c	cloud base portion:
c
c	sgnl: signal
c	rez: lidar resolution = 40 m
c	sp:log-range adjusted signal 
c	rng: range
c	z: range index
c	npts: number of total range points = 100
c	
		real sgnl(600),sgnlnrc(600),rng(600)
		integer npts,itest,icount,ieb,sw,sw2     		
                character*60 lidar,outfile
c
		write(*,'(a\)') '\name of lidar file to look at?'
                read(*,'(a)') lidar
		open (unit=10,file=lidar,status='old',form='formatted')
		write(*,'(a\)') '\name of file to write to?'
		read(*,'(a)') outfile
		open (unit=11,file=outfile,status='new')
c
c	read program parameters
c
		write(*,*) 'npts?'
		read(*,*) npts
		write(*,*) 'rez,zero?'
		read(*,*) rez,zero
		write(*,*) 'total number of traces?'
		read(*,*) ieb
		write(*,*) 'trace number to extract?'
		read(*,*) itest
		write(*,*) 'RiVM data(1) or extinction (0)?'
		read(*,*) sw
c		write(*,*) 'range correct?(1)'
c		read(*,*) sw2
c
		icount=0
		do 2000 ll=1,ieb
		icount=icount+1
		write(*,*)'reading trace:',icount
c
c	set signal registers to zero fill range registers
c
		do 78 i=1,npts
			sgnl(i)=0.00 
			sgnlnrc(i)=0.00 
			rng(i)=i*rez
78		continue
c
c	
c	load lidar data
c	
	if (sw.eq.1) then	
		read(10,209,end=1999) (sgnl(i),i=1,npts)
209		format(100f7.0)
	else if (sw.eq.0) then
		read(10,210,end=1999) (sgnl(i),i=1,npts)
210		format(100e10.3)
	else
	end if
c		
c		check for zeros?
c
1999                if (sw.eq.1) then
                do 79 i=1,npts
                        if (sgnl(i).lt.zero) sgnl(i)=zero
79	continue
        else
	end if
c
c	if RIVM data...remove range correction!
c	RIVM data is range corrected ie.= P(r)*Index of r *some scaling factor
c
                do 769 i=1,npts
                if (sw.eq.1) then
        sgnlnrc(i)=sgnl(i)/((i*1.)**2)
	else
        sgnlnrc(i)=sgnl(i)*((i*1.)**2)
	end if
769     continue
c	
        if (icount.eq.itest) then
	goto 99
	else 
	end if
c
c	next file
c
2000		continue
c...................................................................
c
c	write desired trace data to file
c	
99	do 4445 j=1,npts
	if (sw.eq.1) then 
	write(11,44422) rng(j),sgnl(j),sgnlnrc(j)
	else
	write(11,4442) rng(j),sgnl(j),sgnlnrc(j)
	end if
4445	continue
4442	format(1x,f7.2,1x,e10.4,1x,e10.4)
44422	format(1x,f7.2,1x,e10.4,1x,e10.4)
	write(*,*) 'How steak is done!'
c
c
c..............................................................
c
c	close all files
c
		close(10)
		close(11)
c
                stop
                end

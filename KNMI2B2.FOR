        program knmi2b2
c
c..............................................................
c
c		KNMI2B2.FOR 
c		version 2
c
c 
c	Author: Alexandre Y. Fong, Cirrus Technology Consulting
c				   Copyright 1994
c..............................................................
c
c	[In a Rod Serling voice:] DooDoo-DooDoo
c
c 	Submitted for your approval: a program to calculate cloud 
c	optical extinction and depth 
c	values from and correlate between rivm lidar 
c	and knmi rawidnsonde data (for 1 layer only so far).
c	This version (2B2) allows truthing of data by 
c	comparing the signal just above and below clouds
c	and taking the log of their scattering ratio to obtain
c	a mean value and depth (unit=20) AND reads the *.t files for
c  	saturation flags (unit=22). Saturated traces are NOT inverted. 
c	In addition a non-range corrected backscatter is written
c	to a file (unit=21). Much signal averaging and processing
c	hocus-pocus takes place to make the inversion work...
c
c	In The Twighlight Zone!
c
c...........................................................
c
c	arrays and constants
c
c	cloud base variables:
c
c	sgnl: signal
c	rez: lidar resolution = 40 m
c	sp:log-range adjusted signal 
c	rng: range
c	z: range index
c	npts: number of total range points = 100-200
c	mpt: how many point fit for cloud base algorithm peak search
c	dpt:how many point fit for cloud base algorithm derrivative
c	xovlp: overlap height, clouds below will not register
c	nob: number of traces to average over for cloud base derivation
c	     ussually = 1	
c	nos: number of shots averaged for this trace = 250
c	gain: i.r. gain setting
c	nd: overall neutral density
c	indx: number of cloud layers not to be confused with indx
c		in the subroutine which is something else!
c	ieb: total number of traces to read
c
c..............................................................
c
c	optical variables:
c
c	pr:pressure
c	tmp:temperature
c	rh:relative humidity
c	ntot: number of points in sonde data
c	sigref: reference extinction coefficient
c	iref: inversion extinction reference range index
c	iref2: rayliegh-sonde extinction reference range index
c	k: backscatter/extinction ratio exponent
c	thrush: threshold for diffrence allowed in correlation 
c	algorithm
c	iabove: how many range bins above cloud top to 
c	use as backround for inversion
c..............................................................
c
                real sgnl(200),fit(200)
		real rng(200),mat(9),lambda
                real logsgnl(200),rezkm
                real sp(200),mdcldalt,larry
		real nd,gain,cldthk(9),mdht(9)
		DIMENSION RAYN(50),SIGR(50),BR(50)
		real z0,zt,sig(200),dep,meanx,divlen
		real rez,krez,k,xovlp,dummy,zero,high,conh
                real toll,thrush,diff,it,plusig,avy
                real dpp(3),mnx(3),multi,offset,radd
		real pr(20),tmp(20),dew(20),windvel(20)
		real topdep(9),tmeanx(9),slop,bay,htest
		integer windir(20),rh(20),ht(20),nero,pyt,ls,blt
		integer mdcldx,nob,ieb,ll,iabove,mink,inpts,avlen
		integer itopx,ibotx,isat,baseline,add,seth,asw,rrcor
		integer skipper,botlen,toplen,me,fav,nlog,mixbit
                integer*2 jbase(1,6),jtop(1,6)
		integer*2 jpeak(1,6),jtrace,n,iwk
		integer mode,indx,kl,li,npts,mpt,dpt,tpt,lr
		integer*4 i,l,z,iref,iref2,ntot    
                real info(4)
                character*60 lidar,outfile,outfile1,outfile2,sonde
                character*60 outfile3,outfile4,outfile5,outfile6
		character*60 tea,outfile7
		character sat
c...........................................................
c
c	Advertising... The right to choose!
c
	write(*,*)
	write(*,*)'	     ..................................'
	write(*,*)'	     :	         CIRRUSWARE	      :'
	write(*,*)'	     :		    from	      :'
	write(*,*)'	     :	CIRRUS TECHNOLOGY CONSULTING  :'
        write(*,*)'             :     Copyright 1994, Canada     :'
	write(*,*)'	     :	      (416) 650-5343	      :'	
        write(*,*)'             ..................................'
	write(*,*)
	write(*,*) '		 KNMI2B2'
	write(*,*)
        write(*,*) 'last update/compile:Wed.Oct26th,1994[12:50p.m.]'
	write(*,*)
	write(*,*) '<ENTER> to run...'
	read(*,*)
c...........................................................
c
c	open program parameter, file list input and output files
c
		open (unit=8,file='para.dat',status='old')
		open (unit=14,file='knmi2b.lis',status='old')
c		write(*,'(a\)') '\name of lidar file to look at?'
                read(14,'(a)',end=1999) lidar
c		write(*,'(a\)') '\name of sonde file to use for Rayleigh?'
		read(14,'(a)') sonde
		open (unit=10,file=lidar,status='old',form='formatted')
		open (unit=9,file=sonde,status='old')
c		write(*,'(a\)') '\names of files to write 
c      +	cloud bases, ext.write, correlated, Andre format, 
c      +	truthing and non-range corrected data to?'
		read(14,'(a)') outfile
		read(14,'(a)') outfile1
		read(14,'(a)') outfile2
		read(14,'(a)') outfile3
		read(14,'(a)') outfile4
		read(14,'(a)') outfile5
		read(14,'(a)') outfile6
		open (unit=11,file=outfile,status='new')
		open (unit=16,file=outfile1,status='new')
		open (unit=17,file=outfile2,status='new')
		open (unit=18,file=outfile3,status='new')
		open (unit=19,file=outfile4,status='new')
		open (unit=20,file=outfile5,status='new')
		open (unit=21,file=outfile6,status='new')
c		write(*,'(a\)') '\name of *.t file?'
		read(14,'(a)') tea
		open (unit=22,file=tea,status='old')
c		write(*,'(a\)') '\name of processed signal data file?'
		read(14,'(a)') outfile7
		open (unit=23,file=outfile7,status='new')
c...........................................................
c
c	write(17,'(a\)') 'no. layers:  opt. depth 1,2,3:  mean ext.1,2,3:'
c	write(17,*) 
c...........................................................
c
c	read program parameters
c
		read(8,*)
		read(8,*) npts,mpt,dpt,tpt,inpts,rrcor,lambda
		write(*,*) npts,mpt,dpt,tpt,inpts,rrcor,lambda
		read(8,*)
		read(8,*) rez,zero,nob,high,nlog,mixbit
		write(*,*) rez,zero,nob,high,nlog,mixbit
		read(8,*)
		read(8,*) corrfak,xovlp,ieb,k,seth,asw,conh,avlen
		write(*,*) corrfak,xovlp,ieb,k,seth,asw,conh,avlen
                read(8,*)
		read(8,*) nd,gain,nos,multi,offset,fav
		write(*,*) nd,gain,nos,multi,offset,fav
		read(8,*)
		read(8,*) thrush,ntot,iabove,iwk,isat,baseline,add,plusig
		write(*,*) thrush,ntot,iabove,iwk,isat,baseline,add,plusig
c
c...........................................................
		do 200 ll=1,ieb
c
c	set signal registers to zero fill range registers
c
		do 78 i=1,npts
			sgnl(i)=0.00
			sp(i)=0.00
			rng(i)=i*rez
			sig(i)=0.
c		write(11,*) sgnl(i),rng(i)
78		continue
	do 745 i=1,9
	mat(i)=0.
	topdep(i)=0.
	tmeanx(i)=0.
745 	continue
c.....................................
c
c	
c	load lidar data and check saturation flag
c	
	if (npts.le.100) then
		read(10,209,end=1999) (sgnl(i),i=1,npts)
	else
		read(10,2209,end=1999) (sgnl(i),i=1,npts)
	end if
209		format(100f7.0)
2209		format(200f7.0)
		read(22,365) sat
365		format(13x,a1)
c		write(*,'(a\)') '\lidar data read, saturation flag checked'
c		write(*,209) (sgnl(i),i=1,npts)
c		write(*,365) sat		
c...........................................................
c
		radd=high
c
                do 79 i=1,npts
c
c	remove range correction!
c	RIVM data is range corrected ie.= P(r)*Index of r *some scaling factor
c
	if (rrcor.eq.1) then
	sgnl(i)=sgnl(i)/((i*1.)**2)
	else
	end if
c...........................................................
c
c	check for zero's or negatives -find lowest value!
c
			if (sgnl(i).lt.radd) then
			radd=sgnl(i)
			else 
			end if 
                        if (sgnl(i).lt.zero) then
                        sgnl(i)=zero
                        else
                        end if
c		write(*,'(a\)') '\checked for zeros'
c
c...........................................................
c	
c	take the log of the signal for the cloud base algorithm
c		
                if (sgnl(i).le.0.) then
                logsgnl(i)=0.
		else
                logsgnl(i)=alog(sgnl(i))
		if (logsgnl(i).lt.0.) logsgnl(i)=0
		end if
c...........................................................
c
c		mutiply data by constant factor

c               write(*,*) 'sgnl before multi',sgnl(i)
		sgnl(i)=sgnl(i)*multi
c		write(*,*) 'sgnl after multi',sgnl(i)
79              continue
		if (radd.gt.0.) then 
		write(*,*) 'No negative signal values'
		radd=0.
		else 
		end if
c                write(*,*) 'multi*radd=',(radd*multi)
c...........................................................
c
c		add constant background
c
		if (add.eq.1) then
		do 779 i=1,npts
c		write(*,*) 'sgnl before:', sgnl(i)
		sgnl(i)=sgnl(i)-(multi*radd)
c		write(*,*) 'sgnl after addition:', sgnl(i)
		sgnl(i)=sgnl(i)+offset
c		write(*,*) 'after offset:', sgnl(i)
                if (sgnl(i).lt.0) then
		write(*,*) 'choose offset and multi closer, rerun.' 
		stop
		else
		end if
779		continue
		else
		end if
c...................................................................
c
c	write non-range corrected signal to file
c	
4445	if (npts.le.100) then
	write(21,4442) (sgnl(j),j=1,npts)
	else
	write(21,24442) (sgnl(j),j=1,npts)
	end if
4442	format(100e10.3)
24442	format(200e10.3)
c...................................................................
c
c 	gain values and technical details
c
		info(1)=nos
		info(2)=rez
		info(3)=gain
		info(4)=nd
		jtrace = 1				
c
c
c..............................................................
c
c	* CLOUD BASE SECTION: write all cloud base data to file unit=11 *
c
c	call cloud base subroutine
c
c		write(*,'(a\)') '\calling cloud base subroutine'
c
	if (nlog.eq.1) then
c
c		write(*,*) 'log version'	
c
		call club(npts,logsgnl,INFO,corrfak,xovlp,
     +              	jbase,jpeak,jtop,jtrace,indx,mpt,dpt)
c..............................................................
c
	else
c
c	linear version
c
c	** not effective - not used  for 1064 nm **
c
c		write(*,*) 'linear version'
c
		call club(npts,sgnl,INFO,corrfak,xovlp,
     +              	jbase,jpeak,jtop,jtrace,indx,mpt,dpt)
c
	end if
c
c..............................................................
c
c 	write to file
c
c     		write(11,*) 'number of cloud layers',indx
		if (indx.eq.0) then
		write(*,'(A\)') '\NO CLOUD LAYERS DETECTED:'
	        write(11,747) (mat(i),i=1,9)
		do 7 i=1,3
		dpp(i)=0.
		mnx(i)=0.
		cldthk(i)=0.
		mdht(i)=0.
7		continue		
		pyt=20
		tmp(pyt)=0.
		rh(pyt)=0
c
c	detour to smoothing routines
c
		goto 21115
		else 
		end if
		if (indx.gt.3) then
		write(*,'(A\)') '\MORE THEN 3 LAYERS DETECTED:'
		write(*,'(A\)') '\DISREGARDING ANY LAYERS ABOVE THREE.'
		indx=3
		else 
		end if
c
	li=1
	do 746 i=1,indx
	mat(li+i-1)=rng(jbase(1,i))
	mat(li+i)=rng(jpeak(1,i))
	mat(li+i+1)=rng(jtop(1,i))
	li=li+2
746	continue
	do 7746 i=1,indx	
	cldthk(i)=rng(jtop(1,i)-jbase(1,i))
	mdht(i)=rng(((jtop(1,i)-jbase(1,i))/2)+jbase(1,i))
7746	continue
c	
       write(11,747) (mat(i),i=1,9)
747	format(9f7.0)
c 
c..............................................................
c
21115	if (asw.eq.1) then
c
c	fit an average line in between and above clouds down to confidence
c	height - go from top -> bottom.
c
c	setting up an alternate signal array
c
	do 312 me=1,npts
	fit(me)=sgnl(me)
312	continue
c
	do 6000 ls=npts,nint(conh/rez),-1
c
c	check if you're in a cloud, if so skip
c	until you hit the bottom then start again 
c
c	write(*,*) 'ls:',ls
	skipper=0
	if (indx.ge.1) then
	do 5993 blt=1,indx
	if ((ls.ge.jbase(1,blt)).and.(ls.le.jtop(1,blt))) then
	skipper=1
c	write(*,*) 'in-cloud skip'
	goto 6000
	else 
	end if
5993	continue
	else
	end if
c
c	are there clouds in this profile?
c
	if (indx.ge.1) then
c	write(*,*) 'cloud case:'
	toplen=0
	botlen=0
c
c	if yes, find ls with respect to clouds - ie. "case"
c
	m=indx
5000	if (ls.gt.jtop(1,m)) then
	goto 4999
	else
	m=m-1
	if (m.eq.0) then
	goto 5002
	else
	end if
	goto 5000
	end if
c
c	above top cloud case
c
4999	if (m.eq.indx) then
c	write(*,*) 'above top cloud case'
5001	if (((ls+toplen).lt.npts).and.(toplen.lt.avlen)) then
	toplen=toplen+1
	goto 5001
	else
	end if
c
	else 
c	write(*,*) 'general above cloud case'
c
c	general above cloud case
c

5002    if (((ls+toplen).lt.jbase(1,m+1)).and.(toplen.lt.avlen)) then
	toplen=toplen+1
	goto 5002
	else
	end if
c
	end if
c
c	below lowest cloud case
c
	if (m.eq.0) then
c 	write(*,*) 'below lowest cloud case'
5004	if (((ls-botlen).gt.nint(conh/rez)).and.(botlen.lt.avlen)) then
	botlen=botlen+1
	goto 5004
	else 
	end if
c
	else
c	write(*,*) 'general below cloud case'
c
c	general below cloud case
c	
5003	if (((ls-botlen).gt.jtop(1,m)).and.(botlen.lt.avlen)) then
	botlen=botlen+1
	goto 5003
	else
	end if
c
	end if
c	
	else
c
c	if no clouds, set botlen and toplen to avlen and go from top to
c	bottom limits  
c
c	write(*,*) 'no cloud case'
	botlen=avlen
	toplen=avlen
c
c	do top or bottom limits exceed profile length?
c
5998	if (ls+toplen.gt.npts) then
	toplen=toplen-1
	goto 5998
	else
	end if
5994	if (ls-botlen.lt.nint(conh/rez)) then
	botlen=botlen-1
	goto 5994
	else
	end if
c
c	close if loop for cloud vs. no cloud case   
c
	end if
c
c	SMOOTHING PROCEDURES:
c
c	choose either least squares fit over avlen points on either side 
c	of ls...
c
	if (fav.eq.1) then
c
c        write(*,*) 'line fit chosen'
	if ((botlen+toplen).eq.1) then
	write(*,*) 'no fit - single point'
	goto 6000
	else
	call linfit2(sgnl,rez,npts,ls,botlen,toplen,a,b)
	fit(ls)=a*ls*rez+b
c	write(*,*) 'fit,a,b',fit(ls),a,b
	end if
c
c	...or plain good old fashioned averaging of the avlen
c 	points on either side of avlen 
c	
	else
c
	avy=0.0
	divlen=1.
	do 5999 me=(ls-botlen),(ls+toplen),1
	avy=avy+sgnl(me)
5999	continue
	divlen=toplen+botlen+1.0	
	fit(ls)=avy/divlen	
c
	end if
c
c
c	write(*,*) 'toplen,botlen,divlen',toplen,botlen,divlen
6000	continue
c
c	transfer "fitted" or averaged data back into signal array to invert
c
	do 412 me=1,npts
	sgnl(me)=fit(me)
412	continue
c
c	for no cloud case don't proceed to invert, just write to file.
c
	if (indx.eq.0) then
	goto 1115
	else 
	end if
c
	else
	end if
c..............................................................
c
c	* OPTICAL SECTION:  write all optical data to file unit=16 * 
c
c	NOTE: 
c	load sonde data (look for tabulated data, always starts at pr=1000
c			 and is preceeded by 3 empty lines)
c
c
                if (isat.eq.1) then
		if (sat.eq.'j') then
		write(*,'(A\)') '\TRACE HAS SATURATED DATA - NOT INVERTING:'
		do 367 i=1,3
		dpp(i)=0.
		mnx(i)=0.
		cldthk(i)=0.
		mdht(i)=0.
367		continue		
		pyt=20
		tmp(pyt)=0.
		rh(pyt)=0
		goto 1115
		else 
		end if
		else
		end if
c..............................................................
c
c	read sonde data
c
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
	end if
	goto 56
c..................................................................
c
c	calculate rayliegh value
c
c	** Call rayliegh profile subroutine - modified for use with
c	crosub.for
c
58		call rayshiv(ntot,tmp,pr,sigr,br,rayn,lambda)
c
c..................................................................
c	add aerosol correction
c
c	load model and add to above profile
c	[Editors note: not yet implemented - P. Stammes hesitant
c	to part with personal secret knowledge of the sky.]
c..................................................................
c	** ideally use 1 or 2 range bins above 
c	cloud top as reference extinction
c	value. ie. sigref, iref and iref2 **
c..................................................................
c
c	 CORRELATE AND SET REFERENCE HEIGHT
c
c	go to default reference height or.... 
c
	if (seth.eq.1) then
	iref=nint(conh/rez)
c
c	check if this height is in a cloud
c
	do 7000 lr=1,indx
        if ((iref.ge.jbase(1,lr)).and.(iref.le.jtop(1,lr))) then
c	
c	if so move it down below nearest cloud base
c
	iref=jbase(1,lr)-1
	else 
	end if
7000	continue
	else
c	
c	otherwise use highest cloud top plus headroom.
c
	mink=iabove
c
c	lower reference height to 4000 km if higher
c
70	if ((jtop(1,indx)+mink).gt.npts) then
	mink=mink-1
	goto 70
	else	
	end if
c
	iref=jtop(1,indx)+mink
	end if
c
c	match it with the sonde
c
69	diff=abs(rng(jtop(1,indx)+mink)-ht(1))
	do 984 it=2,ntot
	toll=abs(rng(jtop(1,indx)+mink)-ht(it))
	if (toll.gt.diff) then
	goto 984
	else
	end if
	if (diff.gt.toll) then 
c
c	reference height correlated, match with sonde...
c 
 	diff=toll
	iref2=it
	else
	end if
984	continue
c
c	NOTE:sigref is extinction ref. coeff. in 1/km
c
c	If it's the same point take it -> otherwise:
c	Use line fit to interpolate for a ref. ext. value
c	
	htest=rng(jtop(1,indx)+mink)
	if (ht(iref2).eq.htest) then
	sigref=sigr(iref2)*1000.
	goto 187
	else
	end if
	if (ht(iref2).gt.htest) then
	itopx=iref2
	ibotx=iref2-1
	else
	itopx=iref2+1
	ibotx=iref2
	end if
	slop=(sigr(itopx)-sigr(ibotx))/(ht(itopx)-ht(ibotx))
	bay=sigr(itopx)-slop*ht(itopx)
c
c	plusig is a constant added by the user
c
	sigref=(slop*htest+bay)*1000.+plusig
c	write(*,*) 'ht-top, ht-bot, ref,sigref,matched sigref'
c	write(*,*) ht(itopx),ht(ibotx),htest,sigref,sigr(iref2)*1000.
c....................................................................
c
c	calculate optical extinction
c
c	call inversion subroutine
c
c	write(*,*) ' choosing sides'
c
187	if (iwk.eq.1) then
c
c	use Wolfi's Inversion
c
	write(*,*) 'Warning: MWS Inversion selected!'
c
c	much pain and misery caused by this algorithm -
c	somewhat like the namesake himself.
c
	call invert(sgnl,rez,rng(jtop(1,indx)+mink),
     + sigref,sig,jtop(1,indx)+mink,npts,inpts,k)
c....................................................................
c
	else if (iwk.eq.2) then
c
	write(*,*) 'Warning: Slope method inversion selected!'
c
c	or can use slope method for another opinion
c
	call slopsig(npts,rez,rng,sgnl,sig,sigref)
c
c....................................................................
c
	else 
c	
c	use Klett Inversion
c
	call klett(npts,rng,sgnl,sigref,iref,rez,
     +k,sig,inpts,mink,zero,high,baseline)
c
c
	end if
c
c	write(*,*) 'inverted'
c	write(*,*) calculate optical depth
c
c....................................................................
c
c	call optical depth subroutine
c
	do 1110 i=1,3
        dpp(i)=0.
        mnx(i)=0.
1110	continue
        do 1111 i=1,indx
	z0=rng(jbase(1,i))/1000.
	zt=rng(jtop(1,i))/1000.
	rezkm=rez/1000.
        call opdep(npts,rezkm,z0,zt,sig,dep,meanx)
        dpp(i)=dep
        mnx(i)=meanx
1111    continue
c	write(*,*) 'depthed'
c
c..............................................................
c
c	Truth Inversion Data using ratio method:
c
c	write(*,*)'going to find the truth...'
c
	call scrat(rez,npts,sgnl,indx,jbase,jtop,topdep,tmeanx,tpt,
     +zero)
c	write(*,*) 'truth found'
c..............................................................
c
c
c	correlate to cloud base or mid-cloud or top height
c
	if (mixbit.eq.1) then
	diff=abs(rng(jbase(1,1))-ht(1))
	else if (mixbit.eq.2) then
	diff=abs(rng(jtop(1,1))-ht(1))
	else
	diff=abs(mdht(1)-ht(1))
	end if
	do 7984 it=2,ntot
	if (mixbit.eq.1) then
        toll=abs(rng(jbase(1,1))-ht(it))
	else if (mixbit.eq.2) then
        toll=abs(rng(jtop(1,1))-ht(it))
	else
        toll=abs(mdht(1)-ht(it))
	end if
	if (toll.gt.diff) then
	goto 7984
	else
	end if
	if (diff.gt.toll) then 
 	diff=toll
	pyt=it
	else
	end if
7984	continue
c
c..............................................................
c
c	write optical data to files 16 and 17
c
1115	if (npts.le.100) then
	write(16,1112) (sig(i),i=1,npts)
	else
	write(16,21112) (sig(i),i=1,npts)
	end if
1112	format(100e10.3)
21112	format(200e10.3)
c....................................................
c	write processed-averaged signal to file 23
c
	if (npts.le.100) then
	write(23,1112) (sgnl(i),i=1,npts)
	else
	write(23,21112) (sgnl(i),i=1,npts)
	end if
c.....................................................
        write(17,1114) indx,(dpp(i),i=1,3),(mnx(i),i=1,3)
c
c	truthing data file
c
        write(20,1114) indx,(topdep(i),i=1,3),(tmeanx(i),i=1,3)
c
1114	format(i2,3e10.3,3e10.3)
c.............................................................
c
c	write correlated data to file 18
c
c	mid cld ht., geo. thk, opdep, mean ext., temp, rh
c
	if (mixbit.eq.1) then
	write(18,3) rng(jbase(1,1)),cldthk(1),topdep(1),mnx(1),
     +tmp(pyt),rh(pyt)
	else if (mixbit.eq.2) then
	write(18,3) rng(jtop(1,1)),cldthk(1),topdep(1),mnx(1),
     +tmp(pyt),rh(pyt)
	else
	write(18,3) mdht(1),cldthk(1),topdep(1),mnx(1),tmp(pyt),rh(pyt)
	end if
3	format(f7.0,1x,f7.0,1x,e10.3,1x,e10.3,2x,f8.1,1x,I4)
c
c..............................................................
c
c	write optical extinction data to file 19 in Andre's format  
c
	if (npts.le.100) then
	write(19,7772)  (mat(i),i=1,9),(sig(i),i=1,npts)
	else
	write(19,27772)  (mat(i),i=1,9),(sig(i),i=1,npts)
	end if
7772	format(9f7.0,100e10.3)
27772	format(9f7.0,200e10.3)
c..............................................................
c
c	next file
c
200		continue
c..............................................................
c
1999    write(*,*)'ALL TRACES ANALYZED - TASK COMPLETE.'
	write(*,*)
	write(*,*) '<ENTER> to return to system'
	read(*,*)
c
c	close all files
c
		close(8)
		close(9)
		close(10)
		close(11)
		close(12)
		close(14)
		close(16)
		close(17)
		close(18)
		close(19)
		close(20)
		close(21)
		close(22)
		close(23)
c
                stop
                end
c
c..................................................................
c
c	subroutines
c
c...........................................................

	subroutine klett(npts,range,sgnl,sigref,iref,rez,
     +k,sig,inpts,mink,zero,high,baseline)
c
c		klett.for
c
c		 A simple stable LIDAR aerosol
c		inversion FORTRAN subroutine source code
c		based on J.D. Klett, 15 Jan. 1981/
c		Vol. 20, No: 2/ Applied Optics 
c	
c	npts = number of total points
c	range  = array for range bins
c	sgnl = raw signal array, no correction or logging
c	sigref = refernce extinction value 
c		- calculate from Rayleigh + aerosol model (Shettle & Fenn)
c	iref = index number for reference atmosphere height
c	rez = range bin resolution
c	k = backscatter to extinction ratio constant
c	sig = array with extinction profile
c	mink = number of points above cloud top to ref height
c
c.................................................................
c
c
		real Z,REZ,ZREF,SIGREF,TM1,TM2,SUM,K
		real RANGE(200),SGNL(200)
		real S(200),SIG(200)
		real SIGMIN,SIGMAX
		real ZERO,SREF,div,high
		integer baseline
		INTEGER I,L,IREF,NPTS,mink,inpts,itiph,itipl
c	
c	convert everything to km:
c
	rez=rez/1000.
	do 8 i=1,npts
	range(i)=range(i)/1000.
8 	continue	
c	write(*,*) 'range converted'
c
	sref=0.	
	itipl=inpts
	itiph=inpts
24      if ((iref+itiph).gt.npts) then
	itiph=itiph-1
	goto 24
	else 
	end if
	div=(itiph+itipl)*1.+1.
c	div=itiph*1.+1.
c	write(*,*) 'itiph,itipl,div:',itiph,itipl,div
	do 57 i = (iref-itipl),(iref+itiph)
c	do 57 i = iref,(iref+itiph)
	SREF= SREF +sgnl(i)
57	continue
c 	write(*,*) 'sref=',sgnl(iref),sref
	if (baseline.eq.1) then
	sref=zero
	div=1.
	else 
	end if
	S(iref)=ALOG((RANGE(iref)**2)*(sref/div))
c
c		Logarithmic Range Adjusted Power
c
C
		DO 50 I=1,iref-1
c			write(*,*) 'sgnl',sgnl(i)
			S(I)=ALOG((RANGE(I)**2)*SGNL(I))
50		CONTINUE
c		write(*,*) 'all data logged'
c
c		The Klett inversion
c
888		SUM=0.0
		SIG(IREF)=SIGREF
c		write(*,*) 'calculate terms '
		DO 60 I=IREF-1,1,-1
			TM1=EXP((S(I)-S(IREF))/K)
			TM2=EXP((S(I+1)-S(IREF))/K)
			SUM=SUM+(TM1+TM2)/2.0
c		write(*,*) 'terms calculated-test'
		IF (((1.0/SIGREF)+(2.0*SUM*REZ/K)).EQ.0.) THEN
		SIG(I)=high
		write(*,*) 'Warning: denominator explosion!'
		ELSE
c		write(*,*) 'ratio!'
c		write(*,*) tm1
c		write(*,*) ((1.0/SIGREF)+(2.0*SUM*REZ/K)) 
		SIG(I)=TM1/((1.0/SIGREF)+(2.0*SUM*REZ/K)) 
		END IF
60		CONTINUE
c		write(*,*) 'inversion done'
c
c
c	
c	convert  back to m:
c
	rez=rez*1000.
	do 4 i=1,npts
	range(i)=range(i)*1000.
4 	continue	
		return
		end
c
c.................................................................
c
	subroutine slopsig(npts,rez,rng,sgnl,sig,sig0)
c
c	slopsig.for
c
c	program to invert lidar returns
c	for volume extinction coefficients (1/km)
c	using the slope method
c
c	rez: height resolution in metres
c	rng: array with range
c	sgnl:  "    "   signal
c	sig:   "    "   sigmas that are returned
c	npts: number of data points
c	sig0: reference extinction value at npts
c
c	arrays and variables
c	
	real rng(200),sgnl(200),S(200),sig(200)
	real sig0,rez
	integer npts
c
c	inversion
c
	sig(npts)=sig0
	S(npts)=alog(((rng(npts)/1000.)**2)*sgnl(npts))	
	do 1 i=npts-1,1,-1
	S(i)=alog(((rng(i)/1000.)**2)*sgnl(i))
	sig(i)=(S(i+1)-S(i))/(-2.0*(rez/1000.))
1	continue
c
c	end
c
	return
	end

c
c.................................................................
c
	subroutine rayshiv(ntot,temp,prez,sigr,br,rayn,lambda)
c		rayshiv.for	
c	A variation S.R. Pal's program for generating a number density 
c 	reference extinction and backscatter coefficient 
c	from temperature and pressure known as raysig.for
c
c	MODIFIED FOR USE WITH CROSUB.FOR
c	
c-------RAYSIG.FOR FOR CALCULATING RAYLEIGH NUMBER DENSITY FROM
c       PRESSURE AND TEMPERATURE  AND RAYLEIGH ATTENUATION 
c	X1(I)=ALOG10(X(I))
c	Y1(I)=ALOG10(Y(I))
c 
	DIMENSION RAYN(50),SIGR(50),BR(50)
	REAL PREZ(100),TEMP(100)
	REAL WTM,RGAS,AVOG,OMEGAR
	INTEGER NTOT
C
	WTM=28.9644
	RGAS=8.31432
	AVOG=6.022E+23
 	F=AVOG*100./RGAS
	call crosub(lambda,omegar)
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
c
c....................................................
c
		subroutine crosub(lambda,omega)
c
c	crosub.for
c
c	a subroutine to calculate the Rayleigh
c       backscatter cross-section as a function of wavelength.
c
c	arrays and variables
c
        real omega,lambda,x 
c
c	calculate
c
        if (lambda.le.0.55) then
	x=0.389*lambda+(0.09426/lambda)-0.3228	
	else
	x=0.04
	end if
        omega=(4.02e-28)/(lambda**(4.+x))
c
c       convert to m2
c
        omega=omega/(10**4)
c
	return
	end
c......................................................................
	subroutine scrat(rez,npts,sgnl,nlyr,jbase,jtop,topdep,tmeanx,tpt,
     +zero)
c
c	scrat.for
c
c	Cirrus Technology Consulting 
c	copyright 1994
c
c	Author: Alexandre Y. Fong
c	Last Revision/Compile/Update: June 30, 1994 [1:10 p.m.]
c
c
c.............................................................
c
c	Description:
c
c	This subroutine calculates the optical depth and mean
c	extinction coefficient of cloud signals by taking the
c	ratio of the lidar signal just above and below the clouds
c	de-range correcting dividing it by 2 and logging it.
c
c.............................................................
c
c	Arrays Variables and Constants
c
c	rez: per bin range resolution of data
c	npts: total number of points of data
c	sgnl: array with the signal
c	nlyr: number of cloud layers in this trace
c	jbase(1,*): array with cloud base height indices
c	jtop(1,*): array with cloud top height indices
c	range(*): array with range bins
c	topdep: array with cloud optical depths
c	tmeanx: array with cloud mean optical extinctions
c	tpt: number of points averaged over path
c
	real rez,sgnl(200),topdep(9),tmeanx(9),range(200)
	real rsc,rsq,rsqb,rsqt,opy,prod,topav,basav,div1
	real div2,zero
	integer npts,nlyr,tpt,itt,itb,itest
	integer*2 jbase(1,6),jtop(1,6)
	
c	intialize arrays
c
	do 2 k=1,6
	topdep(k)=0.
	tmeanx(k)=0.
2	continue
	do 3 l=1,npts
c
c	distances in km
c
	range(l)=l*rez/1000.
c	write(*,*) range(l)
3	continue	
c
	if(nlyr.eq.0) then
	goto 9
	else
	end if
	do 1 i=1,nlyr
	basav=0.
	itb=tpt
c	write(*,*)'base-2-lo?'
	if (nlyr.eq.1) then
8	itest=jbase(1,1)-(itb)
	if (itest.le.0) then
	itb=itb-1
	goto 8
	else
	end if
	else
7	itest=jbase(1,i)-itb
	if (itest.lt.jtop(1,i-1)) then 
	itb=itb-1
	goto 7
	else
	end if
	end if
	do 6 m=1,itb+1
c	write(*,*)'m=',m
	basav=basav+sgnl(jbase(1,i)-(m-1))
6	continue
	div1=itb+1
	basav=basav/div1
c	write(*,*)'top-2-hi?'
	topav=0.
	itt=tpt
5	itest=jtop(1,i)+itt
	if (itest.gt.npts) then 
	itt=itt-1
	goto 5
	else
	end if
c	write(*,*)'truth-top...'
	do 4 m=1,itt+1
	topav=topav+sgnl(jtop(1,i)+(m-1))
4	continue
	div2=itt+1
c	write(*,*) 'div2=',div2
c	write(*,*) 'div1=',div1
c	write(*,*)'topav=',topav
c	write(*,*)'looking for truth'
	if (topav.eq.0.) then
	rsc=0.
	else
	topav=topav/div2
	rsc=basav/topav
	end if
	rsqb=(range(jbase(1,i)))**2
	rsqt=(range(jtop(1,i)))**2
	rsq=rsqt/rsqb 
	prod=rsc*rsq
	if (prod.le.0.) then
	opy=0.
	else
	opy=alog(prod)
	end if
	if (opy.lt.0.) then 
	topdep(i)=0.
	else
	topdep(i)=opy/2.
	end if
	thk=range(jtop(1,i))-range(jbase(1,i))
	if (thk.eq.0.) then
	tmeanx=0.
	else
	tmeanx(i)=topdep(i)/thk	
	end if
1	continue	
c	write(*,*) 'truth found!'
9	return
	end

c......................................................................
c
	subroutine opdep(npts,rezkm,z0,zt,sig,dep,meanx)
c
c	a subroutine to calculate the optical depth and mean cloud extinction
c	from an exitinction profile from a given cloud base and top altitude,
c	vertical resolution and...
c
c	npts = number of data points
c	rez = vertical resolution
c	zmin = overlap altitude of LIDAR
c	z0 = cloud base
c	zt = cloud top
c	sig = extinction profile
c	dep = optical depth
c	meanx = mean cloud extinction (as defined by C.M.R. Platt,
c		Journ.Atmos.Sci. Feb 1987)
c
c	arrays,integers and real variables
c
        real rezkm,z0,zt,sig(200),dep,sum,meanx
	integer npts,i,it,i0
c
	i0=ifix(z0/rezkm)
	it=ifix(zt/rezkm)
c
	dep=0.00
	sum=0.00
c
	do 1 i=i0,it
	sum=sig(i)*rezkm
	dep=dep+sum
1	continue
c
        if ((zt-z0).eq.0) then
	meanx =0.
	else
	meanx=dep/(zt-z0)
	end if
c
	return
	end
c..................................................................
c
c	more subroutines
c..................................................................
c
	SUBROUTINE club(N,P,INFO,corrfak,xovlp,
     +              	jbase,jpeak,jtop,jtrace,k,npt,dpt)
C..................................................................
C
c	subroutine to calculate cloudbase and cloudtop for a single trace
c
c	parameters: integer*2 n = number of data points in the trace
c		    real      p(*) = array containing the lidar backscatter 
c		    		     signal
c	            real info(*) = array containing additional information
c					about the lidar trace, like time, date,
c					neutral-density filter, height-resolution
c
c			info(1) = number of shots averaged for this trace (=fshot)
c			info(2) = resolution
c
c			info(3) = gain
c		         
c			info(4) = overall neutral density
c
c		    real corrfak = faktor that determines the sensitivity of 
c				   algorithm. low corrfak tends to trigger on
c				   noise spikes, high corrfak will not recognize 
c				   weak cloud layers. good values seem to be 
c				   for 150 shot averages and the green signal:
c				      corrfak = 1.4 ... 1.8
c					
c				   and the IR signal :
c	                               corrfak = 0.8 ... 1.2
c
c		    real xovlp = overlap height. signal peaks below the overlap height
c				 will not registre clouds.	
c
c		    integer*2 jbase(600,*),jpeak(600,*),jtop(600,*) the 
c				 	cloudbase algorithm will report the
c				        index of the altitudes (real altitude =
c					= index * height resolution) in these two 
c				        dimensional arrays. since only one trace 
c				        is processed by icb, the integer*2 quantity
c				        jtrace tells icb into which column to store 
c				        the bases and tops for the present trace.
c					jbase(jtrace,1) will contain the index of the 
c					base of the lowest cloud, jbase(jtrace,2) the 
c					index of the next higher cloud and so on.
c					before I use these arrays I fill them up with
c					-1s, since icb only writes the bases it  finds 
c					to the array. therefor jbase(jtrace,level) = 0 means
c					that no base for the level = level was found for trace 
c				        jtrace.
c
c		integer npt = number of point search for peak 
c		integer dpt = number of point fit for derivative fit
c
c		     integer*2 jtrace = number of the trace being processed.
c
c		     integer*2 indx = used to give the plotting coordinate on the 
c				      time axis for the time - altitude plots.
c
c
c		
C................................................................
c
c	arguments passed in parameter list
c
	real p(*),xovlp,corrfak
	integer*2 jbase(1,*),jpeak(1,*),jtop(1,*)
        integer*2 jtrace,mode,iwidth,n
	integer*4 npt,dpt
	real INFO(*)
	integer indx
c
c
	INTEGER*4 	ihelp
	INTEGER		izero(400),ipeak(400),itop(200)
	INTEGER		i,j,k,l,ILAM,fshot
 	REAL		PZERO(400),PPEAK(400),ptop(200)
	REAL		PPRM(2050),f1hlp(509),f2hlp(509),f3hlp(2044)
	real		f4hlp(2044)
        integer*4       i1hlp(3),i2hlp(3),i3hlp(4),i4hlp(4)
        real             xres,nd0,fnd
	
c	common /help2/ f1hlp,f2hlp,f3hlp,f4hlp,i1hlp,i2hlp,i3hlp,i4hlp
c     + 			,izero,ipeak,itop,pzero,ppeak,ptop	
c
c
C..................................................................
c
c	get some information
c
	  xres=info(2)
	  gain=info(3)
	  nd0=info(4)
	  fshot=info(1)	
	  fnd=nd0
c
c.....................................................................
c	get the first derivative of the signal using a dpt point 
c	least squares linear fit
c
c	  p(1)=p(2)
c
	call difflin(p,xres,n,pprm,dpt)
c
c
C------------------------------------------------------------------
C  --- Check on Zero Crossing---------------
C------------------------------------------------------------------
c
	K=0
	ibgn=xovlp/xres
c	old
	IEND=n-1
c	new
c	IEND=n
	npt2=npt/2
c
c	*SLOWPITCH*
c	
	DO 120 I=IBGN,IEND
c	original
c	  IF((PPRM(I).LT.0.0).AND.(PPRM(I+1).GE.0.0)) then
c	my way
	  IF((PPRM(I).le.0.0).AND.(PPRM(I+1).gt.0.0)) then
c
c	that's a base candidate
c
	if(k.gt.400) then
	write(*,*) 'stop idiot!',k
	stop
	else
	ENDIF
	    k=k+1
	    pzero(k)=p(i) 
	    izero(k)=i
c
c	if there is a smaller value in the neighbourhood-
c	thats the better base-candidate
c 
	    do 100 j=max0(ibgn,i-npt2),min0(iend,i+npt2)
	      if (p(j) .lt. pzero(k)) then 
		pzero(k)=p(j)
	        izero(k)=j
	      endif
100	    continue      
	  ENDIF
c
c........................................................................
c
	if (k.gt.400) then
	write(*,*) 'STOP STUPID!',k
	stop
	else
	endif   
	    ipeak(k)=izero(k)
	    ppeak(k)=pzero(k)
c........................................................................
c	now look for the next peak
c
	    do 11110 jj=izero(k),iend
	      IF((PPRM(jj).GT.0.0).AND.(PPRM(jj+1).LE.0.0))then
c
c	that's a peak candidate
c
	        ppeak(k)=p(jj+1) 
	        ipeak(k)=jj+1
c
c	if there is a higher value in the neighbourhood-
c	however after the corresponding base candidate
c	thats the better peak-candidate
c 
	        do 110 j=izero(k)+1,min0(iend,jj+1+npt2)
	          if (p(j) .gt. ppeak(k)) then 
		    ppeak(k)=p(j)
	            ipeak(k)=j
	          endif
110	        continue      
	        goto 119
	      ENDIF
11110	    continue
c		new	
119	  if (k .ge. 400) goto 130
c
120	CONTINUE
	goto 135 
c
130	pause 'more than 400 base candidates found'
c
135	IF(K.EQ.0) GOTO 80
c
C
c...................................................................
C	Reject all insignificant peaks and corresponding zero crossins 
C
	  fak=sqrt(87./fshot/xres*(10.**fnd))*1.16
c
	DO 200 J=1,K
 	    REJECT=PPEAK(J)
	    REJECT=SQRT(REJECT)*FAK*corrfak
	    if (izero(j) .lt. nint(4000./xres)) reject=reject*2.0
	    if (izero(j) .lt. nint(2000./xres)) reject=reject*1.5	 
c
c	mark the base and peak candidates to be rejected
c
	  DIFF= PPEAK(J)-PZERO(J)
	  IF (DIFF.le.REJECT) then
c
c	  single peak less than reject might be thrown away
c
	    diff21 = p(ipeak(j)-1) - pzero(j)
	    diff22 = p(ipeak(j)+1) - pzero(j)
	    reject = reject/1.5
c
c	   unless the single peak and a neighbour-point are
c	   larger than a reduced reject
c
	    if ((diff .gt. reject).and.((diff21.gt.reject).or.
     +                                  (diff22.gt.reject))) goto 200 	
	    diff31 = p(ipeak(j)-2) - pzero(j)
	    diff33 = p(ipeak(j)+2) - pzero(j)
	    reject = reject/1.22
c
c	   or unless the single peak and two neighbours are
c	   larger than an even further reduced reject 
 	    
	    if ((diff .gt. reject).and.
     +	        (((diff31.gt.reject).and.(diff21.gt.reject)).or.
     +		 ((diff22.gt.reject).and.(diff33.gt.reject)))) goto 200
c
c	   now it has to be rejected
c 
	    	izero(j)=-1
c
	  endif
c
200	CONTINUE
c
c....................................................................
c	now throw all the bad candidates away
c
	j=1
202	  if (izero(j).eq.-1) then
c	throw the j.th point away, shift all other points left,
c	after that look at the j.th point (now new,because shifted) again	
	    k=k-1
	    call wshift(izero,pzero,j,k,1)
	    call wshift(ipeak,ppeak,j,k,1)
	  else
c	look at the next point
	   j=j+1
	  endif
c
c	the thing is finished when we have looked at the k.th point
c
	if (j.le.k) goto 202
c	
	if (k .le. 0) goto 80
c
c...................................................................
c
C-----------------------------------
C   Determine cloud top
C-----------------------------------
c
	pback=0.
	pback=p(2043)
	j=1
4999        do 5100 i=ipeak(j),n
c
c	    thats the simplest cloud top: signal at top <= signal at bottom 
c
	      if (p(i) .le. pzero(j)) then
	        itop(j)=i
	        ptop(j)=p(i)
cc	        write (*,'(a\)') ' top found for candidate',j,'itop=',i 
cc	        write (8,*) ' top found for candidate',j,'itop=',i 
c
c	    it is better if we can find a top were p(top) <= p(base) range 
c								 corrected
c 
	        if (i*xres .lt. 7000) then
	        do 5110 l=i+1,min0(n,i+nint(600./xres))
  		  rcorr=float(izero(j))/float(l)
		  rcorr=rcorr*rcorr
	          if (p(l) .le. (pzero(j)*rcorr+pback*(1.-rcorr))) then
	            itop(j)=l
	            ptop(j)=p(l)
CC	            write (*,'(a\)') ' range corrected',j,'itop=',l 
cc	            write (8,*) ' even better',j,'itop=',l 
		    goto 5120		              
	          endif
5110		continue	
	        endif 
5120	        j=j+1
	        goto 5000
	      endif
	      if (i .eq. n) then
cc	        write (*,'(a\)') ' no top found for candidate',j
cc	        write (8,*) ' no top found for candidate',j
c
c	throw the candidate away, shift the rest 
	        k=k-1
	        call wshift(izero,pzero,j,k,1)
	        call wshift(ipeak,ppeak,j,k,1)
	      endif
5100	  continue 
5000	if (j .le. k) goto 4999
c
c	now make shure to get the highest top of all bases in the same 
c	cloud
c
	j=1
3003	ishift=0	
	do 4000 i=j+1,k
	  if (izero(i) .le. (itop(j)+nint(100./xres))) then
	    if (itop(i) .gt. itop(j)) then
		itop(j)=itop(i)
	        ptop(j)=ptop(i)
	    endif 	
	    ishift=ishift+1
	  endif
4000	continue
	if (ishift .gt. 0) then
c
c	throw these ishift points away, shift all other points left,
c
cc	write (*,'(a\)') 'ishift=',ishift
cc	write (8,*) 'ishift=',ishift
	    k=k-ishift
	    call wshift(izero,pzero,j+1,k,ishift)
	    call wshift(ipeak,ppeak,j+1,k,ishift)
	    call wshift(itop,ptop,j+1,k,ishift) 
	endif
c
c	look at the next point
c
	   j=j+1
c
c	the thing is finished when we have looked at the k.th point
c
	if (j.le.k) goto 3003
c
c.........................................................................
c	find the maximum signal from the cloud
c
	do 4114 j=1,k
	  peak=ppeak(j)
	  imerk=ipeak(j)
	  do 4224 jj=izero(j),itop(j)
	    if (p(jj) .gt. peak) then
              peak=p(jj)
	      imerk=jj
	    endif
4224	  continue
	  ppeak(j)=peak
	  ipeak(j)=imerk
4114	 continue 	
c
c	now lets write the stuff to the cloudbase file
c
	if (k .gt. 6) stop 'icb cannot handle than more than 6 layers'
	DO 4230 J=1,k
	  jbase(jtrace,j)=izero(j)
	  jpeak(jtrace,j)=ipeak(j)
	  jtop(jtrace,j)=itop(j)
c 	  WRITE(11,370) J,izero(J)*xres,pzero(J),itop(j)*xres,ptop(j),      
c     +			ipeak(j)*xres,ppeak(j)
cc 	  WRITE(11,370) J,izero(J)*xres,pzero(J),itop(j)*xres,ptop(j)
4230	CONTINUE
c
c80     WRITE(11,390)
80	WRITE(*,391)
c80       WRITE(*,390)
C-------------------------------------------------------
C	FORMAT STATEMENTS
C............................................................
 300	FORMAT ('@Time:',i6,'Date:',i6,'lambda',
     +	i5,'ND',f4.1,'gain',f5.1,'NPT',i3,'@') 
 310	FORMAT(I5,F12.5,2E15.5)
 320	FORMAT('             ALL  ZERO CROSSING IN dP/dr' )
 330	FORMAT(' Index  Rbase [m]   Pbase       Rtop [m]    Ptop
     +	Rpeak [m]  Ppeak')
 340	FORMAT(2I6,E12.5,i6,E12.5)
 350	FORMAT('    Main zero crossings and cloud peaks ')
 370	FORMAT(I6,6E12.5)
 380	FORMAT('    Cloud bases of separate layers only ')
 390	FORMAT('   NO CLOUD DETECTION')
 391	FORMAT('   TRACE ANALYZED')
395 	format('@ymin=',e12.5,'ymax=',e12.5,'@')	
C-----------------------------------------------------
 400	CONTINUE
	RETURN
	END
c
c	more SUBROUTINES.....................................
c
c
	subroutine wshift (iarr,parr,j,k,ishift)
c
c	this subroutine shifts in the arrays iarr(integer*2)
c	and parr(real) all data points left. 
c	data point j+ishift is shifted into j
c       .....
c	data point k+ishift is shifted into k
c
c	data points with index < j are not affected.
c	data points with j =< index < j+index are lost.
c
c	! make shure k+ishift does not exceed the array-size ! 	 
c
	integer iarr(*),j,k,ishift
	real parr(*)
c
 	do 200 l=j,k
	  iarr(l)=iarr(l+ishift)
	  parr(l)=parr(l+ishift)
200	continue
c
	return
	end

c
c...................................................................
c
	subroutine difflin(fdat,dx,n,fdiff,iwidth)
c
c	takes the derivative df/dx of the function given by 
c	f(xi)=fdat(i) , xi=i*dx.
c	to get the derivative a straight line is fitted through the
c	data points iwidth to the left and the right of the point 
c	where the derivative is to be taken.
c	fdat(1..n) = array containing the values
c	dx = increase in x value for increasing the index
c	fdiff(1..n) = array returning the derivative
c	iwidth = width of the fit
c
	real fdat(*),fdiff(*),dx
	integer iwidth
        integer*2 n
c
	real fx,fy,fx2,fxy,a,b
	integer iold,i
c
	iold=0
	do 100 i=1,n
	 call linfit(fdat,dx,n,i,iwidth,a,b)
	 fdiff(i)=a
cc	 write (*,*) ' i,fdiff(i)',i,fdiff(i) 
100	continue
c
	return
	end
c
c....................................................................
c
	subroutine linfit(fdat,dx,n,icenter,iwidth,a,b)
c
c	fits a straight line a*i*dx+b through the data points given by
c	the coordinates fdat(i),i*dx. fdat(1..n) is a real-array. 
c	the line is fitted through the data points iwidth to the left
c	and to the right of icenter. 		
c	a,b are the coefficients of the linear fit.
c
c
	real fdat(*),dx,a,b
	integer icenter
	integer iwidth
        integer*2 n
c
c
	ileft=max0(1,icenter-iwidth)
	iright=min0(icenter+iwidth,n)
	di=iright-ileft+1 
	if (di .le. 1) goto 999
c
c	reset the other things
c	
	  fx=0.
	  fx2=0.
	  fy=0.
	  fxy=0.
c
c	evaluate the sum yi = y(ileft)+y(ileft+1)+ .. y(iright)
c	  and the sum xi*yi = dx*(ileft*y(ileft)+ .. +iright*y(iright))
c	 and the other sums 
c
105	do 200 i=ileft,iright
cc	  write (*,*) ' adding point ',i
	  fx=fx+i
	  fx2=fx2+float(i)*i 	
	  fhelp=fdat(i)
	  fy=fy+fhelp
	  fxy=fxy+i*fhelp            
200	continue
c
c	calculate the a and b
c
cc	write (*,*) ' fx=',fx
cc	write (*,*) ' fx2=',fx2
cc	write (*,*) ' fy=',fy
cc	write (*,*) ' fxy=',fxy
	fhelp=(fx*fx-di*fx2)*dx
	if (fhelp .eq. 0.0) fhelp=1.e-15
	a=(fx*fy-di*fxy)/fhelp
	b=(fy-dx*a*fx)/di
c
	return
999	pause 'cannot fit a line with only 1 data point'
	return
	end
c
c.........................................................................
	subroutine linfit2(fdat,dx,npts,icenter,botlen,toplen,a,b)
c
c	***MODIFIED by Cirrus Technology Consulting for KNMI
c	
c	original code by W. Steinbrecht of York-ISTS modified
c	for current application by A.Fong of Cirrus.
c
c	fits a straight line a*i*dx+b through the data points given by
c	the coordinates fdat(i),i*dx. fdat(1..n) is a real-array. 
c	the line is fitted through the data points botlen points below
c	and toplen points above of icenter. 		
c	a,b are the coefficients of the linear fit.
c
c
	real fdat(*),dx,a,b
	integer icenter
	integer toplen,botlen
        integer npts
c
c
	di=botlen+toplen+1
	if (di .le. 1) goto 999
c
c	reset the other things
c	
	  fx=0.
	  fx2=0.
	  fy=0.
	  fxy=0.
c
c	evaluate the sum yi = y(ileft)+y(ileft+1)+ .. y(iright)
c	  and the sum xi*yi = dx*(ileft*y(ileft)+ .. +iright*y(iright))
c	 and the other sums 
c
105	do 200 i=icenter-botlen,icenter+toplen
cc	  write (*,*) ' adding point ',i
	  fx=fx+i*dx
	  fx2=fx2+(i*dx)**2 	
	  fhelp=fdat(i)
	  fy=fy+fhelp
	  fxy=fxy+i*dx*fhelp            
200	continue
c
c	calculate the a and b
c
cc	write (*,*) ' fx=',fx
cc	write (*,*) ' fx2=',fx2
cc	write (*,*) ' fy=',fy
cc	write (*,*) ' fxy=',fxy
	fhelp=di*fx2-(fx**2)
	if (fhelp .eq. 0.0) fhelp=1.e-15
	a=(di*fxy-fx*fy)/fhelp
	b=(fx2*fy-fx*fxy)/fhelp
c	write(*,*) 'fit parameters derived'
c
	return
999	pause 'cannot fit a line with only 1 data point'
	return
	end
c
c.........................................................................
c
	subroutine invert(p,res,rfar,sigfar,sigma,ifar,npts,inpts,k)
c
c
c	***MODIFIED by Cirrus Technology Consulting for KNMI
c
c	original code by W. Steinbrecht of York-ISTS modified
c	for current application by A.Fong of Cirrus.
c
c	does a simple far end (klett) inversion,
c	The standart deviation of the 
c	sigma profile is estimated using the standart deviation of the signal
c	and gaussian error analysis.
c	Signal values less than zero are forced to zero.
c
c	input parameters : p = array containing the backscattered power, 
c			       background subtracted, but not range corrected
c			       and/ or log-ged.
c			   dp = array containing the standart deviation 
c				of the signal 
c			   res = range resolution [m]
c			   sigfar = far end starting sigma [1/km]
c			   rfar = far end distance [m]
c			inpts = number of point averaage of signal above cloud
c
c	output parameters : sigma = array containing the sigma profile [1/km]
c			    dsigma = array containing estimated standart 
c				     deviation of sigma	
c			    ifar = array index corresponding to the far end
c				   range
c
        real p(*),sigma(*)
c        real dsigma(*),dp(*)
        real sigfar,rfar
	integer ifar,inpts
c
	real pff,sum,sum2,pfold,sig2far,div
c
c	
	ifar = nint(rfar/res)
c	write(*,*) 'rfar,ifar',rfar,ifar
	pff = 0.
c
c	do 50 i = ifar-5,ifar+5
c50	pff = pff +p(i)
2	div=abs(inpts)*1.+1.
c	write(*,*) 'div=',div
        if ((ifar+inpts).ge.npts) then
	inpts=inpts-1
	goto 2
	else 
	end if
	do 50 i = ifar,ifar+inpts
50	pff = pff +p(i) 

c	far end range corrected power

C	force pff >= 0.
	if (pff .lt. 0.) pff=0.
51      pff = (pff/div)*ifar*ifar
c
c	terms for error analysis
c	
c        if (p(ifar) .ne. 0.) then
c          dx=dp(ifar)/p(ifar)/sqrt(11.)
c        else
c          dx=dp(ifar)/1.e-10/sqrt(11.)
c        endif
c        dy = ifar*ifar*dp(ifar)/sigfar/sqrt(11.)
c        sumz = ifar*ifar*dp(ifar)*dp(ifar)
c        dz = sumz*4.
c
c	terms for sigma calculation
c
	reskm = res/1000.
	sigma(ifar) = sigfar
	sum2 = 0.
	pr2old = pff
	sig2far = pff/sigfar
c
c        dsigma(ifar) = sigma(ifar)*sqrt(dx*dx+
c     +              (dy*dy+dz)/sig2far/sig2far)
c
	do 200 i = ifar-1,1,-1
c
c	that is to calculate the sigma
c
	   pff = p(i)*i*i
c	force pff >= 0.
	   if (pff .lt. 0.) pff = 0.
	   sum2 = sum2 + (pff+pr2old)*reskm
	   if ((sig2far + sum2).eq.0) then
	   sigma(i)=0.
	   goto 54 
	   else
	   end if
	   sigma(i) = pff/(sig2far + (sum2/k))
54	   pr2old = pff
c
c	   write (*,*) i*res,sigma(i)
c
c	that is to calculate the standart deviation of the sigma
c
c           if (p(i) .ne. 0.) then
c             dx=dp(i)/p(i)
c           else
c             dx=dp(i)/1.e-10
c           endif
c           sumz=sumz+i*i*dp(i)*dp(i)
c           dz = 4.*sumz/(ifar-i+1)
c           dsigma(i) = sigma(i)*sqrt(dx*dx+
c     +              (dy*dy+dz)/(sig2far+sum2)/(sig2far+sum2))
c
200	continue
c
	return
	end
c
c.......................................................................
c

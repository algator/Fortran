
	program sosi
c
c	SOSI.FOR
c
c	program to generate Si-So values
c	and invert them by the Klett method from
c	extinction profiles generated by raygen3.for
c 
c	arrays and constants
c
	real siso(100),range(100),sisf(100)
c        real P(600)
        real SIG(600),SIGk(600)
        real BR(600),KT(600)
        real rez
c        real start
        real k
	real SUM,TM,SUM2,TM2
	integer npts 
	character*60 infile,outfile

c...........................................................
c
c	inputs 
c
	write(*,'(a\)') 'name of ext. input file'
	read(*,'(a)') infile
	open(unit=9,file=infile,status='old')
	write(*,'(a\)') 'name of output file'
	read(*,'(a)') outfile
	open(unit=10,file=outfile,status='new')
	write(*,*) 'number of points?'
	read(*,*) npts
        write(*,*) 'increment?'
        read(*,*) rez
        write(*,*) 'k-ratio?'
        read(*,*) k
c	write(*,*) 'starting value?'
c	read(*,*) start
c...........................................................
c
c	do the what-sosi!
c
	sum=sig(1)
c	tm=-2.0*sum*(rez/1000.)
	siso(1)=0.
        read(9,2,end=3) range(1),sig(1),br(1),kt(1)
	do 1 i=2,npts
        read(9,2,end=3) range(i),sig(i),br(i),kt(i)
2	format(1x,f7.0,1x,e10.3,1x,e10.3,1x,e10.3)
	sum=sum+(sig(i-1)+sig(i))/2.
	tm=-2.0*sum*(rez/1000.)
	siso(i)=alog(sig(i)/sig(1))+tm
1       continue
c.............................................................
c
c	Time to re-Klett on some things...
c
c		The Klett inversion
c
888		sum2=0.
		SIGk(npts)=SIG(npts)
		sum=sig(npts)
c		tm=2.0*sum*(rez/1000.)
		sisf(npts)=0
c		write(*,*) 'calculate terms '
		DO 60 I=npts-1,1,-1
		sum=sum+(sig(i-1)+sig(i))/2.
		tm=-2.0*sum*(rez/1000.)
		sisf(i)=alog(sig(i)/sig(1))+tm
		TM2=EXP(sisf(i)/K)
		SUM2=SUM2+TM2/2.0
c		write(*,*) 'terms calculated-test'
		IF (((1.0/SIGk(npts))+(2.0*SUM*REZ/K)).EQ.0.) THEN
		SIG(I)=1.0e+20
		ELSE
c		write(*,*) ((1.0/SIGk(npts))+(2.0*SUM*REZ/K)) 
		SIGk(I)=TM2/((1.0/SIGk(npts))+(2.0*SUM*REZ/K)) 
		END IF
60		CONTINUE
c		write(*,*) 'inversion done'
c
c
c	
c.............................................................
c
c	Everyday, everyday I write the books...
c
1112    do 26 i=1,npts
	write(10,6) range(i),siso(i),sigk(i),sig(i),sisf(i)
6	format(1x,f7.0,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3)
26       continue
c.............................................................
c
c	that's all folks
c
3	close(9)
	close(10)
	write(*,*) '\DONE!'
	stop
	end
c.............................................................

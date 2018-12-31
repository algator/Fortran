	PROGRAM PGEN2
C	
C	PGEN2.FOR
C
C	a program to generate lidar returns (received power)
c	from extinction profiles.
c
c.............................................................
c
c	arrays and constants
c
	real P(600),SIG(600),RNG(600)
        real BR(600),mick,k(600)
	real SUM,TM,POA,CT2,REZ
	integer npts 
	character*60 infile,outfile
c.............................................................
c
c	inputs and outputs
c
	write(*,'(a\)') 'name of ext. input file'
	read(*,'(a)') infile
	open(unit=9,file=infile,status='old')
	write(*,'(a\)') 'name of output file'
	read(*,'(a)') outfile
	open(unit=10,file=outfile,status='new')
	write(*,*) 'number of data points'
	read(*,*) npts
	write(*,*) 'resolution?'
	read(*,*) rez
	write(*,*) 'boundary condition (start power)?'
	read(*,*) mick
	write(*,*) 'factor?'
	read(*,*) POA
c.............................................................
c	
c	set other constants 
c
	SUM=0.
	CT2=((3.0e+8)*(9.0e-9))/2.
c.............................................................
c
c	do it!
c
	P(1)=mick
	read(9,2,end=3) rng(1),sig(1),br(1),k(1)
	write(10,6) rng(1),p(1)
        do 1 i=2,npts
	read(9,2,end=3) rng(i),sig(i),br(i),k(i)
2	format(1x,f7.0,1x,e10.3,1x,e10.3,1x,e10.3)
	sum=sum+(sig(i-1)+sig(i))/2.
	tm=exp(-2.0*sum*(rez/1000.))
	p(i)=(CT2*k(i)*sig(i)*tm*poa)/((rng(i)/1000.)**2)
	write(10,6) rng(i),p(i)
6	format(1x,f7.0,1x,e10.3)
1       continue
c
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

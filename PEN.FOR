	program pen
c
c	pen.for
c
c	program to calculate theoretical
c	lidar penetration depth limits
c
c	arrays and variables
c
	real sig,rng,io,snr,rez,emin,sum
	integer npts,i
	character*60 filein,fileout
c	
c	intialize 
c
	rng=0.0
	sig=0.0
	rez=0.0
	snr=0.0
	io=0.0
	emin=0.0
	sum=0.0
	npts=0
c
c	input
c
	write(*,*) 'extinction profile file?'
	read(*,'(a)') filein
	write(*,*) 'beam output energy, signal-to-noise ratio?'
	read(*,*) io,snr
	write(*,*) 'number of points, height resolution?'
	read(*,*) npts,rez
	write(*,*) 'output file?'
	read(*,'(a)') fileout
c
c	open files
c
	open (unit=8,file=filein,status='old')
	open (unit=9,file=fileout,status='new')
c	
c	define SNR
c	

c
c	read-calculate-write	
c
	do 2 i=1,npts
	read(8,1) rng,sig
	sum=sum+sig*rez
	emin=io*(snr**2)*exp(-2.0*sum)
	write(9,1) rng,emin
1	format(f7.0,1x,e10.3)
2	continue
c
c	end of program
c
	write(*,*) 'program finished - <ENTER> to resume system'
	read(*,*)
	close(8)
	close(9)
	stop
	end

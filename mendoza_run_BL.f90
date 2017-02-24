!**************************************
!	
!	Mendoza model
!
!**************************************

program mendoza_run
	
	use omp_lib
	use mendoza_tools
	use basic
	use tools
	
	implicit none 
	 
	!**************************************
	! Initialize the parameters
	!**************************************
	
	integer :: model_t = 1					! RE Model = 0 , BL Model = 1

	real :: beta = 0.91, rho = 0.878, sigmay = 0.00663 
	real :: r      = 0.027, kappah = 0.93, kappal = 0.63, alpha  = 0.026
	character(100) :: fileplace
	character(100) :: filename
	
	!**************** Technical parameters  *****************
	real :: uptd	= 0.1 		! Weight on new policy function in update to next iteration
								! change to 0.1 or smaller (w # if policy function fluctuate)
	integer :: outfreq  = 5 	! Display frequency (shows in screen each 20th iteration)
	integer :: iter_tol = 3000  ! Iteration tolerance (max iterations)
	real 	:: gamma  = 2.
	integer :: l
	
	! Numerical tolerance (convergence criterion)
	real :: tol	= 1E-4 	! High value for session, increase for accuracy

	!*********************************************************
	! Grid assignment 
	integer, parameter :: NB   = 100,  NK = 2, NL = 44, NTOTAL = 18, NS = 9! number of nodes
	real :: bmin = -0.7, start, finish
	real :: bmax = -0.1
	
	! Probabilities 
    real :: Fhh 
    real :: Fll 
	
	real, dimension(:), allocatable 	:: FFHHHH,FFLLLL
	real, dimension(:,:,:), allocatable :: BBPP, QQ, EERR, PPRROOBB

    real :: kappa(NB,NTOTAL), B(NB,1), Z(NS), Zprob(NS,NS), prob(NTOTAL,NTOTAL)=0., er(NB,NTOTAL)
    real :: q(NB,NTOTAL), bp(NB,NTOTAL)
    
    real :: BB(NB), RR = 0.027
    
    
    write(*,*) 'Model :'  
    
    if(model_t==1)then 
    	allocate(FFHHHH(NL))
       	allocate(FFLLLL(NL))
    	
    	allocate(BBPP(NB,NTOTAL,NL))
    	allocate(QQ(NB,NTOTAL,NL))
    	allocate(EERR(NB,NTOTAL,NL))
    	allocate(PPRROOBB(NTOTAL,NTOTAL,NL))
    	
    	open (1,file='inputData/VE_Fhh.dat')
	    read (1,*) (FFHHHH(l), l=1,NL )
		open (1,file='inputData/VE_Fll.dat')
	    read (1,*) (FFLLLL(l), l=1,NL )	
	    
	 	write(*,*) 'Bayesian Learning'  	    
		write(*,*) ' '
	else
		allocate(FFLLLL(1))
		allocate(FFHHHH(1))

	    FFHHHH = 0.964
    	FFLLLL = 0.964
    	
    	allocate(BBPP(NB,NTOTAL,1))
    	allocate(QQ(NB,NTOTAL,1))
    	allocate(EERR(NB,NTOTAL,1))
    	allocate(PPRROOBB(NTOTAL,NTOTAL,1))
    
	   	write(*,*) 'Rational Expectactions'  	
	   	write(*,*) ' '    
    endif
    
    
	if(model_t==1)then 
		fileplace = "outData/model_bl/"
	else
		fileplace = "outData/model_re/"
    endif
    
	call cpu_time(start)
	do l=1,size(FFHHHH)
		
		write(*,*) l,' iterations of ', size(FFHHHH)
		
		Fhh = FFHHHH(l)
		Fll = FFLLLL(l)
		
		call mendoza(bp,q,er,Prob,Zprob,Z,BB,kappa,    &
		Fhh,Fll,gamma,beta,rho,sigmay, NS,RR,kappah,kappal,alpha,uptd,outfreq,iter_tol,tol,NB,bmin,bmax)
	
		BBPP(:,:,l) = bp
		QQ(:,:,l)	= q
		EERR(:,:,l) = er
		PPRROOBB(:,:,l) = Prob
		
	enddo 

    call cpu_time(finish)
    
	write(*,*) 'Exec. Time: ', finish-start , ' secs'
	call cpu_time(start)
	call print_binary_0(trim(fileplace) //'NB',real(NB))
	call print_binary_0(trim(fileplace) //'NK',real(NK))
	call print_binary_0(trim(fileplace) //'RR',RR)
	call print_binary_0(trim(fileplace) //'model_t',real(model_t))		
	call print_binary_1(trim(fileplace) //'Z',Z)
	call print_binary_1(trim(fileplace) //'BB',BB)
	call print_binary_2(trim(fileplace) //'Zprob',Zprob)
	call print_binary_2(trim(fileplace) //'Prob',Prob)
	call print_binary_2(trim(fileplace) //'kappa',kappa)
	call print_binary_3(trim(fileplace) //'BBPP',BBPP)
	call print_binary_3(trim(fileplace) //'QQ',QQ)
	call print_binary_3(trim(fileplace) //'EERR',EERR)
	call print_binary_3(trim(fileplace) //'PPRROOBB',PPRROOBB)
	
	OPEN (unit = 1, file=trim(fileplace)//'internal/QQ.txt',status='replace')
	write(1,*) QQ
	CLOSE (1)
	OPEN (unit = 1, file=trim(fileplace)//'internal/PPRROOBB.txt',status='replace')
	write(1,*) PPRROOBB
	CLOSE (1)
	OPEN (unit = 1, file=trim(fileplace)//'internal/EERR.txt',status='replace')
	write(1,*) EERR
	CLOSE (1)
	OPEN (unit = 1, file=trim(fileplace)//'internal/BBPP.txt',status='replace')
	write(1,*) BBPP
	CLOSE (1)
	
	call cpu_time(finish)
	write(*,*) 'Writ. Time: ', finish-start 
end program 
program mendoza_forecast
	
	use tools
	implicit none 
	
    integer ::  model_t = 0, i, s, ST
	integer :: T = 100000, cut, l, TT
	integer, parameter :: NB = 100, NK = 2, NTOTAL = 18, NS = 9, NL = 1 
	real, dimension(:,:), allocatable :: bpSIM, qSIM, erSIM, bindSIM, cSIM, CA_SIM, ySIM, KSIM, bplim_SIM, bpySIM, CAy_SIM, temp2
	
	real, dimension(NB,NTOTAL) :: bp, q, er
	real :: Z(NS), BB(NB), RR, temp, yt(NS), Zprob(NS,NS)
	real :: Prob(Ntotal,Ntotal)
	real, allocatable :: S_index(:), S_state(:,:)
	
	real :: bpyLr, cLR, qLR , var
	real, dimension(:), allocatable:: zbpySIM, zqSIM, zCAy_SIM, zcSIM, zerSIM 
	
    real, dimension(NB,NTOTAL,NL) :: BBPP,QQ,EERR, temp3
    real, dimension(NTOTAL,NTOTAL,NL) :: PPRROOBB
    real, dimension(NB,NTOTAL) ::kappa
	real, allocatable :: BS(:,:), chain(:)
    real :: start,finish

	character(100) :: fileplace
	call cpu_time(start)

	if(model_t==1)then 
		fileplace = "outData/model_bl/"
		write(*,*) ' Bayesian Learning Model'
	else
		fileplace = "outData/model_re/"
		write(*,*) ' Rational Expectation Model'
    endif

	OPEN (unit = 1, file = trim(fileplace)//'model_t.txt', status = 'old')
	READ (1,*) temp
	CLOSE (1)
	model_t=int(temp)
	OPEN (unit = 1, file = trim(fileplace)//'RR.txt', status = 'old')
	READ (1,*) RR
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'Z.txt', status = 'old')
	READ (1,*) Z
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'BB.txt', status = 'old')
	READ (1,*) BB
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'Prob.txt', status = 'old')
	READ (1,*) Prob
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'Zprob.txt', status = 'old')
	READ (1,*) Zprob
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'internal/BBPP.txt', status = 'old')
	READ (1,*) BBPP
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'internal/QQ.txt', status = 'old')
	READ (1,*) QQ
	CLOSE (1)
	OPEN (unit = 1, file = trim(fileplace)//'internal/EERR.txt', status = 'old')
	READ (1,*) EERR
	CLOSE (1)
	
	T = 44
	ST = 30000
	cut=0
	
	bpyLr 	= -0.40
	cLR		= 0.99
	qLR		= 0.50
	
	! Allocates Grids and initial values for simulation 
	
	allocate(bpSIM(ST,T+cut))	! Next period Debt
	allocate(bplim_SIM(ST,T+cut))
	allocate(bpySIM(ST,T+cut))
	allocate(qSIM(ST,T+cut)) 	! Land's price
	allocate(erSIM(ST,T+cut))	! Land's price
	allocate(bindSIM(ST,T+cut))	! Indicator of binding debt constraint
	allocate(cSIM(ST,T+cut))	! Consumption
	allocate(CA_SIM(ST,T+cut))
	allocate(CAy_SIM(ST,T+cut))
	allocate(ySIM(ST,T+cut))	! Output
 	allocate(KSIM(ST,T+cut))	

    allocate(S_index(T+cut))
    allocate(S_state(NTOTAL,T+cut))
    allocate(chain(T+cut))
    
    chain(1:36) = 0.93
    chain(37:44) = 0.63
    
    kSIM = kronecker(reshape(chain,(/1,size(chain,1)/)), ones(ST,1))
    bpSIM(:,1) = -0.363
    ySIM(:,1)  = 1.
    yt = Z
    
    do s = 1,ST
    
    	!call markov(S_index,S_state,Zprob,T+cut+1,1,real((/(i,i=1,size(Zprob,1))/)))
		S_index = markov_chain(Zprob,T+cut+1,1)
		
    	bp    = BBPP(:,:,1)
     	q     = QQ(:,:,1)
     	er    = EERR(:,:,1)
    	
		qSIM(s,1)   = piecewise(BB,q(:,int(S_index(1))), bpSIM(s,1))	
		erSIM(s,1)  = piecewise(BB,er(:,int(S_index(1))), bpSIM(s,1))
		
		do i=2,T+cut
			 ySIM(s,i) = yt(int(S_index(i))) 
			 			 
			 ! Retrieve Policy rules obtained in period t-1, these are changing
    		 ! policy rules according to households dynamic beliefs
     		if(model_t==1) then
     			bp = BBPP(:,:,i)
    			q  = QQ(:,:,i)
    			er = EERR(:,:,i)
    		else
    			bp = BBPP(:,:,1)
   				q  = QQ(:,:,1)
    			er = EERR(:,:,1)
    		endif
    		! Simulate control variables from date t until date T. At date 37 kappa
    		! moves from high securitization regime to low securitization regime
    		if (i<=36 .and. i>1)then 
			    bpSIM(s,i) = piecewise(BB,bp(:,int(S_index(i))),bpSIM(s,i-1))
  				qSIM(s,i)  = piecewise(BB,q(:,int(S_index(i))) ,bpSIM(s,i-1))
   				erSIM(s,i) = piecewise(BB,er(:,int(S_index(i))),bpSIM(s,i))
		    else    
			    bpSIM(s,i) = piecewise(BB,bp(:,size(Zprob,1)+int(S_index(i))),bpSIM(s,i-1))
  				qSIM(s,i)  = piecewise(BB,q(:,size(Zprob,1) +int(S_index(i))),bpSIM(s,i-1))
  				erSIM(s,i) = piecewise(BB,er(:,size(Zprob,1)+int(S_index(i))),bpSIM(s,i))
   		 	endif
    		! Simulate debt limit
  
   	 		bplim_SIM(s,i) = -(1+RR)*kSIM(s,i)*(qSIM(s,i))      
   	 		
   	 	enddo
	enddo

	! Compute consumption simulation
	cSIM  = bpSIM(:,1:T-1) + ySIM(:,1:T-1) - (1/(1+RR))*bpSIM(:,2:T)
 
	! bpSIM = bpSIM(:,1:T);
	CA_SIM  = (bpSIM(:,2:T)/(1+RR) - bpSIM(:,1:T-1))
 
	bpySIM = bpSIM/ySIM
	! Current Account = change in NFA = change in debt
	CAy_SIM = CA_SIM/ySIM(:,2:T)
 
	allocate(zbpySIM(T))
 	allocate(zqSIM(T))
 	allocate(zCAy_SIM(T))
 	allocate(zcSIM(T))
 	allocate(zerSIM(T))
 	
	zbpySIM   = mean(bpySIM)-bpyLR
	zqSIM     = mean(qSIM)/qLR-1
	zCAy_SIM  = mean(CAy_SIM)
	zcSIM     = mean(cSIM)/cLR - 1
	zerSIM    = mean(erSIM) - RR - 1
	
	call print_matrix_dat(trim(fileplace)//'zbpySIM.txt',zbpySIM,T,1) 	
	call print_matrix_dat(trim(fileplace)//'zqSIM.txt',zqSIM,T,1) 	
	call print_matrix_dat(trim(fileplace)//'zCAy_SIM.txt',zCAy_SIM,T,1) 	
	call print_matrix_dat(trim(fileplace)//'zcSIM.txt',zcSIM,T,1) 	
	call print_matrix_dat(trim(fileplace)//'zerSIM.txt',zerSIM,T,1) 		
	call cpu_time(finish)
	write(*,*) 'Exec. Time: ', finish-start , ' secs'
 end program 
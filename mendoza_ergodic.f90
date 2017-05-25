program mendoza_ergodic

	use tools
	implicit none 
	
	integer ::  model_t = 0, i
	integer :: T = 100000, cut, l, TT
	integer, parameter :: NB = 100, NK = 2, NTOTAL = 18, NS = 9, NL = 1
	real, dimension(:,:), allocatable :: bpSIM, qSIM, erSIM, bindSIM, cSIM, CA_SIM, ySIM, KSIM, bplim_SIM, bpySIM, CAy_SIM, temp2
	
	real, dimension(NB,NTOTAL) :: bp, q, er, yt
	real :: Z(NS), BB(NB), RR, temp
	real :: Prob(Ntotal,Ntotal)
	real, allocatable :: S_index(:), S_state(:,:)
	
    real, dimension(NB,NTOTAL,NL) :: BBPP,QQ,EERR, temp3
    real, dimension(NTOTAL,NTOTAL,NL) :: PPRROOBB
    real, dimension(NB,NTOTAL) ::kappa
    real :: start,finish
    
   	character(200) :: fileplace
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
	OPEN (unit = 1, file = trim(fileplace)//'kappa.txt', status = 'old')
	READ (1,*) kappa
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
	OPEN (unit = 1, file = trim(fileplace)//'internal/PPRROOBB.txt', status = 'old')
	READ (1,*) PPRROOBB	
	CLOSE (1)
			 
	cut = int(0.01*T)
	
	if(model_t==1) then
		TT = 44
	else
		TT = 1
	endif

	!*************************************************
	!1. Allocates Grids and initial values for simulation
	! Initialize vectors for simulation (decentralized equilibrium)
	
	allocate(bpSIM(TT,T+cut))	! Next period Debt
	allocate(bplim_SIM(TT,T+cut))
	allocate(bpySIM(TT,T+cut))

	allocate(qSIM(TT,T+cut)) 	! Land's price
	allocate(erSIM(TT,T+cut))	! Land's price
	allocate(bindSIM(TT,T+cut))	! Indicator of binding debt constraint
	allocate(cSIM(TT,T+cut))	! Consumption
	allocate(CA_SIM(TT,T+cut))
	allocate(CAy_SIM(TT,T+cut))
	allocate(ySIM(TT,T+cut))	! Output
 	allocate(KSIM(TT,T+cut))	

    allocate(S_index(T+cut))
    allocate(S_state(NTOTAL,T+cut))
	! Values at date 1
	bpSIM(:,1) = -0.363
 
	! Retrieving output values (this is an exogenous variable) 
 	yt = transpose(spread( (/Z,Z/) , 2, NB))    

 	do l=1,TT
    	! Uses the endowment and securitization joint transition matrix to generate a simulation of T periods.
     	bp    = BBPP(:,:,l)
     	q     = QQ(:,:,l)
     	er    = EERR(:,:,l)
     	Prob  = PPRROOBB(:,:,l)
     	
!  	 	call markov(S_index,S_state, Prob, T+cut+1 , 1 , real((/(i,i=1,NTOTAL)/)))
		S_index = markov_chain(Prob,T+cut+1,1)
     	ySIM(l,:) = yt(1,int(S_index))
     	kSIM(l,:) = kappa(1,int(S_index)) 

	 	do i=2,T+cut

			bpSIM(l,i) = piecewise(BB,bp(:,int(S_index(i))), bpSIM(l,i-1))	
			qSIM(l,i)  = piecewise(BB,q( :,int(S_index(i))), bpSIM(l,i-1))	
			erSIM(l,i) = piecewise(BB,er(:,int(S_index(i))), bpSIM(l,i))
			
			bplim_SIM(l,i) = -(1+RR)*kSIM(l,i)*qSIM(l,i) 
		end do
	end do 
	
	! Compute consumption simulation 
	cSIM = bpSIM(:,1:T+cut-1) + ySIM(:,1:T+cut-1) - (1/(1+RR)*BPSIM(:,2:T+cut))
	    
	CA_SIM = bpSIM(:,2:T+cut)/(1+RR) - bpSIM(:,1:T+cut-1)
	
	bpySIM = bpSIM/ySIM

	CAy_SIM = CA_SIM/ySIM(:,2:T+cut)

	bpySIM = bpySIM(:,cut+1:T+cut)
	qSIM = qSIM(:,cut+1:T+cut)
	erSIM = erSIM(:,cut+1:T+cut)
	cSIM = cSIM(:,cut+1:T+cut-1)
	CAy_SIM = CAy_SIM(:,cut+1:T+cut-1)
	
	call print_matrix_dat(trim(fileplace)//'bpySIM.txt',bpySIM,TT,T+cut)
	call cpu_time(finish)
	write(*,*) 'Exec. Time: ', finish-start , ' secs'

end program 
module mendoza_tools
use omp_lib
implicit none
contains
subroutine mendoza(bp,q,er,Prob,Zprob,Z,BB,kappa,    &
		Fhh,Fll,gamma,beta,rho,sigmay, NS,RR,kappah,kappal,alpha,uptd,outfreq,iter_tol,tol,NB,bmin,bmax)

	use basic
	use tools
		
	! local 
	integer 					:: i,j,k, outfreq, iter_tol, iter, pos
	integer, parameter			:: NK=2, NTOTAL=18
	real						:: Kprob(2,2), tol, kappas(NK), QQ, d2, val=0., inter(NTOTAL)=0., delta=1E-5, error_EE = 1.
	real    					:: cm, cl=0.001, ch=2., gridStepSize
	real, dimension(NB,NTOTAL)	:: c, EE, b, mup, emu, emuq, cbind, bpmax, oldc, oldq, oldbp
	real, dimension(NB,NTOTAL)	:: mupq, yq, yt, eyq, prueba(4,4), prueba2(2,2)=1.
	logical 					:: nowarning = .true. 
	
	! Dummy 
 	real,	 intent(in) :: gamma
	real, 	 intent(in) :: Fhh, Fll, beta, rho, sigmay, kappah, kappal, alpha, uptd, bmin, bmax
	integer, intent(in) :: NB,NS
	
!	integer, intent(out) :: 
	real, 	 intent(out) :: bp(NB,NTOTAL), q(NB,NTOTAL), er(NB,NTOTAL), Prob(NTOTAL,NTOTAL)
	real, 	 intent(out) :: Zprob(NS,NS), Z(NS), BB(NB), kappa(NB,NTOTAL)
	
!	integer, intent(inout) ::
	real, 	 intent(inout) :: RR 
	
 	! 1. Construction of Exogenous Processes
 	!***** Securitization Process 

	Kprob(1,:) = (/Fhh,1-Fhh/)
	Kprob(2,:) = (/1-Fll,Fll/)
 	Kprob = Kprob
 	kappas = (/kappah,kappal/)

	!***** Tauchen Hussey  (calls tauchenhussey subroutine) *****
	call tauchenhussey(Z, Zprob, NS, cm, rho, sigmay)
	Z = exp(Z)

	!***** Joint Transition Matrix ******	
	Prob = (kronecker(ones(NK,NK),Zprob)*kronecker(Kprob,ones(NS,NS)))
	!****** Total Number of Exogenous States (s) ******

	!****** Create debt grid ******
	BB(1) = bmin		!(finer grid for low i to improve accuracy near lower bound)
	do i=2,NB
    	BB(i)=BB(i-1)+(bmax-BB(i-1))/((NB-i+1.)**1.05)
	end do
	
	!****** Securitization degree at each state (b,z) *****
	do i =1,NB
		kappa(i,:) = (/(spread(kappas(1),1,NS)),spread(kappas(2),1,NS)/)
	end do
		
	!****** Endowment at each state (b,z) ******
	yt = transpose(spread( (/Z,Z/) , 2, NB))    

	!****** Debt at each state (b,z) ******
	b = spread(BB,2,NTOTAL)  
	
 	!**********************  Bond price at each state (b,z) *************************
	QQ=1/(1+RR);                   
 	  
	! 2. Decentralized Equilibrium
	! 2.1 Initial values
 
	! Recall we want policy functions for each triplet (b,z,kappa) where b: current
	! debt, is the endogenous state and z and kappa are the exogenous states.
	! Hence our policies are matrices of dimension NB x Ntotal
 
	!****** Initial values for equilibrium objects *****
	
	q = 0.1
	bp = b
	c = 1.
	EE = 0.
	
	cbind = 1.
	bpmax = - kappa*(1+RR)*q
	mup = 1.
	emu = 0.
	emuq = 0.
	
	! 2.2 Time Iteration Loop
	
	iter = 0 
	d2  = 100.
	
	write (*,*) 'DE iter 	Norm'
	
	do while (d2 > tol .and. iter_tol > iter )
	
		oldq = q
		oldc = c
		oldbp = bp
		
		mup = real(c**(-gamma))
		mupq = real((q+yt*alpha)*c**(-gamma))
		yq = real (q+yt*alpha)
		
		emu = 0.
		emuq = 0.
		eyq = 0. 
		do i = 1,NB
		do j = 1,NTOTAL		
		
			!******** Duda: los b(i,j) quedan por fuera de la grilla?	
			if (bp(i,j)>BB(NB) .or. bp(i,j)<BB(1)) then 
				!write (*,*) 'Se sale del rango bp en ' ,i , j , ' con valor ', bp(i,j)  
			end if 
			
			!********************************************************
			
			inter=0.
			!***** Piecewise interpolation by Johannes Brumm *****
			do k = 1,NTOTAL
				inter(k) =  piecewise(BB,mup(:,k),bp(i,j))
			end do	
			emu(i,j) = beta*(1+RR)*sum(inter*Prob(j,:))
			inter=0.	
			!*****************************************************
			!***** Piecewise interpolation by Johannes Brumm *****
			do k = 1,NTOTAL				
				inter(k) =  piecewise(BB,mupq(:,k),bp(i,j))	
			end do
			emuq(i,j) = beta*sum(inter*Prob(j,:))
			inter=0.	
			!*****************************************************
			!***** Piecewise interpolation by Johannes Brumm *****
			do k = 1,NTOTAL
				inter(k) =  piecewise(BB,yq(:,k),bp(i,j))			
			end do
			eyq(i,j) = sum(inter*Prob(j,:))
			!*****************************************************
		end do 
		end do 
		
		if (iter==0) then
			!write(*,*) 'emu', iter	
			!call print_matrix(emu,NB,NTOTAL)
			!write(*,*) 'q', iter	
			!call print_matrix(q,NB,NTOTAL)
			!stop
		end if
		
		!***** Find new constrained values *****
		bpmax = -(1/QQ)*kappa*q
		where (bpmax>bmax) bpmax=bmax
		where (bpmax<bmin) bpmax=bmin
		cbind = b+ yt -QQ*bpmax
		!***************************************
		
		do i = 1,NB
        do j = 1,Ntotal
            ! Calculate Euler Equation Error at maximum feasible consumption cbind(b,z)
            ! this is equation (4) in the handout.
            
            EE(i,j)=(cbind(i,j)**(-gamma))-emu(i,j)
         
            if (EE(i,j)>tol*1E-2)  then           ! If positive, constraint will be binding then:
              
                ! debt will be as big as possible               
                bp(i,j) = bpmax(i,j)     
                ! consumption is solved from budget constraint
                c(i,j)  = cbind(i,j)   
                ! land's price that clears the market
                q(i,j)  = emuq(i,j)/((cbind(i,j)**(-gamma))-kappa(i,j)*EE(i,j))
                ! Expected return of land
                er(i,j) =  eyq(i,j)/q(i,j)
                
            else							! Constraint not binding
               
                ! Define function that calculates the absolute value of Euler
                ! Equation Error for a given consumption
                
             	!********************* Bisection Algorithm ***********************
    			cl=0.01
				ch=2.        
				error_EE = 1. 
	             
				cm = (cl+ch)/2
				
				if ( ((1/cl)**gamma-emu(i,j)) * ((1/ch)**gamma-emu(i,j))>0  ) then 
					write (*,*) 'Warning: Bisection can not find a solution in EE'
					write (*,*) 'cl: ', ((1/cl)**gamma-emu(i,j)) , 'ch: ',  ((1/ch)**gamma-emu(i,j))
				end if
				do while (abs(error_EE)>delta) 
					cm = (cl+ch)/2
					if (  (cl**(-gamma)-emu(i,j))  *  (cm**(-gamma)-emu(i,j))<0  ) then 
						ch = cm
					else 
						cl = cm
					end if
					
					error_EE=cm**(-gamma)-emu(i,j)
				end do

				c(i,j)  = cm
                EE(i,j) = error_EE
               	!********************* Bisection Algorithm ***********************
				!*****************************************************************
				
                ! Solve Euler Equation, get suggested consumption as the argmin
                ! of Euler error                
                ! Solve debt from budget constraint, check if it is within grid bounds
                bp(i,j)=max((1/QQ) *( yt(i,j)+b(i,j)-c(i,j) ) , bmin)
                bp(i,j)=min(bp(i,j),bmax)
                ! Here debt may not be consistent with consumption anymore
                ! but it will be the optimal choice given the grid.
                mup(i,j)  = c(i,j)**(-gamma)
                ! Land's price clears the market
                q(i,j) = emuq(i,j)/(mup(i,j));
                ! Expected return of land
                er(i,j) =  eyq(i,j)/q(i,j);
            end if 
        end do
    	end do    	
    	
    	!===============Check collateral constraint==============================
    	! To make sure c is feasible:
    	c = b + yt - max(QQ*bp , (-kappa)*q)    ! bp/(1+R)>=-kappa*q
 		! bp=(1./Q).*(b+yt-c);
    	!========================================================================

    	iter = iter+1 							! Update iteration counter
    	! Calculate difference between new and old policies
   		
   		d2 = maxval((/maxval(maxval(abs(c-oldc),1),1),maxval(maxval(abs(bp-oldbp),1),1),maxval(maxval(abs(q-oldq),1),1)/))
		
    	! Print results once every (outfreq) iterations
    	if (mod(iter,outfreq) == 0) then
        	write (*,*) iter, '     ', d2
    	end if
   
    	!=====================Updating rules for next iteration==============
    	bp = uptd*bp + (1-uptd)*oldbp;
    	c  = uptd*c  + (1-uptd)*oldc;
    	q  = uptd*q  + (1-uptd)*oldq;
    	!====================================================================
   	end do	
	write (*,*) ' Done! '
    write (*,*) iter, '     ', d2
	
end subroutine mendoza
end module mendoza_tools
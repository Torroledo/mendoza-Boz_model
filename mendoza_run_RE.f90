!**************************************
!	
!	Mendoza model
!
!**************************************

program mendoza_run
	
	use mendoza_tools
	use basic
	
	implicit none 
	 
	!**************************************
	! Initialize the parameters
	!**************************************
	
	integer :: model_t = 0 					! RE Model = 0 , BL Model = 1

	real :: beta = 0.91, rho = 0.878, sigmay = 0.00663
	real :: RR      = 0.027, kappah = 0.93, kappal = 0.63, alpha  = 0.026
	
	!**************** Technical parameters  *****************
	real :: uptd	= 0.1 		! Weight on new policy function in update to next iteration
								! change to 0.1 or smaller (w # if policy function fluctuate)
	integer :: outfreq  = 5 	! Display frequency (shows in screen each 20th iteration)
	integer :: iter_tol = 3000  ! Iteration tolerance (max iterations)
	real 	:: gamma  = 2.
	
	! Numerical tolerance (convergence criterion)
	real :: tol	= 1E-4 	! High value for session, increase for accuracy

	!*********************************************************
	! Grid assignment 
	integer, parameter :: NB   = 20,  NK = 2, NTOTAL = 18, NS     = 9! number of nodes
	real :: bmin = -0.7, start, finish
	real :: bmax = -0.1
	
	! Probabilities 
    real :: Fhh = 0.964
    real :: Fll = 0.964
    
    real :: kappa(NB,NTOTAL), B(NB,1), Z(NS,1), Zprob(NS,NS), prob(NTOTAL,NTOTAL)=0., er(NB,NTOTAL)
    real :: q(NB,NTOTAL), bp(NB,NTOTAL)
    
	call cpu_time(start)
	call mendoza(bp,q,er,Prob,Zprob,Z,B,kappa,    Fhh,Fll,gamma,beta,rho,sigmay,RR,kappah,kappal,alpha,uptd,bmin,bmax)
	call cpu_time(finish)

	write(*,*) 'Time: ', finish -start, 's'	
	call print_matrix_dat('policy_bp.dat',bp, NB, NTOTAL)
end program 
module tools

  use basic
  implicit none

contains
  !***********************************************************************
  ! Types
  !***********************************************************************

  !***********************************************************************
  ! TAUCHENHUSSEY
  !***********************************************************************
  subroutine tauchenhussey(z, prob, N, mu, rho, sigma)
    
    !Local
    integer :: i
    real :: z_m(N,1)
    
    !Dummy 
    integer, intent(in) ::  N
    real, intent(in) :: mu, rho, sigma
    
    real, intent(out) :: z(N), prob(N, N)
	
    prob = transpose(loaddata(1, 'inputData/prob_tauchenhussey.dat', N, N))

    z_m = loaddata(2, 'inputData/z_tauchenhussey.dat', N, 1)
    z = z_m(:,1)
    
  end subroutine tauchenhussey
  
  !***********************************************************************
  ! ROUWENHORST
  !***********************************************************************
  subroutine rouwenhorst(z, prob, N, mu, rho, sigma)
    
    !Local
    integer :: i
    real :: z_m(N,1)
    
    !Dummy 
    integer, intent(in) ::  N
    real, intent(in) :: mu, rho, sigma
    
    real, intent(out) :: z(N), prob(N, N)

    prob = loaddata(1, 'prob.dat', N, N)

    z_m = loaddata(2, 'z.dat', N, 1)
    z = z_m(:,1)
    
  end subroutine rouwenhorst
  
  !***********************************************************************
  ! GRIDMAKE
  !***********************************************************************
  !******** GRIDMAKE2
  function gridmake2(a, a_n, b, b_n) result(grid)

    !Local
    integer :: i, j

    !Dummy
    integer :: a_n, b_n 
    real :: a(a_n), b(b_n)

    real :: grid(a_n*b_n,2)
    
    do i = 1,a_n
       do j = 1,b_n
          grid(i + a_n*(j-1),1) = a(i)
          grid(i + a_n*(j-1),2) = b(j)
       end do
    end do
  end function gridmake2

  !******** GRIDMAKE3
  function gridmake3(a, a_n, b, b_n, c, c_n) result(grid)

    !Local
    integer :: i, j, k 

    !Dummy
    integer, intent(in) :: a_n, b_n, c_n 
    real, intent(in) :: a(a_n), b(b_n), c(c_n)

    real :: grid(a_n*b_n*c_n,3)
    do i = 1,a_n
       do j = 1,b_n
          do k = 1,c_n
             grid(i + a_n*(j-1) + a_n*b_n*(k-1),1) = a(i)
             grid(i + a_n*(j-1) + a_n*b_n*(k-1),2) = b(j)
             grid(i + a_n*(j-1) + a_n*b_n*(k-1),3) = c(k)
          end do
       end do
    end do
  end function gridmake3
  
  !******** GRIDMAKE4
  function gridmake4(a, a_n, b, b_n, c, c_n, d, d_n) result(grid)
    
    !Local
    integer :: i, j, k, l
    
    !Dummy
    integer :: a_n, b_n, c_n, d_n 
    real :: a(a_n), b(b_n), c(c_n), d(d_n)
    
    real :: grid(a_n*b_n*c_n*d_n,4)
    
    do i = 1,a_n
       do j = 1,b_n
          do k = 1,c_n
             do l = 1,d_n
                grid(i + a_n*(j-1) + a_n*b_n*(k-1) + a_n*b_n*c_n*(l-1),1) = a(i)
                grid(i + a_n*(j-1) + a_n*b_n*(k-1) + a_n*b_n*c_n*(l-1),2) = b(j)
                grid(i + a_n*(j-1) + a_n*b_n*(k-1) + a_n*b_n*c_n*(l-1),3) = c(k)
                grid(i + a_n*(j-1) + a_n*b_n*(k-1) + a_n*b_n*c_n*(l-1),4) = d(l)
             end do
          end do
       end do
    end do
  end function gridmake4
     
  !**********************************************************
  ! KRONECKER
  !**********************************************************
  function kronecker(b,c) result(a)
        
    !*******************************************
    ! Computes a = b.kron.c
    !*******************************************
    
    !Local 
    integer :: i, j, k, l
    
    !Dummy 
    integer :: b_dim(2), c_dim(2)
    real   :: b(:,:), c(:,:)
    
    real, allocatable :: a(:,:)
	b_dim = shape(b)
	c_dim = shape(c)
	
    allocate(a(b_dim(1)*c_dim(1),b_dim(2)*c_dim(2)))
    
    do i = 1,b_dim(1)
    	do j = 1,b_dim(2)
			do k = 1,c_dim(1)
				do l = 1,c_dim(2)
                   a(k + (i-1)*c_dim(1), l + (j-1)*c_dim(2)) = b(i,j)*c(k,l)
                end do
            end do
       	end do
    end do
    
  end function kronecker

  !**********************************************************
  ! ERGODIC DISTRIBUTION
  !**********************************************************
  function ergdist(P,n) result(d_result)

    !Local
    integer, parameter :: maxit = 50000
    integer :: i
    real :: tol = 2.22e-16, dnext(n,1), change
    
    !Dummy
    integer :: n
    real :: P(n,n)

    real :: d(n,1), d_result(n)

    do i = 1,n
       d(i,1) = (1.)
    end do
    
    do i = 1,maxit
       dnext = transpose(matmul(transpose(d),P))
       change = maxval(abs(dnext-d))
       if (change<tol)then
          exit
       end if
       d = dnext
    end do
    d_result = d(:,1)
    
  end function ergdist
  
  !**********************************************************
  ! PIECEWISE LINEAR INTERPOLATION WITH EXTRAPOL
  !**********************************************************
  function piecewise(xgrid,ygrid,xp) result(yp)
  
	integer :: xr
	real :: xgrid(:), ygrid(:),xp, val, yp
	
	if(size(xgrid)==size(ygrid))then
	if    (xp>xgrid(size(xgrid))) then
		xr = size(xgrid)
	elseif(xp<xgrid(1)) then 
		xr = 2
	else
		val = minval(pack(xgrid,(xp<xgrid)) )	! value 
		xr   = minloc(abs(xgrid-val),1)			! position 
	endif
	
	yp = ((xgrid(xr)-xp)*ygrid(xr-1) + (xp-xgrid(xr-1))*ygrid(xr))/(xgrid(xr)-xgrid(xr-1))
	else
		write(0,*) ' No same size in grids'
		stop
	endif
  end function piecewise	
  
  !**********************************************************
  ! MARKOV CHAIN
  !**********************************************************
  subroutine markov(chain,state, T,n,s0,V)
  	
  	! Local
  	integer :: i, c, r, k, dim(2)
  	real,dimension(:), allocatable :: X
  	integer,dimension(:), allocatable :: s
  	real, allocatable :: cum(:,:), ppi(:), ones(:,:)
  	 
  	! Dummy 
  	integer,intent(in)	:: n, s0
  	real, intent(in)	:: T(:,:), V(:)
  	
  	real, intent(out), allocatable	:: chain(:), state(:,:)
  	
  	
  	dim = shape(T)
  	r=dim(1)
  	c=dim(2)
	allocate(ones(r,c))
	ones =1.
	allocate(chain(n-1))
  	allocate(state(r,n-1))
  	
  	allocate(X(n-1))
  	
  	do i=1,n-1
!		X(i)=rand()
		call random_number(X(i))
    enddo
!	call random_number(X)
  	allocate(s(r))
  	s=0.
  	s(s0) = 1
  	allocate(cum(r,c))

  	cum = matmul(T,triu(ones))
 
  	allocate(ppi(n+1))	
  	
  	do k=1,size(X)
  		state(:,k) = s
  		ppi = (/0.,matmul(s,cum)/)
  		s=0.
  		 do i = 1,r
  		 	if((X(k)<=ppi(i+1)) .and. (X(k)>ppi(i)))then
  		 		s(i)=1.
  		 	endif
  		 enddo
  	enddo
  	chain=matmul(V,state) 
 		  		

  end subroutine
  
  !**********************************************************
  ! TRIANGULAR UPPER MATRIX
  !**********************************************************
  function triu(matrix)
  
  	integer :: dim(2),i,j
  	real, allocatable :: triu(:,:)
  	real :: matrix(:,:)
  
  	triu = matrix
  	dim=shape(matrix)
  	
  	do i=1,dim(1)
  		do j=1,dim(2)
  			if (i>j)then
  				triu(i,j)=0.
  			endif
  		enddo
  	enddo
  end function 
  
  !**********************************************************
  ! PRINT BINARY FILE
  !**********************************************************
  subroutine print_binary_0(filename,var)
    real, intent(in) :: var
    character (len=*) :: filename
	open(unit=1,file=filename//'.txt',status='replace', action='write')
	write(1,*) var
	close(1)
  end subroutine 
  subroutine print_binary_1(filename,var)
    real, intent(in) :: var(:)
    character (len=*) :: filename
	open(unit=1,file=filename//'.txt',status='replace', action='write')
	write(1,*) var
	close(1)
  end subroutine 
  subroutine print_binary_2(filename,var)
    real, intent(in) :: var(:,:)
    character (len=*) :: filename
	integer :: o
	open(unit=1,file=filename//'.txt',status='replace', action='write')
	write(1,*) var
	close(1)
  end subroutine
  
  subroutine print_binary_3(filename,var)
    integer :: dim(3), i
    character(10) :: num
    real, intent(in) :: var(:,:,:)
    character (len=*) :: filename
    
    dim = shape(var)
    if(dim(3)==1)then   
		open(unit=1,file=filename//'.txt',status='replace', action='write')
		write(1,*) var(:,:,1)
		close(1)
	else 
		do i = 1,dim(3)
			open(unit=1,file=filename//'_'//trim(str(i))//'.txt',status='replace', action='write')
			write(1,*) var(:,:,i)
			close(1)
		enddo
	endif
  end subroutine 
  
  !**********************************************************
  ! LOAD BINARY FILE
  !**********************************************************
  function load_binary_0(filename) result(var)
    real :: var
    character (len=*) :: filename
	OPEN (unit = 1, file = filename, status = 'old')
	READ (1,*) var
	CLOSE (1)
  end function
  !**********************************************************
  ! MEAN
  !**********************************************************
  function mean(var)
    integer :: i, j
    real :: var(:,:)
    real, allocatable :: mean(:)   	
   	allocate(mean(size(var,2)))
   	do j = 1,size(var,2)
  		mean(j) = sum(var(:,j))/size(var,1)
  	enddo 		
  end function 
  !**********************************************************
  ! EASY MARKOV CHAIN
  !**********************************************************
  function markov_chain(T,n,sO) result(chain)
    integer :: i, j, time, n, sO, chain(n)
    real :: T(:,:), x 
 	chain = 1
 	chain(1) = sO
 	
	do time = 2,n 
		! Find current state
		i = chain(time-1)		
		call random_number(x)
!		x = rand()
		do j = 1,size(T,2)
	    	if(sum(T(i,1:j)) < x) then
	    		chain(time) = chain(time) + 1
 		   	endif
 		enddo
 		if(chain(time)>size(T,2))then
			chain(time) = size(T,2)	
		endif
    enddo
  end function 
end module tools
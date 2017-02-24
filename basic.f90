module basic
  
  implicit none

contains

  !**********************************************************
  ! REPVEC
  !**********************************************************
  function repvec(vec, n_vec, n_row, n_col)
        
    !*****************************************
    ! mat is a n_mat square matrix
    !*****************************************
    
    integer :: i, j, k, l, n_vec, n_row, n_col
    real :: vec(n_vec), repvec(n_vec*n_row,n_col)

    repvec = 0.
    do i = 1,n_vec
    	do j = 1,n_row
        	do k = 1,n_col
        		repvec(j, i + (k-1)*n_col) = vec(i)
        	end do
    	end do
    end do
  end function repvec

  !**********************************************************
  ! REPMAT
  !**********************************************************
  function repmat(mat, n_mat, n_row, n_col)
        
    !*****************************************
    ! mat is a n_mat square matrix
    !*****************************************
    
    integer :: i, j, k, l, n_mat, n_row, n_col
    real :: mat(n_mat,n_mat), repmat(n_mat*n_row,n_mat*n_col)

    repmat = 0.
    do i = 1,n_mat
       do j = 1,n_mat
          do k = 1,n_row
             do l = 1,n_col
                repmat(i + (k-1)*n_mat,j + (l-1)*n_mat) = mat(i,j)
             end do
          end do
       end do
    end do
  end function repmat
  
  !**********************************************************
  ! EYE
  !**********************************************************  
  function eye(m) result(matrix)
        
    integer :: i, j, m
    real	:: matrix(m,m)
    
    do i = 1,m
       do j = 1,m
          if( i==j )then
             matrix(i,j) = 1.
          else
             matrix(i,j) = 0.
          end if
       end do
    end do
    
  end function eye

  !**********************************************************
  ! ONES
  !**********************************************************  
  function ones(m,n) result(matrix)
        
    integer :: i, j, m, n
    real    :: matrix(m,n)
    
    do i = 1,m
       do j = 1,n
    	  matrix(i,j) = 1.
       end do
    end do
    
  end function ones
  !**********************************************************
  ! LOAD DATA
  !**********************************************************  
  function loaddata(unit, namefile, rows, cols) result(dataa)

    !Local
    
    !Dummy
    integer, intent(in) :: unit, rows, cols
    character (len=*) :: namefile

    real, dimension(rows,cols) :: dataa

    open(unit, file=namefile, action='read')
    read(unit,*) dataa
    close(unit)
    
  end function loaddata

  !**********************************************************
  ! LINSPACE
  !**********************************************************  
  function linspace(init, endd, points) result(line)

    !Local
    integer :: i 
    
    !Dummy
    integer :: points
    real :: init, endd

    real :: line(points)
    
    do i = 1, points
       line(i) = init + (endd-init)*(i-1)/(points-1)
    end do
  end function linspace

  !**********************************************************
  ! RESHAPE1
  !**********************************************************  
  function reshape1(vini, size) result(v)

    !Local
    integer :: i

    !Dummy
    integer :: size
    real :: vini, v(size)
    
    do i = 1, size
       if (i==1)then
          v(i) = vini  
       else
          v(i) = 0.  
       end if
    end do
  end function reshape1
  
  !**********************************************************
  ! PRINT VECTOR
  !**********************************************************  
  subroutine print_vector(v,n) 
    !Local
    integer :: i
    
    !Dummy
    integer, intent(in) :: n
    real, intent(in) :: v(n)

    do i = 1,n
       write(*,*) v(i)
    end do
  end subroutine print_vector

  !**********************************************************
  ! PRINT MATRIX CONSOLA
  !**********************************************************  
  subroutine print_matrix(v,n,m) 
    !Local
    integer :: i
    
    !Dummy
    integer, intent(in) :: n, m
    real, intent(in) :: v(n,m)

    do i = 1,n
       write(*,*) v(i,:)
    end do
  end subroutine print_matrix

  !**********************************************************
  ! PRINT MATRIX DAT
  !**********************************************************  
  subroutine print_matrix_dat(namefile,v,n,m) 
    !Local
    integer :: i
    
    !Dummy
    integer, intent(in) :: n, m
    real, intent(in) :: v(n,m)
    character (len=*) :: namefile

    open(1, file=namefile, action='write',status='replace')

    do i = 1,n
       write(1,*) v(i,:)
    end do
    close(1)
  end subroutine print_matrix_dat
  
  subroutine print_scalar_dat(namefile,v) 
    !Local
    
    !Dummy
    real, intent(in) :: v
    character (len=*) :: namefile
    
    open(1, file=namefile, action='write',status='replace')
    write(1,*) v
    close(1)
  end subroutine print_scalar_dat
  
  function str(k)
	! "Convert an integer to string."
    integer, intent(in) :: k
    character(20) :: str
    write (str, *) k
    str = adjustl(str)
  end function


end module basic

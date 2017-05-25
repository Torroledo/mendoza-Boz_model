all:
	gfortran -fopenmp mendoza_run_BL.f90 mendoza_tools_BL.f90 basic.f90 tools.f90 def_types.f90
	./a.out #> run_model.txt
ergodic:
	gfortran mendoza_ergodic.f90 basic.f90 tools.f90 def_types.f90
	./a.out
forecast: 
	gfortran mendoza_forecast.f90 basic.f90 tools.f90 def_types.f90
	./a.out
clean: 
	 

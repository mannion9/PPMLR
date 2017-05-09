CF		   = gfortran
FLAGS      = -fdefault-real-8
# debugFLAGS = -O0 -g -Wall  -Wtabs -Wuninitialized -Wextra -Wconversion -fbacktrace -fbounds-check
objects = Exact.o io.o Constants.o Init.o

exe : $(objects)
	@$(CF) $(FLAGS) $(debugFLAGS) -o exe $(objects)
Constants.mod : Constants.o Constants.f90
	@$(CF) $(FLAGS) $(debugFLAGS) -c Constants.f90
Constants.o : Constants.f90
	@$(CF) $(FLAGS) $(debugFLAGS) -c Constants.f90
%.o: %.f90 Constants.mod
	@$(CF) $(FLAGS) $(debugFLAGS) -c $<
clean:
	@$(RM) exe *.o *~ *.mod exe
